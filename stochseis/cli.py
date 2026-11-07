"""Command line interface.

    python -m stochseis.cli run params.json out.npz [--set key=value ...]
    python -m stochseis.cli analyze out.npz --plots outdir
    python -m stochseis.cli sweep sweep.json outdir --nproc 4

A sweep file is JSON with a ``base`` parameter dictionary, a ``grid`` of
parameter lists (full cartesian product) and a ``seeds`` list.  Each
combination runs as an independent process; summary metrics are written
to ``outdir/sweep_summary.json`` and ``.csv``.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

from .analysis import summary_metrics
from .catalog import Catalog
from .params import Params
from .simulator import Simulation


def _parse_sets(items):
    out = {}
    for s in items or []:
        k, v = s.split("=", 1)
        try:
            v = json.loads(v)
        except json.JSONDecodeError:
            pass
        out[k] = v
    return out


def _json_safe(d):
    o = {}
    for k, v in d.items():
        if isinstance(v, (np.floating, np.integer)):
            v = v.item()
        if isinstance(v, float) and not np.isfinite(v):
            v = None
        o[k] = v
    return o


def cmd_run(args):
    p = Params.from_json(args.params) if args.params != "-" else Params()
    p = p.replace(**_parse_sets(args.set))
    sim = Simulation(p)
    cat = sim.run()
    cat.save(args.out)
    if args.summary:
        with open(args.summary, "w") as f:
            json.dump(_json_safe(summary_metrics(cat)), f, indent=2)
    print(f"wrote {args.out}: {len(cat)} events, closure {cat.meta['moment_closure']:.3f}")


def cmd_analyze(args):
    from . import plots

    cat = Catalog.load(args.catalog)
    m = summary_metrics(cat, mmain=args.mmain)
    os.makedirs(args.plots, exist_ok=True)
    stem = os.path.splitext(os.path.basename(args.catalog))[0]
    with open(os.path.join(args.plots, stem + "_summary.json"), "w") as f:
        json.dump(_json_safe(m), f, indent=2)
    plots.summary_figure(cat, os.path.join(args.plots, stem + "_summary.png"), mmain=args.mmain)
    plots.deficit_figure(cat, os.path.join(args.plots, stem + "_deficit.png"))
    plots.time_series_figure(cat, os.path.join(args.plots, stem + "_time_series.png"))
    plots.fmd_regimes_figure(cat, os.path.join(args.plots, stem + "_fmd_regimes.png"))
    keys = ["n_events", "n_fault_events", "max_mag", "b_sub", "b_fault_below_sat",
            "b_fault_above_sat", "slope_m0_L_small", "slope_m0_L_large", "omori_p", "omori_c",
            "bath_mean", "bath_n", "foreshock_fraction", "foreshock_chance",
            "cv_m6.5_declustered", "recurrence_m6.5_yr", "moment_closure", "final_mean_deficit",
            "realised_branching_ratio", "frac_stop_deficit_m6", "frac_surface_rupturing_m6.5"]
    for k in keys:
        if k in m:
            v = m[k]
            print(f"{k:28s} {v:.4g}" if isinstance(v, (float, np.floating)) else f"{k:28s} {v}")


def _sweep_one(job):
    """Run one sweep member; never raises (errors are recorded in the summary).

    The summary is also written to <outdir>/runs/<name>.json as soon as the
    run finishes, so an interrupted sweep can be resumed.
    """
    pdict, out_path, mmain, save = job
    pdict = dict(pdict)
    design_index = pdict.pop("design_index", None)
    t0 = time.time()
    m = {}
    try:
        p = Params.from_dict(pdict)
        cat = Simulation(p).run()
        if save:
            cat.save(out_path)
        try:
            m = summary_metrics(cat, mmain=mmain)
        except Exception as e:  # keep the run even if a metric fails
            m = dict(cat.meta)
            m["error"] = f"summary_metrics: {type(e).__name__}: {e}"
    except Exception as e:
        m = {"error": f"simulation: {type(e).__name__}: {e}"}
    m = _json_safe(m)
    m["file"] = out_path if save else ""
    if design_index is not None:
        m["design_index"] = design_index
    m["wall_seconds"] = time.time() - t0
    for k, v in pdict.items():
        m["p_" + k] = v
    rdir = os.path.join(os.path.dirname(out_path), "runs")
    os.makedirs(rdir, exist_ok=True)
    with open(os.path.join(rdir, os.path.basename(out_path).replace(".npz", ".json")), "w") as f:
        json.dump(m, f)
    return m


def cmd_sweep(args):
    with open(args.sweep) as f:
        spec = json.load(f)
    base = spec.get("base", {})
    grid = spec.get("grid", {})
    seeds = spec.get("seeds", [1])
    mmain = spec.get("mmain", 6.0)
    os.makedirs(args.outdir, exist_ok=True)
    keys = list(grid)
    jobs = []
    if spec.get("design") == "random":
        # random design: n_samples parameter vectors, uniform in "ranges" (log-uniform for keys in
        # "log_ranges"), uniform choice from "choices"; each vector runs for every seed
        rng = np.random.default_rng(spec.get("design_seed", 0))
        for i in range(int(spec["n_samples"])):
            d0 = dict(base)
            for k, (lo, hi) in spec.get("ranges", {}).items():
                d0[k] = float(rng.uniform(lo, hi))
            for k, (lo, hi) in spec.get("log_ranges", {}).items():
                d0[k] = float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
            for k, opts in spec.get("choices", {}).items():
                d0[k] = opts[int(rng.integers(len(opts)))]
            for seed in seeds:
                d = dict(d0)
                d["seed"] = seed
                d["design_index"] = i
                jobs.append((d, os.path.join(args.outdir, f"cat_r{i:04d}_s{seed}.npz"), mmain,
                             not args.no_catalogs))
    if spec.get("design") == "list":
        # explicit parameter vectors (e.g. an accepted ABC ensemble); each runs for every seed
        for i, vec in enumerate(spec["vectors"]):
            for seed in seeds:
                d = dict(base)
                d.update(vec)
                d["seed"] = seed
                d["design_index"] = i
                jobs.append((d, os.path.join(args.outdir, f"cat_v{i:04d}_s{seed}.npz"), mmain,
                             not args.no_catalogs))
    for combo in (itertools.product(*[grid[k] for k in keys]) if spec.get("design") not in ("random", "list") else []):
        for seed in seeds:
            d = dict(base)
            d.update(dict(zip(keys, combo)))
            d["seed"] = seed
            tag = "_".join(f"{k}{v}" for k, v in zip(keys, combo)) + f"_s{seed}"
            tag = tag.replace("/", "-")
            jobs.append((d, os.path.join(args.outdir, f"cat_{tag}.npz"), mmain, not args.no_catalogs))
    # resume: runs whose summary already exists are not repeated
    rdir = os.path.join(args.outdir, "runs")
    done, todo = [], []
    for j in jobs:
        fj = os.path.join(rdir, os.path.basename(j[1]).replace(".npz", ".json"))
        prev = None
        if os.path.exists(fj):
            with open(fj) as f:
                prev = json.load(f)
        if prev is not None and "error" not in prev:
            done.append(prev)
        else:                                 # not run yet, or failed: (re)run
            todo.append(j)
    print(f"{len(jobs)} runs ({len(done)} already done; failed runs are repeated) on {args.nproc} processes")
    t0 = time.time()
    if args.nproc > 1:
        with Pool(args.nproc) as pool:
            rows = done + list(pool.imap_unordered(_sweep_one, todo, chunksize=1))
    else:
        rows = done + [_sweep_one(j) for j in todo]
    n_err = sum(1 for r in rows if "error" in r)
    if n_err:
        print(f"warning: {n_err} runs recorded an error (key 'error' in the summary)")
    with open(os.path.join(args.outdir, "sweep_summary.json"), "w") as f:
        json.dump(rows, f, indent=1)
    cols = sorted({k for r in rows for k in r})
    with open(os.path.join(args.outdir, "sweep_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"done in {time.time() - t0:.0f} s; summary in {args.outdir}/sweep_summary.csv")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="stochseis")
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("run", help="run one simulation")
    a.add_argument("params", help="parameter JSON file, or - for defaults")
    a.add_argument("out", help="output catalog (.npz)")
    a.add_argument("--set", nargs="*", help="override parameters, key=value")
    a.add_argument("--summary", help="also write summary metrics JSON")
    a.set_defaults(func=cmd_run)
    a = sub.add_parser("analyze", help="summary metrics and figures for one catalog")
    a.add_argument("catalog")
    a.add_argument("--plots", default="plots")
    a.add_argument("--mmain", type=float, default=6.0)
    a.set_defaults(func=cmd_analyze)
    a = sub.add_parser("sweep", help="run a parameter sweep in parallel")
    a.add_argument("sweep", help="sweep specification JSON")
    a.add_argument("outdir")
    a.add_argument("--nproc", type=int, default=os.cpu_count() or 1)
    a.add_argument("--no-catalogs", action="store_true",
                   help="keep only the summary table (catalogs are 10 to 40 MB each)")
    a.set_defaults(func=cmd_sweep)
    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
