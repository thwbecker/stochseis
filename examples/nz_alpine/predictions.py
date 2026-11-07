# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Predictions for the model comparison: what stochastic models that fit the constraints produce.

    python3 examples/nz_alpine/predictions.py out_prefix label1 dir1[:filter] [label2 dir2[:filter] ...]
        [--site north=1] [--step 100]

Pools any number of cases (directories of catalogs or design vectors of
them, e.g. all cases of the realization sets) and weights every
realization by its likelihood under the paleoseismic constraints: the
fraction of its synthetic records of the observed length that match the
observed RI and CV at all three sites at once (record_match.py). The
weighted ensemble is the set of purely stochastic models consistent with
the constraints, spread over model classes in proportion to how well
each class matches; its weighted 5 to 95 percent ranges are the
predictions to compare with physics-based simulators.

Metrics per realization (long runs): RI, CV and SED at the three sites;
mean interval and CV of whole-fault M >= 7.5 events; maximum magnitude;
b below width saturation; excess of M >= 7.5 events over b = 1 from M5.5;
Bath gap (M >= 6 mainshocks); fraction of M >= 7.5 events rupturing all
three sites; fraction of central-site M >= 7.5 events also rupturing the
south site.

Outputs
  <out_prefix>.png  weighted distribution of each metric, stacked by case
  <out_prefix>.md   weighted median and 5-95 % range per metric, unweighted range for comparison,
                    and the weight carried by each case
"""
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
warnings.filterwarnings("ignore")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from control_summary import realization_metrics  # noqa: E402
from record_match import SITES, case_files, parse, sample_records, tolerances  # noqa: E402
from stochseis import Catalog  # noqa: E402
from stochseis.sites import load_sites  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data", "nz_alpine")
COLORS = ["C0", "C3", "C2", "C1", "C4", "C5", "C6", "C7", "C8", "C9"]
METRICS = [("south RI", "south RI [yr]"), ("central RI", "central RI [yr]"), ("north RI", "north RI [yr]"),
           ("south CV", "south CV"), ("central CV", "central CV"), ("north CV", "north CV"),
           ("central SED", "central SED [m]"), ("whole-fault RI", "whole-fault M>=7.5 mean interval [yr]"),
           ("whole-fault CV", "whole-fault M>=7.5 CV"), ("max M", "maximum magnitude"),
           ("b", "b below saturation"), ("excess M>=7.5", "M>=7.5 excess over b=1"), ("Bath", "Bath gap (M>=6)"),
           ("full-length fraction", "fraction of M>=7.5 rupturing S, C, N"),
           ("co-rupture C-S", "M>=7.5 at central also at south")]


def wquantile(x, w, q):
    """Weighted quantile; realizations with zero weight do not contribute."""
    keep = w > 0
    x, w = x[keep], w[keep]
    if len(x) == 0:
        return np.nan
    o = np.argsort(x)
    x, w = x[o], w[o]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()          # mid-point cumulative weights
    return float(np.interp(q, c, x))


def realization(path, sites, tol, opt):
    m = realization_metrics(path, sites, opt)
    r = __import__("nz_page").run_stats(path, sites, opt)
    for n in SITES:
        m[f"{n} RI"] = r["sites"][n]["ri"]
        m[f"{n} CV"] = r["sites"][n]["cv"]
    si = r["sys_intervals"]
    m["whole-fault RI"] = si.mean() if len(si) else np.nan
    m["max M"] = float(Catalog.load(path).events["mag"].max())
    R, _ = sample_records(path, sites, tol, opt)
    if len(R):
        ok = np.ones(len(R), bool)
        for j, n in enumerate(SITES):
            ok &= np.abs(R[:, j, 0] - tol[n]["ri"]) <= tol[n]["d_ri"]
            ok &= np.abs(R[:, j, 1] - tol[n]["cv"]) <= tol[n]["d_cv"]
        m["_weight"] = float(ok.mean())
    else:
        m["_weight"] = 0.0
    return m


def main():
    out_prefix, pairs, opt = parse(sys.argv)
    opt.setdefault("pattern_mag", 7.5)
    opt.setdefault("window", 3000.0)
    sites = [s for s in load_sites(os.path.join(DATA, "paleo_sites.csv")) if s["name"] in SITES]
    tol = tolerances({s["name"]: s for s in sites})
    rows = []
    for k, (label, spec) in enumerate(pairs):
        files = case_files(spec)
        if not files:
            sys.exit(f"no catalogs for {spec}")
        for p in files:
            m = realization(p, sites, tol, opt)
            m["_case"] = k
            rows.append(m)
        w = [r["_weight"] for r in rows if r["_case"] == k]
        print(f"{label}: {len(files)} realizations, mean match fraction {np.mean(w):.3f}")
    W = np.array([r["_weight"] for r in rows])
    case = np.array([r["_case"] for r in rows])
    if W.sum() == 0:
        sys.exit("no realization matches the constraints")
    ess = W.sum() ** 2 / np.sum(W ** 2)

    fig, axs = plt.subplots(5, 3, figsize=(16, 18))
    for ax, (key, lab) in zip(axs.ravel(), METRICS):
        x = np.array([r.get(key, np.nan) for r in rows], float)
        ok = np.isfinite(x)
        if ok.sum() < 2:
            ax.axis("off")
            continue
        bins = np.linspace(np.nanpercentile(x[ok], 1), np.nanpercentile(x[ok], 99), 25)
        data = [x[ok & (case == k)] for k in range(len(pairs))]
        wts = [W[ok & (case == k)] for k in range(len(pairs))]
        ax.hist(data, bins=bins, weights=wts, stacked=True, color=[COLORS[k % 10] for k in range(len(pairs))],
                label=[p[0] for p in pairs])
        lo, med, hi = (wquantile(x[ok], W[ok], q) for q in (0.05, 0.5, 0.95))
        for v, ls in ((lo, ":"), (med, "-"), (hi, ":")):
            ax.axvline(v, color="k", ls=ls, lw=1)
        ax.set(xlabel=lab, ylabel="weight", title=f"{med:.3g} [{lo:.3g}-{hi:.3g}]")
    axs[0, 0].legend(fontsize=7)
    fig.suptitle(f"constraint-weighted predictions ({len(rows)} realizations, effective sample size {ess:.0f}); "
                 "weight = fraction of records of the observed length matching RI and CV at all sites")
    fig.tight_layout()
    fig.savefig(out_prefix + ".png", dpi=100)
    plt.close(fig)

    with open(out_prefix + ".md", "w") as fh:
        fh.write(f"{len(rows)} realizations of {len(pairs)} cases; weight = fraction of records of the observed "
                 f"length matching RI and CV at all three sites (record_match.py); effective sample size {ess:.0f}.\n\n")
        fh.write("Weight carried by each case:\n\n| case | realizations | mean match fraction | share of total weight |\n"
                 "|---|---|---|---|\n")
        for k, (label, _) in enumerate(pairs):
            sel = case == k
            fh.write(f"| {label} | {sel.sum()} | {W[sel].mean():.3f} | {100 * W[sel].sum() / W.sum():.0f}% |\n")
        fh.write("\nPredictions (weighted median [5-95 %]) and, for comparison, the unweighted range over all "
                 "realizations:\n\n| metric | constraint-weighted | all realizations |\n|---|---|---|\n")
        for key, lab in METRICS:
            x = np.array([r.get(key, np.nan) for r in rows], float)
            ok = np.isfinite(x)
            if ok.sum() < 2:
                continue
            q = [wquantile(x[ok], W[ok], v) for v in (0.05, 0.5, 0.95)]
            u = np.percentile(x[ok], [5, 50, 95])
            fmt = "{:.0f}" if "RI" in key else "{:.2f}"
            fh.write(f"| {lab} | {fmt.format(q[1])} [{fmt.format(q[0])}-{fmt.format(q[2])}] | "
                     f"{fmt.format(u[1])} [{fmt.format(u[0])}-{fmt.format(u[2])}] |\n")
    print("wrote", out_prefix + ".png", out_prefix + ".md")


if __name__ == "__main__":
    main()
