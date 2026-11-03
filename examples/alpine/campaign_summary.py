"""Summary table of an Alpine campaign over all configurations.

    python3 examples/alpine/campaign_summary.py outdir name1 [name2 ...] [--site north=1 ...] > outdir/campaign_summary.md

Reads <outdir>/abc_<name>/sweep_summary.json for every configuration and
reports, at a 2 m detection threshold: the number of design vectors
accepted on the paleoseismic observables alone and with the seismicity
constraints, the best misfit and its site statistics, and median and
range of the main parameters over the accepted set (with seismicity).
Because all configurations draw the same vectors, it also reports how
many vectors are accepted in each pair of configurations (overlap).

If <outdir>/posterior_<name>/sweep_summary.json exists, the posterior
vectors are rescored with their posterior runs (more seeds, longer
catalogs). Selecting the best of many noisy ABC estimates favours vectors
whose short-run statistics were lucky, so the rescored misfits and site
statistics are the ones to quote.
"""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from abc_scoring import accepted, load_summary, parse_site_thresholds, score  # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)
KEYS = ["stress_drop_mpa", "nucleation_exponent", "branching_ratio", "arrest_scale", "admission_exponent",
        "b_value", "budget_scale", "slip_cv", "hazard_max_area_km2"]


def main():
    argv = list(sys.argv)
    site_th = {}
    while "--site" in argv:
        k = argv.index("--site")
        site_th.update(parse_site_thresholds([argv[k + 1]]))
        argv = argv[:k] + argv[k + 2:]
    outdir = argv[1]
    names = argv[2:]
    if site_th:
        print("Site-specific detection thresholds: " + ", ".join(f"{k} {v:g} m" for k, v in site_th.items()) + "\n")
    res = {}
    for n in names:
        path = os.path.join(outdir, f"abc_{n}", "sweep_summary.json")
        if not os.path.exists(path):
            print(f"<!-- {n}: no summary at {path} -->")
            continue
        rows = load_summary(path)
        sp = score(rows, 2.0, False, site_thresholds=site_th)
        ss = score(rows, 2.0, True, site_thresholds=site_th)
        res[n] = dict(n=len(sp), n_err=sum(1 for r in rows if r.get("error")), acc_p=accepted(sp), acc_s=accepted(ss),
                      best=ss[0] if ss else None)
    if not res:
        sys.exit("no summaries found")
    print("## ABC results by configuration (2 m detection threshold)\n")
    print("| configuration | vectors | failed runs | accepted (paleo) | accepted (paleo + seismicity) | best misfit | "
          "best: south RI / CV / SED | central | north | b | Bath |")
    print("|" + "---|" * 11)
    for n, r in res.items():
        b = r["best"]
        d = b["detail"]
        print(f"| {n} | {r['n']} | {r['n_err']} | {len(r['acc_p'])} | {len(r['acc_s'])} | {b['misfit']:.2f} | "
              + " | ".join(f"{d[s][0]:.0f} / {d[s][1]:.2f} / {d[s][2]:.1f}" for s in ("south", "central", "north"))
              + f" | {d['seismicity'][0]:.2f} | {d['seismicity'][1]:.2f} |")
    print("\n## Accepted parameter ranges (paleo + seismicity): median [min-max]\n")
    print("| parameter | " + " | ".join(res) + " |")
    print("|" + "---|" * (len(res) + 1))
    for k in KEYS:
        cells = []
        for n, r in res.items():
            v = np.array([float(x["params"][k]) for x in r["acc_s"] if x["params"].get(k) is not None])
            cells.append(f"{np.median(v):.3g} [{v.min():.3g}-{v.max():.3g}]" if len(v) else "-")
        print(f"| {k} | " + " | ".join(cells) + " |")
    post = {}
    for n in res:
        path = os.path.join(outdir, f"posterior_{n}", "sweep_summary.json")
        if os.path.exists(path):
            post[n] = score(load_summary(path), 2.0, True, site_thresholds=site_th)
    if post:
        print("\n## Posterior vectors rescored with their posterior runs (paleo + seismicity, 2 m)\n")
        print("ABC misfit: from the ABC runs; posterior misfit: from the posterior runs of the same vector.\n")
        print("| configuration | vectors | still accepted | best posterior misfit | its ABC rank | "
              "best: south RI / CV / SED | central | north | b | Bath | median site CV (south, central, north) |")
        print("|" + "---|" * 11)
        for n, ps in post.items():
            b = ps[0]
            d = b["detail"]
            med = [np.nanmedian([x["detail"][s][1] for x in ps]) for s in ("south", "central", "north")]
            print(f"| {n} | {len(ps)} | {len(accepted(ps))} | {b['misfit']:.2f} | {int(b['key']) + 1} | "
                  + " | ".join(f"{d[s][0]:.0f} / {d[s][1]:.2f} / {d[s][2]:.1f}" for s in ("south", "central", "north"))
                  + f" | {d['seismicity'][0]:.2f} | {d['seismicity'][1]:.2f} | "
                  + ", ".join(f"{m:.2f}" for m in med) + " |")
    print("\n## Overlap of accepted vectors (paleo + seismicity): number accepted in both\n")
    print("| | " + " | ".join(res) + " |")
    print("|" + "---|" * (len(res) + 1))
    for a, ra in res.items():
        ka = {x["key"] for x in ra["acc_s"]}
        print(f"| {a} | " + " | ".join(str(len(ka & {x['key'] for x in rb['acc_s']})) for rb in res.values()) + " |")


if __name__ == "__main__":
    main()
