"""Summary table of an Alpine campaign over all configurations.

    python3 examples/alpine/campaign_summary.py outdir name1 [name2 ...] > outdir/campaign_summary.md

Reads <outdir>/abc_<name>/sweep_summary.json for every configuration and
reports, at a 2 m detection threshold: the number of design vectors
accepted on the paleoseismic observables alone and with the seismicity
constraints, the best misfit and its site statistics, and median and
range of the main parameters over the accepted set (with seismicity).
Because all configurations draw the same vectors, it also reports how
many vectors are accepted in each pair of configurations (overlap).
"""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from abc_scoring import accepted, load_summary, score  # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)
KEYS = ["stress_drop_mpa", "nucleation_exponent", "branching_ratio", "arrest_scale", "admission_exponent",
        "b_value", "budget_scale", "slip_cv", "hazard_max_area_km2"]


def main():
    outdir = sys.argv[1]
    names = sys.argv[2:]
    res = {}
    for n in names:
        path = os.path.join(outdir, f"abc_{n}", "sweep_summary.json")
        if not os.path.exists(path):
            print(f"<!-- {n}: no summary at {path} -->")
            continue
        rows = load_summary(path)
        sp = score(rows, 2.0, False)
        ss = score(rows, 2.0, True)
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
    print("\n## Overlap of accepted vectors (paleo + seismicity): number accepted in both\n")
    print("| | " + " | ".join(res) + " |")
    print("|" + "---|" * (len(res) + 1))
    for a, ra in res.items():
        ka = {x["key"] for x in ra["acc_s"]}
        print(f"| {a} | " + " | ".join(str(len(ka & {x['key'] for x in rb['acc_s']})) for rb in res.values()) + " |")


if __name__ == "__main__":
    main()
