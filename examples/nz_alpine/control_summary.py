# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
"""What each irregularity control does: one figure for all realization sets.

    python3 examples/nz_alpine/control_summary.py out_prefix sets_dir [set names ...] [--site north=1]

sets_dir is the output directory of run_realization_sets.sh (it contains
specs/realizations_<set>_cases.txt and one catalog directory per set).
For every case of every set, each metric is computed per realization;
the figure shows the median and 5-95 percent range over realizations
against the cases of the set (one column per set, one row per metric).

Metrics: window CV at the three sites (records of the observed length,
median per realization), mean RI and SED at the central site, b of grown
ruptures below width saturation, excess of M >= 7.5 events over a b = 1
extrapolation from M >= 5.5, Bath gap (M >= 6 mainshocks), CV of the
intervals between whole-fault M >= 7.5 events, fraction of M >= 7.5
events that rupture all three sites, fraction of M >= 7.5 events at the
central site that also rupture the south site.
Writes <out_prefix>.png and <out_prefix>.md.
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

from nz_page import run_stats as page_stats  # noqa: E402
from spectrum_figures import case_files, parse  # noqa: E402
from stochseis import Catalog  # noqa: E402
from stochseis import analysis as an  # noqa: E402
from stochseis.sites import load_sites  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data", "nz_alpine")
METRICS = [("south window CV", "south CV"), ("central window CV", "central CV"), ("north window CV", "north CV"),
           ("central RI", "central RI [yr]"), ("central SED", "central SED [m]"), ("b", "b below saturation"),
           ("excess M>=7.5", "M>=7.5 excess over b=1"), ("Bath", "Bath gap (M>=6)"),
           ("whole-fault CV", "CV whole-fault M>=7.5"), ("full-length fraction", "fraction of M>=7.5 rupturing S, C, N"),
           ("co-rupture C-S", "M>=7.5 at central also at south")]


def realization_metrics(path, sites, opt):
    r = page_stats(path, sites, opt)
    m = {}
    for s in ("south", "central", "north"):
        w = np.array(r["windows"][s])
        m[f"{s} window CV"] = np.median(w[:, 1]) if len(w) else np.nan
    m["central RI"] = r["sites"]["central"]["ri"]
    m["central SED"] = r["sites"]["central"]["sed"]
    m["excess M>=7.5"] = r["excess75"]
    m["Bath"] = r["bath6"]
    si = r["sys_intervals"]
    m["whole-fault CV"] = si.std() / si.mean() if len(si) > 2 else np.nan
    tot = sum(r["patterns"].values())
    m["full-length fraction"] = r["patterns"]["SCN"] / tot if tot else np.nan
    sc = r["patterns"]["SC"] + r["patterns"]["SCN"]
    cen = sum(v for k, v in r["patterns"].items() if "C" in k)
    m["co-rupture C-S"] = sc / cen if cen else np.nan
    cat = Catalog.load(path)
    m["b"] = an.b_value_regimes(cat)["b_fault_below_sat"]
    return m


def main():
    argv = list(sys.argv)
    out_prefix, sets_dir = argv[1], argv[2]
    names = [a for a in argv[3:] if not a.startswith("--") and "=" not in a]
    _, _, opt = parse([argv[0], out_prefix] + [a for a in argv[3:] if a.startswith("--") or "=" in a])
    sites = load_sites(os.path.join(DATA, "paleo_sites.csv"))
    obs = {s["name"]: s for s in sites}
    if not names:
        names = sorted(f[len("realizations_"):-len("_cases.txt")] for f in os.listdir(os.path.join(sets_dir, "specs"))
                       if f.endswith("_cases.txt"))
    results = {}
    for name in names:
        cases = []
        with open(os.path.join(sets_dir, "specs", f"realizations_{name}_cases.txt")) as fh:
            for line in fh:
                if line.strip():
                    label, vec, _ = line.split()
                    files = case_files(f"{os.path.join(sets_dir, name)}:{vec}")
                    cases.append((label, [realization_metrics(p, sites, opt) for p in files]))
                    print(f"{name} {label}: {len(files)} realizations")
        results[name] = cases

    fig, axs = plt.subplots(len(METRICS), len(results), figsize=(4.5 * len(results), 1.9 * len(METRICS)),
                            squeeze=False, sharey="row")
    for j, (name, cases) in enumerate(results.items()):
        for i, (key, lab) in enumerate(METRICS):
            ax = axs[i, j]
            for k, (label, runs) in enumerate(cases):
                v = np.array([r[key] for r in runs], float)
                v = v[np.isfinite(v)]
                if len(v):
                    med = np.median(v)
                    ax.errorbar([k], [med], yerr=[[med - np.percentile(v, 5)], [np.percentile(v, 95) - med]],
                                fmt="o", color="C0", capsize=4)
            site = key.split()[0]
            if key.endswith("window CV") and site in obs:
                ax.axhline(obs[site]["ri_sigma_yr"] / obs[site]["ri_yr"], color="k", ls="--", lw=1)
            if key == "central RI":
                ax.axhline(obs["central"]["ri_yr"], color="k", ls="--", lw=1)
            if key == "central SED":
                ax.axhline(obs["central"]["sed_mean_m"], color="k", ls="--", lw=1)
            if key in ("b", "excess M>=7.5"):
                ax.axhline(1.0, color="k", ls=":", lw=1)
            if key == "Bath":
                ax.axhline(1.2, color="k", ls=":", lw=1)
            ax.set_xticks(range(len(cases)))
            ax.set_xticklabels([c[0] for c in cases] if i == len(METRICS) - 1 else [], fontsize=8)
            if j == 0:
                ax.set_ylabel(lab, fontsize=8)
            if i == 0:
                ax.set_title(name)
            ax.set_xlim(-0.5, len(cases) - 0.5)
    fig.suptitle("median and 5-95 % over realizations; dashed: observed, dotted: reference values")
    fig.tight_layout()
    fig.savefig(out_prefix + ".png", dpi=110)
    plt.close(fig)

    with open(out_prefix + ".md", "w") as fh:
        fh.write("Median [5-95 %] over realizations.\n\n")
        fh.write("| set | case | " + " | ".join(k for k, _ in METRICS) + " |\n")
        fh.write("|" + "---|" * (len(METRICS) + 2) + "\n")
        for name, cases in results.items():
            for label, runs in cases:
                cells = []
                for key, _ in METRICS:
                    v = np.array([r[key] for r in runs], float)
                    v = v[np.isfinite(v)]
                    fmt = "{:.0f}" if key == "central RI" else "{:.2f}"
                    cells.append((fmt.format(np.median(v)) + f" [{fmt.format(np.percentile(v, 5))}-"
                                  f"{fmt.format(np.percentile(v, 95))}]") if len(v) else "-")
                fh.write(f"| {name} | {label} | " + " | ".join(cells) + " |\n")
    print("wrote", out_prefix + ".png", out_prefix + ".md")


if __name__ == "__main__":
    main()
