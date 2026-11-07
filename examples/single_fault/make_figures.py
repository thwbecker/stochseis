# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Figures for the single-fault parameter overview.

    python3 examples/single_fault/make_figures.py spec.json sweep_dir out_dir

spec.json is the overview spec (specs/parameter_overview.json: design
"list" with case_names and families); sweep_dir the output of the sweep
(catalogs cat_vNNNN_sS.npz and sweep_summary.json). Writes to out_dir:

  parameter_effects.png   one column per family, one row per statistic: median and 5-95 %
                          over realizations against the parameter value
  family_<n>.png          per family: frequency-magnitude distribution, stacked Omori decay after
                          M >= 6 mainshocks, Bath gaps, intervals between M >= 7 events
                          (lines: median over realizations; bands: 5-95 %)
  example_ruptures_<case>.png  slip of the six largest events of the first realization
  diagnostics_summary.png, diagnostics_time_series.png, diagnostics_fmd_regimes.png
                          standard diagnostics of the first realization of the first case
                          (stochseis.plots: summary page, sequence through time, frequency-magnitude
                          distribution by rupture regime)
  table.md                the statistics per case, median [5-95 %]

Statistics: b of grown ruptures below width saturation; excess of M >= 7
events over the b extrapolation from M >= 6; maximum magnitude; Omori p;
mean Bath gap; foreshock fraction (10 days, 10 km); CV of recurrence at
points along the fault (slip >= 1 m); CV of declustered M >= 7 recurrence
on the whole fault; moment closure after burn-in.
"""
import glob
import json
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
warnings.filterwarnings("ignore")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from stochseis import Catalog  # noqa: E402
from stochseis import analysis as an  # noqa: E402
from stochseis.plots import fmd_regimes_figure, largest_events_figure, summary_figure, time_series_figure  # noqa: E402

STATS = [("b_fault_below_sat", "b below saturation"), ("excess_m7", "M>=7 excess over b"),
         ("max_mag", "maximum magnitude"), ("omori_p", "Omori p"), ("bath_mean", "Bath gap"),
         ("foreshock_fraction", "foreshock fraction"), ("site_cv_slip1.0", "CV at points (>=1 m)"),
         ("cv_m7_declustered", "CV whole-fault M>=7"), ("closure_after_burnin", "moment closure")]
MAG = np.arange(4.0, 8.01, 0.05)
TOM = np.logspace(-3, np.log10(365.25), 26)
TOMC = np.sqrt(TOM[:-1] * TOM[1:])
COLORS = ["C0", "C3", "C2"]


def curves(path):
    cat = Catalog.load(path)
    ev = cat.events
    t0 = 0.1 * cat.duration
    m = ev["mag"][ev["t"] >= t0]
    years = cat.duration - t0
    out = {"fmd": np.array([np.sum(m >= x) for x in MAG]) / years}
    delays, fit = an.omori_stack(cat, mmain=6.0, t_start=t0)
    cnt, _ = np.histogram(delays, TOM)
    out["omori"] = cnt / np.diff(TOM) / max(fit["n_mainshocks"], 1)
    out["bath"] = an.bath_gaps(cat, mmain=6.0, t_start=t0)[0]
    tt = np.sort(ev["t"][(ev["mag"] >= 7.0) & (ev["t"] >= t0)])
    out["iv7"] = np.diff(tt)
    return out


def band(ax, x, ys, c, label, log=False):
    ys = np.array(ys, float)
    lo, med, hi = (np.nanpercentile(ys, q, axis=0) for q in (5, 50, 95))
    if log:
        lo, med, hi = (np.where(v > 0, v, np.nan) for v in (lo, med, hi))
    ax.fill_between(x, lo, hi, color=c, alpha=0.2, lw=0)
    ax.plot(x, med, color=c, label=label)


def main():
    spec_file, sweep_dir, out_dir = sys.argv[1:4]
    os.makedirs(out_dir, exist_ok=True)
    spec = json.load(open(spec_file))
    names = spec["case_names"]
    rows = json.load(open(os.path.join(sweep_dir, "sweep_summary.json")))
    stats = {i: [] for i in range(len(names))}
    for r in rows:
        if r.get("error"):
            continue
        i = int(r["design_index"])
        n6, n7, b = r.get("n_m6"), r.get("n_m7"), r.get("b_fault_below_sat")
        r["excess_m7"] = (n7 / n6) / 10 ** (-b) if n6 and n7 is not None and b else np.nan
        stats[i].append(r)
    cur = {}
    for i in range(len(names)):
        files = sorted(glob.glob(os.path.join(sweep_dir, f"cat_v{i:04d}_s*.npz")))
        cur[i] = [curves(p) for p in files]
        print(f"{names[i]}: {len(stats[i])} summaries, {len(files)} catalogs")

    def val(i, key):
        v = np.array([float(r.get(key)) if r.get(key) is not None else np.nan for r in stats[i]])
        return v[np.isfinite(v)]

    fams = list(spec["families"].items())
    fig, axs = plt.subplots(len(STATS), len(fams), figsize=(4.2 * len(fams), 1.7 * len(STATS)), squeeze=False,
                            sharey="row")
    for j, (fname, fam) in enumerate(fams):
        for k, (key, lab) in enumerate(STATS):
            ax = axs[k, j]
            for x, i in enumerate(fam["cases"]):
                v = val(i, key)
                if len(v):
                    med = np.median(v)
                    ax.errorbar([x], [med], yerr=[[med - np.percentile(v, 5)], [np.percentile(v, 95) - med]],
                                fmt="o", color="C0", capsize=4)
            ax.set_xticks(range(len(fam["cases"])))
            ax.set_xticklabels(fam["values"] if k == len(STATS) - 1 else [], fontsize=8)
            ax.set_xlim(-0.5, len(fam["cases"]) - 0.5)
            if j == 0:
                ax.set_ylabel(lab, fontsize=8)
            if k == 0:
                ax.set_title(fname, fontsize=10)
    fig.suptitle("single fault: median and 5-95 % over realizations")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "parameter_effects.png"), dpi=110)
    plt.close(fig)

    for n, (fname, fam) in enumerate(fams):
        fig, axs = plt.subplots(2, 2, figsize=(12, 9))
        for c, (i, v) in zip(COLORS, zip(fam["cases"], fam["values"])):
            cs = cur[i]
            if not cs:
                continue
            band(axs[0, 0], MAG, [x["fmd"] for x in cs], c, f"{fam['parameter']} = {v}", log=True)
            band(axs[0, 1], TOMC, [x["omori"] for x in cs], c, f"{v}", log=True)
            g = np.concatenate([x["bath"] for x in cs])
            if len(g):
                axs[1, 0].hist(g, bins=np.arange(0, 3.01, 0.1), density=True, histtype="step", color=c, lw=1.5,
                               label=f"{v}: mean {g.mean():.2f}")
            iv = np.concatenate([x["iv7"] for x in cs])
            if len(iv) > 2:
                axs[1, 1].hist(iv, bins=np.linspace(0, np.percentile(iv, 98), 30), density=True, histtype="step",
                               color=c, lw=1.5, label=f"{v}: mean {iv.mean():.0f} yr, CV {iv.std() / iv.mean():.2f}")
        ref = np.median([x["fmd"] for x in cur[fam["cases"][0]]], axis=0)
        k6 = np.argmin(np.abs(MAG - 5.0))
        axs[0, 0].semilogy(MAG, ref[k6] * 10 ** (-(MAG - 5.0)), "k:", lw=1, label="b = 1")
        axs[0, 0].set(yscale="log", xlabel="magnitude", ylabel="events >= M per yr", title="frequency-magnitude")
        axs[0, 1].set(xscale="log", yscale="log", xlabel="days after M >= 6 mainshock",
                      ylabel="aftershocks per day per mainshock", title="Omori decay")
        axs[1, 0].axvline(1.2, color="k", ls=":", lw=1)
        axs[1, 0].set(xlabel="M mainshock - M largest aftershock", ylabel="density", title="Bath gap (M >= 6)")
        axs[1, 1].set(xlabel="interval between M >= 7 events [yr]", ylabel="density",
                      title="recurrence of M >= 7 on the fault")
        for ax in axs.ravel():
            ax.legend(fontsize=7)
        fig.suptitle(f"{fname}: lines median, bands 5-95 % over realizations")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"family_{n + 1}.png"), dpi=110)
        plt.close(fig)

    for i in sorted({fams[0][1]["cases"][1]} | {fam["cases"][-1] for _, fam in fams}):
        files = sorted(glob.glob(os.path.join(sweep_dir, f"cat_v{i:04d}_s*.npz")))
        if files:
            tag = names[i].replace(" ", "_")
            largest_events_figure(Catalog.load(files[0]), n=6, path=os.path.join(out_dir, f"example_ruptures_{tag}.png"))

    first = sorted(glob.glob(os.path.join(sweep_dir, "cat_v0000_s*.npz")))
    if first:
        cat = Catalog.load(first[0])
        summary_figure(cat, os.path.join(out_dir, "diagnostics_summary.png"))
        time_series_figure(cat, os.path.join(out_dir, "diagnostics_time_series.png"))
        fmd_regimes_figure(cat, os.path.join(out_dir, "diagnostics_fmd_regimes.png"))

    with open(os.path.join(out_dir, "table.md"), "w") as fh:
        fh.write("Median [5-95 %] over realizations.\n\n| case | runs | " + " | ".join(l for _, l in STATS) + " |\n")
        fh.write("|" + "---|" * (len(STATS) + 2) + "\n")
        for i, nme in enumerate(names):
            cells = []
            for key, _ in STATS:
                v = val(i, key)
                cells.append(f"{np.median(v):.2f} [{np.percentile(v, 5):.2f}-{np.percentile(v, 95):.2f}]"
                             if len(v) else "-")
            fh.write(f"| {nme} | {len(stats[i])} | " + " | ".join(cells) + " |\n")
    print("wrote figures and table to", out_dir)


if __name__ == "__main__":
    main()
