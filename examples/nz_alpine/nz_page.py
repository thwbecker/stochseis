# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Page 2 of the Alpine figure set: site recurrence and rupture sizes.

    python3 examples/nz_alpine/nz_page.py out_prefix label1 dir1[:filter] [label2 dir2[:filter] ...]
        [--threshold 2] [--site north=1] [--pattern-mag 7.5]

A case is a directory of catalogs, one design vector of it (dir:v0001) or
one realization (dir:v0001_s5). One colour per case.

Record-length windows: each synthetic site record is cut into consecutive,
non-overlapping windows with as many intervals as the paleoseismic record
at that site (n_ri - 1: 12 south, 13 central, 3 north), and RI and CV are
computed per window. Their distributions show what records of the observed
length look like, and the fraction of windows at least as regular as the
observed record measures how unusual the observation is in each case.

Panels: window CV at the three sites (observed CV marked); window RI
versus window CV at the central site; site SED versus RI (run means)
with the observations; rate of M >= 7.5 events by the sites they rupture;
magnitude distribution of M >= 6.5 events; rupture-length distribution;
intervals between whole-fault M >= 7.5 events.
Writes <out_prefix>.png and <out_prefix>.md.
"""
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
warnings.filterwarnings("ignore")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from spectrum_figures import PATTERNS, case_files, parse  # noqa: E402
from spectrum_figures import run_stats as pattern_stats  # noqa: E402
from stochseis import Catalog  # noqa: E402
from stochseis.sites import load_sites  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "data", "nz_alpine")
COLORS = ["C0", "C3", "C2", "C1", "C4", "C5"]
LEN = np.logspace(np.log10(20), np.log10(800), 30)


def windows(intervals, n):
    """RI and CV of consecutive non-overlapping windows of n intervals."""
    out = []
    for k in range(0, len(intervals) - n + 1, n):
        w = intervals[k:k + n]
        if n >= 2:
            out.append((w.mean(), w.std(ddof=1) / w.mean()))
    return out


def run_stats(path, sites, opt):
    out = pattern_stats(path, sites, opt)
    cat = Catalog.load(path)
    fe = cat.fault_events()
    fe = fe[fe["t"] >= 0.1 * cat.duration]
    big = fe["mag"] >= 6.5
    out["lengths"] = fe["length"][big] / 1e3
    out["mags_big"] = fe["mag"][big]
    tt = np.sort(fe["t"][fe["mag"] >= 7.5])
    out["sys_intervals"] = np.diff(tt)
    out["windows"] = {}
    for s in sites:
        n = int(s["n_ri"]) - 1
        out["windows"][s["name"]] = windows(out["sites"][s["name"]]["intervals"], n)
    return out


def main():
    out_prefix, pairs, opt = parse(sys.argv)
    sites = load_sites(os.path.join(DATA, "paleo_sites.csv"))
    obs = {s["name"]: s for s in sites}
    cases = []
    for label, spec in pairs:
        files = case_files(spec)
        if not files:
            sys.exit(f"no catalogs for {spec}")
        cases.append((label, [run_stats(p, sites, opt) for p in files]))
        print(f"{label}: {len(files)} catalogs")

    fig, axs = plt.subplots(3, 3, figsize=(17, 14))
    rows = []
    for k, (label, runs) in enumerate(cases):
        c = COLORS[k % len(COLORS)]
        row = dict(label=label, runs=len(runs))
        for j, name in enumerate(("south", "central", "north")):
            w = np.array([x for r in runs for x in r["windows"][name]])
            ax = axs[0, j]
            cv_obs = obs[name]["ri_sigma_yr"] / obs[name]["ri_yr"]
            if len(w):
                frac = np.mean(w[:, 1] <= cv_obs)
                ax.hist(w[:, 1], bins=np.linspace(0, 1, 41), density=True, histtype="step", color=c, lw=1.5,
                        label=f"{label}: {len(w)} windows, {100 * frac:.0f}% <= obs")
                row[name] = (len(w), np.median(w[:, 0]), np.median(w[:, 1]), np.percentile(w[:, 1], 5),
                             np.percentile(w[:, 1], 95), frac)
            else:
                row[name] = (0, np.nan, np.nan, np.nan, np.nan, np.nan)
            if k == 0:
                ax.axvline(cv_obs, color="k", ls="--", label=f"observed {cv_obs:.2f}")
            n = int(obs[name]["n_ri"]) - 1
            ax.set(xlabel=f"CV in windows of {n} intervals", ylabel="density",
                   title=f"{name} site: records of the observed length")
            ax.legend(fontsize=7)
        w = np.array([x for r in runs for x in r["windows"]["central"]])
        if len(w):
            axs[1, 0].plot(w[:, 0], w[:, 1], "o", color=c, ms=3, alpha=0.4, label=label)
        for s in sites:
            axs[1, 1].plot([r["sites"][s["name"]]["ri"] for r in runs], [r["sites"][s["name"]]["sed"] for r in runs],
                           "o", color=c, ms=3, alpha=0.4)
        axs[1, 1].plot([], [], "o", color=c, label=label)
        P = np.array([[r["patterns"][p] for p in PATTERNS] for r in runs])
        width = 0.8 / len(cases)
        x = np.arange(len(PATTERNS)) + (k - (len(cases) - 1) / 2) * width
        axs[1, 2].bar(x, np.median(P, 0), width, color=c, edgecolor="k", label=label)
        row["patterns"] = np.median(P, 0)
        allm = np.concatenate([r["mags_big"] for r in runs])
        axs[2, 0].hist(allm, bins=np.arange(6.5, 8.61, 0.1), density=True, histtype="step", color=c, lw=1.5,
                       label=label)
        ccdf = np.array([[np.mean(r["lengths"] >= L) if len(r["lengths"]) else np.nan for L in LEN] for r in runs])
        axs[2, 1].loglog(LEN, np.nanmedian(ccdf, 0), color=c, label=label)
        if len(runs) > 1:
            axs[2, 1].fill_between(LEN, np.nanpercentile(ccdf, 5, 0), np.nanpercentile(ccdf, 95, 0), color=c,
                                   alpha=0.15, lw=0)
        si = np.concatenate([r["sys_intervals"] for r in runs])
        if len(si):
            axs[2, 2].hist(si, bins=np.arange(0, 1001, 25), density=True, histtype="step", color=c, lw=1.5,
                           label=f"{label}: mean {si.mean():.0f} yr, CV {si.std() / si.mean():.2f}")
        rows.append(row)
    s = obs["central"]
    axs[1, 0].errorbar([s["ri_yr"]], [s["ri_sigma_yr"] / s["ri_yr"]], fmt="ks", ms=8, label="observed")
    axs[1, 0].set(xlabel="window mean interval [yr]", ylabel="window CV",
                  title=f"central site: windows of {int(s['n_ri']) - 1} intervals")
    axs[1, 0].legend(fontsize=7)
    for s in sites:
        axs[1, 1].errorbar([s["ri_yr"]], [s["sed_mean_m"]], yerr=[[s["sed_mean_m"] - s["sed_min_m"]],
                                                                   [s["sed_max_m"] - s["sed_mean_m"]]],
                           fmt="ks", ms=6, capsize=3)
        axs[1, 1].annotate(s["name"], (s["ri_yr"], s["sed_mean_m"]), fontsize=7, xytext=(4, 4),
                           textcoords="offset points")
    axs[1, 1].set(xlabel="mean recurrence interval [yr]", ylabel="mean SED [m]",
                  title="SED versus RI, run means (observed: squares)")
    axs[1, 1].legend(fontsize=7)
    axs[1, 2].set_xticks(range(len(PATTERNS)))
    axs[1, 2].set_xticklabels(PATTERNS)
    axs[1, 2].set(ylabel=f"M >= {opt['pattern_mag']:g} events per 1000 yr (median)",
                  title="sites ruptured together (S south, C central, N north)")
    axs[1, 2].legend(fontsize=7)
    axs[2, 0].set(xlabel="Mw", ylabel="density", title="magnitudes of M >= 6.5 events")
    axs[2, 0].legend(fontsize=7)
    axs[2, 1].set(xlabel="rupture length [km]", ylabel="fraction of M >= 6.5 events longer than L",
                  title="rupture lengths (median and 5-95 % over runs)")
    axs[2, 1].legend(fontsize=7)
    axs[2, 2].set(xlabel="interval [yr]", ylabel="density", title="whole-fault recurrence of M >= 7.5")
    axs[2, 2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_prefix + ".png", dpi=110)
    plt.close(fig)

    with open(out_prefix + ".md", "w") as fh:
        th = ", ".join(f"{k} {v:g} m" for k, v in opt["site"].items())
        fh.write(f"Detection threshold {opt['threshold']:g} m{', ' + th if th else ''}. Windows: consecutive records "
                 "with the observed number of intervals. Columns per site: windows, median window RI, median window "
                 "CV [5-95 %], fraction of windows with CV at or below the observed CV.\n\n")
        fh.write("| case | runs | south | central | north | " + " | ".join(PATTERNS) + " |\n")
        fh.write("|" + "---|" * (5 + len(PATTERNS)) + "\n")
        for r in rows:
            cells = []
            for name in ("south", "central", "north"):
                nw, ri, cv, lo, hi, fr = r[name]
                cells.append(f"{nw}: {ri:.0f} yr, {cv:.2f} [{lo:.2f}-{hi:.2f}], {100 * fr:.0f}%" if nw else "-")
            fh.write(f"| {r['label']} | {r['runs']} | " + " | ".join(cells) + " | "
                     + " | ".join(f"{v:.1f}" for v in r["patterns"]) + " |\n")
        fh.write("| observed | | " + " | ".join(
            f"{obs[n]['ri_yr']:.0f} yr, {obs[n]['ri_sigma_yr'] / obs[n]['ri_yr']:.2f}" for n in ("south", "central", "north"))
            + " |" + " |" * len(PATTERNS) + "\n")
    print("wrote", out_prefix + ".png", out_prefix + ".md")


if __name__ == "__main__":
    main()
