# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Overview figures for posterior ensembles on one or more geometries.

    python3 examples/nz_alpine/overview_figures.py out_prefix label1 dir1 [label2 dir2 ...]

Each dir is the output of a posterior sweep (design "list", catalogs kept,
files cat_vNNNN_sS.npz). Vector 0000 is the best-fitting vector (the
posterior spec is sorted by misfit); all vectors form the ensemble.

Figure 1 (<out_prefix>_seismicity.png), one colour per configuration:
frequency-magnitude distribution, stacked Omori decay after M >= 6
mainshocks, distribution of the Bath gap, moment release relative to
loading over time (best run), distribution of the moment released in
100-yr windows (all runs), and distribution of the intervals between
system-wide M >= 7.5 events (all runs). Lines are the best vector (first
seed), shaded bands the 5 to 95 percent range over all runs.

Figure 2 (<out_prefix>_recurrence.png): mean recurrence interval versus
CV at each paleoseismic site (2 m detection threshold) for every run, with
the observations; SED versus RI; mean interval versus CV of system-wide
M >= 7.5 and M >= 8 events for every run; magnitude distribution of
M >= 6.5 events; distribution of rupture lengths for M >= 6.5;
magnitude versus rupture length and mean slip versus rupture length
(best run).
"""
import glob
import os
import re
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
from stochseis.sites import compare, load_sites  # noqa: E402

THRESHOLD = 2.0
MAG = np.arange(5.0, 8.65, 0.1)
TOM = np.logspace(-3, np.log10(365.25), 26)
TOMC = np.sqrt(TOM[:-1] * TOM[1:])
BATH = np.arange(-0.5, 3.01, 0.1)
LEN = np.logspace(np.log10(20), np.log10(800), 30)
WIN = 100.0                                    # window for moment release [yr]
WREL = np.linspace(0, 4, 41)                   # released / loaded moment per window
SYSI = np.arange(0, 1001, 25)                  # system-wide intervals [yr]
MBIG = np.arange(6.5, 8.61, 0.1)
COLORS = ["C0", "C3", "C2", "C1", "C4"]



def case_files(spec):
    """Catalogs of a case: "dir" (all), "dir:v0001" (one vector) or "dir:v0001_s5" (one realization)."""
    d, _, v = spec.partition(":")
    if "_s" in v:                                   # one realization: exact match (s1 must not match s10)
        return sorted(glob.glob(os.path.join(d, f"cat_{v}.npz")))
    pat = f"cat_{v}_s*.npz" if v else "cat_v*_s*.npz"
    return sorted(glob.glob(os.path.join(d, pat)))

def run_stats(path, sites):
    cat = Catalog.load(path)
    ev = cat.events
    T = cat.duration
    t0 = 0.1 * T
    sel = ev["t"] >= t0
    m = ev["mag"][sel]
    years = T - t0
    out = {"fmd": np.array([np.sum(m >= x) for x in MAG]) / years}
    delays, fit = an.omori_stack(cat, mmain=6.0, t_start=t0)
    cnt, _ = np.histogram(delays, TOM)
    out["omori"] = cnt / np.diff(TOM) / max(fit["n_mainshocks"], 1)
    gaps, _ = an.bath_gaps(cat, mmain=6.0, t_start=t0)
    out["bath"] = gaps
    t = np.linspace(0, T, 400)
    cum = np.interp(t, ev["t"], np.cumsum(ev["m0"]), left=0.0)
    out["moment"] = (t, (cum - cat.fault.moment_rate * t) / cat.fault.moment_rate)   # in years of loading
    edges = np.arange(t0, T + 1e-9, WIN)
    rel, _ = np.histogram(ev["t"], edges, weights=ev["m0"])
    out["window_release"] = rel / (cat.fault.moment_rate * WIN)
    fe = cat.fault_events()
    fe = fe[fe["t"] >= t0]
    big = fe["mag"] >= 6.5
    out["lengths"] = fe["length"][big] / 1e3
    out["mag_len"] = (fe["mag"][big], fe["length"][big] / 1e3)
    out["slip_len"] = (fe["dmean"][big], fe["length"][big] / 1e3)
    out["mags_big"] = fe["mag"][big]
    for mth in (7.5, 8.0):
        tt = np.sort(fe["t"][fe["mag"] >= mth])
        dt = np.diff(tt)
        out[f"sys_{mth:g}"] = dict(intervals=dt, ri=dt.mean() if len(dt) > 2 else np.nan,
                                   cv=dt.std(ddof=1) / dt.mean() if len(dt) > 2 else np.nan)
    out["sites"] = {r["name"]: r for r in compare(cat, sites, threshold=THRESHOLD, t_start=t0)}
    return out


def band(ax, x, ys, color, label, logy=False):
    ys = np.array(ys, float)
    lo = np.nanpercentile(ys, 5, axis=0)
    hi = np.nanpercentile(ys, 95, axis=0)
    if logy:
        lo = np.where(lo > 0, lo, np.nan)
        hi = np.where(hi > 0, hi, np.nan)
    ax.fill_between(x, lo, hi, color=color, alpha=0.2, lw=0)


def main():
    out_prefix = sys.argv[1]
    pairs = list(zip(sys.argv[2::2], sys.argv[3::2]))
    sites = load_sites(os.path.join(os.path.dirname(__file__), "..", "..", "data", "nz_alpine", "paleo_sites.csv"))
    data = {}
    for label, d in pairs:
        files = case_files(d)
        if not files:
            sys.exit(f"no catalogs in {d}")
        best = [f for f in files if re.search(r"cat_v0000_", f)]
        stats = [run_stats(f, sites) for f in files]
        data[label] = dict(stats=stats, best=stats[files.index(best[0])] if best else stats[0], n=len(files))
        print(f"{label}: {len(files)} catalogs")

    # ---------------- figure 1: seismicity --------------------------------
    fig, axs = plt.subplots(2, 3, figsize=(18, 9))
    for k, (label, D) in enumerate(data.items()):
        c = COLORS[k % len(COLORS)]
        st, b = D["stats"], D["best"]
        ax = axs[0, 0]
        band(ax, MAG, [s["fmd"] for s in st], c, label, logy=True)
        ax.semilogy(MAG, np.where(b["fmd"] > 0, b["fmd"], np.nan), color=c, label=f"{label} (best; band: {D['n']} runs)")
        ax = axs[0, 1]
        band(ax, TOMC, [s["omori"] for s in st], c, label, logy=True)
        ax.loglog(TOMC, np.where(b["omori"] > 0, b["omori"], np.nan), color=c, label=label)
        ax = axs[1, 0]
        allg = np.concatenate([s["bath"] for s in st if len(s["bath"])])
        ax.hist(allg, bins=BATH, density=True, histtype="step", color=c, lw=1.5,
                label=f"{label}: all runs, mean {np.mean(allg):.2f}")
        if len(b["bath"]):
            ax.hist(b["bath"], bins=BATH, density=True, histtype="stepfilled", color=c, alpha=0.25,
                    label=f"{label}: best, mean {np.mean(b['bath']):.2f}")
        ax = axs[1, 1]
        t, r = b["moment"]
        ax.plot(t, r, color=c, label=f"{label} (best, first seed)")
        ax = axs[0, 2]
        wr = np.concatenate([s["window_release"] for s in st])
        ax.hist(wr, bins=WREL, density=True, histtype="step", color=c, lw=1.5,
                label=f"{label}: all runs (CV {np.std(wr) / np.mean(wr):.2f})")
        ax = axs[1, 2]
        iv = np.concatenate([s["sys_7.5"]["intervals"] for s in st])
        if len(iv):
            ax.hist(iv, bins=SYSI, density=True, histtype="step", color=c, lw=1.5,
                    label=f"{label}: mean {iv.mean():.0f} yr, CV {iv.std() / iv.mean():.2f}")
    m0 = MAG[0]
    ref = data[next(iter(data))]["best"]["fmd"][0]
    axs[0, 0].semilogy(MAG, ref * 10 ** (-(MAG - m0)), "k:", lw=1, label="b = 1")
    axs[0, 0].set(xlabel="magnitude", ylabel="rate of events >= M [1/yr]", title="frequency-magnitude distribution")
    b0 = data[next(iter(data))]["best"]["omori"]
    k1 = np.argmin(np.abs(TOMC - 1.0))
    if b0[k1] > 0:
        axs[0, 1].loglog(TOMC, b0[k1] * TOMC[k1] / TOMC, "k:", lw=1, label="p = 1")
    axs[0, 1].set(xlabel="time after M >= 6 mainshock [days]", ylabel="aftershock rate per mainshock [1/day]",
                  title="Omori decay (window-selected)")
    axs[1, 0].axvline(1.2, color="k", ls=":", lw=1)
    axs[1, 0].set(xlabel="M mainshock - M largest aftershock", ylabel="density", title="Bath's law (M >= 6 mainshocks)")
    axs[1, 1].axhline(0, color="k", lw=0.5)
    axs[1, 1].set(xlabel="time [yr]", ylabel="released - loaded moment [years of loading]",
                  title="moment release relative to loading")
    axs[0, 2].axvline(1.0, color="k", ls=":", lw=1)
    axs[0, 2].set(xlabel=f"moment released in {WIN:.0f} yr / moment loaded", ylabel="density",
                  title=f"moment release in {WIN:.0f}-yr windows")
    axs[1, 2].set(xlabel="interval between system-wide M >= 7.5 events [yr]", ylabel="density",
                  title="whole-fault recurrence of M >= 7.5")
    for ax in axs.ravel():
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_prefix + "_seismicity.png", dpi=130)
    plt.close(fig)

    # ---------------- figure 2: recurrence and sizes --------------------
    fig, axs = plt.subplots(3, 3, figsize=(17, 14))
    for j, s in enumerate(sites):
        ax = axs[0, j]
        name = s["name"]
        ri_obs, sig = s["ri_yr"], s["ri_sigma_yr"]
        se = sig / np.sqrt(s["n_ri"] - 1)
        ax.errorbar([ri_obs], [sig / ri_obs], xerr=[se], yerr=[0.1], fmt="ks", ms=7, capsize=4, label="observed")
        for k, (label, D) in enumerate(data.items()):
            c = COLORS[k % len(COLORS)]
            ri = [x["sites"][name]["ri"] for x in D["stats"]]
            cv = [x["sites"][name]["cv"] for x in D["stats"]]
            ax.plot(ri, cv, "o", color=c, ms=3, alpha=0.5, label=f"{label} runs")
            ax.plot([D["best"]["sites"][name]["ri"]], [D["best"]["sites"][name]["cv"]], "*", color=c, ms=14,
                    mec="k", label=f"{label} best")
        ax.set(xlabel="mean recurrence interval [yr]", ylabel="CV of intervals",
               title=f"{name} site ({THRESHOLD:g} m threshold)")
        ax.legend(fontsize=7)
    ax = axs[1, 0]
    for j, s in enumerate(sites):
        ax.errorbar([s["ri_yr"]], [s["sed_mean_m"]], yerr=[[s["sed_mean_m"] - s["sed_min_m"]],
                                                            [s["sed_max_m"] - s["sed_mean_m"]]],
                    fmt="ks", ms=6, capsize=3)
        ax.annotate(s["name"], (s["ri_yr"], s["sed_mean_m"]), fontsize=7, xytext=(4, 4), textcoords="offset points")
    for k, (label, D) in enumerate(data.items()):
        c = COLORS[k % len(COLORS)]
        for s in sites:
            ax.plot([x["sites"][s["name"]]["ri"] for x in D["stats"]],
                    [x["sites"][s["name"]]["sed"] for x in D["stats"]], "o", color=c, ms=3, alpha=0.4)
        ax.plot([], [], "o", color=c, label=label)
    ax.set(xlabel="mean recurrence interval [yr]", ylabel="mean single-event displacement [m]",
           title="SED versus RI, all sites (observed: squares)")
    ax.legend(fontsize=7)
    ax = axs[1, 1]
    for k, (label, D) in enumerate(data.items()):
        c = COLORS[k % len(COLORS)]
        for mth, mk in ((7.5, "o"), (8.0, "^")):
            ax.plot([x[f"sys_{mth:g}"]["ri"] for x in D["stats"]], [x[f"sys_{mth:g}"]["cv"] for x in D["stats"]],
                    mk, color=c, ms=4, alpha=0.6, label=f"{label}, M >= {mth:g}")
            bb = D["best"][f"sys_{mth:g}"]
            ax.plot([bb["ri"]], [bb["cv"]], mk, color=c, ms=11, mec="k")
    ax.axhline(1.0, color="k", ls=":", lw=1)
    ax.set(xlabel="mean interval between system-wide events [yr]", ylabel="CV of intervals",
           title="whole-fault recurrence (large symbols: best)")
    ax.legend(fontsize=7)
    ax = axs[1, 2]
    for k, (label, D) in enumerate(data.items()):
        c = COLORS[k % len(COLORS)]
        allm = np.concatenate([x["mags_big"] for x in D["stats"]])
        ax.hist(allm, bins=MBIG, density=True, histtype="step", color=c, lw=1.5, label=f"{label} (all runs)")
        ax.hist(D["best"]["mags_big"], bins=MBIG, density=True, histtype="stepfilled", color=c, alpha=0.2,
                label=f"{label} (best)")
    ax.set(xlabel="magnitude", ylabel="density", title="magnitude distribution of M >= 6.5 events")
    ax.legend(fontsize=7)
    ax = axs[2, 0]
    for k, (label, D) in enumerate(data.items()):
        c = COLORS[k % len(COLORS)]
        ccdf = [np.array([np.mean(x["lengths"] >= L) if len(x["lengths"]) else np.nan for L in LEN]) for x in D["stats"]]
        band(ax, LEN, ccdf, c, label, logy=True)
        bl = D["best"]["lengths"]
        ax.loglog(LEN, [np.mean(bl >= L) if len(bl) else np.nan for L in LEN], color=c, label=f"{label} (best)")
    ax.set(xlabel="rupture length [km]", ylabel="fraction of M >= 6.5 events longer than L",
           title="rupture-length distribution (band: all runs)")
    ax.legend(fontsize=7)
    ax = axs[2, 1]
    for k, (label, D) in enumerate(data.items()):
        c = COLORS[k % len(COLORS)]
        mg, ln = D["best"]["mag_len"]
        ax.semilogx(ln, mg, "o", color=c, ms=3, alpha=0.5, label=f"{label} (best)")
    ax.set(xlabel="rupture length [km]", ylabel="magnitude", title="magnitude versus rupture length (M >= 6.5)")
    ax.legend(fontsize=7)
    ax = axs[2, 2]
    for k, (label, D) in enumerate(data.items()):
        c = COLORS[k % len(COLORS)]
        dm, ln = D["best"]["slip_len"]
        ax.loglog(ln, dm, "o", color=c, ms=3, alpha=0.5, label=f"{label} (best)")
    ax.set(xlabel="rupture length [km]", ylabel="mean slip [m]", title="mean slip versus rupture length (M >= 6.5)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_prefix + "_recurrence.png", dpi=130)
    plt.close(fig)
    print("wrote", out_prefix + "_seismicity.png", out_prefix + "_recurrence.png")


if __name__ == "__main__":
    main()
