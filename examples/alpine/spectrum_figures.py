"""Cases side by side: what kinds of earthquake sequences each case produces.

    python3 examples/alpine/spectrum_figures.py out_prefix label1 dir1[:vNNNN] [label2 dir2[:vNNNN] ...]
        [--threshold 2] [--site north=1] [--window 3000] [--pattern-mag 7.5]

A case is a sweep output directory with catalogs (cat_vNNNN_sS.npz), one
design vector of it (dir:v0001) or one realization (dir:v0001_s5). Use it for

* the realization sets (run_realization_sets.sh: cases that differ in one
  irregularity control), with one case per design vector, and
* the scenarios compatible with the constraints, with one case per
  posterior ensemble (all accepted vectors pooled).

Outputs
  <out_prefix>.png  one row per case: space-time plot of M >= 7 ruptures in a
                    window of the first run (bars: along-strike extent in
                    northing; dots: hypocentre, coloured by Mw); interval
                    distributions at the three sites against the paleoseismic
                    normal distributions; frequency-magnitude distribution
                    against b = 1; frequencies of rupture patterns of large
                    events (which sites rupture together)
  <out_prefix>.md   per case: site RI / CV / SED (median over runs), Bath gap
                    (M >= 6 and M >= 7 mainshocks), excess of M >= 7.5 events
                    over a b = 1 extrapolation from M >= 5.5, rupture-pattern
                    frequencies

Site records use the detection threshold (default 2 m; per site with
--site) and the merge interval the catalogs were run with (site_merge_yr).
Rupture patterns use the same per-site thresholds without merging.
"""
import glob
import os
import sys
import warnings
from math import erf

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
warnings.filterwarnings("ignore")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from stochseis import Catalog  # noqa: E402
from stochseis import analysis as an  # noqa: E402
from stochseis.sites import load_sites, site_patch, site_records  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "nz_alpine")
SITE_COLORS = {"south": "C0", "central": "C1", "north": "C2"}
PATTERNS = ["S", "C", "N", "SC", "CN", "SCN", "none"]
MAG = np.arange(5.0, 8.61, 0.05)


def parse(argv):
    opt = dict(threshold=2.0, site={}, window=3000.0, pattern_mag=7.5)
    rest = [argv[1]]
    k = 2
    while k < len(argv):
        a = argv[k]
        if a == "--threshold":
            opt["threshold"] = float(argv[k + 1])
            k += 2
        elif a == "--site":
            n, v = argv[k + 1].split("=")
            opt["site"][n] = float(v)
            k += 2
        elif a == "--window":
            opt["window"] = float(argv[k + 1])
            k += 2
        elif a == "--pattern-mag":
            opt["pattern_mag"] = float(argv[k + 1])
            k += 2
        else:
            rest.append(a)
            k += 1
    return rest[0], list(zip(rest[1::2], rest[2::2])), opt


def case_files(spec):
    d, _, v = spec.partition(":")
    if "_s" in v:                                   # one realization: exact match (s1 must not match s10)
        return sorted(glob.glob(os.path.join(d, f"cat_{v}.npz")))
    pat = f"cat_{v}_s*.npz" if v else "cat_v*_s*.npz"
    return sorted(glob.glob(os.path.join(d, pat)))


def run_stats(path, sites, opt):
    cat = Catalog.load(path)
    ev = cat.events
    f = cat.fault
    T = cat.duration
    t0 = 0.1 * T
    years = T - t0
    out = {}
    m = ev["mag"][ev["t"] >= t0]
    out["fmd"] = np.array([np.sum(m >= x) for x in MAG]) / years * 1e4
    n55 = np.sum(m >= 5.5)
    n75 = np.sum(m >= 7.5)
    out["excess75"] = (n75 / n55) / 10 ** (-2.0) if n55 else np.nan
    g6, _ = an.bath_gaps(cat, mmain=6.0, t_start=t0)
    g7, _ = an.bath_gaps(cat, mmain=7.0, t_start=t0)
    out["bath6"] = float(np.mean(g6)) if len(g6) else np.nan
    out["bath7"] = float(np.mean(g7)) if len(g7) else np.nan
    merge = getattr(cat.params, "site_merge_yr", 0.0)
    th = {s["name"]: opt["site"].get(s["name"], opt["threshold"]) for s in sites}
    out["sites"] = {}
    member = {}
    for s in sites:
        n = s["name"]
        rec = site_records(cat, [s], threshold=th[n], t_start=t0, merge_yr=merge)[n]
        t = rec["times"]
        dt = np.diff(t)
        out["sites"][n] = dict(intervals=dt, ri=dt.mean() if len(dt) > 1 else np.nan,
                               cv=dt.std(ddof=1) / dt.mean() if len(dt) > 1 else np.nan,
                               sed=rec["sed"].mean() if len(rec["sed"]) else np.nan)
        member[n] = set(site_records(cat, [s], threshold=th[n], t_start=t0, merge_yr=0.0)[n]["events"].tolist())
    big = np.where((ev["t"] >= t0) & (ev["mag"] >= opt["pattern_mag"]))[0]
    counts = dict.fromkeys(PATTERNS, 0)
    for e in big:
        key = "".join(c for c, n in (("S", "south"), ("C", "central"), ("N", "north")) if e in member[n]) or "none"
        if key not in counts:            # S and N without C (should not happen on one fault)
            key = "none"
        counts[key] += 1
    out["patterns"] = {k: v / years * 1000.0 for k, v in counts.items()}      # per 1000 yr
    out["n_big"] = len(big)
    # space-time window: last `window` years of the run
    st = []
    for k, e in enumerate(cat.slip_event):
        if ev["mag"][e] < 7.0 or ev["t"][e] < T - opt["window"]:
            continue
        a, b = cat.slip_ptr[k], cat.slip_ptr[k + 1]
        idx = cat.slip_idx[a:b]
        st.append((ev["t"][e], f.y[idx].min(), f.y[idx].max(), f.y[ev["nuc"][e]], ev["mag"][e]))
    out["spacetime"] = np.array(st) if st else np.zeros((0, 5))
    out["site_y"] = {s["name"]: f.y[site_patch(f, s["x_nztm_m"], s["y_nztm_m"])[0]] for s in sites}
    out["T"] = T
    return out


def main():
    out_prefix, pairs, opt = parse(sys.argv)
    sites = load_sites(os.path.join(DATA, "paleo_sites.csv"))
    cases = []
    for label, spec in pairs:
        files = case_files(spec)
        if not files:
            sys.exit(f"no catalogs for {spec}")
        cases.append((label, [run_stats(p, sites, opt) for p in files]))
        print(f"{label}: {len(files)} catalogs")

    nrow = len(cases)
    fig, axs = plt.subplots(nrow, 4, figsize=(22, 4.2 * nrow), squeeze=False,
                            gridspec_kw=dict(width_ratios=[1.6, 1, 1, 1]))
    for i, (label, runs) in enumerate(cases):
        r0 = runs[0]
        ax = axs[i, 0]
        stt = r0["spacetime"]
        if len(stt):
            ax.vlines(stt[:, 0], stt[:, 1] / 1e3, stt[:, 2] / 1e3, color="0.5", lw=0.8)
            sc = ax.scatter(stt[:, 0], stt[:, 3] / 1e3, c=stt[:, 4], cmap="plasma", vmin=7.0, vmax=8.4, s=16, zorder=3)
            plt.colorbar(sc, ax=ax, label="Mw", shrink=0.8)
        for n, y in r0["site_y"].items():
            ax.axhline(y / 1e3, color=SITE_COLORS[n], ls="--", lw=1)
        ax.set(xlabel="time [yr]", ylabel="northing [km]",
               title=f"{label}: M >= 7 ruptures, last {opt['window']:.0f} yr of the first run")
        ax = axs[i, 1]
        xx = np.linspace(0, 1200, 300)
        for s in sites:
            n = s["name"]
            c = SITE_COLORS[n]
            iv = np.sort(np.concatenate([r["sites"][n]["intervals"] for r in runs]))
            if len(iv):
                ax.plot(iv, np.arange(1, len(iv) + 1) / len(iv), color=c, label=f"{n} (model)")
            ax.plot(xx, [0.5 * (1 + erf((u - s["ri_yr"]) / (s["ri_sigma_yr"] * np.sqrt(2)))) for u in xx], color=c,
                    ls="--", lw=1, label=f"{n} (paleo)")
        ax.set(xlabel="interval [yr]", ylabel="cumulative fraction", title="site intervals (pooled over runs)")
        ax.legend(fontsize=7)
        ax = axs[i, 2]
        F = np.array([r["fmd"] for r in runs])
        ax.fill_between(MAG, np.percentile(F, 5, 0), np.percentile(F, 95, 0), color="C0", alpha=0.2, lw=0)
        ax.semilogy(MAG, np.where(np.median(F, 0) > 0, np.median(F, 0), np.nan), color="C0", label="median")
        ref = np.median(F, 0)[np.argmin(np.abs(MAG - 5.5))]
        ax.semilogy(MAG, ref * 10 ** (-(MAG - 5.5)), "k:", lw=1, label="b = 1 from M5.5")
        ax.set(xlabel="Mw", ylabel="events >= M per 10,000 yr", title="frequency-magnitude distribution",
               ylim=(0.3, None))
        ax.legend(fontsize=7)
        ax = axs[i, 3]
        P = np.array([[r["patterns"][p] for p in PATTERNS] for r in runs])
        ax.bar(range(len(PATTERNS)), np.median(P, 0), color="0.6", edgecolor="k")
        ax.errorbar(range(len(PATTERNS)), np.median(P, 0),
                    yerr=[np.median(P, 0) - np.percentile(P, 5, 0), np.percentile(P, 95, 0) - np.median(P, 0)],
                    fmt="none", ecolor="k", capsize=3)
        ax.set_xticks(range(len(PATTERNS)))
        ax.set_xticklabels(PATTERNS)
        ax.set(ylabel=f"M >= {opt['pattern_mag']:g} events per 1000 yr",
               title="sites ruptured together (S south, C central, N north)")
    fig.tight_layout()
    fig.savefig(out_prefix + ".png", dpi=110)
    plt.close(fig)

    with open(out_prefix + ".md", "w") as fh:
        th_txt = ", ".join(f"{k} {v:g} m" for k, v in opt["site"].items())
        fh.write(f"Detection threshold {opt['threshold']:g} m{', ' + th_txt if th_txt else ''}; medians over runs. "
                 f"Excess: rate of M >= 7.5 relative to a b = 1 extrapolation from M >= 5.5. "
                 f"Patterns: M >= {opt['pattern_mag']:g} events per 1000 yr by sites ruptured.\n\n")
        fh.write("| case | runs | south RI / CV / SED | central | north | Bath M>=6 / M>=7 | excess M>=7.5 | "
                 + " | ".join(PATTERNS) + " |\n")
        fh.write("|" + "---|" * (7 + len(PATTERNS)) + "\n")
        for label, runs in cases:
            cells = []
            for n in ("south", "central", "north"):
                ri = np.nanmedian([r["sites"][n]["ri"] for r in runs])
                cv = np.nanmedian([r["sites"][n]["cv"] for r in runs])
                sed = np.nanmedian([r["sites"][n]["sed"] for r in runs])
                cells.append(f"{ri:.0f} / {cv:.2f} / {sed:.1f}")
            b6 = np.nanmedian([r["bath6"] for r in runs])
            b7 = np.nanmedian([r["bath7"] for r in runs])
            ex = np.nanmedian([r["excess75"] for r in runs])
            pat = [f"{np.median([r['patterns'][p] for r in runs]):.1f}" for p in PATTERNS]
            fh.write(f"| {label} | {len(runs)} | " + " | ".join(cells) + f" | {b6:.2f} / {b7:.2f} | {ex:.1f} | "
                     + " | ".join(pat) + " |\n")
        obs = {s["name"]: s for s in load_sites(os.path.join(DATA, "paleo_sites.csv"))}
        fh.write("| observed | | " + " | ".join(
            f"{obs[n]['ri_yr']:.0f} / {obs[n]['ri_sigma_yr'] / obs[n]['ri_yr']:.2f} / {obs[n]['sed_mean_m']:.1f}"
            for n in ("south", "central", "north")) + " | | |" + " |" * len(PATTERNS) + "\n")
    print("wrote", out_prefix + ".png", out_prefix + ".md")


if __name__ == "__main__":
    main()
