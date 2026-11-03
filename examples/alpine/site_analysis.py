"""Site and system statistics for posterior ensembles, in the style of the
RSQSim intercomparison slides and extended.

    python3 examples/alpine/site_analysis.py out_prefix label1 dir1 [label2 dir2 ...] [--threshold 2]

Each dir holds the catalogs of a posterior sweep (design "list", files
cat_vNNNN_sS.npz; vector 0000 is the best-fitting one). Labels should be
"planar" or "variable_dip" for the RSQSim values in
data/nz_alpine/rsqsim_summary.csv to be shown alongside.

Outputs
  <out_prefix>_spacetime.png   map of slip rate with the sites and a space-time plot of M >= 7
                               ruptures (along-strike extent in northing, hypocentre coloured by Mw)
                               for the best run, one row per geometry
  <out_prefix>_statistics.png  slip, magnitude and inter-event time distributions at the sites,
                               cumulative frequency-magnitude distribution, hypocentre depths,
                               area versus Mw (as in the RSQSim summary), plus RI versus CV and SED
                               for every run with observations and RSQSim values, and the
                               dependence of the mean interval on the detection threshold
  <out_prefix>_extended.png    co-rupture of sites, time- and slip-predictability correlations,
                               interval distributions against the paleoseismic normal PDFs
  <out_prefix>_table.md / .csv per site: best model, ensemble median and 5-95 percent range,
                               observations, RSQSim
Distributions labelled "ensemble" pool all posterior runs.
"""
import csv
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
from stochseis.sites import load_sites, site_patch, site_records  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "data", "nz_alpine")
SITE_COLORS = {"south": "C0", "central": "C1", "north": "C2"}
GEO_STYLE = ["-", ":", "--", "-."]
GEO_COLORS = ["C0", "C3", "C2", "C1"]
MAG = np.arange(5.0, 8.61, 0.05)
THRESHOLDS = [0.5, 1.0, 2.0, 3.0, 4.0]


def parse_args(argv):
    th = 2.0
    if "--threshold" in argv:
        k = argv.index("--threshold")
        th = float(argv[k + 1])
        argv = argv[:k] + argv[k + 2:]
    return argv[1], list(zip(argv[2::2], argv[3::2])), th


def run_stats(path, sites, th):
    cat = Catalog.load(path)
    ev = cat.events
    f = cat.fault
    T = cat.duration
    t0 = 0.1 * T
    years = T - t0
    sel = ev["t"] >= t0
    out = {"years": years, "path": path}
    m = ev["mag"][sel]
    out["fmd"] = np.array([np.sum(m >= x) for x in MAG]) / years * 1e4      # per 10,000 yr
    # catalog events only: grown ruptures at or above the completeness magnitude (deficit-capped
    # ruptures through relaxed patches can have tiny moments and are not part of the catalog)
    fe = sel & (ev["kind"] == 1) & (ev["mag"] >= cat.params.mag_min)
    out["hypo_depth"] = f.z[ev["nuc"][fe]] / 1e3
    out["area_mag"] = (ev["area"][fe], ev["mag"][fe])
    # site records at the detection threshold and at the threshold series
    out["sites"] = {}
    for thr in THRESHOLDS + ([th] if th not in THRESHOLDS else []):
        rec = site_records(cat, sites, threshold=thr, t_start=t0)
        for s in sites:
            r = rec[s["name"]]
            t = r["times"]
            d = {"n": len(t), "times": t, "sed": r["sed"], "mag": ev["mag"][r["events"]] if len(t) else np.zeros(0)}
            if len(t) >= 3:
                dt = np.diff(t)
                d.update(intervals=dt, ri=dt.mean(), ri_sd=dt.std(ddof=1), cv=dt.std(ddof=1) / dt.mean(),
                         sed_mean=r["sed"].mean(), sed_sd=r["sed"].std(ddof=1),
                         corr_tp=np.corrcoef(r["sed"][:-1], dt)[0, 1] if len(dt) > 3 else np.nan,
                         corr_sp=np.corrcoef(dt, r["sed"][1:])[0, 1] if len(dt) > 3 else np.nan)
            else:
                d.update(intervals=np.zeros(0), ri=np.nan, ri_sd=np.nan, cv=np.nan, sed_mean=np.nan, sed_sd=np.nan,
                         corr_tp=np.nan, corr_sp=np.nan)
            out["sites"][(s["name"], thr)] = d
    # co-rupture at the detection threshold
    names = [s["name"] for s in sites]
    M = np.full((len(names), len(names)), np.nan)
    for i, a in enumerate(names):
        ea = set(out["sites"][(a, th)]["times"].round(9))
        for j, b in enumerate(names):
            eb = set(out["sites"][(b, th)]["times"].round(9))
            if ea:
                M[i, j] = len(ea & eb) / len(ea)
    out["corupture"] = M
    # space-time data for M >= 7 (stored slip maps)
    st = []
    for k, e in enumerate(cat.slip_event):
        if ev["mag"][e] < 7.0:
            continue
        a, b = cat.slip_ptr[k], cat.slip_ptr[k + 1]
        idx = cat.slip_idx[a:b]
        st.append((ev["t"][e], f.y[idx].min(), f.y[idx].max(), f.y[ev["nuc"][e]], ev["mag"][e]))
    out["spacetime"] = np.array(st) if st else np.zeros((0, 5))
    out["map"] = (f.x, f.y, f.slip_rate * 1e3, f.row == 0)
    out["site_xy"] = [(s["name"], f.x[site_patch(f, s["x_nztm_m"], s["y_nztm_m"])[0]],
                       f.y[site_patch(f, s["x_nztm_m"], s["y_nztm_m"])[0]]) for s in sites]
    out["duration"] = T
    return out


def load_rsqsim():
    path = os.path.join(DATA, "rsqsim_summary.csv")
    out = {}
    if os.path.exists(path):
        for r in csv.DictReader(open(path)):
            out[(r["geometry"], r["site"])] = {k: float(v) for k, v in r.items() if k not in ("geometry", "site")}
    return out


def pct(v, q):
    v = np.asarray([x for x in v if np.isfinite(x)])
    return np.percentile(v, q) if len(v) else np.nan


def main():
    out_prefix, pairs, th = parse_args(sys.argv)
    sites = load_sites(os.path.join(DATA, "paleo_sites.csv"))
    names = [s["name"] for s in sites]
    rsq = load_rsqsim()
    G = {}
    for label, d in pairs:
        files = sorted(glob.glob(os.path.join(d, "cat_v*_s*.npz")))
        if not files:
            sys.exit(f"no catalogs in {d}")
        runs = [run_stats(f, sites, th) for f in files]
        best = [r for r in runs if re.search(r"cat_v0000_s", r["path"])]
        G[label] = dict(runs=runs, best=best[0] if best else runs[0])
        print(f"{label}: {len(runs)} catalogs, {sum(r['years'] for r in runs):.0f} yr after burn-in")

    # ------------------------------------------------------------ space-time
    fig, axs = plt.subplots(len(G), 2, figsize=(14, 5.5 * len(G)), squeeze=False,
                            gridspec_kw=dict(width_ratios=[1, 1.4]))
    for gi, (label, D) in enumerate(G.items()):
        b = D["best"]
        x, y, v, top = b["map"]
        ax = axs[gi, 0]
        sc = ax.scatter(x / 1e3, y / 1e3, c=v, s=2, cmap="viridis", vmin=0, vmax=30)
        plt.colorbar(sc, ax=ax, label="slip rate [mm/yr]", shrink=0.8)
        for (n, sx, sy) in b["site_xy"]:
            ax.plot(sx / 1e3, sy / 1e3, "o", color=SITE_COLORS[n], mec="k", ms=8)
            ax.axhline(sy / 1e3, color=SITE_COLORS[n], ls="--", lw=1)
            ax.annotate(n, (sx / 1e3, sy / 1e3), xytext=(-40, 5), textcoords="offset points", fontsize=8)
        ax.set(xlabel="easting [km NZTM]", ylabel="northing [km NZTM]", title=f"{label}: slip rate and sites")
        ax.set_aspect("equal")
        ax = axs[gi, 1]
        stt = b["spacetime"]
        T = b["duration"]
        win = min(10000.0, T)
        if len(stt):
            w = stt[:, 0] >= T - win
            stt = stt[w]
            ax.vlines(stt[:, 0], stt[:, 1] / 1e3, stt[:, 2] / 1e3, color="0.5", lw=0.7)
            sc = ax.scatter(stt[:, 0], stt[:, 3] / 1e3, c=stt[:, 4], cmap="plasma", s=18, vmin=7.0, vmax=8.4, zorder=3)
            plt.colorbar(sc, ax=ax, label="Mw", shrink=0.8)
        for (n, sx, sy) in b["site_xy"]:
            ax.axhline(sy / 1e3, color=SITE_COLORS[n], ls="--", lw=1)
        ax.set(xlabel="time [yr]", ylabel="northing [km NZTM]",
               title=f"{label}: M >= 7 ruptures, last {win:,.0f} yr of the best run (bars: extent; dots: hypocentre)")
        ax.set_ylim(axs[gi, 0].get_ylim())
    fig.tight_layout()
    fig.savefig(out_prefix + "_spacetime.png", dpi=120)
    plt.close(fig)

    # ------------------------------------------------------------ statistics
    fig, axs = plt.subplots(3, 3, figsize=(17, 14))
    for gi, (label, D) in enumerate(G.items()):
        ls = GEO_STYLE[gi % len(GEO_STYLE)]
        gc = GEO_COLORS[gi % len(GEO_COLORS)]
        runs = D["runs"]
        for n in names:
            c = SITE_COLORS[n]
            sed = np.concatenate([r["sites"][(n, th)]["sed"] for r in runs])
            mg = np.concatenate([r["sites"][(n, th)]["mag"] for r in runs])
            iv = np.concatenate([r["sites"][(n, th)]["intervals"] for r in runs])
            axs[0, 0].hist(sed, bins=np.arange(0, 15.5, 0.5), density=True, histtype="step", color=c, ls=ls, lw=1.5)
            axs[0, 1].hist(mg, bins=np.arange(6.0, 8.61, 0.1), density=True, histtype="step", color=c, ls=ls, lw=1.5)
            axs[0, 2].hist(iv, bins=np.arange(0, 1050, 40), density=True, histtype="step", color=c, ls=ls, lw=1.5)
        # FMD
        F = np.array([r["fmd"] for r in runs])
        axs[1, 0].fill_between(MAG, np.percentile(F, 5, 0), np.percentile(F, 95, 0), color=gc, alpha=0.2, lw=0)
        axs[1, 0].semilogy(MAG, np.where(D["best"]["fmd"] > 0, D["best"]["fmd"], np.nan), color=gc,
                           label=f"{label}: best (band: {len(runs)} runs)")
        dep = np.concatenate([r["hypo_depth"] for r in runs])
        axs[1, 1].hist(dep, bins=np.arange(0, 21, 1), density=True, orientation="horizontal", histtype="step",
                       color=gc, lw=1.5, label=label)
        A, Mw = D["best"]["area_mag"]
        axs[1, 2].semilogy(Mw, A, ".", color=gc, ms=1.5, alpha=0.4, label=f"{label} (best)")
        # RI versus CV and SED, every run
        for n in names:
            c = SITE_COLORS[n]
            ri = [r["sites"][(n, th)]["ri"] for r in runs]
            cv = [r["sites"][(n, th)]["cv"] for r in runs]
            se = [r["sites"][(n, th)]["sed_mean"] for r in runs]
            mk = "o" if gi == 0 else ("s" if gi == 1 else "^")
            axs[2, 0].plot(ri, cv, mk, color=c, ms=4, alpha=0.5, mfc="none" if gi else c)
            axs[2, 1].plot(ri, se, mk, color=c, ms=4, alpha=0.5, mfc="none" if gi else c)
            q = rsq.get((label, n))
            if q:
                axs[2, 0].plot(q["ri_mean_yr"], q["ri_sd_yr"] / q["ri_mean_yr"], "X", color=c, mec="k", ms=11)
                axs[2, 1].plot(q["ri_mean_yr"], q["slip_mean_m"], "X", color=c, mec="k", ms=11)
            rth = [np.nanmedian([r["sites"][(n, t)]["ri"] for r in runs]) for t in THRESHOLDS]
            lo = [pct([r["sites"][(n, t)]["ri"] for r in runs], 5) for t in THRESHOLDS]
            hi = [pct([r["sites"][(n, t)]["ri"] for r in runs], 95) for t in THRESHOLDS]
            axs[2, 2].plot(THRESHOLDS, rth, ls=ls, color=c, marker="o", ms=4)
            axs[2, 2].fill_between(THRESHOLDS, lo, hi, color=c, alpha=0.08, lw=0)
    # observations
    for s in sites:
        c = SITE_COLORS[s["name"]]
        y0 = {"south": 0.66, "central": 0.69, "north": 0.62}[s["name"]]
        axs[0, 0].errorbar([s["sed_mean_m"]], [y0], xerr=[[s["sed_mean_m"] - s["sed_min_m"]], [s["sed_max_m"] - s["sed_mean_m"]]],
                           fmt="o", color=c, mec="k", ms=7, capsize=3, ls="--")
        xx = np.linspace(0, 1000, 400)
        axs[0, 2].fill_between(xx, np.exp(-0.5 * ((xx - s["ri_yr"]) / s["ri_sigma_yr"]) ** 2)
                               / (s["ri_sigma_yr"] * np.sqrt(2 * np.pi)), color=c, alpha=0.2, lw=0)
        se_ri = s["ri_sigma_yr"] / np.sqrt(s["n_ri"] - 1)
        axs[2, 0].errorbar([s["ri_yr"]], [s["ri_sigma_yr"] / s["ri_yr"]], xerr=[se_ri], yerr=[0.1], fmt="s", color=c,
                           mec="k", ms=9, capsize=4)
        axs[2, 1].errorbar([s["ri_yr"]], [s["sed_mean_m"]], xerr=[se_ri],
                           yerr=[[s["sed_mean_m"] - s["sed_min_m"]], [s["sed_max_m"] - s["sed_mean_m"]]],
                           fmt="s", color=c, mec="k", ms=9, capsize=4)
        axs[2, 2].errorbar([th], [s["ri_yr"]], yerr=[se_ri], fmt="s", color=c, mec="k", ms=8, capsize=4)
    # legends and labels
    from matplotlib.lines import Line2D
    site_h = [Line2D([], [], color=SITE_COLORS[n], lw=2, label=f"{n} site") for n in names]
    geo_h = [Line2D([], [], color="k", ls=GEO_STYLE[i], label=f"{lab} (ensemble)") for i, lab in enumerate(G)]
    axs[0, 0].legend(handles=site_h + geo_h + [Line2D([], [], marker="o", color="k", ls="--", label="paleo SED mean and range")],
                     fontsize=7)
    axs[0, 0].set(xlabel="slip at site [m]", ylabel="probability density", title=f"slip distribution (events >= {th:g} m)")
    axs[0, 1].legend(handles=site_h + geo_h, fontsize=7)
    axs[0, 1].set(xlabel="Mw", ylabel="probability density", title="magnitude of events recorded at the sites")
    axs[0, 2].legend(handles=site_h + geo_h + [plt.Rectangle((0, 0), 1, 1, color="0.7", label="paleo (normal PDF)")],
                     fontsize=7)
    axs[0, 2].set(xlabel="years since last event", ylabel="probability density", title="inter-event time distribution")
    axs[1, 0].legend(fontsize=7)
    axs[1, 0].set(xlabel="Mw", ylabel="cumulative number per 10,000 yr", title="frequency-magnitude distribution")
    axs[1, 1].invert_yaxis()
    axs[1, 1].legend(fontsize=7)
    axs[1, 1].set(xlabel="density", ylabel="hypocentre depth [km]", title="hypocentre depths, M >= mag_min (ensemble)")
    axs[1, 2].legend(fontsize=7, markerscale=6)
    axs[1, 2].set(xlabel="Mw", ylabel="rupture area [m$^2$]", title="area versus magnitude")
    mk_h = [Line2D([], [], marker=["o", "s", "^"][i % 3], color="k", ls="", mfc="k" if i == 0 else "none",
                   label=f"{lab}: one run each") for i, lab in enumerate(G)]
    obs_h = [Line2D([], [], marker="s", color="k", ls="", label="observed"),
             Line2D([], [], marker="X", color="k", ls="", ms=9, label="RSQSim (same geometry)")]
    axs[2, 0].legend(handles=site_h + mk_h + obs_h, fontsize=7)
    axs[2, 0].set(xlabel="mean recurrence interval [yr]", ylabel="CV of intervals", title=f"RI versus CV ({th:g} m threshold)")
    axs[2, 1].legend(handles=site_h + mk_h + obs_h, fontsize=7)
    axs[2, 1].set(xlabel="mean recurrence interval [yr]", ylabel="mean SED [m]", title="RI versus single-event displacement")
    axs[2, 2].legend(handles=site_h + geo_h + [obs_h[0]], fontsize=7)
    axs[2, 2].set(xlabel="detection threshold [m]", ylabel="mean recurrence interval [yr]",
                  title="dependence on the detection threshold (median, 5-95 %)")
    fig.tight_layout()
    fig.savefig(out_prefix + "_statistics.png", dpi=120)
    plt.close(fig)

    # ------------------------------------------------------------ extended
    fig, axs = plt.subplots(2, max(len(G), 2) + 1, figsize=(6 * (max(len(G), 2) + 1), 10), squeeze=False)
    for gi, (label, D) in enumerate(G.items()):
        ax = axs[0, gi]
        Mm = np.nanmedian(np.array([r["corupture"] for r in D["runs"]]), axis=0)
        im = ax.imshow(Mm, vmin=0, vmax=1, cmap="Blues")
        for i in range(len(names)):
            for j in range(len(names)):
                ax.text(j, i, f"{Mm[i, j]:.2f}", ha="center", va="center", fontsize=10)
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names)
        ax.set(title=f"{label}: fraction of row-site events\nthat also rupture the column site (median)")
        plt.colorbar(im, ax=ax, shrink=0.7)
    ax = axs[0, -1]
    width = 0.8 / max(len(G), 1)
    for gi, (label, D) in enumerate(G.items()):
        for k, n in enumerate(names):
            tp = [r["sites"][(n, th)]["corr_tp"] for r in D["runs"]]
            sp = [r["sites"][(n, th)]["corr_sp"] for r in D["runs"]]
            x0 = k + gi * width
            ax.bar(x0 - 0.2, np.nanmedian(tp), width / 2, color=SITE_COLORS[n], alpha=0.9, edgecolor="k",
                   hatch="" if gi == 0 else "//")
            ax.bar(x0 + width / 2 - 0.2, np.nanmedian(sp), width / 2, color=SITE_COLORS[n], alpha=0.4, edgecolor="k",
                   hatch="" if gi == 0 else "//")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names)
    ax.set(ylabel="median correlation over runs",
           title="predictability (median over runs)\ndark: corr(SED, next interval), time-predictable\n"
                 "light: corr(previous interval, SED), slip-predictable; hatched: 2nd geometry")
    for k, s in enumerate(sites):
        ax = axs[1, k] if k < axs.shape[1] else None
        if ax is None:
            continue
        n = s["name"]
        xx = np.linspace(0, 1000, 400)
        from math import erf
        ax.plot(xx, [0.5 * (1 + erf((u - s["ri_yr"]) / (s["ri_sigma_yr"] * np.sqrt(2)))) for u in xx], "k--",
                label="paleo (normal)")
        for gi, (label, D) in enumerate(G.items()):
            iv = np.sort(np.concatenate([r["sites"][(n, th)]["intervals"] for r in D["runs"]]))
            if len(iv):
                ax.plot(iv, np.arange(1, len(iv) + 1) / len(iv), color=GEO_COLORS[gi], ls=GEO_STYLE[gi], label=f"{label} (ensemble)")
            ivb = np.sort(D["best"]["sites"][(n, th)]["intervals"])
            if len(ivb):
                ax.step(ivb, np.arange(1, len(ivb) + 1) / len(ivb), color=GEO_COLORS[gi], lw=0.8, alpha=0.6,
                        where="post", label=f"{label} (best run)")
        ax.set(xlabel="interval [yr]", ylabel="cumulative fraction", title=f"{n} site: interval distribution")
        ax.legend(fontsize=7)
    for k in range(len(sites), axs.shape[1]):
        axs[1, k].axis("off")
    fig.tight_layout()
    fig.savefig(out_prefix + "_extended.png", dpi=120)
    plt.close(fig)

    # ------------------------------------------------------------ table
    rows = []
    for label, D in G.items():
        for s in sites:
            n = s["name"]
            b = D["best"]["sites"][(n, th)]
            ens = [r["sites"][(n, th)] for r in D["runs"]]
            q = rsq.get((label, n), {})
            row = dict(geometry=label, site=n,
                       ri_best=b["ri"], ri_sd_best=b["ri_sd"], cv_best=b["cv"], slip_best=b["sed_mean"],
                       slip_sd_best=b["sed_sd"], n_best=b["n"],
                       ri_ens_med=np.nanmedian([e["ri"] for e in ens]), ri_ens_lo=pct([e["ri"] for e in ens], 5),
                       ri_ens_hi=pct([e["ri"] for e in ens], 95),
                       cv_ens_med=np.nanmedian([e["cv"] for e in ens]), cv_ens_lo=pct([e["cv"] for e in ens], 5),
                       cv_ens_hi=pct([e["cv"] for e in ens], 95),
                       slip_ens_med=np.nanmedian([e["sed_mean"] for e in ens]),
                       slip_ens_lo=pct([e["sed_mean"] for e in ens], 5), slip_ens_hi=pct([e["sed_mean"] for e in ens], 95),
                       ri_obs=s["ri_yr"], ri_sd_obs=s["ri_sigma_yr"], cv_obs=s["ri_sigma_yr"] / s["ri_yr"],
                       slip_obs=s["sed_mean_m"],
                       ri_rsqsim=q.get("ri_mean_yr", np.nan), ri_sd_rsqsim=q.get("ri_sd_yr", np.nan),
                       cv_rsqsim=q.get("ri_sd_yr", np.nan) / q.get("ri_mean_yr", np.nan) if q else np.nan,
                       slip_rsqsim=q.get("slip_mean_m", np.nan), slip_sd_rsqsim=q.get("slip_sd_m", np.nan))
            rows.append(row)
    with open(out_prefix + "_table.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(out_prefix + "_table.md", "w") as fh:
        fh.write(f"Detection threshold {th:g} m. Best: best-fitting vector, first seed. Ensemble: median "
                 "[5-95 %] over all posterior runs.\n\n")
        fh.write("| geometry | site | RI best (SD) | RI ensemble | RI obs (SD) | RI RSQSim (SD) | CV best | CV ensemble "
                 "| CV obs | CV RSQSim | slip best (SD) | slip ensemble | slip obs | slip RSQSim (SD) |\n")
        fh.write("|" + "---|" * 14 + "\n")
        for r in rows:
            fh.write(f"| {r['geometry']} | {r['site']} | {r['ri_best']:.0f} ({r['ri_sd_best']:.0f}) | "
                     f"{r['ri_ens_med']:.0f} [{r['ri_ens_lo']:.0f}-{r['ri_ens_hi']:.0f}] | "
                     f"{r['ri_obs']:.0f} ({r['ri_sd_obs']:.0f}) | {r['ri_rsqsim']:.0f} ({r['ri_sd_rsqsim']:.0f}) | "
                     f"{r['cv_best']:.2f} | {r['cv_ens_med']:.2f} [{r['cv_ens_lo']:.2f}-{r['cv_ens_hi']:.2f}] | "
                     f"{r['cv_obs']:.2f} | {r['cv_rsqsim']:.2f} | {r['slip_best']:.1f} ({r['slip_sd_best']:.1f}) | "
                     f"{r['slip_ens_med']:.1f} [{r['slip_ens_lo']:.1f}-{r['slip_ens_hi']:.1f}] | {r['slip_obs']:.1f} | "
                     f"{r['slip_rsqsim']:.1f} ({r['slip_sd_rsqsim']:.1f}) |\n")
    print("wrote", ", ".join(out_prefix + x for x in ("_spacetime.png", "_statistics.png", "_extended.png",
                                                       "_table.md", "_table.csv")))


if __name__ == "__main__":
    main()
