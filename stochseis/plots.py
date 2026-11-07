# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Diagnostic figures (matplotlib, Agg backend)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import analysis as an  # noqa: E402
from .catalog import Catalog, moment_from_magnitude  # noqa: E402
from .etas import DAYS_PER_YEAR  # noqa: E402


def summary_figure(cat: Catalog, path=None, mmain=6.0, t_start=None):
    """Nine-panel overview: FMD, M0-L, L/W-M0, slip-length, Omori, Bath,
    interevent times, cumulative moment, slip map of the largest event."""
    if t_start is None:
        t_start = 0.1 * cat.duration
    ev = cat.events
    fe = cat.fault_events()
    p = cat.params
    fig, axs = plt.subplots(3, 3, figsize=(15, 12))
    axs = axs.ravel()

    # 1 FMD
    ax = axs[0]
    c, inc, cum = an.fmd(ev["mag"], p.mag_min)
    T = cat.duration
    ax.semilogy(c, cum / T, "k.-", label="cumulative")
    ax.semilogy(c[inc > 0], inc[inc > 0] / T, "o", mfc="none", color="gray", label="incremental")
    reg = an.b_value_regimes(cat)
    a0 = np.log10(cum[0] / T) + p.b_value * c[0]
    ax.semilogy(c, 10 ** (a0 - p.b_value * c), "r--", lw=1, label=f"GR b={p.b_value}")
    ax.axvline(cat.meta["m_patch"], color="b", lw=0.8, ls=":", label="patch magnitude")
    if np.isfinite(reg["m_sat"]):
        ax.axvline(reg["m_sat"], color="g", lw=0.8, ls=":", label="W saturation")
    ax.set_xlabel("magnitude")
    ax.set_ylabel("rate [1/yr]")
    ax.set_title(f"FMD: b sub={reg['b_sub']:.2f}, fault<sat={reg['b_fault_below_sat']:.2f}, "
                 f">sat={reg['b_fault_above_sat']:.2f}", fontsize=9)
    ax.legend(fontsize=7)

    # 2 M0 vs L
    ax = axs[1]
    sel = fe["npatch"] >= 2
    ax.loglog(fe["length"][sel] / 1e3, fe["m0"][sel], ".", ms=2, alpha=0.3, color="C0")
    Lr = np.array([1, 2000.0])
    wf = cat.fault.width / 1e3
    # constant stress-drop guides
    mu = cat.fault.mu.mean()
    ds = p.stress_drop_mpa * 1e6
    m0_small = mu * (np.pi * (Lr * 1e3 / 2) ** 2) * (16 / (7 * np.pi)) * ds * (Lr * 1e3 / 2) / mu
    ax.loglog(Lr, m0_small, "r--", lw=1, label="L^3, circular crack")
    w_eff = cat.meta["w_eff"]
    m0_large = mu * (Lr * 1e3 * wf * 1e3) * (np.pi / 4) * ds * w_eff / mu
    ax.loglog(Lr, m0_large, "g--", lw=1, label="L^1, saturated W (W-model)")
    ax.loglog(Lr, m0_large * (Lr / wf), "b:", lw=1, label="L^2 (slip ~ L, as in data)")
    sc = an.scaling_fits(cat)
    ax.set_xlabel("rupture length [km]")
    ax.set_ylabel("M0 [N m]")
    ax.set_title(f"slopes: small {sc.get('slope_m0_L_small', np.nan):.2f}, "
                 f"large {sc.get('slope_m0_L_large', np.nan):.2f}", fontsize=9)
    ax.set_xlim(0.5, 2 * p.fault_length_km)
    ax.legend(fontsize=7)
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xscale("log")
    ax2.set_xticks([])

    # 3 aspect ratio vs M0 (Fig. 3.43b style)
    ax = axs[2]
    ax.loglog(fe["m0"][sel], (fe["length"] / fe["width"])[sel], ".", ms=2, alpha=0.3)
    for mw in (5, 6, 7, 8):
        ax.axvline(moment_from_magnitude(mw), color="gray", lw=0.5)
        ax.text(moment_from_magnitude(mw), ax.get_ylim()[1] if False else 40, f"M{mw}", fontsize=7)
    ax.axhline(1, color="k", lw=0.5)
    ax.set_xlabel("M0 [N m]")
    ax.set_ylabel("L / W")
    ax.set_title(f"aspect ratio; W_fault = {wf:.0f} km", fontsize=9)

    # 4 mean slip vs length, stress drop
    ax = axs[3]
    ax.loglog(fe["length"][sel] / 1e3, fe["dmean"][sel], ".", ms=2, alpha=0.3)
    ax.set_xlabel("rupture length [km]")
    ax.set_ylabel("mean slip [m]")
    ax.set_title(f"median D/L = {sc.get('median_dmean_over_L', np.nan):.1e}; "
                 f"Blaser et al. (2010) data: 1.3e-5", fontsize=9)

    # 5 Omori
    ax = axs[4]
    delays, fit = an.omori_stack(cat, mmain=mmain, t_start=t_start)
    if len(delays) > 20:
        tc, rate = an.omori_binned(delays)
        ax.loglog(tc, rate / max(1, fit["n_mainshocks"]), "o", ms=3)
        tt = np.logspace(-3, np.log10(365.25), 100)
        ax.loglog(tt, fit["K"] * (tt + fit["c"]) ** (-fit["p"]) / max(1, fit["n_mainshocks"]), "r-")
        ax.set_title(f"Omori stack M>={mmain}: p={fit['p']:.2f}, c={fit['c']:.3f} d, "
                     f"n_main={fit['n_mainshocks']}", fontsize=9)
    else:
        ax.set_title("Omori: too few aftershocks", fontsize=9)
    ax.set_xlabel("time after mainshock [days]")
    ax.set_ylabel("rate per mainshock [1/day]")

    # 6 Bath
    ax = axs[5]
    gaps, none = an.bath_gaps(cat, mmain=mmain, t_start=t_start)
    if len(gaps):
        ax.hist(gaps, bins=np.arange(-0.5, 4.0, 0.2), color="C0", alpha=0.7, label="window selection")
    gg = an.bath_gaps_genealogy(cat, mmain=mmain)
    if len(gg):
        ax.hist(gg, bins=np.arange(-0.5, 4.0, 0.2), histtype="step", color="k", label="genealogy")
    ax.axvline(1.2, color="r", ls="--", lw=1, label="1.2")
    ax.set_xlabel("M_main - M_largest aftershock")
    ax.set_title(f"Bath: mean {gaps.mean() if len(gaps) else np.nan:.2f} (n={len(gaps)}, "
                 f"{none} without aftershock)", fontsize=9)
    ax.legend(fontsize=7)

    # 7 interevent times of large events
    ax = axs[6]
    for mmin, col in ((6.0, "C0"), (6.5, "C1"), (7.0, "C3")):
        idx = an.select_mainshocks(cat, mmin, t_start=t_start)
        if len(idx) > 5:
            dt = np.diff(np.sort(ev["t"][idx]))
            cv = dt.std() / dt.mean()
            ax.hist(dt / dt.mean(), bins=np.linspace(0, 5, 26), histtype="step", color=col,
                    density=True, label=f"M>={mmin}: n={len(idx)}, CV={cv:.2f}")
    x = np.linspace(0, 5, 100)
    ax.plot(x, np.exp(-x), "k:", lw=1, label="Poisson")
    ax.set_xlabel("interevent time / mean (declustered)")
    ax.legend(fontsize=7)
    sr = an.site_recurrence_cv(cat, slip_threshold=0.3, t_start=t_start)
    sr1 = an.site_recurrence_cv(cat, slip_threshold=1.0, t_start=t_start)
    ax.set_title(f"whole-fault recurrence; site CV (slip>=0.3 m) {sr['cv']:.2f}, (>=1 m) {sr1['cv']:.2f}",
                 fontsize=9)

    # 8 cumulative moment
    ax = axs[7]
    t, cum, load = an.cumulative_moment(cat)
    ax.plot(t, cum, "k", label="released")
    ax.plot(t, load, "r--", label="tectonic loading")
    ax.set_xlabel("time [yr]")
    ax.set_ylabel("cumulative moment [N m]")
    ax.set_title(f"closure {cat.meta['moment_closure']:.2f}; final mean deficit "
                 f"{cat.meta.get('final_mean_deficit', np.nan):.2f} m", fontsize=9)
    ax.legend(fontsize=7)

    # 9 slip map of the largest stored event and final deficit
    ax = axs[8]
    plot_largest_slip(cat, ax)

    fig.suptitle(f"{cat.params.fault_length_km:.0f} x {cat.params.fault_width_km:.0f} km fault, "
                 f"{cat.params.slip_rate_mm_yr:.0f} mm/yr, stress drop {cat.params.stress_drop_mpa} MPa "
                 f"(ln-sigma {cat.params.stress_drop_ln_sigma}), b={cat.params.b_value}, "
                 f"rho_c={cat.params.arrest_scale}, n={cat.params.branching_ratio}, "
                 f"seed {cat.params.seed}, {len(ev)} events / {cat.duration:.0f} yr", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    if path:
        fig.savefig(path, dpi=130)
        plt.close(fig)
    return fig


def slip_grid(cat: Catalog, event_index):
    """Slip of a stored event on the (ns, nz) fault grid (nan outside)."""
    f = cat.fault
    out = np.full((f.ns, f.nz), np.nan)
    res = cat.slip_of(event_index)
    if res is None:
        return out
    idx, val = res
    out[f.col[idx], f.row[idx]] = val
    return out


def plot_largest_slip(cat: Catalog, ax=None, event_index=None):
    f = cat.fault
    ev = cat.events
    if event_index is None:
        stored = cat.slip_event
        if len(stored) == 0:
            if ax is not None:
                ax.set_title("no stored slip maps")
            return
        event_index = stored[np.argmax(ev["mag"][stored])]
    g = slip_grid(cat, event_index)
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 2.5))
    extent = (0, f.ns * f.dx / 1e3, f.nz * f.dx / 1e3, 0)
    im = ax.imshow(g.T, extent=extent, aspect="auto", cmap="viridis", origin="upper")
    plt.colorbar(im, ax=ax, label="slip [m]", shrink=0.8)
    ax.set_xlabel("along strike [km]")
    ax.set_ylabel("down dip [km]")
    e = ev[event_index]
    ax.set_title(f"event {event_index}: M{e['mag']:.2f}, t={e['t']:.1f} yr, L={e['length']/1e3:.0f} km, "
                 f"Dmean={e['dmean']:.2f} m, Dmax={e['dmax']:.2f} m, stress drop {e['dsig']/1e6:.2f} MPa",
                 fontsize=8)
    return ax


def deficit_figure(cat: Catalog, path=None):
    """Final deficit on the fault and recovered versus prescribed slip rate."""
    f = cat.fault
    fig, axs = plt.subplots(2, 1, figsize=(12, 5))
    d = cat.final_deficit()
    g = np.full((f.ns, f.nz), np.nan)
    g[f.col, f.row] = d
    extent = (0, f.ns * f.dx / 1e3, f.nz * f.dx / 1e3, 0)
    im = axs[0].imshow(g.T, extent=extent, aspect="auto", cmap="RdBu_r", origin="upper",
                       vmin=-np.nanmax(np.abs(g)), vmax=np.nanmax(np.abs(g)))
    plt.colorbar(im, ax=axs[0], label="final deficit [m]")
    axs[0].set_ylabel("down dip [km]")
    rec, pre = an.slip_rate_recovered(cat)
    g2 = np.full((f.ns, f.nz), np.nan)
    g2[f.col, f.row] = rec / pre
    im = axs[1].imshow(g2.T, extent=extent, aspect="auto", cmap="viridis", origin="upper", vmin=0.5, vmax=1.5)
    plt.colorbar(im, ax=axs[1], label="released / prescribed slip rate")
    axs[1].set_xlabel("along strike [km]")
    axs[1].set_ylabel("down dip [km]")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=130)
        plt.close(fig)
    return fig


def sweep_figure(table, xkey, ykeys, path=None, groupkey=None):
    """Plot summary metrics from a sweep table (list of dicts) against one parameter."""
    import itertools

    fig, axs = plt.subplots(1, len(ykeys), figsize=(4 * len(ykeys), 3.5))
    axs = np.atleast_1d(axs)
    groups = sorted({r.get(groupkey) for r in table}) if groupkey else [None]
    markers = itertools.cycle("osd^v")
    for gval in groups:
        mk = next(markers)
        rows = [r for r in table if (groupkey is None or r.get(groupkey) == gval)]
        x = np.array([r[xkey] for r in rows], float)
        for ax, yk in zip(axs, ykeys):
            y = np.array([r.get(yk, np.nan) for r in rows], float)
            o = np.argsort(x)
            ax.plot(x[o], y[o], mk + "-", ms=4, label=f"{groupkey}={gval}" if groupkey else None)
            ax.set_xlabel(xkey)
            ax.set_ylabel(yk)
    if groupkey:
        axs[0].legend(fontsize=7)
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=130)
        plt.close(fig)
    return fig


def time_series_figure(cat: Catalog, path=None, window_events=2000, mag_mark=6.0, mc=None):
    """Magnitude versus time with a running b-value underneath.

    The running b is the Aki-Utsu estimate over sliding windows of
    ``window_events`` consecutive events with M >= mc (default mag_min),
    plotted at the window's median time with one-sigma bounds; a second
    curve uses only events at or above the patch magnitude (grown
    ruptures), which is where the slip budget can act.
    """
    ev = cat.events
    p = cat.params
    if mc is None:
        mc = p.mag_min
    fig, axs = plt.subplots(3, 1, figsize=(14, 9), sharex=True,
                            gridspec_kw=dict(height_ratios=[2, 1.2, 1]))
    ax = axs[0]
    small = (ev["mag"] < mag_mark) & (ev["mag"] >= mc)
    ax.plot(ev["t"][small], ev["mag"][small], ".", ms=1, color="0.6", rasterized=True)
    big = ev["mag"] >= mag_mark
    ax.vlines(ev["t"][big], mc, ev["mag"][big], color="C3", lw=0.8)
    ax.plot(ev["t"][big], ev["mag"][big], "o", color="C3", ms=3)
    ax.set_ylabel("magnitude")
    ax.set_ylim(mc - 0.1, None)
    ax.set_title(f"{int(np.sum(ev['mag'] >= mc))} events with M >= {mc}; red: M >= {mag_mark}", fontsize=9)

    ax = axs[1]
    for sel, col, lab in ((ev["mag"] >= mc, "k", f"all M >= {mc}"),
                          (ev["mag"] >= cat.meta["m_patch"], "C0",
                           f"fault events M >= {cat.meta['m_patch']:.1f}")):
        idx = np.where(sel)[0]
        n = window_events if col == "k" else max(200, window_events // 10)
        if len(idx) < 2 * n:
            continue
        m_lo = mc if col == "k" else cat.meta["m_patch"]
        starts = np.arange(0, len(idx) - n, max(1, n // 4))
        tt, bb, ss = [], [], []
        for s0 in starts:
            w = idx[s0:s0 + n]
            b, sb, _ = an.b_value_mle(ev["mag"][w], m_lo)
            tt.append(np.median(ev["t"][w]))
            bb.append(b)
            ss.append(sb)
        tt, bb, ss = map(np.array, (tt, bb, ss))
        ax.plot(tt, bb, color=col, lw=1, label=f"{lab} (window {n} events)")
        ax.fill_between(tt, bb - ss, bb + ss, color=col, alpha=0.15)
    ax.axhline(p.b_value, color="r", ls="--", lw=0.8, label="input b")
    ax.set_ylabel("running b")
    ax.set_ylim(0.5, 1.8)
    ax.legend(fontsize=7, loc="upper right")

    ax = axs[2]
    t, cum, load = an.cumulative_moment(cat)
    ax.plot(t, (cum - load) / 1e18, "k", lw=1)
    ax.axhline(0, color="r", ls="--", lw=0.8)
    ax.set_ylabel("released - loaded\nmoment [1e18 N m]")
    ax.set_xlabel("time [yr]")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=130)
        plt.close(fig)
    return fig


def largest_events_figure(cat: Catalog, n=6, path=None):
    """Slip maps of the n largest stored events, with the deficit just before each."""
    ev = cat.events
    stored = cat.slip_event
    order = stored[np.argsort(ev["mag"][stored])[::-1][:n]]
    fig, axs = plt.subplots(len(order), 1, figsize=(13, 2.4 * len(order)))
    axs = np.atleast_1d(axs)
    for ax, e in zip(axs, order):
        plot_largest_slip(cat, ax, event_index=int(e))
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=120)
        plt.close(fig)
    return fig


def fmd_regimes_figure(cat: Catalog, path=None):
    """FMD with the three magnitude ranges in which b is reported."""
    ev = cat.events
    p = cat.params
    reg = an.b_value_regimes(cat)
    c, inc, cum = an.fmd(ev["mag"], p.mag_min)
    T = cat.duration
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(c, cum / T, "k.-", label="cumulative")
    ax.semilogy(c[inc > 0], inc[inc > 0] / T, "o", mfc="none", color="gray", label="incremental")
    mp, ms = cat.meta["m_patch"], reg["m_sat"]
    ax.axvspan(p.mag_min, mp, color="C0", alpha=0.08)
    ax.axvspan(mp, ms, color="C2", alpha=0.08)
    ax.axvspan(ms, ev["mag"].max() + 0.2, color="C3", alpha=0.08)
    ymax = ax.get_ylim()[1]
    ax.text(p.mag_min + 0.05, ymax * 0.3, f"sub-patch events\nGR by construction\nb = {reg['b_sub']:.2f}", fontsize=8)
    ax.text(mp + 0.05, ymax * 0.3, f"grown ruptures,\nL < 2W\nb = {reg['b_fault_below_sat']:.2f}", fontsize=8)
    ax.text(ms + 0.05, ymax * 0.3, f"L > 2W\nb = {reg['b_fault_above_sat']:.2f}", fontsize=8)
    ax.set_xlabel("magnitude")
    ax.set_ylabel("rate [1/yr]")
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=130)
        plt.close(fig)
    return fig
