"""Catalog statistics.

Every function takes a :class:`~stochseis.catalog.Catalog` (or plain
arrays where noted) and returns numbers or arrays; nothing here plots.
:func:`summary_metrics` collects the scalar diagnostics used for parameter
sweeps.

Two views of clustering are available: the *observer* view, which selects
mainshocks and aftershocks with space-time windows exactly as one would
do with a real catalog, and the *genealogy* view, which uses the stored
ETAS parent links.  Comparing the two shows how much of a statistic is
selection.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from .catalog import Catalog
from .etas import DAYS_PER_YEAR


# ----------------------------------------------------------------------
# frequency-magnitude distribution
# ----------------------------------------------------------------------
def fmd(mags, mmin=None, dm=0.1):
    """Incremental and cumulative counts per magnitude bin (bin centres)."""
    mags = np.asarray(mags)
    if mmin is None:
        mmin = np.floor(mags.min() / dm) * dm
    edges = np.arange(mmin - dm / 2, mags.max() + dm, dm)
    inc, _ = np.histogram(mags, edges)
    cum = np.cumsum(inc[::-1])[::-1]
    centres = edges[:-1] + dm / 2
    return centres, inc, cum


def b_value_mle(mags, mc, mmax=None, dm=0.0):
    """Maximum-likelihood b-value for mc <= M (< mmax); returns (b, sigma_b, n).

    dm is the magnitude binning; synthetic magnitudes are continuous, so
    the default is 0 (Aki, 1965).  For binned data use dm = 0.1, which
    applies the Utsu (1966) half-bin correction.  With an upper bound the
    likelihood of the doubly truncated exponential is maximised (Page,
    1968); without it the estimate is the closed-form Aki-Utsu value.
    """
    mags = np.asarray(mags, float)
    lo = mc - dm / 2.0
    sel = mags >= lo
    if mmax is not None:
        sel &= mags < mmax
    m = mags[sel] - lo
    n = len(m)
    if n < 10:
        return np.nan, np.nan, n
    mbar = m.mean()
    if mmax is None:
        beta = 1.0 / mbar
    else:
        from scipy.optimize import brentq

        w = mmax - lo
        if w <= 0 or not (0.0 < mbar < 0.5 * w):
            # no interior maximum (empty window, or data not decreasing within it)
            return (0.0 if mbar >= 0.5 * w > 0 else np.nan), np.nan, n

        def h(x):
            # 1/x - 1/(e^x - 1), evaluated stably; decreases from 1/2 at x = 0 to 0
            x = float(x)
            if x < 1e-4:
                return 0.5 - x / 12.0 + x ** 3 / 720.0
            return 1.0 / x - 1.0 / np.expm1(x)

        def g(x):
            # score of the doubly truncated exponential in x = beta * w
            return w * h(x) - mbar

        try:
            x = brentq(g, 1e-10, 2.0 * w / mbar + 50.0)
        except ValueError:
            return np.nan, np.nan, n
        beta = x / w
    b = beta / np.log(10.0)
    return float(b), float(b / np.sqrt(n)), n


def b_value_regimes(cat: Catalog, dm=0.0):
    """b-values in the sub-patch, fault and width-saturated magnitude regimes."""
    ev = cat.events
    m_patch = cat.meta["m_patch"]
    m_sat = saturation_magnitude(cat)
    out = {}
    out["b_sub"] = (b_value_mle(ev["mag"], cat.params.mag_min, m_patch, dm)[0]
                    if m_patch > cat.params.mag_min + 0.2 else np.nan)
    out["b_fault_below_sat"] = b_value_mle(ev["mag"], m_patch, m_sat, dm)[0]
    out["b_fault_above_sat"] = b_value_mle(ev["mag"], m_sat, None, dm)[0]
    out["b_all"] = b_value_mle(ev["mag"], cat.params.mag_min, None, dm)[0]
    out["m_sat"] = m_sat
    return out


def saturation_magnitude(cat: Catalog, ratio=2.0):
    """Magnitude above which the median rupture aspect ratio exceeds ``ratio``."""
    fe = cat.fault_events()
    if len(fe) < 20:
        return np.nan
    ar = fe["length"] / fe["width"]
    edges = np.arange(np.floor(fe["mag"].min() * 10) / 10, fe["mag"].max() + 0.2, 0.2)
    for lo in edges:
        sel = (fe["mag"] >= lo) & (fe["mag"] < lo + 0.2)
        if sel.sum() >= 5 and np.median(ar[sel]) > ratio:
            return float(lo)
    return float(fe["mag"].max())


# ----------------------------------------------------------------------
# source scaling
# ----------------------------------------------------------------------
def scaling_fits(cat: Catalog):
    """Slopes of log10 M0 vs log10 L below and above width saturation and
    median realised stress-drop-equivalent quantities."""
    fe = cat.fault_events()
    out = {}
    if len(fe) < 20:
        return out
    L = fe["length"]
    W = fe["width"]
    M0 = fe["m0"]
    wf = cat.fault.width
    small = (L < 0.8 * wf) & (fe["npatch"] >= 4)
    large = L > 3.0 * wf
    for name, sel in (("slope_m0_L_small", small), ("slope_m0_L_large", large)):
        if sel.sum() >= 10:
            out[name] = float(np.polyfit(np.log10(L[sel]), np.log10(M0[sel]), 1)[0])
        else:
            out[name] = np.nan
    # mean slip / length, Blaser-style
    out["median_dmean_over_L"] = float(np.median(fe["dmean"] / L))
    out["median_dmean_over_L_large"] = float(np.median((fe["dmean"] / L)[large])) if large.sum() else np.nan
    out["frac_surface_rupturing_m6.5"] = float(np.mean(fe["surface"][fe["mag"] >= 6.5])) if (fe["mag"] >= 6.5).sum() else np.nan
    out["frac_surface_rupturing_m5.5_6.5"] = float(np.mean(fe["surface"][(fe["mag"] >= 5.5) & (fe["mag"] < 6.5)])) if ((fe["mag"] >= 5.5) & (fe["mag"] < 6.5)).sum() else np.nan
    out["max_mag"] = float(fe["mag"].max())
    out["frac_multi_segment"] = multi_segment_fraction(cat)
    return out


def multi_segment_fraction(cat: Catalog, mmin=6.0):
    """Fraction of stored events (M >= mmin) that span more than one section."""
    fe = cat.fault_events()
    sel = np.where((fe["has_slip"] == 1) & (fe["mag"] >= mmin))[0]
    if len(sel) == 0 or cat.fault.n_conn == 0:
        return np.nan
    ev_index = np.where(cat.events["kind"] == 1)[0][sel]
    multi = 0
    for e in ev_index:
        idx, _ = cat.slip_of(e)
        if len(np.unique(cat.fault.section[idx])) > 1:
            multi += 1
    return multi / len(sel)


# ----------------------------------------------------------------------
# space-time windows: mainshocks, aftershocks, foreshocks
# ----------------------------------------------------------------------
def _event_radius(ev, r_min, f_len):
    return np.maximum(r_min, f_len * ev["length"])


def select_mainshocks(cat: Catalog, mmain=6.0, t_window_yr=1.0, r_min_km=10.0, f_len=1.0,
                      t_start=0.0):
    """Indices of events with M >= mmain that have no larger event within
    +/- t_window and within max(r_min, f_len * L) of their centroid."""
    ev = cat.events
    t = ev["t"]
    cand = np.where((ev["mag"] >= mmain) & (t >= t_start))[0]
    R = _event_radius(ev, r_min_km * 1e3, f_len)
    keep = []
    for i in cand:
        lo = np.searchsorted(t, t[i] - t_window_yr)
        hi = np.searchsorted(t, t[i] + t_window_yr, side="right")
        j = np.arange(lo, hi)
        j = j[j != i]
        d = np.hypot(ev["cs"][j] - ev["cs"][i], ev["czd"][j] - ev["czd"][i])
        bigger = (ev["mag"][j] >= ev["mag"][i]) & (d <= np.maximum(R[i], R[j]))
        if not np.any(bigger):
            keep.append(i)
    return np.array(keep, int)


def aftershocks_of(cat: Catalog, i, t_window_yr=1.0, r_min_km=10.0, f_len=1.0, mc=None):
    """Indices of events in (t_i, t_i + window] within the mainshock radius."""
    ev = cat.events
    t = ev["t"]
    R = max(r_min_km * 1e3, f_len * ev["length"][i])
    lo = np.searchsorted(t, t[i], side="right")
    hi = np.searchsorted(t, t[i] + t_window_yr, side="right")
    j = np.arange(lo, hi)
    d = np.hypot(ev["cs"][j] - ev["cs"][i], ev["czd"][j] - ev["czd"][i])
    sel = d <= R
    if mc is not None:
        sel &= ev["mag"][j] >= mc
    return j[sel]


def foreshocks_of(cat: Catalog, i, t_window_yr=30.0 / DAYS_PER_YEAR, r_min_km=10.0, f_len=1.0,
                  mc=None):
    ev = cat.events
    t = ev["t"]
    R = max(r_min_km * 1e3, f_len * ev["length"][i])
    lo = np.searchsorted(t, t[i] - t_window_yr)
    hi = np.searchsorted(t, t[i])
    j = np.arange(lo, hi)
    d = np.hypot(ev["cs"][j] - ev["cs"][i], ev["czd"][j] - ev["czd"][i])
    sel = d <= R
    if mc is not None:
        sel &= ev["mag"][j] >= mc
    return j[sel]


# ----------------------------------------------------------------------
# Omori
# ----------------------------------------------------------------------
def omori_loglik(params, t, t0, t1):
    lnK, lnc, p = params
    K, c = np.exp(lnK), np.exp(lnc)
    if p <= 0.2 or p > 3.0:
        return 1e30
    ll = len(t) * lnK - p * np.sum(np.log(t + c))
    if abs(p - 1.0) < 1e-6:
        integ = K * (np.log(t1 + c) - np.log(t0 + c))
    else:
        integ = K / (1.0 - p) * ((t1 + c) ** (1.0 - p) - (t0 + c) ** (1.0 - p))
    return -(ll - integ)


def fit_omori(times_days, t0_days=0.0, t1_days=365.25):
    """Maximum-likelihood (K, c, p) for a stacked Omori sequence in [t0, t1] days."""
    t = np.asarray(times_days, float)
    t = t[(t > t0_days) & (t <= t1_days)]
    if len(t) < 20:
        return dict(K=np.nan, c=np.nan, p=np.nan, n=len(t))
    x0 = np.array([np.log(len(t) / 10.0), np.log(0.05), 1.1])
    best = None
    for p0 in (0.9, 1.1, 1.3):
        x0[2] = p0
        res = minimize(omori_loglik, x0, args=(t, t0_days, t1_days), method="Nelder-Mead",
                       options=dict(xatol=1e-4, fatol=1e-4, maxiter=4000))
        if best is None or res.fun < best.fun:
            best = res
    lnK, lnc, p = best.x
    return dict(K=float(np.exp(lnK)), c=float(np.exp(lnc)), p=float(p), n=len(t))


def omori_stack(cat: Catalog, mmain=6.0, t_window_yr=1.0, r_min_km=10.0, f_len=1.0, mc=None,
                t_start=0.0):
    """Stacked aftershock delays [days] relative to each mainshock and the fit."""
    mains = select_mainshocks(cat, mmain, t_window_yr, r_min_km, f_len, t_start)
    ev = cat.events
    delays = []
    for i in mains:
        j = aftershocks_of(cat, i, t_window_yr, r_min_km, f_len, mc)
        delays.append((ev["t"][j] - ev["t"][i]) * DAYS_PER_YEAR)
    delays = np.concatenate(delays) if delays else np.zeros(0)
    fit = fit_omori(delays, 0.0, t_window_yr * DAYS_PER_YEAR)
    fit["n_mainshocks"] = len(mains)
    fit["per_mainshock"] = len(delays) / max(1, len(mains))
    return delays, fit


def omori_binned(delays_days, nbins=25, tmin=1e-3, tmax=365.25):
    """Rate per day in log-spaced bins (for plotting)."""
    edges = np.logspace(np.log10(tmin), np.log10(tmax), nbins + 1)
    cnt, _ = np.histogram(delays_days, edges)
    rate = cnt / np.diff(edges)
    centres = np.sqrt(edges[:-1] * edges[1:])
    return centres, rate


# ----------------------------------------------------------------------
# Bath's law
# ----------------------------------------------------------------------
def bath_gaps(cat: Catalog, mmain=6.0, t_window_yr=1.0, r_min_km=10.0, f_len=1.0, t_start=0.0):
    """Observer view: M_main - max aftershock magnitude for each mainshock.

    Returns (gaps, n_without_aftershock).  Mainshocks without any
    aftershock above completeness are excluded from ``gaps``.
    """
    mains = select_mainshocks(cat, mmain, t_window_yr, r_min_km, f_len, t_start)
    ev = cat.events
    gaps = []
    none = 0
    for i in mains:
        j = aftershocks_of(cat, i, t_window_yr, r_min_km, f_len)
        if len(j) == 0:
            none += 1
            continue
        gaps.append(ev["mag"][i] - ev["mag"][j].max())
    return np.array(gaps), none


def bath_gaps_genealogy(cat: Catalog, mmain=6.0):
    """Genealogy view: gap to the largest descendant (all generations)."""
    ev = cat.events
    parent = ev["parent"]
    n = len(ev)
    # root ancestor of each event
    root = np.arange(n)
    order = np.arange(n)
    for i in order:
        if parent[i] >= 0:
            root[i] = root[parent[i]]
    mains = np.where((ev["mag"] >= mmain) & (parent < 0))[0]
    gaps = []
    for i in mains:
        desc = np.where((root == i) & (np.arange(n) != i))[0]
        if len(desc):
            gaps.append(ev["mag"][i] - ev["mag"][desc].max())
    return np.array(gaps)


# ----------------------------------------------------------------------
# foreshocks
# ----------------------------------------------------------------------
def foreshock_fraction(cat: Catalog, mmain=6.0, t_fore_days=10.0, r_km=10.0, mc_fore=None,
                       t_window_yr=1.0, n_random=5, seed=0, t_start=0.0, **win):
    """Fraction of window-selected mainshocks preceded by a foreshock.

    A foreshock is an event with M >= mc_fore (default mag_min) within
    ``t_fore_days`` before the mainshock and within ``r_km`` of the
    mainshock epicentre (its nucleation patch), measured between
    epicentres.  The chance level repeats the test at ``n_random`` random
    times per mainshock at the same epicentre.  This replaced (3 November
    2026) a definition with centroid distances and a rupture-length
    radius, whose chance level approached one for large ruptures.
    """
    ev = cat.events
    f = cat.fault
    if mc_fore is None:
        mc_fore = cat.params.mag_min
    mains = select_mainshocks(cat, mmain, t_window_yr, t_start=t_start, **win)
    if len(mains) == 0:
        return dict(fraction=np.nan, chance=np.nan, n_mainshocks=0)
    t = ev["t"]
    es = f.s[ev["nuc"]]
    ez = f.zd[ev["nuc"]]
    dt = t_fore_days / DAYS_PER_YEAR
    r = r_km * 1e3

    def any_before(t0, s0, z0):
        lo = np.searchsorted(t, t0 - dt)
        hi = np.searchsorted(t, t0)
        j = np.arange(lo, hi)
        j = j[ev["mag"][j] >= mc_fore]
        return bool(np.any(np.hypot(es[j] - s0, ez[j] - z0) <= r))

    hit = np.mean([any_before(t[i], es[i], ez[i]) for i in mains])
    rng = np.random.default_rng(seed)
    t_lo = max(t_start, dt)
    trials = [any_before(rng.uniform(t_lo, cat.duration), es[i], ez[i])
              for i in mains for _ in range(n_random)]
    return dict(fraction=float(hit), chance=float(np.mean(trials)), n_mainshocks=len(mains))


# ----------------------------------------------------------------------
# recurrence
# ----------------------------------------------------------------------
def interevent_cv(cat: Catalog, mmin=6.5, t_start=0.0, declustered=True, **win):
    """Coefficient of variation of interevent times of events with M >= mmin.

    With ``declustered`` only mainshocks (window method) are used.
    """
    ev = cat.events
    if declustered:
        idx = select_mainshocks(cat, mmin, t_start=t_start, **win)
    else:
        idx = np.where((ev["mag"] >= mmin) & (ev["t"] >= t_start))[0]
    t = np.sort(ev["t"][idx])
    if len(t) < 5:
        return dict(cv=np.nan, mean_yr=np.nan, n=len(t))
    dt = np.diff(t)
    return dict(cv=float(dt.std() / dt.mean()), mean_yr=float(dt.mean()), n=len(t))


def site_recurrence_cv(cat: Catalog, n_sites=9, slip_threshold=0.5, t_start=0.0):
    """Mean and spread of the point-recurrence CV over sites spaced along
    strike at mid depth (synthetic paleoseismic trenches)."""
    f = cat.fault
    cvs, trs = [], []
    for frac in np.linspace(0.1, 0.9, n_sites):
        site = int(np.argmin((f.col - int(frac * f.ns)) ** 2 + (f.row - f.nz // 2) ** 2))
        pr = point_recurrence(cat, site, slip_threshold, t_start)
        if np.isfinite(pr["cv"]):
            cvs.append(pr["cv"])
            trs.append(pr["mean_yr"])
    if not cvs:
        return dict(cv=np.nan, cv_std=np.nan, mean_yr=np.nan, n_sites=0)
    return dict(cv=float(np.mean(cvs)), cv_std=float(np.std(cvs)), mean_yr=float(np.mean(trs)),
                n_sites=len(cvs))


def point_recurrence(cat: Catalog, patch, slip_threshold=0.5, t_start=0.0):
    """Times at which a given patch slipped at least ``slip_threshold`` [m]
    in a stored event (a synthetic paleoseismic site), with the CV."""
    times = []
    for k, e in enumerate(cat.slip_event):
        a, b = cat.slip_ptr[k], cat.slip_ptr[k + 1]
        idx = cat.slip_idx[a:b]
        m = idx == patch
        if np.any(m) and cat.slip_val[a:b][m][0] >= slip_threshold and cat.events["t"][e] >= t_start:
            times.append(cat.events["t"][e])
    t = np.array(times)
    if len(t) < 3:
        return dict(times=t, cv=np.nan, mean_yr=np.nan)
    dt = np.diff(np.sort(t))
    return dict(times=t, cv=float(dt.std() / dt.mean()), mean_yr=float(dt.mean()))


# ----------------------------------------------------------------------
# budget
# ----------------------------------------------------------------------
def moment_closure(cat: Catalog, t_start=0.0):
    ev = cat.events
    sel = ev["t"] >= t_start
    loaded = cat.fault.moment_rate * (cat.duration - t_start)
    rel = ev["m0"][sel].sum()
    sub = ev["m0"][sel & (ev["kind"] == 0)].sum()
    return dict(closure=rel / loaded, sub_patch_fraction=sub / max(rel, 1e-30),
                released=rel, loaded=loaded)


def cumulative_moment(cat: Catalog, nt=400):
    ev = cat.events
    t = np.linspace(0, cat.duration, nt)
    cum = np.interp(t, ev["t"], np.cumsum(ev["m0"]), left=0.0)
    return t, cum, cat.fault.moment_rate * t


def slip_rate_recovered(cat: Catalog):
    """Released slip per patch divided by duration, versus prescribed rate."""
    return cat.released / cat.duration, cat.fault.slip_rate


# ----------------------------------------------------------------------
# summary
# ----------------------------------------------------------------------
def summary_metrics(cat: Catalog, mmain=6.0, t_start=None) -> dict:
    """Scalar diagnostics for one catalog (used by sweeps)."""
    if t_start is None:
        t_start = 0.1 * cat.duration
    ev = cat.events
    out = dict(cat.meta)
    out.update(b_value_regimes(cat))
    out.update(scaling_fits(cat))
    delays, fit = omori_stack(cat, mmain=mmain, t_start=t_start)
    out.update({"omori_" + k: v for k, v in fit.items()})
    gaps, none = bath_gaps(cat, mmain=mmain, t_start=t_start)
    out["bath_mean"] = float(gaps.mean()) if len(gaps) else np.nan
    out["bath_median"] = float(np.median(gaps)) if len(gaps) else np.nan
    out["bath_n"] = len(gaps)
    out["bath_no_aftershock"] = none
    fs = foreshock_fraction(cat, mmain=mmain, t_start=t_start)
    out["foreshock_fraction"] = fs["fraction"]
    out["foreshock_chance"] = fs["chance"]
    fs = foreshock_fraction(cat, mmain=mmain, t_start=t_start, t_fore_days=3.0, r_km=5.0)
    out["foreshock_fraction_3d_5km"] = fs["fraction"]
    out["foreshock_chance_3d_5km"] = fs["chance"]
    cv = interevent_cv(cat, mmin=6.5, t_start=t_start)
    out["cv_m6.5_declustered"] = cv["cv"]
    out["recurrence_m6.5_yr"] = cv["mean_yr"]
    cv = interevent_cv(cat, mmin=6.5, t_start=t_start, declustered=False)
    out["cv_m6.5_all"] = cv["cv"]
    cv = interevent_cv(cat, mmin=7.0, t_start=t_start)
    out["cv_m7_declustered"] = cv["cv"]
    for thr in (0.3, 1.0):
        sr = site_recurrence_cv(cat, slip_threshold=thr, t_start=t_start)
        out[f"site_cv_slip{thr}"] = sr["cv"]
        out[f"site_cv_std_slip{thr}"] = sr["cv_std"]
        out[f"site_recurrence_yr_slip{thr}"] = sr["mean_yr"]
    mc = moment_closure(cat, t_start)
    out["closure_after_burnin"] = mc["closure"]
    out["sub_patch_moment_fraction"] = mc["sub_patch_fraction"]
    out["n_m6"] = int(np.sum(ev["mag"] >= 6.0))
    out["n_m7"] = int(np.sum(ev["mag"] >= 7.0))
    if getattr(cat.params, "sites_file", ""):
        from .sites import co_rupture, compare, load_sites

        sites = load_sites(cat.params.sites_file)
        for th in cat.params.site_thresholds_m:
            for r in compare(cat, sites, threshold=th, t_start=t_start):
                for k in ("n", "ri", "ri_z", "cv", "sed"):
                    out[f"site_{r['name']}_{k}_th{th:g}"] = r[k]
        names, M = co_rupture(cat, sites, threshold=cat.params.site_thresholds_m[0], t_start=t_start)
        for i, a in enumerate(names):
            for j, b in enumerate(names):
                if i < j:
                    out[f"site_corupture_{a}_{b}"] = float(M[i, j])
    fe = cat.fault_events()
    if len(fe):
        out["frac_stop_hazard"] = float(np.mean(fe["stop"] == 0))
        out["frac_stop_deficit"] = float(np.mean(fe["stop"] == 1))
        out["frac_stop_boundary"] = float(np.mean(fe["stop"] == 2))
        big = fe["mag"] >= 6.0
        out["frac_stop_deficit_m6"] = float(np.mean(fe["stop"][big] == 1)) if big.sum() else np.nan
    return out
