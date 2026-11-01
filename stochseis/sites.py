"""Synthetic paleoseismic records at sites and comparison with data.

A site is a map position (x, y) close to the fault trace.  It is assigned
the top-row patch nearest in map view; an event is recorded at the site
if its slip on that patch is at least ``threshold`` metres (a detection
threshold for surface rupture), and its single-event displacement (SED)
is that slip.  Only events whose slip was stored (``store_slip_min_mag``)
can be recorded, so the storage magnitude must be low enough that every
event with slip above the threshold is stored; ``site_records`` checks
this and warns.

Observed recurrence data are summarised by a mean interval and a standard
deviation of intervals; ``compare`` returns the synthetic counterparts
and, for the mean interval, a z-score using the standard error of the
observed mean (sigma / sqrt(n - 1) for n events, n - 1 intervals).
"""

from __future__ import annotations

import csv
import warnings

import numpy as np

from .catalog import Catalog


def load_sites(path):
    """Read a site table (CSV with name, x_nztm_m, y_nztm_m, ri_yr, ri_sigma_yr, n_ri, sed_* columns)."""
    out = []
    with open(path) as f:
        for row in csv.DictReader(f):
            d = {}
            for k, v in row.items():
                try:
                    d[k] = float(v)
                except (TypeError, ValueError):
                    d[k] = v
            out.append(d)
    return out


def site_patch(fault, x, y, top_rows=1):
    """Index of the patch in the top ``top_rows`` rows nearest to (x, y) in map view, and its distance."""
    top = np.where(fault.row < top_rows)[0]
    d = np.hypot(fault.x[top] - x, fault.y[top] - y)
    k = np.argmin(d)
    return int(top[k]), float(d[k])


def site_records(cat: Catalog, sites, threshold=1.0, t_start=0.0):
    """Event times and SEDs at each site.  Returns {name: dict(times, sed, patch, dist)}."""
    f = cat.fault
    ev = cat.events
    patches = {}
    for s in sites:
        p, dist = site_patch(f, s["x_nztm_m"], s["y_nztm_m"])
        patches[s["name"]] = (p, dist)
    lookup = {}
    for name, (p, _) in patches.items():
        lookup.setdefault(p, []).append(name)
    rec = {name: ([], []) for name in patches}
    for k, e in enumerate(cat.slip_event):
        if ev["t"][e] < t_start:
            continue
        a, b = cat.slip_ptr[k], cat.slip_ptr[k + 1]
        idx = cat.slip_idx[a:b]
        hit = np.isin(idx, list(lookup))
        if not np.any(hit):
            continue
        for j in np.where(hit)[0]:
            val = cat.slip_val[a + j]
            if val >= threshold:
                for name in lookup[idx[j]]:
                    rec[name][0].append(ev["t"][e])
                    rec[name][1].append(val)
    # storage check: the smallest stored magnitude must be below what a threshold slip needs
    m_store = cat.params.store_slip_min_mag
    mu = float(np.mean(f.mu))
    m0_needed = mu * threshold * (f.dx * f.dx)        # one patch carrying the threshold slip
    m_needed = 2.0 / 3.0 * (np.log10(m0_needed) - 9.05)
    if m_store > m_needed + 1.0:
        warnings.warn(f"store_slip_min_mag = {m_store} may miss events with {threshold} m of slip "
                      f"at a site (smallest such events are around M{m_needed + 1.0:.1f})")
    out = {}
    for name, (times, sed) in rec.items():
        o = np.argsort(times)
        out[name] = dict(times=np.asarray(times)[o], sed=np.asarray(sed)[o],
                         patch=patches[name][0], dist=patches[name][1])
    return out


def record_stats(times, sed):
    t = np.sort(np.asarray(times))
    n = len(t)
    if n < 3:
        return dict(n=n, ri=np.nan, ri_sd=np.nan, cv=np.nan, sed_mean=np.nan, sed_sd=np.nan)
    dt = np.diff(t)
    return dict(n=n, ri=float(dt.mean()), ri_sd=float(dt.std(ddof=1)), cv=float(dt.std(ddof=1) / dt.mean()),
                sed_mean=float(np.mean(sed)), sed_sd=float(np.std(sed)))


def compare(cat: Catalog, sites, threshold=1.0, t_start=None):
    """Synthetic versus observed recurrence and SED at each site."""
    if t_start is None:
        t_start = 0.1 * cat.duration
    rec = site_records(cat, sites, threshold, t_start)
    rows = []
    for s in sites:
        r = rec[s["name"]]
        st = record_stats(r["times"], r["sed"])
        n_obs = s.get("n_ri", np.nan)
        se = s["ri_sigma_yr"] / np.sqrt(max(n_obs - 1, 1)) if np.isfinite(n_obs) else s["ri_sigma_yr"]
        rows.append(dict(name=s["name"], patch=r["patch"], dist_km=r["dist"] / 1e3,
                         s_km=float(cat.fault.s[r["patch"]] / 1e3), n=st["n"],
                         ri=st["ri"], ri_obs=s["ri_yr"], ri_z=(st["ri"] - s["ri_yr"]) / se,
                         cv=st["cv"], cv_obs=s["ri_sigma_yr"] / s["ri_yr"],
                         sed=st["sed_mean"], sed_obs=s["sed_mean_m"],
                         sed_in_range=bool(s["sed_min_m"] <= st["sed_mean"] <= s["sed_max_m"])
                         if np.isfinite(st["sed_mean"]) else False))
    return rows


def co_rupture(cat: Catalog, sites, threshold=1.0, t_start=0.0):
    """Fraction of recorded events at each site that also rupture each other site (matrix)."""
    rec = site_records(cat, sites, threshold, t_start)
    names = [s["name"] for s in sites]
    m = np.full((len(names), len(names)), np.nan)
    for i, a in enumerate(names):
        ta = set(np.round(rec[a]["times"], 9))
        for j, b in enumerate(names):
            tb = set(np.round(rec[b]["times"], 9))
            if ta:
                m[i, j] = len(ta & tb) / len(ta)
    return names, m
