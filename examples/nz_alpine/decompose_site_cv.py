"""Decompose site recurrence variability for a saved catalog (see examples/nz_alpine/README.md, "What controls site CV").

    python3 examples/nz_alpine/decompose_site_cv.py catalog.npz [threshold_m]

Requires a catalog saved with store_slip_min_mag low enough (6.0) and data/nz_alpine/paleo_sites.csv.
"""
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from stochseis import Catalog
from stochseis.sites import load_sites, site_patch
cat = Catalog.load(sys.argv[1]); ev = cat.events; f = cat.fault
th = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
T0 = 0.15 * cat.duration
sites = load_sites(os.path.join(os.path.dirname(__file__), "..", "..", "data", "nz_alpine", "paleo_sites.csv"))
# released slip history per site patch from stored events (M >= store_slip_min_mag) plus sub-threshold info
for st in sites:
    k, _ = site_patch(f, st["x_nztm_m"], st["y_nztm_m"])
    v = f.slip_rate[k]
    ts, sl, mg, dist = [], [], [], []
    for j, e in enumerate(cat.slip_event):
        a, b = cat.slip_ptr[j], cat.slip_ptr[j + 1]
        idx = cat.slip_idx[a:b]
        m = idx == k
        if m.any():
            ts.append(ev["t"][e]); sl.append(cat.slip_val[a:b][m][0]); mg.append(ev["mag"][e])
            dist.append(abs(f.s[ev["nuc"][e]] - f.s[k]) / 1e3)
    ts, sl, mg, dist = map(np.array, (ts, sl, mg, dist))
    o = np.argsort(ts); ts, sl, mg, dist = ts[o], sl[o], mg[o], dist[o]
    # deficit at the site just before each slip episode: v t + d0 - cumulative slip before
    cum = np.concatenate([[0.0], np.cumsum(sl)])[:-1]
    dbefore = v * ts + cat.deficit0[k] - cum   # ignores slip from unstored (M<6) events
    rec = (sl >= th) & (ts >= T0)
    t = ts[rec]; s = sl[rec]; db = dbefore[rec]; m_ = mg[rec]; dd = dist[rec]
    T = np.diff(t)
    frac_rec = sl[ts >= T0][sl[ts >= T0] >= th].sum() / max(sl[ts >= T0].sum(), 1e-9)
    print(f"{st['name']:8s} v {v*1e3:4.1f} mm/yr  n {len(t)}  RI {T.mean():5.0f}  CV(T) {T.std(ddof=1)/T.mean():.2f}  "
          f"SED {s.mean():.1f} CV(SED) {s.std()/s.mean():.2f}  deficit before {db.mean():.1f} CV(D_before) {db.std()/db.mean():.2f}  "
          f"slip/deficit {np.mean(s/np.maximum(db,1e-3)):.2f}")
    print(f"         corr(SED_i, T_next) {np.corrcoef(s[:-1], T)[0,1]:+.2f} (time-predictable)  corr(T_prev, SED_i) {np.corrcoef(T, s[1:])[0,1]:+.2f} (slip-predictable)"
          f"  fraction of site slip in recorded events {frac_rec:.2f}  sub-threshold slipping events {np.sum((sl<th)&(ts>=T0))}"
          f"  recorded events: median M {np.median(m_):.2f}, nucleated within 30 km {np.mean(dd<30):.2f}")
