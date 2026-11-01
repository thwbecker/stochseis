"""Accept parameter vectors whose synthetic site records match the data.

    python3 examples/alpine/abc_accept.py abc1_planar/sweep_summary.json [threshold_m] [tolerance]

For each design vector (averaged over its seeds) and the chosen detection
threshold (default 2 m), the distance to the data is

    sum over sites of ((RI - RI_obs) / sigma_RI)^2 + ((CV - CV_obs) / 0.15)^2
                      + ((SED - SED_obs) / sigma_SED)^2

with sigma_RI the standard error of the observed mean interval, 0.15 a
typical sampling error of a CV from about 10 events, and sigma_SED half
the observed SED range. Vectors with distance below ``tolerance`` per
observable (default 2, i.e. about two standard errors on average) are
accepted. The script prints the best vectors and the accepted parameter
ranges.
"""
import csv
import json
import os
import sys
from collections import defaultdict

import warnings

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

here = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(sys.argv[1]))
th = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
tol = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
sites = list(csv.DictReader(open(os.path.join(here, "..", "..", "data", "nz_alpine", "paleo_sites.csv"))))

params = ["stress_drop_mpa", "stress_drop_ln_sigma", "nucleation_exponent", "branching_ratio",
          "deep_slip_km", "condition_weight", "arrest_scale", "cap_by_deficit"]
groups = defaultdict(list)
for r in rows:
    groups[r.get("design_index", r.get("file"))].append(r)

scored = []
for key, rs in groups.items():
    d2 = 0.0
    nobs = 0
    detail = {}
    for s in sites:
        n = s["name"]
        ri_obs, ri_sig, n_ri = float(s["ri_yr"]), float(s["ri_sigma_yr"]), float(s["n_ri"])
        sed_obs = float(s["sed_mean_m"])
        sed_sig = 0.5 * (float(s["sed_max_m"]) - float(s["sed_min_m"]))
        ri = np.nanmean([np.nan if r.get(f"site_{n}_ri_th{th:g}") is None else r[f"site_{n}_ri_th{th:g}"] for r in rs])
        cv = np.nanmean([np.nan if r.get(f"site_{n}_cv_th{th:g}") is None else r[f"site_{n}_cv_th{th:g}"] for r in rs])
        sed = np.nanmean([np.nan if r.get(f"site_{n}_sed_th{th:g}") is None else r[f"site_{n}_sed_th{th:g}"] for r in rs])
        for val, obs, sig in ((ri, ri_obs, ri_sig / np.sqrt(max(n_ri - 1, 1))), (cv, ri_sig / ri_obs, 0.15),
                              (sed, sed_obs, sed_sig)):
            d2 += ((val - obs) / sig) ** 2 if np.isfinite(val) else 25.0
            nobs += 1
        detail[n] = (ri, cv, sed)
    scored.append((d2 / nobs, key, {p: rs[0].get("p_" + p) for p in params}, detail))
scored.sort(key=lambda x: x[0])

print(f"{len(scored)} design vectors, threshold {th} m; best ten (mean squared standardized misfit):")
for d, key, pv, det in scored[:10]:
    print(f"  {d:6.2f}  " + " ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}" for k, v in pv.items()))
    print("          " + "  ".join(f"{n}: RI {v[0]:.0f} CV {v[1]:.2f} SED {v[2]:.1f}" for n, v in det.items()))
acc = [x for x in scored if x[0] < tol ** 2]
print(f"\naccepted (mean squared misfit < {tol ** 2:g}): {len(acc)}")
if acc:
    for p in params:
        v = np.array([float(x[2][p]) for x in acc if x[2][p] is not None])
        if len(v):
            print(f"  {p:22s} median {np.median(v):7.3g}  range {v.min():7.3g} to {v.max():7.3g}")
