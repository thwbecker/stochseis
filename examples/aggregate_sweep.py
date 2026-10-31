"""Aggregate a sweep summary over seeds: mean and standard deviation per
parameter combination for a set of metrics.

    python examples/aggregate_sweep.py sweep_out/sweep_summary.json \
        p_arrest_scale p_stress_drop_ln_sigma > table.txt
"""
import json
import sys
import warnings
from collections import defaultdict

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)   # all-nan metrics in short runs

METRICS = ["max_mag", "n_m7", "b_fault_below_sat", "b_fault_above_sat", "bath_mean",
           "omori_p", "site_cv_slip0.3", "site_cv_slip1.0", "cv_m7_declustered",
           "slope_m0_L_large", "frac_stop_deficit_m6", "final_mean_deficit",
           "closure_after_burnin"]

rows = json.load(open(sys.argv[1]))
keys = sys.argv[2:]
groups = defaultdict(list)
for r in rows:
    groups[tuple(r[k] for k in keys)].append(r)
print("  ".join(f"{k:>12s}" for k in keys) + "  n  " + "  ".join(f"{m:>22s}" for m in METRICS))
for g in sorted(groups):
    rs = groups[g]
    cells = []
    for m in METRICS:
        v = np.array([np.nan if r.get(m) is None else r[m] for r in rs], float)
        cells.append(f"{np.nanmean(v):10.3g} +- {np.nanstd(v):9.2g}")
    print("  ".join(f"{x:>12}" for x in g) + f"  {len(rs)}  " + "  ".join(cells))
