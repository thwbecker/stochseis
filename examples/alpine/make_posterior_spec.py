"""Write a sweep spec that reruns the best accepted ABC vectors with catalogs kept.

    python3 examples/alpine/make_posterior_spec.py abc4_planar/sweep_summary.json base_spec.json out_spec.json \
        [threshold_m=2] [n_best=20] [n_seeds=5] [duration_yr=20000]

The output spec has design "list": every selected vector is run for every
seed; run it without --no-catalogs so that event catalogs are available
for the intercomparison statistics.
"""
import csv
import json
import os
import sys
from collections import defaultdict

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
summary, base_spec, out_spec = sys.argv[1:4]
th = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
n_best = int(sys.argv[5]) if len(sys.argv) > 5 else 20
n_seeds = int(sys.argv[6]) if len(sys.argv) > 6 else 5
duration = float(sys.argv[7]) if len(sys.argv) > 7 else 20000.0
if summary.endswith(".csv"):
    rows = []
    for r in csv.DictReader(open(summary)):
        d = {}
        for k, v in r.items():
            if v in ("", None):
                d[k] = None
            elif v in ("True", "False"):
                d[k] = v == "True"
            else:
                try:
                    d[k] = float(v)
                except ValueError:
                    d[k] = v
        rows.append(d)
else:
    rows = json.load(open(summary))
sites = list(csv.DictReader(open(os.path.join(here, "..", "..", "data", "nz_alpine", "paleo_sites.csv"))))
spec = json.load(open(base_spec))
varied = list(spec.get("ranges", {})) + list(spec.get("log_ranges", {})) + list(spec.get("choices", {}))
groups = defaultdict(list)
for r in rows:
    groups[r["design_index"]].append(r)
scored = []
for k, rs in groups.items():
    d2, n = 0.0, 0
    for s in sites:
        nm = s["name"]
        obs = ((float(s["ri_yr"]), float(s["ri_sigma_yr"]) / np.sqrt(float(s["n_ri"]) - 1), "ri"),
               (float(s["ri_sigma_yr"]) / float(s["ri_yr"]), 0.15, "cv"),
               (float(s["sed_mean_m"]), 0.5 * (float(s["sed_max_m"]) - float(s["sed_min_m"])), "sed"))
        for o, sig, q in obs:
            v = np.nanmean([np.nan if x.get(f"site_{nm}_{q}_th{th:g}") is None else x[f"site_{nm}_{q}_th{th:g}"] for x in rs])
            d2 += ((v - o) / sig) ** 2 if np.isfinite(v) else 25.0
            n += 1
    scored.append((d2 / n, {p: rs[0]["p_" + p] for p in varied if "p_" + p in rs[0]}))
scored.sort(key=lambda x: x[0])
out = {"_comment": f"Posterior ensemble: best {n_best} vectors of {summary} (threshold {th} m), {n_seeds} seeds, catalogs kept.",
       "design": "list", "base": dict(spec["base"]), "vectors": [v for _, v in scored[:n_best]],
       "misfits": [m for m, _ in scored[:n_best]], "seeds": list(range(101, 101 + n_seeds)), "mmain": spec.get("mmain", 6.0)}
out["base"]["duration_yr"] = duration
json.dump(out, open(out_spec, "w"), indent=1)
print(f"wrote {out_spec}: {len(out['vectors'])} vectors x {n_seeds} seeds, misfits {out['misfits'][0]:.2f} to {out['misfits'][-1]:.2f}")
