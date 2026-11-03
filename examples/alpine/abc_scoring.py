"""Scoring of ABC sweep summaries against the Alpine Fault paleoseismic data.

Shared by abc_accept.py, make_posterior_spec.py and campaign_summary.py.

For each design vector (averaged over its seeds) and a detection
threshold, the misfit is the mean over observables of squared
standardized residuals:

  per site:  RI   with the standard error of the observed mean (sigma / sqrt(n - 1))
             CV   with 0.15 (sampling error of a CV from about 10 events)
             SED  with half the observed SED range
  optional (seismicity=True):
             b-value of grown ruptures below width saturation, 1.0 +- 0.1
             mean Bath gap, 1.2 +- 0.15

A missing value (fewer than three recorded events) counts as a residual
of 5 standard errors.

site_thresholds ({site: threshold_m}) overrides the detection threshold
for individual sites, e.g. {"north": 1.0}: the north record includes
events of about 1 m (observed SED range 1 to 9.4 m). Only thresholds that
the sweep summarized (site_thresholds_m, by default 1, 2, 4 m) are
available.
"""
import csv
import json
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SITES_FILE = os.path.join(HERE, "..", "..", "data", "nz_alpine", "paleo_sites.csv")
PARAMS = ["stress_drop_mpa", "stress_drop_ln_sigma", "nucleation_exponent", "branching_ratio",
          "deep_slip_km", "condition_weight", "arrest_scale", "cap_by_deficit", "b_value",
          "segment_boundary_strength", "budget_scale", "slip_fraction", "slip_cv", "admission_exponent",
          "hazard_max_area_km2", "state_smoothing_km"]
SEISMICITY = (("b_fault_below_sat", 1.0, 0.1), ("bath_mean", 1.2, 0.15))
MISSING = 25.0


def _convert(v):
    if v in ("", None):
        return None
    if v in ("True", "False"):
        return v == "True"
    try:
        return float(v)
    except ValueError:
        return v


def load_summary(path):
    """Rows of a sweep summary, from the JSON or the CSV written by the sweep tool."""
    if path.endswith(".csv"):
        with open(path) as fh:
            return [{k: _convert(v) for k, v in r.items()} for r in csv.DictReader(fh)]
    with open(path) as fh:
        return json.load(fh)


def load_sites(path=SITES_FILE):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def _mean(rows, key):
    v = [r.get(key) for r in rows]
    v = [float(x) for x in v if x is not None and x != ""]
    return float(np.nanmean(v)) if v else np.nan


def parse_site_thresholds(args):
    """['north=1', ...] -> {'north': 1.0}"""
    out = {}
    for a in args:
        k, v = a.split("=")
        out[k] = float(v)
    return out


def score(rows, threshold=2.0, seismicity=False, sites=None, site_thresholds=None):
    """Misfit per design vector, sorted best first.

    Returns a list of dicts with keys misfit, key (design index or file), params (the varied and
    fixed p_ values present), detail ({site: (ri, cv, sed)}, plus "seismicity": (b, bath)), rows.
    """
    sites = sites or load_sites()
    groups = defaultdict(list)
    for r in rows:
        if "error" in r and r["error"]:
            continue
        groups[r.get("design_index", r.get("file"))].append(r)
    out = []
    for key, rs in groups.items():
        d2 = 0.0
        nobs = 0
        detail = {}
        for s in sites:
            n = s["name"]
            th = (site_thresholds or {}).get(n, threshold)
            ri_obs, ri_sig, n_ri = float(s["ri_yr"]), float(s["ri_sigma_yr"]), float(s["n_ri"])
            sed_obs = float(s["sed_mean_m"])
            sed_sig = 0.5 * (float(s["sed_max_m"]) - float(s["sed_min_m"]))
            ri = _mean(rs, f"site_{n}_ri_th{th:g}")
            cv = _mean(rs, f"site_{n}_cv_th{th:g}")
            sed = _mean(rs, f"site_{n}_sed_th{th:g}")
            for val, obs, sig in ((ri, ri_obs, ri_sig / np.sqrt(max(n_ri - 1, 1))), (cv, ri_sig / ri_obs, 0.15),
                                  (sed, sed_obs, sed_sig)):
                d2 += ((val - obs) / sig) ** 2 if np.isfinite(val) else MISSING
                nobs += 1
            detail[n] = (ri, cv, sed)
        if seismicity:
            vals = []
            for k, obs, sig in SEISMICITY:
                val = _mean(rs, k)
                d2 += ((val - obs) / sig) ** 2 if np.isfinite(val) else MISSING
                nobs += 1
                vals.append(val)
            detail["seismicity"] = tuple(vals)
        params = {p: rs[0]["p_" + p] for p in PARAMS if "p_" + p in rs[0]}
        out.append(dict(misfit=d2 / nobs, key=key, params=params, detail=detail, rows=rs))
    out.sort(key=lambda x: x["misfit"])
    return out


def accepted(scored, tolerance=2.0):
    """Design vectors with mean squared standardized misfit below tolerance^2."""
    return [x for x in scored if x["misfit"] < tolerance ** 2]
