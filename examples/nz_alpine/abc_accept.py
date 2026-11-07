# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Accept parameter vectors whose synthetic site records match the data.

    python3 examples/nz_alpine/abc_accept.py summary.(json|csv) [threshold_m] [tolerance] [--seismicity]
        [--site north=1 ...] [--bath-sigma 0.15] [--bath-m7]

Scores every design vector of a sweep summary (averaged over its seeds)
against recurrence interval, CV and single-event displacement at the three
paleoseismic sites for a detection threshold (default 2 m), optionally with
the b-value and Bath gap (see abc_scoring.py for the misfit). Vectors with a
mean squared standardized misfit below tolerance^2 (default tolerance 2) are
accepted. Prints the best ten vectors and the accepted parameter ranges.
Without --seismicity a model can match the site records with unphysical
seismicity (round 3: b about 0, Bath 0.7).
"""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from abc_scoring import accepted, load_summary, parse_site_thresholds, score  # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)


def main():
    seis = "--seismicity" in sys.argv
    argv = [a for a in sys.argv if a != "--seismicity"]
    site_th = {}
    while "--site" in argv:
        k = argv.index("--site")
        site_th.update(parse_site_thresholds([argv[k + 1]]))
        argv = argv[:k] + argv[k + 2:]
    bath_sigma = 0.15
    if "--bath-sigma" in argv:
        k = argv.index("--bath-sigma")
        bath_sigma = float(argv[k + 1])
        argv = argv[:k] + argv[k + 2:]
    bath_key = "bath_mean_m7" if "--bath-m7" in argv else "bath_mean"
    argv = [a for a in argv if a != "--bath-m7"]
    rows = load_summary(argv[1])
    th = float(argv[2]) if len(argv) > 2 else 2.0
    tol = float(argv[3]) if len(argv) > 3 else 2.0
    scored = score(rows, th, seis, site_thresholds=site_th, bath_sigma=bath_sigma, bath_key=bath_key)
    extra = "".join(f", {k} at {v:g} m" for k, v in site_th.items())
    print(f"{len(scored)} design vectors, threshold {th} m{extra}{', with seismicity' if seis else ''}; "
          "best ten (mean squared standardized misfit):")
    for x in scored[:10]:
        pv = x["params"]
        print(f"  {x['misfit']:6.2f}  " + " ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}"
                                                 for k, v in pv.items()))
        print("          " + "  ".join(f"{n}: RI {v[0]:.0f} CV {v[1]:.2f} SED {v[2]:.1f}" if n != "seismicity"
                                        else f"b {v[0]:.2f} Bath {v[1]:.2f}" for n, v in x["detail"].items()))
    acc = accepted(scored, tol)
    print(f"\naccepted (mean squared misfit < {tol ** 2:g}): {len(acc)}")
    if acc:
        for p in acc[0]["params"]:
            v = np.array([float(x["params"][p]) for x in acc if x["params"].get(p) is not None])
            if len(v):
                print(f"  {p:22s} median {np.median(v):7.3g}  range {v.min():7.3g} to {v.max():7.3g}")


if __name__ == "__main__":
    main()
