"""Accept parameter vectors whose synthetic site records match the data.

    python3 examples/alpine/abc_accept.py summary.(json|csv) [threshold_m] [tolerance] [--seismicity]

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
from abc_scoring import accepted, load_summary, score  # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)


def main():
    seis = "--seismicity" in sys.argv
    argv = [a for a in sys.argv if a != "--seismicity"]
    rows = load_summary(argv[1])
    th = float(argv[2]) if len(argv) > 2 else 2.0
    tol = float(argv[3]) if len(argv) > 3 else 2.0
    scored = score(rows, th, seis)
    print(f"{len(scored)} design vectors, threshold {th} m{', with seismicity' if seis else ''}; "
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
