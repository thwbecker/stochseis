"""L1 correction of crack-solution amplitudes for dip, rake and burial.

The spectral operator (unrolled grid, free surface by an image) is exact
for vertical strike-slip faults but not for dipping or dip-slip ruptures
or for ruptures that do not reach the surface (E20 in docs/experiments.md).
The factors below are the ratio of mean slip from half-space triangular
dislocation kernels (cutde; Nikkhoo and Walter, 2015) to mean slip from
the spectral operator, for 60 x 20 km ruptures with uniform stress drop
at 2 km resolution (examples/checks/l1_table.py). Buried means a top
edge at 3 km depth. For oblique slip the strike-slip and dip-slip factors
are weighted by cos^2 and sin^2 of the rake. The correction scales the
whole crack solution; it does not change the shape of the slip
distribution near the surface trace.
"""

from __future__ import annotations

import numpy as np

DIPS = np.array([15.0, 30.0, 45.0, 60.0, 75.0, 90.0])
FACTORS = {
    ("strike", True): np.array([1.907, 1.320, 1.129, 1.042, 1.001, 0.989]),
    ("strike", False): np.array([1.219, 1.134, 1.091, 1.067, 1.055, 1.051]),
    ("dip", True): np.array([1.487, 1.181, 1.140, 1.182, 1.242, 1.269]),
    ("dip", False): np.array([1.182, 1.163, 1.164, 1.171, 1.177, 1.179]),
}


def correction_factor(dip_deg, rake_deg, surface):
    """Half-space / spectral mean-slip ratio for a rupture of the given mean dip and rake."""
    d = float(np.clip(abs(dip_deg), DIPS[0], DIPS[-1]))
    fs = np.interp(d, DIPS, FACTORS[("strike", bool(surface))])
    fd = np.interp(d, DIPS, FACTORS[("dip", bool(surface))])
    r = np.deg2rad(rake_deg)
    return float(np.cos(r) ** 2 * fs + np.sin(r) ** 2 * fd)
