# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Persistent hierarchical (multi-scale) patch structure of a fault.

A quenched structure in the spirit of hierarchical-patch models of
earthquake growth (Ide and Aochi, 2005; Aochi and Ide, 2009): circular
patches at a sequence of scales, radii r_max, r_max / ratio, ... down to
about r_min, each level covering a fraction ``fill`` of the fault (so the
number of patches grows as r^-2), with smaller patches placed inside a
patch of the next larger level with probability ``nest``. The structure
is drawn once per simulation and does not change.

Each patch of the fault grid is assigned the radius of the smallest
hierarchical patch that contains it (the background, outside all patches,
gets ratio x r_max). The simulator uses this scale map in two ways
(stochseis.simulator, stochseis.rupture):

* a size-dependent barrier in rupture admission: a rupture of equivalent
  radius R entering a patch of scale r has its admission ratio divided by
  max(1, r / R)^eta, so ruptures smaller than the patch they try to enter
  need a correspondingly larger deficit (cascade growth from small into
  progressively larger patches);
* optionally, a nucleation preference for small scales: nucleation weights
  are multiplied by (r_min / r)^kappa, normalized to a fault mean of one.

Coordinates are the unrolled along-strike and down-dip coordinates
(fault.s, fault.zd).
"""

from __future__ import annotations

import numpy as np


def hierarchy_radii(rmin, rmax, ratio):
    """Patch radii from largest to smallest: rmax, rmax / ratio, ... while >= rmin (within 1 percent)."""
    if rmin <= 0 or rmax < rmin or ratio <= 1:
        raise ValueError("need 0 < rmin <= rmax and ratio > 1")
    radii = []
    r = float(rmax)
    while r >= 0.99 * rmin:
        radii.append(r)
        r /= ratio
    return np.array(radii)


def hierarchical_patches(L, W, radii, fill, nest, rng):
    """Patch centres per level (largest first): list of arrays (n_k, 2) of (s, z) in metres.

    Level k has max(1, round(fill * L * W / (pi r_k^2))) patches. A patch of
    level k > 0 is placed, with probability ``nest``, uniformly inside a
    randomly chosen patch of level k - 1 (its whole disc inside the parent
    where possible); otherwise uniformly on the fault.
    """
    levels = []
    for k, r in enumerate(radii):
        n = max(1, int(round(fill * L * W / (np.pi * r * r))))
        cen = np.empty((n, 2))
        for j in range(n):
            if k > 0 and rng.random() < nest:
                parent = levels[k - 1][rng.integers(len(levels[k - 1]))]
                rr = max(radii[k - 1] - r, 0.0) * np.sqrt(rng.random())
                th = 2 * np.pi * rng.random()
                cen[j] = (parent[0] + rr * np.cos(th), parent[1] + rr * np.sin(th))
            else:
                cen[j] = (rng.uniform(0, L), rng.uniform(0, W))
        levels.append(cen)
    return levels


def scale_map(s, z, radii, levels, ratio):
    """Radius of the smallest patch containing each point (background: ratio x largest radius)."""
    r_map = np.full(len(s), ratio * radii[0])
    for r, cen in zip(radii, levels):                     # largest first, smaller overwrite
        for c in cen:
            inside = (s - c[0]) ** 2 + (z - c[1]) ** 2 <= r * r
            r_map[inside] = np.minimum(r_map[inside], r)
    return r_map


def build_hierarchy(fault, p):
    """Scale map [m] per fault patch and the patch levels, from the hier_* parameters.

    Uses its own random stream (seed + 104729), so the structure is the same for a
    given seed regardless of other settings, and can be rebuilt for analysis.
    """
    rng = np.random.default_rng(p.seed + 104729)
    radii = hierarchy_radii(p.hier_rmin_km * 1e3, p.hier_rmax_km * 1e3, p.hier_ratio)
    s0, z0 = fault.s.min() - 0.5 * fault.dx, fault.zd.min() - 0.5 * fault.dx
    L = fault.s.max() - s0 + 0.5 * fault.dx
    W = fault.zd.max() - z0 + 0.5 * fault.dx
    levels = hierarchical_patches(L, W, radii, p.hier_fill, p.hier_nest, rng)
    levels = [lv + np.array([s0, z0]) for lv in levels]
    r_map = scale_map(fault.s, fault.zd, radii, levels, p.hier_ratio)
    return r_map, radii, levels


def nucleation_factor(r_map, rmin, kappa):
    """Static nucleation multiplier (rmin / r)^kappa, normalized to a fault mean of one."""
    if kappa == 0:
        return np.ones(len(r_map))
    w = (rmin / r_map) ** kappa
    return w / w.mean()


def barrier(r_patch, cum_area):
    """Size-dependent barrier factor max(1, r_patch / R) for a rupture of area cum_area (R = sqrt(A / pi))."""
    R = np.sqrt(max(cum_area, 1.0) / np.pi)
    return max(1.0, r_patch / R)
