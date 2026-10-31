"""Rupture growth and arrest.

A rupture nucleates on one patch and floods outward through the patch
graph.  Candidate patches are the graph neighbours of already admitted
patches; they are visited in order of path distance from the nucleation
point (ties broken at random), so growth is radial, and only contiguous
sets of patches can rupture.  Two rules act on each candidate.

Admission (local, state dependent).  A candidate i is admitted with
probability

    q_i = 1 - exp(-(rho_i / rho_c)^2),

    state_measure = "deficit": rho_i = deficit_i^+ / implied mean slip(A)
    state_measure = "stress":  rho_i = tau_i^+ / dsig

and is otherwise blocked for this rupture (blocked patches that end up
completely enclosed by the rupture are added afterwards, since they are
forced to slip).  Because only neighbours of
admitted patches are candidates, a region in which q falls below the
site-percolation threshold of the lattice cannot be crossed: the rupture
grows around depleted patches and stops when its frontier is exhausted.
This is the state-dependent arrest; rho_c = 0 disables it.

Size hazard (global, state independent).  Before the admitted area grows
from A_k to A_{k+1} the rupture stops with probability

    p_k = 1 - (A_k / A_{k+1})^b                       (stop_measure = "area")
        = 1 - exp(-(2b/3) (ln M0_{k+1} - ln M0_k))     (stop_measure = "moment")

so that on an unbounded fault N(>A) ~ A^-b: GR with the given b below
width saturation by construction, independent of the patch size.  When
the hazard stops the rupture at k patches, a continuous area is drawn
within [A_k, A_{k+1}) from the same hazard (the rupture stopped somewhere
before reaching the next patch) and returned as ``area_scale`` = A / A_k;
the caller rescales slip by sqrt(area_scale) (crack scaling D ~ sqrt(A)),
which removes the discreteness of magnitudes at small patch counts.

Cost per event: one shortest-path computation for the ordering plus a
heap flood that is linear in the number of admitted patches.
"""

from __future__ import annotations

import heapq

import numpy as np
from scipy.sparse.csgraph import dijkstra

STOP_HAZARD = 0
STOP_STATE = 1
STOP_BOUNDARY = 2

_UNSEEN, _QUEUED, _ADMITTED, _BLOCKED = 0, 1, 2, 3


def _implied_slip(area, w_eff, dsig, mu):
    return min((16.0 / (7.0 * np.pi)) * dsig / mu * np.sqrt(area / np.pi),
               (np.pi / 4.0) * dsig / mu * w_eff)


def grow_rupture(fault, graph, i0, state_values_fn, dsig, rng, p, w_eff, m0_table=None):
    """Grow a rupture from patch ``i0``.

    state_values_fn(idx) returns the state measure (stress [Pa] or deficit
    [m]) of patches idx at the current time.  ``m0_table[k-1]`` is the
    realised moment of a k-patch rupture at unit stress-drop ratio (from
    the pilot ensemble); with the moment hazard it replaces the analytic
    crack law for patch counts it covers, so that the hazard is defined
    on the moments the discretised crack solutions actually produce.
    Returns (patch indices in growth order, distances [m] from the
    nucleation patch, stop code, rho of the last considered patch,
    scale), where scale multiplies the slip to represent the continuous
    stopping size within the last increment (moment scale for the moment
    hazard, sqrt of the area scale for the area hazard).
    """
    dist = dijkstra(graph, directed=True, indices=i0)   # graph is symmetric
    if np.sum(np.isfinite(dist)) == 1:
        return np.array([i0]), dist[[i0]], STOP_BOUNDARY, np.inf, 1.0
    key = dist + rng.random(len(dist)) * (1e-3 * fault.dx)
    indptr, indices = graph.indptr, graph.indices
    area = fault.area
    mu = float(fault.mu[i0])
    b = p.b_value
    rho_c = p.arrest_scale
    use_state = rho_c > 0
    use_moment = p.stop_measure == "moment"
    n_tab = 0 if m0_table is None else len(m0_table)

    def ln_m0(k, a):
        if k <= n_tab:
            return np.log(m0_table[k - 1])
        # analytic law, matched to the table at its last entry
        base = np.log(mu * a * _implied_slip(a, w_eff, dsig, mu))
        if n_tab:
            a_t = n_tab * area[i0]
            base += np.log(m0_table[-1]) - np.log(mu * a_t * _implied_slip(a_t, w_eff, dsig, mu))
        return base
    if use_state:
        reach = np.where(np.isfinite(dist))[0]
        vals = np.zeros(len(dist))
        vals[reach] = np.maximum(np.asarray(state_values_fn(reach), float), 0.0)

    status = np.zeros(len(dist), np.int8)
    status[i0] = _QUEUED
    heap = [(key[i0], i0)]
    admitted = []
    cumA = 0.0
    lnm_prev = None
    code = STOP_BOUNDARY
    rho_last = np.inf
    any_blocked = False
    scale = 1.0
    while heap:
        _, i = heapq.heappop(heap)
        a_new = cumA + area[i]
        k = len(admitted)
        if k > 0:
            # admission
            if use_state:
                if p.state_measure == "stress":
                    rho = vals[i] / dsig
                else:
                    rho = vals[i] / _implied_slip(a_new, w_eff, dsig, mu)
                rho_last = rho
                q = 1.0 - np.exp(-(rho / rho_c) ** 2)
                if rng.random() >= q:
                    status[i] = _BLOCKED
                    any_blocked = True
                    continue
            # size hazard
            if use_moment:
                lnm = ln_m0(k + 1, a_new)
                if lnm_prev is None:
                    lnm_prev = ln_m0(k, cumA)
                beta = 2.0 * b / 3.0
                p_stop = 1.0 - np.exp(-beta * (lnm - lnm_prev))
                if rng.random() < p_stop:
                    code = STOP_HAZARD
                    # continuous stopping moment within [M0_k, M0_k+1): survival (M0_k/M0)^beta
                    u = rng.random()
                    scale = (1.0 - u * p_stop) ** (-1.0 / beta)
                    break
                lnm_prev = lnm
            else:
                p_stop = 1.0 - (cumA / a_new) ** b
                if rng.random() < p_stop:
                    code = STOP_HAZARD
                    u = rng.random()
                    scale = np.sqrt((1.0 - u * p_stop) ** (-1.0 / b))
                    break
        status[i] = _ADMITTED
        admitted.append(i)
        cumA = a_new
        for j in indices[indptr[i]:indptr[i + 1]]:
            if status[j] == _UNSEEN:
                status[j] = _QUEUED
                heapq.heappush(heap, (key[j], j))
    else:
        code = STOP_STATE if any_blocked else STOP_BOUNDARY
    idx = np.array(admitted, int)
    idx, d = fill_holes(fault, idx, dist)
    return idx, d, code, float(rho_last), float(scale)


def fill_holes(fault, idx, dist):
    """Add patches that are completely enclosed by the rupture.

    A patch surrounded by slipping patches is forced to slip, so blocked
    patches that do not connect to the exterior of the rupture outline
    are added.  Returns the extended (idx, dist), with distances of the
    added patches taken from the shortest-path array ``dist``.
    """
    from scipy.ndimage import binary_fill_holes

    if len(idx) < 8:
        return idx, dist[idx]
    col = fault.col[idx]
    row = fault.row[idx]
    c0, c1 = col.min(), col.max()
    r0, r1 = row.min(), row.max()
    mask = np.zeros((c1 - c0 + 1, r1 - r0 + 1), bool)
    mask[col - c0, row - r0] = True
    filled = binary_fill_holes(mask)
    new = np.argwhere(filled & ~mask)
    if len(new) == 0:
        return idx, dist[idx]
    cand = fault.grid_index[new[:, 0] + c0, new[:, 1] + r0]
    cand = cand[(cand >= 0) & np.isfinite(dist[np.maximum(cand, 0)])]
    if len(cand) == 0:
        return idx, dist[idx]
    idx2 = np.concatenate([idx, cand])
    return idx2, dist[idx2]
