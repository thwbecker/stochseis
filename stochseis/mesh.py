# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Triangulated fault surfaces to structured section grids.

A fault surface given as a triangle mesh (for example the STL reference
geometries of the New Zealand Alpine Fault intercomparison) is
parametrised by an along-strike coordinate s and a down-dip coordinate d
and resampled onto a regular (s, d) grid of square patches.  The result
is a :class:`~stochseis.geometry.Fault` whose patches carry their real 3-D
positions, normals, dips, strikes and areas, while rupture growth, random
fields and the spectral crack solutions keep working on the unrolled
(s, d) grid.  This is the same representation that fault-section
databases use (a trace plus a down-dip profile; UCERF3 subsections, the
NZ Community Fault Model), so section-based loaders only have to produce
the (s, d) -> (x, y, z) map that is built here from a mesh.

Assumptions of the parametrisation (checked and reported by
:func:`load_mesh_section`):

* the surface is a single sheet whose top edge is a polyline (the trace,
  at the shallowest depth of the mesh) and whose down-dip lines are close
  to perpendicular to the trace in map view;
* s is the arc length of the map-view projection of a point onto the
  trace; d is the arc length along the surface, measured in the vertical
  plane perpendicular to the trace, from the top edge;
* the base of the surface may lie at different down-dip distances along
  strike (variable dip at fixed depth); grid cells below the base are
  omitted, so columns are ragged at the bottom.

Curvature of the trace and dip changes are represented exactly in the
3-D patch positions but only approximately in the static stiffness,
which is evaluated on the unrolled grid (see docs/networks.md).

Coordinates: x east, y north in the units of the input (converted to m),
z positive down (depth).  STL files with negative z for depth are
converted.
"""

from __future__ import annotations

import struct

import numpy as np
from scipy import sparse
from scipy.interpolate import LinearNDInterpolator
from scipy.sparse.csgraph import connected_components

from .geometry import CONN_DTYPE, Fault, _intra_graph


# ----------------------------------------------------------------------
# reading
# ----------------------------------------------------------------------
def read_stl(path):
    """Triangle vertices (n, 3, 3) from an ASCII or binary STL file."""
    with open(path, "rb") as f:
        head = f.read(80)
        rest = f.read()
    is_ascii = head.lstrip().lower().startswith(b"solid") and b"facet" in rest[:1000]
    if is_ascii:
        v = []
        for line in (head + rest).decode("ascii", errors="ignore").splitlines():
            s = line.split()
            if s and s[0] == "vertex":
                v.append((float(s[1]), float(s[2]), float(s[3])))
        return np.array(v, float).reshape(-1, 3, 3)
    n = struct.unpack("<I", rest[:4])[0]
    rec = np.frombuffer(rest[4:4 + 50 * n], dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)),
                                                            ("a", "<u2")]))
    return rec["v"].astype(float)


def weld(tri, tol):
    """Merge vertices closer than tol; returns (vertices, faces)."""
    pts = tri.reshape(-1, 3)
    key = np.round(pts / tol).astype(np.int64)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    verts = np.zeros((len(uniq), 3))
    np.add.at(verts, inv, pts)
    verts /= np.bincount(inv)[:, None]
    return verts, inv.reshape(-1, 3)


# ----------------------------------------------------------------------
# trace and parametrisation
# ----------------------------------------------------------------------
def _boundary_edges(faces):
    e = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    ue, cnt = np.unique(e, axis=0, return_counts=True)
    return ue[cnt == 1], int(np.sum(cnt > 2))


def _top_trace(verts, faces, ztol):
    """Ordered polyline of the top boundary (vertices at the minimum depth)."""
    bnd, nonmanifold = _boundary_edges(faces)
    ztop = verts[:, 2].min()
    on_top = np.abs(verts[:, 2] - ztop) <= ztol
    e = bnd[on_top[bnd[:, 0]] & on_top[bnd[:, 1]]]
    nodes = np.unique(e)
    idx = {v: k for k, v in enumerate(nodes)}
    g = sparse.coo_matrix((np.ones(len(e)), ([idx[a] for a in e[:, 0]], [idx[b] for b in e[:, 1]])),
                          shape=(len(nodes), len(nodes)))
    g = (g + g.T).tocsr()
    ncomp, lab = connected_components(g, directed=False)
    if ncomp > 1:
        keep = np.argmax(np.bincount(lab))
        sel = lab == keep
    else:
        sel = np.ones(len(nodes), bool)
    deg = np.asarray(g.sum(axis=1)).ravel()
    ends = np.where((deg == 1) & sel)[0]
    start = ends[0] if len(ends) else np.where(sel)[0][0]
    # walk the chain
    order = [start]
    prev = -1
    cur = start
    while True:
        nb = g.indices[g.indptr[cur]:g.indptr[cur + 1]]
        nb = [k for k in nb if k != prev and sel[k]]
        if not nb or nb[0] in order:
            break
        prev, cur = cur, nb[0]
        order.append(cur)
    trace = verts[nodes[np.array(order)]]
    # orient so that x increases from start to end overall (convention only)
    if trace[-1, 0] < trace[0, 0]:
        trace = trace[::-1]
    return trace, ztop, nonmanifold, ncomp


def _project_on_polyline(p2, poly2):
    """Arc length s, signed perpendicular distance h and segment index of map points on a polyline."""
    a = poly2[:-1]
    b = poly2[1:]
    ab = b - a
    seglen = np.linalg.norm(ab, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seglen)])
    best_d = np.full(len(p2), np.inf)
    best_s = np.zeros(len(p2))
    best_h = np.zeros(len(p2))
    for k in range(len(a)):
        t = np.clip(((p2 - a[k]) @ ab[k]) / seglen[k] ** 2, 0.0, 1.0)
        q = a[k] + t[:, None] * ab[k]
        dv = p2 - q
        dist = np.hypot(dv[:, 0], dv[:, 1])
        better = dist < best_d
        if np.any(better):
            best_d[better] = dist[better]
            best_s[better] = cum[k] + t[better] * seglen[k]
            cross = ab[k, 0] * dv[better, 1] - ab[k, 1] * dv[better, 0]
            best_h[better] = np.sign(cross) * dist[better]
    return best_s, best_h, cum[-1]


def _simplify_trace(trace, step):
    """Resample the trace polyline at roughly ``step`` spacing (speeds up projection)."""
    seg = np.linalg.norm(np.diff(trace[:, :2], axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    n = max(2, int(np.ceil(cum[-1] / step)) + 1)
    s = np.linspace(0.0, cum[-1], n)
    return np.column_stack([np.interp(s, cum, trace[:, k]) for k in range(3)])


def parametrise(verts, trace, dx):
    """Along-strike s and down-dip d for every vertex."""
    tr = _simplify_trace(trace, max(dx, 250.0))
    s, h, L = _project_on_polyline(verts[:, :2], tr[:, :2])
    depth = verts[:, 2] - verts[:, 2].min()
    # down-dip arc length: in columns of width dx, integrate sqrt(dh^2 + dz^2) along depth
    d = np.zeros(len(verts))
    col = np.floor(s / dx).astype(int)
    for c in np.unique(col):
        m = np.where(col == c)[0]
        o = m[np.argsort(depth[m])]
        z = depth[o]
        hh = h[o]
        # robust down-dip profile: median offset from the trace in depth bins of half a
        # patch (unstructured meshes have vertices at many depths whose offsets scatter
        # where the dip changes along strike; a per-vertex profile would zig-zag)
        zbin = max(0.5 * dx, 250.0)
        kb = np.floor(z / zbin).astype(int)
        lev = np.unique(kb)
        zg = np.array([np.median(z[kb == k]) for k in lev])
        hg = np.array([np.median(hh[kb == k]) for k in lev])
        # anchor the profile at the top of the surface
        if zg[0] > 0.0:
            top = z <= z.min() + 1e-6
            zg = np.concatenate([[z.min()], zg])
            hg = np.concatenate([[np.median(hh[top])], hg])
        seg = np.hypot(np.diff(hg), np.diff(zg))
        dg = np.concatenate([[0.0], np.cumsum(seg)])
        d[o] = np.interp(z, zg, dg)
    return s, d, h, L


# ----------------------------------------------------------------------
# resampling to a section grid
# ----------------------------------------------------------------------
def load_mesh_section(path, dx, units="auto", weld_tol=1.0, slip_rate=None, rake_deg=180.0,
                      rigidity=3.0e10, section_breaks_km=(), verbose=False):
    """Read a triangulated fault surface and return a :class:`Fault`.

    dx                patch size [m]
    units             "m", "km" or "auto" (km if the horizontal extent is below 1e4)
    slip_rate         callable s[m] -> slip rate [m/yr], a constant [m/yr], or None (1e-3 m/yr placeholder)
    rake_deg          rake used for the stiffness (scalar for the section)
    section_breaks_km along-strike positions that split the fault into sections (segment labels only;
                      the grid stays continuous and ruptures cross freely)
    """
    tri = read_stl(path)
    ext = np.ptp(tri.reshape(-1, 3)[:, :2], axis=0).max()
    if units == "km" or (units == "auto" and ext < 1e4):
        tri = tri * 1e3
    if np.median(tri[..., 2]) < 0:          # depth negative in the file
        tri = tri.copy()
        tri[..., 2] = -tri[..., 2]
    verts, faces = weld(tri, weld_tol)
    trace, ztop, nonmanifold, ncomp = _top_trace(verts, faces, ztol=max(1.0, 1e-3 * dx))
    s, d, h, L = parametrise(verts, trace, dx)

    # base of the surface along strike: maximum d per column of width dx
    ncol = int(np.floor(L / dx))
    col_v = np.clip(np.floor(s / dx).astype(int), 0, ncol - 1)
    dmax = np.zeros(ncol)
    np.maximum.at(dmax, col_v, d)
    # fill empty columns by interpolation
    good = dmax > 0
    if not np.all(good):
        dmax[~good] = np.interp(np.where(~good)[0], np.where(good)[0], dmax[good])
    nrow = int(np.floor(dmax.max() / dx))

    # interpolate positions on the (s, d) grid
    interp = LinearNDInterpolator(np.column_stack([s, d]), verts)
    cols, rows = np.meshgrid(np.arange(ncol), np.arange(nrow), indexing="ij")
    sc = (cols + 0.5) * dx
    dc = (rows + 0.5) * dx
    valid = dc + 0.5 * dx <= dmax[cols] + 0.25 * dx
    pts = interp(sc[valid], dc[valid])
    ok = np.all(np.isfinite(pts), axis=1)
    cv, rv = cols[valid][ok], rows[valid][ok]
    xyz = pts[ok]

    # tangent vectors and areas by finite differences of the map
    def at(sv, dv):
        return interp(sv, dv)
    eps = 0.25 * dx
    s0, d0 = (cv + 0.5) * dx, (rv + 0.5) * dx
    sp, sm = np.minimum(s0 + eps, L - 1e-3), np.maximum(s0 - eps, 1e-3)
    dp, dm = np.minimum(d0 + eps, dmax[cv] - 1e-3), np.maximum(d0 - eps, 1e-3)
    ts = (at(sp, d0) - at(sm, d0)) / (sp - sm)[:, None]       # d(xyz)/ds
    td = (at(s0, dp) - at(s0, dm)) / (dp - dm)[:, None]       # d(xyz)/dd
    bad = ~np.all(np.isfinite(ts), axis=1) | ~np.all(np.isfinite(td), axis=1)
    ts[bad] = np.nan
    td[bad] = np.nan
    nrm = np.cross(ts, td)
    nn = np.linalg.norm(nrm, axis=1)
    # Jacobian of the map: area of the patch relative to dx^2
    jac = nn.copy()
    jac[~np.isfinite(jac)] = np.nanmedian(jac)
    area = dx * dx * jac
    nrm = nrm / nn[:, None]
    nrm[~np.all(np.isfinite(nrm), axis=1)] = np.nanmedian(nrm, axis=0)
    nrm[nrm[:, 2] < 0] *= -1.0            # normals pointing down (z is depth)
    dip = np.degrees(np.arccos(np.clip(nrm[:, 2], -1, 1)))
    dip = 90.0 - np.abs(90.0 - dip)
    strike = np.degrees(np.arctan2(ts[:, 0], ts[:, 1])) % 360.0

    n = len(cv)
    s_patch = (cv + 0.5) * dx
    if slip_rate is None:
        v = np.full(n, 1e-3)
    elif callable(slip_rate):
        v = np.asarray(slip_rate(s_patch), float)
    else:
        v = np.full(n, float(slip_rate))
    breaks = np.sort(np.asarray(section_breaks_km, float) * 1e3)
    section = np.searchsorted(breaks, s_patch)
    intra = _intra_graph(cv, rv, ncol, nrow, dx, n)
    fault = Fault(s_patch, (rv + 0.5) * dx, xyz[:, 0], xyz[:, 1], xyz[:, 2], area, v,
                  np.full(n, rigidity), section, cv, rv, dx, ncol, nrow, ztop,
                  float(np.median(dip)), rake_deg, intra, np.zeros(0, dtype=CONN_DTYPE),
                  (np.arange(ncol) + 0.5) * dx, 5e3)
    fault.dip = dip
    fault.strike = strike
    fault.normal = nrm
    fault.trace = trace
    fault.mesh_info = dict(path=str(path), n_triangles=int(len(tri)), n_vertices=int(len(verts)),
                           nonmanifold_edges=nonmanifold, top_components=int(ncomp),
                           trace_length_m=float(L), mesh_area_m2=float(
                               0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0],
                                                             tri[:, 2] - tri[:, 0]), axis=1).sum()),
                           grid_area_m2=float(area.sum()), n_patches=int(n),
                           width_min_m=float(dmax.min()), width_max_m=float(dmax.max()))
    if verbose:
        print({k: (round(v, 1) if isinstance(v, float) else v) for k, v in fault.mesh_info.items()})
    return fault
