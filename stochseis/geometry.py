"""Fault geometry as a patch graph.

The simulator only talks to the :class:`Fault` interface: per-patch arrays
(position, area, slip rate, rigidity, section id, local along-strike and
down-dip coordinates) plus a weighted adjacency graph used for rupture
growth.  A loader for UCERF3-style fault sections or a stochastic fault
network generator only has to produce the same object.

Stage 0 provides :func:`planar_fault`, a single rectangular fault that can
be cut into along-strike segments separated by step-overs.  Along-strike
distance ``s`` is a global coordinate that includes the gaps, so a rupture
spanning several segments is naturally "unrolled" into one (s, z) plane
for the slip-field construction.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse

from .params import Params

SQRT2 = np.sqrt(2.0)


class Fault:
    """Patch graph.

    Attributes (all arrays have length ``n``):

    s, zd       along-strike and down-dip coordinates of patch centres [m]
    x, y, z     Cartesian coordinates [m], z positive down
    depth       depth of patch centre [m]
    area        patch area [m^2]
    slip_rate   long-term slip rate [m/yr]
    mu          rigidity [Pa]
    section     integer section id
    col, row    grid column (global, counting gap columns) and row of each patch
    """

    def __init__(self, s, zd, x, y, z, area, slip_rate, mu, section, col, row,
                 dx, ns, nz, top_depth, dip_deg, rake_deg, intra_graph,
                 connections, col_s, jump_max):
        self.s = np.asarray(s, float)
        self.zd = np.asarray(zd, float)
        self.x = np.asarray(x, float)
        self.y = np.asarray(y, float)
        self.z = np.asarray(z, float)
        self.depth = self.z
        self.area = np.asarray(area, float)
        self.slip_rate = np.asarray(slip_rate, float)
        self.mu = np.asarray(mu, float)
        self.section = np.asarray(section, int)
        self.col = np.asarray(col, int)
        self.row = np.asarray(row, int)
        self.dx = float(dx)
        self.ns = int(ns)          # number of global columns (incl. gaps)
        self.nz = int(nz)
        self.top_depth = float(top_depth)
        self.dip_deg = float(dip_deg)
        self.rake_deg = float(rake_deg)
        self.intra_graph = intra_graph.tocsr()
        # connections: structured array (a, b, w, gap, conn_id)
        self.connections = connections
        self.col_s = np.asarray(col_s, float)   # s of each global column centre
        self.jump_max = float(jump_max)
        self.n = len(self.s)
        self.reaches_surface = self.top_depth <= 1e-6
        self.width = float(self.nz * self.dx)
        # lookup from (col, row) to patch index; -1 where no patch (gap columns)
        self.grid_index = -np.ones((self.ns, self.nz), dtype=int)
        self.grid_index[self.col, self.row] = np.arange(self.n)
        self._full_graph = None
        self.n_conn = len(np.unique(connections["conn_id"])) if len(connections) else 0
        self.conn_gap = np.zeros(self.n_conn)
        for k in range(self.n_conn):
            self.conn_gap[k] = connections["gap"][connections["conn_id"] == k][0]

    # ------------------------------------------------------------------
    @property
    def slip_dir(self) -> str:
        r = abs(((self.rake_deg + 180.0) % 360.0) - 180.0)
        return "strike" if (r < 45.0 or r > 135.0) else "dip"

    @property
    def moment_rate(self) -> float:
        """Tectonic moment rate sum(mu A v) [N m / yr]."""
        return float(np.sum(self.mu * self.area * self.slip_rate))

    @property
    def total_area(self) -> float:
        return float(np.sum(self.area))

    # ------------------------------------------------------------------
    def graph(self, open_conn=None) -> sparse.csr_matrix:
        """Adjacency graph with the given connections open.

        ``open_conn`` is a boolean array over connection ids; ``None`` opens
        all connections.  Edge weights are centre-to-centre distances [m].
        """
        if len(self.connections) == 0:
            return self.intra_graph
        if open_conn is None or np.all(open_conn):
            if self._full_graph is None:
                self._full_graph = self._with_connections(np.ones(self.n_conn, bool))
            return self._full_graph
        return self._with_connections(open_conn)

    def _with_connections(self, open_conn):
        c = self.connections
        sel = open_conn[c["conn_id"]]
        a, b, w = c["a"][sel], c["b"][sel], c["w"][sel]
        extra = sparse.coo_matrix((np.r_[w, w], (np.r_[a, b], np.r_[b, a])),
                                  shape=(self.n, self.n))
        return (self.intra_graph + extra.tocsr()).tocsr()

    # ------------------------------------------------------------------
    def nearest_patch(self, s, zd):
        """Snap (s, zd) positions [m] to the nearest existing patch."""
        s = np.atleast_1d(np.asarray(s, float))
        zd = np.atleast_1d(np.asarray(zd, float))
        j = np.searchsorted(self.col_s, s)
        j = np.clip(j, 1, len(self.col_s) - 1)
        left = self.col_s[j - 1]
        right = self.col_s[j]
        col = np.where(np.abs(s - left) <= np.abs(right - s), j - 1, j)
        row = np.clip(np.round(zd / self.dx - 0.5).astype(int), 0, self.nz - 1)
        idx = self.grid_index[col, row]
        # gap columns have no patch: move to the nearest column with one
        bad = np.where(idx < 0)[0]
        if len(bad):
            valid_cols = np.unique(self.col)
            for k in bad:
                jj = np.argmin(np.abs(valid_cols - col[k]))
                idx[k] = self.grid_index[valid_cols[jj], row[k]]
        return idx

    # ------------------------------------------------------------------
    def to_dict(self) -> dict:
        d = {k: getattr(self, k) for k in
             ("s", "zd", "x", "y", "z", "area", "slip_rate", "mu", "section",
              "col", "row", "col_s")}
        d["connections"] = self.connections
        d["scalars"] = np.array([self.dx, self.ns, self.nz, self.top_depth,
                                 self.dip_deg, self.rake_deg, self.jump_max])
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "Fault":
        dx, ns, nz, top, dip, rake, jmax = d["scalars"]
        n = len(d["s"])
        intra = _intra_graph(d["col"], d["row"], int(ns), int(nz), dx, n)
        return cls(d["s"], d["zd"], d["x"], d["y"], d["z"], d["area"],
                   d["slip_rate"], d["mu"], d["section"], d["col"], d["row"],
                   dx, ns, nz, top, dip, rake, intra, d["connections"],
                   d["col_s"], jmax)


CONN_DTYPE = np.dtype([("a", int), ("b", int), ("w", float), ("gap", float),
                       ("conn_id", int)])


def _intra_graph(col, row, ns, nz, dx, n):
    """8-connected graph within each contiguous run of columns."""
    grid = -np.ones((ns, nz), dtype=int)
    grid[col, row] = np.arange(n)
    rows_a, rows_b, w = [], [], []
    for dc, dr, wt in ((1, 0, dx), (0, 1, dx), (1, 1, dx * SQRT2), (1, -1, dx * SQRT2)):
        c0 = np.arange(max(0, -dc), ns - max(0, dc))
        r0 = np.arange(max(0, -dr), nz - max(0, dr))
        C0, R0 = np.meshgrid(c0, r0, indexing="ij")
        a = grid[C0, R0]
        b = grid[C0 + dc, R0 + dr]
        ok = (a >= 0) & (b >= 0)
        rows_a.append(a[ok])
        rows_b.append(b[ok])
        w.append(np.full(ok.sum(), wt))
    a = np.concatenate(rows_a)
    b = np.concatenate(rows_b)
    w = np.concatenate(w)
    g = sparse.coo_matrix((np.r_[w, w], (np.r_[a, b], np.r_[b, a])), shape=(n, n))
    return g.tocsr()


def planar_fault(p: Params) -> Fault:
    """Build a rectangular (possibly segmented) planar fault from parameters."""
    dx = p.patch_km * 1e3
    L = p.fault_length_km * 1e3
    W = p.fault_width_km * 1e3
    nz = max(1, int(round(W / dx)))
    ncol_fault = max(1, int(round(L / dx)))
    gaps = sorted([(float(a) * 1e3, float(b) * 1e3) for a, b in p.segment_gaps_km])
    # cut positions in column units
    cuts = [int(round(pos / dx)) for pos, _ in gaps]
    cuts = [c for c in cuts if 0 < c < ncol_fault]
    gap_cols = [max(1, int(round(g / dx))) for _, g in gaps][: len(cuts)]
    # global column layout: fault columns interleaved with gap columns
    col_global = []
    section = []
    col_s = []
    gpos = 0
    sec = 0
    cut_set = dict(zip(cuts, gap_cols))
    for c in range(ncol_fault):
        if c in cut_set:
            gpos += cut_set[c]
            sec += 1
        col_global.append(c + gpos)
        section.append(sec)
    ns = ncol_fault + gpos
    col_s = (np.arange(ns) + 0.5) * dx
    col_global = np.array(col_global)
    section = np.array(section)

    cols = np.repeat(col_global, nz)
    rows = np.tile(np.arange(nz), ncol_fault)
    secs = np.repeat(section, nz)
    s = col_s[cols]
    zd = (rows + 0.5) * dx
    dip = np.deg2rad(p.dip_deg)
    x = s
    y = zd * np.cos(dip)
    z = p.top_depth_km * 1e3 + zd * np.sin(dip)
    n = len(s)
    area = np.full(n, dx * dx)
    v = np.full(n, p.slip_rate_mm_yr * 1e-3)
    mu = np.full(n, p.rigidity_pa)

    intra = _intra_graph(cols, rows, ns, nz, dx, n)

    # connections across gaps: same row, last column of segment k to first of k+1
    grid = -np.ones((ns, nz), dtype=int)
    grid[cols, rows] = np.arange(n)
    conn = []
    jmax = p.jump_max_km * 1e3
    for k, c in enumerate(cuts):
        gap_m = gap_cols[k] * dx
        if gap_m > jmax:
            continue
        ca = col_global[c - 1]
        cb = col_global[c]
        for r in range(nz):
            conn.append((grid[ca, r], grid[cb, r], gap_m + dx, gap_m, k))
    connections = np.array(conn, dtype=CONN_DTYPE) if conn else np.zeros(0, dtype=CONN_DTYPE)

    return Fault(s, zd, x, y, z, area, v, mu, secs, cols, rows, dx, ns, nz,
                 p.top_depth_km * 1e3, p.dip_deg, p.rake_deg, intra,
                 connections, col_s, jmax)
