import numpy as np

from stochseis.mesh import load_mesh_section
from stochseis.sites import record_stats, site_patch


def _write_stl(path, tris):
    with open(path, "w") as f:
        f.write("solid test\n")
        for t in tris:
            f.write("  facet normal 0 0 0\n    outer loop\n")
            for v in t:
                f.write(f"      vertex {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
            f.write("    endloop\n  endfacet\n")
        f.write("endsolid test\n")


def _surface(L, depth, dip_fn, curvature, h=0.5):
    """Triangulated surface in km, depth negative: trace along an arc, dip varying along strike."""
    ns = int(L / h)
    nz = int(depth / h)
    pts = np.zeros((ns + 1, nz + 1, 3))
    for i in range(ns + 1):
        s = i * h
        ang = curvature * s
        x0, y0 = np.sin(ang) / curvature if curvature else s, (1 - np.cos(ang)) / curvature if curvature else 0.0
        tx, ty = np.cos(ang), np.sin(ang)
        nxh, nyh = -ty, tx                         # horizontal normal to the trace
        dip = np.deg2rad(dip_fn(s))
        for j in range(nz + 1):
            z = j * h
            off = z / np.tan(dip)
            pts[i, j] = (x0 + off * nxh + 1000.0, y0 + off * nyh + 5000.0, -z)
    tris = []
    for i in range(ns):
        for j in range(nz):
            a, b, c, d = pts[i, j], pts[i + 1, j], pts[i + 1, j + 1], pts[i, j + 1]
            tris.append([a, b, c])
            tris.append([a, c, d])
    return np.array(tris)


def test_planar_dipping_mesh(tmp_path):
    tris = _surface(60.0, 15.0, lambda s: 50.0, 0.0)
    path = tmp_path / "plane.stl"
    _write_stl(path, tris)
    f = load_mesh_section(path, 1000.0, units="km")
    W = 15.0 / np.sin(np.deg2rad(50.0))
    assert f.ns == 60
    assert f.nz == int(W)
    assert abs(np.median(f.dip) - 50.0) < 0.5
    assert abs(f.area.sum() / 1e6 - f.n) < 1e-6 * f.n + 1.0      # unit Jacobian
    assert f.reaches_surface
    assert abs(f.mesh_info["trace_length_m"] - 60e3) < 10.0


def test_curved_variable_dip_mesh(tmp_path):
    tris = _surface(80.0, 15.0, lambda s: 80.0 if s < 40 else 50.0, 1.0 / 300.0)
    path = tmp_path / "curved.stl"
    _write_stl(path, tris)
    f = load_mesh_section(path, 1000.0, units="km")
    mesh_area = f.mesh_info["mesh_area_m2"]
    # partial cells at the base and the far end are dropped (up to about 5 percent here)
    assert 0.93 < f.area.sum() / mesh_area <= 1.0
    assert abs(np.median(f.area) / 1e6 - 1.0) < 0.01
    south = f.s < 35e3
    north = f.s > 45e3
    assert abs(np.median(f.dip[south]) - 80.0) < 2.0
    assert abs(np.median(f.dip[north]) - 50.0) < 2.0
    # ragged bottom: steep part narrower down dip than the shallow-dipping part
    assert f.row[south].max() < f.row[north].max()
    # site mapping picks a surface patch near the requested point
    k, d = site_patch(f, f.x[f.grid_index[20, 0]], f.y[f.grid_index[20, 0]])
    assert f.row[k] == 0 and d < 1.0


def test_record_stats():
    st = record_stats([0, 100, 200, 300, 400], [5, 5, 5, 5, 5])
    assert st["n"] == 5 and abs(st["ri"] - 100) < 1e-9 and st["cv"] == 0.0


def test_dip_correction():
    from stochseis.dipcorr import correction_factor
    from stochseis.geometry import planar_fault
    from stochseis.params import Params
    from stochseis.slip import build_slip

    assert abs(correction_factor(90.0, 180.0, True) - 0.989) < 1e-9
    assert abs(correction_factor(15.0, 90.0, True) - 1.487) < 1e-9
    f50 = correction_factor(50.0, 160.0, True)
    assert 1.05 < f50 < 1.15
    p0 = Params(fault_length_km=60, fault_width_km=15, dip_deg=50.0, rake_deg=160.0, slip_cv=0.0)
    p1 = p0.replace(dip_correction=True)
    fault = planar_fault(p0)
    idx = np.where(fault.col < 40)[0]
    s0, _ = build_slip(fault, idx, np.ones(len(idx)), 3e6, np.random.default_rng(0), p0)
    s1, _ = build_slip(fault, idx, np.ones(len(idx)), 3e6, np.random.default_rng(0), p1)
    assert abs(s1.mean() / s0.mean() - f50) < 1e-6
