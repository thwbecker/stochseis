"""Tests for the persistent hierarchical patch structure (stochseis/hierarchy.py)."""
import numpy as np

from stochseis import Params, Simulation
from stochseis.geometry import planar_fault
from stochseis.hierarchy import (barrier, build_hierarchy, hierarchical_patches, hierarchy_radii,
                                 nucleation_factor, scale_map)


def test_radii():
    r = hierarchy_radii(1000.0, 40000.0, 2.0)
    assert np.allclose(r, [40000, 20000, 10000, 5000, 2500, 1250])


def test_patch_counts_and_nesting():
    rng = np.random.default_rng(1)
    L, W = 200e3, 15e3
    radii = hierarchy_radii(1000.0, 8000.0, 2.0)
    lv = hierarchical_patches(L, W, radii, 0.3, 1.0, rng)
    n = np.array([len(x) for x in lv])
    expected = np.maximum(1, np.round(0.3 * L * W / (np.pi * radii ** 2)))
    assert np.array_equal(n, expected)                       # number grows as r^-2
    for k in range(1, len(radii)):                           # nest = 1: every centre inside a parent disc
        d = np.min(np.hypot(lv[k][:, None, 0] - lv[k - 1][None, :, 0], lv[k][:, None, 1] - lv[k - 1][None, :, 1]),
                   axis=1)
        assert np.all(d <= radii[k - 1] - radii[k] + 1e-6)
    lv0 = hierarchical_patches(L, W, radii, 0.3, 0.0, np.random.default_rng(1))
    d0 = np.min(np.hypot(lv0[-1][:, None, 0] - lv0[-2][None, :, 0], lv0[-1][:, None, 1] - lv0[-2][None, :, 1]),
                axis=1)
    assert np.mean(d0 <= radii[-2]) < 0.8                    # nest = 0: placed independently


def test_scale_map_and_reproducibility():
    p = Params(fault_length_km=100, fault_width_km=15, seed=3, hier_rmin_km=1.0, hier_rmax_km=8.0)
    f = planar_fault(p)
    r1, radii, levels = build_hierarchy(f, p)
    r2, _, _ = build_hierarchy(f, p)
    assert np.array_equal(r1, r2)
    allowed = set(np.round(radii).tolist()) | {round(2.0 * radii[0])}
    assert set(np.round(np.unique(r1)).tolist()) <= allowed
    r3, _, _ = build_hierarchy(f, p.replace(seed=4))
    assert not np.array_equal(r1, r3)
    # every patch's scale is the smallest disc containing it
    s, z = f.s[::37], f.zd[::37]
    direct = scale_map(s, z, radii, levels, 2.0)
    assert np.array_equal(direct, r1[::37])


def test_barrier_and_nucleation_factor():
    assert barrier(1000.0, np.pi * 2000.0 ** 2) == 1.0       # rupture larger than the patch: no barrier
    assert np.isclose(barrier(8000.0, np.pi * 1000.0 ** 2), 8.0)
    r = np.array([1000.0, 1000.0, 4000.0, 16000.0])
    w = nucleation_factor(r, 1000.0, 2.0)
    assert np.isclose(w.mean(), 1.0) and w[0] > w[2] > w[3]
    assert np.all(nucleation_factor(r, 1000.0, 0.0) == 1.0)


def test_simulation_with_hierarchy():
    base = dict(duration_yr=300, mag_min=3.5, seed=5, fault_length_km=60, fault_width_km=15,
                hier_rmin_km=1.0, hier_rmax_km=8.0, hier_barrier_exponent=1.0)
    sim = Simulation(Params(**base, hier_nucleation_exponent=3.0))
    cat = sim.run()
    assert abs(cat.meta["moment_closure"] - 1.0) < 0.5
    rmin = sim.hier_radii[-1]
    small = np.isclose(sim.hier_scale, rmin)
    nuc = cat.events["nuc"][cat.events["kind"] == 1]
    # nucleation concentrates in the smallest patches
    assert np.mean(small[nuc]) > 1.5 * np.mean(small)
    # option off: no hierarchy objects, no change of interface
    sim0 = Simulation(Params(**{**base, "hier_rmin_km": 0.0}))
    assert sim0.hier_scale is None and sim0.nuc_factor is None
