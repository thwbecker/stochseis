# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
import numpy as np

from stochseis.geometry import planar_fault
from stochseis.params import Params
from stochseis.randomfield import fit_von_karman, gaussian_field, radial_psd
from stochseis.slip import build_slip
from stochseis.stiffness import moment_weighted_stress_drop, stress_drop_field


def test_field_variance_and_spectrum():
    rng = np.random.default_rng(0)
    n = 512
    dx = 1.0
    a = 12.0
    H = 0.75
    vars_ = []
    Hs, As = [], []
    for k in range(6):
        f = gaussian_field(n, n, dx, a, a, H, rng)
        vars_.append(f.var())
        kk, pp = radial_psd(f, dx)
        sel = kk < 0.5 * np.pi / dx       # avoid the grid corner
        h, aa = fit_von_karman(kk[sel], pp[sel])
        Hs.append(h)
        As.append(aa)
    assert abs(np.mean(vars_) - 1.0) < 0.1
    assert abs(np.mean(Hs) - H) < 0.1
    assert abs(np.mean(As) / a - 1.0) < 0.2


def test_anisotropic_field_correlation():
    """Longer correlation along strike must give higher lag correlation along strike."""
    rng = np.random.default_rng(1)
    f = gaussian_field(1024, 256, 1.0, 40.0, 5.0, 0.75, rng)
    lag = 5
    cx = np.corrcoef(f[:-lag, :].ravel(), f[lag:, :].ravel())[0, 1]
    cz = np.corrcoef(f[:, :-lag].ravel(), f[:, lag:].ravel())[0, 1]
    assert cx > cz + 0.2


def test_slip_moment_equals_crack_moment_and_positivity():
    p = Params(fault_length_km=60, fault_width_km=15, patch_km=1.0)
    fault = planar_fault(p)
    rng = np.random.default_rng(2)
    # rectangular rupture 30 km long, full width, reaching the surface
    idx = np.where((fault.col >= 10) & (fault.col < 40))[0]
    dsig = 3e6
    slip, info = build_slip(fault, idx, np.ones(len(idx)), dsig, rng, p)
    assert np.all(slip > 0)
    assert info["surface"]
    # moment of the rough field equals the moment of the smooth crack solution
    assert abs(slip.mean() - info["dmean_crack"]) / info["dmean_crack"] < 1e-9
    # smooth crack: recompute its stress drop and check it is uniform ~ dsig
    p0 = p.replace(slip_cv=0.0)
    slip0, info0 = build_slip(fault, idx, np.ones(len(idx)), dsig, rng, p0)
    nx, nz = 30, fault.nz
    D = np.zeros((nx, nz))
    D[fault.col[idx] - 10, fault.row[idx]] = slip0
    S = stress_drop_field(D, fault.dx, fault.mu[0], "strike", mirror_top=True)
    assert np.median(np.abs(S[D > 0] - dsig) / dsig) < 0.02
    # roughness raises the moment-weighted stress drop above the macroscopic one
    assert info["dsig_m"] > dsig
    # amplitude in the right range for a surface-breaking 15 km wide fault
    # (pi/2 dsig W / mu = 2.4 m is the infinitely long antiplane limit)
    assert 0.8 < info["dmean"] < 2.4


def test_deficit_conditioning_moves_slip():
    p = Params(fault_length_km=60, fault_width_km=15, patch_km=1.0)
    fault = planar_fault(p)
    idx = np.where((fault.col >= 10) & (fault.col < 40))[0]
    deficit = np.where(fault.col[idx] < 25, 0.1, 3.0)   # left half depleted
    means = []
    for seed in range(5):
        slip, _ = build_slip(fault, idx, deficit, 3e6, np.random.default_rng(seed), p)
        left = slip[fault.col[idx] < 25].mean()
        right = slip[fault.col[idx] >= 25].mean()
        means.append(right / left)
    assert np.mean(means) > 5.0


def test_buried_vs_surface_amplitude():
    """Same outline, buried fault (top at 5 km) should carry about half the slip."""
    p_s = Params(fault_length_km=100, fault_width_km=15, patch_km=1.0, top_depth_km=0.0)
    p_b = p_s.replace(top_depth_km=5.0)
    out = {}
    for name, p in (("surf", p_s), ("bur", p_b)):
        fault = planar_fault(p)
        idx = np.where((fault.col >= 20) & (fault.col < 80))[0]
        d = []
        for seed in range(8):
            slip, info = build_slip(fault, idx, np.ones(len(idx)), 3e6,
                                    np.random.default_rng(seed), p)
            d.append(info["dmean"])
        out[name] = np.mean(d)
    ratio = out["surf"] / out["bur"]
    assert 1.4 < ratio < 2.2


def test_single_patch_slip():
    p = Params(fault_length_km=20, fault_width_km=10, patch_km=1.0, top_depth_km=5.0)
    fault = planar_fault(p)
    slip, info = build_slip(fault, np.array([55]), np.ones(1), 3e6, np.random.default_rng(0), p)
    # uniform stress drop on a 1 km square patch: slip of order dsig dx / mu
    assert 0.2 * 3e6 * 1e3 / 3e10 < slip[0] < 2.0 * 3e6 * 1e3 / 3e10
