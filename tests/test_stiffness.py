import numpy as np

from stochseis.stiffness import (antiplane_crack_slip, eshelby_circular_slip,
                                 eshelby_mean_slip, implied_mean_slip,
                                 moment_weighted_stress_drop, stress_drop_field)

MU = 3e10
DSIG = 3e6


def test_eshelby_circular_crack_uniform_stress_drop():
    """Full-space circular crack with Eshelby slip must give a uniform stress drop."""
    dx = 250.0
    a = 10e3
    n = 120
    x = (np.arange(n) - n / 2 + 0.5) * dx
    X, Z = np.meshgrid(x, x, indexing="ij")
    r = np.sqrt(X ** 2 + Z ** 2)
    D = eshelby_circular_slip(r, a, DSIG, MU)
    S = stress_drop_field(D, dx, MU, slip_dir="strike", mirror_top=False, pad_factor=6)
    inner = r < 0.8 * a
    rel = np.abs(S[inner] - DSIG) / DSIG
    assert np.median(rel) < 0.03
    # moment-weighted average over the whole crack close to the target as well
    dsm = moment_weighted_stress_drop(D, S)
    assert abs(dsm - DSIG) / DSIG < 0.05
    # mean slip matches analytic
    assert abs(D[r < a].mean() - eshelby_mean_slip(a, DSIG, MU)) / eshelby_mean_slip(a, DSIG, MU) < 0.03


def test_eshelby_dip_slip_direction():
    """Same test with slip along the dip direction (in-plane / anti-plane swap)."""
    dx = 250.0
    a = 10e3
    n = 120
    x = (np.arange(n) - n / 2 + 0.5) * dx
    X, Z = np.meshgrid(x, x, indexing="ij")
    r = np.sqrt(X ** 2 + Z ** 2)
    D = eshelby_circular_slip(r, a, DSIG, MU)
    S = stress_drop_field(D, dx, MU, slip_dir="dip", mirror_top=False, pad_factor=6)
    inner = r < 0.8 * a
    assert np.median(np.abs(S[inner] - DSIG) / DSIG) < 0.03


def test_antiplane_buried_crack():
    """Long buried strike-slip crack: antiplane solution, uniform stress drop."""
    dx = 250.0
    a = 7.5e3                     # half width, crack width 15 km
    nz = 120
    nx = 4000                     # long in strike: 1000 km
    zc = (np.arange(nz) - nz / 2 + 0.5) * dx
    prof = antiplane_crack_slip(zc, a, DSIG, MU)
    D = np.tile(prof[None, :], (nx, 1))
    S = stress_drop_field(D, dx, MU, slip_dir="strike", mirror_top=False, pad_factor=6)
    mid = nx // 2
    inner = np.abs(zc) < 0.8 * a
    rel = np.abs(S[mid, inner] - DSIG) / DSIG
    assert np.median(rel) < 0.03


def test_antiplane_surface_breaking_with_mirror():
    """Surface-breaking strike-slip fault: mirroring makes the half-width the full width."""
    dx = 250.0
    W = 15e3
    nz = int(W / dx)
    nx = 4000
    zd = (np.arange(nz) + 0.5) * dx
    prof = antiplane_crack_slip(zd, W, DSIG, MU)       # crack centred at the surface
    D = np.tile(prof[None, :], (nx, 1))
    S = stress_drop_field(D, dx, MU, slip_dir="strike", mirror_top=True, pad_factor=6)
    mid = nx // 2
    inner = zd < 0.8 * W
    rel = np.abs(S[mid, inner] - DSIG) / DSIG
    assert np.median(rel) < 0.03
    # the mean slip is then (pi/2) dsig W / mu, twice the buried value
    assert abs(prof.mean() - (np.pi / 2) * DSIG * W / MU) / ((np.pi / 2) * DSIG * W / MU) < 0.02


def test_implied_mean_slip_regimes():
    # small: circular; large: antiplane with effective width
    d_small = implied_mean_slip(np.pi * (1e3) ** 2, 30e3, DSIG, MU)
    assert abs(d_small - eshelby_mean_slip(1e3, DSIG, MU)) < 1e-9
    d_large = implied_mean_slip(1e3 * 1e9, 30e3, DSIG, MU)
    assert abs(d_large - (np.pi / 4) * DSIG / MU * 30e3) < 1e-9
