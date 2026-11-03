"""Slip distribution for a given rupture outline.

Steps (see docs/design.md):

1. unroll the rupture patches into their (s, z) bounding box;
2. solve the static crack problem for uniform stress drop on the outline
   with the spectral stiffness (free-surface mirror if the rupture reaches
   the surface); this smooth solution carries the moment, i.e. the
   stress-drop parameter is the macroscopic one that moment-area scaling
   constrains, and width saturation, end effects and the buried versus
   surface-breaking difference all come out of the elasticity;
3. draw a unit-variance von Karman Gaussian field with Mai & Beroza
   correlation lengths for the box dimensions and exponentiate it to a
   lognormal pattern with prescribed coefficient of variation;
4. weight by the local positive state measure (loading stress or slip
   deficit, see state.py), which conditions where slip concentrates
   (overshoot is allowed);
5. multiply crack solution, roughness and deficit weight, and rescale the
   product so that its moment equals the moment of the crack solution.

The roughness therefore redistributes slip within the outline without
changing the moment.  The realised moment-weighted stress drop of the
rough field (Noda et al., 2013) is returned for diagnostics.

Crack solutions depend only on the outline shape and are cached, which
matters because shells grown from a nucleation point produce many
identical small shapes.
"""

from __future__ import annotations

from collections import OrderedDict

import numpy as np

from .randomfield import gaussian_field, mai_beroza_correlation_lengths
from .stiffness import SpectralStiffness, crack_solution, moment_weighted_stress_drop

_CRACK_CACHE: "OrderedDict[tuple, np.ndarray]" = OrderedDict()
_CRACK_CACHE_MAX = 2000
CRACK_SOLVES = [0]   # counter for diagnostics


def unit_crack(mask, dx, slip_dir, surface):
    """Crack solution for unit dsig/mu on ``mask`` (cached)."""
    key = (mask.shape, bool(surface), slip_dir, dx, mask.tobytes())
    sol = _CRACK_CACHE.get(key)
    if sol is None:
        sol, _ = crack_solution(mask, dx, 1.0, slip_dir=slip_dir, mirror_top=surface, dsig=1.0)
        CRACK_SOLVES[0] += 1
        _CRACK_CACHE[key] = sol
        if len(_CRACK_CACHE) > _CRACK_CACHE_MAX:
            _CRACK_CACHE.popitem(last=False)
    else:
        _CRACK_CACHE.move_to_end(key)
    return sol


def build_slip(fault, idx, condition, dsig, rng, p, diagnostics=True):
    """Return (slip [m] per patch in ``idx``, info dict).

    condition: state measure of the patches in ``idx`` (stress [Pa] or
               deficit [m]; may be negative); only its positive part,
               normalised to mean one, is used;
    dsig:      target macroscopic stress drop [Pa].
    """
    dx = fault.dx
    mu = float(np.mean(fault.mu[idx]))
    col = fault.col[idx]
    row = fault.row[idx]
    c0, c1 = col.min(), col.max()
    r0, r1 = row.min(), row.max()
    nx = c1 - c0 + 1
    nz = r1 - r0 + 1
    ic = col - c0
    ir = row - r0
    mask = np.zeros((nx, nz), bool)
    mask[ic, ir] = True
    L = nx * dx
    W = nz * dx
    surface = bool(fault.reaches_surface and r0 == 0)

    # --- smooth crack solution ----------------------------------------
    # optional deep extension: columns that reach the base of the locked
    # width continue to slip over deep_slip_km below it (Shaw-type
    # argument for L^2 moment scaling); only the locked part is kept, so
    # the extension raises slip and moment in the seismogenic zone
    n_ext = int(round(p.deep_slip_km * 1e3 / dx)) if p.deep_slip_km > 0 else 0
    if n_ext > 0 and r1 == fault.nz - 1:
        mask_x = np.zeros((nx, nz + n_ext), bool)
        mask_x[:, :nz] = mask
        mask_x[:, nz:] = mask[:, -1:]
        Dc = unit_crack(mask_x, dx, fault.slip_dir, surface)[:, :nz] * (dsig / mu)
    else:
        Dc = unit_crack(mask, dx, fault.slip_dir, surface) * (dsig / mu)
    if getattr(p, "dip_correction", False):
        from .dipcorr import correction_factor

        dips = getattr(fault, "dip", None)
        dip = float(np.mean(dips[idx])) if isinstance(dips, np.ndarray) else float(fault.dip_deg)
        rakes = getattr(fault, "rake", None)
        rake = float(np.mean(rakes[idx])) if isinstance(rakes, np.ndarray) else float(fault.rake_deg)
        Dc = Dc * correction_factor(dip, rake, surface)
    dc = Dc[ic, ir]

    # --- roughness ----------------------------------------------------
    if nx * nz > 1 and p.slip_cv > 0:
        ax, az = mai_beroza_correlation_lengths(L, W)
        if p.corr_length_scatter > 0:
            ax *= np.exp(p.corr_length_scatter * rng.standard_normal())
            az *= np.exp(p.corr_length_scatter * rng.standard_normal())
        g = gaussian_field(nx, nz, dx, ax, az, p.hurst, rng)
        sg = np.sqrt(np.log1p(p.slip_cv ** 2))
        pattern = np.exp(sg * g[ic, ir] - 0.5 * sg ** 2)
    else:
        ax, az = mai_beroza_correlation_lengths(L, W)
        pattern = np.ones(len(idx))

    # --- state conditioning -------------------------------------------
    dpos = np.maximum(np.asarray(condition, float), 0.0)
    if p.condition_weight > 0 and np.mean(dpos) > 0:
        w = np.maximum((dpos / np.mean(dpos)) ** p.condition_weight, 1e-3)
    else:
        w = np.ones(len(idx))

    slip = dc * pattern * w
    m_crack = np.sum(dc * fault.area[idx])
    slip *= m_crack / np.sum(slip * fault.area[idx])

    info = dict(L=L, W=W, surface=surface, ax=ax, az=az,
                dmean=float(slip.mean()), dmax=float(slip.max()),
                dmean_crack=float(dc.mean()), dsig_m=np.nan)
    if diagnostics:
        D = np.zeros((nx, nz))
        D[ic, ir] = slip
        op = SpectralStiffness(nx, nz, dx, mu, fault.slip_dir, surface)
        info["dsig_m"] = moment_weighted_stress_drop(D[mask], op.apply(D)[mask])
    return slip, info
