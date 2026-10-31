"""ETAS building blocks.

Each event with magnitude M spawns a Poisson number of direct offspring
with mean

    A(M) = K0 * 10^(alpha (M - M_min)),

delayed by Omori-law times with density ~ (t + c)^-p and displaced from a
slip-weighted random point of the parent rupture by a distance r drawn
from f(r) ~ r (r^2 + d^2)^-q (2-D isotropic power-law kernel) in the
unrolled (s, z) plane.  K0 is set from a target branching ratio using the
expected productivity over the magnitude distribution.

Below the single-patch magnitude the simulator assigns Gutenberg-Richter
magnitudes directly; above it magnitudes emerge from rupture growth, so
the expectation of 10^(alpha (M - M_min)) is evaluated from a pilot sample
of grown ruptures rather than from a GR formula.
"""

from __future__ import annotations

import numpy as np

DAYS_PER_YEAR = 365.25


def omori_times(n, c_yr, p_exp, rng):
    """n random delays [yr] with density proportional to (t + c)^-p, p > 1."""
    u = rng.random(n)
    return c_yr * ((1.0 - u) ** (-1.0 / (p_exp - 1.0)) - 1.0)


def omori_cdf(t, c, p_exp):
    return 1.0 - (1.0 + t / c) ** (1.0 - p_exp)


def spatial_offsets(n, d, q, rng):
    """n (ds, dz) offsets with |r| ~ r (r^2 + d^2)^-q, q > 1, isotropic."""
    u = rng.random(n)
    r = d * np.sqrt((1.0 - u) ** (-1.0 / (q - 1.0)) - 1.0)
    th = 2.0 * np.pi * rng.random(n)
    return r * np.cos(th), r * np.sin(th)


def gr_magnitudes(n, b, mlo, mhi, rng):
    """n magnitudes from GR with slope b truncated to [mlo, mhi)."""
    u = rng.random(n)
    frac = 1.0 - 10.0 ** (-b * (mhi - mlo))
    return mlo - np.log10(1.0 - u * frac) / b


def gr_expected_productivity(b, alpha, mlo, mhi):
    """E[10^(alpha (M - mlo))] for GR(b) truncated to [mlo, mhi)."""
    dm = mhi - mlo
    if dm <= 0:
        return 1.0
    norm = 1.0 - 10.0 ** (-b * dm)
    if abs(b - alpha) < 1e-9:
        return b * np.log(10.0) * dm / norm
    return b / (b - alpha) * (1.0 - 10.0 ** (-(b - alpha) * dm)) / norm


def productivity_constant(n_target, alpha, mlo, mags_fault, p_fault, b, m_patch):
    """K0 such that the mean number of offspring per event is n_target.

    mags_fault: sample of fault-event magnitudes (from a pilot run),
    p_fault:    probability that a nucleation becomes a fault event,
    the remainder being sub-patch GR events in [mlo, m_patch).
    """
    e_fault = np.mean(10.0 ** (alpha * (np.asarray(mags_fault) - mlo))) if len(mags_fault) else 1.0
    e_sub = gr_expected_productivity(b, alpha, mlo, m_patch) if m_patch > mlo else 0.0
    expectation = p_fault * e_fault + (1.0 - p_fault) * e_sub
    return n_target / expectation
