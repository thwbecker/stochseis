# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Von Karman random fields by spectral synthesis.

Power spectral density (Mai & Beroza, 2002, eq. 1), with angular
wavenumbers and separate correlation lengths along strike (ax) and down
dip (az):

    P(kx, kz) = (1 + ax^2 kx^2 + az^2 kz^2)^(-(H + 1)).

The field is white Gaussian noise filtered by sqrt(P) and normalised to
unit variance using the discrete mean of P over the grid, so that small
grids do not inherit the bias of a sample-variance normalisation.  The
fractal dimension of the field is D = 3 - H.
"""

from __future__ import annotations

import numpy as np
from scipy.fft import fft2, ifft2


def von_karman_psd(kx, kz, ax, az, H):
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    return (1.0 + (ax * KX) ** 2 + (az * KZ) ** 2) ** (-(H + 1.0))


def gaussian_field(nx, nz, dx, ax, az, H, rng):
    """Zero-mean, unit-variance Gaussian von Karman field of shape (nx, nz).

    dx, ax, az in the same length unit.
    """
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, dx)
    P = von_karman_psd(kx, kz, ax, az, H)
    noise = rng.standard_normal((nx, nz))
    f = np.real(ifft2(np.sqrt(P) * fft2(noise)))
    return f / np.sqrt(np.mean(P))


def mai_beroza_correlation_lengths(L, W):
    """Mai & Beroza (2002) eq. 4: ax = 2 + L/3, az = 1 + W/3 (km in, m out)."""
    Lk = L / 1e3
    Wk = W / 1e3
    return (2.0 + Lk / 3.0) * 1e3, (1.0 + Wk / 3.0) * 1e3


def radial_psd(field, dx, nbins=40):
    """Circularly averaged power spectrum of a 2-D field (for validation)."""
    nx, nz = field.shape
    F = fft2(field - field.mean())
    P = np.abs(F) ** 2 / (nx * nz)
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, dx)
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    k = np.sqrt(KX ** 2 + KZ ** 2).ravel()
    P = P.ravel()
    kmax = np.pi / dx
    kmin = 2.0 * np.pi / (max(nx, nz) * dx)
    edges = np.logspace(np.log10(kmin), np.log10(kmax), nbins + 1)
    idx = np.digitize(k, edges) - 1
    kc, pc = [], []
    for i in range(nbins):
        m = idx == i
        if m.sum() > 3:
            kc.append(np.exp(np.mean(np.log(k[m]))))
            pc.append(np.mean(P[m]))
    return np.array(kc), np.array(pc)


def fit_von_karman(k, P, a0=None):
    """Least-squares fit of log P = log P0 - (H+1) log(1 + a^2 k^2); returns (H, a)."""
    from scipy.optimize import least_squares

    lk = np.log(P)

    def resid(x):
        lp0, H, la = x
        return lp0 - (H + 1.0) * np.log1p(np.exp(2 * la) * k ** 2) - lk

    if a0 is None:
        a0 = 1.0 / k[len(k) // 3]
    x0 = np.array([lk[0], 0.7, np.log(a0)])
    sol = least_squares(resid, x0)
    return sol.x[1], np.exp(sol.x[2])
