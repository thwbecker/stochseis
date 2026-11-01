"""Static stress change from slip on a planar fault, in the wavenumber domain.

For slip D(x, z) on an infinite plane in a homogeneous full space the
static stress change is tau(k) = -K(k) D(k) with (Andrews, 1980)

    K(k) = mu / (2 k) * [ 2 (lam + mu) / (lam + 2 mu) k_par^2 + k_perp^2 ],

where k_par is the wavenumber component along the slip direction.  For a
Poisson solid (lam = mu) the bracket is 4/3 k_par^2 + k_perp^2.  This is
exact for the planar full-space problem and is evaluated with FFTs on a
zero-padded grid.

Free surface: for a vertical fault that breaks the surface the slip field
is mirrored about z = 0 before applying K.  For the antiplane (strike-slip)
components this image makes the surface traction free exactly; for the
in-plane (dip-slip) components it is an approximation.  See docs/design.md.

Sign convention: stress *drop* is positive where slip releases stress.

:func:`crack_solution` solves the inverse problem K D = dsig on a given
outline (zero slip outside), i.e. the static crack with uniform stress
drop, by conjugate gradients with the spectral operator.  This gives the
smooth, crack-like long-wavelength slip for an arbitrary rupture shape,
including the free-surface image and width saturation, without any
scaling relation.
"""

from __future__ import annotations

import numpy as np
from scipy.fft import irfft2, next_fast_len, rfft2
from scipy.sparse.linalg import LinearOperator, cg


def stiffness_kernel(nx, nz, dx, mu, lam=None, slip_dir="strike"):
    """K(kx, kz) on an (nx, nz) grid with spacing dx [m]; returns Pa/m."""
    if lam is None:
        lam = mu
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, dx)
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    k = np.sqrt(KX ** 2 + KZ ** 2)
    if slip_dir == "strike":
        kpar, kperp = KX, KZ
    else:
        kpar, kperp = KZ, KX
    fac = 2.0 * (lam + mu) / (lam + 2.0 * mu)
    with np.errstate(divide="ignore", invalid="ignore"):
        K = mu / (2.0 * k) * (fac * kpar ** 2 + kperp ** 2)
    K[0, 0] = 0.0
    return K


class SpectralStiffness:
    """Spectral stiffness operator for an (nx, nz) slip grid.

    The grid is zero-padded by ``pad_factor`` in each direction and to at
    least ``min_pad`` points so that very small ruptures (one or a few
    patches) see a well resolved kernel.  Periodic images reduce the
    stress drop by roughly 1 percent at pad_factor 4 when the rupture
    fills the grid.
    """

    def __init__(self, nx, nz, dx, mu, slip_dir="strike", mirror_top=False,
                 pad_factor=4.0, min_pad=64, lam=None):
        self.nx, self.nz = int(nx), int(nz)
        self.mirror = bool(mirror_top)
        self.nzm = 2 * self.nz if self.mirror else self.nz
        self.Nx = next_fast_len(max(int(np.ceil(pad_factor * self.nx)) + 1, min_pad), True)
        self.Nz = next_fast_len(max(int(np.ceil(pad_factor * self.nzm)) + 1, min_pad), True)
        K = stiffness_kernel(self.Nx, self.Nz, dx, mu, lam=lam, slip_dir=slip_dir)
        self.K = K[:, : self.Nz // 2 + 1]           # real-FFT half spectrum
        kmin = np.min(self.K[self.K > 0])
        self.Kinv = np.where(self.K > 0, 1.0 / np.maximum(self.K, kmin), 1.0 / kmin)
        self.shape = (self.Nx, self.Nz)

    def _embed(self, D):
        Dm = np.concatenate([D[:, ::-1], D], axis=1) if self.mirror else D
        Dp = np.zeros((self.Nx, self.Nz))
        Dp[: self.nx, : self.nzm] = Dm
        return Dp

    def _extract(self, S):
        S = S[: self.nx, : self.nzm]
        return S[:, self.nz:] if self.mirror else S

    def apply(self, D):
        """Stress drop field [Pa] for slip D [m]."""
        return self._extract(irfft2(self.K * rfft2(self._embed(D)), s=self.shape))

    def apply_inverse(self, S):
        """Approximate inverse (used as a preconditioner)."""
        return self._extract(irfft2(self.Kinv * rfft2(self._embed(S)), s=self.shape))


def stress_drop_field(D, dx, mu, slip_dir="strike", mirror_top=False,
                      pad_factor=4.0, lam=None, min_pad=64):
    """Stress drop [Pa] on the same grid as the slip D [m] (axis 1 down-dip, row 0 top)."""
    D = np.asarray(D, float)
    op = SpectralStiffness(D.shape[0], D.shape[1], dx, mu, slip_dir, mirror_top,
                           pad_factor, min_pad, lam)
    return op.apply(D)


def moment_weighted_stress_drop(D, S):
    """Noda, Lapusta & Kanamori (2013) moment-weighted average stress drop."""
    w = np.sum(D)
    if w <= 0:
        return np.nan
    return float(np.sum(D * S) / w)


def crack_solution(mask, dx, mu, slip_dir="strike", mirror_top=False,
                   dsig=1.0, rtol=1e-4, maxiter=400):
    """Slip [m] for uniform stress drop ``dsig`` on the cells of ``mask``.

    Solves K D = dsig on the mask (D = 0 elsewhere) with preconditioned
    conjugate gradients.  Returns (D, n_iterations).  The solution is
    linear in dsig / mu, so callers can cache the unit solution.
    """
    mask = np.asarray(mask, bool)
    nx, nz = mask.shape
    op = SpectralStiffness(nx, nz, dx, mu, slip_dir, mirror_top)
    cells = np.where(mask)
    n = len(cells[0])

    def matvec(v):
        D = np.zeros((nx, nz))
        D[cells] = v
        return op.apply(D)[cells]

    def precond(v):
        S = np.zeros((nx, nz))
        S[cells] = v
        return op.apply_inverse(S)[cells]

    A = LinearOperator((n, n), matvec=matvec, dtype=float)
    M = LinearOperator((n, n), matvec=precond, dtype=float)
    rhs = np.full(n, float(dsig))
    x0 = precond(rhs)
    it = [0]

    def cb(xk):
        it[0] += 1

    try:
        x, info = cg(A, rhs, x0=x0, M=M, rtol=rtol, maxiter=maxiter, callback=cb)
    except TypeError:                     # scipy < 1.12 calls the tolerance "tol"
        x, info = cg(A, rhs, x0=x0, M=M, tol=rtol, maxiter=maxiter, callback=cb)
    D = np.zeros((nx, nz))
    D[cells] = x
    return D, it[0]


# ----------------------------------------------------------------------
# analytic reference solutions (used in tests and for the arrest rule)
# ----------------------------------------------------------------------
def eshelby_circular_slip(r, a, dsig, mu, nu=0.25):
    """Slip on a circular crack of radius a with uniform stress drop dsig.

    D(r) = 24 / (7 pi) * dsig / mu * sqrt(a^2 - r^2) for nu = 0.25.
    General nu: D = 8 (1 - nu) / (pi (2 - nu)) * dsig / mu * sqrt(a^2 - r^2),
    equivalent to Eshelby's dsig = 7 pi / 16 * mu * mean(D) / a at nu = 0.25.
    """
    c = 8.0 * (1.0 - nu) / (np.pi * (2.0 - nu))
    r = np.asarray(r, float)
    out = np.zeros_like(r)
    inside = r < a
    out[inside] = c * dsig / mu * np.sqrt(a ** 2 - r[inside] ** 2)
    return out


def eshelby_mean_slip(a, dsig, mu, nu=0.25):
    """Mean slip over a circular crack: (2/3) * D(0) = 16/(7 pi) dsig a / mu at nu = 0.25."""
    c = 8.0 * (1.0 - nu) / (np.pi * (2.0 - nu))
    return 2.0 / 3.0 * c * dsig / mu * a


def antiplane_crack_slip(zc, a, dsig, mu):
    """2-D antiplane (screw) crack of half-width a centred at 0: D = 2 dsig/mu sqrt(a^2 - z^2)."""
    zc = np.asarray(zc, float)
    out = np.zeros_like(zc)
    inside = np.abs(zc) < a
    out[inside] = 2.0 * dsig / mu * np.sqrt(a ** 2 - zc[inside] ** 2)
    return out


def implied_mean_slip(area, w_eff, dsig, mu):
    """Rough mean slip for a rupture of given area under constant stress drop.

    Circular crack scaling below width saturation, antiplane crack scaling
    of the effective width above it.  Used only by the arrest rule to
    judge whether the local slip deficit can support further growth; the
    realised slip comes from the crack solution.
    """
    a_circ = np.sqrt(area / np.pi)
    d_circ = eshelby_mean_slip(a_circ, dsig, mu)
    d_long = (np.pi / 4.0) * dsig / mu * w_eff   # mean of antiplane crack, half-width w_eff/2
    return min(d_circ, d_long)
