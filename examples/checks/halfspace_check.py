"""Compare the FFT/image crack solutions with full elastic half-space kernels.

For a rectangular rupture with uniform stress drop, solve K D = dsig with
(a) the spectral stiffness used by stochseis (full space, free surface by
mirror image) and (b) half-space triangular dislocation kernels from cutde
(Nikkhoo and Walter, 2015), and compare mean and maximum slip.  Requires
``pip install cutde``.

    python examples/checks/halfspace_check.py
"""
import numpy as np
import cutde.halfspace as hs

from stochseis.stiffness import crack_solution

MU = 3e10
NU = 0.25
DSIG = 3e6


def plane(L, W, dx, dip_deg, top):
    """Patch centres and triangles for a plane; x along strike, z up (negative depth)."""
    dip = np.deg2rad(dip_deg)
    ns, nz = int(round(L / dx)), int(round(W / dx))
    def pt(s, d):
        return np.array([s, d * np.cos(dip), -(top + d * np.sin(dip))])
    tris, cen = [], []
    for i in range(ns):
        for j in range(nz):
            s0, s1, d0, d1 = i * dx, (i + 1) * dx, j * dx, (j + 1) * dx
            a, b, c, e = pt(s0, d0), pt(s1, d0), pt(s1, d1), pt(s0, d1)
            tris.append([a, b, c])
            tris.append([a, c, e])
            cen.append(pt(0.5 * (s0 + s1), 0.5 * (d0 + d1)))
    tris = np.array(tris)
    # observe at triangle centroids (the patch centre lies on the shared
    # diagonal, where the triangular-dislocation solution is singular)
    cen = tris.mean(axis=1)
    return tris, cen, ns, nz, dip


def halfspace_matrix(tris, cen, comp, dip):
    """Shear traction in the slip direction at patch centres per unit slip of each patch."""
    strain = np.asarray(hs.strain_matrix(cen, tris, NU))
    if strain.shape[1] != 6:                            # (n_obs, n_tri, 3, 6) layout
        strain = strain.transpose(0, 3, 1, 2)
    # strain: (n_obs, 6, n_tri, 3)
    st = strain.transpose(0, 2, 3, 1).reshape(-1, 6)
    stress = np.asarray(hs.strain_to_stress(st, MU, NU))
    stress = stress.reshape(cen.shape[0], tris.shape[0], 3, 6)[:, :, comp, :]   # (obs, tri, 6)
    sxx, syy, szz, sxy, sxz, syz = [stress[..., k] for k in range(6)]
    n = np.array([0.0, -np.sin(dip), -np.cos(dip)])     # normal to the plane
    if comp == 0:
        v = np.array([1.0, 0.0, 0.0])                   # strike direction
    else:
        v = np.array([0.0, np.cos(dip), -np.sin(dip)])  # down-dip direction
    tx = sxx * n[0] + sxy * n[1] + sxz * n[2]
    ty = sxy * n[0] + syy * n[1] + syz * n[2]
    tz = sxz * n[0] + syz * n[1] + szz * n[2]
    tau = tx * v[0] + ty * v[1] + tz * v[2]             # (obs, tri)
    K = tau[:, 0::2] + tau[:, 1::2]                     # two source triangles per patch
    K = 0.5 * (K[0::2, :] + K[1::2, :])                 # average over the two observation triangles
    if np.mean(np.diag(K)) > 0:                         # make stress drop positive for positive slip
        K = -K
    return -K                                           # stress drop per unit slip


def run(name, L, W, dx, dip_deg, top, slip_dir):
    tris, cen, ns, nz, dip = plane(L, W, dx, dip_deg, top)
    comp = 0 if slip_dir == "strike" else 1
    K = halfspace_matrix(tris, cen, comp, dip)
    D_hs = np.linalg.solve(K, np.full(K.shape[0], DSIG))
    surface = top <= 0
    D_fft, _ = crack_solution(np.ones((ns, nz), bool), dx, MU, slip_dir=slip_dir,
                              mirror_top=surface, dsig=DSIG)
    D_fft_full, _ = crack_solution(np.ones((ns, nz), bool), dx, MU, slip_dir=slip_dir,
                                   mirror_top=False, dsig=DSIG)
    print(f"{name:42s} mean slip half-space {D_hs.mean():5.2f} m | FFT {'image' if surface else 'full '} "
          f"{D_fft.mean():5.2f} m (ratio {D_fft.mean() / D_hs.mean():4.2f}) | FFT no image "
          f"{D_fft_full.mean():5.2f} m | max {D_hs.max():5.2f} / {D_fft.max():5.2f} m")


if __name__ == "__main__":
    run("vertical strike-slip, surface, 60x15 km", 60e3, 15e3, 1e3, 90.0, 0.0, "strike")
    run("vertical strike-slip, buried 5 km, 60x15 km", 60e3, 15e3, 1e3, 90.0, 5e3, "strike")
    run("thrust 15 deg, trench-breaking, 60x40 km", 60e3, 40e3, 2e3, 15.0, 0.0, "dip")
    run("thrust 15 deg, buried 5 km, 60x40 km", 60e3, 40e3, 2e3, 15.0, 5e3, "dip")
    run("thrust 45 deg, surface, 40x20 km", 40e3, 20e3, 1e3, 45.0, 0.0, "dip")
