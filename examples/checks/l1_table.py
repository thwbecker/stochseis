"""Tabulate the L1 dip correction: half-space / spectral-image crack slip.

For rectangular ruptures 60 km long and 20 km wide down dip with uniform
stress drop, compute the mean slip with half-space triangular dislocation
kernels (cutde) and with the spectral operator used by stochseis (free
surface by an image), for a range of dips, pure strike and pure dip slip,
and ruptures that reach the surface or are buried (top at 3 km). Prints a
table that is pasted into stochseis/dipcorr.py. Requires cutde.

    python3 examples/checks/l1_table.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from halfspace_check import DSIG, MU, halfspace_matrix, plane  # noqa: E402

from stochseis.stiffness import crack_solution  # noqa: E402

dips = [15, 30, 45, 60, 75, 90]
L, W, dx = 60e3, 20e3, 2e3
print("dip  slip    surface  buried(3 km)")
table = {}
for slip in ("strike", "dip"):
    for surf in (True, False):
        row = []
        for dip in dips:
            tris, cen, ns, nz, dipr = plane(L, W, dx, float(dip), 0.0 if surf else 3e3)
            comp = 0 if slip == "strike" else 1
            K = halfspace_matrix(tris, cen, comp, dipr)
            D_hs = np.linalg.solve(K, np.full(K.shape[0], DSIG))
            D_fft, _ = crack_solution(np.ones((ns, nz), bool), dx, MU, slip_dir=slip, mirror_top=surf, dsig=DSIG)
            row.append(D_hs.mean() / D_fft.mean())
        table[(slip, surf)] = row
for slip in ("strike", "dip"):
    for k, dip in enumerate(dips):
        print(f"{dip:3d}  {slip:6s}  {table[(slip, True)][k]:.3f}    {table[(slip, False)][k]:.3f}")
print("DIPS =", dips)
for key, row in table.items():
    print(f"{key}: [" + ", ".join(f"{v:.3f}" for v in row) + "]")
