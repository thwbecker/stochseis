"""Fault state: slip deficit, stress and maturity on every patch.

Deficit
    delta_i(t) = v_i t + delta_i(0) - R_i(t), with R_i the slip released
    by past events.  Overshoot (negative deficit) is allowed.

Stress
    The fault is loaded by backslip: the stress on the fault is the static
    stress of a dislocation equal to the deficit (zero outside the modelled
    fault), tau = K delta, computed with the same spectral stiffness that
    gives the crack solutions (free-surface mirror included).  Because K
    is linear, tau(t) = tau_static + tau_rate * t, where tau_rate = K v is
    computed once and tau_static is decremented by K slip_e after each
    event.  Small events are accumulated in a pending-slip buffer whose
    effect on the stress is approximated locally (self-stiffness of a
    patch) until the buffer is flushed with one full FFT, which happens
    whenever a rupture larger than ``flush_patches`` occurs or the buffer
    holds more than ``flush_events`` events.  Positive tau is loading; a past
    rupture leaves negative tau inside its outline and positive stress
    concentrations around it, so arrest and slip conditioning driven by
    tau see both shadows and concentrations.  Units: Pa.

Maturity
    A slowly evolving per-patch variable, here simply the cumulative slip
    released, kept as a hook for state-dependent rules (e.g. smoother
    slip or longer ruptures on mature faults).  It is recorded but not
    used by the Stage 0 rules.

Only the planar single-grid geometry is supported for the stress field at
this stage (the whole fault lives on one (ns, nz) grid, gap columns are
zero).  Other geometries fall back to deficit-based rules.
"""

from __future__ import annotations

import numpy as np

from .stiffness import SpectralStiffness


class FaultState:
    def __init__(self, fault, deficit0, with_stress=True, flush_patches=8, flush_events=200):
        self.fault = fault
        f = fault
        self.deficit0 = np.asarray(deficit0, float)
        self.released = np.zeros(f.n)
        self.released_sum = 0.0
        self.with_stress = bool(with_stress)
        self.flush_patches = int(flush_patches)
        self.flush_events = int(flush_events)
        self.n_flush = 0
        if self.with_stress:
            mu = float(np.mean(f.mu))
            self.op = SpectralStiffness(f.ns, f.nz, f.dx, mu, f.slip_dir,
                                        mirror_top=f.reaches_surface)
            self.tau_rate = self._stress_of(f.slip_rate)          # Pa / yr
            self.tau_static = self._stress_of(self.deficit0)      # Pa at t = 0
            self.pending = np.zeros(f.n)                          # buffered slip [m]
            self.n_pending = 0
            # self-stiffness of one patch (stress on a patch per unit slip on it)
            unit = np.zeros(f.n)
            i0 = int(np.argmin((f.col - f.ns // 2) ** 2 + (f.row - f.nz // 2) ** 2))
            unit[i0] = 1.0
            self.k_self = float(self._stress_of(unit)[i0])

    # ------------------------------------------------------------------
    def _grid(self, values, idx=None):
        f = self.fault
        g = np.zeros((f.ns, f.nz))
        if idx is None:
            g[f.col, f.row] = values
        else:
            g[f.col[idx], f.row[idx]] = values
        return g

    def _stress_of(self, slip_values, idx=None):
        """Stress on all patches [Pa] for a slip (or deficit) distribution."""
        f = self.fault
        S = self.op.apply(self._grid(slip_values, idx))
        return S[f.col, f.row]

    # ------------------------------------------------------------------
    def deficit(self, t, idx=None):
        f = self.fault
        if idx is None:
            return f.slip_rate * t + self.deficit0 - self.released
        return f.slip_rate[idx] * t + self.deficit0[idx] - self.released[idx]

    def mean_deficit(self, t):
        f = self.fault
        return float(np.mean(f.slip_rate) * t + np.mean(self.deficit0) - self.released_sum / f.n)

    def stress(self, t, idx=None):
        """Loading stress from the deficit [Pa]; requires with_stress.

        Pending (not yet flushed) slip is accounted for by the local
        self-stiffness only.
        """
        if idx is None:
            return self.tau_static + self.tau_rate * t - self.k_self * self.pending
        return (self.tau_static[idx] + self.tau_rate[idx] * t
                - self.k_self * self.pending[idx])

    def flush(self):
        """Apply all pending slip to the static stress with one FFT."""
        if self.with_stress and self.n_pending > 0:
            self.tau_static -= self._stress_of(self.pending)
            self.pending[:] = 0.0
            self.n_pending = 0
            self.n_flush += 1

    def maturity(self, idx=None):
        return self.released if idx is None else self.released[idx]

    # ------------------------------------------------------------------
    def apply_event(self, idx, slip):
        """Record released slip and update the static stress."""
        self.released[idx] += slip
        self.released_sum += float(np.sum(slip))
        if self.with_stress:
            self.pending[idx] += slip
            self.n_pending += 1
            if len(idx) > self.flush_patches or self.n_pending >= self.flush_events:
                self.flush()

    def stress_heterogeneity(self, t):
        """(mean, std, fraction negative) of the loading stress at time t."""
        self.flush()
        tau = self.stress(t)
        return float(tau.mean()), float(tau.std()), float(np.mean(tau < 0))
