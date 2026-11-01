# Single-fault experiments

This log records the tests run on a single planar strike-slip fault
(200 x 15 km, surface breaking, 20 mm/yr, 1 km patches, unless stated),
what each was meant to establish, the settings, the results and the
conclusion. Spec files that reproduce the multi-seed experiments with the
current code are in `examples/experiments/`;
`examples/experiments/run_all.sh [nproc] [outdir]` runs all of them and
`examples/aggregate_sweep.py` averages over seeds. Summaries of the runs
quoted below are in `runs/regime_2026-11-01/` and `runs/tests_2026-11-02/`.

Two warnings apply to the whole log.

* All results obtained before 1 November 2026 used a nucleation rate
  that was several times too high (E2). They are listed only where they
  led to a change in the code, and their numbers should not be compared
  with later ones.
* Unless stated otherwise, numbers are means over three seeds of 3000 yr
  with a 10 percent burn-in. Seed-to-seed standard deviations are 0.01
  to 0.1 for the CVs and Bath's gap; single-seed values of the site
  recurrence CV have a standard error of about 0.2.

Metrics: b-values are maximum-likelihood estimates for continuous
magnitudes (E4) in the sub-patch range, between the patch magnitude
(M4.2) and width saturation (about M6.6), and above saturation. Bath's
gap and Omori p use window-selected mainshocks (M >= 6, no larger event
within 1 yr and max(10 km, L)). Site recurrence CV is the coefficient of
variation of the times at which a patch slips more than 0.3 m or 1 m,
averaged over nine sites along strike at mid depth. Whole-fault CV uses
declustered events of M >= 6.5 or M >= 7 anywhere on the fault.

## E1. Numerical building blocks (test suite)

Purpose: verify the pieces before interpreting any statistics.

* Spectral stiffness: uniform stress drop recovered to 3 percent for the
  Eshelby circular crack (strike and dip slip), the long buried
  antiplane crack, and the surface-breaking antiplane crack with the
  free-surface image; the crack solver reproduces Eshelby mean slip to
  1 percent.
* Von Karman fields: unit variance, Hurst exponent and correlation
  length recovered from the spectrum; anisotropy correct.
* Growth: without state arrest N(>=k patches) = 1/k for b = 1 (area
  hazard) to within sampling error; a depleted zone blocks ruptures.
* ETAS pieces: Omori delays pass a Kolmogorov-Smirnov test; GR sampling
  and productivity integrals checked analytically.
* Since 1 November: analytic calibration against a 20,000-event pilot;
  deficit-weighted nucleation, state-only mode and gamma > 1 run and
  conserve their invariants; a regression test checks that the defaults
  give GR across the patch magnitude.

Reproduce: `python -m pytest tests`.

## E2. Rate calibration (1 November)

Purpose: find why the fault ran at a low mean deficit and lacked large
events.

The nucleation rate is set so that the expected moment per nucleation
times the rate equals the moment rate. The expectation had been a Monte
Carlo mean over 300 to 1000 pilot ruptures. Moment survival goes as
M0^(-2/3), so the mean is dominated by rare large events: a 600-event
pilot gave 2.0e15 N m against an analytic value of 2.0e16 N m, and a
30,000-event pilot gave 2.6e16 +- 1.0e16 N m. The rate had therefore
been up to ten times too high. Replaced by the analytic expectation
(`docs/design.md` section 2.5). Conclusion: the earlier shortage of
M >= 6.5 events relative to GR and the low equilibrium deficit were
calibration artefacts.

## E3. Continuity of the frequency-magnitude distribution (31 October)

Purpose: remove a 28 percent step in the rate between sub-patch and
grown events. Causes: one-patch ruptures could carry less moment than
the patch magnitude because of stress-drop scatter, sizes jumped by
whole patches, and a one-cell crack carries 1.7 times the moment of a
resolved crack of the same area. Fix: sub-patch events from the same
size hazard, continuous size within the last increment, and the hazard
evaluated on the realised moment versus patch-count relation. Result:
N(>=M)/GR within a few percent across the patch magnitude; now guarded by
a regression test.

## E4. b-value estimator (31 October)

Purpose: explain a running b of 0.9 when 1 was imposed. The Utsu
half-bin correction had been applied to continuous synthetic
magnitudes (bias about -10 percent) and the upper truncation of the
sub-patch range had been ignored (bias upward). Fixed (continuous Aki
estimator, doubly truncated likelihood). Running b in time then
fluctuates at the sampling level (standard deviation 0.022 for windows
of 2000 events, expected 0.022) and does not correlate with the moment
balance (|r| < 0.2), as expected for a state-independent size hazard.
This result predates E2 and should be repeated.

## E5. Regime: arrest scale, stress-drop scatter, nucleation mode (1 November)

Purpose: find which parameters control the large-event statistics after
the calibration fix. Settings: branching ratio 0.8; uniform nucleation
with the slip cap and global feedback, or deficit-weighted nucleation
without either. Specs: `e5_regime_uniform.json`, `e5_regime_deficit.json`,
`e5_regime_scatter.json`.

| nucleation | rho_c | ln-sigma | max M | Bath | Omori p | site CV 0.3 m | site CV 1 m | M>=7 CV | closure | below mag_min |
|---|---|---|---|---|---|---|---|---|---|---|
| uniform | 0.25 | 0.25 | 7.53 | 1.07 | 0.97 | 0.85 | 0.73 | 0.79 | 0.97 | 0.06 |
| uniform | 1.0 | 0.25 | 7.53 | 1.04 | 0.96 | 0.82 | 0.70 | 0.74 | 0.98 | 0.03 |
| deficit | 0.25 | 0.0 | 7.50 | 1.06 | 0.96 | 0.94 | 0.80 | 0.74 | 1.00 | 0 |
| deficit | 0.25 | 0.25 | 7.59 | 1.05 | 0.97 | 0.87 | 0.76 | 0.80 | 1.01 | 0 |
| deficit | 0.25 | 0.5 | 7.58 | 1.04 | 0.96 | 0.84 | 0.67 | 0.78 | 0.98 | 0 |
| deficit | 1.0 | 0.25 | 7.49 | 1.01 | 0.96 | 0.83 | 0.67 | 0.90 | 0.99 | 0 |

Conclusion: these parameters change recurrence CVs by at most about 0.1.
Deficit weighting closes the budget without feedback or cap and
produces no events below mag_min, so it became the default.

## E6. Branching ratio (1 and 2 November)

Purpose: test whether ETAS clustering controls recurrence regularity and
Bath's gap. Defaults otherwise. Spec: `e6_branching.json`.

| branching ratio | Bath | site CV 1 m | M>=7 CV |
|---|---|---|---|
| 0 | 1.82 | 0.57 | 0.60 |
| 0.3 | 1.40 | 0.52 | 0.71 |
| 0.5 | 1.25 | 0.62 | 0.82 |
| 0.65 | 1.13 | 0.68 | 0.88 |
| 0.8 | 1.05 | 0.76 | 0.80 |

Conclusion: the branching ratio sets Bath's gap and the recurrence CVs
together and in opposite directions. 0.5 gives a gap near the observed
1.2 and a site CV near 0.6; it is the default. Without clustering the
window-selected "aftershocks" are background events and the Omori fit
is meaningless.

## E7. Fault length (2 November)

Purpose: check that Bath and Omori do not depend on geometry, and find
why whole-fault M >= 7 recurrence is close to Poisson on long faults.
Defaults. Spec: `e7_fault_length.json`.

| length | max M | N(M>=7) | Bath | Omori p | site CV 1 m | M>=7 CV | M0-L slope above sat |
|---|---|---|---|---|---|---|---|
| 50 km | 7.09 | 3 | 1.22 | 0.95 | 0.54 | (too few) | |
| 100 km | 7.38 | 17 | 1.23 | 0.95 | 0.53 | 0.38 | 1.48 |
| 200 km | 7.51 | 35 | 1.25 | 0.96 | 0.62 | 0.82 | 1.26 |
| 400 km | 7.79 | 54 | 1.25 | 0.96 | 0.65 | 0.88 | 1.43 |

Conclusion: Bath and Omori p are independent of length. M >= 7 events
recur periodically where they nearly fill the fault; on long faults the
whole-fault sequence is a superposition of nearly independent sequences
along strike and approaches Poisson, while site CVs stay at 0.5 to 0.65.

## E8. Strong state control of arrest (2 November)

Purpose: test whether stronger budget control of rupture size makes
large events more periodic. arrest_scale 2, defaults otherwise. Spec:
`e8_arrest_scale.json`. Result: max M 7.31, b above saturation 1.58,
Bath 1.19, site CV 0.66, M >= 7 CV 0.79, M0-L slope 1.12. Conclusion:
no gain in periodicity; the large-magnitude end is truncated.

## E9. Renewal-like nucleation (2 November)

Purpose: test whether a sharper dependence of nucleation on loading
makes large events more periodic. nucleation_exponent 3 (intensity
proportional to deficit^3), defaults otherwise. Spec:
`e9_nucleation_exponent.json`. Result relative to gamma = 1: M >= 7 CV
0.58 (0.82), site CV 0.54 (0.62), Bath 1.23 (1.25), Omori p 0.95 (0.96),
M0-L slope above saturation 1.89 (1.26), max M 7.50 (7.51). Conclusion:
renewal-like initiation makes large events markedly more regular without
affecting the clustering statistics. It is an option, not the default,
because it adds a parameter that these statistics alone do not fix.

## E10. Constant initiation, state-only termination (2 November)

Purpose: test whether GR emerges when initiation is uniform and constant
in time and rupture termination depends only on the fault state.
Settings: uniform nucleation, no rate feedback, no slip cap, no ETAS,
size_hazard off, so a rupture stops only when no neighbour can carry the
slip implied by its current size. Single seed, 1500 yr. Specs:
`e10_state_only.json`, `e10_state_only_variants.json`.

| initiation rate (x calibrated) | variant | N(>=M)/GR at M5.5 / 6.0 / 6.5 / 7.0 | characteristic magnitude |
|---|---|---|---|
| 1 | rho_c 0.25 | 6.3 / 1.5 / 0.02 / 0.05 | about 5.7 |
| 0.1 | rho_c 0.25 | 11 / 24 / 20 / 0.5 | 6 to 6.5 |
| 0.01 | rho_c 0.25 | 16 / 47 / 104 / 127 | 7 to 7.3 |
| 0.1 | rho_c 2 | 20 / 44 / 7.2 / 0 | about 6 |
| 0.1 | rho_c 0.25, ETAS n = 0.5 | 4.4 / 6.8 / 4.0 / 0.12 | about 6 |

The budget closes in all cases (closure 1.00 +- 0.04). Conclusion: the
size distribution is characteristic, the imposed GR below the patch
scale plus a bump whose magnitude is set by the initiation rate.
Broad, percolation-like admission and clustering narrow the bump but do
not remove it. With smooth, crack-like slip the deficit field does not
carry the scale-invariant heterogeneity that a power law would need;
compare the contrast between characteristic events on uniform faults
and GR on heterogeneous ones in fault models (for example Ben-Zion and
Rice, 1993). The mode is kept for the heterogeneity test below; the
default retains the scale-free size hazard, which gives GR.

## E11. Omori decay (1 November)

Purpose: explain an observer-view Omori p of about 0.95 for an input of
1.1. Single seed, 1500 yr. Direct offspring of M >= 6 mainshocks give p
= 1.09 to 1.11, all generations 0.91 to 0.96, and the window-selected
stack 0.81 to 0.98. Conclusion: secondary triggering (Helmstetter and
Sornette, 2002) and background in the window explain the difference.

## E12. Stress-drop and patch-scale roughness (31 October, pre-E2)

Slip CV 0.6 instead of 0.75 reduced the backslip-stress standard
deviation only from about 20 to 18 MPa; most of the heterogeneity comes
from the deficit pattern. Needs repeating with the current code before
any conclusion.

## Remaining single-fault tests

In order of priority; details in `ROADMAP.md`.

1. Quenched, scale-invariant strength heterogeneity with state-only
   termination (follow-up to E10): does GR emerge, and with what b?
2. Resolution: patch size 0.5 and 2 km with the same physical
   parameters; the patch magnitude moves, and the statistics above it
   should not.
3. Initial-condition and burn-in sensitivity: zero versus random initial
   deficit, and runs of 10,000 yr to check stationarity of the CVs.
4. Segmented fault: step-overs of 1 to 5 km, multi-segment rupture
   fraction versus gap, and recurrence on each segment.
5. Dip-slip and megathrust-like geometry: 30 degree dip, width 100 km,
   buried versus trench-breaking; the free-surface image is only
   approximate for in-plane slip.
6. Moment-length scaling above saturation (L^1.3 to 1.9 versus L^2):
   effective depth of slip below the locked zone.
7. Repeat E4 (time-dependent b) and E12 (stress heterogeneity) with the
   current code; settle the foreshock definition and measure foreshock
   fractions.
8. Longer and more seeds for the defaults, to give error bars on every
   number in `HANDOFF.md`.
