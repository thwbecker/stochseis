# Handoff: stochseis, state as of 5 November 2026

Repository: https://github.com/thwbecker/stochseis. All experiments are
logged in `docs/experiments.md`. Written by Claude
(Anthropic) in a chat session with Thorsten Becker; see README.md.

Purpose: a stochastic seismicity generator on fault patch graphs,
intended as the baseline that physics-based earthquake cycle models have
to beat. Owner: Thorsten Becker. This file records what exists, which
decisions were made and why, what was tried and rejected, and what to do
next. Read `docs/design.md` first for the model itself.

## 1. What exists

* Package `stochseis` (Python, numpy/scipy/matplotlib only), 47 pytest
  tests, described in `README.md` (what the code does, what was
  learned, layout) and `docs/` (design, parameters, analysis, the
  development log `experiments.md`, network design).
* Geometries: planar and segmented faults, and triangulated surfaces
  (STL) resampled onto section grids (the Alpine Fault meshes), with
  segment slip rates and a dip correction of crack amplitudes.
* Rules: budget-calibrated, deficit-weighted nucleation; ETAS
  clustering; rupture flood growth with deficit-dependent admission and
  a scale-free or hybrid size rule; crack slip set by stress drop, with
  roughness, optionally capped at the local deficit; optional quenched
  structure (strength field, segment-boundary barriers, hierarchical
  patches).
* Analysis: FMD and b by regime, source scaling, Omori stacks, Bath,
  foreshocks, recurrence at points and on the whole fault, synthetic
  paleoseismic records, record matching against observations,
  constraint-weighted predictions; CLI with run / analyze / sweep.
* Two worked examples, each with a README, results and `run_test.sh`:
  `examples/single_fault/` (parameter overview and the development
  experiments E1 to E23) and `examples/nz_alpine/` (Alpine Fault
  simulator comparison).

## 2. Decisions and their reasons

1. Slip amplitude is the static crack solution for uniform stress drop on
   the realised outline, times mean-one roughness and state weights.
   Normalising a tapered random field to a moment-weighted stress drop
   was rejected: edge singularities dominate that measure and the mean
   slip came out at half the crack value.
2. Arrest is local admission with contiguity (flood growth); a region in
   which the admission probability falls below the lattice percolation
   threshold cannot be crossed. Rejected before it: a front-averaged
   hazard (ruptures crossed barriers because the average mixes azimuths)
   and admission without contiguity (ruptures percolated as sparse,
   hole-ridden patterns across semi-depleted fault). Enclosed holes are
   filled after growth.
3. The deficit is the state measure for arrest and conditioning; the
   backslip stress tau = K deficit is tracked for diagnostics and is
   available as an option, but rough slip leaves tens of MPa of
   patch-scale stress heterogeneity, and as an arrest measure it let an
   unphysical fault-wide rupture through in an early run.
4. Nucleation rate: calibrated from a pilot ensemble so that expected
   release equals loading, with the background rate scaled by the
   fault-mean deficit. Without feedback the deficit drifts because
   shadowed nucleations release less than the pilot assumes; feedback on
   triggered events was removed because it switched aftershock sequences
   off after large events.
5. Sizes are continuous below the patch (sub-patch events from the same
   hazard, with budget debit) and within the last patch increment, and
   the moment hazard runs on the realised moment versus patch-count
   relation from the pilot. This removed a 28 percent step in the FMD
   at the patch magnitude; see docs/design.md section 2.2.
6. Event mean slip is capped at the mean positive deficit on the rupture
   (`cap_by_deficit`). Without it, once the nucleation rate is high the
   forced release by nucleations in relaxed zones drives the fault to a
   near-zero deficit in which large events cannot percolate.
7. Defaults (provisional, see docs/parameters.md): stress drop 3 MPa with ln-scatter 0.25 (0.5 gave site
   recurrence CVs near 0.8; 0.25 gives about 0.5 to 0.6, the
   paleoseismic range), b = 1, rho_c = 0.25, slip CV 0.75, moment
   hazard (b continuous through width saturation), branching ratio 0.8.

8. Rate calibration for the moment hazard is analytic (design.md
   section 2.5). The Monte Carlo pilot mean it replaces was dominated by
   rare large events and typically several times too small, which made
   the nucleation rate too high and, through the feedback, held the fault
   at a low mean deficit with too few large events. This affected every
   result before 1 November 2026.
9. Deficit-weighted nucleation with productivity-preserving placement of
   triggered events is the default; it closes the budget without the
   global feedback and the slip cap (both kept as options for uniform
   nucleation).
10. Branching ratio 0.5: the reduced regime sweep showed that the ETAS
   branching ratio, not arrest or stress-drop scatter, sets the site
   recurrence CV and Bath's gap, in opposite directions; 0.5 gives about
   0.6 and 1.25.
11. Options added on 2 November: ``nucleation_exponent`` (renewal-like
   sharpening; 3 lowers the M>=7 CV from 0.82 to 0.58) and
   ``size_hazard = False`` (state-only termination; yields characteristic
   size distributions). Neither is a default. The background sampler was
   rewritten (thinning on the fault-mean weight, then weighted patch
   choice) because thinning on the fault maximum was too slow for
   gamma > 1. Test results are in ROADMAP step 1.

## 3. Current behaviour of the defaults (200 x 15 km, 20 mm/yr)

Reference values with error bars from 20 seeds of 10,000 yr (E22 in
docs/experiments.md): b = 1.00 +- 0.005 below width saturation and
1.28 +- 0.07 above, maximum M 7.58, Bath 1.24 +- 0.015, Omori p 0.96 +-
0.01 (observer view), site recurrence CV 0.62 +- 0.03 for slip >= 1 m,
whole-fault M >= 7 CV 0.74 +- 0.05, foreshock fraction 0.33 against a
chance level of 0.15 (10 days, 10 km), closure 1.002 +- 0.008, realised
branching ratio 0.41 for a nominal 0.5. Runtime about 12 min per
10,000 yr on one core.

Results quoted before 1 November 2026 and the reasoning that led to the
earlier defaults were produced with the miscalibrated rate (the figures
from that period were removed from the repository on 5 November 2026).

Note on b-values: before 31 October the estimator applied the Utsu
half-bin correction to continuous synthetic magnitudes (b biased about
10 percent low) and ignored upper truncation in the regime fits (biased
high). Both are fixed in `analysis.b_value_mle`; b numbers in earlier
figures and messages are affected. Running b in time fluctuates at the
sampling-noise level and is uncorrelated with the moment balance,
because the size hazard is state independent; a state-dependent hazard
would be needed to produce stress-dependent b variations.

## 4. Open issues, in order of importance

1. **South-central contrast on the Alpine Fault.** No control tested
   (nucleation clock, size rule, strength heterogeneity, first test of
   segment-boundary barriers) gives a regular central and a less
   regular south record together; set D (`examples/nz_alpine/run_test.sh`)
   is the next test.
2. **Regularity versus aftershock statistics.** Regular recurrence comes
   with characteristic sizes and weak triggering (Bath gaps of 1.3 and
   more); this is a property of the model class to be stated in the
   comparison, not a bug.
3. **Hierarchical patches (E23).** Implemented, first runs only; strong
   barriers need about 3000 yr of spin-up. Decision (5 November): keep
   the barrier inside the admission ratio (no extra parameter) and handle
   spin-up by run length and the analysis burn-in rather than a
   recalibration.
4. **Stress as state.** Patch-scale stress roughness is probably too
   large (scales with slip_cv); needed before a damage or maturity
   variable can be given a role.
5. **M0 versus L above saturation** goes as about L^1.5 unless slip
   extends below the seismogenic depth (`deep_slip_km`, E18).
6. **Kinematics.** Events have onset times but only final slip, so no
   moment-rate functions; a slip-time function with a local rise time
   would give them (relevant for cascade-type observations).
7. Minor: capped events can fall below mag_min (analysis filters on
   mag_min); the free-surface image is exact for antiplane components
   only; sub-patch events carry about 2 percent of the moment.

## 5. Next stages (concrete steps in ROADMAP.md)

Run the two `run_test.sh` scripts; the E23 hierarchy sweep on the single
fault; then Stage 1: a UCERF3 example (fault sections with connections
and jump distances, `docs/networks.md` stage N2 for multi-section
faults); later a fault system represented only probabilistically
(stochastic fault-network generator on the same patch-graph interface).

## 6. Practical notes

* `python -m stochseis.cli run - out.npz --set key=value ...` then
  `analyze out.npz --plots dir`; `Params.from_json` for configs; seeds
  make runs deterministic.
* `Simulation._build_m0_table` and `_pilot` run at construction; setup is
  a few seconds. `slip.CRACK_SOLVES` counts cache misses.
* Hot spots: `dijkstra` per grown rupture and the flood loop in
  `rupture.grow_rupture`; `FaultState` full-grid FFT per flush. A numba
  or C version of the flood loop is the first optimisation if networks
  with 10^5 patches are needed.
* References in docs/design.md were written from memory; verify against
  erdbeben.bib before citing.
