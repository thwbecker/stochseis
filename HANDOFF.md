# Handoff: stochseis, state as of 31 October 2026

Repository: https://github.com/thwbecker/stochseis. Written by Claude
(Anthropic) in a chat session with Thorsten Becker; see README.md.

Purpose: a stochastic seismicity generator on fault patch graphs,
intended as the baseline that physics-based earthquake cycle models have
to beat. Owner: Thorsten Becker. This file records what exists, which
decisions were made and why, what was tried and rejected, and what to do
next. Read `docs/design.md` first for the model itself.

## 1. What exists

* Package `stochseis` (Python, numpy/scipy/matplotlib only), 18 pytest
  tests passing, documented in `README.md`, `docs/design.md`,
  `docs/parameters.md`, `docs/analysis.md`.
* Stage 0 geometry: planar, optionally segmented faults with step-overs,
  on the abstract patch-graph interface that later loaders (UCERF3 or
  NZ CFM sections, stochastic fault networks) must produce.
* Physics: spectral static stiffness (validated against Eshelby and the
  antiplane crack to 1 percent), crack solutions for uniform stress drop
  on arbitrary outlines (preconditioned CG, cached by shape), von Karman
  lognormal roughness, deficit-conditioned slip, rupture flood growth
  with percolation-type arrest, scale-free size hazard, ETAS clustering,
  backslip stress state, rupture kinematics (speed, onset times,
  duration), budget-calibrated nucleation rate with weak global feedback.
* Analysis: FMD and b by magnitude regime, source scaling, Omori MLE
  stacks, Bath (observer and genealogy views), foreshock fraction with
  chance level, whole-fault and site recurrence CVs, moment closure,
  stress heterogeneity, arrest statistics; summary and time-series
  figures; CLI with run / analyze / sweep (multiprocessing).

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

## 3. Current behaviour of the defaults (200 x 15 km, 20 mm/yr, 3000 yr, seed 1)

Defaults as of 1 November 2026: deficit-weighted nucleation, no slip cap,
branching ratio 0.5, stress drop 3 MPa with ln-scatter 0.25, b = 1,
arrest_scale 0.25, slip CV 0.75, moment hazard. Figures in
`runs/defaults_2026-11-01/`.

* FMD on the GR line from M3 to M7.5 (b = 1.00 in the sub-patch range
  and below width saturation, 1.35 above it), maximum M7.52, 27 events
  above M7, none below mag_min.
* Bath 1.25, Omori p 0.96 (observer view), site recurrence CV 0.74
  (slip >= 0.3 m) and 0.62 (>= 1 m), whole-fault M >= 7 CV 0.81,
  closure 0.98, mean deficit 2.1 m, M0 versus L slope above saturation
  1.76.
* Runtime: about 70 s per 3000 years, 330,000 events, one core.
* Multi-seed statistics for neighbouring settings are in ROADMAP step 1
  and `runs/regime_2026-11-01/`.

Results quoted before 1 November 2026, the figures in `runs/sweep1`,
`runs/sweep2` and the top level of `runs/`, and the reasoning that led to
the earlier defaults were produced with the miscalibrated rate.

Note on b-values: before 31 October the estimator applied the Utsu
half-bin correction to continuous synthetic magnitudes (b biased about
10 percent low) and ignored upper truncation in the regime fits (biased
high). Both are fixed in `analysis.b_value_mle`; b numbers in earlier
figures and messages are affected. Running b in time fluctuates at the
sampling-noise level and is uncorrelated with the moment balance,
because the size hazard is state independent; a state-dependent hazard
would be needed to produce stress-dependent b variations.

## 4. Open issues, in order of importance

1. **Finish the regime study** (ROADMAP step 1, remaining items): more
   seeds and longer runs for the chosen defaults, arrest_scale 2, and the
   dependence of site CV and Bath on fault length; the trade-off between
   Bath and site CV through the branching ratio is the main result so far
   and should be checked for robustness.
2. **GR from state alone**: the state-only termination mode gives
   characteristic distributions with smooth crack-like slip; test quenched
   scale-invariant strength heterogeneity (ROADMAP step 3). The question
   of whole-fault M >= 7 CV is resolved: superposition along strike and
   clustering, see ROADMAP.
3. **Stress as state.** Patch-scale stress roughness (std about 18 MPa)
   is probably too large; it scales with slip_cv and the deficit
   conditioning. Options: filter the roughness at the patch scale,
   lower slip_cv, or use a low-pass stress for arrest. Needed before the
   maturity/damage variable (recorded, unused) can be given a role.
4. **M0 versus L above saturation** goes as L^1.2 to 1.8 against the
   observed L^2 (W-model limitation). One extra parameter (effective
   depth extent of slip below the locked zone) would test the fix.
5. **Omori p** (largely resolved): in a 1500-yr run the direct offspring
   of M >= 6 mainshocks recover the input p = 1.10; secondary
   triggering flattens the all-generation stack to 0.96 (the known ETAS
   effect, Helmstetter and Sornette, 2002), and the background inside the
   observer window lowers it further to 0.93 to 0.98. Comparisons with
   data should use the observer view and fit p jointly with a background
   term.
6. **Foreshock statistic**: defined as an event of M >= M_main - 2
   within 30 days and max(10 km, L); chance level is computed, but the
   definition should be settled with the user.
7. Minor: single-core sandbox prevented completing the six-run sweep in
   one go; the sweep command itself works. Capped events can fall below
   mag_min (analysis filters on mag_min). The free-surface image is
   exact for antiplane components only. Sub-patch events carry 2 percent
   of the moment.

## 5. Next stages (from the original plan; concrete steps in ROADMAP.md)

Stage 1: UCERF3 or NZ CFM subsection loader producing the patch graph
(sections, connections with jump distances, slip rates, coupling); the
backslip stress state then needs the general-geometry path (per-section
grids or Okada via interact). Stage 2: megathrust meshes (triangles,
dipping, 100 km widths). Stage 3: recurrence-rate calibration mode
(adjust per-section nucleation weights to paleoseismic recurrence).
Later: stochastic fault-network generator on the same interface, with
slip rates constrained to a transect loading rate.

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
