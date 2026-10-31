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
7. Defaults: stress drop 3 MPa with ln-scatter 0.25 (0.5 gave site
   recurrence CVs near 0.8; 0.25 gives about 0.5 to 0.6, the
   paleoseismic range), b = 1, rho_c = 0.25, slip CV 0.75, moment
   hazard (b continuous through width saturation), branching ratio 0.8.

## 3. Current behaviour of the defaults (200 x 15 km, 20 mm/yr, 1500 yr, seed 1)

* FMD flat relative to GR from M3 to M6 within 3 percent; b = 1.00 /
  0.98 / 1.7 in the sub-patch, grown-below-saturation and
  above-saturation ranges (corrected estimator, see below); above M6.5 the rate falls below GR (0.85 at
  6.5, 0.36 at 7.0) with maximum M7.1 and 11 events above M7.
* Bath 0.96, Omori p about 0.93 (input 1.1), site recurrence CV 0.57
  (slip >= 0.3 m) and 0.49 (>= 1 m), whole-fault M >= 7 CV 0.63, closure
  1.02, mean deficit 0.55 m.
* Runtime: about 60 s per 1500 years, 360 k events, one core.

Note on b-values: before 31 October the estimator applied the Utsu
half-bin correction to continuous synthetic magnitudes (b biased about
10 percent low) and ignored upper truncation in the regime fits (biased
high). Both are fixed in `analysis.b_value_mle`; b numbers in earlier
figures and messages are affected. Running b in time fluctuates at the
sampling-noise level and is uncorrelated with the moment balance,
because the size hazard is state independent; a state-dependent hazard
would be needed to produce stress-dependent b variations.

## 4. Open issues, in order of importance

1. **Large-magnitude end and the budget regime.** With the final scheme
   the fault runs at a lower mean deficit (0.5 m) than before the
   continuity fix (1.7 m) and releases its budget through more frequent
   M6 to M7 events; M7.5 fault-wide events, present earlier, are now
   rare. The percolation threshold (deficit above about 0.73 rho_c times
   the implied slip) acts as a size ceiling that self-organises with the
   budget. Whether this regime or the earlier one is the better baseline
   is a modelling question: redo the rho_c (0.1 to 0.5) by stress-drop
   scatter sweep with the final code, 3000 years or longer and several
   seeds, and look at max magnitude, b above saturation, site CV and the
   FMD ratio to GR together. `examples/sweep_arrest_stressdrop.json` is
   the template; `runs/cv` and `runs/sweep2` hold figures from
   intermediate code versions and should not be compared across versions.
2. **Nucleation location.** Nucleations are uniform in space (background)
   or follow the ETAS kernel (triggered), independent of state. Weighting
   nucleation by positive deficit (background sampling proportional to
   deficit, triggered candidates accepted with probability
   min(1, deficit/delta_ref)) is parameter-free, physically motivated
   (aftershocks concentrate at rupture peripheries and low-slip patches),
   would make the global rate feedback unnecessary, and would remove the
   need for the slip cap. It changes realised ETAS productivity, so K0
   would need recalibration on acceptance. Recommended next step. The
   size of the issue: in the default runs 10 to 20 percent of all
   nucleations are capped to below mag_min because they fall on relaxed
   patches; the analysis filters them out, but they should not exist.
3. **Stress as state.** Patch-scale stress roughness (std about 18 MPa)
   is probably too large; it scales with slip_cv and the deficit
   conditioning. Options: filter the roughness at the patch scale,
   lower slip_cv, or use a low-pass stress for arrest. Needed before the
   maturity/damage variable (recorded, unused) can be given a role.
4. **M0 versus L above saturation** goes as L^1.2 to 1.8 against the
   observed L^2 (W-model limitation). One extra parameter (effective
   depth extent of slip below the locked zone) would test the fix.
5. **Omori p** in the observer view is 0.92 to 0.97 for input 1.1;
   understand how much is selection versus secondary triggering before
   comparing to data.
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
