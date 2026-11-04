# Handoff: stochseis, state as of 31 October 2026

Repository: https://github.com/thwbecker/stochseis. All experiments are
logged in `docs/experiments.md`. Written by Claude
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

## 3. Current behaviour of the defaults (200 x 15 km, 20 mm/yr)

Reference values with error bars from 20 seeds of 10,000 yr (E22 in
docs/experiments.md): b = 1.00 +- 0.005 below width saturation and
1.28 +- 0.07 above, maximum M 7.58, Bath 1.24 +- 0.015, Omori p 0.96 +-
0.01 (observer view), site recurrence CV 0.62 +- 0.03 for slip >= 1 m,
whole-fault M >= 7 CV 0.74 +- 0.05, foreshock fraction 0.33 against a
chance level of 0.15 (10 days, 10 km), closure 1.002 +- 0.008, realised
branching ratio 0.41 for a nominal 0.5. Runtime about 12 min per
10,000 yr on one core.

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
2. **Single-fault tests are essentially complete** (docs/experiments.md
   E1 to E21): resolution independence, stationarity, segmented and
   dipping faults, emergence of GR from state (negative: GR must be
   imposed), deep slip and moment scaling, foreshocks, half-space check
   (exact for vertical faults, up to 35 percent too little slip for
   shallow thrusts), lower stress drop. Remaining runs that need more
   cores are specified in docs/experiments.md. Stage 1 (fault networks)
   is next.
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

## 4b. Alpine Fault intercomparison (started 1 November 2026)

Update after ABC round 4: hybrid sizes (scale-free below a threshold
area, state-only termination above) with slip capped at the local deficit
match mean recurrence, CV and single-event displacement at all three
paleoseismic sites (55 of 100 vectors accepted); details in docs/alpine.md.
Round 3 (slip-predictable) matches the site records but has unphysical
seismicity (b about 0, Bath 0.7) and is rejected once b and Bath are
added to the acceptance (`--seismicity`); round 4 keeps 54 of 100.
Low site CVs come with large Bath gaps (1.5 to 1.9) in this model class.
Variable-dip round 4: same acceptance (61 / 69 of the same 100 vectors);
geometry shows mainly in co-rupture between sites (much higher on the
variable-dip geometry), not in the site statistics. Next: posterior
ensemble with catalogs; L1 dip correction before quoting geometry effects
on slip.

Current step (3 November): L1 dip correction implemented
(stochseis/dipcorr.py); full campaign over five configurations (planar,
variable dip, segmented transition on 0.5 km meshes with 1 km patches;
planar and variable dip on 1 km meshes with 2 km patches) driven by
examples/alpine/run_full_campaign.sh, with RSQSim-style site analysis
(site_analysis.py), overview figures and a cross-configuration summary
(campaign_summary.py). See docs/alpine.md, "Full campaign".
Results (docs/alpine.md, "Full campaign results"): 52 to 65 of 100
vectors accepted in every configuration with strong overlap; accepted
parameters close to the round-5 priors; south and central means matched,
central CV and north RI biased high in all ensembles (RSQSim has the
north bias too); site statistics converged at 2 km and insensitive to
geometry. Posterior rescoring: all 20 posterior vectors stay accepted in
every configuration and the best was ranked 1 or 2 by the ABC. The
product for the intercomparison is the calibrated envelope (docs/alpine.md);
physics models are compared against it, RSQSim being one example
(site_analysis.py --reference). Figures: event statistics independent of
geometry; co-rupture differs by at most 0.07 between geometries (the
round-4 co-rupture contrast is not reproduced); hypocentres at the base;
slip-predictable sites; north RI strongly threshold dependent
(`--site north=1` fits RI and SED but raises the north CV to 0.5).
Campaign rescored with north at 1 m: 64 to 69 accepted, best vectors
unchanged, north RI / SED of the best vectors 362 to 436 yr / 3.5 to
4.2 m; the only remaining systematic misfit is regularity at the central
and north sites (CV 0.23 and 0.25 observed, about 0.4 in the ensembles).
Bath gap too large (1.4 ensembles): controlled by b, area threshold and
branching ratio; round 6 (abc6_bath.json) widens branching (0.2-0.7) and
b (0.9-1.2); new metric bath_mean_m7; M8 through-going ruptures keep gaps
near 2 on a single fault (no loaded area left for large aftershocks).
The paragraph below describes the starting point.

See docs/alpine.md and docs/networks.md. STL geometries load as section
grids (stochseis/mesh.py), paleo-site statistics are computed by
stochseis/sites.py and included in sweep summaries when sites_file is
set. First runs with the defaults give single-event displacements of 3 to
5 m against 7 m observed and recurrence CVs of 0.5 to 0.8 against 0.2 to
0.4. Missing inputs: segment slip rates, truncated and Alpine-Hope
geometries.

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

Round 6 (3 November, coarse): branching 0.2-0.7 and b 0.9-1.2 bring the
Bath gap to about 1.0 (M >= 6) and 1.24-1.38 (M >= 7, best vectors) but
raise site CVs (best central 0.42-0.56) and lower acceptance (41-54):
regular sites and a Bath gap of 1.2 are in tension in this model class.
Option site_merge_yr / site_analysis.py --merge merges site events closer
than the dating resolution (modest effect, CV -0.05 to -0.1).
Decision: regularity over Bath. Reference baseline = round 5 (north at
1 m), Bath gap reported as an outcome; round 6 = sensitivity case.
Merging site events closer than 30 yr on the round-5 fine posteriors:
central and south CVs -0.03 to -0.04 (ensemble medians 0.34 to 0.37),
north unchanged; reference envelope in docs/alpine.md.
Round 7 (regularity priors, coarse): best variable-dip vector matches
central RI/CV/SED (264 yr / 0.23 / 7.0 m) and north RI/SED; south too
regular (0.24 vs 0.40), north CV 0.32 vs 0.25. Cost: strongly
characteristic, partly bimodal size distribution, Bath about 1.7 overall.
site_analysis.py now uses the catalogs' site_merge_yr by default.
