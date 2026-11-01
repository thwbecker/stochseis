# Roadmap

Concrete next steps, in the order proposed in `HANDOFF.md`. Each step
lists what to do, what is already in place, and how to judge the
outcome.

Status 1 November 2026: step 2 is implemented and is the default; a
calibration error that had distorted all earlier statistics at the
large-magnitude end was found and fixed on the way (`docs/design.md`
section 2.5); a reduced version of step 1 has been run and the defaults
were updated from it. Steps 3 to 5 are next.

## 1. Map the budget-control regime

Question: how do the large-magnitude statistics depend on how strongly
the slip budget controls arrest, and which regime is the better
baseline?

* Ready: `examples/sweep_regime.json` (arrest_scale 0.25 to 2 by
  nucleation weighting uniform or deficit by stress-drop scatter 0.25 or
  0.5, three seeds, 3000 yr; 48 runs, about 2 min each on one core). In
  the deficit-weighted runs the slip cap stays on because the sweep grid
  is a full product; with weighting it rarely binds (set
  cap_by_deficit = false in a follow-up if the weighted mode is chosen). Run with
  `python -m stochseis.cli sweep examples/sweep_regime.json sweep_regime --nproc N --no-catalogs`
  and summarise with
  `python examples/aggregate_sweep.py sweep_regime/sweep_summary.json p_arrest_scale p_stress_drop_ln_sigma`.
* Look at jointly: maximum magnitude and number of M >= 7, FMD ratio to
  GR above M6.5, b above saturation, site recurrence CV, Bath, fraction
  of M >= 6 events arrested by the state, mean deficit.
* Outcome: a choice of nucleation weighting, arrest_scale and stress-drop
  scatter with seed-to-seed spread, recorded in `docs/parameters.md`.
  Single-seed site-CV estimates have a standard error of about 0.2, so
  three seeds of 3000 yr is the minimum; more seeds for the finalists.
* Done in reduced form (24 runs instead of 48; no arrest_scale 2 and no
  uniform runs with scatter 0.5). Summaries in `runs/regime_2026-11-01/`:

  | nucleation | rho_c | ln-sigma | n | max M | N(M>=7) | Bath | Omori p | site CV 0.3 m | site CV 1 m | M>=7 CV | closure | below mag_min |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | uniform | 0.25 | 0.25 | 0.8 | 7.53 | 30 | 1.07 | 0.97 | 0.85 | 0.73 | 0.79 | 0.97 | 0.06 |
  | uniform | 1.0 | 0.25 | 0.8 | 7.53 | 31 | 1.04 | 0.96 | 0.82 | 0.70 | 0.74 | 0.98 | 0.03 |
  | deficit | 0.25 | 0.0 | 0.8 | 7.50 | 30 | 1.06 | 0.96 | 0.94 | 0.80 | 0.74 | 1.00 | 0 |
  | deficit | 0.25 | 0.25 | 0.8 | 7.59 | 30 | 1.05 | 0.97 | 0.87 | 0.76 | 0.80 | 1.01 | 0 |
  | deficit | 0.25 | 0.5 | 0.8 | 7.58 | 31 | 1.04 | 0.96 | 0.84 | 0.67 | 0.78 | 0.98 | 0 |
  | deficit | 1.0 | 0.25 | 0.8 | 7.49 | 36 | 1.01 | 0.96 | 0.83 | 0.67 | 0.90 | 0.99 | 0 |
  | deficit | 0.25 | 0.25 | 0.0 | 7.58 | 31 | 1.82 | n/a | 0.66 | 0.57 | 0.60 | 0.99 | 0 |
  | deficit | 0.25 | 0.25 | 0.5 | 7.51 | 35 | 1.25 | 0.96 | 0.73 | 0.62 | 0.82 | 0.98 | 0 |

  Means over seeds 1 to 3, 3000 yr each, 200 x 15 km fault, b = 1.00 below
  width saturation in all runs. Seed-to-seed standard deviations are 0.01
  to 0.1 for the CVs and Bath and are in the JSON files.
* Findings. Arrest scale, stress-drop scatter and the nucleation mode
  change the recurrence CV by at most 0.1. The ETAS branching ratio
  controls it: site CV for slip >= 1 m goes from 0.57 without clustering
  to 0.62 at n = 0.5 and 0.76 at n = 0.8, because triggered large events
  re-rupture sites shortly after a mainshock. Bath's gap moves the other
  way (1.82, 1.25, 1.05 in the window-selection view). n = 0.5 gives both
  a Bath gap near the observed 1.2 and a site CV near 0.6.
* Adopted as defaults: nucleation_weighting = "deficit",
  cap_by_deficit = false, branching_ratio = 0.5; stress-drop scatter 0.25
  and arrest_scale 0.25 kept. A 3000-yr run with these defaults is shown
  in `runs/defaults_2026-11-01/`.
* Remaining for this step: more seeds and longer runs for the finalists,
  arrest_scale 2, and the dependence on fault length (the CV and Bath
  values above are for one fault geometry).

### Further tests, 2 November 2026

Defaults unless stated (deficit-weighted nucleation, branching ratio
0.5), three seeds of 3000 yr per row; summaries in
`runs/tests_2026-11-02/`.

| test | max M | Bath | site CV 1 m | whole-fault M>=7 CV | M0-L slope above sat |
|---|---|---|---|---|---|
| fault length 50 km | 7.09 | 1.22 | 0.54 | (too few) | |
| fault length 100 km | 7.38 | 1.23 | 0.53 | 0.38 | 1.48 |
| fault length 200 km (defaults) | 7.51 | 1.25 | 0.62 | 0.82 | 1.26 |
| fault length 400 km | 7.79 | 1.25 | 0.65 | 0.88 | 1.43 |
| branching 0 / 0.3 / 0.5 / 0.65 | 7.5 to 7.6 | 1.82 / 1.40 / 1.25 / 1.13 | 0.57 / 0.52 / 0.62 / 0.68 | 0.60 / 0.71 / 0.82 / 0.88 | |
| arrest_scale 2 | 7.31 | 1.19 | 0.66 | 0.79 | 1.12 |
| nucleation_exponent 3 | 7.50 | 1.23 | 0.54 | 0.58 | 1.89 |

* Bath and Omori p (0.95) do not depend on fault length.
* Whole-fault M>=7 recurrence is periodic (CV 0.38) when such events
  nearly fill the fault and close to Poisson on long faults: the
  irregularity comes from superposed, nearly independent rupture
  sequences along strike and from ETAS clustering, not from the size
  rule.
* A steeper, renewal-like nucleation dependence (deficit^3) lowers the
  M>=7 CV from 0.82 to 0.58 and the site CV to 0.54 without changing
  Bath or Omori.
* State-only termination (single seed, 1500 yr; `size_hazard = false`,
  constant uniform initiation, no feedback, no cap): the
  frequency-magnitude distribution becomes characteristic, a GR part
  below the patch scale (imposed) plus a bump whose magnitude is set by
  the initiation rate: about M5.7 at the calibrated rate, M6 to 6.5 at
  1/10 of it, M7 to 7.3 at 1/100. Broad admission (arrest_scale 2) and
  ETAS clustering (n = 0.5) narrow the bump but do not remove it. With
  the smooth, crack-like slip of this model the deficit field is not
  heterogeneous enough for a power law to emerge from state alone,
  consistent with the classic contrast between characteristic events on
  uniform faults and GR on heterogeneous ones (e.g. Ben-Zion and Rice,
  1993). A test of quenched, scale-invariant strength heterogeneity is
  the natural follow-up (step 3).

Single-fault tests E13 to E19 (resolution, stationarity, segmentation,
dipping fault, GR from state, deep slip, foreshocks) were run on 3
November; see `docs/experiments.md`, which also lists the few remaining
single-fault items. Stage 1 (step 5 below) is next.

## 2. Nucleation weighted by slip deficit (implemented as an option)

Problem: nucleation is state independent, so 10 to 20 percent of
nucleations fall on relaxed patches and are capped below mag_min by
`cap_by_deficit`, and a global rate feedback is needed to close the
budget.

* Done: `nucleation_weighting = "deficit"`. Background intensity is
  proportional to positive deficit (thinning against an adaptive bound);
  triggered events keep their number and are placed at the most loaded
  of eight kernel candidates, chosen with probability proportional to
  deficit. The first design (thinning triggered events and rescaling
  productivity) was rejected because it suppressed the productivity of
  large events (design.md section 2.6).
* Against the acceptance criteria (1500 yr, seed 1, cap off): events
  below mag_min 0 (uniform: 3 to 8 percent); closure 1.006 after burn-in
  without feedback; FMD continuous at the patch magnitude; Bath 1.05 and
  Omori p 0.95 against 0.97 and 0.98 for uniform. All met.
* The high site CV in the first single run was seed noise; over three
  seeds deficit-weighted and uniform nucleation give the same CVs within
  0.05. Now the default (step 1).

## 3. Stress heterogeneity and a fault-state variable

* New first item: quenched, scale-invariant heterogeneity in strength
  (for example a von Karman field multiplying the slip a patch must
  carry before it is admitted) in the state-only termination mode, to
  test whether GR emerges from heterogeneity plus state rather than from
  the imposed size hazard.

* Measure the backslip-stress standard deviation (currently about
  18 MPa at 1 km patches) as a function of slip_cv, condition_weight and
  patch size; decide on a patch-scale filter of the roughness or a
  low-pass stress for state-dependent rules.
* Give the recorded maturity (cumulative slip) a role only if a clear
  observable motivates it, for example smoother slip or lower stress
  drop on mature faults, or a state-dependent size hazard that would
  produce stress-dependent b-value variations (currently b fluctuates
  only at the sampling-noise level).

## 4. Moment-length scaling above width saturation

M0 grows as L^1.2 to 1.8 against the observed L^2. Test one extra
parameter, an effective depth extent of coseismic slip below the locked
width (Shaw, 2009, and related arguments), implemented by extending the
crack-solution grid below the fault with zero loading there. Compare
mean slip versus length against the Blaser et al. (2010) data.

## 5. Stage 1: fault networks

* Loader for UCERF3 fault subsections (OpenSHA fault model and
  deformation model files; trace, dip, depths, rake, slip rate, aseismic
  fraction) and for the New Zealand community fault model, producing the
  `Fault` patch graph: one grid per section, connections between section
  ends with jump distances (UCERF3 plausibility rules as a default).
* Slip on multi-section ruptures: unroll sections along the rupture path
  into one (s, z) plane for the random field and crack solution, with an
  optional correlation penalty across section boundaries.
* State stress across sections: either neglect inter-section stress
  transfer at first (deficit is the default state measure, so this is
  acceptable) or compute it with Okada kernels, possibly through the
  interact code.
* Geometry generalisation required in `geometry.py` (non-planar
  coordinates), `state.py` (per-section stress grids) and
  `slip.unit_crack` (per-section dip and slip direction).

## 6. Later

* Megathrust meshes (triangulated, dipping, widths near 100 km); the
  KL/lognormal construction of LeVeque et al. (2016) and Melgar et al.
  (2016) on triangles as an alternative to the FFT field.
* Recurrence-rate calibration mode: adjust per-section nucleation
  weights to match paleoseismic recurrence at given sites.
* Stochastic fault-network generator on the same interface, with slip
  rates constrained to a transect loading rate.
* Performance: numba or C for the flood loop in `rupture.grow_rupture`
  if networks with 10^5 or more patches are needed.

## Housekeeping

* Choose a license.
* Check the references in `docs/design.md` against erdbeben.bib.
* Line-by-line human review of the core modules (`stiffness.py`,
  `slip.py`, `rupture.py`, `simulator.py`), which have so far been
  checked only through the tests and the diagnostics.
* The GitHub Actions workflow in `.github/workflows/tests.yml` runs the
  test suite on push.
