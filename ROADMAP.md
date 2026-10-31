# Roadmap

Concrete next steps, in the order proposed in `HANDOFF.md`. Each step
lists what to do, what is already in place, and how to judge the
outcome. Steps 1 and 2 are independent of each other; step 2 changes
the regime that step 1 maps, so if time is short, do step 2 first and
run the sweep once.

## 1. Map the budget-control regime

Question: how do the large-magnitude statistics depend on how strongly
the slip budget controls arrest, and which regime is the better
baseline?

* Ready: `examples/sweep_regime.json` (arrest_scale 0.1 to 1.0 by
  stress-drop scatter 0 to 0.5, three seeds, 3000 yr; 36 runs, about
  2 min each on one core). Run with
  `python -m stochseis.cli sweep examples/sweep_regime.json sweep_regime --nproc N --no-catalogs`
  and summarise with
  `python examples/aggregate_sweep.py sweep_regime/sweep_summary.json p_arrest_scale p_stress_drop_ln_sigma`.
* Look at jointly: maximum magnitude and number of M >= 7, FMD ratio to
  GR above M6.5, b above saturation, site recurrence CV, Bath, fraction
  of M >= 6 events arrested by the state, mean deficit.
* Outcome: a choice of default arrest_scale and stress-drop scatter with
  seed-to-seed spread, recorded in `docs/parameters.md`.

## 2. Nucleation weighted by slip deficit

Problem: nucleation is state independent, so 10 to 20 percent of
nucleations fall on relaxed patches and are capped below mag_min by
`cap_by_deficit`, and a global rate feedback is needed to close the
budget.

* Proposal (no new parameter): sample background nucleation patches with
  probability proportional to positive deficit; accept triggered
  nucleations with probability min(1, deficit / delta_ref) at the
  offspring patch. Recalibrate the ETAS productivity K0 for the
  acceptance rate (measure it in the pilot or iterate once).
* Then test whether `rate_feedback` and `cap_by_deficit` can be switched
  off. Keep both as options for comparison.
* Acceptance: fraction of events below mag_min near zero; moment closure
  within 5 percent over 3000 yr without feedback; FMD continuity at the
  patch magnitude unchanged; Bath and Omori statistics unchanged within
  seed spread; tests extended accordingly.
* Code: `Simulation.run` (background draw and offspring acceptance),
  `Simulation._setup` (K0), `params.py` (switch), tests in
  `tests/test_growth_etas_sim.py`.

## 3. Stress heterogeneity and a fault-state variable

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
