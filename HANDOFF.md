# Handoff: stochseis, state as of 6 November 2026

For whoever continues this work (a person or a model). The top-level
`README.md` states the purpose, philosophy and findings; this file
records where the work stands, why things are as they are, what is open,
and what comes next.

## 1. Where we are

* **Package.** `stochseis` (Python with numpy, scipy, matplotlib; no
  compiled code), 47 pytest tests. Planar, segmented and triangulated
  (STL) faults on 1 to 2 km patch grids; budget-calibrated,
  deficit-weighted nucleation; ETAS clustering; flood growth with
  deficit-dependent admission and a scale-free or hybrid size rule;
  crack slip from stress drop with roughness, optionally capped at the
  local deficit; dip correction; optional quenched structure (strength
  field, segment-boundary barriers, hierarchical patches); analysis of
  seismicity statistics, recurrence, synthetic paleoseismic records,
  record matching and constraint-weighted predictions; CLI for runs,
  analysis and restartable parallel sweeps.
* **Single fault, complete.** `examples/single_fault/`: overview of three
  controls (branching ratio, nucleation clock, size rule; 7 cases x 10
  realizations x 10,000 yr) with figures, diagnostics and example
  ruptures; the development experiments E1 to E22 in
  `docs/experiments.md`.
* **Alpine Fault, complete for the intercomparison.**
  `examples/nz_alpine/`: three geometries, calibration by approximate
  Bayesian computation (rounds 5 to 7, 900 runs), realization sets A to D
  (clock, size rule, strength heterogeneity, segment-boundary barriers;
  50 realizations of 20,000 yr per case), record matching against the
  three paleoseismic records, and constraint-weighted predictions with
  ranges (`results/sets/predictions.md`). What remains for the exercise
  is exporting catalogs once the organisers specify a format.
* **Implemented, not yet explored.** The hierarchical patch structure
  (`hier_*`, `stochseis/hierarchy.py`, E23 in `docs/experiments.md`):
  first single-fault runs only.
* **Next.** UCERF3 (a network of connected fault sections), then fault
  systems represented only probabilistically (section 6).

## 2. What has been established

The numbered findings in `README.md` ("What we have learned") are the
summary. In one paragraph: timing regularity, the shape of the size
distribution and aftershock statistics are controlled separately by the
nucleation clock, the size rule and the branching ratio; regular
recurrence at a point requires slip limited by the accumulated deficit;
Gutenberg-Richter has to be imposed; rupture size is decided by arrest;
mean recurrence intervals are easy to match and regularity is the
discriminating observation, judged on records of the observed length;
regular records come with characteristic sizes and weak triggering; on
the Alpine Fault the regular, characteristic, weakly triggered models
match all three records about one time in six, and no control tested
reproduces the contrast between the central and south records.

## 3. Decisions and their reasons

Model rules:

1. **Slip** is the static crack solution for uniform stress drop on the
   realised rupture outline, times mean-one roughness. Normalising a
   tapered random field to a moment-weighted stress drop was rejected
   (edge singularities dominate; mean slip came out at half the crack
   value). For the Alpine Fault, slip is capped patch by patch at the
   local deficit (`local_cap`): with uncapped crack slip, site recurrence
   CV cannot fall below about 0.4; slip equal to the deficit was rejected
   because magnitude then decouples from rupture size (b about 0).
2. **Arrest** is local admission with contiguity (flood growth); a
   region where the admission probability falls below the percolation
   threshold cannot be crossed. Rejected: a front-averaged hazard
   (ruptures crossed barriers) and admission without contiguity (sparse,
   hole-ridden ruptures). Enclosed holes are filled after growth.
3. **State** for arrest and nucleation is the slip deficit. The backslip
   stress is tracked for diagnostics only: rough slip leaves tens of MPa
   of patch-scale stress heterogeneity.
4. **Rate calibration** is analytic (`docs/design.md` section 2.5). The
   earlier Monte Carlo pilot mean was dominated by rare large events,
   several times too small, and affected every result before 1 November
   2026 (those results and figures have been removed).
5. **Nucleation** is deficit-weighted, with productivity-preserving
   placement of triggered events; this closes the budget without a
   global feedback or a slip cap. The exponent of the weighting is the
   nucleation clock.
6. **Sizes** are continuous below the patch size (sub-patch events with
   budget debit) and through the last patch increment, which removed a 28
   percent step in the frequency-magnitude distribution.
7. **Size rule.** Scale-free stopping by default; the hybrid rule
   (state-only arrest above an area threshold) is used for the Alpine
   Fault because regular records need it.
8. **Branching ratio 0.5** by default (Bath gap about 1.2 on the single
   fault); 0.05 to 0.3 for the Alpine Fault, where regularity was given
   priority over a standard Bath gap (decision by T. Becker).
9. **Hierarchical barriers** act inside the admission ratio (no extra
   parameter); spin-up is handled by run length and the analysis burn-in,
   not by recalibration (fewest parameters, fastest at equal accuracy).

Observation model and comparison (Alpine Fault):

10. Detection threshold 2 m at the south and central sites and 1 m at
    the north site (whose record includes events of about 1 m); events
    closer than 30 yr merged.
11. Models are compared with records of the observed length (the last n
    events before a common end time at all sites), not with long-run
    averages, which overstate irregularity.
12. Predictions weight each realization by the fraction of its records
    that match all three observations; model classes contribute in
    proportion to their fit.

Repository conventions: patches against the GitHub head applied with
`git apply --index`; cautious, neutral wording in documentation; no
reference to individual simulators by name (they are "physics-based
models of the intercomparison").

## 4. Behaviour of the defaults (200 x 15 km, 20 mm/yr)

From 20 seeds of 10,000 yr (E22) and the overview (10 seeds): b 1.00 to
1.01 below width saturation and 1.28 above, maximum M7.6, Bath gap 1.17
to 1.24, Omori p 0.96 to 0.98, CV at points 0.62 to 0.63 (slip of 1 m or
more), whole-fault M >= 7 CV 0.74, foreshock fraction 0.15 with catalogs
complete to M4 and 0.33 with catalogs complete to M3 (chance level about
0.15), budget closure within 1 percent. Runtime per 10,000 yr on one
core: about 1 min with events down to M4, about 12 min down to M3.

## 5. Open issues

1. **South-central contrast on the Alpine Fault.** Not reproduced by the
   clock, size rule, strength heterogeneity or segment-boundary barriers
   (barriers separate the sites but make both irregular). Remaining
   possibility: an intrinsically different southern segment (stress
   drop, clock or loading), which is site-specific tuning.
2. **Regularity versus aftershock statistics** is a property of the model
   class, to be stated in the comparison.
3. **Hierarchical patches (E23)**: sweep not yet run.
4. **Moment-length scaling** above width saturation goes as about L^1.5
   unless slip extends below the seismogenic depth (`deep_slip_km`).
5. **No kinematics**: onset times but only final slip, so no moment-rate
   functions.
6. **Elasticity** is spectral with a dip correction (L0 and L1); there is
   no stress transfer between separate fault sections. Needed for
   networks only if interaction is to be represented physically (L2,
   section 6).
7. Minor: events capped below mag_min (filtered in analysis); the
   free-surface image is exact for antiplane components only; sub-patch
   events carry about 2 percent of the moment.

## 6. Next stages

The design for networks is in `docs/networks.md` (data model, three
levels of elasticity, loading and rake, observables, tuning, stage plan
N1 to N5). N1 (single triangulated surfaces) is done.

### 6.1 UCERF3: a network of connected fault sections

Goal: `examples/ucerf3/` with a README and `run_test.sh`, in the same
form as the two existing examples, showing how ruptures propagate
through a real fault system and what the stochastic model predicts for
multi-fault ruptures and their rates.

Steps:

1. **Multi-section fault container (N2).** Extend `Fault` to hold several
   sections, each a structured grid with its own trace, dip, depth
   range, rake and slip rate, plus connections (patch pairs across
   section ends or branch points with separation, change of strike and
   dip, and an opening probability). Rupture growth already runs on an
   adjacency graph; connection edges are added to it, and a rupture
   crosses an edge with the connection's opening probability. Test first
   on the Alpine Fault split at segment boundaries, which must reproduce
   the single-surface results with connections fully open.
2. **UCERF3 loader (N3).** Read the fault subsections (trace, dip, dip
   direction, upper and lower depth, rake, slip rate, aseismicity;
   typically GeoJSON from OpenSHA) and build sections from traces;
   connection rules as in UCERF3 (jumps up to 5 km, limits on azimuth
   and rake change) as defaults for the opening probabilities.
3. **Calibration targets.** Regional rates from the UCERF3 inversion
   (magnitude-frequency distribution, participation rates per section)
   and the paleoseismic sites used there; the record-matching and
   prediction tools carry over unchanged.
4. **Questions for the example.** How often ruptures jump between
   sections as a function of separation; whether the stochastic model
   reproduces the multi-fault ruptures that the inversion requires; how
   the controls identified here (clock, size rule, clustering, boundary
   barriers) act in a network.
5. **Performance.** Networks have 10^4 to 10^5 patches; the flood loop
   in `rupture.grow_rupture` and the per-rupture shortest paths are the
   hot spots (a numba or C version of the flood loop first).

### 6.2 Fault systems represented only probabilistically (N5)

Goal: a model with no fixed fault map. A network generator draws
sections (traces, dips, depths, slip rates consistent with a regional
loading constraint such as a transect sum) and connections with gap- and
angle-dependent opening probabilities; each simulation realises one
network, and the simulator sees only the fault container, so nothing
downstream changes. Epistemic spread then comes from sampling networks as
well as parameters, and predictions follow from the same constraint
weighting. Natural first test: a generated network whose statistics
resemble a UCERF3 region, compared with the UCERF3 example.

### 6.3 Also pending

* E23 sweep on the single fault (hierarchical patches).
* Intercomparison catalog export from the constrained Alpine ensemble.
* L2 half-space operator (cutde, N4): accurate dip effects and, if
  wanted, stress transfer between sections.
* Housekeeping: a line-by-line human review of the code; check the references in `docs/design.md` against the literature.

## 7. Practical notes

* `python3 -m stochseis.cli run - out.npz --set key=value ...`, then
  `analyze out.npz --plots dir` (summary, deficit, time series and
  regime figures, scalar diagnostics); `sweep spec.json outdir --nproc N`
  for parameter studies (restartable). `Params.from_json` for
  configurations; seeds make runs deterministic.
* Each example has `run_test.sh` (all runs and figures its README needs).
  Their outputs are written into the repository folders; discard them
  (`git checkout` and `git clean` on the result folders) before applying
  a patch that contains the same files.
* Figure scripts take a case as a directory of catalogs, one design
  vector (`dir:v0001`) or one realization (`dir:v0001_s5`).
* Hot spots: `dijkstra` per grown rupture and the flood loop in
  `rupture.grow_rupture`; the full-grid FFT in `FaultState` per flush.
* References in `docs/design.md` were written from memory; verify before
  citing.
