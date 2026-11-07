# stochseis

A stochastic earthquake-sequence generator for faults and fault systems,
intended as a baseline against which physics-based earthquake cycle
simulators can be measured. It produces long synthetic catalogs with
slip distributions from a few rules and a slip budget, without friction,
stress transfer or dynamics. When such a model, tuned to the same
constraints, reproduces an observation, that observation does not by
itself require the physics of a simulator; where it fails, the
simulator has something to show.

## Authorship

The code, tests and documentation in this repository were written by
Claude (Anthropic's AI model, Claude Opus 5.5) in a chat session with
Thorsten Becker (UT Austin), who specified the scientific goals, made the
modelling decisions and reviewed the results. The code has been tested
as described below but has not been reviewed line by line by a human;
please treat results accordingly. References in `docs/design.md` were
written from memory and have not all been checked against the
literature.

## Philosophy

* **A null hypothesis for simulators.** Physics-based earthquake cycle
  simulators produce long catalogs from friction, stress transfer and
  dynamics. stochseis asks how much of what they produce follows from
  simpler ingredients: a slip budget, a rule for where ruptures start, a
  rule for where they stop, and clustering. Features that the stochastic
  model reproduces once it is tuned to the same constraints do not by
  themselves demonstrate the physics; features it cannot reproduce are
  where simulators show something.
* **Few rules, few parameters.** Each rule is kept only if it is needed
  to match an observation, and it is described by as few parameters as
  possible. Where two formulations give the same result, the one with
  fewer parameters and the faster run is kept.
* **Same constraints, same observation model.** Models are compared with
  data the way the data were collected: synthetic paleoseismic records
  use a detection threshold, merge events too close to separate, and
  have the length of the real records.
* **Ensembles, not single runs.** Every statement is checked over many
  realizations; single catalogs are used only as illustrations.
* **A range, not a best fit.** The constraints select a class of models
  rather than one parameter vector; predictions are therefore given as
  the range of behaviour of all models consistent with the constraints,
  weighted by how well they match.

## What the code does

A fault (planar, segmented, or a triangulated surface such as the
Alpine Fault meshes) is divided into patches of 1 to 2 km. On it:

* **Loading and state.** Every patch accumulates slip deficit at its
  long-term slip rate; earthquakes release it. The deficit is the state
  of the fault, and the budget closes: over long runs released moment
  equals loaded moment.
* **Nucleation.** Background earthquakes start at a rate set by the
  budget, preferentially where the deficit is high; the exponent of that
  dependence (`nucleation_exponent`) acts as a renewal clock, from
  Poisson-like to strongly time-predictable.
* **Clustering.** Each event triggers aftershocks (ETAS, Omori decay,
  productivity scaling with magnitude, `branching_ratio`).
* **Rupture growth and arrest.** A rupture grows outward from its
  hypocentre over neighbouring patches. At each step it may stop: by a
  scale-free size hazard (Gutenberg-Richter, `b_value`), and wherever the
  fault is not ready, i.e. where the deficit is small compared with the
  slip the growing rupture would need (`arrest_scale`). With the hybrid
  size rule (`hazard_max_area_km2`), ruptures larger than a threshold
  area stop only where the fault is not ready. Optional quenched
  structure: a strength field, barriers at segment boundaries, and a
  persistent hierarchy of nested patches.
* **Slip.** Each rupture receives a crack-shaped slip distribution set by
  its stress drop, with self-affine roughness, optionally capped patch by
  patch at the local deficit (`slip_mode = "local_cap"`); a dip
  correction accounts for dipping faults.
* **Outputs.** Catalogs with event times, magnitudes, geometry and slip
  distributions, and analysis tools for frequency-magnitude statistics,
  aftershock statistics (Omori, Bath), recurrence at points and on the
  whole fault, synthetic paleoseismic records, and comparisons with
  observations.

`docs/design.md` gives the rules and their justification in detail.

## What we have learned

Single fault (`examples/single_fault/`, 10 realizations of 10,000 yr per
case) and Alpine Fault (`examples/nz_alpine/`, 50 realizations of
20,000 yr per case and 900 calibration runs):

1. **Three controls with separate signatures.**
   * The *nucleation clock* (how strongly nucleation favours accumulated
     slip deficit) sets timing regularity and nothing else: on the single
     fault, CV of recurrence at points from 0.63 to 0.51 with b, Bath gap
     and the size distribution unchanged.
   * The *size rule* (whether large ruptures stop at random or only where
     the fault is not ready) sets the shape of the size distribution,
     from Gutenberg-Richter to a characteristic population of large events
     (3.5 times the Gutenberg-Richter rate at M >= 7), and with it the
     Bath gap and the regularity of large events.
   * *Clustering* (the ETAS branching ratio) sets aftershock and foreshock
     statistics (Bath gap 1.43 to 1.01, foreshock fraction 0.08 to 0.24)
     and, when strong, makes recurrence irregular.
2. **Regular recurrence needs slip limited by the deficit.** With slip set
   by stress drop alone, slip at a point varies from event to event and
   recurrence there stays irregular (CV of 0.4 and more); capping slip at
   the accumulated deficit allows regular recurrence with normal
   seismicity statistics.
3. **Gutenberg-Richter has to be imposed**; ruptures stopped only by the
   fault state give characteristic size distributions, with or without
   quenched heterogeneity.
4. **Final size is decided by arrest.** Nothing at nucleation determines
   how large a rupture becomes, consistent with large earthquakes that
   begin like small ones.
5. **Mean intervals are easy, regularity is hard.** Once the slip budget
   is calibrated, mean recurrence intervals and displacements follow from
   slip rate and stress drop; the coefficient of variation is the
   discriminating observation, and it has to be judged on records of the
   observed length (a dozen events), which are more regular than long
   runs.
6. **Regularity and standard aftershock statistics are in tension.**
   Regular records come with characteristic sizes and weak triggering
   (Bath gaps of 1.3 and more); strong triggering gives standard Bath gaps
   but irregular records.
7. **Alpine Fault.** Models with regular, partly characteristic rupture,
   weak triggering and deficit-capped slip match all three paleoseismic
   records together about one time in six; irregular models about one
   time in a hundred. Geometry (planar, variable dip, segmented) and
   resolution (1 or 2 km) matter little. No control tested (clock, size
   rule, strength heterogeneity, segment-boundary barriers) reproduces
   the contrast between the regular central and the less regular south
   record: the sites either share nearly all large ruptures, or, when
   barriers separate them, both become irregular. The constraint-weighted
   ensemble gives predictions with ranges for the comparison with the
   physics-based simulators of the intercomparison.

## Examples

* `examples/single_fault/`: a single planar strike-slip fault; how a few
  parameters change Gutenberg-Richter statistics, Omori decay, Bath's
  law and recurrence, with spreads from several realizations and example
  slip distributions; the development experiments E1 to E23.
* `examples/nz_alpine/`: the Alpine Fault, New Zealand (three mesh
  geometries, segment slip rates, three paleoseismic sites); how the
  parameters, including barriers at segment boundaries, change the
  behaviour of a long fault, and the model tuned to the paleoseismic
  constraints of the simulator comparison, with predictions.
* Planned: UCERF3 fault sections (a network of connected faults), and
  later a fault system represented only probabilistically
  (`docs/networks.md`).

Each example folder has a `README.md` with the results and figures and a
`run_test.sh` that reproduces them.

## Requirements and installation

Python 3.10 or later, numpy (1.22 or later), scipy (1.10 or later),
matplotlib; pytest for the tests; optionally cutde for the half-space
checks in `examples/checks/`. No compiled code.

    python3 -m pip install --user -e ".[test]"     # or set PYTHONPATH to the repository root
    python3 -m pytest tests                          # 47 tests, about a minute

    # one run with defaults, 1500 years, seed 3, and its figures
    python3 -m stochseis.cli run - cat.npz --set duration_yr=1500 seed=3 --summary cat.json
    python3 -m stochseis.cli analyze cat.npz --plots plots/

    # a sweep (parameters and seeds from a JSON spec) on four processes
    python3 -m stochseis.cli sweep examples/single_fault/specs/parameter_overview.json sweep_out --nproc 4

From Python:

```python
from stochseis import Params, run
from stochseis import analysis, plots

cat = run(Params(duration_yr=1000, seed=1))
m = analysis.summary_metrics(cat)        # dict of scalar diagnostics
plots.summary_figure(cat, "summary.png")
cat.save("cat.npz")
```

A 1500-year catalog on a 200 x 15 km fault at 1 km patches (about
160,000 events above M3) takes under a minute on one core. Sweeps run
their members as independent processes; they can be interrupted and
restarted.

## Layout

```
stochseis/          the package
  params.py         all parameters, JSON I/O
  geometry.py       fault patch graph, planar and segmented builders, segment slip rates
  mesh.py           triangulated surfaces (STL) to patch grids
  stiffness.py      spectral static stiffness, crack solutions
  dipcorr.py        dip and rake correction of crack amplitudes
  randomfield.py    von Karman random fields
  hierarchy.py      persistent hierarchical patch structure
  slip.py           slip distribution of a rupture
  rupture.py        growth, admission and arrest
  state.py          slip deficit and stress
  etas.py           aftershock times, locations, productivity
  simulator.py      event-driven simulation and budget calibration
  catalog.py        catalog container and file format
  sites.py          synthetic paleoseismic records
  analysis.py       statistics
  plots.py          figures
  cli.py            run / analyze / sweep
tests/              pytest suite
docs/               design, parameters, analysis, development log, network design
examples/
  quickstart.py
  checks/           half-space comparison and dip-correction table
  single_fault/     README, run_test.sh, overview spec and figures, experiments E1-E23
  nz_alpine/        README, run_test.sh, scripts, specs, results
data/nz_alpine/     Alpine Fault inputs (sites, segments; meshes not included)
```

## Documentation

* `docs/design.md`: core assumptions, the rules, parameter provenance,
  known limitations and references.
* `docs/parameters.md`: every parameter with its default and role.
* `docs/analysis.md`: the statistics and figures.
* `docs/experiments.md`: the development log of the single-fault
  experiments (E1 to E23), with settings, results and the specs that
  reproduce them.
* `docs/networks.md`: design for 3-D geometries, fault networks, UCERF3
  and probabilistic networks.
* `HANDOFF.md`: decisions and their reasons, rejected alternatives, open
  issues; `ROADMAP.md`: next steps.

## Output

`Catalog.save` writes one `.npz` file with the event table (time,
magnitude, moment, centroid, extent, mean and maximum slip, stress drop,
rupture duration and speed, arrest code, ETAS parent and generation),
the slip distributions of events above `store_slip_min_mag`, the fault
geometry, the released slip per patch, the parameters and run metadata.
`Catalog.load` restores it.

## License

No license has been chosen yet.
