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

Results from the single-fault experiments (`examples/single_fault/`) and
the Alpine Fault application (`examples/nz_alpine/`):

* **Two main axes of complexity.** Timing regularity is controlled by
  the nucleation clock, which leaves the size distribution and the
  aftershock statistics unchanged. The shape of the size distribution
  (characteristic to scale-free) is controlled by the size rule, which
  also changes aftershock productivity, the Bath gap, rupture patterns
  and, through them, regularity. Clustering (the branching ratio) is a
  third control that couples aftershock statistics to regularity: strong
  triggering gives standard Bath gaps but irregular recurrence.
* **Prerequisites.** Whether regular recurrence at a point is possible at
  all depends on how slip is assigned: when slip is set by stress drop
  alone, the slip released at a point varies from event to event and
  point recurrence stays irregular (CV of 0.4 and more); capping slip at
  the accumulated deficit allows regular recurrence with normal
  seismicity statistics. Gutenberg-Richter does not emerge from the
  fault state alone and has to be imposed.
* **Final size is decided by arrest.** Nothing at nucleation determines
  how large a rupture becomes; it grows until it reaches fault that is
  not ready to fail. This matches observations of large earthquakes that
  begin like small ones.
* **Mean recurrence intervals are easy, regularity is hard.** With the
  slip budget calibrated, mean recurrence intervals and displacements at
  paleoseismic sites follow from slip rate and stress drop; the
  coefficient of variation of the intervals is the discriminating
  observation, and it must be judged on synthetic records of the
  observed length (a dozen events), which are more regular than long
  runs.
* **Alpine Fault.** Models with regular, partly characteristic rupture,
  weak triggering and deficit-capped slip match the three paleoseismic
  records about one time in six; irregular models (soft clock,
  scale-free sizes) about one time in a hundred. In most models the
  south and central sites share nearly all large ruptures; barriers at
  the segment boundaries separate them but make both irregular, so the
  observed contrast between a regular central and a less regular south
  record is not reproduced by any control tested. The
  constraint-weighted ensemble gives predictions with ranges for the
  model comparison.

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
