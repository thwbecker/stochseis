# stochseis

A stochastic seismicity generator on fault patch graphs, intended as a
baseline that physics-based earthquake cycle models should be measured
against. Large ruptures carry stress-drop-normalised, rough slip; small
events and clustering follow an ETAS process; a slip budget couples
everything to the long-term loading.

## Authorship

The code, tests and documentation in this repository were written by
Claude (Anthropic's AI model, Claude Opus 5.5) in a chat session with
Thorsten Becker (UT Austin), who specified the scientific goals, made the
modelling decisions and reviewed the results. The design discussion that
led to the current rules, including approaches that were tried and
rejected, is summarised in `HANDOFF.md`. The code has been tested as
described below but has not yet been reviewed line by line by a human;
please treat results accordingly. References in `docs/design.md` were
written from memory and have not all been checked against the
literature.

## Status

This is a first iteration (Stage 0: single planar or segmented faults).
It runs and its tests pass, but the statistics it produces have only
begun to be explored and several modelling choices are provisional. On a
single 200 x 15 km strike-slip fault with default parameters the
frequency-magnitude distribution follows Gutenberg-Richter from M3 to M6
up to the fault-size limit near M7.5 (b = 1 by construction below width
saturation), Bath's gap (about 1.25) and Omori decay emerge, site
recurrence coefficients of variation are about 0.6 for slip of 1 m or
more, and the slip budget closes. These values come from a reduced
parameter study on one fault geometry (ROADMAP step 1). A calibration
error that affected all results before 1 November 2026 has been fixed
(`docs/design.md` section 2.5). Figures in `runs/` come from intermediate
versions of the code and single seeds and should not be compared across
versions.

## Documentation

* `docs/design.md`: core assumptions first, then the rules, parameter
  provenance, known limitations and references.
* `docs/parameters.md`: every parameter with its default and role.
* `docs/analysis.md`: the statistics and figures.
* `docs/networks.md`: design for 3-D geometries, fault networks, UCERF3
  and probabilistic networks; `docs/alpine.md`: the New Zealand Alpine
  Fault intercomparison (inputs, status, plan).
* `docs/experiments.md`: log of the single-fault experiments, with
  settings, results, conclusions and the specs in `examples/experiments/`
  that reproduce them.
* `HANDOFF.md`: decisions and their reasons, rejected alternatives,
  current behaviour, open issues in order of importance.
* `ROADMAP.md`: concrete next steps.

## Requirements

Python 3.10 or later, numpy, scipy, matplotlib; pytest for the tests;
optionally cutde for the half-space check in `examples/checks/`.
No compiled code. A 1500-year catalog on a 200 x 15 km fault at 1 km
patch size (about 160,000 events above M3) takes under a minute on one
core.

## Installation and quick start

```bash
pip install -e .            # or set PYTHONPATH to the repository root
python -m pytest tests      # 26 tests, about 60 s

# one run with defaults, 1500 years, seed 3
python -m stochseis.cli run - cat.npz --set duration_yr=1500 seed=3 --summary cat.json

# metrics and figures
python -m stochseis.cli analyze cat.npz --plots plots/

# parameter sweep on four processes
python -m stochseis.cli sweep examples/sweep_arrest_stressdrop.json sweep_out --nproc 4
```

From Python:

```python
from stochseis import Params, run
from stochseis import analysis, plots

cat = run(Params(duration_yr=1000, seed=1))
m = analysis.summary_metrics(cat)        # dict of scalar diagnostics
plots.summary_figure(cat, "summary.png")
plots.time_series_figure(cat, "time_series.png")
cat.save("cat.npz")
```

## Layout

```
stochseis/
  params.py      all parameters, JSON I/O
  geometry.py    patch-graph fault interface, planar builder, fault factory
  mesh.py        triangulated surfaces (STL) to section grids
  sites.py       synthetic paleoseismic records and comparison with data
  stiffness.py   spectral static stiffness, crack solver, analytic checks
  randomfield.py von Karman fields
  slip.py        slip distribution for a rupture outline
  rupture.py     growth, admission and arrest
  state.py       deficit, backslip stress, maturity
  etas.py        Omori times, spatial kernel, productivity
  simulator.py   event-driven simulation and budget calibration
  catalog.py     catalog container and file format
  analysis.py    statistics
  plots.py       figures
  cli.py         run / analyze / sweep
tests/           pytest suite (analytic crack checks, spectra, growth
                 statistics, ETAS pieces, end-to-end run)
docs/            design.md, parameters.md, analysis.md
examples/        quickstart, sweep specifications, sweep plotting
runs/            figures from development runs (see Status)
```

## Parallelism

Single runs are serial. Sweeps over parameters and seeds run as
independent processes through `multiprocessing.Pool`. The numerical hot
spots are the shortest-path computation and the flood loop per grown
rupture and the FFTs of the stiffness operator; a compiled version of
the flood loop would be the first step if fault networks with many more
patches make it necessary.

## Output

`Catalog.save` writes one `.npz` file with the event table (time,
magnitude, moment, centroid, extent, mean and maximum slip, stress drop,
rupture duration and speed, arrest code, ETAS parent and generation),
the slip distributions of events above `store_slip_min_mag`, the fault
geometry, the released slip per patch, the parameters and run metadata.
`Catalog.load` restores it.

## License

No license has been chosen yet.
