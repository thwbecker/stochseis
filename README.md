# stochseis

A stochastic seismicity generator on fault patch graphs, intended as a
baseline for physics-based earthquake cycle models. Large ruptures carry
stress-drop-normalised, rough slip; small events and clustering follow an
ETAS process; a slip budget couples everything to the long-term loading.
The model is described in `docs/design.md`, which lists the core
assumptions first; parameter provenance is in `docs/parameters.md`, the
analysis tools in `docs/analysis.md`.

This is a first iteration (Stage 0: planar faults). It runs and passes
its tests, but the statistics it produces have only begun to be
explored, and several modelling choices are provisional (see the
limitations section of `docs/design.md`).

## Requirements

Python 3.10 or later, numpy, scipy, matplotlib; pytest for the tests.
No compiled code. A 3000-year catalog on a 200 x 15 km fault at 1 km
patch size (about half a million events above M3) takes one to two
minutes on one core.

## Quick start

```bash
# one run with defaults, 2000 years, seed 3
python -m stochseis.cli run - cat.npz --set duration_yr=2000 seed=3 --summary cat.json

# metrics and figures
python -m stochseis.cli analyze cat.npz --plots plots/

# parameter sweep on four processes (see examples/)
python -m stochseis.cli sweep examples/sweep_arrest_stressdrop.json sweep_out --nproc 4

# tests
python -m pytest tests
```

From Python:

```python
from stochseis import Params, run
from stochseis import analysis, plots

cat = run(Params(duration_yr=1000, seed=1))
m = analysis.summary_metrics(cat)        # dict of scalar diagnostics
plots.summary_figure(cat, "summary.png")
cat.save("cat.npz")
```

## Layout

```
stochseis/
  params.py      all parameters, JSON I/O
  geometry.py    patch-graph fault interface and planar builder
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
examples/        sweep specifications
```

## Parallelism

Single runs are serial. Sweeps over parameters and seeds run as
independent processes through `multiprocessing.Pool`; this is the
intended way to use many cores. Numerical hot spots are the shortest-path
computation per fault event and the FFTs of the stiffness operator; both
are already vectorised and would be the first candidates for a compiled
extension if larger fault systems make it necessary.

## State of the work

`HANDOFF.md` records the decisions, what was tried and rejected, the
behaviour of the current defaults and the open issues, in order of
importance. In short: on a single 200 x 15 km strike-slip fault the
frequency-magnitude distribution follows GR from M3 to M6 within a few
percent with b = 1 by construction, Bath's gap and Omori decay emerge
with plausible values, site recurrence CVs are 0.5 to 0.6, and the slip
budget closes. The large-magnitude end depends on how strongly the
budget controls rupture arrest (the `arrest_scale` parameter), and that
regime choice is the first thing to settle with longer runs. Figures in
`runs/` come from intermediate versions of the code and single seeds.

## Output

`Catalog.save` writes one `.npz` file with the event table (time,
magnitude, moment, centroid, extent, mean and maximum slip, stress drop,
rupture duration and speed, arrest code, ETAS parent and generation),
the slip distributions of events above `store_slip_min_mag`, the fault
geometry, the released slip per patch, the parameters and run metadata.
`Catalog.load` restores it.
