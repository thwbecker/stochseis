# Alpine Fault runs

1. Copy the reference STL files into `data/nz_alpine/` (they are not in
   the repository; see `data/nz_alpine/README.md`).
2. First tuning sweep (24 runs of 5000 yr per geometry, about 6
   core-hours each):

       python3 -m stochseis.cli sweep examples/alpine/tune1_planar.json tune1_planar --nproc N --no-catalogs
       python3 -m stochseis.cli sweep examples/alpine/tune1_variable_dip.json tune1_vardip --nproc N --no-catalogs

   The summaries contain, per site and detection threshold (1, 2, 4 m),
   the number of recorded events, mean recurrence interval, its z-score
   against the observed mean, CV and mean single-event displacement
   (keys `site_<name>_<stat>_th<threshold>`), plus co-rupture fractions.
3. Seed-averaged table, for example:

       python3 examples/aggregate_sweep.py tune1_planar/sweep_summary.json p_stress_drop_mpa p_nucleation_exponent p_branching_ratio

   (the aggregation script lists the generic metrics; the site keys are in
   the JSON and CSV summaries).

## Approximate Bayesian computation (with the segment slip rates)

    python3 -m stochseis.cli sweep examples/alpine/abc1_planar.json abc1_planar --nproc N --no-catalogs
    python3 -m stochseis.cli sweep examples/alpine/abc1_variable_dip.json abc1_vardip --nproc N --no-catalogs
    python3 examples/alpine/abc_accept.py abc1_planar/sweep_summary.json 2

Sweeps can be restarted after an interruption: each run writes its
summary to `<outdir>/runs/` when it finishes, and a rerun of the same
command skips runs that are already there. A run whose simulation or
summary fails records the error in its summary (key `error`) instead of
stopping the sweep.

Each spec is 400 runs (200 parameter vectors x 2 seeds) of 5000 yr at
2 km patches, about 1 minute each on one core, i.e. under an hour on 8
cores. The acceptance script scores each vector against RI, CV and SED at
the three sites for a chosen detection threshold (1, 2 or 4 m) and prints
the best vectors and the accepted ranges. The tune1 specs predate the
segment slip rates and are superseded.

## Round 2 (planar)

    python3 -m stochseis.cli sweep examples/alpine/abc2_planar.json abc2_planar --nproc N --no-catalogs

300 runs of 10,000 yr, about 2 to 4 minutes each on one core (10 to 20
core-hours). Round 1 results and the reasons for the narrowed priors are
in docs/alpine.md.

## Round 3 (planar, slip-predictable mode)

    python3 -m stochseis.cli sweep examples/alpine/abc3_planar_slip_predictable.json abc3_planar --nproc N --no-catalogs
    python3 examples/alpine/abc_accept.py abc3_planar/sweep_summary.json 2

300 runs of 10,000 yr; this mode is faster (about 30 to 60 s per run on
one core, roughly 3 to 5 core-hours in total).

## Round 4 (planar, hybrid sizes with deficit-capped slip)

    python3 -m stochseis.cli sweep examples/alpine/abc4_planar_hybrid.json abc4_planar --nproc N --no-catalogs
    python3 examples/alpine/abc_accept.py abc4_planar/sweep_summary.json 2

300 runs of 10,000 yr, about 20 to 40 s each. For any single run kept
with its catalog, `examples/alpine/decompose_site_cv.py catalog.npz 2`
decomposes the site variability (docs/alpine.md).
