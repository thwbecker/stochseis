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
