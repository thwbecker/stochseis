# Alpine Fault calibration

Scripts and specs for the stochastic baseline of the New Zealand Alpine
Fault intercomparison. Results and their interpretation are in
`docs/alpine.md`; this file lists what is here and how to run it.

Mesh files go in `data/nz_alpine/` (see `data/nz_alpine/README.md`).

## Specs

* `abc5_planar_l1.json`: reference calibration (round 5), base of the
  campaign; the geometry is set per configuration.
* `abc6_bath.json`: sensitivity case with stronger triggering (Bath gap).
* `abc7_regularity.json`: calibration toward regular site records.
* `spectrum_area.json`: one-parameter spectrum (area threshold of the
  hybrid size rule: 100, 488, 2000 km^2), catalogs kept.
* `campaign_configs.json`: configurations (geometry, patch size).

## Scripts

* `run_full_campaign.sh`: ABC, acceptance, posterior ensembles, summary,
  figures for a list of configurations (arguments in its header;
  restartable).
* `abc_scoring.py`: misfit and acceptance (shared); `abc_accept.py`,
  `make_posterior_spec.py`, `campaign_summary.py` use it. Options:
  `--seismicity`, `--site north=1`, `--bath-sigma x`, `--bath-m7`.
* `make_campaign_specs.py`: one ABC spec per configuration.
* `overview_figures.py`: event-statistics and recurrence figures for
  posterior ensembles.
* `site_analysis.py`: site figures and tables (paleoseismic panels,
  co-rupture, predictability; `--merge`, `--reference`).
* `spectrum_figures.py`: cases side by side (space-time, site intervals,
  frequency-magnitude, rupture patterns) for the spectrum and for the
  scenarios compatible with the constraints.
* `decompose_site_cv.py`: decomposition of site recurrence variability
  for one catalog.

## Typical commands

    examples/alpine/run_full_campaign.sh $NR_CPUS alpine_r7 100 3 10000 20 5 20000 \
        "planar_coarse variable_dip_coarse" examples/alpine/abc7_regularity.json
    python3 -m stochseis.cli sweep examples/alpine/spectrum_area.json alpine_spectrum --nproc $NR_CPUS
    python3 examples/alpine/spectrum_figures.py alpine_spectrum/spectrum \
        characteristic alpine_spectrum:v0000 reference alpine_spectrum:v0001 scale-free alpine_spectrum:v0002 \
        --site north=1

A quick end-to-end check (about 10 minutes on one core):

    examples/alpine/run_full_campaign.sh 1 alpine_smoke 2 1 800 1 1 1000
