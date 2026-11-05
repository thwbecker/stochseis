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
* `campaign_configs.json`: configurations (geometry, patch size).
* `realization_sets.json`: three families of cases varying one
  irregularity control each (nucleation clock, size rule, barriers).

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
  frequency-magnitude, rupture patterns), e.g. for the scenarios
  compatible with the constraints.
* `record_match.py`: samples synthetic records of the observed length (the
  last n events before an end time, all sites ending together) and gives
  how often they match the observed RI and CV, per site and jointly.
* `predictions.py`: constraint-weighted predictions: pools cases and
  weights each realization by the fraction of its records that match the
  observations at all sites; weighted medians and 5-95 % ranges of the
  metrics to compare with physics-based simulators.
* `control_summary.py`: one figure of what each irregularity control does
  (metrics against the cases of each realization set).
* `parameter_effects.py`: parameter effects pooled over ABC rounds
  (standardized-effect heatmap, binned dependence, acceptance rates).
* `run_realization_sets.sh`: runs the realization sets (50 realizations
  per case) and draws the three pages for one realization and for the
  ensemble; `make_realization_specs.py` writes their specs, `nz_page.py`
  draws page 2 (record-length windows, rupture patterns and sizes).
* `decompose_site_cv.py`: decomposition of site recurrence variability
  for one catalog.

## Typical commands

    examples/alpine/run_full_campaign.sh $NR_CPUS alpine_r7 100 3 10000 20 5 20000 \
        "planar_coarse variable_dip_coarse" examples/alpine/abc7_regularity.json
    examples/alpine/run_realization_sets.sh $NR_CPUS alpine_sets 50 20000
    python3 examples/alpine/control_summary.py alpine_sets/controls alpine_sets --site north=1
    python3 examples/alpine/record_match.py alpine_sets/match_B \
        characteristic alpine_sets/B_size:v0000 reference alpine_sets/B_size:v0001 scale_free alpine_sets/B_size:v0002 \
        --site north=1
    python3 examples/alpine/parameter_effects.py alpine_effects alpine_full/abc_*/sweep_summary.json \
        alpine_r6/abc_*/sweep_summary.json alpine_r7/abc_*/sweep_summary.json --site north=1

Quick end-to-end checks (a few minutes on one core each):

    examples/alpine/run_realization_sets.sh 1 alpine_sets_smoke 2 3000 A_clock

    examples/alpine/run_full_campaign.sh 1 alpine_smoke 2 1 800 1 1 1000
