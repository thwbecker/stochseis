# Alpine Fault intercomparison: stochastic baseline

Status: 4 November 2026. This document summarizes the calibration of
stochseis to the Alpine Fault paleoseismic constraints for the New
Zealand simulator intercomparison. The chronological development record
(rounds 1 to 7, intermediate analyses) is in the git history of this
file; the run summaries kept in `runs/` are listed in section 8.

## 1. Purpose

Physics-based earthquake simulators are run on prescribed Alpine Fault
geometries with prescribed slip rates and compared with paleoseismic
recurrence and single-event displacement at three sites. stochseis
provides a baseline: a stochastic model with the slip budget, a rupture
size rule, deficit-dependent nucleation and aftershock clustering, but
no fault mechanics, calibrated to the same constraints. The calibrated
ensembles show which observations such a model already reproduces and
which features of the physics-based models go beyond it.

## 2. Inputs

* Geometries (`data/nz_alpine/`, STL in km NZTM; not in the repository):
  planar (755 km, dip 50 degrees, 20 km depth), variable dip (761 km;
  80 degrees in the south, 50 central, 55 to 85 in the north) and
  segmented transition (759 km, the variable-dip geometry with abrupt
  dip changes at the segment ends). Meshes are resampled onto 1 km or
  2 km patch grids (`stochseis/mesh.py`; 98.6 to 99.5 percent of the mesh
  area).
* Segment slip rates and rakes (`data/nz_alpine/segments.csv`), blended
  over 10 km at segment boundaries: 27 mm/yr in the south and centre,
  15, 5.5 and 4.9 mm/yr to the north. The north site lies at the Springs
  Junction boundary (9.4 to 9.9 mm/yr).
* Paleoseismic sites (`data/nz_alpine/paleo_sites.csv`): mean recurrence
  interval (RI), its standard deviation (sigma of individual intervals,
  so CV = sigma / RI) and mean single-event displacement (SED):

  | site | RI [yr] | CV | SED [m] (range) |
  |---|---|---|---|
  | south | 268 | 0.40 | 7.5 (7.0 to 8.0) |
  | central | 249 | 0.23 | 7.1 (4.5 to 8.5) |
  | north | 413 | 0.25 | 4.0 (1.0 to 9.4) |

* Synthetic paleoseismic records (`stochseis/sites.py`): an event is
  recorded at a site if its slip on the nearest surface patch reaches a
  detection threshold, 2 m at the south and central sites and 1 m at the
  north site (whose observed displacements include events of about 1 m);
  recorded events closer than 30 yr are merged, as a paleoseismic record
  cannot separate them (`site_merge_yr`).

## 3. Model class

The calibrated model combines:

* backslip loading with the segment slip rates; the slip deficit is the
  state;
* deficit-weighted nucleation (intensity proportional to deficit to the
  power `nucleation_exponent`) and ETAS aftershocks (`branching_ratio`);
* rupture growth by flooding with deficit-dependent admission
  (`arrest_scale`, `admission_exponent`) and a hybrid size rule: a
  scale-free size hazard (Gutenberg-Richter, `b_value`) for ruptures
  smaller than `hazard_max_area_km2`; larger ruptures end only where
  the fault state cannot support them;
* crack slip set by stress drop, capped patch by patch at the local
  deficit (`slip_mode = "local_cap"`), with roughness (`slip_cv`);
* the L1 dip correction of the crack amplitude (`dip_correction`,
  `stochseis/dipcorr.py`): half-space / spectral mean-slip ratio as a
  function of dip, rake and burial, 1.0 to 1.1 for the Alpine Fault.

How this class was arrived at (rounds 1 to 4, coarse planar geometry):

| slip rule and sizes | result |
|---|---|
| stress-drop slip, scale-free sizes | means matched, but site CV not below about 0.4 (the slip released at a site varies from event to event; CV of intervals equals CV of displacement) and north RI 460 to 1100 yr |
| slip equal to the local deficit (slip-predictable) | site means matched, but magnitude decouples from rupture size: b about 0, Bath gap 0.7; rejected |
| deficit-capped slip, hybrid sizes | all three sites matched with standard seismicity; adopted |

The decomposition of site variability behind the first row is in
`examples/alpine/decompose_site_cv.py`.

## 4. Calibration

Approximate Bayesian computation (`examples/alpine/`): 100 parameter
vectors drawn from the priors, 3 seeds x 10,000 yr each; each vector is
scored by the mean squared standardized misfit of RI, CV and SED at the
three sites (errors: standard error of the observed mean RI, 0.15 for
CV, half the observed SED range) and, with `--seismicity`, the b-value of
grown ruptures (1.0 +- 0.1) and the Bath gap; accepted if the misfit is
below 4 (`abc_scoring.py`). The best 20 vectors are rerun as posterior
ensembles (5 seeds x 20,000 yr, catalogs kept) and rescored with those
runs. `run_full_campaign.sh` runs the whole chain for a set of
configurations (`campaign_configs.json`) and writes the summary tables
and figures (`campaign_summary.py`, `overview_figures.py`,
`site_analysis.py`, `spectrum_figures.py`).

Configurations: planar, variable dip and segmented transition on 1 km
patches ("fine"); planar and variable dip on 2 km patches ("coarse").

Decision on priorities (T. Becker): regularity at the paleoseismic sites
matters more than a Bath gap of 1.2; the Bath gap is reported as an
outcome and not imposed (`--bath-sigma 1`).

## 5. Results

### 5.1 Reference calibration (round 5, `abc5_planar_l1.json`)

Priors: stress drop 5 to 10 MPa, nucleation exponent 2 to 12, branching
ratio 0.05 to 0.3, area threshold 300 to 3000 km^2, arrest scale 0.2 to 1,
admission exponent 2 to 8, b 0.8 to 1.1, initiation rate 0.1 to 1 of the
calibrated rate, slip CV 0.1 to 0.5. All five configurations
(`runs/alpine_full_campaign/`).

* 64 to 69 of 100 vectors accepted in each configuration, 56 to 66 of
  them shared between any two configurations; all posterior vectors stay
  accepted after rescoring, and the best rescored vector was ranked first
  or second by the ABC (misfits 0.49 to 0.69).
* The accepted parameter distributions are close to the priors: within
  this model class the site constraints select the class, not a
  parameter vector.
* Site envelope (fine configurations, events merged at 30 yr, 2 m at all
  sites; median of the configuration medians and union of the 5 to 95
  percent ranges):

  | site | RI [yr] | CV | SED [m] | observed |
  |---|---|---|---|---|
  | south | 261 [221-302] | 0.37 [0.25-0.49] | 6.8 [5.9-8.0] | 268 / 0.40 / 7.5 |
  | central | 269 [231-320] | 0.35 [0.26-0.47] | 6.9 [6.0-8.4] | 249 / 0.23 / 7.1 |
  | north | 546 [424-691] | 0.37 [0.26-0.51] | 4.8 [3.9-6.1] | 413 / 0.25 / 4.0 |

  With the 1 m threshold at the north site its RI is about 440 to 490 yr.
* Event statistics (all geometries alike): Gutenberg-Richter with b of
  about 1 up to M7, an excess between M7.3 and M8.0, maximum M8.3 to 8.4;
  Omori decay with p of about 1; Bath gap 1.43 to 1.47 (weak
  triggering); moment released per 100 yr with a CV of 1.1 to 1.2;
  whole-fault M >= 7.5 events every 171 to 177 yr (CV 0.63 to 0.65).
* Two rupture populations: shorter than about 60 km (M <= 7.4) and
  through-going 400 to 700 km (M7.9 to 8.4); events recorded at the sites
  are almost all M8.0 to 8.4. Co-rupture at 2 m: central with south 0.70
  to 0.76, central with north 0.35 to 0.42.
* Hypocentres concentrate at the base of the fault (deficit accumulates
  where crack slip tapers to zero).
* Site records are slip-predictable: displacement correlates with the
  preceding interval (0.55 to 0.70), much less with the following one
  (0.13 to 0.30).
* Geometry and resolution: site statistics, event statistics and
  co-rupture differ between geometries and between 1 and 2 km patches by
  no more than their sampling errors.

### 5.2 Regular sites (round 7, `abc7_regularity.json`)

Priors as round 5 but area threshold 50 to 500 km^2, nucleation exponent
3 to 30, admission exponent 4 to 12; events merged at 30 yr at run time;
coarse configurations (`runs/alpine_r7/`).

* 58 (planar) and 46 (variable dip) of 100 vectors accepted; accepted
  area thresholds 52 to 489 km^2 (median 180 to 200), stress drops higher
  (median 8.1 to 8.5 MPa) and initiation rates lower (0.26 to 0.28) than
  in round 5.
* Best rescored vector (variable dip): south 262 yr / 0.24 / 7.0 m,
  central 264 yr / 0.23 / 7.0 m, north 417 yr / 0.32 / 4.1 m.
* Ensembles (2 m, merged): central and south RI 242 to 251 yr, CV 0.26
  to 0.30 [0.14-0.37], SED 6.4 to 6.6 m; north CV 0.32 to 0.35.
* The central regularity is reached; the south site becomes too regular
  (observed 0.40); north CV 0.25 is approached but not reached.
* Cost: a characteristic, partly bimodal size distribution
  (Gutenberg-Richter to about M6.8, then a plateau to M8.1; in many runs
  segment-scale M7.6 to 8.0 ruptures, in some almost no M6.5 to 7
  events), Bath gap about 1.7 over all runs, whole-fault M >= 7.5 events
  every 103 to 105 yr (CV 0.80 to 0.84), and lower co-rupture (central
  with south 0.34 to 0.38, south with north 0.02 to 0.04). Records are
  strongly slip-predictable (0.74 to 0.86).

### 5.3 Bath gap (round 6, `abc6_bath.json`)

Branching ratio 0.2 to 0.7 and b 0.9 to 1.2 bring the Bath gap to about
1.0 (M >= 6 mainshocks) and 1.24 to 1.38 (M >= 7, best vectors), but make
sites and the whole-fault sequence less regular (site CV medians 0.41 to
0.45; best central 0.42 to 0.56) and reduce acceptance to 41 to 54 of
100 (`runs/alpine_r6/`). Pairs of large events a few years apart, which
a paleoseismic record counts as one, account for part of this (merging
at 30 yr lowers the CVs by 0.04 to 0.05). Through-going M8 ruptures keep
gaps of about 2 at any branching ratio: on an isolated fault no loaded
area remains for a large aftershock. Regular sites and a Bath gap of 1.2
are in tension in this model class.

### 5.4 One-parameter spectrum: the area threshold

The area threshold of the hybrid size rule is the strongest single
control of site regularity. Holding all other parameters at the best
round-7 variable-dip vector (2 km patches; single 10,000-yr seeds;
2 m, north 1 m, merged at 30 yr):

| case | area threshold | south RI / CV / SED | central | north | Bath (M >= 6) | events >= M5 |
|---|---|---|---|---|---|---|
| characteristic | 100 km^2 | 208 / 0.28 / 5.6 | 219 / 0.24 / 5.9 | 377 / 0.30 / 3.9 | 2.10 | 2,800 |
| reference | 488 km^2 | 257 / 0.24 / 6.9 | 257 / 0.23 / 6.9 | 417 / 0.35 / 4.2 | 1.46 | 7,800 |
| scale-free | 2000 km^2 | 310 / 0.40 / 7.3 | 246 / 0.50 / 5.7 | 465 / 0.44 / 4.6 | 1.15 | 26,900 |

Moving the threshold from 100 to 2000 km^2 takes the fault from
characteristic behaviour (few moderate events, regular sites, large
Bath gap) to nearly scale-free behaviour (many moderate and partial
ruptures, irregular sites, standard Bath gap). Only the middle case meets
all mean constraints with the other parameters fixed.
`spectrum_area.json` reruns the three cases with 5 seeds x 20,000 yr and
catalogs, for `spectrum_figures.py`.

### 5.5 Scenarios compatible with the constraints

`spectrum_figures.py` applied to the posterior ensembles (all accepted
vectors pooled per case) shows the kinds of earthquake sequences that
satisfy the constraints: space-time windows of large ruptures, site
interval distributions against the paleoseismic ones, the
frequency-magnitude distribution, and the rate of M >= 7.5 events by the
combination of sites they rupture (south only, south and central, all
three, and so on). Commands in section 7; results to be added.

## 6. Open points for the comparison with physics-based models

* The observed contrast in regularity between the central (CV 0.23) and
  south (0.40) sites: in every calibrated ensemble the two sites share
  most ruptures and have similar CVs. Reproducing the contrast needs
  partly independent rupture histories at the two sites.
* Regular sites require, in this model class, a characteristic size
  distribution and weak triggering (large Bath gaps). A physics-based
  model that matches the regularity with near-Gutenberg-Richter
  statistics and standard aftershock sequences does something the
  stochastic baseline cannot.
* North site: intervals are matched only with a 1 m detection threshold;
  its regularity (CV 0.25) is not reached by any ensemble.
* Testable predictions of the calibrated baseline: slip-predictable site
  records; hypocentres near the base of the seismogenic zone; depletion
  of M6.5 to 7 events in the regular-site calibrations.

## 7. Reproducing

Mesh files in `data/nz_alpine/` (`alpine_planar_1km_26Jan.stl`,
`variable_dip_1km_26Jan.stl`, `Alpine_full_segmented_transition_05km_d20.stl`).

    # reference calibration, five configurations (about 75 core-hours)
    examples/alpine/run_full_campaign.sh $NR_CPUS alpine_full 100 3 10000 20 5 20000 \
        "planar_fine variable_dip_fine segmented_fine planar_coarse variable_dip_coarse" \
        examples/alpine/abc5_planar_l1.json "--site north=1 --bath-sigma 1"
    # regular sites (coarse, about 10 core-hours; add the fine names for all)
    examples/alpine/run_full_campaign.sh $NR_CPUS alpine_r7 100 3 10000 20 5 20000 \
        "planar_coarse variable_dip_coarse" examples/alpine/abc7_regularity.json "--site north=1 --bath-sigma 1"
    # one-parameter spectrum (15 runs of 20,000 yr)
    python3 -m stochseis.cli sweep examples/alpine/spectrum_area.json alpine_spectrum --nproc $NR_CPUS
    python3 examples/alpine/spectrum_figures.py alpine_spectrum/spectrum \
        characteristic alpine_spectrum:v0000 reference alpine_spectrum:v0001 scale-free alpine_spectrum:v0002 \
        --site north=1
    # scenarios compatible with the constraints
    python3 examples/alpine/spectrum_figures.py alpine_scenarios \
        reference alpine_full/posterior_variable_dip_fine regular alpine_r7/posterior_variable_dip_coarse \
        --site north=1

The round-5 campaign in `runs/alpine_full_campaign/` was run before
merging was available at run time; its merged tables come from
`site_analysis.py --merge 30` and its summaries were rescored with
`--site north=1`.

## 8. Kept run summaries

* `runs/alpine_full_campaign/`: round 5, five configurations (summaries
  at 2 m and with north at 1 m, site tables, merged tables, figures).
* `runs/alpine_r6/`: round 6 (Bath gap), coarse.
* `runs/alpine_r7/`: round 7 (regular sites), coarse, merged tables and
  figures.
