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

## 5. What we learned

The parameter effects below come from two kinds of evidence: the
realization sets, which vary one control at a time around the reference
calibration (section 6.5; `control_summary.py`), and regressions over the
design vectors of the ABC rounds (`parameter_effects.py`; standardized
coefficients, metric standard deviations per parameter standard deviation,
quoted here from round 4, 200 vectors).

| ingredient | what it controls | evidence |
|---|---|---|
| slip rule | whether site regularity and standard seismicity are possible at all: stress-drop slip leaves site CV near 0.4 or above (released slip varies from event to event); slip equal to the deficit decouples magnitude from size (b about 0); deficit-capped slip with the hybrid size rule satisfies both | rounds 1 to 4 |
| stress drop and initiation rate | mean RI and SED at all sites (slip budget: RI is about SED / slip rate) | coefficients +0.7 (stress drop) and -0.3 (initiation rate) on RI and SED |
| nucleation exponent (renewal clock) | regularity of sites and of the whole-fault sequence, with the frequency-magnitude distribution, Omori decay and Bath gap unchanged; higher exponents also shorten intervals | set A: whole-fault CV 0.21 to 0.37, central window CV 0.21 to 0.35; coefficient -0.6 on RI |
| area threshold of the size rule | the size distribution (characteristic to scale-free), aftershock productivity, Bath gap, rupture patterns and whole-fault regularity; site regularity is greatest at intermediate thresholds | set B: Bath 2.16 / 1.34 / 1.13, whole-fault CV 0.30 / 0.21 / 0.76; coefficients +0.65 on site CV, +0.8 on whole-fault CV, -0.6 on the M >= 7 excess |
| arrest scale and admission sharpness | how far ruptures run into partly loaded fault: co-rupture between sites and north-site regularity | coefficients -0.74 on central-south co-rupture, +0.72 on north CV |
| branching ratio | Bath gap and Omori p; above about 0.3 also less regular sites and whole-fault sequences | coefficients -0.5 on Bath, +0.5 on Omori p; round 6 |
| b | b, and the Bath gap | coefficients +0.96 on b, -0.6 on Bath |
| quenched barriers (strength heterogeneity) | moderately less regular sites, shorter north intervals, fewer full-length ruptures; they do not give neighbouring sites different regularity | set C |
| slip roughness, stress-drop scatter | little, once slip is capped by the deficit | coefficients below 0.15 |
| geometry, resolution | little for these observables (with the L0 and L1 elastic operators) | five configurations, paired designs |
| observation model | detection threshold (north RI 420 to 460 yr at 1 m, 540 at 2 m), merging of close events (CV -0.03 to -0.05), record length (windows of 12 to 13 intervals have CVs about 0.07 lower than long runs) | sections 6.1, 6.5 |

Overall:

* The site records constrain the model class more than its parameters:
  in each calibration the accepted parameter distributions are close to
  their priors.
* Regular site records come with a characteristic size distribution and
  weak triggering in this model class; standard aftershock statistics
  (Bath gap 1.2) come with less regular sites.
* The contrast between a regular central record and a less regular south
  record is not produced by any control tested: the two sites share most
  large ruptures in every case.
* Judged by records of the observed length, the reference calibration is
  consistent with the central and north records.

## 6. Results in detail

### 6.1 Reference calibration (round 5, `abc5_planar_l1.json`)

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

### 6.2 Regular sites (round 7, `abc7_regularity.json`)

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

### 6.3 Bath gap (round 6, `abc6_bath.json`)

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

### 6.4 Scenarios compatible with the constraints

Page 2 of the figure set (`nz_page.py`, rupture patterns and sizes) and
`spectrum_figures.py` applied to posterior ensembles (all accepted
vectors pooled per case) show the kinds of earthquake sequences that
satisfy the constraints: space-time windows of large ruptures, site
intervals, the frequency-magnitude distribution, and the rate of M >= 7.5
events by the combination of sites they rupture.

### 6.5 Realization sets: three irregularity controls

The area threshold is not the only control of site regularity. Single
10,000-yr tests around the reference vector (best round-7 variable-dip
calibration; 2 km patches; 1 m threshold, merged at 30 yr):

| change from the reference | central CV | south CV | north CV | Bath | other |
|---|---|---|---|---|---|
| none (nucleation exponent 12.3, area threshold 488 km^2) | 0.23 | 0.24 | 0.35 | 1.46 | |
| nucleation exponent 3 | 0.40 | 0.44 | 0.42 | 1.45 | event count about the same |
| area threshold 100 km^2 | 0.24 | 0.28 | 0.30 | 2.10 | few moderate events |
| area threshold 100, strength heterogeneity ln-sigma 1.0, fault-length correlation | 0.28 | 0.30 | 0.32 | 1.84 | M >= 8 events 42 to 13 |
| area threshold 100, strength heterogeneity ln-sigma 1.0, 50 km correlation | 0.35 | 0.46 | 0.28 | 1.90 | central RI 190, south 272 |

* Nucleation sharpness (a renewal clock) controls regularity about as
  strongly as the area threshold, without changing the Bath gap or the
  number of events.
* The area threshold controls the size distribution (characteristic to
  scale-free) and with it regularity.
* Short-wavelength strength heterogeneity (persistent barriers) is the
  only ingredient found so far that gives neighbouring sites different
  regularity and intervals, the contrast between the central and south
  records that no calibrated ensemble reproduced.

`examples/alpine/realization_sets.json` defines three families around the
reference vector: A_clock (nucleation exponent 3, 6, 12.3), B_size (area
threshold 100, 488, 2000 km^2) and C_barriers (strength ln-sigma 0, 0.5,
1.0 at 50 km correlation length). `run_realization_sets.sh` runs each
case as 50 independent realizations of 20,000 yr (with catalogs) and
draws three pages for the first realization ("single") and for all
realizations ("ensemble"): page 1, frequency-magnitude distribution,
Omori decay, Bath gaps, moment release and whole-fault recurrence
(`overview_figures.py`); page 2, site recurrence in record-length
windows, SED versus RI, rupture patterns, magnitudes, rupture lengths
(`nz_page.py`); page 3, the intercomparison panels (`site_analysis.py`).
Record-length windows cut each synthetic site record into consecutive
windows with as many intervals as the paleoseismic record at that site
(12 south, 13 central, 3 north) and give the fraction of windows at least
as regular as the observed record. The extremes of each family shift the
mean intervals and are meant to show what the control does, not as
alternative fits.

Results (50 realizations of 20,000 yr per case; `runs/alpine_sets/`).
Window statistics: median window CV and the fraction of record-length
windows at least as regular as the observed record (for the south site,
at least as irregular, i.e. CV >= 0.40):

| set | case | south: CV, P(CV >= 0.40) | central: CV, P(CV <= 0.23) | north: CV, P(CV <= 0.25) | median RI south / central / north [yr] |
|---|---|---|---|---|---|
| A_clock | soft clock (exponent 3) | 0.35, 25 % | 0.35, 5 % | 0.27, 46 % | 307 / 307 / 403 |
| A_clock | medium (6) | 0.26, 7 % | 0.27, 28 % | 0.26, 46 % | 286 / 286 / 397 |
| A_clock | reference (12.3) | 0.21, 2 % | 0.21, 63 % | 0.28, 41 % | 267 / 266 / 399 |
| B_size | characteristic (100 km^2) | 0.22, 0 % | 0.22, 59 % | 0.25, 51 % | 221 / 223 / 426 |
| B_size | reference (488 km^2) | 0.20, 1 % | 0.19, 73 % | 0.29, 37 % | 266 / 269 / 405 |
| B_size | scale-free (2000 km^2) | 0.45, 66 % | 0.44, 4 % | 0.34, 30 % | 262 / 263 / 466 |
| C_barriers | reference (ln-sigma 0) | 0.20, 0 % | 0.20, 63 % | 0.27, 45 % | 268 / 267 / 399 |
| C_barriers | weak (0.5) | 0.24, 3 % | 0.23, 48 % | 0.22, 56 % | 259 / 260 / 342 |
| C_barriers | strong (1.0) | 0.29, 12 % | 0.28, 31 % | 0.22, 56 % | 249 / 255 / 318 |

Findings.

1. Records of the observed length are more regular than long runs: the
   reference model's windows have median CVs of 0.19 to 0.21 against 0.26
   to 0.30 from whole runs (short records rarely sample the long
   intervals). Judged by windows, the reference calibration is consistent
   with the central (63 to 73 percent of windows at or below 0.23) and
   north records (37 to 45 percent at or below 0.25). The three reference
   copies (independent seeds) agree within a few percent.
2. The south record is the outlier: in the reference cases at most 2
   percent of windows are as irregular as observed. Every case irregular
   enough for the south site (soft clock, scale-free) makes the central
   record improbable (4 to 5 percent), because the two sites share most
   ruptures. The contrast between the south and central records is the
   most diagnostic observation.
3. Fingerprints of the three controls:
   * nucleation clock (A): only the timing changes; frequency-magnitude
     distribution, Omori decay and Bath gap (1.33 to 1.37) are identical;
     whole-fault M >= 7.5 recurrence CV 0.21 to 0.37, mean intervals 268
     to 306 yr;
   * size rule (B): the frequency-magnitude distribution (plateau from
     M6.5 to near Gutenberg-Richter to M7.5), aftershock productivity,
     Bath gap (2.16, 1.34, 1.13), rupture patterns (the scale-free case has
     single-site and no-site M >= 7.5 ruptures) and whole-fault recurrence
     (CV 0.30, 0.21, 0.76) all change; regularity is greatest at the
     reference, not at the most characteristic case;
   * barriers (C, 50 km correlation): moderate increase of site CV (to
     0.28 to 0.29 at south and central), shorter north intervals (318 to
     342 yr) and fewer south-central ruptures; they do not decouple the
     south and central sites (similar CVs at both), contrary to a single
     earlier test run.

## 7. Open points for the comparison with physics-based models

* The observed contrast in regularity between the central (CV 0.23) and
  south (0.40) sites: in every case and ensemble the two sites share most
  ruptures and have similar CVs; judged by record-length windows, the
  reference calibration matches the central and north records but makes
  the south record improbable (at most 2 percent). Neither the size
  rule, the nucleation clock nor 50 km barriers give the two sites partly
  independent rupture histories.
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

## 8. Reproducing

Mesh files in `data/nz_alpine/` (`alpine_planar_1km_26Jan.stl`,
`variable_dip_1km_26Jan.stl`, `Alpine_full_segmented_transition_05km_d20.stl`).

    # reference calibration, five configurations (about 75 core-hours)
    examples/alpine/run_full_campaign.sh $NR_CPUS alpine_full 100 3 10000 20 5 20000 \
        "planar_fine variable_dip_fine segmented_fine planar_coarse variable_dip_coarse" \
        examples/alpine/abc5_planar_l1.json "--site north=1 --bath-sigma 1"
    # regular sites (coarse, about 10 core-hours; add the fine names for all)
    examples/alpine/run_full_campaign.sh $NR_CPUS alpine_r7 100 3 10000 20 5 20000 \
        "planar_coarse variable_dip_coarse" examples/alpine/abc7_regularity.json "--site north=1 --bath-sigma 1"
    # realization sets A, B, C: 450 runs of 20,000 yr at 2 km (about 9 core-hours)
    examples/alpine/run_realization_sets.sh $NR_CPUS alpine_sets 50 20000
    # what each control does (realization sets) and parameter effects (ABC rounds pooled)
    python3 examples/alpine/control_summary.py alpine_sets/controls alpine_sets --site north=1
    python3 examples/alpine/parameter_effects.py alpine_effects alpine_full/abc_*/sweep_summary.json \
        alpine_r6/abc_*/sweep_summary.json alpine_r7/abc_*/sweep_summary.json --site north=1
    # scenarios compatible with the constraints
    python3 examples/alpine/spectrum_figures.py alpine_scenarios \
        reference alpine_full/posterior_variable_dip_fine regular alpine_r7/posterior_variable_dip_coarse \
        --site north=1

The round-5 campaign in `runs/alpine_full_campaign/` was run before
merging was available at run time; its merged tables come from
`site_analysis.py --merge 30` and its summaries were rescored with
`--site north=1`.

## 9. Kept run summaries

* `runs/alpine_full_campaign/`: round 5, five configurations (summaries
  at 2 m and with north at 1 m, site tables, merged tables, figures).
* `runs/alpine_r6/`: round 6 (Bath gap), coarse.
* `runs/alpine_r7/`: round 7 (regular sites), coarse, merged tables and
  figures.
* `runs/alpine_sets/`: realization sets A, B, C (page-2 tables, page-3
  tables, page-1 and page-2 ensemble figures).
