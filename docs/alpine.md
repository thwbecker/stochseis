# Alpine Fault simulator intercomparison

Goal: contribute a stochastic baseline to the New Zealand simulator
intercomparison (Alpine and Alpine-Hope exercises, August 2025), tuned to
the paleoseismic recurrence data given the prescribed slip rates, and
report the range of stochastic models consistent with those constraints
so that event statistics of physics-based simulators can be compared
against it. Inputs are described in `data/nz_alpine/README.md`; the
geometry machinery in `docs/networks.md`.

## Inputs

* Geometries (STL, km NZTM): planar (a single plane, 755 km long,
  dipping 50 degrees, 26 km wide to 20 km depth) and variable dip
  (curved trace, 761 km; 80 degrees in the south, 50 degrees in the
  central section, 55 to 65 and then about 80 degrees in the north;
  20 km depth throughout), each at 1 and 0.5 km resolution; segmented
  transition (`Alpine_full_segmented_transition_05km_d20.stl`, 0.5 km
  only): a single unstructured sheet (Gmsh, 139,706 triangles, 17,247
  km^2, 759 km trace) with dips of 80 degrees in the south, 50 in the
  central section, 56 to 65 and then about 85 in the north, changing
  over a few kilometres at the segment ends (Jacksons, Kaniere, Springs
  Junction, Tophouse). The Alpine-Hope geometry has not been received.
* Paleoseismic sites (`data/nz_alpine/paleo_sites.csv`):

  | site | along strike (planar / variable dip) | RI [yr] | sigma [yr] | n events | SED [m] (range) |
  |---|---|---|---|---|---|
  | south | 134 / 125 km | 268 | 108 | 13 | 7.5 (7.0 to 8.0) |
  | central | 350 / 342 km | 249 | 57 | 14 | 7.1 (4.5 to 8.5) |
  | north | 540 / 531 km | 413 | 103 | 4 | 4.0 (1.0 to 9.4) |

* Slip rates and rakes per segment (`data/nz_alpine/segments.csv`, from
  "alpine_hope_segments_and_slip_rates.xlsx"): George-Jacksons 27.2 mm/yr
  (rake 180), Jacksons-Kaniere 27.0 (160), Kaniere-Springs Junction 15.0
  (135), Springs Junction-Tophouse 5.5 (160), Wairau 4.9 (180); Hope
  branch segments 10.9 to 17.3. Segments are projected onto the trace
  (`segments_file`, `segment_include = ["Alpine", "Wairau"]`; the straight
  planar trace needs `segment_max_dist_km = 40`), and rates are blended
  linearly over 10 km either side of segment boundaries
  (`segment_taper_km`). The north site lies 1 to 2 km from the Springs
  Junction boundary, so its rate (9.4 mm/yr planar, 9.9 variable dip)
  depends on that taper.
* Stress orientations and GNSS-derived convergence are given for
  alternative loading strategies; backslip loading does not need them.

The "1 sigma" of the recurrence intervals is the standard deviation of
individual intervals (confirmed), so the observed CVs are 0.40 (south),
0.23 (central) and 0.25 (north). These are among the most regular
paleoseismic series known; San Andreas and San Jacinto sites have CVs of
0.56 to 1.1 (compilation in Becker, Tectonic Geodynamics, Fig. 3.34). A
useful baseline must be able to produce both regimes, and which
parameters separate them is itself a result.

## Budget check

If the events recorded at a site carry most of the slip there, the mean
recurrence interval is close to SED / slip rate. With the observed SEDs:
central 7.1 m / 27 mm/yr = 263 yr (observed 249 +- 57), south 7.5 m
needs about 28 mm/yr for 268 yr, north 4.0 m needs about 10 mm/yr for
413 yr. The north site therefore requires a slip rate well below the
central value, which the missing segment table presumably provides.

## Exploration with the segment slip rates (planar geometry, 2 km patches)

Single seeds of 4000 yr; with 10 to 25 recorded events per site the CVs
carry about +-0.15 sampling error. Values are RI [yr] / CV / SED [m] at
a detection threshold of 2 m; observed: south 268 / 0.40 / 7.5,
central 249 / 0.23 / 7.1, north 413 / 0.25 / 4.0.

| run | changes from defaults | south | central | north |
|---|---|---|---|---|
| v0 | none | 245 / 0.72 / 5.2 | 179 / 0.80 / 4.4 | 458 / 0.60 / 3.0 |
| v1 | stress drop 5 MPa | 191 / 0.74 / 4.7 | 272 / 0.99 / 7.3 | 723 / 0.42 / 6.8 |
| v2 | 5 MPa, nucleation exponent 3, branching 0.3 | 226 / 0.64 / 6.2 | 285 / 0.61 / 7.0 | 1268 / 0.06 / 7.1 |
| v3 | v2 plus 20 km deep slip | 310 / 0.43 / 7.8 | 231 / 0.59 / 7.2 | 762 / 0.77 / 4.7 |
| v6 | 4.5 MPa, exponent 6, branching 0.15, deep 20 km, slip weight 2 | 291 / 0.62 / 7.7 | 290 / 0.61 / 6.1 | (< 3 events) |
| v7 | as v6 without slip weight, arrest scale 1 | 279 / 0.36 / 5.3 | 222 / 0.48 / 4.1 | 636 / 0.62 / 3.9 |
| v9 | v6 plus deficit cap | 272 / 0.39 / 6.0 | 277 / 0.51 / 6.9 | (< 3 events) |
| v10 | v9 plus arrest scale 1 | 311 / 0.57 / 5.7 | 331 / 0.12 / 6.8 | 611 / 0.69 / 5.1 |

Findings:

* With the segment rates the south and central budgets are consistent
  with the data; both sites release 85 to 95 percent of their slip in
  events with at least 2 m of slip (v9).
* Stress drop sets the single-event displacement; 4.5 to 5 MPa (or 3 MPa
  plus deep slip) gives the observed 7 m at the 27 mm/yr sites.
* Regularity needs a sharper renewal dependence of nucleation (exponent
  of 6 or more), low clustering (branching ratio 0.15 or less) and slip
  that follows the local deficit (deficit cap, slip weight 2). Strong
  size control (arrest scale 1) gave the most regular central record
  (CV 0.12) but lengthened intervals elsewhere.
* The north site is the hard constraint. At 9.4 mm/yr it is ruptured by
  most through-going ruptures from the south, each with a small share of
  slip (mean 1.4 m in v9), whereas the data require fewer, larger events
  (4 m every 413 yr). This is a statement about how often ruptures pass
  the Kaniere and Springs Junction boundaries, the behaviour Howarth et
  al. (2021) attribute to geometry, and it should be addressed by the
  rupture-passage rules at segment boundaries rather than by stress drop.

The approximate Bayesian computation (`examples/alpine/abc1_*.json`,
200 random parameter vectors over stress drop, its scatter, nucleation
exponent, branching ratio, deep slip, slip weight, arrest scale and the
deficit cap, two seeds each; `examples/alpine/abc_accept.py`) maps the
accepted region on both geometries.

## ABC round 1, planar geometry (run by T. Becker)

200 random parameter vectors x 2 seeds x 5000 yr at 2 km patches
(`examples/alpine/abc1_planar.json`; summary in
`runs/alpine_abc1_planar/`), no failed runs. Scored with
`examples/alpine/abc_accept.py` (CV error 0.15, RI error the standard
error of the observed mean, SED error half the observed range):

* Detection threshold 2 m: 7 vectors accepted (mean squared standardized
  misfit below 4). The best (misfit 2.0; 3.45 MPa, ln-sigma 0.28,
  nucleation exponent 1.7, branching 0.18, deep slip 19 km, slip weight
  1.6, arrest scale 0.26, deficit cap) gives RI / CV / SED of 296 /
  0.49 / 6.9 (south), 284 / 0.49 / 6.5 (central) and 519 / 0.55 / 3.5
  (north). Accepted ranges are broad: stress drop 2.8 to 6.3 MPa, deep
  slip 4 to 30 km, nucleation exponent 1.2 to 7.6, branching 0.18 to
  0.57; slip weight above 1 and the deficit cap are favoured.
* Threshold 1 m: one vector accepted; threshold 4 m: none (intervals
  too long because too few events reach 4 m at the site).

What controls what (linear regression of the 200 design means on the
standardized parameters):

| output | strongest controls (sign) | R^2 |
|---|---|---|
| SED at south and central | stress drop (+), nucleation exponent (-), arrest scale, cap and slip weight | 0.7 |
| RI at south and central | stress drop (+), slip weight (+), nucleation exponent (-), cap (-) | 0.4 |
| CV at south and central | stress drop (+), nucleation exponent (-), cap (-), all weak | 0.3 to 0.4 |
| RI at north | none of the parameters | 0.09 |
| Bath gap, Omori p, foreshock fraction | branching ratio | 0.8 to 0.97 |
| b above saturation, max M | nucleation exponent and arrest scale (+ / -) | 0.7 to 0.9 |

Findings.

1. Mean recurrence and single-event displacement at the two 27 mm/yr
   sites are matched by a band of stress drop and deep slip; they do not
   constrain the clustering parameters.
2. Regularity has a floor. Over the whole prior the lowest CV at the
   central site (mean of two runs) is 0.37 and no vector has both south
   and central below 0.40. With about 17 recorded events per run the
   sampling standard deviation of a CV near 0.4 is about 0.07, and that
   of the observed central CV (13 intervals) about 0.05, so the central
   record (CV 0.23) lies roughly three standard errors below what this
   model family produces; the south record (0.40) is matched. Extra
   single-seed tests with a size-hazard b of 0.7, a nucleation exponent
   of 20 with a branching ratio of 0.05, and stronger segment boundaries
   did not lower the central CV below about 0.5.
3. The north site is not fitted by any vector: intervals of 460 to 1100
   yr against 413, insensitive to every parameter. Ruptures from the
   south pass Springs Junction too readily and leave little slip there.
   A strength factor at segment boundaries
   (`segment_boundary_strength`, new) is included in round 2.
4. Interpretation for the intercomparison: a stochastic model with
   constant-stress-drop slip, a scale-free size hazard and deficit-
   dependent nucleation reproduces the mean slip and timing of the
   Alpine Fault and the regularity of the south record, but not the
   exceptional regularity of the central record. That regularity
   therefore needs more deterministic coupling between accumulated slip
   and event timing and size than this baseline contains, which is what
   physics-based simulators can be asked to explain. The natural
   counter-hypothesis to test next within stochseis is a slip-predictable
   amplitude rule (event slip equal to the local accumulated deficit
   rather than set by stress drop).

Round 2 (`examples/alpine/abc2_planar.json`): priors narrowed to the
region with acceptable RI and SED, plus the size-hazard b-value (0.7 to
1.1) and the segment-boundary strength (1 to 10), 100 vectors x 3 seeds
x 10,000 yr for smaller CV errors. The variable-dip round 1
(`abc1_variable_dip.json`) is still to be run.

## ABC round 2, planar geometry (run by T. Becker)

100 vectors x 3 seeds x 10,000 yr, priors narrowed to the round-1
region plus the size-hazard b-value and the segment-boundary strength
(`examples/alpine/abc2_planar.json`; summary in
`runs/alpine_abc2_planar/`), no failed runs.

* Threshold 2 m: 4 vectors accepted. Posterior ranges: stress drop 3.1
  to 6.4 MPa (trading off against deep slip 13 to 25 km), ln-sigma 0.03
  to 0.24, nucleation exponent 1.8 to 4.8, branching 0.15 to 0.26, slip
  weight 1.1 to 1.5, arrest scale 0.10 to 0.80, b 0.75 to 0.99,
  boundary strength 1.4 to 2.3. Best: RI / CV / SED 254 / 0.60 / 6.4
  (south), 261 / 0.61 / 6.6 (central), 583 / 0.54 / 4.9 (north).
* With three 10,000-yr seeds the seed-to-seed standard deviation of a
  site CV is 0.06, so the design means are known to about 0.035. The
  lowest central CV over the prior is 0.41 (5th percentile 0.45), about
  three combined standard errors above the observed 0.23 +- 0.05; the
  south record (0.40) remains within reach.
* The nucleation exponent is now the strongest control of CV (-0.07 per
  standard deviation of log exponent), but it also shortens intervals
  and reduces SED, which a higher stress drop has to compensate.
* The north interval (488 to 1100 yr, median 750, against 413) is not
  fixed by the boundary strength, which lengthens it further (+62 yr per
  standard deviation): the north site is ruptured too often with small
  slip, so its slip is fragmented into events below the detection
  threshold. The budget at that site (SED / slip rate) is met only if
  almost all slip comes in about 4 m events.
* Bath (1.04 to 1.41) is set by branching ratio and b; the b-value below
  saturation equals the input b; Omori p 0.70 to 1.01.

Conclusion: with stress-drop-controlled slip, the model family cannot
produce both the central regularity and the north interval.

## Slip-predictable alternative (single-seed tests, 1 November)

`slip_mode = "deficit"` (new): each patch in a rupture releases a fixed
fraction (`slip_fraction`) of its positive deficit; rupture extent still
comes from growth and arrest, and nucleation remains deficit-weighted.
With the round-2 parameters (exponent 2.6, branching 0.2, arrest 0.12,
b 0.8), 5000 yr, threshold 2 m (RI / CV / SED):

| initiation rate (budget_scale) | south | central | north |
|---|---|---|---|
| 1 (calibrated) | 230 / 0.45 / 4.8 | 173 / 0.33 / 4.6 | 349 / 0.25 / 2.9 |
| 0.3 | 269 / 0.42 / 7.0 | 199 / 0.47 / 5.1 | 576 / 0.32 / 4.6 |
| 0.1 | 265 / 0.47 / 6.6 | 294 / 0.55 / 8.2 | 409 / 0.48 / 4.2 |
| 0.1, exponent 6 | 207 / 0.58 / 4.9 | 183 / 0.39 / 4.4 | 517 / 0.31 / 3.9 |
| observed | 268 / 0.40 / 7.5 | 249 / 0.23 / 7.1 | 413 / 0.25 / 4.0 |

In this mode the initiation rate sets how much deficit accumulates before
release and therefore SED and RI; at a tenth of the calibrated rate all
three sites, including the north, match in mean interval and SED. The
CVs (0.25 to 0.58, single seeds with errors of about 0.07 to 0.1) are in
the range of the stress-drop runs or lower. Round 3
(`examples/alpine/abc3_planar_slip_predictable.json`) maps this mode
properly. If it matches the regularity as well, the result for the
intercomparison is that the Alpine records require slip-predictable
behaviour, which constant-stress-drop slip with a scale-free size
distribution does not provide.

## What controls site CV (assessment, 1 November)

Method: single-seed 10,000-yr runs on the planar geometry with catalogs
kept, decomposed with `examples/alpine/decompose_site_cv.py`. For each
site and detection threshold (2 m) it reports the CV of intervals T, the
CV of single-event displacements (SED), the CV of the deficit just
before each recorded event (the "failure deficit"), the ratio of slip to
that deficit, and the correlations corr(SED_i, T_next) (time-predictable
behaviour) and corr(T_prev, SED_i) (slip-predictable behaviour).
Single-seed CVs carry about +-0.07.

### 1. The variance budget

For a site that fails when its deficit reaches some level and releases
part of it, CV(T)^2 is roughly CV(failure deficit)^2 plus CV(released
slip)^2 plus a contribution from slip released in events below the
detection threshold (fragmentation). In the best round-2 model (stress-
drop slip):

| site | CV(T) | CV(SED) | CV(failure deficit) | corr(SED_i, T_next) | slip in recorded events |
|---|---|---|---|---|---|
| central | 0.71 | 0.71 | 0.29 | +0.64 | 94 % |
| south | 0.56 | 0.57 | 0.29 | +0.78 | 93 % |
| north | 0.69 | 0.38 | 0.35 | -0.05 | 77 % |

At the 27 mm/yr sites the model is time-predictable: sites fail at a
fairly consistent deficit but release very different amounts of slip, so
the interval CV equals the slip CV. At the north site, a quarter of the
slip is released in events below the threshold.

### 2. Sources of released-slip variability, and what removes them

* Slip roughness at the site (von Karman pattern, CV 0.75) and stress-
  drop scatter: reducing them to 0.2 and 0.05 lowers central CV(SED)
  from 0.71 to 0.44 and CV(T) to 0.46. Not varied in rounds 1 and 2.
* Crack taper when a rupture ends near the site (random rupture extents
  under scale-free sizes): remains in the south (CV 0.62) after
  smoothing.
* Partial release: in stress-drop mode a rupture releases 50 to 80 % of
  the local deficit, with large scatter.
* Releasing the local deficit (slip-predictable modes) removes most of
  this, but the plain slip-predictable mode has a flaw: any small
  rupture that includes a loaded site releases its whole deficit, so
  moderate events reset sites at random times. The new `local_cap` mode
  (stress-drop crack slip capped patch by patch at the local deficit,
  with a high stress drop so that large ruptures reach the cap) keeps
  small events stress-drop limited.

### 3. The remaining floor: variability of the failure deficit

With released slip tied to the deficit, CV(T) approaches CV(failure
deficit), which stays at about 0.3 to 0.4 whatever else is done:

| variant (single seeds) | central CV(T) | south | north | central CV(failure deficit) |
|---|---|---|---|---|
| characteristic (no size hazard), local cap, exponent 3 | 0.35 | 0.35 | 0.41 | 0.36 |
| as above, nucleation exponent 10 | 0.34 | 0.36 | 0.33 | 0.29 |
| as above, deficit smoothed over 20 km | 0.41 | 0.32 | 0.29 | 0.34 |
| as above, deficit smoothed over 60 km | 0.38 | 0.42 | 0.46 | 0.34 |
| sharp admission threshold (exponent 8), local cap, scale-free sizes | 0.58 | 0.80 | 0.70 | 0.41 |

Its cause is along-strike phase heterogeneity: ruptures end at different
places, leave deficit steps, and neighbouring stretches of the fault
drift out of phase, so a site joins ruptures at different points of its
own cycle. Sharper nucleation and admission rules do not remove it, nor
does smoothing the deficit as a stand-in for elastic coupling. The
likely missing ingredient is genuine stress interaction (load transfer
to rupture ends that pulls stretches back into phase), which the level-2
half-space operator of `docs/networks.md` would provide; this is
untested.

### 4. Hybrid sizes: Gutenberg-Richter below, characteristic above

`hazard_max_area_km2` (new) applies the scale-free size hazard only to
ruptures smaller than a given area; larger ruptures stop only by state.
With 700 km^2 (about width-saturated M7), local-cap slip, 8 MPa, slip
CV 0.3, branching 0.1, admission scale 0.5 with exponent 4, at 0.3 of
the calibrated initiation rate (2 m threshold, RI / CV / SED):

| nucleation exponent | south | central | north | FMD |
|---|---|---|---|---|
| 3 | 318 / 0.28 / 8.4 | 286 / 0.43 / 7.5 | 692 / 0.32 / 6.9 | GR (b = 1.00) to M7, x2 at M7.5, x5 at M8 |
| 8 | 235 / 0.38 / 6.1 | 276 / 0.33 / 6.6 | 549 / 0.35 / 4.9 | GR (b = 0.97) to M7, x2.3 at M7.5, x5.4 at M8 |

The hybrid keeps Gutenberg-Richter below M7, gives site CVs of 0.28 to
0.43, and switches the sites from time-predictable to slip-predictable
behaviour (central: corr(T_prev, SED_i) = +0.77, corr(SED_i, T_next) =
+0.15). The north interval is still too long.

### 5. How far from the data is the floor?

The observed CVs come from 13 to 14 intervals. A gamma-renewal estimate
of their sampling error is 0.05 (central) and 0.09 (south), but
bootstrapped uncertainties of Alpine Fault records are larger (0.10 for
Hokuri Creek and 0.13 for Lake Paringa in Becker, Tectonic Geodynamics,
Fig. 3.34). Against those, a model floor of 0.30 to 0.35 is about one
standard error above the central record, not three as stated for round 2.

### 6. Conclusions

1. With stress-drop slip and scale-free sizes (rounds 1 and 2), the site
   CV is dominated by point-slip variability and cannot fall below about
   0.4.
2. Getting closer to the Alpine values needs two changes together: slip
   at a site tied to its accumulated deficit, with low roughness and
   stress-drop scatter; and large ruptures that span the loaded segment
   (characteristic at the segment scale), which the hybrid size rule
   provides while keeping GR below M7.
3. The remaining floor of about 0.3 comes from along-strike phase
   heterogeneity and probably needs elastic interaction to remove.
4. The decomposition gives a testable signature: the stress-drop
   baseline predicts time-predictable records (displacement correlates
   with the following interval), the hybrid predicts slip-predictable
   records (displacement correlates with the preceding interval). Where
   paleoseismic records resolve per-event displacement, this
   discriminates between the two.
5. Round 4 (`examples/alpine/abc4_planar_hybrid.json`) maps the hybrid
   mode with the ABC machinery; round 3 (slip-predictable) is running.

## ABC round 4: hybrid sizes with deficit-capped slip (run by T. Becker)

100 vectors x 3 seeds x 10,000 yr, planar geometry
(`examples/alpine/abc4_planar_hybrid.json`; summary in
`runs/alpine_abc4_planar/`), no failed runs. Slip mode `local_cap`,
hybrid size rule with the area threshold drawn from 300 to 3000 km^2.

Acceptance (mean squared standardized misfit below 4): 55 of 100 at a
2 m detection threshold, 48 at 1 m, 1 at 4 m. Best vector at 2 m (misfit
0.57): stress drop 7.6 MPa, ln-sigma 0.19, slip CV 0.15, nucleation
exponent 4.4, branching 0.27, arrest scale 0.29 with exponent 3.8, b 0.86,
initiation rate 0.89 x calibrated, threshold area 810 km^2:

| site | RI [yr] | CV | SED [m] | observed |
|---|---|---|---|---|
| south | 246 | 0.31 | 6.6 | 268 / 0.40 / 7.5 |
| central | 255 | 0.26 | 6.8 | 249 / 0.23 / 7.1 |
| north | 458 | 0.28 | 4.0 | 413 / 0.25 / 4.0 |

Accepted ranges (2 m): stress drop 5.2 to 10 MPa, ln-sigma up to 0.2,
slip CV 0.12 to 0.50, nucleation exponent 2.6 to 12, branching 0.05 to
0.30, arrest scale 0.23 to 0.99, admission exponent 2 to 8, b 0.80 to
1.10, initiation rate 0.10 to 0.95 x calibrated, threshold area 310 to
2900 km^2. Site CVs over the accepted set: medians 0.38 to 0.40, lowest
0.16 (central) to 0.23 (north); over the whole prior the 5th percentile
of the central CV is 0.23. The CV floor of rounds 1 and 2 is gone.

Controls (regression over the 100 design means, per standard deviation
of each parameter, R^2 0.86 to 0.90 for the south and central sites):

* CV: threshold area (+0.08, the strongest: the smaller the area above
  which ruptures stop only by state, the more regular), arrest scale
  (+0.04), stress drop (+0.03 to 0.04), admission exponent (-0.02).
* RI and SED: stress drop (+), nucleation exponent (-), initiation rate
  (-), b (-).
* North RI: additionally arrest scale (+54 yr per standard deviation).

Event statistics of the accepted ensemble (medians, 5 to 95 percent):
Bath 1.41 (1.16 to 1.88), Omori p 1.01 (0.88 to 1.09), b below
saturation 0.91 (0.80 to 1.08), maximum M 8.31 (8.13 to 8.41), 1190 M>=6
and 139 M>=7 per 10,000 yr with N(M>=7) / N(M>=6) = 0.119, as for GR
with the accepted b (the characteristic excess is above M7.5), whole-fault
M>=7 CV 0.90, realised branching ratio 0.17. Co-rupture at 1 m: 42
percent of central events also rupture the south site, 19 percent the
north site; 7 percent of south events reach the north site.

Interpretation: with slip tied to the local deficit and large ruptures
terminated by fault state, the stochastic model reproduces mean
recurrence, regularity and single-event displacement at all three sites
with Gutenberg-Richter statistics below about M7, a characteristic
excess above M7.5, and standard aftershock statistics. Many parameter
combinations do this (55 of 100 draws), so the paleoseismic data
constrain the model class strongly but individual parameters weakly.
Round 3 (slip-predictable mode) results have not been received.

Next: the same round on the variable-dip geometry
(`examples/alpine/abc4_variable_dip_hybrid.json`) and a posterior
ensemble with catalogs for the intercomparison statistics
(`examples/alpine/posterior4_planar.json`: the 20 best vectors x 5 seeds x
20,000 yr; regenerate it for other summaries with
`examples/alpine/make_posterior_spec.py`).

## ABC round 3: slip-predictable mode (run by T. Becker)

100 vectors x 3 seeds x 10,000 yr, planar geometry, `slip_mode =
"deficit"` (`examples/alpine/abc3_planar_slip_predictable.json`; summary
in `runs/alpine_abc3_planar/`), no failed runs.

On the paleoseismic observables alone, 29 vectors are accepted at 2 m
(22 at 1 m, none at 4 m); site CVs are 0.36 to 0.75 (lowest central CV
0.35), worse than round 4. The seismicity, however, is unphysical: the
b-value of grown ruptures is about 0.03 (5 to 95 percent: 0 to 0.13),
the mean Bath gap 0.71, N(M>=7) / N(M>=6) = 0.066, about five times as
many M>=6 events as in round 4, and a whole-fault M>=7 CV of 1.17. This
is the flaw anticipated in "What controls site CV": when slip equals the
local deficit, a small rupture through a loaded patch releases metres of
slip, so magnitude decouples from rupture size, the frequency-magnitude
distribution flattens, and aftershocks are often as large as their
mainshocks. The plain slip-predictable mode is therefore rejected.

## Rescoring with seismicity constraints

`examples/alpine/abc_accept.py ... --seismicity` adds the b-value of
grown ruptures below width saturation (1.0 +- 0.1) and the mean Bath gap
(1.2 +- 0.15) to the misfit. Accepted vectors at 2 m (mean squared
misfit below 4):

| round | slip rule | paleo only | paleo and seismicity |
|---|---|---|---|
| 2 | stress drop, scale-free sizes | 4 | 5 (site CVs 0.5 to 0.6) |
| 3 | slip-predictable | 29 | 0 |
| 4 | hybrid sizes, deficit-capped slip | 55 | 54 |

Within the round-4 set accepted with seismicity: site CV medians 0.40 to
0.41 (lowest 0.16 central, 0.20 south, 0.23 north), b 0.80 to 1.08, Bath
1.15 to 1.75 (median 1.36). The five vectors with all three site CVs at
or below 0.35 all have Bath gaps of 1.46 to 1.89: in this model class,
Alpine-like regularity comes with weak triggering of large aftershocks.
Whether the Alpine Fault's instrumental seismicity supports a large Bath
gap is a check worth making, and a prediction to compare with the
physics-based simulators.

Conclusion across rounds 1 to 4: of the slip rules tried, only deficit-
capped slip with hybrid sizes (scale-free below a threshold area,
characteristic above) satisfies both the paleoseismic records and basic
seismicity statistics.

## ABC round 4 on the variable-dip geometry (run by T. Becker)

Same 100 parameter vectors (same design seed) x 3 seeds x 10,000 yr on
the variable-dip geometry (`examples/alpine/abc4_variable_dip_hybrid.json`;
summary in `runs/alpine_abc4_variable_dip/`), no failed runs. Accepted at
2 m: 61 vectors on the paleoseismic observables, 69 with the seismicity
constraints (planar: 55 and 54). Best vector (the same as on the planar
geometry): RI / CV / SED 259 / 0.31 / 6.9 (south), 258 / 0.32 / 6.9
(central), 489 / 0.28 / 4.4 (north), b 0.89, Bath 1.33. Accepted ranges
are essentially those of the planar round.

Because the vectors are identical, the geometries can be compared design
by design. Paired mean differences, variable dip minus planar, with
t-statistics over the 100 designs (|t| > 3 marked as systematic):

| statistic | planar median | variable dip median | mean difference | t |
|---|---|---|---|---|
| moment rate [N m/yr] | 1.13e19 | 0.99e19 | -13 % | (geometry) |
| N(M>=6), N(M>=7) per 10,000 yr | 1230, 141 | 1020, 119 | -18 %, -21 % | -9.3, -8.4 |
| co-rupture central with south | 0.43 | 0.60 | +0.13 | 14.5 |
| co-rupture central with north | 0.20 | 0.28 | +0.08 | 14.1 |
| co-rupture south with north | 0.08 | 0.21 | +0.07 | 11.4 |
| north RI [yr] | 609 | 554 | -53 | -5.9 |
| south CV / north CV | 0.41 / 0.41 | 0.44 / 0.45 | +0.04 / +0.03 | 6.5 / 3.5 |
| central CV | 0.41 | 0.40 | +0.005 | 0.7 |
| central SED [m] | 6.7 | 7.2 | +0.27 | 3.5 |
| Bath gap | 1.39 | 1.31 | -0.08 | -9.6 |
| b above saturation | 1.13 | 0.96 | -0.23 | -4.8 |
| median D/L, L > 3W | 1.8e-5 | 1.1e-5 | -0.9e-5 | -13 |
| whole-fault M>=7 CV | 0.96 | 0.87 | -0.06 | -6.9 |

Findings.

* The paleoseismic records are matched equally well on both geometries,
  by the same parameter vectors; per-design site CVs correlate between
  the geometries (r = 0.82 central, 0.83 south, 0.60 north). The site
  data therefore do not discriminate between the two geometries in this
  model.
* Geometry shows most clearly in how ruptures span the fault: on the
  variable-dip geometry many more events rupture two or three sites
  (south-north co-rupture nearly triples), there are fewer but longer
  large events, and the north site is ruptured more often (shorter RI).
  Co-rupture fractions are thus the statistic most sensitive to
  geometry, and a natural one to compare with the physics-based
  simulators and with paleoseismic correlations between sites (Howarth
  et al., 2021).
* Caveat: the elastic operator is level L0 (unrolled grid, free-surface
  image), which ignores the effect of dip on stiffness beyond the change
  in down-dip width (`docs/networks.md`, E20). Differences in slip
  amplitude between the 50 and 80 degree sections are therefore only
  partly represented; the L1 dip correction or the L2 half-space
  operator should be used before quoting geometry effects on slip.

## L1 dip correction and the rerun campaign (round 5)

`dip_correction = true` (new, `stochseis/dipcorr.py`) scales each
crack solution by the ratio of half-space to spectral-image mean slip for
the rupture's mean dip, rake and burial, tabulated once from cutde
solutions for 60 x 20 km ruptures (`examples/checks/l1_table.py`):

| rupture | 15 deg | 30 | 45 | 60 | 75 | 90 |
|---|---|---|---|---|---|---|
| strike slip, surface | 1.91 | 1.32 | 1.13 | 1.04 | 1.00 | 0.99 |
| strike slip, buried (top 3 km) | 1.22 | 1.13 | 1.09 | 1.07 | 1.06 | 1.05 |
| dip slip, surface | 1.49 | 1.18 | 1.14 | 1.18 | 1.24 | 1.27 |
| dip slip, buried | 1.18 | 1.16 | 1.16 | 1.17 | 1.18 | 1.18 |

Oblique slip weights the two by cos^2 and sin^2 of the rake. For the
Alpine Fault (rake 135 to 180, dip 50 to 85 degrees) the factor is 1.0
to 1.1; for shallow thrusts it is much larger. It corrects the moment
and mean slip of a rupture, not the shape of its slip distribution near
the trace, and it adds no cost. The default stays off so that earlier
single-fault results are unchanged; the Alpine specs switch it on.

Round 5 repeats round 4 with the correction and the same design seed, so
every design pairs with round 4 and across geometries. The round-5 spec
(`examples/alpine/abc5_planar_l1.json`) is the base of the full campaign
(below), whose posterior ensembles are drawn with
`examples/alpine/overview_figures.py`:

* `<outdir>/fig_seismicity.png`: frequency-magnitude distribution, Omori
  decay after M >= 6 mainshocks, Bath gap distribution, moment release
  relative to loading; best vector as a line, 5 to 95 percent of all
  posterior runs as a band, one colour per geometry.
* `<outdir>/fig_recurrence.png`: RI versus CV at each site for every
  posterior run with the observations, SED versus RI, the distribution
  of rupture lengths of M >= 6.5 events, and magnitude versus rupture
  length.

`docs/alpine_overview_*_example.png` show the figures for a small test
ensemble (2 vectors x 2 seeds x 3000 yr per geometry) and are not
results. Already visible in them: through-going ruptures of 500 to
700 km form a separate population from ruptures shorter than about
100 km, and system-wide moment release is quasi-periodic with a period
of about 250 yr.

## Segmented-transition geometry and a parametrisation fix (round 5)

The segmented-transition mesh is unstructured: within one along-strike
column, vertices sit at many depths, and where the dip changes along
strike their offsets from the trace scatter. The down-dip coordinate was
built from a per-vertex profile, which zig-zagged and gave down-dip
widths of up to 107 km at 2 km resolution. It now uses the median offset
in depth bins of half a patch. The segmented geometry resamples to 4227
patches at 2 km (98.5 percent of the mesh area, widths 8 to 26 km); the
planar grid is unchanged; the variable-dip grid changes from 4294 to 4260
patches at 2 km (maximum width 26.2 instead of 33.8 km), so a few
transition columns had about 0.8 percent too much area in the
variable-dip runs of rounds 1 to 4. Round 5 uses the corrected grids.

A first run on the segmented geometry with the best round-4 vector and
the L1 correction (single seed, 5000 yr, 2 m threshold): RI / CV / SED
249 / 0.25 / 6.7 (south), 221 / 0.39 / 6.0 (central), 454 / 0.21 / 4.3
(north); maximum M 8.30, Bath 1.43, closure 1.01.

The full campaign (`run_full_campaign.sh`, below) covers all three
geometries.

## Site analysis in the style of the RSQSim comparison

`examples/alpine/site_analysis.py out_prefix label dir [label dir ...]
[--threshold 2]` reproduces the panels of the RSQSim summary for the
intercomparison (planar and varying-dip geometries, 50,000 yr; RSQSim
site values in `data/nz_alpine/rsqsim_summary.csv`) for a posterior
ensemble, and extends them:

* `_spacetime.png`: slip-rate map with the sites and a space-time plot of
  M >= 7 ruptures (along-strike extent in northing, hypocentre coloured by
  Mw) for the best run, one row per geometry.
* `_statistics.png`: the RSQSim panels (slip, magnitude and inter-event
  time distributions at the sites with the paleoseismic SED ranges and
  normal interval PDFs; cumulative frequency-magnitude distribution;
  hypocentre depths; rupture area versus Mw), computed over the whole
  ensemble with the best run and a 5 to 95 percent band where useful;
  plus RI versus CV and RI versus SED for every run with the observations
  and the RSQSim values, and the mean interval as a function of the
  detection threshold (0.5 to 4 m).
* `_extended.png`: co-rupture matrices, time- and slip-predictability
  correlations, and interval distributions against the paleoseismic
  normal CDFs.
* `_table.md` / `_table.csv`: per geometry and site, mean and SD of the
  intervals, CV, mean and SD of SED for the best run, ensemble median and
  5 to 95 percent range, observations, and RSQSim.

The full campaign runs it after the posteriors. For posteriors already
computed:

    python3 examples/alpine/site_analysis.py out/sites planar_fine out/posterior_planar_fine \
        variable_dip_fine out/posterior_variable_dip_fine --threshold 2

Example output for a small test ensemble (round-4 best three vectors with
the L1 correction, 2 seeds x 8000 yr per geometry) is in
`docs/alpine_site_analysis_*_example.*`. These are not the campaign
results, but the comparison with RSQSim already shows where the models
differ:

| | stochseis (test ensemble, median) | RSQSim | observed |
|---|---|---|---|
| central RI / CV / SED, planar | 264 / 0.32 / 7.0 | 271 / 0.13 / 7.3 | 249 / 0.23 / 7.1 |
| south RI / CV / SED, planar | 258 / 0.33 / 6.9 | 260 / 0.19 / 7.1 | 268 / 0.40 / 7.5 |
| north RI / CV / SED, planar | 492 / 0.34 / 4.3 | 640 / 0.45 / 5.8 | 413 / 0.25 / 4.0 |
| central RI / CV / SED, variable dip | 250 / 0.28 / 6.7 | 272 / 0.15 / 7.4 | |
| south RI / CV / SED, variable dip | 235 / 0.32 / 6.3 | 214 / 0.24 / 5.8 | |
| north RI / CV / SED, variable dip | 468 / 0.28 / 4.3 | 543 / 0.36 / 5.4 | |

* Both models match the mean interval and displacement at the south and
  central sites. RSQSim is more periodic than the data (CV 0.13 to 0.24
  against 0.23 and 0.40); the stochastic ensemble brackets the
  observations.
* At the north site both are too slow, RSQSim more so (543 to 640 yr and
  5.4 to 5.8 m against 413 yr and 4.0 m; stochseis about 470 to 490 yr
  and 4.3 m).
* Size distribution: events recorded at the sites are almost all
  through-going M8.1 to 8.4 ruptures in stochseis, against M7.5 to 8.2 in
  RSQSim; the space-time plots show the same, a characteristic
  population of fault-length ruptures in stochseis (the hybrid size rule)
  and more partial M7.5 to 8 ruptures in RSQSim. Both frequency-magnitude
  distributions depart from Gutenberg-Richter above about M7.
* Hypocentres are spread over depth on the planar geometry in stochseis
  but concentrate at the base of the fault on the variable-dip geometry;
  RSQSim nucleates mostly between 13 and 20 km on both.

## Full campaign: three geometries at fine and two at coarse resolution

Configurations (`examples/alpine/campaign_configs.json`):

| configuration | mesh | patches | grid area / mesh area |
|---|---|---|---|
| planar_fine | alpine_planar_1km_26Jan.stl (1 km) | 1 km, 19,630 | 0.995 |
| variable_dip_fine | variable_dip_1km_26Jan.stl (1 km) | 1 km, 17,127 | 0.992 |
| segmented_fine | Alpine_full_segmented_transition_05km_d20.stl (0.5 km) | 1 km, 16,995 | 0.993 |
| planar_coarse | alpine_planar_1km_26Jan.stl (1 km) | 2 km, 4,901 | 0.994 |
| variable_dip_coarse | variable_dip_1km_26Jan.stl (1 km) | 2 km, 4,260 | 0.986 |

"Fine" and "coarse" refer to the simulation patch size. The 1 km meshes
resample to the same 1 km grids as the 0.5 km meshes (19,630 and 17,127
against 19,630 and 17,125 patches), so the original 1 km meshes are used
at both resolutions. All earlier rounds correspond to the coarse
configurations. `alpine_planar.stl` and `variable_dip_1km.stl` (in
metres) are the same meshes as the 1 km files (in km).

`examples/alpine/run_full_campaign.sh [ncpu] [outdir] [n_samples]
[abc_seeds] [abc_duration] [n_best] [post_seeds] [post_duration]
["configurations"]` (defaults: 100 vectors x 3 seeds x 10,000 yr for the
ABC, best 20 vectors x 5 seeds x 20,000 yr for the posteriors, all five
configurations, fine first) does, for each configuration:

1. writes the ABC spec from `abc5_planar_l1.json` (round-5 priors, L1
   correction, design seed 4: identical parameter vectors in every
   configuration) with the configuration's mesh and patch size
   (`make_campaign_specs.py`);
2. runs the ABC sweep and writes acceptance tables with and without the
   seismicity constraints;
3. ranks the vectors with the seismicity constraints and runs the
   posterior ensemble with catalogs;

and then writes `campaign_summary.md` (acceptance counts, best fits,
accepted parameter ranges and overlap of accepted vectors between
configurations; `campaign_summary.py`) and two kinds of figures:

* our overview figures (`overview_figures.py`, `fig_<set>_seismicity.png`
  and `fig_<set>_recurrence.png`): frequency-magnitude distribution,
  Omori decay, Bath gaps, moment release relative to loading, moment
  released per 100-yr window, and whole-fault recurrence of M >= 7.5; then
  RI versus CV at each site for every run, SED versus RI, mean interval
  versus CV of system-wide M >= 7.5 and M >= 8 events, magnitude
  distribution of M >= 6.5 events, rupture-length distribution, magnitude
  and mean slip versus rupture length. Best vector as line or large
  symbol, all posterior runs as band or small symbols, one colour per
  configuration;
* the RSQSim-style site analysis (`site_analysis.py`, `sites_<set>_*`).

Sets: the fine configurations (`fine`), the coarse ones (`coarse`), all
five together (`all`, overview figures only), and fine against coarse for
each geometry (`resolution_planar`, `resolution_variable_dip`). The scoring used by all scripts is in
`abc_scoring.py`.

Cost on one core (measured): 2.6 min per 10,000 yr at 1 km patches,
0.6 min at 2 km. Per fine configuration about 13 core-hours for the ABC
and 9 for the posterior; per coarse configuration about 3 and 2; about
75 core-hours in total.

## Full campaign results (run by T. Becker, 3 November)

Five configurations, 100 vectors x 3 seeds x 10,000 yr each for the ABC
and the best 20 vectors x 5 seeds x 20,000 yr for the posteriors, no
failed runs. Tables and acceptance lists in `runs/alpine_full_campaign/`.

### Acceptance

| configuration | accepted (paleo) | accepted (paleo + seismicity) | best misfit (with seismicity) |
|---|---|---|---|
| planar_fine | 54 | 53 | 0.85 |
| variable_dip_fine | 57 | 59 | 0.64 |
| segmented_fine | 57 | 60 | 0.60 |
| planar_coarse | 52 | 56 | 0.99 |
| variable_dip_coarse | 64 | 65 | 0.48 |

The accepted sets overlap strongly: of the 53 to 65 vectors accepted in
each configuration, 45 to 58 are accepted in any other one. The same few
vectors rank at the top everywhere (for example stress drop 8.0 MPa,
nucleation exponent 3.3, branching 0.29, arrest scale 0.27, area
threshold 990 km^2; and 7.4 MPa, exponent 4.0, branching 0.22, arrest
0.27, area 370 km^2).

The accepted parameter distributions are almost the prior distributions
of this round: the medians (stress drop 7.0 to 7.7 MPa, nucleation
exponent 4.9 to 5.3, branching ratio 0.17 to 0.18, arrest scale 0.67,
admission exponent 4.4 to 4.6, b 0.94 to 0.96, initiation rate 0.32 to
0.40 of calibrated, slip CV 0.32 to 0.34, area threshold 1080 to 1140
km^2) differ from the prior medians (7.5, 4.9, 0.175, 0.6, 5.0, 0.95,
0.32, 0.30, 950) by little, and the ranges cover the priors. The priors
of round 5 were themselves narrowed after rounds 1 to 4, so this says
that within the hybrid, deficit-capped model class the paleoseismic and
seismicity constraints accept a broad region of parameters, and do not
pin down any single parameter beyond the earlier rounds. What the data
select is the model class, not a parameter vector.

### Site statistics of the posterior ensembles (2 m threshold)

Ensemble medians and 5 to 95 percent ranges over all posterior runs:

| configuration | south RI / CV / SED | central RI / CV / SED | north RI / CV / SED |
|---|---|---|---|
| planar_fine | 259 [234-292] / 0.38 [0.27-0.52] / 6.6 | 259 [234-301] / 0.38 [0.27-0.52] / 6.6 | 544 [424-659] / 0.36 / 4.8 |
| variable_dip_fine | 254 [216-288] / 0.40 [0.31-0.53] / 6.5 | 264 [228-312] / 0.39 [0.30-0.52] / 6.7 | 543 [438-683] / 0.38 / 4.9 |
| segmented_fine | 256 [221-287] / 0.41 [0.31-0.54] / 6.7 | 266 [231-314] / 0.38 [0.29-0.53] / 6.8 | 520 [425-654] / 0.40 / 4.7 |
| planar_coarse | 252 [218-291] / 0.37 [0.25-0.51] / 6.6 | 257 [226-295] / 0.39 [0.25-0.50] / 6.6 | 545 [422-657] / 0.37 / 4.7 |
| variable_dip_coarse | 251 [224-276] / 0.41 [0.31-0.53] / 6.6 | 260 [230-309] / 0.40 [0.31-0.52] / 6.8 | 511 [413-642] / 0.39 / 4.7 |
| observed | 268 / 0.40 / 7.5 | 249 / 0.23 / 7.1 | 413 / 0.25 / 4.0 |
| RSQSim planar / variable dip | 260 / 0.19 / 7.1; 214 / 0.24 / 5.8 | 271 / 0.13 / 7.3; 272 / 0.15 / 7.4 | 640 / 0.45 / 5.8; 543 / 0.36 / 5.4 |

1. Mean intervals at the south and central sites are matched in every
   configuration; single-event displacements are 0.4 to 0.9 m below the
   observed means, with the observations inside or at the upper edge of
   the 5 to 95 percent ranges.
2. The south CV (0.40) is matched. The central CV (0.23) lies at or
   below the 5th percentile of every ensemble (0.25 to 0.31), and the
   ensemble medians (0.38 to 0.41) are close to the typical CV of the
   model class found in rounds 4 and 5. For the top-ranked vector the
   central CV changes from the ABC (three seeds) to its first posterior
   run as 0.37 to 0.38 (planar fine), 0.41 to 0.44 (variable dip fine),
   0.27 to 0.37 (segmented), 0.31 to 0.30 (planar coarse) and 0.41 to
   0.37 (variable dip coarse). Only the segmented case shows a clear
   increase, and single-run CVs carry about +-0.07, so selection bias
   (ranking 100 vectors by short-run statistics favours lucky ones) is
   plausible but not demonstrated by these numbers.
   `campaign_summary.py` now rescores the posterior vectors with all
   their posterior runs; those misfits are the ones to quote.
3. The north site is too slow in every configuration (medians 511 to 545
   yr against 413, with the observation at the 5th percentile), with a
   mean SED of 4.7 to 4.9 m against 4.0 m. RSQSim has the same bias, more
   strongly (543 to 640 yr, 5.4 to 5.8 m).
4. Resolution: fine and coarse ensembles agree within about 10 yr in RI,
   0.02 in CV and 0.2 m in SED at every site, so site statistics are
   converged at 2 km patches.
5. Geometry: differences between planar, variable-dip and segmented
   ensembles are of the same size as the resolution differences; the
   site records do not discriminate between the geometries in this
   model, as found in round 4. Geometry effects on co-rupture and on
   the size distribution of system events are in the figures
   (`fig_*`, `sites_*_extended.png`), which are not summarized here.

### Comparison with RSQSim

At the south and central sites both models reproduce the mean intervals
and displacements. RSQSim is more periodic (CV 0.13 to 0.24) than the
stochastic ensembles (0.37 to 0.41): RSQSim matches the central record
(0.23) better and the south record (0.40) worse. The stochastic model's
single-event displacements are more variable at a site (SD 1.7 to 2.6 m
for the best vectors against 1.0 to 1.6 m in RSQSim). At the north site
both models are too slow with too much slip per event, RSQSim more so.
A stochastic model calibrated to the paleoseismic records thus spans the
range between the RSQSim regularity and the more irregular records, and
the central and north records are the ones that discriminate.

## Status (1 November 2026)

* Both geometries load as section grids at 1 km (19,630 and 17,142
  patches; 99.5 and 99.2 percent of the mesh area) and run at about
  1 minute per 600 yr on one core with mag_min 4.5.
* Elastic operator: level L0 of `docs/networks.md`, i.e. the unrolled
  grid with a free-surface image. For a 50 degree dip this gives about
  10 to 15 percent too little slip (E20); the L1 correction is next.
* First 6000-yr runs with the placeholder slip rate: see the results
  section below once complete.

## Plan

1. Obtain the segment slip-rate table and the remaining geometries.
2. L1 dip correction; general rake.
3. Baseline runs on both geometries with the defaults; site statistics
   at detection thresholds of 1, 2 and 3 m.
4. Approximate Bayesian computation over stress-drop median and scatter,
   nucleation exponent, branching ratio and detection threshold: a few
   hundred 5000-yr runs per geometry, accepted where RI, CV and SED at
   the three sites fall within the observational uncertainties.
5. For the accepted ensemble: event statistics for the intercomparison
   (magnitude-frequency, rupture lengths, fraction of ruptures crossing
   section boundaries, co-rupture of sites, whole-fault recurrence).
