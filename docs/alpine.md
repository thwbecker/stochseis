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
  20 km depth throughout), each at 1 and 0.5 km resolution. Truncated
  and Alpine-Hope geometries have not been received.
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
