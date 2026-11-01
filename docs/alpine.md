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

* Slip rates: a segment table is referenced in the exercise description
  but has not been received. The description gives 27 mm/yr of strike
  slip for the central section (Barth et al., 2013) and a rake of about
  160 degrees. Runs so far use 27 mm/yr everywhere as a placeholder.
* Stress orientations and GNSS-derived convergence are given for
  alternative loading strategies; backslip loading does not need them.

## Budget check

If the events recorded at a site carry most of the slip there, the mean
recurrence interval is close to SED / slip rate. With the observed SEDs:
central 7.1 m / 27 mm/yr = 263 yr (observed 249 +- 57), south 7.5 m
needs about 28 mm/yr for 268 yr, north 4.0 m needs about 10 mm/yr for
413 yr. The north site therefore requires a slip rate well below the
central value, which the missing segment table presumably provides.

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
