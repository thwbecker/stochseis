# Synthesis: what purely stochastic fault models produce

Status: 5 November 2026. This document draws together the single-fault
experiments (`docs/experiments.md`, E1 to E22) and the Alpine Fault
calibration (`docs/alpine.md`). It states what the model class produces,
which parameters control its complexity, which Alpine cases are
consistent with the paleoseismic constraints, and how predictions for the
model comparison are obtained. Numbers marked provisional rest on small
ensembles and are to be replaced by the full realization sets.

## 1. The model class

stochseis generates earthquake sequences on a fault from rules without
fault mechanics: a slip budget (backslip loading at prescribed rates; the
slip deficit is the state), deficit-weighted nucleation, ETAS clustering,
rupture growth by flooding with a size rule and deficit-dependent
admission, crack-shaped slip set by stress drop, with roughness, and
static elasticity only through crack solutions (spectral operator with a
dip correction). Everything that physics-based simulators obtain from
friction, stress transfer and dynamics is replaced by a few rules with
tunable parameters, so the model shows which observations such rules
already explain.

## 2. Single fault (Stage 0, planar 200 km strike-slip fault)

With the default parameters (stress-drop slip, scale-free sizes,
branching ratio 0.5; E22, 20 seeds x 10,000 yr):

* Gutenberg-Richter statistics by construction (b = 1.01 below width
  saturation, 1.28 above), maximum M7.6, Omori p 0.96, Bath gap 1.24,
  foreshocks before 33 percent of mainshocks (15 percent by chance),
  moment closure within 1 percent.
* Irregular recurrence: site CV 0.62 to 0.73, whole-fault M >= 7 CV 0.74;
  the Bath gap and site CV move in opposite directions with the
  branching ratio (E6).
* Moment-length scaling follows L^1.5 for long ruptures; L^2 needs slip
  below the seismogenic depth (E18). Results are independent of patch
  size above the patch magnitude (E13) and of the initial state (E14).
* A renewal-like nucleation clock lowers site CV (0.62 to 0.54, E9);
  state-only rupture termination gives characteristic size distributions
  (E10, E17): Gutenberg-Richter does not emerge from the state alone and
  has to be imposed.

## 3. What controls complexity

Across the single fault and the Alpine runs, the parameters fall into
groups that act on different aspects of a sequence:

| control | what it changes | what it leaves alone |
|---|---|---|
| slip rule (stress-drop crack slip, slip equal to the deficit, crack slip capped by the deficit) | whether site regularity is possible: with uncapped crack slip, slip released at a site varies from event to event and site CV stays near 0.4 or above; slip equal to the deficit destroys the magnitude-size relation (b about 0); capped crack slip allows regular sites with normal seismicity | |
| stress drop and initiation rate | mean recurrence interval and displacement (slip budget) | regularity, to first order |
| nucleation clock (exponent of the deficit dependence) | timing regularity: site and whole-fault CV | size distribution, Omori decay, Bath gap |
| size rule (area above which ruptures stop only by the fault state) | the size distribution from characteristic to scale-free, aftershock productivity, Bath gap, rupture patterns, and with them regularity | |
| clustering (branching ratio) | Bath gap and Omori p; above about 0.3 also short intervals between large events | mean intervals |
| admission (arrest scale and sharpness) | how far ruptures run into partly loaded fault: co-rupture between sites | |
| quenched barriers (strength heterogeneity) | moderately less regular sites, fewer full-length ruptures | site-to-site differences in regularity were not produced |
| roughness, stress-drop scatter, geometry, resolution | little, once slip is capped by the deficit | |

The two main axes of complexity are therefore timing regularity (the
nucleation clock) and the shape of the size distribution (the size
rule), with clustering as a third that couples aftershock statistics to
regularity. Only the first can be changed without changing the
seismicity statistics.

## 4. The Alpine Fault example

### 4.1 What fits

* The slip budget makes mean recurrence intervals easy to match once
  stress drop and initiation rate are calibrated; regularity (CV) is the
  discriminating observation.
* Judged against the observations, three model classes were rejected:
  stress-drop slip with scale-free sizes (site CV not below about 0.4),
  slip equal to the deficit (no Gutenberg-Richter), and strong clustering
  (branching ratio above about 0.3: irregular sites). Deficit-capped slip
  with the hybrid size rule fits; geometry (planar, variable dip,
  segmented transition) and resolution (1 and 2 km patches) make little
  difference.
* Paleoseismic records are short (4 to 14 events). Records of that
  length drawn from the models are more regular than long runs (window
  CVs about 0.07 lower), so model CVs must be compared with records of
  the observed length, not with long-run averages.

### 4.2 Which cases match

Fraction of synthetic records of the observed length (last 13, 14 and 4
events before a common end time) that match the observed RI and CV at
all three sites at once (`record_match.py`):

| case | match, all sites, RI and CV | basis |
|---|---|---|
| reference (nucleation exponent 12.3, area threshold 488 km^2) | 16 % | set A, 50 realizations |
| medium clock (exponent 6) | 8.6 % | set A, 50 realizations |
| soft clock (exponent 3) | 1.5 % | set A, 50 realizations |
| characteristic (area threshold 100 km^2) | 16 % | 4 realizations, provisional |
| scale-free (area threshold 2000 km^2) | 2 % | 4 realizations, provisional |
| barriers (strength ln-sigma 0.5, 1.0) | not yet computed | set C |

* Regular, characteristic models (reference, characteristic) match the
  records best, in about one record set in six; their best record sets
  agree with the observations in every number (reference: south 238 yr /
  CV 0.35, central 254 / 0.24, north 443 / 0.26; observed 268 / 0.40,
  249 / 0.23, 413 / 0.25).
* Irregular models match the irregular south record but lose the
  regular central record; regular models do the reverse. In every model
  the south and central sites share more than 90 percent of their large
  ruptures, so the observed contrast between the two records is rare in
  all of them and is the main mismatch.
* The Bath gap of the best-matching models is large (M >= 6 mainshocks:
  1.33 in the reference case, 2.16 in the characteristic case; about 2
  for through-going M8 ruptures in all cases): matching the regular
  records favours weak triggering.

## 5. Predictions for the model comparison

The predictions are the behaviour of all stochastic models consistent
with the constraints, across model classes. `predictions.py` pools the
cases of the realization sets and weights each realization by its
likelihood under the constraints, the fraction of its records of the
observed length that match the observed RI and CV at all three sites.
Model classes then contribute in proportion to how well they match, and
the weighted 5 to 95 percent ranges are the predicted ranges.
`run_realization_sets.sh` writes them as its last step
(`predictions.md`, `predictions.png`).

Predicted quantities: RI, CV and displacement at the sites; mean
interval and CV of whole-fault M >= 7.5 events; maximum magnitude; b
and the excess of large events over Gutenberg-Richter; Bath gap; the
fraction of large ruptures that span all three sites; co-rupture of the
central and south sites. Further statistics on the figure pages
(hypocentre depths, Omori decay, moment release, rupture lengths,
slip-predictability) can be added in the same way.

Provisional example (4 cases x 4 realizations, effective sample size
about 7; characteristic case 64 percent of the weight, reference 24):

| quantity | constraint-weighted median [5-95 %] | all realizations |
|---|---|---|
| central RI [yr] | 226 [215-282] | 275 [219-305] |
| central CV (long run) | 0.22 [0.16-0.40] | 0.26 [0.17-0.42] |
| whole-fault M >= 7.5 interval [yr] | 214 [178-286] | 235 [174-314] |
| whole-fault M >= 7.5 CV | 0.27 [0.17-0.72] | 0.31 [0.18-0.76] |
| maximum magnitude | 8.39 [8.34-8.47] | 8.43 [8.34-8.47] |
| excess of M >= 7.5 over b = 1 | 8.5 [0.7-10.4] | 1.6 [0.7-10.1] |
| Bath gap (M >= 6) | 2.2 [1.1-2.3] | 1.35 [1.1-2.3] |
| fraction of M >= 7.5 rupturing all three sites | 0.50 [0.24-0.74] | 0.59 [0.22-0.78] |

Weighting by the constraints shifts the ensemble toward regular,
characteristic behaviour: fewer moderate-to-large events relative to
Gutenberg-Richter, large Bath gaps, and through-going ruptures in about
half of the large events. The full sets will narrow these ranges and
show how much of the weight the barrier cases carry.

## 6. What a physics-based simulator would add

Beyond the stochastic baseline, a simulator demonstrates something new
if it matches the paleoseismic records together with any of:

* the contrast between a regular central and an irregular south record
  (not produced by any stochastic class tested);
* regular sites with near-Gutenberg-Richter statistics and standard
  aftershock sequences (in the stochastic models, regular sites come with
  a characteristic size distribution and weak triggering);
* a mechanism for the controls that are imposed here: the nucleation
  clock, the size rule and the deficit cap on slip.

Shared predictions of the constraint-weighted stochastic ensemble that a
simulator can be checked against: slip-predictable site records
(displacement correlated with the preceding interval), hypocentres near
the base of the seismogenic zone, and a depletion of M6.5 to 7.5 events
relative to Gutenberg-Richter.
