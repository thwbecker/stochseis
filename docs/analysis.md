# Analysis tools

`stochseis.analysis` computes statistics from a `Catalog`;
`stochseis.plots` draws them; `python -m stochseis.cli analyze` runs
both and prints a short table. `summary_metrics(cat)` returns one flat
dictionary and is what sweeps collect into `sweep_summary.csv`.

Two views of clustering are available throughout: the *observer* view
selects mainshocks and aftershocks with space-time windows (default:
1 yr, max(10 km, rupture length) from the centroid, no larger event in
the window), and the *genealogy* view uses the stored ETAS parent links.
Differences between the two measure the effect of selection.

## Frequency-magnitude distribution

* `fmd`, `b_value_mle`: maximum-likelihood b (Aki, 1965) for continuous
  magnitudes; the Utsu half-bin correction only for binned data (dm > 0),
  and the doubly truncated likelihood (Page, 1968) when an upper bound is
  given. Applying the half-bin correction to continuous synthetic
  magnitudes biases b low by about 10 percent; earlier b values in this
  project carried that bias (running b) or a compensating truncation
  bias (regime b).
* `b_value_regimes`: b in the sub-patch range [mag_min, M_patch), in the
  fault range below width saturation, and above saturation.
  `saturation_magnitude` is the magnitude above which the median aspect
  ratio exceeds 2.

## Source scaling

* `scaling_fits`: slopes of log M0 versus log L for L < 0.8 W and
  L > 3 W (3 and 2 are the observed values below and above saturation;
  a pure W-model gives 3 and 1), the median mean-slip to length ratio
  (Blaser et al. 2010 data give about 1.3e-5), the surface-rupturing
  fraction by magnitude, the multi-segment fraction and the maximum
  magnitude.
* The summary figure shows M0-L, L/W-M0 (as in Tectonic Geodynamics
  Fig. 3.43) and mean slip versus length.

## Clustering

* `omori_stack`, `fit_omori`: aftershock delays stacked over mainshocks
  of M >= mmain, maximum-likelihood Omori (K, c, p) over [0, 1 yr].
* `bath_gaps`, `bath_gaps_genealogy`: M_main minus largest aftershock,
  observer and genealogy views; the summary records the mean, median
  and the number of mainshocks without any aftershock above completeness.
* `foreshock_fraction`: fraction of window-selected mainshocks preceded
  by an event of M >= mag_min within 10 days and 10 km of the mainshock
  epicentre (summary also reports 3 days and 5 km), together with the
  chance level from random times at the same epicentres. Epicentres are
  nucleation patches. An earlier definition used centroid distances and
  a rupture-length radius; its chance level approached one for large
  ruptures.
* `realised_branching_ratio` (metadata): fraction of triggered events.

## Recurrence and budget

* `interevent_cv`: coefficient of variation of interevent times of
  M >= 6.5, declustered (mainshocks only) and raw. 1 is Poisson; the
  deficit control lowers it.
* `point_recurrence`, `site_recurrence_cv`: event times at which a
  chosen patch slipped more than a threshold in a stored event, i.e. a
  synthetic paleoseismic site, and the CV averaged over nine sites along
  strike. On a long fault the whole-fault interevent CV of M 6 to 6.5
  events tends to 1 by superposition of independent segments; site
  recurrence is the quantity comparable to paleoseismic CVs.
* `moment_closure`, `cumulative_moment`, `slip_rate_recovered`: released
  versus loaded moment, cumulative curves, and the ratio of recovered to
  prescribed slip rate per patch (`deficit_figure`).
* Metadata records the fault-mean deficit and the mean, standard
  deviation and negative fraction of the backslip stress over time as a
  measure of stress heterogeneity.

## Arrest statistics

Each fault event stores its arrest code (size hazard, state blocked,
fault boundary) and the state ratio at arrest, so the summary reports the
fraction of events, and of M >= 6 events, that were stopped by the
budget rather than by the scale-free hazard.

## Sweeps

`examples/single_fault/experiments/sweep_arrest_stressdrop.json` is a small template and
`examples/single_fault/experiments/sweep_regime.json` the planned regime sweep (ROADMAP step 1).
`--no-catalogs` keeps only the summary table; `examples/single_fault/experiments/aggregate_sweep.py`
averages metrics over seeds. `plots.sweep_figure`
plots any summary metric against one parameter, grouped by another.
