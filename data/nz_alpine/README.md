# Alpine Fault intercomparison data

Inputs for the New Zealand Alpine Fault simulator intercomparison
(August 2025 exercise), as provided by the organisers.

* `paleo_sites.csv`: single-event displacements (SED) and recurrence
  intervals (RI) at three paleoseismic sites, transcribed from
  "Paleoseismology derived Alpine Fault SED and RI.xlsx" (J. Howarth).
  Sources given in the spreadsheet: Langridge et al. (in review), De
  Pascale et al. (2014), Berryman et al. (2012) for SED; Howarth et al.
  (2021) and Langridge et al. (in review) for RI. Site names here are
  positional labels (south, central, north), not the published site
  names. Columns: n_sed and n_ri are the numbers of events behind the
  SED and RI values; sed_min and sed_max include measurement
  uncertainty. ri_sigma_yr is interpreted in stochseis as the standard
  deviation of individual recurrence intervals (so ri_sigma / ri is the
  CV: 0.40, 0.23, 0.25 south to north); this interpretation needs to be
  confirmed with the data providers.
* Fault geometries (STL, coordinates in km NZTM, depth negative) are not
  stored in the repository because of their size (9 to 39 MB). Expected
  file names: `alpine_planar_1km_26Jan.stl`, `alpine_planar05km.stl`,
  `variable_dip_1km_26Jan.stl`, `variable_dip_05km_improved3.stl`,
  `Alpine_full_segmented_transition_05km_d20.stl`.
  Point `geometry_file` at them.
* Segment slip rates (`alpine_hope_faults_segment_slip_rates`,
  `alpine_segments.pdf`) have not been received yet; runs so far use a
  placeholder (see docs/alpine.md).

* `rsqsim_summary.csv`: site statistics of the RSQSim simulations of the
  same exercise (planar and varying-dip geometries, 50,000 yr), transcribed
  from the slides "RSqsim_results.pdf" provided by the organisers (sites 1
  to 3 = south, central, north). An example of a reference table for
  `examples/alpine/site_analysis.py --reference`; other models' results
  can be added in the same format.
