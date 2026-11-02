# Parameters

All parameters live in `stochseis.params.Params`; a run is reproducible
from its JSON file and seed. Units in the interface: km, mm/yr, MPa,
years, days (Omori c). Section numbers refer to `docs/design.md`.

## Fault geometry (planar builder)

| name | default | meaning |
|---|---|---|
| fault_length_km, fault_width_km | 200, 15 | rupture plane dimensions |
| patch_km | 1.0 | patch size; sets the single-patch magnitude (about M4.2 at 3 MPa) |
| dip_deg, top_depth_km, rake_deg | 90, 0, 180 | top_depth 0 means the free-surface image is used; rake selects strike or dip slip for the stiffness |
| slip_rate_mm_yr | 20 | long-term loading; multiplied by patch area and rigidity it fixes the moment budget |
| rigidity_pa | 3e10 | |
| segment_gaps_km | [] | list of [position, gap]: cuts and step-overs |
| jump_max_km | 5 | step-overs up to this width can be crossed, with probability 1 - gap/jump_max (UCERF3 uses a 5 km hard cutoff) |

## Slip amplitude and roughness (A2, A3)

| name | default | meaning |
|---|---|---|
| stress_drop_mpa | 3 | median macroscopic stress drop of the crack solution |
| stress_drop_ln_sigma | 0.25 | lognormal scatter (0.25 is a factor 1.28 at one sigma; 0.5 gave site-recurrence CVs near 0.8); 0 for constant stress drop |
| deep_slip_km | 0 | crack solutions of ruptures reaching the base of the locked width extend this far below it (E18) |
| hurst | 0.75 | von Karman H (Mai and Beroza, 2002) |
| slip_cv | 0.75 | coefficient of variation of the lognormal roughness pattern |
| corr_length_scatter | 0 | ln-scatter on the Mai and Beroza correlation lengths |
| condition_weight | 1 | exponent of the state weighting of slip inside the rupture; 0 switches it off |
| cap_by_deficit | False | mean slip of an event cannot exceed the mean positive deficit on its rupture (roughness still overshoots locally); prevents forced release in relaxed zones |

## Growth, arrest, kinematics (A4, A5, A6)

| name | default | meaning |
|---|---|---|
| b_value | 1.0 | size hazard per unit ln area, equal to the GR b below width saturation |
| stop_measure | "moment" | "moment" keeps b continuous through width saturation, "area" raises it by a factor 1.5 |
| size_hazard | True | False: no scale-free size hazard, ruptures stop only when no neighbour is admitted (state-only termination; gives characteristic size distributions, see ROADMAP) |
| segments_file, segment_include, segment_max_dist_km, segment_taper_km | "", [], 10, 10 | segment table with end points, rake and slip rate; projected onto mesh traces (docs/alpine.md) |
| slip_mode, slip_fraction | "stress_drop", 1 | "deficit": each patch releases slip_fraction of its positive deficit; "local_cap": stress-drop crack slip capped patch by patch at slip_fraction of the deficit (docs/alpine.md) |
| hazard_max_area_km2 | -1 | hybrid sizes: scale-free size hazard only below this area, state-only termination above |
| admission_exponent | 2 | m in the admission probability 1 - exp(-(rho/rho_c)^m) |
| state_smoothing_km | 0 | Gaussian smoothing of the deficit seen by admission and nucleation |
| segment_boundary_strength, segment_boundary_width_km | 1, 4 | strength factor of patches near segment boundaries, dividing the admission ratio (rupture-passage control) |
| sites_file, site_thresholds_m | "", [1, 2, 4] | paleoseismic sites and detection thresholds for the summary statistics |
| strength_ln_sigma, strength_hurst, strength_corr_km | 0, 0.25, -1 | quenched strength heterogeneity dividing the admission ratio (E17); off by default |
| state_measure | "deficit" | state used for admission and slip conditioning; "stress" uses the backslip stress (noisy, see design.md) |
| arrest_scale | 0.25 | rho_c of the admission rule; 0 disables the budget control of ruptures. Ruptures percolate where deficit exceeds about 0.73 rho_c times the implied slip; 0.5 confines the budget release to M < 7 on the test fault |
| track_stress | True | compute the backslip stress field for heterogeneity diagnostics |
| vs_km_s, rupture_speed_frac | 3.5, 0.8 | shear-wave speed and sub-Rayleigh rupture speed |
| supershear_prob, supershear_speed_frac | 0.1, 1.5 | fraction and speed of supershear ruptures |

## Nucleation and clustering (A7)

| name | default | meaning |
|---|---|---|
| mag_min | 3.0 | completeness of the synthetic catalog |
| etas | True | False gives stationary Poisson nucleation |
| branching_ratio | 0.5 | nominal mean number of direct offspring per event (realised about 0.41 with the defaults, E22); controls both Bath's gap and site recurrence CV |
| alpha | -1 (= b) | productivity exponent |
| omori_p, omori_c_days | 1.1, 0.01 | temporal kernel |
| spatial_d_km, spatial_q | 1.0, 1.5 | distance kernel r (r^2 + d^2)^-q |
| nucleation_rate_per_km2_yr | -1 (auto) | fault-event nucleation rate; negative means budget calibration |
| budget_scale | 1.0 | multiplier on the calibrated rate |
| rate_feedback | True | background rate proportional to fault-mean deficit (uniform nucleation only; ignored with deficit weighting) |
| nucleation_exponent | 1.0 | gamma in the deficit weight (deficit / delta_ref)^gamma; 3 makes large events markedly more periodic (ROADMAP) |
| nucleation_weighting | "deficit" | "deficit": local nucleation intensity proportional to positive deficit, offspring placed preferentially on loaded patches; makes rate_feedback and cap_by_deficit unnecessary (design.md section 2.6) |

## Run control

| name | default | meaning |
|---|---|---|
| duration_yr | 3000 | |
| seed | 1 | numpy Generator seed; runs are deterministic |
| store_slip_min_mag | 6.0 | slip maps stored for events at or above this |
| pilot_events | 2000 | pilot ruptures; used for the calibration only with the area hazard, otherwise for diagnostics (the moment hazard is calibrated analytically) |
| initial_deficit | "random" | smooth random field in [0, fault-wide slip]; "zero" |
| verbose | False | |

## How the defaults were chosen

After the calibration fix (design.md section 2.5), a reduced regime sweep
(ROADMAP step 1, three seeds of 3000 yr per setting) showed that
arrest_scale, stress-drop scatter and the nucleation mode change site
recurrence CVs by at most 0.1, while the branching ratio controls both
the CV and Bath's gap. branching_ratio = 0.5 gives Bath 1.25 and a site
CV of 0.62 for slip >= 1 m; 0.8 gives 1.05 and 0.76. Deficit-weighted
nucleation is the default because it closes the budget without the
global feedback and the slip cap and produces no events below mag_min.
stress_drop_ln_sigma = 0.25 and arrest_scale = 0.25 are kept; they are
not well constrained by these statistics.

## What is effectively free

For single-fault experiments the quantities worth sweeping are
stress_drop_mpa and its scatter, b_value, arrest_scale, branching_ratio
and the two spatial-kernel numbers, plus stop_measure and state_measure
as discrete choices. Everything else is either geometry, a convention
(patch size, completeness) or taken from the literature.
