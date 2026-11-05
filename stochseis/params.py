"""Model parameters.

All tunable quantities of a simulation live in one flat dataclass so that a
run is fully described by one JSON file.  Parameters are grouped by comment
blocks below; see docs/parameters.md for a discussion of which ones are
genuinely free and which are fixed from published constraints.

Units: km for lengths in the interface, mm/yr for slip rates, MPa for stress
drops, years for times, days for the Omori c constant.  Internally the code
uses SI units (m, Pa, yr).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields


@dataclass
class Params:
    # ------------------------------------------------------------------
    # Fault geometry (planar builder; see geometry.py)
    # ------------------------------------------------------------------
    fault_length_km: float = 200.0
    fault_width_km: float = 15.0
    patch_km: float = 1.0
    dip_deg: float = 90.0
    top_depth_km: float = 0.0
    rake_deg: float = 180.0
    slip_rate_mm_yr: float = 20.0
    rigidity_pa: float = 3.0e10
    # geometry from a file (triangulated surface, STL); empty = planar builder
    geometry_file: str = ""
    geometry_units: str = "auto"          # "m", "km" or "auto"
    # piecewise-linear slip rate along strike, [[s_km, rate_mm_yr], ...]; empty = slip_rate_mm_yr everywhere
    slip_rate_profile: list = field(default_factory=list)
    section_breaks_km: list = field(default_factory=list)   # segment labels along strike (mesh geometries)
    segments_file: str = ""               # segment table with end points, rake, slip rate (overrides the profile)
    segment_taper_km: float = 10.0        # linear blending of rates across segment boundaries
    segment_max_dist_km: float = 10.0     # segments farther than this from the trace are ignored
    segment_include: list = field(default_factory=list)     # name prefixes of segments to use; empty = all
    segment_boundary_strength: float = 1.0  # strength factor of patches within segment_boundary_width_km of a segment boundary
    segment_boundary_width_km: float = 4.0
    # paleoseismic sites for summary statistics (CSV, see data/nz_alpine/paleo_sites.csv); empty = none
    sites_file: str = ""
    site_thresholds_m: list = field(default_factory=lambda: [1.0, 2.0, 4.0])
    site_merge_yr: float = 0.0            # merge recorded site events closer than this (paleoseismic dating resolution)
    # list of [position_km, gap_km]: the fault is cut at position_km along
    # strike and the two segments are separated by a step-over of gap_km
    segment_gaps_km: list = field(default_factory=list)
    # connections with gap <= jump_max_km are open with probability
    # 1 - gap/jump_max_km for each event (UCERF3 uses a 5 km hard cutoff)
    jump_max_km: float = 5.0

    # ------------------------------------------------------------------
    # Slip amplitude and roughness
    # ------------------------------------------------------------------
    stress_drop_mpa: float = 3.0          # median of lognormal stress drop
    stress_drop_ln_sigma: float = 0.25    # scatter in ln(stress drop); 0 = constant
    dip_correction: bool = False          # L1: scale crack slip by the half-space / spectral ratio for dip, rake, burial
    deep_slip_km: float = 0.0             # crack solutions extend this far below the locked width where ruptures reach its base
    hurst: float = 0.75                   # von Karman Hurst exponent (Mai & Beroza 2002)
    slip_cv: float = 0.75                 # coefficient of variation of the lognormal slip pattern
    corr_length_scatter: float = 0.0      # ln-scatter on correlation lengths; 0 = Mai & Beroza mean
    condition_weight: float = 1.0         # exponent on the state measure used to condition the slip pattern; 0 = none
    slip_mode: str = "stress_drop"        # "stress_drop": crack slip; "deficit": slip = slip_fraction x local deficit; "local_cap": crack slip capped per patch at slip_fraction x local deficit
    slip_fraction: float = 1.0            # fraction of the local positive deficit released in slip_mode "deficit"
    cap_by_deficit: bool = False          # mean slip of an event cannot exceed the mean positive deficit on its rupture

    # ------------------------------------------------------------------
    # Rupture growth, arrest and kinematics
    # ------------------------------------------------------------------
    b_value: float = 1.0                  # GR b of the size hazard (2b/3 per ln M0 or b per ln A, see stop_measure)
    hazard_max_area_km2: float = -1.0     # size hazard acts only below this rupture area; above it ruptures stop by state only (negative = always)
    size_hazard: bool = True              # False: no scale-free size hazard; ruptures stop only by state (admission)
    stop_measure: str = "moment"          # "moment": hazard 2b/3 per ln M0 (b continuous); "area": hazard b per ln A
    strength_ln_sigma: float = 0.0        # quenched strength heterogeneity: ln-std of a lognormal von Karman field (0 = none)
    strength_hurst: float = 0.25          # Hurst exponent of the strength field
    strength_corr_km: float = -1.0        # correlation length of the strength field; negative = fault length
    # persistent hierarchical patch structure (stochseis/hierarchy.py); hier_rmin_km = 0 switches it off
    hier_rmin_km: float = 0.0             # smallest patch radius
    hier_rmax_km: float = 40.0            # largest patch radius
    hier_ratio: float = 2.0               # radius ratio between successive levels (background: ratio x largest)
    hier_fill: float = 0.3                # area fraction covered by each level
    hier_nest: float = 0.8                # probability that a patch lies inside a patch of the next larger level
    hier_barrier_exponent: float = 1.0    # eta: admission ratio divided by max(1, r_patch / R_rupture)^eta
    hier_nucleation_exponent: float = 0.0  # kappa: nucleation weight multiplied by (r_min / r_patch)^kappa
    state_measure: str = "deficit"        # "deficit" or "stress" (tau = K deficit) drives arrest and conditioning
    track_stress: bool = True             # compute the backslip stress field for heterogeneity diagnostics
    state_smoothing_km: float = 0.0       # Gaussian smoothing (sd) of the deficit seen by admission and nucleation; mimics elastic coupling
    admission_exponent: float = 2.0       # m in q = 1 - exp(-(rho/rho_c)^m); larger = sharper deficit threshold
    arrest_scale: float = 0.25            # rho_c: patch admitted with prob 1 - exp(-(rho/rho_c)^2); 0 disables state-dependent arrest
    vs_km_s: float = 3.5                  # shear wave speed used for rupture kinematics
    rupture_speed_frac: float = 0.8       # sub-Rayleigh rupture speed as a fraction of vs
    supershear_prob: float = 0.1          # probability that a rupture runs supershear
    supershear_speed_frac: float = 1.5    # supershear rupture speed as a fraction of vs (sqrt(2) vs .. vp)

    # ------------------------------------------------------------------
    # Nucleation and ETAS clustering
    # ------------------------------------------------------------------
    mag_min: float = 3.0                  # smallest simulated magnitude (completeness of the synthetic catalog)
    etas: bool = True                     # if False, nucleation is a stationary Poisson process
    branching_ratio: float = 0.5          # mean number of direct offspring per event
    alpha: float = -1.0                   # productivity exponent; negative means alpha = b_value
    omori_p: float = 1.1
    omori_c_days: float = 0.01
    spatial_d_km: float = 1.0             # scale of the (r^2 + d^2)^-q offspring distance kernel
    spatial_q: float = 1.5
    nucleation_rate_per_km2_yr: float = -1.0  # fault-event nucleation rate; negative = set from moment budget
    budget_scale: float = 1.0             # multiplies the budget-derived nucleation rate
    rate_feedback: bool = True            # scale background nucleation rate by fault-mean deficit / reference deficit
    nucleation_exponent: float = 1.0      # gamma: deficit-weighted intensity ~ (deficit/delta_ref)^gamma; >1 is more renewal-like
    nucleation_weighting: str = "deficit" # "uniform" or "deficit": local nucleation intensity proportional to positive deficit

    # ------------------------------------------------------------------
    # Run control
    # ------------------------------------------------------------------
    duration_yr: float = 3000.0
    seed: int = 1
    store_slip_min_mag: float = 6.0       # slip distributions are stored for events at or above this magnitude
    pilot_events: int = 2000              # growths used to calibrate the nucleation rate from the budget
    initial_deficit: str = "random"       # "random" (uniform in [0, one large-event slip]) or "zero"
    verbose: bool = False

    # ------------------------------------------------------------------
    def to_json(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(asdict(self), f, indent=2)

    @classmethod
    def from_json(cls, path: str) -> "Params":
        with open(path) as f:
            return cls.from_dict(json.load(f))

    @classmethod
    def from_dict(cls, d: dict) -> "Params":
        names = {f.name for f in fields(cls)}
        unknown = set(d) - names
        if unknown:
            raise KeyError(f"unknown parameter(s): {sorted(unknown)}")
        return cls(**d)

    def replace(self, **kw) -> "Params":
        d = asdict(self)
        for k, v in kw.items():
            if k not in d:
                raise KeyError(f"unknown parameter: {k}")
            d[k] = v
        return Params(**d)

    def as_dict(self) -> dict:
        return asdict(self)

    @property
    def alpha_eff(self) -> float:
        return self.b_value if self.alpha < 0 else self.alpha
