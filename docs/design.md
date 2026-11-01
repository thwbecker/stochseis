# stochseis: model description

This document states the assumptions of the generator in order of
importance, then the rules that follow from them, then what is fixed
from published work and what remains free. The code is meant as a
baseline that physics-based earthquake cycle models should be measured
against, so the aim is the smallest set of rules that still reproduces
the first-order statistics of seismicity, with everything else left to
emerge.

## 1. Core assumptions

A1. **Elasticity.** The fault is a planar surface in a homogeneous,
isotropic, linear elastic half space (Poisson solid). Static stress
changes from slip follow from the exact full-space stiffness of a planar
dislocation distribution (Andrews, 1980); the free surface is included
by an image for faults that reach it. There is no viscoelasticity and no
off-fault inelasticity.

A2. **Constant loading, constant stress drop.** Each patch is loaded at
a prescribed long-term slip rate. The size of a rupture given its
outline is set by a macroscopic static stress drop drawn from a
lognormal distribution with a median of a few MPa and a scatter of about
a factor of two to three, independent of magnitude (Kanamori and
Anderson, 1975; Allmann and Shearer, 2009; the compilation in Becker,
Tectonic Geodynamics, Fig. 3.43). The slip of a rupture is therefore the
static crack solution for uniform stress drop on its outline. Width
saturation at the brittle-ductile transition and the factor of about two
between buried and surface-breaking ruptures follow from this, not from
scaling laws.

A3. **Self-similar roughness.** Within the outline, slip is the crack
solution multiplied by a lognormal von Karman random field with Hurst
exponent H = 0.75 and correlation lengths proportional to rupture length
and width (Mai and Beroza, 2002). Roughness redistributes slip but does
not change the moment.

A4. **Causality and kinematics.** Nothing propagates faster than the
P-wave speed. Rupture fronts travel at about 0.8 of the shear-wave
speed, with a minority of events (default 10 percent) supershear at
about 1.5 times the shear-wave speed (Bao et al., 2022, found about 14
percent of large shallow strike-slip events supershear). Each patch has
an onset time equal to its path distance from the nucleation point
divided by the rupture speed; triggering of other events is counted from
that onset, never from before it. Stage 0 uses no waveforms; the
kinematics only order events in time and define rupture durations.

A5. **Scale-free size hazard.** Rupture arrest, in the absence of any
information about the fault state, is a constant hazard per unit
logarithm of rupture size. This is the minimal rule that produces a
power-law size distribution; its single parameter is the Gutenberg-
Richter b-value (Gutenberg and Richter, 1944). Whether the hazard is per
unit log-area or per unit log-moment decides whether b changes at width
saturation (Pacheco et al., 1992) or not; both are available.

A6. **Slip budget as the only state.** The fault state is the slip
deficit, deficit = (slip rate x time) - released slip, together with the
loading stress it implies through A1, tau = K deficit (backslip loading;
Savage, 1983), and a maturity variable (cumulative slip, recorded as a
hook, unused by the Stage 0 rules). Ruptures are conditioned on this
state in two ways only: a patch whose state cannot support the slip
about to be imposed is unlikely to be admitted into the rupture, and
within the rupture slip concentrates where the state is high. Overshoot
is allowed; the budget closes in the long-term mean.

A7. **Clustering is epidemic.** Nucleations occur at all times as a
stationary Poisson background plus Omori-law triggering by every event
(ETAS; Ogata, 1988). The productivity exponent equals b, so that the
Bath gap is an outcome rather than a constraint (Helmstetter and
Sornette, 2003). The background rate is not free: it is set so that the
expected moment release equals the tectonic moment rate, with a weak
global feedback on the fault-average deficit that keeps the balance
stable.

Everything not listed above (magnitude distribution of fault events and
its curvature, Bath's law, recurrence statistics, foreshock fractions,
moment-area scaling, stress heterogeneity, aftershock productivity of
individual sequences) is emergent and is what the analysis tools
measure.

## 2. Rules

### 2.1 Geometry

The fault is a graph of patches with position, area, slip rate,
rigidity, section id, local along-strike and down-dip coordinates, and
a weighted adjacency (8-connected within a section; step-overs are
connections with their own distance). Stage 0 builds planar, optionally
segmented, faults. The simulator never sees anything but this graph, so
UCERF3-style fault sections or a stochastic fault-network generator only
have to produce the same object.

### 2.2 Nucleation

Background nucleations are uniform over the fault area in space and
Poisson in time with rate nu_bg = nu_tot (1 - n), where n is the ETAS
branching ratio. Every event with magnitude M triggers a Poisson number
of offspring with mean K0 10^(alpha (M - M_min)), delayed by
(t + c)^-p from the onset time of a slip-weighted random patch of the
parent, and displaced from it by a 2-D power-law kernel
(r^2 + d^2)^-q. K0 follows from the target branching ratio and the
expected productivity over the realised magnitude distribution.

Every nucleation is a rupture governed by the same size hazard from the
area A_min at which a rupture has magnitude M_min for its drawn stress
drop. With probability (A_min / A_patch)^b it reaches the patch size and
is grown on the graph; otherwise it is a sub-patch event whose area is
drawn from the hazard within [A_min, A_patch), whose moment follows the
one-patch solution scaled as A^(3/2), and whose equivalent slip is
debited on the nucleation patch. Ruptures stopped by the hazard at k
patches receive a continuous size within the last increment. The
frequency-magnitude distribution is therefore continuous across the
patch size by construction. Because a one-cell boxcar carries 1.7 times
the moment of a resolved crack of the same area, the moment hazard is
evaluated on the realised moment versus patch-count relation measured
in the pilot ensemble rather than on the analytic crack law; without
this correction the rate of grown ruptures sits 28 percent below the
continuation of the sub-patch GR line.

### 2.3 Growth and arrest

The rupture floods outward from the nucleation patch: candidates are the
graph neighbours of already admitted patches, visited in order of path
distance from the nucleation point, so only contiguous sets can rupture.
Each candidate is admitted with probability q = 1 - exp(-(rho / rho_c)^2),
where rho is the positive state measure of the patch (deficit relative
to the mean slip implied by the current rupture size and stress drop, or
loading stress relative to the stress drop), and is otherwise blocked.
A region in which q falls below the site-percolation threshold of the
8-connected lattice (about 0.41, i.e. rho below about 0.73 rho_c) cannot
be crossed; the rupture grows around depleted patches and stops when its
frontier is exhausted. Blocked patches that end up fully enclosed are
added, since they are forced to slip. Independently, the admitted
rupture stops with probability 1 - (A_k / A_{k+1})^b at each increment
of its area (or the moment-based equivalent).

In addition, the mean slip of any event is capped at the mean positive
deficit on its rupture (``cap_by_deficit``); roughness and conditioning
still redistribute slip, so local overshoot remains. Without the cap,
events nucleating in relaxed zones release their full stress-drop slip
regardless of the budget, which drives the fault to a depleted
equilibrium once the nucleation rate is high.

With rho_c = 0 the admission rule is off and the size distribution on an
unbounded fault is N(>A) proportional to A^-b exactly. With
``size_hazard = False`` the size hazard is off instead and ruptures stop
only when no neighbouring patch is admitted; this state-only mode
produces characteristic rather than GR size distributions (ROADMAP,
tests of 2 November 2026). rho_c therefore
controls how strongly the slip budget shapes ruptures: small values give
a pure epidemic GR model with ruptures independent of the budget, large
values give budget-controlled, quasi-periodic behaviour. Two earlier
formulations were tried and rejected: a front-averaged arrest hazard
(smooths over barriers, so ruptures crossed them) and admission without
the contiguity requirement (ruptures percolated through semi-depleted
fault as sparse, hole-ridden patterns).

### 2.4 Slip

For the admitted outline the static crack problem K D = dsig is solved
with a spectral stiffness operator and preconditioned conjugate
gradients (a few iterations; solutions are cached by outline shape).
The free surface is included by mirroring when the rupture reaches the
top of a fault with zero top depth. The crack solution is multiplied by
a lognormal von Karman pattern of mean one and by the normalised
positive state measure, and the product is rescaled to the moment of the
crack solution. The realised moment-weighted stress drop of the rough
field (Noda et al., 2013) is larger than the macroscopic value and is
recorded as a diagnostic.

### 2.5 Rate calibration

The total nucleation rate is chosen so that the expected moment release
per nucleation, absent any state-dependent arrest, times the rate equals
the tectonic moment rate. For the moment-based size hazard the expected
moment is computed analytically: sub-patch sizes follow a truncated
power law in moment; a grown rupture survives from k to k+1 patches with
probability (M0_k / M0_{k+1})^(2b/3), using the realised moment versus
patch count from the pilot, stops at a continuous moment inside that
increment, and stops at the fault boundary otherwise; the stress-drop
distribution is integrated by quadrature. The same calculation gives the
expected ETAS productivity that sets K0. Before 1 November 2026 both
expectations were Monte Carlo means over 10^3 pilot ruptures. Because
the moment distribution is heavy-tailed (survival M0^(-2/3)), such means
are dominated by events that a pilot of that size rarely contains; they
were typically several times, in one case ten times, too small, so the
nucleation rate was too high and the feedback below held the fault at a
low mean deficit (0.5 m instead of about 3 m). The resulting shortage of
M >= 6.5 events relative to GR, attributed earlier to the arrest rule,
was largely this calibration error. The pilot is still used for the
area-based hazard and for diagnostics.

### 2.6 State update and feedback

After each event the released slip is added per patch and the backslip
stress is updated (one FFT, buffered for small events). Two ways of
letting the state feed back on nucleation are available.

Uniform nucleation (``nucleation_weighting = "uniform"``): the
background rate is multiplied by g(t) = fault-mean deficit / delta_ref
with delta_ref = half the slip of a fault-wide event; triggered events
are not modulated, so aftershock sequences continue when the background
is quiet. Without this feedback the deficit drifts, because nucleations
in the shadow of a recent rupture release less moment than the
calibration assumes. Event mean slip is capped by the available deficit
(``cap_by_deficit``), and nucleations on relaxed patches therefore give
events below mag_min (3 to 8 percent with the corrected calibration).

Deficit-weighted nucleation (default, ``nucleation_weighting = "deficit"``): the
local background intensity is proportional to w_i = (max(deficit_i, 0) /
delta_ref)^gamma (gamma = ``nucleation_exponent``, 1 by default),
implemented by thinning on the fault-mean weight against an adaptive
bound and choosing the patch with probability proportional to w_i. The fault mean of w replaces g(t). Triggered events keep
their number, so ETAS productivity is unchanged, but each offspring is
placed at one of eight candidate locations from the spatial kernel,
chosen with probability proportional to w. Neither the global feedback
nor the slip cap is needed, and no events fall below mag_min. A first
version thinned triggered events at their location and rescaled
productivity by the global acceptance; because acceptance inside a fresh
rupture is far below average, that halved the direct aftershocks of
M >= 6 events and raised Bath's gap to about 1.3, and it was replaced.

## 3. Parameters and their provenance

Fixed from published work (treated as data, not tuned):

| quantity | value | source |
|---|---|---|
| Hurst exponent H | 0.75 | Mai and Beroza (2002), median of 44 slip models |
| correlation lengths | a_x = 2 + L/3 km, a_z = 1 + W/3 km | Mai and Beroza (2002), eq. 4 |
| slip CV | 0.75 | Graves and Pitarka (2010) use 0.85; LeVeque et al. (2016) 0.75 |
| jump distance | 5 km | UCERF3 plausibility filter (Field et al., 2014) |
| rupture speed | 0.8 v_s | typical sub-Rayleigh values, e.g. Heaton (1990) |
| supershear fraction and speed | 0.1, 1.5 v_s | Bao et al. (2022) |
| Omori c, p | 0.01 d, 1.1 | typical California ETAS fits (Ogata, 1988; Helmstetter et al., 2006) |
| alpha | = b | needed for Bath's law to emerge (Helmstetter and Sornette, 2003) |

Free, physically meaningful:

| parameter | default | role |
|---|---|---|
| stress_drop_mpa, stress_drop_ln_sigma | 3, 0.25 | amplitude and scatter of slip |
| b_value | 1.0 | size hazard, GR slope |
| arrest_scale (rho_c) | 0.25 | strength of budget control on ruptures |
| branching_ratio | 0.5 | clustering strength; sets Bath's gap and site recurrence CV together |
| spatial_d_km, spatial_q | 1, 1.5 | aftershock distance kernel |

Choices that are conventions rather than parameters: patch size (1 km),
minimum magnitude (3.0), and the initial deficit (a smooth random field,
marginally uniform in [0, slip of a fault-wide event]); the first tenth
of a run should be treated as burn-in.

## 4. Known limitations of Stage 0

* With constant stress drop and a saturating width, mean slip saturates
  for long ruptures and M0 grows as L rather than the observed L^2
  (Scholz, 1982; Blaser et al., 2010). The honest version is kept; an
  effective depth extent of slip below the locked zone would be the one
  extra parameter to test.
* The free-surface image is exact for the antiplane components only.
* The backslip stress of rough slip is heterogeneous at the patch scale
  (tens of MPa); as an arrest measure it is noisy, which is why the
  deficit is the default state measure and stress is a diagnostic. The
  stress standard deviation scales with slip_cv and is probably too
  large at the default value; this needs a decision.
* Aftershocks are placed on the fault plane; off-fault seismicity is not
  represented.
* Sub-patch events do not debit the budget (their moment is a few
  percent of the total for b = 1 and a patch magnitude two to three units
  below the maximum).
* Single-fault results depend on the fault length through the maximum
  magnitude; a fault network changes this.
* The analytic calibration assumes every patch of the fault is reachable
  from every nucleation point; on segmented faults with step-overs that
  are open only some of the time it overestimates the boundary-limited
  sizes slightly.

## 5. References

Allmann, B. P., and P. M. Shearer (2009), Global variations of stress
drop for moderate and large earthquakes, J. Geophys. Res., 114, B01310.

Andrews, D. J. (1980), A stochastic fault model: 1. Static case,
J. Geophys. Res., 85, 3867-3877.

Bao, H., L. Xu, L. Meng, J.-P. Ampuero, L. Gao and H. Zhang (2022),
Global frequency of oceanic and continental supershear earthquakes,
Nature Geosci., 15, 942-949.

Blaser, L., F. Kruger, M. Ohrnberger and F. Scherbaum (2010), Scaling
relations of earthquake source parameter estimates with special focus on
subduction environment, Bull. Seismol. Soc. Am., 100, 2914-2926.

Field, E. H., et al. (2014), Uniform California Earthquake Rupture
Forecast, version 3 (UCERF3): The time-independent model, Bull. Seismol.
Soc. Am., 104, 1122-1180.

Graves, R. W., and A. Pitarka (2010), Broadband ground-motion simulation
using a hybrid approach, Bull. Seismol. Soc. Am., 100, 2095-2123.

Gutenberg, B., and C. F. Richter (1944), Frequency of earthquakes in
California, Bull. Seismol. Soc. Am., 34, 185-188.

Heaton, T. H. (1990), Evidence for and implications of self-healing
pulses of slip in earthquake rupture, Phys. Earth Planet. Inter., 64,
1-20.

Helmstetter, A., and D. Sornette (2003), Bath's law derived from the
Gutenberg-Richter law and from aftershock properties, Geophys. Res.
Lett., 30, 2069.

Helmstetter, A., Y. Y. Kagan and D. D. Jackson (2006), Comparison of
short-term and time-independent earthquake forecast models for southern
California, Bull. Seismol. Soc. Am., 96, 90-106.

Kanamori, H., and D. L. Anderson (1975), Theoretical basis of some
empirical relations in seismology, Bull. Seismol. Soc. Am., 65,
1073-1095.

LeVeque, R. J., K. Waagan, F. I. Gonzalez, D. Rim and G. Lin (2016),
Generating random earthquake events for probabilistic tsunami hazard
assessment, Pure Appl. Geophys., 173, 3671-3692.

Mai, P. M., and G. C. Beroza (2002), A spatial random field model to
characterize complexity in earthquake slip, J. Geophys. Res., 107, 2308.

Melgar, D., R. J. LeVeque, D. S. Dreger and R. M. Allen (2016),
Kinematic rupture scenarios and synthetic displacement data: An example
application to the Cascadia subduction zone, J. Geophys. Res. Solid
Earth, 121, 6658-6674.

Noda, H., N. Lapusta and H. Kanamori (2013), Comparison of average
stress drop measures for ruptures with heterogeneous stress change and
implications for earthquake physics, Geophys. J. Int., 193, 1691-1712.

Ogata, Y. (1988), Statistical models for earthquake occurrences and
residual analysis for point processes, J. Am. Stat. Assoc., 83, 9-27.

Pacheco, J. F., C. H. Scholz and L. R. Sykes (1992), Changes in
frequency-size relationship from small to large earthquakes, Nature,
355, 71-73.

Savage, J. C. (1983), A dislocation model of strain accumulation and
release at a subduction zone, J. Geophys. Res., 88, 4984-4996.

Scholz, C. H. (1982), Scaling laws for large earthquakes: Consequences
for physical models, Bull. Seismol. Soc. Am., 72, 1-14.

Page and volume numbers above were written from memory and should be
checked against erdbeben.bib before being cited.
