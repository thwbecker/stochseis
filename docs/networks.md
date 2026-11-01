# 3-D geometries and fault networks: design

Status: design agreed in outline on 1 November 2026; step N1 is
implemented (`stochseis/mesh.py`), the rest is planned. This document
explains how stochseis moves from single planar faults to the New
Zealand Alpine Fault reference geometries, then to UCERF3-type fault
systems and finally to probabilistic fault networks, without changing
the simulation rules established on single faults
(`docs/design.md`, `docs/experiments.md`).

## 1. Principle: sections are structured grids with a 3-D embedding

Everything the simulator does (rupture growth on a graph, von Karman
slip fields, spectral crack solutions, the deficit and stress state)
works on a regular grid of square patches indexed by an along-strike
coordinate s and a down-dip coordinate d. A *section* keeps that grid
and adds, per patch, its real 3-D position, normal, dip, strike, area,
slip-rate vector and rigidity. The grid may be ragged at the bottom
(variable dip at fixed seismogenic depth) and has gaps where a section is
absent.

Why grids rather than raw triangles:

* fault-section databases already describe faults this way (a trace, a
  dip or down-dip profile, upper and lower seismogenic depths: UCERF3
  subsections, the NZ Community Fault Model), so loaders for them only
  have to build the (s, d) -> (x, y, z) map;
* random slip fields with prescribed spectra are FFTs on a grid; on
  arbitrary triangles they need a Karhunen-Loeve expansion of an N x N
  covariance (LeVeque et al., 2016; Melgar et al., 2016), which is
  affordable per rupture but not for 20,000-patch faults;
* the spectral crack solver and the backslip state are fast on grids and
  exact for vertical strike-slip faults (E20).

Triangle meshes are converted to sections (N1). Geometries that are not
a single sheet with a top edge (the truncated Alpine model, intersecting
or branching faults) are split into several sections joined by
connections (N2).

## 2. Data model

`Fault` (current class) becomes the network container:

* per patch: s, d, col, row (local to the section), section id, x, y, z,
  normal, dip, strike, area, slip rate, rake, rigidity;
* per section: name, grid shape, trace, top and bottom depth, whether it
  reaches the surface, parent fault (for UCERF3 parent sections);
* connections: patch pairs across section ends or branch points, with 3-D
  separation, change of strike and of dip, and an opening probability
  (or a rule to compute it per rupture);
* a global adjacency graph (within-section 8-neighbour edges plus
  connection edges) used by rupture growth.

The planar segmented fault of Stage 0 is the special case of sections
laid end to end on one unrolled grid with gap columns; it stays as is.

## 3. Elasticity: three levels

L0 (now): spectral operator on the unrolled grid of the rupture, free
surface by an image. Exact for vertical planar strike-slip; for dipping
faults it gives 13 to 35 percent too little slip, most near a shallowly
dipping surface trace (E20). Curvature and dip changes along strike are
ignored by the operator (they enter only through patch areas and
positions).

L1 (next, cheap): a correction factor for the crack-solution amplitude
as a function of dip, rake and whether the rupture reaches the surface,
tabulated once from half-space solutions (`examples/checks/
halfspace_check.py` extended to a grid of dips). This fixes the mean slip
of dipping ruptures to a few percent but not the slip distribution near
the trace.

L2 (for networks and for publication-quality dipping faults): a
half-space stiffness for the whole network from triangular dislocations
(cutde; Nikkhoo and Walter, 2015) at the section-grid resolution, stored
as a hierarchical low-rank matrix (cutde's ACA routines) so that 20,000
to 100,000 patches fit in memory. Crack solutions of large ruptures then
use conjugate gradients with this operator and the L0 spectral operator
as preconditioner; small ruptures stay on L0. The same operator gives
stress transfer between sections, which L0 cannot.

Multi-section ruptures: on L0 a rupture along a chain of sections is
unrolled along its path (sections concatenated in s, with the step-over
as a gap), which is what the segmented planar fault does already. A
rupture that takes a branch has no single unrolling; L0 then treats each
branch as its own unrolled grid with the shared length entering the
correlation lengths, and L2 is the reference.

## 4. Loading and rake

Backslip loading with a slip-rate vector per patch: rate and rake from a
segment table, interpolated along strike within a section (the
`slip_rate_profile` parameter implements the scalar part). The deficit
stays a scalar measured along the local rake. The spectral kernel needs
the angle between slip and strike (it currently distinguishes only
strike from dip slip); generalising it to any rake is a one-line change
(the kernel depends on the component of the wavenumber along the slip
direction).

## 5. Observables at sites and for model comparison

`stochseis/sites.py` turns catalogs into synthetic paleoseismic records:
a site maps to the nearest surface patch, an event is recorded if its
slip there exceeds a detection threshold, and the record gives the mean
recurrence interval, its CV and the single-event displacement. Further
observables for comparison with physics-based simulators: co-rupture
fractions between sites (Howarth et al., 2021, report spatio-temporal
clustering on the Alpine Fault), the fraction of ruptures crossing
segment boundaries, rupture-length and magnitude distributions, and
hypocentre locations.

## 6. Tuning and the range of consistent models

The paleoseismic data constrain the model through a budget identity: if
detected events release most of the slip at a site, the mean interval is
about the mean single-event displacement divided by the slip rate
(central Alpine Fault: 7.1 m / 27 mm/yr = 263 yr against 249 +- 57 yr
observed). In stochseis the displacement per event is set by the stress
drop and the rupture width, and the regularity by the nucleation
exponent and the branching ratio (E6, E9), so the free parameters with
the most leverage are the stress-drop median and scatter, the nucleation
exponent, the branching ratio and the detection threshold.

Procedure: draw parameter vectors from broad priors (Latin hypercube or
random), run a catalog of several thousand years for each, compute the
site statistics, and keep the parameter vectors whose statistics fall
within the observational uncertainties (approximate Bayesian
computation by rejection). The accepted set is the "range of stochastic
models consistent with the constraints"; its event statistics are what
we compare with the physics-based simulators. A few hundred runs of
5000 yr on the 755 km Alpine Fault at 1 km patches cost about 8 minutes
each on one core, i.e. one to two days on a 16-core machine; an emulator
can follow if more is needed.

## 7. Stage plan

N1 (done): triangulated single-sheet surfaces to a section grid
(`mesh.load_mesh_section`), per-patch 3-D geometry, slip-rate profile
along strike, paleo-site observables. Alpine planar and variable-dip
geometries run.

N2: multi-section Fault container and connections; split meshes into
sections (truncated Alpine, Alpine plus Hope); section-from-trace
builder (trace, dip, depths) for database inputs; general-rake kernel;
L1 dip correction.

N3: UCERF3 loader (fault subsections from OpenSHA, typically GeoJSON:
trace, dip, dip direction, upper and lower depth, rake, slip rate,
aseismicity) and UCERF3 connection rules (jumps up to 5 km, azimuth and
rake changes) as defaults for opening probabilities; NZ CFM loader by the
same route.

N4: L2 half-space operator (cutde, hierarchical) for crack solutions and
inter-section stress transfer.

N5: probabilistic networks. A network generator draws sections (traces,
dips, depths, slip rates consistent with a regional loading constraint
such as a transect sum) and connections with gap- and angle-dependent
probabilities; each simulation realises one network. Nothing downstream
changes because the simulator sees only the Fault container; epistemic
spread then comes from sampling networks as well as parameters.
