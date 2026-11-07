# Roadmap

State as of 5 November 2026. Earlier steps (budget-control regime,
deficit-weighted nucleation, stress state, moment-length scaling, the
Alpine calibration rounds) are complete and documented in
`docs/experiments.md` and `examples/nz_alpine/README.md`; their history is
in git.

## Done (5 November 2026)

1. Single-fault parameter overview (`examples/single_fault/run_test.sh`;
   branching ratio, nucleation clock, size rule; 10 realizations of
   10,000 yr per case).
2. Alpine realization sets A to D, control summary, pooled
   constraint-weighted predictions, parameter effects over the ABC rounds
   (`examples/nz_alpine/run_test.sh`).

## Next

3. UCERF3 (`examples/ucerf3/`): multi-section fault container with
   connections (`docs/networks.md` stage N2), UCERF3 subsection loader and
   connection rules (N3), calibration against UCERF3 regional and
   section participation rates; see HANDOFF section 6.1.
4. Fault systems represented only probabilistically (N5): a network
   generator for sections and connections constrained by a regional
   loading rate, networks sampled along with parameters; see HANDOFF
   section 6.2.

## Also pending

5. E23 (`docs/experiments.md`): sweep of the hierarchical patch option on
   the single fault (barrier and nucleation exponents, radius range,
   nesting; at least 3000 yr burn-in).
6. Intercomparison submission: catalog export from the constrained
   Alpine ensemble once the organisers specify a format.
7. Alpine-Hope geometry (needs the Hope mesh; follows from N2); the L2
   half-space operator (cutde, N4) for accurate dip effects and stress
   transfer between sections.

## Housekeeping

* Line-by-line review of the code by a human; check the references in
  `docs/design.md` against the literature.
