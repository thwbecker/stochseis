# Roadmap

State as of 5 November 2026. Earlier steps (budget-control regime,
deficit-weighted nucleation, stress state, moment-length scaling, the
Alpine calibration rounds) are complete and documented in
`docs/experiments.md` and `examples/nz_alpine/README.md`; their history is
in git.

## Now

1. `examples/single_fault/run_test.sh`: the single-fault parameter
   overview (branching ratio, nucleation clock, size rule; 10
   realizations of 10,000 yr per case), figures for
   `examples/single_fault/README.md`.
2. `examples/nz_alpine/run_test.sh`: realization set D (segment-boundary
   barriers), the control summary of sets A to D, the pooled
   constraint-weighted predictions, and parameter effects over the ABC
   rounds; then update `examples/nz_alpine/README.md` sections 3 and 4.3.

## Next

3. E23 (`docs/experiments.md`): sweep of the hierarchical patch option on
   the single fault (barrier and nucleation exponents, radius range,
   nesting; at least 3000 yr burn-in).
4. Intercomparison submission: catalog export from the constrained
   ensemble once the organisers specify a format.
5. Alpine-Hope geometry (needs the Hope mesh and multi-section faults,
   `docs/networks.md` stage N2); the L2 half-space operator (cutde) for
   accurate dip effects and inter-section stress transfer.

## Later

6. `examples/ucerf3/`: UCERF3 fault sections as a network of connected
   faults (jump distances, connection probabilities, section slip rates).
7. A fault system represented only probabilistically: a stochastic
   fault-network generator on the same patch-graph interface, with slip
   rates constrained by a regional loading rate.

## Housekeeping

* Choose a license.
* Line-by-line review of the code by a human; check the references in
  `docs/design.md` against the literature.
