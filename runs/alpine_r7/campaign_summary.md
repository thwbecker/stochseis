Site-specific detection thresholds: north 1 m

## ABC results by configuration (2 m detection threshold)

| configuration | vectors | failed runs | accepted (paleo) | accepted (paleo + seismicity) | best misfit | best: south RI / CV / SED | central | north | b | Bath |
|---|---|---|---|---|---|---|---|---|---|---|
| planar_coarse | 100 | 0 | 52 | 58 | 0.56 | 263 / 0.21 / 7.1 | 263 / 0.20 / 7.0 | 478 / 0.33 / 4.5 | 1.12 | 1.31 |
| variable_dip_coarse | 100 | 0 | 40 | 46 | 0.51 | 261 / 0.24 / 6.9 | 269 / 0.19 / 7.1 | 388 / 0.31 / 3.8 | 1.11 | 1.36 |

## Accepted parameter ranges (paleo + seismicity): median [min-max]

| parameter | planar_coarse | variable_dip_coarse |
|---|---|---|
| stress_drop_mpa | 8.11 [5.19-9.99] | 8.52 [6.25-9.99] |
| nucleation_exponent | 8.53 [3.22-26.7] | 8.53 [3.24-26.7] |
| branching_ratio | 0.166 [0.0535-0.293] | 0.169 [0.0535-0.293] |
| arrest_scale | 0.514 [0.203-0.993] | 0.65 [0.203-0.993] |
| admission_exponent | 7.6 [4.03-12] | 7.3 [4.03-12] |
| b_value | 0.941 [0.802-1.1] | 0.935 [0.802-1.1] |
| budget_scale | 0.281 [0.108-0.893] | 0.263 [0.108-0.893] |
| slip_cv | 0.331 [0.106-0.475] | 0.333 [0.106-0.467] |
| hazard_max_area_km2 | 200 [58.4-489] | 180 [52.4-489] |

## Posterior vectors rescored with their posterior runs (paleo + seismicity, 2 m)

ABC misfit: from the ABC runs; posterior misfit: from the posterior runs of the same vector.

| configuration | vectors | still accepted | best posterior misfit | its ABC rank | best: south RI / CV / SED | central | north | b | Bath | median site CV (south, central, north) |
|---|---|---|---|---|---|---|---|---|---|---|
| planar_coarse | 20 | 20 | 0.60 | 4 | 262 / 0.13 / 7.1 | 258 / 0.17 / 6.9 | 466 / 0.33 / 4.4 | 1.05 | 2.13 | 0.27, 0.26, 0.35 |
| variable_dip_coarse | 20 | 20 | 0.39 | 1 | 262 / 0.24 / 7.0 | 264 / 0.23 / 7.0 | 417 / 0.32 / 4.1 | 1.10 | 1.35 | 0.31, 0.28, 0.38 |

## Overlap of accepted vectors (paleo + seismicity): number accepted in both

| | planar_coarse | variable_dip_coarse |
|---|---|---|
| planar_coarse | 58 | 43 |
| variable_dip_coarse | 43 | 46 |
