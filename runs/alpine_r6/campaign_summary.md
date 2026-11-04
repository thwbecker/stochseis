Site-specific detection thresholds: north 1 m

## ABC results by configuration (2 m detection threshold)

| configuration | vectors | failed runs | accepted (paleo) | accepted (paleo + seismicity) | best misfit | best: south RI / CV / SED | central | north | b | Bath |
|---|---|---|---|---|---|---|---|---|---|---|
| planar_coarse | 100 | 0 | 47 | 41 | 1.28 | 290 / 0.55 / 7.3 | 277 / 0.56 / 7.0 | 515 / 0.35 / 4.7 | 0.99 | 1.35 |
| variable_dip_coarse | 100 | 0 | 53 | 54 | 1.03 | 278 / 0.51 / 7.1 | 275 / 0.43 / 7.3 | 442 / 0.59 / 4.3 | 0.97 | 1.29 |

## Accepted parameter ranges (paleo + seismicity): median [min-max]

| parameter | planar_coarse | variable_dip_coarse |
|---|---|---|
| stress_drop_mpa | 7.49 [5.08-9.93] | 7.75 [5.08-9.93] |
| nucleation_exponent | 5.49 [2.11-11.8] | 5.42 [2.11-11.8] |
| branching_ratio | 0.428 [0.208-0.691] | 0.434 [0.208-0.691] |
| arrest_scale | 0.663 [0.204-0.997] | 0.65 [0.204-0.99] |
| admission_exponent | 5.24 [2.58-7.92] | 5.27 [2.43-7.92] |
| b_value | 1.08 [0.913-1.19] | 1.08 [0.913-1.2] |
| budget_scale | 0.393 [0.114-0.984] | 0.349 [0.103-0.984] |
| slip_cv | 0.297 [0.106-0.497] | 0.288 [0.104-0.498] |
| hazard_max_area_km2 | 1.13e+03 [368-2.97e+03] | 953 [300-2.97e+03] |

## Posterior vectors rescored with their posterior runs (paleo + seismicity, 2 m)

ABC misfit: from the ABC runs; posterior misfit: from the posterior runs of the same vector.

| configuration | vectors | still accepted | best posterior misfit | its ABC rank | best: south RI / CV / SED | central | north | b | Bath | median site CV (south, central, north) |
|---|---|---|---|---|---|---|---|---|---|---|
| planar_coarse | 20 | 20 | 1.46 | 8 | 242 / 0.39 / 6.1 | 241 / 0.42 / 6.1 | 418 / 0.39 / 3.8 | 1.01 | 1.52 | 0.40, 0.41, 0.55 |
| variable_dip_coarse | 20 | 20 | 1.06 | 2 | 263 / 0.57 / 6.7 | 267 / 0.52 / 6.8 | 441 / 0.47 / 4.3 | 0.99 | 1.26 | 0.41, 0.43, 0.52 |

## Overlap of accepted vectors (paleo + seismicity): number accepted in both

| | planar_coarse | variable_dip_coarse |
|---|---|---|
| planar_coarse | 41 | 39 |
| variable_dip_coarse | 39 | 54 |
