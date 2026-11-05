Site-specific detection thresholds: north 1 m

## ABC results by configuration (2 m detection threshold)

| configuration | vectors | failed runs | accepted (paleo) | accepted (paleo + seismicity) | best misfit | best: south RI / CV / SED | central | north | b | Bath |
|---|---|---|---|---|---|---|---|---|---|---|
| planar_fine | 100 | 0 | 68 | 64 | 0.74 | 263 / 0.42 / 7.1 | 278 / 0.37 / 7.2 | 421 / 0.43 / 4.0 | 0.95 | 1.38 |
| variable_dip_fine | 100 | 0 | 64 | 68 | 0.55 | 265 / 0.41 / 7.0 | 268 / 0.41 / 7.0 | 402 / 0.36 / 3.8 | 0.96 | 1.38 |
| segmented_fine | 100 | 0 | 65 | 64 | 0.59 | 261 / 0.29 / 7.0 | 268 / 0.27 / 7.1 | 385 / 0.36 / 3.7 | 1.07 | 1.42 |
| planar_coarse | 100 | 0 | 63 | 64 | 0.68 | 265 / 0.41 / 7.1 | 284 / 0.37 / 7.5 | 435 / 0.38 / 4.0 | 1.01 | 1.22 |
| variable_dip_coarse | 100 | 0 | 66 | 69 | 0.49 | 264 / 0.46 / 7.0 | 271 / 0.41 / 7.1 | 402 / 0.37 / 3.9 | 1.02 | 1.18 |

## Accepted parameter ranges (paleo + seismicity): median [min-max]

| parameter | planar_fine | variable_dip_fine | segmented_fine | planar_coarse | variable_dip_coarse |
|---|---|---|---|---|---|
| stress_drop_mpa | 7.31 [5.3-9.95] | 7.81 [5.3-9.95] | 7.81 [5.37-9.95] | 7.31 [5.18-9.95] | 7.61 [5.3-9.95] |
| nucleation_exponent | 5.19 [2.16-12] | 4.97 [2.16-12] | 4.97 [2.16-12] | 5.19 [2.24-12] | 5.04 [2.16-12] |
| branching_ratio | 0.17 [0.0506-0.296] | 0.176 [0.0506-0.299] | 0.17 [0.0506-0.299] | 0.17 [0.0506-0.296] | 0.179 [0.0506-0.299] |
| arrest_scale | 0.693 [0.229-0.999] | 0.689 [0.229-0.999] | 0.685 [0.229-0.999] | 0.685 [0.229-0.999] | 0.686 [0.229-0.999] |
| admission_exponent | 4.48 [2-7.93] | 4.48 [2-7.93] | 4.46 [2-7.93] | 4.54 [2-7.93] | 4.49 [2-7.93] |
| b_value | 0.954 [0.802-1.1] | 0.954 [0.802-1.1] | 0.954 [0.802-1.1] | 0.955 [0.802-1.1] | 0.953 [0.802-1.1] |
| budget_scale | 0.334 [0.101-0.952] | 0.344 [0.101-0.952] | 0.357 [0.101-0.952] | 0.357 [0.101-0.952] | 0.339 [0.101-0.952] |
| slip_cv | 0.312 [0.117-0.495] | 0.319 [0.117-0.495] | 0.312 [0.117-0.495] | 0.311 [0.117-0.495] | 0.307 [0.117-0.495] |
| hazard_max_area_km2 | 1.09e+03 [308-2.92e+03] | 1.13e+03 [302-2.92e+03] | 1.13e+03 [302-2.92e+03] | 1.13e+03 [308-2.82e+03] | 1.14e+03 [302-2.92e+03] |

## Posterior vectors rescored with their posterior runs (paleo + seismicity, 2 m)

ABC misfit: from the ABC runs; posterior misfit: from the posterior runs of the same vector.

| configuration | vectors | still accepted | best posterior misfit | its ABC rank | best: south RI / CV / SED | central | north | b | Bath | median site CV (south, central, north) |
|---|---|---|---|---|---|---|---|---|---|---|
| planar_fine | 20 | 20 | 0.66 | 1 | 266 / 0.40 / 7.0 | 274 / 0.37 / 7.2 | 436 / 0.39 / 4.2 | 0.96 | 1.38 | 0.38, 0.37, 0.43 |
| variable_dip_fine | 20 | 20 | 0.65 | 1 | 252 / 0.41 / 6.6 | 247 / 0.43 / 6.5 | 417 / 0.36 / 3.9 | 0.96 | 1.37 | 0.41, 0.40, 0.41 |
| segmented_fine | 20 | 20 | 0.64 | 1 | 250 / 0.33 / 6.7 | 252 / 0.34 / 6.7 | 380 / 0.32 / 3.7 | 1.07 | 1.43 | 0.41, 0.39, 0.45 |
| planar_coarse | 20 | 20 | 0.69 | 2 | 281 / 0.40 / 7.4 | 284 / 0.41 / 7.5 | 428 / 0.40 / 4.0 | 1.02 | 1.21 | 0.40, 0.39, 0.41 |
| variable_dip_coarse | 20 | 20 | 0.49 | 2 | 260 / 0.33 / 7.0 | 261 / 0.32 / 6.9 | 362 / 0.38 / 3.5 | 1.10 | 1.29 | 0.41, 0.41, 0.44 |

## Overlap of accepted vectors (paleo + seismicity): number accepted in both

| | planar_fine | variable_dip_fine | segmented_fine | planar_coarse | variable_dip_coarse |
|---|---|---|---|---|---|
| planar_fine | 64 | 57 | 56 | 59 | 60 |
| variable_dip_fine | 57 | 68 | 63 | 57 | 66 |
| segmented_fine | 56 | 63 | 64 | 56 | 64 |
| planar_coarse | 59 | 57 | 56 | 64 | 60 |
| variable_dip_coarse | 60 | 66 | 64 | 60 | 69 |
