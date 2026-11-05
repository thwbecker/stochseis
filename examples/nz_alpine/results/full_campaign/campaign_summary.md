## ABC results by configuration (2 m detection threshold)

| configuration | vectors | failed runs | accepted (paleo) | accepted (paleo + seismicity) | best misfit | best: south RI / CV / SED | central | north | b | Bath |
|---|---|---|---|---|---|---|---|---|---|---|
| planar_fine | 100 | 0 | 54 | 53 | 0.85 | 263 / 0.42 / 7.1 | 278 / 0.37 / 7.2 | 500 / 0.36 / 4.5 | 0.95 | 1.38 |
| variable_dip_fine | 100 | 0 | 57 | 59 | 0.64 | 265 / 0.41 / 7.0 | 268 / 0.41 / 7.0 | 466 / 0.38 / 4.2 | 0.96 | 1.38 |
| segmented_fine | 100 | 0 | 57 | 60 | 0.60 | 261 / 0.29 / 7.0 | 268 / 0.27 / 7.1 | 465 / 0.28 / 4.1 | 1.07 | 1.42 |
| planar_coarse | 100 | 0 | 52 | 56 | 0.99 | 254 / 0.32 / 6.8 | 250 / 0.31 / 6.7 | 484 / 0.36 / 4.2 | 0.89 | 1.53 |
| variable_dip_coarse | 100 | 0 | 64 | 65 | 0.48 | 264 / 0.46 / 7.0 | 271 / 0.41 / 7.1 | 442 / 0.33 / 4.2 | 1.02 | 1.18 |

## Accepted parameter ranges (paleo + seismicity): median [min-max]

| parameter | planar_fine | variable_dip_fine | segmented_fine | planar_coarse | variable_dip_coarse |
|---|---|---|---|---|---|
| stress_drop_mpa | 6.95 [5.3-9.95] | 7.43 [5.3-9.95] | 7.74 [5.37-9.95] | 7.31 [5.18-9.95] | 7.4 [5.18-9.95] |
| nucleation_exponent | 5.04 [2.24-12] | 4.91 [2.16-12] | 4.97 [2.16-12] | 5.31 [2.58-12] | 5.12 [2.24-12] |
| branching_ratio | 0.166 [0.0506-0.296] | 0.18 [0.0506-0.299] | 0.177 [0.0506-0.299] | 0.173 [0.0506-0.296] | 0.18 [0.0506-0.299] |
| arrest_scale | 0.674 [0.229-0.997] | 0.674 [0.229-0.999] | 0.673 [0.229-0.999] | 0.673 [0.229-0.999] | 0.674 [0.229-0.999] |
| admission_exponent | 4.48 [2-7.93] | 4.48 [2-7.93] | 4.36 [2-7.93] | 4.46 [2-7.93] | 4.6 [2-7.93] |
| b_value | 0.958 [0.802-1.1] | 0.935 [0.802-1.1] | 0.954 [0.802-1.1] | 0.944 [0.802-1.1] | 0.935 [0.802-1.1] |
| budget_scale | 0.365 [0.101-0.952] | 0.321 [0.101-0.952] | 0.36 [0.101-0.952] | 0.4 [0.101-0.952] | 0.365 [0.101-0.952] |
| slip_cv | 0.336 [0.117-0.495] | 0.326 [0.117-0.495] | 0.319 [0.117-0.495] | 0.339 [0.117-0.495] | 0.32 [0.117-0.495] |
| hazard_max_area_km2 | 1.08e+03 [308-2.82e+03] | 1.08e+03 [302-2.92e+03] | 1.11e+03 [302-2.92e+03] | 1.08e+03 [308-2.82e+03] | 1.14e+03 [302-2.92e+03] |

## Posterior vectors rescored with their posterior runs (paleo + seismicity, 2 m)

ABC misfit: from the ABC runs; posterior misfit: from the posterior runs of the same vector.

| configuration | vectors | still accepted | best posterior misfit | its ABC rank | best: south RI / CV / SED | central | north | b | Bath | median site CV (south, central, north) |
|---|---|---|---|---|---|---|---|---|---|---|
| planar_fine | 20 | 20 | 0.80 | 1 | 266 / 0.40 / 7.0 | 274 / 0.37 / 7.2 | 497 / 0.35 / 4.5 | 0.96 | 1.38 | 0.38, 0.37, 0.39 |
| variable_dip_fine | 20 | 20 | 0.69 | 2 | 253 / 0.33 / 6.8 | 257 / 0.32 / 6.8 | 471 / 0.31 / 4.1 | 1.08 | 1.44 | 0.41, 0.40, 0.37 |
| segmented_fine | 20 | 20 | 0.60 | 1 | 250 / 0.33 / 6.7 | 252 / 0.34 / 6.7 | 421 / 0.30 / 3.9 | 1.07 | 1.43 | 0.41, 0.39, 0.40 |
| planar_coarse | 20 | 20 | 0.85 | 2 | 281 / 0.40 / 7.4 | 284 / 0.41 / 7.5 | 502 / 0.37 / 4.3 | 1.02 | 1.21 | 0.40, 0.39, 0.37 |
| variable_dip_coarse | 20 | 20 | 0.41 | 2 | 260 / 0.33 / 7.0 | 261 / 0.32 / 6.9 | 430 / 0.35 / 3.9 | 1.10 | 1.29 | 0.41, 0.41, 0.39 |

## Overlap of accepted vectors (paleo + seismicity): number accepted in both

| | planar_fine | variable_dip_fine | segmented_fine | planar_coarse | variable_dip_coarse |
|---|---|---|---|---|---|
| planar_fine | 53 | 46 | 46 | 49 | 50 |
| variable_dip_fine | 46 | 59 | 56 | 45 | 57 |
| segmented_fine | 46 | 56 | 60 | 47 | 58 |
| planar_coarse | 49 | 45 | 47 | 56 | 52 |
| variable_dip_coarse | 50 | 57 | 58 | 52 | 65 |
