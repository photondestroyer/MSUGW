# Support table: one-page summary (Task 1.7)

Mode stage3. Cells in scope: 128156.

## Sample masks

| mask | counties | cohort_cells | rainfed_cells | cohort_cell_months | county_years | dry_county_years | dry_episodes | dry_episodes_multiyear | deficit_peak_min | deficit_peak_max | S_min | S_max | area_ha |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mask_sm | 33 | 7483 | 16445 | 54394 | 244 | 74 | 11 | 7 | 80.2 | 269 | -1.47 | 1.15 | 4.85e+05 |
| mask_peak | 37 | 7784 | 16728 | 68683 | 332 | 100 | 12 | 7 | 80.2 | 269 | -1.47 | 1.15 | 5.04e+05 |
| mask_common | 33 | 7483 | 16445 | 54394 | 244 | 74 | 11 | 7 | 80.2 | 269 | -1.47 | 1.15 | 4.85e+05 |
| mask_timing | 37 | 7784 | 16728 | 336018 | 333 | 100 | 12 | 7 | 80.2 | 269 | -1.71 | 1.92 | 5.04e+05 |

## Exclusions

| reason | cells |
|---|---|
| not rainfed: irrigated in 2011-2015 | 34861 |
| not rainfed: not cropland | 29218 |
| retained | 24512 |
| outside the aquifer polygon | 19768 |
| irrigated, not cohort: irrigated < 4 of 2011-2015 | 9757 |
| irrigated, not cohort: corn < 2 of 2013-2015 | 6244 |
| not rainfed: irregular cell | 1426 |
| not rainfed: no irrigation map | 1264 |
| no month with >= 4 ET models | 443 |
| irrigated, not cohort: irregular cell | 377 |
| irrigated, not cohort: surface water | 258 |
| irrigated, not cohort: irrigation maps disagree | 25 |
| no measured well within 8 km in any year | 3 |

## Area

| role | retained | area_ha |
|---|---|---|
| cohort | False | 3.62e+03 |
| cohort | True | 5.04e+05 |
| excluded | False | 6.66e+06 |
| rainfed | False | 2.5e+04 |
| rainfed | True | 1.08e+06 |