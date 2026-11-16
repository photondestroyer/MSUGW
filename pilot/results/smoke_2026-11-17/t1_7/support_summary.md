# Support table: one-page summary (Task 1.7)

Mode smoke. Cells in scope: 2304.

## Sample masks

| mask | counties | cohort_cells | rainfed_cells | cohort_cell_months | county_years | dry_county_years | dry_episodes | dry_episodes_multiyear | deficit_peak_min | deficit_peak_max | S_min | S_max | area_ha |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mask_sm | 1 | 480 | 297 | 3840 | 8 | 2 | 2 | 2 | 114 | 228 | -0.33 | 0.581 | 3.12e+04 |
| mask_peak | 1 | 480 | 297 | 4320 | 9 | 3 | 3 | 2 | 114 | 228 | -0.33 | 0.581 | 3.12e+04 |
| mask_common | 1 | 480 | 297 | 3840 | 8 | 2 | 2 | 2 | 114 | 228 | -0.33 | 0.581 | 3.12e+04 |
| mask_timing | 1 | 480 | 297 | 20539 | 9 | 3 | 3 | 2 | 114 | 228 | -1.25 | 1.54 | 3.12e+04 |

## Exclusions

| reason | cells |
|---|---|
| retained | 777 |
| irrigated, not cohort: corn < 2 of 2013-2015 | 456 |
| not rainfed: irrigated in 2011-2015 | 429 |
| irrigated, not cohort: irrigated < 4 of 2011-2015 | 414 |
| not rainfed: not cropland | 213 |
| irrigated, not cohort: surface water | 14 |
| irrigated, not cohort: irrigation maps disagree | 1 |

## Area

| role | retained | area_ha |
|---|---|---|
| cohort | True | 3.12e+04 |
| excluded | False | 9.92e+04 |
| rainfed | True | 1.93e+04 |