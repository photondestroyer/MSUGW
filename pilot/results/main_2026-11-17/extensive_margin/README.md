# Extensive margin (secondary; specified before running)

Extensive margin: irrigated area against aquifer thickness (the analogue of Mieno et al. 2024, Fig. 3).

SPECIFIED 2026-11-18, BEFORE ANY RESULT OF IT WAS SEEN. Not in the frozen decision log: a secondary analysis.
It reads no ET (the blinding rule covers ET and surface temperature only).
Outcome: a cell-year is irrigated when its LANID irrigated share meets D04 (definitions.is_irrigated_year),
2016-2020 (LANID ends in 2020); unknown where LANID covers too little of the cell (D06), and then left out.
Thickness: cell-year A of Check 2.4 (frozen main run); poor-thickness cell-years (D25) left out.
Baseline cropland: regular cells inside the aquifer, not served by surface water (D09), cropland in every year
2011-2015 (the D07 cropland rule).
  S1  share of baseline cropland irrigated, by 10 m bins of A: (a) mean of county means (between places, the
      Mieno Fig. 3 analogue); (b) bin effects with county x year effects and root-zone storage (within county).
  S2  cohort cells (irrigated in 2011-2015): share still irrigated, by A bin as in S1, and the within-county slope
      per 10 m of A (county x year effects, root-zone storage, corn years; A error corrected as in main_models).
  S3  cohort cells: does irrigation drop more on thin aquifer in dry seasons? Cell and county x year effects,
      A, A x seasonal deficit, controls x deficit (the D30 form); contrast = beta[A x deficit] x (A_thin - A_thick)
      x (deficit q90 - median), in percentage points.
Expected signs: irrigated share rises with A (S1, S2 positive slopes); S3 contrast negative.
Evidence: a 95% interval excluding zero (county-clustered; t with G-1 df below 50 counties).
Run:  PILOT_RUN_MODE=main python -m pilot.phase3.extensive_margin


## S1_between

| A_bin | mean | std | counties | half95 |
|---|---|---|---|---|
| 0-10 | 18.1 | 22.3 | 34 | 7.49 |
| 10-20 | 21.8 | 17.8 | 37 | 5.74 |
| 20-30 | 25.5 | 17.9 | 37 | 5.78 |
| 30-40 | 32.1 | 19.4 | 35 | 6.41 |
| 40-50 | 36.5 | 21.4 | 31 | 7.53 |
| 50-60 | 40.8 | 23.6 | 30 | 8.46 |
| 60-70 | 50.4 | 25.8 | 23 | 10.6 |
| 70-80 | 46.9 | 23.8 | 16 | 11.7 |
| 80-100 | 51.5 | 31.8 | 10 | 19.7 |
| >=100 | 63.1 | 24.3 | 6 | 19.5 |

## S1_within

| A_bin | effect_pp | lo_pp | hi_pp | n |
|---|---|---|---|---|
| 0-10 | 0 | 0 | 0 | 42818 |
| 10-20 | 6.03 | 3.31 | 8.76 | 61274 |
| 20-30 | 9.8 | 5.68 | 13.9 | 53606 |
| 30-40 | 15.1 | 8.82 | 21.4 | 40896 |
| 40-50 | 16.4 | 8.64 | 24.2 | 27659 |
| 50-60 | 21.2 | 10.4 | 32 | 15151 |
| 60-70 | 27 | 16.1 | 37.8 | 9379 |
| 70-80 | 27.1 | 17.4 | 36.7 | 6012 |
| 80-100 | 25.5 | 11.8 | 39.1 | 5292 |
| >=100 | 17.8 | -3.85 | 39.4 | 5674 |

## S2_between

| A_bin | mean | std | counties | half95 |
|---|---|---|---|---|
| 0-10 | 87.8 | 11.4 | 31 | 4.02 |
| 10-20 | 88.3 | 5.93 | 33 | 2.02 |
| 20-30 | 90.1 | 6.62 | 34 | 2.22 |
| 30-40 | 90.1 | 6.78 | 31 | 2.39 |
| 40-50 | 91.2 | 5.67 | 29 | 2.06 |
| 50-60 | 91.4 | 7.89 | 27 | 2.98 |
| 60-70 | 91.9 | 9.87 | 21 | 4.22 |
| 70-80 | 78.3 | 34 | 14 | 17.8 |
| 80-100 | 93.8 | 5.22 | 8 | 3.62 |
| >=100 | 96.3 | 3.77 | 6 | 3.02 |

## S2_within

| A_bin | effect_pp | lo_pp | hi_pp | n |
|---|---|---|---|---|
| 0-10 | 0 | 0 | 0 | 1615 |
| 10-20 | 3.78 | 0.1 | 7.46 | 4961 |
| 20-30 | 4.98 | 1.47 | 8.49 | 5557 |
| 30-40 | 6.21 | 2.15 | 10.3 | 6588 |
| 40-50 | 6.04 | 2.12 | 9.95 | 5350 |
| 50-60 | 5.59 | 0.81 | 10.4 | 3317 |
| 60-70 | 7.03 | 2.52 | 11.6 | 2961 |
| 70-80 | 6.03 | 1.1 | 10.9 | 2219 |
| 80-100 | 6.73 | 2.12 | 11.3 | 1893 |
| >=100 | 9.03 | 4.71 | 13.3 | 2181 |

## Summary

```
{
 "freeze": "pilot-freeze-2026-11-17",
 "S1_cell_years": 267761,
 "S1_cells": 55429,
 "S1_counties": 37,
 "S2_cell_years": 36642,
 "S2_cells": 7540,
 "S2_slope_per_10m": {
  "corrected": {
   "pp_per_10m": 0.8862490220993985,
   "lo": 0.32605731371878843,
   "hi": 1.4464407304800084
  },
  "uncorrected": {
   "pp_per_10m": 0.49727947366062475,
   "lo": 0.18295265247818954,
   "hi": 0.81160629484306
  }
 },
 "S3_drought_contrast_pp": {
  "value_mm": 0.999416752183351,
  "se_mm": 2.488488105237053,
  "lo_mm": -4.047471045559123,
  "hi_mm": 6.046304549925826,
  "A_thin_m": 18.070421575120108,
  "A_thick_m": 77.62691372442625,
  "deficit_q90_mm": 715.9075403386933,
  "deficit_median_mm": 538.9184256262405,
  "clusters": 37,
  "n": 36642,
  "cells": 7540
 }
}
```
