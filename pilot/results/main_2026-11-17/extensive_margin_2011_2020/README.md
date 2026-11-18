# Extensive margin 2011-2020 (cohort from 2006-2010; secondary, specified before running)

Extensive margin, long window (pilot Section 11 fallback): a cohort labelled from 2006-2010, outcome 2011-2020.

SPECIFIED 2026-11-18, BEFORE ANY RESULT OF IT WAS SEEN. Secondary analysis; not in the frozen decision log; no ET.
Why: the 2016-2020 window has 5 map years and no severe drought. LANID covers 1997-2020 and thickness is kriged
for 2011-2024, so 2011-2020 gives 10 years including the 2011-2012 drought, with no new data.
Cohort F: regular cells inside the aquifer, not served by surface water (D09), LANID irrigated (D04 rule, D06
coverage) in at least COHORT_MIN_IRR_YEARS (4) of 2006-2010, LANID known in all five years.
  shortcut: no corn rule and no AIM-HPA agreement check (the outcome is irrigation of any crop); add both if this
  becomes a primary analysis.
Outcome: irrigated per LANID (D04), 2011-2020; unknown cell-years left out. Thickness: cell-year A of Check 2.4
(frozen main run), poor cell-years (D25) left out; thirds re-cut for cohort F by the D26 rule on its 2011-2020 means.
Deficit: seasonal (May-Sep) reference ET minus rain, D01, county-area-weighted gridMET, 2011-2020 (the Task 1.1
computation, kept here so the frozen module stays untouched).
Control: root-zone storage. Corn years 2013-2015 are NOT a control here: they fall inside the outcome window and are
partly an outcome of irrigation.
  S1  share of baseline cropland (the 2016-2020 analysis's population) irrigated, 2011-2020, by A bin.
  S2  cohort F still irrigated, by A bin, and within-county slope per 10 m (A error corrected).
  S3  cohort F: thin minus thick change in irrigated (pp) in a dry season (q90 vs median deficit), cell and
      county x year effects, A error corrected (D30 form).
  S4  (display) share of cohort F irrigated each year 2011-2020, by thickness third, with the deficit.
Expected signs: S1, S2 slopes positive; S3 contrast negative. Evidence: 95% interval excluding zero.
Run:  PILOT_RUN_MODE=main python -m pilot.phase3.extensive_margin_2011_2020


## S1_between

| A_bin | mean | std | counties | half95 |
|---|---|---|---|---|
| 0-10 | 20 | 25.8 | 36 | 8.42 |
| 10-20 | 21.8 | 18.3 | 37 | 5.89 |
| 20-30 | 25.3 | 18 | 37 | 5.8 |
| 30-40 | 32.2 | 19.9 | 35 | 6.58 |
| 40-50 | 36.1 | 20.7 | 32 | 7.19 |
| 50-60 | 41.8 | 24.4 | 30 | 8.71 |
| 60-70 | 50.1 | 27.7 | 25 | 10.9 |
| 70-80 | 55.7 | 23.2 | 17 | 11.1 |
| 80-100 | 59 | 22.7 | 10 | 14.1 |
| >=100 | 67.6 | 21.4 | 7 | 15.9 |

## S1_within

| A_bin | effect_pp | lo_pp | hi_pp | n |
|---|---|---|---|---|
| 0-10 | 0 | 0 | 0 | 82987 |
| 10-20 | 5.8 | 3.27 | 8.33 | 121243 |
| 20-30 | 9.33 | 5.4 | 13.3 | 106050 |
| 30-40 | 14.3 | 8.37 | 20.3 | 81090 |
| 40-50 | 16.6 | 9.14 | 24.2 | 55825 |
| 50-60 | 21.2 | 10.9 | 31.4 | 31342 |
| 60-70 | 27.3 | 16.2 | 38.3 | 18967 |
| 70-80 | 27.1 | 16.7 | 37.5 | 12695 |
| 80-100 | 26.2 | 13.4 | 39.1 | 11654 |
| >=100 | 18.4 | -3.74 | 40.6 | 11878 |

## S2_between

| A_bin | mean | std | counties | half95 |
|---|---|---|---|---|
| 0-10 | 75.6 | 16.3 | 35 | 5.39 |
| 10-20 | 79.1 | 12.4 | 36 | 4.04 |
| 20-30 | 83.8 | 8.02 | 35 | 2.66 |
| 30-40 | 83.8 | 8.45 | 32 | 2.93 |
| 40-50 | 87.1 | 8.95 | 30 | 3.2 |
| 50-60 | 88 | 8.09 | 27 | 3.05 |
| 60-70 | 88.7 | 8.64 | 22 | 3.61 |
| 70-80 | 86.9 | 13.9 | 16 | 6.82 |
| 80-100 | 94.9 | 3.66 | 10 | 2.27 |
| >=100 | 93.7 | 5.47 | 7 | 4.05 |

## S2_within

| A_bin | effect_pp | lo_pp | hi_pp | n |
|---|---|---|---|---|
| 0-10 | 0 | 0 | 0 | 6737 |
| 10-20 | 5.17 | 2.87 | 7.48 | 18902 |
| 20-30 | 6.14 | 3.03 | 9.26 | 21836 |
| 30-40 | 9.86 | 6.25 | 13.5 | 23461 |
| 40-50 | 12.1 | 8.09 | 16.1 | 19416 |
| 50-60 | 13.4 | 8.07 | 18.6 | 13722 |
| 60-70 | 15.2 | 9.34 | 21 | 10442 |
| 70-80 | 16.5 | 10.7 | 22.4 | 7092 |
| 80-100 | 18.2 | 12.1 | 24.3 | 6596 |
| >=100 | 17.2 | 12.1 | 22.3 | 6994 |

## S4_by_year

| year | thin | middle | thick |
|---|---|---|---|
| 2011 | 85.4 | 89.3 | 93.6 |
| 2012 | 81.2 | 86.9 | 89.2 |
| 2013 | 81.9 | 88.7 | 91 |
| 2014 | 83.7 | 89.1 | 92.5 |
| 2015 | 84.1 | 89 | 91.3 |
| 2016 | 81.6 | 87.8 | 91.4 |
| 2017 | 76.5 | 80.6 | 85.5 |
| 2018 | 81.2 | 87.9 | 91.6 |
| 2019 | 81 | 86.1 | 89.3 |
| 2020 | 81.4 | 88.1 | 91.4 |

## Summary

```
{
 "freeze": "pilot-freeze-2026-11-17",
 "years": "2011-2020",
 "cohort_from": "2006-2010",
 "cohort_f_cells": 14132,
 "cohort_f_cell_years": 135198,
 "thirds_cuts_m": [
  28.347176270360542,
  49.97252693919637
 ],
 "S1_cell_years": 533731,
 "S1_counties": 37,
 "S2_slope_per_10m": {
  "corrected": {
   "pp_per_10m": 2.8021059485600706,
   "lo": 1.1128898294502025,
   "hi": 4.491322067669938
  },
  "uncorrected": {
   "pp_per_10m": 1.5322280993737836,
   "lo": 0.6085426816452664,
   "hi": 2.4559135171023008
  }
 },
 "S3_drought_contrast_pp": {
  "value_mm": 1.4382237797014727,
  "se_mm": 2.0238805460121925,
  "lo_mm": -2.666396214366893,
  "hi_mm": 5.542843773769839,
  "A_thin_m": 17.763723585613892,
  "A_thick_m": 75.01647010680976,
  "deficit_q90_mm": 1019.3415493826994,
  "deficit_median_mm": 591.5700238269401,
  "clusters": 37,
  "n": 135198,
  "cells": 13869
 }
}
```
