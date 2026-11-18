# Extensive margin 2011-2020 (cohort from 2006-2010; secondary, specified before running)

Extensive margin, long window (pilot Section 11 fallback): a cohort labelled from 2006-2010, outcome 2011-2020.

SPECIFIED 2026-11-18, BEFORE ANY RESULT OF IT WAS SEEN. Secondary analysis; not in the frozen decision log; no ET.
Why: the 2016-2020 window has 5 map years and no severe drought. LANID covers 1997-2020 and thickness is kriged
for 2011-2024, so 2011-2020 gives 10 years including the 2011-2012 drought, with no new data.
Cohort F: regular cells inside the aquifer, not served by surface water (D09), LANID irrigated (D04 rule, D06
coverage) in at least COHORT_MIN_IRR_YEARS (4) of 2006-2010, LANID known in all five years, and LANID and
AIM-HPA agreeing in at least D06 min_agree_years of 2006-2010 (definitions.maps_agree; AIM-HPA 2006-2010 read from
"pilot test/AIM HPA"; this rule added 2026-11-18 before the second run, replacing the no-agreement shortcut).
  shortcut: no corn rule (the outcome is irrigation of any crop); add it if this becomes a primary analysis.
S4 also draws AIM-HPA (2011-2017) for the same cells, dashed, as a check on single-year map artefacts.
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
| 0-10 | 76 | 16 | 35 | 5.29 |
| 10-20 | 79.2 | 12.4 | 36 | 4.05 |
| 20-30 | 84 | 7.8 | 35 | 2.59 |
| 30-40 | 84 | 8.49 | 32 | 2.94 |
| 40-50 | 87.2 | 8.96 | 30 | 3.21 |
| 50-60 | 88.7 | 7.5 | 27 | 2.83 |
| 60-70 | 88.8 | 8.61 | 22 | 3.6 |
| 70-80 | 88.8 | 11.9 | 16 | 5.85 |
| 80-100 | 95 | 3.6 | 10 | 2.23 |
| >=100 | 93.7 | 5.48 | 7 | 4.06 |

## S2_within

| A_bin | effect_pp | lo_pp | hi_pp | n |
|---|---|---|---|---|
| 0-10 | 0 | 0 | 0 | 6692 |
| 10-20 | 4.98 | 2.71 | 7.25 | 18758 |
| 20-30 | 6 | 2.9 | 9.11 | 21628 |
| 30-40 | 9.7 | 6.13 | 13.3 | 23168 |
| 40-50 | 11.8 | 7.86 | 15.7 | 19191 |
| 50-60 | 13.2 | 7.93 | 18.4 | 13563 |
| 60-70 | 14.9 | 9.11 | 20.7 | 10373 |
| 70-80 | 16.3 | 10.6 | 22.1 | 7055 |
| 80-100 | 17.9 | 11.8 | 24 | 6583 |
| >=100 | 16.7 | 11.7 | 21.8 | 6974 |

## S4_by_year

| year | thin | middle | thick |
|---|---|---|---|
| 2011 | 85.7 | 89.7 | 93.8 |
| 2012 | 81.5 | 87.2 | 89.5 |
| 2013 | 82.1 | 88.9 | 91.1 |
| 2014 | 83.9 | 89.4 | 92.6 |
| 2015 | 84.3 | 89.1 | 91.4 |
| 2016 | 81.8 | 88.2 | 91.5 |
| 2017 | 76.6 | 80.9 | 85.6 |
| 2018 | 81.4 | 88.1 | 91.9 |
| 2019 | 81.2 | 86.3 | 89.5 |
| 2020 | 81.7 | 88.2 | 91.5 |

## S4_by_year_AIM

| year | thin | middle | thick |
|---|---|---|---|
| 2011 | 85.6 | 87.5 | 91.1 |
| 2012 | 79.8 | 83.7 | 87.6 |
| 2013 | 77.6 | 84.2 | 88.3 |
| 2014 | 79.9 | 85.2 | 89.8 |
| 2015 | 83.4 | 87.7 | 90.8 |
| 2016 | 88.4 | 91.2 | 94.5 |
| 2017 | 84 | 88.8 | 90.9 |

## Summary

```
{
 "freeze": "pilot-freeze-2026-11-17",
 "years": "2011-2020",
 "cohort_from": "2006-2010",
 "cohort_f_cells": 14007,
 "cohort_f_cell_years": 133985,
 "thirds_cuts_m": [
  28.32580862519749,
  50.02632290139953
 ],
 "S1_cell_years": 533731,
 "S1_counties": 37,
 "S2_slope_per_10m": {
  "corrected": {
   "pp_per_10m": 2.739977932740019,
   "lo": 1.055521108311886,
   "hi": 4.424434757168152
  },
  "uncorrected": {
   "pp_per_10m": 1.5014875045291822,
   "lo": 0.5784177076609571,
   "hi": 2.4245573013974075
  }
 },
 "S3_drought_contrast_pp": {
  "value_mm": 1.403289412444602,
  "se_mm": 2.0209878985065974,
  "lo_mm": -2.695464020570715,
  "hi_mm": 5.502042845459919,
  "A_thin_m": 17.750245866534257,
  "A_thick_m": 75.13492369858454,
  "deficit_q90_mm": 1019.3415493826994,
  "deficit_median_mm": 591.5700238269401,
  "clusters": 37,
  "n": 133985,
  "cells": 13745
 }
}
```
