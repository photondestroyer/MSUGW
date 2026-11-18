# EXPLORATORY: T1 by month (not pre-registered)

EXPLORATORY (not pre-registered): the T1 timing model with its thickness terms split by calendar month.

Frozen T1 (D30) has A x peak-month indicator and A x exposure x peak indicator. Here the same monthly data, fixed
effects (cell x month, county x month x year), controls and thickness-error correction (main_models) are used,
with one A term and one A x monthly-deficit term PER MONTH (May-September):
    A x month_m            within a cell, how month m's ET moves with that cell's own A from year to year
    A x deficit x month_m  in month m, how the ET response to that month's deficit differs with A
For each month: contrast_m = beta[A x deficit x month_m] * (A_thin - A_thick) * (deficit q90_m - median_m),
the thin-minus-thick ET difference (mm) at a high deficit of that month against its median.
Reads the frozen main-mode outputs; writes only to <out>/explore/t1_by_month (never into a frozen task folder).
Run:  PILOT_RUN_MODE=main python -m pilot.phase3.explore_t1_by_month


| month | n | beta_A | lo_A | hi_A | contrast_value_mm | contrast_lo_mm | contrast_hi_mm |
|---|---|---|---|---|---|---|---|
| 5 | 56633 | 0.1 | 0.056 | 0.143 | 5.43 | 3.31 | 7.55 |
| 6 | 65353 | -0.07 | -0.141 | 0.001 | -1.88 | -3.76 | -0.011 |
| 7 | 63662 | 0.024 | -0.078 | 0.126 | 0.356 | -2.07 | 2.78 |
| 8 | 64075 | -0.079 | -0.174 | 0.015 | -2.19 | -4.52 | 0.143 |
| 9 | 64391 | -0.074 | -0.194 | 0.046 | -1.24 | -3.14 | 0.659 |