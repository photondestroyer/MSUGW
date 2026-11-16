# Processing test log (Task 1.6)

Mode: smoke. Counties: [np.int64(20081)]. Cells: 75 ({'cohort': 50, 'rainfed': 25}).

No feasibility conclusion is drawn from this sample (pilot Section 1).

## Checks

| check | passes | detail |
|---|---|---|
| coverage rows only for cells in scope | True | 2304 cells |
| ET rows only for cohort or rainfed cells | True |  |
| no cohort ET in dry months (gate) | True |  |
| every ET row has a county and a dry flag | True |  |
| one row per cell-month | True |  |
| ensemble ET within 0-400 mm/month | True | min 13.11781206171108, max 218.39316239316238 |
| dry flags cover every study county-month | True |  |
| blinding audit of the whole output root | True | 133 ET tables; [] |

## Seconds per task (last completed run)

```
{
 "t0_1": 51.9,
 "t0_4": 37.6,
 "counties": 5.7,
 "t0_2": 72.2,
 "t1_1": 6.2,
 "t1_2": 40.3,
 "t1_3": 5.7,
 "dry": 0.1,
 "t1_5": 7.8,
 "sample": 0.1,
 "t1_4": 568.1,
 "et_cells": 994.0,
 "t1_6": 3.2,
 "t1_inv": 0.1,
 "t1_7": 1.3,
 "rep": 17.3,
 "c2_1": 0.1,
 "c2_2": 0.1,
 "c2_3": 0.1,
 "c2_4": 28.6,
 "c2_5": 1.6,
 "c2_6": 2.4,
 "c2_7": 4.8,
 "c2_8": 2.5,
 "c2_9": 0.1,
 "memo": 0.0,
 "paper": 1.5
}
```