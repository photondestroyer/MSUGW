# Processing test log (Task 1.6)

Mode: stage2. Counties: [np.int64(20055), np.int64(20069), np.int64(20175), np.int64(20185)]. Cells: 300 ({'cohort': 200, 'rainfed': 100}).

No feasibility conclusion is drawn from this sample (pilot Section 1).

## Checks

| check | passes | detail |
|---|---|---|
| coverage rows only for cells in scope | True | 300 cells |
| ET rows only for cohort or rainfed cells | True |  |
| no cohort ET in dry months (gate) | True |  |
| every ET row has a county and a dry flag | True |  |
| one row per cell-month | True |  |
| ensemble ET within 0-400 mm/month | True | min 14.224066390041493, max 223.0154711673699 |
| dry flags cover every study county-month | True |  |
| blinding audit of the whole output root | True | 28 ET tables; [] |

## Seconds per task (last completed run)

```
{
 "t0_1": 47.2,
 "t0_4": 49.3,
 "counties": 7.5,
 "t0_2": 490.9,
 "t1_1": 2.9,
 "t1_2": 28.9,
 "t1_3": 13.7,
 "dry": 0.1,
 "t1_5": 8.1,
 "sample": 0.5,
 "t1_4": 469.3,
 "et_cells": 829.6
}
```