# Processing test log (Task 1.6)

Mode: stage3. Counties: [np.int64(20055), np.int64(20069), np.int64(20175), np.int64(20185)]. Cells: 300 ({'cohort': 200, 'rainfed': 100}).

No feasibility conclusion is drawn from this sample (pilot Section 1).

## Checks

| check | passes | detail |
|---|---|---|
| coverage rows only for cells in scope | True | 125809 cells |
| ET rows only for cohort or rainfed cells | True |  |
| no cohort ET in dry months (gate) | True |  |
| every ET row has a county and a dry flag | True |  |
| one row per cell-month | True |  |
| ensemble ET within 0-400 mm/month | True | min 4.670623145400594, max 229.44154057771664 |
| dry flags cover every study county-month | True |  |
| blinding audit of the whole output root | True | 262 ET tables; [] |

## Seconds per task (last completed run)

```
{
 "t0_1": 32.2,
 "t0_4": 46.6,
 "counties": 5.2,
 "t0_2": 330.8,
 "t1_1": 4.6,
 "t1_2": 27.4,
 "t1_3": 2.5,
 "dry": 0.1,
 "t1_5": 8.1,
 "sample": 0.5,
 "t1_4": 177.1,
 "et_cells": 201.8
}
```