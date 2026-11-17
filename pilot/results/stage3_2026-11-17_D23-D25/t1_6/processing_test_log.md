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
| blinding audit of the whole output root | True | 344 ET tables; [] |

## Seconds per task (last completed run)

```
{
 "t0_1": 32.2,
 "t0_4": 46.6,
 "counties": 5.2,
 "t0_2": 172.1,
 "t1_1": 4.9,
 "t1_2": 23.0,
 "t1_3": 1.5,
 "dry": 0.0,
 "t1_5": 5.1,
 "sample": 0.4,
 "t1_4": 90.6,
 "et_cells": 245.6,
 "t1_6": 13.6,
 "t1_inv": 0.3,
 "t1_7": 51.4,
 "rep": 17.8,
 "c2_1": 0.1,
 "c2_2": 0.1,
 "c2_3": 1.7,
 "c2_4": 316.1,
 "c2_5": 3.9,
 "c2_6": 1162.1,
 "c2_7": 10.0,
 "c2_8": 12.0,
 "c2_9": 3.0,
 "memo": 0.0,
 "paper": 20.3,
 "fig_thick": 4.8,
 "t0_3": 2032.5
}
```