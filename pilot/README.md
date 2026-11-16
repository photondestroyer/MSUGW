# pilot: the groundwater-ET pilot (Kansas, 2016-2024)

Code for `GW_ET_High_Plains_pilot.docx`. Design: `docs/superpowers/specs/2026-11-14-pilot-design.md`.
Definitions: `docs/pilot/DEFINITIONS.md`. Decision log: `pilot/decisions.yaml` (`docs/pilot/DECISIONS.md`).

## The physical idea, and where each piece is computed

| Idea | Code |
|---|---|
| Wells limit irrigation where the aquifer is thin -> thickness A with its error | Task 1.5 `phase1/t1_5_wells.py`, Check 2.4 `phase2/c2_4_thickness.py` |
| The shortfall shows in the peak-demand month -> peak month from weather only | Task 1.1 `phase1/t1_1_weather.py` |
| A drought signal irrigation cannot contaminate -> S over irrigation-free SMAP footprints | Tasks 1.2-1.3 `phase1/t1_2_footprints.py`, `t1_3_smap_dates.py` |
| Rainfed cells give the baseline of ET against soil moisture | Task 0.2 `phase0/t0_2_cells.py`, Check 2.5 `phase2/c2_5_sm_axis.py` |

## Rules, and how the code keeps them

1. **Blinding** - `blinding.py`. Dry flags come from weather and S only. `blinding.gate()` is the only way ET reaches
   disk: rainfed rows, and cohort rows of non-dry months; dry cohort ET is dropped in memory. `python -m pilot audit`
   re-reads every table and fails on a leak. No surface temperature is extracted.
2. **One ET release** - OpenET v2.1 only (`config.DATASETS["openet_tiles"]`, made by
   `Colab notebooks/OpenET_v2_1_Kansas_pilot.ipynb`).
3. **Decision log** - `decisions.yaml`. A task on real data (stage2, stage3) stops while any decision it uses is
   `proposed`. Values are part of each task's signature; every output's sidecar records the log's MD5 and the commit.
4. **Targets are reported, not enforced** - checks write observed value next to target; nothing stops on a miss.

## Commands

```
python -m pilot decisions                      # the log; then: decisions fix D01 D02 ... --by "<name>"
python -m pilot status                         # every task: done / to run / blocked by a proposed decision
python -m pilot run phase0                     # or a task id: run 1.1, run 2.4, run all
python -m pilot audit                          # blinding audit
python -m pilot test                           # definitions + unit + blinding (seconds); test e2e = smoke chain
python -m pilot preflight                      # which inputs exist in this mode
python -m pilot sync-inputs --dry-run          # what Colab needs on Drive
python -m pilot r-env --dry-run                # the isolated R environment for Task 0.3 (then --create)
```

Modes (`PILOT_RUN_MODE`): `test` (unit tests), `smoke` (one county, real data, a code test - never a result),
`stage2` (the processing-test counties of D27), `stage3` (all Kansas). Environment: `PILOT_RUN_ENV=local|colab`.

## How to run it (PowerShell, from `H:\MSU_GWB\datasets`)

The mode is an environment variable. Outputs go to `H:\MSU_GWB\pilot_out\` (`_smoke\`, `_stage2\`, or the root
for stage3), one folder per task and signature; a finished task is skipped on the next run.

```
# 1. code test on one county (Haskell). Runs with proposed decisions; its outputs are never results.
$env:PILOT_RUN_MODE = "smoke";  python -m pilot run all;  python -m pilot audit

# 2. confirm the decision log (once). Read it first: python -m pilot decisions
python -m pilot decisions fix all --by "<your name>"

# 3. Stage 2: the processing-test counties (D27), then read t1_6\...\processing_test_log.md
$env:PILOT_RUN_MODE = "stage2";  python -m pilot run all;  python -m pilot audit

# 4. Stage 3: every eligible Kansas cell
$env:PILOT_RUN_MODE = "stage3";  python -m pilot run all;  python -m pilot audit
```

`python -m pilot status` shows each task as done / to run / blocked by a proposed decision, and the last error.
`run all` leaves out Task 0.3 (the Mieno benchmark), which needs the R environment: `python -m pilot r-env`.
A run that stops (time budget, a Drive read error) is started again with the same command and continues.

## Order of work

`counties` -> 0.1, 0.4, 0.2 -> 1.1, 1.5, 1.2, 1.3, dry, sample -> 1.4, et_cells, 1.6, 1.7 -> 2.1-2.9 -> memo, paper, fig_thick.
`fig_thick` draws the ET version of Mieno et al. (2024) Figures 3-5 on simulated ET (real design, no ET value
read); the real figures are made after the freeze, as the blinding rule requires.
On disk: OpenET v2.1 Kansas (99 tile-years), LANID 1997-2020 and the LANID 2018-2020 CONUS maps.
Still waiting: the unmasked SMAP extraction (`Colab notebooks/SMAP_L3E_Kansas_pilot.ipynb`; until then Tasks 1.2-1.3
use a file that ends in June 2022 and holds data only inside the aquifer polygon, flagged `smap_source=fallback`);
WIMAS rates (Check 2.7 sheet).

Irrigation maps (D06): LANID is primary - it decides which cell-years are irrigated; AIM-HPA is the secondary
check (agreement in at least 3 of 2011-2015). Reference footprints (D12) are still labelled from AIM-HPA
2011-2015, because the LANID file holds labels only inside the aquifer polygon; they are checked against the
LANID 2018-2020 CONUS maps.

## Layout

Every task module starts with the pilot's own words, the physical idea, its decisions and deliverables, then
numbered `STEP` banners. Shared code: `definitions.py` (every term), `common/` (grids, cells, zonal statistics,
readers, statistics, report writer), `core/` (runner, provenance, I/O; copied from `paper1/core`).
