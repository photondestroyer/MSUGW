# Plan: audit and correction of the residual-ET training (XGBoost + differentiable LSTM-HBV)

## Context

Research Plan §6.1 defines the buffering signal as `Bi(t) = Fi(t) − F̂i(P, SMrz, VPD, Rn, T, PFT, season)`,
where `Fi` is the ET **anomaly**. Two models estimate `F̂`:

- XGBoost: `XGB model/gw_buffering_xgb_et.py` (+ `buffering_monthly_maps_timeseries.py`), TBI §11.15–11.16.
- Differentiable LSTM-parameterised HBV, full aquifer: `G:\MSU_GWB\delta_gwb\{data_full,train_full,buffering_full,hbv,dpl,config}.py`
  (outside this repo, not under git), TBI §11.18–11.19.

You asked me to check seven bug classes, fix them, look for more, and write the architecture and changes to a
new `.md`. I read the three documents, both codebases, the run logs, the data audit and the file headers.
Every class you listed is present. The headline numbers in TBI §11.15 and §11.19 (holdout R² 0.654,
τ = +1.62, "null placebo", the LSTM's flat WTD gradient) do not survive the audit.

Neither script can run today: the Earth Engine files they read (`GRIDMET_Merged_Ogallala.nc`,
`DROUGHT_Merged_Ogallala.nc`) are gone, and the lab files that replace them use `lat/lon` dims, the native
287×233 grid and packed int16.

**Your decisions (2026-10-07)**
- Forcing: `GRIDMET_Ogallala_lab_latest.nc`. I checked it: complete, 1979-01-01 → 2026-10-05. The zarr fallback is not needed.
- SPI: `GRIDMET_DROUGHT_Ogallala_lab_colab.nc`. I checked it: complete, has `spi90d`, 73 pentads in every year including 2020. The lab's SPI-90 is used as stored; no recomputation as the primary index.
- Split: 5-fold rotation over half-year blocks; each fold 70 train : 20 test : 10 validation.
- Run scope: fix and test only. No full training. The `.md` carries run instructions.
- Drought conditions: taken from published definitions, several indices, with checks that each definition is used correctly.

## Findings

Status: **R** = reproduced numerically, **F** = read in a file header or log, **C** = read in the code, not run.

| # | Class | Where | Bug | Status |
|---|---|---|---|---|
| L1 | Leakage / SPI | XGB l.375 | Old drought file has zero steps in 2020. "Nearest step" lookup gives 19 dekads of 2020 the 2019-12-31 value and 17 the **2021-01-05** value; 8 of 15 growing-season dekads use a future value. 2020 was then picked as the "driest holdout year" from that index. | R |
| L2 | Leakage | XGB l.456–458 | `Bi` is predicted on all rows: train-year residuals are in-sample, holdout residuals out-of-sample, and every effect pools both. Placebo years 2017 and 2019 are both training years, so the "null placebo" is partly in-sample shrinkage. | C |
| L3 | Leakage / anomaly | XGB l.481–486 | Dry-down climatology is `groupby("doyk")` over **all pixels and all years**: a spatial anomaly, not a per-pixel temporal one, and it includes holdout years. | C |
| L4 | Leakage | `data_full.py` l.413 | The LSTM's static inputs include **WTD** and irrigation (`irrig_frac`, `gfsad`). The exposure variable is inside the baseline model, so the residual cannot show a WTD gradient (reported τ ≈ +0.02). | C |
| L5 | Leakage | `train_full.py` l.359 | "Static" θ is re-inferred from whichever window is simulated: 730 days in training, 1,096 days of test-period forcing at test. | C |
| L6 | Leakage | `data_full.py` l.327 | Land-cover mode taken over 2001–2024, which includes validation and test years. | C |
| L7 | Leakage | both | No embargo between adjacent periods of different splits, although features reach back 90 days. | C |
| S1 | SPI | docs, both | `spi90d` is treated as SPI-3M. It is a rank-based index on a 1981–2016 baseline with about 36 distinct values, bounded at ±2.09, posted every 5 days. In May–Sep 2012, 22 % of aquifer values sit exactly at −2.09. | R |
| S2 | SPI | both | Threshold −0.8 (XGB) vs −1.0 (LSTM), neither cited; nearest step (XGB) vs last preceding step (LSTM). | C |
| S3 | SPI | lab file, catalogue | The lab's category layer cuts at −0.7 / −1.2 / −1.5; the US Drought Monitor table cuts at −0.8 / −1.3 / −1.6. The stored levels −0.71 and −1.28 land in different classes. The catalogue says EDDI shares SPI's sign; it does not (corr. with SPEI −0.80; positive = dry). | R |
| C1 | Crop | GFSAD file, XGB l.240 | File legend is invented. Real classes: 1–2 irrigated cropland, 3–5 rainfed cropland. `crop = isin([1,2])` labels irrigated cropland "rainfed-crop" and sends real rainfed cropland (35 % of the box) to "natural". Product is nominal 2010, labelled 2019. | F |
| C2 | Crop | XGB l.229, `data_full.py` l.352 | Irrigation map: 0 none, 1 low-to-medium, 2 high. XGB counts only class 1 as irrigated, so **high-irrigation cells are non-irrigated**; the LSTM counts `>0`, so 99 % of dry cells are irrigated and rainfed n = 0. Logged shares: irrigated 0.78, rainfed crop 0.02. | F |
| C3 | Crop | both | Each 4 km cell takes one nearest 500 m / 1 km pixel. `crop_frac` is the share of *years* that pixel was cropland, not the cropland share of the cell. | C |
| C4 | Crop | `data_full.py` l.414 | Class codes fed as numbers (`gfsad/5`, `hsg/4`). | C |
| T1 | Split | both | XGB 3:1:2 years (50:17:33). LSTM 5:1:2 years nominal, less in test because SSEBop ends May 2022. Neither is stratified by a climate index. | C |
| A1 | Anomaly | XGB | Target is raw ET with day-of-year features; R² is mostly the seasonal cycle and was never compared with a climatology baseline. | C |
| A2 | Anomaly | `train_full.py` l.206–219 | Soil-moisture loss: observations standardised by train statistics, model output by its **own window** mean and SD, so interannual soil-moisture anomalies are invisible to the loss. Soil moisture is not deseasonalised. | C |
| A3 | Anomaly | `data_full.py` l.444 | ET anomaly SD floor is 1 mm, equal to SSEBop's integer step, from 5 years of data; winter dekads dominate the loss. | C |
| A4 | Anomaly | `train_full.py` l.296 | Trained on anomaly loss, checkpoint chosen on raw RMSE (picked epoch 1). | F |
| R1 | Rainfall | XGB, `data_full.py` | Rainfall enters as raw mm; climatologies come from 5–6 years although gridMET starts in 1979. Daily rainfall is z-scored against a 5-year smoothed mean. Flagged source days (`srad` = 0 on 2014-12-22, `pr` 2013-12-30 = 12-31) are not masked. | C |
| B1 | Batch | XGB l.342–358 | Antecedent days are set to NaN before the cumulative sum. **P90 is wrong in 63 of 216 dekads** (on average 47 % of the true value), P30 in 21. | R |
| B2 | Batch | XGB l.371 | `SMrz_lag1` is NaN for the first dekad of each batch; 7 dekads are dropped, 3 of them in the growing season. | R |
| B3 | Batch | `train_full.py` l.121–133 | `--resume` reloads the **best** checkpoint (epoch 1) and no optimiser state, so each continuation restarted from epoch 1. Train loss stayed at 2.0; anomaly RMSE 1.68–1.93 where climatology scores 1.0. | F |
| B4 | Batch | `train_full.py` | 32 optimiser steps per epoch, one random window shared by a micro-batch, under 10 % of cells per epoch. | C |
| B5 | Batch | both | The monthly script duplicates the panel builder with B1/B2. XGB uses every second pixel, not the full aquifer. Both read files that no longer exist; XGB opens with `mask_and_scale=False`, which on the packed lab files returns raw integers. | F |
| U1 | UQ | both | Only a normal-approximation CI over pixel-dekads; nothing on data processing, climatology, split or prediction. | C |
| U2 | UQ | both | ET is point-sampled to 4 km (bilinear in XGB, nearest in LSTM) with no spread; the Earth Engine exports of SSEBop, land cover, irrigation and WTD skipped 22 % of native columns. | F |
| U3 | UQ | XGB l.757 | Tower ET sums whatever half-hours exist and keeps `LE` up to 14,980 W m⁻² and leaked −9999. | C |
| U4 | Analysis | XGB l.567–582 | Volume check compares residual ET summed over six years of dry dekads with depletion accumulated since ~1950. "Two orders of magnitude" is a time-base mismatch. | C |

## Approach

The two models disagree on data handling in ways that are themselves bugs (S2, C2, U2). The shared steps go
into one package that both models import.

**New package `residual_et/` in this repo.** Old scripts and `delta_gwb/` stay untouched, so earlier results remain reproducible.

| Module | Content | Fixes |
|---|---|---|
| `config.py` | Paths, windows, thresholds with their citations, seeds, baselines; one place | S2, B5 |
| `sources.py` | One adapter per dataset. Decodes packed integers, applies the lab `_qc` flags, maps `lat/lon`. Raises on a missing variable, a gap in a time axis, or a period not covered. | L1, R1, B5 |
| `grid.py` | Analysis grid = native gridMET 1/24° (287×233); aquifer cells = `roi_fraction ≥ 0.5` (27,333); area aggregation returning mean, SD and count per cell | U2 |
| `landcover.py` | Correct GFSAD, irrigation and IGBP legends; per-cell area fractions; year-matched MCD12Q1; management classes incl. "mixed/uncertain" (Plan §6.2) | C1–C4, L6 |
| `climate.py` | Per-cell dekad climatologies with standard errors; anomalies; standardised rainfall anomalies (gamma fit with a zero-rain mass) for 30 and 90 days ending the day before each dekad | L3, A1, A3, R1 |
| `drought.py` | Named, cited drought definitions over the lab indices; strict "last pentad before the dekad" lookup; definition checks | L1, S1–S3 |
| `features.py` | Year-chunked feature store with a carried antecedent buffer, so results do not depend on chunk size; all aquifer cells; written once as `.npy` + manifest and read by both models | B1, B2, B5 |
| `splits.py` | Half-year blocks (growing Apr–Sep, dormant Oct–Mar), stratified dry / normal / wet, 5 folds of 70:20:10, embargo set from the measured autocorrelation of ET anomalies | T1, L7 |
| `uq.py` | Conformal intervals, spatio-temporal block bootstrap, error budget | U1 |
| `train_xgb.py` | Cross-fitted model on ET anomalies; mean + 5/50/95 % quantiles; out-of-fold `Bi` for every row | L2, A1 |
| `train_dlstm.py`, `hbv.py`, `dpl.py` | Copied from `delta_gwb` and corrected: no exposure variables in static inputs, θ from a fixed climatology sequence, relative-saturation soil-moisture loss, checkpoint on the training metric, resume from the last state with optimiser, several windows per step, seed ensemble | L4, L5, A2, A4, B3, B4 |
| `buffering.py` | One dry-down rule and one set of strata for both models; within-pixel dry-minus-non-dry contrast; placebo on out-of-fold residuals; per-year volume check; screened tower validation; bootstrap CIs | S2, U3, U4 |
| `tests/` | See Verification | |

**Drought conditions from the literature** (all references confirmed to exist in OpenAlex)

| Definition | Rule | Source |
|---|---|---|
| Moderate meteorological drought | SPI-90 ≤ −1.0 | McKee, Doesken & Kleist 1993 |
| US Drought Monitor D1 or worse; classes D0–D4 | SPI or SPEI ≤ −0.8; class cuts −0.5 / −0.8 / −1.3 / −1.6 / −2.0 | Svoboda et al. 2002; USDM classification table (fetched) |
| Water-balance drought | SPEI-90 ≤ −1.0 | Vicente-Serrano et al. 2010 |
| High evaporative demand | EDDI-30 above its 80th percentile, positive = dry | Hobbins et al. 2016 |
| Flash drought | Root-zone soil moisture falls from ≥ 40th to ≤ 20th percentile within 20 days | Ford & Labosier 2017; Otkin et al. 2018 |
| Compound dry-down (Plan §6.1) | Precipitation deficit + soil-moisture decline + high VPD; your function below | Research Plan |

The lab indices are rank-based (Farahmand & AghaKouchak 2015) on the lab's grid (Abatzoglou 2013). The McKee,
EDDI and flash-drought thresholds above are from my memory of the papers; I will confirm each against a
readable source while implementing and mark in the `.md` any I could not confirm. Every result is reported
under each definition, so sensitivity to the definition is visible.

**Other design points**

- *Climatology baselines:* gridMET variables on 1985–2014; SSEBop ET on 2003–2014, which lies outside every split; SMAP from the training blocks of each fold (no earlier record exists).
- *Target:* ET anomaly in mm/dekad. Skill is reported against the per-cell climatology baseline.
- *Model inputs:* climate, soil moisture, soils, land-cover fractions, season. WTD and irrigation only stratify results, in both models.
- *Uncertainty, by step:* regridding (within-cell SD, count), observation (SSEBop vs screened towers; integer step), climatology (standard error), rainfall-anomaly fit (bootstrap), split (spread over folds), prediction (conformal quantiles for XGBoost; seed ensemble + conformal for the LSTM), drought definition (spread over the table above), effects (block bootstrap). Combined per row into σ(Bi).
- *Old 4 km layers:* PAW is re-aggregated from the 120 m gSSURGO stack; `dWL_*` (post-hoc screen only) is remapped by nearest cell and flagged.

**Your contributions (Learning mode).** Two short functions carry the judgement calls; I stub each with `TODO(human)` and ask one at a time:
1. `landcover.classify_management()`: thresholds turning area fractions into irrigated / rainfed crop / natural / mixed.
2. `buffering.is_drydown()`: how the cited index conditions, soil-moisture and VPD anomalies combine into the compound dry-down flag.

**Documentation.** New `RESIDUAL_ET_MODELS.md` at the repo root: architecture of both models, data flow, the
finding table with evidence, each change and its reason, drought definitions and their checks, split and
embargo, the uncertainty budget, assumptions, **how to run** (prepare → XGBoost → LSTM → buffering; the two
Python environments; resume; expected outputs), what was and was not run, and which TBI §11.15–11.19 numbers
are superseded. One line added to `README.md`. `TBI_DOCUMENTATION.md` has your uncommitted edits and is left alone.

## Files

- New: `residual_et/` (modules above), `RESIDUAL_ET_MODELS.md`.
- Edited: `README.md` (one line).
- Read only: `merged_datasets/**`, `G:\MSU_GWB\delta_gwb\**`, old `XGB model/*.py`.
- Reused: WTD block-mean logic (`delta_gwb/data_full.py::_wtd_cells`), HBV and ParamNet (`delta_gwb/hbv.py`, `dpl.py`), map helper (`gw_buffering_xgb_et.py::plot_grid`), retry-open (`data_tile.py::_open`).
- Environments: preparation and XGBoost in the system Python 3.13 (xgboost 3.4.1, scipy); LSTM in conda `SSM` (torch 2.5.1, RTX 3050 Ti 4 GB), which reads the prepared `.npy` store and needs no scipy.

## Verification (fix and test only)

1. Unit tests (`residual_et/tests/`):
   - features identical for chunk sizes 1, 7 and 32 dekads (B1, B2);
   - a rain impulse enters and leaves the 30- and 90-day windows on the right days; no feature uses a day on or after the dekad start, except the concurrent dekad means the Plan allows (L1);
   - packed-integer decoding matches the lab's scale and offset; the loader raises on a time gap such as the old 2020 hole (L1, B5);
   - drought definitions: lab `long_name` states the expected accumulation window; index date is strictly before the dekad; sign checks (SPI rises with rainfall, EDDI positive = dry); class cuts equal the cited table; lab category layer vs USDM classes reported; in-house 30/90-day rainfall anomaly rank-correlates with lab `spi30d`/`spi90d` and has mean 0, SD 1 on the baseline (S1–S3);
   - split: ratios within tolerance, every block tested once, no train row inside an embargo, each fold holds all three climate classes (T1, L7);
   - land cover: fractions sum to 1; irrigated share far below the old 0.78 (C1, C2);
   - no WTD or irrigation column in either model's input matrix (L4);
   - conformal intervals reach nominal coverage on synthetic held-out data (U1);
   - resume continues from the last epoch with optimiser state (B3).
2. Smoke runs on a small tile and short period: prepare → XGBoost (few rounds, 5 folds) → LSTM (2 epochs, resume once) → buffering. These prove the code path, not the science.
3. No full training. All result numbers in the `.md` are marked pending, with the commands to produce them.

## Open items I will not fix here

- No crop-type (CDL) or field-scale irrigation (LANID) layer exists in the repo; the Plan requires them. Documented as a gap.
- SSEBop, irrigation, GFSAD, HSG, MCD12Q1 and WTD files still carry the Earth Engine column subsampling; quantified in the uncertainty budget, not repaired.
- HBV has no irrigation input (TBI §11.18 verdict). A model-structure question, not a bug.
- The pilot-tile code (`train.py`, `buffering.py`, `data_tile.py`) is left as legacy.
