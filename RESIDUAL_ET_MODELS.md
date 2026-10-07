# Residual-ET models: audit, corrected architecture and how to run

**Scope.** The two models that estimate expected evapotranspiration (ET) for the buffering residual of
Research Plan §6.1, `Bi(t) = Fi(t) − F̂i(P, SMrz, VPD, Rn, T, PFT, season)`, where `Fi` is the ET anomaly:

- XGBoost (was `XGB model/gw_buffering_xgb_et.py`, TBI §11.15–11.16);
- a differentiable HBV water balance whose parameters come from an LSTM (was
  `G:\MSU_GWB\delta_gwb\train_full.py`, TBI §11.18–11.19).

**Date.** 2026-10-07. **Code.** `residual_et/` in this repository. The old scripts and `delta_gwb/` are
untouched.

## 1. Status

| | |
|---|---|
| Audit | Done. 29 findings, every one of the seven classes you listed is present (§2). |
| Corrected code | Written: one shared package used by both models (§3–§8). |
| Tests | 41 pass: 34 in the system Python (28 synthetic, 6 on the real files), 7 in conda `SSM` (§10). |
| Smoke runs | Prepare → XGBoost → LSTM (train, resume, infer) → buffering ran end to end on a 144-cell tile. They prove the code path only. |
| Tower comparison | Run on the real data (§8.3). It does not depend on the feature store. |
| Full training | **Not run**, as agreed. No scientific number in this file comes from the new code. Commands are in §9. |
| Waiting on you | `landcover.classify_management()` (thresholds for management classes). `buffering.is_drydown()` holds a provisional rule. |

**The numbers in TBI §11.15, §11.16 and §11.19 should not be used.** §2 explains why for each; the main ones:

| Old statement | Why it does not hold |
|---|---|
| XGBoost holdout R² 0.654 | Target was raw ET with day-of-year inputs; most of that is the seasonal cycle (A1). 27 % of dekads had a truncated 90-day rainfall input (B1). |
| "Holdout = two driest growing seasons, 2020 and 2021" | The drought file had no 2020. Its 2020 values were copied from 2019-12-31 and 2021-01-05 (L1). |
| Mean dry-down `Bi` +0.25 ± 0.11; τ = +1.62 | Residuals of training years were in-sample (L2); the dry-down rule used a spatial anomaly (L3); irrigation classes were wrong (C1, C2); the interval ignored clustering (U1). |
| Wet-year placebo −0.02 ("null, as required") | Both placebo years were training years (L2). |
| "Two orders of magnitude smaller than depletion, so no mined label" | Six years of dry dekads were compared with depletion since about 1950 (U4). |
| LSTM-HBV: flat water-table gradient, τ ≈ +0.02 | Water-table depth and irrigation were inputs of the parameter network (L4). |
| LSTM-HBV: test NSE 0.30, anomaly RMSE 1.93 | The network did not train: each resume restarted from epoch 1 (B3). |

## 2. Findings

Status: **R** reproduced numerically · **F** read in a file header or a log · **C** read in the code, not executed.
Line numbers refer to `XGB model/gw_buffering_xgb_et.py` (XGB) and `delta_gwb/` (LSTM).

### 2.1 Data leakage

| # | Where | Finding | Status | Fix |
|---|---|---|---|---|
| L1 | XGB l.375 | `DROUGHT_Merged_Ogallala.nc` had zero steps in 2020. The "nearest step" lookup gave 19 dekads of 2020 the 2019-12-31 value and 17 the 2021-01-05 value; 8 of 15 growing-season dekads used a future value. | R | `drought.lookup_preceding`: last pentad strictly before the dekad; raises on a gap. |
| L2 | XGB l.456–458 | `Bi` was predicted for all rows. On the smoke tile the same fold has RMSE 6.0 on rows it was fitted to and 11.0 on held-out rows. Placebo years 2017 and 2019 were training years. | C, R | 5-fold cross-fitting: every residual is out-of-fold (§5, §6.1). |
| L3 | XGB l.481–486 | Dry-down climatology was one mean per dekad over all pixels and all years, so climatologically dry pixels were always "below normal", and holdout years fed the baseline. | C | Per-cell baselines that end before the model window (§4.3). |
| L4 | `data_full.py` l.413 | Static inputs of the parameter network included water-table depth, `irrig_frac` and the GFSAD irrigation class. | C | `config.EXPOSURE_VARS` are barred from both models by `assert_no_exposure`. |
| L5 | `train_full.py` l.359 | "Static" parameters were re-inferred from the simulated window: 730 days in training, 1,096 days of test-period forcing at test. | C | Parameters come from fixed inputs (§6.2). |
| L6 | `data_full.py` l.327 | Land-cover mode over 2001–2024 includes validation and test years. | C | Land cover averaged over 2010–2014 only. |
| L7 | both | No embargo, although inputs reach back 90 days and ET anomalies are autocorrelated. | C | Embargo length from the measured autocorrelation (§5). |

### 2.2 SPI definition

| # | Finding | Status | Fix |
|---|---|---|---|
| S1 | `spi90d` was treated as "SPI-3M". It is a rank-based index on a 1981–2016 reference with 36 distinct values, bounded at ±2.09, posted every 5 days, "aggregated over the last 90 days". In May–Sep 2012, 22 % of aquifer values sit exactly at −2.09; on the smoke tile, 12 % of 2015–2022. | R | Described and checked at prepare time (§7). A continuous rainfall anomaly is computed in-house as a model input and cross-check. |
| S2 | Threshold −0.8 in XGBoost and −1.0 in the LSTM, neither cited; nearest step vs last preceding step. | C | One table of cited definitions for both models (§7). |
| S3 | The lab's `category` layer cuts at −0.7 / −1.2 / −1.5; the US Drought Monitor table cuts at −0.8 / −1.3 / −1.6. The stored levels −0.71 and −1.28 fall in different classes (4.4 % of values on the smoke tile). The catalogue says EDDI shares SPI's sign; its correlation with SPEI is −0.77 to −0.80, so positive EDDI is dry. | R | Classes are computed from the value with the published cuts; the difference is reported. |

Per your instruction the lab's stored `spi90d` is the drought index (`GRIDMET_DROUGHT_Ogallala_lab_colab.nc`:
complete, 73 pentads in every year including 2020).

### 2.3 Crop and irrigation classification

| # | Where | Finding | Status | Fix |
|---|---|---|---|---|
| C1 | GFSAD file; XGB l.240 | The file's legend is invented. Published classes: 0 non-cropland, 1–2 irrigated cropland, 3 rainfed cropland, 4–5 rainfed with minor / very minor fragments. `crop = isin([1, 2])` called irrigated cropland "rainfed-crop" and sent real rainfed cropland to "natural". The product is nominal 2010, labelled 2019. | F | Published legend in `landcover.py`. |
| C2 | XGB l.229; `data_full.py` l.352 | Irrigation map: 0 none, 1 low-to-medium, 2 high. XGBoost counted only class 1, so high-irrigation cells were non-irrigated. The LSTM counted every class above 0: 99 % of dry cells irrigated, rainfed n = 0. Logged shares were irrigated 0.78, rainfed crop 0.02. | F | Both classes kept as frequencies; fill is missing, not "not irrigated". |
| C3 | both | One nearest 500 m or 1 km pixel stood for a 4 km cell. `crop_frac` was the share of years that pixel was cropland. | C | Area fractions per cell. |
| C4 | `data_full.py` l.414 | Class codes fed as numbers (`gfsad/5`, `hsg/4`). | C | Fractions per class. |

### 2.4 Train / validation / test split

| # | Finding | Status | Fix |
|---|---|---|---|
| T1 | XGBoost 3:1:2 years (50:17:33). LSTM 5:1:2 years nominal, less in test because SSEBop ends in May 2022. No climate stratification. The XGBoost holdout was chosen from the broken 2020 index. | C | 5 folds of 70:20:10 over half-year blocks, stratified by drought class (§5). |

### 2.5 Seasonality and anomaly

| # | Where | Finding | Status | Fix |
|---|---|---|---|---|
| A1 | XGB | Target was raw ET with day-of-year inputs; skill was never compared with climatology. | C | Target is the anomaly; skill is reported against climatology. |
| A2 | `train_full.py` l.206–219 | Soil-moisture loss standardised the observations with training statistics and the model output with its own window, so a dry year looked like any other. Soil moisture was not deseasonalised. | C | Relative saturation with a fixed per-cell range (§6.2). |
| A3 | `data_full.py` l.444 | ET anomaly SD from 5 years with a 1 mm floor, the size of SSEBop's integer step; winter dekads dominated the loss. | C | 12-year baseline, 3 mm floor. |
| A4 | `train_full.py` l.296 | Trained on the anomaly loss, checkpoint chosen on raw RMSE (it picked epoch 1). | F | Selection on the trained quantity. |

### 2.6 Rainfall and climatology

| # | Finding | Status | Fix |
|---|---|---|---|
| R1 | Rainfall entered as raw mm. Baselines were 5–6 years although gridMET starts in 1979. Daily rainfall was z-scored against a 5-year smoothed mean. Days the lab flags (duplicate fields on 2011-01-01 and 2013-12-31, `srad` = 0 on 2014-12-22) were used as data. | C, F | 30-year baseline 1985–2014; standardised 30- and 90-day rainfall; lab flags applied (§4). |

### 2.7 Batch processing

| # | Where | Finding | Status | Fix |
|---|---|---|---|---|
| B1 | XGB l.342–358 | Look-back days were set to NaN before the cumulative sum. P90 was wrong in 63 of 216 dekads (on average 47 % of the true value), P30 in 21. | R | Each chunk reads its own look-back; result is independent of chunk size. |
| B2 | XGB l.371 | `SMrz_lag1` was NaN for the first dekad of every batch; 7 dekads dropped (2016-01-01, 2016-11-21, 2017-10-11, 2018-09-01, 2019-07-21, 2020-06-11, 2021-05-01). | R | Lag is a shift of the finished series; a missing input no longer drops a row. |
| B3 | `train_full.py` l.121–133 | `--resume` reloaded the best checkpoint (epoch 1) and no optimiser state. Train loss stayed at 2.0; anomaly RMSE 1.68–1.93 where climatology scores 1.0. | F | `last.pt` holds weights, optimiser and random state. |
| B4 | `train_full.py` | 32 optimiser steps per epoch, one window per micro-batch, under 10 % of cells per epoch. | C | Every cell every epoch; 4 windows per step. |
| B5 | both | The monthly script copied the panel builder with B1 and B2. XGBoost used every second pixel. Both read files that no longer exist; XGBoost opened gridMET with `mask_and_scale=False`, which on the packed lab files returns raw integers. | F | One builder, all aquifer cells, explicit decoding. |

### 2.8 Uncertainty and analysis

| # | Finding | Status | Fix |
|---|---|---|---|
| U1 | The only interval was 1.96·SD/√n over pixel-dekads. | C | Error budget per step and clustered intervals (§8). |
| U2 | ET was point-sampled to 4 km (bilinear in XGBoost, nearest in the LSTM). The Earth Engine exports of SSEBop, land cover, irrigation, soils and WTD skipped 22 % of native columns. | F | Area mean with within-cell SD and pixel count. The skipped columns cannot be restored here. |
| U3 | Tower ET summed whatever half-hours existed and kept `LE` up to 14,980 W m⁻² and leaked −9999. | C | Screened, gap-aware comparison (§8.3). |
| U4 | The volume check compared six years of dry-dekad residuals with depletion since about 1950. | C | Both sides in mm per year. |

## 3. Architecture

```
merged_datasets/ (read-only)
   │
   ▼  residual_et.features            system Python 3.13 (scipy, netCDF4)
residual_et_store/<profile>/          .npy arrays + manifest.json, definition_checks.json, splits.json
   │
   ├─► residual_et.train_xgb          system Python (xgboost)   → <profile>/xgb/
   ├─► residual_et.train_dlstm        conda SSM (torch)         → <profile>/dlstm/
   │
   ▼  residual_et.buffering           system Python             → <model>/buffering.json
residual_et.towers                    independent of the store  → residual_et_store/towers/
```

| Module | Role |
|---|---|
| `config.py` | Paths, periods, thresholds with citations, seeds, model settings. The only place they are defined. |
| `sources.py` | Readers that decode explicitly, apply the lab's QC flags, and raise on a missing variable, a time gap, an unfinished file or values outside physical ranges. |
| `grid.py` | Analysis grid and area aggregation (mean, within-cell SD, pixel count; class fractions). |
| `landcover.py` | Published legends; cropland, irrigation, land-cover and soil-group fractions; management classes. |
| `climate.py` | Dekad calendar, per-cell climatologies with standard errors, anomalies, standardised rainfall. |
| `drought.py` | Cited drought definitions, the strict time lookup, the definition checks. |
| `features.py` | Prepare stage; resumable (`axes → et → met → smap → drought → static → splits`). |
| `splits.py` | Blocks, climate classes, folds, embargo. |
| `uq.py` | Conformal intervals, clustered bootstrap, error budget. |
| `train_xgb.py` | Cross-fitted XGBoost. |
| `hbv.py`, `dpl.py`, `train_dlstm.py` | Water balance, parameter network, training / resume / inference. |
| `buffering.py` | Dry-down statistics, strata, placebo, volume check. Same code for both models. |
| `towers.py` | SSEBop against flux towers. |
| `tests/` | `run.py` (runner), `test_core.py`, `test_data.py`, `test_torch.py`, `synth.py`. |

Two profiles exist. `full` is the aquifer. `smoke` is a 0.5° tile in south-west Kansas (144 cells) with a
10-year rainfall baseline; it exists to prove that the code runs.

## 4. Data processing

### 4.1 Sources and grid

| Input | File | Handling |
|---|---|---|
| Daily meteorology | `GRIDMET_Ogallala_lab_latest.nc` (1979-01-01 → 2026-10-05, complete) | Packed int16 decoded with the file's scale and offset; records the lab flags are set to NaN (see below). |
| Drought indices | `GRIDMET_DROUGHT_Ogallala_lab_colab.nc` (1980-01-05 → 2026-10-02, complete) | Value at the last pentad strictly before each dekad (lag 1–5 days). |
| ET | `MODIS_ET_SSEBop_Merged_Ogallala.nc` | Area mean of the 1 km pixels in each cell (median 15 pixels), with SD and count. |
| Root-zone soil moisture | `SPL4SMGP_Ogallala_FULL.nc` | Nearest 9 km EASE cell (the product is coarser than the grid); daily mean of the 3-hourly values. |
| Water-table depth | `HRES-WTD_2015_Ogallala.nc` | Area mean and SD of the 30 m pixels (about 26,000 per cell). Stratification only. |
| Soil water storage | `GSSURGO_Ogallala_all_attrs_120m.nc`, `rootznaws` | Area mean and SD after projecting the 120 m pixels to the grid. |
| Land cover | `MCD12Q1_Merged_Ogallala.nc`, `LC_Type1` | Area fractions of ten vegetation groups, mean of 2010–2014. |
| Cropland | `GFSAD1000_V1_2019_Ogallala.nc` | Fractions: irrigated (1–2), rainfed (3), fragments (4–5), non-cropland (0). Stratification only. |
| Irrigation | `Global_irrigation_Area_Merged_Ogallala.nc` | Share of 2001–2015 in class 1 and in class 2 at the nearest 9 km pixel. Stratification only. |
| Soil group | `Hydrologic_Soil_Group_250m_2019_Ogallala.nc` | Fractions A–D; dual codes 14/24/34 go to their drained group. |
| Water-level change | `dWL_*_4km.nc` | Nearest cell of the older stretched grid. Used only in the volume check. |

**Grid.** The lab's native gridMET grid, 287 × 233 at 1/24°. A cell is an aquifer cell when at least half of
it lies inside `hp_bound2010` (`roi_fraction ≥ 0.5`): 27,333 cells. The old "4 km grid" (287 × 181) was an Earth
Engine resampling with 22 % of columns dropped.

**Lab QC flags.** A record is dropped when the lab marks it `no_data`, when it repeats the previous record
without being a constant field, or when a field that cannot be constant is constant (`srad`, `tmmn`, `tmmx`,
`vpd`, `etr`, `pet`). A constant rainfall field is a dry day and is kept. In 1985–2022 this removes rainfall and
temperature on 2011-01-01 and 2013-12-31 and radiation on 2013-12-31 and 2014-12-22. Dekad sums are scaled to
the full dekad; a dekad needs 80 % of its days. The LSTM's daily forcing takes the day-of-year normal on those
days, because the water balance cannot step through a NaN. The identical 31 December minimum temperatures of
2017–2019 are not flagged by the lab and remain.

### 4.2 Time base

Dekads start on the 1st, 11th and 21st, as in SSEBop: 258 dekads from 2015-04-01 to 2022-05-21 (SMAP starts
2015-03-31; SSEBop ends with the 2022-05-21 dekad). The LSTM's daily forcing starts 2013-04-01, two years of
spin-up.

### 4.3 Seasonality and anomalies

Every anomaly is per cell and per dekad-of-year, pooled over the neighbouring dekad on each side for the SD.

| Variable | Baseline | Why |
|---|---|---|
| Rainfall, VPD, radiation, temperature, reference ET | 1985–2014 (30 years) | Ends before the model window, so no split can leak into it. |
| SSEBop ET | 2003–2014 (12 years) | Same reason; SSEBop starts in 2003. |
| SMAP soil moisture | Training dekads of each fold | SMAP has no earlier record. Computed inside each fold. |

The standard error of each climatological mean is stored (`*_clim_se`), using the number of years, because
pooled dekads of one year are not independent.

ET in 2015–2022 is not centred on the 2003–2014 normal (+1.6 mm/dekad on the smoke tile). A constant offset of
a cell cancels in the within-pixel contrast of §8.2, which is the statistic to prefer for that reason.

### 4.4 Rainfall and climatology

Three rainfall inputs reach the models, all relative to the 30-year baseline:

- `P_dek_anom`: dekad rainfall minus its normal.
- `SPI30_g`, `SPI90_g`: rainfall summed over the 30 and 90 days that end the day before the dekad, transformed
  to a standard-normal score with a gamma distribution plus a probability mass at zero, fitted per cell and
  dekad-of-year on the baseline (Thom's closed-form estimator; zero rainfall takes the centre of the zero mass,
  Stagge et al. 2015). Scores are continuous to ±3.7, so they do not saturate where the lab index does.
- Their bootstrap SD over resampled baseline years (`SPI30_se`, `SPI90_se`).

The lab `spi30d` / `spi90d` are not model inputs. They define drought conditions (§7) and the climate class of
the split (§5).

## 5. Split

**Unit.** A half-year block: growing season April–September, dormant season October–March. A dry-down and
its recovery stay in one block. All cells of a block share its role, so no split holds a neighbouring pixel
of the same date. 15 blocks: G2015 … G2022 (the last one has 6 dekads) and D2015 … D2021.

**Climate index.** Aquifer-mean lab SPI-90 averaged over the block. Blocks are ranked within their season
and cut into three classes of five: dry, normal, wet.

**Folds.** Each of the five folds tests one dry, one normal and one wet block (20 %), validates on three
half-blocks, one per class (10 %), and trains on the rest (70 %). Every block is tested exactly once.

**Embargo.** Training dekads next to a validation or test dekad are removed. The length is the first lag at
which the median autocorrelation of the ET anomaly falls below 0.2, limited to 2–6 dekads.

Realised on the smoke tile (the full aquifer will give different classes, since the index is its own mean):

| Fold | Train | Embargoed | Validation | Test | Test blocks (class) |
|---|---|---|---|---|---|
| 0 | 0.647 | 0.085 | 0.105 | 0.163 | G2016 (wet), D2020 (normal), G2022 (dry) |
| 1 | 0.593 | 0.093 | 0.105 | 0.209 | D2016 (normal), D2018 (wet), G2020 (dry) |
| 2 | 0.632 | 0.054 | 0.105 | 0.209 | G2015 (wet), D2017 (normal), D2021 (dry) |
| 3 | 0.609 | 0.078 | 0.105 | 0.209 | D2015 (wet), G2019 (normal), D2019 (dry) |
| 4 | 0.632 | 0.078 | 0.081 | 0.209 | G2017 (wet), G2018 (dry), G2021 (normal) |

Train plus embargoed is 0.69–0.73. The ET-anomaly autocorrelation there is 0.25, 0.18, 0.19, 0.06 at lags
1–4, giving a 2-dekad embargo. `splits.json` in the store holds the realised split of any run.

XGBoost runs all five folds, so every dekad has a residual. The LSTM is trained per fold; one fold gives
residuals for its three test blocks.

## 6. Models

Both models see the same kinds of input: climate, soil moisture, soils, land cover, season. Neither sees
water-table depth, irrigation, cropland-irrigation class, water-level change, or the ET climatology. The
last one is excluded because a high summer normal identifies an irrigated cell.

### 6.1 XGBoost

| | |
|---|---|
| Target | ET anomaly, mm per dekad |
| Rows | One per cell and dekad with an observed target. Missing inputs stay NaN. |
| Inputs (25) | `P_dek_anom`, `SPI30_g`, `SPI90_g`, `VPD_anom`, `SRAD_anom`, `TMEAN_anom`, `ETR_anom`, `SMrz_anom`, `SMrz_lag1_anom`; `doy_sin`, `doy_cos`; `PAW`, `aridity`, `P_annual`, `T_annual`; ten land-cover fractions |
| Mean model | `reg:squarederror`, depth 6, learning rate 0.05, subsample 0.8, column sample 0.7, L2 1.0, min child weight 10, `hist`, 256 bins, up to 800 rounds, early stopping 50 on the validation dekads |
| Interval model | `reg:quantileerror` at 0.05 / 0.50 / 0.95, same settings, then widened on the validation dekads so that 90 % of them are covered (conformalised quantile regression, Romano et al. 2019) |
| Batching | Rows are streamed 12 dekads at a time into a `QuantileDMatrix`; the full table is never built |
| Output | `oof_pred`, `oof_q50`, `oof_lo`, `oof_hi`, `oof_fold`, `Bi`, `sigma_*`, `metrics.json`, one model file per fold |

`metrics.json` reports, per fold and pooled: RMSE, the RMSE of climatology (anomaly = 0), skill against
climatology, anomaly correlation, 90 % coverage, interval width, and the in-sample fit on training dekads
for comparison.

### 6.2 Differentiable LSTM-HBV

**Water balance (`hbv.py`).** Unchanged from the pilot: daily snow, soil, upper and lower subsurface stores;
ET = `min(Ep · min((Ss/(FC·LP))^γ, 1), Ss)`; mass-conservative (tested). It has no groundwater-to-ET pathway
and no irrigation input. Fifteen parameters with the pilot's literature bounds (`config.HBV_BOUNDS`).

**Parameter network (`dpl.py`).** LSTM, 2 layers of 256 units, 21 inputs, linear head to 15 parameters,
sigmoid into the bounds; 815,887 weights. Its inputs are fixed per cell:

- 18 static attributes, standardised: `PAW`, `aridity`, `P_annual`, `T_annual`, four soil-group fractions,
  ten land-cover fractions;
- a 366-step sequence of the cell's day-of-year normals of rainfall, temperature and reference ET
  (1985–2014, 31-day smoothing).

So a cell has one parameter set for all periods. This departs from Feng et al. (2022), who feed the forcing
of the simulated window; that made the parameters depend on the test period here.

**Loss.** `RMSE(z_ET) + 0.25 · RMSE(relative saturation)`

- `z_ET = (model dekad ET − observed dekad ET) / max(SD_clim, 3 mm)`, SD from the 2003–2014 baseline.
- Relative saturation: model `Ss/FC` against `(SM − p1)/(p99 − p1)` clipped to [0, 1], where p1 and p99 are the
  cell's 1st and 99th soil-moisture percentiles on training days. A dry year stays visible.

**Training.** Windows of 730 days, loss on the last 365, only on training dekads and training days. Each
optimiser step uses 4 windows of 64 cells; an epoch visits every cell once. Adam, learning rate 1e-3,
gradient clipping at 1. Water balance on the CPU, network on the GPU.

**Validation and selection.** 512 fixed cells, one continuous simulation from 2013-04-01, the same loss
expression on validation dekads. `best.pt` is the epoch with the lowest value.

**Resume.** `last.pt` stores weights, optimiser state, epoch, best value and the random state. A two-epoch
run and a one-epoch run resumed for one epoch give identical weights (tested).

**Inference.** Continuous simulation of every cell; expected anomaly = dekad ET sum − normal. A
split-conformal half-width from validation residuals gives the 90 % interval; several members
(`--member`) add their spread. Residuals are written for the fold's test dekads only.

## 7. Drought conditions and definition checks

| Name in `config.DROUGHT_DEFS` | Rule | Source | Threshold checked against |
|---|---|---|---|
| `mckee_moderate` | SPI-90 ≤ −1.0 | McKee, Doesken & Kleist (1993) | My recollection of the paper; not re-read |
| `usdm_d1` | SPI-90 ≤ −0.8 | Svoboda et al. (2002); US Drought Monitor classification table | The table on droughtmonitor.unl.edu, fetched 2026-10-07 |
| `spei_moderate` | SPEI-90 ≤ −1.0 | Vicente-Serrano, Beguería & López-Moreno (2010) | Recollection; not re-read |
| `eddi_high` | EDDI-30 ≥ 0.84 (80th percentile) | Hobbins et al. (2016) | Sign verified on the data; the 80 % bound is a recollection |
| `flash_sm` | Soil-moisture percentile ≥ 40 then ≤ 20 within two dekads | Ford & Labosier (2017); Otkin et al. (2018) | Recollection; not re-read |
| `compound_drydown` | `buffering.is_drydown()` | Research Plan §6.1 | Provisional (below) |

US Drought Monitor classes used for labels: D0 ≤ −0.5, D1 ≤ −0.8, D2 ≤ −1.3, D3 ≤ −1.6, D4 ≤ −2.0. All seven
references were confirmed to exist in OpenAlex. Please check the four recalled thresholds against the papers
before publication.

`compound_drydown` currently reads: SPI-90 at D1 or worse, soil moisture below its normal, VPD above its normal.
It is the old XGBoost rule with the cited D1 bound. Every statistic is reported under all six conditions, so
the dependence on the definition is visible.

**Checks.** Hard checks stop the prepare stage; the rest is written to `definition_checks.json`.

| Check | Kind | Result on the smoke tile |
|---|---|---|
| The file's `long_name` states the window each rule assumes ("aggregated over the last N days") | hard | passes for all five indices |
| Each rule's index has the window the rule names | hard | passes |
| Class cuts equal the published US Drought Monitor table | hard | passes |
| Index date strictly before the dekad; no gap over 6 days | hard | lag 1–5 days |
| Signs: SPI rises with rainfall; SPEI follows SPI; EDDI rises with reference ET and opposes SPEI | hard | +0.89, +0.88, +0.36, −0.77 |
| Lab `spi30d` / `spi90d` rank-correlate with the in-house 30 / 90-day score (≥ 0.7 required) | hard | 0.90, 0.94 |
| Distinct levels, bounds, share at the bounds | report | 36 levels, ±2.09, 12.4 % at the minimum |
| Lab category cuts vs US Drought Monitor cuts | report | 4.4 % classified differently (levels −1.28, −0.71) |
| The old drought time axis is refused | test | passes: "no drought index within 6 days before 2020-01-11" |

## 8. Uncertainty

### 8.1 By step

| Step | Quantity | Where |
|---|---|---|
| Regridding ET | Within-cell SD and pixel count; SE = SD/√n | `ET_sd`, `ET_n` |
| SSEBop storage | Whole millimetres: rounding SD 1/√12 mm per pixel | `uq.observation_sigma` |
| Regridding statics | Within-cell SD of WTD and PAW | `WTD_sd`, `PAW_sd` |
| Climatology | Standard error of each normal | `*_clim_se` |
| Rainfall score | Bootstrap SD over baseline years | `SPI30_se`, `SPI90_se` |
| Split | Spread of skill over the five folds | `metrics.json` → `fold_spread` |
| Prediction | Conformal 90 % interval; member spread for the LSTM | `oof_lo`, `oof_hi`, `sigma_ensemble` |
| Drought definition | Each statistic under six conditions | `buffering.json` → `definitions` |
| Effects | Two-way cluster bootstrap: 0.5° spatial blocks and half-year time blocks | `uq.block_bootstrap` |
| Observation | SSEBop against towers | §8.3 |

**Per-row σ(Bi).** `sigma_Bi = sqrt(sigma_model² + sigma_obs² + sigma_clim²)`, where `sigma_model` is the SD
of a normal with the width of the conformal interval. This is conservative: the interval is calibrated on
observed anomalies, so it already contains the errors of the target. `sigma_model` alone is the calibrated
term. The components are assumed independent, which is not strictly true in sparse or cloudy scenes.

**Not quantified.** The SSEBop retrieval error per pixel (no such layer exists; §8.3 bounds it), the SMAP L4
error (the file carries none), and the effect of the columns Earth Engine skipped.

### 8.2 Effect estimates

For each condition, on growing-season (May–September) dekads with an out-of-fold residual:

- mean `Bi` with a clustered 95 % interval, and the old-style interval beside it for comparison;
- within-pixel contrast: mean over cells of (mean `Bi` in dry dekads − mean `Bi` in the cell's other dekads);
- by WTD tercile and by management class;
- τ = shallow − deep, overall and within each management class;
- share of dry-down residuals above 1.645 σ(Bi).

Plus the wet placebo (SPI-90 ≥ +0.8), on out-of-fold residuals like everything else, and the volume check in
mm per year: residual ET of 2017–2019 against the 2017→2019 water-level change × Sy / 2, with the long-term
rate (change since about 1950 × Sy / 69) for scale.

An interval is reported only when the data span at least three spatial and three time blocks. On the smoke
tile the clustered half-widths are about ten times the old ones (for example ±7 against ±0.64 mm/dekad).

### 8.3 SSEBop against flux towers (run on the real data)

`python -m residual_et.towers`. Latent heat flux outside −100 … 800 W m⁻² and sentinels are removed; a dekad
needs 80 % valid records; ET = mean flux × dekad length / 2.45 MJ kg⁻¹; SSEBop is the mean of the 3 × 3
pixels at the tower. 21 towers inside the SSEBop box, 2003–2021.

| | n (site-dekads) | SSEBop − tower | RMSE | Unbiased RMSE | r |
|---|---|---|---|---|---|
| All dekads | 2,797 | +1.4 | 12.0 | 11.9 | 0.74 |
| May–September | 1,205 | +6.2 | 15.9 | 14.6 | 0.56 |

Units mm per dekad. The old figure was RMSE 16.06 over 1,820 dekads without screening. These include the
mismatch between a tower footprint and 9 km², and tower energy-balance closure, so they are an upper bound on
the SSEBop error. Most of the towers lie outside the aquifer polygon.

## 9. How to run

Two environments. Everything is run from the repository root.

```bash
# system Python 3.13: numpy, scipy, netCDF4, pyproj, xgboost, pandas
python -m residual_et.tests.run core data
```
```bash
# conda env SSM: numpy, torch
C:/Users/AlienX/anaconda3/envs/SSM/python.exe -m residual_et.tests.run torch
```

**1. Prepare the store** (once; resumable; `--stage <name> --force` redoes one stage).

```bash
python -m residual_et.features --profile full
```

Then read `residual_et_store/full/manifest.json`, `definition_checks.json` and `splits.json`. The smoke tile
took 78 s. For the aquifer I estimate 15–30 minutes, mostly reading, and a store of 2–3 GB; the full profile
has not been run, so both are estimates.

**2. XGBoost, all five folds.**

```bash
python -m residual_et.train_xgb --profile full
```

Output in `residual_et_store/full/xgb/`. About 7 million rows; the smoke tile took 16 s.

**3. LSTM-HBV, one fold and member at a time.**

```bash
C:/Users/AlienX/anaconda3/envs/SSM/python.exe -m residual_et.train_dlstm --profile full --fold 0
```
```bash
C:/Users/AlienX/anaconda3/envs/SSM/python.exe -m residual_et.train_dlstm --profile full --fold 0 --resume --epochs 10
```
```bash
C:/Users/AlienX/anaconda3/envs/SSM/python.exe -m residual_et.train_dlstm --profile full --fold 0 --infer --members 0
```

`--epochs` is the number of additional epochs (default 30). `--member 1`, `--member 2` train further ensemble
members; list them in `--members` at inference. On the smoke tile a window of 64 cells took about 2 s, which
puts an epoch of the full aquifer near 15 minutes and 30 epochs near 7–8 hours per fold and member. Watch
`train_log.csv`: `train_loss` and `val_loss` should fall. A flat `train_et` near 1.5–2 means the network is
not learning.

**4. Buffering statistics** (after `classify_management()` is written).

```bash
python -m residual_et.buffering --profile full --model xgb
```
```bash
python -m residual_et.buffering --profile full --model dlstm
```

**5. Tower comparison** (already run; rerun if the tower file changes).

```bash
python -m residual_et.towers
```

`--profile smoke` on any command gives the quick tile. Set `RESIDUAL_ET_STORE` to move the store.

## 10. What was verified

| Test group | Count | What it shows |
|---|---|---|
| Batch processing | 3 | Inputs identical for chunk sizes 1, 7, 32, 1000; antecedent rainfall equals a brute-force sum; the old cumulative sum is truncated |
| Leakage in time | 2 | A rain spike enters and leaves the 30- and 90-day windows on the right days; the lookup is strictly before the dekad and refuses a missing year |
| Sources | 2 | Packed integers decode; lab flags are respected (a dry day is kept, a duplicate dropped); gaps and unfinished files are refused |
| Definitions | 4 | US Drought Monitor classes; each rule names an index of its own window and a source; a flipped index is caught; flash drought needs a fast fall |
| Rainfall, climatology | 3 | The rainfall score has mean 0 and SD 1 on its baseline and is monotonic; climatology is per cell with a correct standard error; dekad calendar |
| Split | 4 | 70:20:10 within tolerance; every block tested once; every fold tests and validates all three classes; no training dekad inside the embargo; embargo follows the data's memory |
| Regridding, land cover | 2 | Mean, SD and count are exact; fractions sum to 1 with the published legend |
| Model inputs | 1 | Exposure variables are refused |
| Uncertainty | 2 | Conformal interval reaches 90 %; the clustered interval is over ten times the old one on clustered data |
| Cross-fitting | 4 | Every residual is out-of-fold and from the right fold; a noise target gets no out-of-fold skill; the same target looks skilful in-sample; a missing lag drops no row |
| Towers | 1 | Screening and the gap rule |
| Real files | 6 | Lab drought windows and cadence; index signs in 2012; gridMET decodes to Kelvin; the old drought axis is refused; GFSAD classes; SSEBop cadence and integer storage |
| LSTM-HBV | 7 | Mass balance; no exposure input; parameters identical for every period; soil-moisture loss sees a dry year; loss uses only dekads of its role and ignores test targets; resume is exact; inference writes test residuals only |

**Smoke results (code path only; 144 cells, 10-year rainfall baseline).** XGBoost: out-of-fold RMSE 9.9
against 11.1 mm/dekad for climatology, 90 % interval covering 87.7 %. LSTM-HBV: 9 optimiser steps, so its
numbers mean nothing. One real bug surfaced only here: on the GPU the saved random state was loaded onto
CUDA and resume failed. It is fixed.

**Not verified.** The full profile has not been run, so memory and run time at 27,333 cells are estimates.
The LSTM has not been trained to convergence, so I do not know whether the corrected model beats climatology;
the pilot did not, and the missing irrigation input is still missing. `buffering.py` ran only with a stand-in
management rule. Nothing compares old and new residual maps, because the old inputs no longer exist.

## 11. Assumptions and open items

1. **Management classes.** `landcover.classify_management()` is unwritten. On the smoke tile GFSAD gives
   irrigated 0.17, rainfed cropland 0.48, fragments 0.24; the irrigation map is class 1 in 94 % of years and
   class 2 in under 1 %, so its class 1 does not discriminate there.
2. **GFSAD is nominal 2010 and the irrigation map ends in 2015.** The analysis runs to 2022.
3. **No crop-type or field-scale irrigation layer** (USDA CDL, LANID) is in the repository. The Research
   Plan requires them before any irrigation-subsidy label.
4. **Earth Engine column subsampling** remains in SSEBop, MCD12Q1, GFSAD, the irrigation map, soil groups
   and WTD. Coordinates are correct; 22 % of native columns are absent.
5. **HBV has no irrigation input.** The LSTM-HBV residual over irrigated land contains all irrigation water.
   A model-structure question, unchanged.
6. **Embargo and blocks.** Half-year blocks and a 2–6 dekad embargo are choices. Whole-year blocks would be
   cleaner in time but allow only about 71:14:14.
7. **ET baseline 2003–2014** assumes the normal is stable into 2015–2022. Use the within-pixel contrast.
8. **Lab SPI saturates** at −2.09 in the driest fifth of 2012 and in 2022; it cannot rank severity there.
9. **Soil-moisture loss.** Comparing `Ss/FC` with SMAP relative saturation assumes both span wilting to
   saturation. The 0.25 weight is unchanged from the pilot.
10. **Specific yield 0.15**, the user constant of TBI §11.3, in the volume check.
11. **Pilot-tile code** (`delta_gwb/train.py`, `buffering.py`, `data_tile.py`) is legacy and uncorrected.
12. **`MERGED_DATASETS_CATALOG.md`** still describes the removed Earth Engine gridMET and drought files.

## 12. References

- Abatzoglou, J. T. (2013). Development of gridded surface meteorological data for ecological applications and modelling. *Int. J. Climatol.* 33. doi:10.1002/joc.3413
- Farahmand, A. & AghaKouchak, A. (2015). A generalized framework for deriving nonparametric standardized drought indicators. *Adv. Water Resour.* 76. doi:10.1016/j.advwatres.2014.11.012
- Feng, D., Liu, J., Lawson, K. & Shen, C. (2022). Differentiable, learnable, regionalized process-based models. *Water Resour. Res.* 58. (not re-verified here; cited by the pilot code)
- Ford, T. W. & Labosier, C. F. (2017). Meteorological conditions associated with the onset of flash drought in the Eastern United States. *Agric. For. Meteorol.* 247. doi:10.1016/j.agrformet.2017.08.031
- Hobbins, M. T. et al. (2016). The Evaporative Demand Drought Index. Part I. *J. Hydrometeor.* 17. doi:10.1175/JHM-D-15-0121.1
- McKee, T. B., Doesken, N. J. & Kleist, J. (1993). The relationship of drought frequency and duration to time scales. *8th Conf. on Applied Climatology.*
- Nagaraj, D. et al. (2021). A new dataset of global irrigation areas from 2001 to 2015. *Adv. Water Resour.* 152, 103910.
- Otkin, J. A. et al. (2018). Flash droughts: a review and assessment of the challenges imposed by rapid-onset droughts in the United States. *BAMS* 99. doi:10.1175/BAMS-D-17-0149.1
- Romano, Y., Patterson, E. & Candès, E. (2019). Conformalized quantile regression. *NeurIPS.* (not re-verified here)
- Stagge, J. H. et al. (2015). Candidate distributions for climatological drought indices (SPI and SPEI). *Int. J. Climatol.* 35. (not re-verified here)
- Svoboda, M. et al. (2002). The Drought Monitor. *BAMS* 83. doi:10.1175/1520-0477-83.8.1181
- Vicente-Serrano, S. M., Beguería, S. & López-Moreno, J. I. (2010). A multiscalar drought index sensitive to global warming: the SPEI. *J. Climate* 23. doi:10.1175/2009JCLI2909.1
