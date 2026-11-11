# paper1 configuration reference

Generated from `config.py` by `python -m paper1 docs` - do not edit by hand (a test checks it is current).

Changing a setting listed under *stages* gives those stages (and every stage after them) a new signature,
so they recompute into a new folder; outputs made with the old value stay where they are.


## where and how

| setting | default | meaning | stages |
|---|---|---|---|
| `RUN_ENV` | `os.environ.get("P1_RUN_ENV", "colab")` | "colab" \| "local" | - |
| `RUN_MODE` | `os.environ.get("P1_RUN_MODE", "test")` | "test" \| "focus" \| "full" | - |
| `COLAB_MYDRIVE` | `"/content/drive/MyDrive"` |  | - |
| `COLAB_DATA_ROOT` | `f"{COLAB_MYDRIVE}/MSUGWB"` | where the dataset notebooks write | - |
| `COLAB_SCRATCH` | `"/content/p1_scratch"` | fast local disk of the Colab VM | - |
| `LOCAL_DATA_ROOT` | `os.environ.get("P1_LOCAL_DATA_ROOT", "G:/MSU_GWB/drive_mirror/MSUGWB")` |  | - |
| `LOCAL_SCRATCH` | `os.environ.get("P1_LOCAL_SCRATCH", "G:/MSU_GWB/p1_scratch")` |  | - |
| `REPO_ROOT` | `os.path.dirname(os.path.dirname(os.path.abspath(__file__)))` | G:/MSU_GWB/datasets locally | - |
| `LOCAL_MERGED` | `os.path.join(REPO_ROOT, "merged_datasets")` |  | - |
| `OUTPUT_FOLDER` | `"P1_rain_only"` | outputs: <data root>/P1_rain_only[/_test\|/_focus]/<stage>/<signature>/ | - |
| `SHARED_FOLDER` | `"P1_rain_only/_shared"` | downloads reused by every mode (NHDPlus HR subsets, counties, DEM tiles) | - |
| `SESSION_BUDGET_H` | `11.0` | a stage stops cleanly (and resumes next session) after this many hours | - |
| `LOCK_STALE_H` | `13.0` | a lock older than this is considered abandoned | - |
| `CODE_IN_SIGNATURE` | `True` | a stage's output folder changes when its code changes (safer); False freezes it | - |

## which tiles / years

| setting | default | meaning | stages |
|---|---|---|---|
| `TEST_TILES` | `["z14e02n29"]` | south-west Nebraska (Republican / Platte / Sandhills edge) | tiles |
| `TEST_YEARS` | `[2015, 2016]` | calendar years of 30 m data used in test mode (Landsat v2 on Drive starts in 2015) | blocks, classes |
| `TEST_WINDOW_PX` | `1024` | test mode reads a square of this many 30 m pixels (multiple of 16) | tiles |
| `TEST_WINDOW_CENTER` | `None` | (row, col) in frame pixels; None -> the frame centre | tiles |
| `TEST_SPACE_GROUP_DEG` | `0.05` | spatial CV group size in test mode (the window is ~30 km) | folds |
| `TEST_SPACE_FOLDS` | `2` |  | folds |
| `TEST_YEAR_GROUPS` | `[[2015, 2016]]` | one year group: no year hold-out in test mode | folds |
| `FOCUS_POINTS` | `{ "FL1_central_platte": (-98.9, 40.75), "FL2_sandhills": (-101.0, 4...` |  | tiles |
| `AGG_YEAR_FIRST` | `2007` | blocks are also aggregated for 2007 (pre-season of WY2008) | blocks |
| `DRIVER_MONTH_FIRST` | `"2005-01"` | monthly drivers kept from here (24-month history of WY2008) | drivers |

## source datasets (registry)

| setting | default | meaning | stages |
|---|---|---|---|
| `SOURCE_SUBFOLDER` | `""` | "" = the finished products; "/_test" would read the notebooks' test outputs | tiles |
| `ROI_LABEL` | `"Ogallala"` |  | - |
| `LANDSAT_TAG` | `"LANDSAT_C2L2_indices_30m_monthly_tiles_v2_latest"` | on Drive 2026-11-11: 2015-2026 (2008-2014 to come) | tiles, blocks |
| `OPENET_TAG` | `"OPENET_v2_0_monthly_30m_tiles"` | release 2.0 (1999-2024); the v2_1 folder is empty | blocks |
| `CDL_TAG` | `"USDA_NASS_CDL_30m_tiles"` |  | classes |
| `NWI_TAG` | `"NWI_wetlands_30m_tiles"` |  | static |
| `DATASETS` | `{ "landsat": dict(kind="tiles", tag=LANDSAT_TAG, per_year=True, var...` |  | - |
| `LOCATE_BY_NAME` | `True` | if a "file" is not at its registry path, search the data root (depth 3) by name | - |

## 30 m tile geometry (from the notebooks)

| setting | default | meaning | stages |
|---|---|---|---|
| `PIXEL_M` | `30.0` |  | tiles |
| `TILE_PX` | `5120` |  | tiles |
| `LATTICE_ORIGIN_M` | `15.0` | Landsat pixel edges sit at 15 m mod 30 | tiles |
| `BLOCK_PX` | `16` | 16 x 16 pixels = 480 m blocks, nested in the tile lattice | tiles |
| `BLOCK_MIN_ROI_FRACTION` | `0.5` | a block is part of the analysis when >= half its pixels are inside the aquifer | tiles |
| `LANDSAT_FILL` | `-32768` |  | blocks |
| `NOBS_FILL` | `255` |  | blocks |

## P1.0 block aggregates (plan v3)

| setting | default | meaning | stages |
|---|---|---|---|
| `HIST_BINS` | `16` | histogram of the 30 m values per block-month (32 is the plan's fallback) | blocks |
| `HIST_EDGES` | `dict(ET=(0.0, 240.0), NIRv=(-0.05, 0.55), NDMI=(-0.45, 0.65))` | outer bins also take values beyond | blocks |
| `OPENET_COUNT_VARS` | `["count_ssebop"]` | clear-scene count per pixel-month: mean of these OpenET model counts | blocks |
| `OBS_REGIMES` | `[("L5+L7", "2005-01"), ("L7", "2012-01"), ("L7+L8", "2013-04"), ("L...` |  | train, verify_p12 |
| `P10_VALID_SHARE_MIN` | `0.70` | P1.0 check: valid block-months Apr-Sep 2012 per class x decile | verify_p11 |

## Step 0

| setting | default | meaning | stages |
|---|---|---|---|
| `DRIVER_VARS` | `["pr", "pet", "vpd", "srad", "tmmn", "tmmx"]` | lab gridMET daily -> monthly | drivers |
| `DRIVER_MIN_DAY_FRACTION` | `0.9` | a month needs this share of usable days (sums scaled to the full month) | drivers |
| `RAIN_WINDOWS_D` | `[30, 60, 90, 180, 365]` | rain summed over [month start - w, month start - 1] | drivers |
| `ARIDITY_BASE_YEARS` | `(1985, 2014)` | aridity index P / PET climatology (same baseline as residual_et) | drivers |
| `QC_CONSTANT_IS_BAD` | `("srad", "tmmn", "tmmx", "vpd", "pet")` | rain keeps constant (all-dry) days | drivers |
| `BUCKET_CAPACITY_BINS_MM` | `[25, 50, 75, 100, 150, 200, 250, 300]` |  | drivers, static |
| `BUCKET_LP` | `0.6` | AET = PET * min(1, S / (LP * capacity))   (HBV soil routine) | drivers |
| `BUCKET_SPINUP_FROM` | `"1979-01-01"` |  | drivers |
| `AWS_FOR_BUCKET` | `"aws0_150"` |  | - |
| `SPI_VAR` | `"spi90d"` | lab drought product variable used for the drought class (stratum only) | drivers |
| `STATIC_SOIL_VARS` | `["aws0_100", "aws0_150", "rootznaws", "nccpi3all", "soc0_30"]` |  | static |
| `WTD_SHALLOW_M` | `[2, 5, 10]` | S0.5 extras: fraction of the block shallower than these depths | static |
| `NHD_GDB_URL` | `"https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHDPlu...` |  | hand |
| `TNM_API` | `"https://tnmaccess.nationalmap.gov/api/v1/products"` |  | - |
| `HU4_CODES` | `None` | None -> asked from the TNM Access API for the aquifer bounding box | hand |
| `STREAM_FCODES` | `[46000, 46003, 46006, 55800, 33400]` | stream/river (unspecified, intermittent, perennial), artificial path, connector | hand |
| `INCLUDE_EPHEMERAL` | `False` | FCode 46007 | hand |
| `MIN_DRAINAGE_KM2` | `10.0` | flowlines draining less than this (VAA TotDASqKm) are not channels | hand |
| `WATERBODY_FTYPES` | `[390, 436]` | lake/pond, reservoir; playas (361) and swamp/marsh (466) are not drainage | hand |
| `BURN_DEPTH_M` | `10.0` |  | hand |
| `HU8_BUFFER_M` | `2000.0` |  | hand |
| `CANAL_FTYPES` | `[336]` | NHD canal / ditch; kept as their own layer for the canal buffer (never channels) | hand |
| `DEM_COG_URL` | `"https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_...` |  | hand |
| `COUNTIES_URL` | `"https://www2.census.gov/geo/tiger/GENZ2023/shp/cb_2023_us_county_5...` |  | counties |
| `AQUIFER_STATE_FIPS` | `["08", "20", "31", "35", "40", "46", "48", "56"]` | CO KS NE NM OK SD TX WY | counties |

## P1.1 rules (block-year)

| setting | default | meaning | stages |
|---|---|---|---|
| `WTD_MIN_M` | `30.0` | strict depth rung: used wherever DEPTH_ACCEPTED has no entry | reference, verify_p12 |
| `WTD_SENSITIVITY_M` | `[15.0, 60.0]` | rung bits WTD15 / WTD60 (WTD10 and WTD30 come from DEPTH_RUNGS_M) | reference |
| `DEPTH_RUNGS_M` | `[60.0, 30.0, 15.0, 10.0]` | ladder fixed in advance (plan v3, P1.1) | reference, verify_p11 |
| `DEPTH_BY_FUNCTIONAL` | `dict(winter_cereal=15.0, summer_crop=15.0, other_annual=15.0, peren...` |  | reference, verify_p12 |
| `WOODY_DEPTH_M` | `30.0` |  | reference |
| `WOODY_NOISE_PX` | `3` | a block with >= this many shrub/forest pixels has woody cover | reference |
| `DEPTH_ACCEPTED` | `{}` | {"<functional>\|<decile>": rung m} copied from the ring-test report after review | reference, verify_p11 |
| `WTD_MIN_COVERAGE` | `0.95` | share of the block's aquifer pixels with a valid HRES-WTD cell | reference |
| `HAND_MIN_M` | `10.0` |  | static, reference |
| `HAND_PASS_FRACTION` | `0.9` | share of the block's pixels with HAND > HAND_MIN_M | reference |
| `IRR_BUFFER_M` | `800.0` | no LANID ever-irrigated pixel within this distance | static, reference, verify_p11 |
| `LANID_FREQ_MIN_YEARS` | `3` | sensitivity set: exclude only pixels irrigated in at least this many years | static, reference |
| `CANAL_BUFFER_M` | `800.0` | no NHD canal / ditch within this distance | static, reference |
| `FIRE_PERIMETERS` | `None` | slot: MTBS perimeters (path); None -> FIRE passes and is reported as not assessed | reference |
| `IRRIGATION_REGISTRIES` | `[]` | future: [{"name": "KS_WIMAS", "path": ..., "lon": ..., "lat": ...}] | reference |
| `WETLAND_NOISE_PX` | `3` | NWI wetland/riparian pixels tolerated in a block | reference |
| `NWI_MISSING_POLICY` | `"fail"` | block without an NWI tile: "fail" \| "pass" (only for rehearsals) | reference, verify_p11 |
| `TEST_NWI_MISSING_POLICY` | `"pass"` | test mode only (user, 2026-11-11: no NWI tile on Drive yet); reported as not applied | reference, verify_p11 |
| `WATER_NOISE_PX` | `3` | open-water pixel-months tolerated in a block-year ... | reference |
| `WATER_MAX_MONTHS` | `1` | ... and months with any open water | reference |
| `WATER_MIN_NOBS` | `2` | an open-water pixel-month needs at least this many observations | blocks |
| `GDE_SENSITIVITY_MAX` | `0.10` | GDE is NOT a main rule (user, 2026-11-04); sensitivity set only | reference |
| `CLASS_PURITY_MIN` | `0.8` | dominant CDL group share of the block | reference |
| `DRYLAND_SHARE_MIN` | `0.5` | crop is dryland-eligible in a county-year when this share lies outside LANID | dryland |
| `DRYLAND_MIN_PIXELS` | `1000` |  | dryland |
| `MIN_VALID_MONTHS` | `6` | a block-year needs this many months with all three targets | reference |
| `TARGET_MIN_PIXEL_FRACTION` | `0.5` | a block-month target needs this share of the block's aquifer pixels | reference |
| `MIN_REF_BLOCKS_PER_STRATUM` | `500` |  | verify_p11 |
| `V2_TOLERANCE` | `0.10` | water-year ET <= (P + storage change) * (1 + tolerance) | verify_p11, train |
| `V2_PASS_SHARE` | `0.95` |  | verify_p11 |
| `V4_SMD_MAX` | `0.25` | representativeness: \|standardised mean difference\| per decile | verify_p11 |
| `V4_VARS` | `["aws0_150", "rootznaws", "nccpi3all", "elevation", "P_annual", "PE...` |  | verify_p11 |

## gates (plan v3, section 6)

| setting | default | meaning | stages |
|---|---|---|---|
| `GATE_OVERRIDE` | `[]` | gate names (e.g. "G1a") whose failure may not halt focus/full runs | - |

## P1.2 model

| setting | default | meaning | stages |
|---|---|---|---|
| `TARGETS` | `["ET", "NIRv", "NDMI"]` |  | samples |
| `QUANTILES` | `[0.05, 0.5, 0.95]` |  | train |
| `ALPHA` | `0.10` |  | train |
| `DYN_INPUTS` | `["P", "PET", "VPD", "SRAD", "TMIN", "TMAX", "BUCKET", "P90d", "P365d"]` |  | train |
| `STATIC_INPUTS` | `["aws0_100", "aws0_150", "rootznaws", "nccpi3all", "soc0_30", "elev...` |  | samples, train |
| `MEAS_INPUTS` | `["oli_share", "nobs_mean", "obs_day_mean", "ET_scenes_mean"]` | how a value was observed (rule R2 holds) | samples, train, verify_p12 |
| `SPECIES_INPUT` | `True` | CDL group one-hot as a secondary class input (functional group is primary) | train |
| `EXPOSURE_VARS` | `["WTD", "wtd_min", "wtd_mean", "wtd_p05", "HAND", "hand_frac", "han...` |  | train |
| `PRESEASON_MONTHS` | `3` |  | train |
| `SPACE_GROUP_DEG` | `0.5` |  | folds |
| `SPACE_FOLDS` | `5` |  | folds |
| `YEAR_GROUPS` | `[list(range(2008, 2014)), list(range(2014, 2020)), list(range(2020,...` |  | folds |
| `EMBARGO_YEARS` | `1` |  | folds |
| `INNER_SPLIT` | `(0.70, 0.15, 0.15)` | train / validation (early stopping) / calibration, by spatial group | folds |
| `CAL_YEAR_HOLDOUT` | `False` | True: also withhold one training year group for calibration | train |
| `CAL_STRATA` | `["target", "functional", "aridity_tercile", "season", "drought", "r...` |  | train |
| `CAL_MERGE_ORDER` | `["regime", "aridity_tercile", "drought", "functional", "season"]` | dropped in this order | train |
| `CAL_MIN_N` | `1000` | calibration block-months a stratum needs before it gets its own margin | train |
| `TEST_CAL_MIN_N` | `30` |  | train |
| `SAMPLE_BLOCKS_PER_STRATUM` | `1500` | blocks sampled per aridity decile x class (all their reference years) | samples |
| `TEST_SAMPLE_BLOCKS_PER_STRATUM` | `10**9` |  | samples |
| `M2_WINDOWS_MONTHS` | `[1, 2, 3, 4, 6, 9, 12]` |  | train |
| `M2_MIN_SAMPLES` | `50` |  | train |
| `M0_MIN_YEARS` | `3` |  | train |
| `RUN_ABLATION` | `True` | also train fold (0, 0) without pre-season inputs (leakage check) | train |
| `RUN_MEAS_ABLATION` | `True` | ... and without measurement inputs | train |
| `RUN_RATIO_MODEL` | `True` | ... and with ET / PET as the ET target (sensitivity; no mass penalty) | train |
| `CLASS_FLIP` | `dict(share=0.05, seed=11)` | class-flip test: relabel this share of test blocks | train |
| `CLASS_CONFUSION` | `dict(winter_cereal="fallow", fallow="winter_cereal", summer_crop="o...` |  | train |
| `M1` | `dict(hidden=64, tcn_levels=3, kernel=3, static_hidden=64, dropout=0...` |  | train |
| `TEST_M1` | `dict(max_epochs=3, patience=2, batch=128)` |  | train |
| `PLANT` | `dict(share=0.10, ET_mm=20.0, NIRv=0.02, NDMI=0.04, seed=7)` | synthetic subsidy for verify (3) | train |
| `PASS` | `dict(coverage=(0.88, 0.92), plant_bias=0.10, fpr_max=0.07, fpr_nomi...` |  | train, verify_p12 |
| `TEST_PASS_STRATUM_MIN_N` | `30` |  | verify_p12 |
| `DROUGHT_TEST_YEARS` | `dict(d2012=[2011, 2012, 2013], d2022=[2022])` | unseen droughts of gate G1b (water years) | verify_p12 |
| `DRIFT_YEARS` | `((2008, 2017), (2018, 2025))` | LANID drift test: before / after LANID ends | verify_p12 |
