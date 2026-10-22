# Catalogue of the Ogallala datasets, 22 October 2026

This file is a new catalogue. It does not replace `MERGED_DATASETS_CATALOG.md` (24 September), which is kept unchanged.
That older catalogue still describes files that are gone (`GRIDMET_Merged`, `DROUGHT_Merged`, `ameriflux.nc`, `GLO_DEM_ogallala.tif`)
and misses about twenty files that arrived or were rebuilt since.

## 0. How this catalogue was made, and what it does not cover

| | |
|---|---|
| Folder read | `G:\MSU_GWB\datasets\merged_datasets\`: **58 files, 42 NetCDF, 42.5 GB** (top level, `USGS data\derived_usgs\`, `USGS data\derived_usgs\dsi_fit_analysis\`, `USGS_HPA_native\`) |
| Every file | opened read-only on 22 October 2026. Headers, all global attributes, every variable (dtype, dimensions, shape, units, fill, scale), coordinate ranges and time axes were read for all of them |
| Values | read as a sample, never a full scan: three time slices per gridded file (first, middle, last) with a spatial stride, about 400,000 cells per slice. "Sampled" in this file means that. Anything marked *(audit)* comes from `data audit/DATA_AUDIT_REPORT_2026-10-20.md`, which read the flagged files end to end; *(old)* means it was carried over from the 24 September catalogue and not re-checked |
| Redundancy | single-year subsets, mesh files and the wells files were compared directly (section 2). Cross-grid supersessions were judged from metadata |
| Notebooks | the Landsat, ESPA, CDL, OpenET, GEDI, LANID and other Colab notebooks were read for the datasets they produce but whose output is **not in the folder** (section 4) |
| Not done | no end-to-end scan of any file in this session; nothing was modified, moved or deleted; the Colab notebooks were not run; values of the not-on-disk datasets are therefore described from the notebooks, not from data |

Status tags used throughout:

| Tag | Meaning |
|---|---|
| **CURRENT** | built by the current notebooks from the provider's own files; follows the project rules (native grid, no manipulation, flags kept as variables) |
| **REPROCESSED** | rebuilt after the 20 October audit |
| **OLD-FLAW** | made before the audit fixes and carries a defect that was found and fixed in the notebook, but the file itself was not rebuilt |
| **PARTIAL** | export unfinished or a test slice |
| **REDUNDANT** | duplicated by another file in the folder (section 2 gives the proof) |
| **SUPERSEDED** | a newer file holds the same quantity on a better grid, but the two are not bit-identical |
| **NOT-DATA** | log, report, figure or leftover, not a dataset |

## 1. What is in the folder at a glance

| # | File | Size | Variables ("n total" = every variable in the file, including coordinates and provenance; a bare number = data variables only) | Time | Tag |
|---|---|---|---|---|---|
| **ET and energy flux** | | | | | |
| 1 | `ECHOET_Ogallala_UNH_R1.nc` | 2.8 GB | ET | hourly, 2018-01-02 to 2022-12-30 | CURRENT |
| 2 | `ECO_L3T_JET_Ogallala_ET (1).nc` | 0.97 GB | 6 ET layers | 1,200 overpasses, 2018-07-29 to 2020-08-27 | PARTIAL, deviates from the 70 m native rule |
| 3 | `MODIS_ET_SSEBop_Merged_Ogallala.nc` | 0.27 GB | et | 699 dekads, 2003-01-01 to 2022-05-21 | OLD-FLAW |
| 4 | `modis_et_v5_dekadal_2020_Ogallala.nc` | 0.11 GB | et | 36 dekads of 2020 | **REDUNDANT** |
| **SIF** | | | | | |
| 5 | `OCO2_L2_Lite_SIF_Ogallala.nc` | 1.06 GB | 111 total | 3,780 days, 2014-09-06 to 2025-09-02 | PARTIAL |
| 6 | `OCO3_L2_Lite_SIF_Ogallala.nc` | 0.57 GB | 110 total | 2,118 days, 2019-08-06 to 2026-08-30 | CURRENT |
| 7 | `TROPOSIF_L2B_v2.1_Ogallala.nc` | 1.87 GB | 40 total | 1,186 days, 2018-05-01 to 2021-12-31 | CURRENT |
| 8 | `TROPOMI_SIF740_Caltech_Ogallala.nc` | 0.43 GB | 34 total | 626 days, 2017-11-23 to 2019-10-31 | PARTIAL |
| 9 | `GOSIF_Ogallala_8day_monthly.nc` | 0.10 GB | sif, sif_monthly, n_monthly | 8-day 2000-02-26 to 2025-12-19 | CURRENT (older script) |
| 10 | `GOSIF_Ogallala_mesh.npz` | 4 KB | grid sidecar | | CURRENT |
| **Climate** | | | | | |
| 11 | `GRIDMET_Ogallala_lab_latest.nc` | 9.7 GB | 18 variables | daily, 1979-01-01 to 2026-10-05 | CURRENT |
| 12 | `GRIDMET_DROUGHT_Ogallala_lab_colab.nc` | 3.2 GB | 87 total (26 indices, 26 categories, 26 qc) | pentads, 1980-01-05 to 2026-10-02 | CURRENT (provider defect in `eddi2y`) |
| **Soil moisture, water storage, water table, soils** | | | | | |
| 13 | `SPL4SMGP_Ogallala_FULL.nc` | 0.93 GB | 5 data | 3-hourly, 2015-03-31 to 2024-05-09 | OLD-FLAW (notebook fixed 10-20, file not rebuilt) |
| 14 | `SPL3SMP_E_VOD_Ogallala_FULL.nc` | 0.19 GB | 12 data | daily, 2015-03-31 to 2026-09-14 | OLD-FLAW (same) |
| 15 | `GRACE_Merged_Ogallala.nc` | 0.16 MB | lwe_thickness, uncertainty | monthly, 2002-04-03 to 2024-09-30 | OLD-FLAW |
| 16 | `HRES-WTD_2015_Ogallala.nc` | 7.4 GB | b1 | static | OLD-FLAW |
| 17 | `Hydrologic_Soil_Group_250m_2019_Ogallala.nc` | 25 MB | b1 | static | OLD-FLAW |
| 18 | `USGS data/derived_usgs/GSSURGO_Ogallala_all_attrs_120m.nc` (+ `.report.csv`, `gssurgo_iso19115_summary.json`) | 0.94 GB | 61 total | static | CURRENT |
| **Groundwater and USGS products** | | | | | |
| 19 | `HPA_wells_water_level_timeseries.nc` | 52 MB | 40 total | 1900 to 2026, 2.27 M measurements | CURRENT |
| 20 | `USGS_HPA_native/USGS_HPA_wells.nc` | 9.7 MB | 34 total | 1901 to 2019, 24,417 records | CURRENT |
| 21 | `USGS_HPA_native/USGS_HPA_water_level_change_500m.nc` | 4.5 MB | 9 total | static | CURRENT |
| 22 | `USGS_HPA_native/USGS_HPA_hydraulic_conductivity_polygons.gpkg` | 2.8 MB | 787 polygons | static | CURRENT |
| 23 | `derived_usgs/GW_wells_HPA_CF.nc` | 1.8 MB | 19 total | 1901 to 2019, 20,580 records | **REDUNDANT** |
| 24-28 | `derived_usgs/aquifer_mask_4km.nc`, `dWL_2017_to_2019_ft_4km.nc`, `dWL_predev_to_2019_ft_4km.nc`, `K_hydraulic_conductivity_mday_4km.nc`, `PAW_rootznaws_4km.nc` | 0.06 to 0.22 MB each | 1 each | static | **SUPERSEDED** (old 287 x 181 grid) |
| 29-33 | the same five with `_reprocessed_20261020` | 0.03 to 0.12 MB each | 1 each | static | REPROCESSED |
| 34 | `derived_usgs/dsi_fit_analysis/dsi_run_fits.nc` (+ 2 csv, 2 png) | 5.5 MB | 84 data (86 total) | static | OLD-FLAW (old grid, old drought file) |
| 35 | `derived_usgs/preprocess_log.txt`, `preprocess_log_reprocessed_20261020.txt` | 1 KB | | | NOT-DATA |
| **Vegetation indices** | | | | | |
| 36 | `MOD13A3_Merged_Ogallala.nc` | 4.0 GB | 11 | monthly, 2000-02 to 2026-06 | OLD-FLAW |
| 37 | `MOD13A3_2000_Ogallala_reprocessed_20261021.nc` | 1.0 GB | 11 | 11 months of 2000 | REPROCESSED, PARTIAL (one year) |
| 38 | `LANDSAT_C2L2_indices_500m_p031r031_Ogallala.nc` | 3.7 GB | 18 indices + 2 counts | 1,102 scenes, 1982-11-18 to 2023-10-03 | PARTIAL (1 of 59 path/rows) |
| **Land cover, cropland, irrigation** | | | | | |
| 39 | `MCD12Q1_v061_500m_Ogallala.nc` | 0.13 GB | 13 | yearly, 2001 to 2024 | CURRENT |
| 40 | `MCD12Q1_Merged_Ogallala.nc` | 0.22 GB | 13 | yearly, 2001 to 2024 | **SUPERSEDED** |
| 41 | `GFSAD1000_V1_2019_Ogallala.nc` (+ `.cf.txt`) | 2.9 MB | landcover | static | OLD-FLAW |
| 42 | `Global_irrigation_Area_Merged_Ogallala.nc` | 0.06 MB | classification | yearly, 2001 to 2015 | OLD-FLAW |
| 43 | `Global_irrigation_Area_2001_Ogallala.nc` | 0.10 MB | classification | 2001 | **REDUNDANT** |
| 44 | `LANID_v3_30m_1997-2017_Ogallala.nc` | 0.67 GB | irrigated + 7 derived + 1 coarse | yearly, 1997 to 2017 | CURRENT |
| **Topography, ecosystems** | | | | | |
| 45 | `GLO30_DEM_Ogallala_1arcsec.tif` | 1.3 GB | elevation | static | CURRENT |
| 46 | `GDE_30arcsec.tif` | 0.49 GB | 5 bands, global | static | **REDUNDANT** for this project (global source raster) |
| 47 | `derived_usgs/GDE_30arcsec_Ogallala_CF.nc` | 3 MB | 5 | static | CURRENT |
| **Other** | | | | | |
| 48 | `GOSIF_ROI_mesh.npz` (in `derived_usgs/`) | 4 KB | | | **REDUNDANT** |
| 49 | `AmeriFlux_NEON_data_readme.pdf` | 31 KB | | | NOT-DATA, orphan: the tower file it documents was removed |
| 50 | `Unconfirmed 191249.crdownload` | 260 MB | | | **NOT-DATA**, an unfinished browser download |

Datasets that exist only as notebooks and whose output is **not** in the folder (Landsat 30 m, Landsat ESPA 4 km, CDL, OpenET, GEDI, ECOSTRESS 70 m) are in section 4.

## 2. Redundancy register

"Bit-identical" means every value of the overlapping records was compared in this session.

| File | Redundant with | Evidence | What to do |
|---|---|---|---|
| `modis_et_v5_dekadal_2020_Ogallala.nc` (113 MB) | the 2020 records of `MODIS_ET_SSEBop_Merged_Ogallala.nc` | all 36 dekads matched by date; **0 differing cells**; same x and y | safe to remove |
| `Global_irrigation_Area_2001_Ogallala.nc` | first record of `Global_irrigation_Area_Merged_Ogallala.nc` | **bit-identical**, same grid | safe to remove |
| `USGS data/derived_usgs/GOSIF_ROI_mesh.npz` | `GOSIF_Ogallala_mesh.npz` | mask, lat and lon **equal** | remove; it also sits in the wrong folder |
| `derived_usgs/GW_wells_HPA_CF.nc` | `USGS_HPA_native/USGS_HPA_wells.nc` | same 8,222 wells (lat and lon equal). The native file keeps every column of the three USGS shapefiles, flags repeats, and invents no dates; the CF file has fewer columns, 20,580 records (native: 24,417, of which 20,651 are not repeats) and fills missing dates with a mid-season nominal date | use the native file |
| `derived_usgs/aquifer_mask_4km.nc`, `dWL_2017_to_2019_ft_4km.nc`, `dWL_predev_to_2019_ft_4km.nc`, `K_hydraulic_conductivity_mday_4km.nc`, `PAW_rootznaws_4km.nc`, `preprocess_log.txt` | their `_reprocessed_20261020` versions | same quantities, but on the old Earth Engine grid (287 x 181, dx 0.0536 deg) that does not line up with the gridMET lab grid (287 x 233) *(audit)*. Not bit-comparable | use the reprocessed files; `residual_et/features.py` still remaps the old ones |
| `MCD12Q1_Merged_Ogallala.nc` | `MCD12Q1_v061_500m_Ogallala.nc` | same product (MCD12Q1 v061), same 24 years 2001 to 2024, same 13 variables. The old file is on a regridded EPSG:4326 grid with 1.286x too wide longitude steps; the new one is on MODIS's native grid and carries class legends. Not bit-comparable | use the new file |
| `MOD13A3_Merged_Ogallala.nc` (4.0 GB) | will be redundant with the full-record `MOD13A3_v061_1km_Ogallala.nc` once `MODIS.ipynb` is run | today the only 2000 to 2026 record, so keep it. Its `MOD13A3_2000_..._reprocessed` sibling holds only Feb-Dec 2000 | keep until the full file exists, then remove both |
| `GDE_30arcsec.tif` (global, 42,600 x 13,200) | `GDE_30arcsec_Ogallala_CF.nc` | the NetCDF is the same five bands clipped to the aquifer, masked and decoded (fractions divided by 1e8) | keep the TIF only as the raw source |
| `ECO_L3T_JET_Ogallala_ET (1).nc` | an older copy, `G:\ECO_L3T_JET_Ogallala_ET.nc` (357 MB), exists outside this folder | the " (1)" in the name is a re-download; the sizes differ, so they are not identical | decide which one to keep; neither follows the 70 m native rule |
| `Unconfirmed 191249.crdownload` | nothing | starts with the HDF5 signature, so it is an unfinished NetCDF-4 or HDF5 download; its content cannot be identified | delete or finish the download |
| `GFSAD1000_V1_2019_Ogallala.nc.cf.txt` | | an IOOS compliance report, not data | keep as a log |
| the 26 `<index>_category` variables inside `GRIDMET_DROUGHT_Ogallala_lab_colab.nc` | the 26 index variables beside them | the category is a fixed binning of the index value (provider's 11 classes), so it carries no new information; it is roughly a third of the 3.2 GB | keep (provider layer), but do not feed both to a model |

**Overlaps that are not duplicates** (different products of the same quantity; choose one on purpose):

* TROPOSIF L2B and the Caltech TROPOMI file are two retrievals from the same TROPOMI radiances; they overlap from 2018-05 to 2019-10. GOSIF is a model built from OCO-2 SIF, MODIS and meteorology, so it is not independent of OCO-2.
* ECHO-ET is a random forest trained on ECOSTRESS ET; SSEBop is also one of the OpenET models. Using them together in one model repeats information.
* Irrigation: Global irrigation (9 km, 2001 to 2015), GFSAD (1 km cropland), LANID (30 m, 1997 to 2017) and the planned CDL cover the same question at different scales.
* Landsat greenness: MOD13A3 (1 km), Landsat 500 m (one path/row), Landsat 30 m and ESPA 4 km (not on disk) describe the same vegetation signal.
* HSG and GSSURGO: `rootznaws` inside GSSURGO is the source of `PAW_rootznaws_4km`.

Outside the folder (not part of this catalogue): `G:\ECHOET_Ogallala.nc` (2.2 GB) and two `.bak-pre-timefix` files of ECHO-ET are older versions of item 1.

## 3. Datasets on disk, with their variables

Fill and scale conventions are stated per file. Where a value says *sampled*, it comes from the three-slice read described in section 0.

### 3.1 ET and energy flux

#### `ECHOET_Ogallala_UNH_R1.nc`: ECHO-ET hourly ET, University of New Hampshire (CURRENT)
* Product: ECHO-ET hourly ET, local-solar-time version; random forest trained on ERA5 hourly reanalysis, MODIS VI and ECOSTRESS 70 m ET, bias-corrected with 194 eddy-covariance sites (Bu et al. 2026, Nature Communications). Producer's own check: r2 0.95, RMSE 29 W m-2.
* Grid: 240 x 194, 0.05 deg, lat 31.725 to 43.675, lon -105.925 to -96.275, clipped to the aquifer polygon. The same lattice as GOSIF.
* Time: 43,745 hourly records (hours since 2000-01-01, **local solar time, not UTC**), 2018-01-02 to 2022-12-30. Step sizes 1, 2, 3, 4 and 7 hours occur, so there are gaps. The time axis was repaired on 2026-09-29 (it had been written in seconds).
* Variables: `ET` (float32, W m-2, valid_range 0 to 1000, fill -9999; the producer clipped negatives to 0; sampled 40.7 % valid = the polygon share); `roi_mask`; `lat`, `lon` (2-D); `crs`.

#### `ECO_L3T_JET_Ogallala_ET (1).nc`: ECOSTRESS ECO_L3T_JET v002 mosaics (PARTIAL)
* Grid: 1194 x 969, 0.01 deg (lat 31.74 to 43.67, lon -105.93 to -96.25), EPSG:4326. One record per ISS overpass: 1,200 records on 335 UTC days, 2018-07-29 to 2020-08-27 (v002 runs to the present).
* **It breaks the project's rule for ECOSTRESS.** Per its own `comment`, the 70 m MGRS tiles were averaged onto a 0.01 deg grid and cloud, water and polygon were masked. The rule recorded for the notebook is native 70 m per tile, QA kept as variables, nothing masked.
* Variables (float32, NaN fill): `ETdaily` (mm day-1), `PTJPLSMinst`, `STICinst`, `MOD16inst`, `BESSinst` (instantaneous ET of four models, W m-2), `ETinstUncertainty` (W m-2); `roi_mask`; `crs`. Sampled valid share is only 1.1 % because each overpass covers a narrow swath.

#### `MODIS_ET_SSEBop_Merged_Ogallala.nc`: SSEBop MODIS v5 dekadal ET (OLD-FLAW)
* Source: USGS SSEBop MODIS v5 (`projects/usgs-ssebop/modis_et_v5_dekadal`) through Earth Engine.
* Grid: 1235 x 779, EPSG:4326, lon -105.915 to -96.259. **x steps are 1.286x the native step (22 % of columns skipped)**; fixed in the notebook, file not rebuilt *(audit)*.
* Time: 699 dekads, 2003-01-01 to 2022-05-21 (the asset ends there).
* Variable: `et` (int16, fill -32768, **no units or scale attribute**; sampled 0 to 106). The unit is millimetres; whether a dekad value is a 10-day total is not stated by the provider. A value of 109 mm only fits a dekad total, so it very probably is one *(notebook finding)*. About a third of aquifer values are stored 0 and not masked.
* Replacement planned: `MODIS_ET_SSEBop_v5_1km_Ogallala.nc` from `ET_&_Energy_Thermal_Flux.ipynb` (1235 x 1000 cells, time = 1st, 11th, 21st of each month).

### 3.2 Solar-induced fluorescence

#### `OCO2_L2_Lite_SIF_Ogallala.nc` and `OCO3_L2_Lite_SIF_Ogallala.nc`: NASA OCO-2 and OCO-3 Lite SIF soundings
* Everything NASA publishes, at native footprint resolution, for every sounding whose centre or any footprint corner lies in the aquifer bounding box. No gridding, scaling, quality filtering or masking. Source: NASA GES DISC, `OCO2_L2_Lite_SIF` versions 11r and 11.2r (later version wins for a day; `collection_version` says which), `OCO3_L2_Lite_SIF` 11r.
* OCO-2: 6,892,940 soundings, 3,780 days, 2014-09-06 to 2025-09-02, `processing_status = in_progress` (506 of 818 version-11.2r days are in). NASA has no file for 2017-07-31 to 2017-09-18. OCO-3: 3,624,762 soundings, 2,118 days, 2019-08-06 to 2026-08-30, complete (NASA has no file for 2023-11-13 to 2024-07-15).
* Layout: a `sounding` axis (point data), a `day` axis (one row per source file) and a `sequence` axis.
* Variables, by group (111 in OCO-2, 110 in OCO-3):
  * **Indexing and provenance**: `day_index`, `source_index`, `collection_version`, `in_polygon` (centre inside the polygon), `footprint_in_polygon`, `sequence_row`, and per day `day_time`, `day_source_file`, `day_collection_version`, `day_build`, `day_production_stamp`, `day_source_etag`, `day_source_last_modified`, `day_source_size`, `day_n_soundings_file`, `day_n_soundings`, `day_n_in_polygon`, `day_n_sequences`, `day_first_sounding`, `day_first_sequence`, `day_global_*` (provider attributes), `sequence_day_index`, `sequence_local_index`.
  * **Time and geometry**: `Delta_Time`, `Geolocation_time_tai93` (true TAI), `Latitude`, `Longitude`, `Latitude_Corners`, `Longitude_Corners`, `SZA`, `VZA`, `SAz`, `VAz`, `Geolocation_*` (the same quantities again, bit-identical, plus altitude and footprint vertices).
  * **SIF**: `SIF_740nm`, `SIF_Uncertainty_740nm`, `Daily_SIF_740nm`, `Daily_SIF_757nm`, `Daily_SIF_771nm`; `Science_SIF_757nm`, `Science_SIF_771nm`, each also as `_Unadjusted_`, `_Relative_`, `_Unadjusted_Relative_` and `_Uncertainty_`; `Science_continuum_radiance_757nm/771nm`, `Science_daily_correction_factor`, `Science_sounding_land_fraction`, `Science_IGBP_index`.
  * **Quality**: `Quality_Flag`, `Science_sounding_qual_flag` (11r days only), `SimplyGoodOrBadQualityFlag`, `Science_sounding_l1b_quality_flag`, `Diagnostics_detailed_quality_bitflag_best/good` and their `Diagnostics_quality_bitflag_thresholds_*` (11.2r and OCO-3 only; all fill for the other version).
  * **Cloud and meteorology**: `Cloud_surface_albedo_abp`, `Cloud_cloud_flag_abp`, `Cloud_delta_pressure_abp`, `Cloud_co2_ratio`, `Cloud_o2_ratio`; `Meteo_surface_pressure`, `Meteo_specific_humidity`, `Meteo_vapor_pressure_deficit`, `Meteo_temperature_skin`, `Meteo_temperature_two_meter`, `Meteo_wind_speed`.
  * **Instrument metadata**: `Metadata_OrbitId`, `Metadata_SoundingId`, `Metadata_FootprintId`, `Metadata_MeasurementMode`, `Metadata_CollectionLabel`, `Metadata_BuildId`, `Sequences_*`.
  * **Per-day offset tables** (day x 227 signal bins x 8 footprints): `Offset_signal_histogram_*`, `Offset_SIF_Mean/Median/Relative_Mean/Relative_Median/Relative_SDev_757nm and _771nm`.
* **Defect (provider, audit):** `Science_SIF_Unadjusted_771nm` holds an undeclared **-4.637617e12** in 7,271 OCO-2 and 274,164 OCO-3 soundings, mostly with `Quality_Flag` 0, so quality screening keeps them. Mask it before use. Similar ±1e12 values sit in the `Offset_*` tables; `Meteo_vapor_pressure_deficit` goes to -324 Pa (OCO-2) and -499 Pa (OCO-3).

#### `TROPOSIF_L2B_v2.1_Ogallala.nc`: ESA TROPOSIF (TROPOMI) L2B v2.1
* Source: SRON open-access FTP (`TROPOSIF_L2B_<date>.nc`), native footprints (about 3.5 x 5.5 km), every sounding in the bounding box. 27,333,141 soundings, 1,186 days, 2018-05-01 to 2021-12-31 (the end of the v2.1 archive), complete; some days are absent at the provider (for example 2019-06-15).
* Variables (40): `SIF_743`, `SIF_Corr_743`, `SIF_ERROR_743` (743 to 758 nm window), `SIF_735`, `SIF_Corr_735`, `SIF_ERROR_735` (735 to 758 nm window), all in mW m-2 sr-1 nm-1; `latitude`, `longitude`, `GEOLOCATIONS_latitude_bounds/longitude_bounds` (4 corners), `GEOLOCATIONS_viewing_zenith_angle`, `GEOLOCATIONS_solar_zenith_angle`, `GEOLOCATIONS_relative_azimuth_angle`; `DETAILED_RESULTS_TOA_RFL` (7 windows), `DETAILED_RESULTS_WVL_RFL`, `DETAILED_RESULTS_Mean_TOA_RAD_743/735`; `INPUT_DATA_cloud_fraction_L2`, `INPUT_DATA_LC_MASK` (16 IGBP-like classes); time as `reference_time` (per day), `delta_time` (ms) and a computed `sounding_time`; and the same provenance/index variables as OCO (`day_*`, `in_polygon`, `footprint_in_polygon`, `day_retrieval_settings`).
* The provider declares no fill; the netCDF default 9.96921e+36 occurs (about 0.4 % of soundings in the 735 nm variables) and is declared as `_FillValue`.

#### `TROPOMI_SIF740_Caltech_Ogallala.nc`: Caltech TROPOMI SIF 740 nm (Koehler et al.) (PARTIAL)
* Native ungridded soundings. 9,659,110 soundings, 626 days, 2017-11-23 to 2019-10-31, `processing_status = in_progress` (the provider has monthly archives to 2022-12).
* Variables (34): `sif` (SIF 740 nm), `sif_err`, `sif_relative` (percent of continuum), `dcSIF` (day-length corrected), `daily_correction_factor`, `NIR` (reflected radiance 743 to 758.5 nm), `cloud_fraction` (from VIIRS), `sza`, `vza`, `saa`, `raa`, `phase_angle`, `lat`, `lon`, `lat_bnds`, `lon_bnds`, `TIME`, plus the provenance variables. Fill -999.
* **Caution:** the notebook was changed afterwards to drop damaged rows (a Caltech file of 2020-02-13 holds 109,267 rows with the wrong time and position). The file on disk stops in 2019, before that day, but it was written with the older code, so it has no `day_n_damaged_soundings` column.

#### `GOSIF_Ogallala_8day_monthly.nc` and `GOSIF_Ogallala_mesh.npz`: GOSIF v2 (Li and Xiao 2019)
* Grid 240 x 194, 0.05 deg, masked to the aquifer polygon (sampled 40.8 % finite). Time: 1,183 eight-day composites (start dates 2000-02-26 to 2025-12-19) and 311 calendar months (2000-02 to 2025-12). Five composites are absent *(audit)*.
* Variables: `sif` (8-day), `sif_monthly` (mean of the contributing composites), both in W m-2 um-1 sr-1, float32 NaN fill, small negatives are retrieval noise (sampled -0.016 to 0.19); `n_monthly` (uint8, number of composites in the month; 0, not fill, outside the footprint); `time`, `time_monthly`, `lat`, `lon`, `crs`.
* Sidecar `GOSIF_Ogallala_mesh.npz`: `bounds`, `H`, `W`, `res`, `shp_sig`, `mask` (240 x 194), `lat`, `lon`.

### 3.3 Climate

#### `GRIDMET_Ogallala_lab_latest.nc`: gridMET daily, from the Climatology Lab's own files (CURRENT)
* Grid: **287 x 233 cells of 1/24 deg** (lat 31.733 to 43.650, lon -105.933 to -96.267), the lab's native grid; the whole bounding box is stored and `roi_mask` (cells touched by the polygon) and `roi_fraction` (share of the cell inside it) describe the aquifer. This is the project's analysis grid.
* Time: 17,445 days, 1979-01-01 to 2026-10-05. Variable `status` (int8 per day): 0 permanent (17,145 days), 1 provisional (270), 2 early (30), 3 not stated.
* Values are the lab's **packed int16** with `scale_factor`, `add_offset` and `_FillValue` 32767; read with auto-scaling. Each variable also has a `<name>_qc` flag per day (constant field, identical to previous record, no data).

| Variable | Meaning | Unit | Scale, offset | Sampled range |
|---|---|---|---|---|
| `pr` | daily precipitation | mm | 0.1 | 0 to 17.8 |
| `tmmn`, `tmmx` | daily minimum, maximum temperature | K | 0.1, offsets 210 and 220 | 238 to 291, 245 to 304 |
| `vpd` | mean vapour pressure deficit | kPa | 0.01 | 0.02 to 2.14 |
| `srad` | mean downward shortwave | W m-2 | 0.1 | 75 to 215 |
| `sph` | mean specific humidity | kg kg-1 | 1e-5 | 0.0002 to 0.012 |
| `rmax`, `rmin` | daily max, min relative humidity | % | 0.1 | 22 to 100, 8 to 73 |
| `vs` | mean wind speed at 10 m | m s-1 | 0.1 | 1.1 to 14.2 |
| `th` | mean wind direction | degrees | 1 | 0 to 360 |
| `pet` | reference ET, short grass (GEE calls it `eto`) | mm | 0.1 | 0.2 to 7 |
| `etr` | reference ET, alfalfa | mm | 0.1 | 0.3 to 11.2 |
| `bi` | burning index (BI-G) | 1 | 1 | 0 to 87 |
| `erc` | energy release component (ERC-G) | 1 | 1 | 4 to 67 |
| `fm100`, `fm1000` | 100-hour, 1000-hour dead fuel moisture | % | 0.1 | 7 to 22, 10 to 26 |
| `fm1`, `fm10` | 1-hour, 10-hour dead fuel moisture | % | 0.1 | 3 to 13, 2 to 35; **start in early 2020**, empty before |

* Defects *(audit, provider)*: placeholder days exist in the lab's own files (`th` = 359 everywhere on many 31 December; `th`, `vs` or `srad` = 0 everywhere on a few days such as 2012-03-01, 2014-12-22, 2025-05-10; `tmmn`/`tmmx` identical across several 31 Decembers). Nov-Dec of 2020, 2022, 2024, 2025 are still provisional at the lab.

#### `GRIDMET_DROUGHT_Ogallala_lab_colab.nc`: gridMET drought indices (CURRENT)
* Same grid and conventions as gridMET above. 3,413 pentads (every 5 days), 1980-01-05 to 2026-10-02, complete. 87 variables.
* Indices, each as packed int16 (scale 0.01, offset -5, or -15 for PDSI and Z; fill 32767), **26 in all**: `spi`, `spei`, `eddi` at 14d, 30d, 90d, 180d, 270d, 1y, 2y, 5y; `pdsi`; `z` (Palmer Z).
* For each index also `<index>_category` (int16, the lab's 11 drought classes aligned with the US Drought Monitor) and `<index>_qc` (uint8 per pentad); plus `status`, `roi_mask`, `roi_fraction`.
* Sampled: the standardized indices stop at **±2.09** (the lab's own cap, about 36 distinct values); `pdsi` -9.1 to 7.1, `z` -2.9 to 19.4. The 2-year indices are empty for the first 73 pentads (1980) and the 5-year ones for the first 292 (1980-83).
* **Defect (provider, audit): `eddi2y` repeats 1981-83 as 1984-86** (219 fields bit-identical to the pentads three years earlier, confirmed on the lab's OPeNDAP). Treat `eddi2y` for 1984-86 as missing.

### 3.4 Soil moisture, water storage, water table and soils

#### `SPL4SMGP_Ogallala_FULL.nc`: SMAP L4 soil moisture (OLD-FLAW)
* SMAP Level-4 global geophysical, EASE-Grid 2.0 9 km, 134 x 104, 2-D lat/lon, 3-hourly, 26,623 steps, 2015-03-31 to 2024-05-09 (NASA now has data to 2026-10).
* Variables: `sm_rootzone` (m3 m-3, sampled 0.03 to 0.40), `sm_rootzone_pctl` (percentile, only ~10 % valid, valid where it rained), `land_evapotranspiration_flux` and `baseflow_flux` (model fluxes, kg m-2 s-1 as NASA gives them, not observations), `depth_to_water_table_from_surface` (**100 % fill**, NASA fills it only for peatland), `roi_mask`, `lat`, `lon`, `crs`. Fill -9999.
* The notebook (`SMAP.ipynb`) was fixed on 10-20 and 10-22 (highest build per step, output audit, no longer stops at 2024); the file was not rebuilt.

#### `SPL3SMP_E_VOD_Ogallala_FULL.nc`: SMAP L3 enhanced soil moisture and VOD (OLD-FLAW)
* SPL3SMP_E, DCA algorithm, 9 km, 134 x 104, 4,100 daily records between 2015-03-31 and 2026-09-14 (about 85 days of that span are absent). AM = descending, PM = ascending.
* Variables: `soil_moisture_am/pm` (m3 m-3), `vod_am/pm` (L-band, dimensionless), `retrieval_qual_flag_am/pm` and `surface_flag_am/pm` (uint16 bit fields; bit 0 of the retrieval flag = recommended quality), `soil_moisture_error_am/pm` and `vod_error_am/pm` (**all fill**, NASA's layers are empty), `roi_mask`, `lat`, `lon`. Fill -9999.
* Sampled: AM valid 19.5 %, PM only 1 % (PM is far sparser); `lat`/`lon` hold -9999 in 921 cells and `vod` has about 20 values near -1e-12 *(audit)*. Apply the flags before use.

#### `GRACE_Merged_Ogallala.nc`: GRACE and GRACE-FO mascons (OLD-FLAW)
* NASA GRACE mascon V04 (Earth Engine copy), 24 x 15 cells of 0.5 deg, 238 months (irregular, 10 to 374 days apart), 2002-04-03 to 2024-09-30. The native window has 19 columns; this file has 15 (the longitude flaw).
* Variables: `lwe_thickness` (liquid water equivalent anomaly, cm; sampled -34 to 18), `uncertainty` (cm; 1.6 to 17.5). No unit attributes in the file. Basin-scale constraint only.

#### `HRES-WTD_2015_Ogallala.nc`: modelled water-table depth (OLD-FLAW)
* 52,397 x 33,030 (about 30 m), 7.4 GB, one variable `b1` (float32, metres, NaN fill, "Water Table Depth"). The time label 2015-01-01 is a placeholder: the product (Ma et al. 2026) is a **long-term mean, 1914 to 2023** *(notebook finding)*. Static. 0 to 223 m, mean 23 m *(old)*.
* About 80 % of aquifer cells are nearest-neighbour from coarser state images, and the grid matches none of the source grids. Replacement planned: `HRES_WTD_Ogallala.nc` (49,820 x 40,383 with an added `source_image` band).

#### `Hydrologic_Soil_Group_250m_2019_Ogallala.nc`: HiHydroSoil v2 hydrologic soil group (OLD-FLAW)
* 5,237 x 3,302, 250 m, `b1` (uint8, fill 255): 1 = A, 2 = B, 3 = C, 4 = D; sampled 99.4 % valid, median 3. Dual codes 14, 24, 34 (A/D, B/D, C/D, inferred from the digits) cover about 0.5 % of cells *(old)*. The year label 2019 is a placeholder (Earth Engine date 2020-05-01).

#### `derived_usgs/GSSURGO_Ogallala_all_attrs_120m.nc` (+ `.report.csv`, `gssurgo_iso19115_summary.json`): gSSURGO soil attributes (CURRENT)
* 120 m NAD83 Conus Albers (EPSG:5070), 11,112 x 6,568, static, 61 variables, clipped to the aquifer (31.6 M cells inside). The CSV lists n_valid, minimum, median and maximum per column.
* Variables: `mukey` (map unit key); available water storage `aws0_5`, `aws5_20`, `aws20_50`, `aws50_100`, `aws100_150`, `aws150_999`, and the cumulative `aws0_20`, `aws0_30`, `aws0_100`, `aws0_150`, `aws0_999` (cm); component thickness `tk*a` and `tk*s` (cm) for the same depth intervals; soil organic carbon `soc*` (g m-2) for the same intervals; `musumcpcta`, `musumcpcts`, `musumcpct` (percent of map unit composition covered); `pctearthmc`; National Commodity Crop Productivity Indices `nccpi3corn`, `nccpi3soy`, `nccpi3cot`, `nccpi3sg`, `nccpi3all`; `rootznemc`, `rootznaws` (root zone depth and available water, cm; median 213), `droughty` (0/1), `pwsl1pomu` (potential wetland soil landscape; 13.7 M valid, the sparsest). Float32 with fill 9.97e+36 (mukey -9999).

### 3.5 Groundwater and USGS products

#### `HPA_wells_water_level_timeseries.nc`: all depth-to-water measurements (CURRENT)
* Point time series: **79,485 wells, 2,269,791 measurements**, three sources one after another on both axes, selected by a `source` axis (first row and count given): USGS Water Data API (47,849 wells, 1,196,552 records, 1903 to 2026), Texas Water Development Board (9,267 wells, 223,761 records, 1900 to 2026), Nebraska CSD (22,369 wells, 849,478 records, 1900 to 2025). Kansas, Colorado, Oklahoma, New Mexico, Wyoming and South Dakota state databases are **not** included.
* Station variables: `well_id`, `source`, `station_source`, `agency`, `well_name`, `state`, `county`, `aquifer`, `usgs_site_no`, `lat`, `lon`, `well_depth` (ft), `in_polygon`, `well_group` (wells of different sources within 25 m).
* Observation variables: `station_index`, `time` (days since 1900-01-01), `year`, `depth_to_water` (ft below land surface), `status`, `approval`, `measuring_agency`, `method`, `record_source` (coded, with `_text` lookup tables), `duplicate_of_other_source`. Source tables: `source_name`, `source_long_name`, `source_url`, `source_selection`, `source_citation`, `source_first_*`, `source_*_count`, `source_first_year`, `source_last_year`.
* No quality filter is applied: read `status`. *(Audit)* 17 records of one USGS well give 2,977 ft (probably water-level altitude stored as depth); 10 TWDB records have the date 1899-12-30; one Nebraska record is dated 3785; 2 Nebraska depths are -99; 10 wells have depth 99999; 842 exact in-source duplicates; 495,833 records carry only year and season or no date; 1,088 wells are outside the aquifer polygon; USGS thins out after 2017.

#### `USGS_HPA_native/USGS_HPA_wells.nc`: USGS campaign wells (CURRENT)
* 8,222 wells, 24,417 water-level records, 1901 to 2019. The three USGS shapefiles (F02 primary network, F03 supplemental incl. a 1980 epoch, F05 the 2017-19 network), every column, values unchanged, in feet. `lat`, `lon` in NAD83.
* Variables: `sitebadge`, `usgs_id`, `station_na`, `state`, `county`, `lat`, `lon`, `well_depth`; per-file flags `in_f02`, `in_f03`, `in_f05`; water-level changes `f02_deltapd_19`, `f02_dpd_19src`, `f03_deltapd_15/16/17/18`, `f03_dpd_19est`, `f03_dpd_19src`, `f03_delta80_19`, `f03_PD80_CNTR`, `f03_PD80_RANGE`, `f05_delta17_19`; records `station_index`, `source_file`, `campaign`, `time`, `date_precision`, `date_string`, `time_of_day`, `water_level` (ft below land surface), `use_code`, `repeat_of_other_file` (keep the 0s to count each measurement once), `same_measurement_as_other_campaign`. Fill -9999. Read `use_code` before using an early level: only code 1 is a true predevelopment level.

#### `USGS_HPA_native/USGS_HPA_water_level_change_500m.nc` and `..._hydraulic_conductivity_polygons.gpkg`: USGS maps (CURRENT)
* 2,659 x 1,569 cells of 500 m, NAD83 / Conus Albers, the lattice of the Landsat 500 m grid. Variables: `dWL_predev_to_2019` and `dWL_2017_to_2019` (ft; negative = decline; fill -3.4e38), `aquifer_mask` (hp_bound2010), `K_lower`, `K_upper` (ft day-1 class limits, fill -1), `K_limits_conflict` (1 for the 9 polygons whose RANGE text disagrees with MAJOR1/MINOR1).
* These are contour-based maps, not continuous surfaces: exactly 0 over 59 % (predevelopment to 2019) and 74 % (2017 to 2019) of the aquifer.
* The GeoPackage holds 787 polygons (OFR 98-548, Cederstrand and Becker 1998) with columns `ArcIds`, `AREA`, `PERIMETER`, `COND#`, `COND-ID`, `RANGE`, `MAJOR1`, `MINOR1`, `K_lower_ftday`, `K_upper_ftday`, `limits_conflict`, EPSG:5070.

#### `derived_usgs/*_reprocessed_20261020.nc`: the five 4 km products on the gridMET lab grid (REPROCESSED)
* 287 x 233, the lab grid cell for cell (built by `preprocess_usgs_gw.py` and `process_gssurgo_paw.py` from the native files). `aquifer_mask` (uint8; 27,291 cells inside, 40.8 %), `dWL_predev_to_2019` (ft; valid 27,278 cells, -265 to 54, mean -14.8), `dWL_2017_to_2019` (ft; -19.8 to 12.8, mean 0.12), `K_mday` (m day-1, class midpoint; 26,692 valid cells, 3.8 to 137, mean 19.2; 37 of 787 polygons have no class), `root_zone_available_water_storage` (cm, area mean of `rootznaws`). All float32 NaN fill, `lat_bnds`, `lon_bnds`, `crs`.

#### `derived_usgs/GW_wells_HPA_CF.nc`, old 4 km files and `dsi_fit_analysis/` (see section 2)
* `dsi_run_fits.nc` (287 x 181, old grid): 12 `best_dist_<indicator>_<DD|DS|IDD>` maps (AIC-winning distribution: exponential, lognormal, gamma, Weibull) and 72 fitted parameter maps, for the indicators `spi90d`, `spi180d`, `spei90d`, `spei180d` with threshold -1 (run characteristics DD, DS, IDD; first and last events discarded; DD and IDD in pentads). Companion `aic_winner_shares.csv`, `param_summary.csv` and two figures. Computed 2026-08-26 from the removed drought file, so its input differs from the current drought file.

### 3.6 Vegetation indices

#### `MOD13A3_Merged_Ogallala.nc`: MODIS monthly vegetation indices (OLD-FLAW)
* Source: Earth Engine `MODIS/061/MOD13A3`, 1 km, EPSG:4326 grid 1432 x 903 (longitude flaw), 317 months, 2000-02-01 to 2026-06-01.
* Variables: `NDVI`, `EVI` (int16, scale 0.0001, fill -3000), `sur_refl_b01` (red), `b02` (NIR), `b03` (blue), `b07` (SWIR) (scale 0.0001; NASA's real fill is -1000, the file declares -3000, and `valid_range` hides it when auto-masking), `SummaryQA` (0 good, 1 marginal, 2 snow/ice, 3 cloudy), `DetailedQA` (uint16 bit field), `ViewZenith`, `RelativeAzimuth` (scale 0.01 deg), `SolarZenith` (**stored raw, no scale attribute**).
* **Defects:** sampled `RelativeAzimuth` shows only **18 % valid** (the old year-merge hid about 86 % of it); `SolarZenith` must be multiplied by 0.01 by hand.

#### `MOD13A3_2000_Ogallala_reprocessed_20261021.nc` (REPROCESSED, PARTIAL)
* The same 11 variables for **Feb-Dec 2000 only**, on MODIS's native sinusoidal grid (1431 x 2454, sphere radius 6371007.181 m, `SR-ORG:6974`), with 2-D lat/lon, correct scales for all angles, `SummaryQA` flag legend. Sampled 99.9 % valid. It is a one-year result of the fixed export; the full record comes from `MODIS.ipynb`.

#### `LANDSAT_C2L2_indices_500m_p031r031_Ogallala.nc`: per-scene Landsat indices at 500 m (PARTIAL)
* Landsat Collection 2 Level-2 (Landsat 4, 5, 7, 8, 9), WRS-2 **path 31 row 31 only** (one of the 59 path/rows over the aquifer), 1,102 scenes from 1982-11-18 to 2023-10-03 (`processing_status = in_progress`; the file stops there). One record per scene, no compositing.
* Grid: the USGS High Plains 500 m grid (NAD83 / Conus Albers), 508 x 767, x -597,958 to -214,958, y 1,966,008 to 2,219,508 m, lat 40.5 to 42.9 N. Each cell is the mean of the usable 30 m pixels whose centre lies in it (about 278 pixels per full cell). Values stored as float32 holding half-precision numbers.
* Pixels left out: fill, dilated cloud, cirrus, cloud, cloud shadow, snow (`QA_PIXEL`); saturated (`QA_RADSAT`); any reflective band outside USGS's valid range 7273 to 43636.
* Variables (all float32, NaN fill, unitless): `NDVI`, `EVI`, `EVI2`, `SAVI` (L = 0.5), `MSAVI`, `OSAVI`, `GNDVI`, `NIRv`, `kNDVI`, `NDMI`, `NMDI`, `MSI`, `GVMI`, `NDWI`, `MNDWI`, `AWEInsh`, `AWEIsh`, `WI2015`; `n_clear` and `n_valid` (uint16, pixels averaged and pixels of the scene in the cell, 65535 = outside); per scene `satellite`, `cloud_cover`, `sun_elevation`, `clear_cell_fraction`, `covered_cell_fraction`, `source_scene`, `source_product_id`, `source_asset_updated`; `roi_mask`, 2-D `lat`, `lon`, `crs`.
* Sampled 25 % of cells valid (the aquifer share of this path/row); EVI and MSI have rare extremes where a denominator is near zero (audit: EVI -22 to 31, MSI up to 733 in under 0.03 % of cells). Landsat 7 is included through 2023.

### 3.7 Land cover, cropland and irrigation

#### `MCD12Q1_v061_500m_Ogallala.nc`: MODIS land cover, native grid (CURRENT)
* MCD12Q1 v061 (Terra + Aqua), 500 m (463.3 m) sinusoidal grid, 2,861 x 3,743, yearly 2001-01-01 to 2024-01-01 (24 records, `time_bnds` give the year). Cells outside the polygon hold 255; `roi_mask`; 2-D `lat`, `lon`; `source_images`, `source_update_time` record the Earth Engine image of each year.
* Variables (uint8, fill 255, each with its legend as `flag_values`/`flag_meanings`): `LC_Type1` (IGBP, 1-17), `LC_Type2` (UMD, 0-15), `LC_Type3` (LAI, 0-10), `LC_Type4` (BGC, 0-8), `LC_Type5` (plant functional types, 0-11), `LC_Prop1` (FAO-LCCS land cover), `LC_Prop2` (land use), `LC_Prop3` (surface hydrology), their `_Assessment` confidences (percent), `QC` (0-9), `LW` (1 water, 2 land, from MOD44W). Sampled valid 19.8 % (the polygon share of the window). 255 also occurs as data in `LC_Prop2_Assessment` (0.16 %).

#### `MCD12Q1_Merged_Ogallala.nc` (SUPERSEDED)
* The same 13 variables on a 2,864 x 1,806 regridded EPSG:4326 grid, with no legends and no units. Sampled `LC_Type1` 1-17, median 10 (grassland). See section 2.

#### `GFSAD1000_V1_2019_Ogallala.nc` (OLD-FLAW)
* USGS GFSAD1000 cropland extent, 1 km, 1,335 x 842, one record. `landcover` (int16, fill -32768): 0 non-cropland, 1 cropland, 2 cropland fallow, 3 water, 4 urban, 5 natural vegetation, 6 to 9 other; sampled classes 0 to 5, median 0. The date label 2019 is a placeholder: Earth Engine dates the layer 2010-01-01 (the nominal year of the data). `flag_values` has 10 entries but `flag_meanings` 7 (the compliance report flags this). A CF compliance report sits beside it (`.cf.txt`).

#### `Global_irrigation_Area_Merged_Ogallala.nc` (OLD-FLAW) and `..._2001_...nc` (REDUNDANT)
* Nagaraj et al. 2021 global irrigation maps, 144 x 91 (about 9 km), 15 annual records 2001 to 2015. `classification` (int32, fill -2,147,483,648): **0 none, 1 low-to-medium irrigation, 2 high irrigation** (the meanings come from the source's Zenodo record 4659476, not from the file, which carries no legend). Sampled counts of the 13,104 cells: class 0 falls from 6,156 (2001) to 3,527 (2015), class 1 rises from 6,622 to 9,202, class 2 stays at 326 to 506 (2 to 4 %).
* One column of cells (105.08 to 104.99 W, 21 inside the aquifer) is masked by Earth Engine but holds stored values here.

#### `LANID_v3_30m_1997-2017_Ogallala.nc`: Landsat-based Irrigation Dataset (CURRENT)
* LANID-US v3.0 (Xie, Gibbs and Lark 2021; Zenodo 5548555; CC-BY-4.0), on LANID's own 30 m NAD83 Conus Albers grid (EPSG:5070), 44,322 x 26,145, x -805,770 to -21,450 m, y 983,190 to 2,312,820 m; cells outside the polygon hold the fill. Complete. 2018 to 2020 USGS maps are not included.
* Root: `irrigated` (uint8 [year, y, x]: 1 irrigated, 0 not, 255 outside the polygon; 21 years 1997 to 2017), `roi_mask`, `time` (year), `crs`.
* Group `derived` (same grid, **open with `group="derived"`**): `irrFreq` (years irrigated, 1 to 21), `maxIrrExt` (irrigated at least once), `irrFreqChange` (years irrigated 2008-17 minus 1998-2007), `formerIrr` (not irrigated 2015-17 but at least 3 times before), `irrWholePeriod` (irrigated in both 1997-99 and 2015-17), `intermittentIrr` (the same with irrFreq <= 18), `irrFreqPasture_West` (pasture/hay irrigation frequency, western US).
* Group `derived_irrAreaGain`: `irrAreaGain` (float64, hectares, irrigation gain 1998-2007 to 2008-2017) on a 6 km grid, 225 x 133, with its own `x`, `y`, `roi_mask`.

### 3.8 Topography and ecosystems

#### `GLO30_DEM_Ogallala_1arcsec.tif`: Copernicus GLO-30 DEM (CURRENT, rebuilt)
* Copernicus DEM GLO-30 (AWS Open Data tiles), native 1 arc-second, 34,786 x 42,913 (0.000278 deg), EPSG:4326, lon -105.921 to -96.258, lat 31.743 to 43.664, clipped to the aquifer polygon. One float32 band, "elevation above the EGM2008 geoid (m)", **nodata -9999 now declared**. Statistics tags: minimum 330 m, maximum 2,395 m, mean 982 m, 100 % valid. Cell centres lie on whole arc-seconds (GDAL places the edges half a cell away). This rebuild removed the zero hole that the old file had in eastern Nebraska.

#### `GDE_30arcsec.tif` (REDUNDANT, global) and `derived_usgs/GDE_30arcsec_Ogallala_CF.nc` (CURRENT)
* Rohde et al. (2024) groundwater-dependent ecosystem map, 30 arc-second aggregates. The TIF is global (42,600 x 13,200, five uint32 bands, nodata 4,294,967,295): `GDE_sqm`, `GDE_frac_AA`, `GDE_frac_GA`, `AA_sqm`, `AA_frac_GA` (fractions stored as 1e8). Its band statistics tags are stale (mean -9999).
* The NetCDF holds the same five layers for the aquifer, decoded to float32 (m2; fractions 0 to 1), 1,431 x 1,160 cells of 0.00833 deg (its summary text says 0.01 deg; the coordinates say 0.00833), 9.3 % of cells valid *(old)*, static, CC BY 4.0.

## 4. Datasets generated by Colab notebooks whose output is not in the folder

None of these outputs is in `merged_datasets`, and none was found elsewhere under `G:\MSU_GWB`; they live on Google Drive (`MyDrive/MSUGWB/...`) or have not been run. Everything below is taken from the notebooks and from the decisions recorded while building them, so it states what the notebook writes, not what has been measured on real output. Where a real measurement exists it is named as such.

### 4.1 Landsat vegetation and water indices at 30 m, monthly: `Landsat_Indices_30m_Monthly_Ogallala.ipynb`
* **What it is:** the main Landsat vegetation-index product. Landsat 5, 7, 8 and 9 Collection 2 Level-2 surface reflectance (Tier 1), **monthly median per 30 m pixel** from 2005-01 (`START_MONTH`) to the latest month that ended 45 days ago. Nothing is resampled; no gap filling, no BRDF or terrain correction.
* **Output:** `MyDrive/MSUGWB/LANDSAT_C2L2_indices_30m_monthly_tiles_v2/<year>/<tile>/LANDSAT_C2L2_indices_30m_monthly_tiles_v2_<tile>_<year>_Ogallala.nc`. One file per tile and year, zlib level 2, chunked month x 512 x 512.
* **Grid:** Landsat's own 30 m lattice of each UTM zone (the aquifer lies in zones 13 and 14 and a corner of 15), cut into tiles of 5120 x 5120 pixels (153.6 km), named like `z14e02n29`. Zones are not merged: along the 13/14 boundary both hold the same ground. Only pixels whose centre lies inside the aquifer are kept (`roi_mask`). A pixel takes every pass of its zone that sees it, so neighbouring paths join without a seam in the grid; the values still step across path overlaps because the passes fall on different days (`n_obs` shows where).
* **Size (my earlier measurement, version 1):** about 26 bytes per pixel-month for 18 indices, 506 million aquifer pixels, so roughly **3.5 TB** for 2005 to now; version 2 adds four layers and will be somewhat larger. This is why it is not on disk.
* **Variables:**

| Group | Variables | Type |
|---|---|---|
| Greenness (9) | `NDVI`, `EVI`, `EVI2`, `SAVI` (L = 0.5), `MSAVI`, `OSAVI`, `GNDVI`, `NIRv`, `kNDVI` | int16, scale 1e-4 |
| Canopy water (4) | `NDMI`, `NMDI`, `MSI` (scale 1e-3), `GVMI` | int16, scale 1e-4 (MSI 1e-3) |
| Open water (5) | `NDWI`, `MNDWI`, `AWEInsh` (scale 2.5e-4, Feyisa 2014 form `4(G-S1) - (0.25 N + 2.75 S2)`), `AWEIsh` (scale 2e-4), `WI2015` (scale 1e-2) | int16 |
| Observation support | `n_obs` (observations behind the monthly median, 0 = none), `n_obs_oli` (how many are Landsat 8/9), `obs_day` (mean day of month of the observations), `NDVI_range`, `NDMI_range` (largest minus smallest of the month's observations) | uint8 / int16 |
| Per month | `n_scenes`, `n_scenes_by_satellite`, `source_scenes` (Earth Engine ids), `satellite`, `time`, `time_bnds` | |
| Grid | `x`, `y`, `crs`, `roi_mask` | |

* **Screening that defines the values:** an observation is used unless it is fill, dilated cloud, cirrus, cloud, cloud shadow or snow (`QA_PIXEL`), saturated (`QA_RADSAT`), outside USGS's valid reflectance range, or hazy (aerosol "high" for Landsat 8/9, `SR_ATMOS_OPACITY` > 0.3 for Landsat 5/7); Landsat 7 after 2016-12-31 is excluded. One satellite pass counts once. Water is kept.
* **Switch to know:** `HARMONIZE = True` adjusts Landsat 5/7 to Landsat 8/9 with Roy et al. (2016) coefficients, which Google now calls outdated and "not recommended" for Collection 2; the choice is written into every file and files of the two kinds must not share a folder.
* **Companion:** `Landsat_Indices_30m_Monthly_Maps_Ogallala.ipynb` averages one chosen month to 4 km (mean of the valid 30 m pixel values, moved to the 4 km EPSG:5070 grid by 480 m block centres) and draws three indices per figure; six figures for 2020-07 are in `landsat indices plot/`. Output `MyDrive/MSUGWB/LANDSAT_indices_4km_maps/<YYYY-MM>/` (`*.png`, `*_4km_<YYYY-MM>.npz`).
* **Tested:** offline against a stand-in service with fault scenarios, independent recomputation, and small real Earth Engine rehearsals (about 0.9 GB of pixels). Version 2 was **not** rehearsed against the real service; no full tile-year has run in Colab.

### 4.2 Landsat indices at 4 km from USGS ESPA: `Landsat_ESPA_Indices_4km_Ogallala.ipynb`
* ESPA computes the indices per scene on USGS servers and resizes them (bilinear) to the gridMET lab grid (233 x 287 cells of 1/24 deg). One record per scene (Landsat 4, 5, 7, 8, 9, Collection 2 Level-2), one file per WRS-2 path/row: `MyDrive/MSUGWB/LANDSAT_ESPA_indices_4km/LANDSAT_ESPA_indices_4km_p<path>r<row>_Ogallala.nc`.
* Variables: `ndvi`, `evi`, `savi`, `msavi`, `ndmi` (int16 as delivered by ESPA, scale 0.0001, fill -9999, 20000 = saturated), `qa_pixel`, `qa_radsat`, `roi_mask`, `roi_fraction`. `nbr`, `nbr2`, `ndsi` can be added in the settings.
* **Clouds are not removed** (ESPA's index files carry no screening) and the quality bands are nearest-neighbour from one 30 m pixel. What bilinear does at a 133x reduction is not documented. The ESPA side was **tested only against a stand-in** written from the API documentation; no real order has been placed. Treat as unverified.

### 4.3 USDA NASS Cropland Data Layer at 30 m: `USDA_NASS_CDL_30m_Ogallala.ipynb`
* **Output:** `MyDrive/MSUGWB/USDA_NASS_CDL_30m_tiles/<year>/<tile>/USDA_NASS_CDL_30m_tiles_<tile>_<year>_Ogallala.nc`, one small file per tile and year, 2008 to the latest year NASS has published (2025). Same 30 m UTM tiles as the Landsat and OpenET files (their `frames.json` is reused), so all three line up pixel for pixel. A tile pixel takes the CDL pixel that contains its centre (nearest neighbour, computed exactly with pyproj); no value is averaged or filled. Pixels outside the polygon hold 255.
* **Variables** (all uint8, fill 255; no time axis, the year is the file):

| Variable | Content | Source by year |
|---|---|---|
| `cropland` | crop-type class code, with NASS's class names and colours as attributes (class 0 = Background) | NASS, 2008 to latest. From 2024 NASS makes the CDL at 10 m; the notebook uses NASS's own 30 m copy |
| `confidence` | percent confidence (NASS: not a measure of accuracy) | NASS 2017 to 2023; Earth Engine 2008 to 2016; not available for 2024 on |
| `cultivated` | cultivated / non-cultivated (cultivated if classed as a crop in at least two of the last five years) | NASS 2025; Earth Engine 2013 to 2023; not available 2008 to 2012 and 2024 |
| `cropland_source`, `confidence_source`, `cultivated_source` | where the layer came from: 0 none, 1 NASS, 2 Earth Engine | |
| `roi_mask`, `crs`, `x`, `y` | grid | |

* **Cost to run:** NASS national files, about 2 GB per year for crop type and about 6.5 GB for each confidence year; deleted after each year. **Tested:** with a stand-in NASS server and real Earth Engine windows (every aquifer pixel equal to a direct lookup, resumed download, rerun skips finished files). **Not tested:** NASS's real national files.
* An older, different CDL export exists outside this folder at `G:\MSU_GWB\codes\Ogallala_Phase0\01_gee_rasters\01j_usda_cdl\`.

### 4.4 OpenET monthly evapotranspiration at 30 m: `OpenET_Monthly_30m_Ogallala.ipynb`
* **Output:** `MyDrive/MSUGWB/OPENET_v2_1_monthly_30m_tiles/<year>/<tile>/OPENET_v2_1_monthly_30m_tiles_<tile>_<year>_Ogallala.nc` (a `v2_0` switch writes to its own folder; **do not mix releases**). One file per tile and year, monthly records from 1999-10 (DisALEXI from 2001-01) to the latest month OpenET has (v2.1: to 2026-09, DisALEXI one month behind). Same 30 m tiles as Landsat and CDL.
* **Values are OpenET's own integers, untouched** (unit mm per month, scale 1); no scaling, screening or resampling. No data = -32768 (also outside the polygon). v2.1 stores no-data as 0, so the notebook reads OpenET's mask separately and marks exactly the masked pixels. A month is written only when every chosen collection has it.
* **Variables (18):**

| Variable | Content | Type |
|---|---|---|
| `et_ensemble_mad` | OpenET ensemble ET: mean of the models left after median-absolute-deviation outlier removal | int16 |
| `et_ensemble_mad_min`, `et_ensemble_mad_max` | smallest, largest model value among those kept | int16 |
| `et_ensemble_mad_count` | number of models kept | int16 |
| `et_ensemble_mad_index` | bit mask of the models used | int16 |
| `et_ensemble_sam` | simple arithmetic mean of all models | int16 |
| `et_disalexi`, `et_eemetric`, `et_geesebal`, `et_ptjpl`, `et_sims`, `et_ssebop` | monthly ET of each model (ALEXI/DisALEXI, eeMETRIC, geeSEBAL, PT-JPL, SIMS, SSEBop) | int16 |
| `count_disalexi`, `count_eemetric`, `count_geesebal`, `count_ptjpl`, `count_sims`, `count_ssebop` | cloud-free Landsat observations behind the month (0 = filled from neighbouring dates) | uint8, fill 255 |
| per file | `n_images`, `source_images`, `time`, `time_bnds`, `roi_mask`, `crs` | |

* **Daily OpenET is not produced:** it exists only through the OpenET API (400 requests a month), and the aquifer needs on the order of 100,000.
* **Tested** with real Earth Engine on 256-pixel frames (2 tiles; every pixel of every variable identical to an independent reading, both releases). **Not tested:** a whole tile, a long run, Colab. A month already in a file is never refetched, so a later OpenET reissue is not picked up unless the file is removed.
* The old OpenET cell in `ET_&_Energy_Thermal_Flux.ipynb` was removed because it took its grid from a UTM zone 10 image while the aquifer is in zones 13 and 14.

### 4.5 Other notebook-defined datasets
* **GEDI L4B gridded biomass, 1 km** (`GEDI_L4B_Biomass_1km_Ogallala.ipynb`): NASA GEDI L4B v2.1, one map for April 2019 to March 2023 (not a time series), on the product's own grid, to `GEDI_L4B_v2_1_1km_Ogallala.nc`. Ten layers named `MU` (mean aboveground biomass, Mg/ha), `SE` (standard error), `PE` (standard error in percent), `V1`, `V2` (variance components), `QF` (quality flag), `NC` (number of tracks), `NS` (number of shots), `MI` (mode of inference), `PS` (prediction stratum). The layer meanings are from memory of the user guide and are **unverified**; nothing was run against real NASA files.
* **ECOSTRESS ECO_L3T_JET at 70 m** (`ECOSTRESS_JET_Ogallala.ipynb`): one file per MGRS tile, native 70 m, v002 `ETdaily`, `cloud`, `water` only, nothing masked; about 250 GB download and 140 to 160 GB output. Not on disk; the file `ECO_L3T_JET_Ogallala_ET (1).nc` is a different, aggregated product (item 2).
* **Files the current notebooks will write, replacing old flawed ones** (none run in Colab yet): `MODIS_ET_SSEBop_v5_1km_Ogallala.nc` (replaces item 3), `MOD13A3_v061_1km_Ogallala.nc` (items 36 and 37), `GFSAD1000_V1_1km_Ogallala.nc` (item 41), `Global_irrigation_Area_9km_Ogallala.nc` (items 42 and 43), `HiHydroSoil_v2_HSG_250m_Ogallala.nc` (item 17), `HRES_WTD_Ogallala.nc` (item 16), `GRACE_MASCON_RL06.3Mv04_Ogallala.nc` (item 15), and new `*_reprocessed_<date>.nc` SMAP L4 and L3 files (items 13 and 14).

## 5. Grids, time coverage and fill values

| Grid | Files |
|---|---|
| gridMET lab, 1/24 deg, 287 x 233 | gridMET, drought, the five `_reprocessed_20261020` USGS products |
| 0.05 deg, 240 x 194 (identical lattice) | ECHO-ET, GOSIF |
| 0.01 deg, 1194 x 969 | ECOSTRESS mosaic |
| EASE-Grid 2.0 9 km, 134 x 104 | SMAP L4, SMAP L3 |
| MODIS sinusoidal, 463 m / 926 m | MCD12Q1 v061, MOD13A3 2000 |
| USGS 500 m, EPSG:5070 (same lattice) | Landsat 500 m, USGS native water-level change and K |
| 120 m / 30 m Albers, EPSG:5070 | gSSURGO, LANID |
| Old Earth Engine regridded EPSG:4326 (x steps 1.286x too wide) | SSEBop, GRACE, irrigation, GFSAD, HSG, MOD13A3 merged, MCD12Q1 merged, HRES-WTD |
| Old 287 x 181 Earth Engine grid | the five superseded 4 km USGS files, DSI fits |
| 30 m UTM tiles | Landsat 30 m, CDL, OpenET (not on disk) |
| Point soundings or wells | OCO-2, OCO-3, TROPOSIF, TROPOMI Caltech, wells |

| Quantity | Record | Gaps or limits |
|---|---|---|
| gridMET | 1979-01-01 to 2026-10-05 daily | `fm1`, `fm10` only from early 2020 |
| Drought | 1980-01-05 to 2026-10-02 pentadal | `eddi2y` 1984-86 invalid |
| MOD13A3 | 2000-02 to 2026-06 monthly | old file flawed |
| GOSIF | 2000-02-26 to 2025-12-19 | 5 composites absent |
| GRACE | 2002-04 to 2024-09 | 0.5 deg |
| SSEBop | 2003-01 to 2022-05 | asset ends there |
| Landsat 500 m | 1982-11 to 2023-10 | one path/row |
| SMAP L4 / L3 | 2015-03 to 2024-05 / 2026-09 | L3 PM sparse |
| ECHO-ET | 2018-01 to 2022-12 hourly | scattered gaps |
| TROPOSIF / Caltech | 2018-05 to 2021-12 / 2017-11 to 2019-10 | Caltech partial |
| OCO-2 / OCO-3 | 2014-09 to 2025-09 / 2019-08 to 2026-08 | OCO-2 partial |
| Wells | 1900 to 2026 | KS, CO, OK, NM, WY, SD not included |

Fill conventions in use: NaN (float files, most new products); 32767 on packed gridMET/drought int16; -9999 (ECHO-ET, SMAP, USGS wells, GLO DEM); -3000, -1000, -10000, -4000, 65535 (MODIS fields, each different); 255 (land cover, HSG, LANID, CDL); -32768 (SSEBop, GFSAD, OpenET); 9.96921e+36 (netCDF default, TROPOSIF, GSSURGO); 2^32 - 1 (GDE TIF). Mask before any statistic.

## 6. Priority actions the catalogue points to

1. Remove the proven duplicates: `modis_et_v5_dekadal_2020_Ogallala.nc`, `Global_irrigation_Area_2001_Ogallala.nc`, `derived_usgs/GOSIF_ROI_mesh.npz`, `Unconfirmed 191249.crdownload`.
2. Run the rebuild notebooks listed in section 4.5 so the eight OLD-FLAW files stop being used; until then do not combine them cell by cell with gridMET.
3. Finish or retire the PARTIAL files: OCO-2 (stops 2025-09-02), Caltech TROPOMI (stops 2019-10-31), Landsat 500 m (1 of 59 path/rows), ECOSTRESS mosaic.
4. Mask `Science_SIF_Unadjusted_771nm` = -4.637617e12 in OCO-2 and OCO-3, treat `eddi2y` 1984-86 as missing, and read `use_code` and `status` in the two well files.
5. Decide the Landsat 30 m `HARMONIZE` setting before starting the multi-terabyte run, and rehearse version 2 against the real Earth Engine first, because it has not been.
6. Update `README.md`, `residual_et/config.py` and the scripts that still open removed files; the 20 October audit lists ten.
