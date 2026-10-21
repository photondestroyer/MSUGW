# Data audit of `merged_datasets/`, 20 October 2026

Re-run of the 6 October audit (`DATA_AUDIT_REPORT.md`, kept unchanged next to this file). Since then five
NetCDF files were removed and twelve files added (eleven NetCDF and a GeoPackage); gridMET, the drought
indices and OCO-2 were rebuilt from their providers' own files. All files were opened read-only.

## Scope and method

The folder holds **49 files, 53.7 GB**: 34 NetCDF, 3 GeoTIFF, 1 GeoPackage, and side files (npz, csv, png,
json, txt, pdf). As you asked on 20 October, only the files the previous audit flagged got the full read;
everything else got header, time-axis and targeted checks.

| What | Files | How |
|---|---|---|
| Header dump | all 34 NetCDF | `dump_meta.py`: every group, attribute, chunking and filter |
| Time axes and grids | all 34 NetCDF | `coord_checks.py`: cadence, gaps, order, duplicates, GeoTransform agreement, grid identity |
| **Full read of every value** | the 17 files flagged on 6 Oct, or the files that replaced them: 289 variables, 95.6 GB uncompressed | `scan_all.py`, lazy xarray/dask reads one ~16 MB block at a time, raw values (no automatic mask/scale), packed integers unpacked for the physical-range checks; `perstep.py` turns the per-time-step counts and CRCs into anomaly lists |
| Full read | `GLO_DEM_ogallala.tif`, `GDE_30arcsec.tif`, the GeoPackage | `scan_rasters.py`, plus an exact count of DEM zeros inside the aquifer polygon (`targeted/dem_mask.py`) |
| Targeted checks | new files: OCO-3, HPA wells, Landsat 500 m, ECOSTRESS, the USGS native grids and wells | `targeted/*.py`: index tables, flags, odd records, and values on every 10th record, all through lazy reads |
| Bit comparison | redundant files | `dupcheck.py`: single-year subsets in full; unfinished gridMET/drought copies on every 25th record |
| Not re-read | GSSURGO, GW wells CF, PAW, DSI fits (unchanged since 26 Aug, clean on 6 Oct); `Hydrologic_Soil_Group_250m.tif` (12.5-gigapixel global source raster) | header only |

The 17 fully read files: `GRIDMET_Ogallala_lab_latest`, `GRIDMET_DROUGHT_Ogallala_lab_colab`,
`OCO2_L2_Lite_SIF_Ogallala` (replacing the flagged gridMET, drought and OCO-2 files), `SPL4SMGP`,
`SPL3SMP_E_VOD`, `GOSIF`, `ECHOET`, `GFSAD1000`, `MOD13A3`, both `Global_irrigation_Area` files,
`Hydrologic_Soil_Group_250m_2019`, `GDE_30arcsec_Ogallala_CF`, `aquifer_mask_4km`, `K_hydraulic_conductivity_mday_4km`
and both `dWL_*_4km`. Files flagged on 6 Oct only for their longitude spacing (SSEBop and its 2020 subset,
GRACE, MCD12Q1, HRES-WTD) got the coordinate check, which is where that problem shows.

## 0. Bottom line

* **No file is physically corrupt.** All 34 NetCDF files open; the 289 fully read variables, the DEM, the GDE
  raster and the GeoPackage read end to end with zero errors.
* **Fixed since 6 October** by the rebuilt files: the missing year 2020 of the drought indices, the unsorted
  OCO-2 file, the Earth Engine garbage in gridMET wind direction, and the subsampled gridMET/drought grid.
  `ameriflux.nc` was removed rather than fixed.
* **New content defects (section 2):** a 370 km² hole of zero elevation inside the aquifer in the DEM; a
  provider-side repeat of three years in `eddi2y`; an undeclared `-4.6e12` value in both OCO files that the
  quality flag does not screen; a few impossible records in the wells file.
* **Unfinished exports sit next to finished ones (section 3):** OCO-2 stops at 2025-09-02, Landsat holds 1 of
  59 path/rows, the ECOSTRESS file is a partial 0.01° mosaic, and two unfinished gridMET/drought copies
  (10.7 GB) duplicate the complete files.
* **The derived 4 km USGS grids no longer line up with gridMET (section 4).**
* **Everything flagged on 6 October in files that were not rebuilt is still open (section 5).**
* **Documentation and ten scripts still describe or open files that are gone (section 7).**

Items marked *provider* are in the source as distributed and were not introduced by your pipeline.

## 1. What changed since 6 October

| Removed | Added |
|---|---|
| `GRIDMET_Merged_Ogallala.nc`, `GRIDMET_2016_Ogallala.nc` | `GRIDMET_Ogallala_lab_latest.nc` (complete), `GRIDMET_Ogallala_lab_drive.nc` (unfinished copy) |
| `DROUGHT_Merged_Ogallala.nc` | `GRIDMET_DROUGHT_Ogallala_lab_colab.nc` (complete), `GRIDMET_DROUGHT_Ogallala_lab.nc` (unfinished copy) |
| `OCO2_L2_Lite_SIF_Ogallala_FULL.nc` | `OCO2_L2_Lite_SIF_Ogallala.nc`, `OCO3_L2_Lite_SIF_Ogallala.nc` |
| `ameriflux.nc` | `HPA_wells_water_level_timeseries.nc` |
| | `LANDSAT_C2L2_indices_500m_p031r031_Ogallala.nc`, `ECO_L3T_JET_Ogallala_ET (1).nc` |
| | `USGS_HPA_native/`: water-level change and K on the 500 m grid, wells, K polygons (GeoPackage) |

The three GeoTIFFs were already present on 6 October but were not audited then; they are covered here for the
first time.

### Status of the 6 October findings

| 6 Oct finding | Now |
|---|---|
| 1.1 Drought: 2020 missing, empty scenes in 2017 and 2019 | **Fixed.** 3,413 pentads, 1980-01-05 → 2026-10-02, no gaps. The only empty fields are the first 73 pentads of the 2-year indices (1980) and the first 292 of the 5-year indices (1980-83), which is when those windows have too little history. The ±2.09 cap is the lab's own (*provider*). |
| 1.2 OCO-2: unsorted, August 2017 missing, wrong span | **Fixed.** Soundings sorted by time, sounding IDs unique, index tables consistent. The 2017-07-30 → 2017-09-19 gap is NASA's (no daily files exist). New problem: the export is unfinished (section 3). |
| 1.3 gridMET: `th` = -1,799,820, placeholder and duplicated days | Earth Engine garbage **gone** (`th` 0-360). Placeholder and duplicated days **remain; they are in the lab's own files** (*provider*): `th` = 359 everywhere on 23 Dec 31sts (1981-2010; 6 Oct said to 2023), `th` = 0 everywhere on 2012-03-01, 2012-12-31 and 2025-05-10, `vs` = 0 on 2025-05-10, `srad` = 0 on 2014-12-22; `tmmn`/`tmmx` identical on Dec 31 of 2017, 2018 and 2019 and on 2010-12-31/2011-01-01; `srad`, `pr`, `tmmn`, `tmmx` identical on 2013-12-30/31. **Correction to the 7 Oct addendum:** `vs` 2012-02-29 is a copy of 2011-03-01 in the lab's file too (all 66,871 cells, 99 distinct values); only `th` differs. The six `pr` values above 600 mm (max 644.9 mm, 2015-10-23 near 31.9°N 96.6°W) lie outside the aquifer polygon and match the remnants of Hurricane Patricia, so they are probably real (not checked against gauges). |
| 1.4 `ameriflux.nc` garbage and no metadata | File **removed**; no tower data remain in the folder. `residual_et/config.py` and the README still refer to it. |
| 1.5 SMAP L4 ends 2024-05-09; `depth_to_water_table_from_surface` all fill; `sm_rootzone_pctl` 8.8 % valid | **Unchanged** (file not rebuilt). |
| 1.5 SMAP L3: four `*_error_*` variables all fill; lat/lon -9999 in 921 cells; `vod` slightly below 0 | **Unchanged**: error variables 0 % valid, lat/lon 93.4 % valid, 242 / 449 `vod` values around -1e-10. |
| 1.5 K has 427 cells outside the aquifer mask; GOSIF and ECHO-ET masks differ | **Unchanged**: 427 cells; 80 cells only in GOSIF, 76 only in ECHO-ET. |
| 1.5 GOSIF composites missing | 5 single 8-day composites are absent: 2009-01-17, 2009-03-06, 2009-10-24, 2025-07-20, 2025-11-09 (6 Oct said 6). `n_monthly` holds 0, not fill, outside the footprint. |
| 2 Longitude spacing 1.286× native (22 % of columns skipped) | **Fixed for gridMET and drought** (native 1/24° grid, 233 columns). **Remains** in SSEBop (+2020), GRACE, both irrigation files, HSG, GFSAD, MOD13A3, MCD12Q1 and HRES-WTD (not rebuilt). |
| 3 Metadata (GFSAD legend and year, MOD13A3 `SolarZenith`, SMAP flags, HSG dual classes, irrigation legend, derived-file coordinates and paths, missing units) | **All unchanged** (section 5). The OCO-2 `comment` problem went with the old file. |

## 2. New content defects

### 2.1 `GLO_DEM_ogallala.tif`: a hole of zeros inside the aquifer
* 58.9 % of the 1.49 billion pixels are exactly 0 and the file declares **no nodata value**, so the zeros count
  as elevations. 878.4 million of them lie outside the aquifer polygon (the clip fill).
* **523,074 zero pixels lie inside the aquifer polygon**, 499,891 of them more than 3 pixels from its
  boundary. 520,796 form one block within 97.14-96.26°W, 40.64-42.38°N in eastern Nebraska (about 370 km²);
  the ground just west of it is about 363 m. This is missing data, not terrain. The other 2,278 lie elsewhere (not located further).
* The catalogue reads the 0 m minimum as "low plains near sea-level datum"; the lowest ground in the High
  Plains is about 300 m. Any mean, slope or elevation covariate that includes the zeros is wrong. The
  sidecar `.aux.xml` statistics (mean 402.8 m) include them.
* Stored uncompressed (6.0 GB).

### 2.2 `GRIDMET_DROUGHT_..._lab_colab.nc`: `eddi2y` repeats 1981-83 as 1984-86 (*provider*)
Every pentad of 1984-1986 is bit-for-bit the pentad three years earlier (219 fields). Two pairs were read
again from the lab's OPeNDAP server: the repeat is in the lab's file, and the server's values equal this
file's. No other index has repeated fields. The file's `_qc` flags only compare a record with the one
before, so they do not catch this. Treat `eddi2y` for 1984-86 as missing.

### 2.3 OCO-2 and OCO-3: an undeclared `-4.637617e12` (*provider*)
`Science_SIF_Unadjusted_771nm` holds this one value in 7,271 OCO-2 soundings (all version 11.2r, 2024-04-03
onward) and 274,164 OCO-3 soundings (7.6 %, the whole record). It is not the declared `_FillValue`, and most
of these soundings have `Quality_Flag` = 0 ("best"), so quality screening keeps them. Similar values (±1e12)
sit in the `Offset_SIF_Mean/Median_*` tables, as noted when the notebook was built.

Smaller OCO items: OCO-2 `SIF_740nm` has 1,185 values outside [-3, 8] (min -5.3, max 21.6), OCO-3 27;
`Meteo_vapor_pressure_deficit` goes down to -324 Pa (OCO-2) and -499 Pa (OCO-3). Version-specific variables
are all fill for the other version, as designed: `Science_sounding_qual_flag` exists for 11r days only;
`SimplyGoodOrBadQualityFlag`, `Science_sounding_l1b_quality_flag` and the `Diagnostics_*` bit flags for
11.2r only.

### 2.4 `HPA_wells_water_level_timeseries.nc`: impossible records
2,269,791 measurements at 79,485 wells (USGS, Texas TWDB, Nebraska CSD). Index tables are consistent and
every well has measurements. Records to fix or mask:

| Records | What |
|---|---|
| 17 | USGS well `410733099573201` (Nebraska, 320 ft deep): depth to water 2,976.7-2,980.3 ft. Probably water-level altitude stored as depth |
| 10 | TWDB well `1142315`: year 2024 but date 1899-12-30 (raw -2, an Access/Excel zero date) |
| 1 | NECSD well `424057102411601`: year 3785 (date 3785-01-01) |
| 2 | NECSD: depth exactly -99 (a missing code, not declared) |
| 10 | `well_depth` = 99999 (undeclared missing code); one more well 5,098 ft deep |
| 842 | exact duplicates within a source (same well, date, depth and year) |

Also note: 495,833 records have no date (494,566 Nebraska records with only year and season, 1,267 Texas),
2,163 depths are negative (water above land surface, mostly status "Static"), and 1,088 wells lie outside
the aquifer polygon (`in_polygon` = 0), mainly in Texas, South Dakota and New Mexico. 22,006 records are
marked as repeats of another source. The USGS source thins out after 2017, as documented when the file was
built.

## 3. Unfinished, partial or redundant files

| File | State |
|---|---|
| `OCO2_L2_Lite_SIF_Ogallala.nc` | `processing_status = in_progress`, stops at **2025-09-02**. All 3,274 version-11r days are present; 506 of the 818 version-11.2r daily files NASA listed on 7 Oct (to 2026-07-27) are in. OCO-3 is complete: 2,118 days, 2019-08-06 → 2026-08-30. |
| `LANDSAT_C2L2_indices_500m_p031r031_Ogallala.nc` | `in_progress`: **one WRS-2 path/row** (31/31) of the 59 with aquifer cells, 1,102 scenes from 1982-11-18 to **2023-10-03**. Values are sound (n_clear ≤ n_valid; no index where no pixel was used); under 0.03 % of cell means fall outside usual index ranges (EVI -22 to 31, MSI up to 733), where a denominator is near zero. Landsat 7 is included through 2023. |
| `ECO_L3T_JET_Ogallala_ET (1).nc` | 1,200 records, **2018-07-29 → 2020-08-27** (ECOSTRESS v002 runs to the present). Its own `comment` says per-overpass mosaics "reprojected onto a shared EPSG:4326 grid (average resampling), cloud/water + polygon masked" at 0.01°. That is what the decisions recorded for the ECOSTRESS notebook on 4 Oct ruled out (native 70 m tiles, QA kept as variables, nothing masked). No script in the repo writes this title; the closest is `local data processing/ecostress_et_extractor.py`. The ` (1)` suffix is a browser's duplicate-download name. Records are per granule time (median step 50 s), not per overpass. Values plausible (ETdaily 0-6.2 mm/day). |
| `GRIDMET_Ogallala_lab_drive.nc` (7.9 GB) | Unfinished copy (`in_progress`, 1979 → 2018-12-31, no `fm1`/`fm10`). Identical to `_lab_latest` on every compared record. **Still read by** `plotting codes/plot_gw_clusters_yearly.py` and `plot_gw_shallow_deep_clusters.py`. |
| `GRIDMET_DROUGHT_Ogallala_lab.nc` (2.8 GB) | Unfinished copy (`in_progress`, lacks `spi90d` and `z`). Identical to `_lab_colab` on every compared record. |
| `GDE_30arcsec.tif`, `Hydrologic_Soil_Group_250m.tif` | Global source rasters (not clipped); clipped NetCDF versions exist (`derived_usgs/GDE_30arcsec_Ogallala_CF.nc`, `Hydrologic_Soil_Group_250m_2019_Ogallala.nc`). The GDE TIFF's band statistics tags say mean = -9999 (stale). |

## 4. Grids

* gridMET and drought now use the lab's native grid: 287 × 233 cells of 1/24°, x/y equal to the GeoTransform.
* **The derived 4 km USGS products** (`aquifer_mask`, `dWL_*`, `K_*`, `PAW_*`, the DSI fits) are still on the
  old Earth Engine grid: 287 × 181, dx 0.05358°. Their rows are 0.0071° (~780 m) off the new gridMET rows and
  none of their columns falls within 0.001° of a gridMET column. They cannot be combined cell by cell with the
  new gridMET/drought files. `residual_et/features.py` already remaps them (`_remap_old_4km`); other code has to
  do the same, or the products should be rebuilt from `USGS_HPA_native/`.
* The Landsat 500 m grid and the USGS native 500 m grid (EPSG:5070) lie on the same lattice (origins differ by
  whole cells). ECHO-ET and GOSIF share an identical 0.05° grid. ECOSTRESS is on its own 0.01° grid.
* The README's "gridded holdings share EPSG:4326; the ~4 km gridMET grid is the common analysis grid" no
  longer describes the folder.

## 5. Still open from 6 October (files not rebuilt)

All re-checked on 20 October against the headers and, for the flagged files, every value:
* **GFSAD1000:** invented legend (Water, Urban, ... in place of the official irrigation/rainfed classes),
  `valid_range` 0-9 for data 0-5, labelled 2019 for a nominal-2010 product.
* **MOD13A3:** `SolarZenith` has no `scale_factor`, `units` or `valid_range` (raw 901-8030, hundredths of a
  degree); `RelativeAzimuth` is ±36° after its 0.01 scale and only 13.8 % valid.
* **SMAP L3 flags:** `flag_values = [0 1]` with 10 (or 6) `flag_masks` and 9 (or 6) meanings; `surface_flag` is
  mostly 2047. **SMAP L4:** no `grid_mapping`.
* **HSG:** `valid_range = [1 4]` but 79,411 cells hold the dual classes 14, 24, 34; netCDF4/xarray hide them.
* **Irrigation:** no class legend (Nagaraj et al. 2021: 0 none or very little, 1 low-medium, 2 high).
* **Derived 4 km files:** `_FillValue = NaN` on x/y, no CRS, `aquifer_mask:flag_values` is a string, `dWL_*`
  `source` points into `AppData\Local\Temp\opencode\` (still present today, but a temp folder), the GDE CF `source` into the old `high_plains_quifer\` folder.
* **No `units`/`long_name`** on the data variables of SSEBop, GRACE, SMAP L4 and the irrigation files.
* **Longitude subsampling** in the Earth Engine exports listed in section 1.

## 6. Smaller items
* ECHO-ET (hourly, local solar time) runs 2018-01-02 00:00 → 2022-12-30 23:00 and lacks 31 hours in 8 gaps
  (the longest 9 hours, 2022-08-10).
* gridMET `fm1` and `fm10` start on 2020-01-01 (14 % of records, *provider*). 300 records are not final
  (`status` provisional 270, early 30): the last 60 days (2026-08-07 → 2026-10-05) and also Nov 2 - Dec 31 of
  2020, 2022, 2024 and 2025, which the lab's yearly files still mark provisional/early (*provider*). Select
  `status = 0` for reproducible work. Drought files carry no status (all `not_stated`).
* USGS native water-level change is exactly 0 over 59 % (predevelopment-2019) and 74.5 % (2017-19) of the
  aquifer, as the file's own comment explains (contour-based maps). The K polygons (787, valid, no overlaps)
  include 37 with unknown K (`RANGE = -1`) and 9 whose RANGE text disagrees with their class limits (1,530
  cells at 500 m), and leave 7,446 km² (1.6 %) of the aquifer uncovered.
* The two GOSIF mesh files (`GOSIF_Ogallala_mesh.npz`, `derived_usgs/GOSIF_ROI_mesh.npz`) hold the same
  mask, lat, lon and bounds in different layouts.

## 7. Documentation and code that point at files that are gone

* **`MERGED_DATASETS_CATALOG.md`** describes seven files that are not in the folder (`DROUGHT_Merged`,
  `GRIDMET_2016`, `GRIDMET_Merged`, `OCO2_..._FULL`, `SPL4SMGP_Ogallala_latest`, `ameriflux.nc`,
  `merged_ameriflux.nc`) and none of the files added since (gridMET/drought lab, OCO-2 Lite, OCO-3, ECOSTRESS,
  Landsat, HPA wells, USGS native), nor the older ECHO-ET file. Its §3.1 heading has been reduced to `###da`, so the SSEBop
  merged file has no entry. Its DEM entry misreads the zeros (2.1).
* **`README.md`** lists AmeriFlux tower records (gone), says all files are "standalone CF 1.8 compliant" (the
  new files declare CF-1.9; `aquifer_mask`, `K` and both `dWL` 4 km files declare no `Conventions`) and that everything is at native resolution
  (the ECOSTRESS file is averaged to 0.01°).
* **Scripts that open removed files:** `XGB model/gw_buffering_xgb_et.py`, `XGB model/buffering_monthly_maps_timeseries.py`,
  `temporal buffering/run_tbi_temporal_maps.py`, `run_tbi_v2_usgs.py`, `run_tbi_v3_usgs.py`,
  `temporal_buffering_index.py`, `local data processing/preprocess_usgs_gw.py`, `process_gssurgo_paw.py`,
  `plotting codes/check_and_plot_gridmet.py` (all GRIDMET_Merged / DROUGHT_Merged / merged_ameriflux), and
  `residual_et/config.py` (`ameriflux.nc`). Two plotting scripts read the unfinished `_lab_drive` copy (section 3).

## 8. Suggested fixes, in order of value
1. **DEM:** re-export the hole (97.14-96.26°W, 40.64-42.38°N) from GLO-30, declare a nodata value for the
   outside fill, and compress the file.
2. **OCO:** mask `-4.637617e12` in `Science_SIF_Unadjusted_771nm` (both files) and the ±1e12 `Offset` values
   before any use, or list them in `source_missing_values`.
3. **Drought:** treat `eddi2y` 1984-86 as missing and tell the Climatology Lab.
4. **Wells:** mask or correct the records in 2.4 and drop the 842 exact duplicates.
5. **Finish the exports:** OCO-2 from 2025-09-03; Landsat (58 path/rows to go); ECOSTRESS with the per-tile
   notebook, then retire the ` (1)` file.
6. **Remove or archive** `GRIDMET_Ogallala_lab_drive.nc` and `GRIDMET_DROUGHT_Ogallala_lab.nc` (10.7 GB) after
   pointing the two plotting scripts at `_lab_latest`.
7. **Rebuild the 4 km USGS products** on the lab gridMET grid (from `USGS_HPA_native/`), or document the
   remap every user must do.
8. **Update the ten scripts** (including `residual_et/config.py`) to the current file names.
9. **Rewrite the catalogue and README** (section 7).
10. The remaining 6 October fixes for files not rebuilt (section 5).

## 9. Checked and fine
* Every read in this audit completed without a decompression or HDF5 error.
* Time axes strictly increasing with no duplicates in all gridded files; gridMET complete 1979-01-01 →
  2026-10-05; drought complete 1980-01-05 → 2026-10-02; GRACE and SMAP gaps match known mission outages.
* OCO-2 and OCO-3: soundings sorted, IDs unique, daily counts add up to the sounding axis, day blocks contiguous.
* x/y equal each file's GeoTransform; the derived 4 km grids are identical to one another; ECHO-ET and GOSIF
  coordinates identical.
* Bit-identical: SSEBop 2020 subset vs merged, irrigation 2001 vs merged (all records); the two unfinished
  gridMET/drought copies vs the complete files (every 25th record).
* gridMET and drought values within physical bounds apart from the items above, and no undeclared missing
  codes in either; SMAP, GOSIF, ECHO-ET, MOD13A3 VI bands and the derived 4 km files within plausible ranges.
* ECOSTRESS: no value outside `roi_mask`; every sampled record holds data. Landsat: no value where no pixel
  was used.
* GeoPackage: 787 valid polygons, no overlaps.

## 10. Reproducing this audit
```bash
cd "G:/MSU_GWB/datasets/data audit"
bash run.sh dump_meta.py
bash run.sh coord_checks.py
AUDIT_TAG=flagged AUDIT_WORKERS=2 bash run.sh scan_all.py <file-name fragments>
bash run.sh summarize.py flagged
bash run.sh perstep.py
bash run.sh scan_rasters.py GLO_DEM GDE_30arcsec.tif gpkg
bash run.sh dupcheck.py
bash run.sh targeted/dem_mask.py      # and the other targeted/*.py checks
```
Results go to `results/<AUDIT_RUN>/` (default `2026-10-20`). `run.sh` adds the `pyexpat` fix the `MSUGW` env
needs and caps BLAS/OpenMP threads at one per process (with 20 threads, importing numpy alone reserved
0.6 GB). The first version of the updated scripts kept each block's sample as a NumPy view, which held the
whole block in memory: up to 9 GB per scan worker and 11 GB in the raster scan, enough to exhaust the PC's
memory next to the other open applications. Samples are now copied. Expect about 0.4 GB per scan worker
and 32 minutes for the flagged set with two workers.
