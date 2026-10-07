# Data audit of `merged_datasets/` (28 NetCDF files, 6 October 2026)

Scope: every `.nc` under `merged_datasets/` (top level and `USGS data/derived_usgs/`), compared with the
Colab notebook (or local script) that produced it. All files were opened read-only.

Method (scripts in this folder):

| Script | What it does |
|---|---|
| `dump_meta.py` | Full header dump (attributes, dims, chunking, filters) of every file |
| `scan_all.py` | **Reads every byte of all 328 variables** with auto mask/scale off, so decompression errors, raw sentinels and NaN/fill conventions are visible. Records min/max/mean, percentiles, fill/NaN/Inf counts, range violations, per-timestep valid counts, per-timestep CRC32, per-pixel valid map |
| `coord_checks.py` | Time-axis cadence/gaps/duplicates; x/y regularity against each file's GeoTransform |
| `perstep.py`, `summarize.py` | Turn the scan output into anomaly lists (empty steps, duplicated slabs, level shifts) |
| `dupcheck.py` | Bit-for-bit comparison of the redundant single-year files with their merged parents |

Run with the `MSUGW` conda env through `run.sh` (needs `import pyexpat` on `PYTHONPATH`, see memory notes).

## 0. Bottom line

* **No file is physically corrupt.** All 28 open, and all 328 variables decompressed end to end with zero
  HDF5/read errors.
* **Five files have real content defects** (missing data, garbage values, unsorted records), listed in section 1.
* **One systematic grid defect** affects eight files (section 2).
* **Metadata is wrong or missing in most files** (section 3). Two legends are actively misleading.
* Things I checked that are fine are in section 6, so you do not have to re-check them.

Items marked *source?* look like they come from the upstream product and not from your pipeline. I did not have the
raw sources on disk, so I did not compare against them.

## 1. Content defects (data missing, wrong or out of order)

### 1.1 `DROUGHT_Merged_Ogallala.nc`: the whole of 2020 is missing
* Time axis jumps 2019-12-31 → 2021-01-05 (371 days). 2020 has zero of its 73 pentads. The header says
  "Merged 42 yearly files" for 1984-2026, which is 43 years.
* Notebook `GridMet & drought indices.ipynb` loops over years and counts failures (`n_fail`) but still merges whatever
  succeeded, so a failed year disappears without stopping the run.
* 14 of 24 SPEI/EDDI indices also contain entirely-NaN scenes (2-4 each): eddi180d/270d/90d 2017-03..05,
  eddi1y/2y 2017-01..03, eddi30d 2017-07/08, spei14d..270d/30d/90d 2019-12, spei2y/5y 2019-09/10.
  Pixels have no other holes.
* The ±2.09 cap on every SPI/SPEI/EDDI is **not** a production winsorisation: the notebook contains no clipping. It
  is in the source asset (more than 1 % of values sit exactly at each cap).

### 1.2 `OCO2_L2_Lite_SIF_Ogallala_FULL.nc`: unsorted, month missing, wrong stated span
* Time is not monotonic. Rows 0..146,690 cover 2017-09-06 → 2018-02-03, then time jumps back to 2015-02-05 and runs to
  2026-07-27. Both blocks are internally sorted. The two blocks cover mostly different days (67 vs 2,087; 9 shared days,
  no shared soundings), i.e. a catch-up run was appended out of order.
* **August 2017 has no soundings at all** (the only empty month).
* True span is **2015-02-05 → 2026-07-27**; the catalogue's "6 Sep 2017" start is just the first row. The notebook's
  `START_DATE` is 2014-09-06, so Sep 2014 - Jan 2015 is also absent.
* Anyone doing `ds.isel(obs=slice(a,b))` as a time slice, or `searchsorted` on `time`, will get wrong answers.
  Sort by `time` before use.
* 103 `sif_740` values are outside [-3, 8] (max 16.4, min -4.6) even though `quality_flag` is 0 everywhere. Outliers, not fill.

### 1.3 `GRIDMET_Merged_Ogallala.nc`: finite garbage and placeholder days (*source?*)
* **`th` (wind direction): 3,869 cells on 334 days (2011-2014)** hold exactly `-1,799,820` or `-1,809,819`
  (= -9999 × 180 / × 181). They are finite, so they are counted as valid. Scattered over 3,698 positions, 2-63 per day.
* **Constant fields:** `th` = 359° everywhere on 23 different Dec 31s (1981-2023); `srad` = 0 everywhere on 2014-12-22;
  `th` and `vs` = 0 everywhere on 2025-05-10.
* **Duplicated days (impossible naturally):** `tmmn`/`tmmx` Dec 31 of 2017, 2018 and 2019 are identical to each other
  (and quantised to 0.1 K); `tmmn`/`tmmx` 2010-12-31 = 2011-01-01; `srad`, `pr` 2013-12-30 = 2013-12-31; `th`, `vs`
  2011-03-01 = 2012-02-29.
* Outliers: `pr` max 644.9 mm on 2015-10-23 at 31.89°N 96.62°W (four values above 600 mm); Palmer `z` peaks at 26.5 in
  July 2023 (12,432 values outside ±12).
* Tiny negative zeros: `bi`, `erc` have ~-4e-6 where 0 is meant (9 % of `bi` cells). Harmless but trips `>= 0` checks.

### 1.4 `ameriflux.nc`: garbage values, leaked sentinels, no metadata
* `GPP_uStar_f`: 99 values above 1e6 (up to 1.4e37) plus 51 Inf. `Reco_uStar`: 141 values above 1e6 plus 60 Inf. Sites
  CA-Man (1996), US-Prr (2012), US-Uaf (2005), US-Fcr, US-A03, CA-Cbo.
* `P`: two values of ±9999.2 at US-BRG (2018) that do not equal the declared fill -9999.0, so nothing masks them.
* Unscreened extremes: `H` 2,051 values outside [-300, 1200] (max 73,210 W m-2); `LE` 786 (max 14,980); `Tsoil` 55
  (min -1346 °C); `rH` up to 123.7 %.
* Missing data are encoded **both** as NaN (12.9 M for NEE) and as -9999 (70,781 in LE), while only -9999 is declared.
* True span is **1994 → 2021-11-02**, not 2005 as the catalogue says. `time` agrees exactly with Year/DoY/Hour for all
  32.6 M records, and each site is one contiguous block.
* The only attribute in the file is `completed_batches = 001,002`. No notebook or script in the repo creates this file,
  so it cannot be traced or regenerated.

### 1.5 Smaller content items
* `SPL4SMGP_Ogallala_FULL.nc` ends 2024-05-09 although the notebook loop runs "until now" (the VOD file runs to
  2026-09). `depth_to_water_table_from_surface` is 100 % fill in every step; `sm_rootzone_pctl` is 8.8 % valid
  (only about half of the ROI cells, and only 3 % of cells at every step).
* `SPL3SMP_E_VOD...nc`: `soil_moisture_error_*` and `vod_error_*` (4 variables) are 100 % fill. `lat`/`lon` hold -9999
  in 921 cells (all outside the ROI). The same grid in the L4 file has no holes and is identical elsewhere, so the
  coordinates can be restored from it. `vod_*` has 242/449 values slightly below `valid_min = 0` (about -1e-10);
  netCDF4/xarray will mask them silently.
* `K_hydraulic_conductivity_mday_4km.nc` has 427 valid cells outside `aquifer_mask_4km` (every other 4 km product is
  masked to it). GOSIF and ECHO-ET share a 0.05° grid but their polygon masks differ by 79 / 75 cells (0.4 %).
* `GOSIF`: `n_monthly` is not masked outside the polygon (counts 0-4 everywhere). 6 of 1,189 eight-day composites are absent (known).

## 2. Longitude spacing: 22 % of native columns never made it into the files

`build_grid_params()` (identical in the gridMET, MODIS, ET and cropland notebooks) widens the x step by
`1/cos(latitude)` for degree-based CRSs, `scale_x = scale_y / cos_lat`. For sources that are natively on a square
degree grid this makes pixels 1.286x wider than native in x only, so nearest-neighbour export skips columns.

| File | native dx | file dx | file dx / native | native columns spanned → used | skipped |
|---|---|---|---|---|---|
| GRIDMET, GRIDMET_2016, DROUGHT | 0.04167° | 0.05358° | 1.286 | 232 → 181 | 22.0 % |
| MODIS_ET_SSEBop (+2020) | 0.00965° | 0.01241° | 1.286 | 1002 → 779 | 22.3 % |
| GRACE | 0.5° | 0.6439° | 1.288 | 19 → 15 | 21.1 % |
| Global irrigation (both) | 0.0833° | 0.1072° | 1.286 | 117 → 91 | 22.2 % |
| Hydrologic soil group | 0.00228° | 0.00293° | 1.286 | 4246 → 3302 | 22.2 % |
| GFSAD1000 | 0.00893° | 0.01148° | 1.286 | 1082 → 842 | 22.2 % |

Rows (y) are native in all of them. MOD13A3, MCD12Q1 and HRES-WTD show the same 1.286 ratio. Their sources are not on
a degree grid (or have no usable native transform in the file), so I could not count skipped columns. This contradicts the
README ("native resolution with minimum changes"). The four derived 4 km products and the DSI fits inherit the gridMET grid.
The coordinates in the files are correct, so nothing is mis-located; the data are subsampled.

## 3. Metadata problems

### Wrong (actively misleading)
* **`GFSAD1000_V1_2019_Ogallala.nc`** legend is invented. Official classes (GEE catalogue): 0 non-croplands, 1 irrigation
  major, 2 irrigation minor, 3 rainfed, 4 rainfed minor fragments, 5 rainfed very minor fragments. The file says
  1 Cropland, 2 Cropland_fallow, **3 Water, 4 Urban, 5 Natural_vegetation**, and `valid_range`/`flag_values` run 0-9
  (data hold only 0-5). The product is nominal **2010**, not 2019 as the title and time coordinate say. The notebook
  comment admits the legend was a placeholder.
* **`MOD13A3_Merged_Ogallala.nc` `SolarZenith`** has no `scale_factor`, `units` or `valid_range`. Cause: the notebook's
  `BAND_METADATA` key is `SunZenith`, the GEE band is `SolarZenith`. Raw values 901-8030 are hundredths of a degree.
  `RelativeAzimuth` spans exactly ±3600 raw and is only 14 % valid (NDVI is 99.99 %); with the declared 0.01 scale that is
  ±36°, which is not a relative azimuth. Check the MOD13 user guide for the scale before using it.
* **`SPL3SMP_E_VOD` flags:** `flag_values = [0 1]`, `flag_masks` is a string of 10 numbers (or 6 for surface flag) and
  there are 9 (or 6) `flag_meanings`; counts do not match. 52 % of ROI cell-days in `surface_flag_*` are 2047 (all
  11 bits set), which means "no observation", not 11 simultaneous conditions.

### Missing or non-conforming
* No `units`/`long_name` on any data variable in GRIDMET, DROUGHT, SSEBop, GRACE, SMAP L4, OCO-2 geometry, Global
  irrigation; **none at all in `ameriflux.nc`** (not even `Conventions`, time units, or `title`).
* Irrigation has no class legend. The source (Nagaraj et al. 2021) defines 0 = none/very little irrigation,
  1 = low-to-medium (≤ 2000 ha per 86 km²), 2 = high (> 2000 ha per 86 km²).
* HSG `valid_range = [1 4]`, `flag_meanings = A B C D`, but the file holds 14 (1 cell), 24 (4,578), 34 (74,832) which
  are dual-class codes. netCDF4/xarray mask those 79,411 cells silently.
* Stale GEE leftovers: GRIDMET `date_range`/`period_mapping` end in 2022 (data run to 2026), titles say only the first
  year ("... (1979)", "(2003)", "(2000)"), `crs_transform` attributes describe the global source grid.
* Derived 4 km files (`aquifer_mask`, `dWL_*`, `K_*`): coordinate variables x/y carry `_FillValue = NaN` (not allowed on
  coordinates), no CRS/grid_mapping, no `Conventions`, `aquifer_mask:flag_values` is the *string* "0 outside, 1 inside".
  `dWL_*` `source` points to `C:\Users\AlienX\AppData\Local\Temp\opencode\...`, which will not exist later.
  `GDE_30arcsec_Ogallala_CF.nc` `source` points to the old folder `high_plains_quifer\`.
* OCO-2 `comment` mentions `sif_757_corrected`, no such variable exists (`sif_757_daily`).
* `SPL4SMGP`: data variables have no `grid_mapping`, no x/y; EASE-2 `crs` exists but is not linked.

## 4. Corrections to `MERGED_DATASETS_CATALOG.md`
* OCO-2 span (2015-02-05 → 2026-07-27, unsorted), AmeriFlux span (1994-2021), DROUGHT "complete and consistent" (2020 missing),
  "winsorisation" (it is in the source), "AmeriFlux site ids space-padded" (NUL-padded), "SMAP depth to water table is a
  critical gap" (a variable that is entirely empty over the Ogallala; check whether the SMAP name refers to peat-only data
  before treating it as a missing product).

## 5. Suggested fixes, in order of value
1. Re-export DROUGHT 2020 and re-merge; re-check the failed-year handling in the merge cell.
2. Fix `build_grid_params()` for degree CRSs (`scale_x = scale_y`) and re-export the eight affected products, or
   document the subsampling prominently.
3. Re-export GFSAD with the correct legend and year; fix the `SolarZenith` key; add the irrigation legend.
4. OCO-2: sort by `time`, append Aug 2017 and Sep 2014-Jan 2015 if available.
5. gridMET: replace the 3,869 `th` cells and the placeholder days with NaN (or re-export and compare to source); decide
   on the 2013-12-30/31 etc. duplicate days after comparing with the raw gridMET files.
6. AmeriFlux: mask |x| > 1e6, Inf and ±9999.2; add units/Conventions; keep the generating script in the repo.
7. Add `units`/`long_name` to the GEE-exported files; fix derived-file coordinate/CRS metadata.

## 6. Checked and fine
* All 328 variables read end to end with no errors (including the 29 GB gridMET file and the 7.4 GB HRES-WTD grid).
* `GRIDMET_2016` vs gridMET 2016 slice, `modis_et_v5_dekadal_2020` vs SSEBop 2020 slice, `Global_irrigation_2001` vs
  merged step 0: **bit-identical** (all variables, all steps).
* Time axes: strictly increasing, no duplicates, cadence as expected in all files except OCO-2; GRACE and SMAP gaps match
  known mission outages; GOSIF resets at year ends are normal; ECHO-ET axis passes a diurnal test (peak at hour 12 =
  local solar noon).
* x/y strictly regular and equal to each file's GeoTransform; all 4 km products share one identical grid.
* No flipped axes: gridMET `tmmx` vs latitude r = -0.85, `pr` vs longitude r = +0.85.
* PAW, dWL, mask agree (21,259 / 21,250 / 21,259). GOSIF and ECHO-ET share lat/lon exactly.
* gSSURGO (58 layers), GDE, HRES-WTD, SSEBop, MOD13A3 VI bands, MCD12Q1, SMAP L4 `sm_rootzone`, GOSIF, GRACE, GW wells:
  values within plausible ranges, no sentinel leakage.

## 7. Addendum, 7 October 2026: where the gridMET defects come from

Section 1.3 marked the gridMET oddities as "source?". Checked against Earth Engine and against the Climatology Lab's
own files (`https://www.northwestknowledge.net/metdata/data/`):

| Defect in `GRIDMET_Merged_Ogallala.nc` | Earth Engine asset | Lab's own file |
|---|---|---|
| `th` = -1,799,820 / -1,809,819 (2011-2014) | present | **absent** (0-360) |
| `th`, `vs` 2011-03-01 = 2012-02-29 | present | **absent** (the days differ) |
| `th` constant 359 on 31 Dec (e.g. 1981) | present | present |
| `tmmn` identical on 31 Dec 2017, 2018, 2019 | present | present |
| `tmmn` 2010-12-31 = 2011-01-01; `srad`, `pr` 2013-12-30 = 12-31 | present | present |
| `srad` = 0 on 2014-12-22; `th`, `vs` = 0 on 2025-05-10 | present | present |
| `pr` 644.9 mm on 2015-10-23 | present | present |

Earth Engine's older years are also a different, unquantised version of the data than the lab's current files.

Drought indices: Earth Engine has all years (2020 and 1980-83 were lost only by the old export), no placeholder or
duplicated scenes, and agrees with the lab's files where both have data. But 40 scene-bands are fully masked on
Earth Engine (EDDI in Feb-Aug 2017, SPEI in Dec 2019-Jan 2020) that the lab's files do contain. The +-2.09 bound
with about 36 distinct values per field is in the lab's files too.

Both datasets are now exported from the lab's files by `Colab notebooks/GridMET_Lab_Ogallala.ipynb`, which fixes the
grid subsampling of section 2 (index window on the native grid) and flags, without altering, the lab-side records.
