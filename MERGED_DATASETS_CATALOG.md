# Catalog of the Merged Datasets Repository — Ogallala Aquifer

## 1. Purpose and scope

This document provides a meticulous, file-by-file inventory of the data holdings
in `G:\MSU_GWB\datasets\merged_datasets\` (24 entries as of 24 September 2026;
top-level files only, plus explicitly noted derived products). For each file
it records every variable present, the temporal extent of the data, and an
assessment of data quality. Files are grouped into nine thematic categories,
namely: (1) ET and Energy/Thermal Flux; (2) Carbon Flux and Photosynthesis;
(3) Climate and Meteorological Forcing; (4) Hydrology, Soil Moisture, Soils and
Groundwater; (5) Vegetation Indices and Phenology; (6) Land Cover and Surface
Type; (7) Topography and DEMs; (8) Cropland and Agricultural Management; and
(9) In-situ Flux and Eco Networks. Only holdings physically present in the
folder are documented.

## 2. Methods and conventions

All files were opened read-only. No source file was modified. NetCDF files
were inspected through metadata headers in full, and data values were examined
exclusively through strided samples (a maximum of eight time steps spread
across each record and a spatial stride limiting each read to approximately
250,000 cells), so that statistics reported below are estimates subject to
sampling variation. GeoTIFF files were inspected through raster profiles with
band statistics computed on a decimated read.

Two conventions apply throughout. First, most NetCDF files in this repository
carry Goddard Earth Engine export metadata but no physical-unit attributes.
Where units are stated, they are inferred from the source product specification
and flagged [UNIT-ASSUMPTION]. Second, each file uses its own missing-value
convention (NaN, -9999.0, -32768, 255, or 2^32-1); the convention in force is
stated per file, and all completeness fractions reported below exclude fill
values.

## 3. Category 1 — ET and Energy/Thermal Flux

### 3.1 MODIS_ET_SSEBop_Merged_Ogallala.nc (270.7 MB)

Source product: MODIS SSEBop evapotranspiration. Grid: 1235 by 779
(approximately 1 km), EPSG:4326. Temporal extent: 1 January 2003 to
21 May 2022, 699 steps at a median interval of 10 days (dekadal).

| Variable | Dimensions | Type | Observed range (sample) | Remarks |
|---|---|---|---|---|
| et | time, y, x | int16 | 0 to 106 | Fill value -32768. Units are not recorded in the file [UNIT-ASSUMPTION: mm per dekad as stored]. |

Quality: 100 percent of sampled cells contain valid (non-fill) data. Values
are non-negative throughout the sample, consistent with an evapotranspiration
flux. The record terminates in May 2022 and therefore does not cover the most
recent drought years in full.

### 3.2 modis_et_v5_dekadal_2020_Ogallala.nc (113.3 MB)

A single-year extract on the same 1235 by 779 grid: 1 January to
21 December 2020, 36 steps at a 10-day interval. Single variable `et`
(int16, fill -32768, observed 0 to 109). Quality is equivalent to Section 3.1.
This file is redundant with the corresponding year of the merged file above
and appears to be a processing intermediate; the merged file supersedes it for
multi-year analysis.

## 4. Category 2 — Carbon Flux and Photosynthesis

### 4.1 OCO2_L2_Lite_SIF_Ogallala_FULL.nc (129.7 MB)

Source product: OCO-2 Level-2 Lite solar-induced fluorescence soundings.
Structure: 5,323,544 discrete soundings (one-dimensional `obs` dimension, no
gridded axes). Temporal extent: 6 September 2017 to 27 July 2026.

| Variable | Type | Observed range (sample) | Remarks |
|---|---|---|---|
| time | float64 | — | Seconds since 1970-01-01 (Unix epoch), proleptic Gregorian |
| lat | float32 | degrees north | Fill -9999.0 |
| lon | float32 | degrees east | Fill -9999.0 |
| sif_740 | float32 | -0.30 to 3.11 | W m-2 sr-1 um-1; instantaneous 740 nm retrieval |
| sif_740_unc | float32 | 0.29 to 0.57 | Per-sounding uncertainty; same units |
| sif_757_daily | float32 | -0.15 to 0.43 | Daily-corrected 757 nm channel |
| sif_771_daily | float32 | -0.11 to 0.78 | Daily-corrected 771 nm channel |
| sif_740_daily | float32 | -0.10 to 1.19 | Daily-corrected 740 nm channel |
| solar_zenith | float32 | 26.8 to 59.5 | Degrees |
| sensor_zenith | float32 | 0.1 to 33.3 | Degrees |
| rel_azimuth | float32 | 11.2 to 179.8 | Degrees |
| quality_flag | int8 | 0 throughout sample | Fill -1; sampled soundings are uniformly best-quality |

Quality: soundings carry full per-retrieval uncertainty and geometry metadata,
and the sampled quality flags are uniformly zero. As point soundings rather
than a gridded time series, the file requires binning or geostatistical
aggregation before pixel-scale analysis, and its coverage is irregular in
space and time. Negative SIF values within the uncertainty envelope are
normal retrieval noise, not corruption.

### 4.2 GOSIF_Ogallala_8day_monthly.nc (96.5 MB) with GOSIF_Ogallala_mesh.npz (3.9 KB)

Source product: GOSIF v2 solar-induced fluorescence, 8-day global composites
(Li and Xiao 2019), clipped to the Ogallala aquifer and aggregated to calendar
months by `process_gosif_ogallala.py`. Grid: 240 by 194 (0.05 degree,
EPSG:4326), longitude -105.95 to -96.25, latitude 31.70 to 43.70, masked to
the hp_bound2010 polygon (NaN outside it). Temporal extent: 1,183 eight-day
steps from 26 February 2000 to 19 December 2025, with 311 calendar months
from February 2000 to December 2025. Six of the 1,189 source scenes are
absent from the time axis (unreadable or skipped during production); the
record is otherwise continuous at the nominal 8-day cadence.

| Variable | Dimensions | Type | Observed range (sample) | Remarks |
|---|---|---|---|---|
| time | time (1183) | float64 | 2000-02-26..2025-12-19 | Days since 1970-01-01; 8-day composite start dates |
| time_monthly | time_monthly (311) | float64 | 2000-02-01..2025-12-01 | First day of each calendar month |
| lat, lon | 240 / 194 | float64 | — | Cell centers, degrees |
| sif | time, lat, lon | float32 | -0.018 to 0.71 | W m-2 um-1 sr-1; 40.8 percent finite, matching the in-polygon share |
| sif_monthly | time_monthly, lat, lon | float32 | -0.041 to 0.77 | Same units; 40.8 percent finite; cell_methods time_monthly mean |
| n_monthly | time_monthly, lat, lon | uint8 | Counts 0-4 | Contributing 8-day composites per month |
| crs | scalar | int32 | — | latitude_longitude, EPSG:4326 |

Global attributes record CF-1.8 conventions, title, summary, source,
reference (doi.org/10.3390/rs11050517), creation history, geospatial bounds
and resolution. The companion mesh sidecar stores the window bounds, grid
dimensions, resolution, shapefile signature, polygon mask, and lat/lon
vectors needed to reproduce the extraction geometry. An earlier mesh file
(`GOSIF_ROI_mesh.npz`) exists under `USGS data/derived_usgs/`; the top-level
mesh matches the current product.

Quality: finite fractions coincide with the aquifer mask share, monthly
means are consistent with the 8-day series, and small negative values are
normal retrieval noise. This is the repository's gridded photosynthesis
record and the appropriate carbon-leg input for multi-sensor analysis.

## 5. Category 3 — Climate and Meteorological Forcing

### 5.1 GRIDMET_Merged_Ogallala.nc (29.27 GB)

Source product: gridMET daily surface meteorology. Grid: 287 by 181
(approximately 4 km), EPSG:4326; this is the de facto common analysis grid of
the repository. Temporal extent: 1 January 1979 to 1 August 2026, 17,380
daily steps. All variables are float32 with NaN fill; sampled completeness is
100 percent throughout.

| Variable | Interpretation [UNIT-ASSUMPTION] | Observed sample range |
|---|---|---|
| pr | Daily precipitation, mm | 0 to 47.5, mean 1.73 |
| tmmn | Daily minimum temperature, Kelvin | 238.1 to 297.8 |
| tmmx | Daily maximum temperature, Kelvin | 245.2 to 312.9 |
| vpd | Vapor pressure deficit, kPa | 0 to 3.73 |
| srad | Downward shortwave radiation, W m-2 | 29.1 to 345.8 |
| sph | Specific humidity, kg kg-1 | 0.00021 to 0.01764 |
| rmax, rmin | Relative humidity bounds, percent | 22.5 to 100 |
| vs | Wind speed, m s-1 | 0.96 to 16.0 |
| th | Wind direction, degrees | 0 to 360 |
| eto, etr | Grass/alfalfa reference ET, mm | 0.1 to 10.4 / 0.1 to 15.6 |
| bi | Burning index | 0 to 85 |
| erc | Energy release component | 0.05 to 83 |
| fm100, fm1000 | 100/1000-hour fuel moisture | 5.4 to 24.9 / 7.6 to 26.8 |

Quality: excellent. Complete records, physically plausible ranges, and
temperatures correctly stored in Kelvin (conversion required before use in
degree-Celsius formulations). At 29.3 GB this is the largest file in the
repository; all access must remain lazy (per-time-slice reads).

### 5.2 GRIDMET_2016_Ogallala.nc (1.54 GB)

A single-year extract (1 January to 31 December 2016, 366 daily steps) on the
identical 287 by 181 grid carrying the same 16 variables with equivalent
ranges. This is a processing intermediate; the merged file in Section 5.1
supersedes it.

### 5.3 DROUGHT_Merged_Ogallala.nc (4.61 GB)

Source family: gridMET-derived standardized drought indicators. Grid:
287 by 181 (approximately 4 km). Temporal extent: 5 January 1984 to
29 July 2026, 3,035 steps at a median interval of 5 days. All 26 indicator
variables are float32 with NaN fill and 100 percent sampled completeness.

| Variable group | Members | Observed range (sample) |
|---|---|---|
| SPI (precipitation-only) | spi14d, spi30d, spi90d, spi180d, spi270d, spi1y, spi2y, spi5y | Capped at ±2.09 |
| SPEI (precipitation minus evaporative demand) | spei14d through spei5y, same scales | Capped at ±2.09 |
| EDDI (evaporative demand) | eddi14d through eddi5y, same scales | Capped at ±2.09 |
| Palmer indices | pdsi | -9.85 to 7.22 |
| Palmer Z-index | z | -4.63 to 9.70 |

Quality: complete and internally consistent; sample means lie near zero as
expected of standardized anomalies. Two observations require attention.
First, every SPI, SPEI and EDDI variable is hard-capped at ±2.09, evidently a
winsorization applied during production; distribution tails beyond two
standard deviations are therefore truncated, which compresses extreme-drought
(and extreme-wet) signals. Second, EDDI shares the sign convention of the
other indices in this file as stored; users must confirm directionality
against the EDDI literature before interpreting it, since published EDDI
conventionally takes positive values for drier conditions.

## 6. Category 4 — Hydrology, Soil Moisture, Soils and Groundwater

### 6.1 SPL4SMGP_Ogallala_FULL.nc (931.9 MB)

Source product: SMAP Level-4 Soil Moisture Geophysical (land-model
assimilation). Grid: 134 by 104 EASE cells (approximately 9 km) with
two-dimensional latitude/longitude coordinates; an `roi_mask` variable marks
the 40.2 percent of cells inside the region of interest. Temporal resolution
is 3-hourly. Temporal extent: 31 March 2015 to 9 May 2024 (26,623 steps).
This is the sole SMAP Level-4 holding in the folder.

| Variable | Observed valid range (fill-excluded) | Remarks |
|---|---|---|
| sm_rootzone | 0.030 to 0.459 m3 m-3 | Root-zone soil moisture; 37.0 percent valid (fill -9999.0 elsewhere, concentrated outside the ROI) |
| sm_rootzone_pctl | 0 to 99.99 | Percentile transform; only 7.5 percent valid — a sparse product that cannot support dense analysis |
| land_evapotranspiration_flux | -0.00009 to 0.00030 | Model-assimilated flux, not an independent observation |
| baseflow_flux | 0 to 0.0000095 | Model-assimilated flux, near-zero throughout the sample |
| depth_to_water_table_from_surface | No valid data | 0.0 percent valid: the variable is entirely fill over the sampled record. This is a critical gap — no dynamic water-table series is available from this product |
| roi_mask | 0/1 flags | 40.2 percent inside ROI |

Quality: the root-zone moisture record itself is sound within its valid
domain, but three qualifications apply. First, effective coverage is the ROI
fraction (approximately 37 percent of stored cells). Second, the percentile
layer is too sparse for continuous use. Third, all flux and water-table-depth
variables are model-assimilated outputs of the SMAP Level-4 system rather
than independent observations, and the water-table variable is entirely
absent; none of them may be treated as ground truth for well dynamics.

### 6.2 GRACE_Merged_Ogallala.nc (0.16 MB)

Source product: GRACE/GRACE-FO mascon total water storage. Grid: 24 by 15
(approximately 50 km). Temporal extent: 3 April 2002 to 30 September 2024,
238 monthly steps.

| Variable | Type | Observed range (sample) | Remarks |
|---|---|---|---|
| lwe_thickness | float64 | -34.3 to 23.2 | Liquid-water-equivalent thickness anomaly [UNIT-ASSUMPTION: cm] |
| uncertainty | float64 | 1.3 to 17.5 | Formal mascon uncertainty, same units |

Quality: complete records with uncertainties of the expected magnitude.
Resolution (~50 km, 360 cells) restricts legitimate use to basin-scale
storage constraints and screening, consistent with the research plan; the
product cannot resolve pixel-scale depletion.

### 6.3 HRES-WTD_2015_Ogallala.nc (7.42 GB)

Source product: high-resolution modeled water-table depth (Ma et al. style
CONUS product), static snapshot dated 1 January 2015. Grid: 52,397 by 33,030
(approximately 30 m), EPSG:4326 — the finest grid in the repository and the
second-largest file. Single variable `b1` (float32, NaN fill, units meters,
long name "Water Table Depth").

Quality: 99.998 percent valid in the sampled window; observed range 0 to
222.8 m with a mean of 23.2 m, physically plausible for the High Plains. The
file is static: it provides the spatial depth axis and distance-to-threshold
information only, and must not be interpreted as a time series. Its 7.4 GB
size requires strided or windowed access at all times.

### 6.4 Hydrologic_Soil_Group_250m_2019_Ogallala.nc (25.3 MB) and Hydrologic_Soil_Group_250m.tif (153.2 MB)

Source product: HiHydroSoil hydrologic soil groups. The NetCDF (5,237 by
3,302, static 2019) carries variable `b1` (uint8, fill 255, valid range 1-4,
flag meanings A/B/C/D). The GeoTIFF is the global 250 m source raster
(158,159 by 79,080, uint8, nodata 255).

Quality: within the Ogallala bounding box, 99.84 percent of cells are valid
and coded 1-4 with the expected distribution (group C dominant at roughly
86 percent, followed by D, B and A). Approximately 0.5 percent of in-ROI
cells carry dual-group codes 14, 24 and 34 (drained/undrained complexes),
which any reclassification scheme must handle explicitly rather than silently
dropping. Global samples of the source TIFF contain additional out-of-range
codes (up to 34 observed), consistent with non-soil or miscellaneous classes
outside the region of interest.

### 6.5 SPL3SMP_E_VOD_Ogallala_FULL.nc (193.1 MB)

Source product: SMAP Level-3 enhanced vegetation optical depth and surface
soil moisture (DCA algorithm), AM (descending) and PM (ascending) passes.
Grid: 134 by 104. Temporal extent: 31 March 2015 to 14 September 2026,
4,100 daily steps.

| Variable | Observed valid range (fill-excluded) | Remarks |
|---|---|---|
| soil_moisture_am / _pm | to 0.30 / 0.34 m3 m-3 | DCA surface soil moisture; fill -9999.0 |
| vod_am / vod_pm | to 0.40 / 0.37, dimensionless | L-band vegetation optical depth; fill -9999.0 |
| soil_moisture_error_am/pm, vod_error_am/pm | entirely fill | Uncertainty layers contain no data in this extract |
| retrieval_qual_flag_am/pm | uint16 bit field | Bit 0 = recommended quality; quality filtering is mandatory before use |
| surface_flag_am/pm | uint16 bit field | Open water, precipitation, snow, frozen ground, slope, urban |
| roi_mask | 0/1 | 40.2 percent inside ROI |

Quality: the geophysical layers (soil moisture, VOD) are present with
plausible ranges, but the accompanying uncertainty layers are entirely fill
and therefore unusable, and valid coverage is limited to the ROI fraction.
The AM/PM pair permits diurnal-difference diagnostics. Users must apply the
retrieval quality flags; unfiltered use will admit frozen-ground, snow and
open-water contamination.

### 6.6 GDE_30arcsec.tif (490.7 MB)

Source product: global groundwater-dependent ecosystem mapping (Rohde et
al.), 30-arcsecond aggregates. Grid: 42,600 by 13,200 global (EPSG:4326),
five uint32 bands with nodata 2^32-1: `GDE_sqm` (GDE area, m2, observed
0-853,510), `GDE_frac_AA` and `GDE_frac_GA` (fractions scaled by 10^8),
`AA_sqm` (analyzed area, m2, observed 691-854,795) and `AA_frac_GA`
(analyzed fraction scaled by 10^8).

Quality: within the Ogallala bounding box approximately 19.8 percent of
cells carry data; the remainder lies outside the GDE model domain (nodata by
design, principally the eastern sub-humid portion of the aquifer). Fraction
bands attain exactly 10^8, confirming the documented scaling. This is a
static contextual layer for stratifying groundwater-dependent systems, not a
time series.

### 6.7 GDE_30arcsec_Ogallala_CF.nc (3.0 MB; under `USGS data/derived_usgs/`)

Derived product: the five GDE aggregate layers of Section 6.6, windowed to
the aquifer bounding box, masked to the hp_bound2010 polygon, and decoded
from stored integers to physical values. Grid: 1,431 by 1,160 (native
30-arcsecond lattice preserved). Static (no time dimension). CF-1.8 with a
`crs` grid-mapping variable (latitude_longitude, EPSG:4326); licensed
CC BY 4.0 with citation of Rohde et al. (2024) required.

| Variable | Dimensions | Type | Observed range | Remarks |
|---|---|---|---|---|
| GDE_sqm | y, x | float32 | 0 to 723,600 m2 | GDE area within grid cell |
| GDE_frac_AA | y, x | float32 | 0 to 1, valid_range recorded | Fraction of analyzed area that is a GDE |
| GDE_frac_GA | y, x | float32 | 0 to 1, valid_range recorded | Fraction of grid-cell area that is a GDE |
| AA_sqm | y, x | float32 | 692 to 729,500 m2 | Analyzed area within grid cell |
| AA_frac_GA | y, x | float32 | 0.001 to 1, valid_range recorded | Fraction of grid-cell area analyzed |
| crs | scalar | int32 | — | Grid mapping |

Quality: 9.3 percent of cells valid (in-polygon analyzed cells); all
fraction layers lie inside [0, 1]. This analysis-ready form is preferable to
the source TIFF for any work on these grids, as scaling, masking and
projection metadata are already resolved. Provenance attributes record the
former source paths; folder reorganization since production does not affect
content.

## 7. Category 5 — Vegetation Indices and Phenology

### 7.1 MOD13A3_Merged_Ogallala.nc (4.05 GB)

Source product: MODIS monthly vegetation indices. Grid: 1,432 by 903
(approximately 1 km). Temporal extent: 1 February 2000 to 1 June 2026,
317 monthly steps.

| Variable | Type/scale | Observed range (sample) | Remarks |
|---|---|---|---|
| NDVI | int16 x0.0001, fill -3000, valid -2000..10000 | 99.995 percent in range | Primary greenness record |
| EVI | int16 x0.0001, fill -3000, valid -2000..10000 | 99.995 percent in range | Primary greenness record |
| sur_refl_b01 (red) | int16 x0.0001, fill -3000, valid 0..10000 | 99.995 percent in range | |
| sur_refl_b02 (NIR) | int16 x0.0001, fill -3000, valid 0..10000 | 99.995 percent in range | |
| sur_refl_b03 (blue) | int16 x0.0001, fill -3000, valid 0..10000 | 99.995 percent in range | |
| sur_refl_b07 (SWIR) | int16 x0.0001, fill -3000, valid 0..10000 | 99.995 percent in range | |
| SummaryQA | int16, fill -1, values -1..3 | Mean 0.09 | Reliability flag; predominantly good quality |
| DetailedQA | uint16 bit field, fill 65535 | — | Per-pixel VI quality indicators |
| RelativeAzimuth | int16 x0.01 degrees, fill -4000 | — | Observation geometry |
| SolarZenith | int16, fill -32768 | — | Observation geometry (no scale attribute recorded) |
| ViewZenith | int16 x0.01 degrees, fill -10000 | — | Observation geometry |

Quality: excellent. Effectively all reflectance and index values fall inside
their documented valid ranges after fill exclusion, and the reliability flag
is near zero (good quality) across the sample. Scale factors are recorded in
file attributes and must be applied (values are stored as integers). This is
the repository's long-record greenness source for context, region screening
and reversibility analysis; per the research plan it must not be treated as
an independent carbon signal.

## 8. Category 6 — Land Cover and Surface Type

### 8.1 MCD12Q1_Merged_Ogallala.nc (217.7 MB)

Source product: MODIS annual land cover. Grid: 2,864 by 1,806
(approximately 500 m). Temporal extent: 2001 to 2024, 24 annual steps.
All variables are uint8 with fill value 255.

| Variable | Observed range (sample) | Remarks |
|---|---|---|
| LC_Type1 | 1-17 | Primary IGBP classification; mean 10.3 (grassland-dominated sample) |
| LC_Type2 | 0-15 | UMD scheme |
| LC_Type3 | 0-10 | LAI scheme |
| LC_Type4 | 0-8 | BGC scheme |
| LC_Type5 | 0-11 | PFT scheme |
| LC_Prop1/2/3 | 1-51 | Sub-pixel proportions |
| LC_Prop1/2/3_Assessment | 34-255 | Confidence assessments (34 minimum observed; 255 is fill) |
| LW | 1-2 | Land/water mask |
| QC | 0-9 | Quality control, mean 0.04 (overwhelmingly good quality) |

Quality: complete records with near-zero QC flags. Five parallel
classification schemes allow plant-functional-type stratification under
alternative legends; LC_Type1 (IGBP) is the recommended default for
consistency with the wider literature.

## 9. Category 7 — Topography and DEMs

### 9.1 GLO_DEM_ogallala.tif (5.97 GB)

Source product: GLO-30 digital elevation model extract, pre-clipped to the
Ogallala bounding box. Grid: 34,785 by 42,913 (1/3,600 degree, approximately
30 m), EPSG:4326, longitude -105.92 to -96.26, latitude 31.74 to 43.66.
Single float32 band; no band description recorded.

| Property | Value |
|---|---|
| Nodata tag | None recorded |
| Sampled range (decimated read) | 0.00 to 2,206.12 m, mean 404.17 m |
| Sampled completeness | 100 percent finite |

Quality: relief is physically plausible for the High Plains (low plains
near sea-level datum rising westward past 2,000 m). One caveat governs use:
with no nodata tag recorded, the sampled minimum of exactly 0.00 m cannot
be distinguished from a potential fill sentinel; the source documentation
should be consulted before masking on zero, and any derived slope or
elevation covariate should be validated against an independent DEM. At
6.0 GB this is the second-largest file in the repository; windowed or
strided access is mandatory.

## 10. Category 8 — Cropland and Agricultural Management

### 10.1 GFSAD1000_V1_2019_Ogallala.nc (2.9 MB)

Source product: GFSAD cropland extent classification, static 2019. Grid:
1,335 by 842 (approximately 1 km). Single variable `landcover` (int16,
fill -32768, valid range 0-9) with a documented legend: 0 Non-cropland,
1 Cropland, 2 Cropland_fallow, 3 Water, 4 Urban, 5 Natural_vegetation,
6-9 Other classes.

Quality: complete records; observed classes in the sample are 0 through 5
(mean 1.5). This is the primary static cropland mask for separating managed
from natural systems.

### 10.2 Global_irrigation_Area_Merged_Ogallala.nc (59.6 KB) and Global_irrigation_Area_2001_Ogallala.nc (98.2 KB)

Source product: global irrigated-area mapping. Grid: 144 by 91
(approximately 8 km). The merged file spans 2001 to 2015 (15 annual steps);
the 2001 file is a single-year extract on the same grid. Single variable
`classification` (int32, fill -2,147,483,648) taking values 0, 1 and 2.

Quality: complete records. A material deficiency is that the files carry no
class legend: the meanings of codes 0, 1 and 2 are undocumented in file
metadata and must be confirmed against the source publication before any
irrigation stratification is performed. The coarse 8 km grid additionally
limits the product to approximate stratification duty.

Users requiring sub-kilometre crop-type detail should note that the finest
agricultural-management grids in this folder are the 1 km GFSAD mask above
and the 8 km irrigation grids; no 30 m crop-type layer is present.

## 11. Category 9 — In-situ Flux and Eco Networks

### 11.1 ameriflux.nc (3.11 GB)

Source product: AmeriFlux/NEON eddy-covariance tower records compiled to a
single table. Structure: 290 sites by 32,628,912 half-hourly records
(one-dimensional `record` dimension with `site_index` linkage and a 290-entry
`site_lookup` table of space-padded 32-character site identifiers).

Temporal extent: records span 2005 (earliest Unix timestamp 1,115,974,800)
through 2021 (Year variable 2005-2021; Hour 0-23.5 confirms half-hourly
resolution; DayOfYear 1-326 observed).

| Variable | Observed range (sample) | Remarks |
|---|---|---|
| NEE | -8.44 to 2.17 | Net ecosystem exchange; fill -9999.0 |
| LE | -1.37 to 94.6 | Latent heat flux; fill -9999.0; 65.9 percent valid excluding fill |
| H | -18.6 to 197.0 | Sensible heat flux; fill -9999.0 |
| Rg | 0 to 716.9 | Incoming shortwave; fill -9999.0 |
| Rg_QC, Tair_QC, rH_QC, VPD_QC | 0-1 | Quality flags |
| Tair | 8.7 to 30.2 | Air temperature; fill -9999.0 |
| Tsoil | 1.3 to 24.2 | Soil temperature; fill -9999.0 |
| rH | 15.9 to 93.2 | Relative humidity; fill -9999.0 |
| VPD | 0.99 to 36.1 | Vapor pressure deficit; fill -9999.0 |
| P | 0 to 0.073 | Precipitation; fill -9999.0 |
| P_QC | 0-1 | Quality flag (mean 0.5) |
| NEE_uStar_f, NEE_uStar_fqc | — | u*-filtered NEE and its flag (0-3) |
| Reco_uStar | 0.83 to 4.57 | Ecosystem respiration estimate |
| GPP_uStar_f, GPP_uStar_fqc | 0.18 to 13.0 / flag 0-3 | Gross primary productivity and flag |
| Year, DoY, Hour | 2005-2021 / 1-326 / 0-23.5 | Time coordinates (fill -9999) |
| time | Unix seconds (float64) | No CF time units attribute recorded |
| site_lookup, site_id, site_index | — | 290 sites, indices 0-289; identifiers are space-padded fixed strings |

A companion document, `AmeriFlux_NEON_data_readme.pdf` (30.8 KB), is present
in the same folder and should be consulted for usage policy and variable
definitions.

Quality: the core fluxes and meteorology are present at expected magnitudes
with documented fill conventions; approximately two-thirds of sampled LE
records are valid. Two cautions apply. First, the `time` variable lacks a
units attribute (values are verifiably Unix seconds); any reader must assert
this rather than assume it. Second, a repository note: an identically sized
file previously inventoried as `merged_ameriflux.nc` (3,106,230,064 bytes)
no longer exists under that name following the September 2026 repository
reorganization, while `ameriflux.nc` carries exactly that byte count; the two
are consistent with a pure rename, but byte-level equivalence was not
directly verified and should be confirmed by checksum before the old name is
cited.

In-situ validation from tower records is therefore restricted to the
AmeriFlux/NEON network holdings above.

## 12. Repository-level observations

**Renames and supersessions.** The September 2026 reorganization renamed
several holdings: `merged_ameriflux.nc` to `ameriflux.nc` (byte count
identical); the `high_plains_quifer/` folder to `HPA_polygon/`; derived
products consolidated under `merged_datasets/USGS data/derived_usgs/`; and
`GDE_30arcsec.tif` moved into `merged_datasets/`. The shorter SMAP Level-4
extract (`SPL4SMGP_Ogallala_latest.nc`, ending January 2024) has been removed;
`SPL4SMGP_Ogallala_FULL.nc` (ending May 2024) is the sole Level-4 holding.
`GRIDMET_2016_Ogallala.nc` and `modis_et_v5_dekadal_2020_Ogallala.nc`
are single-year extracts redundant with their merged parents. New arrivals
since the previous revision are `GLO_DEM_ogallala.tif`,
`GOSIF_Ogallala_8day_monthly.nc` and `GOSIF_Ogallala_mesh.npz`.

**Grids.** Native resolutions span five orders of magnitude: 30 m (water-table
depth, elevation model), 250 m (soils), 500 m (land cover), ~1 km (SSEBop ET,
MODIS VI, GFSAD), 0.05 degree (GOSIF SIF), ~4 km (GRIDMET, drought indices —
the de facto common analysis grid), ~8 km (irrigation), ~9 km EASE (SMAP),
~50 km (GRACE), point soundings (OCO-2) and point towers (AmeriFlux). All
gridded products share EPSG:4326. Harmonization to the 4 km grid with
nearest-neighbour treatment of categorical fields is required before joint
analysis.

**Temporal coverage.** The longest records are GRIDMET (1979-2026, daily) and
drought indices (1984-2026, pentadal). The satellite era constrains joint
analysis: GOSIF SIF and MODIS VI from 2000, SSEBop ET 2003-2022, GRACE
2002-2024, SMAP 2015-2024/2026, OCO-2 SIF 2017-2026, AmeriFlux records
2005-2021. The maximum overlap of all dynamic sensors is therefore
approximately 2017-2021.

**Fill conventions.** Five conventions coexist (NaN; -9999.0; -32768/-3000/
-4000/-10000 with per-variable scales; 255; 2^32-1). Several files store
scaled integers whose scale factors are recorded in attributes (MODIS VI
products) while others store unscaled integers (irrigation, land cover).
Masking must always precede statistics; unmasked means in this catalog's
source tables are dominated by fill values and are not reported as data.

## 13. Data-quality summary

| File | Completeness (fill-excluded, sampled) | Range plausibility | Verdict |
|---|---|---|---|
| DROUGHT_Merged_Ogallala.nc | 100% | Indices capped ±2.09 (production winsorization); PDSI/Z physical | Sound; note tail truncation |
| GRIDMET_Merged_Ogallala.nc | 100% | All variables physical; Kelvin temperatures | Sound; convert temperature units |
| GRIDMET_2016_Ogallala.nc | 100% | As above | Sound; redundant extract |
| SSEBop merged + v5-2020 | 100% in sample | 0-109, non-negative | Sound; record ends May 2022 |
| SPL4SMGP FULL | 37% (ROI-limited); pctl 7.5% | SM 0.03-0.46 physical | Usable within ROI; pctl sparse |
| GOSIF 8-day/monthly | 40.8% (polygon-masked) | SIF -0.04..0.77, small negatives are noise | Sound; 6 of 1189 source dates absent |
| GLO_DEM tif | 100% in decimated sample | 0-2206 m, mean 404 m | Sound; nodata tag absent, verify zero handling |
| SMAP depth-to-water-table | 0% (entirely fill) | — | Missing; critical gap |
| SMAP fluxes (ET, baseflow) | 37% | Near-zero magnitudes | Model-assimilated; screening only |
| OCO-2 SIF | Full record, flags all best-quality | Physical incl. normal negative noise | Sound; requires gridding |
| MOD13A3 NDVI/EVI | 99.995% in valid range | Integer-scaled, QA~0 | Excellent |
| MCD12Q1 | Complete; QC~0 | Classes 1-17 observed | Sound |
| HRES-WTD 2015 | 99.998% | 0-222.8 m | Sound; strictly static |
| HSG nc + tif | 99.8% in ROI; classes 1-4 + 0.5% dual codes | — | Sound; handle dual codes |
| GRACE | 100% | ±35 with uncertainties 1-17 | Sound at basin scale only |
| GFSAD 2019 | Complete | Classes 0-5 observed | Sound; static |
| Global irrigation (both) | Complete | Codes 0/1/2, legend absent | Usable only after legend confirmed |
| GDE tif | ~20% in ROI (domain-limited) | Fractions attain documented 10^8 scaling | Sound within model domain |
| ameriflux.nc | LE 66% valid excl. fill | Physical magnitudes | Sound; time units unattributed |
| SPL3SMP VOD | Geophysical layers present; error layers entirely fill | SM to 0.34, VOD to 0.40 | Usable with flag filtering; uncertainties missing |


### 13.1 USGS derived 4 km products: aquifer_mask, water-level change, hydraulic conductivity
All files in this subsection share the 287 by 181 GRIDMET-matched 4 km EPSG:4326 grid. Provenance is recorded in preprocess_log.txt, which must be cited with the files.

aquifer_mask_4km.nc (0.06 MB): Variable aquifer_mask, uint8, (y, x), flag_values 0 outside, 1 inside, long_name High Plains aquifer boundary (hp_bound2010). Valid 51,947 of 51,947 (100 percent of grid defined), mean 0.4092, inside pixels 21,259 (40.9 percent of grid). This mask defines the valid domain for all 4 km derived products. Quality is definitional; accuracy inherits the hp_bound2010 polygon.

dWL_2017_to_2019_ft_4km.nc (0.21 MB): Variable dWL_2017_to_2019, float32, (y, x), FillValue NaN, long_name Mapped water-level change, 2017 to 2019, units feet. Source 500 m USGS TIFF hp_wlc1719t reprojected by averaging to 4 km. Valid 21,250 of 51,947 (40.9 percent), matching the mask (9 pixels fewer than mask due to source nodata). Sampled statistics: range -18.50 to 12.36 ft, mean 0.12 ft (log reports mean 0.11 ft, NaN 56.2 percent pre-masking). Negative denotes decline.

dWL_predev_to_2019_ft_4km.nc (0.21 MB): Variable dWL_predev_to_2019, float32, (y, x), FillValue NaN, long_name Mapped water-level change, predevelopment (approximately 1950) to 2019, units feet. Source 500 m TIFF hp_wlcpd19t. Valid 21,250 of 51,947 (40.9 percent). Range -264.97 to 53.73 ft, mean -13.81 ft (sampled mean -14.80 ft on masked grid; difference reflects masking and floating-point aggregation). NaN 56.2 percent pre-masking per log.

K_hydraulic_conductivity_mday_4km.nc (0.21 MB): Variable K_mday, float32, (y, x), FillValue NaN, long_name Horizontal hydraulic conductivity (zone midpoint), units m per day. Source USGS Open-File Report 98-548 digital hydraulic conductivity map (Cederstrand and Becker 1998), E00 PAL zones, 787 polygons. Zone range strings in ft per day converted by interval midpoint multiplied by 0.3048, rasterized by area-descending paint. Valid 21,342 of 51,947 (41.1 percent), masked to 20,915 inside aquifer per log. Range 0.00 to 137.16 m per day, mean 19.01 m per day, 90th percentile 45.72 m per day. Classes include -1, 0 to 25, 10 to 25, 25 to 50, 50 to 100, 100 to 200, 200 to 300, 300 to 400, 400 to 500 ft per day.

Quality of data: High derivation transparency; original USGS zips, TIFFs and E00 files unmodified per global attributes. Water-level change grids are suitable for regional decline analysis but are smoothed by 500 m to 4 km averaging and should not be interpreted at well scale. Hydraulic conductivity is a coarse zonation midpoint approximation suitable for regional modeling only; within-zone heterogeneity is not represented and the -1 class requires special handling per source documentation.

4.7 GW_wells_HPA_CF.nc
Source datasets: USGS High Plains water-level program shapefiles: F02 primary predevelopment to 2019 network, F03 supplemental network including 1980 epoch, F05 2017 to 2019 network. Native coordinates EPSG:5070 with NAD83 latitude and longitude attributes (approximately 1 m from WGS84). File size 1.72 MB. Conventions CF-1.8, featureType point. License: USGS data are public domain. History notes sentinels -9999 and -999 converted to NaN and mid-season nominal dates substituted where measurement-date blanks occurred in some epochs (documented approximation).

Spatial extent and resolution: Point observations. Dimensions: station 8,222, obs 20,580. Variables lat (degrees_north), lon (degrees_east), float64. Geospatial bounds latitude 31.8636 to 43.5156, longitude -105.1836 to -96.2837.

Temporal extent and resolution: Variable time, float64, units days since 1900-01-01, calendar gregorian, dimensions obs 20,580. Range 1901-07-25 to 2019-06-15. Sampling is irregular campaign-based; variables station_index (obs) and campaign (8-character labels; pd denotes predevelopment approximately 1950, otherwise calendar year; 1980 epoch only in F03) link measurements to stations.

Variables present:

Identifiers: site_badge (24 characters), usgs_id (32), station_name (32), state (4), county (6), source_dataset (8), all S1 character arrays, dimensions station by string length.
station, obs, station_index, integer indices.
well_depth_ft, float32 per station, total well depth below land surface, feet. Valid 6,924 of 8,222 (84.2 percent), range 7.0 to 900.0 ft, mean 205.74 ft.
water_level_ft, float32 per observation, depth to water below land surface, feet, positive down. Valid 20,580 of 20,580 (100 percent), range -1.9 to 587.4 ft, mean 106.88 ft. Negative values near -1.9 ft denote artesian or measurement datum effects and warrant site-level review.
dWL_pd19_primary_ft per station (F02 deltapd_19, negative equals decline), valid 2,741 of 8,222 (33.3 percent), range -264.97 to 85.8 ft, mean -14.59 ft.
dWL_pd19_supplemental_ft per station (F03 dpd_19est), companion estimate.
dWL_1980_2019_ft per station (F03 delta80_19), valid 292 of 8,222 (3.6 percent), range -188.2 to 38.77 ft, mean -48.84 ft.
dWL_2017_2019_ft per station (F05 delta17_19), valid 7,195 of 8,222 (87.5 percent), range -27.23 to 33.68 ft, mean 0.23 ft.
lat, lon, time as above.
Quality of data: Authoritative in-situ groundwater record for the aquifer. Completeness varies by epoch by design (different campaign networks); the 1980 to 2019 change is available for only 3.6 percent of stations and must not be treated as representative without weighting. Temporal approximation for blank dates and NAD83 to WGS84 datum offset are documented and negligible at 4 km analysis scale but relevant for site-scale work. Appropriate for calibration and validation of gridded water-level change, storage and model results.

4.8 GSSURGO soils: GSSURGO_Ogallala_all_attrs_120m.nc with report.csv and iso19115_summary.json
Source datasets: USDA-NRCS gSSURGO CONUS geodatabase, raster MURASTER_30m with table Valu1, joined by mapunit key (mukey) and clipped to hp_bound2010. Grid uniform approximately 120 m Albers (overview factor 4, effective resolution 120 m, GeoTransform 120 m, origin -807735, 2314755). Coordinate system NAD83 Conus Albers, EPSG:5070. Dimensions y 11,112 by x 6,568 (72.98 million cells). Pixels inside aquifer 31,631,581. File size 894.70 MB.

Spatial and temporal extent: Static soil survey compilation. No time dimension. Geospatial bounds longitude -106.0435 to -96.2071, latitude 31.5789 to 43.8243.

Variables present (61 total):

x, y, crs (Albers Equal Area, EPSG:5070).
mukey, int32, map unit key, FillValue -9999.
Available water storage (cm water): aws0_5, aws5_20, aws20_50, aws50_100, aws100_150, aws150_999, aws0_20, aws0_30, aws0_100, aws0_150, aws0_999.
Thickness attributes (cm or as documented): tk0_5a, tk5_20a, tk20_50a, tk50_100a, tk100_150a, tk150_999a, tk0_20a, tk0_30a, tk0_100a, tk0_150a, tk0_999a and sand-associated tk0_5s, tk5_20s, tk20_50s, tk50_100s, tk100_150s, tk150_999s, tk0_20s, tk0_30s, tk0_100s, tk0_150s, tk0_999s.
Soil organic carbon stocks (g m-2 per depth interval): soc0_5, soc5_20, soc20_50, soc50_100, soc100_150, soc150_999, soc0_20, soc0_30, soc0_100, soc0_150, soc0_999.
Composition and indices (dimensionless unless noted): musumcpcta, musumcpcts, musumcpct, pctearthmc, nccpi3corn, nccpi3soy, nccpi3cot, nccpi3sg, nccpi3all (National Commodity Crop Productivity Indices), rootznemc (cm), rootznaws (cm), droughty (0/1), pwsl1pomu.
All float32 data variables use FillValue approximately 9.96921e+36. Units follow SSURGO conventions as recorded in file comments; users must consult SSURGO documentation for precise horizon definitions.

Accompanying files: GSSURGO_Ogallala_all_attrs_120m.report.csv provides per-column n_valid, minimum, median and maximum (for example aws0_5 n_valid 31,460,908, range 0 to 20.0 cm, median 7.5 cm; rootznaws n_valid 31,390,882, range 7 to 600, median 213; nccpi3all median 0.423). gssurgo_iso19115_summary.json records ISO 19115 titles, CONUS bounding box, scale denominators and abstract snippet confirming the most detailed NRCS survey level, field-verified and cartographically compiled.

Quality of data: Comprehensive and high-resolution for regional soils analysis. Valid counts of 13.7M to 31.6M per column reflect non-soil, urban, water and survey-gap masking; pwsl1pomu is sparsest (13,748,622 valid). Vintage and mapping scale vary by survey area; values are static and do not capture management-induced change. Appropriate for available water, carbon, productivity and root-zone characterization with explicit masking.

4.9 PAW_rootznaws_4km.nc
Source and derivation: gSSURGO MURASTER_30m mukey grid joined to Valu1.rootznaws, averaged onto the GRIDMET-matched 287 by 181 EPSG:4326 grid and masked to the aquifer boundary. Created by process_gssurgo_paw.py. File size 0.09 MB.

Spatial and temporal extent: 4 km grid, 287 by 181. Single time step days since 2019-07-01 (value 2019-07-01), static for analysis purposes.

Variables present:

root_zone_available_water_storage, float32, (time, y, x), FillValue approximately 9.96921e+36, long_name Plant available water stored in the soil root zone (PAW stock), units cm, cell_methods area mean, valid pixels 21,259, median samples per cell 1,502. Sampled statistics: valid 21,259 of 51,947 (40.9 percent), range 30.58 to 369.1 cm, mean 205.2 cm.
Quality of data: High aggregation robustness due to large samples per 4 km cell. Represents real PAW for soil buffering indices. Static; within-cell heterogeneity is averaged out. Suitable for direct use on the 4 km analysis grid without further resampling.


## 14. Recommendations

1. Confirm the irrigation classification legend and the AmeriFlux time units
   in writing; both are currently inferred.
2. Checksum-verify `ameriflux.nc` against the former `merged_ameriflux.nc`
   before citing continuity.
3. Treat the entirely-missing SMAP depth-to-water-table series and SMAP
   uncertainty layers as the binding constraints on decline/recovery
   analysis; no buffering attribution requiring observed water-table dynamics
   can proceed until observed well time series are secured, noting that the
   campaign-based well compilation under `USGS data/derived_usgs/` provides
   only sparse temporal snapshots.
4. Apply MODIS scale factors, SMAP/VOD quality-flag filtering, and irrigation
   class confirmation as mandatory preprocessing steps, and retain the ~4 km
   GRIDMET grid as the harmonization target with nearest-neighbour handling
   of all categorical fields.
5. Consult the GLO DEM source documentation to establish the missing-data
   sentinel (if any) before masking on zero elevation, and use windowed
   access for the 6.0 GB raster at all times.

