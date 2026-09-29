#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ecostress_et_gpp_extractor.py — ECOSTRESS ET + GPP ROI gridded pipeline (earthaccess)
=====================================================================================

Processes ECOSTRESS evapotranspiration (ET) **and** gross primary production (GPP)
over the High Plains / Ogallala aquifer polygon (``HPA_polygon/hp_bound2010.*``)
for the WHOLE mission duration via NASA Earthdata (``earthaccess``), and streams
results into **two separate CF-1.8 compliant NetCDF files** — one for ET, one for GPP.

Colab-notebooks review (10 notebooks in ``Colab notebooks/``)
--------------------------------------------------------------
* 7 notebooks use Google Earth Engine only (``ee.Authenticate`` + ``xee``/``geemap``):
  ``Cropland_&_Ag_Mgmt``, ``ET_&_Energy_Thermal_Flux`` (SSEBop/OpenET),
  ``GridMet & drought indices``, ``MODIS``, ``hydrology_&_soils``,
  ``grace_extractor``, ``GOSIF_processor`` (local GeoTIFF, no earthaccess).
* Only 3 notebooks use the NASA earthdata module (``earthaccess``):
  ``SMAP.ipynb`` (SPL4SMGP L4 soil moisture, HDF5 global EASE-2 grid, 3-hourly),
  ``SMAP_VOD.ipynb`` (SPL3SMP_E L3 enhanced, HDF5, daily AM/PM), and
  ``OCO SIF.ipynb`` (OCO2_L2_Lite_SIF point soundings).
* ``SMAP.ipynb`` is therefore the template. This script mirrors its architecture
  1:1 — login, ROI dissolve/simplify, batched search/download, per-file
  validation + targeted redownload, resume-from-output, periodic flush/sync,
  MD5-verified Drive copy, runtime checkpoints, logging/gc/cleanup — and marks
  every ECOSTRESS-specific deviation with ``[ECO-ADAPT]``.

SMAP -> ECOSTRESS adaptation notes  [ECO-ADAPT]
-----------------------------------------------
1. Products: ``SPL4SMGP`` (single HDF5 global EASE-2 grid, 3-hourly, fill -9999)
   becomes TWO tiled COG products:
     ET  = ``ECO_L3T_JET.002`` "ECOSTRESS Tiled ET Instantaneous and Daytime
           L3 Global 70 m V002" (DOI 10.5067/ECOSTRESS/ECO_L3T_JET.002,
           2018-07-10 to Present, 12 COGs per acquisition: ETdaily,
           PTJPLSMinst, STICinst, MOD16inst, BESSinst, ETinstUncertainty,
           PTJPLSMcanopy/interception/soil, STICcanopy, cloud, water).
           V1 swath product ECO3ETPTJPL.001 (HDF5 + ECO1BGEO pairing,
           swath-to-grid required) was deprecated 2025-09-30 — this script
           targets V2 tiled COGs (georeferencing embedded, no GEO pairing).
     GPP = ``ECO_L4T_WUE.002`` "ECOSTRESS Tiled Water Use Efficiency
           Instantaneous L4 Global 70 m V002" (4 COGs per acquisition: WUE,
           GPP, cloud, water). GPP is the BESS-JPL estimate
           (units umol m-2 s-1); WUE = BESS GPP / PT-JPL-SM transpiration
           (units g C kg-1 H2O). Prior to Oct-2024 the WUE field was not
           reliably populated for some scenes (see LP DAAC known issues).
   Set ``ET_VERSION``/``GPP_VERSION`` to ``"003"`` to use Collection 3 when
   available (ECO_L3T_JET.003 / ECO_L4T_WUE.003).
2. Outputs: SMAP builds ONE fixed mask on the global EASE-2 grid and appends
   every granule as one timestep. ECOSTRESS tiles are per-tile UTM (MGRS,
   109.8 km, 70 m, ~1568x1568 px) spanning many zones over the Ogallala —
   there is no single native grid, and a wall-to-wall 70 m stack of the whole
   aquifer x whole mission is infeasible (~1e8 px/timestep). Hence this script
   mosaics each overpass onto ONE shared EPSG:4326 regular grid
   (``TARGET_RES_DEG``, default 0.01 deg ~1 km; set 0.0416667 for ~4 km
   GRIDMET-match, the repo common analysis grid) and writes TWO files:
     ``ECO_L3T_JET_Ogallala_ET.nc``   (ET science bands, time,y,x)
     ``ECO_L4T_WUE_Ogallala_GPP.nc``  (GPP + WUE, time,y,x)
   Each file is standalone CF-1.8 (own time axis, lat/lon, roi_mask, crs).
   One timestep = one ISS overpass datetime (all intersecting tiles mosaicked,
   overlap-averaged). This matches ``merged_datasets/`` conventions
   (EPSG:4326, CF-1.8, polygon-masked, minimum resampling).
3. Search: SMAP searches by time only (global product). ECOSTRESS MUST add
   ``bounding_box=ROI.bounds`` (W,S,E,N), otherwise CMR pages through
   millions of global tiles.
4. Time: SMAP is regular 3-hourly (resume = last_time + 3 h). ECOSTRESS is
   irregular ISS overpasses (day+night, revisit days, gaps: Sep-Dec 2018 MSU
   anomaly, Mar-May 2019 MSU, Feb-2020 safehold, 2019-2023 3-band mode,
   noisy May-Jul 2025, etc.). Datetime is parsed from the filename
   (``..._YYYYMMDDTHHMMSS_...``); resume = skip datetimes already in output;
   batches advance by calendar days (``BATCH_DAYS``).
5. Filenames: ``SMAP_L4_SM_gph_YYYYMMDDTHHMM_V....h5`` becomes
   ``ECOv002_L3T_JET_<orbit>_<scene>_<MGRS-tile>_<YYYYMMDDTHHMMSS>_<build>_<iter>_<band>.tif``
   (and ``ECOv002_L4T_WUE_..._<band>.tif`` with bands WUE/GPP/cloud/water).
6. Validation: ``h5py`` structure check becomes lenient Tier-1 ``rasterio``
   open + CRS/shape/nearest-sample finite-pixel check (``check_cog``) with a
   logged reason. Range + constant-plane screening happens AFTER reprojection
   (``reproject_tile_to_common`` NaN-masks out-of-range values;
   ``mosaic_overpass`` drops all-NaN and flat planes) — judging the full 70 m
   tile strictly misfires on healthy tiles that are mostly fill outside the
   small ROI overlap. Redownloads re-fetch the ONE matching granule from the
   batch results (``index_results_by_filename``), never a blind time-window
   re-search (which pulls 700+ files per retry and returns the same bytes).
7. Masking: fixed EASE-2 bbox+polygon slice becomes per-overpass:
   reproject each tile COG onto the common grid (average resampling for
   float science bands, nearest for masks), accumulate sum+count for
   overlap averaging, apply mosaicked cloud/water QA mask, then apply the
   static polygon ``roi_mask`` (outside = NaN). Polygon is rasterized once
   with ``rasterio.features.rasterize``.
8. Fill: SMAP ``-9999.0`` becomes ``NaN`` (float32 COGs use NaN; uint8
   cloud/water masks use 255). NetCDF science vars use NaN _FillValue.
9. Batch sizing: SMAP ``BATCH_DATES=150-175`` d (8 files/day global).
   ECOSTRESS yields tiles x bands per overpass, so ``BATCH_DAYS=30`` default.
10. Prior work ``local data processing/ecostress_et_extractor.py`` writes a
    polygon-mean TIME SERIES (one row per COG: CSV + obs-vector NC, all bands
    mixed). This script instead writes GRIDDED mosaicked stacks with SEPARATE
    ET vs GPP files as requested. The two approaches are complementary
    (lightweight stats vs analysis-ready grids).

Usage — Google Colab (no command line, no argparse; just Run cells in order)
--------------------------------------------------------------------------
CELL 1 — install + mount Drive + stage shapefile (run once)::

    !pip install -q earthaccess rasterio netCDF4 geopandas shapely tqdm
    from google.colab import drive; drive.mount('/content/drive')
    import zipfile
    with zipfile.ZipFile('/content/drive/MyDrive/high_plains_quifer.zip') as z:
        z.extractall('/content/ogallala_shp')

CELL 2 — paste this whole file into a cell (or upload it and ``%run`` it).
         Edit the 4 settings in SETTINGS below (RUN_WHICH, versions, resolution).

CELL 3 — run (same cell, auto-runs on ``Run all``)::

    run_pipeline()

Local use is identical — just ``python ecostress_et_gpp_extractor.py`` after
editing SETTINGS. There are no CLI flags.

Outputs (in LOCAL_NC_DIR, mirrored to DRIVE_OUTPUT with MD5 verification)::

    ECO_L3T_JET_Ogallala_ET.nc    — ETdaily, PTJPLSMinst, STICinst, MOD16inst,
                                    BESSinst, ETinstUncertainty (time,y,x)
    ECO_L4T_WUE_Ogallala_GPP.nc   — GPP, WUE (time,y,x)

Re-running is safe (resume skips datetimes already stored).
"""

import os
import re
import gc
import glob
import time
import shutil
import logging
import hashlib
import warnings
from datetime import datetime, timezone, timedelta
from collections import defaultdict

import numpy as np
import netCDF4 as nc4
import geopandas as gpd

try:
    import earthaccess
except ImportError as exc:  # pragma: no cover
    raise ImportError("Please install earthaccess: pip install earthaccess") from exc

try:
    import rasterio
    from rasterio.warp import reproject, Resampling
    from rasterio.transform import from_origin
    from rasterio.features import rasterize
except ImportError as exc:  # pragma: no cover
    raise ImportError("Please install rasterio: pip install rasterio") from exc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("ecostress_et_gpp")

# ============================================================================
# SETTINGS — EDIT THESE 4 LINES ONLY, then Run all (no argparse, no CLI)
# ============================================================================

RUN_WHICH = "both"   # "both" -> ET + GPP (2 files); "ET" -> ET file only; "GPP" -> GPP file only
ET_VERSION = "002"   # set "003" for Collection 3 (ECO_L3T_JET.003)
GPP_VERSION = "002"  # set "003" for Collection 3 (ECO_L4T_WUE.003)
TARGET_RES_DEG = 0.01  # ~1 km EPSG:4326 grid; use 0.0416667 for ~4 km GRIDMET-match

# ============================================================================
# CONFIGURATION  (mirrors SMAP.ipynb config block — do not edit below)
# ============================================================================

# [ECO-ADAPT] two products instead of one SHORT_NAME
ET_SHORT_NAME = "ECO_L3T_JET"
GPP_SHORT_NAME = "ECO_L4T_WUE"

SHAPEFILE_DIR = "/content/ogallala_shp"          # same polygon as SMAP pipeline
DRIVE_OUTPUT = "/content/drive/MyDrive/MSUGWB/ECOSTRESS_gridded"
LOCAL_RAW_DIR = "/content/ecostress_raw"
LOCAL_NC_DIR = "/content/ecostress_nc"

# Local fallback: plain workstation (no /content) -> paths relative to repo.
if not os.path.isdir("/content"):
    _HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    _BASE = os.path.dirname(_HERE) if os.path.basename(_HERE).lower().startswith("local") else _HERE
    SHAPEFILE_DIR = os.path.join(_BASE, "HPA_polygon")
    DRIVE_OUTPUT = os.path.join(_HERE, "ecostress_out")
    LOCAL_RAW_DIR = os.path.join(_HERE, "ecostress_raw")
    LOCAL_NC_DIR = os.path.join(_HERE, "ecostress_nc")

ROI_LABEL = "Ogallala"

# [ECO-ADAPT] science-band subsets (keep small for whole-duration runs).
# Full JET set also includes PTJPLSMcanopy/interception/soil + STICcanopy;
# add them to ET_BANDS if partition analysis is needed.
ET_BANDS = [
    "ETdaily",
    "PTJPLSMinst",
    "STICinst",
    "MOD16inst",
    "BESSinst",
    "ETinstUncertainty",
]
GPP_BANDS = ["GPP", "WUE"]
QA_BANDS = ("cloud", "water")  # sibling COGs used for masking only

ET_META = {
    # Units per LP DAAC ECO_L3T_JET v002 user guide / product spec.
    "ETdaily": {"long_name": "Daily evapotranspiration (daytime total)", "units": "mm day-1"},
    "PTJPLSMinst": {"long_name": "Instantaneous ET, PT-JPL-SM model", "units": "W m-2"},
    "STICinst": {"long_name": "Instantaneous ET, STIC model", "units": "W m-2"},
    "MOD16inst": {"long_name": "Instantaneous ET, MOD16 model", "units": "W m-2"},
    "BESSinst": {"long_name": "Instantaneous ET, BESS model", "units": "W m-2"},
    "ETinstUncertainty": {"long_name": "Instantaneous ET uncertainty", "units": "W m-2"},
}
GPP_META = {
    # Units confirmed in ECO_L3T/L4T v002/v003 user guides, Table 6.
    "GPP": {"long_name": "Gross primary production, BESS-JPL", "units": "umol m-2 s-1"},
    "WUE": {"long_name": "Water use efficiency (BESS GPP / PT-JPL-SM transpiration)", "units": "g C kg-1 H2O"},
}

# Plausible physical ranges used for corruption / gibberish screening.
# Values outside these bounds are treated as corrupt and set to NaN (never written).
# Bounds are intentionally generous so real extremes survive; obvious fill-errors
# (e.g. 1e20, -1e10, constant-plane artifacts) are caught.
BAND_VALID_RANGE = {
    "ETdaily": (-1.0, 30.0),             # mm day-1 (small negatives = retrieval noise)
    "PTJPLSMinst": (-50.0, 1200.0),      # W m-2 (night can be slightly negative)
    "STICinst": (-50.0, 1200.0),
    "MOD16inst": (-50.0, 1200.0),
    "BESSinst": (-50.0, 1200.0),
    "ETinstUncertainty": (0.0, 1000.0),
    "GPP": (0.0, 100.0),                 # umol m-2 s-1
    "WUE": (0.0, 30.0),                  # g C kg-1 H2O
}
# A reprojected plane is dropped as a flat artifact only when it has at least
# this many finite pixels sharing one value (tiny footprints can be legitimately uniform).
PLANE_CONST_MIN_PX = 200
# A plane is rejected as gibberish when >99.5% of its finite pixels share one value
# (flat-plane decode error) — real ET/GPP fields always vary.
GIBBERISH_CONST_FRAC = 0.995

ET_OUTPUT = "ECO_L3T_JET_Ogallala_ET.nc"
GPP_OUTPUT = "ECO_L4T_WUE_Ogallala_GPP.nc"

MISSION_START = datetime(2018, 7, 10)  # [ECO-ADAPT] ECO_L3T_JET temporal extent start
TIME_UNITS = "seconds since 2000-01-01 00:00:00 UTC"  # CF-compliant, second precision for ISS times
TIME_CALENDAR = "proleptic_gregorian"

ROI_SIMPLIFY_DEG = 0.01
APPLY_QA_MASK = True       # mask science pixels where mosaicked cloud/water == 1
QA_COVER_FRAC = 0.5        # mosaicked QA >= 0.5 counts as cloudy/water

BATCH_DAYS = 30            # [ECO-ADAPT] small: tiles x bands per overpass
FLUSH_EVERY_STEPS = 10     # nc sync + gc every N overpass timesteps
BACKUP_EVERY_STEPS = 50    # MD5-verified Drive backup every N timesteps
CHECKPOINT_HOURS = 1
DOWNLOAD_THREADS = 6
MAX_RETRY = 3              # [ECO-ADAPT] mirrors SMAP MAX_RETRIES redownload logic
COMPRESS_LEVEL = 3

# ECOv002_L3T_JET_40665_017_51KVB_20250907T091510_0713_01_STICinst.tif
# ECOv002_L4T_WUE_40665_017_14TML_20220801T143210_0710_01_GPP.tif
ECO_FILENAME_RE = re.compile(
    r"ECOv\d+_(L[34]T_(?:JET|WUE))_(\d+)_(\d+)_([A-Z0-9]+)_(\d{8}T\d{6})_\d+_\d+_([A-Za-z0-9]+)\.tif$"
)


# ============================================================================
# FILENAME PARSING  (SMAP parse_dt analog)
# ============================================================================

def parse_ecostress_filename(filename):
    """Extract (product, orbit, scene, tile, datetime, band) from a V2/V3 tiled COG name.

    Returns dict or None if the pattern does not match (caller skips with warning).
    """
    m = ECO_FILENAME_RE.search(os.path.basename(filename))
    if not m:
        return None
    product, orbit, scene, tile, dt_str, band = m.groups()
    try:
        dt = datetime.strptime(dt_str, "%Y%m%dT%H%M%S")
    except ValueError:
        return None
    return {"product": product, "orbit": orbit, "scene": scene,
            "tile": tile, "datetime": dt, "band": band}


def timestep_key(dt):
    """Canonical per-overpass key (one NetCDF timestep per ISS overpass)."""
    return dt.strftime("%Y%m%dT%H%M%S")


# ============================================================================
# ROI + COMMON GRID  (SMAP load_roi/build_mask analogs)
# ============================================================================

def load_roi_geometry(shapefile_dir, simplify_deg=ROI_SIMPLIFY_DEG):
    """Load HPA polygon, reproject to EPSG:4326, dissolve, simplify."""
    shp_files = glob.glob(os.path.join(shapefile_dir, "**", "*.shp"), recursive=True)
    if not shp_files:
        raise FileNotFoundError("No .shp found in {}".format(shapefile_dir))
    log.info("Shapefile: %s", shp_files[0])
    gdf = gpd.read_file(shp_files[0])
    if gdf.crs is None or gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(epsg=4326)
    roi = (gdf.geometry.union_all() if hasattr(gdf.geometry, "union_all")
           else gdf.geometry.unary_union)
    roi_simple = roi.simplify(simplify_deg, preserve_topology=True)
    log.info("ROI bounds (EPSG:4326): %s",
             tuple(round(v, 4) for v in roi_simple.bounds))
    return roi_simple


def define_common_grid(roi_4326, res=TARGET_RES_DEG):
    """Snap ROI bbox to a regular EPSG:4326 grid; return (lon1d, lat1d, transform).

    Latitude is north-first (descending) to match rasterio/GDAL convention.
    """
    minx, miny, maxx, maxy = roi_4326.bounds
    minx = float(np.floor(minx / res) * res)
    miny = float(np.floor(miny / res) * res)
    maxx = float(np.ceil(maxx / res) * res)
    maxy = float(np.ceil(maxy / res) * res)
    nx = int(round((maxx - minx) / res)) + 1
    ny = int(round((maxy - miny) / res)) + 1
    lon = minx + np.arange(nx, dtype=np.float64) * res
    lat = maxy - np.arange(ny, dtype=np.float64) * res
    transform = from_origin(minx - res / 2.0, maxy + res / 2.0, res, res)
    log.info("Common grid: %d (y) x %d (x) at %.5f deg | lon %.4f..%.4f lat %.4f..%.4f",
             ny, nx, res, lon[0], lon[-1], lat[-1], lat[0])
    return lon, lat, transform


def rasterize_roi_mask(roi_4326, lon, lat, transform):
    """Burn polygon onto the common grid (1=inside, 0=outside)."""
    mask = rasterize(
        [(g, 1) for g in getattr(roi_4326, "geoms", [roi_4326])],
        out_shape=(len(lat), len(lon)),
        transform=transform,
        fill=0,
        dtype="uint8",
        all_touched=True,
    ).astype(np.int8)
    log.info("ROI mask: %d inside / %d total (%.1f %%)",
             int(mask.sum()), mask.size, 100.0 * mask.sum() / mask.size)
    return mask


# ============================================================================
# CF-1.8 NETCDF CREATION  (SMAP create_nc analog; one file per product)
# ============================================================================

def _stamp_global_attrs(ds, short_name, version, bands, lon, lat, res, references):
    ds.Conventions = "CF-1.8"
    ds.title = "{} v{} — {} ROI gridded subset ({} deg)".format(
        short_name, version, ROI_LABEL, res)
    ds.source = "NASA ECOSTRESS {} v{} via earthaccess (LP DAAC)".format(short_name, version)
    ds.institution = "NASA JPL / LP DAAC"
    ds.history = "Created {}".format(datetime.now(timezone.utc).isoformat())
    ds.references = references
    ds.comment = ("Per-overpass mosaics of 70 m MGRS tiles reprojected onto a shared "
                  "EPSG:4326 grid (average resampling), cloud/water + polygon masked, "
                  "NaN fill. One timestep per ISS overpass datetime (irregular). "
                  "Bands: {}.".format(", ".join(bands)))
    ds.roi = ROI_LABEL
    ds.geospatial_lat_min = float(lat.min())
    ds.geospatial_lat_max = float(lat.max())
    ds.geospatial_lon_min = float(lon.min())
    ds.geospatial_lon_max = float(lon.max())
    ds.geospatial_bounds_crs = "EPSG:4326"
    ds.grid_resolution = "{} degrees".format(res)


def create_output_nc(path, short_name, version, bands, band_meta, lon, lat, roi_mask,
                      references, res=TARGET_RES_DEG):
    """Create a fresh CF-1.8 gridded file with unlimited time axis."""
    if os.path.exists(path):
        return
    ny, nx = len(lat), len(lon)
    ds = nc4.Dataset(path, "w", format="NETCDF4")
    ds.createDimension("time", None)
    ds.createDimension("y", ny)
    ds.createDimension("x", nx)

    t = ds.createVariable("time", "f8", ("time",))
    t.units = TIME_UNITS
    t.calendar = TIME_CALENDAR
    t.axis = "T"
    t.standard_name = "time"
    t.long_name = "ECOSTRESS overpass time (from filename, UTC)"

    y = ds.createVariable("lat", "f4", ("y",))
    y.units = "degrees_north"
    y.standard_name = "latitude"
    y.axis = "Y"
    y.long_name = "latitude of grid cell center"
    y[:] = lat.astype(np.float32)

    x = ds.createVariable("lon", "f4", ("x",))
    x.units = "degrees_east"
    x.standard_name = "longitude"
    x.axis = "X"
    x.long_name = "longitude of grid cell center"
    x[:] = lon.astype(np.float32)

    m = ds.createVariable("roi_mask", "i1", ("y", "x"), zlib=True, complevel=1)
    m.long_name = "Region of Interest mask (1=inside {}, 0=outside)".format(ROI_LABEL)
    m.flag_values = np.array([0, 1], dtype="i1")
    m.flag_meanings = "outside_roi inside_roi"
    m[:] = roi_mask

    chunk = (1, min(512, ny), min(512, nx))
    for band in bands:
        meta = band_meta.get(band, {})
        v = ds.createVariable(band, "f4", ("time", "y", "x"),
                              fill_value=np.float32(np.nan),
                              zlib=True, complevel=COMPRESS_LEVEL,
                              chunksizes=chunk)
        v.long_name = meta.get("long_name", band)
        v.units = meta.get("units", "1")
        v.coordinates = "lat lon"
        v.grid_mapping = "crs"
        v.cell_methods = "time: point (ISS overpass mosaic)"

    crs = ds.createVariable("crs", "i4")
    crs.grid_mapping_name = "latitude_longitude"
    crs.semi_major_axis = 6378137.0
    crs.inverse_flattening = 298.257223563
    crs.crs_wkt = ('GEOGCRS["WGS 84",DATUM["World Geodetic System 1984",'
                   'ELLIPSOID["WGS 84",6378137,298.257223563]],CS[ellipsoidal,2],'
                   'AXIS["lat",north],AXIS["lon",east],UNIT["degree",0.0174532925199433]]')

    _stamp_global_attrs(ds, short_name, version, bands, lon, lat, res, references)
    ds.close()
    log.info("Created %s (%d bands, %d x %d grid)", os.path.basename(path), len(bands), ny, nx)


def get_resume_datetimes(path):
    """Return set of timestep_keys already stored + max datetime (or empty/None)."""
    keys, max_dt = set(), None
    if not os.path.exists(path):
        return keys, max_dt
    try:
        with nc4.Dataset(path, "r") as ds:
            if len(ds.variables["time"]) == 0:
                return keys, max_dt
            units = ds.variables["time"].units
            cal = getattr(ds.variables["time"], "calendar", TIME_CALENDAR)
            vals = ds.variables["time"][:]
            for v in vals:
                dt = nc4.num2date(v, units, cal)
                dt = datetime(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
                keys.add(timestep_key(dt))
                if max_dt is None or dt > max_dt:
                    max_dt = dt
        log.info("Resume %s: %d timesteps (max %s)", os.path.basename(path), len(keys), max_dt)
    except Exception as e:
        log.warning("Cannot read resume point from %s: %s", path, e)
    return keys, max_dt


# ============================================================================
# SEARCH + DOWNLOAD  (SMAP pattern + [ECO-ADAPT] bounding_box filter)
# ============================================================================

def search_granules(short_name, version, start_dt, end_dt, bbox):
    """CMR search restricted to ROI bbox — REQUIRED for tiled ECOSTRESS.

    bbox = (minx, miny, maxx, maxy) == (W, S, E, N).
    """
    results = earthaccess.search_data(
        short_name=short_name,
        version=version,
        bounding_box=(bbox[0], bbox[1], bbox[2], bbox[3]),
        temporal=(start_dt.strftime("%Y-%m-%dT%H:%M:%SZ"),
                  end_dt.strftime("%Y-%m-%dT%H:%M:%SZ")),
    )
    log.info("Search %s %s: %s -> %s returned %d granules",
             short_name, version, start_dt.date(), end_dt.date(), len(results))
    return results


def download_granules(results, local_dir, threads=DOWNLOAD_THREADS):
    os.makedirs(local_dir, exist_ok=True)
    files = earthaccess.download(results, local_path=local_dir, threads=threads)
    return [str(f) for f in files]


def index_results_by_filename(results):
    """Map remote basename -> granule for single-file redownloads.

    Lets a bad-file retry re-fetch exactly ONE granule instead of re-searching
    a time window (which re-downloads the whole batch — hundreds of files per
    retry — and usually returns the same bad bytes, because earthaccess skips
    files already present locally).
    """
    index = {}
    for g in results:
        try:
            links = g.data_links()
        except Exception:
            continue
        for url in links or []:
            base = os.path.basename(url.split("?")[0])
            if base.endswith(".tif"):
                index.setdefault(base, g)
    return index


def check_cog(path):
    """Tier-1 integrity check: is this file worth reprojecting? Returns (ok, reason).

    Verifies the file is a readable single-band GeoTIFF with a CRS and non-zero
    dims, and that a strided NEAREST-neighbour sample holds at least a handful
    of finite (non-fill) pixels. (SMAP validate_h5 analog.)

    Deliberately LENIENT about content: physical-range and constant-plane
    screening happen AFTER reprojection onto the common grid (see
    reproject_tile_to_common / mosaic_overpass), because a healthy tile can be
    mostly fill outside our small ROI overlap — judging the full 70 m tile
    strictly is what previously flagged whole overpasses "corrupt". Truncated
    downloads fail here with an "unreadable" reason and are the files actually
    worth redownloading.
    """
    try:
        with rasterio.open(path) as src:
            if src.count < 1 or src.width < 1 or src.height < 1:
                return False, "empty dims"
            if src.crs is None:
                return False, "missing CRS"
            try:
                _ = src.read(1, window=((0, min(2, src.height)), (0, min(2, src.width))))
            except Exception as e:
                return False, "unreadable (%s)" % e
            # Strided sample for the finite-pixel check. MUST be nearest: the
            # default resampling smears NaN across the sample and makes healthy
            # but sparse tiles look all-NaN (mass false "corrupt").
            h, w = src.height, src.width
            stride_y = max(1, h // 500)
            stride_x = max(1, w // 500)
            sample = src.read(1,
                              out_shape=(max(1, h // stride_y), max(1, w // stride_x)),
                              resampling=Resampling.nearest)
            sample = sample.astype(np.float32)
            nodata = src.nodata
            info = parse_ecostress_filename(path)
            band = info["band"] if info else None
            if band in QA_BANDS:
                sample[sample == 255] = np.nan
            else:
                if nodata is not None and not (isinstance(nodata, float) and np.isnan(nodata)):
                    sample[sample == nodata] = np.nan
                sample[~np.isfinite(sample)] = np.nan
            n_fin = int(np.isfinite(sample).sum())
            if n_fin == 0:
                return False, "empty (0 finite pixels in sample)"
            if n_fin < 8:
                return False, "empty (only %d finite px in sample)" % n_fin
        return True, ""
    except Exception as e:
        return False, "unreadable (%s)" % e


def validate_cog(path):
    """Bool wrapper around check_cog (kept for compatibility)."""
    ok, _ = check_cog(path)
    return ok


# ============================================================================
# MOSAICKING  [ECO-ADAPT core: many UTM tiles -> one EPSG:4326 timestep]
# ============================================================================

def reproject_tile_to_common(cog_path, dst_transform, dst_width, dst_height, is_mask=False,
                             band=None):
    """Reproject one tile COG onto the full common grid; return float32 array.

    Applies corruption screening: nodata -> NaN, out-of-range science values
    -> NaN (per BAND_VALID_RANGE), strict {0,1} enforcement for QA masks.
    All-NaN sources skip the warp and return a NaN plane directly.
    """
    if band is None:
        info = parse_ecostress_filename(cog_path)
        band = info["band"] if info else None
    with rasterio.open(cog_path) as src:
        arr = src.read(1).astype(np.float32)
        nodata = src.nodata
        if is_mask or band in QA_BANDS:
            arr[arr == 255] = np.nan
            if nodata is not None and nodata != 255 and not (isinstance(nodata, float) and np.isnan(nodata)):
                arr[arr == nodata] = np.nan
            # Strict mask: only 0/1 survive; gibberish mask codes -> NaN
            bad = np.isfinite(arr) & (arr != 0) & (arr != 1)
            arr[bad] = np.nan
        else:
            if nodata is not None and not (isinstance(nodata, float) and np.isnan(nodata)):
                arr[arr == nodata] = np.nan
            arr[~np.isfinite(arr)] = np.nan
            if band in BAND_VALID_RANGE:
                vmin, vmax = BAND_VALID_RANGE[band]
                arr[(arr < vmin) | (arr > vmax)] = np.nan
        if np.isfinite(arr).sum() == 0:
            return np.full((dst_height, dst_width), np.nan, dtype=np.float32)
        dst = np.full((dst_height, dst_width), np.nan, dtype=np.float32)
        reproject(arr, dst,
                  src_transform=src.transform, src_crs=src.crs,
                  dst_transform=dst_transform, dst_crs="EPSG:4326",
                  resampling=Resampling.nearest if (is_mask or band in QA_BANDS) else Resampling.average,
                  dst_nodata=np.nan)
    # Post-warp screen: a tile that missed the common grid is all-NaN -> keep as NaN plane
    return dst


def mosaic_overpass(file_list, science_bands, dst_transform, dst_width, dst_height,
                     roi_mask, apply_qa=APPLY_QA_MASK):
    """Mosaic all tile COGs sharing one overpass datetime.

    Returns {band: 2D mosaicked array (NaN outside footprint/ROI)} or {} if empty.
    Overlaps are averaged (NaN-aware sum/count). QA siblings (cloud/water) are
    mosaicked with max (conservative: cloudy if any tile says cloudy).
    """
    by_band = defaultdict(list)
    for f in file_list:
        info = parse_ecostress_filename(f)
        if info is None:
            continue
        by_band[info["band"]].append(f)

    # Mosaic QA first (nearest-neighbor, max combine)
    qa_mosaic = None
    if apply_qa:
        qa_parts = []
        for qb in QA_BANDS:
            for f in by_band.get(qb, []):
                try:
                    plane = reproject_tile_to_common(
                        f, dst_transform, dst_width, dst_height, is_mask=True)
                    if np.isfinite(plane).sum() == 0:
                        continue  # tile missed grid or was all-fill -> ignore
                    qa_parts.append(plane)
                except Exception as e:
                    log.warning("QA reproject failed %s: %s", os.path.basename(f), e)
        if qa_parts:
            stack = np.stack(qa_parts, axis=0)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                with np.errstate(invalid="ignore"):
                    qa_mosaic = np.nanmax(stack, axis=0)  # 1=bad anywhere
            del stack, qa_parts
            gc.collect()

    out = {}
    for band in science_bands:
        parts = []
        for f in by_band.get(band, []):
            try:
                plane = reproject_tile_to_common(
                    f, dst_transform, dst_width, dst_height, is_mask=False)
                fin = plane[np.isfinite(plane)]
                if fin.size == 0:
                    continue  # corrupt/empty/off-footprint tile -> ignore
                # Flat-artifact screen (post-reproject, where it belongs: full-tile
                # screening cannot tell fill-outside-ROI from real corruption).
                if fin.size >= PLANE_CONST_MIN_PX:
                    vals, counts = np.unique(fin, return_counts=True)
                    if counts.max() / fin.size > GIBBERISH_CONST_FRAC:
                        log.warning("Dropping constant-plane tile %s (v=%s)",
                                    os.path.basename(f), vals[np.argmax(counts)])
                        continue
                parts.append(plane)
            except Exception as e:
                log.warning("Reproject failed %s: %s", os.path.basename(f), e)
        if not parts:
            continue
        stack = np.stack(parts, axis=0)
        # Guard: every column all-NaN -> nanmean warns "Mean of empty slice".
        # Suppress the warning AND verify afterwards that anything survived.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            with np.errstate(invalid="ignore"):
                mosaic = np.nanmean(stack, axis=0).astype(np.float32)
        del stack, parts
        gc.collect()
        if qa_mosaic is not None and mosaic.shape == qa_mosaic.shape:
            mosaic[qa_mosaic >= QA_COVER_FRAC] = np.nan
        mosaic[roi_mask == 0] = np.nan  # static polygon mask (repo convention)
        if np.all(~np.isfinite(mosaic)):
            log.info("Band %s: no valid pixels after masking (skipped)", band)
            continue
        out[band] = mosaic
    return out


def append_timestep(nc_path, dt, band_arrays):
    """Append one overpass timestep. Skips all-NaN bands so empty cells never fill the file.

    Returns True if at least one band was written, False if the timestep was
    entirely empty (caller marks it done WITHOUT appending — no gibberish/NaN
    timesteps accumulate in the outputs).
    """
    clean = {}
    for band, arr in band_arrays.items():
        a = np.asarray(arr, dtype=np.float32)
        if a.size == 0 or np.isfinite(a).sum() == 0:
            log.info("Band %s: all-NaN, not written for %s", band, timestep_key(dt))
            continue
        clean[band] = a
    if not clean:
        return False
    ds = nc4.Dataset(nc_path, "a")
    try:
        t_idx = len(ds.variables["time"])
        ds.variables["time"][t_idx] = nc4.date2num(dt, TIME_UNITS, TIME_CALENDAR)
        for band, arr in clean.items():
            if band in ds.variables:
                ds.variables[band][t_idx, :, :] = arr
            else:
                log.warning("Band %s not in %s — skipped", band, os.path.basename(nc_path))
        ds.sync()
    finally:
        ds.close()
    return True


# ============================================================================
# DRIVE COPY  (identical to SMAP robust_drive_copy)
# ============================================================================

def get_md5(filepath, chunk_bytes=10 * 1024 * 1024):
    h = hashlib.md5()
    with open(filepath, "rb") as f:
        while True:
            chunk = f.read(chunk_bytes)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def robust_drive_copy(src, dst, max_retries=3):
    src_md5 = get_md5(src)
    for attempt in range(1, max_retries + 1):
        try:
            os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
            with open(src, "rb") as fs, open(dst, "wb") as fd:
                shutil.copyfileobj(fs, fd, length=10 * 1024 * 1024)
                fd.flush()
                os.fsync(fd.fileno())
            time.sleep(2)
            if os.path.getsize(src) != os.path.getsize(dst):
                raise IOError("Size mismatch")
            if src_md5 != get_md5(dst):
                raise IOError("MD5 mismatch")
            log.info("Drive copy verified (%.1f MiB): %s",
                     os.path.getsize(dst) / 1024 ** 2, os.path.basename(dst))
            return True
        except Exception as e:
            log.error("Drive copy attempt %d failed: %s", attempt, e)
            if os.path.exists(dst):
                try:
                    os.remove(dst)
                except OSError:
                    pass
            if attempt < max_retries:
                time.sleep(10 * attempt)
            else:
                raise RuntimeError("Drive copy failed: {}".format(e))


# ============================================================================
# MAIN PIPELINE (one product -> one CF file; run twice for ET + GPP)
# ============================================================================

def run_product(short_name, version, science_bands, band_meta, output_name,
                references, roi, bbox, lon, lat, transform, roi_mask,
                raw_root, nc_path):
    """Batched whole-mission download -> mosaic -> append loop for ONE product."""
    raw_dir = os.path.join(raw_root, short_name)
    os.makedirs(raw_dir, exist_ok=True)

    # Resume: restore from Drive if missing locally (SMAP pattern)
    drive_nc = os.path.join(DRIVE_OUTPUT, output_name)
    if not os.path.exists(nc_path) and os.path.exists(drive_nc):
        log.info("Restoring %s from Drive for resume...", output_name)
        shutil.copy2(drive_nc, nc_path)
    if not os.path.exists(nc_path):
        create_output_nc(nc_path, short_name, version, science_bands, band_meta,
                          lon, lat, roi_mask, references)

    done_keys, max_dt = get_resume_datetimes(nc_path)
    start_dt = (max_dt + timedelta(seconds=1)) if max_dt else MISSION_START
    log.info("Product %s v%s start: %s (mission %s -> now)",
             short_name, version, start_dt, MISSION_START.date())

    ny, nx = len(lat), len(lon)
    script_t0 = time.time()
    next_checkpoint = script_t0 + CHECKPOINT_HOURS * 3600
    n_checkpoint = 0
    steps_since_backup = 0
    total_steps = 0
    today = datetime.now()
    current = start_dt

    while current < today:
        batch_end = min(current + timedelta(days=BATCH_DAYS), today)
        log.info("=" * 50)
        log.info("BATCH %s %s: %s -> %s", short_name, version, current.date(), batch_end.date())
        log.info("=" * 50)

        try:
            results = search_granules(short_name, version, current, batch_end, bbox)
        except Exception as e:
            log.error("Search failed: %s", e)
            current = batch_end
            continue
        if not results:
            current = batch_end
            continue

        try:
            downloaded = download_granules(results, raw_dir)
        except Exception as e:
            log.error("Download failed: %s", e)
            current = batch_end
            continue

        # Keep only target + QA bands; group by overpass datetime
        groups = defaultdict(list)
        for f in downloaded:
            if not f.endswith(".tif"):
                continue
            info = parse_ecostress_filename(f)
            if info is None:
                log.warning("Cannot parse filename, skipping: %s", os.path.basename(f))
                continue
            if info["band"] not in set(science_bands) | set(QA_BANDS):
                continue
            groups[info["datetime"]].append(f)

        new_dts = sorted(dt for dt in groups if timestep_key(dt) not in done_keys)
        log.info("%d overpasses (%d new) in batch.", len(groups), len(new_dts))
        granule_index = index_results_by_filename(results)

        for dt in new_dts:
            files = groups[dt]
            # Per-file integrity check with SINGLE-granule redownload.
            valid_files = []
            for f in files:
                ok, reason = check_cog(f)
                if ok:
                    valid_files.append(f)
                    continue
                base = os.path.basename(f)
                log.warning("Bad COG (%s), redownloading single granule: %s", reason, base)
                try:
                    os.remove(f)
                except OSError:
                    pass
                granule = granule_index.get(base)
                if granule is None:
                    # Never blind re-search a time window here: it pulls the
                    # whole batch again (700+ files per retry in the logs).
                    log.error("No granule match for %s in batch results — skipping.", base)
                    continue
                # "empty" content is usually real (commissioning fill / night /
                # off-swath tile), not a transfer error: one confirmation
                # re-fetch, then skip. "unreadable" gets full retries.
                max_attempts = 1 if reason.startswith("empty") else MAX_RETRY
                redone = False
                for attempt in range(1, max_attempts + 1):
                    try:
                        rfiles = earthaccess.download([granule], local_path=raw_dir, threads=1)
                        got = None
                        for rf in rfiles:
                            rf = str(rf)
                            if os.path.basename(rf) == base:
                                got = rf
                                break
                        if got is None:
                            log.error("Redownload attempt %d: no file returned for %s",
                                      attempt, base)
                            continue
                        ok2, reason2 = check_cog(got)
                        if ok2:
                            valid_files.append(got)
                            redone = True
                            log.info("Redownload OK: %s", base)
                            break
                        log.warning("Redownload attempt %d still bad (%s): %s",
                                    attempt, reason2, base)
                        try:
                            os.remove(got)
                        except OSError:
                            pass
                    except Exception as e:
                        log.error("Redownload attempt %d failed: %s", attempt, e)
                    time.sleep(5 * attempt)
                if not redone:
                    log.error("Skipping bad file after %d retries: %s", max_attempts, f)

            if not valid_files:
                continue
            try:
                mosaics = mosaic_overpass(valid_files, science_bands, transform, nx, ny,
                                          roi_mask, apply_qa=APPLY_QA_MASK)
            except Exception as e:
                log.error("Mosaic failed for %s: %s", dt, e)
                mosaics = {}
            for f in valid_files:
                try:
                    if os.path.exists(f):
                        os.remove(f)
                except OSError:
                    pass
            if not mosaics:
                # Tiles were all corrupt/empty/off-footprint: mark done WITHOUT
                # appending so resume never retries it, but no NaN timestep pollutes the file.
                log.info("No valid pixels for %s — marking done, nothing appended.",
                         timestep_key(dt))
                done_keys.add(timestep_key(dt))
                continue

            try:
                wrote = append_timestep(nc_path, dt, mosaics)
                # Mark done either way: re-downloading an empty overpass forever
                # is worse than skipping it; append_timestep already refused to
                # write all-NaN bands so no gibberish enters the file.
                done_keys.add(timestep_key(dt))
                if wrote:
                    total_steps += 1
                    steps_since_backup += 1
                    log.info("Appended %s (%s) bands=%s", dt.strftime("%Y-%m-%dT%H:%M:%SZ"),
                             timestep_key(dt), sorted(mosaics.keys()))
                else:
                    log.info("Empty mosaic for %s — marking done, nothing appended.",
                             timestep_key(dt))
            except Exception as e:
                log.error("Append failed for %s: %s", dt, e)
            finally:
                del mosaics
                gc.collect()

            if total_steps and total_steps % FLUSH_EVERY_STEPS == 0:
                gc.collect()
            if steps_since_backup >= BACKUP_EVERY_STEPS:
                try:
                    robust_drive_copy(nc_path, drive_nc)
                except Exception as e:
                    log.error("Backup copy failed: %s", e)
                steps_since_backup = 0
            if time.time() >= next_checkpoint:
                n_checkpoint += 1
                cp = os.path.join(DRIVE_OUTPUT, "ECOSTRESS_{}_{}hr.nc".format(
                    short_name, n_checkpoint * CHECKPOINT_HOURS))
                try:
                    with nc4.Dataset(nc_path, "a") as _ds:
                        _ds.sync()
                    robust_drive_copy(nc_path, cp)
                except Exception as e:
                    log.error("Checkpoint copy failed: %s", e)
                next_checkpoint = time.time() + CHECKPOINT_HOURS * 3600

        for f in glob.glob(os.path.join(raw_dir, "*.tif")):
            try:
                os.remove(f)
            except OSError:
                pass
        current = batch_end
        gc.collect()

    log.info("Backing up %s to Drive...", output_name)
    try:
        robust_drive_copy(nc_path, drive_nc)
    except Exception as e:
        log.error("Final backup failed: %s", e)
    log.info("Product %s complete — %d new timesteps this run.", short_name, total_steps)
    return total_steps


def run_pipeline(products=None):
    """Run the whole-mission pipeline. products=None honours SETTINGS RUN_WHICH."""
    if products is None:
        products = _wanted_products()
    log.info("=" * 64)
    log.info("  ECOSTRESS gridded ET+GPP extraction — %s (products=%s)", ROI_LABEL, products)
    log.info("=" * 64)

    earthaccess.login(strategy="interactive")
    roi = load_roi_geometry(SHAPEFILE_DIR)
    bbox = roi.bounds  # (minx, miny, maxx, maxy) == (W, S, E, N)

    for d in (LOCAL_RAW_DIR, LOCAL_NC_DIR, DRIVE_OUTPUT):
        os.makedirs(d, exist_ok=True)

    lon, lat, transform = define_common_grid(roi, res=TARGET_RES_DEG)
    roi_mask = rasterize_roi_mask(roi, lon, lat, transform)

    results = {}
    if "ET" in products:
        results["ET"] = run_product(
            ET_SHORT_NAME, ET_VERSION, ET_BANDS, ET_META, ET_OUTPUT,
            "https://doi.org/10.5067/ECOSTRESS/ECO_L3T_JET.002",
            roi, bbox, lon, lat, transform, roi_mask,
            LOCAL_RAW_DIR, os.path.join(LOCAL_NC_DIR, ET_OUTPUT))
    if "GPP" in products:
        results["GPP"] = run_product(
            GPP_SHORT_NAME, GPP_VERSION, GPP_BANDS, GPP_META, GPP_OUTPUT,
            "https://doi.org/10.5067/ECOSTRESS/ECO_L4T_WUE.002",
            roi, bbox, lon, lat, transform, roi_mask,
            LOCAL_RAW_DIR, os.path.join(LOCAL_NC_DIR, GPP_OUTPUT))

    log.info("=" * 64)
    log.info("  PIPELINE COMPLETE — %s", results)
    log.info("  ET : %s", os.path.join(DRIVE_OUTPUT, ET_OUTPUT))
    log.info("  GPP: %s", os.path.join(DRIVE_OUTPUT, GPP_OUTPUT))
    log.info("=" * 64)
    return results


# ============================================================================
# RUN — Colab: just press Run all. Local: python ecostress_et_gpp_extractor.py
# Settings come from SETTINGS at the top (RUN_WHICH / versions / resolution).
# ============================================================================

def _wanted_products():
    key = str(RUN_WHICH).strip().lower()
    if key == "both":
        return ("ET", "GPP")
    if key in ("et", "gpp"):
        return (key.upper(),)
    raise ValueError('RUN_WHICH must be "both", "ET" or "GPP" — got %r' % RUN_WHICH)


if __name__ == "__main__":
    # Running as a script (Colab "Run all" or python file.py) launches the pipeline.
    # In a notebook, calling run_pipeline() manually does the same thing.
    run_pipeline(products=_wanted_products())
else:
    # Imported as a module: do NOT auto-run; user calls run_pipeline() in next cell.
    pass
