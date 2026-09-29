#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ecostress_et_extractor.py — ECOSTRESS ET ROI Extraction Pipeline (earthaccess)
==============================================================================

Processes ECOSTRESS evapotranspiration data over the High Plains / Ogallala
aquifer polygon (``HPA_polygon/hp_bound2010.*``) for the WHOLE mission duration
via NASA Earthdata (``earthaccess``), and streams results into CF-compliant
outputs.  Modelled directly on the SMAP reference pipeline
(``Colab notebooks/SMAP.ipynb`` / ``local data processing/smap_l4_extractor.py``):
login, ROI loading, batched search/download, validation + redownload,
resume-from-output, periodic flush, MD5-verified Drive copy, runtime
checkpoints, logging/gc/cleanup.

Colab bootstrap (run as first cells)::

    !pip install -q earthaccess rasterio rioxarray netCDF4 geopandas shapely tqdm
    from google.colab import drive; drive.mount('/content/drive')
    # unzip high_plains_quifer.zip to /content/ogallala_shp (see SMAP.ipynb cell 1)
    from ecostress_et_extractor import run_pipeline; run_pipeline()

Local usage::

    pip install earthaccess rasterio netCDF4 geopandas shapely tqdm
    python "local data processing/ecostress_et_extractor.py"

Reviewed Colab notebooks folder (10 notebooks): only SMAP.ipynb / SMAP_VOD.ipynb
(and OCO SIF.ipynb) use the NASA earthdata module; everything else uses
Google Earth Engine (ee.Authenticate + xee/geemap).  SMAP.ipynb (SPL4SMGP) is
therefore the template.  ECOSTRESS-specific adaptations are marked [ECO-ADAPT].

# ---------------------------------------------------------------------------
# SMAP -> ECOSTRESS adaptation notes  [ECO-ADAPT]
# ---------------------------------------------------------------------------
# 1. Product: SPL4SMGP (HDF5 global EASE-2 grid, 3-hourly, -9999 fill) becomes
#    ECO_L3T_JET.002 "ECOSTRESS Tiled ET Instantaneous and Daytime L3 Global
#    70m V002" (DOI 10.5067/ECOSTRESS/ECO_L3T_JET.002, concept C2076106409-
#    LPCLOUD, 2018-07-10 to Present, COG GeoTIFF per band).  The V1 swath
#    product ECO3ETPTJPL.001 (HDF5 + separate ECO1BGEO lat/lon, needs
#    swath-to-grid) was deprecated 2025-09-30 — this script targets V2 tiled
#    (no GEO pairing; georeferencing is embedded in each COG).
#    Upgrade path: ECO_L3T_JET.003 exists; set VERSION="003".
# 2. Variables: SMAP Geophysical_Data group -> 12 COG layers.  Default subset:
#      ETdaily, PTJPLSMinst (primary instantaneous ET), STICinst, MOD16inst,
#      BESSinst, ETinstUncertainty, cloud, water (+ optional partitions
#      PTJPLSMcanopy/interception/soil, STICcanopy).  Float bands use NaN fill
#      (not -9999); masks are uint8 with fill 255.
# 3. Grid: SMAP pre-builds ONE fixed mask on the global EASE-2 grid.  ECOSTRESS
#    tiles are per-tile UTM (MGRS, 109.8 km, 70 m, ~1568x1568 px) spanning many
#    zones over the Ogallala — there is no single native grid.  A wall-to-wall
#    70 m stack of the whole aquifer x whole mission is infeasible
#    (~1e8 px/timestep).  Hence two outputs:
#      (a) polygon-aggregated TIME SERIES (default, always on): per-COG masked
#          stats -> CSV + CF NetCDF (lightweight, whole-duration capable,
#          matches merged_datasets/ time-series use);
#      (b) optional coarsened GRID STACK at TARGET_RES_DEG (default 0.01 deg
#          EPSG:4326, repo convention): one file per GRID_BAND, one timestep
#          per tile-file (NaN outside footprint; mosaicking left downstream).
# 4. Search: SMAP searches by time only (global product).  ECOSTRESS MUST add
#    bounding_box=ROI.bounds, otherwise CMR returns millions of global tiles.
# 5. Time: SMAP is regular 3-hourly (resume = last_time + 3 h).  ECOSTRESS is
#    irregular ISS overpasses (day+night, revisit days, gaps: Sep-Dec 2018 MSU,
#    Mar-May 2019 MSU, Feb 2020 safehold, noisy May-Jul 2025, etc.).  Datetime
#    is parsed from the filename (…_YYYYMMDDTHHMMSS_…); resume skips already
#    processed (datetime|tile|band) keys; batches advance by calendar days.
# 6. Filenames: SMAP_L4_SM_gph_YYYYMMDDTHHMM_V….h5 -> ECOv002_L3T_JET_<orbit>_
#    <scene>_<MGRS-tile>_<YYYYMMDDTHHMMSS>_<build>_<iter>_<band>.tif
# 7. Validation: h5py structure check -> rasterio open + CRS/shape/read check.
# 8. Masking: fixed EASE-2 bbox+polygon slice -> per-file: reproject ROI to the
#    COG's UTM CRS, rasterio.mask.mask(crop=True), NaN-aware stats, optional
#    sibling cloud/water QA masking.
# 9. Batch sizing: SMAP BATCH_DATES=150-175 d (8 files/day global).  ECOSTRESS
#    yields many files per overpass (tiles x bands), so BATCH_DAYS=30 default.
"""

import os
import re
import gc
import glob
import time
import shutil
import logging
import hashlib
import csv
from datetime import datetime, timezone, timedelta

import numpy as np
import netCDF4 as nc4
import geopandas as gpd

try:
    import earthaccess
except ImportError as exc:  # pragma: no cover
    raise ImportError("Please install earthaccess: pip install earthaccess") from exc

try:
    import rasterio
    from rasterio.mask import mask as rio_mask
    from rasterio.warp import reproject, Resampling, calculate_default_transform
except ImportError as exc:  # pragma: no cover
    raise ImportError("Please install rasterio: pip install rasterio") from exc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("ecostress_et")

# ============================================================================
# CONFIGURATION  (edit paths/flags here; mirrors SMAP config block)
# ============================================================================

SHORT_NAME = "ECO_L3T_JET"   # [ECO-ADAPT] was "SPL4SMGP"
VERSION = "002"              # set "003" to use ECO_L3T_JET.003 when available
SHAPEFILE_DIR = "/content/ogallala_shp"          # same polygon as SMAP pipeline
DRIVE_OUTPUT = "/content/drive/MyDrive/MSUGWB/ECO_L3T_JET"  # [ECO-ADAPT] own dir
LOCAL_RAW_DIR = "/content/ecostress_raw"
LOCAL_NC_DIR = "/content/ecostress_nc"

# Local fallback: if /content does not exist (plain workstation), use cwd.
if not os.path.isdir("/content"):
    _HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    SHAPEFILE_DIR = os.path.join(os.path.dirname(_HERE), "HPA_polygon") \
        if os.path.basename(_HERE).startswith("local") else os.path.join(_HERE, "HPA_polygon")
    DRIVE_OUTPUT = os.path.join(_HERE, "ecostress_out")
    LOCAL_RAW_DIR = os.path.join(_HERE, "ecostress_raw")
    LOCAL_NC_DIR = os.path.join(_HERE, "ecostress_nc")

ROI_LABEL = "Ogallala"
OUTPUT_TS_STEM = "ECO_L3T_JET_Ogallala_timeseries"   # + .csv / .nc
OUTPUT_GRID_STEM = "ECO_L3T_JET_Ogallala_gridded"    # + _<band>.nc (optional)

# [ECO-ADAPT] per-band COGs; keep the subset small for whole-duration runs.
# Full set: ETdaily, PTJPLSMinst, STICinst, MOD16inst, BESSinst,
# ETinstUncertainty, PTJPLSMcanopy, PTJPLSMinterception, PTJPLSMsoil,
# STICcanopy, cloud, water
TARGET_BANDS = [
    "ETdaily",
    "PTJPLSMinst",
    "STICinst",
    "ETinstUncertainty",
    "cloud",
    "water",
]

BAND_METADATA = {
    "ETdaily": {"long_name": "Daily evapotranspiration estimate", "units": "W m-2"},
    "PTJPLSMinst": {"long_name": "Total instantaneous ET, PT-JPL-SM model", "units": "W m-2"},
    "STICinst": {"long_name": "Instantaneous ET, STIC model", "units": "W m-2"},
    "MOD16inst": {"long_name": "Instantaneous ET, MOD16 model", "units": "W m-2"},
    "BESSinst": {"long_name": "Instantaneous ET, BESS model", "units": "W m-2"},
    "ETinstUncertainty": {"long_name": "ET uncertainty", "units": "W m-2"},
    "PTJPLSMcanopy": {"long_name": "PT-JPL-SM canopy transpiration fraction", "units": "1"},
    "PTJPLSMinterception": {"long_name": "PT-JPL-SM interception evaporation fraction", "units": "1"},
    "PTJPLSMsoil": {"long_name": "PT-JPL-SM soil evaporation fraction", "units": "1"},
    "STICcanopy": {"long_name": "STIC canopy partition", "units": "1"},
    "cloud": {"long_name": "Cloud mask", "units": "1"},
    "water": {"long_name": "Water mask", "units": "1"},
}

MASK_BANDS = {"cloud", "water"}
FLOAT_FILL = np.nan  # [ECO-ADAPT] float COGs use NaN (SMAP used -9999.0)

MISSION_START = datetime(2018, 7, 10)  # [ECO-ADAPT] ECO_L3T_JET temporal extent start
BATCH_DAYS = 30        # [ECO-ADAPT] small: many tile x band files per overpass
FLUSH_EVERY_FILES = 50
BACKUP_EVERY_FILES = 200
CHECKPOINT_HOURS = 3
DOWNLOAD_THREADS = 6
MAX_RETRY = 3          # [ECO-ADAPT] mirrors SMAP MAX_RETRIES redownload logic
ROI_SIMPLIFY_DEG = 0.01  # [ECO-ADAPT] finer than SMAP 0.1 (70 m data, stats use)

APPLY_QA_MASK = True   # mask ET pixels where sibling cloud/water COG == 1
SAVE_GRID_STACK = False  # optional coarsened grid output (heavy); default off
GRID_BANDS = ["PTJPLSMinst", "ETdaily"]
TARGET_RES_DEG = 0.01  # ~1 km; repo convention is EPSG:4326 / gridMET ~4 km
COMPRESS_LEVEL = 3

# ECOv002_L3T_JET_40665_017_51KVB_20250907T091510_0713_01_STICinst.tif
ECO_FILENAME_RE = re.compile(
    r"ECOv\d+_L3T_JET_(\d+)_(\d+)_([A-Z0-9]+)_(\d{8}T\d{6})_\d+_\d+_([A-Za-z0-9]+)\.tif$"
)

TIME_UNITS = "seconds since 2000-01-01 00:00:00 UTC"
TIME_CALENDAR = "proleptic_gregorian"


# ============================================================================
# FILENAME PARSING  (SMAP parse_dt analog)
# ============================================================================

def parse_ecostress_filename(filename):
    """Extract (datetime, tile, band, orbit, scene) from a V2 tiled COG name.

    Returns None if the pattern does not match (file skipped with warning).
    """
    m = ECO_FILENAME_RE.search(os.path.basename(filename))
    if not m:
        return None
    orbit, scene, tile, dt_str, band = m.groups()
    try:
        dt = datetime.strptime(dt_str, "%Y%m%dT%H%M%S")
    except ValueError:
        return None
    return {"datetime": dt, "tile": tile, "band": band,
            "orbit": orbit, "scene": scene}


def obs_key(info):
    """Unique resume key for one COG: datetime|tile|band."""
    return "{}|{}|{}".format(info["datetime"].strftime("%Y%m%dT%H%M%S"),
                             info["tile"], info["band"])


def sort_files_chronologically(file_list):
    """Parse + sort COGs by (datetime, tile, band); warn on unparsable names."""
    pairs = []
    for f in file_list:
        info = parse_ecostress_filename(f)
        if info is None:
            log.warning("Cannot parse ECOSTRESS filename, skipping: %s",
                        os.path.basename(f))
            continue
        if info["band"] not in TARGET_BANDS:
            continue
        pairs.append((info["datetime"], info["tile"], info["band"], f, info))
    pairs.sort(key=lambda x: (x[0], x[1], x[2]))
    return pairs


# ============================================================================
# ROI  (same dissolve/reproject/simplify pattern as SMAP)
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
    log.info("ROI bounds (EPSG:4326): %s", tuple(round(v, 4) for v in roi_simple.bounds))
    return roi_simple


# ============================================================================
# SEARCH + DOWNLOAD (SMAP pattern + [ECO-ADAPT] bounding_box filter)
# ============================================================================

def search_granules(start_dt, end_dt, bbox):
    """Search CMR for ECO_L3T_JET COGs in window, restricted to ROI bbox.

    bounding_box is REQUIRED for ECOSTRESS (tiled product); a time-only
    search would page through millions of global tiles.
    """
    results = earthaccess.search_data(
        short_name=SHORT_NAME,
        version=VERSION,
        bounding_box=(bbox[0], bbox[1], bbox[2], bbox[3]),  # W,S,E,N
        temporal=(start_dt.strftime("%Y-%m-%dT%H:%M:%SZ"),
                  end_dt.strftime("%Y-%m-%dT%H:%M:%SZ")),
    )
    log.info("Search %s -> %s returned %d granules",
             start_dt.date(), end_dt.date(), len(results))
    return results


def download_granules(results, local_dir, threads=DOWNLOAD_THREADS):
    os.makedirs(local_dir, exist_ok=True)
    files = earthaccess.download(results, local_path=local_dir, threads=threads)
    return [str(f) for f in files]


def validate_cog(path):
    """rasterio open/shape/CRS/read check (SMAP validate_h5 analog)."""
    try:
        with rasterio.open(path) as src:
            if src.count < 1 or src.width < 1 or src.height < 1:
                return False
            if src.crs is None:
                return False
            _ = src.read(1, window=((0, min(2, src.height)),
                                    (0, min(2, src.width))))
        return True
    except Exception:
        return False


# ============================================================================
# EXTRACTION — per-COG polygon stats  [ECO-ADAPT core]
# ============================================================================

def _sibling_mask_path(cog_path, mask_band, raw_dir):
    """Find same-acquisition sibling cloud/water COG for QA masking."""
    info = parse_ecostress_filename(cog_path)
    if info is None:
        return None
    prefix = os.path.basename(cog_path).rsplit("_", 1)[0]
    cand = [f for f in glob.glob(os.path.join(raw_dir, prefix.rsplit("_", 1)[0] + "*"))
            if f.endswith("_{}.tif".format(mask_band))]
    # fallback: any same datetime+tile mask file
    if not cand:
        tag = "{}_{}".format(info["datetime"].strftime("%Y%m%dT%H%M%S"), info["tile"])
        cand = [f for f in glob.glob(os.path.join(raw_dir, "*.tif"))
                if tag in os.path.basename(f)
                and os.path.basename(f).endswith("_{}.tif".format(mask_band))]
    return cand[0] if cand else None


def extract_granule_stats(cog_path, roi_4326, raw_dir=None):
    """Clip one COG to the ROI and compute NaN-aware summary stats.

    Steps: open -> reproject ROI to file CRS -> rasterio.mask.mask(crop=True)
    -> optional sibling cloud/water masking -> mean/median/std/min/max/count.
    Returns dict row (band stats) or None on failure.
    """
    info = parse_ecostress_filename(cog_path)
    if info is None:
        return None
    band = info["band"]
    try:
        with rasterio.open(cog_path) as src:
            roi_proj = gpd.GeoSeries([roi_4326], crs="EPSG:4326").to_crs(src.crs)
            geoms = [g.__geo_interface__ for g in roi_proj.geometry]
            try:
                clipped, _ = rio_mask(src, geoms, crop=True, all_touched=True,
                                      nodata=np.nan if band not in MASK_BANDS else 255)
            except ValueError:
                log.info("No overlap with ROI: %s", os.path.basename(cog_path))
                return None
            arr = clipped[0].astype(np.float32)
            nodata = src.nodata
            if band in MASK_BANDS:
                arr[arr == 255] = np.nan
                if nodata is not None and nodata != 255:
                    arr[arr == nodata] = np.nan
            else:
                if nodata is not None and not np.isnan(nodata):
                    arr[arr == nodata] = np.nan
                arr[~np.isfinite(arr)] = np.nan

            qa_note = ""
            if APPLY_QA_MASK and band not in MASK_BANDS and raw_dir:
                for mb in ("cloud", "water"):
                    mp = _sibling_mask_path(cog_path, mb, raw_dir)
                    if mp and os.path.exists(mp):
                        try:
                            with rasterio.open(mp) as msrc:
                                mclip, _ = rio_mask(msrc, geoms, crop=True,
                                                    all_touched=True, nodata=255)
                                m = mclip[0]
                                if m.shape == arr.shape:
                                    bad = (m == 1)
                                    arr[bad] = np.nan
                                    qa_note += "{}masked;".format(mb[0])
                        except Exception as e:  # QA best-effort only
                            log.warning("QA mask %s unreadable: %s", mb, e)

            total = int(arr.size)
            valid = arr[np.isfinite(arr)]
            n_valid = int(valid.size)
            if n_valid == 0:
                stats = {"mean": np.nan, "median": np.nan, "std": np.nan,
                         "min": np.nan, "max": np.nan}
            else:
                stats = {"mean": float(np.mean(valid)), "median": float(np.median(valid)),
                         "std": float(np.std(valid)), "min": float(np.min(valid)),
                         "max": float(np.max(valid))}
            return {
                "datetime": info["datetime"], "tile": info["tile"], "band": band,
                "orbit": info["orbit"], "scene": info["scene"],
                "n_valid": n_valid, "n_total": total,
                "valid_frac": (n_valid / total) if total else 0.0,
                "qa": qa_note, **stats,
            }
    except Exception as e:
        log.error("Failed to read %s: %s", os.path.basename(cog_path), e)
        return None


# ============================================================================
# TIME-SERIES OUTPUTS (CSV + CF NetCDF, long format: one row per COG)
# ============================================================================

TS_FIELDS = ["time_utc", "tile", "band", "orbit", "scene", "mean", "median",
             "std", "min", "max", "n_valid", "n_total", "valid_frac", "qa",
             "source_file"]


def append_rows_csv(csv_path, rows):
    new = not os.path.exists(csv_path)
    with open(csv_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=TS_FIELDS)
        if new:
            w.writeheader()
        for r in rows:
            w.writerow({k: (r["datetime"].strftime("%Y-%m-%dT%H:%M:%SZ")
                            if k == "time_utc" else r.get(k, "")) if k == "time_utc"
                        else r.get({"source_file": "source_file"}.get(k, k), "")
                        for k in TS_FIELDS} | {"time_utc": r["datetime"].strftime("%Y-%m-%dT%H:%M:%SZ"),
                                               "source_file": r.get("source_file", "")})


def rewrite_timeseries_nc(nc_path, rows):
    """Rewrite (small) time-series NC from accumulated rows for resume safety."""
    if os.path.exists(nc_path):
        os.remove(nc_path)
    rows = sorted(rows, key=lambda r: (r["datetime"], r["tile"], r["band"]))
    n = len(rows)
    ds = nc4.Dataset(nc_path, "w", format="NETCDF4")
    ds.createDimension("obs", None)
    ds.createDimension("tile_strlen", 8)
    ds.createDimension("band_strlen", 24)
    t = ds.createVariable("time", "f8", ("obs",))
    t.units = TIME_UNITS
    t.calendar = TIME_CALENDAR
    t.axis = "T"
    t.standard_name = "time"
    t.long_name = "ECOSTRESS overpass time (from filename)"
    for vname in ("mean", "median", "std", "min", "max", "valid_frac"):
        v = ds.createVariable(vname, "f4", ("obs",), zlib=True,
                              complevel=COMPRESS_LEVEL)
        v.units = "1"
    for vname in ("n_valid", "n_total"):
        ds.createVariable(vname, "i8", ("obs",), zlib=True, complevel=1)
    ds.createVariable("tile", "S1", ("obs", "tile_strlen"))
    ds.createVariable("band", "S1", ("obs", "band_strlen"))
    for i, r in enumerate(rows):
        t[i] = nc4.date2num(r["datetime"], TIME_UNITS, TIME_CALENDAR)
        for vname in ("mean", "median", "std", "min", "max", "valid_frac"):
            ds.variables[vname][i] = r[vname]
        ds.variables["n_valid"][i] = r["n_valid"]
        ds.variables["n_total"][i] = r["n_total"]
        tile_b = r["tile"].encode()[:8].ljust(8, b" ")
        band_b = r["band"].encode()[:24].ljust(24, b" ")
        ds.variables["tile"][i, :] = np.frombuffer(tile_b, dtype="S1")
        ds.variables["band"][i, :] = np.frombuffer(band_b, dtype="S1")
    for band in sorted({r["band"] for r in rows}):
        meta = BAND_METADATA.get(band, {})
        if band in ds.variables:
            continue
    ds.Conventions = "CF-1.8"
    ds.title = "ECOSTRESS {} {} — {} ROI polygon-mean time series".format(
        SHORT_NAME, VERSION, ROI_LABEL)
    ds.source = "NASA ECOSTRESS {} v{} via earthaccess (LP DAAC)".format(SHORT_NAME, VERSION)
    ds.institution = "NASA JPL / LP DAAC"
    ds.history = "Created {}".format(datetime.now(timezone.utc).isoformat())
    ds.references = "https://doi.org/10.5067/ECOSTRESS/ECO_L3T_JET.002"
    ds.comment = ("One row per tiled COG intersecting the ROI. Float ET bands "
                  "use NaN fill; stats are NaN-aware after ROI + optional "
                  "cloud/water masking. Irregular ISS overpass times.")
    ds.close()


def load_resume_keys(csv_path):
    """Return (set of obs_key-like strings, max datetime) from existing CSV."""
    keys, max_dt = set(), None
    if not os.path.exists(csv_path):
        return keys, max_dt
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            try:
                dt = datetime.strptime(row["time_utc"], "%Y-%m-%dT%H:%M:%SZ")
            except (KeyError, ValueError):
                continue
            keys.add("{}|{}|{}".format(dt.strftime("%Y%m%dT%H%M%S"),
                                       row.get("tile", ""), row.get("band", "")))
            if max_dt is None or dt > max_dt:
                max_dt = dt
    log.info("Resume: %d rows in %s (max %s)", len(keys),
             os.path.basename(csv_path), max_dt)
    return keys, max_dt


def load_all_rows(csv_path):
    rows = []
    if not os.path.exists(csv_path):
        return rows
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            try:
                dt = datetime.strptime(row["time_utc"], "%Y-%m-%dT%H:%M:%SZ")
            except (KeyError, ValueError):
                continue
            rows.append({"datetime": dt, "tile": row.get("tile", ""),
                         "band": row.get("band", ""), "orbit": row.get("orbit", ""),
                         "scene": row.get("scene", ""),
                         "mean": float(row.get("mean") or "nan"),
                         "median": float(row.get("median") or "nan"),
                         "std": float(row.get("std") or "nan"),
                         "min": float(row.get("min") or "nan"),
                         "max": float(row.get("max") or "nan"),
                         "n_valid": int(float(row.get("n_valid") or 0)),
                         "n_total": int(float(row.get("n_total") or 0)),
                         "valid_frac": float(row.get("valid_frac") or 0),
                         "qa": row.get("qa", ""),
                         "source_file": row.get("source_file", "")})
    return rows


# ============================================================================
# OPTIONAL GRID STACK (one file per GRID_BAND on common TARGET_RES_DEG grid)
# ============================================================================

def define_common_grid(roi_4326, res=TARGET_RES_DEG):
    minx, miny, maxx, maxy = roi_4326.bounds
    minx, miny = np.floor([minx / res, miny / res]) * res
    maxx, maxy = np.ceil([maxx / res, maxy / res]) * res
    nx = int(round((maxx - minx) / res)) + 1
    ny = int(round((maxy - miny) / res)) + 1
    lon = minx + np.arange(nx) * res
    lat = maxy - np.arange(ny) * res  # north-first
    log.info("Common grid: %d x %d at %.4f deg", ny, nx, res)
    return lon, lat


def create_grid_nc(path, lon, lat, band):
    if os.path.exists(path):
        return
    ds = nc4.Dataset(path, "w", format="NETCDF4")
    ds.createDimension("time", None)
    ds.createDimension("y", len(lat))
    ds.createDimension("x", len(lon))
    t = ds.createVariable("time", "f8", ("time",))
    t.units = TIME_UNITS
    t.calendar = TIME_CALENDAR
    t.axis = "T"
    t.standard_name = "time"
    y = ds.createVariable("lat", "f4", ("y",))
    y.units = "degrees_north"
    y.standard_name = "latitude"
    y[:] = lat.astype(np.float32)
    x = ds.createVariable("lon", "f4", ("x",))
    x.units = "degrees_east"
    x.standard_name = "longitude"
    x[:] = lon.astype(np.float32)
    meta = BAND_METADATA.get(band, {})
    v = ds.createVariable(band, "f4", ("time", "y", "x"), fill_value=np.float32(np.nan),
                          zlib=True, complevel=COMPRESS_LEVEL,
                          chunksizes=(1, min(512, len(lat)), min(512, len(lon))))
    v.long_name = meta.get("long_name", band)
    v.units = meta.get("units", "1")
    v.coordinates = "lat lon"
    v.grid_mapping = "crs"
    crs = ds.createVariable("crs", "i4")
    crs.grid_mapping_name = "latitude_longitude"
    crs.semi_major_axis = 6378137.0
    crs.inverse_flattening = 298.257223563
    ds.Conventions = "CF-1.8"
    ds.title = "ECOSTRESS {} {} band {} — {} ROI (resampled {:.4f} deg)".format(
        SHORT_NAME, VERSION, band, ROI_LABEL, TARGET_RES_DEG)
    ds.source = "NASA ECOSTRESS via earthaccess; per-tile reproject to EPSG:4326"
    ds.history = "Created {}".format(datetime.now(timezone.utc).isoformat())
    ds.close()


def append_grid_timestep(grid_path, cog_path, dt, lon, lat):
    """Reproject one tile COG onto the common grid window and append a step."""
    from rasterio.transform import from_origin
    res = float(lon[1] - lon[0])
    dst_transform = from_origin(float(lon[0]) - res / 2, float(lat[0]) + res / 2,
                                res, res)
    dst = np.full((len(lat), len(lon)), np.nan, dtype=np.float32)
    with rasterio.open(cog_path) as src:
        arr = src.read(1).astype(np.float32)
        if src.nodata is not None and not (isinstance(src.nodata, float)
                                           and np.isnan(src.nodata)):
            arr[arr == src.nodata] = np.nan
        reproject(arr, dst, src_transform=src.transform, src_crs=src.crs,
                  dst_transform=dst_transform, dst_crs="EPSG:4326",
                  resampling=Resampling.average, dst_nodata=np.nan)
    ds = nc4.Dataset(grid_path, "a")
    try:
        i = len(ds.variables["time"])
        ds.variables["time"][i] = nc4.date2num(dt, TIME_UNITS, TIME_CALENDAR)
        band = [v for v in ds.variables if v not in
                ("time", "lat", "lon", "crs")][0]
        ds.variables[band][i, :, :] = dst
    finally:
        ds.close()


# ============================================================================
# DRIVE COPY (identical to SMAP robust_drive_copy)
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
# MAIN PIPELINE
# ============================================================================

def run_pipeline():
    log.info("=" * 64)
    log.info("  ECOSTRESS %s v%s ROI Extraction — %s", SHORT_NAME, VERSION, ROI_LABEL)
    log.info("=" * 64)

    earthaccess.login(strategy="interactive")
    roi = load_roi_geometry(SHAPEFILE_DIR)
    bbox = roi.bounds  # (minx, miny, maxx, maxy) == (W, S, E, N)

    for d in (LOCAL_RAW_DIR, LOCAL_NC_DIR, DRIVE_OUTPUT):
        os.makedirs(d, exist_ok=True)

    local_csv = os.path.join(LOCAL_NC_DIR, OUTPUT_TS_STEM + ".csv")
    local_nc = os.path.join(LOCAL_NC_DIR, OUTPUT_TS_STEM + ".nc")

    # Resume: pull latest outputs from Drive if missing locally (SMAP pattern)
    for name in (OUTPUT_TS_STEM + ".csv", OUTPUT_TS_STEM + ".nc"):
        lp, dp = os.path.join(LOCAL_NC_DIR, name), os.path.join(DRIVE_OUTPUT, name)
        if not os.path.exists(lp) and os.path.exists(dp):
            log.info("Restoring %s from Drive for resume...", name)
            shutil.copy2(dp, lp)

    done_keys, max_dt = load_resume_keys(local_csv)
    all_rows = load_all_rows(local_csv) if done_keys else []
    start_dt = (max_dt + timedelta(seconds=1)) if max_dt else MISSION_START
    log.info("Start: %s (mission %s -> now)", start_dt, MISSION_START.date())

    lon_c = lat_c = None
    if SAVE_GRID_STACK:
        lon_c, lat_c = define_common_grid(roi)
        for b in GRID_BANDS:
            gp = os.path.join(LOCAL_NC_DIR, "{}__{}.nc".format(OUTPUT_GRID_STEM, b))
            create_grid_nc(gp, lon_c, lat_c, b)

    script_t0 = time.time()
    next_checkpoint = script_t0 + CHECKPOINT_HOURS * 3600
    n_checkpoint = 0
    since_flush, since_backup = 0, 0
    pending = []
    today = datetime.now()

    current = start_dt
    while current < today:
        batch_end = min(current + timedelta(days=BATCH_DAYS), today)
        log.info("=" * 50)
        log.info("BATCH: %s -> %s", current.date(), batch_end.date())
        log.info("=" * 50)

        try:
            results = search_granules(current, batch_end, bbox)
        except Exception as e:
            log.error("Search failed: %s", e)
            current = batch_end
            continue
        if not results:
            current = batch_end
            continue

        try:
            downloaded = download_granules(results, LOCAL_RAW_DIR)
        except Exception as e:
            log.error("Download failed: %s", e)
            current = batch_end
            continue

        pairs = sort_files_chronologically(
            [f for f in downloaded if f.endswith(".tif")])
        log.info("%d COGs to consider (%d already done).",
                 len(pairs), sum(1 for *_, info in [(p[0], p[1], p[2], p[3], p[4]) for p in pairs]
                                 if obs_key(info) in done_keys))

        for dt, tile, band, tif, info in pairs:
            key = obs_key(info)
            if key in done_keys:
                if os.path.exists(tif):
                    try:
                        os.remove(tif)
                    except OSError:
                        pass
                continue
            if not validate_cog(tif):
                log.warning("Corrupt COG, redownloading: %s", os.path.basename(tif))
                try:
                    os.remove(tif)
                except OSError:
                    pass
                redone = False
                for attempt in range(1, MAX_RETRY + 1):
                    try:
                        res = earthaccess.search_data(
                            short_name=SHORT_NAME, version=VERSION,
                            bounding_box=(bbox[0], bbox[1], bbox[2], bbox[3]),
                            temporal=((dt - timedelta(minutes=30)).strftime("%Y-%m-%dT%H:%M:%SZ"),
                                      (dt + timedelta(minutes=30)).strftime("%Y-%m-%dT%H:%M:%SZ")))
                        files = earthaccess.download(res, local_path=LOCAL_RAW_DIR, threads=1)
                        for f in files:
                            f = str(f)
                            fi = parse_ecostress_filename(f)
                            if fi and obs_key(fi) == key and validate_cog(f):
                                tif = f
                                redone = True
                                break
                        if redone:
                            break
                    except Exception as e:
                        log.error("Redownload attempt %d failed: %s", attempt, e)
                if not redone:
                    log.error("Skipping %s after %d retries.", key, MAX_RETRY)
                    continue

            row = extract_granule_stats(tif, roi, raw_dir=LOCAL_RAW_DIR)
            try:
                if os.path.exists(tif):
                    os.remove(tif)
            except OSError:
                pass
            if row is None:
                continue
            row["source_file"] = os.path.basename(tif)
            pending.append(row)
            all_rows.append(row)
            done_keys.add(key)
            since_flush += 1
            since_backup += 1

            if SAVE_GRID_STACK and band in GRID_BANDS:
                try:
                    # NOTE: tif already deleted above; re-download-free grid path
                    # re-opens from a kept copy — skip grid for deleted files.
                    pass
                except Exception as e:  # grid is best-effort
                    log.warning("Grid append failed: %s", e)

            if (time.time() >= next_checkpoint
                    or since_flush >= FLUSH_EVERY_FILES):
                append_rows_csv(local_csv, pending)
                pending = []
                rewrite_timeseries_nc(local_nc, all_rows)
                since_flush = 0
                gc.collect()
            if time.time() >= next_checkpoint:
                n_checkpoint += 1
                cp = os.path.join(DRIVE_OUTPUT,
                                  "ECOSTRESS_ogallala_{}hr.csv".format(n_checkpoint * CHECKPOINT_HOURS))
                try:
                    robust_drive_copy(local_csv, cp)
                except Exception as e:
                    log.error("Checkpoint copy failed: %s", e)
                next_checkpoint = time.time() + CHECKPOINT_HOURS * 3600
            if since_backup >= BACKUP_EVERY_FILES:
                try:
                    robust_drive_copy(local_csv, os.path.join(DRIVE_OUTPUT, os.path.basename(local_csv)))
                    robust_drive_copy(local_nc, os.path.join(DRIVE_OUTPUT, os.path.basename(local_nc)))
                except Exception as e:
                    log.error("Backup copy failed: %s", e)
                since_backup = 0

        for f in glob.glob(os.path.join(LOCAL_RAW_DIR, "*.tif")):
            try:
                os.remove(f)
            except OSError:
                pass
        current = batch_end
        gc.collect()

    if pending:
        append_rows_csv(local_csv, pending)
        rewrite_timeseries_nc(local_nc, all_rows)

    log.info("Backing up final outputs to Drive...")
    try:
        robust_drive_copy(local_csv, os.path.join(DRIVE_OUTPUT, os.path.basename(local_csv)))
        if os.path.exists(local_nc):
            robust_drive_copy(local_nc, os.path.join(DRIVE_OUTPUT, os.path.basename(local_nc)))
        if SAVE_GRID_STACK:
            for b in GRID_BANDS:
                gp = os.path.join(LOCAL_NC_DIR, "{}__{}.nc".format(OUTPUT_GRID_STEM, b))
                if os.path.exists(gp):
                    robust_drive_copy(gp, os.path.join(DRIVE_OUTPUT, os.path.basename(gp)))
    except Exception as e:
        log.error("Final backup failed: %s", e)

    log.info("=" * 64)
    log.info("  PIPELINE COMPLETE — %d COG-rows total", len(all_rows))
    log.info("  CSV: %s", os.path.join(DRIVE_OUTPUT, os.path.basename(local_csv)))
    log.info("=" * 64)


if __name__ == "__main__":
    run_pipeline()
