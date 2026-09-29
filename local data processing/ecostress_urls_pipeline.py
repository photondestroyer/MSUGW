#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ECHO-ET (ECOSTRESS ET) URL-list pipeline — chunked download / ROI filter / CF NetCDF.

NO earthaccess, NO login. Data come from the public UNH directory
(plain HTTP, h5ai listing)::

    https://data.globalecology.unh.edu/data/ECHO-ET/Local_solar_time/YYYY/MM/ECHO-ET_hourly_LST_YYYYMMDDHH.tif

Layout verified 2026: ``Local_solar_time/<year>/<month>/`` folders holding
hourly global GeoTIFFs (``ECHO-ET_hourly_LST_YYYYMMDDHH.tif``, ~17 MB each,
24/day, ~44,000 files over 2018-2022). One file = one hourly timestep, so no
tile-mosaicking is needed: each file is windowed to the ROI bbox, masked to
the ``.shp`` polygon, and appended to the ``.nc`` time axis.

What it does
------------
1. Crawls the site (h5ai JSON + HTML fallback) and writes ``urls.txt`` (every
   file found) and ``tif_urls.txt`` (``.tif`` files only, ~44k lines). Both
   are copied to Google Drive. The raw listing is cached as JSON on Drive so
   re-runs skip the crawl (plus HEAD size lookup, one-time cost).
2. Downloads 10,000 files at a time into scratch space (resumable ``.part``
   files via HTTP Range). Scratch free space is checked before every download
   sub-batch and at least every 30 minutes: below 40 GB, downloading pauses
   and the pipeline processes what is already on scratch.
3. Verifies every raw file (TIFF magic bytes, no HTML error page, expected
   size, rasterio open + sample read, CRS/dims/finite-pixel check). Failures
   are deleted and redownloaded individually (up to MAX_RETRIES, backoff);
   every corrupt event is appended to ``corrupted_files.txt``.
4. Filters each file to the ``.shp`` polygon and appends one hourly timestep
   to a CF-1.8 ``.nc`` file written with compression level 1.
5. Deletes the raw files from scratch, appends their names to the
   processed-files log, and repeats the cycle for the next 10,000.
6. Every 1 hour the ``.nc`` file (plus the url lists and logs) is copied to
   Google Drive as backup.
7. Resume: on restart the processed-files log is restored from Drive and any
   basename already logged is never downloaded again; the ``.nc`` time axis
   is a second guard against duplicate timesteps.

Units warning
-------------
Filenames carry the tag ``LST`` but the product folder is an ET product.
Because the server documents no units, the ``.nc`` variable stores raw values
with ``units = "1"`` and a flagging comment — VERIFY the physical units
before analysis. At startup the pipeline logs sample min/max/mean so you can
tell at a glance (Kelvin ~250-330  vs  mm hr-1 ~0-2). If the values prove to
be Kelvin LST, set VAR_ATTRS below to
``{"long_name": "Land surface temperature", "units": "K"}``.

Colab use (no argparse, no CLI — just Run cells in order)
----------------------------------------------------------
CELL 1 — setup (run once)::

    !pip install -q netCDF4 rasterio geopandas shapely tqdm beautifulsoup4 requests
    from google.colab import drive; drive.mount('/content/drive')
    import zipfile
    with zipfile.ZipFile('/content/drive/MyDrive/high_plains_quifer.zip') as z:
        z.extractall('/content/ogallala_shp')

CELL 2 — paste this whole file into a cell (or upload it and ``%run`` it),
         edit SETTINGS below, then Run all. ``run_pipeline()`` launches.
         Set TEST_MODE = True first: only 100 files are processed.

Local use is identical: ``python echo_et_urls_pipeline.py`` after editing
SETTINGS. There are no command-line flags.
"""

import os
import re
import gc
import glob
import json
import time
import shutil
import logging
import hashlib
import warnings
from datetime import datetime, timezone, timedelta
from collections import defaultdict
from urllib.parse import urljoin, urlsplit, quote, unquote

import requests

try:
    import numpy as np
except ImportError as exc:
    raise ImportError("Please install numpy: pip install numpy") from exc
import netCDF4 as nc4
import geopandas as gpd

try:
    import rasterio
    from rasterio.warp import reproject, Resampling
    from rasterio.windows import from_bounds
    from rasterio.transform import from_origin
    from rasterio.features import rasterize
except ImportError as exc:
    raise ImportError("Please install rasterio: pip install rasterio") from exc

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None  # progress bars become plain loops

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("echo_et")

# ============================================================================
# SETTINGS — EDIT THESE, then Run all (no argparse, no CLI)
# ============================================================================

TEST_MODE = True    # True -> only TEST_N files are processed (dry run); False -> full ~44k
TEST_N = 100        # files processed when TEST_MODE is True

VAR_NAME = "et"     # .nc variable name for the single data band
VAR_ATTRS = {       # SEE "Units warning" above — confirm before publishing analysis
    "long_name": "ECHO-ET hourly estimate (source filename tag: LST; verify units)",
    "units": "1",
}
TARGET_RES_DEG = 0.01  # fallback common-grid resolution IF files are not regular lat/lon

CHUNK_FILES = 10000             # files per download/process cycle
SUB_BATCH = 100                 # files between scratch space checks
MIN_SCRATCH_GB = 40.0           # pause downloading when scratch free space falls below this
SPACE_CHECK_MIN = 30            # re-check free space at least this often (minutes)
BACKUP_SECONDS = 3600           # copy .nc (+ lists/logs) to Drive this often (1 hour)
MAX_RETRIES = 5                 # per-file download/verify attempts
CHUNK = 8 * 1024 * 1024         # streaming chunk (bytes)
MAX_DEPTH = 4                   # crawl depth below BASE_URL (year/month fits easily)
COMPRESS_LEVEL = 1              # CF .nc compression level (as specified)

# ============================================================================
# CONFIGURATION — do not edit below
# ============================================================================

BASE_URL = "https://data.globalecology.unh.edu/data/ECHO-ET/Local_solar_time/"
SHAPEFILE_DIR = "/content/ogallala_shp"
DRIVE_OUT_DIR = "/content/drive/MyDrive/ECHOET"
SCRATCH_DIR = "/content/echo_scratch"            # intermediate raw-file download area
LOCAL_NC_DIR = "/content/echo_nc"

# Local fallback: plain workstation (no /content) -> paths relative to repo.
if not os.path.isdir("/content"):
    _HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    _BASE = os.path.dirname(_HERE) if os.path.basename(_HERE).lower().startswith("local") else _HERE
    SHAPEFILE_DIR = os.path.join(_BASE, "HPA_polygon")
    DRIVE_OUT_DIR = os.path.join(_HERE, "echoet_out")
    SCRATCH_DIR = os.path.join(_HERE, "echo_scratch")
    LOCAL_NC_DIR = os.path.join(_HERE, "echo_nc")

ROI_LABEL = "Ogallala"

OUTPUT_NC = "ECHO_ET_Ogallala_hourly.nc"
URLS_TXT = "urls.txt"                 # ALL files found on server
TIF_TXT = "tif_urls.txt"              # .tif files only
LIST_CACHE = "_filelist.json"         # raw crawl cache (name + size per file)
CORRUPT_LOG = "corrupted_files.txt"   # every corrupt/failed download event (JSON lines)
FAILED_LOG = "failed_files.txt"       # files failing all retries
PROCESSED_TXT = "echo_et_processed.txt"  # resume record: basenames fully processed

MISSION_START_HINT = datetime(2018, 1, 1)
TIME_UNITS = "seconds since 2000-01-01 00:00:00 UTC"
TIME_CALENDAR = "proleptic_gregorian"
ROI_SIMPLIFY_DEG = 0.01
PLANE_CONST_MIN_PX = 200
GIBBERISH_CONST_FRAC = 0.995

# ECHO-ET_hourly_LST_2022010100.tif -> 2022-01-01 00:00
ECHO_FILENAME_RE = re.compile(r"(\d{8})(\d{2})\.tif$", re.IGNORECASE)
# Generic fallback: first YYYYMMDD(HH[MMSS]) found anywhere in the basename
GENERIC_DATE_RE = re.compile(r"(19|20)\d{2}[_-]?(0[1-9]|1[0-2])[_-]?(0[1-9]|[12]\d|3[01])")

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "Mozilla/5.0 ColabDownloader/1.0"})


# ============================================================================
# UTILITIES
# ============================================================================

def human(n):
    """Human-readable byte count."""
    n = float(n)
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return "%.1f %s" % (n, u)
        n /= 1024
    return "%.1f PB" % n


def safe_url(u):
    """Percent-encode a URL path without breaking its structure."""
    p = urlsplit(u)
    q = "?%s" % p.query if p.query else ""
    return "%s://%s%s%s" % (p.scheme, p.netloc, quote(unquote(p.path), safe="/"), q)


def disk_free_gb(path):
    return shutil.disk_usage(path).free / 1024 ** 3


def space_ok(path):
    free = disk_free_gb(path)
    log.info("Scratch free space: %.1f GB (floor %.0f GB)", free, MIN_SCRATCH_GB)
    return free >= MIN_SCRATCH_GB


def pause_until_space(path):
    """Block, re-checking every SPACE_CHECK_MIN minutes, until space is free."""
    while not space_ok(path):
        log.warning("Scratch below %.0f GB — pausing %d min before re-check.",
                    MIN_SCRATCH_GB, SPACE_CHECK_MIN)
        time.sleep(SPACE_CHECK_MIN * 60)


def parse_echo_datetime(basename):
    """Datetime for one hourly file. Primary: ECHO-ET_hourly_LST_YYYYMMDDHH.

    Fallback: first YYYYMMDD (+HH if trailing) anywhere in the name.
    Returns None when no date is recoverable (caller skips the file).
    """
    m = ECHO_FILENAME_RE.search(basename)
    if m:
        try:
            return datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H")
        except ValueError:
            pass
    m = GENERIC_DATE_RE.search(basename)
    if m:
        dofinans = m.group(0)
        digits = re.sub(r"\D", "", basename[m.end():m.end() + 7])
        hh = digits[:2] if len(digits) >= 2 and digits[:2].isdigit() and int(digits[:2]) < 24 else "12"
        try:
            return datetime.strptime(re.sub(r"\D", "", dofinans) + hh, "%Y%m%d%H")
        except ValueError:
            return None
    return None


def timestep_key(dt):
    return dt.strftime("%Y%m%dT%H%M")


# --- server crawl (h5ai JSON + HTML fallback) ---

def h5ai_json(url):
    try:
        r = SESSION.get(safe_url(url) + "?json", timeout=60)
        if r.status_code != 200:
            return []
        data = r.json()
    except Exception:
        return []
    items = []

    def walk(o):
        if isinstance(o, dict):
            if o.get("href") and o.get("type") in ("file", "folder"):
                items.append({"href": o["href"], "is_dir": o["type"] == "folder",
                              "size": o.get("size")})
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(data)
    return items


def html_listing(url):
    try:
        from bs4 import BeautifulSoup
    except ImportError:
        raise ImportError("Please install beautifulsoup4: pip install beautifulsoup4")
    try:
        r = SESSION.get(safe_url(url), timeout=60)
        r.raise_for_status()
    except Exception:
        return []
    soup = BeautifulSoup(r.text, "html.parser")
    items = []
    for a in soup.find_all("a", href=True):
        href = a["href"]
        if href.startswith(("#", "javascript", "mailto:")):
            continue
        absu = safe_url(urljoin(url, href))
        if not absu.startswith(BASE_URL):
            continue
        items.append({"href": absu, "is_dir": href.endswith("/"), "size": None})
    return items


def crawl(url, depth=0):
    """Recursively list files below *url* (skips doc/md5 sidecars)."""
    url = url if url.endswith("/") else url + "/"
    files, dirs = [], []
    for it in (h5ai_json(url) or html_listing(url)):
        name = unquote(it["href"].rstrip("/").rsplit("/", 1)[-1])
        if it["is_dir"] and not name.startswith("."):
            dirs.append(it["href"])
        elif not it["is_dir"] and not name.lower().endswith((".html", ".htm", ".md5")):
            files.append({"url": it["href"], "name": name, "size": it["size"]})
    if depth < MAX_DEPTH:
        for d in dirs:
            time.sleep(0.2)  # be gentle with the server
            files += crawl(d, depth + 1)
    return files


# --- download with resume + verification (GEDI pattern) ---

def download(url, dest, expected_size=None):
    """Stream *url* to *dest* (resumes .part via HTTP Range)."""
    url, tmp = safe_url(url), dest + ".part"
    pos = os.path.getsize(tmp) if os.path.exists(tmp) else 0
    headers = {"Range": "bytes=%d-" % pos} if pos else {}
    with SESSION.get(url, headers=headers, stream=True, timeout=(30, 180)) as r:
        if r.status_code == 416:  # .part already complete
            os.replace(tmp, dest)
            return
        r.raise_for_status()
        resume = (r.status_code == 206 and pos > 0)
        if not resume:
            pos = 0
        total = int(r.headers.get("Content-Length", 0)) + (pos if resume else 0) \
            or expected_size
        bar = tqdm(total=total, initial=pos, unit="B", unit_scale=True, unit_divisor=1024,
                   desc=os.path.basename(dest)[:45], leave=False) if tqdm else None
        try:
            with open(tmp, "ab" if resume else "wb") as fh:
                for chunk in r.iter_content(chunk_size=CHUNK):
                    if chunk:
                        fh.write(chunk)
                        if bar is not None:
                            bar.update(len(chunk))
        finally:
            if bar is not None:
                bar.close()
    if expected_size is not None and os.path.getsize(tmp) != expected_size:
        raise IOError("size mismatch: got %d, expected %s" % (os.path.getsize(tmp), expected_size))
    os.replace(tmp, dest)


def verify_tif(path, expected_size=None):
    """Transfer-integrity check: size, TIFF magic bytes, no HTML error page, open."""
    if expected_size is not None and os.path.getsize(path) != expected_size:
        return False, "size mismatch"
    try:
        with open(path, "rb") as fh:
            head = fh.read(512)
        if b"<html" in head.lower():
            return False, "got HTML error page instead of data"
        if head[:2] not in (b"II", b"MM"):
            return False, "not a valid TIFF (bad magic bytes)"
    except Exception as e:
        return False, str(e)
    try:
        with rasterio.open(path) as src:
            _ = src.bounds
            src.read(1, out_shape=(1, min(src.height, 64) or 1, min(src.width, 64) or 1))
    except Exception as e:
        return False, "rasterio read failed: %s" % e
    return True, "ok"


def check_cog(path):
    """Content check: CRS present, dims valid, a few finite pixels (nearest sample)."""
    try:
        with rasterio.open(path) as src:
            if src.count < 1 or src.width < 1 or src.height < 1:
                return False, "empty dims"
            if src.crs is None:
                return False, "missing CRS"
            h, w = src.height, src.width
            sample = src.read(1,
                              out_shape=(max(1, h // max(1, h // 500)), max(1, w // max(1, w // 500))),
                              resampling=Resampling.nearest)
            sample = sample.astype(np.float32)
            nodata = src.nodata
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


def log_corrupt(path, entry):
    """Append one JSON line immediately (survives crashes)."""
    with open(path, "a") as f:
        f.write(json.dumps(entry) + "\n")


# --- ROI + CF grid ---

def load_roi_geometry(shapefile_dir, simplify_deg=ROI_SIMPLIFY_DEG):
    shp_files = glob.glob(os.path.join(shapefile_dir, "**", "*.shp"), recursive=True)
    if not shp_files:
        raise FileNotFoundError("No .shp found in {}".format(shapefile_dir))
    log.info("Shapefile: %s", shp_files[0])
    gdf = gpd.read_file(shp_files[0])
    if gdf.crs is None or gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(epsg=4326)
    roi = (gdf.geometry.union_all() if hasattr(gdf.geometry, "union_all")
           else gdf.geometry.unary_union)
    return roi.simplify(simplify_deg, preserve_topology=True)


def build_window(sample_path, roi_4326):
    """One-time geometry: windowed read of the global grid + polygon mask.

    Returns dict(win, win_transform, height, width, lon, lat, mask_2d, mode)
    with mode "windowed" for regular geographic grids, else "reproject" onto
    a shared TARGET_RES_DEG EPSG:4326 grid (UTM-tiled style inputs).
    """
    import pyproj
    with rasterio.open(sample_path) as src:
        crs = src.crs
        is_geo = False
        try:
            is_geo = pyproj.CRS.from_user_input(crs).is_geographic
        except Exception:
            pass
        if is_geo and src.transform.a == -src.transform.e and src.transform.b == 0 \
                and src.transform.d == 0:
            # Regular geographic grid: slice the ROI bbox directly (fast path).
            win = from_bounds(*roi_4326.bounds, transform=src.transform)
            win = win.round_offsets().round_lengths()
            h, w = int(win.height), int(win.width)
            tr = src.window_transform(win)
            xs = tr.c + (np.arange(w, dtype=np.float64) + 0.5) * tr.a
            ys = tr.f + (np.arange(h, dtype=np.float64) + 0.5) * tr.e
            mask = rasterize(
                [(g, 1) for g in getattr(roi_4326, "geoms", [roi_4326])],
                out_shape=(h, w), transform=tr, fill=0, dtype="uint8",
                all_touched=True).astype(np.int8)
            log.info("Windowed mode: %d x %d subset, %d px in ROI (%.1f %%)",
                     h, w, int(mask.sum()), 100.0 * mask.sum() / mask.size)
            # Sample stats so the operator can judge the units question.
            sub = src.read(1, window=win).astype(np.float32)
            if src.nodata is not None and not (isinstance(src.nodata, float)
                                               and np.isnan(src.nodata)):
                sub[sub == src.nodata] = np.nan
            fin = sub[np.isfinite(sub)]
            if fin.size:
                log.info("Sample-value probe (%s): n=%d min=%.4g max=%.4g mean=%.4g "
                         "(K~250-330 suggests LST; ~0-2 suggests mm/hr ET)",
                         os.path.basename(sample_path), fin.size,
                         float(fin.min()), float(fin.max()), float(fin.mean()))
            return {"mode": "windowed", "win": win, "lon": xs, "lat": ys,
                    "mask_2d": mask, "height": h, "width": w, "transform": tr}
    # Fallback: reprojected common grid (non-geographic / rotated inputs).
    minx, miny, maxx, maxy = roi_4326.bounds
    res = TARGET_RES_DEG
    minx = float(np.floor(minx / res) * res)
    miny = float(np.floor(miny / res) * res)
    maxx = float(np.ceil(maxx / res) * res)
    maxy = float(np.ceil(maxy / res) * res)
    nx = int(round((maxx - minx) / res)) + 1
    ny = int(round((maxy - miny) / res)) + 1
    lon = minx + np.arange(nx, dtype=np.float64) * res
    lat = maxy - np.arange(ny, dtype=np.float64) * res
    tr = from_origin(minx - res / 2.0, maxy + res / 2.0, res, res)
    mask = rasterize(
        [(g, 1) for g in getattr(roi_4326, "geoms", [roi_4326])],
        out_shape=(ny, nx), transform=tr, fill=0, dtype="uint8",
        all_touched=True).astype(np.int8)
    log.info("Reproject mode: %d x %d common grid", ny, nx)
    return {"mode": "reproject", "lon": lon, "lat": lat,
            "mask_2d": mask, "height": ny, "width": nx, "transform": tr}


def extract_filtered(tif_path, grid):
    """Read one file, keep ROI pixels, return 2D float32 (NaN outside ROI) or None."""
    arr = None
    if grid["mode"] == "windowed":
        with rasterio.open(tif_path) as src:
            a = src.read(1, window=grid["win"]).astype(np.float32)
            if src.nodata is not None and not (isinstance(src.nodata, float)
                                               and np.isnan(src.nodata)):
                a[a == src.nodata] = np.nan
            a[~np.isfinite(a)] = np.nan
            if a.shape != (grid["height"], grid["width"]):
                dst = np.full((grid["height"], grid["width"]), np.nan, dtype=np.float32)
                h = min(a.shape[0], grid["height"])
                w = min(a.shape[1], grid["width"])
                dst[:h, :w] = a[:h, :w]
                a = dst
            arr = a
    else:
        with rasterio.open(tif_path) as src:
            a = src.read(1).astype(np.float32)
            if src.nodata is not None and not (isinstance(src.nodata, float)
                                               and np.isnan(src.nodata)):
                a[a == src.nodata] = np.nan
            a[~np.isfinite(a)] = np.nan
            dst = np.full((grid["height"], grid["width"]), np.nan, dtype=np.float32)
            reproject(a, dst, src_transform=src.transform, src_crs=src.crs,
                      dst_transform=grid["transform"], dst_crs="EPSG:4326",
                      resampling=Resampling.average, dst_nodata=np.nan)
            arr = dst
    arr[grid["mask_2d"] == 0] = np.nan
    fin = arr[np.isfinite(arr)]
    if fin.size == 0:
        return None
    if fin.size >= PLANE_CONST_MIN_PX:
        vals, counts = np.unique(fin, return_counts=True)
        if counts.max() / fin.size > GIBBERISH_CONST_FRAC:
            log.warning("Dropping constant-plane file %s (v=%s)",
                        os.path.basename(tif_path), vals[np.argmax(counts)])
            return None
    return arr


def create_output_nc(path, lon, lat):
    """Fresh CF-1.8 file, compression level 1 throughout."""
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
    t.long_name = "ECHO-ET file time (UTC, from filename)"
    y = ds.createVariable("lat", "f4", ("y",), zlib=True, complevel=COMPRESS_LEVEL)
    y.units = "degrees_north"
    y.standard_name = "latitude"
    y.axis = "Y"
    y[:] = lat.astype(np.float32)
    x = ds.createVariable("lon", "f4", ("x",), zlib=True, complevel=COMPRESS_LEVEL)
    x.units = "degrees_east"
    x.standard_name = "longitude"
    x.axis = "X"
    x[:] = lon.astype(np.float32)
    v = ds.createVariable(VAR_NAME, "f4", ("time", "y", "x"),
                          fill_value=np.float32(np.nan),
                          zlib=True, complevel=COMPRESS_LEVEL,
                          chunksizes=(1, min(512, ny), min(512, nx)))
    v.long_name = VAR_ATTRS["long_name"]
    v.units = VAR_ATTRS["units"]
    v.coordinates = "lat lon"
    v.grid_mapping = "crs"
    v.cell_methods = "time: point (hourly file)"
    crs = ds.createVariable("crs", "i4")
    crs.grid_mapping_name = "latitude_longitude"
    crs.semi_major_axis = 6378137.0
    crs.inverse_flattening = 298.257223563
    ds.Conventions = "CF-1.8"
    ds.title = "ECHO-ET hourly — {} ROI subset".format(ROI_LABEL)
    ds.source = "UNH ECHO-ET Local_solar_time (%s)" % BASE_URL
    ds.institution = "University of New Hampshire / ECOSTRESS"
    ds.history = "Created {}".format(datetime.now(timezone.utc).isoformat())
    ds.comment = ("Hourly files windowed to the ROI bbox, polygon-masked, NaN fill. "
                  "Stored values are raw file values; VERIFY physical units before "
                  "analysis (server filenames carry the tag LST).")
    ds.roi = ROI_LABEL
    ds.close()
    log.info("Created %s (%d x %d grid, complevel=%d)",
             os.path.basename(path), ny, nx, COMPRESS_LEVEL)


def get_resume_datetimes(path):
    keys, max_dt = set(), None
    if not os.path.exists(path):
        return keys, max_dt
    try:
        with nc4.Dataset(path, "r") as ds:
            if len(ds.variables["time"]) == 0:
                return keys, max_dt
            units = ds.variables["time"].units
            cal = getattr(ds.variables["time"], "calendar", TIME_CALENDAR)
            for v in ds.variables["time"][:]:
                dt = nc4.num2date(v, units, cal)
                dt = datetime(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
                keys.add(timestep_key(dt))
                if max_dt is None or dt > max_dt:
                    max_dt = dt
        log.info("Resume %s: %d timesteps (max %s)", os.path.basename(path), len(keys), max_dt)
    except Exception as e:
        log.warning("Cannot read resume point from %s: %s", path, e)
    return keys, max_dt


def append_timestep(nc_path, dt, arr):
    """Append one hourly timestep; all-NaN arrays are refused. Returns bool."""
    a = np.asarray(arr, dtype=np.float32)
    if a.size == 0 or np.isfinite(a).sum() == 0:
        return False
    ds = nc4.Dataset(nc_path, "a")
    try:
        t_idx = len(ds.variables["time"])
        ds.variables["time"][t_idx] = nc4.date2num(dt, TIME_UNITS, TIME_CALENDAR)
        ds.variables[VAR_NAME][t_idx, :, :] = a
        ds.sync()
    finally:
        ds.close()
    return True


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


# --- resume-log helpers ---

def load_processed(path):
    done = set()
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line:
                    done.add(line)
        log.info("Resume log: %d files already processed", len(done))
    return done


def mark_processed(path, basenames):
    if not basenames:
        return
    if isinstance(basenames, str):
        basenames = [basenames]
    with open(path, "a") as f:
        for b in basenames:
            f.write(b + "\n")
        f.flush()
        os.fsync(f.fileno())


# ============================================================================
# PIPELINE
# ============================================================================

def run_pipeline():
    log.info("=" * 64)
    log.info("  ECHO-ET URL-list pipeline — %s (TEST_MODE=%s)", ROI_LABEL, TEST_MODE)
    log.info("=" * 64)

    for d in (SCRATCH_DIR, LOCAL_NC_DIR, DRIVE_OUT_DIR):
        os.makedirs(d, exist_ok=True)

    local_nc = os.path.join(LOCAL_NC_DIR, OUTPUT_NC)
    local_urls = os.path.join(LOCAL_NC_DIR, URLS_TXT)
    local_tifs = os.path.join(LOCAL_NC_DIR, TIF_TXT)
    local_cache = os.path.join(LOCAL_NC_DIR, LIST_CACHE)
    local_corrupt = os.path.join(LOCAL_NC_DIR, CORRUPT_LOG)
    local_failed = os.path.join(LOCAL_NC_DIR, FAILED_LOG)
    local_processed = os.path.join(LOCAL_NC_DIR, PROCESSED_TXT)

    # Restore resume artifacts from Drive when missing locally.
    for name, local in ((OUTPUT_NC, local_nc), (URLS_TXT, local_urls),
                        (TIF_TXT, local_tifs), (LIST_CACHE, local_cache),
                        (CORRUPT_LOG, local_corrupt), (FAILED_LOG, local_failed),
                        (PROCESSED_TXT, local_processed)):
        remote = os.path.join(DRIVE_OUT_DIR, name)
        if not os.path.exists(local) and os.path.exists(remote):
            log.info("Restoring %s from Drive for resume...", name)
            shutil.copy2(remote, local)

    roi = load_roi_geometry(SHAPEFILE_DIR)

    # ── STAGE 1: urls.txt + tif_urls.txt ──────────────────────────────────
    if os.path.exists(local_cache) and os.path.getsize(local_cache) > 0:
        all_files = json.load(open(local_cache))
        log.info("Loaded cached listing: %d files", len(all_files))
    else:
        log.info("Crawling server (one-time, ~minutes)...")
        all_files = crawl(BASE_URL)
        if tqdm is not None:
            bar = tqdm(all_files, desc="HEAD size lookup", unit="file")
        else:
            bar = all_files
        for x in bar:
            if x["size"] is None:
                try:
                    h = SESSION.head(safe_url(x["url"]), allow_redirects=True, timeout=60)
                    x["size"] = int(h.headers["Content-Length"]) if "Content-Length" in h.headers else None
                except Exception:
                    x["size"] = None
        if tqdm is not None:
            bar.close()
        json.dump(all_files, open(local_cache, "w"))
    with open(local_urls, "w") as f:
        f.write("\n".join(x["url"] for x in all_files))
    log.info("%d total URLs -> %s", len(all_files), local_urls)

    tif_files = [x for x in all_files if x["name"].lower().endswith(".tif")]
    tif_files.sort(key=lambda x: x["name"])
    with open(local_tifs, "w") as f:
        f.write("\n".join(x["url"] for x in tif_files))
    total_bytes = sum(x["size"] or 0 for x in tif_files)
    log.info("%d .tif files (%s) -> %s", len(tif_files), human(total_bytes), local_tifs)
    for name, local in ((URLS_TXT, local_urls), (TIF_TXT, local_tifs), (LIST_CACHE, local_cache)):
        try:
            robust_drive_copy(local, os.path.join(DRIVE_OUT_DIR, name))
        except Exception as e:
            log.error("%s Drive backup failed: %s", name, e)

    if TEST_MODE:
        tif_files = tif_files[:TEST_N]
        log.info("TEST_MODE: processing first %d files only", len(tif_files))
    log.info("Sample: %s", [x["name"] for x in tif_files[:5]])
    log.info("Scratch free: %s | total to fetch: %s",
             human(shutil.disk_usage(SCRATCH_DIR).free), human(total_bytes))

    processed = load_processed(local_processed)
    todo = [x for x in tif_files if x["name"] not in processed]
    log.info("Files: %d total, %d already processed, %d to do",
             len(tif_files), len(tif_files) - len(todo), len(todo))

    # One-time ROI geometry from the first available file (download it if needed).
    grid = None
    if todo:
        probe = todo[0]
        probe_local = os.path.join(SCRATCH_DIR, probe["name"])
        if not os.path.exists(probe_local):
            download(probe["url"], probe_local, probe["size"])
        grid = build_window(probe_local, roi)
    if grid is None:
        log.info("Nothing to do — all files already processed.")
        return {}
    if not os.path.exists(local_nc):
        create_output_nc(local_nc, grid["lon"], grid["lat"])
    done_keys, _ = get_resume_datetimes(local_nc)

    total_steps = 0
    last_backup = time.time()
    last_space_check = 0.0

    def backup_if_due(force=False):
        nonlocal last_backup
        if not force and (time.time() - last_backup) < BACKUP_SECONDS:
            return
        log.info("1-hour backup: copying outputs to Drive...")
        try:
            with nc4.Dataset(local_nc, "a") as _ds:
                _ds.sync()
            for name, local in ((OUTPUT_NC, local_nc), (PROCESSED_TXT, local_processed),
                                (CORRUPT_LOG, local_corrupt), (FAILED_LOG, local_failed)):
                if os.path.exists(local):
                    robust_drive_copy(local, os.path.join(DRIVE_OUT_DIR, name))
        except Exception as e:
            log.error("Hourly backup failed: %s", e)
        last_backup = time.time()

    failed = []
    if os.path.exists(local_failed):
        failed = [l.strip() for l in open(local_failed) if l.strip()]

    # ── STAGE 2: 10,000-file download / process / delete cycles ───────────
    pos = 0
    while pos < len(todo):
        chunk = todo[pos:pos + CHUNK_FILES]
        pos += CHUNK_FILES
        log.info("=" * 50)
        log.info("CHUNK: %d files", len(chunk))
        log.info("=" * 50)

        # — download phase (space-guarded, resumable) —
        staged = []
        pending = list(chunk)
        if not space_ok(SCRATCH_DIR):
            pause_until_space(SCRATCH_DIR)
        while pending:
            if (time.time() - last_space_check) >= SPACE_CHECK_MIN * 60:
                last_space_check = time.time()
                if not space_ok(SCRATCH_DIR):
                    log.warning("Space floor hit — pausing download, processing staged files.")
                    break
            sub = pending[:SUB_BATCH]
            pending = pending[SUB_BATCH:]
            for x in sub:
                dest = os.path.join(SCRATCH_DIR, x["name"])
                if os.path.exists(dest):
                    staged.append(x)
                    continue
                ok = False
                for attempt in range(1, MAX_RETRIES + 1):
                    try:
                        if os.path.exists(dest):
                            os.remove(dest)
                        download(x["url"], dest, x["size"])
                        good, why = verify_tif(dest, x["size"])
                        if not good:
                            raise IOError("corrupt: %s" % why)
                        good, why = check_cog(dest)
                        if not good:
                            raise IOError("content check: %s" % why)
                        ok = True
                        break
                    except Exception as e:
                        msg = "attempt %d/%d ERROR: %s" % (attempt, MAX_RETRIES, e)
                        log.warning("  [%s] %s", x["name"], msg)
                        entry = {"url": x["url"], "attempt": attempt, "reason": str(e)}
                        log_corrupt(local_corrupt, entry)
                        time.sleep(min(2 ** attempt, 30))
                if ok:
                    staged.append(x)
                else:
                    failed.append(x["url"])
                    with open(local_failed, "w") as f:
                        f.write("\n".join(failed))
            gc.collect()
        if pending:
            log.warning("%d files not downloaded this cycle (space) — resume next run.", len(pending))

        # — filter-to-ROI + append phase (one timestep per file) —
        staged.sort(key=lambda x: x["name"])
        for x in staged:
            dest = os.path.join(SCRATCH_DIR, x["name"])
            if not os.path.exists(dest):
                continue
            dt = parse_echo_datetime(x["name"])
            if dt is None:
                log.warning("No date in filename, skipping: %s", x["name"])
                continue
            key = timestep_key(dt)
            if key in done_keys or x["name"] in processed:
                try:
                    os.remove(dest)
                except OSError:
                    pass
                if x["name"] not in processed:
                    mark_processed(local_processed, x["name"])
                    processed.add(x["name"])
                done_keys.add(key)
                continue
            try:
                arr = extract_filtered(dest, grid)
            except Exception as e:
                log.error("Extract failed for %s: %s", x["name"], e)
                arr = None
            try:
                if os.path.exists(dest):
                    os.remove(dest)
            except OSError:
                pass
            if arr is None:
                log.info("No valid ROI pixels for %s — marking done, nothing appended.", key)
                mark_processed(local_processed, x["name"])
                processed.add(x["name"])
                done_keys.add(key)
                continue
            try:
                if append_timestep(local_nc, dt, arr):
                    total_steps += 1
                    if total_steps % 25 == 0:
                        log.info("Appended %d timesteps (latest %s).", total_steps, key)
                done_keys.add(key)
                mark_processed(local_processed, x["name"])
                processed.add(x["name"])
            except Exception as e:
                log.error("Append failed for %s: %s", x["name"], e)
            finally:
                del arr
                gc.collect()
            backup_if_due()

        # — end-of-chunk scratch sweep —
        left = [f for f in glob.glob(os.path.join(SCRATCH_DIR, "*.tif"))
                if os.path.basename(f) in processed]
        for f in left:
            try:
                os.remove(f)
            except OSError:
                pass
        n_left = len(glob.glob(os.path.join(SCRATCH_DIR, "*.tif*")))
        log.info("Chunk done. Steps appended: %d. Scratch leftovers: %d.", total_steps, n_left)
        backup_if_due()
        gc.collect()

    backup_if_due(force=True)
    n_in_nc = len(nc4.Dataset(local_nc, "r").variables["time"])
    log.info("=" * 64)
    log.info("  PIPELINE COMPLETE — %d new steps, %d total in .nc | TEST_MODE=%s",
             total_steps, n_in_nc, TEST_MODE)
    log.info("  NC  : %s", os.path.join(DRIVE_OUT_DIR, OUTPUT_NC))
    log.info("  URLS: %s , %s", os.path.join(DRIVE_OUT_DIR, URLS_TXT),
             os.path.join(DRIVE_OUT_DIR, TIF_TXT))
    log.info("  LOG : %s", os.path.join(DRIVE_OUT_DIR, PROCESSED_TXT))
    log.info("=" * 64)
    return {"steps": total_steps, "total": n_in_nc}


if __name__ == "__main__":
    # Running as a script (Colab "Run all" or python file.py) launches the pipeline.
    # In a notebook, calling run_pipeline() manually does the same thing.
    run_pipeline()
else:
    # Imported as a module: do NOT auto-run; user calls run_pipeline() in next cell.
    pass
