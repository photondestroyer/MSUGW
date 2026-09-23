#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
smap_l4_extractor.py — SMAP L4 (SPL4SMGP) ROI Extraction Pipeline for Google Colab
====================================================================================

Downloads SPL4SMGP HDF5 granules from NASA Earthdata, extracts 5 target variables
within the Ogallala aquifer ROI on the native EASE-2 9km grid, and streams them
into a single CF-compliant NetCDF file on Google Drive.

Architecture:
  1. Builds a one-time spatial mask from shapefile + EASE-2 grid coordinates.
  2. Downloads granules in batches of 175 dates (1400 files) via earthaccess (4 threads).
  3. Extracts target variables at the ROI subset, applies mask, appends to NetCDF.
  4. Flushes to disk every 50 dates. Copies to Drive after each batch.
  5. Full resume capability: reads last timestamp from existing output file.
  6. Handles corrupt granules with automatic re-download.

Native grid: EASE-2 9km cylindrical equal-area (1624 × 3856 global).
Cadence: 3-hourly (8 granules per day).

Usage on Colab:
  Cell 1: !pip install -q earthaccess h5py netCDF4 geopandas shapely tqdm
  Cell 2: from google.colab import drive; drive.mount('/content/drive')
  Cell 3: from smap_l4_extractor import run_pipeline; run_pipeline()
"""

# ==========================================
# Install (run as first Colab cell)
# ==========================================
# !pip install -q earthaccess h5py netCDF4 geopandas shapely tqdm

import os
import re
import gc
import glob
import time
import shutil
import logging
import hashlib
from datetime import datetime, timezone, timedelta
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import h5py
import netCDF4 as nc4
import geopandas as gpd
from shapely.geometry import Point
from tqdm.auto import tqdm

try:
    import earthaccess
except ImportError:
    raise ImportError("Please install earthaccess: !pip install -q earthaccess")

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('smap_l4')

# ╔════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION                                                     ║
# ╚════════════════════════════════════════════════════════════════════╝

SHORT_NAME      = "SPL4SMGP"
SHAPEFILE_DIR   = "/content/ogallala_shp"
DRIVE_OUTPUT    = "/content/drive/MyDrive/MSUGWB/SPL4SMGP"
LOCAL_RAW_DIR   = "/content/smap_raw"
LOCAL_NC_DIR    = "/content/smap_nc"

OUTPUT_FILENAME = "SPL4SMGP_Ogallala.nc"
ROI_LABEL       = "Ogallala"

# Variables to extract from Geophysical_Data group
TARGET_VARS = [
    "baseflow_flux",
    "depth_to_water_table_from_surface",
    "land_evapotranspiration_flux",
    "sm_rootzone",
    "sm_rootzone_pctl",
]

FILL_VALUE      = np.float32(-9999.0)
COMPRESS_LEVEL  = 3
BATCH_DATES     = 175       # Download 175 dates at a time
FLUSH_EVERY     = 50        # Flush to disk every 50 dates (400 granules)
DOWNLOAD_THREADS = 4
MAX_RETRY       = 2         # Max re-download attempts for corrupt files
ROI_SIMPLIFY_DEG = 0.1      # Coarsen shapefile to 0.1 degree

# Regex for SMAP L4 filenames: SMAP_L4_SM_gph_YYYYMMDDT\d{4}0_V...
SMAP_FILENAME_RE = re.compile(
    r'SMAP_L4_SM_gph_(\d{8})T(\d{4})0_V.*\.h5$'
)

# ╔════════════════════════════════════════════════════════════════════╗
# ║  FILENAME PARSING                                                  ║
# ╚════════════════════════════════════════════════════════════════════╝

def parse_smap_datetime(filename):
    """
    Extract datetime from SMAP L4 filename using regex.

    Filename pattern: SMAP_L4_SM_gph_YYYYMMDDT\\d{4}0_Vv8010_001.h5
    Example: SMAP_L4_SM_gph_20170101T013000_Vv8010_001.h5
             -> datetime(2017, 1, 1, 1, 30)

    Returns datetime or None if pattern doesn't match.
    """
    basename = os.path.basename(filename)
    m = SMAP_FILENAME_RE.search(basename)
    if not m:
        return None
    date_str = m.group(1)  # YYYYMMDD
    time_str = m.group(2)  # HHMM
    try:
        dt = datetime.strptime(date_str + time_str, '%Y%m%d%H%M')
        return dt
    except ValueError:
        return None


def sort_files_chronologically(file_list):
    """Sort HDF5 files chronologically by parsed datetime from filename."""
    pairs = []
    for f in file_list:
        dt = parse_smap_datetime(f)
        if dt is not None:
            pairs.append((dt, f))
        else:
            log.warning('Cannot parse datetime from filename, skipping: %s',
                        os.path.basename(f))
    pairs.sort(key=lambda x: x[0])
    return pairs


# ╔════════════════════════════════════════════════════════════════════╗
# ║  ROI MASK BUILDING                                                 ║
# ╚════════════════════════════════════════════════════════════════════╝

def load_roi_geometry(shapefile_dir, simplify_deg=0.1):
    """Load shapefile, reproject to EPSG:4326, simplify to given tolerance."""
    shp_files = glob.glob(os.path.join(shapefile_dir, '**', '*.shp'), recursive=True)
    if not shp_files:
        raise FileNotFoundError(f'No .shp found in {shapefile_dir}')
    gdf = gpd.read_file(shp_files[0])
    if gdf.crs is None or gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(epsg=4326)
    roi = (gdf.geometry.union_all() if hasattr(gdf.geometry, 'union_all')
           else gdf.geometry.unary_union)
    roi_simple = roi.simplify(simplify_deg, preserve_topology=True)
    log.info('ROI loaded and simplified to %.2f° tolerance. Bounds: %s',
             simplify_deg, roi_simple.bounds)
    return roi_simple


def build_roi_mask(sample_h5_path, roi_geometry):
    """
    Build a boolean mask and row/col slices for the ROI on the EASE-2 grid.

    Steps:
      1. Read cell_lat / cell_lon from the HDF5 file root.
      2. Compute bounding box mask (fast, vectorized).
      3. Within the bbox, test each point against the simplified ROI polygon.
      4. Return the minimal row/col slice and the 2D boolean mask within that slice.

    Returns:
        (row_slice, col_slice, mask_2d, lat_subset, lon_subset)
    """
    log.info('Building ROI mask from EASE-2 grid...')
    with h5py.File(sample_h5_path, 'r') as f:
        cell_lat = f['cell_lat'][:]  # (1624, 3856)
        cell_lon = f['cell_lon'][:]  # (1624, 3856)

    minx, miny, maxx, maxy = roi_geometry.bounds

    # Fast bounding box filter
    bbox_mask = (
        (cell_lat >= miny) & (cell_lat <= maxy) &
        (cell_lon >= minx) & (cell_lon <= maxx)
    )

    # Find the minimal row/col range containing all bbox pixels
    rows_any = np.any(bbox_mask, axis=1)
    cols_any = np.any(bbox_mask, axis=0)
    row_indices = np.where(rows_any)[0]
    col_indices = np.where(cols_any)[0]

    if len(row_indices) == 0 or len(col_indices) == 0:
        raise ValueError('No EASE-2 grid cells found within ROI bounding box!')

    r0, r1 = int(row_indices[0]), int(row_indices[-1]) + 1
    c0, c1 = int(col_indices[0]), int(col_indices[-1]) + 1

    row_slice = slice(r0, r1)
    col_slice = slice(c0, c1)

    lat_sub = cell_lat[row_slice, col_slice]
    lon_sub = cell_lon[row_slice, col_slice]
    bbox_sub = bbox_mask[row_slice, col_slice]

    # Fine-grained polygon containment test within bbox subset
    # Use prepared geometry for speed
    from shapely.prepared import prep
    roi_prep = prep(roi_geometry)

    ny, nx = lat_sub.shape
    mask_2d = np.zeros((ny, nx), dtype=bool)

    log.info('Testing %d candidate pixels against ROI polygon...', int(bbox_sub.sum()))
    for iy in range(ny):
        for ix in range(nx):
            if bbox_sub[iy, ix]:
                pt = Point(float(lon_sub[iy, ix]), float(lat_sub[iy, ix]))
                if roi_prep.contains(pt):
                    mask_2d[iy, ix] = True

    n_in = int(mask_2d.sum())
    log.info('ROI mask: %d × %d subset, %d pixels inside polygon '
             '(rows %d:%d, cols %d:%d)',
             ny, nx, n_in, r0, r1, c0, c1)

    if n_in == 0:
        raise ValueError('No grid cells fall within the ROI polygon!')

    return row_slice, col_slice, mask_2d, lat_sub, lon_sub


# ╔════════════════════════════════════════════════════════════════════╗
# ║  NETCDF CREATION & RESUME                                         ║
# ╚════════════════════════════════════════════════════════════════════╝

def create_output_netcdf(filepath, lat_2d, lon_2d, mask_2d,
                         row_slice, col_slice):
    """
    Create the pre-allocated CF-compliant NetCDF4 output file.

    Dimensions: time (unlimited), y, x
    Coordinates: lat(y,x), lon(y,x), time
    Variables: 5 target vars (time, y, x), float32, fill=-9999.0, zlib complevel=3
    """
    ny, nx = lat_2d.shape

    ds = nc4.Dataset(filepath, 'w', format='NETCDF4')

    # Dimensions
    ds.createDimension('time', None)  # unlimited
    ds.createDimension('y', ny)
    ds.createDimension('x', nx)

    # Time coordinate
    t_var = ds.createVariable('time', 'f8', ('time',))
    t_var.units = 'hours since 2000-01-01 00:00:00 UTC'
    t_var.calendar = 'proleptic_gregorian'
    t_var.axis = 'T'
    t_var.standard_name = 'time'
    t_var.long_name = 'time of observation'

    # 2D lat/lon coordinates (native EASE-2 grid, NOT interpolated)
    lat_var = ds.createVariable('lat', 'f4', ('y', 'x'),
                                zlib=True, complevel=1)
    lat_var.units = 'degrees_north'
    lat_var.standard_name = 'latitude'
    lat_var.long_name = 'latitude of grid cell center'
    lat_var[:] = lat_2d

    lon_var = ds.createVariable('lon', 'f4', ('y', 'x'),
                                zlib=True, complevel=1)
    lon_var.units = 'degrees_east'
    lon_var.standard_name = 'longitude'
    lon_var.long_name = 'longitude of grid cell center'
    lon_var[:] = lon_2d

    # ROI mask (static, useful for downstream users)
    m_var = ds.createVariable('roi_mask', 'i1', ('y', 'x'),
                              zlib=True, complevel=1)
    m_var.long_name = 'Region of Interest mask (1=inside ROI, 0=outside)'
    m_var.flag_values = np.array([0, 1], dtype='i1')
    m_var.flag_meanings = 'outside_roi inside_roi'
    m_var[:] = mask_2d.astype(np.int8)

    # Data variables
    var_metadata = {
        'baseflow_flux': {
            'long_name': 'Baseflow flux',
            'units': 'kg m-2 s-1',
        },
        'depth_to_water_table_from_surface': {
            'long_name': 'Depth to water table from surface',
            'units': 'm',
        },
        'land_evapotranspiration_flux': {
            'long_name': 'Land evapotranspiration flux',
            'units': 'kg m-2 s-1',
        },
        'sm_rootzone': {
            'long_name': 'Root zone soil moisture',
            'units': 'm3 m-3',
        },
        'sm_rootzone_pctl': {
            'long_name': 'Root zone soil moisture percentile',
            'units': '%',
        },
    }

    for vname in TARGET_VARS:
        v = ds.createVariable(
            vname, 'f4', ('time', 'y', 'x'),
            fill_value=FILL_VALUE,
            zlib=True,
            complevel=COMPRESS_LEVEL,
            chunksizes=(1, min(256, ny), min(256, nx)),
        )
        meta = var_metadata.get(vname, {})
        v.long_name = meta.get('long_name', vname)
        v.units = meta.get('units', 'unknown')
        v.grid_mapping = 'crs'
        v.coordinates = 'lat lon'

    # CRS / grid mapping for EASE-2
    crs_var = ds.createVariable('crs', 'i4')
    crs_var.grid_mapping_name = 'lambert_cylindrical_equal_area'
    crs_var.longitude_of_central_meridian = 0.0
    crs_var.standard_parallel = 30.0
    crs_var.false_easting = 0.0
    crs_var.false_northing = 0.0
    crs_var.long_name = 'EASE-2 Global Cylindrical Equal-Area Projection'
    crs_var.crs_wkt = (
        'PROJCRS["EASE-Grid 2.0 Global",'
        'BASEGEOGCRS["WGS 84",DATUM["World Geodetic System 1984",'
        'ELLIPSOID["WGS 84",6378137,298.257223563]]],'
        'CONVERSION["EASE-Grid 2.0 Global",'
        'METHOD["Lambert Cylindrical Equal Area"],'
        'PARAMETER["Latitude of 1st standard parallel",30],'
        'PARAMETER["Longitude of natural origin",0],'
        'PARAMETER["False easting",0],'
        'PARAMETER["False northing",0]],'
        'CS[Cartesian,2],AXIS["easting",east],AXIS["northing",north],'
        'UNIT["metre",1]]'
    )

    # Global attributes
    ds.Conventions = 'CF-1.8'
    ds.title = f'SMAP L4 SM (SPL4SMGP) — {ROI_LABEL} ROI Subset'
    ds.source = 'NASA SMAP Level 4 Surface and Root Zone Soil Moisture (SPL4SMGP)'
    ds.institution = 'NASA GMAO'
    ds.history = f'Created {datetime.now(timezone.utc).isoformat()}'
    ds.references = 'https://nsidc.org/data/spl4smgp'
    ds.comment = (
        f'Spatial subset of SPL4SMGP on native EASE-2 9km grid. '
        f'ROI: {ROI_LABEL}. '
        f'Row slice: {row_slice.start}:{row_slice.stop}, '
        f'Col slice: {col_slice.start}:{col_slice.stop}. '
        f'No interpolation or regridding applied.'
    )

    ds.close()
    log.info('Created output NetCDF: %s (%d×%d grid, %d variables)',
             filepath, ny, nx, len(TARGET_VARS))


def get_resume_point(filepath):
    """
    Read the last timestamp from an existing output file for resume.

    Returns datetime or None if file doesn't exist or has no time steps.
    """
    if not os.path.exists(filepath):
        return None, 0
    try:
        ds = nc4.Dataset(filepath, 'r')
        t_var = ds.variables['time']
        n_times = len(t_var)
        if n_times == 0:
            ds.close()
            return None, 0
        last_val = float(t_var[-1])
        units = t_var.units
        calendar = getattr(t_var, 'calendar', 'proleptic_gregorian')
        last_dt = nc4.num2date(last_val, units, calendar)
        # Convert cftime to standard datetime
        last_dt = datetime(last_dt.year, last_dt.month, last_dt.day,
                           last_dt.hour, last_dt.minute, last_dt.second)
        ds.close()
        log.info('Resume point: last timestamp = %s (index %d)',
                 last_dt.isoformat(), n_times - 1)
        return last_dt, n_times
    except Exception as e:
        log.warning('Could not read resume point: %s', e)
        return None, 0


# ╔════════════════════════════════════════════════════════════════════╗
# ║  DOWNLOAD                                                         ║
# ╚════════════════════════════════════════════════════════════════════╝

def search_granules(start_date, end_date):
    """Search earthaccess for SPL4SMGP granules in the date range."""
    results = earthaccess.search_data(
        short_name=SHORT_NAME,
        temporal=(start_date.strftime('%Y-%m-%d'),
                  end_date.strftime('%Y-%m-%d')),
    )
    log.info('Search returned %d granules for %s to %s',
             len(results), start_date.date(), end_date.date())
    return results


def download_granules(results, local_dir, threads=4):
    """Download granules to local_dir using earthaccess with N threads."""
    os.makedirs(local_dir, exist_ok=True)
    files = earthaccess.download(results, local_path=local_dir,
                                 threads=threads)
    return [str(f) for f in files]


def validate_h5(filepath):
    """Check if an HDF5 file is readable and has the expected structure."""
    try:
        with h5py.File(filepath, 'r') as f:
            _ = f['Geophysical_Data'][TARGET_VARS[0]].shape
        return True
    except Exception:
        return False


# ╔════════════════════════════════════════════════════════════════════╗
# ║  EXTRACTION                                                        ║
# ╚════════════════════════════════════════════════════════════════════╝

def extract_granule(h5_path, row_slice, col_slice, mask_2d):
    """
    Read target variables from one HDF5 granule at the ROI subset.

    Returns dict of {var_name: 2d_array} or None if corrupt.
    Out-of-ROI pixels are set to FILL_VALUE.
    """
    try:
        with h5py.File(h5_path, 'r') as f:
            gd = f['Geophysical_Data']
            data = {}
            for vname in TARGET_VARS:
                arr = gd[vname][row_slice, col_slice].astype(np.float32)
                # Apply ROI mask: set pixels outside ROI to fill value
                arr[~mask_2d] = FILL_VALUE
                data[vname] = arr
        return data
    except Exception as e:
        log.error('Failed to read %s: %s', os.path.basename(h5_path), e)
        return None


def append_timestep(nc_path, dt, data_dict):
    """
    Append one time step to the output NetCDF file.

    dt: datetime of this granule
    data_dict: {var_name: 2D array (y, x)}
    """
    ds = nc4.Dataset(nc_path, 'a')
    try:
        t_var = ds.variables['time']
        t_idx = len(t_var)

        # Convert datetime to numeric
        t_num = nc4.date2num(dt, t_var.units, t_var.calendar)
        t_var[t_idx] = t_num

        for vname, arr in data_dict.items():
            ds.variables[vname][t_idx, :, :] = arr
    finally:
        ds.close()


# ╔════════════════════════════════════════════════════════════════════╗
# ║  DRIVE COPY                                                        ║
# ╚════════════════════════════════════════════════════════════════════╝

def get_md5(filepath, chunk_bytes=10 * 1024 * 1024):
    """Calculate MD5 checksum."""
    h = hashlib.md5()
    with open(filepath, 'rb') as f:
        while True:
            chunk = f.read(chunk_bytes)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def robust_drive_copy(src, dst, max_retries=3):
    """Copy file to Drive with MD5 verification."""
    src_md5 = get_md5(src)
    for attempt in range(1, max_retries + 1):
        try:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with open(src, 'rb') as fs, open(dst, 'wb') as fd:
                shutil.copyfileobj(fs, fd, length=10 * 1024 * 1024)
                fd.flush()
                os.fsync(fd.fileno())
            try:
                os.sync()
            except AttributeError:
                pass
            time.sleep(2)
            if os.path.getsize(src) != os.path.getsize(dst):
                raise IOError('Size mismatch')
            if src_md5 != get_md5(dst):
                raise IOError('MD5 mismatch')
            log.info('Drive copy verified (%.1f MiB)',
                     os.path.getsize(dst) / 1024**2)
            return True
        except Exception as e:
            log.error('Drive copy attempt %d failed: %s', attempt, e)
            if os.path.exists(dst):
                try:
                    os.remove(dst)
                except OSError:
                    pass
            if attempt < max_retries:
                time.sleep(10 * attempt)
            else:
                raise RuntimeError(f'Drive copy failed: {e}')


# ╔════════════════════════════════════════════════════════════════════╗
# ║  MAIN PIPELINE                                                     ║
# ╚════════════════════════════════════════════════════════════════════╝

def run_pipeline():
    log.info('=' * 64)
    log.info('  SMAP L4 (SPL4SMGP) ROI Extraction Pipeline')
    log.info('=' * 64)

    # ── 1. Authenticate ──────────────────────────────────────────────
    earthaccess.login(strategy="interactive")
    log.info('NASA Earthdata authenticated.')

    # ── 2. Determine operational date range ──────────────────────────
    # SPL4SMGP starts 2015-03-31. Search for the latest available.
    MISSION_START = datetime(2015, 3, 31)
    today = datetime.now()

    # ── 3. Build ROI mask ────────────────────────────────────────────
    roi_geom = load_roi_geometry(SHAPEFILE_DIR, simplify_deg=ROI_SIMPLIFY_DEG)

    # We need one sample HDF5 to read the EASE-2 grid coordinates.
    # Download exactly 1 granule to get the grid.
    log.info('Downloading 1 sample granule for grid coordinates...')
    sample_results = earthaccess.search_data(
        short_name=SHORT_NAME,
        temporal=('2020-01-01', '2020-01-01'),
        count=1,
    )
    if not sample_results:
        raise RuntimeError('Cannot find any SPL4SMGP granules!')

    os.makedirs(LOCAL_RAW_DIR, exist_ok=True)
    sample_files = earthaccess.download(sample_results, local_path=LOCAL_RAW_DIR)
    sample_h5 = str(sample_files[0])

    row_slice, col_slice, mask_2d, lat_sub, lon_sub = build_roi_mask(
        sample_h5, roi_geom)

    # Clean up sample file
    os.remove(sample_h5)

    # ── 4. Setup output NetCDF ───────────────────────────────────────
    os.makedirs(LOCAL_NC_DIR, exist_ok=True)
    os.makedirs(DRIVE_OUTPUT, exist_ok=True)

    local_nc = os.path.join(LOCAL_NC_DIR, OUTPUT_FILENAME)
    drive_nc = os.path.join(DRIVE_OUTPUT, OUTPUT_FILENAME)

    # Check for resume: copy from Drive if exists locally missing
    if not os.path.exists(local_nc) and os.path.exists(drive_nc):
        log.info('Copying existing output from Drive for resume...')
        shutil.copy2(drive_nc, local_nc)

    resume_dt, resume_idx = get_resume_point(local_nc)

    if resume_dt is None:
        # Create fresh output file
        create_output_netcdf(local_nc, lat_sub, lon_sub, mask_2d,
                             row_slice, col_slice)
        search_start = MISSION_START
    else:
        # Resume from after the last processed timestamp
        search_start = resume_dt + timedelta(hours=3)
        log.info('Resuming from %s (index %d)', search_start.isoformat(), resume_idx)

    # ── 5. Main batch loop ───────────────────────────────────────────
    current_start = search_start
    total_appended = 0
    batch_num = 0

    while current_start < today:
        batch_num += 1
        # Calculate end date for this batch (175 days)
        batch_end = min(current_start + timedelta(days=BATCH_DATES),
                        today)

        log.info('=' * 50)
        log.info('BATCH %d: %s to %s', batch_num,
                 current_start.date(), batch_end.date())
        log.info('=' * 50)

        # Search for granules
        results = search_granules(current_start, batch_end)
        if not results:
            log.info('No granules found for this date range. Moving on.')
            current_start = batch_end + timedelta(days=1)
            continue

        # Download with 4 threads
        log.info('Downloading %d granules (4 threads)...', len(results))
        downloaded = download_granules(results, LOCAL_RAW_DIR,
                                       threads=DOWNLOAD_THREADS)

        # Sort chronologically
        sorted_pairs = sort_files_chronologically(downloaded)
        log.info('Sorted %d valid granules chronologically.', len(sorted_pairs))

        # Filter out already-processed timestamps (for resume safety)
        if resume_dt is not None:
            sorted_pairs = [(dt, f) for dt, f in sorted_pairs if dt > resume_dt]
            log.info('After resume filter: %d granules to process.', len(sorted_pairs))

        # Process granules
        dates_since_flush = 0
        last_date_str = None

        pbar = tqdm(sorted_pairs, desc=f'Batch {batch_num} extraction',
                    unit='granule')

        for dt, h5_path in pbar:
            pbar.set_postfix_str(dt.strftime('%Y-%m-%d %H:%M'))

            # Validate file
            if not validate_h5(h5_path):
                log.warning('Corrupt granule: %s — attempting re-download...',
                            os.path.basename(h5_path))
                os.remove(h5_path)

                # Re-download this specific granule
                retry_results = earthaccess.search_data(
                    short_name=SHORT_NAME,
                    temporal=(dt.strftime('%Y-%m-%d'), dt.strftime('%Y-%m-%d')),
                )
                # Find the matching granule by time
                redownloaded = False
                if retry_results:
                    retry_files = earthaccess.download(retry_results,
                                                       local_path=LOCAL_RAW_DIR)
                    for rf in retry_files:
                        rf_dt = parse_smap_datetime(str(rf))
                        if rf_dt == dt:
                            h5_path = str(rf)
                            if validate_h5(h5_path):
                                redownloaded = True
                            break

                if not redownloaded:
                    log.error('Re-download failed for %s — SKIPPING.',
                              dt.isoformat())
                    continue

            # Extract and append
            data = extract_granule(h5_path, row_slice, col_slice, mask_2d)
            if data is None:
                log.error('Extraction failed for %s — SKIPPING.',
                          os.path.basename(h5_path))
                continue

            append_timestep(local_nc, dt, data)
            total_appended += 1

            # Track dates for flush logic
            current_date_str = dt.strftime('%Y%m%d')
            if current_date_str != last_date_str:
                dates_since_flush += 1
                last_date_str = current_date_str

            # Flush every FLUSH_EVERY dates
            if dates_since_flush >= FLUSH_EVERY:
                log.info('Flushing after %d dates (%d total timesteps)...',
                         dates_since_flush, total_appended)
                # Sync is handled by open/close in append_timestep
                dates_since_flush = 0
                gc.collect()

            # Update resume point for next iteration
            resume_dt = dt

            del data
            gc.collect()

        # End of batch: clean up raw files
        log.info('Cleaning up raw HDF5 files...')
        for _, h5_path in sorted_pairs:
            if os.path.exists(h5_path):
                try:
                    os.remove(h5_path)
                except OSError:
                    pass
        # Also clean any remaining files in raw dir
        for f in glob.glob(os.path.join(LOCAL_RAW_DIR, '*.h5')):
            try:
                os.remove(f)
            except OSError:
                pass

        # Copy updated NC to Drive
        log.info('Copying updated NetCDF to Drive...')
        try:
            robust_drive_copy(local_nc, drive_nc)
        except Exception as e:
            log.error('Drive copy failed: %s (continuing pipeline)', e)

        # Move to next batch
        current_start = batch_end + timedelta(days=1)

        gc.collect()

    # ── 6. Final summary ─────────────────────────────────────────────
    final_size = os.path.getsize(local_nc) / 1024**2 if os.path.exists(local_nc) else 0

    log.info('=' * 64)
    log.info('  PIPELINE COMPLETE')
    log.info('  Total timesteps appended: %d', total_appended)
    log.info('  Output file: %s (%.1f MiB)', drive_nc, final_size)
    log.info('=' * 64)


if __name__ == '__main__':
    run_pipeline()
