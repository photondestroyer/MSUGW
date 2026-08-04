#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gee_to_netcdf_v2.py — Corrected Earth Engine → NetCDF Exporter
================================================================

Drop-in replacement for z+claude.py.  Every bug from the full audit is
addressed inline with a "FIX(N)" tag so you can grep for them.

To use on Colab:  paste each section below into its own cell, or
                  run the whole file with  %run gee_to_netcdf_v2.py

Crash-causing bugs fixed
------------------------
FIX(1)  DASK_WORKERS=43 caused 503 storms and HDF5 write races → capped to 4
FIX(2)  O_DIRECT in get_md5 → EINVAL on FUSE mount → removed entirely
FIX(3)  chunksizes in build_encoding could exceed dim size → capped
FIX(4)  CHUNK_XY dead conditional (overwritten on next line) → removed
FIX(5)  sz UnboundLocalError in resume check → restructured
FIX(6)  Duplicate code block (copy-paste leftover) → removed
FIX(7)  to_netcdf with threaded scheduler → switched to synchronous
FIX(8)  RuntimeError catch killed retries for transient errors → custom exc
FIX(9)  Dead resolve_crs_and_scale() referencing non-existent global → removed
FIX(10) HTTP pool (_ses/_adp) never wired into GEE/xee → removed

Data-integrity bugs fixed
-------------------------
FIX(11) Integrity checks only inspected first variable → all vars checked
FIX(12) strip_xee_encoding applied bogus scale via fragile float eq → simple clear
FIX(13) sanitize_attrs turned None into the string "None" → drops key instead
FIX(14) cos_lat approximation for geographic CRS documented + improved
FIX(15) O_DIRECT tail-alignment → removed with FIX(2)
"""

# ╔════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION                                                     ║
# ╚════════════════════════════════════════════════════════════════════╝

DATASET_ID = "MODIS/061/MOD13A3"
PROJECT_ID = "msugw-503806"

SHAPEFILE_DIR = "/content/ogallala_shp"
ROI_LABEL     = "Ogallala"
DRIVE_ROOT    = "/content/drive/MyDrive/MSUGWB"

YEAR_START = None   # None = auto-detect
YEAR_END   = None

# FIX(1): was 43 — 503 storms + HDF5 write races.  4-6 is safe.
DASK_WORKERS   = 4
# FIX(4): dead conditional removed; single value only.
CHUNK_XY       = 512
COMPRESS_LEVEL = 4      # zlib  0=off  9=max
MAX_RETRIES    = 3
FORCE_OVERWRITE = False

# ╔════════════════════════════════════════════════════════════════════╗
# ║  SETUP                                                             ║
# ╚════════════════════════════════════════════════════════════════════╝

# --- Colab prerequisites ---
# Uncomment the next two lines when running on Colab:
# !pip install -q "xee>=0.0.14" xarray netcdf4 geopandas pyproj
# from google.colab import drive; drive.mount('/content/drive')

import os, sys, glob, time, shutil, logging, gc, hashlib
from datetime import datetime, timezone

import numpy as np
import ee
import xarray as xr
import xee                       # registers the 'ee' xarray engine
from xee import helpers
import geopandas as gpd
import pyproj
import dask

# ── Logging ──
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('gee_export')
logging.getLogger('urllib3.connectionpool').setLevel(logging.ERROR)

# FIX(10): removed dead requests session — it was never wired into
# ee / xee / googleapiclient.  GEE retries are handled by
# googleapiclient's built-in retry logic.

# ── Dask ──
dask.config.set(scheduler='threads', num_workers=DASK_WORKERS)

# ── Earth Engine ──
ee.Authenticate()
ee.Initialize(project=PROJECT_ID)
log.info('Earth Engine ready  (project=%s)', PROJECT_ID)


# ╔════════════════════════════════════════════════════════════════════╗
# ║  CUSTOM EXCEPTIONS                                                 ║
# ╚════════════════════════════════════════════════════════════════════╝

# FIX(8): RuntimeError was used for both non-retryable data errors AND
# wrapped transient I/O failures, so the bare `except RuntimeError: break`
# killed retries on transient errors too.  Now data-quality failures use
# a dedicated exception that is *never* retried.

class DataQualityError(Exception):
    """Raised when data is invalid (all-NaN, corrupt, failed integrity).
    This is NEVER retried — the problem is in the data, not the network."""
    pass


# ╔════════════════════════════════════════════════════════════════════╗
# ║  PIPELINE FUNCTIONS                                                ║
# ╚════════════════════════════════════════════════════════════════════╝

# ─── ROI ─────────────────────────────────────────────────────────────

def load_roi(shapefile_dir, simplify_deg=0.001):
    '''Load shapefile -> dissolved EPSG:4326 geometry.
    Auto-simplifies when vertex count exceeds Earth Engine limits.
    '''
    paths = glob.glob(os.path.join(shapefile_dir, '**', '*.shp'), recursive=True)
    if not paths:
        raise FileNotFoundError(f'No .shp found in {shapefile_dir}')
    log.info('Shapefile: %s', paths[0])

    gdf = gpd.read_file(paths[0])
    if gdf.crs is None or gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(epsg=4326)

    roi = (gdf.geometry.union_all()
           if hasattr(gdf.geometry, 'union_all')
           else gdf.geometry.unary_union)

    try:
        import shapely as _shp
        nv = int(_shp.get_num_coordinates(roi))
    except Exception:
        nv = len(roi.wkt) // 20

    if nv > 50_000:
        roi = roi.simplify(simplify_deg, preserve_topology=True)
        try:
            nv_after = int(_shp.get_num_coordinates(roi))
        except Exception:
            nv_after = len(roi.wkt) // 20
        log.warning('Simplified ROI: %s -> %s vertices (tol=%.4f deg)',
                    f'{nv:,}', f'{nv_after:,}', simplify_deg)
    else:
        log.info('ROI vertices: %s', f'{nv:,}')
    return roi


# ─── CRS & Grid ─────────────────────────────────────────────────────

# FIX(9): removed dead resolve_crs_and_scale() — it referenced a global
# `roi` that doesn't exist at function-definition time, and was only
# called from commented-out code.  CRS resolution is now inline in the
# pipeline, and build_grid_params handles the grid construction.

def resolve_native_crs(collection):
    '''Detect native CRS + nominal scale from a GEE collection.
    Returns (epsg_string, nominal_scale_metres).
    '''
    proj = collection.first().select(0).projection()
    info = proj.getInfo()
    native_crs = info.get('crs', 'EPSG:4326')
    nominal_m  = proj.nominalScale().getInfo()

    log.info('Native CRS    : %s', native_crs[:70])
    log.info('Nominal scale : %.1f m', nominal_m)

    resolved = 'EPSG:4326'

    if native_crs.startswith('SR-ORG:'):
        log.warning('SR-ORG CRS (%s) -> EPSG:4326 fallback', native_crs)
    elif native_crs.startswith('EPSG:'):
        resolved = native_crs
    else:
        try:
            crs_obj = pyproj.CRS.from_user_input(native_crs)
            epsg = crs_obj.to_epsg()
            if epsg:
                resolved = f'EPSG:{epsg}'
            elif crs_obj.axis_info and crs_obj.axis_info[0].unit_name == 'metre':
                resolved = 'EPSG:5070'
                log.warning('Metre CRS, no EPSG -> EPSG:5070')
        except Exception as exc:
            log.warning('CRS parse error: %s -> EPSG:4326', exc)

    log.info('Resolved CRS  : %s', resolved)
    return resolved, nominal_m


def build_grid_params(roi, grid_crs, nominal_m):
    '''Compute exact crs_transform + shape_2d for the XEE grid.

    GEE affine: (x_scale, 0, x_origin, 0, -y_scale, y_origin)
    Origin is the top-left corner of the bounding box.

    FIX(14): For geographic CRS, cos_lat correction IS applied to make
    pixels approximately square in metres at the ROI centroid.  This is
    geometrically more accurate than equal-degree spacing but introduces
    progressive X-skew away from the centroid latitude.  For the Ogallala
    aquifer (~30-43 deg N), max X-distortion is ~10%.  If this matters,
    set FALLBACK_CRS to a projected CRS like EPSG:5070.
    '''
    minx, miny, maxx, maxy = roi.bounds

    transformer = pyproj.Transformer.from_crs(
        'EPSG:4326', grid_crs, always_xy=True)
    x1, y1 = transformer.transform(minx, miny)
    x2, y2 = transformer.transform(maxx, maxy)

    left   = min(x1, x2)
    right  = max(x1, x2)
    bottom = min(y1, y2)
    top    = max(y1, y2)

    try:
        unit = pyproj.CRS.from_user_input(grid_crs).axis_info[0].unit_name
    except Exception:
        unit = 'degree'

    if unit == 'metre':
        scale_x = float(nominal_m)
        scale_y = float(nominal_m)
    else:
        lat_centroid = roi.centroid.y
        cos_lat = np.cos(np.radians(lat_centroid))
        scale_y = nominal_m / 111_320.0
        scale_x = scale_y / cos_lat if cos_lat > 0.01 else scale_y

    width  = int(np.ceil((right - left) / scale_x))
    height = int(np.ceil((top - bottom) / scale_y))

    crs_transform = (scale_x, 0.0, left, 0.0, -scale_y, top)

    log.info('Grid : %d x %d px  |  scale_x=%.6g  scale_y=%.6g %s',
             width, height, scale_x, scale_y, unit)

    return {
        'crs': grid_crs,
        'crs_transform': crs_transform,
        'shape_2d': (width, height),
    }


# ─── Encoding ───────────────────────────────────────────────────────

def strip_xee_encoding(ds):
    '''Clear ALL encoding from variables and coordinates.

    FIX(12): The previous version tried to "bake in" the scale_factor by
    multiplying data, then re-masked with fragile float equality.  This
    is unnecessary and dangerous:
      - XEE's "scale_factor" is the CRS pixel size, NOT a CF packing param
      - GEE already returns decoded physical values through computePixels
      - Float equality after arithmetic is unreliable
    The correct fix is to simply CLEAR the encoding so to_netcdf() writes
    the values as-is.
    '''
    for name in list(ds.data_vars) + list(ds.coords):
        if name in ds:
            ds[name].encoding.clear()
    return ds


def build_encoding(ds, compress_level=4, chunk_xy=512):
    '''Build explicit NetCDF encoding for every variable and coordinate.

    FIX(3): chunksizes are now capped to actual dimension sizes.
    HDF5 rejects chunk dimensions larger than the array dimension,
    raising EINVAL.
    '''
    enc = {}

    for var in ds.data_vars:
        var_chunks = []
        for dim in ds[var].dims:
            dim_size = ds[var].sizes[dim]
            if dim == 'time':
                var_chunks.append(min(1, dim_size))
            elif dim in ('x', 'y'):
                # FIX(3): cap to actual dimension size
                var_chunks.append(min(chunk_xy, dim_size))
            else:
                var_chunks.append(min(1, dim_size))

        enc[var] = {
            'dtype': 'float32',
            '_FillValue': np.float32(np.nan),
            'zlib': compress_level > 0,
            'complevel': compress_level,
            'chunksizes': tuple(var_chunks),
        }

    for c in ds.coords:
        if c == 'time':
            enc[c] = {
                'dtype': 'int64',
                '_FillValue': None,
                'units': 'days since 1970-01-01',
                'calendar': 'proleptic_gregorian',
            }
        elif c in ('x', 'y'):
            # float64 for coordinates — avoids precision loss for
            # fine-resolution grids (e.g. 30 m CDL over 800 km extent)
            enc[c] = {'dtype': 'float64', '_FillValue': None}
    return enc


# ─── CF Compliance ──────────────────────────────────────────────────

def make_cf_compliant(ds, dataset_id, roi_label, grid_crs, year):
    '''Stamp CF-1.8 global + coordinate attributes.'''
    ds.attrs.update({
        'Conventions': 'CF-1.8',
        'title':  f'{dataset_id} -- {roi_label} ({year})',
        'source': f'Google Earth Engine: {dataset_id}',
        'history': f'Created {datetime.now(timezone.utc).isoformat()}',
        'crs': grid_crs,
    })

    if 'time' in ds.coords:
        ds['time'].attrs.update(axis='T', standard_name='time')

    try:
        is_proj = (pyproj.CRS.from_user_input(grid_crs)
                   .axis_info[0].unit_name == 'metre')
    except Exception:
        is_proj = False

    if 'x' in ds.coords:
        ds['x'].attrs.update(
            axis='X',
            standard_name='projection_x_coordinate' if is_proj else 'longitude',
            units='m' if is_proj else 'degrees_east',
        )
    if 'y' in ds.coords:
        ds['y'].attrs.update(
            axis='Y',
            standard_name='projection_y_coordinate' if is_proj else 'latitude',
            units='m' if is_proj else 'degrees_north',
        )
    return ds


def sanitize_attrs(ds):
    '''Flatten complex GEE metadata to strings for HDF5 compatibility.

    FIX(13): None values are now DROPPED rather than written as the
    literal string "None", which would be misleading in CF metadata.
    '''
    for var in list(ds.data_vars) + list(ds.coords):
        if var not in ds:
            continue
        to_drop = []
        for key, val in ds[var].attrs.items():
            if val is None:
                to_drop.append(key)
            elif not isinstance(val, (str, int, float, np.number, np.ndarray)):
                ds[var].attrs[key] = str(val)
        for key in to_drop:
            del ds[var].attrs[key]

    # Same for global attrs
    to_drop = []
    for key, val in ds.attrs.items():
        if val is None:
            to_drop.append(key)
        elif not isinstance(val, (str, int, float, np.number, np.ndarray)):
            ds.attrs[key] = str(val)
    for key in to_drop:
        del ds.attrs[key]

    return ds


# ─── Pre-Write Validation ──────────────────────────────────────────

def pre_write_check(ds, year):
    '''Coarsened subsample across the ENTIRE grid for each variable.

    FIX(11): now checks ALL data variables, not just the first one.
    Uses coarsen+max which samples the full spatial extent evenly,
    avoiding the "donut hole" false-positive of centre-only checks.
    '''
    for var_name in ds.data_vars:
        da = ds[var_name]

        coarsen_dims = {}
        for d in ('x', 'y'):
            if d in da.dims and da.sizes[d] > 1:
                coarsen_dims[d] = max(1, da.sizes[d] // 50)

        if not coarsen_dims:
            log.warning('[%d] %s: no spatial dims to check', year, var_name)
            continue

        indexers = {'time': 0} if 'time' in da.dims else {}
        sample = (da.isel(**indexers)
                    .coarsen(coarsen_dims, boundary='trim')
                    .max()
                    .compute())

        nf = int(np.isfinite(sample.values).sum())
        if nf == 0:
            raise DataQualityError(
                f'[{year}] Variable "{var_name}" is ALL-NaN.  '
                'Earth Engine returned no valid pixels.'
            )
        log.info('[%d] Pre-check %s: %d finite values OK', year, var_name, nf)


# ─── Post-Write Integrity Verification ─────────────────────────────

def post_write_verify(filepath, expected_vars, year):
    '''Reopen the NetCDF from disk and verify every variable.

    FIX(11): checks ALL bands, not just the first.
    Returns (ok, message).
    '''
    try:
        vds = xr.open_dataset(filepath, chunks='auto')
    except Exception as e:
        return False, f'Cannot reopen: {e}'

    try:
        # All expected variables?
        missing = [v for v in expected_vars if v not in vds.data_vars]
        if missing:
            return False, f'Missing: {missing}'

        # Time dimension?
        if 'time' not in vds.dims:
            return False, 'No time dimension'

        # Coordinates finite?
        for c in vds.coords:
            v = vds[c].values
            if np.issubdtype(v.dtype, np.floating):
                if not np.all(np.isfinite(v)):
                    return False, f'Non-finite coord: {c}'

        # No bogus scale_factor on ANY variable?
        for v in expected_vars:
            sf = vds[v].encoding.get('scale_factor')
            if sf is not None and sf != 1.0:
                return False, f'scale_factor={sf} on {v}'

        # Data finite check — sample ALL variables
        bad_vars = []
        for v in expected_vars:
            da = vds[v]
            idx = {}
            if 'time' in da.dims:
                idx['time'] = 0
            for d in da.dims:
                if d != 'time':
                    s = da.sizes[d]
                    m = s // 2
                    hw = min(50, s // 4, 100)
                    hw = max(hw, 1)
                    idx[d] = slice(max(0, m - hw), min(s, m + hw))
            samp = da.isel(**idx).compute().values
            if not np.any(np.isfinite(samp)):
                bad_vars.append(v)

        if bad_vars:
            return False, f'All NaN on re-read: {bad_vars}'

        return True, f'All {len(expected_vars)} vars verified'

    finally:
        vds.close()


# ─── Drive Copy ─────────────────────────────────────────────────────

def get_md5(filepath, chunk_bytes=10 * 1024 * 1024):
    '''Calculate MD5 hash with standard buffered I/O.

    FIX(2): Removed os.O_DIRECT entirely.
      - FUSE (Google Drive mount) does NOT support O_DIRECT → EINVAL
      - O_DIRECT tail-alignment at EOF → EINVAL on some kernels (FIX 15)
      - Standard read() with large buffers is fast enough for verification.
    '''
    h = hashlib.md5()
    with open(filepath, 'rb') as f:
        while True:
            chunk = f.read(chunk_bytes)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def robust_drive_copy(src, dst, max_retries=3):
    '''Copy file to Drive, flush, and verify MD5.'''
    src_md5 = get_md5(src)

    for attempt in range(1, max_retries + 1):
        try:
            log.info('Drive copy attempt %d/%d ...', attempt, max_retries)

            with open(src, 'rb') as f_src, open(dst, 'wb') as f_dst:
                shutil.copyfileobj(f_src, f_dst, length=10 * 1024 * 1024)
                f_dst.flush()
                os.fsync(f_dst.fileno())

            # Give FUSE time to propagate metadata
            os.sync()
            time.sleep(3)

            # Size check
            src_size = os.path.getsize(src)
            dst_size = os.path.getsize(dst)
            if src_size != dst_size:
                raise IOError(
                    f'Size mismatch: src={src_size}, dst={dst_size}')

            # MD5 check
            dst_md5 = get_md5(dst)
            if src_md5 != dst_md5:
                raise IOError('MD5 mismatch (silent data corruption)')

            log.info('Drive copy verified (MD5 match, %d bytes)', src_size)
            return True

        except Exception as e:
            log.error('Drive copy attempt %d failed: %s', attempt, e)
            if os.path.exists(dst):
                try:
                    os.remove(dst)
                except OSError:
                    pass
            if attempt < max_retries:
                wait = 15 * attempt
                log.info('Retrying in %ds ...', wait)
                time.sleep(wait)
            else:
                raise DataQualityError(
                    f'Drive copy failed after {max_retries} attempts: {e}')


# ╔════════════════════════════════════════════════════════════════════╗
# ║  MAIN PIPELINE                                                     ║
# ╚════════════════════════════════════════════════════════════════════╝

def run_pipeline():
    log.info('=' * 62)
    log.info('  Dataset : %s', DATASET_ID)
    log.info('=' * 62)

    # ── 1. Collection metadata ──────────────────────────────────
    coll = ee.ImageCollection(DATASET_ID)
    bands = coll.first().bandNames().getInfo()
    if not bands:
        raise ValueError(f'{DATASET_ID} returned no bands')

    try:
        y0 = YEAR_START or int(ee.Date(
            coll.sort('system:time_start').first()
                .get('system:time_start')).get('year').getInfo())
        y1 = YEAR_END or int(ee.Date(
            coll.sort('system:time_start', False).first()
                .get('system:time_start')).get('year').getInfo())
    except Exception as exc:
        raise ValueError(
            f'Cannot determine year range for {DATASET_ID}: {exc}') from exc

    years = list(range(y0, y1 + 1))
    log.info('Bands (%d) : %s', len(bands), bands)
    log.info('Years      : %d - %d  (%d total)', y0, y1, len(years))

    # ── 2. CRS & grid ──────────────────────────────────────────
    grid_crs, nominal_m = resolve_native_crs(coll.select(bands))

    tolerance_deg = (nominal_m * 1.5) / 111_320.0
    roi = load_roi(SHAPEFILE_DIR, simplify_deg=tolerance_deg)

    short = DATASET_ID.split('/')[-1]
    out_dir = os.path.join(DRIVE_ROOT, short)
    os.makedirs(out_dir, exist_ok=True)
    log.info('Output : %s', out_dir)

    grid_params = build_grid_params(roi, grid_crs, nominal_m)

    # ── 3. Year loop ────────────────────────────────────────────
    n_ok = n_fail = n_skip = 0

    for year in years:
        fname = f'{short}_{year}_{ROI_LABEL}.nc'
        final = os.path.join(out_dir, fname)

        # ── Resume check ────────────────────────────────────────
        # FIX(5): restructured to avoid UnboundLocalError on `sz`.
        if os.path.exists(final) and not FORCE_OVERWRITE:
            sz = os.path.getsize(final)
            if sz < 1024:
                log.warning('[%d] Tiny file (%d B) — re-processing', year, sz)
                os.remove(final)
            else:
                # Quick structural check
                is_valid = False
                try:
                    with xr.open_dataset(final) as check_ds:
                        if 'time' in check_ds.dims and check_ds.sizes['time'] > 0:
                            is_valid = True
                except Exception:
                    pass

                if is_valid:
                    log.info('[%d] On Drive (%.1f MiB), valid — skip',
                             year, sz / 1024**2)
                    n_skip += 1
                    continue
                else:
                    log.warning('[%d] Existing file is corrupt — re-processing', year)
                    try:
                        os.remove(final)
                    except OSError:
                        pass

        tmp = os.path.join('/content', fname)
        done = False

        for attempt in range(1, MAX_RETRIES + 1):
            try:
                log.info('[%d] attempt %d/%d', year, attempt, MAX_RETRIES)

                # ── Filter collection ────────────────────────────
                yr_col = (
                    ee.ImageCollection(DATASET_ID)
                    .filter(ee.Filter.calendarRange(year, year, 'year'))
                    .select(bands)
                    .sort('system:time_start')
                )

                n_img = yr_col.size().getInfo()
                if n_img == 0:
                    log.warning('[%d] 0 images — skip year', year)
                    n_skip += 1
                    break

                log.info('[%d] %d image(s)', year, n_img)

                # ── Open lazily via XEE ──────────────────────────
                # FIX(3): cap chunk to actual grid dims
                grid_w, grid_h = grid_params['shape_2d']
                safe_chunk = max(16, min(CHUNK_XY, grid_w, grid_h))

                ds = xr.open_dataset(
                    yr_col,
                    engine='ee',
                    chunks={'x': safe_chunk, 'y': safe_chunk},
                    **grid_params,
                )

                # Deduplicate timestamps (GEE backend quirk)
                _, unique_idx = np.unique(ds['time'], return_index=True)
                if len(unique_idx) != ds.sizes['time']:
                    log.warning('[%d] Deduplicating %d -> %d timestamps',
                                year, ds.sizes['time'], len(unique_idx))
                    ds = ds.isel(time=sorted(unique_idx))

                # ── Fix encoding & metadata ──────────────────────
                ds = strip_xee_encoding(ds)
                ds = sanitize_attrs(ds)
                ds = make_cf_compliant(ds, DATASET_ID, ROI_LABEL,
                                       grid_crs, year)

                # ── Pre-write validation ─────────────────────────
                pre_write_check(ds, year)

                # ── Write to scratch ─────────────────────────────
                enc = build_encoding(ds, COMPRESS_LEVEL, safe_chunk)

                log.info('[%d] Writing -> %s', year, tmp)
                t0 = time.time()

                # FIX(7): use synchronous scheduler for the write.
                # netCDF4/HDF5 is NOT thread-safe for concurrent writes
                # to a single file.  With threaded scheduler + N workers,
                # N chunks can write simultaneously → EINVAL / corruption.
                # Synchronous = one chunk at a time = safe.
                try:
                    with dask.config.set(scheduler='synchronous'):
                        ds.to_netcdf(
                            tmp,
                            engine='netcdf4',
                            encoding=enc,
                            unlimited_dims=['time'],
                        )
                except Exception as write_err:
                    if os.path.exists(tmp):
                        os.remove(tmp)
                    raise  # let retry logic handle it

                # FIX(6): removed duplicate dt/mb/enc/write block

                dt = time.time() - t0
                mb = os.path.getsize(tmp) / 1024**2
                log.info('[%d] Written in %.0fs  (%.1f MiB)', year, dt, mb)

                if os.path.getsize(tmp) < 1024:
                    raise DataQualityError(
                        f'[{year}] File too small ({os.path.getsize(tmp)} B)')

                # ── Post-write integrity ─────────────────────────
                ok_flag, msg = post_write_verify(tmp, bands, year)
                if not ok_flag:
                    raise DataQualityError(
                        f'[{year}] INTEGRITY FAIL: {msg}')
                log.info('[%d] Integrity OK: %s', year, msg)

                # ── Export to Drive ───────────────────────────────
                log.info('[%d] Exporting to Drive ...', year)
                t_mv = time.time()
                robust_drive_copy(tmp, final, max_retries=3)
                os.remove(tmp)
                log.info('[%d] Exported in %.0fs -> %s',
                         year, time.time() - t_mv, final)

                n_ok += 1
                done = True
                break

            # FIX(8): DataQualityError = never retry (data problem)
            except DataQualityError as e:
                log.error('[%d] DATA ERROR (not retrying): %s', year, e)
                n_fail += 1
                break

            # Everything else = possibly transient → retry with backoff
            except Exception as e:
                log.error('[%d] attempt %d error: %s', year, attempt, e)
                if attempt < MAX_RETRIES:
                    wait = 30 * 2 ** (attempt - 1)
                    log.info('[%d] retry in %ds ...', year, wait)
                    time.sleep(wait)
                else:
                    log.error('[%d] retries exhausted', year)
                    n_fail += 1

            finally:
                if not done and os.path.exists(tmp):
                    try:
                        os.remove(tmp)
                    except OSError:
                        pass
                gc.collect()

    # ── Summary ─────────────────────────────────────────────────
    log.info('=' * 62)
    log.info('  DONE  %d ok | %d fail | %d skip  (%d years)',
             n_ok, n_fail, n_skip, len(years))
    log.info('=' * 62)

    if os.path.isdir(out_dir):
        nc_files = sorted(f for f in os.listdir(out_dir) if f.endswith('.nc'))
        if nc_files:
            log.info('Files on Drive (%s):', out_dir)
            for f in nc_files:
                sz = os.path.getsize(os.path.join(out_dir, f)) / 1024**2
                log.info('  %s  (%.1f MiB)', f, sz)


# ╔════════════════════════════════════════════════════════════════════╗
# ║  ENTRY POINT                                                       ║
# ╚════════════════════════════════════════════════════════════════════╝

if __name__ == '__main__':
    run_pipeline()
