!pip install -q xee xarray netcdf4 geopandas geemap pyproj dask shapely
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
high_res_v2.py — Optimized Earth Engine → CF-Compliant NetCDF Exporter
===================================================================

Optimizations:
- Threaded scheduler for parallel downloads
- High-volume EE API endpoint
- Larger spatial chunks (2048)
- Month-by-month processing with lazy concatenation
- Lightweight center-sample pre-write and post-write checks
"""

DATASET_ID = "projects/openet/assets/ensemble/conus/gridmet/monthly/v2_1"
PROJECT_ID = "msugw-503806"

SHAPEFILE_DIR = "/content/ogallala_shp"
ROI_LABEL     = "Ogallala"
DRIVE_ROOT    = "/content/drive/MyDrive/MSUGWB"

YEAR_START = None
YEAR_END   = None

DASK_WORKERS    = 4
CHUNK_XY        = 2048
COMPRESS_LEVEL  = 0
MAX_RETRIES     = 3
FORCE_OVERWRITE = False

BAND_METADATA = {}

# ╔════════════════════════════════════════════════════════════════════╗
# ║  SETUP                                                             ║
# ╚════════════════════════════════════════════════════════════════════╝

import os, sys, glob, time, shutil, logging, gc, hashlib, json
from datetime import datetime, timezone
import numpy as np
import ee
import xarray as xr
import xee
import geopandas as gpd
import pyproj
import dask


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('gee_export')
logging.getLogger('urllib3.connectionpool').setLevel(logging.ERROR)

dask.config.set(scheduler='threads', num_workers=DASK_WORKERS)

ee.Authenticate()
ee.Initialize(project=PROJECT_ID, opt_url='https://earthengine-highvolume.googleapis.com')
log.info('Earth Engine ready (high-volume) (project=%s)', PROJECT_ID)

# ╔════════════════════════════════════════════════════════════════════╗
# ║  EXCEPTIONS                                                        ║
# ╚════════════════════════════════════════════════════════════════════╝

class DataQualityError(Exception):
    """Non-retryable data quality error."""
    pass

# ╔════════════════════════════════════════════════════════════════════╗
# ║  PIPELINE FUNCTIONS                                                ║
# ╚════════════════════════════════════════════════════════════════════╝

def load_roi(shapefile_dir, simplify_deg=0.001):
    '''Load shapefile -> dissolved EPSG:4326 geometry.'''
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

def resolve_native_crs(collection):
    '''Detect native CRS + nominal scale.  Returns (epsg_str, metres).'''
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
    '''Compute crs_transform + shape_2d for the actual ROI subset.'''
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

def detect_band_types(collection, bands):
    '''Auto-detect native dtypes from GEE ee.Image.bandTypes().'''
    detected = {}
    try:
        img = collection.first()
        bt = img.bandTypes().getInfo()
        for band_name in bands:
            if band_name not in bt:
                continue
            info = bt[band_name]
            precision = info.get('precision', 'float')
            mn = info.get('min', 0)
            mx = info.get('max', 0)

            if precision == 'int':
                if mn >= 0 and mx <= 255:
                    dtype = 'uint8'
                elif mn >= 0 and mx <= 65535:
                    dtype = 'uint16'
                elif mn >= -32768 and mx <= 32767:
                    dtype = 'int16'
                else:
                    dtype = 'int32'
            elif precision == 'double':
                dtype = 'float64'
            else:
                dtype = 'float32'

            detected[band_name] = {'dtype': dtype}
            log.info('  Band %-20s : %s (GEE: %s, range [%s, %s])',
                     band_name, dtype, precision, mn, mx)
    except Exception as e:
        log.warning('Cannot auto-detect band types: %s', e)
    return detected

def merge_band_info(bands, user_metadata, auto_detected):
    '''Merge user-provided BAND_METADATA with auto-detected types.'''
    result = {}
    warned = False
    for band in bands:
        info = {}
        if band in auto_detected:
            info.update(auto_detected[band])
        if band in user_metadata:
            info.update(user_metadata[band])
        if 'dtype' not in info:
            info['dtype'] = 'float32'
        if '_FillValue' not in info:
            dt = np.dtype(info['dtype'])
            if np.issubdtype(dt, np.unsignedinteger):
                info['_FillValue'] = int(np.iinfo(dt).max)
            elif np.issubdtype(dt, np.signedinteger):
                info['_FillValue'] = int(np.iinfo(dt).min)
            else:
                info['_FillValue'] = None
        if 'scale_factor' not in info and not warned:
            if band not in user_metadata:
                log.warning('BAND_METADATA missing for "%s" (and possibly others).', band)
                warned = True
        result[band] = info
    return result

def strip_xee_encoding(ds):
    '''Clear ALL XEE encoding (the bogus CRS-derived scale_factor).'''
    for name in list(ds.data_vars) + list(ds.coords):
        if name in ds:
            ds[name].encoding.clear()
    return ds

def apply_band_metadata(ds, band_info):
    '''Write scale_factor, add_offset, valid_range, units, long_name as attrs.'''
    for var in list(ds.data_vars):
        if var not in band_info:
            continue
        info = band_info[var]
        if 'scale_factor' in info:
            ds[var].attrs['scale_factor'] = np.float64(info['scale_factor'])
        if 'add_offset' in info:
            ds[var].attrs['add_offset'] = np.float64(info['add_offset'])
        if 'valid_range' in info:
            ds[var].attrs['valid_range'] = np.array(
                info['valid_range'], dtype=info.get('dtype', 'float32'))
        if 'long_name' in info:
            ds[var].attrs['long_name'] = info['long_name']
        if 'units' in info:
            ds[var].attrs['units'] = info['units']
        if 'flag_meanings' in info:
            ds[var].attrs['flag_meanings'] = info['flag_meanings']
    return ds

def prepare_for_write(ds, band_info):
    '''Cast variables to their target dtype, filling NaN → _FillValue.'''
    for var in list(ds.data_vars):
        if var not in band_info:
            continue
        info = band_info[var]
        target_dtype = np.dtype(info['dtype'])
        fill_val = info.get('_FillValue')
        if np.issubdtype(target_dtype, np.integer):
            if fill_val is None:
                fill_val = int(np.iinfo(target_dtype).min)
            ds[var] = ds[var].fillna(fill_val).astype(target_dtype)
        else:
            ds[var] = ds[var].astype(target_dtype)
    return ds

def build_encoding(ds, band_info, compress_level=4, chunk_xy=2048):
    '''Build NetCDF encoding dict with per-band dtype and _FillValue.'''
    enc = {}
    for var in ds.data_vars:
        info = band_info.get(var, {})
        dtype = info.get('dtype', 'float32')
        fill = info.get('_FillValue')
        if fill is None and np.issubdtype(np.dtype(dtype), np.floating):
            fill = np.float32(np.nan) if dtype == 'float32' else np.float64(np.nan)
        if fill is not None:
            try:
                fill = np.dtype(dtype).type(fill)
            except (ValueError, OverflowError):
                pass
        var_chunks = []
        for dim in ds[var].dims:
            dim_size = ds[var].sizes[dim]
            if dim == 'time':
                var_chunks.append(min(1, dim_size))
            elif dim in ('x', 'y'):
                var_chunks.append(min(chunk_xy, dim_size))
            else:
                var_chunks.append(min(1, dim_size))
        enc[var] = {
            'dtype': dtype,
            '_FillValue': fill,
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
            enc[c] = {'dtype': 'float64', '_FillValue': None}
    return enc

def add_grid_mapping(ds, grid_crs, crs_transform, grid_shape):
    '''Add a CF-compliant grid_mapping variable (spatial_ref).'''
    crs_obj = pyproj.CRS.from_user_input(grid_crs)
    cf_attrs = crs_obj.to_cf()
    wkt = crs_obj.to_wkt()
    cf_attrs['crs_wkt'] = wkt
    cf_attrs['spatial_ref'] = wkt
    geo_transform = (
        f'{crs_transform[2]} {crs_transform[0]} {crs_transform[1]} '
        f'{crs_transform[5]} {crs_transform[3]} {crs_transform[4]}'
    )
    cf_attrs['GeoTransform'] = geo_transform
    width, height = grid_shape
    cf_attrs['grid_width']  = int(width)
    cf_attrs['grid_height'] = int(height)
    ds['spatial_ref'] = xr.DataArray(data=np.int32(0), attrs=cf_attrs)
    for var in list(ds.data_vars):
        if var != 'spatial_ref':
            ds[var].attrs['grid_mapping'] = 'spatial_ref'
    return ds

def make_cf_compliant(ds, dataset_id, roi_label, grid_crs, crs_transform, grid_shape, year):
    '''Stamp CF-1.8 global + coordinate attributes + grid_mapping.'''
    width, height = grid_shape
    ds.attrs.update({
        'Conventions': 'CF-1.8',
        'title':  f'{dataset_id} -- {roi_label} ({year})',
        'source': f'Google Earth Engine: {dataset_id}',
        'history': f'Created {datetime.now(timezone.utc).isoformat()}',
        'crs': grid_crs,
        'geospatial_bounds_crs': grid_crs,
        'grid_width': int(width),
        'grid_height': int(height),
    })
    if 'time' in ds.coords:
        ds['time'].attrs.update(axis='T', standard_name='time')
    try:
        is_proj = (pyproj.CRS.from_user_input(grid_crs).axis_info[0].unit_name == 'metre')
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
    ds = add_grid_mapping(ds, grid_crs, crs_transform, grid_shape)
    return ds

def sanitize_attrs(ds):
    '''Flatten complex GEE metadata; drop None values.'''
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
    to_drop = []
    for key, val in ds.attrs.items():
        if val is None:
            to_drop.append(key)
        elif not isinstance(val, (str, int, float, np.number, np.ndarray)):
            ds.attrs[key] = str(val)
    for key in to_drop:
        del ds.attrs[key]
    return ds

def pre_write_check(ds, year, month=None):
    '''Lightweight pre-check: 5x5 center sample.'''
    label = f'[{year}-{month:02d}]' if month else f'[{year}]'
    for var_name in ds.data_vars:
        if var_name == 'spatial_ref':
            continue
        da = ds[var_name]
        ny, nx = da.sizes.get('y',1), da.sizes.get('x',1)
        indexers = {'time': 0} if 'time' in da.dims else {}
        indexers['y'] = slice(max(0, ny//2 - 2), min(ny, ny//2 + 3))
        indexers['x'] = slice(max(0, nx//2 - 2), min(nx, nx//2 + 3))
        
        sample = da.isel(**indexers).compute()
        vals = sample.values
        
        if np.issubdtype(vals.dtype, np.floating):
            nf = int(np.isfinite(vals).sum())
        else:
            fill = da.attrs.get('_FillValue', da.encoding.get('_FillValue', None))
            if fill is not None:
                nf = int((vals != fill).sum())
            else:
                nf = int(vals.size)
                
        if nf == 0:
            raise DataQualityError(f'{label} Variable "{var_name}" center sample all-fill/NaN.')
        log.info('%s Pre-check %s: %d valid center pixels OK', label, var_name, nf)

def post_write_verify(filepath, expected_vars, band_info, year, month=None):
    '''Lightweight post-check: opens file, checks metadata and 5x5 center sample.'''
    label = f'[{year}-{month:02d}]' if month else f'[{year}]'
    try:
        vds = xr.open_dataset(filepath, chunks='auto', mask_and_scale=False)
    except Exception as e:
        return False, f'Cannot reopen: {e}'

    try:
        missing = [v for v in expected_vars if v not in vds.data_vars]
        if missing:
            return False, f'Missing vars: {missing}'
        if 'time' not in vds.dims:
            return False, 'No time dimension'
        if 'spatial_ref' not in vds:
            return False, 'Missing spatial_ref grid_mapping variable'
            
        for c in ('x', 'y'):
            if c in vds.coords:
                vals = vds[c].values
                if np.issubdtype(vals.dtype, np.floating):
                    if not np.all(np.isfinite(vals)):
                        return False, f'Non-finite coord: {c}'
                        
        bad_vars = []
        for v in expected_vars:
            da = vds[v]
            info = band_info.get(v, {})
            fill_val = info.get('_FillValue')
            
            ny, nx = da.sizes.get('y', 1), da.sizes.get('x', 1)
            indexers = {'time': 0} if 'time' in da.dims else {}
            indexers['y'] = slice(max(0, ny//2 - 2), min(ny, ny//2 + 3))
            indexers['x'] = slice(max(0, nx//2 - 2), min(nx, nx//2 + 3))
            
            samp = da.isel(**indexers).compute().values
            
            if np.issubdtype(samp.dtype, np.floating):
                valid = np.any(np.isfinite(samp))
            else:
                if fill_val is not None:
                    valid = np.any(samp != fill_val)
                else:
                    valid = samp.size > 0
                    
            if not valid:
                bad_vars.append(v)
                
        if bad_vars:
            return False, f'No valid data in center on re-read: {bad_vars}'
            
        return True, f'All {len(expected_vars)} vars verified'
    finally:
        vds.close()

def get_md5(filepath, chunk_bytes=10 * 1024 * 1024):
    '''MD5 with standard buffered I/O (no O_DIRECT — FUSE incompatible).'''
    h = hashlib.md5()
    with open(filepath, 'rb') as f:
        while True:
            chunk = f.read(chunk_bytes)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()

def robust_drive_copy(src, dst, max_retries=3):
    '''Copy to Drive, flush, verify MD5.'''
    src_md5 = get_md5(src)
    for attempt in range(1, max_retries + 1):
        try:
            log.info('Drive copy attempt %d/%d ...', attempt, max_retries)
            with open(src, 'rb') as f_src, open(dst, 'wb') as f_dst:
                shutil.copyfileobj(f_src, f_dst, length=10 * 1024 * 1024)
                f_dst.flush()
                os.fsync(f_dst.fileno())
            try:
                os.sync()
            except AttributeError:
                pass  # os.sync() not available on Windows
            time.sleep(3)
            src_size = os.path.getsize(src)
            dst_size = os.path.getsize(dst)
            if src_size != dst_size:
                raise IOError(f'Size mismatch: src={src_size}, dst={dst_size}')
            dst_md5 = get_md5(dst)
            if src_md5 != dst_md5:
                raise IOError('MD5 mismatch')
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
                raise DataQualityError(f'Drive copy failed after {max_retries} attempts: {e}')

# ╔════════════════════════════════════════════════════════════════════╗
# ║  MAIN PIPELINE                                                     ║
# ╚════════════════════════════════════════════════════════════════════╝

def run_pipeline():
    log.info('=' * 62)
    log.info('  Dataset : %s', DATASET_ID)
    log.info('=' * 62)

    coll = ee.ImageCollection(DATASET_ID)
    bands = coll.first().bandNames().getInfo()
    if not bands:
        raise ValueError(f'{DATASET_ID} returned no bands')

    try:
        y0 = YEAR_START or int(ee.Date(coll.sort('system:time_start').first().get('system:time_start')).get('year').getInfo())
        y1 = YEAR_END or int(ee.Date(coll.sort('system:time_start', False).first().get('system:time_start')).get('year').getInfo())
    except Exception as exc:
        raise ValueError(f'Cannot determine year range: {exc}') from exc

    years = list(range(y0, y1 + 1))
    log.info('Bands (%d) : %s', len(bands), bands)
    log.info('Years      : %d - %d  (%d total)', y0, y1, len(years))

    grid_crs, nominal_m = resolve_native_crs(coll.select(bands))
    tolerance_deg = (nominal_m * 1.5) / 111_320.0
    roi = load_roi(SHAPEFILE_DIR, simplify_deg=tolerance_deg)

    # Use bounding box for ee.Geometry — the full polygon (120K+ vertices)
    # exceeds GEE's geometry limit.  crs_transform + shape_2d already
    # defines the exact spatial extent, so a bbox is sufficient for
    # server-side clipping.
    minx, miny, maxx, maxy = roi.bounds
    ee_roi = ee.Geometry.Rectangle([minx, miny, maxx, maxy])

    short = DATASET_ID.split('/')[-1]
    out_dir = os.path.join(DRIVE_ROOT, short)
    os.makedirs(out_dir, exist_ok=True)
    log.info('Output : %s', out_dir)

    grid_params = build_grid_params(roi, grid_crs, nominal_m)
    crs_transform = grid_params['crs_transform']
    grid_shape    = grid_params['shape_2d']

    log.info('Detecting band types from GEE ...')
    auto_types = detect_band_types(coll.select(bands), bands)
    band_info  = merge_band_info(bands, BAND_METADATA, auto_types)

    log.info('Band info summary:')
    for b in bands:
        info = band_info[b]
        sf = info.get('scale_factor', '—')
        fv = info.get('_FillValue', '—')
        log.info('  %-20s  dtype=%-7s  scale=%-8s  fill=%s', b, info['dtype'], sf, fv)

    n_ok = n_fail = n_skip = 0

    for year in years:
        fname = f'{short}_{year}_{ROI_LABEL}.nc'
        final = os.path.join(out_dir, fname)

        if os.path.exists(final) and not FORCE_OVERWRITE:
            sz = os.path.getsize(final)
            if sz < 1024:
                log.warning('[%d] Tiny file (%d B) — re-processing', year, sz)
                os.remove(final)
            else:
                is_valid = False
                try:
                    with xr.open_dataset(final, mask_and_scale=False) as chk:
                        if ('time' in chk.dims and chk.sizes['time'] > 0 and 'spatial_ref' in chk):
                            is_valid = True
                except Exception:
                    pass

                if is_valid:
                    log.info('[%d] On Drive (%.1f MiB), valid — skip', year, sz / 1024**2)
                    n_skip += 1
                    continue
                else:
                    log.warning('[%d] Existing file missing metadata — re-processing', year)
                    try:
                        os.remove(final)
                    except OSError:
                        pass

        month_files = []
        for month in range(1, 13):
            tmp_month = os.path.join('/content', f'{short}_{year}_{month:02d}_tmp.nc')
            
            if os.path.exists(tmp_month):
                log.info('[%d-%02d] Found existing temp file %s, skipping download', year, month, tmp_month)
                month_files.append(tmp_month)
                continue
                
            done = False
            for attempt in range(1, MAX_RETRIES + 1):
                try:
                    log.info('[%d-%02d] attempt %d/%d', year, month, attempt, MAX_RETRIES)
                    month_col = (
                        ee.ImageCollection(DATASET_ID)
                        .filter(ee.Filter.calendarRange(year, year, 'year'))
                        .filter(ee.Filter.calendarRange(month, month, 'month'))
                        .select(bands)
                        .sort('system:time_start')
                    )

                    n_img = month_col.size().getInfo()
                    if n_img == 0:
                        log.info('[%d-%02d] 0 images — skip month', year, month)
                        done = True
                        break

                    log.info('[%d-%02d] %d image(s)', year, month, n_img)

                    grid_w, grid_h = grid_shape
                    safe_chunk = max(16, min(CHUNK_XY, grid_w, grid_h))

                    ds = xr.open_dataset(
                        month_col,
                        engine='ee',
                        geometry=ee_roi,
                        chunks={'x': safe_chunk, 'y': safe_chunk},
                        **grid_params,
                    )

                    _, unique_idx = np.unique(ds['time'], return_index=True)
                    if len(unique_idx) != ds.sizes['time']:
                        ds = ds.isel(time=sorted(unique_idx))

                    ds = strip_xee_encoding(ds)
                    ds = apply_band_metadata(ds, band_info)
                    ds = prepare_for_write(ds, band_info)
                    ds = sanitize_attrs(ds)
                    ds = make_cf_compliant(ds, DATASET_ID, ROI_LABEL, grid_crs, crs_transform, grid_shape, year)
                    
                    pre_write_check(ds, year, month)
                    
                    enc = build_encoding(ds, band_info, COMPRESS_LEVEL, safe_chunk)

                    log.info('[%d-%02d] Writing -> %s', year, month, tmp_month)
                    t0 = time.time()

                    with dask.config.set(scheduler='threads', num_workers=DASK_WORKERS):
                        ds.to_netcdf(
                            tmp_month,
                            engine='netcdf4',
                            encoding=enc,
                            unlimited_dims=['time'],
                        )

                    dt = time.time() - t0
                    mb = os.path.getsize(tmp_month) / 1024**2
                    log.info('[%d-%02d] Written in %.0fs  (%.1f MiB)', year, month, dt, mb)

                    if os.path.getsize(tmp_month) < 1024:
                        raise DataQualityError(f'[{year}-{month:02d}] File too small')

                    ok_flag, msg = post_write_verify(tmp_month, bands, band_info, year, month)
                    if not ok_flag:
                        raise DataQualityError(f'[{year}-{month:02d}] INTEGRITY FAIL: {msg}')
                        
                    log.info('[%d-%02d] Integrity OK: %s', year, month, msg)
                    
                    month_files.append(tmp_month)
                    done = True
                    break

                except DataQualityError as e:
                    log.error('[%d-%02d] DATA ERROR (not retrying): %s', year, month, e)
                    break
                except Exception as e:
                    log.error('[%d-%02d] attempt %d error: %s', year, month, attempt, e)
                    if attempt < MAX_RETRIES:
                        wait = 30 * 2 ** (attempt - 1)
                        log.info('[%d-%02d] retry in %ds ...', year, month, wait)
                        time.sleep(wait)
                    else:
                        log.error('[%d-%02d] retries exhausted', year, month)
                finally:
                    if not done and os.path.exists(tmp_month):
                        try:
                            os.remove(tmp_month)
                        except OSError:
                            pass
                    gc.collect()

        if len(month_files) == 0:
            log.warning('[%d] No monthly files successfully downloaded, skipping concatenation.', year)
            n_fail += 1
            continue

        log.info('[%d] Concatenating %d monthly files...', year, len(month_files))
        tmp_yearly = os.path.join('/content', fname)
        try:
            ds_yearly = xr.open_mfdataset(
                month_files, concat_dim='time', combine='nested',
                mask_and_scale=False,  # keep raw packed values
            )
            # Re-apply encoding so the yearly file preserves dtypes,
            # _FillValue, chunksizes, and compression from the monthly files.
            yearly_enc = build_encoding(ds_yearly, band_info, COMPRESS_LEVEL,
                                        max(16, min(CHUNK_XY, *grid_shape)))
            log.info('[%d] Writing yearly file...', year)
            ds_yearly.to_netcdf(
                tmp_yearly, engine='netcdf4', encoding=yearly_enc,
                unlimited_dims=['time'],
            )
            ds_yearly.close()
            
            log.info('[%d] Exporting to Drive...', year)
            t_mv = time.time()
            robust_drive_copy(tmp_yearly, final, max_retries=3)
            os.remove(tmp_yearly)
            log.info('[%d] Exported in %.0fs -> %s', year, time.time() - t_mv, final)
            
            for f in month_files:
                os.remove(f)
                
            n_ok += 1
        except Exception as e:
            log.error('[%d] Failed during concatenation or drive copy: %s', year, e)
            n_fail += 1
            if os.path.exists(tmp_yearly):
                os.remove(tmp_yearly)

    log.info('=' * 62)
    log.info('  DONE  %d ok | %d fail | %d skip  (%d years)', n_ok, n_fail, n_skip, len(years))
    log.info('=' * 62)

    if os.path.isdir(out_dir):
        nc_files = sorted(f for f in os.listdir(out_dir) if f.endswith('.nc'))
        if nc_files:
            log.info('Files on Drive (%s):', out_dir)
            for f in nc_files:
                sz = os.path.getsize(os.path.join(out_dir, f)) / 1024**2
                log.info('  %s  (%.1f MiB)', f, sz)


if __name__ == '__main__':
    run_pipeline()
