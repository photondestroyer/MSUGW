"""
GOSIF v2 8-day global GeoTIFFs -> Ogallala-clipped, CF-1.8 NetCDF.

Reads every GOSIF_YYYYDDD.tif.gz in SRC_DIR through rasterio's /vsigzip/
handler (no extraction to disk), windows each scene to the hp_bound2010
bbox snapped to the native 0.05-degree grid, masks cells outside the
aquifer polygon, applies scale/offset + fill handling, and appends the
result to an unlimited time axis. Monthly means are derived afterwards by
streaming the finished 8-day series back in 50-date chunks.

Parallelism (dask, threaded): per-file reads are independent, so each
50-date chunk is fetched concurrently across N_WORKERS threads while the
single main thread keeps NetCDF appends strictly serial and in date order.
Progress: an overall tqdm bar over all dates to ingest plus a tqdm bar for
the monthly pass. Every fetch is bracketed by timestamped heartbeat logs
(chunk size, worker count, files/sec) so a stall is always localizable;
dask's own display thread is deliberately NOT used (it fights tqdm/terminals
for output and can mask a stall as a freeze). If reads ever stall, rerun
with --workers 1 (serial fallback).

Mesh reuse: the ROI window, polygon mask and lat/lon mesh derive only from
the regular 0.05-degree grid + shapefile, so they are computed once, cached
to a small .npz sidecar (validated against the shapefile signature on load),
and reused for every data file and every resumed run - never recomputed
per date.

Memory design: only one 50-date chunk (50 x H x W float32, ~10 MB for the
Ogallala window) plus small monthly accumulators are ever resident. The
global 7200 x 3600 scenes are never fully loaded.

Resume: dates already present in OUT_PATH (matched by calendar date) are
skipped; monthly layers are always recomputed from the complete 8-day
series at the end of a run. Re-running the script is idempotent.

Usage (SSM env: rasterio, netCDF4, geopandas, tqdm + dask recommended):
    conda run -n SSM python G:\\MSU_GWB\\datasets\\process_gosif_ogallala.py
Outputs:
    <OUT_DIR>/GOSIF_Ogallala_8day_monthly.nc (+ GOSIF_Ogallala_mesh.npz)
"""

import glob
import os
import sys
import time
from datetime import date, timedelta

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window, from_bounds
from rasterio.features import rasterize
import geopandas as gpd
import netCDF4 as nc

try:
    from tqdm import tqdm
except ImportError:  # graceful fallback: silent no-op bar, script still runs
    class tqdm:  # noqa: D101
        def __init__(self, *a, **k):
            self.total = k.get("total", 0)
            self.n = 0

        def update(self, n=1):
            self.n += n

        def set_description(self, *a, **k):
            pass

        def close(self):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

try:
    import dask
    _HAS_DASK = True
except ImportError:  # fallback: serial reads (slower, same results)
    _HAS_DASK = False

# --------------------------------------------------------------------------
# Config (edit here; everything else derives from these + the first file)
# --------------------------------------------------------------------------
SRC_DIR = r"G:\MSU_GWB\New folder"
SHP_PATH = r"G:\MSU_GWB\datasets\HPA_polygon\hp_bound2010.shp"
OUT_DIR = r"G:\MSU_GWB\datasets\merged_datasets"
OUT_NAME = "GOSIF_Ogallala_8day_monthly.nc"
CHUNK_DATES = 100          # dates processed (and NetCDF-synced) per batch
N_WORKERS = 4             # dask threads for parallel per-file reads in a chunk
# (4 is deliberately gentle: 8 concurrent gzip readers stalled drive I/O
#  during testing. Raise to 8 with --workers 8 once a chunk streams cleanly.)
MESH_NAME = "GOSIF_Ogallala_mesh.npz"  # sidecar caching the precomputed ROI mesh
RES = 0.05                # GOSIF native resolution, degrees (asserted per file)
SCALE = 0.0001            # per GOSIF v2 tags
FILL_CODES = (32766, 32767)   # snow/ice, water (per GOSIF v2 tags)
UNITS = "W m-2 um-1 sr-1"     # per GOSIF v2 tags
COMPRESSION = 3
# --------------------------------------------------------------------------

OUT_PATH = os.path.join(OUT_DIR, OUT_NAME)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def decode_time(time_var):
    """Version-proof decode of 'days since 1970-01-01' (avoids
    netCDF4/cftime quirks across environments)."""
    assert str(time_var.units) == "days since 1970-01-01", \
        f"unexpected time units: {time_var.units}"
    return pd.to_datetime(np.asarray(time_var[:], dtype="float64"),
                          unit="D", origin="1970-01-01")


def parse_gosif_date(fname):
    """GOSIF_YYYYDDD.tif.gz -> datetime.date. Raises ValueError if malformed."""
    base = os.path.basename(fname)
    if not (base.startswith("GOSIF_") and base.endswith(".tif.gz")):
        raise ValueError(f"unexpected filename: {base}")
    core = base[len("GOSIF_"):-len(".tif.gz")]
    year, doy = int(core[:4]), int(core[4:])
    if not (1 <= doy <= 366):
        raise ValueError(f"day-of-year out of range in {base}")
    return date(year, 1, 1) + timedelta(days=doy - 1)


def expected_8day_grid(year):
    """Nominal GOSIF 8-day grid DOYs (informational gap check only)."""
    doys = list(range(1, 362, 8))
    if date(year, 12, 31).timetuple().tm_yday == 366:
        doys.append(361)  # already included; leap years: same 46 slots
    return sorted(set(doys))


def vsi(p):
    return "/vsigzip/" + p.replace("\\", "/")


def read_one(f, col_off, row_off, w, h):
    """Read + decode one GOSIF file's ROI window. Module-level so dask
    threaded workers can call it; opens its own dataset handle (thread-safe).
    Returns (ok: bool, image_or_error: np.ndarray | str). Image is float32
    SIF already scaled, polygon-masked and NaN-filled; on failure ok=False
    and the payload is the error message (caller logs + skips, like before).
    Uses module globals FILL_CODES/SCALE/MASK_ARR (read-only sharing)."""
    import rasterio as _rio
    import numpy as _np
    try:
        with _rio.open(vsi(f)) as src:
            if (src.width, src.height) != (7200, 3600):
                return False, f"unexpected shape {(src.width, src.height)}"
            a = src.read(1, window=Window(col_off, row_off, w, h)
                         ).astype("float32")
    except Exception as e:  # noqa: BLE001 - per-file faults must not kill the chunk
        return False, f"unreadable ({e})"
    a[(a == FILL_CODES[0]) | (a == FILL_CODES[1])] = _np.nan
    a = a * _np.float32(SCALE)
    a[~MASK_ARR] = _np.nan
    return True, a


def main():
    try:  # unbuffered stdout so every log line appears instantly, even piped
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass
    t_start = time.time()
    log("process_gosif_ogallala starting "
        f"(dask={'yes' if _HAS_DASK else 'NO - serial fallback'}, "
        f"workers={N_WORKERS if _HAS_DASK else 1})")
    log(f"source dir: {SRC_DIR}")
    files = sorted(glob.glob(os.path.join(SRC_DIR, "GOSIF_*.tif.gz")))
    if not files:
        raise SystemExit("no GOSIF_*.tif.gz files found - nothing to do")
    log(f"found {len(files)} input files")

    # ---- parse + audit dates ----
    dated = []
    bad_names = []
    for f in files:
        try:
            dated.append((parse_gosif_date(f), f))
        except ValueError as e:
            bad_names.append(str(e))
    for b in bad_names:
        log(f"WARNING skipping file: {b}")
    dated.sort()
    by_year = {}
    for d, f in dated:
        by_year.setdefault(d.year, []).append(d)
    for y in sorted(by_year):
        have = sorted({dd.timetuple().tm_yday for dd in by_year[y]})
        want = expected_8day_grid(y)
        missing = [x for x in want if x not in have]
        log(f"year {y}: {len(have)} dates" +
            (f", missing DOYs vs nominal grid: {missing}" if missing else ""))

    # ---- reference grid from first file (strict) ----
    # (vsi() helper is module-level; reused here)
    with rasterio.open(vsi(dated[0][1])) as src:
        assert (src.width, src.height) == (7200, 3600), \
            f"unexpected raster size {(src.width, src.height)}"
        resx, resy = abs(src.transform.a), abs(src.transform.e)
        assert abs(resx - RES) < 1e-9 and abs(resy - RES) < 1e-9, \
            f"unexpected resolution {(resx, resy)}"
        assert src.count == 1 and src.dtypes[0] in ("int16", "uint16"), \
            f"unexpected bands/dtype {(src.count, src.dtypes)}"
        assert str(src.crs) == "EPSG:4326", f"unexpected CRS {src.crs}"
        transform = src.transform
    log("reference grid OK: global 0.05-deg EPSG:4326, 1 x int16 band")

    # ---- ROI mesh: load cached or compute once, reuse for every file ----
    # The window, polygon mask and lat/lon mesh depend only on the regular
    # 0.05-degree grid + shapefile, so they are built a single time and cached
    # in a sidecar .npz (validated against the shapefile signature on load).
    mesh_path = os.path.join(os.path.dirname(OUT_PATH), MESH_NAME)
    shp_sig = (os.path.getmtime(SHP_PATH), os.path.getsize(SHP_PATH))
    gdf = gpd.read_file(SHP_PATH).to_crs("EPSG:4326")
    log(f"shapefile: {len(gdf)} polygon(s), reprojected to EPSG:4326")
    minx, miny, maxx, maxy = (float(v) for v in gdf.total_bounds)
    west = float(np.floor(minx / RES) * RES)
    east = float(np.ceil(maxx / RES) * RES)
    south = float(np.floor(miny / RES) * RES)
    north = float(np.ceil(maxy / RES) * RES)
    mesh = None
    if os.path.exists(mesh_path):
        try:
            cached = np.load(mesh_path, allow_pickle=False)
            if (tuple(cached["bounds"]) == (west, south, east, north)
                    and tuple(cached["shp_sig"]) == tuple(shp_sig)
                    and float(cached["res"]) == RES):
                H, W = int(cached["H"]), int(cached["W"])
                mask = cached["mask"].astype(bool)
                lat, lon = cached["lat"], cached["lon"]
                mesh = (H, W)
                log(f"mesh cache HIT: {mesh_path} "
                    f"({W} x {H}, mask {100 * mask.mean():.1f}% inside)")
            else:
                log("mesh cache STALE (bounds/shapefile changed) - rebuilding")
        except Exception as e:  # noqa: BLE001 - corrupt cache must not be fatal
            log(f"mesh cache unreadable ({e}) - rebuilding")
    if mesh is None:
        win = from_bounds(west, south, east, north, transform)
        win = win.round_offsets().round_lengths().intersection(
            Window(0, 0, 7200, 3600))
        win = Window(int(win.col_off), int(win.row_off),
                     int(win.width), int(win.height))
        H, W = win.height, win.width
        wtrans = rasterio.windows.transform(win, transform)
        log(f"ROI window: {W} x {H} px, lon [{west:.2f},{east:.2f}], "
            f"lat [{south:.2f},{north:.2f}]")
        mask = rasterize(
            [(geom, 1) for geom in gdf.geometry],
            out_shape=(H, W), transform=wtrans, fill=0,
            all_touched=False, dtype="uint8").astype(bool)
        log(f"polygon mask: {mask.sum()} / {mask.size} cells inside "
            f"({100 * mask.mean():.1f}%)")
        lon = west + (np.arange(W) + 0.5) * RES
        lat = north - (np.arange(H) + 0.5) * RES  # decreasing N->S, raster order
        np.savez_compressed(
            mesh_path, bounds=np.array([west, south, east, north]),
            H=np.array(H), W=np.array(W), res=np.array(RES),
            shp_sig=np.array(shp_sig, dtype="float64"),
            mask=mask.astype("uint8"), lat=lat, lon=lon)
        log(f"mesh cache saved: {mesh_path}")
    else:
        win = Window(int(round((west + 180.0) / RES)),
                     int(round((90.0 - north) / RES)), W, H)
    globals()["MASK_ARR"] = mask  # read-only share for dask worker threads
    win_tpl = (int(win.col_off), int(win.row_off), int(win.width),
               int(win.height))

    # ---- open/create output (resume-aware) ----
    os.makedirs(OUT_DIR, exist_ok=True)
    existing_dates = set()
    if os.path.exists(OUT_PATH):
        with nc.Dataset(OUT_PATH, "a") as ds:
            for vname in ("sif", "lat", "lon"):
                if vname not in ds.variables:
                    raise SystemExit(
                        f"{OUT_PATH} exists but lacks '{vname}' - refusing to "
                        f"append to an incompatible file; move it aside first")
            if tuple(ds["sif"].dimensions) != ("time", "lat", "lon"):
                raise SystemExit("existing file has unexpected sif dimensions")
            if ds["lat"].shape != (H,) or ds["lon"].shape != (W,) \
                    or not np.allclose(ds["lat"][:], lat) \
                    or not np.allclose(ds["lon"][:], lon):
                raise SystemExit("existing file grid != current ROI grid; "
                                 "move it aside first")
            existing_dates = {d.strftime("%Y-%m-%d") for d in
                                decode_time(ds["time"])}
            log(f"resuming: {len(existing_dates)} dates already present")
    else:
        with nc.Dataset(OUT_PATH, "w", format="NETCDF4") as ds:
            ds.createDimension("time", None)
            ds.createDimension("time_monthly", None)
            ds.createDimension("lat", H)
            ds.createDimension("lon", W)
            tv = ds.createVariable("time", "f8", ("time",))
            tv.units = "days since 1970-01-01"
            tv.calendar = "standard"
            tv.standard_name = "time"
            tv.axis = "T"
            tv.long_name = "start date of 8-day composite"
            tm = ds.createVariable("time_monthly", "f8", ("time_monthly",))
            tm.units = "days since 1970-01-01"
            tm.calendar = "standard"
            tm.standard_name = "time"
            tm.axis = "T"
            tm.long_name = "first day of calendar month"
            lav = ds.createVariable("lat", "f8", ("lat",))
            lav.units = "degrees_north"
            lav.standard_name = "latitude"
            lav.axis = "Y"
            lav[:] = lat
            lov = ds.createVariable("lon", "f8", ("lon",))
            lov.units = "degrees_east"
            lov.standard_name = "longitude"
            lov.axis = "X"
            lov[:] = lon
            sv = ds.createVariable("sif", "f4", ("time", "lat", "lon"),
                                   zlib=True, complevel=COMPRESSION,
                                   chunksizes=(1, H, W),
                                   fill_value=np.float32("nan"))
            sv.long_name = ("Solar-induced chlorophyll fluorescence, "
                            "GOSIF v2 8-day composite")
            sv.units = UNITS
            sv.grid_mapping = "crs"
            sv.coordinates = "time lat lon"
            mv = ds.createVariable("sif_monthly", "f4",
                                   ("time_monthly", "lat", "lon"),
                                   zlib=True, complevel=COMPRESSION,
                                   chunksizes=(1, H, W),
                                   fill_value=np.float32("nan"))
            mv.long_name = ("Solar-induced chlorophyll fluorescence, "
                            "calendar-month mean of contributing 8-day composites")
            mv.units = UNITS
            mv.grid_mapping = "crs"
            mv.coordinates = "time_monthly lat lon"
            mv.cell_methods = "time_monthly: mean over contributing 8-day composites"
            nv = ds.createVariable("n_monthly", "u1",
                                   ("time_monthly", "lat", "lon"),
                                   zlib=True, complevel=COMPRESSION,
                                   chunksizes=(1, H, W),
                                   fill_value=np.uint8(255))
            nv.long_name = ("number of 8-day composites contributing to "
                            "sif_monthly (0 = no data that month)")
            nv.grid_mapping = "crs"
            nv.coordinates = "time_monthly lat lon"
            crs = ds.createVariable("crs", "i4", ())
            crs.grid_mapping_name = "latitude_longitude"
            crs.longitude_of_prime_meridian = 0.0
            crs.semi_major_axis = 6378137.0
            crs.inverse_flattening = 298.257223563
            crs.crs_wkt = ('GEOGCRS["WGS 84",DATUM["World Geodetic System 1984",'
                           'ELLIPSOID["WGS 84",6378137,298.257223563]],'
                           'PRIMEM["Greenwich",0],CS[ellipsoidal,2]]')
            crs.epsg_code = "EPSG:4326"
            ds.Conventions = "CF-1.8"
            ds.title = ("GOSIF v2 solar-induced fluorescence, 8-day and monthly, "
                        "clipped to the Ogallala (High Plains) aquifer")
            ds.summary = ("Per-8-day GOSIF SIF composites windowed to the "
                          "hp_bound2010 bbox on the native 0.05-degree grid, "
                          "masked to the aquifer polygon (NaN outside it), "
                          "plus calendar-month means derived from them.")
            ds.source = ("GOSIF v2 8-day global GeoTIFFs (Li, X. & Xiao, J. "
                         "2019, Remote Sensing 11, 517; FORMAT_VERSION GOSIF "
                         "v2 compatible encoding)")
            ds.references = "https://doi.org/10.3390/rs11050517"
            ds.history = (f"Created {time.strftime('%Y-%m-%dT%H:%M:%S')} by "
                          "process_gosif_ogallala.py (SSM env)")
            ds.geospatial_lat_min, ds.geospatial_lat_max = float(south), float(north)
            ds.geospatial_lon_min, ds.geospatial_lon_max = float(west), float(east)
            ds.geospatial_lat_resolution = ds.geospatial_lon_resolution = RES
            ds.time_coverage_start = str(min(d for d, _ in dated))
            ds.time_coverage_end = str(max(d for d, _ in dated))
        log(f"created {OUT_PATH}")

    todo = [(d, f) for d, f in dated if d.strftime("%Y-%m-%d") not in existing_dates]
    _maxn = globals().get("_MAX_DATES")
    if _maxn:
        todo = todo[:_maxn]
        log(f"SMOKE MODE: capped at first {_maxn} dates")
    log(f"{len(todo)} dates to ingest ({len(dated) - len(todo)} already present)")

    # ---- ingest in CHUNK_DATES-date batches (dask-parallel reads) ----
    # Per-file reads are independent, so each chunk is fetched concurrently
    # across N_WORKERS threads; the single main thread then appends results
    # to the NetCDF serially and strictly in date order (resume-safe).
    # Heartbeat logs bracket every fetch so a stall is always localizable.
    # Calm GDAL for many concurrent reads: no directory scans, no unbounded
    # VSI cache growth across 1000+ files, no nested decompression threads
    # (parallelism lives at file level via the dask workers). Entered once,
    # intentionally process-lifetime (closes on exit).
    if not _HAS_DASK:
        log("WARNING dask not installed - serial reads "
            "(pip install dask for parallel ingest)")
    rio_env = rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                           VSI_CACHE="FALSE", GDAL_NUM_THREADS="1")
    rio_env.__enter__()
    n_new = 0
    skipped = []
    n_chunks = max(1, (len(todo) + CHUNK_DATES - 1) // CHUNK_DATES)
    bar = tqdm(total=len(todo), unit="dates", desc="ingest",
               disable=(len(todo) == 0))
    for ci in range(0, len(todo), CHUNK_DATES):
        batch = todo[ci:ci + CHUNK_DATES]
        nw = N_WORKERS if _HAS_DASK else 1
        log(f"chunk {ci // CHUNK_DATES + 1}/{n_chunks}: fetching "
            f"{len(batch)} files with {nw} worker(s)...")
        t_fetch = time.time()
        if _HAS_DASK:
            tasks = [dask.delayed(read_one)(f, *win_tpl) for _, f in batch]
            results = dask.compute(
                *tasks, scheduler="threads", num_workers=N_WORKERS)
        else:
            results = [read_one(f, *win_tpl) for _, f in batch]
        dt_fetch = time.time() - t_fetch
        log(f"chunk {ci // CHUNK_DATES + 1}/{n_chunks}: fetched in "
            f"{dt_fetch:.1f}s "
            f"({len(batch) / max(dt_fetch, 1e-6):.1f} files/s)")
        arrs, tnums = [], []
        for (d, f), (ok, payload) in zip(batch, results):
            if not ok:
                log(f"WARNING {os.path.basename(f)}: {payload} - skipped")
                skipped.append(os.path.basename(f))
                continue
            arrs.append(payload)
            tnums.append((d - date(1970, 1, 1)).days)
        if not arrs:
            log(f"chunk {ci // CHUNK_DATES + 1}/{n_chunks}: no readable "
                f"dates, skipping write")
            continue
        with nc.Dataset(OUT_PATH, "a") as ds:
            t0 = ds.dimensions["time"].size
            n = len(arrs)
            ds["time"][t0:t0 + n] = np.array(tnums, dtype="float64")
            ds["sif"][t0:t0 + n, :, :] = np.stack(arrs).astype("float32")
            ds.sync()
        n_new += len(arrs)
        vf = float(np.mean([np.isfinite(a).mean() for a in arrs]))
        log(f"chunk {ci // CHUNK_DATES + 1}/{n_chunks}: wrote {len(arrs)} "
            f"dates (mean in-polygon valid frac {vf:.3f})")
        bar.update(len(arrs))
        del arrs, results
    bar.close()
    log(f"ingest done: {n_new} new dates, {len(skipped)} skipped "
        f"({(time.time() - t_start) / 60:.1f} min)")
    if skipped:
        log(f"skipped files: {skipped[:10]}"
            f"{' ...' if len(skipped) > 10 else ''}")

    # ---- monthly aggregation: stream the finished 8-day series ----
    log("deriving monthly means (streaming 8-day series back in chunks)")
    with nc.Dataset(OUT_PATH, "a") as ds:
        nt = ds.dimensions["time"].size
        if nt == 0:
            raise SystemExit(
                "TIME DIMENSION IS EMPTY (0 dates): no data were ingested - "
                "refusing to write an empty product. Check that input files "
                "were readable (see warnings above).")
        tall = decode_time(ds["time"])
        ym = np.array([(t.year, t.month) for t in tall])
        months = sorted({tuple(x) for x in ym.tolist()})
        nmon = len(months)
        # Monthly layers are rewritten fully on every run (idempotent across
        # resumes). Unlimited dims auto-extend on indexed write, so no
        # explicit resize call is needed (and none is used for portability).
        msum = np.zeros((nmon, H, W), dtype="float64")
        mcnt = np.zeros((nmon, H, W), dtype="int32")
        midx_of = {m: i for i, m in enumerate(months)}
        msteps = list(range(0, nt, CHUNK_DATES))
        for ci in tqdm(msteps, unit="chunk", desc="monthly"):
            sl = slice(ci, min(ci + CHUNK_DATES, nt))
            block = ds["sif"][sl, :, :].astype("float64")  # (<=50, H, W)
            for k, t in enumerate(tall[sl]):
                j = midx_of[(t.year, t.month)]
                v = block[k]
                ok = np.isfinite(v)
                msum[j][ok] += v[ok]
                mcnt[j][ok] += 1
            del block
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(mcnt > 0, msum / np.maximum(mcnt, 1), np.nan)
        ds["sif_monthly"][:, :, :] = mean.astype("float32")
        ds["n_monthly"][:, :, :] = np.minimum(mcnt, 255).astype("uint8")
        firsts = [date(y, m, 1) for (y, m) in months]
        ds["time_monthly"][:] = np.array(
            [(d - date(1970, 1, 1)).days for d in firsts], dtype="float64")
        ds.sync()
    log(f"monthly done: {nmon} months")

    # ---- final verification (reopen read-only) ----
    with nc.Dataset(OUT_PATH, "r") as ds:
        nt = ds.dimensions["time"].size
        nm = ds.dimensions["time_monthly"].size
        assert nt > 0, "verification failed: time is empty"
        assert ds["sif"].dimensions == ("time", "lat", "lon")
        assert ds["sif_monthly"].dimensions == ("time_monthly", "lat", "lon")
        assert ds["sif"].units == UNITS
        tall = decode_time(ds["time"])
        s = ds["sif"][0, :, :]
        v = np.asarray(s, dtype="float64")
        vf = float(np.isfinite(v).mean())
        log(f"VERIFY OK: time={nt} ({tall.min().date()}..{tall.max().date()}), "
            f"time_monthly={nm}, first-date in-polygon valid frac={vf:.3f}, "
            f"first-date valid range "
            f"[{np.nanmin(v):.4f}, {np.nanmax(v):.4f}] {UNITS}")
    log(f"DONE in {(time.time() - t_start) / 60:.1f} min -> {OUT_PATH}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="GOSIF v2 8-day -> Ogallala-clipped CF NetCDF "
                    "(50-date chunks, resume-safe)")
    ap.add_argument("--max-dates", type=int, default=None,
                    help="process only the first N dates (smoke test; "
                         "resume with a full run afterwards)")
    ap.add_argument("--out", type=str, default=None,
                    help="override output .nc path (default: derived_usgs)")
    ap.add_argument("--workers", type=int, default=None,
                    help=f"dask reader threads per chunk (default: {N_WORKERS}; "
                         f"use --workers 1 for serial reads if parallel "
                         f"reads ever stall on your drive)")
    args = ap.parse_args()
    if args.out:
        globals()["OUT_PATH"] = args.out
    if args.max_dates is not None:
        globals()["_MAX_DATES"] = args.max_dates
    if args.workers is not None:
        globals()["N_WORKERS"] = max(1, args.workers)
    main()
