"""
GOSIF v2 8-day global GeoTIFFs -> Ogallala-clipped, CF-1.8 NetCDF.

Reads every GOSIF_YYYYDDD.tif.gz in SRC_DIR through rasterio's /vsigzip/
handler (no extraction to disk), windows each scene to the hp_bound2010
bbox snapped to the native 0.05-degree grid, masks cells outside the
aquifer polygon, applies scale/offset + fill handling, and appends the
result to an unlimited time axis. Monthly means are derived afterwards by
streaming the finished 8-day series back in 50-date chunks.

Parallelism: per-scene work (gunzip -> windowed read -> mask/scale) runs
on a local dask.distributed cluster (one process per worker; scenes are
independent and per-scene results are ~0.2 MB). A one-chunk prefetch
pipeline keeps NetCDF compression overlapped with decompression. Falls
back to a ThreadPoolExecutor (or fully serial) if dask is unavailable.

Mesh reuse: window, polygon mask and lon/lat are a pure function of
(shapefile, RES, global grid) on this regular grid. They are computed
once and cached to GOSIF_ROI_mesh.npz (rebuilt only when the shapefile
mtime/size or grid signature changes), then loaded lazily exactly once
per worker and reused for every scene.

Memory design: ~2 x 50-date chunks (50 x H x W float32, ~25 MB for the
Ogallala window) are resident in the parent, plus one scene per worker.
Global 7200 x 36000 scenes are never fully loaded (except by --bench,
one scene at a time).

Resume: dates already present in OUT_PATH (matched by calendar date) are
skipped; monthly layers are always recomputed at the end of a run.
Re-running is idempotent.

Usage (SSM env: rasterio, netCDF4, geopandas, numpy, pandas, tqdm;
dask.distributed recommended):
    conda install -n SSM -c conda-forge tqdm dask distributed
    conda run -n SSM python G:\\MSU_GWB\\datasets\\process_gosif_ogallala.py \
        [--workers 8] [--bench 3] [--max-dates N] [--out path.nc]
Outputs:
    <OUT_DIR>/GOSIF_Ogallala_8day_monthly.nc
    <OUT_DIR>/GOSIF_ROI_mesh.npz            (reusable ROI mesh cache)
"""

import argparse
import collections
import glob
import os
import time
from datetime import date, timedelta
from types import SimpleNamespace

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window, from_bounds
from rasterio.features import rasterize
import geopandas as gpd
import netCDF4 as nc
from tqdm import tqdm

# --------------------------------------------------------------------------
# Config (edit here; everything else derives from these + the first file)
# --------------------------------------------------------------------------
SRC_DIR = r"G:\MSU_GWB\New folder"
SHP_PATH = r"G:\MSU_GWB\datasets\HPA_polygon\hp_bound2010.shp"
OUT_DIR = r"G:\MSU_GWB\datasets\merged_datasets"
OUT_NAME = "GOSIF_Ogallala_8day_monthly.nc"
ROI_NAME = "GOSIF_ROI_mesh.npz"
CHUNK_DATES = 50          # dates processed (and NetCDF-synced) per batch
RES = 0.05                # GOSIF native resolution, degrees (asserted per file)
GRID_W, GRID_H = 7200, 3600
SCALE = 0.0001            # per GOSIF v2 tags
FILL_CODES = (32766, 32767)   # snow/ice, water (per GOSIF v2 tags)
UNITS = "W m-2 um-1 sr-1"     # per GOSIF v2 tags
COMPRESSION = 4
# --------------------------------------------------------------------------

OUT_PATH = os.path.join(OUT_DIR, OUT_NAME)
ROI_PATH = os.path.join(OUT_DIR, ROI_NAME)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def vsi(p):
    return "/vsigzip/" + p.replace("\\", "/")


def decode_time(time_var):
    """Version-proof decode of 'days since 1970-01-01'."""
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
        doys.append(361)
    return sorted(set(doys))


# --------------------------------------------------------------------------
# ROI mesh: precompute once, cache to disk, reuse for every scene / worker
# --------------------------------------------------------------------------
_ROI = None


def _roi_signature():
    st = os.stat(SHP_PATH)
    return (f"{os.path.abspath(SHP_PATH)}|{int(st.st_mtime)}|{st.st_size}"
            f"|{RES}|{GRID_W}x{GRID_H}")


def _load_roi():
    """Read the cached mesh file. Self-contained: no loop variables, no
    external state beyond ROI_PATH."""
    with np.load(ROI_PATH, allow_pickle=False) as z:
        win = Window(int(z["win"][0]), int(z["win"][1]),
                     int(z["win"][2]), int(z["win"][3]))
        bounds = tuple(float(v) for v in z["bounds"])
        return SimpleNamespace(win=win, mask=z["mask"], lon=z["lon"],
                               lat=z["lat"], bounds=bounds)


def get_roi():
    """Per-process lazy cache: the parent and each dask worker each load
    the mesh from disk exactly once, then reuse it for every scene."""
    global _ROI
    if _ROI is None:
        _ROI = _load_roi()
    return _ROI


def ensure_roi_cache(transform):
    """Build window + polygon mask + lon/lat once and cache to disk;
    later runs and all workers are plain cache hits."""
    sig = _roi_signature()
    if os.path.exists(ROI_PATH):
        try:
            with np.load(ROI_PATH, allow_pickle=False) as z:
                if str(z["sig"][0]) == sig:
                    log(f"ROI mesh cache hit: {ROI_PATH}")
                    return
        except Exception as e:
            log(f"ROI cache unreadable ({e}) - rebuilding")
    gdf = gpd.read_file(SHP_PATH).to_crs("EPSG:4326")
    log(f"shapefile: {len(gdf)} polygon(s), reprojected to EPSG:4326")
    minx, miny, maxx, maxy = gdf.total_bounds
    west = np.floor(minx / RES) * RES
    east = np.ceil(maxx / RES) * RES
    south = np.floor(miny / RES) * RES
    north = np.ceil(maxy / RES) * RES
    win = from_bounds(west, south, east, north, transform)
    win = win.round_offsets().round_lengths().intersection(
        Window(0, 0, GRID_W, GRID_H))
    win = Window(int(win.col_off), int(win.row_off),
                 int(win.width), int(win.height))
    H, W = int(win.height), int(win.width)
    wtrans = rasterio.windows.transform(win, transform)
    mask = rasterize(
        [(geom, 1) for geom in gdf.geometry],
        out_shape=(H, W), transform=wtrans, fill=0,
        all_touched=False, dtype="uint8").astype(bool)
    lon = west + (np.arange(W) + 0.5) * RES
    lat = north - (np.arange(H) + 0.5) * RES   # decreasing N->S, raster order
    os.makedirs(OUT_DIR, exist_ok=True)
    tmp = ROI_PATH + ".building"
    with open(tmp, "wb") as fh:
        np.savez_compressed(
            fh, mask=mask, lon=lon, lat=lat,
            win=np.asarray([win.col_off, win.row_off,
                            win.width, win.height], dtype="int64"),
            bounds=np.asarray([west, south, east, north]),
            sig=np.asarray([sig]))
    os.replace(tmp, ROI_PATH)
    log(f"ROI mesh built & cached: {W} x {H} px, "
        f"{int(mask.sum())}/{mask.size} cells in polygon "
        f"({100 * mask.mean():.1f}%) -> {ROI_PATH}")# --------------------------------------------------------------------------
# Per-scene worker (runs in dask processes / pool threads / parent)
# --------------------------------------------------------------------------
def process_one(task):
    """Open one .tif.gz via /vsigzip/, window to the ROI, apply fill/scale/
    mask. Returns (date, float32 HxW | None, meta). Catches everything so
    one bad scene never kills a batch."""
    d, path = task
    try:
        roi = get_roi()
        t0 = time.perf_counter()
        with rasterio.open(vsi(path)) as src:
            if (src.width, src.height) != (GRID_W, GRID_H):
                return d, None, {"warn": f"unexpected shape "
                                         f"{(src.width, src.height)}"}
            a = src.read(1, window=roi.win, out_dtype="float32")
        read_ms = (time.perf_counter() - t0) * 1e3
        if a.shape != roi.mask.shape:
            return d, None, {"warn": f"window/mask shape mismatch "
                                     f"{a.shape} vs {roi.mask.shape}"}
        bad = (a == FILL_CODES[0]) | (a == FILL_CODES[1])
        a[bad | ~roi.mask] = np.nan      # fills + outside-polygon, one pass
        a *= np.float32(SCALE)
        return d, a, {"read_ms": read_ms,
                      "valid_frac": float(np.isfinite(a).mean())}
    except Exception as e:
        return d, None, {"warn": f"unreadable ({e})"}


def run_bench(paths, roi):
    """Efficiency verification: is a windowed /vsigzip/ read actually cheaper
    than a full read? Usually barely - gzip has no seek index, so GDAL must
    inflate the stream from byte 0 up to the last ROI row either way. That
    is precisely why worker parallelism across files is the right lever."""
    log(f"bench: windowed vs full read via /vsigzip/, single worker "
        f"({len(paths)} file(s))")
    ratios = []
    for p in paths:
        with rasterio.open(vsi(p)) as src:
            t0 = time.perf_counter()
            wa = src.read(1, window=roi.win)
            t_win = (time.perf_counter() - t0) * 1e3
        del wa
        with rasterio.open(vsi(p)) as src:
            t0 = time.perf_counter()
            fa = src.read(1)
            t_full = (time.perf_counter() - t0) * 1e3
        del fa
        ratios.append(t_win / t_full)
        log(f"  {os.path.basename(p)}: windowed {t_win:6.0f} ms | "
            f"full {t_full:6.0f} ms | ratio {t_win / t_full:.0%}")
    r = float(np.mean(ratios))
    if r > 0.25:
        log(f"mean ratio {r:.0%}: windowed reads still pay most of the gzip "
            f"inflate cost -> per-scene cost is ~O(file), single-core in "
            f"GDAL; --workers N gives ~Nx throughput (disk permitting). "
            f"Masking/scaling is <1 ms and irrelevant.")
    else:
        log(f"mean ratio {r:.0%}: window reads are cheap for this ROI; "
            f"ingest is closer to I/O-optimal already.")


# --------------------------------------------------------------------------
def main(args):
    t_start = time.time()
    out_path = args.out or OUT_PATH
    workers = args.workers
    log(f"source dir: {SRC_DIR}")
    files = sorted(glob.glob(os.path.join(SRC_DIR, "GOSIF_*.tif.gz")))
    if not files:
        raise SystemExit("no GOSIF_*.tif.gz files found - nothing to do")
    log(f"found {len(files)} input files")

    # ---- parse + audit dates ----
    dated, bad_names = [], []
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
    with rasterio.open(vsi(dated[0][1])) as src:
        assert (src.width, src.height) == (GRID_W, GRID_H), \
            f"unexpected raster size {(src.width, src.height)}"
        resx, resy = abs(src.transform.a), abs(src.transform.e)
        assert abs(resx - RES) < 1e-9 and abs(resy - RES) < 1e-9, \
            f"unexpected resolution {(resx, resy)}"
        assert src.count == 1 and src.dtypes[0] in ("int16", "uint16"), \
            f"unexpected bands/dtype {(src.count, src.dtypes)}"
        assert str(src.crs) == "EPSG:4326", f"unexpected CRS {src.crs}"
        transform = src.transform
    log("reference grid OK: global 0.05-deg EPSG:4326, 1 x int16 band")

    # ---- ROI mesh (cached; reused by parent and all workers) ----
    ensure_roi_cache(transform)
    roi = get_roi()
    H, W = roi.mask.shape
    west, south, east, north = roi.bounds
    log(f"ROI window: {W} x {H} px, lon [{west:.2f},{east:.2f}], "
        f"lat [{south:.2f},{north:.2f}]")

    # ---- optional single-worker I/O bench (verification only) ----
    if args.bench:
        run_bench([f for _, f in dated[:args.bench]], roi)
        return

    # ---- open/create output (resume-aware) ----
    os.makedirs(OUT_DIR, exist_ok=True)
    existing_dates = set()
    if os.path.exists(out_path):
        with nc.Dataset(out_path, "a") as ds:
            for vname in ("sif", "lat", "lon"):
                if vname not in ds.variables:
                    raise SystemExit(
                        f"{out_path} exists but lacks '{vname}' - refusing "
                        f"to append to an incompatible file; move it aside")
            if tuple(ds["sif"].dimensions) != ("time", "lat", "lon"):
                raise SystemExit("existing file has unexpected sif dimensions")
            if ds["lat"].shape != (H,) or ds["lon"].shape != (W,) \
                    or not np.allclose(ds["lat"][:], roi.lat) \
                    or not np.allclose(ds["lon"][:], roi.lon):
                raise SystemExit("existing file grid != current ROI grid; "
                                 "move it aside first")
            existing_dates = {d.strftime("%Y-%m-%d") for d in
                              decode_time(ds["time"])}
            log(f"resuming: {len(existing_dates)} dates already present")
    else:
        with nc.Dataset(out_path, "w", format="NETCDF4") as ds:
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
            lav[:] = roi.lat
            lov = ds.createVariable("lon", "f8", ("lon",))
            lov.units = "degrees_east"
            lov.standard_name = "longitude"
            lov.axis = "X"
            lov[:] = roi.lon
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
                          "process_gosif_ogallala.py (SSM env, dask-parallel "
                          "ingest, cached ROI mesh)")
            ds.geospatial_lat_min, ds.geospatial_lat_max = float(south), float(north)
            ds.geospatial_lon_min, ds.geospatial_lon_max = float(west), float(east)
            ds.geospatial_lat_resolution = ds.geospatial_lon_resolution = RES
            ds.time_coverage_start = str(min(d for d, _ in dated))
            ds.time_coverage_end = str(max(d for d, _ in dated))
        log(f"created {out_path}")

    todo = [(d, f) for d, f in dated
            if d.strftime("%Y-%m-%d") not in existing_dates]
    if args.max_dates:
        todo = todo[:args.max_dates]
        log(f"SMOKE MODE: capped at first {args.max_dates} dates")
    log(f"{len(todo)} dates to ingest "
        f"({len(dated) - len(todo)} already present)")

    # ---- ingest: parallel per-scene workers + tqdm, chunked NetCDF writes ----
    n_new = 0
    stats = {"read_ms": [], "vf": [], "warn": 0}
    if todo:
        t_ing = time.time()
        bar = tqdm(total=len(todo), desc="ingest 8-day", unit="file",
                   dynamic_ncols=True)

        def write_chunk(results):
            nonlocal n_new
            arrs, tnums = [], []
            for d, a, meta in results:
                if a is None:
                    stats["warn"] += 1
                    tqdm.write(f"[{time.strftime('%H:%M:%S')}] WARNING {d}: "
                               f"{meta['warn']} - scene skipped")
                    continue
                arrs.append(a)
                tnums.append((d - date(1970, 1, 1)).days)
                stats["read_ms"].append(meta["read_ms"])
                stats["vf"].append(meta["valid_frac"])
            if arrs:
                with nc.Dataset(out_path, "a") as ds:
                    t0 = ds.dimensions["time"].size
                    ds["time"][t0:t0 + len(arrs)] = \
                        np.asarray(tnums, dtype="float64")
                    ds["sif"][t0:t0 + len(arrs)] = np.stack(arrs)
                    ds.sync()
                n_new += len(arrs)
            bar.update(len(results))
            if stats["read_ms"]:
                bar.set_postfix(ok=n_new, skip=stats["warn"],
                                ms=f"{np.mean(stats['read_ms'][-200:]):.0f}",
                                vf=f"{np.mean(stats['vf'][-200:]):.3f}")

        client = cluster = None
        if workers > 1:
            try:
                from dask.distributed import Client, LocalCluster
                cluster = LocalCluster(n_workers=workers,
                                       threads_per_worker=1,
                                       processes=True, silence_logs=True)
                client = Client(cluster)
                client.run(get_roi)          # warm mesh cache once per worker
                log(f"dask LocalCluster: {workers} workers | "
                    f"dashboard: {client.dashboard_link}")
            except Exception as e:
                if client is not None:
                    try:
                        client.close()
                        cluster.close()
                    except Exception:
                        pass
                client = cluster = None
                log(f"WARNING dask backend unavailable "
                    f"({type(e).__name__}: {e}) - thread-pool fallback")

        try:
            if client is not None:
                # Rolling pipeline: keep <=2 chunks in flight so NetCDF
                # compression of chunk i overlaps gunzip of chunk i+1.
                # gather() preserves order -> time axis stays ascending.
                queued, nxt = collections.deque(), 0

                def submit_more():
                    nonlocal nxt
                    while nxt < len(todo) and len(queued) < 2:
                        queued.append(client.map(
                            process_one, todo[nxt:nxt + CHUNK_DATES]))
                        nxt += CHUNK_DATES

                submit_more()
                while queued:
                    write_chunk(client.gather(queued.popleft()))
                    submit_more()
            elif workers > 1:
                from concurrent.futures import ThreadPoolExecutor
                with ThreadPoolExecutor(max_workers=workers) as ex:
                    for ci in range(0, len(todo), CHUNK_DATES):
                        write_chunk(list(ex.map(
                            process_one, todo[ci:ci + CHUNK_DATES])))
            else:
                for ci in range(0, len(todo), CHUNK_DATES):
                    write_chunk([process_one(t)
                                 for t in todo[ci:ci + CHUNK_DATES]])
        finally:
            bar.close()
            if client is not None:
                client.close()
                cluster.close()

        ms = np.asarray(stats["read_ms"])
        wall = time.time() - t_ing
        if ms.size:
            log(f"ingest perf: {ms.size} scenes, mean {ms.mean():.0f} ms/scene, "
                f"p95 {np.percentile(ms, 95):.0f} ms, "
                f"{ms.size / wall:.1f} scenes/s end-to-end, "
                f"mean valid frac {np.mean(stats['vf']):.3f}")
    else:
        log("no new dates to ingest")
    log(f"ingest done: {n_new} new dates, {stats['warn']} skipped "
        f"({(time.time() - t_start) / 60:.1f} min)")

    # ---- monthly aggregation: stream the finished 8-day series ----
    log("deriving monthly means (streaming 8-day series back in chunks)")
    with nc.Dataset(out_path, "a") as ds:
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
        msum = np.zeros((nmon, H, W), dtype="float64")
        mcnt = np.zeros((nmon, H, W), dtype="int32")
        midx_of = {m: i for i, m in enumerate(months)}
        with tqdm(total=nt, desc="monthly means", unit="file",
                  dynamic_ncols=True) as mbar:
            for ci in range(0, nt, CHUNK_DATES):
                sl = slice(ci, min(ci + CHUNK_DATES, nt))
                block = np.ma.filled(ds["sif"][sl, :, :], np.nan
                                     ).astype("float32")
                okb = np.isfinite(block)
                for k in range(block.shape[0]):
                    t = tall[ci + k]
                    j = midx_of[(t.year, t.month)]
                    m = okb[k]
                    msum[j][m] += block[k][m]
                    mcnt[j][m] += 1
                mbar.update(block.shape[0])
                del block, okb
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
    with nc.Dataset(out_path, "r") as ds:
        nt = ds.dimensions["time"].size
        nm = ds.dimensions["time_monthly"].size
        assert nt > 0, "verification failed: time is empty"
        assert ds["sif"].dimensions == ("time", "lat", "lon")
        assert ds["sif_monthly"].dimensions == ("time_monthly", "lat", "lon")
        assert ds["sif"].units == UNITS
        tall = decode_time(ds["time"])
        v = np.asarray(np.ma.filled(ds["sif"][0, :, :], np.nan),
                       dtype="float64")
        vf = float(np.isfinite(v).mean())
        log(f"VERIFY OK: time={nt} ({tall.min().date()}..{tall.max().date()}), "
            f"time_monthly={nm}, first-date in-polygon valid frac={vf:.3f}, "
            f"first-date valid range "
            f"[{np.nanmin(v):.4f}, {np.nanmax(v):.4f}] {UNITS}")
    log(f"DONE in {(time.time() - t_start) / 60:.1f} min -> {out_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="GOSIF v2 8-day -> Ogallala-clipped CF NetCDF "
                    "(dask-parallel, tqdm, cached ROI mesh, resume-safe)")
    ap.add_argument("--max-dates", type=int, default=None,
                    help="process only the first N dates (smoke test)")
    ap.add_argument("--out", type=str, default=None,
                    help="override output .nc path")
    ap.add_argument("--workers", type=int,
                    default=max(1, min(os.cpu_count() or 4, 8)),
                    help="parallel scene readers (default min(8, cores)); "
                         "~4 on a spinning disk, more on SSD/NVMe")
    ap.add_argument("--bench", type=int, default=0, metavar="N",
                    help="benchmark windowed vs full /vsigzip/ reads on the "
                         "first N files, then exit (verification only)")
    main(ap.parse_args())