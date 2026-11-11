"""Source identity checks: is each input the file the pipeline was written for?

Used by `python -m paper1 preflight` and by the provenance tests. Each check returns a list of problems
(empty = fine). Checks read headers and coordinates only (never whole variables).
"""
import json
import os

import numpy as np

from . import config as C
from .core import paths, ncio

LANDSAT_SCALES = dict(NIRv=1e-4, NDMI=1e-4, NDVI=1e-4, MNDWI=1e-4, AWEIsh=2e-4)


def _nc_vars(path, want, group=None):
    probs = []
    with ncio.open_ro(path) as ds:
        g = ds.groups[group] if group else ds
        miss = [v for v in want if v not in g.variables]
        if miss:
            probs.append(f"{os.path.basename(path)}: variables missing {miss}")
    return probs


def check_gridmet(p):
    probs = _nc_vars(p, C.DATASETS["gridmet"]["vars"])
    with ncio.open_ro(p) as ds:
        if getattr(ds, "processing_status", "") != "complete":
            probs.append("gridmet: processing_status is not 'complete'")
        lat, lon = ds["lat"][:], ds["lon"][:]
        if not (lat[0] > lat[-1] and lon[0] < lon[-1]):
            probs.append("gridmet: expected lat descending, lon ascending")
        if abs(abs(lat[1] - lat[0]) - 1 / 24) > 1e-6:
            probs.append("gridmet: not the 1/24 degree lab grid")
    return probs


def check_hres(p):
    probs = _nc_vars(p, ["b1", "source_image"])
    with ncio.open_ro(p) as ds:
        yn = "lat" if "lat" in ds.variables else "y"
        xn = "lon" if "lon" in ds.variables else "x"
        dy = abs(float(ds[yn][1] - ds[yn][0]))
        dx = abs(float(ds[xn][1] - ds[xn][0]))
        if dx / dy > 1.2:
            probs.append(f"hres_wtd: x/y step ratio {dx / dy:.3f} - this is the stretched-longitude export, refused")
    return probs


def check_dem(p):
    import rasterio
    probs = []
    with rasterio.open(p) as r:
        if r.crs is None or r.crs.to_epsg() != 4326:
            probs.append(f"dem: CRS {r.crs} is not EPSG:4326")
        if abs(r.res[0] * 3600 - 1) > 1e-6 or abs(r.res[1] * 3600 - 1) > 1e-6:
            probs.append(f"dem: cell {r.res} is not 1 arc-second")
        if r.nodata is None:
            probs.append("dem: no nodata declared (the old QGIS file had a zero hole; use the rebuilt file)")
    return probs


def check_lanid(p):
    probs = _nc_vars(p, ["irrigated", "roi_mask"]) + _nc_vars(p, ["maxIrrExt"], group="derived")
    with ncio.open_ro(p) as ds:
        if getattr(ds, "processing_status", "") != "complete":
            probs.append("lanid: processing_status is not 'complete'")
        x = ds["x"][:]
        if abs(abs(float(x[1] - x[0])) - 30) > 1e-6:
            probs.append("lanid: not a 30 m grid")
    return probs


def check_gssurgo(p):
    probs = _nc_vars(p, C.DATASETS["gssurgo"]["vars"])
    with ncio.open_ro(p) as ds:
        if abs(abs(float(ds["x"][1] - ds["x"][0])) - 120) > 1e-6:
            probs.append("gssurgo: not the 120 m stack")
    return probs


def check_gde(p):
    probs = _nc_vars(p, C.DATASETS["gde"]["vars"])
    with ncio.open_ro(p) as ds:
        if abs(abs(float(ds["x"][1] - ds["x"][0])) * 120 - 1) > 1e-6:
            probs.append("gde: not the 30 arc-second grid")
    return probs


FILE_CHECKS = dict(gridmet=check_gridmet, hres_wtd=check_hres, dem=check_dem, lanid=check_lanid,
                   gssurgo=check_gssurgo, gde=check_gde)


def check_tile_file(key, p, frame=None):
    """Landsat / OpenET / CDL / NWI tile file: variables, types, scale factors, status and frame."""
    probs = []
    with ncio.open_ro(p) as ds:
        miss = [v for v in C.DATASETS[key]["vars"] if v not in ds.variables]
        if miss:
            return [f"{os.path.basename(p)}: variables missing {miss}"]
        st = getattr(ds, "processing_status", "")
        if st not in ("complete", "partial"):
            probs.append(f"{os.path.basename(p)}: processing_status {st!r}")
        if key == "landsat":
            for v, s in LANDSAT_SCALES.items():
                var = ds[v]
                if var.dtype != np.int16 or abs(float(getattr(var, "scale_factor", 0)) - s) > 1e-9:
                    probs.append(f"{os.path.basename(p)}: {v} is {var.dtype} scale {getattr(var, 'scale_factor', None)}, "
                                 f"expected int16 scale {s}")
            if ds["n_obs"].dtype != np.uint8:
                probs.append(f"{os.path.basename(p)}: n_obs is not uint8")
        if key == "openet":
            var = ds["et_ensemble_mad"]
            if var.dtype != np.int16 or float(ncio.fill_of(var, 0)) != C.LANDSAT_FILL:
                probs.append(f"{os.path.basename(p)}: et_ensemble_mad is {var.dtype} fill {ncio.fill_of(var)}")
        if key == "cdl" and ds["cropland"].dtype != np.uint8:
            probs.append(f"{os.path.basename(p)}: cropland is not uint8")
        if frame is not None and "frame" in ds.ncattrs():
            F = json.loads(ds.getncattr("frame"))
            if any(abs(float(F[c]) - float(frame[c])) > 1e-6 for c in ("x0", "y0", "w", "h")):
                probs.append(f"{os.path.basename(p)}: frame differs from frames.json")
    return probs


def inventory(tiles, years):
    """Which tile files exist for the run, with sizes."""
    rows = []
    for key in ("landsat", "openet", "cdl", "nwi"):
        spec = C.DATASETS[key]
        for t in tiles:
            for y in (years if spec["per_year"] else [None]):
                if key == "cdl" and y is not None and y < C.YEAR_FIRST:
                    continue
                p = paths.tile_path(key, t, y)
                rows.append(dict(dataset=key, tile=t, year=y, path=p, exists=os.path.exists(p),
                                 size=os.path.getsize(p) if os.path.exists(p) else None))
    return rows
