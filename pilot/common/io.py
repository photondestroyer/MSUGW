"""Readers for every pilot input. Paths come from config.src(key); nothing here writes to a source.

All NetCDF files are opened read-only. gridMET and SMAP are read with netCDF4's automatic mask and scale
(their scale_factor / _FillValue attributes are standard); OpenET and CDL tiles are read raw and decoded
explicitly because their fills differ by variable.
"""
import os
import glob

import numpy as np

from .. import config as C
from ..core.errors import SourceError
from . import grids


def path(key, must_exist=True):
    p = C.src(key)
    if must_exist and not os.path.exists(p):
        raise SourceError(f"{key}: not found at {p}")
    return p


def exists(key):
    return os.path.exists(C.src(key))


# ------------------------------------------------------------------------------------------- vectors and tables
def vector(key):
    import geopandas as gpd
    p = path(key)
    if p.endswith(".zip"):
        return gpd.read_file(f"zip://{p}")
    return gpd.read_file(p)


def table(key, **kw):
    import pandas as pd
    p = path(key)
    if p.endswith(".txt") or p.endswith(".csv"):
        kw.setdefault("low_memory", False)
        kw.setdefault("encoding", "latin-1")
        return pd.read_csv(p, **kw)
    raise ValueError(f"{key}: not a table ({p})")


def counties():
    """Census 2023 1:500k counties of all states (GEOID, STATEFP, NAME, geometry)."""
    g = vector("counties_zip")
    return g[["GEOID", "STATEFP", "NAME", "geometry"]].copy()


# ------------------------------------------------------------------------------------------- gridMET lab
def gridmet():
    import netCDF4 as nc
    ds = nc.Dataset(path("gridmet"), "r")
    ds.set_auto_maskandscale(True)
    return ds


def gridmet_axes(ds):
    import pandas as pd
    import netCDF4 as nc
    t = ds["time"]
    dates = pd.DatetimeIndex(nc.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"),
                                         only_use_cftime_datetimes=False))
    return np.asarray(ds["lat"][:], float), np.asarray(ds["lon"][:], float), dates


def gridmet_read(ds, var, i0, i1, r0, r1, c0, c1):
    """Days i0:i1 of var over the lat/lon box rows r0:r1, cols c0:c1 -> float array, NaN at fill."""
    a = ds[var][i0:i1, r0:r1, c0:c1]
    return np.ma.filled(a.astype("float64"), np.nan)


# ------------------------------------------------------------------------------------------- SMAP
def smap_source():
    """(path, label) of the SMAP L3E file to use: the unmasked Kansas extraction, else a fallback."""
    for key in C.SMAP_ORDER:
        p = C.src(key)
        if os.path.isdir(p):                       # the pilot folder: newest dated file, checkpoints ignored
            files = sorted(f for f in glob.glob(f"{p}/SPL3SMP_E_Kansas_pilot_*.nc") if "_ckpt_" not in f)
            if not files:
                continue
            p = files[-1]
        if os.path.exists(p):
            return p, ("pilot" if key == "smap_pilot" else f"fallback:{key}")
    raise SourceError("no SMAP L3E file found (run Colab notebooks/SMAP_L3E_Kansas_pilot.ipynb)")


def smap_open(p):
    import netCDF4 as nc
    ds = nc.Dataset(p, "r")
    ds.set_auto_maskandscale(True)
    return ds


def smap_dates(ds):
    import pandas as pd
    import netCDF4 as nc
    t = ds["time"]
    return pd.DatetimeIndex(nc.num2date(t[:], t.units, only_use_cftime_datetimes=False)).normalize()


# ------------------------------------------------------------------------------------------- irrigation, soils, DEM
def aim_path(year):
    p = f"{path('aim_dir')}/{year}_AIM-HPA_rse_finalMapBinary_20190909.tif"
    if not os.path.exists(p):
        raise SourceError(f"AIM-HPA {year}: not found at {p}")
    return p


AIM_IRRIGATED, AIM_NOT_IRRIGATED = 1, 2          # 0 = outside the study area (the file's nodata tag says 255)


def lanid():
    import netCDF4 as nc
    ds = nc.Dataset(path("lanid"), "r")
    ds.set_auto_maskandscale(False)
    return ds


def lanid_transform(ds):
    from rasterio.transform import Affine
    gt = [float(v) for v in str(ds["crs"].GeoTransform).split()]
    return Affine(gt[1], gt[2], gt[0], gt[4], gt[5], gt[3])


def lanid_1820_paths():
    """{year: path} of the LANID 2018-2020 rasters as released (lanid2018.tif ...); {} when not on disk."""
    d = C.src("lanid_1820_dir")
    out = {}
    for p in sorted(glob.glob(f"{d}/lanid[0-9][0-9][0-9][0-9].tif")) if os.path.isdir(d) else []:
        out[int(os.path.basename(p)[5:9])] = p.replace(os.sep, "/")
    return out


def gssurgo():
    import netCDF4 as nc
    ds = nc.Dataset(path("gssurgo"), "r")
    ds.set_auto_maskandscale(True)
    return ds


def gssurgo_transform(ds):
    from rasterio.transform import Affine
    txt = str(ds["crs"].comment)                   # "(120.0, 0.0, -807735.0, 0.0, -120.0, 2314755.0)": Affine order
    a = [float(v) for v in txt[txt.index("(") + 1:txt.index(")")].split(",")]
    tr = Affine(*a)
    x = np.asarray(ds["x"][:2], float)             # cross-check with the file's own pixel-centre axis
    if abs(x[0] - (tr.c + tr.a / 2)) > 1e-3:
        raise SourceError(f"gSSURGO transform {a} does not match its x axis {x[:2]}")
    return tr


def sample_raster(key_or_path, x, y, crs_xy):
    """Values of band 1 at points (x, y) given in crs_xy; NaN at nodata / outside."""
    import rasterio
    from pyproj import Transformer
    p = path(key_or_path) if key_or_path in C.DATASETS else key_or_path
    with rasterio.open(p) as r:
        tx, ty = Transformer.from_crs(crs_xy, r.crs, always_xy=True).transform(np.asarray(x), np.asarray(y))
        vals = np.array([v[0] for v in r.sample(zip(tx, ty))], float)
        nd = r.nodata
        if nd is not None:
            vals[np.isclose(vals, nd)] = np.nan
        h, w = r.height, r.width
        from rasterio.transform import rowcol
        rr, cc = rowcol(r.transform, np.asarray(tx), np.asarray(ty))
        outside = (np.asarray(rr) < 0) | (np.asarray(rr) >= h) | (np.asarray(cc) < 0) | (np.asarray(cc) >= w)
        vals[outside] = np.nan
    return vals


def well_dates(t_days_since_1900):
    """(year, month, known) of well-file times (days since 1900-01-01). Missing dates and impossible ones (the
    data audit found a record dated 3785 and others in 1899) are 'not known' (year 0, month 0)."""
    import pandas as pd
    t = np.asarray(t_days_since_1900, float)
    known = np.isfinite(t) & (t >= 0) & (t <= 60000)            # 1900-01-01 .. 2064
    d = pd.to_datetime("1900-01-01") + pd.to_timedelta(np.where(known, t, 0), unit="D")
    return np.where(known, d.year, 0), np.where(known, d.month, 0), known


# ------------------------------------------------------------------------------------------- 30 m tiles
def cdl_file(tile, year):
    return grids.tile_file(path("cdl_tiles"), C.CDL_TAG, tile, year)


def openet_file(tile, year):
    folder = path("openet_tiles")
    tag = os.path.basename(folder.rstrip("/")).replace("_KS_pilot", "")    # file names keep the release tag
    return grids.tile_file(folder, tag, tile, year)


def tile_frames(kind):
    return grids.load_frames(path("cdl_tiles" if kind == "cdl" else "openet_tiles"))


def open_tile(p):
    import netCDF4 as nc
    ds = nc.Dataset(p, "r")
    ds.set_auto_maskandscale(False)
    if str(getattr(ds, "processing_status", "complete")) != "complete":
        ds.close()
        raise SourceError(f"{p}: processing_status is not 'complete'")
    return ds


def tile_frame_of(ds, kind_frames, tile):
    """The frame of a tile file: from frames.json, checked against the file's own x/y axes."""
    F = kind_frames[tile]
    x, y = ds["x"][:], ds["y"][:]
    if abs(float(x[0]) - (F["x0"] + C.PIXEL_M / 2)) > 1e-3 or abs(float(y[0]) - (F["y0"] - C.PIXEL_M / 2)) > 1e-3 \
            or len(x) != F["w"] or len(y) != F["h"]:
        raise SourceError(f"{tile}: file axes do not match frames.json")
    return F


def openet_months(ds):
    from ..core.ncio import months_of
    return months_of(ds)
