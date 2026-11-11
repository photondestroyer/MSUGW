"""A small synthetic data root in the EXACT formats of the real inputs, for end-to-end tests.

One tile (z14e02n29, the real tile geometry) with a 160 x 160 pixel frame (10 x 10 blocks) in south-west
Nebraska; two years (2012, 2013) of Landsat / OpenET / CDL tile files; an NWI tile; and around it a lab
gridMET extract (5 x 6 cells, 2004-2013 daily, packed int16 with _qc), HRES-WTD, GLO-30 DEM (GeoTIFF,
1 arc-second, centres on whole arc-seconds), LANID (EPSG:5070 30 m with the derived group), gSSURGO
(EPSG:5070 120 m), GDE (30 arc-seconds), the aquifer polygon and two counties. A valley runs north-south
through the frame (shallow water table, low HAND, a pond, irrigated corn on one side), so the P1.1 rules
have something to exclude. Plan v3 adds: Landsat n_obs_oli / obs_day, OpenET model min / max / count and a
model clear-scene count, the lab drought file (spi90d every 5 days; 2012 summer dry, 2013 summer wet), LANID
derived/irrFreq (a patch irrigated in only 2 years) and an east-west canal across the north of the frame.
"""
import json
import os
import zipfile

import numpy as np

TILE = "z14e02n29"
ZONE, LEFT, TOP = 14, 15.0 + 2 * 153600.0, 15.0 + 30 * 153600.0
R0 = C0 = 1600
W = H = 160
X0, Y0 = LEFT + C0 * 30.0, TOP - R0 * 30.0
YEARS = [2012, 2013]
SEED = 3
VERSION = 2                      # bump when the synthetic inputs change: a cached root of another version is rebuilt
LANDSAT_TAG = "LANDSAT_C2L2_indices_30m_monthly_tiles_v2"
OPENET_TAG = "OPENET_v2_1_monthly_30m_tiles"
CDL_TAG = "USDA_NASS_CDL_30m_tiles"
NWI_TAG = "NWI_wetlands_30m_tiles"


def _nc(path):
    import netCDF4 as nc
    os.makedirs(os.path.dirname(path), exist_ok=True)
    ds = nc.Dataset(path, "w", format="NETCDF4")
    ds.set_auto_maskandscale(False)
    return ds


def pixel_grid():
    xs = X0 + (np.arange(W) + 0.5) * 30
    ys = Y0 - (np.arange(H) + 0.5) * 30
    return xs, ys


def lonlat_bounds(margin=0.02):
    from pyproj import Transformer
    tr = Transformer.from_crs(32614, 4326, always_xy=True)
    xs = [X0, X0 + W * 30, X0, X0 + W * 30]
    ys = [Y0, Y0, Y0 - H * 30, Y0 - H * 30]
    lo, la = tr.transform(xs, ys)
    return min(lo) - margin, min(la) - margin, max(lo) + margin, max(la) + margin


def valley_mask_utm(X, Y):
    """A north-south valley 6 pixels wide at frame column ~60."""
    c = (X - X0) / 30
    return (c > 57) & (c < 63)


def valley_of_lonlat(lon, lat):
    from pyproj import Transformer
    X, Y = Transformer.from_crs(4326, 32614, always_xy=True).transform(lon, lat)
    return valley_mask_utm(X, Y), X, Y


def frames_json(root, tag):
    d = f"{root}/{tag}"
    os.makedirs(d, exist_ok=True)
    F = {TILE: dict(crs="EPSG:32614", x0=X0, y0=Y0, w=W, h=H, n_roi=W * H - 300), "z14e02n28": None}
    with open(f"{d}/frames.json", "w") as f:
        json.dump(dict(made_with=dict(pixel=30.0, tile_px=5120), created="synthetic", frames=F), f)


def roi_mask():
    m = np.ones((H, W), "int8")
    m[:10, :30] = 0                       # a corner outside the aquifer (300 pixels)
    return m


def tile_common(ds, year=None):
    ds.createDimension("y", H)
    ds.createDimension("x", W)
    xs, ys = pixel_grid()
    ds.createVariable("x", "f8", ("x",))[:] = xs
    ds.createVariable("y", "f8", ("y",))[:] = ys
    cr = ds.createVariable("crs", "i4")
    cr.spatial_ref = "EPSG:32614"
    rv = ds.createVariable("roi_mask", "i1", ("y", "x"))
    rv[:] = roi_mask()
    ds.frame = json.dumps(dict(crs="EPSG:32614", x0=X0, y0=Y0, w=W, h=H, n_roi=W * H - 300))
    ds.processing_status = "complete"
    ds.tile = TILE
    if year is not None:
        ds.createDimension("time", None)
        tv = ds.createVariable("time", "f8", ("time",))
        tv.units, tv.calendar = "days since 1970-01-01 00:00:00", "proleptic_gregorian"
        tv[:] = [(np.datetime64(f"{year}-{m:02d}-01") - np.datetime64("1970-01-01")).astype(int) for m in range(1, 13)]


def season(m):
    return np.exp(-0.5 * ((m - 7) / 1.8) ** 2)


def landsat(root, year, rng):
    p = f"{root}/{LANDSAT_TAG}/{year}/{TILE}/{LANDSAT_TAG}_{TILE}_{year}_Ogallala.nc"
    ds = _nc(p)
    tile_common(ds, year)
    xs, ys = pixel_grid()
    X, Y = np.meshgrid(xs, ys)
    val = valley_mask_utm(X, Y)
    pond = val & ((Y0 - Y) / 30 > 100) & ((Y0 - Y) / 30 < 106)
    corn = ((X - X0) / 30 > 70) & ((X - X0) / 30 < 110) & ((Y0 - Y) / 30 > 40) & ((Y0 - Y) / 30 < 80)
    scales = dict(NIRv=1e-4, NDMI=1e-4, NDVI=1e-4, MNDWI=1e-4, AWEIsh=2e-4)
    vars_ = {}
    for v, s in scales.items():
        x = ds.createVariable(v, "i2", ("time", "y", "x"), fill_value=-32768)
        x.set_auto_maskandscale(False)                     # values below are the stored integers
        x.scale_factor = np.float32(s)
        x.units = "1"
        vars_[v] = x
    nob = ds.createVariable("n_obs", "u1", ("time", "y", "x"), fill_value=255)
    oli = ds.createVariable("n_obs_oli", "u1", ("time", "y", "x"), fill_value=255)
    oday = ds.createVariable("obs_day", "u1", ("time", "y", "x"), fill_value=255)
    roi = roi_mask() == 1
    for m in range(12):
        sz = season(m + 1)
        rain_boost = 1.0 if year == 2013 else 0.7          # 2012 was a drought year
        nirv = 0.05 + 0.20 * sz * rain_boost + 0.15 * sz * (val | corn) + rng.normal(0, 0.01, (H, W))
        ndmi = -0.10 + 0.25 * sz * rain_boost + 0.20 * sz * (val | corn) + rng.normal(0, 0.01, (H, W))
        ndvi = 0.20 + 0.45 * sz * rain_boost + 0.20 * sz * (val | corn) + rng.normal(0, 0.01, (H, W))
        mndwi = np.where(pond, 0.3, -0.4) + rng.normal(0, 0.02, (H, W))
        awei = np.where(pond, 0.5, -1.2) + rng.normal(0, 0.05, (H, W))
        n = rng.integers(2, 6, (H, W)).astype("u1")
        n[rng.random((H, W)) < 0.1] = 0                     # cloudy pixel-months: no observation
        n[~roi] = 255
        obs = (n > 0) & (n != 255)
        for v, a in dict(NIRv=nirv, NDMI=ndmi, NDVI=ndvi, MNDWI=mndwi, AWEIsh=awei).items():
            raw = np.round(a / scales[v]).astype("i2")
            raw[~obs] = -32768
            vars_[v][m] = raw
        nob[m] = n
        o = np.where(obs, (n // 2) if (year > 2013 or (year == 2013 and m >= 3)) else 0, 255).astype("u1")
        oli[m] = o
        oday[m] = np.where(obs, rng.integers(3, 27, (H, W)), 255).astype("u1")
    ds.close()
    return p


def openet(root, year, rng):
    p = f"{root}/{OPENET_TAG}/{year}/{TILE}/{OPENET_TAG}_{TILE}_{year}_Ogallala.nc"
    ds = _nc(p)
    tile_common(ds, year)
    xs, ys = pixel_grid()
    X, Y = np.meshgrid(xs, ys)
    val = valley_mask_utm(X, Y)
    corn = ((X - X0) / 30 > 70) & ((X - X0) / 30 < 110) & ((Y0 - Y) / 30 > 40) & ((Y0 - Y) / 30 < 80)
    v = ds.createVariable("et_ensemble_mad", "i2", ("time", "y", "x"), fill_value=-32768)
    v.units = "mm"
    vmin = ds.createVariable("et_ensemble_mad_min", "i2", ("time", "y", "x"), fill_value=-32768)
    vmax = ds.createVariable("et_ensemble_mad_max", "i2", ("time", "y", "x"), fill_value=-32768)
    vcnt = ds.createVariable("et_ensemble_mad_count", "u1", ("time", "y", "x"), fill_value=255)
    vsc = ds.createVariable("count_ssebop", "u1", ("time", "y", "x"), fill_value=255)
    roi = roi_mask() == 1
    for m in range(12):
        sz = season(m + 1)
        et = 15 + 80 * sz * (0.7 if year == 2012 else 1.0) + 70 * sz * (val | corn) + rng.normal(0, 3, (H, W))
        raw = np.round(et).astype("i2")
        raw[~roi] = -32768
        raw[rng.random((H, W)) < 0.02] = -32768
        v[m] = raw
        bad = raw == -32768
        vmin[m] = np.where(bad, -32768, raw - rng.integers(3, 15, (H, W))).astype("i2")
        vmax[m] = np.where(bad, -32768, raw + rng.integers(3, 15, (H, W))).astype("i2")
        vcnt[m] = np.where(bad, 255, rng.integers(3, 7, (H, W))).astype("u1")
        vsc[m] = np.where(bad, 255, rng.integers(1, 5, (H, W))).astype("u1")
    ds.close()
    return p


def cdl(root, year):
    p = f"{root}/{CDL_TAG}/{year}/{TILE}/{CDL_TAG}_{TILE}_{year}_Ogallala.nc"
    ds = _nc(p)
    tile_common(ds)
    xs, ys = pixel_grid()
    X, Y = np.meshgrid(xs, ys)
    c = np.full((H, W), 176, "u1")                         # grassland/pasture
    c[:, 112:] = 24                                         # winter wheat on the east side (dryland)
    corn = ((X - X0) / 30 > 70) & ((X - X0) / 30 < 110) & ((Y0 - Y) / 30 > 40) & ((Y0 - Y) / 30 < 80)
    c[corn] = 1
    c[roi_mask() == 0] = 255
    v = ds.createVariable("cropland", "u1", ("y", "x"), fill_value=255)
    v[:] = c
    ds.close()
    return p


def nwi(root):
    p = f"{root}/{NWI_TAG}/{TILE}/{NWI_TAG}_{TILE}_Ogallala.nc"
    ds = _nc(p)
    tile_common(ds)
    xs, ys = pixel_grid()
    X, Y = np.meshgrid(xs, ys)
    pond = valley_mask_utm(X, Y) & ((Y0 - Y) / 30 > 100) & ((Y0 - Y) / 30 < 106)
    roi = roi_mask() == 1
    sysv = np.where(pond, 5, 0).astype("u1")
    rip = valley_mask_utm(X, Y).astype("u1")
    for name, a in (("wetland_system", sysv), ("riparian", rip), ("nwi_mapped", np.ones((H, W), "u1"))):
        a = a.copy()
        a[~roi] = 255
        ds.createVariable(name, "u1", ("y", "x"), fill_value=255)[:] = a
    ds.close()
    return p


def gridmet(root, rng):
    p = f"{root}/GRIDMET_LAB/GRIDMET_Ogallala_lab.nc"
    lo0, la0, lo1, la1 = lonlat_bounds(0.1)
    lat_all = 43.65 - np.arange(287) / 24
    lon_all = -105.933333333 + np.arange(233) / 24
    li = np.nonzero((lat_all <= la1) & (lat_all >= la0))[0]
    lj = np.nonzero((lon_all >= lo0) & (lon_all <= lo1))[0]
    lat, lon = lat_all[li[0]:li[-1] + 1], lon_all[lj[0]:lj[-1] + 1]
    days = np.arange(np.datetime64("2004-01-01"), np.datetime64("2014-01-01"))
    ds = _nc(p)
    ds.createDimension("time", days.size)
    ds.createDimension("lat", lat.size)
    ds.createDimension("lon", lon.size)
    tv = ds.createVariable("time", "f8", ("time",))
    tv.units, tv.calendar = "days since 1900-01-01", "standard"
    tv[:] = (days - np.datetime64("1900-01-01")).astype(int)
    ds.createVariable("lat", "f8", ("lat",))[:] = lat
    ds.createVariable("lon", "f8", ("lon",))[:] = lon
    ds.createVariable("roi_fraction", "f4", ("lat", "lon"))[:] = np.ones((lat.size, lon.size), "f4")
    doy = (days - days.astype("datetime64[Y]")).astype(int)
    shp = (days.size, lat.size, lon.size)
    sz = np.exp(-0.5 * ((doy - 190) / 45.0) ** 2)[:, None, None]
    yr = days.astype("datetime64[Y]").astype(int) + 1970
    dry = np.where(yr == 2012, 0.4, 1.0)[:, None, None]
    pr = np.where(rng.random(shp) < 0.25, rng.gamma(0.8, 6.0, shp), 0.0) * (0.5 + sz) * dry
    vals = dict(pr=pr, pet=1 + 6 * sz + rng.normal(0, 0.3, shp), vpd=0.3 + 1.5 * sz + rng.normal(0, 0.05, shp),
                srad=120 + 200 * sz + rng.normal(0, 10, shp), tmmn=263 + 25 * sz + rng.normal(0, 1, shp),
                tmmx=275 + 25 * sz + rng.normal(0, 1, shp))
    packing = dict(pr=(0.1, 0), pet=(0.1, 0), vpd=(0.01, 0), srad=(0.1, 0), tmmn=(0.1, 210), tmmx=(0.1, 220))
    for name, a in vals.items():
        s, o = packing[name]
        v = ds.createVariable(name, "i2", ("time", "lat", "lon"), fill_value=32767)
        v.set_auto_maskandscale(False)
        v.scale_factor, v.add_offset = np.float32(s), np.float32(o)
        raw = np.round((np.maximum(a, 0 if name in ("pr", "pet", "vpd", "srad") else a) - o) / s).astype("i2")
        v[:] = raw
        q = ds.createVariable(name + "_qc", "u1", ("time",))
        qq = np.zeros(days.size, "u1")
        qq[100] = 4                                         # one no-data day per variable
        q[:] = qq
        if name == "pr":
            v[100] = np.full((lat.size, lon.size), 32767, "i2")
    ds.processing_status = "complete"
    ds.close()
    return p


def drought(root, rng):
    """Lab drought product on the synthetic gridMET cells: spi90d every 5 days, packed int16 0.01, with _qc."""
    p = f"{root}/GRIDMET_LAB/GRIDMET_DROUGHT_Ogallala_lab.nc"
    import netCDF4 as nc
    with nc.Dataset(f"{root}/GRIDMET_LAB/GRIDMET_Ogallala_lab.nc") as g:
        lat, lon = np.asarray(g["lat"][:]), np.asarray(g["lon"][:])
    days = np.arange(np.datetime64("2004-01-03"), np.datetime64("2014-01-01"), 5)
    ds = _nc(p)
    ds.createDimension("time", days.size)
    ds.createDimension("lat", lat.size)
    ds.createDimension("lon", lon.size)
    tv = ds.createVariable("time", "f8", ("time",))
    tv.units, tv.calendar = "days since 1900-01-01 00:00:00", "standard"
    tv[:] = (days - np.datetime64("1900-01-01")).astype(int)
    ds.createVariable("lat", "f8", ("lat",))[:] = lat
    ds.createVariable("lon", "f8", ("lon",))[:] = lon
    yr = days.astype("datetime64[Y]").astype(int) + 1970
    mo = days.astype("datetime64[M]").astype(int) % 12 + 1
    base = np.where((yr == 2012) & (mo >= 5) & (mo <= 10), -1.5, np.where((yr == 2013) & (mo >= 5) & (mo <= 9), 1.2, 0.0))
    spi = base[:, None, None] + rng.normal(0, 0.1, (days.size, lat.size, lon.size))
    v = ds.createVariable("spi90d", "i2", ("time", "lat", "lon"), fill_value=32767)
    v.set_auto_maskandscale(False)
    v.scale_factor = np.float32(0.01)
    v[:] = np.round(spi / 0.01).astype("i2")
    ds.createVariable("spi90d_qc", "u1", ("time",))[:] = np.zeros(days.size, "u1")
    ds.processing_status = "complete"
    ds.close()
    return p


def hres(root):
    p = f"{root}/HRES_WTD/HRES_WTD_Ogallala.nc"
    lo0, la0, lo1, la1 = lonlat_bounds(0.01)
    step = 0.000239
    lat = np.arange(la1, la0, -step)
    lon = np.arange(lo0, lo1, step)
    LO, LA = np.meshgrid(lon, lat)
    val, X, Y = valley_of_lonlat(LO, LA)
    depth = np.where(val, 3.0, 45.0 + 10 * ((X - X0) / (W * 30.0))).astype("f4")
    ds = _nc(p)
    ds.createDimension("time", None)
    ds.createDimension("lat", lat.size)
    ds.createDimension("lon", lon.size)
    ds.createVariable("lat", "f8", ("lat",))[:] = lat
    ds.createVariable("lon", "f8", ("lon",))[:] = lon
    b1 = ds.createVariable("b1", "f4", ("time", "lat", "lon"), fill_value=np.finfo("f4").max)
    b1[0] = depth
    si = ds.createVariable("source_image", "i2", ("time", "lat", "lon"), fill_value=0)
    si[0] = np.ones(depth.shape, "i2")
    ds.processing_status = "complete"
    ds.close()
    return p


def dem(root):
    import rasterio
    from rasterio.transform import from_origin
    p = f"{root}/GLO_DEM/GLO30_DEM_Ogallala_1arcsec.tif"
    os.makedirs(os.path.dirname(p), exist_ok=True)
    lo0, la0, lo1, la1 = lonlat_bounds(0.05)
    top, bot = int(np.ceil(la1 * 3600)), int(np.floor(la0 * 3600))
    left, right = int(np.floor(lo0 * 3600)), int(np.ceil(lo1 * 3600))
    lat = (top - np.arange(top - bot + 1)) / 3600
    lon = (left + np.arange(right - left + 1)) / 3600
    LO, LA = np.meshgrid(lon, lat)
    val, X, Y = valley_of_lonlat(LO, LA)
    z = 900 + 0.004 * (X - X0) + 0.002 * (Y - Y0) + 25 * np.abs((X - X0) / 30 - 60) ** 0.5
    z = np.where(val, z - 15, z).astype("f4")
    tr = from_origin((left - 0.5) / 3600, (top + 0.5) / 3600, 1 / 3600, 1 / 3600)
    with rasterio.open(p, "w", driver="GTiff", height=z.shape[0], width=z.shape[1], count=1, dtype="float32",
                       crs="EPSG:4326", transform=tr, nodata=-9999) as r:
        r.write(z, 1)
    return p


def _albers_grid(step, margin):
    from pyproj import Transformer
    tr = Transformer.from_crs(32614, 5070, always_xy=True)
    xs = [X0 - margin, X0 + W * 30 + margin, X0 - margin, X0 + W * 30 + margin]
    ys = [Y0 + margin, Y0 + margin, Y0 - H * 30 - margin, Y0 - H * 30 - margin]
    ax, ay = tr.transform(xs, ys)
    x0 = np.floor(min(ax) / step) * step + (15.0 if step == 30 else 15.0)
    y1 = np.ceil(max(ay) / step) * step + (15.0 if step == 30 else 15.0)
    x = np.arange(x0, max(ax), step)
    y = np.arange(y1, min(ay), -step)
    return x, y


def lanid(root):
    p = f"{root}/LANID/LANID_v3_30m_1997-2017_Ogallala.nc"
    from pyproj import Transformer
    x, y = _albers_grid(30, 1500)
    X, Y = np.meshgrid(x, y)
    ux, uy = Transformer.from_crs(5070, 32614, always_xy=True).transform(X, Y)
    corn = ((ux - X0) / 30 > 70) & ((ux - X0) / 30 < 110) & ((Y0 - uy) / 30 > 40) & ((Y0 - uy) / 30 < 80)
    inside = ((ux - X0) / 30 > -40) & ((Y0 - uy) / 30 > -40)          # LANID coverage ends west/north of the frame
    ds = _nc(p)
    ds.createDimension("time", 21)
    ds.createDimension("y", y.size)
    ds.createDimension("x", x.size)
    ds.createVariable("time", "i4", ("time",))[:] = np.arange(1997, 2018)
    ds.createVariable("x", "f8", ("x",))[:] = x
    ds.createVariable("y", "f8", ("y",))[:] = y
    ds.createVariable("crs", "i4")
    ds.createVariable("roi_mask", "i1", ("y", "x"))[:] = inside.astype("i1")
    irr = np.where(inside, corn.astype("u1"), 255).astype("u1")
    v = ds.createVariable("irrigated", "u1", ("time", "y", "x"), fill_value=255)
    for k in range(21):
        v[k] = irr
    g = ds.createGroup("derived")
    m = g.createVariable("maxIrrExt", "u1", ("y", "x"), fill_value=255)
    occasional = ((ux - X0) / 30 > 130) & ((ux - X0) / 30 < 140) & ((Y0 - uy) / 30 > 120) & ((Y0 - uy) / 30 < 130)
    m[:] = np.where(inside, (corn | occasional).astype("u1"), 255).astype("u1")
    f = g.createVariable("irrFreq", "u1", ("y", "x"), fill_value=255)
    f[:] = np.where(inside, np.where(corn, 21, np.where(occasional, 2, 0)), 255).astype("u1")
    ds.processing_status = "complete"
    ds.close()
    return p


def gssurgo(root):
    p = f"{root}/USGS data/derived_usgs/GSSURGO_Ogallala_all_attrs_120m.nc"
    x, y = _albers_grid(120, 600)
    ds = _nc(p)
    ds.createDimension("y", y.size)
    ds.createDimension("x", x.size)
    ds.createVariable("x", "f8", ("x",))[:] = x
    ds.createVariable("y", "f8", ("y",))[:] = y
    ds.createVariable("crs", "i4")
    rng = np.random.default_rng(1)
    for v, base in dict(aws0_100=15.0, aws0_150=20.0, rootznaws=18.0, nccpi3all=0.5, soc0_30=4000.0).items():
        a = ds.createVariable(v, "f4", ("y", "x"), fill_value=np.float32(9.96921e36))
        a[:] = (base * (1 + 0.1 * rng.standard_normal((y.size, x.size)))).astype("f4")
    ds.close()
    return p


def gde(root):
    p = f"{root}/USGS data/derived_usgs/GDE_30arcsec_Ogallala_CF.nc"
    lo0, la0, lo1, la1 = lonlat_bounds(0.03)
    lat = np.arange(np.ceil(la1 * 120) / 120 - 1 / 240, la0, -1 / 120)
    lon = np.arange(np.floor(lo0 * 120) / 120 + 1 / 240, lo1, 1 / 120)
    ds = _nc(p)
    ds.createDimension("y", lat.size)
    ds.createDimension("x", lon.size)
    ds.createVariable("y", "f8", ("y",))[:] = lat
    ds.createVariable("x", "f8", ("x",))[:] = lon
    LO, LA = np.meshgrid(lon, lat)
    val, _, _ = valley_of_lonlat(LO, LA)
    ds.createVariable("GDE_frac_AA", "f4", ("y", "x"), fill_value=np.float32(9.96921e36))[:] = np.where(val, 0.6, 0.02)
    ds.createVariable("AA_frac_GA", "f4", ("y", "x"), fill_value=np.float32(9.96921e36))[:] = np.ones(LO.shape)
    ds.close()
    return p


def aquifer(root):
    import geopandas as gpd
    from shapely.geometry import box
    d = f"{root}/aquifer"
    os.makedirs(d, exist_ok=True)
    g = gpd.GeoDataFrame(dict(AQUIFER=["synthetic"]), geometry=[box(X0 - 5000, Y0 - H * 30 - 5000, X0 + W * 30 + 5000,
                                                                     Y0 + 5000)], crs=32614).to_crs(5070)
    p = f"{d}/hp_bound2010.shp"
    g.to_file(p)
    return p


def counties_zip(root):
    import geopandas as gpd
    from shapely.geometry import box
    d = f"{root}/_downloads"
    os.makedirs(d, exist_ok=True)
    xm = X0 + W * 30 / 2
    geoms = [box(X0 - 20000, Y0 - H * 30 - 20000, xm, Y0 + 20000), box(xm, Y0 - H * 30 - 20000, X0 + W * 30 + 20000, Y0 + 20000)]
    g = gpd.GeoDataFrame(dict(GEOID=["31001", "31003"], STATEFP=["31", "31"], NAME=["West", "East"]), geometry=geoms,
                         crs=32614).to_crs(4269)
    shp = f"{d}/cb_2023_us_county_500k.shp"
    g.to_file(shp)
    z = f"{d}/cb_2023_us_county_500k.zip"
    with zipfile.ZipFile(z, "w") as zz:
        for ext in ("shp", "shx", "dbf", "prj", "cpg"):
            f = f"{d}/cb_2023_us_county_500k.{ext}"
            if os.path.exists(f):
                zz.write(f, os.path.basename(f))
    return z


def nhd_gpkg(root):
    """Channels along the valley axis, no waterbodies, one HU8 polygon covering the frame."""
    import geopandas as gpd
    from shapely.geometry import LineString, box
    d = f"{root}/_downloads"
    os.makedirs(d, exist_ok=True)
    xc = X0 + 60 * 30
    line = LineString([(xc, Y0 + 3000), (xc, Y0 - H * 30 - 3000)])
    p = f"{d}/nhd_0000.gpkg"
    if os.path.exists(p):
        os.remove(p)
    gpd.GeoDataFrame(dict(NHDPlusID=[1.0], FCode=[46006], FType=[460], TotDASqKm=[50.0]), geometry=[line],
                     crs=32614).to_crs(4326).to_file(p, layer="channels")
    gpd.GeoDataFrame(dict(NHDPlusID=[], FType=[], FCode=[]), geometry=[], crs=4326).to_file(p, layer="waters")
    canal = LineString([(X0 + 90 * 30, Y0 - 12 * 30), (X0 + 150 * 30, Y0 - 12 * 30)])
    gpd.GeoDataFrame(dict(NHDPlusID=[2.0], FCode=[33600], FType=[336]), geometry=[canal],
                     crs=32614).to_crs(4326).to_file(p, layer="canals")
    gpd.GeoDataFrame(dict(HUC8=["10250000"], Name=["Synthetic"]),
                     geometry=[box(X0 - 1500, Y0 - H * 30 - 1500, X0 + W * 30 + 1500, Y0 + 1500)],
                     crs=32614).to_crs(4326).to_file(p, layer="hu8")
    return p


def build(root):
    """Write every synthetic input under root (idempotent). Returns {name: path}."""
    marker = f"{root}/_synthetic_complete.json"
    if os.path.exists(marker):
        got = json.load(open(marker))
        if got.get("_version") == VERSION:
            return got
        import shutil
        shutil.rmtree(root, ignore_errors=True)            # older synthetic inputs and every output made from them
    rng = np.random.default_rng(SEED)
    out = {}
    for tag in (LANDSAT_TAG, OPENET_TAG, CDL_TAG, NWI_TAG):
        frames_json(root, tag)
    for y in YEARS:
        out[f"landsat_{y}"] = landsat(root, y, rng)
        out[f"openet_{y}"] = openet(root, y, rng)
        out[f"cdl_{y}"] = cdl(root, y)
    out["nwi"] = nwi(root)
    out["gridmet"] = gridmet(root, rng)
    out["drought"] = drought(root, rng)
    out["hres_wtd"] = hres(root)
    out["dem"] = dem(root)
    out["lanid"] = lanid(root)
    out["gssurgo"] = gssurgo(root)
    out["gde"] = gde(root)
    out["aquifer"] = aquifer(root)
    out["counties_zip"] = counties_zip(root)
    out["nhd_gpkg"] = nhd_gpkg(root)
    out["_version"] = VERSION
    with open(marker, "w") as f:
        json.dump(out, f, indent=1)
    return out
