"""S0.6: height above nearest drainage (HAND) from Copernicus GLO-30 and NHDPlus HR channels.

1. NHDPlus HR vector geodatabases (one per HU4) listed by the TNM Access API for the area processed
   (or config.HU4_CODES), downloaded to scratch, reduced to an aquifer subset GeoPackage on Drive:
   channels  = NHDFlowline with FCode in STREAM_FCODES (+ 46007 when INCLUDE_EPHEMERAL) and
               VAA TotDASqKm >= MIN_DRAINAGE_KM2; canals/ditches (336), pipelines (428) and underground
               conduits (420) never qualify because their FCodes are not listed
   waters    = NHDWaterbody with FType in WATERBODY_FTYPES (lake/pond, reservoir; no playas, no marsh)
   canals    = NHDFlowline with FType in CANAL_FTYPES (canal / ditch; plan v3 canal buffer of P1.1, never channels)
   hu8       = WBDHU8 polygons
2. Per HU8 (+ HU8_BUFFER_M): DEM window on the 1-arc-second lattice of GLO-30 (cell centres on whole
   arc-seconds) from the public COG tiles (unclipped, so windows may cross the aquifer edge; the Drive
   file is clipped to the polygon), channels and water rasterised (all touched), channel cells lowered
   by BURN_DEPTH_M, pits and depressions filled, flats resolved, D8 directions (pysheds 0.5) computed
   with the TRUE metric cell size (dx = 1/3600 deg * 111320 m * cos(lat), dy = 1/3600 deg * 110574 m;
   pysheds takes D8 spacing from the affine), HAND = unburned elevation - elevation of the channel cell
   the D8 path reaches. Cells inside the HU8 are kept; a path that leaves the window gives no HAND.
3. Output per HU8: hu8/HAND_<hu8>.nc  (lat, lon) hand int16 0.1 m (fill -32768), fdir uint8 (D8 code,
   0 none), channel uint8, inside uint8, facc float32 (km2, D8 flow accumulation: area draining through the
   cell from inside the HU8 window, so a measure of local convergence; drainage entering from upstream HU8s
   is not counted); hu8_index.gpkg (polygon, code, file) and canals.gpkg (all canal / ditch lines).
"""
import glob
import json
import os
import re
import shutil
import zipfile

import numpy as np

from .. import config as C
from ..core import paths, ncio, net
from ..core.errors import SourceError

DIRMAP = (64, 128, 1, 2, 4, 8, 16, 32)
ARCSEC = 1.0 / 3600.0
# D8 code -> (drow, dcol) in pysheds' dirmap order: N, NE, E, SE, S, SW, W, NW
D8_STEP = {64: (-1, 0), 128: (-1, 1), 1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1)}


# ---------------------------------------------------------------- NHDPlus HR
def hu4_codes(bbox, log=None):
    if C.HU4_CODES:
        return list(C.HU4_CODES)
    import requests
    params = dict(datasets="National Hydrography Dataset Plus High Resolution (NHDPlus HR)",
                  bbox=",".join(f"{v:.4f}" for v in bbox), max=500)
    r = requests.get(C.TNM_API, params=params, timeout=120, headers={"User-Agent": net.USER_AGENT})
    r.raise_for_status()
    codes = set()
    for it in r.json().get("items", []):
        m = re.search(r"NHDPLUS_H_(\d{4})_HU4", it.get("downloadURL", "") or "")
        if m:
            codes.add(m.group(1))
    if not codes:
        raise SourceError(f"TNM Access API returned no NHDPlus HR HU4 for bbox {bbox}")
    if log:
        log.info("NHDPlus HR HU4s for %s: %s", bbox, sorted(codes))
    return sorted(codes)


def _layer(gdb, want):
    import pyogrio
    names = [n for n, _ in pyogrio.list_layers(gdb)]
    for n in names:
        if n.lower() == want.lower():
            return n
    raise SourceError(f"{gdb}: layer {want} not found (layers: {names})")


def nhd_subset(hu4, log=None):
    """<shared>/nhdplus_hr/<hu4>.gpkg with layers channels, waters, hu8 (EPSG:4326, 2-D)."""
    import geopandas as gpd
    import pyogrio
    import shapely
    d = paths.shared_dir("nhdplus_hr")
    sig = dict(STREAM_FCODES=C.STREAM_FCODES, INCLUDE_EPHEMERAL=C.INCLUDE_EPHEMERAL,
               MIN_DRAINAGE_KM2=C.MIN_DRAINAGE_KM2, WATERBODY_FTYPES=C.WATERBODY_FTYPES, CANAL_FTYPES=C.CANAL_FTYPES)
    import hashlib
    tag = hashlib.md5(json.dumps(sig, sort_keys=True).encode()).hexdigest()[:8]
    out = f"{d}/{hu4}_{tag}.gpkg"
    meta_p = out + ".json"
    if os.path.exists(out) and os.path.exists(meta_p):
        return out, json.load(open(meta_p))
    raw = paths.scratch_dir("nhd_raw")
    url = C.NHD_GDB_URL.format(hu4=hu4)
    z, meta = net.download(url, f"{raw}/{os.path.basename(url)}", log=log)
    ex = f"{raw}/{hu4}"
    if not glob.glob(f"{ex}/*.gdb"):
        with zipfile.ZipFile(z) as zz:
            zz.extractall(ex)
    gdbs = glob.glob(f"{ex}/*.gdb") or glob.glob(f"{ex}/*/*.gdb")
    if not gdbs:
        raise SourceError(f"{z}: no .gdb inside")
    gdb = gdbs[0]
    fl = pyogrio.read_dataframe(gdb, layer=_layer(gdb, "NHDFlowline"), columns=["NHDPlusID", "FCode", "FType"])
    vaa = pyogrio.read_dataframe(gdb, layer=_layer(gdb, "NHDPlusFlowlineVAA"), columns=["NHDPlusID", "TotDASqKm"],
                                 read_geometry=False)
    fl = fl.merge(vaa, on="NHDPlusID", how="left")
    fcodes = list(C.STREAM_FCODES) + ([46007] if C.INCLUDE_EPHEMERAL else [])
    keep = fl["FCode"].isin(fcodes) & (fl["TotDASqKm"].fillna(0) >= C.MIN_DRAINAGE_KM2)
    ch = fl[keep].copy()
    cn = fl[fl["FType"].isin(C.CANAL_FTYPES)].copy()
    wb = pyogrio.read_dataframe(gdb, layer=_layer(gdb, "NHDWaterbody"), columns=["NHDPlusID", "FType", "FCode"])
    wb = wb[wb["FType"].isin(C.WATERBODY_FTYPES)].copy()
    hu8 = pyogrio.read_dataframe(gdb, layer=_layer(gdb, "WBDHU8"), columns=["HUC8", "Name"])
    tmp = f"{raw}/{os.path.basename(out)}"
    if os.path.exists(tmp):
        os.remove(tmp)
    for name, g in (("channels", ch), ("waters", wb), ("canals", cn), ("hu8", hu8)):
        g = g.set_geometry(shapely.force_2d(g.geometry.values)).to_crs(4326)
        g.to_file(tmp, layer=name, driver="GPKG")
    from ..core.drive import deliver
    deliver(tmp, out)
    meta = dict(hu4=hu4, source=meta, n_flowlines=int(len(fl)), n_channels=int(len(ch)), n_waters=int(len(wb)),
                n_canals=int(len(cn)), n_hu8=int(len(hu8)), rules=sig)
    with open(meta_p, "w") as f:
        json.dump(meta, f, indent=1)
    shutil.rmtree(ex, ignore_errors=True)
    os.remove(z)
    if log:
        log.info("NHDPlus HR %s: %d of %d flowlines are channels, %d waterbodies, %d HU8",
                 hu4, len(ch), len(fl), len(wb), len(hu8))
    return out, meta


# ---------------------------------------------------------------- DEM tiles
def dem_tile(lat_s, lon_w, log=None):
    """Path of the 1-degree COG tile whose south edge is lat_s and west edge lon_w (None if none exists)."""
    ns, ew = ("N" if lat_s >= 0 else "S"), ("E" if lon_w >= 0 else "W")
    url = C.DEM_COG_URL.format(ns=ns, lat=abs(lat_s), ew=ew, lon=abs(lon_w))
    d = paths.shared_dir("glo30_tiles")
    dest = f"{d}/{os.path.basename(url)}"
    try:
        p, _ = net.download(url, dest, log=log)
    except FileNotFoundError:
        return None
    return p


def dem_window(top, bottom, left, right, log=None):
    """GLO-30 elevations for cell centres lat = k/3600 (k = top .. bottom, descending) and lon = m/3600
    (m = left .. right). Returns float64 array, NaN where no tile."""
    import rasterio
    H, W = top - bottom + 1, right - left + 1
    out = np.full((H, W), np.nan)
    used = []
    for lat_s in range(int(np.floor(bottom * ARCSEC)) - 1, int(np.floor(top * ARCSEC)) + 1):
        for lon_w in range(int(np.floor(left * ARCSEC)) - 1, int(np.floor(right * ARCSEC)) + 1):
            p = dem_tile(lat_s, lon_w, log)
            if p is None:
                continue
            with rasterio.open(p) as r:
                t = r.transform
                k0 = (t.f + t.e * 0.5) * 3600          # lattice index (arc-seconds) of the centre of row 0
                m0 = (t.c + t.a * 0.5) * 3600          # ... and of column 0
                if abs(k0 - round(k0)) > 1e-3 or abs(m0 - round(m0)) > 1e-3 or abs(abs(t.e) * 3600 - 1) > 1e-6:
                    raise SourceError(f"{p}: not on the whole-arc-second lattice")
                k0, m0 = int(round(k0)), int(round(m0))
                if abs(t.a * 3600 - 1) > 1e-6:
                    continue                                   # tiles north of 50 deg are wider; not needed here
                r0, r1 = max(k0 - top, 0), min(k0 - bottom + 1, r.height)
                c0, c1 = max(left - m0, 0), min(right - m0 + 1, r.width)
                if r1 <= r0 or c1 <= c0:
                    continue
                a = r.read(1, window=((r0, r1), (c0, c1))).astype("float64")
                if r.nodata is not None:
                    a[a == r.nodata] = np.nan
                rr, cc = r0 - k0 + top, m0 + c0 - left   # window row/column of the first cell read
                out[rr:rr + (r1 - r0), cc:cc + (c1 - c0)] = a
                used.append(os.path.basename(p))
    return out, sorted(set(used))


# ---------------------------------------------------------------- HAND for one HU8
def compute_hand(dem, channel, lat_c):
    """dem (H, W) unburned elevations on the 1" lattice, channel bool -> hand, fdir (pysheds)."""
    from affine import Affine
    from pysheds.grid import Grid
    from pysheds.sview import Raster, ViewFinder
    dxm = ARCSEC * 111320.0 * np.cos(np.radians(lat_c))
    dym = ARCSEC * 110574.0
    bad = ~np.isfinite(dem)
    z = np.where(bad, np.nanmax(dem) + 1000.0, dem)
    vf = ViewFinder(affine=Affine(dxm, 0, 0, 0, -dym, 0), shape=dem.shape, nodata=np.nan)
    grid = Grid(viewfinder=vf)
    burned = Raster(np.where(channel, z - C.BURN_DEPTH_M, z), viewfinder=vf)
    filled = grid.fill_pits(burned)
    filled = grid.fill_depressions(filled)
    flat = grid.resolve_flats(filled)
    fdir = grid.flowdir(flat, dirmap=DIRMAP)
    vf_mask = ViewFinder(affine=vf.affine, shape=dem.shape, nodata=False)    # a bool raster cannot hold NaN
    hand = grid.compute_hand(fdir, Raster(z, viewfinder=vf), Raster(channel.astype(bool), viewfinder=vf_mask),
                             dirmap=DIRMAP, nodata_out=np.nan)
    hand = np.asarray(hand, "float64")
    fd = np.asarray(fdir)
    fd = np.where(np.isin(fd, DIRMAP), fd, 0).astype("uint8")
    hand[bad] = np.nan
    return hand, fd


def d8_accumulation(fdir, valid):
    """Number of valid cells draining through each cell (itself included) along D8 codes (DIRMAP).
    A cell whose code points outside the array, or code 0, is an outlet. Cells on a cycle stay at their
    own count (resolved flats have none). Pure numpy: processed in topological order, frontier by frontier."""
    H, W = fdir.shape
    n = H * W
    f = np.asarray(fdir).ravel()
    dr, dc = np.zeros(n, "int64"), np.zeros(n, "int64")
    for code, (a, b) in D8_STEP.items():
        m = f == code
        dr[m], dc[m] = a, b
    i = np.arange(n)
    r, c = i // W + dr, i % W + dc
    has = (f > 0) & (r >= 0) & (r < H) & (c >= 0) & (c < W)
    down = np.where(has, r * W + c, -1)
    acc = np.asarray(valid, bool).ravel().astype("float64")
    indeg = np.bincount(down[has], minlength=n).astype("int64")
    front = np.nonzero(indeg == 0)[0]
    while front.size:
        d = down[front]
        m = d >= 0
        src, dst = front[m], d[m]
        np.add.at(acc, dst, acc[src])
        np.subtract.at(indeg, dst, 1)
        dst = np.unique(dst)
        front = dst[indeg[dst] == 0]
    return acc.reshape(H, W)


def hu8_window(geom4326):
    import geopandas as gpd
    b = gpd.GeoSeries([geom4326], crs=4326).to_crs(5070).buffer(C.HU8_BUFFER_M).to_crs(4326).total_bounds
    left, bottom, right, top = b
    return (int(np.ceil(top * 3600)), int(np.floor(bottom * 3600)), int(np.floor(left * 3600)),
            int(np.ceil(right * 3600)))


def build_hu8(code, geom, channels, waters, out_path, scratch, prov, log=None):
    from rasterio import features
    from rasterio.transform import from_origin
    top, bottom, left, right = hu8_window(geom)
    dem, used = dem_window(top, bottom, left, right, log)
    H, W = dem.shape
    tr = from_origin((left - 0.5) * ARCSEC, (top + 0.5) * ARCSEC, ARCSEC, ARCSEC)
    shapes = [(gm, 1) for gm in channels.geometry if gm is not None and not gm.is_empty]
    shapes += [(gm, 1) for gm in waters.geometry if gm is not None and not gm.is_empty]
    channel = features.rasterize(shapes, out_shape=(H, W), transform=tr, fill=0, dtype="uint8",
                                 all_touched=True).astype(bool) if shapes else np.zeros((H, W), bool)
    inside = features.rasterize([(geom, 1)], out_shape=(H, W), transform=tr, fill=0, dtype="uint8",
                                all_touched=False).astype(bool)
    lat_c = (top + bottom) / 2 * ARCSEC
    hand, fdir = compute_hand(dem, channel, lat_c)
    cell_km2 = (ARCSEC * 111320.0 * np.cos(np.radians(lat_c))) * (ARCSEC * 110574.0) / 1e6
    facc = d8_accumulation(fdir, np.isfinite(dem)) * cell_km2
    keep = inside & np.isfinite(hand)
    hv = np.where(keep, np.clip(np.round(hand * 10), -32767, 32767), -32768).astype("int16")
    lat = (top - np.arange(H)) * ARCSEC
    lon = (left + np.arange(W)) * ARCSEC
    variables = dict(
        lat=dict(dims=("lat",), data=lat, dtype="f8", attrs=dict(units="degrees_north", long_name="cell centre")),
        lon=dict(dims=("lon",), data=lon, dtype="f8", attrs=dict(units="degrees_east", long_name="cell centre")),
        hand=dict(dims=("lat", "lon"), data=hv, dtype="i2", fill=-32768,
                  attrs=dict(units="m", scale_factor=np.float32(0.1), long_name="height above nearest drainage")),
        fdir=dict(dims=("lat", "lon"), data=np.where(inside, fdir, 0).astype("u1"), dtype="u1",
                  attrs=dict(long_name="D8 flow direction (pysheds dirmap 64 128 1 2 4 8 16 32 = N NE E SE S SW W NW)")),
        channel=dict(dims=("lat", "lon"), data=channel.astype("u1"), dtype="u1"),
        facc=dict(dims=("lat", "lon"), data=np.where(inside, facc, -9999).astype("f4"), dtype="f4", fill=np.float32(-9999),
                  attrs=dict(units="km2", long_name="D8 flow accumulation inside the HU8 window (local convergence)")),
        inside=dict(dims=("lat", "lon"), data=inside.astype("u1"), dtype="u1"),
        dem=dict(dims=("lat", "lon"), data=np.where(np.isfinite(dem), dem, -9999).astype("f4"), dtype="f4",
                 fill=np.float32(-9999), attrs=dict(units="m", long_name="GLO-30 elevation (unburned)")),
    )
    n_in = int(inside.sum())
    stats = dict(hu8=code, cells_inside=n_in, cells_with_hand=int(keep.sum()), dem_tiles=used,
                 channel_cells_inside=int((channel & inside).sum()),
                 no_hand_share=float(1 - keep.sum() / max(n_in, 1)))
    prov = dict(prov, extra=stats)
    ncio.write_nc(out_path, dict(lat=H, lon=W), variables, dict(hu8=code, title=f"HAND for HU8 {code}"), prov, scratch)
    return stats


def run(ctx):
    """Stage `hand`."""
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import box
    from . import tiles as TL
    from pyproj import Transformer
    T, G = TL.load_table(ctx.dep("tiles"))
    boxes = []
    for t, g in G.items():
        tr = Transformer.from_crs(g.epsg, 4326, always_xy=True)
        x0 = g.left + g.bc0 * C.BLOCK_PX * C.PIXEL_M
        y1 = g.top - g.br0 * C.BLOCK_PX * C.PIXEL_M
        x1 = x0 + g.nbx * C.BLOCK_PX * C.PIXEL_M
        y0 = y1 - g.nby * C.BLOCK_PX * C.PIXEL_M
        lo, la = tr.transform([x0, x1, x0, x1], [y0, y0, y1, y1])
        boxes.append(box(min(lo), min(la), max(lo), max(la)))
    region = gpd.GeoSeries(boxes, crs=4326).union_all()
    aq = gpd.read_file(paths.aquifer_shapefile()).to_crs(4326).union_all()
    region = region.intersection(aq.buffer(0.01))
    codes = hu4_codes(region.bounds, ctx.log)
    subsets = {}
    for hu4 in codes:
        ctx.check_budget()
        subsets[hu4] = nhd_subset(hu4, ctx.log)
        ctx.unit(f"hu4/{hu4}", "done", **{k: v for k, v in subsets[hu4][1].items() if k.startswith("n_")})
    hu8 = pd.concat([gpd.read_file(p, layer="hu8") for p, _ in subsets.values()], ignore_index=True)
    hu8 = hu8[hu8.intersects(region)].drop_duplicates("HUC8").reset_index(drop=True)
    chans = pd.concat([gpd.read_file(p, layer="channels") for p, _ in subsets.values()], ignore_index=True)
    write_canals(ctx, [p for p, _ in subsets.values()])
    waters = pd.concat([gpd.read_file(p, layer="waters") for p, _ in subsets.values()], ignore_index=True)
    os.makedirs(f"{ctx.out}/hu8", exist_ok=True)
    ci, wi = chans.sindex, waters.sindex
    results = []
    for _, row in hu8.iterrows():
        ctx.check_budget()
        code = row["HUC8"]
        p = f"{ctx.out}/hu8/HAND_{code}.nc"
        if ncio.is_complete(p, ctx.sig):
            continue
        top, bottom, left, right = hu8_window(row.geometry)
        win = box(left * ARCSEC, bottom * ARCSEC, right * ARCSEC, top * ARCSEC)
        ch = chans.iloc[ci.query(win, predicate="intersects")]
        wb = waters.iloc[wi.query(win, predicate="intersects")]
        prov = ctx.record({f"nhd_{k}": v[1]["source"] for k, v in subsets.items()}, dict(hu8=code))
        st = build_hu8(code, row.geometry, ch, wb, p, ctx.scratch, prov, ctx.log)
        results.append(st)
        ctx.unit(f"hu8/{code}", "done", **{k: v for k, v in st.items() if k != "dem_tiles"})
        ctx.log.info("HU8 %s: %d cells, %.1f %% without HAND", code, st["cells_inside"], 100 * st["no_hand_share"])
    idx = hu8.rename(columns={"HUC8": "hu8"})[["hu8", "Name", "geometry"]].copy()
    idx["file"] = [f"hu8/HAND_{c}.nc" for c in idx["hu8"]]
    tmp = f"{ctx.scratch}/hu8_index.gpkg"
    if os.path.exists(tmp):
        os.remove(tmp)
    idx.to_file(tmp, driver="GPKG")
    from ..core.drive import deliver
    deliver(tmp, f"{ctx.out}/hu8_index.gpkg")
    ctx.complete(dict(hu4=codes, hu8=int(len(idx)), computed_now=len(results)))


def write_canals(ctx, gpkgs):
    """All canal / ditch lines of the HU4 subsets -> <out>/canals.gpkg (EPSG:4326; may be empty)."""
    import geopandas as gpd
    import pandas as pd
    import pyogrio
    parts = []
    for p in gpkgs:
        if "canals" in [n for n, _ in pyogrio.list_layers(p)]:
            parts.append(gpd.read_file(p, layer="canals")[["geometry"]])
        else:
            ctx.log.warning("%s has no canals layer (subset made before plan v3?): no canal buffer there", p)
    g = (pd.concat(parts, ignore_index=True) if parts else gpd.GeoDataFrame(geometry=[], crs=4326))
    g = gpd.GeoDataFrame(g, geometry="geometry", crs=4326)
    tmp = f"{ctx.scratch}/canals.gpkg"
    if os.path.exists(tmp):
        os.remove(tmp)
    g.to_file(tmp, layer="canals", driver="GPKG")
    from ..core.drive import deliver
    deliver(tmp, f"{ctx.out}/canals.gpkg")
    return len(g)


def open_canals(hand_dir):
    import geopandas as gpd
    p = f"{hand_dir}/canals.gpkg"
    if not os.path.exists(p):
        return None
    return gpd.read_file(p, layer="canals")


# ---------------------------------------------------------------- readers
def open_hand(hand_dir, hu8):
    ds = ncio.open_ro(f"{hand_dir}/hu8/HAND_{hu8}.nc")
    lat, lon = np.asarray(ds["lat"][:]), np.asarray(ds["lon"][:])
    return dict(ds=ds, lat0=lat[0], lon0=lon[0], H=lat.size, W=lon.size)


def lookup(h, lon, lat, var="hand"):
    """Value of `var` (HAND in m, or facc in km2) of the 1" cells containing the points; NaN outside / none."""
    r = np.rint((h["lat0"] - lat) * 3600).astype("int64")
    c = np.rint((lon - h["lon0"]) * 3600).astype("int64")
    ok = (r >= 0) & (r < h["H"]) & (c >= 0) & (c < h["W"])
    out = np.full(np.shape(lon), np.nan)
    if ok.any():
        r0, r1, c0, c1 = r[ok].min(), r[ok].max() + 1, c[ok].min(), c[ok].max() + 1
        if var not in h["ds"].variables:
            return out
        var = h["ds"][var]
        win = np.asarray(var[r0:r1, c0:c1])
        out[ok] = ncio.decode(win, var)[r[ok] - r0, c[ok] - c0]
    return out
