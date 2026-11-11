"""Per-block static layers (S0.5 exposure, S0.6 HAND statistics, soils, terrain, LANID, NWI, GDE, county).

Every value is computed from the 30 m PIXEL CENTRES of the block (aquifer pixels only): each centre is
projected exactly (pyproj) into the source grid and takes the source cell that contains it - no
resampling, no GDAL warps. Per block (480 m):
  soils      mean of gSSURGO aws0_100, aws0_150, rootznaws (cm), nccpi3all, soc0_30 (EPSG:5070, 120 m)
  elevation  mean GLO-30 (m)
  WTD        HRES-WTD min, mean, p05, SD, fraction shallower than 2/5/10 m, coverage (valid share)
  HAND       fraction > HAND_MIN_M, min, p10, mean, coverage (from the `hand` stage)
  GDE        max of GDE_frac_AA over every 30" cell the block touches, mean at pixel centres, assessed share
  LANID      distance (m) from the block's pixels to the nearest ever-irrigated pixel (derived/maxIrrExt),
             and whether every pixel within IRR_BUFFER_M is assessed (inside LANID's polygon)
  NWI        wetland / riparian / mapped pixel counts (NWI 30 m tiles; absent tile -> nwi_available 0)
  plan v3    lanid3_dist_m: distance to pixels irrigated in >= LANID_FREQ_MIN_YEARS years (derived/irrFreq);
             canal_dist_m: distance (m) from the block's pixels to the nearest NHD canal / ditch (the `hand`
             stage's canals.gpkg; inf when none within CANAL_BUFFER_M + margin);
             facc_max_km2, facc_mean_km2: D8 flow accumulation at the block's pixels (from the HAND files)
  county     county of the block centre; bucket bin nearest to aws0_150
Also written: lanid_<tile>.npz with the per-pixel ever-irrigated and assessed masks of the block window
(bit-packed), used by the `classes` stage.
"""
import os

import numpy as np

from .. import config as C
from ..core import paths, ncio
from . import tiles as TL

STRIP_BLOCKS = 16                     # block rows per strip (256 pixel rows)


# ---------------------------------------------------------------- generic grid lookups
class RegularGrid:
    """A 1-D regular coordinate pair (centres) with lookups of the cell containing a point."""

    def __init__(self, ycoords, xcoords):
        self.y, self.x = np.asarray(ycoords, "float64"), np.asarray(xcoords, "float64")
        self.dy, self.dx = self.y[1] - self.y[0], self.x[1] - self.x[0]
        for c, d in ((self.y, self.dy), (self.x, self.dx)):
            if np.max(np.abs(np.diff(c) - d)) > 1e-6 * abs(d):
                raise ValueError("coordinate axis is not regular")

    def index(self, Y, X):
        r = np.floor((Y - (self.y[0] - self.dy / 2)) / self.dy).astype("int64")
        c = np.floor((X - (self.x[0] - self.dx / 2)) / self.dx).astype("int64")
        ok = (r >= 0) & (r < self.y.size) & (c >= 0) & (c < self.x.size)
        return r, c, ok


def gather(var, r, c, ok, t=None, fill=None):
    """Values of var at (r, c) read through one window read; outside -> fill (raw dtype)."""
    out = np.full(r.shape, fill if fill is not None else 0, dtype=var.dtype)
    if not ok.any():
        return out
    r0, r1, c0, c1 = r[ok].min(), r[ok].max() + 1, c[ok].min(), c[ok].max() + 1
    sl = (slice(r0, r1), slice(c0, c1)) if t is None else (t, slice(r0, r1), slice(c0, c1))
    win = np.asarray(var[sl])
    out[ok] = win[r[ok] - r0, c[ok] - c0]
    return out


def reduce_blocks(a, roi, nbx, how):
    """a, roi: (strip_rows, nbx*16) -> per block statistic over finite values at aquifer pixels."""
    nby = a.shape[0] // C.BLOCK_PX
    b = TL.blockify(a.astype("float64"), nby, nbx)
    m = TL.blockify(roi, nby, nbx) & np.isfinite(b)
    n = m.sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        if how == "mean":
            v = np.where(m, b, 0).sum(-1) / n
        elif how == "min":
            v = np.where(m, b, np.inf).min(-1)
        elif how == "max":
            v = np.where(m, b, -np.inf).max(-1)
        elif how == "sd":
            mu = np.where(m, b, 0).sum(-1) / n
            v = np.sqrt(np.where(m, (b - mu[..., None]) ** 2, 0).sum(-1) / n)
        elif how.startswith("p"):
            q = float(how[1:]) / 100
            bb = np.where(m, b, np.nan)
            v = np.nanquantile(bb, q, axis=-1) if n.max() > 0 else np.full(n.shape, np.nan)
        elif how == "n":
            return n.astype("float64")
        else:
            raise ValueError(how)
    v = np.where(n > 0, v, np.nan)
    return v


# ---------------------------------------------------------------- per-source block statistics
def soils_dem_wtd(g, roi_win, lonlat_strips, xy5070_strips, log=None):
    import rasterio
    out = {k: np.full((g.nby, g.nbx), np.nan, "float32") for k in
           C.STATIC_SOIL_VARS + ["elevation", "wtd_min", "wtd_mean", "wtd_p05", "wtd_sd", "wtd_coverage"]
           + [f"wtd_lt{d}m" for d in C.WTD_SHALLOW_M]}
    inputs = {}
    gs = paths.dataset_path("gssurgo")
    dm = paths.dataset_path("dem")
    wt = paths.dataset_path("hres_wtd")
    inputs.update(gssurgo=gs, dem=dm, hres_wtd=wt)
    with ncio.open_ro(gs) as dss, ncio.open_ro(wt) as dsw, rasterio.open(dm) as rd:
        gsg = RegularGrid(dss["y"][:], dss["x"][:])
        ylat = "lat" if "lat" in dsw.variables else "y"
        xlon = "lon" if "lon" in dsw.variables else "x"
        wg = RegularGrid(dsw[ylat][:], dsw[xlon][:])
        b1, srcimg = dsw["b1"], dsw["source_image"]
        b1_fill = ncio.fill_of(b1, np.finfo("f4").max)
        inv = ~rd.transform
        dem_nodata = rd.nodata
        for (ra, rb), (lon, lat), (X5, Y5) in zip(_strip_bounds(g), lonlat_strips(), xy5070_strips()):
            roi = roi_win[ra:rb]
            bs = slice(ra // C.BLOCK_PX, rb // C.BLOCK_PX)
            r, c, ok = gsg.index(Y5, X5)
            for v in C.STATIC_SOIL_VARS:
                var = dss[v]
                raw = gather(var, r, c, ok, fill=ncio.fill_of(var))
                out[v][bs] = reduce_blocks(ncio.decode(raw, var), roi, g.nbx, "mean")
            # DEM: pixel centre -> cell containing it (rasterio affine, floor)
            cc, rr = inv * (lon, lat)
            rr, cc = np.floor(rr).astype("int64"), np.floor(cc).astype("int64")
            okd = (rr >= 0) & (rr < rd.height) & (cc >= 0) & (cc < rd.width)
            if okd.any():
                r0, r1, c0, c1 = rr[okd].min(), rr[okd].max() + 1, cc[okd].min(), cc[okd].max() + 1
                win = rd.read(1, window=((r0, r1), (c0, c1))).astype("float64")
                z = np.full(rr.shape, np.nan)
                z[okd] = win[rr[okd] - r0, cc[okd] - c0]
                if dem_nodata is not None:
                    z[z == dem_nodata] = np.nan
                out["elevation"][bs] = reduce_blocks(z, roi, g.nbx, "mean")
            # HRES-WTD
            r, c, ok = wg.index(lat, lon)
            t = 0 if b1.ndim == 3 else None
            w = gather(b1, r, c, ok, t=t, fill=b1_fill).astype("float64")
            s = gather(srcimg, r, c, ok, t=t, fill=0)
            w[(w == np.float32(b1_fill)) | (s == 0) | ~ok | ~np.isfinite(w)] = np.nan
            out["wtd_min"][bs] = reduce_blocks(w, roi, g.nbx, "min")
            out["wtd_mean"][bs] = reduce_blocks(w, roi, g.nbx, "mean")
            out["wtd_p05"][bs] = reduce_blocks(w, roi, g.nbx, "p5")
            out["wtd_sd"][bs] = reduce_blocks(w, roi, g.nbx, "sd")
            nroi = TL.blockify(roi, rb // C.BLOCK_PX - ra // C.BLOCK_PX, g.nbx).sum(-1)
            nval = reduce_blocks(w, roi, g.nbx, "n")
            with np.errstate(invalid="ignore", divide="ignore"):
                out["wtd_coverage"][bs] = np.where(nroi > 0, nval / nroi, np.nan)
                for d in C.WTD_SHALLOW_M:
                    out[f"wtd_lt{d}m"][bs] = reduce_blocks(np.where(np.isfinite(w), (w < d).astype(float), np.nan),
                                                           roi, g.nbx, "mean")
    return out, inputs


def gde_stats(g, roi_win, lonlat_strips):
    p = paths.dataset_path("gde")
    out = dict(gde_max=np.full((g.nby, g.nbx), np.nan, "float32"), gde_mean=np.full((g.nby, g.nbx), np.nan, "float32"),
               gde_assessed=np.full((g.nby, g.nbx), np.nan, "float32"))
    from pyproj import Transformer
    tr = Transformer.from_crs(g.epsg, 4326, always_xy=True)
    with ncio.open_ro(p) as ds:
        grid = RegularGrid(ds["y"][:], ds["x"][:])
        fa, aa = ds["GDE_frac_AA"], ds["AA_frac_GA"]
        FA = ncio.decode(np.asarray(fa[:]), fa)
        AA = ncio.decode(np.asarray(aa[:]), aa)
    for (ra, rb), (lon, lat) in zip(_strip_bounds(g), lonlat_strips()):
        bs = slice(ra // C.BLOCK_PX, rb // C.BLOCK_PX)
        roi = roi_win[ra:rb]
        r, c, ok = grid.index(lat, lon)
        v = np.where(ok, FA[np.clip(r, 0, FA.shape[0] - 1), np.clip(c, 0, FA.shape[1] - 1)], np.nan)
        a = np.where(ok, AA[np.clip(r, 0, AA.shape[0] - 1), np.clip(c, 0, AA.shape[1] - 1)], np.nan)
        out["gde_mean"][bs] = reduce_blocks(v, roi, g.nbx, "mean")
        out["gde_assessed"][bs] = reduce_blocks(a, roi, g.nbx, "mean")
    # max over every 30" cell the block touches: block corners -> lon/lat bounding box
    xs, ys = g.block_centres()
    half = C.BLOCK_PX * C.PIXEL_M / 2
    X, Y = np.meshgrid(xs, ys)
    lons, lats = [], []
    for dx in (-half, half):
        for dy in (-half, half):
            lo, la = tr.transform(X + dx, Y + dy)
            lons.append(lo)
            lats.append(la)
    lo_min, lo_max = np.min(lons, 0), np.max(lons, 0)
    la_min, la_max = np.min(lats, 0), np.max(lats, 0)
    r_a, c_a, _ = grid.index(la_max, lo_min)
    r_b, c_b, _ = grid.index(la_min, lo_max)
    rmin, rmax = np.minimum(r_a, r_b), np.maximum(r_a, r_b)
    cmin, cmax = np.minimum(c_a, c_b), np.maximum(c_a, c_b)
    gm = np.full(X.shape, np.nan, "float32")
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            a0, a1 = max(rmin[i, j], 0), min(rmax[i, j] + 1, FA.shape[0])
            b0, b1 = max(cmin[i, j], 0), min(cmax[i, j] + 1, FA.shape[1])
            if a1 > a0 and b1 > b0:
                sub = FA[a0:a1, b0:b1]
                if np.isfinite(sub).any():
                    gm[i, j] = np.nanmax(sub)
    out["gde_max"] = gm
    return out, {"gde": p}


def lanid_distance(g, roi_win, log=None):
    """Distance to ever-irrigated pixels on the tile's own 30 m lattice (block window + margin)."""
    from pyproj import Transformer
    from scipy import ndimage
    p = paths.dataset_path("lanid")
    M = int(np.ceil(C.IRR_BUFFER_M / C.PIXEL_M)) + 1
    H, W = g.nby * C.BLOCK_PX + 2 * M, g.nbx * C.BLOCK_PX + 2 * M
    tr = Transformer.from_crs(g.epsg, 5070, always_xy=True)
    irr = np.zeros((H, W), bool)
    irr3 = np.zeros((H, W), bool)
    unk = np.zeros((H, W), bool)
    have_freq = True
    tc = g.bc0 * C.BLOCK_PX - M + np.arange(W)
    xs = g.left + (tc + 0.5) * C.PIXEL_M
    with ncio.open_ro(p) as ds:
        grid = RegularGrid(ds["y"][:], ds["x"][:])
        var = ds.groups["derived"]["maxIrrExt"]
        fill = ncio.fill_of(var, 255)
        fvar = ds.groups["derived"].variables.get("irrFreq")
        have_freq = fvar is not None
        ffill = ncio.fill_of(fvar, 255) if have_freq else None
        step = 256
        for a in range(0, H, step):
            b = min(a + step, H)
            tr_rows = g.br0 * C.BLOCK_PX - M + np.arange(a, b)
            ys = g.top - (tr_rows + 0.5) * C.PIXEL_M
            XX, YY = np.meshgrid(xs, ys)
            X5, Y5 = tr.transform(XX, YY)
            r, c, ok = grid.index(Y5, X5)
            v = gather(var, r, c, ok, fill=fill)
            irr[a:b] = (v == 1) & ok
            unk[a:b] = (v == fill) | ~ok
            if have_freq:
                fq = gather(fvar, r, c, ok, fill=ffill)
                irr3[a:b] = (fq != ffill) & (fq >= C.LANID_FREQ_MIN_YEARS) & ok
    d_irr = ndimage.distance_transform_edt(~irr) * C.PIXEL_M if irr.any() else np.full(irr.shape, np.inf)
    d_unk = ndimage.distance_transform_edt(~unk) * C.PIXEL_M if unk.any() else np.full(unk.shape, np.inf)
    inner = (slice(M, H - M), slice(M, W - M))
    di, du = d_irr[inner], d_unk[inner]
    out = dict(lanid_dist_m=_block_min(di, roi_win, g), lanid_unassessed_dist_m=_block_min(du, roi_win, g))
    if have_freq:
        d3 = ndimage.distance_transform_edt(~irr3) * C.PIXEL_M if irr3.any() else np.full(irr3.shape, np.inf)
        out["lanid3_dist_m"] = _block_min(d3[inner], roi_win, g)
    else:
        if log:
            log.warning("%s: LANID has no derived/irrFreq - lanid3_dist_m left NaN (IRR3 sensitivity fails)", p)
        out["lanid3_dist_m"] = np.full((g.nby, g.nbx), np.nan, "float32")
    pix = dict(ever=np.packbits(irr[inner]), assessed=np.packbits(~unk[inner]), shape=np.array(irr[inner].shape))
    return out, pix, {"lanid": p}


def canal_distance(g, roi_win, canals):
    """Distance from each block's aquifer pixels to the nearest canal / ditch line, on the tile's 30 m lattice
    (block window + a CANAL_BUFFER_M margin; lines rasterised with all_touched). inf where none is that close."""
    from rasterio import features
    from rasterio.transform import from_origin
    from scipy import ndimage
    M = int(np.ceil(C.CANAL_BUFFER_M / C.PIXEL_M)) + 1
    H, W = g.nby * C.BLOCK_PX + 2 * M, g.nbx * C.BLOCK_PX + 2 * M
    if canals is None or len(canals) == 0:
        return np.full((g.nby, g.nbx), np.inf, "float32")
    left = g.left + (g.bc0 * C.BLOCK_PX - M) * C.PIXEL_M
    top = g.top - (g.br0 * C.BLOCK_PX - M) * C.PIXEL_M
    tr = from_origin(left, top, C.PIXEL_M, C.PIXEL_M)
    geoms = [gm for gm in canals.to_crs(g.epsg).geometry if gm is not None and not gm.is_empty]
    if not geoms:
        return np.full((g.nby, g.nbx), np.inf, "float32")
    line = features.rasterize(((gm, 1) for gm in geoms), out_shape=(H, W), transform=tr, fill=0, dtype="uint8",
                              all_touched=True).astype(bool)
    if not line.any():
        return np.full((g.nby, g.nbx), np.inf, "float32")
    d = ndimage.distance_transform_edt(~line) * C.PIXEL_M
    return _block_min(d[M:H - M, M:W - M], roi_win, g)


def _block_min(a, roi, g):
    b = TL.blockify(a, g.nby, g.nbx)
    m = TL.blockify(roi, g.nby, g.nbx)
    v = np.where(m, b, np.inf).min(-1)
    return np.where(m.any(-1), v, np.nan).astype("float32")


def nwi_counts(g, roi_win):
    p = paths.tile_path("nwi", g.tile)
    z = np.zeros((g.nby, g.nbx), "uint16")
    if not os.path.exists(p):
        return dict(nwi_available=np.zeros((g.nby, g.nbx), "uint8"), nwi_wetland_px=z, nwi_riparian_px=z.copy(),
                    nwi_mapped_px=z.copy()), {}
    with ncio.open_ro(p) as ds:
        rows, cols = g.frame_rows(), g.frame_cols()
        sysv = ncio.read_window(ds["wetland_system"], rows, cols, 255)
        rip = ncio.read_window(ds["riparian"], rows, cols, 255)
        mp = ncio.read_window(ds["nwi_mapped"], rows, cols, 255)
    cnt = lambda m: TL.blockify((m & roi_win).astype("int16"), g.nby, g.nbx).sum(-1).astype("uint16")
    return dict(nwi_available=np.ones((g.nby, g.nbx), "uint8"),
                nwi_wetland_px=cnt((sysv > 0) & (sysv != 255)), nwi_riparian_px=cnt(rip == 1),
                nwi_mapped_px=cnt(mp == 1)), {"nwi": p}


def hand_stats(g, roi_win, lonlat_strips, hand_dir):
    """HAND at pixel centres: HU8 of each pixel (polygon), then that HU8's HAND raster."""
    import geopandas as gpd
    from rasterio import features
    from rasterio.transform import from_origin
    from . import hand as HD
    keys = ["hand_frac", "hand_min", "hand_p10", "hand_mean", "hand_coverage", "facc_max_km2", "facc_mean_km2"]
    out = {k: np.full((g.nby, g.nbx), np.nan, "float32") for k in keys}
    idx = gpd.read_file(f"{hand_dir}/hu8_index.gpkg").to_crs(g.epsg)
    cache = {}
    for (ra, rb), (lon, lat) in zip(_strip_bounds(g), lonlat_strips()):
        bs = slice(ra // C.BLOCK_PX, rb // C.BLOCK_PX)
        roi = roi_win[ra:rb]
        left = g.left + g.bc0 * C.BLOCK_PX * C.PIXEL_M
        top = g.top - (g.br0 * C.BLOCK_PX + ra) * C.PIXEL_M
        shape = (rb - ra, g.nbx * C.BLOCK_PX)
        tr = from_origin(left, top, C.PIXEL_M, C.PIXEL_M)
        which = features.rasterize(((geom, i + 1) for i, geom in enumerate(idx.geometry)), out_shape=shape,
                                   transform=tr, fill=0, dtype="int32", all_touched=False)
        h = np.full(shape, np.nan)
        fa = np.full(shape, np.nan)
        for k in np.unique(which[which > 0]):
            hu8 = idx.iloc[k - 1]["hu8"]
            if hu8 not in cache:
                cache[hu8] = HD.open_hand(hand_dir, hu8)
            sel = which == k
            h[sel] = HD.lookup(cache[hu8], lon[sel], lat[sel])
            fa[sel] = HD.lookup(cache[hu8], lon[sel], lat[sel], var="facc")
        valid = np.where(np.isfinite(h), 1.0, np.nan)
        out["hand_frac"][bs] = reduce_blocks(np.where(np.isfinite(h), (h > C.HAND_MIN_M).astype(float), np.nan),
                                             roi, g.nbx, "mean")
        out["hand_min"][bs] = reduce_blocks(h, roi, g.nbx, "min")
        out["hand_p10"][bs] = reduce_blocks(h, roi, g.nbx, "p10")
        out["hand_mean"][bs] = reduce_blocks(h, roi, g.nbx, "mean")
        out["facc_max_km2"][bs] = reduce_blocks(fa, roi, g.nbx, "max")
        out["facc_mean_km2"][bs] = reduce_blocks(fa, roi, g.nbx, "mean")
        nroi = TL.blockify(roi, rb // C.BLOCK_PX - ra // C.BLOCK_PX, g.nbx).sum(-1)
        with np.errstate(invalid="ignore", divide="ignore"):
            out["hand_coverage"][bs] = np.where(nroi > 0, reduce_blocks(valid, roi, g.nbx, "n") / nroi, np.nan)
    for v in cache.values():
        v["ds"].close()
    return out


def county_of_blocks(T, sl, counties):
    import geopandas as gpd
    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(T["lon"][sl], T["lat"][sl]), crs=4326)
    j = gpd.sjoin(pts, counties.to_crs(4326)[["GEOID", "geometry"]], how="left", predicate="within")
    j = j[~j.index.duplicated(keep="first")]
    geoid = j["GEOID"].reindex(range(len(pts)))
    lut = {gid: i for i, gid in enumerate(counties["GEOID"])}
    return np.array([lut.get(x, -1) if isinstance(x, str) else -1 for x in geoid], "int32")


# ---------------------------------------------------------------- strips of pixel coordinates
def _strip_bounds(g):
    H = g.nby * C.BLOCK_PX
    step = STRIP_BLOCKS * C.BLOCK_PX
    return [(a, min(a + step, H)) for a in range(0, H, step)]


def _coords(g, epsg):
    from pyproj import Transformer
    tr = Transformer.from_crs(g.epsg, epsg, always_xy=True)

    def gen():
        for ra, rb in _strip_bounds(g):
            xs, ys = g.pixel_centres((ra, rb))
            XX, YY = np.meshgrid(xs, ys)
            yield tr.transform(XX, YY)
    return gen


def roi_window(g):
    src = TL.first_source_file("landsat", g.tile, C.agg_years())
    with ncio.open_ro(src) as ds:
        return ncio.read_window(ds["roi_mask"], g.frame_rows(), g.frame_cols(), 0) == 1, src


def build_tile(g, T, sl, counties, hand_dir, log=None):
    roi, roi_src = roi_window(g)
    lonlat = _coords(g, 4326)
    xy5070 = _coords(g, 5070)
    V, inputs = {}, {"landsat_roi": roi_src}
    a, i = soils_dem_wtd(g, roi, lonlat, xy5070, log)
    V.update(a)
    inputs.update(i)
    a, i = gde_stats(g, roi, lonlat)
    V.update(a)
    inputs.update(i)
    a, pix, i = lanid_distance(g, roi, log)
    V.update(a)
    inputs.update(i)
    a, i = nwi_counts(g, roi)
    V.update(a)
    inputs.update(i)
    V.update(hand_stats(g, roi, lonlat, hand_dir))
    from . import hand as HD
    V["canal_dist_m"] = canal_distance(g, roi, HD.open_canals(hand_dir))
    V["county"] = county_of_blocks(T, sl, counties).reshape(g.nby, g.nbx)
    aws_mm = V["aws0_150"].astype("float64") * 10.0
    caps = np.asarray(C.BUCKET_CAPACITY_BINS_MM, "float64")
    V["bucket_bin"] = np.where(np.isfinite(aws_mm), np.abs(aws_mm[..., None] - caps).argmin(-1), -1).astype("int8")
    V["roi_px"] = TL.blockify(roi.astype("int16"), g.nby, g.nbx).sum(-1).astype("uint16")
    return V, pix, inputs


def run(ctx):
    """Stage `static`."""
    from ..core import provenance as P
    from . import counties as CO
    T, G = TL.load_table(ctx.dep("tiles"))
    counties = CO.load(ctx.dep("counties"))
    hand_dir = ctx.dep("hand")
    tiles = list(T["tiles"])
    for ti, t in enumerate(tiles):
        ctx.check_budget()
        g = G[t]
        p = f"{ctx.out}/{t}/P1_static_{t}.nc"
        if ncio.is_complete(p, ctx.sig) and os.path.exists(f"{ctx.out}/{t}/lanid_{t}.npz"):
            continue
        V, pix, inputs = build_tile(g, T, TL.tile_slice(T, ti), counties, hand_dir, ctx.log)
        os.makedirs(f"{ctx.out}/{t}", exist_ok=True)
        np.savez_compressed(f"{ctx.out}/{t}/lanid_{t}.npz", **pix)
        variables = {k: dict(dims=("by", "bx"), data=v, dtype=str(v.dtype)) for k, v in V.items()}
        variables["county_geoid"] = dict(dims=("county",), data=np.array(counties["GEOID"], dtype=object), dtype="str")
        prov = ctx.record({k: P.fingerprint(v) for k, v in inputs.items()}, dict(tile=t, window=g.to_dict()))
        ncio.write_nc(p, dict(by=g.nby, bx=g.nbx, county=len(counties)), variables,
                      dict(tile=t, title="Per-block static layers (paper1 Step 0)"), prov, ctx.scratch)
        ctx.unit(t, "done", nwi=int(V["nwi_available"].max()))
        ctx.log.info("%s: static layers written", t)
    ctx.complete(dict(tiles=len(tiles)))


def open_tile(stage_dir, tile):
    with ncio.open_ro(f"{stage_dir}/{tile}/P1_static_{tile}.nc") as ds:
        return {k: np.asarray(ds[k][...]) for k in ds.variables if k != "county_geoid"}


def open_lanid(stage_dir, tile):
    with np.load(f"{stage_dir}/{tile}/lanid_{tile}.npz") as z:
        shp = tuple(z["shape"])
        n = shp[0] * shp[1]
        return (np.unpackbits(z["ever"])[:n].reshape(shp).astype(bool),
                np.unpackbits(z["assessed"])[:n].reshape(shp).astype(bool))


def open_county_geoid(stage_dir, tile):
    """County GEOIDs (index = the `county` values of the tile's static file)."""
    with ncio.open_ro(f"{stage_dir}/{tile}/P1_static_{tile}.nc") as ds:
        return [str(x) for x in np.asarray(ds["county_geoid"][...])]
