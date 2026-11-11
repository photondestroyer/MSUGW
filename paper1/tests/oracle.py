"""Independent re-computations used by the provenance tests.

Deliberately written WITHOUT importing paper1 (an AST test enforces it), from the file formats alone:
scale factors and fill values come from each file's attributes, months from its time variable, pixel
positions from its x / y coordinate arrays, and the tile corner from the tile-key formula of the Landsat
notebook. If the pipeline and this module disagree, one of them is wrong.
"""
import json

import netCDF4
import numpy as np
from pyproj import Transformer

TILE_M = 5120 * 30.0


def open_raw(path):
    ds = netCDF4.Dataset(path)
    ds.set_auto_maskandscale(False)
    return ds


def decoded(var, raw):
    a = raw.astype("float64")
    s, o = float(getattr(var, "scale_factor", 1.0)), float(getattr(var, "add_offset", 0.0))
    a = a * s + o
    if "_FillValue" in var.ncattrs():
        a[raw == var.getncattr("_FillValue")] = np.nan
    return a


def tile_corner(key):
    zone, i, j = int(key[1:3]), int(key[4:6]), int(key[7:9])
    return zone, 15.0 + i * TILE_M, 15.0 + (j + 1) * TILE_M


def block_of_lonlat(key, lon, lat):
    """(br, bc) of the 480 m block (aligned to the tile corner) holding a point."""
    zone, left, top = tile_corner(key)
    x, y = Transformer.from_crs(4326, 32600 + zone, always_xy=True).transform(lon, lat)
    return int((top - y) // 480), int((x - left) // 480)


def block_pixel_centres(key, br, bc):
    zone, left, top = tile_corner(key)
    xs = left + (bc * 16 + np.arange(16) + 0.5) * 30
    ys = top - (br * 16 + np.arange(16) + 0.5) * 30
    return np.meshgrid(xs, ys), zone


def file_rows_cols(ds, xs, ys):
    """Row / column of each pixel centre in a tile file, from its coordinate arrays (-1 = not in the file)."""
    fx, fy = np.asarray(ds["x"][:]), np.asarray(ds["y"][:])
    c = np.array([int(np.argmin(np.abs(fx - x))) if fx.min() - 15 <= x <= fx.max() + 15 else -1 for x in xs])
    r = np.array([int(np.argmin(np.abs(fy - y))) if fy.min() - 15 <= y <= fy.max() + 15 else -1 for y in ys])
    return r, c


def record_of_month(ds, year, month):
    t = netCDF4.num2date(ds["time"][:], ds["time"].units, only_use_cftime_datetimes=False)
    hits = [k for k, d in enumerate(t) if d.year == year and d.month == month]
    return hits[0] if hits else None


def tile_block_values(path, key, br, bc, year, month, var, nobs="n_obs"):
    """Decoded values of one variable over one block's pixels in one month (block found from coordinates)."""
    (X, Y), _ = block_pixel_centres(key, br, bc)
    with open_raw(path) as ds:
        rec = record_of_month(ds, year, month)
        if rec is None:
            return None
        r, c = file_rows_cols(ds, X[0], Y[:, 0])
        rr, cc = r[r >= 0], c[c >= 0]
        if rr.size == 0 or cc.size == 0:
            return np.array([])
        v = ds[var]
        a = decoded(v, np.asarray(v[rec, rr.min():rr.max() + 1, cc.min():cc.max() + 1]))
        if nobs in ds.variables:
            n = np.asarray(ds[nobs][rec, rr.min():rr.max() + 1, cc.min():cc.max() + 1])
            a[(n == 0) | (n == 255)] = np.nan
        return a.ravel()


def lab_month(path, var, iy, ix, year, month):
    """Monthly sum (pr, pet) or mean of a lab gridMET daily variable, with the lab QC bits applied."""
    with open_raw(path) as ds:
        t = netCDF4.num2date(ds["time"][:], ds["time"].units, only_use_cftime_datetimes=False)
        idx = [k for k, d in enumerate(t) if d.year == year and d.month == month]
        v = ds[var]
        a = decoded(v, np.asarray(v[idx[0]:idx[-1] + 1, iy, ix]))
        if var in ("tmmn", "tmmx"):
            a = a - 273.15
        q = np.asarray(ds[var + "_qc"][idx[0]:idx[-1] + 1]).astype(int)
        const, ident, nodata = q & 1 > 0, q & 2 > 0, q & 4 > 0
        bad = nodata | (ident & ~const)
        if var in ("srad", "tmmn", "tmmx", "vpd", "pet"):
            bad |= const
        a[bad] = np.nan
        good = np.isfinite(a).sum()
        if good < 0.9 * len(idx):
            return np.nan
        return np.nansum(a) * len(idx) / good if var in ("pr", "pet") else np.nanmean(a)


def lab_window_sum(path, iy, ix, first_day, n_days):
    """Rain over n_days ending the day before first_day."""
    with open_raw(path) as ds:
        t = netCDF4.num2date(ds["time"][:], ds["time"].units, only_use_cftime_datetimes=False)
        days = np.array([np.datetime64(d.strftime("%Y-%m-%d")) for d in t])
        k = int(np.searchsorted(days, np.datetime64(first_day)))
        v = ds["pr"]
        a = decoded(v, np.asarray(v[k - n_days:k, iy, ix]))
        q = np.asarray(ds["pr_qc"][k - n_days:k]).astype(int)
        a[(q & 4 > 0) | ((q & 2 > 0) & ~(q & 1 > 0))] = np.nan
        good = np.isfinite(a).sum()
        return np.nansum(a) * n_days / good if good >= 0.9 * n_days else np.nan


def roi_of_block(landsat_path, key, br, bc):
    (X, Y), _ = block_pixel_centres(key, br, bc)
    with open_raw(landsat_path) as ds:
        r, c = file_rows_cols(ds, X[0], Y[:, 0])
        m = np.zeros((16, 16), bool)
        for i in range(16):
            for j in range(16):
                if r[i] >= 0 and c[j] >= 0:
                    m[i, j] = int(ds["roi_mask"][r[i], c[j]]) == 1
    return m


def wtd_block(hres_path, key, br, bc, roi):
    """min / mean of HRES-WTD at the block's aquifer pixel centres (cell containing each centre)."""
    (X, Y), zone = block_pixel_centres(key, br, bc)
    lon, lat = Transformer.from_crs(32600 + zone, 4326, always_xy=True).transform(X[roi], Y[roi])
    with open_raw(hres_path) as ds:
        yn = "lat" if "lat" in ds.variables else "y"
        xn = "lon" if "lon" in ds.variables else "x"
        la, lo = np.asarray(ds[yn][:]), np.asarray(ds[xn][:])
        vals = []
        for a, b in zip(lat, lon):
            i, j = int(np.argmin(np.abs(la - a))), int(np.argmin(np.abs(lo - b)))
            v = float(ds["b1"][0, i, j] if ds["b1"].ndim == 3 else ds["b1"][i, j])
            s = int(ds["source_image"][0, i, j] if ds["source_image"].ndim == 3 else ds["source_image"][i, j])
            fill = float(ds["b1"].getncattr("_FillValue"))
            vals.append(np.nan if (v == np.float32(fill) or s == 0 or not np.isfinite(v)) else v)
    v = np.array(vals)
    return np.nanmin(v), np.nanmean(v), np.isfinite(v).mean()


def dem_block_mean(dem_path, key, br, bc, roi):
    import rasterio
    (X, Y), zone = block_pixel_centres(key, br, bc)
    lon, lat = Transformer.from_crs(32600 + zone, 4326, always_xy=True).transform(X[roi], Y[roi])
    with rasterio.open(dem_path) as r:
        z = [r.read(1, window=((row, row + 1), (col, col + 1)))[0, 0] for row, col in
             (r.index(a, b) for a, b in zip(lon, lat))]
        z = np.array(z, "float64")
        z[z == r.nodata] = np.nan
    return np.nanmean(z)


def lanid_block_distance(lanid_path, key, br, bc, roi, reach_px=30):
    """Distance (m) from the block's aquifer pixels to the nearest tile pixel whose centre falls in an
    ever-irrigated LANID cell (derived/maxIrrExt == 1), searching reach_px pixels around the block."""
    zone, left, top = tile_corner(key)
    rows = br * 16 + np.arange(-reach_px, 16 + reach_px)
    cols = bc * 16 + np.arange(-reach_px, 16 + reach_px)
    xs, ys = left + (cols + 0.5) * 30, top - (rows + 0.5) * 30
    XX, YY = np.meshgrid(xs, ys)
    ax, ay = Transformer.from_crs(32600 + zone, 5070, always_xy=True).transform(XX, YY)
    with open_raw(lanid_path) as ds:
        lx, ly = np.asarray(ds["x"][:]), np.asarray(ds["y"][:])
        ci = np.round((ax - lx[0]) / (lx[1] - lx[0])).astype(int)
        ri = np.round((ay - ly[0]) / (ly[1] - ly[0])).astype(int)
        ok = (ci >= 0) & (ci < lx.size) & (ri >= 0) & (ri < ly.size)
        v = ds.groups["derived"]["maxIrrExt"]
        r0, r1, c0, c1 = ri[ok].min(), ri[ok].max() + 1, ci[ok].min(), ci[ok].max() + 1
        win = np.asarray(v[r0:r1, c0:c1])
        irr = np.zeros(XX.shape, bool)
        irr[ok] = win[ri[ok] - r0, ci[ok] - c0] == 1
    inner = (slice(reach_px, reach_px + 16), slice(reach_px, reach_px + 16))
    bx, by = XX[inner][roi], YY[inner][roi]
    ix, iy = XX[irr], YY[irr]
    if ix.size == 0 or bx.size == 0:
        return np.inf
    d = np.sqrt((bx[:, None] - ix[None]) ** 2 + (by[:, None] - iy[None]) ** 2)
    return float(d.min())


def hand_trace(hand_path, i, j, max_steps=200000):
    """Follow the stored D8 directions from cell (i, j) to a channel cell; HAND = dem(i,j) - dem(channel)."""
    step = {64: (-1, 0), 128: (-1, 1), 1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1)}
    with open_raw(hand_path) as ds:
        fd = np.asarray(ds["fdir"][:])
        ch = np.asarray(ds["channel"][:]).astype(bool)
        dem = np.asarray(ds["dem"][:]).astype("float64")
        h = decoded(ds["hand"], np.asarray(ds["hand"][i, j]).reshape(1))[0]
    r, c = i, j
    for _ in range(max_steps):
        if ch[r, c]:
            return h, dem[i, j] - dem[r, c]
        code = int(fd[r, c])
        if code not in step:
            return h, np.nan
        r, c = r + step[code][0], c + step[code][1]
        if not (0 <= r < fd.shape[0] and 0 <= c < fd.shape[1]):
            return h, np.nan
    return h, np.nan


def frame_of(path):
    with open_raw(path) as ds:
        return json.loads(ds.getncattr("frame")) if "frame" in ds.ncattrs() else None
