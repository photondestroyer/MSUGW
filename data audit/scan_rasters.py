"""Full-read scan of the GeoTIFFs and GeoPackages under merged_datasets (read-only).

Per band: every pixel is read in windows; counts of nodata/NaN/zero, min/max/mean, percentiles of a
systematic sample, the value histogram for small-integer bands, and the same statistics restricted to
the aquifer's bounding box for global rasters. Exact zeros in a band with no declared nodata are
sampled and tested against the aquifer polygon (a zero outside the polygon is clip fill, not data).
"""
import json
import math
import os
import sys
import time
from collections import Counter

import geopandas as gpd
import numpy as np
import pyogrio
import rasterio
from rasterio.windows import Window, from_bounds
from shapely import contains_xy
from shapely.ops import unary_union

from audit_paths import OUT, rel, raster_files

POLY = gpd.read_file("G:/MSU_GWB/datasets/HPA_polygon/hp_bound2010.shp")
TARGET_PX = 8 * 1024 * 1024         # pixels per window
os.environ.setdefault("GDAL_CACHEMAX", "256")   # MB; GDAL's default block cache grew to ~11 GB on the global HSG raster


def windows(r, area=None):
    """Block-aligned windows of about TARGET_PX pixels, split in both directions for very wide rasters."""
    full = Window(0, 0, r.width, r.height) if area is None else area
    bh, bw = r.block_shapes[0]
    c0, c1 = int(full.col_off), int(full.col_off + full.width)
    r0_, r1 = int(full.row_off), int(full.row_off + full.height)
    width = c1 - c0
    if width * bh <= TARGET_PX:          # whole rows fit: strips
        cols = width
        rows = max(bh, (TARGET_PX // width) // bh * bh)
    else:                                # tiles of one block-row height group
        rows = bh if bh > 1 else max(1, TARGET_PX // 4096)
        cols = max(bw, (TARGET_PX // rows) // bw * bw)
    for r0 in range(r0_, r1, rows):
        for cc in range(c0, c1, cols):
            yield Window(cc, r0, min(cols, c1 - cc), min(rows, r1 - r0))


def band_stats(r, b, area=None, zero_probe=False, poly=None):
    nodata = r.nodatavals[b - 1]
    is_int = np.dtype(r.dtypes[b - 1]).kind in "iu"
    n = n_nod = n_nan = n_zero = n_valid = 0
    mn, mx, s1, s2 = math.inf, -math.inf, 0.0, 0.0
    hist = Counter()
    samples = []
    zero_pts = []
    rng = np.random.default_rng(0)
    for w in windows(r, area):
        a = r.read(b, window=w)
        n += a.size
        bad = np.zeros(a.shape, bool)
        if nodata is not None and not (isinstance(nodata, float) and math.isnan(nodata)):
            nodm = a == np.asarray(nodata).astype(a.dtype)
            n_nod += int(nodm.sum())
            bad |= nodm
        if a.dtype.kind == "f":
            nanm = np.isnan(a)
            n_nan += int(nanm.sum())
            bad |= nanm
        zm = (a == 0) & ~bad
        nz = int(zm.sum())
        n_zero += nz
        valid = ~bad
        nv = int(valid.sum())
        n_valid += nv
        if nv:
            vals = a[valid]
            mn, mx = min(mn, float(vals.min())), max(mx, float(vals.max()))
            vf = vals.astype(np.float64)
            s1 += float(vf.sum())
            s2 += float((vf * vf).sum())
            samples.append(vf[:: max(1, nv // 20000)][:20000].copy())   # copy: a view keeps the block alive
            if is_int and len(hist) <= 300:
                u, c = np.unique(vals, return_counts=True)
                hist.update(dict(zip(u.tolist(), c.tolist())))
        if zero_probe and nz:
            rr, cc = np.nonzero(zm)
            k = rng.choice(rr.size, size=min(200, rr.size), replace=False)
            xs, ys = rasterio.transform.xy(r.transform, rr[k] + w.row_off, cc[k] + w.col_off)
            zero_pts.append(np.column_stack([xs, ys]))
    out = {"n": n, "nodata_value": nodata, "n_nodata": n_nod, "n_nan": n_nan, "n_zero": n_zero, "n_valid": n_valid,
           "min": None if not n_valid else mn, "max": None if not n_valid else mx,
           "mean": s1 / n_valid if n_valid else None,
           "std": math.sqrt(max(0.0, s2 / n_valid - (s1 / n_valid) ** 2)) if n_valid else None}
    if samples:
        s = np.concatenate(samples)
        out["pcts"] = {str(p): float(np.percentile(s, p)) for p in (0.1, 1, 5, 25, 50, 75, 95, 99, 99.9)}
    if is_int and len(hist) <= 300:
        out["histogram"] = {str(k): v for k, v in sorted(hist.items())}
    if zero_pts and poly is not None:
        pts = np.concatenate(zero_pts)
        inside = contains_xy(poly, pts[:, 0], pts[:, 1])
        out["zero_probe"] = {"n_sampled": int(pts.shape[0]), "n_inside_polygon": int(inside.sum()),
                             "inside_examples": pts[inside][:5].tolist()}
    return out


def scan_tif(path):
    t0 = time.time()
    res = {"file": rel(path), "size_bytes": os.path.getsize(path)}
    with rasterio.open(path) as r:
        res.update(driver=r.driver, count=r.count, dtypes=list(r.dtypes), width=r.width, height=r.height,
                   crs=str(r.crs), transform=list(r.transform)[:6], res=list(r.res), bounds=list(r.bounds),
                   nodata=r.nodatavals, compression=str(r.compression), block_shapes=[list(x) for x in r.block_shapes[:1]],
                   descriptions=list(r.descriptions), tags=r.tags(), band_tags=[r.tags(i) for i in range(1, r.count + 1)],
                   overviews=r.overviews(1))
        poly = unary_union(POLY.to_crs(r.crs).geometry.values)
        minx, miny, maxx, maxy = poly.bounds
        is_global = (r.bounds.right - r.bounds.left) > 2 * (maxx - minx)
        res["global_extent"] = bool(is_global)
        no_nodata = all(v is None for v in r.nodatavals)
        res["bands"] = []
        for b in range(1, r.count + 1):
            bs = {"band": b, "all": band_stats(r, b, zero_probe=no_nodata and not is_global, poly=poly)}
            if is_global:
                area = from_bounds(minx, miny, maxx, maxy, r.transform).round_offsets().round_lengths()
                bs["aquifer_bbox"] = band_stats(r, b, area=area, zero_probe=True, poly=poly)
            res["bands"].append(bs)
            print(f"  band {b} done {time.time() - t0:.0f}s", flush=True)
    res["secs"] = round(time.time() - t0, 1)
    return res


def scan_gpkg(path):
    res = {"file": rel(path), "size_bytes": os.path.getsize(path), "layers": []}
    for name, gtype in pyogrio.list_layers(path):
        g = gpd.read_file(path, layer=name)
        info = {"layer": name, "geometry_type": gtype, "crs": str(g.crs), "n_features": len(g),
                "n_invalid": int((~g.geometry.is_valid).sum()), "n_empty": int(g.geometry.is_empty.sum()),
                "total_area_km2": float(g.geometry.area.sum() / 1e6), "bounds": list(g.total_bounds)}
        poly = POLY.to_crs(g.crs)
        info["aquifer_area_km2"] = float(poly.geometry.area.sum() / 1e6)
        u = unary_union(g.geometry.buffer(0).values)
        info["union_area_km2"] = float(u.area / 1e6)
        info["overlap_area_km2"] = info["total_area_km2"] - info["union_area_km2"]
        info["aquifer_not_covered_km2"] = float(unary_union(poly.geometry.values).difference(u).area / 1e6)
        cols = {}
        for c in g.columns:
            if c == "geometry":
                continue
            s = g[c]
            d = {"dtype": str(s.dtype), "n_null": int(s.isna().sum()), "n_unique": int(s.nunique())}
            if s.dtype.kind in "iuf":
                d.update(min=float(s.min()), max=float(s.max()))
            if s.nunique() <= 25:
                d["values"] = {str(k): int(v) for k, v in s.value_counts(dropna=False).items()}
            cols[c] = d
        info["columns"] = cols
        res["layers"].append(info)
    return res


def main():
    only = sys.argv[1:]
    out = []
    for f in raster_files(only):
        print(rel(f), flush=True)
        out.append(scan_gpkg(f) if f.endswith(".gpkg") else scan_tif(f))
        json.dump(out, open(OUT + f"/scan_rasters{'_' + '_'.join(only) if only else ''}.json", "w"), indent=1, default=str)
    print("DONE")


if __name__ == "__main__":
    main()
