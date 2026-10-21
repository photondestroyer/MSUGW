"""Time-axis and coordinate-grid checks for every .nc (cheap; coordinate variables only).

Time steps are computed numerically from the units (only the first/last values are decoded), so the
6.9 M-row sounding times of the OCO files and the per-record well times are handled too.
"""
import json
import sys
from collections import Counter

import netCDF4
import numpy as np

from audit_paths import OUT, nc_files, rel

UNIT_HOURS = {"second": 1 / 3600, "minute": 1 / 60, "hour": 1, "day": 24}
out = {}
grids = {}


def unit_hours(units):
    u = units.split(" since ")[0].strip().lower().rstrip("s")
    return UNIT_HOURS.get(u)


def decode(v, vals):
    try:
        return [str(x) for x in netCDF4.num2date(vals, v.units, getattr(v, "calendar", "standard"))]
    except Exception as e:  # noqa: BLE001
        return f"decode failed: {e}"


for f in nc_files(sys.argv[1:]):
    name = rel(f)
    rep = {}
    ds = netCDF4.Dataset(f)
    for k, v in ds.variables.items():
        if v.ndim != 1:
            continue
        units = str(getattr(v, "units", ""))
        isdim = k in ds.dimensions and v.dimensions == (k,)
        is_time = " since " in units
        if not (isdim or is_time or k in ("lat", "lon", "x", "y")):
            continue
        if v.dtype == str or v.dtype.kind in "SUO":
            continue
        v.set_auto_maskandscale(False)
        a = np.asarray(v[:])
        fill = getattr(v, "_FillValue", None)
        good = np.ones(a.shape, bool) if fill is None else (a != fill)
        if a.dtype.kind == "f":
            good &= np.isfinite(a)
        r = {"n": int(a.size), "dim": v.dimensions[0], "dtype": str(a.dtype), "n_fill_or_nan": int((~good).sum())}
        g = a[good].astype(np.float64)
        if g.size:
            r["min"], r["max"] = float(g.min()), float(g.max())
        d = np.diff(g)
        if d.size:
            r["monotonic_inc"] = bool((d > 0).all())
            r["monotonic_dec"] = bool((d < 0).all())
            r["n_decreasing_steps"] = int((d < 0).sum())
            r["n_zero_steps"] = int((d == 0).sum())
            r["n_duplicates"] = int(g.size - len(np.unique(g)))
            if k in ("x", "y", "lat", "lon"):
                r["step_counts_top"] = Counter(np.round(d, 9).tolist()).most_common(4)
                r["step_min"], r["step_max"], r["step_mean"] = float(d.min()), float(d.max()), float(d.mean())
                r["step_rel_spread"] = float((d.max() - d.min()) / abs(d.mean())) if d.mean() else None
        if is_time and g.size:
            r["units"] = units
            r["first_last_decoded"] = decode(v, np.array([g[0], g[-1]]))
            r["min_max_decoded"] = decode(v, np.array([g.min(), g.max()]))
            uh = unit_hours(units)
            if uh and g.size > 1 and isdim:
                dh = np.round(d * uh, 3)
                r["step_hours_counts"] = Counter(dh.tolist()).most_common(8)
                med = float(np.median(dh))
                big = np.where(dh > med * 1.5)[0]
                r["median_step_hours"] = med
                r["n_gaps_gt1.5x_median"] = int(big.size)
                big_sorted = big[np.argsort(-dh[big])][:15]
                dec = decode(v, np.concatenate([g[big_sorted], g[big_sorted + 1]])) if big.size else []
                if isinstance(dec, list):
                    nb = len(big_sorted)
                    r["largest_gaps"] = [(dec[i], dec[i + nb], float(dh[j])) for i, j in enumerate(big_sorted)]
        rep[k] = r
        if k in ("x", "y", "lat", "lon"):
            grids.setdefault(name, {})[k] = a
    for gm in ("spatial_ref", "crs"):
        if gm in ds.variables and "GeoTransform" in ds.variables[gm].ncattrs():
            gt = [float(s) for s in ds.variables[gm].getncattr("GeoTransform").split()]
            rep["_geotransform"] = gt
            xs = grids.get(name, {}).get("x", grids.get(name, {}).get("lon"))
            ys = grids.get(name, {}).get("y", grids.get(name, {}).get("lat"))
            if xs is not None and ys is not None:
                rep["_x_vs_gt_maxdev_px"] = float(np.abs(xs - (gt[0] + (np.arange(xs.size) + 0.5) * gt[1])).max() / abs(gt[1]))
                rep["_y_vs_gt_maxdev_px"] = float(np.abs(ys - (gt[3] + (np.arange(ys.size) + 0.5) * gt[5])).max() / abs(gt[5]))
    ds.close()
    out[name] = rep
    print("##", name)
    for k, r in rep.items():
        print("   ", k, r if k.startswith("_") else {kk: vv for kk, vv in r.items() if kk != "dtype"})

print("\n== grid equality among files with identical x/y (or lat/lon) sizes")
by_shape = {}
for name, g in grids.items():
    xs = g.get("x", g.get("lon"))
    ys = g.get("y", g.get("lat"))
    if xs is None or ys is None or xs.ndim != 1:
        continue
    by_shape.setdefault((ys.size, xs.size), []).append((name, xs, ys))
for shp, lst in by_shape.items():
    if len(lst) < 2:
        continue
    ref = lst[0]
    for name, xs, ys in lst[1:]:
        print(shp, name, "vs", ref[0], "| x maxdiff", float(np.abs(xs - ref[1]).max()), "| y maxdiff", float(np.abs(ys - ref[2]).max()))
json.dump(out, open(OUT + "/coord_checks.json", "w"), indent=1, default=str)
