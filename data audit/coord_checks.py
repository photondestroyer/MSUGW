"""Time-axis and coordinate-grid checks for every .nc (cheap; coordinate variables only)."""
import glob
import json
import os
from collections import Counter

import netCDF4
import numpy as np

ROOT = "G:/MSU_GWB/datasets/merged_datasets"
files = sorted(glob.glob(ROOT + "/*.nc") + glob.glob(ROOT + "/USGS data/derived_usgs/**/*.nc", recursive=True))
out = {}
grids = {}


def decode(ds, v, vals):
    try:
        return netCDF4.num2date(vals, v.units, getattr(v, "calendar", "standard"), only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    except Exception:  # noqa: BLE001
        try:
            return netCDF4.num2date(vals, v.units, getattr(v, "calendar", "standard"))
        except Exception as e:  # noqa: BLE001
            return f"decode failed: {e}"


for f in files:
    name = os.path.basename(f)
    rep = {}
    ds = netCDF4.Dataset(f)
    for k, v in ds.variables.items():
        if v.ndim != 1:
            continue
        isdim = k in ds.dimensions and v.dimensions == (k,)
        if not (isdim or k in ("time", "lat", "lon", "x", "y")):
            continue
        if v.dtype.kind in "SUO":
            continue
        v.set_auto_maskandscale(False)
        a = v[:]
        r = {"n": int(a.size), "dtype": str(a.dtype), "nan": int(np.isnan(a).sum()) if a.dtype.kind == "f" else 0,
             "min": float(a.min()), "max": float(a.max())}
        d = np.diff(a.astype(np.float64))
        if d.size:
            r["monotonic_inc"] = bool((d > 0).all())
            r["monotonic_dec"] = bool((d < 0).all())
            r["n_nonpositive_steps"] = int((d <= 0).sum()) if not (d < 0).all() else 0
            r["n_duplicates"] = int(a.size - len(np.unique(a)))
            uniq = Counter(np.round(d, 6).tolist())
            r["step_counts_top"] = uniq.most_common(8)
            if k in ("x", "y", "lat", "lon"):
                r["step_min"], r["step_max"], r["step_mean"] = float(d.min()), float(d.max()), float(d.mean())
                r["step_rel_spread"] = float((d.max() - d.min()) / abs(d.mean())) if d.mean() else None
        if k == "time" or "since" in str(getattr(v, "units", "")):
            dt = decode(ds, v, a)
            if not isinstance(dt, str):
                r["first"], r["last"] = str(dt[0]), str(dt[-1])
                if a.size > 1:
                    dd = np.diff(np.array([np.datetime64(x.replace(tzinfo=None) if hasattr(x, "tzinfo") else x) for x in dt]).astype("datetime64[h]")).astype(np.int64)
                    r["step_hours_counts"] = Counter(dd.tolist()).most_common(8)
                    big = np.where(dd > np.median(dd) * 1.5)[0]
                    r["n_gaps_gt1.5x_median"] = int(big.size)
                    r["gaps_sample"] = [(str(dt[i]), str(dt[i + 1]), int(dd[i])) for i in big[:12]]
            else:
                r["decode"] = dt
            r["units"] = getattr(v, "units", None)
        rep[k] = r
        if k in ("x", "y", "lat", "lon") and v.ndim == 1:
            grids.setdefault(name, {})[k] = a
    # compare against GeoTransform attrs of grid_mapping variable
    for gm in ("spatial_ref", "crs"):
        if gm in ds.variables and "GeoTransform" in ds.variables[gm].ncattrs():
            gt = [float(s) for s in ds.variables[gm].getncattr("GeoTransform").split()]
            rep["_geotransform"] = gt
            if "x" in rep and "y" in rep:
                x = grids[name]["x"]
                y = grids[name]["y"]
                px, py = gt[1], gt[5]
                rep["_x_vs_gt_maxdev_px"] = float(np.abs(x - (gt[0] + (np.arange(x.size) + 0.5) * px)).max() / abs(px))
                rep["_y_vs_gt_maxdev_px"] = float(np.abs(y - (gt[3] + (np.arange(y.size) + 0.5) * py)).max() / abs(py))
    ds.close()
    out[name] = rep
    print(name)
    for k, r in rep.items():
        if k.startswith("_"):
            print("   ", k, r)
        else:
            show = {kk: vv for kk, vv in r.items() if kk not in ("dtype",)}
            print("   ", k, show)

# 4 km grid identity across files
print("\n== grid equality (x,y) among the 287x181 files")
ref = None
for name, g in grids.items():
    if "x" in g and "y" in g and g["x"].size == 181 and g["y"].size == 287:
        if ref is None:
            ref = (name, g)
        else:
            print(name, "x equal:", bool(np.array_equal(g["x"], ref[1]["x"])), "maxdiff", float(np.abs(g["x"] - ref[1]["x"]).max()),
                  "| y equal:", bool(np.array_equal(g["y"], ref[1]["y"])), "maxdiff", float(np.abs(g["y"] - ref[1]["y"]).max()), "vs", ref[0])
json.dump(out, open(os.path.dirname(os.path.abspath(__file__)) + "/results/coord_checks.json", "w"), indent=1, default=str)
