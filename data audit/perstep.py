"""Anomaly list from the per-timestep scan output: empty/partial steps, duplicated slabs, level shifts."""
import glob
import os
import sys

import netCDF4
import numpy as np

from audit_paths import OUT, ROOT

cache = {}


def dates(fn):
    if fn in cache:
        return cache[fn]
    ds = netCDF4.Dataset(f"{ROOT}/{fn}")
    tv = ds["time"] if "time" in ds.variables else ds["day_time"]
    d = netCDF4.num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"))
    cache[fn] = d
    ds.close()
    return d


only = sys.argv[1:]
for p in sorted(glob.glob(OUT + "/perstep/*.npz")):
    base = os.path.basename(p)[:-4]
    fn_tag, var = base.rsplit("__", 1)
    fn = fn_tag.replace("__", "/")
    if only and not any(o in fn for o in only):
        continue
    z = np.load(p)
    n, mean, crc = z["n"], z["mean"], z["crc"]
    if len(n) < 3:
        continue
    try:
        d = dates(fn)
    except Exception:  # noqa: BLE001
        d = None

    def lab(i):
        return str(d[i])[:10] if d is not None and len(d) == len(n) else str(i)

    out = []
    empty = np.where(n == 0)[0]
    if len(empty):
        out.append(f"allEmptySteps={len(empty)} first={lab(empty[0])} last={lab(empty[-1])} e.g.{[lab(i) for i in empty[:6]]}")
    med = np.median(n[n > 0]) if (n > 0).any() else 0
    part = np.where((n > 0) & (n < 0.5 * med))[0]
    if len(part):
        out.append(f"partial(<50%median)={len(part)} e.g.{[lab(i) for i in part[:4]]}")
    nz = n > 0
    idx = np.where(nz)[0]
    vals, cnt = np.unique(crc[nz], return_counts=True)
    dup_vals = vals[cnt > 1]
    if len(dup_vals):
        groups = [[lab(i) for i in idx[crc[nz] == v]] for v in dup_vals]
        out.append(f"duplicateSlabGroups={len(dup_vals)} slabsInGroups={int(cnt[cnt > 1].sum())} e.g.{groups[:6]}")
    consec = int(((crc[1:] == crc[:-1]) & nz[1:]).sum())
    if consec:
        out.append(f"consecutiveIdentical={consec}")
    m = mean[np.isfinite(mean)]
    if len(m) > 10:
        dm = np.diff(mean)
        s = np.nanstd(dm)
        big = np.where(np.abs(dm) > 8 * s)[0] if s > 0 else []
        if len(big) and len(big) < 40:
            out.append(f"meanJumps>8sd={len(big)} e.g.{[lab(i + 1) for i in big[:4]]}")
    vm = z["vmap"]
    if vm.ndim == 2 and vm.size > 1:
        out.append(f"pixelValidMap: allSteps={(vm == vm.max()).mean():.3f} never={(vm == 0).mean():.3f}")
    print(f"{fn[:40]:40s} {var[:30]:30s} steps={len(n)} px/step={int(np.nanmax(n))} " + (" | ".join(out) if out else "ok"))
