"""Bit-for-bit comparison of redundant files with the files that contain them (lazy xarray reads).

Single-year subsets are compared on every record; unfinished copies of the gridMET/drought exports on
every STEP-th record (and the last), which is enough to show they are prefixes of the complete files.
"""
import sys

import numpy as np
import xarray as xr

from audit_paths import ROOT

STEP = 25


def compare(small, big, step=1):
    a = xr.open_dataset(f"{ROOT}/{small}", chunks={}, decode_cf=False)
    b = xr.open_dataset(f"{ROOT}/{big}", chunks={}, decode_cf=False)
    ta = xr.decode_cf(a[["time"]])["time"].values
    tb = xr.decode_cf(b[["time"]])["time"].values
    pos = {t: i for i, t in enumerate(tb)}
    ia = [i for i, t in enumerate(ta) if t in pos][::step]
    ib = [pos[ta[i]] for i in ia]
    print(f"== {small} vs {big}: {len(ta)} records, {sum(t in pos for t in ta)} found in the larger file, {len(ia)} compared")
    for k in ("x", "y", "lat", "lon"):
        if k in a.variables and k in b.variables:
            print(f"   {k} identical: {bool(np.array_equal(a[k].values, b[k].values))}")
    for k, v in a.data_vars.items():
        if v.ndim != 3 or k not in b:
            continue
        x, y = v.isel(time=ia).values, b[k].isel(time=ib).values
        neq = ~((x == y) | (np.isnan(x) & np.isnan(y))) if x.dtype.kind == "f" else (x != y)
        print(f"   {k:16s} differing cells: {int(neq.sum()):>10,}")
    missing = [k for k in a.data_vars if k not in b]
    if missing:
        print("   only in the smaller file:", missing)


if __name__ == "__main__":
    pairs = [("modis_et_v5_dekadal_2020_Ogallala.nc", "MODIS_ET_SSEBop_Merged_Ogallala.nc", 1),
             ("Global_irrigation_Area_2001_Ogallala.nc", "Global_irrigation_Area_Merged_Ogallala.nc", 1),
             ("GRIDMET_Ogallala_lab_drive.nc", "GRIDMET_Ogallala_lab_latest.nc", STEP),
             ("GRIDMET_DROUGHT_Ogallala_lab.nc", "GRIDMET_DROUGHT_Ogallala_lab_colab.nc", STEP)]
    only = sys.argv[1:]
    for s, b, st in pairs:
        if not only or any(o in s for o in only):
            compare(s, b, st)
