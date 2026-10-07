"""Full-read integrity scan of every data variable in merged_datasets (read-only).

Reads each variable end to end with auto mask/scale OFF (so raw sentinels and decompression
errors are visible), in blocks of ~BLOCK_MB. Per variable it records value statistics, fill/NaN
counts, sentinel leakage, range violations, and (for time-stacked variables) per-timestep valid
counts, means and CRC32 of each raw slab, plus a per-pixel valid-count map.
"""
import glob
import json
import math
import os
import sys
import time
import traceback
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import netCDF4

ROOT = "G:/MSU_GWB/datasets/merged_datasets"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = HERE + "/results"
os.makedirs(OUT + "/perstep", exist_ok=True)
BLOCK_MB = 24
SENTINELS = [-9999.0, -999.0, -99.0, -32768.0, -32767.0, 32767.0, 65535.0, 65534.0, 255.0, 9.96921e36, 1e20, -1e20]

# physical plausibility bounds (loose, to catch gross errors rather than outliers)
PHYS = {
    "pr": (0, 600), "tmmn": (200, 320), "tmmx": (210, 335), "vpd": (0, 15), "srad": (0, 450), "sph": (0, 0.04),
    "rmax": (0, 100), "rmin": (0, 100), "vs": (0, 40), "th": (0, 360), "eto": (0, 20), "etr": (0, 30),
    "bi": (0, 400), "erc": (0, 250), "fm100": (0, 40), "fm1000": (0, 45),
    "pdsi": (-12, 12), "z": (-12, 12),
    "lwe_thickness": (-100, 100), "uncertainty": (0, 100),
    "et": (0, 1500),                        # SSEBop dekadal mm (raw int16)
    "ET": (0, 1000),                        # ECHO-ET W m-2
    "sm_rootzone": (0, 0.7), "sm_rootzone_pctl": (0, 100), "land_evapotranspiration_flux": (-1e-3, 1e-2),
    "baseflow_flux": (0, 1e-3), "depth_to_water_table_from_surface": (-5, 500),
    "soil_moisture_am": (0, 0.8), "soil_moisture_pm": (0, 0.8), "vod_am": (0, 6), "vod_pm": (0, 6),
    "sif_740": (-3, 8), "sif_757_daily": (-3, 8), "sif_771_daily": (-3, 8), "sif_740_daily": (-3, 8),
    "solar_zenith": (0, 90), "sensor_zenith": (0, 90), "rel_azimuth": (0, 360),
    "sif": (-0.5, 2), "sif_monthly": (-0.5, 2),
    "NDVI": (-2000, 10000), "EVI": (-2000, 10000), "sur_refl_b01": (0, 10000), "sur_refl_b02": (0, 10000),
    "sur_refl_b03": (0, 10000), "sur_refl_b07": (0, 10000), "SummaryQA": (0, 3),
    "b1": (0, 1000),
    "NEE": (-100, 100), "LE": (-300, 1200), "H": (-300, 1200), "Rg": (-50, 1500), "Tair": (-60, 60),
    "Tsoil": (-60, 70), "rH": (0, 105), "VPD": (0, 100), "P": (0, 200), "Reco_uStar": (-5, 100),
    "GPP_uStar_f": (-20, 100), "NEE_uStar_f": (-100, 100),
    "Year": (1990, 2026), "DoY": (1, 366), "Hour": (0, 24),
}


def slabs(shape, itemsize):
    """Yield index tuples that cover the variable in blocks of ~BLOCK_MB; mode tells how it was split."""
    nbytes_lead = int(np.prod(shape[1:], dtype=np.int64)) * itemsize if len(shape) > 1 else itemsize
    budget = BLOCK_MB * 1024 * 1024
    if len(shape) == 1:
        step = max(1, budget // itemsize)
        for i in range(0, shape[0], step):
            yield (slice(i, min(i + step, shape[0])),), "axis0"
    elif nbytes_lead <= budget:
        step = max(1, budget // nbytes_lead)
        for i in range(0, shape[0], step):
            yield (slice(i, min(i + step, shape[0])),), "axis0"
    else:
        row_bytes = nbytes_lead // shape[1]
        step = max(1, budget // row_bytes)
        for i in range(shape[0]):
            for r in range(0, shape[1], step):
                yield (i, slice(r, min(r + step, shape[1]))), "rows"


def scan_var(path, name):
    t0 = time.time()
    ds = netCDF4.Dataset(path)
    v = ds[name]
    v.set_auto_maskandscale(False)
    res = {"file": os.path.basename(path), "var": name, "dtype": str(v.dtype), "shape": list(v.shape), "dims": list(v.dimensions)}
    attrs = {k: v.getncattr(k) for k in v.ncattrs()}
    fill = attrs.get("_FillValue", None)
    mval = attrs.get("missing_value", None)
    vr = attrs.get("valid_range", None)
    vmin = attrs.get("valid_min", None)
    vmax = attrs.get("valid_max", None)
    res["fill"] = None if fill is None else float(fill)
    res["missing_value"] = None if mval is None else float(mval)
    if v.dtype.kind in "SUO":
        try:
            a = v[: min(v.shape[0], 5000)]
            res["char_sample"] = [bytes(r).decode("latin-1").strip() for r in a[:5]] if a.ndim == 2 else None
        except Exception as e:  # noqa: BLE001
            res["error"] = str(e)
        ds.close()
        return res
    kind = v.dtype.kind
    phys = PHYS.get(name)
    n_total = n_fill = n_nan = n_inf = n_valid = 0
    n_sent = {str(s): 0 for s in SENTINELS}
    n_range = n_phys = n_vr = 0
    vmin_seen, vmax_seen = math.inf, -math.inf
    s1 = s2 = 0.0
    samples = []
    errors = []
    perstep_n, perstep_mean, perstep_crc, perstep_min, perstep_max = [], [], [], [], []
    vmap = None
    is_axis0 = None
    for idx, mode in slabs(tuple(v.shape), v.dtype.itemsize):
        try:
            a = v[idx]
        except Exception as e:  # noqa: BLE001
            errors.append({"index": str(idx), "error": f"{type(e).__name__}: {str(e)[:200]}"})
            if mode == "axis0" and isinstance(idx[0], slice):
                for k in range(idx[0].start, idx[0].stop):
                    perstep_n.append(-1)
                    perstep_mean.append(float("nan"))
                    perstep_crc.append(0)
                    perstep_min.append(float("nan"))
                    perstep_max.append(float("nan"))
            continue
        a = np.asarray(a)
        n_total += a.size
        if kind == "f":
            nanm = np.isnan(a)
            infm = np.isinf(a)
            if fill is not None and not np.isnan(fill):
                fillm = a == np.asarray(fill, dtype=a.dtype)
            else:
                fillm = nanm
            n_nan += int(nanm.sum())
            n_inf += int(infm.sum())
            n_fill += int(fillm.sum())
            valid = ~(fillm | nanm | infm)
        else:
            fillm = (a == np.asarray(fill, dtype=a.dtype)) if fill is not None else np.zeros(a.shape, bool)
            n_fill += int(fillm.sum())
            valid = ~fillm
        nv = int(valid.sum())
        n_valid += nv
        if nv:
            vals = a[valid]
            vmin_seen = min(vmin_seen, float(vals.min()))
            vmax_seen = max(vmax_seen, float(vals.max()))
            vf = vals.astype(np.float64)
            s1 += float(vf.sum())
            s2 += float((vf * vf).sum())
            for s in SENTINELS:
                n_sent[str(s)] += int((vals == s).sum())
            if phys is not None:
                n_phys += int(((vals < phys[0]) | (vals > phys[1])).sum())
            if vr is not None:
                n_vr += int(((vals < vr[0]) | (vals > vr[1])).sum())
            if vmin is not None:
                n_range += int((vals < vmin).sum())
            if vmax is not None:
                n_range += int((vals > vmax).sum())
            k = max(1, nv // 50000)
            samples.append(vals[::k][:50000].astype(np.float64))
        if mode == "axis0" and a.ndim >= 3 and isinstance(idx[0], slice):
            is_axis0 = True
            n = a.shape[0]
            vm = valid.reshape(n, -1)
            cnt = vm.sum(1)
            perstep_n.extend(cnt.tolist())
            af = np.where(vm, a.reshape(n, -1).astype(np.float64), np.nan)
            with np.errstate(all="ignore"):
                perstep_mean.extend(np.nanmean(af, axis=1).tolist())
                perstep_min.extend(np.nanmin(af, axis=1).tolist())
                perstep_max.extend(np.nanmax(af, axis=1).tolist())
            perstep_crc.extend(zlib.crc32(np.ascontiguousarray(a[i]).tobytes()) for i in range(n))
            vc = valid.sum(0).astype(np.int32)
            vmap = vc if vmap is None else vmap + vc
    res.update(n_total=n_total, n_fill=n_fill, n_nan=n_nan, n_inf=n_inf, n_valid=n_valid,
               valid_frac=(n_valid / n_total) if n_total else None,
               min=None if n_valid == 0 else vmin_seen, max=None if n_valid == 0 else vmax_seen,
               mean=(s1 / n_valid) if n_valid else None,
               std=(math.sqrt(max(0.0, s2 / n_valid - (s1 / n_valid) ** 2)) if n_valid else None),
               sentinels={k: c for k, c in n_sent.items() if c},
               n_outside_valid_range_attr=n_vr, n_outside_valid_minmax_attr=n_range,
               phys_bounds=phys, n_outside_phys=n_phys, errors=errors)
    if samples:
        s = np.concatenate(samples)
        res["pcts"] = {str(p): float(np.percentile(s, p)) for p in (0.1, 1, 5, 25, 50, 75, 95, 99, 99.9)}
        res["n_unique_sample"] = int(len(np.unique(s)))
    if perstep_n:
        np.savez_compressed(f"{OUT}/perstep/{os.path.basename(path)}__{name}.npz", n=np.array(perstep_n), mean=np.array(perstep_mean),
                            mn=np.array(perstep_min), mx=np.array(perstep_max), crc=np.array(perstep_crc, dtype=np.uint32),
                            vmap=vmap if vmap is not None else np.zeros(1))
    res["secs"] = round(time.time() - t0, 1)
    ds.close()
    return res


def task_list(files):
    tasks = []
    for f in files:
        try:
            ds = netCDF4.Dataset(f)
        except Exception as e:  # noqa: BLE001
            print("OPEN FAIL", f, e)
            continue
        for k, v in ds.variables.items():
            if v.ndim == 1 and k in ds.dimensions and v.dimensions == (k,):
                continue          # coordinate variables handled by coord_checks
            if v.ndim == 0:
                continue
            nbytes = int(np.prod(v.shape, dtype=np.int64)) * v.dtype.itemsize
            tasks.append((nbytes, f, k))
        ds.close()
    tasks.sort(reverse=True)
    return tasks


def main():
    only = sys.argv[1:]
    files = sorted(glob.glob(ROOT + "/*.nc") + glob.glob(ROOT + "/USGS data/derived_usgs/**/*.nc", recursive=True))
    if only:
        files = [f for f in files if any(o in f for o in only)]
    tasks = task_list(files)
    print(len(tasks), "variables to scan", flush=True)
    results = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(scan_var, f, k): (f, k) for _, f, k in tasks}
        for fu in as_completed(futs):
            f, k = futs[fu]
            try:
                r = fu.result()
            except Exception as e:  # noqa: BLE001
                r = {"file": os.path.basename(f), "var": k, "fatal": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()[-600:]}
            results.append(r)
            print(f"[{time.time() - t0:7.0f}s] {r['file']}::{r['var']} "
                  f"valid={r.get('valid_frac')} errs={len(r.get('errors', []))} {r.get('fatal', '')}", flush=True)
            json.dump(results, open(OUT + f"/scan_{'_'.join(only) if only else 'all'}.json", "w"), indent=1, default=str)
    print("DONE", time.time() - t0)


if __name__ == "__main__":
    main()
