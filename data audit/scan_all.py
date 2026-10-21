"""Full-read integrity scan of every data variable in merged_datasets (read-only).

Reads each variable end to end with auto mask/scale OFF (so raw sentinels and decompression
errors are visible), in blocks of ~BLOCK_MB. Per variable it records value statistics, fill/NaN
counts, sentinel leakage, range violations, and (for time-stacked variables) per-timestep valid
counts, means and CRC32 of each raw slab, plus a per-pixel valid-count map.

Packed variables (scale_factor/add_offset) are checked against valid_range in packed units, as CF
defines it, and against the physical bounds below after unpacking. Variables in subgroups are
scanned too (names are '/group/var').
"""
import json
import math
import os
import sys
import time
import traceback
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import dask
import netCDF4
import psutil
import xarray as xr

from audit_paths import OUT, ROOT, get_var, nc_files, rel, walk

os.makedirs(OUT + "/perstep", exist_ok=True)
dask.config.set(scheduler="synchronous")
BLOCK_MB = 16
SENTINELS = [-9999.0, -999.0, -99.0, -32768.0, -32767.0, 32767.0, 65535.0, 65534.0, 255.0, 9.96921e36, 1e20, -1e20, -9e30]

# physical plausibility bounds in unpacked units (loose, to catch gross errors rather than outliers)
PHYS = {
    # gridMET (lab files)
    "pr": (0, 600), "tmmn": (200, 320), "tmmx": (210, 335), "vpd": (0, 15), "srad": (0, 450), "sph": (0, 0.04),
    "rmax": (0, 100), "rmin": (0, 100), "vs": (0, 40), "th": (0, 360), "pet": (0, 20), "etr": (0, 30),
    "bi": (0, 400), "erc": (0, 250), "fm1": (0, 60), "fm10": (0, 60), "fm100": (0, 40), "fm1000": (0, 45),
    # drought indices
    "pdsi": (-12, 12), "z": (-12, 12),
    **{f"{p}{w}": (-5, 5) for p in ("spi", "spei", "eddi") for w in ("14d", "30d", "90d", "180d", "270d", "1y", "2y", "5y")},
    # other gridded products
    "lwe_thickness": (-100, 100), "uncertainty": (0, 100),
    "et": (0, 1500),                        # SSEBop dekadal mm (raw int16)
    "ET": (0, 1000),                        # ECHO-ET W m-2
    "ETdaily": (0, 15), "PTJPLSMinst": (-100, 1000), "STICinst": (-100, 1000), "MOD16inst": (-100, 1000),
    "BESSinst": (-100, 1000), "ETinstUncertainty": (0, 1000),
    "sm_rootzone": (0, 0.7), "sm_rootzone_pctl": (0, 100), "land_evapotranspiration_flux": (-1e-3, 1e-2),
    "baseflow_flux": (0, 1e-3), "depth_to_water_table_from_surface": (-5, 500),
    "soil_moisture_am": (0, 0.8), "soil_moisture_pm": (0, 0.8), "vod_am": (0, 6), "vod_pm": (0, 6),
    "sif": (-0.5, 2), "sif_monthly": (-0.5, 2),
    "NDVI": (-2000, 10000), "EVI": (-2000, 10000), "sur_refl_b01": (0, 10000), "sur_refl_b02": (0, 10000),
    "sur_refl_b03": (0, 10000), "sur_refl_b07": (0, 10000), "SummaryQA": (0, 3),
    "b1": (0, 1000),
    # OCO-2 / OCO-3 Lite SIF (W m-2 sr-1 um-1, degrees, Pa, K)
    "SIF_740nm": (-3, 8), "Daily_SIF_740nm": (-3, 8), "Daily_SIF_757nm": (-3, 8), "Daily_SIF_771nm": (-3, 8),
    "Science_SIF_757nm": (-3, 8), "Science_SIF_771nm": (-3, 8),
    "SZA": (0, 90), "VZA": (0, 90), "Latitude": (31, 44), "Longitude": (-106.5, -95.5),
    "Meteo_surface_pressure": (5e4, 1.1e5), "Meteo_temperature_skin": (220, 345), "Meteo_temperature_two_meter": (220, 330),
    # Landsat 500 m cell means (dimensionless indices; MSI, AWEI and WI2015 are unbounded ratios/sums)
    "EVI2": (-1, 1.5), "SAVI": (-1.5, 1.5), "MSAVI": (-1, 1), "OSAVI": (-1.5, 1.5), "GNDVI": (-1, 1), "NIRv": (-1, 1),
    "kNDVI": (-1, 1), "NDMI": (-1, 1), "NMDI": (-2, 3), "NDWI": (-1, 1), "MNDWI": (-1, 1),
    # wells (ft below land surface)
    "depth_to_water": (-50, 1500), "water_level": (-50, 1500), "well_depth": (0, 3000),
}
# Landsat NDVI/EVI are true ratios in the Landsat file but scaled integers in MOD13A3
PHYS_BY_FILE = {"LANDSAT_C2L2": {"NDVI": (-1, 1), "EVI": (-1, 1.5)}}


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


def phys_bounds(path, name):
    base = name.rsplit("/", 1)[-1]
    for key, table in PHYS_BY_FILE.items():
        if key in os.path.basename(path) and base in table:
            return table[base]
    return PHYS.get(base)


def scan_var(path, name):
    t0 = time.time()
    ds = netCDF4.Dataset(path)
    v = get_var(ds, name)
    v.set_auto_maskandscale(False)
    res = {"file": rel(path), "var": name, "dtype": str(v.dtype), "shape": list(v.shape), "dims": list(v.dimensions)}
    attrs = {k: v.getncattr(k) for k in v.ncattrs()}
    if v.dtype == str or v.dtype.kind in "SUO":
        try:
            a = np.asarray(v[: min(v.shape[0], 200000)], dtype=object).ravel()
            s = [x if isinstance(x, str) else bytes(x).decode("latin-1") for x in a]
            res["n_total"] = int(np.prod(v.shape))
            res["n_read"] = len(s)
            res["n_empty"] = sum(1 for x in s if not x.strip())
            res["n_unique_read"] = len(set(s))
            res["char_sample"] = s[:5]
        except Exception as e:  # noqa: BLE001
            res["error"] = str(e)
        ds.close()
        return res
    ds.close()
    # numeric data: lazy xarray/dask view of the variable, raw (undecoded) values; only one block is
    # ever computed at a time, on dask's synchronous scheduler so memory stays at one block
    gpath, base = (name.rsplit("/", 1) if "/" in name.strip("/") else (None, name))
    xds = xr.open_dataset(path, group=gpath, decode_cf=False, mask_and_scale=False, decode_times=False,
                          chunks={}, engine="netcdf4")
    v = xds[base]
    fill = attrs.get("_FillValue", None)
    mval = attrs.get("missing_value", None)
    vr = attrs.get("valid_range", None)
    vmin = attrs.get("valid_min", None)
    vmax = attrs.get("valid_max", None)
    scale = float(attrs.get("scale_factor", 1.0))
    offset = float(attrs.get("add_offset", 0.0))
    packed = "scale_factor" in attrs or "add_offset" in attrs
    res["fill"] = None if fill is None else float(fill)
    res["missing_value"] = None if mval is None else float(np.ravel(mval)[0])
    res["packed"] = [scale, offset] if packed else None
    kind = v.dtype.kind
    phys = phys_bounds(path, name)
    n_total = n_fill = n_nan = n_inf = n_valid = n_mval = n_huge = n_near9999 = 0
    n_sent = {str(s): 0 for s in SENTINELS}
    n_range = n_phys = n_vr = n_neg_zero_ish = 0
    vmin_seen, vmax_seen = math.inf, -math.inf
    s1 = s2 = 0.0
    samples = []
    errors = []
    perstep_n, perstep_mean, perstep_crc, perstep_min, perstep_max = [], [], [], [], []
    vmap = None
    for idx, mode in slabs(tuple(v.shape), v.dtype.itemsize):
        try:
            a = v[idx].values      # computes this block only
        except MemoryError:
            raise                 # the PC ran short of memory: not a file defect, the variable is retried
        except Exception as e:  # noqa: BLE001
            errors.append({"index": str(idx), "error": f"{type(e).__name__}: {str(e)[:200]}"})
            if mode == "axis0" and isinstance(idx[0], slice):
                for _ in range(idx[0].start, idx[0].stop):
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
            bad = fillm | nanm | infm
        else:
            fillm = (a == np.asarray(fill).astype(a.dtype)) if fill is not None else np.zeros(a.shape, bool)
            n_fill += int(fillm.sum())
            bad = fillm
        if mval is not None:
            mm = np.isin(a, np.ravel(np.asarray(mval)).astype(a.dtype)) & ~bad
            n_mval += int(mm.sum())
            bad = bad | mm
        valid = ~bad
        nv = int(valid.sum())
        n_valid += nv
        if nv:
            vals = a[valid]
            vmin_seen = min(vmin_seen, float(vals.min()))
            vmax_seen = max(vmax_seen, float(vals.max()))
            vf = vals.astype(np.float64)
            if packed:
                vf = vf * scale + offset
            s1 += float(vf.sum())
            s2 += float((vf * vf).sum())
            for s in SENTINELS:
                n_sent[str(s)] += int((vals == s).sum())
            if kind == "f":
                n_huge += int((np.abs(vf) > 1e20).sum())
                n_near9999 += int(((np.abs(np.abs(vf) - 9999) < 1) & (np.abs(vf) != 9999)).sum())
                n_neg_zero_ish += int(((vf < 0) & (vf > -1e-5)).sum())
            if phys is not None:
                n_phys += int(((vf < phys[0]) | (vf > phys[1])).sum())
            if vr is not None:
                n_vr += int(((vals < vr[0]) | (vals > vr[1])).sum())
            if vmin is not None:
                n_range += int((vals < vmin).sum())
            if vmax is not None:
                n_range += int((vals > vmax).sum())
            k = max(1, nv // 50000)
            samples.append(vf[::k][:50000].copy())   # copy: a view would keep the whole block alive
        if mode == "axis0" and a.ndim >= 2 and isinstance(idx[0], slice) and v.dims[0] in ("time", "day"):
            n = a.shape[0]
            vm = valid.reshape(n, -1)
            cnt = vm.sum(1)
            perstep_n.extend(cnt.tolist())
            af = np.where(vm, a.reshape(n, -1).astype(np.float64) * scale + offset, np.nan)
            with np.errstate(all="ignore"):
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    perstep_mean.extend(np.nanmean(af, axis=1).tolist())
                    perstep_min.extend(np.nanmin(af, axis=1).tolist())
                    perstep_max.extend(np.nanmax(af, axis=1).tolist())
            perstep_crc.extend(zlib.crc32(np.ascontiguousarray(a[i]).tobytes()) for i in range(n))
            vc = valid.sum(0).astype(np.int32)
            vmap = vc if vmap is None else vmap + vc
    res.update(n_total=n_total, n_fill=n_fill, n_missing_value=n_mval, n_nan=n_nan, n_inf=n_inf, n_valid=n_valid,
               valid_frac=(n_valid / n_total) if n_total else None,
               raw_min=None if n_valid == 0 else vmin_seen, raw_max=None if n_valid == 0 else vmax_seen,
               min=None if n_valid == 0 else (vmin_seen * scale + offset if packed else vmin_seen),
               max=None if n_valid == 0 else (vmax_seen * scale + offset if packed else vmax_seen),
               mean=(s1 / n_valid) if n_valid else None,
               std=(math.sqrt(max(0.0, s2 / n_valid - (s1 / n_valid) ** 2)) if n_valid else None),
               sentinels={k: c for k, c in n_sent.items() if c},
               n_abs_gt_1e20=n_huge, n_near_9999_not_fill=n_near9999, n_tiny_negative=n_neg_zero_ish,
               n_outside_valid_range_attr=n_vr, n_outside_valid_minmax_attr=n_range,
               phys_bounds=phys, n_outside_phys=n_phys, errors=errors)
    if samples:
        s = np.concatenate(samples)
        res["pcts"] = {str(p): float(np.percentile(s, p)) for p in (0.1, 1, 5, 25, 50, 75, 95, 99, 99.9)}
        res["n_unique_sample"] = int(len(np.unique(s)))
    if perstep_n:
        tag = rel(path).replace("/", "__") + "__" + name.strip("/").replace("/", "~")
        np.savez_compressed(f"{OUT}/perstep/{tag}.npz", n=np.array(perstep_n), mean=np.array(perstep_mean),
                            mn=np.array(perstep_min), mx=np.array(perstep_max), crc=np.array(perstep_crc, dtype=np.uint32),
                            vmap=vmap if vmap is not None else np.zeros(1))
    res["secs"] = round(time.time() - t0, 1)
    res["peak_private_gb"] = round(psutil.Process().memory_info().peak_pagefile / 2**30, 2)
    xds.close()
    return res


def task_list(files):
    tasks = []
    for f in files:
        try:
            ds = netCDF4.Dataset(f)
        except Exception as e:  # noqa: BLE001
            print("OPEN FAIL", f, e)
            continue
        for gpath, g in walk(ds):
            for k, v in g.variables.items():
                if v.ndim == 1 and k in g.dimensions and v.dimensions == (k,):
                    continue          # coordinate variables handled by coord_checks
                if v.ndim == 0:
                    continue
                itemsize = 8 if v.dtype == str else v.dtype.itemsize
                nbytes = int(np.prod(v.shape, dtype=np.int64)) * itemsize
                name = k if gpath == "/" else f"{gpath}/{k}"
                tasks.append((nbytes, f, name))
        ds.close()
    tasks.sort(reverse=True)
    return tasks


def main():
    only = sys.argv[1:]
    tasks = task_list(nc_files(only))
    print(len(tasks), "variables to scan", flush=True)
    results = []
    t0 = time.time()
    out_json = OUT + f"/scan_{os.environ.get('AUDIT_TAG') or ('_'.join(only) if only else 'all')}.json"
    with ProcessPoolExecutor(max_workers=int(os.environ.get("AUDIT_WORKERS", 4))) as ex:
        futs = {ex.submit(scan_var, f, k): (f, k) for _, f, k in tasks}
        for i, fu in enumerate(as_completed(futs)):
            f, k = futs[fu]
            try:
                r = fu.result()
            except Exception as e:  # noqa: BLE001
                r = {"file": rel(f), "var": k, "fatal": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()[-600:]}
            results.append(r)
            print(f"[{time.time() - t0:7.0f}s] {i + 1}/{len(tasks)} {r['file']}::{r['var']} "
                  f"valid={r.get('valid_frac')} errs={len(r.get('errors', []))} {r.get('fatal', '')}", flush=True)
            if i % 20 == 0 or i + 1 == len(tasks):
                json.dump(results, open(out_json, "w"), indent=1, default=str)
    # variables that died for lack of memory are read again one at a time
    retry = [(r["file"], r["var"]) for r in results if "MemoryError" in r.get("fatal", "")]
    if retry:
        print(len(retry), "variables retried one at a time after MemoryError", flush=True)
        results = [r for r in results if "MemoryError" not in r.get("fatal", "")]
        for f, k in retry:
            try:
                r = scan_var(f"{ROOT}/{f}", k)
            except Exception as e:  # noqa: BLE001
                r = {"file": f, "var": k, "fatal": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()[-600:]}
            r["retried"] = True
            results.append(r)
            print(f"[{time.time() - t0:7.0f}s] retry {f}::{k} valid={r.get('valid_frac')} {r.get('fatal', '')}", flush=True)
    json.dump(results, open(out_json, "w"), indent=1, default=str)
    print("DONE", time.time() - t0)


if __name__ == "__main__":
    main()
