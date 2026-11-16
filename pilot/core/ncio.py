"""NetCDF reading and writing.

Reading: files are opened with automatic masking/scaling OFF and decoded explicitly from their own
attributes (scale_factor, add_offset, _FillValue) - vendored from residual_et/sources.py (open_ro,
decode), which was audited on 2026-10-07; tests/test_unit_core.py checks the copies still agree.

Writing: every output is written on scratch, read back and compared with what was computed, then moved
(locally) or copied with size + MD5 check (Drive). A file is complete only when its global attribute
processing_status == "complete" and its provenance signature equals the stage's signature.
"""
import json
import os
import time as _time

import numpy as np

from .errors import SourceError, DataIntegrityError
from . import provenance as P


def open_ro(path, tries=4):
    import netCDF4 as nc
    last = None
    for k in range(tries):
        try:
            ds = nc.Dataset(str(path), "r")
            ds.set_auto_maskandscale(False)
            return ds
        except OSError as e:
            last = e
            _time.sleep(2 * (k + 1))
    raise SourceError(f"cannot open {path}: {last}")


def decode(raw, var):
    """Packed integers -> float32 with NaN at the fill value (same arithmetic as residual_et.sources.decode)."""
    out = raw.astype("float32")
    scale = float(getattr(var, "scale_factor", 1.0))
    offset = float(getattr(var, "add_offset", 0.0))
    if scale != 1.0 or offset != 0.0:
        out = out * np.float32(scale) + np.float32(offset)
    if "_FillValue" in var.ncattrs():
        fill = float(var.getncattr("_FillValue"))
        if not np.isnan(fill):
            out[raw == raw.dtype.type(fill)] = np.nan
    return out


def fill_of(var, default=None):
    return var.getncattr("_FillValue") if "_FillValue" in var.ncattrs() else default


def months_of(ds, name="time"):
    """Calendar month (1-12) and year of every record of a monthly file, from its own time variable."""
    import netCDF4 as nc
    v = ds[name]
    t = nc.num2date(np.asarray(v[:], "float64"), v.units, calendar=getattr(v, "calendar", "standard"),
                    only_use_cftime_datetimes=False)
    return np.array([d.year for d in t]), np.array([d.month for d in t])


def json_attr(path, name):
    with open_ro(path) as ds:
        return json.loads(ds.getncattr(name)) if name in ds.ncattrs() else None


def is_complete(path, sig):
    """True when path is a finished output of the stage with this signature."""
    if not os.path.exists(path):
        return False
    try:
        with open_ro(path, tries=1) as ds:
            if getattr(ds, "processing_status", "") != "complete":
                return False
            prov = json.loads(getattr(ds, "provenance", "{}"))
            return prov.get("signature") == sig
    except Exception:
        return False


def write_nc(path, dims, variables, attrs, provenance, scratch, complevel=2):
    """Write one CF-style NetCDF4 file and deliver it to `path`.

    dims        {name: size}
    variables   {name: dict(dims=tuple, data=array, dtype=str, fill=None|number, attrs=dict)}
    attrs       global attributes (strings / numbers)
    provenance  dict from provenance.record(); stored as JSON in the global attribute 'provenance'
    """
    import netCDF4 as nc
    from . import drive
    os.makedirs(scratch, exist_ok=True)
    tmp = os.path.join(scratch, os.path.basename(path) + f".{os.getpid()}.tmp")
    if os.path.exists(tmp):
        os.remove(tmp)
    with nc.Dataset(tmp, "w", format="NETCDF4") as ds:
        for d, n in dims.items():
            ds.createDimension(d, n)
        for name, v in variables.items():
            dtype = v.get("dtype", str(np.asarray(v["data"]).dtype))
            kw = dict(zlib=complevel > 0, complevel=max(complevel, 1), shuffle=True) if v["dims"] else {}
            var = ds.createVariable(name, dtype, v["dims"], fill_value=v.get("fill"), **kw)
            var.set_auto_maskandscale(False)
            for a, val in (v.get("attrs") or {}).items():
                var.setncattr(a, val)
            data = np.asarray(v["data"])
            if dtype == "str":
                for i, s in enumerate(data):
                    var[i] = str(s)
            else:
                var[...] = data.astype(dtype, copy=False)
        for a, val in attrs.items():
            ds.setncattr(a, val)
        ds.setncattr("Conventions", "CF-1.9")
        ds.setncattr("provenance", P.dumps(provenance))
        ds.setncattr("history", f"{P.utcnow()} created by pilot stage {provenance.get('stage')}")
        ds.setncattr("processing_status", "complete")
    # read back and compare with what was computed
    with open_ro(tmp) as ds:
        for name, v in variables.items():
            if v.get("dtype") == "str":
                continue
            got = np.asarray(ds[name][...])
            want = np.asarray(v["data"]).astype(v.get("dtype", str(np.asarray(v["data"]).dtype)), copy=False)
            if got.shape != want.shape or not np.array_equal(got, want, equal_nan=got.dtype.kind == "f"):
                raise DataIntegrityError(f"{os.path.basename(path)}: variable {name} differs after read-back")
    drive.deliver(tmp, path)
    return path


def read_window(var, rows, cols, fill, t=None):
    """Read var[t?, rows, cols] where rows/cols are (start, stop) in the FILE's pixel coordinates and may
    extend beyond it; the part outside the file is filled with `fill`."""
    ny, nx = var.shape[-2], var.shape[-1]
    r0, r1 = rows
    c0, c1 = cols
    out = np.full((r1 - r0, c1 - c0), fill, dtype=var.dtype)
    a0, a1 = max(r0, 0), min(r1, ny)
    b0, b1 = max(c0, 0), min(c1, nx)
    if a1 > a0 and b1 > b0:
        sl = (slice(a0, a1), slice(b0, b1)) if t is None else (t, slice(a0, a1), slice(b0, b1))
        out[a0 - r0:a1 - r0, b0 - c0:b1 - c0] = np.asarray(var[sl])
    return out
