"""Read-only source adapters with integrity guards.

Every reader decodes explicitly (scale, offset, fill), checks the time axis and
fails loudly instead of returning something plausible. The old scripts opened
gridMET with mask_and_scale=False and looked up the "nearest" drought step;
on the present files that yields raw packed integers and, on the old drought
file, values from the wrong year.
"""
import time as _time

import numpy as np
import netCDF4 as nc

from . import config as C


class SourceError(RuntimeError):
    """A source file does not satisfy what the pipeline assumes about it."""


def open_ro(path, tries=4):
    """Open read-only with retries (drive hiccups have been seen on G:)."""
    last = None
    for k in range(tries):
        try:
            ds = nc.Dataset(str(path), "r")
            ds.set_auto_maskandscale(False)     # decoding is done explicitly below
            return ds
        except OSError as e:
            last = e
            _time.sleep(2 * (k + 1))
    raise SourceError(f"cannot open {path}: {last}")


def nc_days(ds, name="time"):
    """Time coordinate as datetime64[D] (and the sub-daily datetime64[s])."""
    v = ds[name]
    cal = getattr(v, "calendar", "standard")
    t = nc.num2date(np.asarray(v[:], dtype="float64"), v.units, calendar=cal,
                    only_use_cftime_datetimes=False)
    t = np.asarray(t, dtype="datetime64[s]")
    return t.astype("datetime64[D]"), t


def assert_cadence(days, max_step_days, what):
    """Raise if the axis is not strictly increasing or has a gap."""
    step = np.diff(days).astype(int)
    if step.size and step.min() <= 0:
        raise SourceError(f"{what}: time axis is not strictly increasing")
    if step.size and step.max() > max_step_days:
        i = int(step.argmax())
        raise SourceError(f"{what}: gap of {int(step[i])} days between {days[i]} and "
                          f"{days[i + 1]} (allowed {max_step_days})")


def assert_covers(days, t0, t1, what):
    t0, t1 = np.datetime64(t0, "D"), np.datetime64(t1, "D")
    if days[0] > t0 or days[-1] < t1 - np.timedelta64(1, "D"):
        raise SourceError(f"{what}: covers {days[0]}..{days[-1]}, need {t0}..{t1} (exclusive)")


def decode(raw, var):
    """Packed integers -> float32 with NaN at the fill value."""
    out = raw.astype("float32")
    scale = float(getattr(var, "scale_factor", 1.0))
    offset = float(getattr(var, "add_offset", 0.0))
    if scale != 1.0 or offset != 0.0:
        out = out * np.float32(scale) + np.float32(offset)
    if "_FillValue" in var.ncattrs():
        fill = float(var.getncattr("_FillValue"))
        if not np.isnan(fill):                  # a NaN fill is already NaN
            out[raw == raw.dtype.type(fill)] = np.nan
    return out


class LabCube:
    """A Climatology Lab extract written by GridMET_Lab_Ogallala.ipynb
    (dims time/lat/lon, packed int16, per-record <var>_qc flags)."""

    def __init__(self, path, step_days):
        self.path = path
        self.ds = open_ro(path)
        try:
            for d in ("time", "lat", "lon"):
                if d not in self.ds.dimensions:
                    raise SourceError(f"{path.name}: dimension {d!r} missing "
                                      f"(found {list(self.ds.dimensions)})")
            status = getattr(self.ds, "processing_status", "unknown")
            if status != "complete":
                raise SourceError(f"{path.name}: processing_status is {status!r}, not 'complete'")
            self.days, _ = nc_days(self.ds)
            assert_cadence(self.days, step_days, path.name)
            self.lat = np.asarray(self.ds["lat"][:], dtype="float64")
            self.lon = np.asarray(self.ds["lon"][:], dtype="float64")
            if not (self.lat[0] > self.lat[-1] and self.lon[0] < self.lon[-1]):
                raise SourceError(f"{path.name}: expected lat descending and lon ascending")
        except Exception:
            self.ds.close()
            raise
        self.qc_masked = {}                     # var -> number of records masked so far

    def close(self):
        self.ds.close()

    def require(self, names):
        missing = [n for n in names if n not in self.ds.variables]
        if missing:
            raise SourceError(f"{self.path.name}: variables missing: {missing}")

    def index(self, t0, t1):
        """Record slice for [t0, t1) (dates)."""
        i0 = int(np.searchsorted(self.days, np.datetime64(t0, "D")))
        i1 = int(np.searchsorted(self.days, np.datetime64(t1, "D")))
        return i0, i1

    def bad_records(self, name, i0, i1):
        """Records the lab's own QC marks as unusable for this variable."""
        qn = name + "_qc"
        if qn not in self.ds.variables:
            return np.zeros(i1 - i0, bool)
        q = np.asarray(self.ds[qn][i0:i1]).astype(int)
        const, ident, nodata = (q & C.QC_CONSTANT) > 0, (q & C.QC_IDENTICAL) > 0, (q & C.QC_NODATA) > 0
        bad = nodata | (ident & ~const)          # a repeated non-trivial field is a duplicate
        if name in C.QC_CONSTANT_IS_BAD:
            bad |= const
        return bad

    def read(self, name, i0, i1, ysl=slice(None), xsl=slice(None), qc=True):
        """Decoded float32 (T, h, w). Lab-flagged records become NaN when qc=True."""
        self.require([name])
        var = self.ds[name]
        out = decode(np.asarray(var[i0:i1, ysl, xsl]), var)
        if qc:
            bad = self.bad_records(name, i0, i1)
            if bad.any():
                out[bad] = np.nan
                self.qc_masked[name] = self.qc_masked.get(name, 0) + int(bad.sum())
        return out

    def long_name(self, name):
        return str(getattr(self.ds[name], "long_name", ""))


class Raster:
    """An Earth Engine export: dims (time, y, x), y descending, x ascending."""

    def __init__(self, path, var):
        self.path, self.name = path, var
        self.ds = open_ro(path)
        if var not in self.ds.variables:
            raise SourceError(f"{path.name}: variable {var!r} missing")
        self.var = self.ds[var]
        self.lat = np.asarray(self.ds["y"][:], dtype="float64")
        self.lon = np.asarray(self.ds["x"][:], dtype="float64")
        if not (self.lat[0] > self.lat[-1] and self.lon[0] < self.lon[-1]):
            raise SourceError(f"{path.name}: expected y descending and x ascending")
        self.days, _ = nc_days(self.ds)

    def close(self):
        self.ds.close()

    def window(self, grid, pad=1):
        """Row/column slices of this raster that cover the analysis window."""
        le, lo = grid.lat_edges, grid.lon_edges
        rows = np.where((self.lat <= le[grid.ysl.start]) & (self.lat >= le[grid.ysl.stop]))[0]
        cols = np.where((self.lon >= lo[grid.xsl.start]) & (self.lon <= lo[grid.xsl.stop]))[0]
        if rows.size == 0 or cols.size == 0:
            raise SourceError(f"{self.path.name}: does not overlap the analysis window")
        r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + 1 + pad, self.lat.size)
        c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + 1 + pad, self.lon.size)
        return slice(int(r0), int(r1)), slice(int(c0), int(c1))

    def read(self, i0, i1, ysl, xsl):
        """Decoded float32 (T, h, w) with NaN at fill."""
        return decode(np.asarray(self.var[i0:i1, ysl, xsl]), self.var)

    def read_classes(self, i, ysl, xsl):
        """Integer class codes for one time step; fill -> -1."""
        raw = np.asarray(self.var[i, ysl, xsl])
        out = raw.astype("int32")
        if "_FillValue" in self.var.ncattrs():
            out[raw == raw.dtype.type(self.var.getncattr("_FillValue"))] = -1
        return out


class SmapL4:
    """SPL4SMGP subset: 3-hourly sm_rootzone on a 2-D EASE grid."""

    def __init__(self, path):
        self.path = path
        self.ds = open_ro(path)
        self.var = self.ds["sm_rootzone"]
        self.lat = np.asarray(self.ds["lat"][:], dtype="float64")
        self.lon = np.asarray(self.ds["lon"][:], dtype="float64")
        self.days, self.stamps = nc_days(self.ds)

    def close(self):
        self.ds.close()

    def daily(self, t0, t1, flat_cells, slab=1024):
        """Daily means (ndays, ncell) for the SMAP cells in `flat_cells`
        (indices into the flattened EASE grid). Days without data are NaN."""
        t0, t1 = np.datetime64(t0, "D"), np.datetime64(t1, "D")
        days = np.arange(t0, t1)
        out_sum = np.zeros((days.size, flat_cells.size), "float64")
        out_n = np.zeros((days.size, flat_cells.size), "int32")
        i0 = int(np.searchsorted(self.days, t0))
        i1 = int(np.searchsorted(self.days, t1))
        for s in range(i0, i1, slab):
            e = min(s + slab, i1)
            a = decode(np.asarray(self.var[s:e]), self.var).reshape(e - s, -1)[:, flat_cells]
            a[(a < C.RANGES["sm"][0]) | (a > C.RANGES["sm"][1])] = np.nan
            d = (self.days[s:e] - t0).astype(int)
            ok = np.isfinite(a)
            np.add.at(out_sum, d, np.where(ok, a, 0.0))
            np.add.at(out_n, d, ok.astype("int32"))
        with np.errstate(invalid="ignore", divide="ignore"):
            out = np.where(out_n > 0, out_sum / np.maximum(out_n, 1), np.nan).astype("float32")
        return days, out


def check_range(name, arr, key):
    """Raise if finite values leave the plausible range; return (min, max)."""
    lo, hi = C.RANGES[key]
    a = arr[np.isfinite(arr)]
    if a.size == 0:
        raise SourceError(f"{name}: no finite values")
    mn, mx = float(a.min()), float(a.max())
    if mn < lo - 1e-6 or mx > hi + 1e-6:
        raise SourceError(f"{name}: values [{mn:.3f}, {mx:.3f}] outside plausible [{lo}, {hi}]")
    return mn, mx
