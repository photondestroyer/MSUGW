"""S0.3 (the part P1.1/P1.2 use): monthly drivers per lab gridMET cell, rain windows, aridity, bucket.

Source: the Climatology Lab extract (GridMET_Lab_Ogallala.ipynb), daily, packed int16 with per-record
<var>_qc flags. Decoding and the QC rule are those of residual_et/sources.py (LabCube):
    bad day = no_data | (identical & not constant); for srad/tmmn/tmmx/vpd/pet also constant.
Monthly value: sum (pr, pet; scaled to the full month) or mean (vpd, srad, tmin, tmax) of the usable
days, NaN when fewer than DRIVER_MIN_DAY_FRACTION of the days are usable.
Rain windows P<w>d: rain over [month start - w, month start - 1] (ends the day BEFORE the month).
Aridity index: mean annual P / mean annual PET over ARIDITY_BASE_YEARS; deciles over aquifer cells
(lab roi_fraction >= 0.5).
Bucket (one layer, HBV soil routine without snow): S = min(S + P, cap); AET = min(PET * min(1, S/(LP cap)), S);
S -= AET; run daily from BUCKET_SPINUP_FROM for each capacity bin; output = end-of-month S / cap.
Cells computed: those holding an aquifer block or with lab roi_fraction > 0.
SPI90 (plan v3, drought class of P1.2 calibration and gate G1b): SPI_VAR of the lab drought product (same
1/24 degree grid, records every 5-6 days); the monthly value is the record with the latest date inside the
month, NaN when the month has none or that record is flagged. A stratum label only, never a model input.
"""
import calendar

import numpy as np

from .. import config as C
from ..core import paths, ncio
from ..core.errors import SourceError
from . import tiles as TL

UNITS = dict(P="mm month-1", PET="mm month-1", VPD="kPa", SRAD="W m-2", TMIN="degC", TMAX="degC")
SRC = dict(P="pr", PET="pet", VPD="vpd", SRAD="srad", TMIN="tmmn", TMAX="tmmx")


def bad_records(ds, name, i0, i1):
    qn = name + "_qc"
    if qn not in ds.variables:
        return np.zeros(i1 - i0, bool)
    q = np.asarray(ds[qn][i0:i1]).astype(int)
    const, ident, nodata = (q & C.QC_CONSTANT) > 0, (q & C.QC_IDENTICAL) > 0, (q & C.QC_NODATA) > 0
    bad = nodata | (ident & ~const)
    if name in C.QC_CONSTANT_IS_BAD:
        bad |= const
    return bad


def read_days(ds, name, i0, i1, iy, ix):
    """Decoded daily values (days, cells) with lab-flagged days set to NaN."""
    var = ds[name]
    y0, y1, x0, x1 = iy.min(), iy.max() + 1, ix.min(), ix.max() + 1
    raw = np.asarray(var[i0:i1, y0:y1, x0:x1])
    v = ncio.decode(raw, var)[:, iy - y0, ix - x0]
    v[bad_records(ds, name, i0, i1)] = np.nan
    if name in ("tmmn", "tmmx"):
        v = v - np.float32(273.15)
    return v


def month_starts(first_day, last_day):
    """First days of every complete month inside [first_day, last_day]."""
    y, m = first_day.astype(object).year, first_day.astype(object).month
    if first_day.astype(object).day != 1:
        m += 1
        if m == 13:
            y, m = y + 1, 1
    out = []
    while True:
        s = np.datetime64(f"{y:04d}-{m:02d}-01")
        e = s + np.timedelta64(calendar.monthrange(y, m)[1], "D")
        if e - np.timedelta64(1, "D") > last_day:
            break
        out.append(s)
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return np.array(out, dtype="datetime64[D]")


def bucket_step(S, P, E, cap, lp=C.BUCKET_LP):
    S = np.minimum(S + P, cap)
    aet = np.minimum(E * np.minimum(1.0, S / (lp * cap)), S)
    return S - aet, aet


def cells_needed(T, roi_fraction):
    need = np.zeros(roi_fraction.shape, bool)
    need |= roi_fraction > 0
    ok = T["valid"] & (T["lab_iy"] >= 0)
    need[T["lab_iy"][ok], T["lab_ix"][ok]] = True
    iy, ix = np.nonzero(need)
    return iy, ix


def compute(path, iy, ix, log=None):
    """Everything for the cells (iy, ix). Returns dict of arrays."""
    with ncio.open_ro(path) as ds:
        if getattr(ds, "processing_status", "") != "complete":
            raise SourceError(f"{path}: processing_status is not 'complete'")
        import netCDF4 as nc
        tv = ds["time"]
        days = np.asarray(nc.num2date(np.asarray(tv[:], "float64"), tv.units,
                                      calendar=getattr(tv, "calendar", "standard"),
                                      only_use_cftime_datetimes=False), dtype="datetime64[D]")
        step = np.diff(days).astype(int)
        if step.min() != 1 or step.max() != 1:
            raise SourceError(f"{path}: daily axis has gaps or repeats")
        starts = month_starts(days[0], days[-1])
        keep_from = np.datetime64(C.DRIVER_MONTH_FIRST + "-01")
        kept = starts[starts >= keep_from]
        nm, nc_ = kept.size, iy.size
        out = {k: np.full((nm, nc_), np.nan, "float32") for k in SRC}
        for w in C.RAIN_WINDOWS_D:
            out[f"P{w}d"] = np.full((nm, nc_), np.nan, "float32")
        caps = np.asarray(C.BUCKET_CAPACITY_BINS_MM, "float64")[:, None]
        out["BUCKET"] = np.full((nm, caps.size, nc_), np.nan, "float32")
        b0, b1 = C.ARIDITY_BASE_YEARS
        p_base = np.zeros(nc_)
        e_base = np.zeros(nc_)
        n_base = np.zeros(nc_)
        S = np.repeat(0.5 * caps, nc_, axis=1)
        spin = np.datetime64(C.BUCKET_SPINUP_FROM)
        prev_pr = np.full((0, nc_), np.nan, "float32")
        bucket_gaps = 0
        kidx = {d: k for k, d in enumerate(kept)}
        # chunks are whole calendar years, so no month is split between two chunks
        years_ = days.astype("datetime64[Y]")
        bounds = [0] + [int(b) for b in np.nonzero(years_[1:] != years_[:-1])[0] + 1] + [days.size]
        for i, j in zip(bounds[:-1], bounds[1:]):
            D = {k: read_days(ds, SRC[k], i, j, iy, ix) for k in SRC}
            cat = np.concatenate([prev_pr, D["P"]], 0)
            off = prev_pr.shape[0]
            # months that start inside this chunk and end inside it
            ms = starts[(starts >= days[i]) & (starts <= days[j - 1])]
            for s in ms:
                a = int((s - days[0]).astype(int))
                y, m = s.astype(object).year, s.astype(object).month
                nd = calendar.monthrange(y, m)[1]
                if a + nd > j:
                    continue                                   # the file ends inside this month
                la = a - i
                agg = {}
                for k in SRC:
                    v = D[k][la:la + nd]
                    good = np.isfinite(v).sum(0)
                    ok = good >= C.DRIVER_MIN_DAY_FRACTION * nd
                    with np.errstate(invalid="ignore", divide="ignore"):
                        if k in ("P", "PET"):
                            val = np.nansum(v, 0) * nd / good
                        else:
                            val = np.nanmean(v, 0) if v.size else np.full(nc_, np.nan)
                    agg[k] = np.where(ok, val, np.nan)
                if b0 <= y <= b1:
                    ok = np.isfinite(agg["P"]) & np.isfinite(agg["PET"])
                    p_base += np.where(ok, agg["P"], 0)
                    e_base += np.where(ok, agg["PET"], 0)
                    n_base += ok
                if s in kidx:
                    k_ = kidx[s]
                    for k in SRC:
                        out[k][k_] = agg[k]
                    sl = off + la
                    for w in C.RAIN_WINDOWS_D:
                        if sl - w < 0:
                            continue
                        win = cat[sl - w:sl]
                        good = np.isfinite(win).sum(0)
                        with np.errstate(invalid="ignore", divide="ignore"):
                            val = np.nansum(win, 0) * w / good
                        out[f"P{w}d"][k_] = np.where(good >= C.DRIVER_MIN_DAY_FRACTION * w, val, np.nan)
            # bucket, day by day; end-of-month storage recorded for kept months
            for d in range(i, j):
                if days[d] < spin:
                    continue
                Pd, Ed = D["P"][d - i], D["PET"][d - i]
                gap = ~np.isfinite(Pd) | ~np.isfinite(Ed)
                bucket_gaps += int(gap.sum())
                S, _ = bucket_step(S, np.where(np.isfinite(Pd), Pd, 0.0), np.where(np.isfinite(Ed), Ed, 0.0), caps)
                nxt = days[d] + np.timedelta64(1, "D")
                if nxt.astype(object).day == 1:
                    ms_ = np.datetime64(f"{days[d].astype(object).year:04d}-{days[d].astype(object).month:02d}-01")
                    if ms_ in kidx:
                        out["BUCKET"][kidx[ms_]] = (S / caps).astype("float32")
            prev_pr = cat[-365:].copy()
            if log:
                log.info("drivers: days %s .. %s done", days[i], days[j - 1])
    with np.errstate(invalid="ignore", divide="ignore"):
        complete = n_base >= 0.95 * 12 * (b1 - b0 + 1)
        ai = np.where(complete, p_base / e_base, np.nan).astype("float32")
    return dict(months=kept, aridity=ai, bucket_gap_cell_days=bucket_gaps, **out)


def _match(axis, values, path, name):
    """Index of each value on a coordinate axis (exact to 1e-6 degree), or SourceError."""
    values = np.asarray(values, "float64")
    k = np.abs(axis[None, :] - values[:, None]).argmin(1)
    if values.size and np.abs(axis[k] - values).max() > 1e-6:
        raise SourceError(f"{path}: {name} of the gridMET cells are not on its grid (a different grid?)")
    return k


def spi_monthly(path, lat_c, lon_c, months, log=None):
    """SPI of the lab drought file per (month, cell): the last record dated inside the month. Cells are found by
    their centre coordinates, so the file may cover a different box of the same lab grid. months: datetime64[M]."""
    import netCDF4 as nc
    with ncio.open_ro(path) as ds:
        if getattr(ds, "processing_status", "") != "complete":
            raise SourceError(f"{path}: processing_status is not 'complete'")
        if C.SPI_VAR not in ds.variables:
            raise SourceError(f"{path}: no variable {C.SPI_VAR}")
        iy = _match(np.asarray(ds["lat"][:], "float64"), lat_c, path, "latitudes")
        ix = _match(np.asarray(ds["lon"][:], "float64"), lon_c, path, "longitudes")
        tv = ds["time"]
        days = np.asarray(nc.num2date(np.asarray(tv[:], "float64"), tv.units, calendar=getattr(tv, "calendar", "standard"),
                                      only_use_cftime_datetimes=False), dtype="datetime64[D]")
        rec_month = days.astype("datetime64[M]")
        # latest record of each month: records are in time order, so the last index per month wins
        last = {}
        for r, m in enumerate(rec_month):
            last[m] = r
        want = np.array([last.get(m, -1) for m in months])
        out = np.full((months.size, iy.size), np.nan, "float32")
        recs = np.unique(want[want >= 0])
        if recs.size:
            var = ds[C.SPI_VAR]
            y0, y1, x0, x1 = iy.min(), iy.max() + 1, ix.min(), ix.max() + 1
            raw = np.asarray(var[recs, y0:y1, x0:x1])
            v = ncio.decode(raw, var)[:, iy - y0, ix - x0]
            bad = bad_records(ds, C.SPI_VAR, 0, days.size)[recs]
            v[bad] = np.nan
            pos = {r: k for k, r in enumerate(recs)}
            for i, r in enumerate(want):
                if r >= 0:
                    out[i] = v[pos[r]]
    if log:
        log.info("SPI: %d of %d months have a record", int((want >= 0).sum()), months.size)
    return out


def drought_class(spi):
    """0 dry (SPI <= DROUGHT_DRY_SPI), 1 between or unknown, 2 wet (SPI >= DROUGHT_WET_SPI)."""
    s = np.asarray(spi, "float64")
    with np.errstate(invalid="ignore"):
        return np.where(s <= C.DROUGHT_DRY_SPI, 0, np.where(s >= C.DROUGHT_WET_SPI, 2, 1)).astype("int8")


def aridity_deciles(ai, aquifer):
    edges = np.nanquantile(ai[aquifer & np.isfinite(ai)], np.linspace(0, 1, 11)[1:-1])
    dec = np.digitize(ai, edges) + 1
    return np.where(np.isfinite(ai), dec, 0).astype("int8"), edges


def run(ctx):
    """Stage `drivers`."""
    from ..core import provenance as P
    T, _ = TL.load_table(ctx.dep("tiles"))
    src = paths.dataset_path("gridmet")
    path = paths.stage_to_scratch(src, ctx.log)
    with ncio.open_ro(path) as ds:
        roi_fraction = np.asarray(ds["roi_fraction"][:], "float32")
        lat = np.asarray(ds["lat"][:], "float64")
        lon = np.asarray(ds["lon"][:], "float64")
    iy, ix = cells_needed(T, roi_fraction)
    R = compute(path, iy, ix, ctx.log)
    dec, edges = aridity_deciles(R["aridity"], roi_fraction[iy, ix] >= 0.5)
    months = R.pop("months")
    dsrc = paths.dataset_path("drought")
    R["SPI90"] = spi_monthly(dsrc, lat[iy], lon[ix], months.astype("datetime64[M]"), ctx.log)
    nm, ncell = months.size, iy.size
    t = (months - np.datetime64("1970-01-01")).astype("int64").astype("float64")
    variables = dict(
        time=dict(dims=("time",), data=t, dtype="f8", attrs=dict(units="days since 1970-01-01", calendar="standard",
                                                                    long_name="first day of month")),
        cell_iy=dict(dims=("cell",), data=iy.astype("int16"), dtype="i2", attrs=dict(long_name="lab grid row")),
        cell_ix=dict(dims=("cell",), data=ix.astype("int16"), dtype="i2", attrs=dict(long_name="lab grid column")),
        cell_lat=dict(dims=("cell",), data=lat[iy], dtype="f8"), cell_lon=dict(dims=("cell",), data=lon[ix], dtype="f8"),
        roi_fraction=dict(dims=("cell",), data=roi_fraction[iy, ix], dtype="f4"),
        aridity=dict(dims=("cell",), data=R.pop("aridity"), dtype="f4",
                     attrs=dict(long_name=f"P / PET {C.ARIDITY_BASE_YEARS}")),
        aridity_decile=dict(dims=("cell",), data=dec, dtype="i1",
                            attrs=dict(long_name="decile over aquifer cells (1 driest .. 10 wettest); 0 unknown",
                                       edges=" ".join(f"{e:.4f}" for e in edges))),
        capacity_mm=dict(dims=("bin",), data=np.asarray(C.BUCKET_CAPACITY_BINS_MM, "f4"), dtype="f4"),
        lab_shape=dict(dims=("two",), data=np.array(roi_fraction.shape, "i4"), dtype="i4",
                       attrs=dict(long_name="rows, columns of the lab grid")),
    )
    gaps = R.pop("bucket_gap_cell_days")
    for k, v in R.items():
        dims = ("time", "bin", "cell") if k == "BUCKET" else ("time", "cell")
        variables[k] = dict(dims=dims, data=v, dtype="f4",
                            attrs=dict(units=UNITS.get(k, "mm" if k.startswith("P") and k != "SPI90" else "1")))
    variables["SPI90"]["attrs"]["long_name"] = f"{C.SPI_VAR} of the lab drought product, last record of the month"
    prov = ctx.record({"gridmet": P.fingerprint(src), "drought": P.fingerprint(dsrc)},
                      dict(cells=int(ncell), bucket_gap_cell_days=gaps))
    ncio.write_nc(f"{ctx.out}/drivers_monthly.nc", dict(time=nm, cell=ncell, bin=len(C.BUCKET_CAPACITY_BINS_MM), two=2),
                  variables, dict(title="Monthly drivers per lab gridMET cell (paper1 S0.3)"), prov, ctx.scratch)
    ctx.complete(dict(cells=int(ncell), months=int(nm), first=str(months[0]), last=str(months[-1]),
                      aridity_edges=[float(e) for e in edges], bucket_gap_cell_days=gaps))


def load(stage_dir):
    with ncio.open_ro(f"{stage_dir}/drivers_monthly.nc") as ds:
        D = {k: np.asarray(ds[k][...]) for k in ds.variables}
    D["months"] = (np.datetime64("1970-01-01") + D["time"].astype("int64")).astype("datetime64[M]")
    return D


def cell_lookup(D):
    """(lab_iy, lab_ix) -> cell index (-1 where not computed)."""
    lut = np.full(tuple(int(v) for v in D["lab_shape"]), -1, "int32")
    lut[D["cell_iy"], D["cell_ix"]] = np.arange(D["cell_iy"].size)
    return lut
