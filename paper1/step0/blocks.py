"""S0.1 (part): 30 m Landsat / OpenET monthly values -> 480 m block statistics, one file per tile-year.

Per block and calendar month:
    <IDX>_mean, <IDX>_sd, <IDX>_n   for NIRv, NDMI, NDVI  (valid = not fill and n_obs > 0)
    water_px                        pixels with MNDWI > 0 and AWEIsh > 0 and n_obs >= WATER_MIN_NOBS
    ET_mean, ET_sd, ET_n            OpenET et_ensemble_mad (mm per month, stored integers, fill -32768)
    roi_px                          aquifer pixels of the block (from roi_mask)
P1.0 (plan v3) adds, per block and month:
    <T>_hist (month, bin, by, bx)   counts of the valid 30 m values of ET, NIRv, NDMI in HIST_BINS bins between
                                    HIST_EDGES[T] (values beyond the edges fall in the outer bins), so the
                                    distribution inside a block survives the aggregation
    oli_share                       sum of n_obs_oli / sum of n_obs over observed pixels (Landsat 8/9 share)
    obs_day_mean                    mean obs_day (mean day of the month of the observations) over observed pixels
    ET_models_mean                  mean et_ensemble_mad_count (OpenET models kept after outlier removal)
    ET_min_mean, ET_max_mean        mean of the smallest / largest model value (mm)
    ET_scenes_mean                  mean clear-scene count of the OpenET models in OPENET_COUNT_VARS
Only the variables used are read (8 of the 23 Landsat variables, 4 + len(OPENET_COUNT_VARS) of the OpenET
variables), one month at a time, inside the block window. Values are decoded from each file's own
scale_factor / _FillValue. SD is the population SD of the valid pixels. Months a file does not hold stay
NaN / 0. A support variable missing from a source file stays NaN and is named in the attribute missing_vars.
"""
import os

import numpy as np

from .. import config as C
from ..core import paths, ncio
from ..core.errors import SourceError
from . import tiles as TL

INDICES = ["NIRv", "NDMI", "NDVI"]
HIST_TARGETS = ["ET", "NIRv", "NDMI"]
COUNT_FILL = 255
SUPPORT = ["oli_share", "obs_day_mean", "ET_models_mean", "ET_min_mean", "ET_max_mean", "ET_scenes_mean"]


def block_stats(a, nby, nbx):
    """mean, SD, count of the finite values of a (nby*16, nbx*16) array, per block."""
    b = TL.blockify(a, nby, nbx).astype("float64")
    ok = np.isfinite(b)
    n = ok.sum(-1)
    s = np.where(ok, b, 0.0).sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s / n
        var = np.where(ok, (b - mean[..., None]) ** 2, 0.0).sum(-1) / n
    mean[n == 0] = np.nan
    var[n == 0] = np.nan
    return mean.astype("float32"), np.sqrt(var).astype("float32"), n.astype("uint16")


def hist_edges(name):
    lo, hi = C.HIST_EDGES[name]
    return np.linspace(lo, hi, C.HIST_BINS + 1)


def block_hist(a, nby, nbx, name):
    """Counts of the finite values of a (nby*16, nbx*16) per block in the bins of hist_edges(name):
    (HIST_BINS, nby, nbx) uint16. Values below the first / above the last edge go to the outer bins."""
    b = TL.blockify(a, nby, nbx).astype("float64")
    ok = np.isfinite(b)
    e = hist_edges(name)
    k = np.clip(np.searchsorted(e, np.where(ok, b, e[0]), side="right") - 1, 0, C.HIST_BINS - 1)
    blk = np.broadcast_to(np.arange(nby * nbx).reshape(nby, nbx, 1), b.shape)
    c = np.bincount((blk[ok] * C.HIST_BINS + k[ok]).ravel(), minlength=nby * nbx * C.HIST_BINS)
    return np.moveaxis(c.reshape(nby, nbx, C.HIST_BINS), -1, 0).astype("uint16")


def _counts(ds, name, rows, cols, rec, valid):
    """uint8 count variable as float32 (NaN at its fill or where not valid); None when the file lacks it."""
    if name not in ds.variables:
        return None
    var = ds[name]
    fill = ncio.fill_of(var, COUNT_FILL)
    raw = ncio.read_window(var, rows, cols, fill, t=rec)
    v = raw.astype("float32")
    v[(raw == raw.dtype.type(fill)) | ~valid] = np.nan
    return v


def block_ratio(num, den, nby, nbx):
    """sum(num) / sum(den) per block over pixels where both are finite."""
    a, b = TL.blockify(num, nby, nbx).astype("float64"), TL.blockify(den, nby, nbx).astype("float64")
    ok = np.isfinite(a) & np.isfinite(b)
    sa, sb = np.where(ok, a, 0).sum(-1), np.where(ok, b, 0).sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(sb > 0, sa / sb, np.nan).astype("float32")


def landsat_month(ds, rec, g):
    """Decoded NIRv/NDMI/NDVI (NaN where invalid), the water mask and n_obs for one record."""
    rows, cols = g.frame_rows(), g.frame_cols()
    nobs = ncio.read_window(ds["n_obs"], rows, cols, C.NOBS_FILL, t=rec).astype("int16")
    obs = (nobs > 0) & (nobs != C.NOBS_FILL)
    out = {}
    for name in INDICES + ["MNDWI", "AWEIsh"]:
        var = ds[name]
        raw = ncio.read_window(var, rows, cols, ncio.fill_of(var, C.LANDSAT_FILL), t=rec)
        v = ncio.decode(raw, var)
        v[~obs] = np.nan
        out[name] = v
    with np.errstate(invalid="ignore"):
        water = (out["MNDWI"] > 0) & (out["AWEIsh"] > 0) & (nobs >= C.WATER_MIN_NOBS) & obs
    out["n_obs_oli"] = _counts(ds, "n_obs_oli", rows, cols, rec, obs)
    out["obs_day"] = _counts(ds, "obs_day", rows, cols, rec, obs)
    return out, water, nobs


def openet_month(ds, rec, g, support=False):
    """Decoded ET (mm); with support=True also {name: array or None} of the observation-support variables."""
    import warnings
    rows, cols = g.frame_rows(), g.frame_cols()

    def dec(name):
        if name not in ds.variables:
            return None
        var = ds[name]
        return ncio.decode(ncio.read_window(var, rows, cols, ncio.fill_of(var, C.LANDSAT_FILL), t=rec), var)
    et = dec("et_ensemble_mad")
    if not support:
        return et
    valid = np.isfinite(et)
    sup = dict(ET_min=dec("et_ensemble_mad_min"), ET_max=dec("et_ensemble_mad_max"),
               ET_models=_counts(ds, "et_ensemble_mad_count", rows, cols, rec, valid))
    for k in ("ET_min", "ET_max"):
        if sup[k] is not None:
            sup[k][~valid] = np.nan
    cs = [c for c in (_counts(ds, n, rows, cols, rec, valid) for n in C.OPENET_COUNT_VARS) if c is not None]
    sup["ET_scenes"] = None
    if cs:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)          # pixels where no model counted -> NaN
            sup["ET_scenes"] = np.nanmean(np.stack(cs), 0)
    return et, sup


def build_tile_year(g, year):
    """All block statistics of one tile-year as {name: (dims, array)}; also the input paths used."""
    nb = (g.nby, g.nbx)
    shp = (12,) + nb
    V = {}
    for nm in INDICES + ["ET"]:
        V[f"{nm}_mean"] = np.full(shp, np.nan, "float32")
        V[f"{nm}_sd"] = np.full(shp, np.nan, "float32")
        V[f"{nm}_n"] = np.zeros(shp, "uint16")
    V["water_px"] = np.zeros(shp, "uint16")
    V["nobs_mean"] = np.full(shp, np.nan, "float32")
    for nm in HIST_TARGETS:
        V[f"{nm}_hist"] = np.zeros((12, C.HIST_BINS) + nb, "uint16")
    for k in SUPPORT:
        V[k] = np.full(shp, np.nan, "float32")
    missing = set()
    roi_px = None
    inputs, months_found = {}, {"landsat": [], "openet": []}
    lp = paths.tile_path("landsat", g.tile, year)
    if os.path.exists(lp):
        inputs["landsat"] = lp
        with ncio.open_ro(lp) as ds:
            roi = ncio.read_window(ds["roi_mask"], g.frame_rows(), g.frame_cols(), 0)
            roi_px = TL.blockify((roi == 1).astype("int16"), *nb).sum(-1).astype("uint16")
            yrs, mons = ncio.months_of(ds)
            for rec, (yy, mm) in enumerate(zip(yrs, mons)):
                if yy != year:
                    raise SourceError(f"{lp}: record {rec} is {yy}-{mm:02d}, not in {year}")
                vals, water, nobs = landsat_month(ds, rec, g)
                for nm in INDICES:
                    m, s, n = block_stats(vals[nm], *nb)
                    V[f"{nm}_mean"][mm - 1], V[f"{nm}_sd"][mm - 1], V[f"{nm}_n"][mm - 1] = m, s, n
                    if nm in HIST_TARGETS:
                        V[f"{nm}_hist"][mm - 1] = block_hist(vals[nm], *nb, nm)
                V["water_px"][mm - 1] = TL.blockify(water.astype("int16"), *nb).sum(-1)
                nb_obs = np.where((nobs > 0) & (nobs != C.NOBS_FILL), nobs, np.nan).astype("float32")
                V["nobs_mean"][mm - 1] = block_stats(nb_obs, *nb)[0]
                if vals["n_obs_oli"] is None:
                    missing.add("n_obs_oli")
                else:
                    V["oli_share"][mm - 1] = block_ratio(vals["n_obs_oli"], nb_obs, *nb)
                if vals["obs_day"] is None:
                    missing.add("obs_day")
                else:
                    V["obs_day_mean"][mm - 1] = block_stats(vals["obs_day"], *nb)[0]
                months_found["landsat"].append(int(mm))
    op = paths.tile_path("openet", g.tile, year)
    if os.path.exists(op):
        inputs["openet"] = op
        with ncio.open_ro(op) as ds:
            yrs, mons = ncio.months_of(ds)
            if roi_px is None:
                roi = ncio.read_window(ds["roi_mask"], g.frame_rows(), g.frame_cols(), 0)
                roi_px = TL.blockify((roi == 1).astype("int16"), *nb).sum(-1).astype("uint16")
            for rec, (yy, mm) in enumerate(zip(yrs, mons)):
                if yy != year:
                    raise SourceError(f"{op}: record {rec} is {yy}-{mm:02d}, not in {year}")
                et, sup = openet_month(ds, rec, g, support=True)
                m, s, n = block_stats(et, *nb)
                V["ET_mean"][mm - 1], V["ET_sd"][mm - 1], V["ET_n"][mm - 1] = m, s, n
                V["ET_hist"][mm - 1] = block_hist(et, *nb, "ET")
                for k, a in sup.items():
                    if a is None:
                        missing.add(k)
                    else:
                        V[f"{k}_mean"][mm - 1] = block_stats(a, *nb)[0]
                months_found["openet"].append(int(mm))
    if roi_px is None:
        roi_px = np.zeros(nb, "uint16")
    V["roi_px"] = roi_px
    months_found["missing_vars"] = sorted(missing)
    return V, inputs, months_found


def out_path(ctx, tile, year):
    return f"{ctx.out}/{year}/{tile}/P1_blocks_{tile}_{year}.nc"


def run(ctx):
    """Stage `blocks`."""
    from ..core import provenance as P
    T, G = TL.load_table(ctx.dep("tiles"))
    n_done = n_skip = 0
    missing = []
    for t in sorted(G):
        g = G[t]
        for year in C.agg_years():
            ctx.check_budget()
            p = out_path(ctx, t, year)
            if ncio.is_complete(p, ctx.sig):
                n_skip += 1
                continue
            V, inputs, found = build_tile_year(g, year)
            if not inputs:
                missing.append(f"{t}/{year}")
                ctx.unit(f"{t}/{year}", "missing", note="no Landsat and no OpenET file")
                continue
            dims_of = {4: ("month", "bin", "by", "bx"), 3: ("month", "by", "bx"), 2: ("by", "bx")}
            variables = {k: dict(dims=dims_of[v.ndim], data=v, dtype=str(v.dtype), fill=None) for k, v in V.items()}
            for nm in HIST_TARGETS:
                variables[f"{nm}_hist"]["attrs"] = dict(bin_edges=hist_edges(nm).astype("f8"),
                                                        comment="outer bins also count values beyond the edges")
            variables["month"] = dict(dims=("month",), data=np.arange(1, 13, dtype="int8"), dtype="int8",
                                      attrs=dict(long_name="calendar month"))
            prov = ctx.record({k: P.fingerprint(v) for k, v in inputs.items()},
                              dict(tile=t, year=year, months_found=found, window=g.to_dict()))
            ncio.write_nc(p, dict(month=12, bin=C.HIST_BINS, by=g.nby, bx=g.nbx), variables,
                          dict(tile=t, year=year, title="480 m block statistics of 30 m Landsat indices and OpenET",
                               months_landsat=" ".join(map(str, found["landsat"])),
                               months_openet=" ".join(map(str, found["openet"])),
                               missing_vars=" ".join(found["missing_vars"]) or "none"), prov, ctx.scratch)
            n_done += 1
            ctx.unit(f"{t}/{year}", "done", landsat_months=len(found["landsat"]), openet_months=len(found["openet"]),
                     missing_vars=found["missing_vars"])
            ctx.log.info("%s %d: Landsat months %s, OpenET months %s", t, year, found["landsat"], found["openet"])
    ctx.complete(dict(written=n_done, already_complete=n_skip, missing=missing))


def open_tile_year(stage_dir, tile, year):
    """Decoded block statistics of one tile-year (dict of arrays) or None when absent."""
    p = f"{stage_dir}/{year}/{tile}/P1_blocks_{tile}_{year}.nc"
    if not os.path.exists(p):
        return None
    with ncio.open_ro(p) as ds:
        return {k: np.asarray(ds[k][...]) for k in ds.variables if k != "month"}
