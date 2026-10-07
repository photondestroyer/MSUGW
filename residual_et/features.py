"""Prepare stage: sources -> feature store (run once, in the system Python).

    python -m residual_et.features --profile full

Stages are cached (rerun skips finished ones; --force redoes them):
    axes -> et -> met -> smap -> drought -> static -> splits -> manifest

What changed against the old panel builders (findings B1, B2, B5, L1, R1, U2):
  * antecedent rainfall is computed from the daily series with its own
    look-back in every chunk, so the result does not depend on chunk size;
  * the soil-moisture lag is a shift of the finished series, not of a batch;
  * all aquifer cells are used, not every second one;
  * ET is the area mean of the 1 km pixels in each cell, with spread and count;
  * anomalies use per-cell baselines that end before the model window.
"""
import argparse
import os
import time

import numpy as np

from . import config as C
from . import climate as K
from . import drought as D
from . import grid as G
from . import landcover as L
from . import sources as S
from . import splits as SP
from .store import Store

DAY = np.timedelta64(1, "D")
MET_VARS = ("pr", "tmmn", "tmmx", "vpd", "srad", "etr")
MET_DEKAD = ("P_dek", "ETR", "VPD", "SRAD", "TMEAN") + tuple(f"P{w}" for w in C.ANTECEDENT_WINDOWS)
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.0f}s] {msg}", flush=True)


# ------------------------------------------------------------------ meteorology
class LabMetReader:
    """Daily gridMET for the aquifer cells, decoded and QC-masked."""

    def __init__(self, grid):
        self.cube = S.LabCube(C.SRC["gridmet"], step_days=1)
        self.cube.require(MET_VARS)
        self.grid = grid
        self.wy, self.wx = grid.window_cells()

    def read(self, t0, t1):
        i0, i1 = self.cube.index(t0, t1)
        out = {}
        for v in MET_VARS:
            a = self.cube.read(v, i0, i1, self.grid.ysl, self.grid.xsl)[:, self.wy, self.wx]
            if np.isfinite(a).any():
                S.check_range(f"gridMET {v}", a, v)
            out[v] = a
        return self.cube.days[i0:i1], out


def build_met(reader, starts, out, chunk=36, on_days=None):
    """Dekad sums/means and antecedent rainfall for `starts`, written into the
    arrays of `out` (keys MET_DEKAD). Each chunk reads its own look-back, so
    the result is identical for every chunk size."""
    starts = np.asarray(starts, "datetime64[D]")
    ends = K.dekad_end(starts)
    pad = max(C.ANTECEDENT_WINDOWS)
    for k0 in range(0, len(starts), chunk):
        k1 = min(k0 + chunk, len(starts))
        days, met = reader.read(starts[k0] - pad * DAY, ends[k1 - 1])
        tmean = (met["tmmn"] + met["tmmx"]) / 2 - 273.15
        s = starts[k0:k1]
        out["P_dek"][k0:k1] = K.aggregate_dekads(days, met["pr"], s, "sum")
        out["ETR"][k0:k1] = K.aggregate_dekads(days, met["etr"], s, "sum")
        out["VPD"][k0:k1] = K.aggregate_dekads(days, met["vpd"], s, "mean")
        out["SRAD"][k0:k1] = K.aggregate_dekads(days, met["srad"], s, "mean")
        out["TMEAN"][k0:k1] = K.aggregate_dekads(days, tmean, s, "mean")
        for w in C.ANTECEDENT_WINDOWS:
            out[f"P{w}"][k0:k1] = K.trailing_sum(days, met["pr"], s, w)
        if on_days is not None:
            inside = (days >= starts[k0]) & (days < ends[k1 - 1])
            on_days(days[inside], met["pr"][inside], tmean[inside], met["etr"][inside])


class DailySink:
    """Collects, in the same pass, the daily forcing of the LSTM and the
    day-of-year normals (baseline years only) that parameterise it."""

    def __init__(self, store, n, profile):
        self.store = store
        self.days = np.arange(np.datetime64(profile.forcing_t0, "D"), np.datetime64(profile.model_t1, "D"))
        self.arr = {k: store.create(f"{k}_daily", (self.days.size, n)) for k in ("P", "T", "EP")}
        self.base = profile.base_met
        self.s = {k: np.zeros((366, n)) for k in ("P", "T", "EP")}
        self.c = {k: np.zeros((366, n)) for k in ("P", "T", "EP")}

    def __call__(self, days, p, t, ep):
        year = days.astype("datetime64[Y]").astype(int) + 1970
        data = dict(P=p, T=t, EP=ep)
        b = (year >= self.base[0]) & (year <= self.base[1])
        if b.any():
            dd = K.doy(days[b])
            for k, x in data.items():
                ok = np.isfinite(x[b])
                np.add.at(self.s[k], dd, np.where(ok, x[b], 0.0))
                np.add.at(self.c[k], dd, ok.astype("float64"))
        f = (days >= self.days[0]) & (days <= self.days[-1])
        if f.any():
            i = (days[f] - self.days[0]).astype(int)
            for k, x in data.items():
                self.arr[k][i] = x[f]

    def finalize(self, win=31):
        """Smoothed day-of-year normals; gaps in the daily forcing take the normal
        of that day (HBV cannot run through a NaN). Returns the number filled."""
        filled = {}
        dd = K.doy(self.days)
        for k in ("P", "T", "EP"):
            with np.errstate(invalid="ignore", divide="ignore"):
                clim = self.s[k] / self.c[k]
            clim[365] = np.where(np.isfinite(clim[365]), clim[365], clim[364])
            h = win // 2
            padded = np.vstack([clim[-h:], clim, clim[:h]])
            cs = np.vstack([np.zeros((1, clim.shape[1])), np.cumsum(padded, axis=0)])
            clim = ((cs[win:] - cs[:-win]) / win).astype("float32")
            self.store.save(f"{k}_dclim", clim)
            a = self.arr[k]
            n = 0
            for i0 in range(0, a.shape[0], 512):
                blk = np.asarray(a[i0:i0 + 512])
                bad = ~np.isfinite(blk)
                if bad.any():
                    blk[bad] = clim[dd[i0:i0 + 512]][bad]
                    a[i0:i0 + 512] = blk
                    n += int(bad.sum())
            a.flush()
            filled[k] = n
        self.store.save("day", self.days)
        return filled


def stage_met(profile, grid, store, starts):
    """Baseline climatologies, anomalies, standardised rainfall, daily forcing."""
    all_starts = K.dekad_starts(f"{profile.base_met[0]}-01-01", profile.model_t1)
    n = grid.n
    scratch = {k: store.create(f"_all_{k}", (all_starts.size, n)) for k in MET_DEKAD}
    reader = LabMetReader(grid)
    S.assert_covers(reader.cube.days, all_starts[0] - max(C.ANTECEDENT_WINDOWS) * DAY,
                    profile.model_t1, "gridMET")
    sink = DailySink(store, n, profile)
    build_met(reader, all_starts, scratch, chunk=36, on_days=sink)
    qc_masked = dict(reader.cube.qc_masked)
    reader.cube.close()
    filled = sink.finalize()
    log(f"met: {all_starts.size} dekads read; lab-QC records masked {qc_masked}; "
        f"daily gaps filled with normals {filled}")

    year = all_starts.astype("datetime64[Y]").astype(int) + 1970
    dek = K.dek36(all_starts)
    base = (year >= profile.base_met[0]) & (year <= profile.base_met[1])
    km = np.searchsorted(all_starts, starts[0])
    assert (all_starts[km:km + starts.size] == starts).all(), "model dekads misaligned"
    mdek = dek[km:km + starts.size]
    clim = {}
    for name, target in (("P_dek", "P_dek_anom"), ("VPD", "VPD_anom"), ("SRAD", "SRAD_anom"),
                         ("TMEAN", "TMEAN_anom"), ("ETR", "ETR_anom")):
        cl = K.DekadClimatology(n)
        a = np.asarray(scratch[name])
        cl.add(dek[base], a[base])
        mean, sd, se, _ = cl.finalize(pool=1)
        clim[name] = mean
        store.save(f"{name}_clim", mean)
        store.save(f"{name}_clim_se", se)
        store.save(target, K.anomaly(a[km:km + starts.size], mdek, mean))
        if name == "P_dek":
            store.save("P_dek", a[km:km + starts.size])
    nyr = profile.base_met[1] - profile.base_met[0] + 1
    for w in C.ANTECEDENT_WINDOWS:
        a = np.asarray(scratch[f"P{w}"])
        samples = a[base].reshape(nyr, 36, n)
        score, se = K.standardize_rainfall(samples, a[km:km + starts.size], mdek, pool=1,
                                           n_boot=profile.n_boot)
        store.save(f"P{w}", a[km:km + starts.size])
        store.save(f"SPI{w}_g", score)
        store.save(f"SPI{w}_se", se if se is not None else np.full(score.shape, np.nan, "float32"))
    p_annual = clim["P_dek"].sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        aridity = np.clip(clim["ETR"].sum(axis=0) / np.maximum(p_annual, 1e-3), 0, C.ARIDITY_CAP)
    store.save("P_annual", p_annual.astype("float32"))
    store.save("aridity", aridity.astype("float32"))
    store.save("T_annual", clim["TMEAN"].mean(axis=0).astype("float32"))
    a = samples = None
    for k in MET_DEKAD:
        del scratch[k]
        try:
            os.remove(store.path / f"_all_{k}.npy")
        except OSError:
            pass                                    # Windows may keep the memmap open; harmless
    return dict(qc_masked=qc_masked, daily_gaps_filled=filled, base_years=list(profile.base_met),
                n_baseline_dekads=int(base.sum()))


# ------------------------------------------------------------------ evapotranspiration
def stage_et(profile, grid, store, starts):
    """SSEBop dekadal ET: area mean, within-cell SD and pixel count per cell;
    climatology from the pre-window baseline; anomaly."""
    r = S.Raster(C.SRC["ssebop"], "et")
    try:
        S.assert_cadence(r.days, 11, "SSEBop")
        dom = (r.days - r.days.astype("datetime64[M]").astype("datetime64[D]")).astype(int)
        if not np.isin(dom, (0, 10, 20)).all():
            raise S.SourceError("SSEBop: a time stamp is not a dekad start (1st, 11th, 21st)")
        pos = np.searchsorted(r.days, starts)
        if pos.max() >= r.days.size or not (r.days[pos] == starts).all():
            raise S.SourceError(f"SSEBop does not hold every dekad of {starts[0]}..{starts[-1]}")
        ysl, xsl = r.window(grid)
        cid = G.cell_ids(r.lat[ysl], r.lon[xsl], grid)
        year = r.days.astype("datetime64[Y]").astype(int) + 1970
        base = np.where((year >= profile.base_et[0]) & (year <= profile.base_et[1]))[0]
        if base.size != 36 * (profile.base_et[1] - profile.base_et[0] + 1):
            raise S.SourceError("SSEBop baseline years are incomplete")
        n = grid.n
        et = np.full((starts.size, n), np.nan, "float32")
        et_sd = np.full((starts.size, n), np.nan, "float32")
        et_n = np.zeros((starts.size, n), "int16")
        clim = K.DekadClimatology(n)
        need = np.union1d(base, pos)
        is_model = np.zeros(r.days.size, bool)
        is_model[pos] = True
        row = {int(p): k for k, p in enumerate(pos)}
        in_base = np.zeros(r.days.size, bool)
        in_base[base] = True
        for a in range(0, need.size, 36):
            idx = need[a:a + 36]
            blk = r.read(int(idx[0]), int(idx[-1]) + 1, ysl, xsl)
            for i in idx:
                mean, sd, cnt = G.block_stats(blk[i - idx[0]], cid, n)
                if in_base[i]:
                    clim.add(K.dek36(r.days[i:i + 1]), mean[None, :])
                if is_model[i]:
                    et[row[int(i)]], et_sd[row[int(i)]], et_n[row[int(i)]] = mean, sd, cnt
        S.check_range("SSEBop ET", et, "et_dekad")
    finally:
        r.close()
    cmean, csd, cse, cyears = clim.finalize(pool=1)
    dek = K.dek36(starts)
    for name, arr in (("ET", et), ("ET_sd", et_sd), ("ET_n", et_n), ("ET_clim", cmean),
                      ("ET_clim_sd", csd), ("ET_clim_se", cse),
                      ("ET_anom", K.anomaly(et, dek, cmean))):
        store.save(name, arr)
    log(f"ET: {starts.size} dekads, median pixels per cell {int(np.median(et_n[et_n > 0]))}, "
        f"baseline years {int(cyears.max())}")
    return dict(median_pixels_per_cell=int(np.median(et_n[et_n > 0])),
                cells_never_observed=int((et_n.max(axis=0) == 0).sum()),
                base_years=list(profile.base_et))


# ------------------------------------------------------------------ soil moisture
def stage_smap(profile, grid, store, starts):
    """SMAP L4 root-zone soil moisture: nearest EASE cell (the product is coarser
    than the grid), daily means, dekad means and the previous-dekad lag."""
    sm = S.SmapL4(C.SRC["smap"])
    try:
        idx, dist = G.nearest_source(sm.lat, sm.lon, grid.cell_lat, grid.cell_lon)
        uniq, inv = np.unique(idx, return_inverse=True)
        days = np.arange(np.datetime64(profile.forcing_t0, "D"), np.datetime64(profile.model_t1, "D"))
        out = store.create("SM_daily", (days.size, grid.n))
        y0 = max(int(str(days[0])[:4]), int(str(sm.days[0])[:4]))
        y1 = int(str(days[-1])[:4])
        for y in range(y0, y1 + 1):
            a = max(np.datetime64(f"{y}-01-01"), days[0], sm.days[0])
            b = min(np.datetime64(f"{y + 1}-01-01"), days[-1] + DAY, sm.days[-1] + DAY)
            if a >= b:
                continue
            d, vals = sm.daily(a, b, uniq)
            out[(d - days[0]).astype(int)] = vals[:, inv]
        out.flush()
    finally:
        sm.close()
    prev = starts[0] - DAY
    prev = K.dekad_starts(str(prev.astype("datetime64[M]")) + "-01", str(starts[0]))[-1]
    ext = np.concatenate([[prev], starts])
    dk = K.aggregate_dekads(days, out, ext, "mean", min_fraction=0.5)
    store.save("SMrz", dk[1:])
    store.save("SMrz_lag1", dk[:-1])                 # a shift of the finished series: no batch seam
    log(f"SMAP: {uniq.size} EASE cells feed {grid.n} cells (max distance {float(dist.max()):.3f} deg); "
        f"dekads without soil moisture: {int(np.isnan(dk[1:]).all(axis=1).sum())}")
    return dict(smap_cells_used=int(uniq.size), max_distance_deg=float(dist.max()),
                dekads_without_sm=int(np.isnan(dk[1:]).all(axis=1).sum()))


# ------------------------------------------------------------------ drought indices
def stage_drought(profile, grid, store, starts):
    """Lab indices at the last pentad strictly before each dekad, plus the
    definition checks (hard checks raise; the rest is reported)."""
    cube = S.LabCube(C.SRC["drought"], step_days=C.PENTAD_MAX_GAP_DAYS)
    checks = {}
    try:
        cube.require(C.LAB_INDICES)
        checks["long_names"] = D.check_long_names(cube)
        checks["usdm_table"] = D.check_usdm_table()
        j = D.lookup_preceding(cube.days, starts)
        lag = (starts - cube.days[j]).astype(int)
        checks["lookup"] = dict(rule="last pentad strictly before the dekad start",
                                min_lag_days=int(lag.min()), max_lag_days=int(lag.max()))
        wy, wx = grid.window_cells()
        idx = {}
        for name in C.LAB_INDICES:
            a = cube.read(name, int(j.min()), int(j.max()) + 1, grid.ysl, grid.xsl)[:, wy, wx]
            S.check_range(f"lab {name}", a, "index")
            idx[name] = a[j - j.min()]
            store.save(name, idx[name])
        checks["qc_masked_records"] = dict(cube.qc_masked)
    finally:
        cube.close()
    checks["index_description"] = {k: D.describe_index(v) for k, v in idx.items()}
    checks["lab_category_vs_usdm"] = D.compare_class_schemes(idx["spi90d"])
    p30 = np.asarray(store.arr("SPI30_g"))
    checks["signs"] = D.check_signs(idx, p30, np.asarray(store.arr("ETR_anom")))
    checks["rank_agreement_with_inhouse"] = dict(
        spi30d_vs_SPI30_g=D.rank_agreement(idx["spi30d"], p30),
        spi90d_vs_SPI90_g=D.rank_agreement(idx["spi90d"], np.asarray(store.arr("SPI90_g"))))
    for k, v in checks["rank_agreement_with_inhouse"].items():
        if np.isfinite(v) and v < 0.7:
            raise S.SourceError(f"definition check failed: {k} rank correlation {v:.2f} < 0.7; "
                                f"the stored index and the rainfall window disagree")
    checks["definitions"] = C.DROUGHT_DEFS
    store.write_json("definition_checks", checks)
    log(f"drought: lookup lag {checks['lookup']['min_lag_days']}-{checks['lookup']['max_lag_days']} d; "
        f"rank agreement {checks['rank_agreement_with_inhouse']}")
    return checks


# ------------------------------------------------------------------ statics
def _wtd(grid):
    r = S.Raster(C.SRC["wtd"], "b1")
    try:
        ysl, xsl = r.window(grid)
        acc, lut = G.Accumulator(grid.n), grid.lut()
        cy, cx = r.var.chunking()[1:]
        for r0 in range((ysl.start // cy) * cy, ysl.stop, cy):       # one stored chunk per read
            a0, a1 = max(r0, ysl.start), min(r0 + cy, ysl.stop)
            for c0 in range((xsl.start // cx) * cx, xsl.stop, cx):
                b0, b1 = max(c0, xsl.start), min(c0 + cx, xsl.stop)
                cid = G.cell_ids(r.lat[a0:a1], r.lon[b0:b1], grid, lut)
                if (cid >= 0).any():
                    acc.add(cid, r.read(0, 1, slice(a0, a1), slice(b0, b1))[0])
    finally:
        r.close()
    return acc.result()


def _paw(grid):
    """Root-zone available water storage (cm) from the 120 m gSSURGO stack (EPSG:5070)."""
    from pyproj import Transformer
    tr = Transformer.from_crs(5070, 4326, always_xy=True)
    ds = S.open_ro(C.SRC["gssurgo"])
    try:
        v = ds["rootznaws"]
        x = np.asarray(ds["x"][:], "float64")
        y = np.asarray(ds["y"][:], "float64")
        acc, lut = G.Accumulator(grid.n), grid.lut()
        cy, cx = v.chunking()
        for r0 in range(0, y.size, cy):
            for c0 in range(0, x.size, cx):
                xx, yy = np.meshgrid(x[c0:c0 + cx], y[r0:r0 + cy])
                lon, lat = tr.transform(xx, yy)
                cid = G.cell_ids_points(lat, lon, grid, lut)
                if (cid >= 0).any():
                    acc.add(cid, S.decode(np.asarray(v[r0:r0 + cy, c0:c0 + cx]), v))
    finally:
        ds.close()
    return acc.result()


def _remap_old_4km(path, var, grid):
    """Nearest cell of an older product on the stretched 287x181 grid (post-hoc screens only)."""
    ds = S.open_ro(path)
    try:
        y = np.asarray(ds["y"][:], "float64")
        x = np.asarray(ds["x"][:], "float64")
        a = S.decode(np.asarray(ds[var][:]), ds[var])
    finally:
        ds.close()
    a = a.reshape(y.size, x.size)
    iy = np.abs(y[:, None] - grid.cell_lat[None, :]).argmin(axis=0)
    ix = np.abs(x[:, None] - grid.cell_lon[None, :]).argmin(axis=0)
    return a[iy, ix].astype("float32")


def stage_static(profile, grid, store, starts):
    out = {}
    out["WTD"], out["WTD_sd"], wtd_n = _wtd(grid)
    log(f"static: WTD done (median {np.nanmedian(out['WTD']):.1f} m, median pixels {int(np.median(wtd_n))})")
    out["PAW"], out["PAW_sd"], paw_n = _paw(grid)
    log(f"static: PAW done (median {np.nanmedian(out['PAW']):.0f} cm)")
    for fn in (L.gfsad_fractions, L.igbp_fractions, L.hsg_fractions, L.irrigation_frequency):
        out.update(fn(grid))
    out["dWL_long"] = _remap_old_4km(C.SRC["dwl_long"], "dWL_predev_to_2019", grid)
    out["dWL_recent"] = _remap_old_4km(C.SRC["dwl_recent"], "dWL_2017_to_2019", grid)
    for k, v in out.items():
        store.save(k, v)
    summary = dict(wtd_median_m=float(np.nanmedian(out["WTD"])), paw_median_cm=float(np.nanmedian(out["PAW"])),
                   mean_f_irrigated=float(np.nanmean(out["f_irrigated"])),
                   mean_f_rainfed_crop=float(np.nanmean(out["f_rainfed_crop"])),
                   mean_f_crop_fragments=float(np.nanmean(out["f_crop_fragments"])),
                   mean_f_natural=float(np.nanmean(out["f_natural"])),
                   mean_gir_high_freq=float(np.nanmean(out["gir_high_freq"])),
                   mean_gir_low_freq=float(np.nanmean(out["gir_low_freq"])),
                   igbp_years=[int(y) for y in out["igbp_years"]],
                   cells_without_wtd=int(np.isnan(out["WTD"]).sum()),
                   cells_without_paw=int(np.isnan(out["PAW"]).sum()))
    log(f"static: GFSAD shares irrigated {summary['mean_f_irrigated']:.2f}, rainfed crop "
        f"{summary['mean_f_rainfed_crop']:.2f}, fragments {summary['mean_f_crop_fragments']:.2f}; "
        f"MCD12Q1 natural {summary['mean_f_natural']:.2f}")
    return summary


# ------------------------------------------------------------------ splits
def stage_splits(profile, grid, store, starts):
    spi = np.asarray(store.arr("spi90d"))
    dryness = np.nanmean(spi, axis=1)                    # aquifer-mean SPI-90 per dekad
    embargo, acf = SP.embargo_from_acf(np.asarray(store.arr("ET_anom")))
    folds = SP.make_folds(starts, dryness, embargo)
    SP.check(folds)
    np.savez(store.path / "folds.npz", split=folds["split"], block=folds["block"],
             block_class=folds["block_class"], block_test_fold=folds["block_test_fold"],
             block_season=folds["block_season"], block_dryness=folds["block_dryness"],
             block_label=np.array(folds["block_label"]))
    summary = dict(embargo_dekads=embargo, et_anomaly_acf_by_lag=acf, folds=SP.summarize(folds),
                   blocks=[dict(label=l, cls=int(c), dryness=float(d), test_fold=int(f))
                           for l, c, d, f in zip(folds["block_label"], folds["block_class"],
                                                 folds["block_dryness"], folds["block_test_fold"])],
                   validation_halves=folds["val_halves"])
    store.write_json("splits", summary)
    s0 = summary["folds"][0]
    log(f"splits: embargo {embargo} dekads; fold 0 train {s0['train']:.2f} (+{s0['embargo']:.2f} "
        f"embargoed) / test {s0['test']:.2f} / val {s0['val']:.2f}")
    return summary


# ------------------------------------------------------------------ driver
def stage_axes(profile, grid, store, starts):
    store.save("dekad_start", starts)
    store.save("dek36", K.dek36(starts))
    store.save("cell_lat", grid.cell_lat)
    store.save("cell_lon", grid.cell_lon)
    store.save("cell_iy", grid.iy)
    store.save("cell_ix", grid.ix)
    store.save("space_block", grid.space_blocks())
    return dict(n_cells=grid.n, n_dekads=int(starts.size), grid_shape=list(grid.frac.shape),
                first_dekad=str(starts[0]), last_dekad=str(starts[-1]))


STAGES = (("axes", stage_axes), ("et", stage_et), ("met", stage_met), ("smap", stage_smap),
          ("drought", stage_drought), ("static", stage_static), ("splits", stage_splits))


def prepare(profile_name, force=False, only=None):
    profile = C.get_profile(profile_name)
    store = Store(profile.store)
    grid = G.load_grid(profile)
    starts = K.dekad_starts(profile.model_t0, profile.model_t1)
    log(f"profile {profile.name}: {grid.n} aquifer cells, {starts.size} dekads "
        f"({starts[0]}..{starts[-1]}), store {store.path}")
    report = store.read_json("manifest") if (store.path / "manifest.json").exists() else {}
    for name, fn in STAGES:
        if only and name != only:
            continue
        if store.done(name) and not force:
            log(f"stage {name}: cached")
            continue
        report[name] = fn(profile, grid, store, starts)
        store.mark(name)
        store.write_json("manifest", report)
    if not only:
        store.require_complete()
        report["sources"] = {k: dict(file=p.name, bytes=p.stat().st_size,
                                     modified=time.strftime("%Y-%m-%d %H:%M", time.localtime(p.stat().st_mtime)))
                             for k, p in C.SRC.items() if p.exists()}
        report["profile"] = dict(name=profile.name, bbox=profile.bbox, model=[profile.model_t0, profile.model_t1],
                                 forcing_t0=profile.forcing_t0, base_met=list(profile.base_met),
                                 base_et=list(profile.base_et))
        store.write_json("manifest", report)
        log("store complete")
    return store


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build the residual-ET feature store.")
    ap.add_argument("--profile", default="full", choices=sorted(C.PROFILES))
    ap.add_argument("--force", action="store_true", help="redo finished stages")
    ap.add_argument("--stage", default=None, choices=[n for n, _ in STAGES])
    a = ap.parse_args()
    prepare(a.profile, a.force, a.stage)
