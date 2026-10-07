"""Tests on synthetic data (system Python: numpy, scipy, netCDF4, xgboost).
Each test names the finding it guards against."""
import tempfile
from pathlib import Path

import numpy as np

from .. import config as C
from .. import climate as K
from .. import drought as D
from .. import features as F
from .. import grid as G
from .. import landcover as L
from .. import sources as S
from .. import splits as SP
from .. import uq as U

DAY = np.timedelta64(1, "D")


# ------------------------------------------------------------------ helpers
class FakeMet:
    """Deterministic daily 'gridMET' held in memory."""

    def __init__(self, t0="2014-06-01", t1="2018-01-01", n=5, seed=1):
        rng = np.random.default_rng(seed)
        self.days = np.arange(np.datetime64(t0), np.datetime64(t1))
        m = self.days.size
        self.data = dict(pr=rng.gamma(0.3, 6.0, (m, n)).astype("float32"),
                         tmmn=(280 + rng.normal(0, 5, (m, n))).astype("float32"),
                         tmmx=(292 + rng.normal(0, 5, (m, n))).astype("float32"),
                         vpd=rng.uniform(0, 3, (m, n)).astype("float32"),
                         srad=rng.uniform(50, 350, (m, n)).astype("float32"),
                         etr=rng.uniform(0, 12, (m, n)).astype("float32"))

    def read(self, t0, t1):
        i0, i1 = np.searchsorted(self.days, t0), np.searchsorted(self.days, t1)
        return self.days[i0:i1], {k: v[i0:i1] for k, v in self.data.items()}


def _build(reader, starts, chunk):
    out = {k: np.full((starts.size, 5), np.nan, "float32") for k in F.MET_DEKAD}
    F.build_met(reader, starts, out, chunk=chunk)
    return out


# ------------------------------------------------------------------ batch processing (B1, B2)
def test_features_do_not_depend_on_chunk_size():
    reader = FakeMet()
    starts = K.dekad_starts("2015-01-01", "2017-12-01")
    ref = _build(reader, starts, 1)
    for chunk in (7, 32, 1000):
        got = _build(reader, starts, chunk)
        for k in F.MET_DEKAD:
            assert np.allclose(ref[k], got[k], equal_nan=True), f"{k} differs for chunk size {chunk}"
    assert np.isfinite(ref["P90"]).all(), "antecedent rainfall has gaps although the data are complete"


def test_antecedent_rainfall_matches_brute_force():
    reader = FakeMet()
    starts = K.dekad_starts("2015-01-01", "2016-01-01")
    got = _build(reader, starts, 32)
    for w in C.ANTECEDENT_WINDOWS:
        for k, s in enumerate(starts):
            i = np.searchsorted(reader.days, s)
            assert np.allclose(got[f"P{w}"][k], reader.data["pr"][i - w:i].sum(axis=0), rtol=1e-4)


def test_old_batch_code_was_truncated():
    """Documents B1: the old cumulative sum zeroed the look-back days."""
    reader = FakeMet()
    starts = K.dekad_starts("2015-01-01", "2015-12-01")
    days, met = reader.read(starts[0] - 95 * DAY, K.dekad_end(starts)[-1])
    inb = days >= starts[0]
    cs = np.cumsum(np.where(inb[:, None], met["pr"], 0.0), axis=0)        # what the old script did
    di = int((starts[0] - days[0]) / DAY)
    old_p90_first = cs[di] - cs[di - 90]
    true_p90 = met["pr"][di - 89:di + 1].sum(axis=0)
    assert old_p90_first.sum() < 0.3 * true_p90.sum()


# ------------------------------------------------------------------ leakage in time (L1)
def test_rain_enters_and_leaves_the_window_on_the_right_days():
    days = np.arange(np.datetime64("2016-01-01"), np.datetime64("2016-12-31"))
    rain = np.zeros((days.size, 1), "float32")
    spike = np.datetime64("2016-06-15")
    rain[days == spike] = 50.0
    targets = np.arange(np.datetime64("2016-06-01"), np.datetime64("2016-10-15"))
    for w in (30, 90):
        got = K.trailing_sum(days, rain, targets, w)[:, 0] > 0
        want = (targets > spike) & (targets <= spike + w * DAY)          # window = [t-w, t-1]
        assert (got == want).all(), f"{w}-day window is misaligned"
    assert K.trailing_sum(days, rain, np.array([spike]), 30)[0, 0] == 0, "the target day itself leaked in"


def test_index_lookup_is_strictly_before_and_refuses_gaps():
    idx = np.arange(np.datetime64("2019-01-05"), np.datetime64("2022-01-01"), 5)
    starts = K.dekad_starts("2019-02-01", "2021-12-01")
    j = D.lookup_preceding(idx, starts)
    lag = (starts - idx[j]).astype(int)
    assert lag.min() >= 1 and lag.max() <= 5
    same_day = np.array([idx[10]])
    assert idx[D.lookup_preceding(idx, same_day)[0]] == idx[9], "an index dated on the start day was used"
    holed = idx[(idx < np.datetime64("2020-01-01")) | (idx > np.datetime64("2020-12-31"))]
    try:
        D.lookup_preceding(holed, starts)
    except S.SourceError as e:
        assert "2020" in str(e)
    else:
        raise AssertionError("a year without drought index was accepted (old nearest-step behaviour)")


# ------------------------------------------------------------------ source decoding (B5)
def _lab_file(path, days, status="complete"):
    import netCDF4 as nc
    ds = nc.Dataset(path, "w")
    ds.processing_status = status
    for d, n in (("time", len(days)), ("lat", 3), ("lon", 4)):
        ds.createDimension(d, n)
    t = ds.createVariable("time", "f8", ("time",))
    t.units, t.calendar = "days since 1900-01-01 00:00:00", "standard"
    t[:] = (np.asarray(days, "datetime64[D]") - np.datetime64("1900-01-01")).astype(int)
    ds.createVariable("lat", "f8", ("lat",))[:] = [40.0, 39.9, 39.8]
    ds.createVariable("lon", "f8", ("lon",))[:] = [-101.0, -100.9, -100.8, -100.7]
    rng = np.random.default_rng(0)
    raw = {}
    for name, scale, off in (("pr", 0.1, 0.0), ("srad", 0.1, 0.0), ("tmmn", 0.1, 210.0)):
        v = ds.createVariable(name, "i2", ("time", "lat", "lon"), fill_value=32767)
        v.scale_factor, v.add_offset = scale, off
        v.set_auto_maskandscale(False)
        a = rng.integers(1, 900, (len(days), 3, 4)).astype("i2")
        if name == "pr":
            a[2] = 0                       # a dry day: constant field, natural
            a[5] = a[4]                    # a repeated wet field: a duplicate
        if name == "srad":
            a[3] = 0                       # constant radiation: a placeholder
        a[0, 0, 0] = 32767
        v[:] = a
        raw[name] = a
        q = ds.createVariable(name + "_qc", "u1", ("time",))
        flags = np.zeros(len(days), "u1")
        if name == "pr":
            flags[2], flags[5] = C.QC_CONSTANT, C.QC_IDENTICAL
        if name == "srad":
            flags[3] = C.QC_CONSTANT
        q[:] = flags
    ds.close()
    return raw


def test_packed_integers_are_decoded_and_lab_flags_respected():
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        p = Path(tmp) / "lab.nc"
        days = np.arange(np.datetime64("2015-01-01"), np.datetime64("2015-01-11"))
        raw = _lab_file(p, days)
        cube = S.LabCube(p, step_days=1)
        tm = cube.read("tmmn", 0, 10)
        assert np.isnan(tm[0, 0, 0]) and np.allclose(tm[1], raw["tmmn"][1] * 0.1 + 210.0, atol=1e-3)
        pr = cube.read("pr", 0, 10)
        assert np.isfinite(pr[2]).all(), "a dry day was thrown away"
        assert np.isnan(pr[5]).all(), "a duplicated rainfall field was kept"
        assert np.isnan(cube.read("srad", 0, 10)[3]).all(), "a constant radiation field was kept"
        assert cube.qc_masked == {"pr": 1, "srad": 1}
        cube.close()


def test_source_guards_refuse_gaps_and_unfinished_files():
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        days = np.arange(np.datetime64("2015-01-01"), np.datetime64("2015-01-21"))
        holed = np.concatenate([days[:8], days[12:]])
        for name, d, status in (("gap.nc", holed, "complete"), ("wip.nc", days, "in_progress")):
            p = Path(tmp) / name
            _lab_file(p, d, status)
            try:
                S.LabCube(p, step_days=1).close()
            except S.SourceError:
                continue
            raise AssertionError(f"{name} was accepted")


# ------------------------------------------------------------------ drought definitions (S1-S3)
def test_usdm_classes_follow_the_published_table():
    D.check_usdm_table()
    v = np.array([-0.49, -0.5, -0.79, -0.8, -1.29, -1.3, -1.59, -1.6, -1.99, -2.0, -2.09, 1.0])
    assert D.usdm_class(v).tolist() == [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0]
    # the two stored levels that the lab's category layer classifies differently
    assert D.usdm_class(np.array([-0.71]))[0] == 1 and D.lab_class(np.array([-0.71]))[0] == 2
    assert D.usdm_class(np.array([-1.28]))[0] == 2 and D.lab_class(np.array([-1.28]))[0] == 3


def test_every_rule_names_an_index_of_its_own_window_and_a_source():
    for name, d in C.DROUGHT_DEFS.items():
        digits = int("".join(ch for ch in d["index"] if ch.isdigit()))
        assert digits == d["window_days"], name
        assert d["index"] in C.LAB_INDICES and d["cite"] and d["confirmed"], name
        assert (d["op"] == "le") == (d["threshold"] < 0), f"{name}: threshold sign and operator disagree"


def test_sign_checks_catch_a_flipped_index():
    rng = np.random.default_rng(0)
    rain = rng.normal(0, 1, (40, 30))
    etr = rng.normal(0, 1, (40, 30))
    good = dict(spi30d=rain, spei30d=rain - 0.5 * etr, eddi30d=etr)
    assert D.check_signs(good, rain, etr)["ok"]
    try:
        D.check_signs(dict(good, eddi30d=-etr), rain, etr)
    except S.SourceError:
        return
    raise AssertionError("an EDDI with the wrong sign passed")


def test_flash_drought_needs_a_fast_fall():
    pct = np.array([[60.0], [45.0], [18.0], [15.0], [50.0], [35.0], [10.0]])
    got = D.flash_drought(pct)[:, 0].tolist()
    assert got == [False, False, True, True, False, False, True]      # 60/45 -> 18, 45 -> 15, 50 -> 10


# ------------------------------------------------------------------ rainfall and climatology (R1, L3)
def test_standardised_rainfall_is_standard_on_its_baseline():
    rng = np.random.default_rng(3)
    base = rng.gamma(2.0, 20.0, (40, 36, 60))
    base[:, :, :5] *= rng.random((40, 36, 5)) > 0.3                     # some cells with dry spells
    vals = base.reshape(-1, 60)
    dek = np.tile(np.arange(36), 40)
    z, se = K.standardize_rainfall(base, vals, dek, pool=1, n_boot=4)
    assert abs(np.nanmean(z)) < 0.05 and 0.9 < np.nanstd(z) < 1.1
    assert np.isfinite(z).all() and np.isfinite(se).all() and (se >= 0).all()
    order = np.argsort(vals[:, 10])
    assert (np.diff(z[order, 10][dek[order] == 0]) >= -1e-6).all(), "score is not monotonic in rainfall"


def test_climatology_is_per_cell_and_reports_its_error():
    dek = np.tile(np.arange(36), 12)
    level = np.array([1.0, 50.0])                                       # a dry cell and a wet cell
    x = level[None, :] + np.random.default_rng(0).normal(0, 1, (dek.size, 2))
    cl = K.DekadClimatology(2)
    cl.add(dek, x)
    mean, sd, se, years = cl.finalize(pool=1)
    assert np.allclose(mean, level[None, :], atol=0.6) and (years == 12).all()
    assert np.allclose(se, sd / np.sqrt(12), rtol=1e-5)
    anom = K.anomaly(x, dek, mean)
    assert abs(anom[:, 0].mean()) < 0.1 and abs(anom[:, 1].mean()) < 0.1, "a cell's own mean was not removed"


def test_dekad_calendar():
    s = K.dekad_starts("2016-01-01", "2017-01-01")
    assert s.size == 36 and (K.dek36(s) == np.arange(36)).all()
    assert K.dekad_end(np.array(["2016-02-21"], "datetime64[D]"))[0] == np.datetime64("2016-03-01")
    days = np.arange(np.datetime64("2016-01-01"), np.datetime64("2016-03-01"))
    tot = K.aggregate_dekads(days, np.ones((days.size, 1), "float32"), s[:6], "sum")[:, 0]
    assert tot.tolist() == [10, 10, 11, 10, 10, 9]


# ------------------------------------------------------------------ split (T1, L7)
def _folds(embargo=3, seed=C.SEED):
    starts = K.dekad_starts("2015-04-01", "2022-06-01")
    dry = np.random.default_rng(5).normal(0, 1, starts.size)
    return starts, SP.make_folds(starts, dry, embargo, seed=seed)


def test_split_is_70_20_10_and_every_block_is_tested_once():
    starts, f = _folds()
    assert SP.check(f)
    assert len(f["block_label"]) == 15 and f["block_label"][0] == "G2015" and f["block_label"][-1] == "G2022"
    for s in SP.summarize(f):
        assert 0.15 <= s["test"] <= 0.24, s
        assert 0.08 <= s["val"] <= 0.12, s
        assert 0.66 <= s["train_assigned"] <= 0.75, s
    mean = {k: np.mean([s[k] for s in SP.summarize(f)]) for k in ("train_assigned", "test", "val")}
    assert abs(mean["test"] - 0.20) < 0.01 and abs(mean["val"] - 0.10) < 0.015


def test_each_fold_tests_and_validates_every_climate_class():
    starts, f = _folds()
    for k, s in enumerate(f["split"]):
        for role in (SP.TEST, SP.VAL):
            cls = set(f["block_class"][np.unique(f["block"][s == role])].tolist())
            assert cls == {0, 1, 2}, f"fold {k} role {role} holds classes {cls}"
        assert (np.bincount(f["block_class"]) == 5).all()


def test_no_training_dekad_inside_the_embargo():
    for e in (0, 2, 4):
        starts, f = _folds(embargo=e)
        for s in f["split"]:
            held = np.where((s == SP.VAL) | (s == SP.TEST))[0]
            train = np.where(s == SP.TRAIN)[0]
            assert np.abs(train[:, None] - held[None, :]).min() > e


def test_embargo_follows_the_memory_of_the_data():
    rng = np.random.default_rng(0)
    white = rng.normal(0, 1, (250, 40))
    red = np.zeros_like(white)
    for t in range(1, 250):
        red[t] = 0.9 * red[t - 1] + white[t]
    assert SP.embargo_from_acf(white)[0] == C.EMBARGO_RANGE[0]
    assert SP.embargo_from_acf(red)[0] == C.EMBARGO_RANGE[1]


# ------------------------------------------------------------------ regridding and land cover (U2, C1-C4)
def _toy_grid():
    lat = np.array([40.0, 39.9, 39.8])
    lon = np.array([-101.0, -100.9, -100.8, -100.7])
    frac = np.ones((3, 4), "float32")
    frac[2, 3] = 0.2                                                    # mostly outside the aquifer
    return G.make_grid(lat, lon, frac)


def test_area_aggregation_returns_mean_spread_and_count():
    g = _toy_grid()
    assert g.n == 11
    fine_lat = 40.05 - 0.0125 - 0.025 * np.arange(12)                   # 4 x 4 pixels per cell
    fine_lon = -101.05 + 0.0125 + 0.025 * np.arange(16)
    cid = G.cell_ids(fine_lat, fine_lon, g)
    assert (cid[8:, 12:] == -1).all() and (np.bincount(cid[cid >= 0]) == 16).all()
    vals = np.repeat(np.repeat(np.arange(12).reshape(3, 4), 4, 0), 4, 1).astype("float32")
    vals[0, 0] = np.nan
    mean, sd, n = G.block_stats(vals + np.tile([[0, 1], [1, 0]], (6, 8)), cid, g.n)
    assert np.allclose(mean[1:], np.arange(12)[1:11] + 0.5) and np.allclose(sd[1:], 0.5)
    assert n[0] == 15 and (n[1:] == 16).all()


def test_class_fractions_sum_to_one_and_use_the_real_legend():
    g = _toy_grid()
    fine_lat = 40.05 - 0.0125 - 0.025 * np.arange(12)
    fine_lon = -101.05 + 0.0125 + 0.025 * np.arange(16)
    cid = G.cell_ids(fine_lat, fine_lon, g)
    cls = np.random.default_rng(0).integers(0, 6, (12, 16))
    frac, n = G.class_fractions(cls, cid, g.n, [(k,) for k in range(6)])
    assert np.allclose(frac.sum(axis=0), 1.0, atol=1e-6)
    irrigated = frac[list(L.GFSAD_IRRIGATED)].sum(axis=0)
    rainfed = frac[list(L.GFSAD_RAINFED)].sum(axis=0)
    fragments = frac[list(L.GFSAD_FRAGMENTS)].sum(axis=0)
    assert np.allclose(irrigated + rainfed + fragments + frac[0], 1.0, atol=1e-6)
    assert L.GFSAD_IRRIGATED == (1, 2) and L.GFSAD_RAINFED == (3,) and L.GFSAD_FRAGMENTS == (4, 5)
    assert "rainfed" in L.GFSAD_LEGEND[3] and "high" in L.GIR_LEGEND[2]


def test_exposure_variables_cannot_be_model_inputs():
    C.assert_no_exposure(C.XGB_DYNAMIC + C.SEASON_FEATURES + C.STATIC_FEATURES)
    C.assert_no_exposure(C.DLSTM_STATIC)
    for bad in ("WTD", "f_irrigated", "gir_high_freq", "ET_clim"):
        try:
            C.assert_no_exposure(C.STATIC_FEATURES + (bad,))
        except ValueError:
            continue
        raise AssertionError(f"{bad} was accepted as a model input")


# ------------------------------------------------------------------ uncertainty (U1)
def test_conformal_interval_reaches_nominal_coverage():
    rng = np.random.default_rng(1)
    x = rng.uniform(0.5, 3, 20000)
    y = rng.normal(0, x)
    lo, hi = -0.8 * x, 0.8 * x                                          # a band that is too narrow
    cal, test = slice(0, 6000), slice(6000, None)
    assert U.coverage(y[test], lo[test], hi[test]) < 0.7
    m = U.cqr_margin(y[cal], lo[cal], hi[cal], alpha=0.10)
    cov = U.coverage(y[test], lo[test] - m, hi[test] + m)
    assert 0.88 <= cov <= 0.92, cov


def test_clustered_interval_is_wider_than_the_old_one():
    rng = np.random.default_rng(2)
    kd, n = 120, 300
    tblock, sblock = np.arange(kd) // 12, np.arange(n) // 30
    values = (rng.normal(0, 1, (10, 1))[tblock] + rng.normal(0, 1, (1, 10))[:, sblock]
              + rng.normal(0, 0.3, (kd, n)))
    mask = np.ones((kd, n), bool)
    b = U.block_bootstrap(values, mask, sblock, tblock, n_boot=300)
    old = U.naive_ci(values, mask)
    assert (b["ci_high"] - b["ci_low"]) / 2 > 10 * old["half_width"]
    assert b["ci_low"] < 0 < b["ci_high"]
    d = U.block_bootstrap(values + 2.0 * (sblock < 5)[None, :], mask & (sblock < 5)[None, :], sblock, tblock,
                          n_boot=200, mask_b=mask & (sblock >= 5)[None, :])
    assert d["ci_low"] < 2.0 < d["ci_high"] or abs(d["estimate"] - 2.0) < 1.5


# ------------------------------------------------------------------ cross-fitting (L2)
def _xgb_run(noise):
    from .. import train_xgb as TX
    from . import synth
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        keep = C.STORE
        C.STORE = Path(tmp)
        C.PROFILES["_test"] = C.Profile("_test", None)
        try:
            st, folds = synth.make_store(Path(tmp) / "_test", n_cells=24, noise_target=noise)
            m = TX.run("_test", rounds=60, threads=2)
            oof_fold = np.load(Path(tmp) / "_test" / "xgb" / "oof_fold.npy")
            bi = np.load(Path(tmp) / "_test" / "xgb" / "Bi.npy")
        finally:
            C.STORE = keep
            C.PROFILES.pop("_test", None)
    return m, oof_fold, bi, folds


def test_every_residual_is_out_of_fold():
    m, oof_fold, bi, folds = _xgb_run(noise=False)
    assert (oof_fold >= 0).all() and np.isfinite(bi).all()
    assert (oof_fold == folds["block_test_fold"][folds["block"]]).all(), "a dekad was predicted by the wrong fold"
    assert m["pooled_out_of_fold"]["skill_vs_climatology"] > 0.2, "the synthetic signal was not learned"
    assert 0.80 <= m["pooled_coverage90"] <= 0.97, m["pooled_coverage90"]
    C.assert_no_exposure(m["features"])


def test_noise_target_gives_no_out_of_fold_skill():
    """If residuals were in-sample, a model could 'explain' pure noise."""
    m, _, _, _ = _xgb_run(noise=True)
    assert m["pooled_out_of_fold"]["skill_vs_climatology"] < 0.03


def test_in_sample_residuals_would_have_looked_like_skill():
    """Documents L2: on a pure-noise target a model scores well on the rows it
    was fitted to and not at all on a held-out block."""
    import xgboost as xgb
    from .. import train_xgb as TX
    from . import synth
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        st, folds = synth.make_store(Path(tmp) / "s", n_cells=24, noise_target=True)
        d = TX.Design(st)
        s = folds["split"][0]
        d.set_fold(s)
        rows = lambda role: [np.concatenate(x) for x in zip(*[d.block(k)[:2] for k in np.where(s == role)[0]])]
        (xtr, ytr), (xte, yte) = rows(SP.TRAIN), rows(SP.TEST)
        bst = xgb.train(dict(C.XGB_PARAMS, objective="reg:squarederror", nthread=2),
                        xgb.DMatrix(xtr, ytr), 80)
        skill = lambda x, y: 1 - np.mean((y - bst.predict(xgb.DMatrix(x))) ** 2) / np.mean(y ** 2)
        assert skill(xtr, ytr) > 0.15 and skill(xte, yte) < 0.03, (skill(xtr, ytr), skill(xte, yte))
        del d


def test_a_missing_lag_does_not_drop_the_dekad():
    from .. import train_xgb as TX
    from . import synth
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        st, folds = synth.make_store(Path(tmp) / "s", n_cells=12)
        d = TX.Design(st)
        d.set_fold(folds["split"][0])
        X, y, cells = d.block(0)                                        # first dekad: no previous soil moisture
        assert X.shape[0] == 12 and np.isnan(X[:, d.names.index("SMrz_lag1_anom")]).all()
        del d


# ------------------------------------------------------------------ tower comparison (U3)
def test_tower_dekads_are_screened_and_gap_aware():
    from .. import towers as TW
    starts = K.dekad_starts("2016-06-01", "2016-07-01")
    t = (np.arange(np.datetime64("2016-06-01T00:00"), np.datetime64("2016-07-01T00:00"), np.timedelta64(30, "m"))
         - np.datetime64("1970-01-01T00:00")).astype("timedelta64[s]").astype("float64")
    le = np.full(t.size, 100.0)
    le[5], le[6], le[7] = 14980.0, -9999.0, np.inf                    # values the old code kept or summed
    le[480:480 + 200] = np.nan                                        # a 4-day gap in the second dekad
    clean = TW.screen_le(le)
    assert np.isnan(clean[5:8]).all()
    et, cov = TW.dekad_et(t, clean, starts)
    full = 100.0 * 10 * 86400 / TW.LATENT_HEAT                        # 35.3 mm for a complete dekad
    assert abs(et[0] - full) < 0.01 and abs(et[2] - full) < 0.01
    assert np.isnan(et[1]) and cov[1] < TW.MIN_COVERAGE, "a dekad with a long gap was accepted"
    le2 = clean.copy()
    le2[480:480 + 60] = np.nan                                        # a short gap: mean flux, not a short sum
    le2[480 + 60:480 + 200] = 100.0
    assert abs(TW.dekad_et(t, le2, starts)[0][1] - full) < 0.01
