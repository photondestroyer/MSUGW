"""Unit tests of P1.2 pieces that need no torch (calendars, features, folds, baselines, conformal)."""
import numpy as np

from .. import config as C
from .. import cdl_groups as CG
from ..p1_2 import calendars as CA, features as FE, folds as FO, baselines as BL
from .util import override


def test_preseason_months_are_strictly_before_the_target():
    t = np.arange(2010 * 12, 2014 * 12)
    for start, n in [(4, 6), (10, 9), (1, 12), (12, 3)]:
        d = CA.derived(t, start, n)
        assert (d["preseason"] < t[:, None]).all()
        assert ((d["start"] % 12 + 1) == start).all() and (d["start"] <= t).all() and (t - d["start"] < 12).all()
        assert (d["in_season"] == ((t - d["start"]) < n)).all()


def test_season_window_contract():
    """Checks the TODO(human) function once it is written."""
    summer = np.array([0.05, 0.05, 0.06, 0.09, 0.15, 0.24, 0.30, 0.27, 0.18, 0.10, 0.06, 0.05])
    try:
        s, n = CA.season_window(summer)
    except NotImplementedError:
        return "SKIP calendars.season_window not written yet (Learn by Doing)"
    assert 1 <= s <= 12 and 1 <= n <= 12
    months = [(s - 1 + k) % 12 for k in range(n)]
    assert 6 in months and 0 not in months                     # July in season, January not
    wheat = np.array([0.12, 0.12, 0.16, 0.25, 0.33, 0.28, 0.10, 0.07, 0.07, 0.10, 0.14, 0.13])
    s2, n2 = CA.season_window(wheat)
    assert 4 in [(s2 - 1 + k) % 12 for k in range(n2)]         # May in season
    s3, n3 = CA.season_window(np.where(np.arange(12) == 3, np.nan, summer))
    assert 1 <= s3 <= 12 and 1 <= n3 <= 12                     # a missing month does not break it


def test_exposure_variables_are_refused():
    for bad in (["WTD"], ["hand_frac"], ["lanid_dist_m"], ["gde_max"], ["lat"], ["M0"]):
        try:
            FE.assert_no_exposure(bad + ["P"])
            raise AssertionError(f"{bad} accepted")
        except ValueError:
            pass
    names = FE.input_names()
    assert not set(names) & set(C.EXPOSURE_VARS)


def test_folds_hold_out_space_and_years_with_embargo():
    rng = np.random.default_rng(0)
    lon, lat = rng.uniform(-103, -99, 3000), rng.uniform(37, 41, 3000)
    ar = rng.integers(1, 11, 3000)
    with override(RUN_MODE="full"):
        F = FO.make(lon, lat, ar)
        _, gk = FO.fold_of(F, lon, lat)
        wy = rng.choice(np.arange(2008, 2026), lon.size)
        for k in range(F["folds"]):
            for g in range(len(F["year_groups"])):
                r = FO.roles(F, k, g, gk, wy)
                test_groups = set(gk[r == "test"])
                for role in ("train", "val", "cal"):
                    assert not (set(gk[r == role]) & test_groups)
                    years = set(wy[r == role])
                    tg = F["year_groups"][g]
                    assert not years & set(range(min(tg) - 1, max(tg) + 2))
                assert not (set(gk[r == "val"]) & set(gk[r == "train"]))
                assert not (set(gk[r == "cal"]) & set(gk[r == "train"]))
        sizes = np.bincount(list(F["group_fold"].values()))
        assert sizes.max() - sizes.min() <= 1


def test_conformal_margins_restore_coverage():
    rng = np.random.default_rng(1)
    y = rng.normal(0, 1, 20000)
    lo, hi = np.full(y.size, -0.8), np.full(y.size, 0.8)               # too narrow on purpose
    z = np.zeros(y.size, int)
    lv = BL.strata_levels(0, rng.integers(0, 3, y.size), rng.integers(1, 11, y.size), rng.random(y.size) < 0.5, z + 1, z)
    m = BL.fit_margins(y[:10000], lo[:10000], hi[:10000], lv[:, :10000], min_n=500)
    l2, h2 = BL.apply_margins(lo[10000:], hi[10000:], lv[:, 10000:], m)
    cov = BL.coverage(y[10000:], l2, h2)
    assert 0.88 <= cov <= 0.92, cov
    assert BL.margins_from_json(BL.margins_to_json(m)) == m


def _toy(nb=40, years=(2010, 2011, 2012, 2013)):
    rng = np.random.default_rng(2)
    months = np.arange(np.datetime64("2007-01"), np.datetime64("2014-01"))
    nm, ncell = months.size, 3
    D = dict(months=months, lab_shape=np.array([2, 2]))
    for v in C.DYN_INPUTS:
        D[v] = rng.gamma(2, 20, (nm, ncell)).astype("float32")
    D["BUCKET"] = rng.random((nm, len(C.BUCKET_CAPACITY_BINS_MM), ncell)).astype("float32")
    ya = np.array(years)
    BS = dict(years_agg=ya, years=ya, cell=rng.integers(0, ncell, nb), bucket_bin=rng.integers(0, 8, nb),
              aridity=rng.integers(1, 11, nb), static=rng.random((nb, len(C.STATIC_INPUTS))).astype("float32"),
              REF=np.ones((nb, ya.size), "u1"), GRP=np.full((nb, ya.size), CG.NAMES.index("grass_pasture"), "u1"),
              block_id=np.arange(nb), lon=rng.uniform(-101, -100, nb), lat=rng.uniform(40, 41, nb))
    Y = np.empty((nb, ya.size, 12, 3), "float32")
    for k in range(ya.size):
        for m in range(12):
            P3 = D["P"][np.searchsorted(months, np.datetime64(f"{ya[k]}-{m + 1:02d}")) - 2:
                        np.searchsorted(months, np.datetime64(f"{ya[k]}-{m + 1:02d}")) + 1].sum(0)
            Y[:, k, m, 0] = 5 + 0.3 * P3[BS["cell"]]
            Y[:, k, m, 1] = 0.1 + 0.001 * P3[BS["cell"]]
            Y[:, k, m, 2] = -0.1 + 0.002 * P3[BS["cell"]]
    BS["Y"] = Y
    BS["REF"][0, 1] = 0
    return BS, D


def test_features_shapes_masks_and_no_future_preseason():
    BS, D = _toy()
    cal = {g: (4, 6) for g in range(CG.N)}
    i = np.array([0, 1, 2])
    X = FE.build(BS, D, i, np.array([2012, 2012, 2013]), cal)
    assert X["dyn"].shape == (3, 24, len(C.DYN_INPUTS)) and X["pm"].shape == (3, 12, FE.n_pm())
    assert X["stat"].shape == (3, FE.n_stat()) and X["y"].shape == (3, 12, 3) and X["pm"].shape[-1] == FE.n_pm()
    X2 = FE.build(BS, D, np.array([0]), np.array([2012]), cal)
    assert not X2["ymask"][0, :3].any()                        # Oct-Dec 2011: block 0 not reference in 2011
    assert X2["ymask"][0, 3:].all()
    # pre-season of an April target is January-March of the same year, taken from Y
    jan_mar = BS["Y"][0, 2, 0:3, 1].mean()
    assert abs(X2["pm"][0, 6, FE.PM_NAMES.index("pre_NIRv")] - jan_mar) < 1e-6
    sc = FE.fit_scaler(X)
    Xs = FE.apply_scaler(X, sc)
    assert np.isfinite(Xs["dyn"]).all() and np.isfinite(Xs["pm"]).all()


def test_m2_recovers_a_rain_use_response():
    BS, D = _toy(nb=400)
    cal = {g: (4, 6) for g in range(CG.N)}
    i = np.arange(400)
    X = FE.build(BS, D, i, np.full(400, 2012), cal)
    with override(M2_MIN_SAMPLES=20):
        model = BL.m2_fit(X, X["ymask"].any(-1))
        p = BL.m2_predict(model, X)
    ok = X["ymask"][..., 0]
    err = np.abs(p[..., 0, 1][ok] - X["y"][..., 0][ok])
    assert np.median(err) < 0.5, np.median(err)
    ws = {v[1] for k, v in model.items() if k[0] == 0}
    assert 3 in ws                                           # the 3-month window that made the data is found


# ---------------------------------------------------------------- plan v3
def test_small_strata_borrow_the_margin_of_the_next_merge_level():
    rng = np.random.default_rng(3)
    n = 6000
    func = np.where(np.arange(n) < 5900, 0, 1)                         # functional group 1 has only 100 months
    y = rng.normal(0, np.where(func == 1, 3.0, 1.0))                   # ... and is much noisier
    lo, hi = np.full(n, -1.0), np.full(n, 1.0)
    lv = BL.strata_levels(0, func, np.full(n, 5), np.ones(n, bool), np.ones(n, int), np.zeros(n, int))
    m = BL.fit_margins(y, lo, hi, lv, min_n=1000)
    used = BL.level_of(lv, m)
    assert (used[func == 0] == 0).all()                                 # its own full stratum
    assert (used[func == 1] > 0).all()                                  # merged upward
    order = [d for d in C.CAL_MERGE_ORDER]
    assert order.index("regime") == 0 and order.index("aridity_tercile") == 1   # plan: regime first, then tercile
    l0 = BL.strata_levels(0, 0, 5, True, 1, 0)
    l1 = BL.strata_levels(1, 0, 5, True, 1, 0)
    assert not set(l0.ravel()) & set(l1.ravel())                       # targets never share a stratum


def test_regimes_follow_the_landsat_record():
    t = lambda y, m: y * 12 + m - 1
    names = [r[0] for r in C.OBS_REGIMES]
    got = [names[i] for i in FE.regime_of([t(2011, 12), t(2012, 6), t(2013, 3), t(2013, 4), t(2016, 12), t(2017, 1),
                                           t(2022, 1)])]
    assert got == ["L5+L7", "L7", "L7", "L7+L8", "L7+L8", "L8", "L8+L9"]


def test_measurement_inputs_class_inputs_and_ratio_target():
    BS, D = _toy()
    nb = BS["Y"].shape[0]
    BS["MEAS"] = np.random.default_rng(4).random((nb, BS["Y"].shape[1], 12, len(C.MEAS_INPUTS))).astype("float32")
    BS["GRP"][:] = CG.NAMES.index("winter_wheat")
    BS["GRP_PREV"] = np.full_like(BS["GRP"], CG.NAMES.index("fallow"))
    D["SPI90"] = np.full(D["P"].shape, -1.0, "float32")
    cal = {g: (4, 6) for g in range(CG.N)}
    X = FE.build(BS, D, np.array([0, 1]), np.array([2012, 2012]), cal)
    j = FE.pm_names().index(C.MEAS_INPUTS[0])
    assert abs(X["pm"][0, 3, j] - BS["MEAS"][0, 2, 0, 0]) < 1e-6        # January 2012 of block 0 (year index 2)
    ns = len(C.STATIC_INPUTS)
    names = FE.stat_names()
    assert X["stat"][0, names.index("func_winter_cereal")] == 1 and X["stat"][0, names.index("prev_func_fallow")] == 1
    assert (X["drought"] == 0).all() and (X["func"] == CG.FUNC_NAMES.index("winter_cereal")).all()
    Xn = FE.build(BS, D, np.array([0]), np.array([2012]), cal, no_meas=True)
    assert np.isnan(Xn["pm"][..., j]).all()
    Xr = FE.build(BS, D, np.array([0]), np.array([2012]), cal, ratio=True)
    ok = np.isfinite(Xr["y"][0, :, 0])
    assert np.allclose(Xr["y"][0, ok, 0] * Xr["pet"][0, ok], X["y"][0, ok, 0], rtol=1e-5)
    FE.set_class(X, np.array([1]), CG.FUNC_NAMES.index("fallow"))
    assert X["stat"][1, names.index("func_fallow")] == 1 and X["stat"][1, names.index("func_winter_cereal")] == 0
    assert X["stat"][1, names.index("group_fallow")] == 1 and X["func"][1] == CG.FUNC_NAMES.index("fallow")
    sc = FE.fit_scaler(X)
    assert np.isfinite(FE.apply_scaler(X, sc)["pm"]).all()
    for bad in (["canal_dist_m"], ["facc_max_km2"], ["state"], ["SPI90"], ["reference_ring"]):
        try:
            FE.assert_no_exposure(bad)
            raise AssertionError(f"{bad} accepted as an input")
        except ValueError:
            pass
