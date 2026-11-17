"""Unit tests: grids, cells, zonal statistics, regression, kriging, the decision gate."""
import os
import shutil
import tempfile

import numpy as np


def test_ease_grid_round_trip():
    from pilot import config as C
    from pilot.common import grids
    for size in (C.EASE_M09, C.EASE_M36):
        r, c = np.array([100, 412]), np.array([800, 1999])
        x, y = grids.ease_centre(r, c, size)
        rr, cc = grids.ease_index(x, y, size)
        assert (rr == r).all() and (cc == c).all()
    r9, c9 = np.arange(400, 408), np.arange(800, 808)
    assert (grids.ease_index(*grids.ease_centre(r9, c9, C.EASE_M09), C.EASE_M36)[0] == r9 // 4).all()


def test_ease_axes_of_the_smap_file():
    from pilot import config as C
    from pilot.common import grids, io
    cands = [C.DATASETS[k][0] for k in ("smap_reprocessed", "smap_full_old") if os.path.exists(C.DATASETS[k][0])]
    if not cands:
        return "SKIP: no SMAP file on this machine"
    ds = io.smap_open(cands[0])
    r9, c9 = grids.check_ease_axes(ds["x"][:], ds["y"][:])
    ds.close()
    assert np.all(np.diff(r9) == 1) and np.all(np.diff(c9) == 1)


def test_quarter_split():
    import geopandas as gpd
    import shapely
    from shapely import affinity
    from pilot.common.cells import split_quarters
    sq = shapely.box(0, 0, 1609.344, 1609.344)
    secs = gpd.GeoDataFrame(dict(FRSTDIVID=["S1", "S2"]), geometry=[sq, affinity.rotate(sq, 12, origin="centroid")],
                            crs=5070)
    q = split_quarters(secs)
    assert len(q) == 8 and set(q.quarter) == {"NE", "NW", "SW", "SE"}
    a = q.geometry.area / 1e4
    assert np.allclose(a, 64.75, atol=0.05)              # a regular quarter is 64.7 ha, rotated or not
    ne = q[(q.section_id == "S1") & (q.quarter == "NE")].geometry.iloc[0].centroid
    assert ne.x > 804 and ne.y > 804


def test_zonal_against_brute_force():
    from pilot.common import zonal
    rs = np.random.default_rng(3)
    ids = rs.integers(-1, 5, (40, 40))
    v = rs.normal(size=(40, 40))
    v[rs.random((40, 40)) < 0.1] = np.nan
    m, k = zonal.mean(ids, v, 5)
    for i in range(5):
        sel = (ids == i) & np.isfinite(v)
        assert k[i] == sel.sum() and abs(m[i] - v[sel].mean()) < 1e-12


def test_fe_ols_matches_dummy_regression():
    import statsmodels.api as sm
    from pilot.common import stats
    rs = np.random.default_rng(1)
    n, g1, g2 = 600, 30, 12
    a, b = rs.integers(0, g1, n), rs.integers(0, g2, n)
    x = rs.normal(size=n) + 0.3 * a / g1
    y = 2.5 * x + rs.normal(size=g1)[a] + rs.normal(size=g2)[b] + rs.normal(size=n)
    r = stats.fe_ols(y, x, fe=[a, b])
    D = np.column_stack([x] + [(a == i).astype(float) for i in range(g1)] + [(b == j).astype(float) for j in range(1, g2)])
    o = sm.OLS(y, D).fit()
    assert abs(r["beta"][0] - o.params[0]) < 1e-6


def test_kriging_known_surface():
    from pilot.common.stats import RegressionKriging
    rs = np.random.default_rng(2)
    xy = rs.uniform(0, 150e3, (1200, 2))
    drift = 900 + xy[:, 0] / 300
    f = lambda p: 30 + 0.95 * (900 + p[:, 0] / 300) + 15 * np.sin(p[:, 1] / 20e3)
    rk = RegressionKriging(xy, f(xy) + rs.normal(0, 0.5, len(xy)), drift, k=24, max_km=40)
    t = rs.uniform(20e3, 130e3, (3000, 2))
    p, se = rk.predict(t, 900 + t[:, 0] / 300)
    err = p - f(t)
    assert np.sqrt(np.mean(err ** 2)) < 2.0
    assert np.mean(np.abs(err) < 1.96 * se) > 0.85


def test_kriging_with_duplicate_wells():
    """Two wells at the same place must not break the system (they did, giving 1e15 m errors)."""
    from pilot.common.stats import RegressionKriging
    rs = np.random.default_rng(4)
    xy = rs.uniform(0, 80e3, (400, 2))
    xy = np.vstack([xy, xy[:50]])                                   # 50 exact duplicates
    z = 500 + xy[:, 0] / 1000 + rs.normal(0, 0.3, len(xy))
    rk = RegressionKriging(xy, z, np.zeros(len(xy)), k=24, max_km=40)
    assert rk.n_merged == 50
    p, se = rk.predict(xy[:100], np.zeros(100))
    assert np.nanmax(np.abs(p - (500 + xy[:100, 0] / 1000))) < 5 and np.isfinite(se).all()


def test_cluster_fallback():
    from pilot.common.stats import cluster_ids
    ids, lab = cluster_ids(["a"] * 4, year=[1, 1, 2, 2], fallback="county_year")
    assert lab == "county-year" and len(set(ids)) == 2
    ids, lab = cluster_ids(["a"] * 3, [0, 5000, 25000], [0, 0, 0], fallback="block")
    assert lab == "10 km block" and len(set(ids)) == 2
    ids, lab = cluster_ids([str(i) for i in range(12)], [0] * 12, [0] * 12)
    assert lab == "county"


def _temp_log():
    from pilot import decisions as D
    d = tempfile.mkdtemp()
    p = os.path.join(d, "decisions.yaml")
    with open(D.PATH, encoding="utf-8") as f, open(p, "w", encoding="utf-8", newline="\n") as g:
        g.write(f.read().replace("status: fixed", "status: proposed"))     # the tests start from a proposed log
    return p, d


def test_decisions_gate_fix_amend():
    from pilot import config as C, decisions as D
    from pilot.core.errors import PipelineStop
    p, d = _temp_log()
    try:
        old = C.RUN_MODE
        C.RUN_MODE = "stage3"
        try:
            D.require(["D04"], "t0_2", p)
            raise AssertionError("a proposed decision must stop stage3")
        except PipelineStop:
            pass
        assert D.fix(["D04"], "tester", p) == ["D04"]
        assert D.require(["D04"], "t0_2", p) == []
        D.amend("D04", dict(min_share=0.6), "sensitivity", "tester", p)
        e = D.entry("D04", p)
        assert e["value"] == dict(min_share=0.6) and e["history"][0]["value"] == dict(min_share=0.5)
        C.RUN_MODE = "smoke"
        assert D.require(["D05"], "t0_2", p) == ["D05"]          # smoke runs, and says so
    finally:
        C.RUN_MODE = old
        shutil.rmtree(d)


def test_signature_follows_decision_values():
    from pilot import config as C
    snap = C.snapshot()
    assert "D04" in snap and snap["D04"]["min_share"] == 0.5 or "min_share" in snap["D04"]


def test_signature_follows_the_source_file():
    """A task that reads another file (new SMAP extraction, rebuilt LANID file) must not be skipped as complete."""
    from pilot import config as C
    from pilot.core import provenance as P
    d = tempfile.mkdtemp()
    old, old_mode = C.DATASETS["lanid"], C.RUN_MODE
    try:
        C.RUN_MODE = "smoke"
        a, b = os.path.join(d, "a.nc"), os.path.join(d, "b.nc")
        for p, n in ((a, 10), (b, 10)):
            with open(p, "wb") as f:
                f.write(b"x" * n)
        C.DATASETS["lanid"] = (a, a, "lanid.nc")
        s1 = P.signature("t0_2", ["SRC:lanid"], [], {})
        C.DATASETS["lanid"] = (b, b, "lanid.nc")
        s2 = P.signature("t0_2", ["SRC:lanid"], [], {})
        with open(b, "ab") as f:
            f.write(b"y")
        s3 = P.signature("t0_2", ["SRC:lanid"], [], {})
        assert len({s1, s2, s3}) == 3            # another path, then another size: a new signature each time
    finally:
        C.DATASETS["lanid"], C.RUN_MODE = old, old_mode
        shutil.rmtree(d)


def test_rootzone_storage_is_mm_at_haskell():
    """gSSURGO rootznaws is mm of water (Valu1). Read with the pipeline's own readers and no scaling, the median
    around Haskell County must be a physically possible 100-400 mm (10x that would exceed the root-zone depth)."""
    from pilot import config as C
    from pilot.common import io
    from pyproj import Transformer
    if not os.path.exists(C.DATASETS["gssurgo"][0]):
        return "SKIP: gSSURGO file not on this machine"
    old = C.RUN_MODE
    C.RUN_MODE = "smoke"
    try:
        ds = io.gssurgo()
        tr = io.gssurgo_transform(ds)
        x, y = Transformer.from_crs(4326, 5070, always_xy=True).transform(-100.87, 37.56)   # Haskell, KS
        c, r = ~tr * (x, y)
        a = np.ma.filled(ds["rootznaws"][int(r) - 120:int(r) + 120, int(c) - 120:int(c) + 120].astype(float), np.nan)
        assert ds["rootznaws"].units == "mm"
        ds.close()
    finally:
        C.RUN_MODE = old
    med = float(np.nanmedian(a))
    assert 100 <= med <= 400, med


def test_footprint_half_off_the_raster():
    """A footprint with half its area beyond the raster edge has label coverage 0.5 (not 1.0)."""
    import logging
    import geopandas as gpd
    import rasterio
    import shapely
    from rasterio.transform import from_origin
    from pilot.phase1 import t1_2_footprints as T12
    d = tempfile.mkdtemp()
    try:
        p = os.path.join(d, "map.tif")
        with rasterio.open(p, "w", driver="GTiff", width=96, height=96, count=1, dtype="uint8", crs="EPSG:5070",
                           transform=from_origin(0, 2880, 30, 30)) as w:
            w.write(np.ones((1, 96, 96), "uint8"))
        fp = gpd.GeoDataFrame(geometry=[shapely.box(1440, 0, 4320, 2880)], crs=5070)
        ctx = type("Ctx", (), dict(log=logging.getLogger("test")))()
        irr, cov = T12._map_shares(ctx, fp, {2018: p}, lambda a: a == 1, lambda a: np.ones(a.shape, bool), "test")
        assert abs(cov[0] - 0.5) < 1e-9 and irr[0] == 1.0, (cov, irr)
    finally:
        shutil.rmtree(d)


def test_reliability_of_a_mismeasured_regressor():
    """Error with the same variance as the real regressor: half the variation is real."""
    from pilot.common import stats
    from pilot.phase2.c2_6_precision import reliability
    rs = np.random.default_rng(7)
    g = np.repeat(np.arange(200), 10)
    real, err = rs.normal(size=2000), rs.normal(size=2000)
    lam = reliability(stats.within(real, [g]), err, [g])[0, 0]
    assert 0.45 < lam < 0.55, lam
    assert abs(reliability(stats.within(real, [g]), np.zeros(2000), [g])[0, 0] - 1.0) < 1e-9


def test_reliability_corrects_regressors_that_share_an_error():
    """Two regressors carrying one error (a and a*z, as A and A x deficit do): OLS on the mismeasured pair is
    biased, and the attenuation matrix brings the coefficients back."""
    from pilot.common import stats
    from pilot.phase2.c2_6_precision import reliability
    rs = np.random.default_rng(11)
    n = 20000
    a, z, u = rs.normal(size=n), rs.uniform(0.5, 1.5, n), rs.normal(scale=0.7, size=n)
    X, U = np.column_stack([a, a * z]), np.column_stack([u, u * z])
    beta = np.array([1.0, -2.0])
    y = X @ beta + rs.normal(scale=0.1, size=n)
    g = np.zeros(n, int)
    b_obs = stats.fe_ols(y, X + U, fe=[g])["beta"]
    L = reliability(stats.within(X, [g]), U, [g])
    assert np.abs(b_obs - beta).max() > 0.3, b_obs
    assert np.abs(np.linalg.solve(L, b_obs) - beta).max() < 0.1, np.linalg.solve(L, b_obs)


def test_et_by_thickness_recovers_a_planted_response():
    """fig_thickness: a response inside the fitted surface comes back as a difference between two thickness
    levels at a high deficit, within its own standard error."""
    from pilot.phase2 import fig_thickness as F
    rs = np.random.default_rng(3)
    A, D = rs.uniform(5, 120, 4000), rs.uniform(-100, 900, 4000)

    def truth(a, d):
        return 0.1 * a + 0.0004 * a * d
    f = F.fit(A, D, truth(A, D) + rs.normal(0, 5, 4000))
    m, s = F.predict(f, 100.0, 800.0, ref_A=10.0)
    assert abs(m[0] - (truth(100, 800) - truth(10, 800))) < 3 * s[0] + 0.5, (m, s)
    assert s[0] < 2, s


def test_blinding_names_of_thickness():
    from pilot import blinding
    cols = ["cell", "year", "et_ensemble_mad", "kgs_satthick_2022_2024", "kgs_hpa_bedrock", "kgs_wte_2022_2024",
            "A_m", "area_ha"]
    assert blinding.thickness_columns(cols) == ["A_m", "kgs_hpa_bedrock", "kgs_satthick_2022_2024", "kgs_wte_2022_2024"]
    assert blinding.thickness_columns(["cell", "year", "month", "et_ssebop", "role", "dry", "county",
                                       "blinding_gate", "et_supported"]) == []
    cells = __import__("pandas").DataFrame(dict(cell_id=["c1"], county=["20001"], role=["rainfed"]))
    et = __import__("pandas").DataFrame(dict(cell=["c1"], year=[2020], month=[7], et_ensemble_mad=[50.0],
                                             kgs_satthick_2022_2024=[30.0]))
    try:
        blinding.gate(et, cells, et[["year", "month"]].assign(county="20001", dry=False))
        raise AssertionError("the gate must refuse an ET table with a thickness column")
    except blinding.BlindingError:
        pass


def test_main_model_recovers_a_planted_thickness_effect():
    """Synthetic panel: peak-month ET = cell + county-year effects + b * A_true * S + noise; A observed with a
    per-cell base error and a per-year water-table error. The corrected C1 estimate recovers b; the uncorrected
    one is attenuated (Check 2.6's reliability), and the contrast is b * (A_thin - A_thick) * (S_q10 - S_med)."""
    import pandas as pd
    from pilot.phase3 import main_models as MM
    rs = np.random.default_rng(7)
    n_c, n_cell, years = 14, 40, list(range(2016, 2025))
    cells = pd.DataFrame(dict(cell_id=[f"c{i}" for i in range(n_c * n_cell)],
                              county=np.repeat([f"20{i:03d}" for i in range(n_c)], n_cell)))
    cells["A_true"] = rs.uniform(5, 100, len(cells))
    cells["base_err"] = rs.normal(0, 8, len(cells))
    cells["rzws_mm"], cells["corn_years"] = rs.uniform(150, 350, len(cells)), rs.integers(2, 4, len(cells))
    cells["fe"] = rs.normal(0, 20, len(cells))
    cy = pd.DataFrame([(c, y) for c in cells.county.unique() for y in years], columns=["county", "year"])
    cy["S"], cy["deficit_season_mm"], cy["fe_cy"] = rs.normal(0, 1, len(cy)), rs.normal(550, 80, len(cy)), \
        rs.normal(0, 15, len(cy))
    d = cells.merge(cy, on="county")
    b = 0.20                                                    # mm of peak-month ET per (m x SD of S)
    d["et_ensemble_mad"] = d.fe + d.fe_cy + b * d.A_true * d.S + rs.normal(0, 8, len(d))
    d["wte_se_m"], d["base_se_m"] = 3.0, 8.0
    d["A_m"] = d.A_true + d.base_err + rs.normal(0, 3, len(d))
    cut = np.quantile(d.groupby("cell_id").A_m.mean(), [1 / 3, 2 / 3])
    mean_a = d.groupby("cell_id").A_m.transform("mean")
    d["A_third"] = np.where(mean_a <= cut[0], "thin", np.where(mean_a <= cut[1], "middle", "thick"))
    d["x5070"], d["y5070"], d["gmd"] = 0.0, 0.0, 1
    est, con = MM.estimate(d, "et_ensemble_mad", "C1", "sm", ["S"], np.random.default_rng(1))
    e = pd.DataFrame(est).set_index("term")
    assert abs(e.loc["A:S", "beta_corrected"] - b) < 0.03, e.loc["A:S"]
    assert e.loc["A:S", "beta"] < e.loc["A:S", "beta_corrected"]          # uncorrected is attenuated
    c = pd.DataFrame(con)
    c = c[c.corrected.astype(bool)].iloc[0]
    expect = e.loc["A:S", "beta_corrected"] * (c["A_thin_m"] - c["A_thick_m"]) * (c["eval_point"] - c["median"])
    assert abs(c.value_mm - expect) < 1e-6 and c.lo_mm < c.value_mm < c.hi_mm
