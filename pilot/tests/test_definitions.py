"""The pilot's definitions: worked examples for every term, the registry, and "never re-implemented"."""
import os
import re

import numpy as np

from pilot import definitions as DEF

D04 = dict(min_share=0.5)
D05 = dict(min_share=0.5)
D06 = dict(primary="LANID", check="AIM-HPA", min_agree_years=3, min_map_coverage=0.9)
D07 = dict(max_irr_share=0.02, min_cropland_share=0.5)
D08 = dict(dominant_min_share=0.5, dominant_min_years=2, rzws_tol_mm=25)
D12 = dict(max_irr_share=0.05, min_label_coverage=0.9)
D13 = dict(radius_km=100, fallback_km=150, min_footprints=3, min_share_counties=0.8)
D14 = dict(rain_tol_frac=0.2, rzws_tol_frac=0.3)
D20 = dict(season_pct=70, s_pct=30, missing_s_is_dry=True)
D25 = dict(max_se_m=15, third_z=1.96, max_well_km=8)


def _cohort(primary, corn, sw=False, regular=True, check=None):
    """primary: yearly shares of the map that decides "irrigated" (D06); check: the other map (same if None)."""
    primary = np.array([primary], float)
    check = primary if check is None else np.array([check], float)
    ok, why = DEF.is_cohort(primary, check, np.array([corn], float), np.array([sw]), np.array([regular]), D04, D05,
                            D06)
    return bool(ok[0]), why[0]


def test_cohort_worked_examples():
    # irrigated 4 of 5 years, corn 2 of 3 -> cohort
    assert _cohort([0.9, 0.9, 0.9, 0.9, 0.0], [0.0, 0.8, 0.8]) == (True, "")
    # irrigated only 3 of 5 years
    assert _cohort([0.9, 0.9, 0.9, 0.0, 0.0], [0.8, 0.8, 0.8])[1] == "irrigated < 4 of 2011-2015"
    # corn only 1 of 3
    assert _cohort([0.9] * 5, [0.8, 0.1, 0.1])[1] == "corn < 2 of 2013-2015"
    # served by surface water
    assert _cohort([0.9] * 5, [0.8] * 3, sw=True)[1] == "surface water"
    # irregular quarter
    assert _cohort([0.9] * 5, [0.8] * 3, regular=False)[1] == "irregular cell"
    # maps agree in only 2 of 5 years
    assert _cohort([0.9] * 5, [0.8] * 3, check=[0.9, 0.9, 0.0, 0.0, 0.0])[1] == "irrigation maps disagree"
    # the PRIMARY map decides "irrigated": irrigated in the check map only is not cohort ...
    assert _cohort([0.9, 0.9, 0.9, 0.0, 0.0], [0.8] * 3, check=[0.9] * 5)[1] == "irrigated < 4 of 2011-2015"
    # ... and the check map may miss one year without removing a cell the primary map calls irrigated
    assert _cohort([0.9] * 5, [0.8] * 3, check=[0.9, 0.9, 0.9, 0.9, 0.0]) == (True, "")
    # the primary map unknown in a year (cell cut by the map edge) -> no label
    assert _cohort([0.9, 0.9, 0.9, 0.9, np.nan], [0.8] * 3, check=[0.9] * 5)[1] == "no irrigation map"


def test_map_share_and_roles():
    # 100 pixels, 95 labelled, 57 irrigated: known, share 0.6
    s, cov = DEF.map_share([57], [95], [100], D06)
    assert abs(s[0] - 0.6) < 1e-12 and cov[0] == 0.95
    # only 40% of the cell labelled (cell cut by the aquifer boundary): unknown, whatever the sliver says
    s, cov = DEF.map_share([40], [40], [100], D06)
    assert np.isnan(s[0]) and cov[0] == 0.4
    assert DEF.map_roles(D06) == ("lanid", "aim")
    assert DEF.map_roles(dict(D06, primary="AIM-HPA", check="LANID")) == ("aim", "lanid")
    try:
        DEF.map_roles(dict(D06, check="LANID"))
        raise AssertionError("the same map cannot be primary and check")
    except ValueError:
        pass
    from pilot import decisions as D
    assert D.value("D06")["primary"] == "LANID" and D.value("D06")["check"] == "AIM-HPA"
    # exactly at the threshold counts (>=)
    assert _cohort([0.5] * 5, [0.5] * 3) == (True, "")


def test_cohort_is_fixed_by_construction():
    """Later years are never an input: a cell that stops irrigating in 2019 cannot leave the cohort."""
    import inspect
    sig = inspect.signature(DEF.is_cohort)
    assert not any("201" in p for p in sig.parameters)
    assert _cohort([0.9] * 5, [0.8] * 3) == (True, "")


def test_rainfed_worked_examples():
    ok, why = DEF.is_rainfed(np.array([[0.0] * 5]), np.array([[0.01] * 5]), np.array([[0.9] * 5]), np.array([True]), D07)
    assert ok[0] and why[0] == ""
    ok, why = DEF.is_rainfed(np.array([[0.0, 0, 0, 0.3, 0]]), np.array([[0.0] * 5]), np.array([[0.9] * 5]),
                             np.array([True]), D07)
    assert why[0] == "irrigated in 2011-2015"
    ok, why = DEF.is_rainfed(np.array([[0.0] * 5]), np.array([[np.nan] * 5]), np.array([[0.9] * 5]), np.array([True]), D07)
    assert why[0] == "no irrigation map"                 # both maps must be known
    ok, why = DEF.is_rainfed(np.array([[0.0] * 5]), np.array([[0.0] * 5]), np.array([[0.2] * 5]), np.array([True]), D07)
    assert why[0] == "not cropland"
    assert DEF.irrigated_later(np.array([[[0.0, 0.0], [0.3, 0.0]]]), D07)[0]


def test_crop_class_and_matching():
    names = ["corn", "wheat"]
    shares = np.array([[[0.8, 0.1], [0.1, 0.8], [0.9, 0.0]], [[0.3, 0.3], [0.2, 0.2], [0.4, 0.4]]])
    dom = DEF.dominant_class(shares, names, D08)
    assert list(dom[0]) == ["corn", "wheat", "corn"] and list(dom[1]) == ["", "", ""]
    cls = DEF.crop_class(dom, D08)
    assert list(cls) == ["corn", "none"]
    assert list(DEF.crop_class(np.array([["corn", "wheat", "sorghum"]], dtype=object), D08)) == ["rotation"]
    m, k = DEF.match_rainfed(["A", "A", "A", "B"], ["corn", "corn", "corn", "corn"], [100, 120, 200, 100],
                             [True, False, False, False], [False, True, True, True], D08)
    assert list(m) == [False, True, False, False]       # same county+crop within 25 mm only
    assert list(k) == [1, 0, 0, 0]


def test_deficit_peak_month():
    assert DEF.deficit(150, 40) == 110
    pm = DEF.peak_month(np.array([100, 90, 50]), np.array([120, 90, np.nan]), np.array([110, 80, 60]))
    assert list(pm) == [7, 6, 0]                         # tie -> earlier; NaN -> undefined


def test_dry_and_episodes():
    x = np.arange(1, 11, dtype=float)
    dry, thr = DEF.season_dry(x, D20)
    assert abs(thr - np.percentile(x, 70)) < 1e-9 and dry.sum() == 3
    S = np.array([-2, -1, 0, 1, 2, np.nan])
    d, t = DEF.month_dry(S, D20)
    assert d[-1] and d[0] and not d[-2]                  # missing S counts as dry
    year = np.array([2020, 2020, 2020, 2021, 2022])
    cty = np.array(["a", "b", "c", "a", "a"])
    drys = np.array([True, True, True, True, True])
    nb = {"a": {"b"}, "b": {"a"}, "c": set()}
    ep, epm = DEF.dry_episodes(year, cty, drys, nb)
    assert ep[0] == ep[1] != ep[2]                       # a-b neighbours, c alone
    assert len(set(ep)) == 4
    assert epm[0] == epm[3] == epm[4] == epm[1]          # a dry 2020-2022 links the years


def test_thickness():
    A, se = DEF.thickness(300.0, 3.0, 250.0, 4.0)
    assert A == 50 and abs(se - 5.0) < 1e-12
    cuts = DEF.thickness_thirds(np.arange(1, 10), dict(quantiles=[1 / 3, 2 / 3]))
    assert list(DEF.assign_third([1, 5, 9], cuts)) == ["thin", "middle", "thick"]
    cuts = [30.0, 50.0]
    # ok thin; SE above 15 m; well too far; A not positive; thin that could be thick (10+1.96*14 > 50);
    # thick that could be thin (55-1.96*14 < 30); a middle cell with the same wide SE is not third-tested
    poor = DEF.thickness_poor([10, 60, 60, -1, 25, 55, 40], [5, 16, 5, 1, 14, 14, 14],
                              [2, 2, 9, 2, 2, 2, 2], ["thin", "thick", "thick", "thin", "thin", "thick", "middle"],
                              cuts, D25)
    assert list(poor) == [False, True, True, True, True, True, False]
    assert not DEF.thickness_poor([10], [9.0], [2], ["thin"], cuts, D25)[0]     # 10 m +- 17.6 stays below 50


def test_supported_month_and_footprint():
    ok, n, t = DEF.is_supported_month([1, 3, 5, 7, 9, 12, 15, 22], 30)
    assert ok and n == 8 and t
    ok, n, t = DEF.is_supported_month([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 30)
    assert not ok and not t                              # nothing in the last third
    ok, n, t = DEF.is_supported_month([1, 1, 3, 12, 15, 22, 25, 28], 31)
    assert not ok and n == 7                             # a date counts once
    v = np.full(16, np.nan)
    v[:8] = 0.2
    assert DEF.footprint_value(v, dict(min_nested_9km=8)) == 0.2
    v[7] = np.nan
    assert np.isnan(DEF.footprint_value(v, dict(min_nested_9km=8)))
    S, n = DEF.soil_moisture_state([1.0, -1.0, np.nan], [True, True, True], dict(weights="equal"),
                                   dict(min_supported_footprints=2))
    assert S == 0.0 and n == 2
    S, n = DEF.soil_moisture_state([1.0, -1.0], [True, False], dict(weights="equal"), dict(min_supported_footprints=2))
    assert np.isnan(S)


def test_reference_footprint_and_radius():
    ok, why = DEF.is_reference_footprint([0.01, 0.1, 0.01, 0.01, 0.01], [1, 1, 0.5, 1, 1], [50, 50, 50, 150, 50],
                                         [1.0, 1.0, 1.0, 1.0, 1.5], [1, 1, 1, 1, 1], [True] * 5, 100, D12, D14)
    assert list(why) == ["", "irrigated >= 5%", "irrigation map coverage", "beyond radius", "rainfall not similar"]
    r, share = DEF.footprint_radius({100: np.array([3, 3, 2, 3, 3]), 150: np.array([5] * 5)}, D13)
    assert r == 100 and share == 0.8
    r, _ = DEF.footprint_radius({100: np.array([3, 2, 2, 3, 3]), 150: np.array([5] * 5)}, D13)
    assert r == 150


def test_et_support_and_masks():
    assert list(DEF.et_model_ok([0.85, 0.79, np.nan], dict(min_valid_share=0.8))) == [True, False, False]
    assert list(DEF.et_supported([4, 3], dict(min_models=4))) == [True, False]
    m = DEF.sample_masks(["cohort", "cohort", "rainfed", "cohort"], [True, True, True, True],
                         [True, False, True, True], [True, True, False, False], [True, True, True, False])
    assert list(m["mask_sm"]) == [True, False, True, False]
    assert list(m["mask_peak"]) == [True, True, True, False]
    assert list(m["mask_timing"]) == [True, True, True, False]


def test_contrasts():
    assert DEF.contrast_sm(2.0, 20.0, 80.0, -1.0, 0.0) == 120.0
    assert DEF.contrast_peak(-0.1, 20.0, 80.0, 200.0, 150.0) == 300.0


def test_registry_complete():
    for term, wording, funcs, dids, task in DEF.TERMS:
        for f in funcs:
            assert callable(getattr(DEF, f)), f"{term}: {f} missing"
            assert getattr(DEF, f).__doc__, f"{term}: {f} has no docstring"
    from pilot import decisions as D
    for _, _, _, dids, _ in DEF.TERMS:
        for did in dids:
            D.entry(did)


THRESHOLD_KEYS = ["min_share", "max_irr_share", "min_cropland_share", "season_pct", "s_pct", "min_valid_share",
                  "min_models", "max_se_m", "third_z", "max_well_km", "min_nested_9km", "min_label_coverage",
                  "min_supported_footprints", "rzws_tol_mm", "dominant_min_share", "min_agree_years",
                  "min_map_coverage",
                  "rain_tol_frac", "rzws_tol_frac"]


def test_no_task_reimplements_a_rule():
    """Decision thresholds may be read only in definitions.py: a task that reads one re-implements a rule."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    bad = []
    for sub in ("phase0", "phase1", "phase2"):
        for f in os.listdir(os.path.join(root, sub)):
            if f.endswith(".py"):
                txt = open(os.path.join(root, sub, f), encoding="utf-8").read()
                for k in THRESHOLD_KEYS:
                    if re.search(rf"\[[\"']{k}[\"']\]", txt):
                        bad.append(f"{sub}/{f}: {k}")
    assert not bad, bad


def test_docs_current():
    from pilot import docs_gen
    assert docs_gen.current(), "docs/pilot is out of date: run python -m pilot docs"
