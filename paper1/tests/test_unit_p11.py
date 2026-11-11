"""Unit tests of the P1.1 rules on hand-made blocks."""
import numpy as np

from .. import config as C
from .. import cdl_groups as CG
from ..p1_1 import rules as R


def _blocks(n=8):
    """n blocks in a row that pass every rule; tests then break one rule per block."""
    S = dict(wtd_min=np.full((1, n), 50.0), wtd_coverage=np.ones((1, n)), hand_frac=np.ones((1, n)),
             hand_coverage=np.ones((1, n)), nwi_available=np.ones((1, n), "u1"), nwi_wetland_px=np.zeros((1, n), "u2"),
             nwi_riparian_px=np.zeros((1, n), "u2"), lanid_dist_m=np.full((1, n), 5000.0),
             lanid_unassessed_dist_m=np.full((1, n), 5000.0), county=np.zeros((1, n), "i4"),
             lanid3_dist_m=np.full((1, n), 5000.0), canal_dist_m=np.full((1, n), np.inf),
             roi_px=np.full((1, n), 256, "u2"), gde_max=np.full((1, n), 0.01), gde_assessed=np.ones((1, n)))
    B = dict(water_px=np.zeros((12, 1, n), "u2"), ET_n=np.full((12, 1, n), 256, "u2"),
             NIRv_n=np.full((12, 1, n), 256, "u2"), NDMI_n=np.full((12, 1, n), 256, "u2"))
    cdl = np.zeros((CG.N, 1, n), "int64")
    cdl[CG.NAMES.index("grass_pasture")] = 256
    return S, B, dict(cdl_n=cdl)


def test_each_rule_fails_alone():
    S, B, Cn = _blocks(10)
    S["wtd_min"][0, 1] = 20.0                                        # WTD
    S["hand_frac"][0, 2] = 0.5                                       # HAND
    S["nwi_wetland_px"][0, 3] = 10                                   # NWI
    B["water_px"][6, 0, 4] = 5                                       # WATER
    S["lanid_dist_m"][0, 5] = 300.0                                  # IRR
    Cn["cdl_n"][:, 0, 6] = 0
    Cn["cdl_n"][CG.NAMES.index("developed"), 0, 6] = 256             # CLASS (never eligible)
    B["ET_n"][:, 0, 7] = 10                                          # DATA
    S["lanid_unassessed_dist_m"][0, 8] = 500.0                       # IRR (LANID does not cover the buffer)
    S["canal_dist_m"][0, 9] = 400.0                                  # CANAL
    valid = np.ones((1, 10), bool)
    E = R.evaluate(S, B, Cn, valid, valid, 2012, {})
    expect = [None, "WTD", "HAND", "NWI", "WATER", "IRR", "CLASS", "DATA", "IRR", "CANAL"]
    for k, rule in enumerate(expect):
        fails = R.explain(E["flags"][0, k])
        assert fails == ([] if rule is None else [rule]), (k, fails)
    assert list(E["reference"][0]) == [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]


def test_wtd_sensitivity_bits_and_gde_is_not_a_main_rule():
    S, B, Cn = _blocks(3)
    S["wtd_min"][0] = [20.0, 40.0, 70.0]
    S["gde_max"][0, 2] = 0.9
    v = np.ones((1, 3), bool)
    E = R.evaluate(S, B, Cn, v, v, 2012, {})
    f = E["flags"][0]
    assert [(int(x) >> R.BITS["WTD15"]) & 1 for x in f] == [1, 1, 1]
    assert [(int(x) >> R.BITS["WTD60"]) & 1 for x in f] == [0, 0, 1]
    assert list(E["reference"][0]) == [0, 1, 1]
    assert list(E["reference_gde"][0]) == [0, 1, 0]


def test_crop_needs_county_dryland_eligibility_and_purity():
    S, B, Cn = _blocks(3)
    ww = CG.NAMES.index("winter_wheat")
    Cn["cdl_n"][:] = 0
    Cn["cdl_n"][ww] = 256
    Cn["cdl_n"][ww, 0, 2] = 150
    Cn["cdl_n"][CG.NAMES.index("corn"), 0, 2] = 106                  # purity 0.59
    S["county"][0] = [0, 1, 0]
    dry = {(0, 2013, ww): True, (1, 2013, ww): False}
    v = np.ones((1, 3), bool)
    E = R.evaluate(S, B, Cn, v, v, 2013, dry)
    assert list(E["reference"][0]) == [1, 0, 0]
    assert R.run_on_set(E["flags"], E["group"]).sum() == 0


def test_run_on_set_is_crop_failing_only_hand():
    S, B, Cn = _blocks(2)
    ww = CG.NAMES.index("winter_wheat")
    Cn["cdl_n"][:] = 0
    Cn["cdl_n"][ww] = 256
    S["hand_frac"][0] = [0.2, 0.2]
    S["wtd_min"][0, 1] = 5.0
    v = np.ones((1, 2), bool)
    E = R.evaluate(S, B, Cn, v, v, 2013, {(0, 2013, ww): True})
    assert list(R.run_on_set(E["flags"], E["group"])[0]) == [True, False]


def test_main_mask_matches_doc_and_registry_slot_refuses():
    assert R.MAIN_MASK == sum(1 << R.BITS[b] for b in R.MAIN)
    import copy
    old = copy.deepcopy(C.IRRIGATION_REGISTRIES)
    C.IRRIGATION_REGISTRIES = [dict(name="KS_WIMAS")]
    try:
        S, B, Cn = _blocks(1)
        v = np.ones((1, 1), bool)
        R.evaluate(S, B, Cn, v, v, 2012, {})
        raise AssertionError("a configured registry without loader was ignored")
    except NotImplementedError:
        pass
    finally:
        C.IRRIGATION_REGISTRIES = old


# ---------------------------------------------------------------- plan v3: depth ladder, ring, sensitivity sets
def _with_group(Cn, name, cols, woody=0):
    for c in cols:
        Cn["cdl_n"][:, 0, c] = 0
        Cn["cdl_n"][CG.NAMES.index(name), 0, c] = 256 - woody
        Cn["cdl_n"][CG.NAMES.index("shrub"), 0, c] += woody


def test_class_rung_makes_a_ring_until_the_stratum_is_accepted():
    from .util import override
    S, B, Cn = _blocks(5)
    S["wtd_min"][0] = [20.0, 20.0, 20.0, 20.0, 12.0]
    _with_group(Cn, "grass_pasture", [0, 1, 4])                         # herbaceous: class rung 15 m
    _with_group(Cn, "grass_pasture", [2], woody=5)                       # woody cover: 30 m
    _with_group(Cn, "hay_alfalfa", [3])                                   # hay: 30 m
    dry = {(0, 2012, CG.NAMES.index("hay_alfalfa")): True}
    v = np.ones((1, 5), bool)
    dec = np.array([[5, 6, 5, 5, 5]])
    E = R.evaluate(S, B, Cn, v, v, 2012, dry, decile=dec)
    assert list(E["depth_class_m"][0]) == [15, 15, 30, 30, 15]
    assert list(E["depth_used_m"][0]) == [30, 30, 30, 30, 30]             # nothing accepted: strict rung
    assert list(E["reference"][0]) == [0, 0, 0, 0, 0]
    assert list(E["reference_ring"][0]) == [1, 1, 0, 0, 0]                # 12 m fails even the class rung
    with override(DEPTH_ACCEPTED={"perennial_grass|5": 15.0}):
        E2 = R.evaluate(S, B, Cn, v, v, 2012, dry, decile=dec)
    assert list(E2["reference"][0]) == [1, 0, 0, 0, 0]                    # only decile 5, never woody / hay
    assert list(E2["reference_ring"][0]) == [0, 1, 0, 0, 0]
    assert list(E2["depth_used_m"][0]) == [15, 30, 30, 30, 15]
    f = E2["flags"][0]
    assert [(int(x) >> R.BITS["WTD10"]) & 1 for x in f] == [1, 1, 1, 1, 1]
    assert [(int(x) >> R.BITS["WTD30"]) & 1 for x in f] == [0, 0, 0, 0, 0]


def test_malformed_depth_acceptance_is_refused():
    from .util import override
    S, B, Cn = _blocks(1)
    v = np.ones((1, 1), bool)
    for bad in ({"grass|5": 15.0}, {"perennial_grass|11": 15.0}, {"perennial_grass|5": 20.0}, {"none|3": 15.0}):
        with override(DEPTH_ACCEPTED=bad):
            assert R.check_depth_accepted(), bad
            try:
                R.evaluate(S, B, Cn, v, v, 2012, {})
                raise AssertionError(f"accepted {bad}")
            except ValueError:
                pass


def test_irr3_sensitivity_ignores_occasional_labels_and_fire_slot_refuses():
    from .util import override
    S, B, Cn = _blocks(2)
    S["lanid_dist_m"][0, 1] = 300.0                                       # a pixel labelled once or twice nearby
    v = np.ones((1, 2), bool)
    E = R.evaluate(S, B, Cn, v, v, 2012, {})
    assert list(E["reference"][0]) == [1, 0] and list(E["reference_irr3"][0]) == [1, 1]
    with override(FIRE_PERIMETERS="mtbs.shp"):
        try:
            R.evaluate(S, B, Cn, v, v, 2012, {})
            raise AssertionError("a configured fire layer without loader was ignored")
        except NotImplementedError:
            pass


def test_functional_groups_cover_every_cdl_group():
    assert set(CG.FUNCTIONAL_OF) == set(CG.NAMES)
    assert set(C.DEPTH_BY_FUNCTIONAL) | {"none"} == set(CG.FUNC_NAMES)
    for n in CG.NAMES:
        if CG.ELIGIBILITY[n] == "never":
            assert CG.FUNCTIONAL_OF[n] == "none", n
