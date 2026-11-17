"""The blinding rule: the gate keeps rainfed and non-dry cohort ET only; the audit catches leaks."""
import os
import shutil
import tempfile

import numpy as np
import pandas as pd


def _fixture():
    cells = pd.DataFrame(dict(cell_id=["c1", "c2", "c3"], county=["20001", "20001", "20003"],
                              role=["cohort", "rainfed", "excluded"]))
    dry = pd.DataFrame(dict(county=["20001", "20001", "20003"], year=[2020, 2020, 2020], month=[7, 8, 7],
                            dry=[True, False, False]))
    et = pd.DataFrame(dict(cell=["c1", "c1", "c2", "c2", "c3"], year=2020, month=[7, 8, 7, 8, 7],
                           et_ensemble_mad=[100.0, 90.0, 60.0, 50.0, 10.0]))
    return cells, dry, et


def test_gate_keeps_only_what_is_allowed():
    from pilot import blinding
    cells, dry, et = _fixture()
    out, dropped = blinding.gate(et, cells, dry)
    got = set(zip(out.cell, out.month))
    assert got == {("c1", 8), ("c2", 7), ("c2", 8)}       # cohort dry July dropped; excluded cell dropped
    assert dropped["cohort_dry_dropped"] == 1 and dropped["other_dropped"] == 1
    assert not ((out.role == "cohort") & out.dry).any()


def test_gate_unknown_dry_is_dry():
    from pilot import blinding
    cells, dry, et = _fixture()
    out, _ = blinding.gate(et, cells, dry[dry.month != 8])
    assert ("c1", 8) not in set(zip(out.cell, out.month))  # no flag -> treated as dry -> hidden


def test_gate_refuses_thickness_and_other_states():
    from pilot import blinding
    cells, dry, et = _fixture()
    try:
        blinding.gate(et.assign(A_m=50.0), cells, dry)
        raise AssertionError("thickness next to ET must be refused")
    except blinding.BlindingError:
        pass
    try:
        blinding.gate(et, cells.assign(county=["48001", "20001", "20003"]), dry)
        raise AssertionError("replication-region cells must be refused")
    except blinding.BlindingError:
        pass


def test_audit_catches_a_leak():
    from pilot import blinding
    cells, dry, et = _fixture()
    d = tempfile.mkdtemp()
    try:
        good, _ = blinding.gate(et, cells, dry)
        good.to_parquet(os.path.join(d, "good.parquet"))
        ok, problems, n = blinding.audit(d)
        assert ok and n == 1, problems
        leak = et.merge(cells, left_on="cell", right_on="cell_id").merge(dry, on=["county", "year", "month"], how="left")
        leak.to_csv(os.path.join(d, "leak.csv"), index=False)
        ok, problems, n = blinding.audit(d)
        assert not ok and any("cohort rows in dry" in p for p in problems)
        pd.DataFrame(dict(cell=["c1"], et_ensemble_mad=[1.0])).to_csv(os.path.join(d, "raw.csv"), index=False)
        pd.DataFrame(dict(file=["x"], et_hits=[3])).to_csv(os.path.join(d, "not_et.csv"), index=False)
        ok, problems, _ = blinding.audit(d)
        assert any("did not pass the gate" in p for p in problems)
        assert not any("not_et.csv" in p for p in problems)          # a count column named et_* is not ET
    finally:
        shutil.rmtree(d)


def test_dry_table():
    from pilot import blinding
    from pilot import config as C
    w = pd.DataFrame(dict(county=["20001"] * 9, year=C.PANEL_YEARS, deficit_season_mm=np.arange(9) * 100.0))
    s = pd.DataFrame(dict(county="20001", year=2016, month=[5, 6], S=[-2.0, 2.0]))
    out, thr = blinding.dry_table(w, s, dict(season_pct=70, s_pct=30, missing_s_is_dry=True))
    assert out.loc[(out.year == 2024)].dry_season.all()                 # high deficit -> dry season
    assert out.loc[(out.year == 2016) & (out.month == 7)].dry.all()     # missing S -> dry month


# ------------------------------------------------------------------------------------------- the freeze
def _decision_log_copy(make_fixed):
    """A temporary copy of the decision log, every entry fixed (or the first fixed one made proposed)."""
    from pilot import decisions as D
    d = tempfile.mkdtemp()
    p = os.path.join(d, "decisions.yaml")
    with open(D.PATH, encoding="utf-8") as f:
        txt = f.read()
    txt = txt.replace("status: proposed", "status: fixed") if make_fixed else txt.replace("status: fixed",
                                                                                          "status: proposed", 1)
    with open(p, "w", encoding="utf-8") as f:
        f.write(txt)
    return d, p


def _stops(fn, *a, **k):
    from pilot.core.errors import PipelineStop
    try:
        fn(*a, **k)
    except PipelineStop:
        return True
    return False


def test_freeze_record_needs_fixed_decisions_and_detects_changes():
    import json
    from pilot import freeze
    d, p = _decision_log_copy(make_fixed=False)
    try:
        assert _stops(freeze.make_record, "me", "sup", "abc", decisions_path=p)      # a decision still proposed
    finally:
        shutil.rmtree(d)
    d, p = _decision_log_copy(make_fixed=True)
    try:
        rec = freeze.make_record("me", "sup", "abc", decisions_path=p, today="2026-12-01")
        assert rec["tag"] == "pilot-freeze-2026-12-01" and rec["history"] == []
        fp = os.path.join(d, "freeze.json")
        with open(fp, "w") as f:
            json.dump(rec, f)
        assert freeze.check(fp, p)["tag"] == rec["tag"]
        assert _stops(freeze.make_record, "me", "sup", "abd", previous=rec, decisions_path=p)   # no deviation
        dev = freeze.make_record("me", "sup", "abd", previous=rec, deviation="why", decisions_path=p,
                                 today="2026-12-02")
        assert dev["tag"].endswith("-dev1") and dev["history"][0]["tag"] == rec["tag"]
        with open(p, "a", encoding="utf-8") as f:                                    # log edited after the freeze
            f.write("\n# edited\n")
        assert _stops(freeze.check, fp, p)
        assert _stops(freeze.check, os.path.join(d, "missing.json"), p)
    finally:
        shutil.rmtree(d)


def test_gate_unblinds_only_in_main_mode_after_a_freeze():
    from pilot import blinding, config as C, freeze
    cells, dry, et = _fixture()
    mode, chk = C.RUN_MODE, freeze.check
    try:
        try:
            blinding.gate(et, cells, dry, unblind=True)                              # not main mode
            assert False, "unblinded outside main mode"
        except blinding.BlindingError:
            pass
        C.RUN_MODE = "main"
        freeze.check = lambda *a, **k: (_ for _ in ()).throw(blinding.PipelineStop("no freeze"))
        assert _stops(blinding.gate, et, cells, dry, unblind=True)                   # main mode, no record
        freeze.check = lambda *a, **k: dict(tag="pilot-freeze-x")
        out, dropped = blinding.gate(et, cells, dry, unblind=True)
        assert set(zip(out.cell, out.month)) == {("c1", 7), ("c1", 8), ("c2", 7), ("c2", 8)}
        assert dropped["cohort_dry_dropped"] == 0 and out.blinding_gate.str.startswith(blinding.AFTER_FREEZE).all()
    finally:
        C.RUN_MODE, freeze.check = mode, chk
