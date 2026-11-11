"""End-to-end run of every stage up to P1.1 (+ samples, folds) on the synthetic data root.

Downloads are replaced by local synthetic files (counties, NHDPlus HR subset, DEM tiles); everything else -
reading, aggregation, rules, provenance, resume - is the real code. The provenance tests are then run
against the synthetic outputs, so they are themselves tested.
"""
import json
import os

import numpy as np

from .. import config as C
from .util import synthetic, stage_dir

STAGES = ["verify_p11", "samples", "folds"]


def _run():
    from ..pipeline import runner
    runner().run(STAGES)


def test_e2e_runs_step0_p11_samples_folds():
    with synthetic():
        _run()
        for s in ["tiles", "blocks", "drivers", "counties", "hand", "static", "classes", "dryland", "reference",
                  "verify_p11", "samples", "folds"]:
            assert stage_dir(s), f"stage {s} not complete"


def test_e2e_rules_exclude_what_the_synthetic_landscape_puts_there():
    from ..step0 import tiles as TL, static as ST
    from ..p1_1 import build as RB, rules as RU
    from . import synth
    with synthetic():
        _run()
        T, G = TL.load_table(stage_dir("tiles"))
        S = ST.open_tile(stage_dir("static"), synth.TILE)
        R = RB.open_tile(stage_dir("reference"), synth.TILE)
        f = R["flags"][0]
        col_valley = (synth.C0 + 60) // 16 - G[synth.TILE].bc0          # block column holding the valley
        bit = lambda name: ((f >> RU.BITS[name]) & 1).astype(bool)
        assert not bit("WTD")[:, col_valley].any(), "valley blocks must fail the depth-to-water rule"
        assert not bit("NWI")[:, col_valley].any(), "valley blocks hold riparian NWI pixels"
        assert S["lanid_dist_m"].min() == 0 and not bit("IRR").all(), "blocks on the irrigated corn fail IRR"
        assert R["reference"].sum() > 0, "some blocks must qualify"
        rep = json.load(open(f"{stage_dir('verify_p11')}/reports/verify_p11.json"))
        assert rep["reference"]["V1"]["total_distinct"] > 0
        assert os.path.exists(f"{stage_dir('verify_p11')}/reports/reference_2012.png")


def test_e2e_rerun_recomputes_nothing():
    from ..core import provenance as P
    with synthetic():
        _run()
        led = P.Ledger(f"{C.out_root()}/ledger.jsonl")
        n0 = len(led.read())
        f = sorted(os.listdir(f"{stage_dir('blocks')}/2012"))[0]
        p = f"{stage_dir('blocks')}/2012/{f}/P1_blocks_{f}_2012.nc"
        m0 = os.path.getmtime(p)
        _run()
        assert os.path.getmtime(p) == m0
        assert len(led.read()) == n0                   # complete stages are skipped without ledger events


def test_e2e_a_changed_setting_gets_a_new_folder():
    from .util import override
    from ..pipeline import runner
    with synthetic():
        _run()
        old = runner().sig("reference")
        with override(WTD_MIN_M=31.0):
            assert runner().sig("reference") != old
            assert runner().sig("blocks") == runner().sig("blocks")
        assert runner().sig("reference") == old


def test_e2e_provenance_tests_pass_on_synthetic_outputs():
    from . import test_provenance as TP
    with synthetic():
        _run()
        ran = 0
        for name in [n for n in dir(TP) if n.startswith("test_")]:
            out = getattr(TP, name)()
            if isinstance(out, str) and out.startswith("SKIP"):
                assert "pysheds" in out or "D8" in out, f"{name} skipped on synthetic data: {out}"
            else:
                ran += 1
        assert ran >= 8
