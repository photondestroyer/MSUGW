"""Unit tests of the core layer (synthetic, seconds)."""
import ast
import os
import tempfile
import time

import numpy as np

from .. import config as C
from ..core import provenance as P, ncio
from .util import override


def test_vendored_copies_agree_with_residual_et():
    try:
        from residual_et import sources as RS, uq as RU
    except Exception as e:                              # residual_et is not deployed to Colab
        return f"SKIP residual_et not importable ({type(e).__name__})"
    from ..p1_2 import baselines as BL
    import netCDF4 as nc
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "x.nc")
        with nc.Dataset(p, "w") as ds:
            ds.createDimension("n", 5)
            v = ds.createVariable("v", "i2", ("n",), fill_value=32767)
            v.scale_factor, v.add_offset = np.float32(0.1), np.float32(210)
            v.set_auto_maskandscale(False)
            v[:] = np.array([0, 1, 32767, -5, 100], "i2")
        with ncio.open_ro(p) as ds:
            raw = np.asarray(ds["v"][:])
            a = ncio.decode(raw, ds["v"])
            b = RS.decode(raw, ds["v"])
    assert np.array_equal(a, b, equal_nan=True)
    rng = np.random.default_rng(0)
    y, lo, hi = rng.normal(size=500), rng.normal(-1, 0.1, 500), rng.normal(1, 0.1, 500)
    assert BL.cqr_margin(y, lo, hi) == RU.cqr_margin(y, lo, hi)
    assert BL.conformal_quantile(y) == RU.conformal_quantile(y)


def test_fingerprint_identity_ignores_mtime_but_not_content():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "f.bin")
        data = np.random.default_rng(1).bytes(6 * 1024 * 1024)
        open(p, "wb").write(data)
        with override(LOCAL_SCRATCH=d, RUN_ENV="local"):
            a = P.fingerprint(p, full=True)
            os.utime(p, (time.time() - 1000, time.time() - 1000))
            b = P.fingerprint(p, full=True)
            assert a["mtime"] != b["mtime"] and P.same_file(a, b)
            open(p, "r+b").write(b"X")
            c = P.fingerprint(p, full=False)
            assert not P.same_file(a, c)


def test_signature_follows_config_and_upstream():
    a = P.signature("x", ["WTD_MIN_M"], ["config.py"], {"up": "11111111"})
    with override(WTD_MIN_M=31.0):
        b = P.signature("x", ["WTD_MIN_M"], ["config.py"], {"up": "11111111"})
    c = P.signature("x", ["WTD_MIN_M"], ["config.py"], {"up": "22222222"})
    assert a != b and a != c and len(a) == 8


def test_write_nc_reads_back_and_marks_complete():
    with tempfile.TemporaryDirectory() as d:
        with override(RUN_ENV="local", LOCAL_SCRATCH=d):
            p = os.path.join(d, "out", "a.nc")
            arr = np.arange(12, dtype="float32").reshape(3, 4)
            arr[0, 0] = np.nan
            ncio.write_nc(p, dict(y=3, x=4), dict(v=dict(dims=("y", "x"), data=arr, dtype="f4")), dict(t="1"),
                          dict(stage="s", signature="abcd1234"), os.path.join(d, "s"))
            assert ncio.is_complete(p, "abcd1234") and not ncio.is_complete(p, "ffffffff")
            with ncio.open_ro(p) as ds:
                assert np.array_equal(np.asarray(ds["v"][:]), arr, equal_nan=True)


def test_read_window_pads_outside_the_file():
    class V:
        shape = (4, 5)
        dtype = np.dtype("int16")

        def __getitem__(self, sl):
            return np.arange(20, dtype="int16").reshape(4, 5)[sl]
    out = ncio.read_window(V(), (-2, 3), (3, 7), -1)
    assert out.shape == (5, 4)
    assert (out[:2] == -1).all() and (out[:, 2:] == -1).all()
    assert out[2, 0] == 3 and out[4, 1] == 14


def test_config_snapshot_is_json_and_check_refuses_bad_modes():
    import json
    json.dumps(P._jsonable(C.snapshot()))
    with override(RUN_MODE="nope"):
        try:
            C.check()
            raise AssertionError("bad RUN_MODE accepted")
        except ValueError:
            pass


def test_oracle_is_independent_of_the_pipeline():
    """Value-tracing oracles must not import paper1 code (they would share its bugs)."""
    src = open(os.path.join(os.path.dirname(__file__), "oracle.py"), encoding="utf-8").read()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.ImportFrom):
            assert node.level == 0 and not (node.module or "").startswith("paper1"), f"oracle imports {node.module}"
        if isinstance(node, ast.Import):
            assert not any(a.name.startswith("paper1") for a in node.names)


def test_config_md_is_current():
    from ..docs_gen import render_config_md
    p = os.path.join(os.path.dirname(os.path.dirname(__file__)), "docs", "CONFIG.md")
    if not os.path.exists(p):
        return "SKIP docs/CONFIG.md not generated yet (python -m paper1 docs)"
    assert open(p, encoding="utf-8").read() == render_config_md(), "docs/CONFIG.md is stale: run python -m paper1 docs"


def test_a_failed_gate_halts_dependents_in_full_mode_unless_overridden():
    import json
    from ..core.stages import Runner
    from ..core.errors import PipelineStop
    root = tempfile.mkdtemp()
    stages = [dict(name="a", deps=[], keys=[], modules=[], func="x:y", gate="G1a"),
              dict(name="b", deps=["a"], gated_by=["a"], keys=[], modules=[], func="x:y")]
    with override(RUN_ENV="local", LOCAL_DATA_ROOT=root, LOCAL_SCRATCH=f"{root}/s", GATE_OVERRIDE=[]):
        for mode in ("full", "test"):
            with override(RUN_MODE=mode):
                r = Runner(stages)
                d = f"{C.out_root()}/a/{r.sig('a')}"
                os.makedirs(d, exist_ok=True)
                with open(f"{d}/_complete.json", "w") as f:
                    json.dump(dict(summary=dict(gate=dict(name="G1a", passes=False, reasons=["V2 fails"]))), f)
                st = r.gate_states(r.stages["b"])
                assert st["G1a"]["passes"] is False and st["G1a"]["overridden"] is False
                if mode == "full":
                    try:
                        r.check_gates(r.stages["b"])
                        raise AssertionError("a failed gate did not halt in full mode")
                    except PipelineStop as e:
                        assert "G1a" in str(e) and "GATE_OVERRIDE" in str(e)
                    with override(GATE_OVERRIDE=["G1a"]):
                        r.check_gates(r.stages["b"])                    # overridden: runs
                        assert r.gate_states(r.stages["b"])["G1a"]["overridden"]
                else:
                    r.check_gates(r.stages["b"])                        # test mode: reported only
