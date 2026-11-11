"""Tests of M1 and of the P1.2 stages (need torch: Colab, or locally the conda env SSM)."""
import contextlib
import os
import tempfile

import numpy as np
import torch

from .. import config as C
from ..p1_2 import m1 as M1, calendars as CA
from .util import synthetic, stage_dir, patched


def _toy(n=64, seed=0):
    rng = np.random.default_rng(seed)
    X = dict(dyn=rng.normal(size=(n, 24, 9)).astype("float32"), stat=rng.normal(size=(n, 12)).astype("float32"),
             pm=rng.normal(size=(n, 12, 7)).astype("float32"), y_std=rng.normal(size=(n, 12, 3)).astype("float32"),
             ymask=rng.random((n, 12, 3)) < 0.8, bound=np.full(n, 400.0, "float32"))
    X["y_std"] += X["dyn"][:, 12:, :3]
    return X


def test_quantiles_never_cross():
    m = M1.make(9, 12, 7, dict(C.M1))
    X = _toy()
    q = m(*(torch.as_tensor(X[k]) for k in ("dyn", "stat", "pm")))
    assert q.shape == (64, 12, 3, 3)
    assert (q[..., 0] <= q[..., 1]).all() and (q[..., 1] <= q[..., 2]).all()


def test_loss_ignores_masked_targets():
    X = _toy()
    q = torch.zeros(64, 12, 3, 3)
    y, m = torch.as_tensor(X["y_std"]), torch.as_tensor(X["ymask"])
    a = M1.pinball(q, y, m)
    y2 = y.clone()
    y2[~m] += 100.0
    assert torch.allclose(a, M1.pinball(q, y2, m))


def test_mass_penalty_acts_only_above_the_bound():
    q = torch.zeros(4, 12, 3, 3)
    q[..., 0, 1] = torch.tensor([0.0, 1.0, 2.0, 3.0])[:, None]      # ET median in sd units
    bound = torch.full((4,), 12 * 10.0 * 1.0)
    p_low = M1.mass_penalty(q[:1], bound[:1], et_mean=10.0, et_sd=1.0)
    p_high = M1.mass_penalty(q[3:], bound[3:], et_mean=10.0, et_sd=1.0)
    assert float(p_low) == 0.0 and float(p_high) > 0.0


def test_resume_gives_the_same_model():
    X, V = _toy(seed=1), _toy(seed=2)
    sc = {"y": (np.zeros(3, "float32"), np.ones(3, "float32"))}
    s = dict(C.M1, max_epochs=2, patience=5, batch=16, amp=False)
    with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
        m1, _ = M1.train(M1.make(9, 12, 7, s), X, V, sc, s, a, device="cpu")
        s1 = dict(s, max_epochs=1)
        M1.train(M1.make(9, 12, 7, s), X, V, sc, s1, b, device="cpu")
        m2, _ = M1.train(M1.make(9, 12, 7, s), X, V, sc, s, b, device="cpu")       # resumes epoch 2
        for p, q in zip(m1.state_dict().values(), m2.state_dict().values()):
            assert torch.allclose(p, q, atol=1e-6)


def _stand_in_season_window(curve):
    """Used ONLY by this test while calendars.season_window is not written: months above the mid-range."""
    c = np.asarray(curve, "float64")
    if not np.isfinite(c).any():
        return 4, 7
    c = np.where(np.isfinite(c), c, np.nanmin(c))
    on = c >= c.min() + 0.5 * (c.max() - c.min())
    start = int(np.argmax(on & ~np.roll(on, 1))) + 1 if (on & ~np.roll(on, 1)).any() else 1
    return start, int(on.sum())


def test_e2e_train_predict_verify_on_synthetic():
    with synthetic():
        if not stage_dir("folds"):
            return "SKIP run `python -m paper1 test e2e` (system Python / pysheds venv) first"
        ctx = contextlib.nullcontext()
        try:
            CA.season_window(np.linspace(0, 1, 12))
        except NotImplementedError:
            ctx = patched(CA, "season_window", _stand_in_season_window)
            print("      (calendars.season_window not written yet: a stand-in is used for this test only)")
        with ctx:
            from ..pipeline import runner
            runner().run(["verify_p12"])
        d = stage_dir("verify_p12")
        assert d and os.path.exists(f"{d}/reports/verify_p12.md")
        import json
        rep = json.load(open(f"{d}/reports/verify_p12.json"))
        for t in C.TARGETS:
            assert rep["skill"][t]["n"] > 0
            assert 0.0 <= rep["skill"][t]["coverage"] <= 1.0
        pdir = stage_dir("predict")
        assert any(f.endswith(".nc") for _, _, fs in os.walk(pdir) for f in fs)
