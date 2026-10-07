"""Tests of the differentiable model (conda env SSM: numpy + torch)."""
import copy
import gc
import tempfile
from pathlib import Path

import numpy as np
import torch

from .. import config as C
from .. import splits as SP
from .. import train_dlstm as T
from ..hbv import HBV
from . import synth

SMALL = dict(hidden=8, layers=1, win=420, warm=200, batch_cells=6, windows_per_step=2, val_cells=6)


class _small_cfg:
    """Shrinks the network and the windows for the duration of a test."""

    def __enter__(self):
        self.keep = copy.deepcopy(C.DLSTM)
        C.DLSTM.update(SMALL)

    def __exit__(self, *exc):
        C.DLSTM.clear()
        C.DLSTM.update(self.keep)


class _tmp_profile:
    def __init__(self, **kw):
        self.kw = kw

    def __enter__(self):
        self.tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.keep = C.STORE
        C.STORE = Path(self.tmp.name)
        C.PROFILES["_test"] = C.Profile("_test", None)
        self.store, self.folds = synth.make_store(C.STORE / "_test", **self.kw)
        return self

    def __exit__(self, *exc):
        C.STORE = self.keep
        C.PROFILES.pop("_test", None)
        gc.collect()
        self.tmp.cleanup()


def test_hbv_conserves_mass():
    torch.manual_seed(0)
    hbv = HBV(C.HBV_BOUNDS)
    B, n = 16, 500
    P = torch.rand(B, n) * 12 * (torch.rand(B, n) > 0.7)
    T = 10 + 15 * torch.sin(torch.arange(n) / 58.0).repeat(B, 1) + torch.randn(B, n)
    Ep = 3 + 2 * torch.rand(B, n)
    out = hbv(P, T, Ep, hbv.scale(torch.randn(B, len(C.HBV_BOUNDS))))
    assert out["mb_resid"].abs().max() < 0.05, float(out["mb_resid"].abs().max())
    assert (out["ET"] >= 0).all() and (out["ET"] <= Ep + 1e-5).all()


def test_static_inputs_hold_no_exposure_variable():
    with _small_cfg(), _tmp_profile(n_cells=12) as p:
        data = T.Data(p.store, 0)
        assert data.A.shape[1] == len(C.DLSTM_STATIC)
        for bad in C.EXPOSURE_VARS:
            assert bad not in C.DLSTM_STATIC
        del data


def test_theta_is_the_same_for_every_period():
    """L5: parameters come from fixed inputs, so no window can change them."""
    with _small_cfg(), _tmp_profile(n_cells=12) as p:
        torch.manual_seed(1)
        data = T.Data(p.store, 0)
        model = T.Model(data.A.shape[1], torch.device("cpu"))
        model.gA.eval()
        cells = np.arange(12)
        with torch.no_grad():
            th = model.theta(data, cells)
            a = model.simulate(data, cells, 0, 900, th)["ET"]
            b = model.simulate(data, cells, 0, 1400, model.theta(data, cells))["ET"]
        assert torch.equal(th, model.theta(data, cells).detach())
        assert torch.allclose(a, b[:, :900], atol=1e-5), "a longer run changed the earlier days"
        del data, model


def test_soil_moisture_loss_sees_a_dry_year():
    """A2: lowering observed soil moisture in the loss window must change the loss.
    The old in-window standardisation removed exactly this signal."""
    with _small_cfg(), _tmp_profile(n_cells=12) as p:
        torch.manual_seed(2)
        data = T.Data(p.store, 0)
        model = T.Model(data.A.shape[1], torch.device("cpu"))
        cells = np.arange(12)
        i0 = int(T.valid_starts(data, C.DLSTM["win"], C.DLSTM["warm"])[40])
        lo, hi = i0 + C.DLSTM["warm"], i0 + C.DLSTM["win"]
        with torch.no_grad():
            th = model.theta(data, cells)
            out = model.simulate(data, cells, i0, hi, th)
            base = T.sums(out, th, model, data, cells, i0, lo, hi, SP.TRAIN)
            data.SM = np.array(data.SM) - 0.06
            dry = T.sums(out, th, model, data, cells, i0, lo, hi, SP.TRAIN)
        assert base[3] > 0 and abs(float(dry[2]) - float(base[2])) > 1e-3 * max(float(base[2]), 1e-6)
        del data, model


def test_loss_uses_only_dekads_of_its_role():
    with _small_cfg(), _tmp_profile(n_cells=12) as p:
        data = T.Data(p.store, 0)
        model = T.Model(data.A.shape[1], torch.device("cpu"))
        cells = np.arange(12)
        with torch.no_grad():
            th = model.theta(data, cells)
            out = model.simulate(data, cells, 0, data.n_days, th)
            n = {r: T.sums(out, th, model, data, cells, 0, 0, data.n_days, r)[1] for r in
                 (SP.TRAIN, SP.VAL, SP.TEST)}
        for role, code in (("train", SP.TRAIN), ("val", SP.VAL), ("test", SP.TEST)):
            assert n[code] == int((data.split == code).sum()) * 12, role
        before = T.sums(out, th, model, data, cells, 0, 0, data.n_days, SP.TRAIN)[0]
        data.et = np.array(data.et)
        data.et[data.split == SP.TEST] += 100.0                       # corrupt the test targets
        after = T.sums(out, th, model, data, cells, 0, 0, data.n_days, SP.TRAIN)[0]
        assert torch.equal(before, after), "test targets reached the training loss"
        del data, model


def _weights(path):
    return torch.load(path, map_location="cpu", weights_only=False)["gA"]


def test_resume_continues_from_the_last_epoch():
    """B3: 2 epochs in one go == 1 epoch + resume for 1 epoch, weight for weight."""
    with _small_cfg(), _tmp_profile(n_cells=12) as p:
        T.train("_test", 0, member=0, epochs=2, device="cpu", steps=2)
        full = _weights(C.STORE / "_test" / "dlstm" / "fold0" / "m0" / "last.pt")
        T.train("_test", 0, member=1, epochs=1, device="cpu", steps=2)
        d = C.STORE / "_test" / "dlstm" / "fold0"
        # member 1 uses another seed; rerun member 0 in two parts instead
        import shutil
        shutil.rmtree(d / "m0")
        T.train("_test", 0, member=0, epochs=1, device="cpu", steps=2)
        T.train("_test", 0, member=0, epochs=1, device="cpu", steps=2, resume=True)
        parts = _weights(d / "m0" / "last.pt")
        for k in full:
            assert torch.allclose(full[k], parts[k], atol=1e-6), f"{k} differs after resume"
        log = (d / "m0" / "train_log.csv").read_text().strip().splitlines()
        assert [ln.split(",")[0] for ln in log[1:]] == ["0", "1"], "epoch numbering restarted or repeated"
        other = _weights(d / "m1" / "last.pt")
        assert any(not torch.allclose(full[k], other[k]) for k in full), "ensemble members are identical"


def test_inference_writes_test_residuals_only():
    with _small_cfg(), _tmp_profile(n_cells=12) as p:
        T.train("_test", 0, member=0, epochs=1, device="cpu", steps=2)
        res = T.infer("_test", 0, members=(0,), device="cpu")
        d = C.STORE / "_test" / "dlstm"
        bi, fold = np.load(d / "Bi.npy"), np.load(d / "oof_fold.npy")
        te = p.folds["split"][0] == SP.TEST
        assert (fold[te] == 0).all() and (fold[~te] == -1).all()
        assert np.isfinite(bi[te]).all() and np.isnan(bi[~te]).all()
        assert np.isfinite(res["test_rmse"]) and 0 <= res["coverage90"] <= 1
