"""Differentiable LSTM-parameterised HBV: train, resume, infer (conda env SSM).

    python -m residual_et.train_dlstm --profile full --fold 0               # train
    python -m residual_et.train_dlstm --profile full --fold 0 --resume      # continue
    python -m residual_et.train_dlstm --profile full --fold 0 --infer       # residuals

Corrections against delta_gwb/train_full.py:
  L4  static inputs hold no water-table depth and no irrigation descriptor;
  L5  theta comes from fixed inputs (attributes + day-of-year normals), so it is
      the same in training, validation and test;
  A2  soil moisture is compared as relative saturation with a fixed per-cell
      range from training days, not standardised inside each window;
  A3  the ET anomaly SD has a 3 mm floor and a 12-year baseline;
  A4  the checkpoint is chosen on the quantity that is trained;
  B3  --resume restores the LAST weights, the optimiser and the random state;
  B4  every epoch visits every cell, with several windows per optimiser step.
Reads only the prepared store; needs numpy and torch.
"""
import argparse
import json
import os
import time

import numpy as np
import torch

from . import config as C
from . import climate as K
from . import splits as SP
from .dpl import ParamNet
from .hbv import HBV
from .store import Store

CFG = C.DLSTM
Z90 = 1.6449


class Data:
    """Everything the model needs, for one fold."""

    def __init__(self, store, fold):
        store.require_complete()
        self.days, self.starts = store.days, store.dekad_start
        self.P, self.T, self.EP = (store.arr(n) for n in ("P_daily", "T_daily", "EP_daily"))
        self.SM = store.arr("SM_daily")
        self.n_days, self.n_cells = self.P.shape
        self.et = np.asarray(store.arr("ET"))
        self.dek = np.asarray(store.arr("dek36")).astype(int)
        self.cm = np.asarray(store.arr("ET_clim"))
        self.cs = np.maximum(np.asarray(store.arr("ET_clim_sd")), CFG["et_sd_floor"])
        self.split = np.load(store.path / "folds.npz")["split"][fold]
        self.d0 = (self.starts - self.days[0]).astype(int)                 # first day of each dekad
        self.d1 = (K.dekad_end(self.starts) - self.days[0]).astype(int)    # one past its last day
        self.day_role = np.full(self.n_days, SP.EMBARGO, "int8")           # spin-up days carry no loss
        for k in range(self.starts.size):
            self.day_role[self.d0[k]:self.d1[k]] = self.split[k]
        C.assert_no_exposure(C.DLSTM_STATIC)
        A = np.column_stack([np.asarray(store.arr(n), "float64") for n in C.DLSTM_STATIC])
        mu, sd = np.nanmean(A, axis=0), np.nanstd(A, axis=0)
        A = (A - mu) / np.where(sd > 1e-9, sd, 1.0)
        self.A = np.where(np.isfinite(A), A, 0.0).astype("float32")
        seq = []
        for n in ("P_dclim", "T_dclim", "EP_dclim"):
            x = np.asarray(store.arr(n), "float64")
            seq.append((x - x.mean()) / max(x.std(), 1e-9))
        self.clim_seq = np.stack(seq, axis=-1).astype("float32")           # (366, N, 3), fixed
        self.sm_lo, self.sm_hi = self._sm_range()

    def _sm_range(self, chunk=4096):
        """Per-cell soil-moisture range from TRAINING days only."""
        tr = np.where(self.day_role == SP.TRAIN)[0]
        lo = np.full(self.n_cells, np.nan, "float32")
        hi = np.full(self.n_cells, np.nan, "float32")
        for c0 in range(0, self.n_cells, chunk):
            blk = np.asarray(self.SM[:, c0:c0 + chunk])[tr]
            ok = np.isfinite(blk).sum(axis=0) >= 30
            if ok.any():
                with np.errstate(invalid="ignore"):
                    q = np.nanpercentile(blk[:, ok], CFG["sm_pct"], axis=0)
                lo[c0:c0 + chunk][ok], hi[c0:c0 + chunk][ok] = q[0], q[1]
        return lo, hi

    def forcing(self, cells, i0, i1):
        f = lambda a: torch.from_numpy(np.ascontiguousarray(np.asarray(a[i0:i1])[:, cells].T))
        return f(self.P), f(self.T), f(self.EP)


class Model:
    def __init__(self, n_attr, device, dropout=CFG["dropout"]):
        self.dev = device
        self.hbv = HBV(C.HBV_BOUNDS)
        self.gA = ParamNet(n_attr, len(C.HBV_BOUNDS), CFG["hidden"], CFG["layers"], dropout).to(device)
        self.i_fc = self.hbv.names.index("FC")

    def theta(self, data, cells):
        """Static parameters from fixed inputs; identical for every simulated period."""
        dyn = torch.from_numpy(np.ascontiguousarray(data.clim_seq[:, cells].transpose(1, 0, 2)))
        A = torch.from_numpy(data.A[cells])
        return self.hbv.scale(self.gA(dyn.to(self.dev), A.to(self.dev)).cpu())

    def simulate(self, data, cells, i0, i1, theta):
        return self.hbv(*data.forcing(cells, i0, i1), theta)


def sums(out, theta, model, data, cells, i0, lo, hi, role):
    """Squared-error sums of one simulation that started on day i0, over the
    days [lo, hi) and the dekads lying fully inside them, for one role."""
    ks = np.where((data.d0 >= lo) & (data.d1 <= hi) & (data.split == role))[0]
    et_ss, et_n = torch.zeros(()), 0
    if ks.size:
        cols = (np.concatenate([np.arange(data.d0[k], data.d1[k]) for k in ks]) - i0).astype("int64")
        pos = np.concatenate([np.full(data.d1[k] - data.d0[k], j) for j, k in enumerate(ks)]).astype("int64")
        mod = torch.zeros(len(cells), ks.size).index_add_(1, torch.from_numpy(pos),
                                                          out["ET"][:, torch.from_numpy(cols)])
        obs = torch.from_numpy(data.et[ks][:, cells].T)
        cs = torch.from_numpy(data.cs[data.dek[ks]][:, cells].T)
        v = torch.isfinite(obs) & torch.isfinite(cs)
        if v.any():
            et_ss = (((mod - torch.nan_to_num(obs)) / torch.nan_to_num(cs, nan=1.0))[v] ** 2).sum()
            et_n = int(v.sum())
    dcols = (np.where(data.day_role[lo:hi] == role)[0] + lo).astype("int64")
    sm_ss, sm_n = torch.zeros(()), 0
    if dcols.size:
        # column selection first: only the batch's cells are copied out of the memmap
        obs = np.asarray(data.SM[dcols.min():dcols.max() + 1][:, cells])[dcols - dcols.min()].T
        rng_ = np.maximum(data.sm_hi[cells] - data.sm_lo[cells], 1e-3)[:, None]
        rs = torch.from_numpy(np.clip((obs - data.sm_lo[cells][:, None]) / rng_, 0.0, 1.0).astype("float32"))
        mod = out["Ss"][:, torch.from_numpy(dcols - i0)] / theta[:, model.i_fc].unsqueeze(1)
        v = torch.isfinite(rs)
        if v.any():
            sm_ss = ((mod - torch.nan_to_num(rs))[v] ** 2).sum()
            sm_n = int(v.sum())
    return et_ss, et_n, sm_ss, sm_n


def combine(et_ss, et_n, sm_ss, sm_n):
    """The trained quantity: RMSE of the standardised ET anomaly plus the weighted
    RMSE of relative saturation. Validation uses the same expression."""
    et = torch.sqrt(et_ss / max(et_n, 1) + 1e-8)
    sm = torch.sqrt(sm_ss / max(sm_n, 1) + 1e-8)
    return et + CFG["sm_weight"] * sm, et, sm


def evaluate(model, data, cells, role, batch=1024):
    """Continuous simulation from the first forcing day (two years of spin-up)."""
    model.gA.eval()
    tot = [torch.zeros(()), 0, torch.zeros(()), 0]
    with torch.no_grad():
        for c0 in range(0, len(cells), batch):
            cb = cells[c0:c0 + batch]
            th = model.theta(data, cb)
            out = model.simulate(data, cb, 0, data.n_days, th)
            for j, x in enumerate(sums(out, th, model, data, cb, 0, 0, data.n_days, role)):
                tot[j] = tot[j] + x
    loss, et, sm = combine(*tot)
    return float(loss), float(et), float(sm)


def valid_starts(data, win, warm, min_dekads=3):
    """Window starts whose loss part holds at least `min_dekads` training dekads."""
    tr = data.split == SP.TRAIN
    starts = np.arange(0, data.n_days - win + 1)
    lo, hi = starts + warm, starts + win
    n = np.array([(tr & (data.d0 >= a) & (data.d1 <= b)).sum() for a, b in zip(lo, hi)])
    return starts[n >= min_dekads]


def save_atomic(obj, path):
    tmp = str(path) + ".tmp"
    torch.save(obj, tmp)
    os.replace(tmp, path)


def train(profile_name, fold, member=0, epochs=None, resume=False, device=None, steps=None):
    t0 = time.time()
    profile = C.get_profile(profile_name)
    store = Store(profile.store)
    out = profile.store / "dlstm" / f"fold{fold}" / f"m{member}"
    out.mkdir(parents=True, exist_ok=True)
    seed = C.SEED + 1000 * member + fold
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    dev = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
    data = Data(store, fold)
    model = Model(data.A.shape[1], dev)
    opt = torch.optim.Adam(model.gA.parameters(), lr=CFG["lr"])
    win, warm, B, Kw = CFG["win"], CFG["warm"], CFG["batch_cells"], CFG["windows_per_step"]
    starts = valid_starts(data, win, warm)
    if starts.size == 0:
        raise RuntimeError("no training window holds three training dekads; check the split")
    val_cells = np.random.default_rng(7).permutation(data.n_cells)[:CFG["val_cells"]]
    epoch0, best, logp, lastp = 0, float("inf"), out / "train_log.csv", out / "last.pt"
    header = "epoch,train_loss,train_et,train_sm,val_loss,val_et,val_sm,seconds\n"
    if resume:
        if not lastp.exists():
            raise FileNotFoundError(f"--resume: {lastp} not found")
        ck = torch.load(lastp, map_location=dev, weights_only=False)
        model.gA.load_state_dict(ck["gA"])
        opt.load_state_dict(ck["opt"])
        rng.bit_generator.state = ck["np_rng"]
        torch.set_rng_state(ck["torch_rng"].cpu())      # map_location moved it to the GPU
        epoch0, best = ck["epoch"] + 1, ck["best"]
        if logp.read_text().splitlines()[0] + "\n" != header:
            raise RuntimeError(f"{logp} was written by a different version; start a fresh run")
        print(f"resumed after epoch {ck['epoch']} (best val {best:.4f})", flush=True)
    else:
        logp.write_text(header)
    n_epochs = epochs or CFG["epochs"]
    steps = steps or int(np.ceil(data.n_cells / (B * Kw)))
    print(f"fold {fold} member {member}: {data.n_cells} cells, {starts.size} window starts, "
          f"{steps} steps/epoch, device {dev}", flush=True)
    for ep in range(epoch0, epoch0 + n_epochs):
        model.gA.train()
        perm = rng.permutation(data.n_cells)
        acc, nacc, ptr = np.zeros(3), 0, 0
        for _ in range(steps):
            opt.zero_grad()
            for _w in range(Kw):
                cb = np.sort(perm[np.arange(ptr, ptr + B) % data.n_cells])
                ptr += B
                i0 = int(rng.choice(starts))
                th = model.theta(data, cb)
                o = model.simulate(data, cb, i0, i0 + win, th)
                loss, et, sm = combine(*sums(o, th, model, data, cb, i0, i0 + warm, i0 + win, SP.TRAIN))
                (loss / Kw).backward()
                acc += [float(loss), float(et), float(sm)]
                nacc += 1
            torch.nn.utils.clip_grad_norm_(model.gA.parameters(), CFG["grad_clip"])
            if all(torch.isfinite(p.grad).all() for p in model.gA.parameters() if p.grad is not None):
                opt.step()
        v = evaluate(model, data, val_cells, SP.VAL)
        tr = acc / max(nacc, 1)
        with open(logp, "a") as fh:
            fh.write(f"{ep},{tr[0]:.5f},{tr[1]:.5f},{tr[2]:.5f},{v[0]:.5f},{v[1]:.5f},{v[2]:.5f},"
                     f"{time.time() - t0:.0f}\n")
        if np.isfinite(v[0]) and v[0] < best:                 # selection on the trained quantity
            best = v[0]
            save_atomic({"gA": model.gA.state_dict(), "epoch": ep, "val": v[0]}, out / "best.pt")
        save_atomic({"gA": model.gA.state_dict(), "opt": opt.state_dict(), "epoch": ep, "best": best,
                     "np_rng": rng.bit_generator.state, "torch_rng": torch.get_rng_state()}, lastp)
        print(f"epoch {ep:02d} train {tr[0]:.4f} (et {tr[1]:.4f}, sm {tr[2]:.4f}) "
              f"val {v[0]:.4f} (et {v[1]:.4f}, sm {v[2]:.4f}) best {best:.4f} "
              f"[{time.time() - t0:.0f}s]", flush=True)
    (out / "run.json").write_text(json.dumps(dict(
        fold=fold, member=member, seed=seed, epochs_done=epoch0 + n_epochs, best_val=best,
        static_inputs=list(C.DLSTM_STATIC), cfg=CFG, hbv_bounds=C.HBV_BOUNDS), indent=1))
    return best


def infer(profile_name, fold, members=(0,), device=None, batch=1024):
    """Dekadal expected-ET anomaly for every cell from the best checkpoint(s),
    a split-conformal interval calibrated on validation dekads, and the residual
    on this fold's TEST dekads only."""
    profile = C.get_profile(profile_name)
    store = Store(profile.store)
    dev = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
    data = Data(store, fold)
    Kn, N = data.starts.size, data.n_cells
    runs, thetas = [], []
    for m in members:
        ck = torch.load(profile.store / "dlstm" / f"fold{fold}" / f"m{m}" / "best.pt",
                        map_location=dev, weights_only=False)
        model = Model(data.A.shape[1], dev)
        model.gA.load_state_dict(ck["gA"])
        model.gA.eval()
        et = np.full((Kn, N), np.nan, "float32")
        th_all = np.zeros((N, len(C.HBV_BOUNDS)), "float32")
        with torch.no_grad():
            for c0 in range(0, N, batch):
                cb = np.arange(c0, min(c0 + batch, N))
                th = model.theta(data, cb)
                o = model.simulate(data, cb, 0, data.n_days, th)["ET"].numpy()
                for k in range(Kn):
                    et[k, cb] = o[:, data.d0[k]:data.d1[k]].sum(axis=1)
                th_all[cb] = th.numpy()
        runs.append(et - data.cm[data.dek])
        thetas.append(th_all)
    runs = np.stack(runs)
    pred = runs.mean(axis=0)
    ens_sd = runs.std(axis=0) if len(members) > 1 else np.zeros_like(pred)
    y = np.asarray(store.arr("ET_anom"))
    from . import uq as U
    va, te = data.split == SP.VAL, data.split == SP.TEST
    margin = U.abs_margin((y - pred)[va])
    lo, hi = pred - margin, pred + margin
    out = Store(profile.store / "dlstm")
    names = ("oof_pred", "oof_lo", "oof_hi", "sigma_ensemble")
    cur = {n: (np.array(out.arr(n)) if out.has(n) else np.full((Kn, N), np.nan, "float32")) for n in names}
    tested = np.array(out.arr("oof_fold")) if out.has("oof_fold") else np.full(Kn, -1, "int8")
    for n, a in zip(names, (pred, lo, hi, ens_sd)):
        cur[n][te] = a[te]
    tested[te] = fold
    for n in names:
        out.save(n, cur[n])
    sig_model = U.combine_sigma(U.interval_to_sigma(cur["oof_lo"], cur["oof_hi"]), cur["sigma_ensemble"])
    sig_obs = U.observation_sigma(np.asarray(store.arr("ET_sd")), np.asarray(store.arr("ET_n")))
    sig_clim = np.asarray(store.arr("ET_clim_se"))[data.dek]
    out.save("oof_fold", tested)
    out.save("Bi", (y - cur["oof_pred"]).astype("float32"))
    for n, a in (("sigma_model", sig_model), ("sigma_obs", sig_obs), ("sigma_clim", sig_clim),
                 ("sigma_Bi", U.combine_sigma(sig_model, sig_obs, sig_clim))):
        out.save(n, np.where(np.isfinite(cur["oof_pred"]), a, np.nan).astype("float32"))
    out.save(f"theta_fold{fold}", np.stack(thetas).mean(axis=0))
    r = (y - pred)[te]
    ok = np.isfinite(r)
    res = dict(model="dlstm", fold=int(fold), members=list(members), n_test_dekads=int(te.sum()),
               test_rmse=float(np.sqrt(np.mean(r[ok] ** 2))) if ok.any() else float("nan"),
               test_rmse_climatology=float(np.sqrt(np.nanmean(y[te] ** 2))),
               conformal_margin=float(margin), coverage90=U.coverage(y[te], lo[te], hi[te]),
               folds_assembled=sorted(int(f) for f in np.unique(tested) if f >= 0))
    res["skill_vs_climatology"] = float(1 - (res["test_rmse"] / res["test_rmse_climatology"]) ** 2)
    out.write_json(f"metrics_fold{fold}", res)
    print(json.dumps(res), flush=True)
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Differentiable LSTM-HBV expected-ET model.")
    ap.add_argument("--profile", default="full", choices=sorted(C.PROFILES))
    ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--member", type=int, default=0, help="ensemble member (changes the seed)")
    ap.add_argument("--members", type=int, nargs="*", default=[0], help="members to average in --infer")
    ap.add_argument("--epochs", type=int, default=None)
    ap.add_argument("--steps", type=int, default=None, help="optimiser steps per epoch (default: all cells)")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--infer", action="store_true")
    ap.add_argument("--device", default=None, choices=["cuda", "cpu"])
    a = ap.parse_args()
    if a.infer:
        infer(a.profile, a.fold, tuple(a.members), a.device)
    else:
        train(a.profile, a.fold, a.member, a.epochs, a.resume, a.device, a.steps)
