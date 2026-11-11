"""M1: multi-output probabilistic sequence model of rain-only vegetation (torch).

dyn (B, 24, n_dyn) -> causal temporal convolution (TCN, dilations 1, 2, 4, ...) -> the 12 output steps
stat (B, n_stat)   -> MLP, broadcast to the 12 output months
pm  (B, 12, n_pm)  -> per-month features (calendar, class clock, pre-season state)
head -> per output month and target three quantiles (5, 50, 95 %), non-crossing by construction:
        q50 = a, q05 = a - softplus(b), q95 = a + softplus(c)
Loss = pinball loss over valid (masked) targets, in standardised units
     + mass_weight * mean(relu(sum_ET_q50 - bound * (1 + V2_TOLERANCE)) / 100 mm)^2
       (water-year ET of the median may not exceed the water the block could have used: P - storage change)
Only this module imports torch.
"""
import math

import torch
import torch.nn as nn
import torch.nn.functional as F

from .. import config as C


class CausalConv(nn.Module):
    def __init__(self, cin, cout, k, d, dropout):
        super().__init__()
        self.pad = (k - 1) * d
        self.conv = nn.Conv1d(cin, cout, k, dilation=d)
        self.res = nn.Conv1d(cin, cout, 1) if cin != cout else nn.Identity()
        self.drop = nn.Dropout(dropout)

    def forward(self, x):                       # x (B, C, T)
        y = self.conv(F.pad(x, (self.pad, 0)))
        return F.gelu(self.drop(y)) + self.res(x)


class M1(nn.Module):
    def __init__(self, n_dyn, n_stat, n_pm, hidden=64, tcn_levels=3, kernel=3, static_hidden=64, dropout=0.1,
                 n_targets=3, **_):
        super().__init__()
        layers, c = [], n_dyn
        for lvl in range(tcn_levels):
            layers.append(CausalConv(c, hidden, kernel, 2 ** lvl, dropout))
            c = hidden
        self.tcn = nn.Sequential(*layers)
        self.static = nn.Sequential(nn.Linear(n_stat, static_hidden), nn.GELU(), nn.Linear(static_hidden, static_hidden))
        self.head = nn.Sequential(nn.Linear(hidden + static_hidden + n_pm, hidden), nn.GELU(), nn.Dropout(dropout),
                                  nn.Linear(hidden, n_targets * 3))
        self.n_targets = n_targets

    def forward(self, dyn, stat, pm):
        h = self.tcn(dyn.transpose(1, 2)).transpose(1, 2)[:, -12:, :]
        s = self.static(stat)[:, None, :].expand(-1, 12, -1)
        z = self.head(torch.cat([h, s, pm], -1)).view(dyn.shape[0], 12, self.n_targets, 3)
        mid = z[..., 0]
        lo = mid - F.softplus(z[..., 1])
        hi = mid + F.softplus(z[..., 2])
        return torch.stack([lo, mid, hi], -1)    # (B, 12, targets, 3) in standardised units


def pinball(q, y, mask, quantiles=C.QUANTILES):
    """Mean pinball loss over masked entries. q (B,12,T,3), y (B,12,T), mask (B,12,T) bool."""
    tot, n = 0.0, mask.sum().clamp(min=1)
    for k, tau in enumerate(quantiles):
        e = y - q[..., k]
        tot = tot + (torch.maximum(tau * e, (tau - 1) * e) * mask).sum()
    return tot / (n * len(quantiles))


def mass_penalty(q, bound, et_mean, et_sd):
    """q50 of ET (target 0) summed over the 12 months vs the water available in the water year (mm)."""
    et = q[:, :, 0, 1] * et_sd + et_mean
    excess = F.relu(et.sum(1) - bound * (1 + C.V2_TOLERANCE)) / 100.0
    return (excess ** 2).mean()


def make(n_dyn, n_stat, n_pm, settings):
    torch.manual_seed(settings["seed"])
    return M1(n_dyn, n_stat, n_pm, **settings)


def to_tensors(X, idx, device):
    t = lambda a, dt=torch.float32: torch.as_tensor(a[idx], dtype=dt, device=device)
    return (t(X["dyn"]), t(X["stat"]), t(X["pm"]), t(X["y_std"]), t(X["ymask"], torch.bool), t(X["bound"]))


def predict(model, X, device, batch=2048):
    model.eval()
    out = []
    with torch.no_grad():
        for s in range(0, X["dyn"].shape[0], batch):
            idx = slice(s, s + batch)
            dyn, stat, pm = (torch.as_tensor(X[k][idx], dtype=torch.float32, device=device) for k in ("dyn", "stat", "pm"))
            out.append(model(dyn, stat, pm).float().cpu())
    return torch.cat(out).numpy() if out else None


def train(model, Xtr, Xva, sc, settings, ckpt_dir, log=None, device=None):
    """Train with early stopping on validation pinball loss; resumes from ckpt_dir/last.pt."""
    import os
    import numpy as np
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=settings["lr"], weight_decay=settings["weight_decay"])
    use_amp = settings.get("amp", True) and device == "cuda"
    scaler = torch.amp.GradScaler("cuda", enabled=use_amp)
    et_mean, et_sd = float(sc["y"][0][0]), float(sc["y"][1][0])
    state = dict(epoch=0, best=math.inf, bad=0, history=[])
    last, best = f"{ckpt_dir}/last.pt", f"{ckpt_dir}/best.pt"
    if os.path.exists(last):
        ck = torch.load(last, map_location=device, weights_only=False)
        model.load_state_dict(ck["model"])
        opt.load_state_dict(ck["opt"])
        state = ck["state"]
        torch.set_rng_state(ck["rng_cpu"])
        np.random.set_state(ck["rng_np"])
        if log:
            log.info("resumed at epoch %d (best %.4f)", state["epoch"], state["best"])
    n = Xtr["dyn"].shape[0]
    while state["epoch"] < settings["max_epochs"] and state["bad"] < settings["patience"]:
        model.train()
        perm = np.random.default_rng(settings["seed"] + state["epoch"]).permutation(n)   # same order after a resume
        tl = 0.0
        for s in range(0, n, settings["batch"]):
            idx = perm[s:s + settings["batch"]]
            dyn, stat, pm, y, m, bound = to_tensors(Xtr, idx, device)
            with torch.autocast(device_type="cuda", enabled=use_amp):
                q = model(dyn, stat, pm)
            q = q.float()
            loss = pinball(q, y, m) + settings["mass_weight"] * mass_penalty(q, bound, et_mean, et_sd)
            opt.zero_grad(set_to_none=True)
            scaler.scale(loss).backward()
            scaler.step(opt)
            scaler.update()
            tl += float(loss) * len(idx)
        vl = evaluate(model, Xva, device)
        state["epoch"] += 1
        state["history"].append(dict(epoch=state["epoch"], train=tl / max(n, 1), val=vl))
        if vl < state["best"] - 1e-5:
            state["best"], state["bad"] = vl, 0
            torch.save(model.state_dict(), best + ".tmp")
            os.replace(best + ".tmp", best)
        else:
            state["bad"] += 1
        torch.save(dict(model=model.state_dict(), opt=opt.state_dict(), state=state, rng_cpu=torch.get_rng_state(),
                        rng_np=np.random.get_state()), last + ".tmp")
        os.replace(last + ".tmp", last)
        if log:
            log.info("epoch %d  train %.4f  val %.4f  (best %.4f, %d bad)", state["epoch"], tl / max(n, 1), vl,
                     state["best"], state["bad"])
    model.load_state_dict(torch.load(best, map_location=device, weights_only=True))
    return model, state


def evaluate(model, X, device, batch=4096):
    model.eval()
    tot, cnt = 0.0, 0
    with torch.no_grad():
        for s in range(0, X["dyn"].shape[0], batch):
            idx = slice(s, s + batch)
            dyn, stat, pm, y, m, _ = (torch.as_tensor(X[k][idx], dtype=dt, device=device) for k, dt in
                                      (("dyn", torch.float32), ("stat", torch.float32), ("pm", torch.float32),
                                       ("y_std", torch.float32), ("ymask", torch.bool), ("bound", torch.float32)))
            q = model(dyn, stat, pm)
            k = int(m.sum())
            tot += float(pinball(q, y, m)) * max(k, 1)
            cnt += max(k, 1)
    return tot / max(cnt, 1)
