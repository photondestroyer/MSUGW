"""Differentiable time-discrete HBV backbone (Feng et al. 2022, Table 1).

Copied from delta_gwb/hbv.py (2026-09-24). The water balance is unchanged; the
unused `warmup` argument was removed. Spin-up is handled by the caller, which
excludes the first days of a run from the loss.

States per cell: Sp (snowpack), Sliq (snow liquid), Ss (soil), Suz (upper
subsurface), Slz (lower subsurface). Fluxes: ET, Q0, Q1, Q2, routed Q.
The model has no groundwater-to-ET pathway and no irrigation input: whatever
ET those sustain is meant to appear in the residual.
"""
import torch
import torch.nn as nn


class HBV(nn.Module):
    def __init__(self, bounds, uh_len=60):
        super().__init__()
        self.names = list(bounds.keys())
        self.register_buffer("lo", torch.tensor([bounds[k][0] for k in self.names]))
        self.register_buffer("hi", torch.tensor([bounds[k][1] for k in self.names]))
        self.uh_len = uh_len

    def scale(self, raw):
        """Unconstrained (..., n_param) -> physical bounds."""
        return self.lo + (self.hi - self.lo) * torch.sigmoid(raw)

    def uh_kernel(self, a, tau):
        """Gamma unit hydrograph, shape (B, L), normalised to sum 1."""
        t = torch.arange(1, self.uh_len + 1, device=a.device, dtype=a.dtype)
        a, tau = a.unsqueeze(1), tau.unsqueeze(1)
        k = 1.0 / (torch.exp(torch.lgamma(a)) * tau.pow(a)) * t.pow(a - 1) * torch.exp(-t / tau)
        return k / k.sum(dim=1, keepdim=True).clamp_min(1e-12)

    def forward(self, P, T, Ep, theta):
        """P, T, Ep: (B, T) in mm/d, degC, mm/d. theta: (B, n_param), already scaled.
        Returns (B, T) tensors ET, Q, Q0, Q1, Q2, Ss, Suz, Slz, Sp and the
        cumulative mass-balance residual `mb_resid` (B,), ~0 if conservative."""
        B, Tp = P.shape
        dev, dt = P.device, P.dtype
        i = {k: j for j, k in enumerate(self.names)}
        g = lambda k: theta[:, i[k]].unsqueeze(1)
        FC, LP, BETA, GAMMA, PERC = g("FC"), g("LP"), g("BETA"), g("GAMMA"), g("PERC")
        K0, K1, K2, UZL = g("K0"), g("K1"), g("K2"), g("UZL")
        TT, DD, CWH, CFR = g("TT"), g("DD"), g("CWH"), g("CFR")
        zero = lambda: torch.zeros(B, 1, device=dev, dtype=dt)
        Sp, Sliq, Ss, Suz, Slz = zero(), zero(), zero(), zero(), zero()
        outs = {k: [] for k in ["ET", "Q0", "Q1", "Q2", "Ss", "Suz", "Slz", "Sp"]}
        cumP, cumOut = zero(), zero()
        for t in range(Tp):
            p = P[:, t:t + 1].clamp_min(0.0)
            tt = T[:, t:t + 1]
            ep = Ep[:, t:t + 1].clamp_min(0.0)
            # snow
            is_snow = (tt < TT).to(dt)
            psnow, prain = p * is_snow, p * (1.0 - is_snow)
            smelt = torch.minimum((DD * (tt - TT)).clamp_min(0.0), Sp)
            rfz = torch.minimum((DD * CFR * (TT - tt)).clamp_min(0.0), Sliq)
            Sp = Sp + psnow + rfz - smelt
            Sliq = Sliq + smelt - rfz
            isnow = (Sliq - CWH * Sp).clamp_min(0.0)
            Sliq = Sliq - isnow
            # soil. Bases are clamped at 1e-6: at Ss = 0 with an exponent below 1
            # the power-rule gradient is infinite and poisons the batch.
            Ss = Ss + prain + isnow
            W = torch.minimum((Ss / FC).clamp_min(1e-6) ** BETA, torch.ones_like(Ss))
            peff = W * (prain + isnow)
            eta = torch.minimum((Ss / (FC * LP)).clamp_min(1e-6) ** GAMMA, torch.ones_like(Ss))
            et = torch.minimum(eta * ep, Ss)
            Ss = Ss - peff - et
            ex = (Ss - FC).clamp_min(0.0)
            Ss = Ss - ex
            # upper zone
            Suz = Suz + peff + ex
            perc = torch.minimum(PERC, Suz)
            q0 = torch.minimum(K0 * (Suz - UZL).clamp_min(0.0), Suz)
            q1 = torch.minimum(K1 * Suz, Suz - q0)
            Suz = Suz - perc - q0 - q1
            # lower zone
            Slz = Slz + perc
            q2 = torch.minimum(K2 * Slz, Slz)
            Slz = Slz - q2
            cumP = cumP + p
            cumOut = cumOut + et + q0 + q1 + q2
            for k, v in (("ET", et), ("Q0", q0), ("Q1", q1), ("Q2", q2), ("Ss", Ss),
                         ("Suz", Suz), ("Slz", Slz), ("Sp", Sp)):
                outs[k].append(v)
        out = {k: torch.cat(v, dim=1) for k, v in outs.items()}
        q = out["Q0"] + out["Q1"] + out["Q2"]
        k = self.uh_kernel(theta[:, i["ROUT_A"]], theta[:, i["ROUT_TAU"]])
        L = self.uh_len
        uw = torch.nn.functional.pad(q, (L - 1, 0)).unfold(1, L, 1)       # (B, T, L)
        out["Q"] = (uw * k.flip(1).unsqueeze(1)).sum(dim=-1)              # causal routing
        s_end = (out["Sp"][:, -1] + out["Ss"][:, -1] + out["Suz"][:, -1] + out["Slz"][:, -1]
                 + Sliq.squeeze(1))
        out["mb_resid"] = cumP.squeeze(1) - cumOut.squeeze(1) - s_end
        return out
