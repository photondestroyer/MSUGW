"""Buffering analysis of out-of-fold residuals; identical for both models.

    python -m residual_et.buffering --profile full --model xgb

Bi = observed ET anomaly - expected ET anomaly, in mm/dekad, only where the
expectation came from a model that never saw that period. Positive Bi during a
dry-down is ET that weather and soil moisture do not explain: a CANDIDATE
subsurface or management signal (Research Plan 6.1), not groundwater by itself.

Changes against the two old analyses (findings L2, L3, S2, U1, U3, U4):
  * one dry-down rule and one set of strata for both models;
  * anomalies are per cell; drought conditions come from cited definitions and
    every statistic is reported under each of them;
  * confidence intervals are clustered in space and time;
  * the placebo uses out-of-fold residuals, like the main estimate;
  * the volume check compares like with like (mm per year on both sides).
"""
import argparse

import numpy as np

from . import config as C
from . import climate as K
from . import drought as D
from . import landcover as L
from . import uq as U
from .store import Store


def is_drydown(indices, sm_anom, sm_trend, vpd_anom):
    """Compound dry-down flag (Research Plan 6.1: precipitation deficit,
    root-zone soil-moisture decline, VPD and standardised indices).

    PROVISIONAL rule, kept from the old XGBoost script with the cited
    US Drought Monitor D1 bound in place of the uncited -0.8:
    SPI-90 at D1 or worse AND soil moisture below normal AND VPD above normal.
    indices: dict of (K, N) lab indices; the other arguments are (K, N) too."""
    with np.errstate(invalid="ignore"):
        return D.flag("usdm_d1", indices) & (sm_anom < 0) & (vpd_anom > 0)


def conditions(store):
    """Every drought condition as a boolean (K, N) field, plus the wet placebo."""
    idx = {n: np.asarray(store.arr(n)) for n in C.LAB_INDICES}
    dek = np.asarray(store.arr("dek36")).astype(int)
    sm = np.asarray(store.arr("SMrz"))
    # Event labels may use the whole record: they select rows, they are not model inputs.
    cl = K.DekadClimatology(sm.shape[1])
    cl.add(dek, sm)
    sm_mean, _, _, _ = cl.finalize(pool=1)
    sm_anom = sm - sm_mean[dek]
    sm_trend = np.vstack([np.full((1, sm.shape[1]), np.nan, "float32"), np.diff(sm, axis=0)])
    vpd_anom = np.asarray(store.arr("VPD_anom"))
    out = {name: D.flag(name, idx) for name in C.DROUGHT_DEFS}
    out["flash_sm"] = D.flash_drought(K.percentile_of(sm, dek, sm, dek, pool=1))
    out["compound_drydown"] = is_drydown(idx, sm_anom, sm_trend, vpd_anom)
    with np.errstate(invalid="ignore"):
        wet = idx["spi90d"] >= C.WET_THRESHOLD
    return out, wet, D.usdm_class(idx["spi90d"])


def strata(store):
    """Exposure and management strata. Never model inputs."""
    wtd = np.asarray(store.arr("WTD"))
    cuts = np.nanpercentile(wtd, [100 / 3, 200 / 3])
    wtd_class = np.where(np.isfinite(wtd), np.digitize(wtd, cuts), -1)       # 0 shallow .. 2 deep
    mgmt = L.classify_management(*(np.asarray(store.arr(n)) for n in
                                   ("f_irrigated", "f_rainfed_crop", "f_natural", "gir_high_freq")))
    return wtd_class, cuts, np.asarray(mgmt, "int8")


def within_pixel_contrast(bi, dry, base, sblock, n_boot=C.N_BOOT_EFFECT, seed=C.SEED, min_n=2):
    """Mean over cells of (mean Bi in dry dekads - mean Bi in the other dekads of
    the same cell). A constant bias of a cell cancels. Spatial-block bootstrap."""
    ok = np.isfinite(bi)
    a, b = dry & base & ok, ~dry & base & ok
    na, nb = a.sum(0), b.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.where(a, bi, 0).sum(0) / na - np.where(b, bi, 0).sum(0) / nb
    use = (na >= min_n) & (nb >= min_n) & np.isfinite(d)
    if use.sum() < 5:
        return dict(estimate=float("nan"), n_cells=int(use.sum()))
    d, sb = d[use], sblock[use]
    _, sb = np.unique(sb, return_inverse=True)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        w = np.bincount(rng.integers(0, sb.max() + 1, sb.max() + 1), minlength=sb.max() + 1)[sb]
        boots[i] = np.sum(w * d) / max(np.sum(w), 1)
    lo, hi = np.nanquantile(boots, [0.025, 0.975])
    return dict(estimate=float(d.mean()), ci_low=float(lo), ci_high=float(hi), n_cells=int(use.sum()))


def volume_check(bi, dry, grow, starts, store):
    """Residual ET per year against water-level change per year, in the same
    unit. 2017-2019 uses the 2017->2019 change map; the long-term rate uses
    predevelopment (~1950) -> 2019. Sy is the user constant."""
    year = starts.astype("datetime64[Y]").astype(int) + 1970
    sel = np.isin(year, (2017, 2018, 2019))[:, None] & dry & grow & np.isfinite(bi)
    extra = np.where(sel, bi, 0).sum(0) / 3.0                                     # mm per year
    n = sel.sum(0)
    mm_per_ft = 304.8
    recent = -np.asarray(store.arr("dWL_recent")) * C.SPECIFIC_YIELD * mm_per_ft / 2.0
    longterm = -np.asarray(store.arr("dWL_long")) * C.SPECIFIC_YIELD * mm_per_ft / 69.0
    ok = (n >= 3) & np.isfinite(recent)
    if ok.sum() < 10:
        return dict(n_cells=int(ok.sum()), note="too few cells with three dry dekads in 2017-2019")
    return dict(n_cells=int(ok.sum()),
                median_residual_ET_mm_per_year=float(np.median(extra[ok])),
                median_depletion_2017_2019_mm_per_year=float(np.median(recent[ok])),
                median_depletion_longterm_mm_per_year=float(np.nanmedian(longterm[ok])),
                correlation_with_recent_depletion=float(np.corrcoef(extra[ok], recent[ok])[0, 1]),
                note="both sides are mm of water per year; depletion = -dWL * Sy")


def run(profile_name, model):
    profile = C.get_profile(profile_name)
    store = Store(profile.store)
    res = Store(profile.store / model)
    bi = np.asarray(res.arr("Bi"))
    sigma = np.asarray(res.arr("sigma_Bi"))
    tested = np.asarray(res.arr("oof_fold"))
    starts = store.dekad_start
    month = (starts.astype("datetime64[M]") - starts.astype("datetime64[Y]").astype("datetime64[M]")).astype(int) + 1
    grow = np.isin(month, C.ANALYSIS_MONTHS)[:, None] & (tested >= 0)[:, None]
    fz = np.load(store.path / "folds.npz")
    tblock, sblock = fz["block"], np.asarray(store.arr("space_block"))
    conds, wet, usdm = conditions(store)
    wtd_class, wtd_cuts, mgmt = strata(store)
    shallow, deep = (wtd_class == 0)[None, :], (wtd_class == 2)[None, :]
    out = dict(model=model, profile=profile.name, n_dekads_with_residual=int((tested >= 0).sum()),
               folds_present=sorted(int(f) for f in np.unique(tested) if f >= 0),
               wtd_tercile_cuts_m=[float(c) for c in wtd_cuts],
               management_shares={L.MGMT_LABELS[k]: float(np.mean(mgmt == k)) for k in L.MGMT_LABELS},
               usdm_class_shares={n: float(np.mean(usdm == k)) for k, n in
                                  enumerate(("none", "D0", "D1", "D2", "D3", "D4"))},
               definitions={})
    boot = lambda a, b=None: U.block_bootstrap(bi, a, sblock, tblock, mask_b=b)
    for name, dry in conds.items():
        dg = dry & grow
        d = dict(share_of_growing_season=float(dry[grow[:, 0]].mean()) if grow.any() else float("nan"),
                 mean_Bi=boot(dg), old_style_ci=U.naive_ci(bi, dg),
                 within_pixel_dry_minus_other=within_pixel_contrast(bi, dry, grow, sblock),
                 share_exceeding_1p645_sigma=float(np.nanmean((bi > 1.645 * sigma)[dg & np.isfinite(bi)]))
                 if (dg & np.isfinite(bi)).any() else float("nan"),
                 by_wtd={lab: boot(dg & (wtd_class == k)[None, :]) for k, lab in
                         enumerate(("shallow", "middle", "deep"))},
                 by_management={L.MGMT_LABELS[k]: boot(dg & (mgmt == k)[None, :]) for k in L.MGMT_LABELS},
                 tau_shallow_minus_deep=boot(dg & shallow, dg & deep),
                 tau_within_management={L.MGMT_LABELS[k]: boot(dg & shallow & (mgmt == k)[None, :],
                                                               dg & deep & (mgmt == k)[None, :])
                                        for k in L.MGMT_LABELS})
        out["definitions"][name] = d
    wg = wet & grow
    out["placebo_wet"] = dict(rule=f"SPI-90 >= {C.WET_THRESHOLD}", mean_Bi=boot(wg),
                              tau_shallow_minus_deep=boot(wg & shallow, wg & deep))
    cd = conds["compound_drydown"]
    out["volume_check"] = volume_check(bi, cd, grow, starts, store)
    pix = np.where((cd & grow & np.isfinite(bi)).sum(0) >= 2,
                   np.nansum(np.where(cd & grow, bi, np.nan), axis=0)
                   / np.maximum((cd & grow & np.isfinite(bi)).sum(0), 1), np.nan)
    res.save("Bi_drydown_mean", pix.astype("float32"))
    res.write_json("buffering", out)
    m = out["definitions"]["compound_drydown"]["mean_Bi"]
    t = out["definitions"]["compound_drydown"]["tau_shallow_minus_deep"]
    print(f"{model}: compound dry-down mean Bi {m['estimate']:+.2f} [{m['ci_low']:+.2f}, {m['ci_high']:+.2f}] "
          f"mm/dekad (n={m['n']}); tau {t['estimate']:+.2f} [{t['ci_low']:+.2f}, {t['ci_high']:+.2f}]; "
          f"placebo tau {out['placebo_wet']['tau_shallow_minus_deep']['estimate']:+.2f}", flush=True)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Buffering statistics from out-of-fold residuals.")
    ap.add_argument("--profile", default="full", choices=sorted(C.PROFILES))
    ap.add_argument("--model", default="xgb", choices=["xgb", "dlstm"])
    a = ap.parse_args()
    run(a.profile, a.model)
