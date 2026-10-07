"""A small synthetic store with a known signal, for tests that need the whole
schema (trainers, splits) without touching the real data."""
import numpy as np

from .. import climate as K
from .. import splits as SP
from ..store import Store, DYNAMIC, STATIC


def bucket(P, EP, fc=200.0, lp=0.6):
    """One-store water balance in numpy: the 'truth' the synthetic ET follows."""
    S = np.full(P.shape[1], 0.5 * fc)
    et = np.zeros_like(P)
    sm = np.zeros_like(P)
    for t in range(P.shape[0]):
        S = np.minimum(S + P[t], fc)
        e = np.minimum(EP[t] * np.minimum(S / (lp * fc), 1.0), S)
        S = S - e
        et[t], sm[t] = e, S / fc
    return et, sm


def make_store(path, n_cells=24, seed=0, t0="2015-04-01", t1="2022-06-01", forcing_t0="2013-04-01",
               noise_target=False):
    """Writes a complete store. With noise_target=True the ET anomaly is pure
    noise, unrelated to any feature (used to show that residuals are out-of-fold)."""
    rng = np.random.default_rng(seed)
    st = Store(path)
    starts = K.dekad_starts(t0, t1)
    kn, dek = starts.size, K.dek36(starts)
    days = np.arange(np.datetime64(forcing_t0, "D"), np.datetime64(t1, "D"))
    seas = np.sin(2 * np.pi * (K.doy(days) - 100) / 365.25)[:, None]
    wet = rng.normal(0, 1, (days.size // 30 + 2, 1)).repeat(30, axis=0)[:days.size]      # slow wet/dry spells
    P = (rng.gamma(0.25, 8.0, (days.size, n_cells)) * np.clip(1 + 0.5 * seas + 0.5 * wet, 0.05, None)).astype("float32")
    T = (12 + 12 * seas + rng.normal(0, 2, (days.size, n_cells))).astype("float32")
    EP = np.clip(3.5 + 3 * seas + rng.normal(0, 0.4, (days.size, n_cells)), 0.2, None).astype("float32")
    et_d, sm_d = bucket(P.astype("float64"), EP.astype("float64"))
    sm_d = (0.08 + 0.25 * sm_d).astype("float32")
    et = K.aggregate_dekads(days, et_d, starts, "sum")
    cl = K.DekadClimatology(n_cells)
    cl.add(dek, et)
    cm, cs, cse, _ = cl.finalize(pool=1)
    et_anom = et - cm[dek]
    sm = K.aggregate_dekads(days, sm_d, starts, "mean")
    feats = {}
    for name, how in (("P_dek", "sum"), ("VPD_anom", None), ("SRAD_anom", None), ("TMEAN_anom", None),
                      ("ETR_anom", None)):
        if how:
            feats[name] = K.aggregate_dekads(days, P, starts, how)
        else:
            feats[name] = rng.normal(0, 1, (kn, n_cells)).astype("float32")
    p30 = K.trailing_sum(days, P, starts, 30)
    p90 = K.trailing_sum(days, P, starts, 90)
    z = lambda a: ((a - np.nanmean(a, axis=0)) / np.nanstd(a, axis=0)).astype("float32")
    if noise_target:
        et_anom = rng.normal(0, 5, (kn, n_cells)).astype("float32")
        et = (cm[dek] + et_anom).astype("float32")
    dyn = dict(ET=et, ET_sd=np.full((kn, n_cells), 2.0, "float32"), ET_n=np.full((kn, n_cells), 12, "int16"),
               ET_anom=et_anom.astype("float32"), P_dek=feats["P_dek"], P_dek_anom=z(feats["P_dek"]),
               P30=p30, P90=p90, SPI30_g=z(p30), SPI90_g=z(p90),
               SPI30_se=np.full((kn, n_cells), 0.1, "float32"), SPI90_se=np.full((kn, n_cells), 0.1, "float32"),
               VPD_anom=feats["VPD_anom"], SRAD_anom=feats["SRAD_anom"], TMEAN_anom=feats["TMEAN_anom"],
               ETR_anom=feats["ETR_anom"], SMrz=sm, SMrz_lag1=np.vstack([np.full((1, n_cells), np.nan, "float32"), sm[:-1]]),
               spi30d=np.clip(z(p30), -2.09, 2.09), spi90d=np.clip(z(p90), -2.09, 2.09),
               spei30d=np.clip(z(p30), -2.09, 2.09), spei90d=np.clip(z(p90), -2.09, 2.09),
               eddi30d=np.clip(-z(p30), -2.09, 2.09))
    assert set(dyn) == set(DYNAMIC), set(DYNAMIC) ^ set(dyn)
    for k, v in dyn.items():
        st.save(k, v)
    for k, v in (("ET_clim", cm), ("ET_clim_sd", np.maximum(cs, 0.5)), ("ET_clim_se", cse)):
        st.save(k, v)
    for k in STATIC:
        st.save(k, rng.uniform(0, 1, n_cells).astype("float32"))
    st.save("WTD", rng.uniform(2, 80, n_cells).astype("float32"))
    for k, v in (("P_daily", P), ("T_daily", T), ("EP_daily", EP), ("SM_daily", sm_d)):
        st.save(k, v)
    dd = K.doy(days)
    for k, v in (("P_dclim", P), ("T_dclim", T), ("EP_dclim", EP)):
        c = np.stack([v[dd == d].mean(axis=0) if (dd == d).any() else v.mean(axis=0) for d in range(366)])
        st.save(k, c.astype("float32"))
    lat = 36 + 0.04 * (np.arange(n_cells) // 6)
    lon = -101 + 0.04 * (np.arange(n_cells) % 6)
    for k, v in (("dekad_start", starts), ("dek36", dek), ("day", days), ("cell_lat", lat), ("cell_lon", lon),
                 ("cell_iy", np.arange(n_cells) // 6), ("cell_ix", np.arange(n_cells) % 6),
                 ("space_block", (np.arange(n_cells) // 6).astype("int32"))):
        st.save(k, v)
    folds = SP.make_folds(starts, np.nanmean(dyn["spi90d"], axis=1), embargo=2)
    SP.check(folds)
    np.savez(st.path / "folds.npz", split=folds["split"], block=folds["block"],
             block_class=folds["block_class"], block_test_fold=folds["block_test_fold"],
             block_season=folds["block_season"], block_dryness=folds["block_dryness"],
             block_label=np.array(folds["block_label"]))
    st.require_complete()
    return st, folds
