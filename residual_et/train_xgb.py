"""Cross-fitted XGBoost model of the expected ET anomaly.

    python -m residual_et.train_xgb --profile full

For each of the 5 folds the model is fitted on that fold's training dekads,
stopped early and calibrated on its validation dekads, and used to predict only
its test dekads. Every dekad is tested once, so every residual

    Bi = ET anomaly (observed) - ET anomaly (expected)

comes from a model that never saw that period (finding L2). The target is the
anomaly, so skill is measured against climatology, not against the seasonal
cycle (finding A1). Rows are streamed dekad by dekad, so the full aquifer
never has to sit in memory as one table.
"""
import argparse
import json
import os
import time

import numpy as np
import xgboost as xgb

from . import config as C
from . import climate as K
from . import splits as SP
from . import uq as U
from .store import Store

SM_SD_FLOOR = 0.005          # m3/m3


class Design:
    """Feature rows, one dekad at a time, read from the store."""

    def __init__(self, store):
        self.dek = np.asarray(store.arr("dek36")).astype(int)
        self.y = store.arr("ET_anom")
        self.et_n = store.arr("ET_n")
        self.dyn = {n: store.arr(n) for n in C.XGB_DYNAMIC if not n.startswith("SMrz")}
        self.sm = np.asarray(store.arr("SMrz"))
        self.sm_lag = np.asarray(store.arr("SMrz_lag1"))
        self.static = np.column_stack([np.asarray(store.arr(n), "float32") for n in C.STATIC_FEATURES])
        self.names = list(C.XGB_DYNAMIC) + list(C.SEASON_FEATURES) + list(C.STATIC_FEATURES)
        C.assert_no_exposure(self.names)
        self.n_cells = self.static.shape[0]
        self.sm_mean = self.sm_sd = None

    def set_fold(self, split_row):
        """Soil-moisture climatology from this fold's TRAINING dekads only.
        (SMAP starts in 2015, so no earlier baseline exists.)"""
        tr = split_row == SP.TRAIN
        cl = K.DekadClimatology(self.n_cells)
        cl.add(self.dek[tr], self.sm[tr])
        self.sm_mean, sd, _, _ = cl.finalize(pool=1)
        self.sm_sd = np.maximum(sd, SM_SD_FLOOR)

    def block(self, k):
        """X (rows, features), y (rows,), cell index of each row, for dekad k.
        Rows need an observed target; missing features stay NaN (XGBoost routes
        them itself), so no dekad is dropped because one input is absent."""
        y = np.asarray(self.y[k])
        ok = np.isfinite(y) & (np.asarray(self.et_n[k]) > 0)
        q = self.dek[k]
        cols = []
        for name in C.XGB_DYNAMIC:
            if name == "SMrz_anom":
                v = (self.sm[k] - self.sm_mean[q]) / self.sm_sd[q]
            elif name == "SMrz_lag1_anom":
                v = (self.sm_lag[k] - self.sm_mean[(q - 1) % 36]) / self.sm_sd[(q - 1) % 36]
            else:
                v = np.asarray(self.dyn[name][k])
            cols.append(v)
        ang = 2 * np.pi * (q + 0.5) / 36.0
        cols += [np.full(self.n_cells, np.sin(ang)), np.full(self.n_cells, np.cos(ang))]
        X = np.column_stack(cols + [self.static]).astype("float32")
        return X[ok], y[ok].astype("float32"), np.where(ok)[0]


class DekadIter(xgb.DataIter):
    """Feeds XGBoost a few dekads at a time."""

    def __init__(self, design, dekads, batch=C.XGB_BATCH_DEKADS):
        self.design, self.dekads, self.batch, self.i = design, [int(k) for k in dekads], batch, 0
        super().__init__()

    def reset(self):
        self.i = 0

    def next(self, input_data):
        while self.i < len(self.dekads):
            ks = self.dekads[self.i:self.i + self.batch]
            self.i += self.batch
            parts = [self.design.block(k)[:2] for k in ks]
            X = np.vstack([p[0] for p in parts])
            if X.shape[0]:
                input_data(data=X, label=np.concatenate([p[1] for p in parts]))
                return True
        return False


def _fit(params, dtrain, dval, rounds):
    return xgb.train(params, dtrain, rounds, evals=[(dval, "val")],
                     early_stopping_rounds=C.XGB_EARLY_STOP, verbose_eval=False)


def _predict(bst, X):
    return bst.inplace_predict(X, iteration_range=(0, bst.best_iteration + 1))


def _scores(y, p):
    """Skill of an anomaly prediction. Climatology predicts zero anomaly."""
    ok = np.isfinite(y) & np.isfinite(p)
    y, p = y[ok], p[ok]
    if y.size < 3:
        return dict(n=int(y.size))
    mse, mse_clim = float(np.mean((y - p) ** 2)), float(np.mean(y ** 2))
    return dict(n=int(y.size), rmse=float(np.sqrt(mse)), mae=float(np.mean(np.abs(y - p))),
                rmse_climatology=float(np.sqrt(mse_clim)),
                skill_vs_climatology=float(1 - mse / mse_clim) if mse_clim > 0 else float("nan"),
                anomaly_correlation=float(np.corrcoef(y, p)[0, 1]), bias=float(np.mean(p - y)))


def run(profile_name, rounds=C.XGB_ROUNDS, folds=None, threads=None):
    t0 = time.time()
    profile = C.get_profile(profile_name)
    store = Store(profile.store)
    store.require_complete()
    out = Store(profile.store / "xgb")
    out.path.mkdir(parents=True, exist_ok=True)
    y_all = np.asarray(store.arr("ET_anom"))
    fz = np.load(store.path / "folds.npz")
    split = fz["split"]
    starts = store.dekad_start
    month = (starts.astype("datetime64[M]") - starts.astype("datetime64[Y]").astype("datetime64[M]")).astype(int) + 1
    grow = np.isin(month, C.ANALYSIS_MONTHS)
    design = Design(store)
    K_, N = split.shape[1], design.n_cells
    threads = threads or min(8, os.cpu_count() or 4)
    base = dict(C.XGB_PARAMS, seed=C.SEED, nthread=threads)
    pred = np.full((K_, N), np.nan, "float32")
    q50, lo, hi = pred.copy(), pred.copy(), pred.copy()
    tested_in = np.full(K_, -1, "int8")
    fold_metrics, gains = [], []
    for f in (range(split.shape[0]) if folds is None else folds):
        s = split[f]
        tr, va, te = (np.where(s == code)[0] for code in (SP.TRAIN, SP.VAL, SP.TEST))
        design.set_fold(s)
        dtrain = xgb.QuantileDMatrix(DekadIter(design, tr), max_bin=C.XGB_PARAMS["max_bin"])
        dval = xgb.QuantileDMatrix(DekadIter(design, va), ref=dtrain)
        m_mean = _fit(dict(base, objective="reg:squarederror"), dtrain, dval, rounds)
        m_q = _fit(dict(base, objective="reg:quantileerror",
                        quantile_alpha=np.array(C.XGB_QUANTILES)), dtrain, dval, rounds)
        # calibrate the 90 % interval on validation (conformalised quantile regression)
        yv, lv, hv = [], [], []
        for k in va:
            X, y, _ = design.block(k)
            if X.shape[0]:
                qq = _predict(m_q, X)
                yv.append(y), lv.append(qq[:, 0]), hv.append(qq[:, -1])
        margin = U.cqr_margin(np.concatenate(yv), np.concatenate(lv), np.concatenate(hv)) if yv else 0.0
        # out-of-fold prediction of the test dekads
        for k in te:
            X, _, cells = design.block(k)
            if X.shape[0]:
                qq = _predict(m_q, X)
                pred[k, cells] = _predict(m_mean, X)
                q50[k, cells], lo[k, cells], hi[k, cells] = qq[:, 1], qq[:, 0] - margin, qq[:, -1] + margin
            tested_in[k] = f
        m = dict(fold=int(f), n_train_dekads=int(tr.size), n_val_dekads=int(va.size),
                 n_test_dekads=int(te.size), n_embargo_dekads=int((s == SP.EMBARGO).sum()),
                 best_iteration=int(m_mean.best_iteration), conformal_margin=float(margin),
                 test=_scores(y_all[te], pred[te]),
                 test_growing_season=_scores(y_all[te][grow[te]], pred[te][grow[te]]),
                 coverage90=U.coverage(y_all[te], lo[te], hi[te]), mean_width=U.mean_width(lo[te], hi[te]))
        # in-sample fit on a few training dekads: shows why residuals must be out-of-fold
        Xs, ys = zip(*[design.block(k)[:2] for k in tr[:: max(1, tr.size // 12)]])
        m["train_in_sample"] = _scores(np.concatenate(ys), _predict(m_mean, np.vstack(Xs)))
        fold_metrics.append(m)
        g = m_mean.get_score(importance_type="gain")
        gains.append({design.names[int(k[1:])]: v for k, v in g.items()})
        m_mean.save_model(str(out.path / f"model_fold{f}.json"))
        m_q.save_model(str(out.path / f"qmodel_fold{f}.json"))
        print(f"fold {f}: test RMSE {m['test'].get('rmse', float('nan')):.2f} mm/dekad "
              f"(climatology {m['test'].get('rmse_climatology', float('nan')):.2f}), skill "
              f"{m['test'].get('skill_vs_climatology', float('nan')):+.3f}, coverage "
              f"{m['coverage90']:.2f} [{time.time() - t0:.0f}s]", flush=True)
        del dtrain, dval

    bi = (y_all - pred).astype("float32")
    sig_model = U.interval_to_sigma(lo, hi)
    sig_obs = U.observation_sigma(np.asarray(store.arr("ET_sd")), np.asarray(store.arr("ET_n")))
    sig_clim = np.asarray(store.arr("ET_clim_se"))[np.asarray(store.arr("dek36")).astype(int)]
    for name, arr in (("oof_pred", pred), ("oof_q50", q50), ("oof_lo", lo), ("oof_hi", hi),
                      ("oof_fold", tested_in), ("Bi", bi), ("sigma_model", sig_model),
                      ("sigma_obs", sig_obs), ("sigma_clim", sig_clim),
                      ("sigma_Bi", U.combine_sigma(sig_model, sig_obs, sig_clim))):
        out.save(name, arr)
    done = tested_in >= 0
    keys = sorted({k for g in gains for k in g})
    summary = dict(
        model="xgboost", profile=profile.name, features=design.names, n_cells=int(N),
        n_dekads=int(K_), dekads_with_out_of_fold_residual=int(done.sum()),
        pooled_out_of_fold=_scores(y_all[done], pred[done]),
        pooled_growing_season=_scores(y_all[done & grow], pred[done & grow]),
        pooled_coverage90=U.coverage(y_all[done], lo[done], hi[done]),
        fold_spread=dict(
            rmse=[m["test"].get("rmse") for m in fold_metrics],
            skill=[m["test"].get("skill_vs_climatology") for m in fold_metrics]),
        folds=fold_metrics,
        mean_gain={k: float(np.mean([g.get(k, 0.0) for g in gains])) for k in keys},
        params=dict(C.XGB_PARAMS, rounds=rounds, early_stopping=C.XGB_EARLY_STOP,
                    quantiles=list(C.XGB_QUANTILES), alpha=C.ALPHA),
        seconds=round(time.time() - t0, 1))
    out.write_json("metrics", summary)
    print(f"out-of-fold: {json.dumps(summary['pooled_out_of_fold'])}", flush=True)
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Cross-fitted XGBoost expected-ET-anomaly model.")
    ap.add_argument("--profile", default="full", choices=sorted(C.PROFILES))
    ap.add_argument("--rounds", type=int, default=C.XGB_ROUNDS)
    ap.add_argument("--folds", type=int, nargs="*", default=None, help="subset of folds (default all)")
    ap.add_argument("--threads", type=int, default=None)
    a = ap.parse_args()
    run(a.profile, a.rounds, a.folds, a.threads)
