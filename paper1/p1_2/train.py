"""Stage `train` (P1.2): for every fold (k, g) - calendars, scalers, M1, conformal margins, M2, M0,
out-of-fold test predictions, planted-subsidy and class-flip tests; and on fold (0, 0) the variants of plan v3:
no pre-season inputs (RUN_ABLATION), no measurement inputs (RUN_MEAS_ABLATION) and ET / PET as the ET target
(RUN_RATIO_MODEL, no mass penalty).

Everything a fold learns comes from its training part only:
    calendars  <- NIRv of training reference block-months
    scalers    <- training items
    M1 weights <- training items (early stopping on validation items)
    margins    <- calibration items (CQR, per CAL_STRATA with the merge levels of baselines.py)
    M2         <- training items
    M0         <- the test block's own reference years outside the test year group +- embargo
Output per fold: <out>/k<k>_g<g>/ with m1/best.pt, last.pt, fold.json (calendars, scaler, margins, history,
plant tests, counts) and oof.npz (test items, observations, M1 / M2 / M0 predictions, in physical units).
Runs where torch exists (Colab, or locally the conda env SSM); it never imports scipy.
"""
import json
import os

import numpy as np

from .. import config as C
from .. import cdl_groups as CG
from ..core.drive import write_text
from ..step0 import drivers as DR
from . import blockset as BKS, folds as FO, features as FE, calendars as CA, baselines as BL

MAX_TRAIN_ITEMS = 1_500_000


def all_items(BS):
    """Every (block, water year) with a computed driver cell and >= 1 reference calendar year touched."""
    yrs = np.asarray(BS["years"]).astype(int)
    ii, ww = [], []
    for wy in yrs:
        r = np.searchsorted(yrs, wy)
        refy = BS["REF"][:, r] == 1
        if wy - 1 in yrs:
            refy |= BS["REF"][:, np.searchsorted(yrs, wy - 1)] == 1
        sel = np.nonzero(refy & (BS["cell"] >= 0))[0]
        ii.append(sel)
        ww.append(np.full(sel.size, wy))
    return np.concatenate(ii), np.concatenate(ww)


def fold_calendars(BS, train_blocks, train_years):
    yrs = np.asarray(BS["years"]).astype(int)
    ya = np.asarray(BS["years_agg"]).astype(int)
    nirv, month, grp = [], [], []
    for y in train_years:
        if y not in yrs or y not in ya:
            continue
        r, k = np.searchsorted(yrs, y), np.searchsorted(ya, y)
        b = train_blocks[BS["REF"][train_blocks, r] == 1]
        if b.size == 0:
            continue
        v = BS["Y"][b, k, :, 1]                                    # NIRv (n, 12)
        nirv.append(v.ravel())
        month.append(np.tile(np.arange(1, 13), b.size))
        grp.append(np.repeat(BS["GRP"][b, r], 12))
    if not nirv:
        return CA.fit({})
    nirv, month, grp = map(np.concatenate, (nirv, month, grp))
    return CA.fit(CA.curves_from_training(nirv, month, grp, np.ones(nirv.size, bool)))


def to_physical(q, sc):
    """Standardised (B,12,T,3) -> physical units."""
    return q * sc["y"][1][None, None, :, None] + sc["y"][0][None, None, :, None]


def levels_of(X, t):
    """Calibration strata keys (levels, B, 12) of target t for the items of X (baselines.strata_levels)."""
    return BL.strata_levels(t, X["func"][:, None], X["aridity"][:, None], X["in_season"], X["drought"], X["regime"])


def calibrate(pred, X, margins):
    out = pred.copy()
    for t in range(3):
        lo, hi = BL.apply_margins(pred[..., t, 0], pred[..., t, 2], levels_of(X, t), margins)
        out[..., t, 0], out[..., t, 2] = lo, hi
    return out


def fit_margins(pred, X):
    ys, los, his, keys = [], [], [], []
    for t in range(3):
        m = X["ymask"][..., t]
        ys.append(X["y"][..., t][m])
        los.append(pred[..., t, 0][m])
        his.append(pred[..., t, 2][m])
        keys.append(levels_of(X, t)[:, m])
    return BL.fit_margins(np.concatenate(ys), np.concatenate(los), np.concatenate(his), np.concatenate(keys, 1))


def class_flip(model, X_test, pred_test, items, margins, sc, device):
    """Plan v3 class-flip test: relabel CLASS_FLIP['share'] of the held-out blocks with the functional group they
    are most often confused with (CLASS_CONFUSION) and measure how often their in-season block-months then
    exceed the calibrated 95 % quantile (must stay <= PASS['fpr_max']); the same blocks unflipped alongside."""
    from . import m1 as M1
    import copy
    rng = np.random.default_rng(C.CLASS_FLIP["seed"])
    blocks = np.unique(items[0])
    chosen = set(rng.choice(blocks, max(1, int(round(C.CLASS_FLIP["share"] * blocks.size))), replace=False).tolist())
    rows = np.array([i in chosen for i in items[0]])
    Xf = {k: (v.copy() if isinstance(v, np.ndarray) else copy.copy(v)) for k, v in X_test.items()}
    flipped = np.zeros(rows.size, bool)
    for r in np.nonzero(rows)[0]:
        old = CG.FUNC_NAMES[int(Xf["func"][r])]
        new = C.CLASS_CONFUSION.get(old)
        if new:
            FE.set_class(Xf, np.array([r]), CG.FUNC_NAMES.index(new))
            flipped[r] = True
    if not flipped.any():
        return None
    qf = calibrate(to_physical(M1.predict(model, FE.apply_scaler(Xf, sc), device), sc), Xf, margins)
    qo = calibrate(pred_test, X_test, margins)
    out = {}
    for t, name in enumerate(C.TARGETS):
        m = X_test["ymask"][..., t] & X_test["in_season"] & flipped[:, None]
        if m.sum() == 0:
            out[name] = None
            continue
        y = X_test["y"][..., t][m]
        rate_f = float(np.mean(y > qf[..., t, 2][m]))
        rate_o = float(np.mean(y > qo[..., t, 2][m]))
        out[name] = dict(n=int(m.sum()), rate_flipped=rate_f, rate_original=rate_o,
                         passes=rate_f <= C.PASS["fpr_max"],
                         mean_abs_change_q50=float(np.mean(np.abs(qf[..., t, 1][m] - qo[..., t, 1][m]))))
    return out


def planted(model, BS, D, X_test, pred_test, items, cal, sc, device):
    """Verify (3): plant a subsidy of known size into 10 % of held-out reference blocks."""
    import torch  # noqa: F401
    from . import m1 as M1
    rng = np.random.default_rng(C.PLANT["seed"])
    blocks = np.unique(items[0])
    chosen = set(rng.choice(blocks, max(1, int(round(C.PLANT["share"] * blocks.size))), replace=False).tolist())
    is_p = np.array([i in chosen for i in items[0]])
    delta = np.array([C.PLANT["ET_mm"], C.PLANT["NIRv"], C.PLANT["NDMI"]], "float32")
    res = {}
    base_res = (X_test["y"] - pred_test[..., 1])                           # (B, 12, T) observed - median
    for variant in ("season_only", "with_preseason"):
        if variant == "season_only":
            q50p = pred_test[..., 1]
        else:
            BSp = dict(BS)
            BSp["Y"] = BS["Y"].copy()
            bl = np.array(sorted(chosen))
            BSp["Y"][bl] = BSp["Y"][bl] + delta                            # every month of the planted blocks
            Xp = FE.apply_scaler(FE.build(BSp, D, items[0], items[1], cal), sc)
            q50p = to_physical(M1.predict(model, Xp, device), sc)[..., 1]
        out = {}
        for t, name in enumerate(C.TARGETS):
            m_season = X_test["ymask"][..., t] & X_test["in_season"]
            yp = X_test["y"][..., t] + delta[t]
            planted_m = m_season & is_p[:, None]
            control_m = m_season & ~is_p[:, None]
            if planted_m.sum() == 0 or control_m.sum() == 0:
                out[name] = None
                continue
            est = float(np.mean((yp - q50p[..., t])[planted_m]) - np.mean(base_res[..., t][control_m]))
            out[name] = dict(planted=float(delta[t]), estimated=est, bias=est / float(delta[t]) - 1,
                             passes=abs(est / float(delta[t]) - 1) < C.PASS["plant_bias"], n=int(planted_m.sum()))
        res[variant] = out
    return res


VARIANTS = dict(main=dict(), nopre=dict(no_preseason=True), nomeas=dict(no_meas=True), ratio=dict(ratio=True))


def run_fold(ctx, BS, D, F, k, g, items, gkeys, device, variant="main"):
    from . import m1 as M1
    ablation = variant != "main"
    tag = f"k{k}_g{g}" + ("" if variant == "main" else f"_{variant}")
    fdir = f"{ctx.out}/{tag}"
    os.makedirs(f"{fdir}/m1", exist_ok=True)
    if os.path.exists(f"{fdir}/_done.json"):
        return json.load(open(f"{fdir}/_done.json"))
    fold = np.array([F["group_fold"].get(int(x), -1) for x in gkeys])
    pool = gkeys[(fold != k) & np.isin(items[1], FO.training_years(F, g))]
    inner = FO.inner_split(pool, F["seed"] + 1000 + k)
    role = FO.roles(F, k, g, gkeys, items[1], inner)
    sets = {r: (items[0][role == r], items[1][role == r]) for r in ("train", "val", "cal", "test")}
    if sets["train"][0].size > MAX_TRAIN_ITEMS:
        sel = np.random.default_rng(k * 10 + g).choice(sets["train"][0].size, MAX_TRAIN_ITEMS, replace=False)
        sets["train"] = (sets["train"][0][sel], sets["train"][1][sel])
    for r, (i, _) in sets.items():
        if i.size == 0:
            raise RuntimeError(f"fold {tag}: no {r} items - too few space groups or years for this design "
                               f"(see docs/CONFIG.md: SPACE_FOLDS, SPACE_GROUP_DEG, YEAR_GROUPS)")
    train_years = FO.training_years(F, g)
    cal = fold_calendars(BS, np.unique(sets["train"][0]), train_years)
    X = {r: FE.build(BS, D, i, w, cal, **VARIANTS[variant]) for r, (i, w) in sets.items()}
    sc = FE.fit_scaler(X["train"])
    Xs = {r: FE.apply_scaler(x, sc) for r, x in X.items()}
    s = C.m1_settings()
    if variant == "ratio":
        s["mass_weight"] = 0.0                     # the water-year bound is in mm; the ratio target has no mm
    model = M1.make(Xs["train"]["dyn"].shape[-1], Xs["train"]["stat"].shape[-1], Xs["train"]["pm"].shape[-1], s)
    model, state = M1.train(model, Xs["train"], Xs["val"], sc, s, f"{fdir}/m1", ctx.log, device)
    pred = {r: to_physical(M1.predict(model, Xs[r], device), sc) for r in ("cal", "test")}
    margins = fit_margins(pred["cal"], X["cal"])
    m1_test = calibrate(pred["test"], X["test"], margins)
    lv = np.stack([BL.level_of(levels_of(X["test"], t), margins) for t in range(3)], -1)
    out = dict(tag=tag, k=k, g=g, ablation=ablation, variant=variant,
               counts={r: int(i.size) for r, (i, _) in sets.items()}, inner=inner,
               calendars=CA.to_json(cal), scaler=FE.scaler_to_json(sc), margins=BL.margins_to_json(margins),
               history=state["history"], device=device, n_stat=FE.n_stat(), n_pm=FE.n_pm(),
               input_names=FE.input_names())
    if not ablation:
        m2 = BL.m2_fit(X["train"], X["train"]["ymask"].any(-1))
        m2_test = BL.m2_predict(m2, X["test"])
        m0_test = BL.m0(BS, sets["test"][0], X["test"]["t_out"], set(F["year_groups"][g]), BS["years"])
        out["plant"] = planted(model, BS, D, X["test"], pred["test"], sets["test"], cal, sc, device)
        out["class_flip"] = class_flip(model, X["test"], pred["test"], sets["test"], margins, sc, device)
    else:
        m2_test = m0_test = np.full(m1_test.shape, np.nan, "float32")
    xt = X["test"]
    np.savez_compressed(f"{fdir}/oof.npz", i=sets["test"][0], wy=sets["test"][1], t_out=xt["t_out"],
                        y=xt["y"], ymask=xt["ymask"], in_season=xt["in_season"],
                        aridity=xt["aridity"], group=xt["group"], m1=m1_test.astype("float32"),
                        m1_raw=pred["test"].astype("float32"), m2=m2_test, m0=m0_test,
                        block_id=BS["block_id"][sets["test"][0]], bound=xt["bound"], func=xt["func"],
                        func_prev=xt["func_prev"], spi=xt["spi"], drought=xt["drought"], regime=xt["regime"],
                        pet=xt["pet"], cal_level=lv, ratio=np.array(variant == "ratio"))
    write_text(f"{fdir}/fold.json", json.dumps(out, indent=1, default=float))
    write_text(f"{fdir}/_done.json", json.dumps(dict(tag=tag, counts=out["counts"]), indent=1))
    return out


def run(ctx):
    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    BS = BKS.load(ctx.dep("samples"))
    D = DR.load(ctx.dep("drivers"))
    F = FO.load(ctx.dep("folds"))
    items = all_items(BS)
    _, gk = FO.fold_of(F, BS["lon"][items[0]], BS["lat"][items[0]])
    ctx.log.info("train: %d block-water-year items from %d sampled blocks, device %s", items[0].size,
                 BS["block_id"].size, device)
    results = []
    for k in range(F["folds"]):
        for g in range(len(F["year_groups"])):
            ctx.check_budget()
            r = run_fold(ctx, BS, D, F, k, g, items, gk, device)
            results.append(dict(tag=r["tag"], counts=r["counts"]))
            ctx.unit(f"k{k}_g{g}", "done", **r["counts"])
    for flag, variant in ((C.RUN_ABLATION, "nopre"), (C.RUN_MEAS_ABLATION, "nomeas"), (C.RUN_RATIO_MODEL, "ratio")):
        if flag:
            ctx.check_budget()
            r = run_fold(ctx, BS, D, F, 0, 0, items, gk, device, variant=variant)
            results.append(dict(tag=r["tag"], counts=r["counts"]))
    ctx.complete(dict(folds=results, device=device))


def load_fold(train_dir, k, g):
    fdir = f"{train_dir}/k{k}_g{g}"
    meta = json.load(open(f"{fdir}/fold.json"))
    with np.load(f"{fdir}/oof.npz") as z:
        oof = {key: z[key] for key in z.files}
    return meta, oof
