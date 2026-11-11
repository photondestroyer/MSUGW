"""Stage `predict` (P1.2 output): rain-only expectation for EVERY valid block and month, scale F.

A block-month in water year wy is predicted by model (k, g) with k = the block's space fold and g = the
year group of wy - the same model whose test predictions the verification judged. Intervals carry the
conformal margins of that model. Calendar-year files: January-September come from water year y,
October-December from water year y + 1.
Output per tile-year: predict/<year>/<tile>/P1_expectation_<tile>_<year>.nc, dims (month, by, bx):
    <TARGET>_q05, _q50, _q95   int16 scaled (ET 0.1 mm, NIRv/NDMI 1e-4), fill -32768
    model_k, model_g           which model predicted (by, bx) / (month)
    reference                  P1.1 reference flag of the block-year
"""
import json
import os

import numpy as np

from .. import config as C
from ..core import ncio
from ..step0 import tiles as TL, drivers as DR
from . import blockset as BKS, folds as FO, features as FE, calendars as CA, baselines as BL

SCALE = dict(ET=0.1, NIRv=1e-4, NDMI=1e-4)


def load_models(train_dir, F, device):
    import torch
    from . import m1 as M1
    models = {}
    for k in range(F["folds"]):
        for g in range(len(F["year_groups"])):
            fdir = f"{train_dir}/k{k}_g{g}"
            meta = json.load(open(f"{fdir}/fold.json"))
            sc = FE.scaler_from_json(meta["scaler"])
            cal = CA.from_json(meta["calendars"])
            s = C.m1_settings()
            n_dyn = len(C.DYN_INPUTS)
            if meta.get("input_names") and meta["input_names"] != FE.input_names():
                raise RuntimeError(f"{fdir}: the model was trained with other inputs than features.py builds now")
            model = M1.make(n_dyn, meta.get("n_stat", FE.n_stat()), meta.get("n_pm", FE.n_pm()), s)
            model.load_state_dict(torch.load(f"{fdir}/m1/best.pt", map_location=device, weights_only=True))
            model.to(device)
            models[(k, g)] = dict(model=model, sc=sc, cal=cal, margins=BL.margins_from_json(meta["margins"]))
    return models


def predict_tile(T, G, ti, dirs, D, lut, F, models, device):
    from . import m1 as M1
    from .train import to_physical, calibrate
    BS = BKS.tile_blockset(T, G, ti, dirs, D, lut)
    BS["years_agg"], BS["years"] = np.array(C.agg_years()), np.array(C.years())
    ok = BS["valid"] & (BS["cell"] >= 0)
    k_of, _ = FO.fold_of(F, BS["lon"], BS["lat"])
    n = BS["block_id"].size
    years = C.years()
    wys = list(range(years[0], years[-1] + 2))
    dmax = FE.abs_month(D["months"][-1])
    res = {}                                                     # wy -> (n, 12, 3, 3) physical, calibrated
    for wy in wys:
        if 12 * wy + 8 > dmax:
            continue                                             # drivers end before this water year does
        g = int(FO.year_group_of(F, wy))
        P = np.full((n, 12, 3, 3), np.nan, "float32")
        for k in range(F["folds"]):
            sel = np.nonzero(ok & (k_of == k))[0]
            if sel.size == 0:
                continue
            M = models[(k, g)]
            X = FE.build(BS, D, sel, np.full(sel.size, wy), M["cal"])
            q = to_physical(M1.predict(M["model"], FE.apply_scaler(X, M["sc"]), device), M["sc"])
            P[sel] = calibrate(q, X, M["margins"])
        res[wy] = (P, g)
    return BS, res, k_of


def run(ctx):
    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    T, G = TL.load_table(ctx.dep("tiles"))
    D = DR.load(ctx.dep("drivers"))
    lut = DR.cell_lookup(D)
    F = FO.load(ctx.dep("folds"))
    dirs = dict(static=ctx.dep("static"), reference=ctx.dep("reference"), blocks=ctx.dep("blocks"))
    models = load_models(ctx.dep("train"), F, device)
    years = C.years()
    for ti, t in enumerate(T["tiles"]):
        ctx.check_budget()
        g_ = G[t]
        if all(ncio.is_complete(f"{ctx.out}/{y}/{t}/P1_expectation_{t}_{y}.nc", ctx.sig) for y in years):
            continue
        BS, res, k_of = predict_tile(T, G, ti, dirs, D, lut, F, models, device)
        for yi, y in enumerate(years):
            p = f"{ctx.out}/{y}/{t}/P1_expectation_{t}_{y}.nc"
            out = np.full((12, BS["block_id"].size, 3, 3), np.nan, "float32")
            mg = np.full(12, -1, "int8")
            if y in res:                                          # Jan-Sep from water year y (output months 3..11)
                out[0:9] = np.moveaxis(res[y][0][:, 3:12], 1, 0)
                mg[0:9] = res[y][1]
            if y + 1 in res:                                      # Oct-Dec from water year y + 1 (output months 0..2)
                out[9:12] = np.moveaxis(res[y + 1][0][:, 0:3], 1, 0)
                mg[9:12] = res[y + 1][1]
            variables = {}
            for j, nm in enumerate(C.TARGETS):
                for qi, qn in enumerate(("q05", "q50", "q95")):
                    a = out[:, :, j, qi].reshape(12, g_.nby, g_.nbx)
                    v = np.where(np.isfinite(a), np.clip(np.round(a / SCALE[nm]), -32767, 32767), -32768).astype("int16")
                    variables[f"{nm}_{qn}"] = dict(dims=("month", "by", "bx"), data=v, dtype="i2", fill=-32768,
                                                   attrs=dict(scale_factor=np.float32(SCALE[nm]),
                                                              units="mm month-1" if nm == "ET" else "1"))
            variables["model_k"] = dict(dims=("by", "bx"), data=k_of.reshape(g_.nby, g_.nbx).astype("int8"), dtype="i1")
            variables["model_g"] = dict(dims=("month",), data=mg, dtype="i1")
            variables["reference"] = dict(dims=("by", "bx"), data=BS["REF"][:, yi].reshape(g_.nby, g_.nbx), dtype="u1")
            prov = ctx.record({}, dict(tile=t, year=y))
            ncio.write_nc(p, dict(month=12, by=g_.nby, bx=g_.nbx), variables,
                          dict(tile=t, year=y, title="P1.2 rain-only expectation (calibrated 5/50/95 % quantiles)"),
                          prov, ctx.scratch)
        ctx.unit(t, "done")
    ctx.complete(dict(tiles=len(T["tiles"]), years=years))


def open_tile_year(stage_dir, tile, year):
    p = f"{stage_dir}/{year}/{tile}/P1_expectation_{tile}_{year}.nc"
    if not os.path.exists(p):
        return None
    with ncio.open_ro(p) as ds:
        out = {}
        for k in ds.variables:
            v = ds[k]
            out[k] = ncio.decode(np.asarray(v[...]), v) if k.split("_")[-1] in ("q05", "q50", "q95") else np.asarray(v[...])
        return out
