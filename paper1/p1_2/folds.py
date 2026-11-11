"""Spatial x year folds for P1.2 (plan rule R3), defined for EVERY block, not only the sample.

Space: blocks are grouped by SPACE_GROUP_DEG cells of their centre (0.5 deg ~ 50 km). Groups are sorted by
mean aridity, and each run of SPACE_FOLDS consecutive groups is dealt to the folds in a random order, so
every fold spans the aridity gradient.
Years: water years (Oct-Sep, labelled by the year in which they end) are split into contiguous
YEAR_GROUPS. Model (k, g) is trained on blocks outside space fold k and water years outside group g and
outside EMBARGO_YEARS around it, and predicts blocks of fold k in water years of g. A block-month is
therefore always predicted by the same model, whether the block is a reference block or not.
Inside the training part, space groups are split INNER_SPLIT into train / validation (early stopping) /
calibration (conformal margins); with CAL_YEAR_HOLDOUT the calibration years are also withheld.
"""
import json

import numpy as np

from .. import config as C


def group_key(lon, lat, deg=None):
    deg = deg or C.space_group_deg()
    return (np.floor(np.asarray(lon) / deg).astype("int64") * 100000 + np.floor(np.asarray(lat) / deg).astype("int64"))


def make(lon, lat, aridity, seed=42):
    K = C.space_folds()
    gk = group_key(lon, lat)
    groups = np.unique(gk)
    mean_ar = np.array([np.nanmean(np.where(aridity[gk == g] > 0, aridity[gk == g], np.nan)) for g in groups])
    order = groups[np.argsort(np.nan_to_num(mean_ar, nan=99), kind="stable")]
    rng = np.random.default_rng(seed)
    fold = {}
    for s in range(0, order.size, K):
        run_ = order[s:s + K]
        for g, k in zip(run_, rng.permutation(K)[:run_.size]):
            fold[int(g)] = int(k)
    inner = {k: inner_split([g for g, f in fold.items() if f != k], seed + 1000 + k) for k in range(K)}
    return dict(deg=C.space_group_deg(), folds=K, group_fold=fold, inner=inner, year_groups=C.year_groups(),
                embargo_years=C.EMBARGO_YEARS, seed=seed)


def fold_of(F, lon, lat):
    gk = group_key(lon, lat, F["deg"])
    return np.array([F["group_fold"].get(int(g), -1) for g in gk], "int16"), gk


def year_group_of(F, wy):
    """Index of the year group holding water year wy (outside every group: the nearest one)."""
    wy = np.asarray(wy)
    out = np.empty(wy.shape, "int16")
    mids = np.array([np.mean(g) for g in F["year_groups"]])
    for i, y in np.ndenumerate(wy):
        hit = [k for k, g in enumerate(F["year_groups"]) if int(y) in g]
        out[i] = hit[0] if hit else int(np.argmin(np.abs(mids - y)))
    return out


def training_years(F, g):
    test = set(F["year_groups"][g])
    lo, hi = min(test) - F["embargo_years"], max(test) + F["embargo_years"]
    allyears = sorted({y for grp in F["year_groups"] for y in grp})
    return [y for y in allyears if not (lo <= y <= hi)] if len(F["year_groups"]) > 1 else allyears


def inner_split(groups, seed):
    """Split space groups into train / validation / calibration (INNER_SPLIT), at least one group each when
    there are three or more."""
    pool = np.array(sorted(set(int(x) for x in groups)))
    perm = np.random.default_rng(seed).permutation(pool)
    n_tr = int(round(C.INNER_SPLIT[0] * pool.size))
    n_va = int(round(C.INNER_SPLIT[1] * pool.size))
    if pool.size >= 3:
        n_va = max(n_va, 1)
        n_tr = min(max(n_tr, 1), pool.size - n_va - 1)
    return dict(train=perm[:n_tr].tolist(), val=perm[n_tr:n_tr + n_va].tolist(), cal=perm[n_tr + n_va:].tolist())


def roles(F, k, g, group_keys, wy, inner=None):
    """Role of each (block, water year) item for model (k, g): 'train' 'val' 'cal' 'test' or ''.
    inner: the train/val/cal split of the training space groups; default = the one drawn in `make` over
    all groups. The train stage passes a split drawn over the groups that actually hold items."""
    gk = np.asarray(group_keys)
    fold = np.array([F["group_fold"].get(int(x), -1) for x in gk])
    inner = inner or F["inner"][k]
    tr_years = set(training_years(F, g))
    test_years = set(F["year_groups"][g])
    out = np.full(gk.shape, "", dtype=object)
    in_tr_years = np.isin(wy, sorted(tr_years))
    out[(fold == k) & np.isin(wy, sorted(test_years))] = "test"
    for role in ("train", "val", "cal"):
        m = np.isin(gk, inner[role]) & (fold != k) & in_tr_years
        out[m & (out == "")] = role
    return out


def run(ctx):
    """Stage `folds`: computed over ALL valid blocks of the run (so predict can use the same folds)."""
    from ..step0 import tiles as TL, drivers as DR
    T, _ = TL.load_table(ctx.dep("tiles"))
    D = DR.load(ctx.dep("drivers"))
    lut = DR.cell_lookup(D)
    ok = T["valid"] & (T["lab_iy"] >= 0)
    cell = lut[np.maximum(T["lab_iy"], 0), np.maximum(T["lab_ix"], 0)]
    ar = np.where(ok & (cell >= 0), D["aridity_decile"][np.maximum(cell, 0)], 0)
    F = make(T["lon"][ok], T["lat"][ok], ar[ok])
    with open(f"{ctx.out}/folds.json", "w") as f:
        json.dump(F, f, indent=1)
    sizes = np.bincount(np.array(list(F["group_fold"].values())), minlength=F["folds"]).tolist()
    ctx.complete(dict(space_groups=len(F["group_fold"]), groups_per_fold=sizes, year_groups=F["year_groups"]))


def load(stage_dir):
    with open(f"{stage_dir}/folds.json") as f:
        F = json.load(f)
    F["group_fold"] = {int(k): v for k, v in F["group_fold"].items()}
    F["inner"] = {int(k): v for k, v in F["inner"].items()}
    return F
