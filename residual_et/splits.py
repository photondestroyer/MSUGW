"""Time-blocked, climate-stratified 70:20:10 folds with an embargo.

Unit of the split: a half-year block (growing Apr-Sep, dormant Oct-Mar), so a
dry-down and its recovery stay together. All cells of a block share its role,
so no split ever holds a neighbouring pixel of the same date.

Each of the 5 folds tests one dry, one normal and one wet block (20 %),
validates on three half-blocks, again one per climate class (10 %), and trains
on the rest (70 %). Every block is tested exactly once, which gives an
out-of-fold residual for every dekad.
"""
import numpy as np

from . import config as C

TRAIN, VAL, TEST, EMBARGO = C.SPLIT_TRAIN, C.SPLIT_VAL, C.SPLIT_TEST, C.SPLIT_EMBARGO


def block_of(starts):
    """Chronological block number of each dekad, and per block: season
    (1 growing, 0 dormant) and a label such as 'G2018' or 'D2018' (Oct 2018-Mar 2019)."""
    s = np.asarray(starts, "datetime64[D]")
    year = s.astype("datetime64[Y]").astype(int) + 1970
    month = (s.astype("datetime64[M]") - s.astype("datetime64[Y]").astype("datetime64[M]")).astype(int) + 1
    grow = np.isin(month, C.BLOCK_GROW_MONTHS)
    byear = np.where(grow | (month >= 10), year, year - 1)
    key = byear * 2 + np.where(grow, 0, 1)            # growing season precedes the dormant one
    uniq, block = np.unique(key, return_inverse=True)
    season = (uniq % 2 == 0).astype("int8")
    labels = [("G" if g else "D") + str(int(k // 2)) for k, g in zip(uniq, season)]
    return block.astype("int32"), season, labels


def climate_classes(dryness, season):
    """0 dry, 1 normal, 2 wet. Blocks are ranked within their season (a dry
    winter and a dry summer are both 'dry'), then cut into three equal groups."""
    pct = np.empty(len(dryness), "float64")
    for g in (0, 1):
        m = season == g
        r = np.argsort(np.argsort(dryness[m], kind="stable"), kind="stable")
        pct[m] = (r + 0.5) / m.sum()
    order = np.argsort(pct, kind="stable")
    cls = np.empty(len(dryness), "int8")
    cls[order] = (np.arange(len(dryness)) * 3) // len(dryness)
    return cls


def assign_test_fold(cls, n_folds=C.N_FOLDS, seed=C.SEED):
    """Fold in which each block is tested. Blocks of one climate class are dealt
    to different folds, so every fold tests each class about equally."""
    rng = np.random.default_rng(seed)
    fold = np.empty(len(cls), "int8")
    ptr = 0
    for c in (0, 1, 2):
        idx = np.where(cls == c)[0]
        idx = idx[rng.permutation(idx.size)]
        for b in idx:
            fold[b] = ptr % n_folds
            ptr += 1
    return fold


def make_folds(starts, dryness_by_dekad, embargo, n_folds=C.N_FOLDS, seed=C.SEED):
    """Role of every dekad in every fold.

    starts            (K,) dekad start dates
    dryness_by_dekad  (K,) climate index (aquifer-mean SPI-90; lower = drier)
    embargo           dekads removed from training next to validation or test
    Returns a dict with `split` (n_folds, K) of TRAIN / VAL / TEST / EMBARGO."""
    block, season, labels = block_of(starts)
    nb = len(labels)
    dry = np.array([np.nanmean(dryness_by_dekad[block == b]) for b in range(nb)])
    cls = climate_classes(dry, season)
    tfold = assign_test_fold(cls, n_folds, seed)
    K = len(starts)
    split = np.full((n_folds, K), TRAIN, "int8")
    val_blocks = []
    for f in range(n_folds):
        split[f, np.isin(block, np.where(tfold == f)[0])] = TEST
        chosen = []
        for c in (0, 1, 2):
            cand = np.where((cls == c) & (tfold != f))[0]
            if cand.size == 0:
                continue
            nxt = cand[tfold[cand] == (f + 1) % n_folds]      # rotate, so validation moves with the fold
            b = int(nxt[0] if nxt.size else cand[0])
            dk = np.where(block == b)[0]
            half = dk[:dk.size // 2] if (f + c) % 2 == 0 else dk[dk.size // 2:]
            split[f, half] = VAL
            chosen.append(labels[b] + ("a" if (f + c) % 2 == 0 else "b"))
        val_blocks.append(chosen)
        # embargo: training dekads within `embargo` of a validation or test dekad
        held = split[f] != TRAIN
        near = np.zeros(K, bool)
        for lag in range(1, embargo + 1):
            near[lag:] |= held[:-lag]
            near[:-lag] |= held[lag:]
        split[f, (split[f] == TRAIN) & near] = EMBARGO
    return dict(split=split, block=block, block_label=labels, block_season=season,
                block_dryness=dry, block_class=cls, block_test_fold=tfold,
                val_halves=val_blocks, embargo=int(embargo))


def embargo_from_acf(anom, threshold=C.EMBARGO_ACF, bounds=C.EMBARGO_RANGE, max_lag=9):
    """Embargo length from the data: the first lag (in dekads) at which the
    median temporal autocorrelation of the ET anomaly falls below `threshold`."""
    a = np.asarray(anom, "float64")
    med = []
    for lag in range(1, max_lag + 1):
        x, y = a[:-lag], a[lag:]
        ok = np.isfinite(x) & np.isfinite(y)
        n = ok.sum(axis=0)
        xm = np.where(ok, x, 0).sum(0) / np.maximum(n, 1)
        ym = np.where(ok, y, 0).sum(0) / np.maximum(n, 1)
        cov = (np.where(ok, (x - xm) * (y - ym), 0)).sum(0)
        vx = (np.where(ok, (x - xm) ** 2, 0)).sum(0)
        vy = (np.where(ok, (y - ym) ** 2, 0)).sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            r = cov / np.sqrt(vx * vy)
        med.append(float(np.nanmedian(r[n >= 20])) if (n >= 20).any() else np.nan)
    lag = next((k + 1 for k, r in enumerate(med) if np.isfinite(r) and r < threshold), max_lag)
    return int(np.clip(lag, *bounds)), med


def summarize(folds):
    """Realised shares per fold, as assigned (embargo counted separately)."""
    out = []
    for f, s in enumerate(folds["split"]):
        k = s.size
        out.append(dict(fold=f, train=float((s == TRAIN).mean()), val=float((s == VAL).mean()),
                        test=float((s == TEST).mean()), embargo=float((s == EMBARGO).mean()),
                        train_assigned=float(((s == TRAIN) | (s == EMBARGO)).mean()), n_dekads=int(k)))
    return out


def check(folds):
    """Invariants that make the split leakage-free. Raises on violation."""
    split, block = folds["split"], folds["block"]
    tested = (split == TEST).sum(axis=0)
    if not (tested == 1).all():
        raise AssertionError("every dekad must be tested in exactly one fold")
    e = folds["embargo"]
    for f, s in enumerate(split):
        for b in np.unique(block[s == TEST]):
            if not (s[block == b] == TEST).all():
                raise AssertionError(f"fold {f}: a test block is split across roles")
        held = np.where(s != TRAIN)[0]
        held = held[s[held] != EMBARGO]
        train = np.where(s == TRAIN)[0]
        if train.size and held.size and e > 0:
            gap = np.abs(train[:, None] - held[None, :]).min()
            if gap <= e:
                raise AssertionError(f"fold {f}: a training dekad lies within the {e}-dekad embargo")
    return True
