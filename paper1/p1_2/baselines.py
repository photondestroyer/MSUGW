"""M0 (climatology floor) and M2 (ecological benchmark) of the plan, plus conformal margins.

M0  for a block and calendar month: median and 5/95 % quantiles of the block's OWN reference values of
    that month in water years outside the test year group and its embargo (>= M0_MIN_YEARS years).
    Under spatial hold-out a test block never appears in training, so its own other years are the only
    climatology it has; that is the comparison the plan asks for ("anomaly RMSE").
M2  rain-use-efficiency / residual-trend benchmark (Evans & Geerken 2004; Burrell et al. 2017): per class
    group x aridity decile x calendar month, target = a + b * rain over the last w months (incl. the month),
    w chosen from M2_WINDOWS_MONTHS by the highest R^2 on the training items of the stratum; intervals
    from the training residual quantiles. Strata with < M2_MIN_SAMPLES fall back to group x month, then month.
Conformal: conformalised quantile regression margins (Romano et al. 2019), vendored from residual_et/uq.py.
Calibration strata (plan v3): CAL_STRATA = target x functional group x aridity tercile x season x drought class
(SPI-90) x observation regime. A margin is fitted at every merge level: level 0 uses all dimensions, level k drops
the first k names of CAL_MERGE_ORDER (regime, then tercile, then drought, functional, season), the last level is
the target alone. A stratum gets its own margin when it holds >= cal_min_n() calibration block-months (the
target-only level always does); a block-month takes the margin of the finest level that has one.
"""
import numpy as np

from .. import config as C


# ------------------------------------------------------------------ conformal (vendored from residual_et.uq)
def conformal_quantile(scores, alpha=C.ALPHA):
    s = np.sort(np.asarray(scores, "float64")[np.isfinite(scores)])
    n = s.size
    if n == 0:
        return float("nan")
    k = int(np.ceil((n + 1) * (1 - alpha))) - 1
    return float(s[min(max(k, 0), n - 1)])


def cqr_margin(y, lo, hi, alpha=C.ALPHA):
    return conformal_quantile(np.maximum(lo - y, y - hi), alpha)


def coverage(y, lo, hi):
    ok = np.isfinite(y) & np.isfinite(lo) & np.isfinite(hi)
    return float(np.mean((y[ok] >= lo[ok]) & (y[ok] <= hi[ok]))) if ok.any() else float("nan")


DIMS = ["target", "functional", "aridity_tercile", "season", "drought", "regime"]


def tercile(aridity):
    a = np.asarray(aridity)
    return np.where(a <= 3, 0, np.where(a <= 7, 1, 2))


def strata_levels(target, func, aridity, season, drought, regime):
    """Keys of every merge level, shape (levels, *broadcast shape). Key = sum over kept dimensions of
    (value + 1) * 100 ** position, so keys of different levels never collide."""
    vals = dict(target=target, functional=func, aridity_tercile=tercile(aridity), season=np.asarray(season).astype(int),
                drought=drought, regime=regime)
    shp = np.broadcast_shapes(*[np.shape(v) for v in vals.values()])
    vals = {k: np.broadcast_to(np.asarray(v).astype("int64"), shp) for k, v in vals.items()}
    used = [d for d in DIMS if d in C.CAL_STRATA or d == "target"]
    dropped, levels = set(), []
    for d in [None] + [x for x in C.CAL_MERGE_ORDER if x in used and x != "target"]:
        if d:
            dropped.add(d)
        key = np.zeros(shp, "int64")
        for pos, name in enumerate(DIMS):
            if name in used and name not in dropped:
                key += (vals[name] + 1) * 100 ** pos
        levels.append(key)
    return np.stack(levels)


def fit_margins(y, lo, hi, levels, min_n=None):
    """{key: margin} for every stratum of every level with >= min_n block-months (last level: any)."""
    min_n = C.cal_min_n() if min_n is None else min_n
    out = {}
    L = levels.shape[0]
    for li in range(L):
        keys = levels[li]
        u, inv, cnt = np.unique(keys, return_inverse=True, return_counts=True)
        for j, k in enumerate(u):
            if cnt[j] >= min_n or (li == L - 1 and cnt[j] > 0):
                m = inv == j
                out[int(k)] = cqr_margin(y[m], lo[m], hi[m])
    return out


def margin_of(levels, margins):
    """Margin per element from the finest level that has one (0 where none - no calibration data at all)."""
    m = np.full(levels.shape[1:], np.nan)
    for li in range(levels.shape[0]):
        todo = np.isnan(m)
        if not todo.any():
            break
        u, inv = np.unique(levels[li][todo], return_inverse=True)
        val = np.array([margins.get(int(k), np.nan) for k in u])
        m[todo] = val[inv]
    return np.nan_to_num(m)


def level_of(levels, margins):
    """Index of the level whose margin each element uses (levels count when none)."""
    lv = np.full(levels.shape[1:], levels.shape[0], "int16")
    for li in range(levels.shape[0] - 1, -1, -1):
        u, inv = np.unique(levels[li], return_inverse=True)
        has = np.array([int(k) in margins for k in u])[inv].reshape(lv.shape)
        lv[has] = li
    return lv


def margins_to_json(m):
    return {str(k): v for k, v in m.items()}


def margins_from_json(d):
    return {int(k): v for k, v in d.items() if k != "_all"}


def apply_margins(lo, hi, levels, margins):
    m = margin_of(levels, margins)
    return lo - m, hi + m


# ------------------------------------------------------------------ M0
def m0(BS, items_i, t_out, test_wy_set, yrs):
    """Median/quantiles per (item, output month) from the block's own reference values in other water years."""
    yrs = np.asarray(yrs).astype(int)
    ya = np.asarray(BS["years_agg"]).astype(int)
    out = np.full(t_out.shape + (3, 3), np.nan, "float32")                 # (B, 12, target, q)
    lo_e, hi_e = min(test_wy_set) - C.EMBARGO_YEARS, max(test_wy_set) + C.EMBARGO_YEARS
    cache = {}
    for b in range(items_i.size):
        i = int(items_i[b])
        if i not in cache:
            Y = BS["Y"][i]                                                  # (nya, 12, 3)
            vals = {}
            for m in range(12):
                rows = []
                for k, y in enumerate(ya):
                    wy = y + 1 if m >= 9 else y                              # water year of (y, month m)
                    if lo_e <= wy <= hi_e:
                        continue
                    r = np.searchsorted(yrs, y)
                    if r < yrs.size and yrs[r] == y and BS["REF"][i, r] == 1:
                        rows.append(Y[k, m])
                if len(rows) >= C.M0_MIN_YEARS:
                    a = np.array(rows)
                    with np.errstate(invalid="ignore"):
                        vals[m] = np.stack([np.nanquantile(a, q, axis=0) for q in C.QUANTILES], -1)
            cache[i] = vals
        for j in range(12):
            m = int(t_out[b, j] % 12)
            if m in cache[i]:
                out[b, j] = cache[i][m]
    return out


# ------------------------------------------------------------------ M2
def m2_fit(X, ok):
    """X from features.build (unscaled); ok (B, 12) usable training item-months. Returns a model dict."""
    P = X["dyn"][..., C.DYN_INPUTS.index("P")]                              # (B, 24)
    acc = {w: np.stack([P[:, 12 + j - w + 1:12 + j + 1].sum(1) for j in range(12)], 1) for w in C.M2_WINDOWS_MONTHS}
    month = (X["t_out"] % 12)
    grp = np.broadcast_to(X["group"][:, None], month.shape)
    dec = np.broadcast_to(X["aridity"][:, None], month.shape)
    model = {}
    for t in range(3):
        y = X["y"][..., t]
        okt = ok & np.isfinite(y)
        for level, keyf in (("gdm", lambda: grp * 10000 + dec * 100 + month), ("gm", lambda: grp * 100 + month),
                            ("m", lambda: month)):
            keys = keyf()
            for k in np.unique(keys[okt]):
                m = okt & (keys == k)
                if m.sum() < C.M2_MIN_SAMPLES:
                    continue
                best = None
                for w, A in acc.items():
                    x = A[m]
                    if np.nanstd(x) < 1e-6:
                        continue
                    b, a = np.polyfit(x, y[m], 1)
                    res = y[m] - (a + b * x)
                    r2 = 1 - np.var(res) / max(np.var(y[m]), 1e-12)
                    if best is None or r2 > best[0]:
                        best = (r2, w, a, b, np.quantile(res, C.QUANTILES))
                if best:
                    model[(t, level, int(k))] = best
    return model


def _acc(P):
    """Rain over the last w months incl. the output month: {w: (B, 12)}."""
    return {w: np.stack([P[:, 12 + j - w + 1:12 + j + 1].sum(1) for j in range(12)], 1) for w in C.M2_WINDOWS_MONTHS}


def m2_predict(model, X):
    P = X["dyn"][..., C.DYN_INPUTS.index("P")]
    acc = _acc(P)
    wlist = list(C.M2_WINDOWS_MONTHS)
    A = np.stack([acc[w] for w in wlist])                                   # (nw, B, 12)
    month = (X["t_out"] % 12)
    grp = np.broadcast_to(X["group"][:, None], month.shape)
    dec = np.broadcast_to(X["aridity"][:, None], month.shape)
    keys = {"gdm": grp * 10000 + dec * 100 + month, "gm": grp * 100 + month, "m": month}
    out = np.full(month.shape + (3, 3), np.nan, "float32")
    for t in range(3):
        done = np.zeros(month.shape, bool)
        for level in ("gdm", "gm", "m"):
            u, inv = np.unique(keys[level], return_inverse=True)
            inv = inv.reshape(month.shape)
            par = [model.get((t, level, int(k))) for k in u]
            have = np.array([p is not None for p in par])
            sel = have[inv] & ~done
            if not sel.any():
                continue
            wi = np.array([wlist.index(p[1]) if p else 0 for p in par])[inv]
            a = np.array([p[2] if p else 0.0 for p in par])[inv]
            b = np.array([p[3] if p else 0.0 for p in par])[inv]
            q = np.array([p[4] if p else [0.0] * 3 for p in par])[inv]       # (B, 12, 3)
            x = np.take_along_axis(A, wi[None], 0)[0]
            pred = (a + b * x)[..., None] + q
            out[..., t, :][sel] = pred[sel]
            done |= sel
    return out
