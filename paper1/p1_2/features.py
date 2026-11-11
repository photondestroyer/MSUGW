"""Model inputs for block-water-year items (built per fold, never stored).

Item = (block i of a block set, water year wy). Water year wy = October of wy-1 .. September of wy.
Absolute month index = year * 12 + (month - 1); output months t = 12*wy - 3 .. 12*wy + 8.
    dyn     (24, len(DYN_INPUTS))  drivers of the 24 months ending in September of wy, from the block's
                                    lab cell (BUCKET from the block's capacity bin)
    stat    stat_names(): STATIC_INPUTS, one-hot FUNCTIONAL group of water year wy (the class input, plan v3),
            one-hot CDL group (secondary "species" input, SPECIES_INPUT), one-hot functional group of the
            previous calendar year (a fallow year stores water the next crop uses)
    pm      per output month, pm_names(): sin/cos of the calendar month, season clock, in-season flag (class
            calendar), pre-season NIRv and NDMI (mean of the PRESEASON_MONTHS before the season start; all
            months strictly before the target month), a flag saying whether any pre-season month was valid,
            and the measurement inputs MEAS_INPUTS of that block-month (how the value was observed: Landsat 8/9
            share, observations, mean observation day, OpenET clear scenes - rule R2 holds)
    y       ET, NIRv, NDMI of the 12 output months; ymask = target valid AND the month's calendar year is a
            reference year of the block (loss and evaluation use only those). ratio=True replaces ET by ET / PET
            of the month (the plan's sensitivity model: OpenET multiplies by gridMET reference ET).
    bound   water the block could have used in the water year: P - storage change (mm), for the mass penalty
Stratum labels (never inputs): func, func_prev, drought (0 dry / 1 between / 2 wet, from SPI90), regime
(OBS_REGIMES index of the output month), spi.
EXPOSURE_VARS can never be among the inputs: names are checked against config.EXPOSURE_VARS.
"""
import warnings

import numpy as np

from .. import config as C
from .. import cdl_groups as CG
from . import calendars as CA

PM_BASE = ["month_sin", "month_cos", "season_clock", "in_season", "pre_NIRv", "pre_NDMI", "pre_valid"]


def pm_names():
    return PM_BASE + list(C.MEAS_INPUTS)


def pm_scaled():
    """Indices of pm columns standardised by the scaler (pre-season state and measurement inputs)."""
    return [4, 5] + list(range(len(PM_BASE), len(pm_names())))


def stat_names():
    names = list(C.STATIC_INPUTS) + [f"func_{n}" for n in CG.FUNC_NAMES]
    if C.SPECIES_INPUT:
        names += [f"group_{n}" for n in CG.NAMES]
    return names + [f"prev_func_{n}" for n in CG.FUNC_NAMES]


def n_stat():
    return len(stat_names())


def n_pm():
    return len(pm_names())


# kept for callers of revision 2
PM_NAMES = pm_names()


def assert_no_exposure(names):
    bad = sorted(set(names) & set(C.EXPOSURE_VARS))
    if bad:
        raise ValueError(f"exposure variables may not be model inputs: {bad}")


def input_names():
    names = list(C.DYN_INPUTS) + stat_names() + pm_names()
    assert_no_exposure(names)
    return names


def abs_month(dt64m):
    y = dt64m.astype("datetime64[Y]").astype(int) + 1970
    m = dt64m.astype(int) % 12
    return y * 12 + m


def regime_of(t_abs):
    """OBS_REGIMES index of absolute months (year * 12 + month - 1)."""
    starts = np.array([abs_month(np.datetime64(s, "M")) for _, s in C.OBS_REGIMES])
    return np.clip(np.searchsorted(starts, np.asarray(t_abs), side="right") - 1, 0, len(starts) - 1).astype("int8")


def _year_index(years, y):
    k = np.searchsorted(years, y)
    ok = (k < years.size) & (years[np.minimum(k, years.size - 1)] == y)
    return k, ok


def set_class(X, rows, func):
    """Overwrite the functional and species one-hots of items `rows` with functional group `func` (class-flip
    test); the species becomes the first CDL group of that functional group."""
    ns = len(C.STATIC_INPUTS)
    st = X["stat"]
    st[rows, ns:ns + CG.NF] = 0
    st[rows, ns + func] = 1
    if C.SPECIES_INPUT:
        st[rows, ns + CG.NF:ns + CG.NF + CG.N] = 0
        g = [i for i, n in enumerate(CG.NAMES) if CG.FUNCTIONAL_OF[n] == CG.FUNC_NAMES[func]][0]
        st[rows, ns + CG.NF + g] = 1
    X["func"][rows] = func


def build(BS, D, items_i, items_wy, cal, no_preseason=False, no_meas=False, ratio=False):
    input_names()
    items_i = np.asarray(items_i, "int64")
    items_wy = np.asarray(items_wy, "int64")
    B = items_i.size
    ya = np.asarray(BS["years_agg"]).astype(int)
    yrs = np.asarray(BS["years"]).astype(int)
    d0 = abs_month(D["months"][0])
    t_hist = 12 * items_wy[:, None] - 15 + np.arange(24)[None, :]          # (B, 24)
    t_out = t_hist[:, 12:]                                                 # (B, 12)
    di = t_hist - d0
    if di.min() < 0 or di.max() >= D["months"].size:
        raise ValueError("driver cube does not cover the requested water years")
    cell = BS["cell"][items_i][:, None]
    dyn = np.empty((B, 24, len(C.DYN_INPUTS)), "float32")
    for j, v in enumerate(C.DYN_INPUTS):
        if v == "BUCKET":
            b = np.maximum(BS["bucket_bin"][items_i], 0)[:, None]
            dyn[:, :, j] = D["BUCKET"][di, b, cell]
        else:
            dyn[:, :, j] = D[v][di, cell]
    # targets of the output months and their reference status
    yo = t_out // 12
    mo = t_out % 12
    k, k_ok = _year_index(ya, yo)
    y = np.full((B, 12, 3), np.nan, "float32")
    ii = np.broadcast_to(items_i[:, None], yo.shape)
    y[k_ok] = BS["Y"][ii[k_ok], k[k_ok], mo[k_ok], :]
    meas = np.full((B, 12, len(C.MEAS_INPUTS)), np.nan, "float32")
    if "MEAS" in BS and len(C.MEAS_INPUTS):
        meas[k_ok] = BS["MEAS"][ii[k_ok], k[k_ok], mo[k_ok], :]
    r, r_ok = _year_index(yrs, yo)
    ref = np.zeros(yo.shape, bool)
    ref[r_ok] = BS["REF"][ii[r_ok], r[r_ok]] == 1
    pet_out = D["PET"][di[:, 12:], cell].astype("float32")
    if ratio:
        with np.errstate(invalid="ignore", divide="ignore"):
            y[..., 0] = np.where(pet_out > 1.0, y[..., 0] / pet_out, np.nan)
    ymask = np.isfinite(y) & ref[..., None]
    # class of the water year and its calendar
    rw, rw_ok = _year_index(yrs, items_wy)
    grp = np.full(B, 255, "int64")
    grp[rw_ok] = BS["GRP"][items_i[rw_ok], rw[rw_ok]]
    # a water year after the last processed calendar year (October-December of the last year) keeps the
    # class of that last year
    rp, rp_ok = _year_index(yrs, items_wy - 1)
    rp_ok &= ~rw_ok
    grp[rp_ok] = BS["GRP"][items_i[rp_ok], rp[rp_ok]]
    grp[grp == 255] = CG.OTHER
    func = CG.functional_of(grp).astype("int64")
    gprev = np.full(B, 255, "int64")
    if "GRP_PREV" in BS:
        gprev[rw_ok] = BS["GRP_PREV"][items_i[rw_ok], rw[rw_ok]]
        gprev[rp_ok] = BS["GRP"][items_i[rp_ok], rp[rp_ok]]
    func_prev = CG.functional_of(gprev).astype("int64")
    start = np.array([cal[int(g)][0] for g in grp])[:, None]
    nmon = np.array([cal[int(g)][1] for g in grp])[:, None]
    der = CA.derived(t_out, start, nmon)
    # pre-season NIRv / NDMI from the block's own earlier months
    pre = der["preseason"]                                                  # (B, 12, P)
    py, pm_ = pre // 12, pre % 12
    pk, pk_ok = _year_index(ya, py)
    vals = np.full(pre.shape + (2,), np.nan, "float32")
    iii = np.broadcast_to(items_i[:, None, None], pre.shape)
    vals[pk_ok] = BS["Y"][iii[pk_ok], pk[pk_ok], pm_[pk_ok], 1:3]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)           # all-NaN pre-season windows -> NaN
        pre_mean = np.nanmean(vals, axis=2) if vals.size else np.full((B, 12, 2), np.nan)
    pre_valid = np.isfinite(pre_mean).all(-1)
    if no_preseason:
        pre_mean[:] = np.nan
        pre_valid[:] = False
    if no_meas:
        meas[:] = np.nan
    ang = 2 * np.pi * mo / 12.0
    pm = np.concatenate([np.stack([np.sin(ang), np.cos(ang), der["clock"], der["in_season"].astype("float32"),
                                   pre_mean[..., 0], pre_mean[..., 1], pre_valid.astype("float32")], -1),
                         meas], -1).astype("float32")
    onehot = lambda idx, n: np.eye(n, dtype="float32")[np.asarray(idx)]
    parts = [BS["static"][items_i].astype("float32"), onehot(func, CG.NF)]
    if C.SPECIES_INPUT:
        parts.append(onehot(grp, CG.N))
    parts.append(onehot(func_prev, CG.NF))
    stat = np.concatenate(parts, 1)
    # water available in the water year: P - storage change (end Sep wy minus end Sep wy-1)
    caps = np.asarray(C.BUCKET_CAPACITY_BINS_MM, "float64")
    b = np.maximum(BS["bucket_bin"][items_i], 0)
    P_wy = D["P"][di[:, 12:], cell].sum(1)
    s_end = D["BUCKET"][di[:, -1], b, cell[:, 0]]
    s_start = D["BUCKET"][di[:, 11], b, cell[:, 0]]
    bound = (P_wy - (s_end - s_start) * caps[b]).astype("float32")
    spi = D["SPI90"][di[:, 12:], cell].astype("float32") if "SPI90" in D else np.full((B, 12), np.nan, "float32")
    from ..step0.drivers import drought_class
    return dict(dyn=dyn, stat=stat, pm=pm, y=y, ymask=ymask, bound=bound, t_out=t_out, in_season=der["in_season"],
                group=grp, func=func, func_prev=func_prev, aridity=BS["aridity"][items_i], i=items_i, wy=items_wy,
                spi=spi, drought=drought_class(spi), regime=regime_of(t_out), pet=pet_out, ratio=ratio)


# ------------------------------------------------------------------ scaling (fit on training items only)
def fit_scaler(X):
    def ms(a, axes):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)      # an all-missing input gets mean 0, SD 1
            mu = np.nanmean(a, axis=axes)
            sd = np.nanstd(a, axis=axes)
        sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
        return np.nan_to_num(mu).astype("float32"), sd.astype("float32")
    sc = {}
    sc["dyn"] = ms(X["dyn"], (0, 1))
    sc["stat"] = ms(X["stat"][:, :len(C.STATIC_INPUTS)], 0)
    sc["pm"] = ms(X["pm"][..., pm_scaled()], (0, 1))
    y = np.where(X["ymask"], X["y"], np.nan)
    sc["y"] = ms(y, (0, 1))
    return sc


def apply_scaler(X, sc):
    out = dict(X)
    out["dyn"] = np.nan_to_num((X["dyn"] - sc["dyn"][0]) / sc["dyn"][1]).astype("float32")
    st = X["stat"].copy()
    ns = len(C.STATIC_INPUTS)
    st[:, :ns] = (st[:, :ns] - sc["stat"][0]) / sc["stat"][1]
    out["stat"] = np.nan_to_num(st).astype("float32")
    pm = X["pm"].copy()
    idx = pm_scaled()
    pm[..., idx] = (pm[..., idx] - sc["pm"][0]) / sc["pm"][1]
    out["pm"] = np.nan_to_num(pm).astype("float32")
    out["y_std"] = np.nan_to_num((X["y"] - sc["y"][0]) / sc["y"][1]).astype("float32")
    return out


def scaler_to_json(sc):
    return {k: [v[0].tolist(), v[1].tolist()] for k, v in sc.items()}


def scaler_from_json(d):
    return {k: (np.asarray(v[0], "float32"), np.asarray(v[1], "float32")) for k, v in d.items()}
