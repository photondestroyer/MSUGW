"""Class-specific growing-season calendars and the inputs derived from them (user decision: class-specific).

For each CDL class group, the calendar is read from the climatological NIRv curve of REFERENCE blocks of
that group in the TRAINING part of a fold (never from test blocks or test years):
    curve[m] = median NIRv of training reference block-months with calendar month m (m = 0 .. 11)
    (start_month, n_months) = season_window(curve)
From a calendar, for any target month t (absolute month index):
    a          = latest month <= t whose calendar month is start_month (the season that contains t, or
                 the last one that began before it)
    pre-season = months a-3 .. a-1 (PRESEASON_MONTHS), all strictly before t
    clock      = (t - a) / n_months, in_season = (t - a) < n_months
A group with too few training reference block-months borrows the calendar of all natural groups pooled.
"""
import numpy as np

from .. import config as C
from .. import cdl_groups as CG

MIN_SAMPLES = 200


def season_window(curve):
    """Growing season of one class from its 12-month climatological NIRv curve.

    curve : float array of 12 values, curve[0] = January ... curve[11] = December (NaN = no data).
    Returns (start_month, n_months): start_month in 1..12, n_months in 1..12. The season may wrap
    over the new year (e.g. winter wheat: start 10, n_months 9 = October .. June).
    """
    # TODO(human): choose the growing-season months from the curve (see the Learn by Doing request).
    raise NotImplementedError("calendars.season_window is not written yet (TODO(human))")


def fit(curves_by_group):
    """{group: curve or None} -> {group: (start, n)}; groups without a curve use the pooled natural one."""
    out = {}
    pooled = curves_by_group.get("_natural")
    for gi in range(CG.N):
        c = curves_by_group.get(gi)
        if c is None:
            c = pooled
        out[gi] = season_window(np.asarray(c, "float64")) if c is not None else (4, 7)
    return out


def curves_from_training(nirv, month, group, ok):
    """Median NIRv per calendar month and group from training reference block-months.
    nirv, month (1-12), group: 1-D arrays of block-months; ok: bool mask of usable training block-months."""
    curves = {}
    for key, sel in [(gi, group == gi) for gi in range(CG.N)] + [("_natural", np.isin(group, CG.NATURAL_GROUPS))]:
        m = ok & sel & np.isfinite(nirv)
        if m.sum() < MIN_SAMPLES:
            continue
        curves[key] = np.array([np.median(nirv[m & (month == k + 1)]) if (m & (month == k + 1)).any() else np.nan
                                for k in range(12)])
    return curves


def latest_start(t_abs, start_month):
    """Absolute index of the latest month <= t_abs whose calendar month is start_month.
    Absolute month index = year * 12 + (month - 1)."""
    cm = t_abs % 12 + 1
    back = (cm - start_month) % 12
    return t_abs - back


def derived(t_abs, start_month, n_months):
    """Arrays for target months t_abs (any shape): season start a, clock, in_season, pre-season months."""
    t_abs = np.asarray(t_abs)
    a = latest_start(t_abs, np.asarray(start_month))
    since = t_abs - a
    clock = since / np.asarray(n_months, "float64")
    in_season = since < np.asarray(n_months)
    pre = np.stack([a - k for k in range(C.PRESEASON_MONTHS, 0, -1)], -1)   # a-3, a-2, a-1
    return dict(start=a, clock=clock.astype("float32"), in_season=in_season, preseason=pre)


def to_json(cal):
    return {CG.NAMES[g]: dict(start_month=int(s), n_months=int(n)) for g, (s, n) in cal.items()}


def from_json(d):
    return {CG.NAMES.index(k): (v["start_month"], v["n_months"]) for k, v in d.items()}
