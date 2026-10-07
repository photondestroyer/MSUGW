"""Drought conditions from published definitions, with checks that each
definition is applied as published.

The stored lab indices are rank-based standardised indices (Farahmand &
AghaKouchak 2015) on a 1981-2016 reference: about 36 distinct values bounded at
+-2.09, posted every 5 days, each "aggregated over the last N days".
"""
import numpy as np

from . import config as C
from .sources import SourceError


# ------------------------------------------------------------------ time lookup
def lookup_preceding(index_days, starts, max_gap_days=C.PENTAD_MAX_GAP_DAYS):
    """Record of the last index date STRICTLY BEFORE each dekad start.

    An index dated d aggregates the days up to and including d, so it is known
    before a dekad that starts on d+1 or later. Raises on a gap: the old
    "nearest step" lookup silently reached across a missing year."""
    index_days = np.asarray(index_days, "datetime64[D]")
    starts = np.asarray(starts, "datetime64[D]")
    j = np.searchsorted(index_days, starts, side="left") - 1
    if (j < 0).any():
        raise SourceError(f"drought index starts {index_days[0]}, after dekad {starts[j < 0][0]}")
    gap = (starts - index_days[j]).astype(int)
    if gap.max() > max_gap_days:
        k = int(np.argmax(gap > max_gap_days))          # the first dekad left without an index
        raise SourceError(f"no drought index within {max_gap_days} days before {starts[k]} "
                          f"(last one is {index_days[j[k]]}, {int(gap[k])} days earlier)")
    return j


# ------------------------------------------------------------------ classes
def usdm_class(v):
    """US Drought Monitor class of an SPI/SPEI value: 0 none, 1..5 = D0..D4."""
    v = np.asarray(v, "float32")
    out = np.zeros(v.shape, "int8")
    for k, (_, cut) in enumerate(C.USDM_CUTS, start=1):
        out[v <= cut] = k
    out[~np.isfinite(v)] = -1
    return out


def lab_class(v):
    """Same scale, using the cut points written in the lab file's category attribute."""
    v = np.asarray(v, "float32")
    out = np.zeros(v.shape, "int8")
    cuts = C.LAB_CATEGORY_CUTS
    out[v <= cuts[0]] = 1
    for k in range(1, 4):
        out[v <= cuts[k]] = k + 1
    out[v < cuts[4]] = 5                    # the lab's top class is strictly below -2
    out[~np.isfinite(v)] = -1
    return out


def flag(name, indices):
    """Boolean field for one named definition in config.DROUGHT_DEFS."""
    d = C.DROUGHT_DEFS[name]
    v = indices[d["index"]]
    with np.errstate(invalid="ignore"):
        return (v <= d["threshold"]) if d["op"] == "le" else (v >= d["threshold"])


def flash_drought(sm_pct):
    """Ford & Labosier (2017): soil-moisture percentile at or above the start
    level and at or below the end level within the following two dekads
    (20 days). sm_pct (K, N) in dekad order; True on the dekad the decline ends."""
    f = C.FLASH
    out = np.zeros(sm_pct.shape, bool)
    with np.errstate(invalid="ignore"):
        low = sm_pct <= f["end_pct"]
        for lag in range(1, f["within_dekads"] + 1):
            out[lag:] |= low[lag:] & (sm_pct[:-lag] >= f["start_pct"])
    return out


# ------------------------------------------------------------------ definition checks
def check_long_names(cube):
    """The file must itself state the accumulation window each rule assumes."""
    res = {}
    for name in C.LAB_INDICES:
        days = int("".join(ch for ch in name if ch.isdigit()))
        ln = cube.long_name(name)
        ok = f"aggregated over the last {days} days" in ln
        res[name] = dict(ok=ok, long_name=ln, expected_window_days=days)
        if not ok:
            raise SourceError(f"{name}: long_name {ln!r} does not state a {days}-day window")
    for rule, d in C.DROUGHT_DEFS.items():
        digits = int("".join(ch for ch in d["index"] if ch.isdigit()))
        if digits != d["window_days"]:
            raise SourceError(f"{rule}: index {d['index']} is not a {d['window_days']}-day index")
    return res


def check_usdm_table():
    """Class cuts in config must equal the published US Drought Monitor table."""
    published = (("D0", -0.5), ("D1", -0.8), ("D2", -1.3), ("D3", -1.6), ("D4", -2.0))
    if tuple(C.USDM_CUTS) != published:
        raise SourceError(f"USDM cut points {C.USDM_CUTS} differ from the published {published}")
    if C.DROUGHT_DEFS["usdm_d1"]["threshold"] != dict(published)["D1"]:
        raise SourceError("usdm_d1 threshold is not the published D1 bound")
    return dict(ok=True, cuts=dict(published))


def check_signs(indices, p30_anom, etr_anom):
    """Sign conventions, tested on data: SPI rises with rainfall; SPEI follows
    SPI; EDDI rises with evaporative demand and is opposite to SPEI."""
    def r(a, b):
        m = np.isfinite(a) & np.isfinite(b)
        return float(np.corrcoef(a[m], b[m])[0, 1]) if m.sum() > 50 else float("nan")
    res = dict(spi30_vs_rain30=r(indices["spi30d"], p30_anom),
               spei30_vs_spi30=r(indices["spei30d"], indices["spi30d"]),
               eddi30_vs_etr=r(indices["eddi30d"], etr_anom),
               eddi30_vs_spei30=r(indices["eddi30d"], indices["spei30d"]))
    expect = dict(spi30_vs_rain30=1, spei30_vs_spi30=1, eddi30_vs_etr=1, eddi30_vs_spei30=-1)
    for k, sign in expect.items():
        if np.isfinite(res[k]) and np.sign(res[k]) != sign:
            raise SourceError(f"sign check failed: corr {k} = {res[k]:+.2f}, expected sign {sign:+d}")
    res["ok"] = True
    return res


def describe_index(v):
    """What the stored index looks like: discrete levels, bounds, saturation."""
    a = v[np.isfinite(v)]
    lv = np.unique(np.round(a, 2))
    return dict(n_levels=int(lv.size), minimum=float(lv[0]), maximum=float(lv[-1]),
                share_at_minimum=float(np.mean(a <= lv[0] + 1e-6)),
                share_at_maximum=float(np.mean(a >= lv[-1] - 1e-6)),
                share_le_minus1=float(np.mean(a <= -1.0)), share_le_minus08=float(np.mean(a <= -0.8)))


def compare_class_schemes(v):
    """Share of values that the lab's category cuts and the USDM cuts classify differently."""
    a = v[np.isfinite(v)]
    u, lab = usdm_class(a), lab_class(a)
    diff = u != lab
    lv = np.unique(np.round(a[diff], 2))
    return dict(share_classified_differently=float(diff.mean()), levels_affected=lv.tolist()[:12])


def rank_agreement(a, b):
    """Spearman rank correlation of two fields (finite pairs only)."""
    from scipy.stats import spearmanr
    m = np.isfinite(a) & np.isfinite(b)
    return float(spearmanr(a[m], b[m]).statistic) if m.sum() > 50 else float("nan")
