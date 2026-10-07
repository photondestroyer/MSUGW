"""Dekads, climatologies, anomalies and standardised rainfall anomalies.

Rules that fix the old code:
  * every anomaly is per cell and per dekad-of-year, from a baseline that
    contains no validation or test data;
  * antecedent windows end the day BEFORE the dekad starts and are computed
    from the daily series, never from a truncated batch;
  * rainfall enters the models relative to a 30-year climatology.

Only `standardize_rainfall` needs scipy; the rest runs in the LSTM environment.
"""
import numpy as np

from . import config as C

DAY = np.timedelta64(1, "D")


# ------------------------------------------------------------------ dekad calendar
def dek36(days):
    """Dekad-of-year 0..35 (three per month: days 1-10, 11-20, 21-end)."""
    d = np.asarray(days, "datetime64[D]")
    m = d.astype("datetime64[M]")
    month0 = (m - d.astype("datetime64[Y]").astype("datetime64[M]")).astype(int)
    dom = (d - m.astype("datetime64[D]")).astype(int)          # 0-based day of month
    return (month0 * 3 + np.minimum(dom // 10, 2)).astype("int16")


def dekad_starts(t0, t1):
    """All dekad start dates in [t0, t1)."""
    months = np.arange(np.datetime64(t0, "M"), np.datetime64(t1, "M") + 1)
    starts = (months[:, None].astype("datetime64[D]") + np.array([0, 10, 20])).ravel()
    return starts[(starts >= np.datetime64(t0, "D")) & (starts < np.datetime64(t1, "D"))]


def dekad_end(starts):
    """Exclusive end of each dekad (= start of the next one)."""
    s = np.asarray(starts, "datetime64[D]")
    m = s.astype("datetime64[M]")
    dom = (s - m.astype("datetime64[D]")).astype(int)
    nxt = (m + 1).astype("datetime64[D]")
    return np.where(dom == 20, nxt, s + 10 * DAY)


def doy(days):
    """0-based day of year (0..365)."""
    d = np.asarray(days, "datetime64[D]")
    return (d - d.astype("datetime64[Y]").astype("datetime64[D]")).astype(int)


# ------------------------------------------------------------------ daily -> dekad
def aggregate_dekads(days, arr, starts, how, min_fraction=C.MIN_VALID_FRACTION):
    """Dekad mean or sum of a daily (T, N) series. A dekad with fewer than
    `min_fraction` valid days is NaN; sums are scaled to the full dekad length."""
    days = np.asarray(days, "datetime64[D]")
    ends = dekad_end(starts)
    out = np.full((len(starts), arr.shape[1]), np.nan, "float32")
    for k, (a, b) in enumerate(zip(starts, ends)):
        i0, i1 = np.searchsorted(days, a), np.searchsorted(days, b)
        length = int((b - a) / DAY)
        if i1 - i0 != length:
            continue                                   # dekad not fully inside the daily series
        blk = arr[i0:i1]
        ok = np.isfinite(blk)
        cnt = ok.sum(axis=0)
        tot = np.where(ok, blk, 0.0).sum(axis=0, dtype="float64")
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = tot / cnt
        val = mean * length if how == "sum" else mean
        out[k] = np.where(cnt >= min_fraction * length, val, np.nan)
    return out


def trailing_sum(days, arr, targets, window, min_fraction=C.MIN_VALID_FRACTION):
    """Sum of the `window` days that END THE DAY BEFORE each target date.
    NaN when the window is not fully inside `days` or too few days are valid."""
    days = np.asarray(days, "datetime64[D]")
    ok = np.isfinite(arr)
    cs = np.vstack([np.zeros((1, arr.shape[1])), np.cumsum(np.where(ok, arr, 0.0), axis=0, dtype="float64")])
    cn = np.vstack([np.zeros((1, arr.shape[1]), "int64"), np.cumsum(ok, axis=0)])
    out = np.full((len(targets), arr.shape[1]), np.nan, "float32")
    for k, t in enumerate(np.asarray(targets, "datetime64[D]")):
        i1 = int(np.searchsorted(days, t))             # index of the target day = one past the window
        i0 = i1 - window
        if i0 < 0 or i1 > len(days) or days[i0] != t - window * DAY:
            continue
        n = cn[i1] - cn[i0]
        with np.errstate(invalid="ignore", divide="ignore"):
            val = (cs[i1] - cs[i0]) * window / n
        out[k] = np.where(n >= min_fraction * window, val, np.nan)
    return out


# ------------------------------------------------------------------ climatology
class DekadClimatology:
    """Streaming per-cell statistics for each dekad-of-year."""

    def __init__(self, n):
        self.s = np.zeros((36, n), "float64")
        self.ss = np.zeros((36, n), "float64")
        self.c = np.zeros((36, n), "float64")

    def add(self, dek, values):
        """dek (T,), values (T, N)."""
        ok = np.isfinite(values)
        v = np.where(ok, values, 0.0).astype("float64")
        np.add.at(self.s, dek, v)
        np.add.at(self.ss, dek, v * v)
        np.add.at(self.c, dek, ok.astype("float64"))

    def finalize(self, pool=1, sd_floor=0.0):
        """mean, SD, standard error of the mean, number of years.
        `pool` neighbouring dekads on each side sharpen the SD estimate; the
        standard error uses the number of YEARS, because pooled dekads of one
        year are not independent samples."""
        q = np.arange(36)
        s, ss, c = self.s.copy(), self.ss.copy(), self.c.copy()
        for k in range(1, pool + 1):
            for sh in (-k, k):
                s += self.s[(q + sh) % 36]
                ss += self.ss[(q + sh) % 36]
                c += self.c[(q + sh) % 36]
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(c > 0, s / c, np.nan)
            var = np.where(c > 1, (ss - c * mean ** 2) / (c - 1), np.nan)
        sd = np.sqrt(np.maximum(var, 0))
        years = self.c
        with np.errstate(invalid="ignore", divide="ignore"):
            se = np.where(years > 0, sd / np.sqrt(years), np.nan)
        return (mean.astype("float32"), np.maximum(sd, sd_floor).astype("float32"),
                se.astype("float32"), years.astype("int16"))


def anomaly(values, dek, mean):
    return (values - mean[dek]).astype("float32")


def std_anomaly(values, dek, mean, sd, floor=1e-6):
    return ((values - mean[dek]) / np.maximum(sd[dek], floor)).astype("float32")


# ------------------------------------------------------------------ rainfall anomaly
def fit_gamma_zero(x, min_positive=8):
    """Two-parameter gamma with a probability mass at zero, per column of x (S, N).
    Thom's closed-form maximum-likelihood approximation, so the fit is vectorised
    over cells. Returns shape, scale, number of zeros, number of samples."""
    ok = np.isfinite(x)
    pos = ok & (x > 0)
    n = ok.sum(axis=0).astype("float64")
    npos = pos.sum(axis=0).astype("float64")
    xs = np.where(pos, x, 1.0).astype("float64")
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(pos, xs, 0.0).sum(axis=0) / npos
        mlog = np.where(pos, np.log(xs), 0.0).sum(axis=0) / npos
        a = np.maximum(np.log(mean) - mlog, 1e-6)
        shape = (1.0 + np.sqrt(1.0 + 4.0 * a / 3.0)) / (4.0 * a)
        scale = mean / shape
    bad = npos < min_positive
    shape[bad], scale[bad] = np.nan, np.nan
    return shape, scale, n - npos, n


def gamma_to_normal(x, shape, scale, nzero, n, p_floor=1e-4):
    """Standard-normal score of x under the fitted mixed distribution."""
    from scipy.special import gammainc, ndtri
    q0 = nzero / n
    with np.errstate(invalid="ignore", divide="ignore"):
        g = gammainc(shape, np.maximum(x, 0.0) / scale)
        p = q0 + (1.0 - q0) * g
        p0 = (nzero + 1.0) / (2.0 * (n + 1.0))        # centre of the zero mass (Stagge et al. 2015)
        p = np.where(x <= 0, p0, p)
    return ndtri(np.clip(p, p_floor, 1.0 - p_floor)).astype("float32")


def standardize_rainfall(base, values, dek, pool=1, n_boot=0, seed=C.SEED):
    """Standardised rainfall anomaly ("in-house SPI") of accumulated rainfall.

    base   (Y, 36, N)  accumulations at each dekad start in the baseline years
    values (K, N)      accumulations to transform; dek (K,) their dekad-of-year
    Fits use the baseline only, pooled over +-`pool` dekads. Returns the score
    and, when n_boot > 0, its bootstrap SD over resampled baseline years."""
    out = np.full(values.shape, np.nan, "float32")
    acc = np.zeros(values.shape, "float64")
    acc2 = np.zeros(values.shape, "float64")
    rng = np.random.default_rng(seed)
    years = base.shape[0]
    draws = [rng.integers(0, years, years) for _ in range(n_boot)]
    for q in range(36):
        rows = np.where(dek == q)[0]
        if rows.size == 0:
            continue
        qs = [(q + s) % 36 for s in range(-pool, pool + 1)]
        samp = base[:, qs, :]                           # (Y, 2*pool+1, N)
        fit = fit_gamma_zero(samp.reshape(-1, samp.shape[-1]))
        out[rows] = gamma_to_normal(values[rows], *fit)
        for d in draws:
            z = gamma_to_normal(values[rows], *fit_gamma_zero(samp[d].reshape(-1, samp.shape[-1])))
            acc[rows] += z
            acc2[rows] += z.astype("float64") ** 2
    if n_boot == 0:
        return out, None
    with np.errstate(invalid="ignore"):
        se = np.sqrt(np.maximum(acc2 / n_boot - (acc / n_boot) ** 2, 0.0))
    return out, se.astype("float32")


def percentile_of(values, dek, sample, sample_dek, pool=1):
    """Empirical percentile (0-100) of each value within the cell's own
    dekad-of-year sample (pooled over +-`pool` dekads)."""
    out = np.full(values.shape, np.nan, "float32")
    for q in range(36):
        rows = np.where(dek == q)[0]
        if rows.size == 0:
            continue
        near = np.isin(sample_dek, [(q + s) % 36 for s in range(-pool, pool + 1)])
        ref = sample[near]                               # (S, N)
        n = np.isfinite(ref).sum(axis=0)
        for r in rows:
            v = values[r]
            with np.errstate(invalid="ignore", divide="ignore"):
                pct = 100.0 * (np.nansum(ref <= v[None, :], axis=0)) / n
            out[r] = np.where((n >= 5) & np.isfinite(v), pct, np.nan)
    return out
