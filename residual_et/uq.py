"""Uncertainty: prediction intervals, clustered confidence intervals and the
per-row error budget of the residual.

The old analysis reported 1.96 * SD / sqrt(n) over pixel-dekads. Neighbouring
pixels share a 9 km soil-moisture cell and a weather system, and consecutive
dekads share a drought, so n was overstated by orders of magnitude.
"""
import numpy as np

from . import config as C


# ------------------------------------------------------------------ conformal prediction
def conformal_quantile(scores, alpha=C.ALPHA):
    """Finite-sample (1 - alpha) quantile of calibration scores."""
    s = np.sort(np.asarray(scores, "float64")[np.isfinite(scores)])
    n = s.size
    if n == 0:
        return float("nan")
    k = int(np.ceil((n + 1) * (1 - alpha))) - 1
    return float(s[min(max(k, 0), n - 1)])


def cqr_margin(y, lo, hi, alpha=C.ALPHA):
    """Conformalised quantile regression (Romano, Patterson & Candes 2019):
    the amount by which [lo, hi] must be widened (or may be narrowed) so that
    it covers 1 - alpha of held-out targets."""
    return conformal_quantile(np.maximum(lo - y, y - hi), alpha)


def abs_margin(resid, alpha=C.ALPHA):
    """Split-conformal half-width from absolute calibration residuals."""
    return conformal_quantile(np.abs(resid), alpha)


def coverage(y, lo, hi):
    ok = np.isfinite(y) & np.isfinite(lo) & np.isfinite(hi)
    return float(np.mean((y[ok] >= lo[ok]) & (y[ok] <= hi[ok]))) if ok.any() else float("nan")


def mean_width(lo, hi):
    ok = np.isfinite(lo) & np.isfinite(hi)
    return float(np.mean(hi[ok] - lo[ok])) if ok.any() else float("nan")


def interval_to_sigma(lo, hi, alpha=C.ALPHA):
    """SD of a normal distribution with the same central (1 - alpha) width."""
    z = {0.10: 1.6449, 0.05: 1.9600, 0.20: 1.2816}.get(round(alpha, 2))
    if z is None:
        from scipy.special import ndtri
        z = float(ndtri(1 - alpha / 2))
    return (hi - lo) / (2.0 * z)


# ------------------------------------------------------------------ error budget
def combine_sigma(*components):
    """Root-sum-of-squares of independent error components (NaN components are skipped).
    Independence is an assumption: regridding, climatology and model errors come
    from different data, but all three grow together in sparse or cloudy scenes."""
    tot = None
    for c in components:
        c2 = np.where(np.isfinite(c), np.asarray(c, "float64") ** 2, 0.0)
        tot = c2 if tot is None else tot + c2
    return np.sqrt(tot).astype("float32")


def observation_sigma(et_sd, et_n, step=C.SSEBOP_STEP_MM):
    """Uncertainty of a cell-mean ET value from the data handling alone:
    sampling of the cell by the available 1 km pixels (SD / sqrt(n)) and the
    whole-millimetre storage of SSEBop (uniform rounding error, step / sqrt(12))."""
    with np.errstate(invalid="ignore", divide="ignore"):
        samp = et_sd / np.sqrt(np.maximum(et_n, 1))
    return combine_sigma(samp, np.full(np.shape(et_sd), step / np.sqrt(12.0)) / np.sqrt(np.maximum(et_n, 1)))


# ------------------------------------------------------------------ clustered bootstrap
def _group_sums(values, mask, sblock, tblock, ns, nt):
    v = np.where(mask, values, 0.0).astype("float64")
    S = np.zeros((ns, nt))
    Cn = np.zeros((ns, nt))
    tb = np.broadcast_to(tblock[:, None], values.shape)
    sb = np.broadcast_to(sblock[None, :], values.shape)
    np.add.at(S, (sb[mask], tb[mask]), v[mask])
    np.add.at(Cn, (sb[mask], tb[mask]), 1.0)
    return S, Cn


def block_bootstrap(values, mask, sblock, tblock, n_boot=C.N_BOOT_EFFECT, seed=C.SEED,
                    mask_b=None, alpha=0.05):
    """Mean of values[mask] (or mean[mask] - mean[mask_b]) with a two-way
    cluster bootstrap: spatial blocks and time blocks are resampled with
    replacement, so correlation in space and in time both widen the interval.

    values (K, N); mask (K, N) bool; sblock (N,) and tblock (K,) integer ids."""
    mask = mask & np.isfinite(values)
    _, sb = np.unique(sblock, return_inverse=True)
    _, tb = np.unique(tblock, return_inverse=True)
    ns, nt = sb.max() + 1, tb.max() + 1
    Sa, Ca = _group_sums(values, mask, sb, tb, ns, nt)
    if mask_b is not None:
        mask_b = mask_b & np.isfinite(values)
        Sb, Cb = _group_sums(values, mask_b, sb, tb, ns, nt)

    def mean(w, S, Cn):
        den = (w * Cn).sum()
        return (w * S).sum() / den if den > 0 else np.nan      # a draw without data is not a zero

    def stat(ws, wt):
        w = ws[:, None] * wt[None, :]
        a = mean(w, Sa, Ca)
        return a if mask_b is None else a - mean(w, Sb, Cb)

    est = stat(np.ones(ns), np.ones(nt))
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for b in range(n_boot):
        ws = np.bincount(rng.integers(0, ns, ns), minlength=ns).astype("float64")
        wt = np.bincount(rng.integers(0, nt, nt), minlength=nt).astype("float64")
        boots[b] = stat(ws, wt)
    valid = np.isfinite(boots)
    n = int(mask.sum())
    res = dict(estimate=float(est), ci_low=float("nan"), ci_high=float("nan"), se=float("nan"), n=n,
               n_space_blocks=int((Ca.sum(1) > 0).sum()), n_time_blocks=int((Ca.sum(0) > 0).sum()),
               n_valid_draws=int(valid.sum()))
    # With data in fewer than three blocks of either kind an interval would be an artefact.
    if valid.sum() >= 0.5 * n_boot and min(res["n_space_blocks"], res["n_time_blocks"]) >= 3:
        lo, hi = np.quantile(boots[valid], [alpha / 2, 1 - alpha / 2])
        res.update(ci_low=float(lo), ci_high=float(hi), se=float(boots[valid].std()))
    return res


def naive_ci(values, mask):
    """The old interval, for comparison only: 1.96 * SD / sqrt(n)."""
    v = values[mask & np.isfinite(values)]
    if v.size < 2:
        return dict(estimate=float("nan"), half_width=float("nan"), n=int(v.size))
    return dict(estimate=float(v.mean()), half_width=float(1.96 * v.std(ddof=1) / np.sqrt(v.size)),
                n=int(v.size))
