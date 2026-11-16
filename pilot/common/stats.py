"""Statistics used by the checks: fixed-effects OLS with cluster-robust SE, regression kriging, percentiles.

fe_ols(y, X, fe, cluster)   y on X after sweeping out any number of fixed effects (alternating projections),
                            CR1 cluster-robust covariance. Returns dict(beta, se, cov, n, df_resid, resid).
within(X, fe)               the same sweep for a matrix: what variation of X is left after the FE (Check 2.9)
RegressionKriging           z ~ drift (OLS) + ordinary kriging of the residual with an exponential variogram,
                            local neighbourhoods; prediction and kriging variance (Check 2.4, D23)
"""
import numpy as np


# ------------------------------------------------------------------------------------------- fixed effects
def _codes(g):
    _, inv = np.unique(np.asarray(g), return_inverse=True)
    return inv


def within(X, fe, tol=1e-10, max_iter=1000):
    """Sweep fixed effects (list of group arrays) out of X (n,) or (n,k) by alternating projections."""
    X = np.asarray(X, float)
    one = X.ndim == 1
    Z = X[:, None].copy() if one else X.copy()
    codes = [_codes(g) for g in fe]
    sizes = [np.bincount(c) for c in codes]
    for _ in range(max_iter):
        before = Z.copy()
        for c, s in zip(codes, sizes):
            for j in range(Z.shape[1]):
                Z[:, j] -= (np.bincount(c, weights=Z[:, j]) / s)[c]
        if np.max(np.abs(Z - before)) < tol * (1 + np.max(np.abs(before))):
            break
    return Z[:, 0] if one else Z


def fe_ols(y, X, fe=(), cluster=None, names=None):
    """OLS of y on X with fixed effects swept out; CR1 cluster-robust (or HC1 when cluster is None)."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    y, X = y[ok], X[ok]
    fe = [np.asarray(g)[ok] for g in fe]
    if fe:
        yw, Xw = within(y, fe), within(X, fe)
        n_fe = sum(len(np.unique(g)) for g in fe) - (len(fe) - 1)
    else:
        yw, Xw = y - y.mean(), X - X.mean(axis=0)
        n_fe = 1
    n, k = Xw.shape
    XtX = Xw.T @ Xw
    XtX_inv = np.linalg.pinv(XtX)
    beta = XtX_inv @ (Xw.T @ yw)
    e = yw - Xw @ beta
    dof = max(n - k - n_fe, 1)
    if cluster is None:
        meat = (Xw * e[:, None]).T @ (Xw * e[:, None])
        cov = XtX_inv @ meat @ XtX_inv * n / dof
        G = n
    else:
        g = _codes(np.asarray(cluster)[ok])
        G = g.max() + 1
        S = np.zeros((G, k))
        for j in range(k):
            S[:, j] = np.bincount(g, weights=Xw[:, j] * e, minlength=G)
        meat = S.T @ S
        cov = XtX_inv @ meat @ XtX_inv * (G / max(G - 1, 1)) * ((n - 1) / dof)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    return dict(beta=beta, se=se, cov=cov, n=int(n), n_clusters=int(G), df_resid=int(dof), resid=e,
                names=list(names) if names is not None else [f"x{j}" for j in range(k)], rows=np.where(ok)[0])


def cluster_ids(county, x=None, y=None, year=None, fallback="block", min_counties=10, block_m=10000.0):
    """Clusters for robust SE: counties when there are at least `min_counties`; otherwise (a single county gives
    one cluster and a meaningless zero SE) county-years (fallback="county_year", for county-level regressors such
    as S) or 10 km blocks (fallback="block"). Returns (ids, label)."""
    import pandas as pd
    county = np.asarray(county).astype(str)
    if len(np.unique(county)) >= min_counties:
        return county, "county"
    if fallback == "county_year":
        return (pd.Series(county) + "_" + pd.Series(np.asarray(year)).astype(str)).to_numpy(), "county-year"
    b = (pd.Series((np.asarray(x) // block_m).astype(int)).astype(str) + "_" +
         pd.Series((np.asarray(y) // block_m).astype(int)).astype(str)).to_numpy()
    return b, f"{int(block_m / 1000)} km block"


def ci95(beta, se, n_clusters=None):
    """95% interval; t with G-1 dof when clusters are few."""
    from scipy import stats
    q = stats.t.ppf(0.975, max(n_clusters - 1, 1)) if n_clusters and n_clusters < 50 else 1.959964
    return beta - q * se, beta + q * se


# ------------------------------------------------------------------------------------------- kriging
def exp_variogram(h, nugget, sill, rng):
    return nugget + sill * (1.0 - np.exp(-np.asarray(h, float) / rng))


def empirical_variogram(xy, r, max_km=60.0, n_bins=15, max_pairs=400000, seed=0):
    """Binned semivariance of residuals r at points xy (metres)."""
    rs = np.random.default_rng(seed)
    n = len(r)
    m = min(max_pairs, n * (n - 1) // 2)
    i = rs.integers(0, n, m)
    j = rs.integers(0, n, m)
    keep = i != j
    i, j = i[keep], j[keep]
    h = np.hypot(*(xy[i] - xy[j]).T) / 1000.0
    g = 0.5 * (r[i] - r[j]) ** 2
    edges = np.linspace(0, max_km, n_bins + 1)
    b = np.digitize(h, edges) - 1
    ok = (b >= 0) & (b < n_bins)
    cnt = np.bincount(b[ok], minlength=n_bins)
    gam = np.bincount(b[ok], weights=g[ok], minlength=n_bins) / np.maximum(cnt, 1)
    mid = (edges[:-1] + edges[1:]) / 2
    return mid[cnt > 0], gam[cnt > 0], cnt[cnt > 0]


def fit_exponential(h, gam, cnt):
    from scipy.optimize import least_squares
    var = float(np.average(gam, weights=cnt))
    x0 = [0.1 * var, 0.9 * var, max(h.max() / 3, 1.0)]

    def res(p):
        return np.sqrt(cnt) * (exp_variogram(h, *p) - gam)
    fit = least_squares(res, x0, bounds=([0, 1e-9, 0.5], [np.inf, np.inf, 500.0]))
    return dict(nugget=float(fit.x[0]), sill=float(fit.x[1]), range_km=float(fit.x[2]))


class RegressionKriging:
    """z = b0 + b1 * drift + e;  e kriged (ordinary, exponential variogram, k nearest within max_km).
    Prediction variance = kriging variance of e + variance of the drift prediction."""

    def __init__(self, xy, z, drift, k=32, max_km=40.0, seed=0):
        from scipy.spatial import cKDTree
        xy, z, drift = np.asarray(xy, float), np.asarray(z, float), np.asarray(drift, float)
        # wells at the same place (within 1 m) are one data point: duplicates make the kriging system singular
        key = np.round(xy).astype(np.int64)
        _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
        inv = inv.ravel()
        self.n_merged = int(len(z) - len(cnt))
        self.xy = np.column_stack([np.bincount(inv, xy[:, j]) / cnt for j in range(2)])
        z = np.bincount(inv, z) / cnt
        drift = np.bincount(inv, drift) / cnt
        X = np.column_stack([np.ones(len(z)), drift])
        self.coef, *_ = np.linalg.lstsq(X, z, rcond=None)
        self.r = z - X @ self.coef
        p = X.shape[1]
        s2 = float(self.r @ self.r / max(len(z) - p, 1))
        self.cov_coef = s2 * np.linalg.pinv(X.T @ X)
        h, g, c = empirical_variogram(self.xy, self.r, seed=seed)
        self.vario = fit_exponential(h, g, c) if len(h) >= 3 else dict(nugget=s2, sill=1e-9, range_km=1.0)
        # a small nugget keeps the system well conditioned when two wells are very close
        self.vario["nugget"] = max(self.vario["nugget"], 1e-4 * (self.vario["sill"] + self.vario["nugget"]) + 1e-6)
        self.tree = cKDTree(self.xy)
        self.k, self.max_km = min(k, len(z)), max_km

    def _gamma(self, h_km):
        return exp_variogram(h_km, **{"nugget": self.vario["nugget"], "sill": self.vario["sill"],
                                       "rng": self.vario["range_km"]})

    def _solve_batch(self, dist, idx):
        """Ordinary kriging of the residual for targets that all have k neighbours (batched solve)."""
        P = self.xy[idx]                                                    # (B, k, 2)
        H = np.sqrt(((P[:, :, None, :] - P[:, None, :, :]) ** 2).sum(-1)) / 1000.0
        B, k = idx.shape
        A = np.ones((B, k + 1, k + 1))
        A[:, :k, :k] = self._gamma(H)
        A[:, np.arange(k), np.arange(k)] = 0.0
        A[:, k, k] = 0.0
        b = np.ones((B, k + 1))
        b[:, :k] = self._gamma(dist / 1000.0)
        w = np.linalg.solve(A, b[..., None])[..., 0]
        res = (w[:, :k] * self.r[idx]).sum(axis=1)
        kvar = (w * b).sum(axis=1)
        return res, kvar

    def predict(self, xy, drift, batch=4000):
        xy = np.asarray(xy, float)
        X0 = np.column_stack([np.ones(len(xy)), np.asarray(drift, float)])
        trend = X0 @ self.coef
        trend_var = np.einsum("ij,jk,ik->i", X0, self.cov_coef, X0)
        dist, idx = self.tree.query(xy, k=self.k)
        if self.k == 1:
            dist, idx = dist[:, None], idx[:, None]
        res = np.full(len(xy), np.nan)
        kvar = np.full(len(xy), np.nan)
        c_inf = self.vario["nugget"] + self.vario["sill"]
        full = (dist <= self.max_km * 1000.0).all(axis=1)                  # all k neighbours within max_km
        fi = np.where(full)[0]
        for s in range(0, len(fi), batch):
            t = fi[s:s + batch]
            try:
                res[t], kvar[t] = self._solve_batch(dist[t], idx[t])
            except np.linalg.LinAlgError:
                full[t] = False
        for t in np.where(~full)[0]:                                        # fewer neighbours: one by one
            ok = dist[t] <= self.max_km * 1000.0
            if ok.sum() < 3:
                continue
            nb = idx[t][ok]
            P = self.xy[nb]
            m = len(nb)
            H = np.hypot(*(P[:, None, :] - P[None, :, :]).transpose(2, 0, 1)) / 1000.0
            A = np.ones((m + 1, m + 1))
            A[:m, :m] = self._gamma(H)
            np.fill_diagonal(A[:m, :m], 0.0)
            A[m, m] = 0.0
            b = np.ones(m + 1)
            b[:m] = self._gamma(dist[t][ok] / 1000.0)
            try:
                w = np.linalg.solve(A, b)
            except np.linalg.LinAlgError:
                w = np.linalg.lstsq(A, b, rcond=None)[0]
            res[t] = w[:m] @ self.r[nb]
            kvar[t] = float(w @ b)
        kvar = np.minimum(kvar, c_inf)
        return trend + res, np.sqrt(np.clip(kvar, 0, None) + trend_var)


def nearest_km(xy_from, xy_to):
    """Distance (km) from each point of xy_from to the nearest point of xy_to (inf when xy_to is empty)."""
    from scipy.spatial import cKDTree
    if len(xy_to) == 0:
        return np.full(len(xy_from), np.inf)
    d, _ = cKDTree(np.asarray(xy_to, float)).query(np.asarray(xy_from, float))
    return d / 1000.0
