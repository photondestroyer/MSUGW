"""Zonal statistics over cell-id rasters with np.bincount (no per-cell loops).

A cell-id raster holds, for every source pixel, the integer index of the cell whose polygon contains the
pixel centre (-1 = no cell). Every statistic of every task is then one bincount over (ids, values).
"""
import numpy as np


def _ok(ids, *arrays):
    m = ids >= 0
    for a in arrays:
        if a is not None:
            m &= np.isfinite(a) if np.issubdtype(np.asarray(a).dtype, np.floating) else True
    return m


def count(ids, n, mask=None):
    """Number of pixels per cell (where mask, if given)."""
    ids = np.asarray(ids).ravel()
    m = ids >= 0
    if mask is not None:
        m &= np.asarray(mask, bool).ravel()
    return np.bincount(ids[m], minlength=n)[:n]


def total(ids, values, n):
    """Sum of finite values per cell, and number of finite values."""
    ids, v = np.asarray(ids).ravel(), np.asarray(values, "float64").ravel()
    m = (ids >= 0) & np.isfinite(v)
    return np.bincount(ids[m], weights=v[m], minlength=n)[:n], np.bincount(ids[m], minlength=n)[:n]


def mean(ids, values, n):
    s, k = total(ids, values, n)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(k > 0, s / k, np.nan), k


class Accumulator:
    """Sums and counts per cell across many windows (strips of a big raster)."""

    def __init__(self, n, names):
        self.n = n
        self.sum = {k: np.zeros(n) for k in names}
        self.cnt = {k: np.zeros(n, "int64") for k in names}

    def add(self, name, ids, values=None, mask=None):
        ids = np.asarray(ids).ravel()
        m = ids >= 0
        if mask is not None:
            m &= np.asarray(mask, bool).ravel()
        if values is None:
            self.cnt[name] += np.bincount(ids[m], minlength=self.n)[:self.n]
            return
        v = np.asarray(values, "float64").ravel()
        m &= np.isfinite(v)
        self.sum[name] += np.bincount(ids[m], weights=v[m], minlength=self.n)[:self.n]
        self.cnt[name] += np.bincount(ids[m], minlength=self.n)[:self.n]

    def mean(self, name):
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(self.cnt[name] > 0, self.sum[name] / self.cnt[name], np.nan)

    def share(self, num, den):
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(self.cnt[den] > 0, self.cnt[num] / self.cnt[den], np.nan)
