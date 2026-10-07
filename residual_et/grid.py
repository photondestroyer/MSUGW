"""Analysis grid and area aggregation.

The analysis grid is the lab's native gridMET grid (1/24 degree). Finer
sources are aggregated to it by area and return the within-cell spread and the
number of contributing pixels, so the regridding step carries an uncertainty.
The old scripts took one point per cell (bilinear or nearest).
"""
from dataclasses import dataclass

import numpy as np

from . import config as C
from . import sources as S


@dataclass
class Grid:
    lat: np.ndarray        # (H,) cell centres, descending
    lon: np.ndarray        # (W,) cell centres, ascending
    frac: np.ndarray       # (H, W) share of each cell inside the aquifer polygon
    iy: np.ndarray         # (N,) row of each aquifer cell
    ix: np.ndarray         # (N,) column of each aquifer cell
    ysl: slice             # bounding window of the cells, in grid rows
    xsl: slice

    @property
    def n(self):
        return int(self.iy.size)

    @property
    def dlat(self):
        return float(self.lat[0] - self.lat[1])

    @property
    def dlon(self):
        return float(self.lon[1] - self.lon[0])

    @property
    def lat_edges(self):       # (H+1,) descending
        return np.append(self.lat + self.dlat / 2, self.lat[-1] - self.dlat / 2)

    @property
    def lon_edges(self):       # (W+1,) ascending
        return np.append(self.lon - self.dlon / 2, self.lon[-1] + self.dlon / 2)

    @property
    def cell_lat(self):
        return self.lat[self.iy]

    @property
    def cell_lon(self):
        return self.lon[self.ix]

    def lut(self):
        """(H, W) int32: aquifer-cell number, -1 elsewhere."""
        out = np.full(self.frac.shape, -1, "int32")
        out[self.iy, self.ix] = np.arange(self.n, dtype="int32")
        return out

    def window_cells(self):
        """Position of each cell inside the bounding window (for hyperslab reads)."""
        return self.iy - self.ysl.start, self.ix - self.xsl.start

    def to_map(self, values, fill=np.nan):
        out = np.full(self.frac.shape, fill, dtype="float32")
        out[self.iy, self.ix] = values
        return out

    def space_blocks(self, deg=C.SPACE_BLOCK_DEG):
        """Integer id of the lat/lon box each cell falls in (bootstrap clusters)."""
        by = np.floor(self.cell_lat / deg).astype(int)
        bx = np.floor(self.cell_lon / deg).astype(int)
        _, inv = np.unique(by * 100000 + bx, return_inverse=True)
        return inv.astype("int32")


def make_grid(lat, lon, frac, bbox=None, min_fraction=C.ROI_MIN_FRACTION):
    keep = frac >= min_fraction
    if bbox is not None:
        lon0, lon1, lat0, lat1 = bbox
        keep &= ((lat[:, None] >= lat0) & (lat[:, None] <= lat1)
                 & (lon[None, :] >= lon0) & (lon[None, :] <= lon1))
    iy, ix = np.where(keep)
    if iy.size == 0:
        raise ValueError("no aquifer cells selected; check the bounding box")
    return Grid(lat, lon, frac, iy.astype("int32"), ix.astype("int32"),
                slice(int(iy.min()), int(iy.max()) + 1), slice(int(ix.min()), int(ix.max()) + 1))


def load_grid(profile):
    """Grid, cells and aquifer fraction straight from the lab gridMET file."""
    ds = S.open_ro(C.SRC["gridmet"])
    try:
        lat = np.asarray(ds["lat"][:], dtype="float64")
        lon = np.asarray(ds["lon"][:], dtype="float64")
        frac = np.asarray(ds["roi_fraction"][:], dtype="float32")
    finally:
        ds.close()
    return make_grid(lat, lon, frac, profile.bbox)


# ------------------------------------------------------------------ aggregation
def rows_cols(src_lat, src_lon, grid):
    """Grid row of each source latitude and grid column of each source longitude (-1 = outside)."""
    le, lo = grid.lat_edges, grid.lon_edges
    ri = np.floor((le[0] - np.asarray(src_lat, "float64")) / grid.dlat).astype("int64")
    ci = np.floor((np.asarray(src_lon, "float64") - lo[0]) / grid.dlon).astype("int64")
    ri[(ri < 0) | (ri >= grid.lat.size)] = -1
    ci[(ci < 0) | (ci >= grid.lon.size)] = -1
    return ri, ci


def cell_ids(src_lat, src_lon, grid, lut=None):
    """(ny, nx) aquifer-cell number of every pixel of a regular lat/lon raster (-1 = none)."""
    lut = grid.lut() if lut is None else lut
    ri, ci = rows_cols(src_lat, src_lon, grid)
    out = np.full((ri.size, ci.size), -1, "int32")
    ok_r, ok_c = ri >= 0, ci >= 0
    out[np.ix_(ok_r, ok_c)] = lut[np.ix_(ri[ok_r], ci[ok_c])]
    return out


def cell_ids_points(lat, lon, grid, lut=None):
    """Aquifer-cell number of arbitrary points (same shape as lat/lon)."""
    lut = grid.lut() if lut is None else lut
    le, lo = grid.lat_edges, grid.lon_edges
    ri = np.floor((le[0] - lat) / grid.dlat).astype("int64")
    ci = np.floor((lon - lo[0]) / grid.dlon).astype("int64")
    ok = (ri >= 0) & (ri < grid.lat.size) & (ci >= 0) & (ci < grid.lon.size)
    out = np.full(np.shape(lat), -1, "int32")
    out[ok] = lut[ri[ok], ci[ok]]
    return out


class Accumulator:
    """Streaming per-cell sum, sum of squares and count."""

    def __init__(self, n):
        self.n = n
        self.s = np.zeros(n, "float64")
        self.ss = np.zeros(n, "float64")
        self.c = np.zeros(n, "float64")

    def add(self, cid, values):
        cid = np.ravel(cid)
        v = np.ravel(values).astype("float64")
        ok = (cid >= 0) & np.isfinite(v)
        if not ok.any():
            return
        cid, v = cid[ok], v[ok]
        self.s += np.bincount(cid, weights=v, minlength=self.n)
        self.ss += np.bincount(cid, weights=v * v, minlength=self.n)
        self.c += np.bincount(cid, minlength=self.n)

    def result(self):
        """mean, within-cell SD (population), count; NaN where no pixel contributed."""
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(self.c > 0, self.s / self.c, np.nan)
            var = np.where(self.c > 0, self.ss / self.c - mean ** 2, np.nan)
        return (mean.astype("float32"), np.sqrt(np.maximum(var, 0)).astype("float32"),
                self.c.astype("int32"))


def block_stats(values, cid, n):
    """mean, SD, count per cell for one 2-D field."""
    acc = Accumulator(n)
    acc.add(cid, values)
    return acc.result()


def class_fractions(classes, cid, n, class_values):
    """Area share of each class per cell: (len(class_values), n), plus valid-pixel count.
    Pixels whose code is not in `class_values` are counted in the denominator only if >= 0."""
    cid = np.ravel(cid)
    cls = np.ravel(classes)
    ok = (cid >= 0) & (cls >= 0)
    cid, cls = cid[ok], cls[ok]
    total = np.bincount(cid, minlength=n).astype("float64")
    out = np.zeros((len(class_values), n), "float32")
    for k, cv in enumerate(class_values):
        m = np.isin(cls, cv)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[k] = np.where(total > 0, np.bincount(cid[m], minlength=n) / total, np.nan)
    return out, total.astype("int32")


def nearest_source(src_lat, src_lon, cell_lat, cell_lon):
    """Index of the nearest source point (flattened) for each cell, and the distance in degrees
    of latitude-equivalent. Used only for sources coarser than the grid (SMAP, irrigation)."""
    from scipy.spatial import cKDTree
    k = np.cos(np.deg2rad(np.mean(cell_lat)))
    tree = cKDTree(np.column_stack([np.ravel(src_lat), np.ravel(src_lon) * k]))
    dist, idx = tree.query(np.column_stack([cell_lat, cell_lon * k]))
    return idx.astype("int64"), dist.astype("float32")
