"""Cells (quarter sections): geometry, the run's scope, and cell-id rasters on every source grid.

split_quarters(sections)   section polygons -> 4 quarters each, split through the section centre along the
                           section's own axes (from its minimum rotated rectangle), labelled NE NW SW SE
scope(cells)               which cells a run uses: smoke = one county, stage2 = the processing-test sample
                           (Task 1.6 selection), stage3 / test = all
strips(cells, transform, shape)   cell-id raster of a big grid (EPSG:5070 rasters) window by window
tile_ids(cells, F, cache)  cell-id raster of one 30 m UTM tile frame (cached as .npz)
The pixel -> cell rule everywhere: the cell whose polygon contains the PIXEL CENTRE (all_touched=False).
"""
import hashlib
import math
import os

import numpy as np

from .. import config as C
from . import grids

QUARTERS = ("NE", "NW", "SW", "SE")


def split_quarters(sections, id_col="FRSTDIVID"):
    """GeoDataFrame of sections (projected CRS) -> GeoDataFrame of quarters: cell_id, section_id, quarter."""
    import geopandas as gpd
    import shapely
    from shapely import affinity
    ids, secs, qs, geoms = [], [], [], []
    for sid, g in zip(sections[id_col].astype(str), sections.geometry):
        if g is None or g.is_empty:
            continue
        if not g.is_valid:
            g = g.buffer(0)
        rect = g.minimum_rotated_rectangle
        xs, ys = rect.exterior.coords.xy
        ang = math.degrees(math.atan2(ys[1] - ys[0], xs[1] - xs[0]))
        ang = ((ang + 45.0) % 90.0) - 45.0                     # the section's own axes, within +-45 deg
        c = g.centroid
        r = affinity.rotate(g, -ang, origin=c)
        x0, y0, x1, y1 = r.bounds
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2                   # quarter corners at the side midpoints
        boxes = dict(NE=(cx, cy, x1, y1), NW=(x0, cy, cx, y1), SW=(x0, y0, cx, cy), SE=(cx, y0, x1, cy))
        for q in QUARTERS:
            part = r.intersection(shapely.box(*boxes[q]))
            if part.is_empty or part.area <= 0:
                continue
            ids.append(f"{sid}{q}")
            secs.append(sid)
            qs.append(q)
            geoms.append(affinity.rotate(part, ang, origin=c))
    out = gpd.GeoDataFrame(dict(cell_id=ids, section_id=secs, quarter=qs), geometry=geoms, crs=sections.crs)
    return out.sort_values("cell_id").reset_index(drop=True)


def cells_signature(cells):
    """Short hash of the cell set (ids + rounded centroids): names the cell-id raster cache."""
    h = hashlib.md5()
    h.update("\n".join(cells["cell_id"].astype(str)).encode())
    h.update(np.round(cells["x5070"].to_numpy(), 0).tobytes())
    return h.hexdigest()[:10]


# ------------------------------------------------------------------------------------------- scope of a run
def scope(cells, sample_ids=None):
    """Boolean mask of the cells a run uses (cells: DataFrame with county, cell_id)."""
    if C.RUN_MODE == "smoke":
        return (cells["county"] == C.SMOKE_COUNTY).to_numpy()
    if C.RUN_MODE == "stage2":
        if sample_ids is None:
            raise RuntimeError("stage2 needs the processing-test sample (Task 1.6 selection)")
        return cells["cell_id"].isin(set(sample_ids)).to_numpy()
    return np.ones(len(cells), bool)


# ------------------------------------------------------------------------------------------- big grids (EPSG:5070)
def strips(cells, transform, shape, strip_rows=1024, crs="EPSG:5070"):
    """Yield (r0, r1, c0, c1, ids) windows of a big raster grid covering the cells' bounding box; ids is the
    int32 cell-id raster of the window (index = position of the cell in `cells`)."""
    from rasterio import features, windows
    g = cells.to_crs(crs) if str(cells.crs) != str(crs) else cells
    x0, y0, x1, y1 = g.total_bounds
    win = windows.from_bounds(x0, y0, x1, y1, transform=transform).round_offsets().round_lengths()
    r_a, c_a = max(int(win.row_off) - 1, 0), max(int(win.col_off) - 1, 0)
    r_b = min(int(win.row_off + win.height) + 2, shape[0])
    c_b = min(int(win.col_off + win.width) + 2, shape[1])
    sidx = g.sindex
    geoms = g.geometry.values
    for r0 in range(r_a, r_b, strip_rows):
        r1 = min(r0 + strip_rows, r_b)
        wt = windows.transform(windows.Window(c_a, r0, c_b - c_a, r1 - r0), transform)
        bx0, by1 = wt * (0, 0)
        bx1, by0 = wt * (c_b - c_a, r1 - r0)
        hit = sidx.query(__import__("shapely").box(bx0, by0, bx1, by1))
        if len(hit) == 0:
            continue
        ids = features.rasterize(((geoms[i], int(i)) for i in hit), out_shape=(r1 - r0, c_b - c_a),
                                 transform=wt, fill=-1, dtype="int32", all_touched=False)
        yield r0, r1, c_a, c_b, ids


# ------------------------------------------------------------------------------------------- 30 m UTM tiles
def tile_ids(cells, F, cache_dir=None, sig=None):
    """Cell-id raster (h, w) int32 of one tile frame F (dict crs, x0, y0, w, h). Cached when cache_dir given."""
    path = None
    if cache_dir:
        path = os.path.join(cache_dir, f"ids_{sig}_{grids.frame_key(F)}.npz")
        if os.path.exists(path):
            return np.load(path)["ids"]
    from rasterio import features
    import shapely
    g = cells.to_crs(F["crs"])
    tr = grids.frame_transform(F)
    x0, y0 = float(F["x0"]), float(F["y0"])
    box = shapely.box(x0, y0 - F["h"] * C.PIXEL_M, x0 + F["w"] * C.PIXEL_M, y0)
    hit = g.sindex.query(box)
    ids = np.full((int(F["h"]), int(F["w"])), -1, "int32")
    if len(hit):
        geoms = g.geometry.values
        ids = features.rasterize(((geoms[i], int(i)) for i in hit), out_shape=ids.shape, transform=tr, fill=-1,
                                 dtype="int32", all_touched=False)
    if path:
        os.makedirs(cache_dir, exist_ok=True)
        tmp = path + ".tmp.npz"
        np.savez_compressed(tmp, ids=ids)
        os.replace(tmp, path)
    return ids
