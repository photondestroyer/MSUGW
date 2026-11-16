"""Grid arithmetic: EASE-Grid 2.0 (SMAP), the 30 m UTM tile lattice (CDL, OpenET), the gridMET lab grid.

EASE-Grid 2.0 global (EPSG:6933): column c, row r of a grid with cell size s have centre
    x = X_MIN + (c + 0.5) s,   y = Y_MAX - (r + 0.5) s.
The 9 km grid nests exactly 4 x 4 in the 36 km grid: r36 = r9 // 4, c36 = c9 // 4.

Tile key z{zone}e{i:02d}n{j:02d}: EPSG:326{zone}; tile corner left = 15 + i*153600, top = 15 + (j+1)*153600.
A product's frame (frames.json) is the tile cropped to its aquifer pixels: pixel (0,0) of a file sits at
(x0, y0) = upper-left corner, and spans w x h pixels of 30 m.
"""
import json
import math
import os

import numpy as np

from .. import config as C

TILE_M = C.TILE_PX * C.PIXEL_M


# ---------------------------------------------------------------- EASE-Grid 2.0
def ease_index(x, y, size):
    """(row, col) of the EASE cell holding EPSG:6933 point(s) (x, y) on a grid of cell size `size`."""
    col = np.floor((np.asarray(x, float) - C.EASE_X_MIN) / size).astype(int)
    row = np.floor((C.EASE_Y_MAX - np.asarray(y, float)) / size).astype(int)
    return row, col


def ease_centre(row, col, size):
    return C.EASE_X_MIN + (np.asarray(col) + 0.5) * size, C.EASE_Y_MAX - (np.asarray(row) + 0.5) * size


def ease_box(row, col, size):
    """(left, bottom, right, top) of one EASE cell in EPSG:6933."""
    left = C.EASE_X_MIN + col * size
    top = C.EASE_Y_MAX - row * size
    return left, top - size, left + size, top


def fp_id(r36, c36):
    """Footprint id (36 km EASE cell): 'R{row}C{col}'."""
    return np.char.add(np.char.add("R", np.asarray(r36).astype(str)), np.char.add("C", np.asarray(c36).astype(str)))


def check_ease_axes(x, y, size=C.EASE_M09, tol=1.0):
    """Assert a file's own x/y axes are EASE cell centres of the given size (within tol metres)."""
    fx = (np.asarray(x, float) - C.EASE_X_MIN) / size - 0.5
    fy = (C.EASE_Y_MAX - np.asarray(y, float)) / size - 0.5
    if np.abs(fx - np.round(fx)).max() * size > tol or np.abs(fy - np.round(fy)).max() * size > tol:
        raise ValueError("SMAP x/y are not EASE-Grid 2.0 cell centres of the expected size")
    return np.round(fy).astype(int), np.round(fx).astype(int)


# ---------------------------------------------------------------- UTM 30 m tiles
def parse_tile(key):
    if not (len(key) == 9 and key[0] == "z" and key[3] == "e" and key[6] == "n"):
        raise ValueError(f"not a tile key: {key!r}")
    return int(key[1:3]), int(key[4:6]), int(key[7:9])


def tile_corner(key):
    zone, i, j = parse_tile(key)
    return zone, C.LATTICE_ORIGIN_M + i * TILE_M, C.LATTICE_ORIGIN_M + (j + 1) * TILE_M


def epsg_of_zone(zone):
    return 32600 + zone


def load_frames(folder):
    p = os.path.join(folder, "frames.json")
    with open(p) as f:
        return {k: v for k, v in json.load(f)["frames"].items() if v}


def frame_transform(F):
    """rasterio Affine of a frame (upper-left corner x0, y0; 30 m pixels)."""
    from rasterio.transform import Affine
    return Affine(C.PIXEL_M, 0.0, float(F["x0"]), 0.0, -C.PIXEL_M, float(F["y0"]))


def frame_key(F):
    return f"{F['crs'].split(':')[-1]}_{int(F['x0'])}_{int(F['y0'])}_{int(F['w'])}_{int(F['h'])}"


def tile_file(folder, tag, tile, year):
    return f"{folder}/{year}/{tile}/{tag}_{tile}_{year}_{C.ROI_LABEL}.nc"


# ---------------------------------------------------------------- regular lat/lon grids (gridMET lab)
def cell_edges(centres):
    c = np.asarray(centres, float)
    step = np.diff(c)
    return np.concatenate([[c[0] - step[0] / 2], (c[:-1] + c[1:]) / 2, [c[-1] + step[-1] / 2]])


def latlon_cell_polygons(lat, lon, rows, cols):
    """shapely boxes of grid cells (rows, cols index lat, lon)."""
    import shapely
    le, oe = cell_edges(lat), cell_edges(lon)
    lo0, lo1 = oe[cols], oe[cols + 1]
    la0, la1 = np.minimum(le[rows], le[rows + 1]), np.maximum(le[rows], le[rows + 1])
    return shapely.box(lo0, la0, lo1, la1)


def haversine_km(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = (np.radians(np.asarray(v, float)) for v in (lon1, lat1, lon2, lat2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 6371.0088 * 2 * np.arcsin(np.sqrt(a))


def days_in_month(year, month):
    import calendar
    return calendar.monthrange(int(year), int(month))[1]


def km(x):
    return math.floor(x / 1000.0)
