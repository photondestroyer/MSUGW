"""S0.1 (part): the 480 m block lattice nested in the 30 m tile lattice, and the block table.

Tile key z{zone}e{i:02d}n{j:02d} (Landsat 30 m notebook): EPSG:326{zone}, tile corner
    left = 15 + i * 153600 m,   top = 15 + (j + 1) * 153600 m      (lattice origin 15 m mod 30)
A tile's FRAME (frames.json) is the tile cropped to the bounding box of its aquifer pixels, so a file's
pixel (0, 0) is NOT the tile corner: frame offset r0 = (top - y0) / 30, c0 = (x0 - left) / 30.
Blocks are 16 x 16 pixels aligned to the TILE corner (5120 px = 320 blocks), so they line up across
years, products and frames. Block (br, bc) of a tile covers tile rows br*16 .. br*16+15.

Blocks along the 13/14 (and 14/15) zone boundary exist in both zones. `primary` marks the copy that is
counted: the copy in the zone whose 6-degree band holds the block centre, unless that zone has no valid
block there among the selected tiles, in which case the other copy is primary.
"""
import json
import math
import os
from dataclasses import dataclass, asdict

import numpy as np

from .. import config as C
from ..core import paths, ncio
from ..core.errors import SourceError, DataIntegrityError

TILE_M = C.TILE_PX * C.PIXEL_M
NB_TILE = C.TILE_PX // C.BLOCK_PX


def parse_tile(key):
    if not (len(key) == 9 and key[0] == "z" and key[3] == "e" and key[6] == "n"):
        raise ValueError(f"not a tile key: {key!r}")
    return int(key[1:3]), int(key[4:6]), int(key[7:9])


def tile_key(zone, i, j):
    return f"z{zone:02d}e{i:02d}n{j:02d}"


def tile_corner(key):
    zone, i, j = parse_tile(key)
    return zone, C.LATTICE_ORIGIN_M + i * TILE_M, C.LATTICE_ORIGIN_M + (j + 1) * TILE_M


def epsg_of_zone(zone):
    return 32600 + zone


def band_zone(lon):
    return np.floor((np.asarray(lon) + 180.0) / 6.0).astype(int) + 1


@dataclass
class TileGeom:
    tile: str
    zone: int
    epsg: int
    left: float
    top: float
    x0: float
    y0: float
    w: int
    h: int
    r0: int          # frame offset in tile pixels
    c0: int
    br0: int         # block window (tile block coordinates)
    bc0: int
    nby: int
    nbx: int

    # pixel window of the block window, in FRAME pixel coordinates (may extend beyond the frame)
    def frame_rows(self):
        a = self.br0 * C.BLOCK_PX - self.r0
        return a, a + self.nby * C.BLOCK_PX

    def frame_cols(self):
        a = self.bc0 * C.BLOCK_PX - self.c0
        return a, a + self.nbx * C.BLOCK_PX

    def pixel_centres(self, rows=None):
        """UTM x (cols) and y (rows) of the pixel centres of the block window (or a row range of it)."""
        r_a, r_b = rows if rows else (0, self.nby * C.BLOCK_PX)
        tr = self.br0 * C.BLOCK_PX + np.arange(r_a, r_b)
        tc = self.bc0 * C.BLOCK_PX + np.arange(self.nbx * C.BLOCK_PX)
        return self.left + (tc + 0.5) * C.PIXEL_M, self.top - (tr + 0.5) * C.PIXEL_M

    def block_centres(self):
        br = self.br0 + np.arange(self.nby)
        bc = self.bc0 + np.arange(self.nbx)
        half = C.BLOCK_PX * C.PIXEL_M / 2
        return self.left + bc * C.BLOCK_PX * C.PIXEL_M + half, self.top - br * C.BLOCK_PX * C.PIXEL_M - half

    def to_dict(self):
        return asdict(self)


def geom_from_frame(key, F):
    zone, left, top = tile_corner(key)
    if int(str(F["crs"]).split(":")[-1]) != epsg_of_zone(zone):
        raise SourceError(f"{key}: frame CRS {F['crs']} is not EPSG:{epsg_of_zone(zone)}")
    c0f = (F["x0"] - left) / C.PIXEL_M
    r0f = (top - F["y0"]) / C.PIXEL_M
    if abs(c0f - round(c0f)) > 1e-6 or abs(r0f - round(r0f)) > 1e-6:
        raise SourceError(f"{key}: frame origin is not on the tile's 30 m lattice ({F['x0']}, {F['y0']})")
    r0, c0 = int(round(r0f)), int(round(c0f))
    if not (0 <= r0 and 0 <= c0 and r0 + F["h"] <= C.TILE_PX and c0 + F["w"] <= C.TILE_PX):
        raise SourceError(f"{key}: frame lies outside its tile")
    br0, bc0 = r0 // C.BLOCK_PX, c0 // C.BLOCK_PX
    nby = -(-(r0 + F["h"]) // C.BLOCK_PX) - br0
    nbx = -(-(c0 + F["w"]) // C.BLOCK_PX) - bc0
    return TileGeom(key, zone, epsg_of_zone(zone), left, top, float(F["x0"]), float(F["y0"]), int(F["w"]),
                    int(F["h"]), r0, c0, br0, bc0, nby, nbx)


def test_window(g):
    """Restrict a geometry to the test square (TEST_WINDOW_PX), snapped to blocks and kept inside the frame."""
    nb = C.TEST_WINDOW_PX // C.BLOCK_PX
    if C.TEST_WINDOW_CENTER:
        cr, cc = C.TEST_WINDOW_CENTER
        cr, cc = (g.r0 + cr) // C.BLOCK_PX, (g.c0 + cc) // C.BLOCK_PX
    else:
        cr, cc = g.br0 + g.nby // 2, g.bc0 + g.nbx // 2
    br0 = int(min(max(cr - nb // 2, g.br0), max(g.br0 + g.nby - nb, g.br0)))
    bc0 = int(min(max(cc - nb // 2, g.bc0), max(g.bc0 + g.nbx - nb, g.bc0)))
    g.br0, g.bc0, g.nby, g.nbx = br0, bc0, min(nb, g.nby), min(nb, g.nbx)
    return g


def load_frames(key="landsat"):
    p = paths.frames_path(key)
    if not os.path.exists(p):
        raise SourceError(f"{key}: frames.json not found at {p}")
    with open(p) as f:
        return json.load(f)["frames"]


def tiles_for_points(frames, points):
    from pyproj import Transformer
    out = []
    for name, (lon, lat) in points.items():
        z = int(band_zone(lon))
        x, y = Transformer.from_crs(4326, epsg_of_zone(z), always_xy=True).transform(lon, lat)
        k = tile_key(z, int((x - C.LATTICE_ORIGIN_M) // TILE_M), int((y - C.LATTICE_ORIGIN_M) // TILE_M))
        if frames.get(k):
            out.append(k)
    return sorted(set(out))


def selected_tiles(frames):
    if C.RUN_MODE == "test":
        missing = [t for t in C.TEST_TILES if not frames.get(t)]
        if missing:
            raise SourceError(f"test tiles without a frame in frames.json: {missing}")
        return list(C.TEST_TILES)
    if C.RUN_MODE == "focus":
        return tiles_for_points(frames, C.FOCUS_POINTS)
    return sorted(k for k, F in frames.items() if F)


def geometries(frames=None):
    frames = frames or load_frames()
    out = {}
    for k in selected_tiles(frames):
        g = geom_from_frame(k, frames[k])
        out[k] = test_window(g) if C.RUN_MODE == "test" else g
    return out


def check_frames_agree(frames, tiles):
    """OpenET / CDL / NWI products must use the same frames as Landsat for the tiles used."""
    report = {}
    for key in ("openet", "cdl", "nwi"):
        p = paths.frames_path(key)
        if not os.path.exists(p):
            report[key] = "no frames.json"
            continue
        other = json.load(open(p))["frames"]
        bad = [t for t in tiles if other.get(t) and any(
            abs(float(other[t][c]) - float(frames[t][c])) > 1e-6 for c in ("x0", "y0", "w", "h"))]
        if bad:
            raise DataIntegrityError(f"{key} frames differ from Landsat frames for tiles {bad}")
        report[key] = "identical"
    return report


def blockify(a, nby, nbx):
    """(nby*16, nbx*16) -> (nby, nbx, 256)."""
    b = C.BLOCK_PX
    return a.reshape(nby, b, nbx, b).transpose(0, 2, 1, 3).reshape(nby, nbx, b * b)


def first_source_file(key, tile, years):
    for y in years:
        p = paths.tile_path(key, tile, y)
        if os.path.exists(p):
            return p
    return None


def lab_cell(lat, lon, lab_lat, lab_lon):
    """Index of the lab gridMET cell containing each point (-1 outside the grid)."""
    dlat, dlon = lab_lat[0] - lab_lat[1], lab_lon[1] - lab_lon[0]
    iy = np.rint((lab_lat[0] - lat) / dlat).astype(int)
    ix = np.rint((lon - lab_lon[0]) / dlon).astype(int)
    ok = (iy >= 0) & (iy < lab_lat.size) & (ix >= 0) & (ix < lab_lon.size)
    iy, ix = np.where(ok, iy, -1), np.where(ok, ix, -1)
    return iy, ix


def block_id(zone, i, j, br, bc):
    return ((((np.int64(zone) * 100 + i) * 100 + j) * 1000 + br) * 1000 + bc).astype("int64")


def run(ctx):
    """Stage `tiles`: block table of the selected tiles."""
    from pyproj import Transformer
    frames = load_frames()
    geos = geometries(frames)
    agree = check_frames_agree(frames, list(geos))
    gm = paths.dataset_path("gridmet")
    with ncio.open_ro(gm) as ds:
        lab_lat = np.asarray(ds["lat"][:], "float64")
        lab_lon = np.asarray(ds["lon"][:], "float64")
        roi_frac = np.asarray(ds["roi_fraction"][:], "float32")
    cols = {k: [] for k in ("block_id", "tile_idx", "zone", "br", "bc", "x", "y", "lon", "lat", "roi_px",
                            "roi_frac", "lab_iy", "lab_ix", "lab_roi_fraction", "band_zone")}
    tiles = sorted(geos)
    inputs = {"gridmet": gm}
    for ti, t in enumerate(tiles):
        g = geos[t]
        src = first_source_file("landsat", t, C.agg_years())
        if src is None:
            raise SourceError(f"{t}: no Landsat file for any of {C.agg_years()} (needed for roi_mask)")
        inputs[f"landsat_roi_{t}"] = src
        with ncio.open_ro(src) as ds:
            roi = ncio.read_window(ds["roi_mask"], g.frame_rows(), g.frame_cols(), 0)
        npx = blockify((roi == 1).astype("int16"), g.nby, g.nbx).sum(-1)
        xs, ys = g.block_centres()
        X, Y = np.meshgrid(xs, ys)
        lon, lat = Transformer.from_crs(g.epsg, 4326, always_xy=True).transform(X, Y)
        iy, ix = lab_cell(lat, lon, lab_lat, lab_lon)
        BR, BC = np.meshgrid(g.br0 + np.arange(g.nby), g.bc0 + np.arange(g.nbx), indexing="ij")
        z, i, j = parse_tile(t)
        vals = dict(block_id=block_id(z, i, j, BR, BC), tile_idx=np.full(BR.shape, ti), zone=np.full(BR.shape, z),
                    br=BR, bc=BC, x=X, y=Y, lon=lon, lat=lat, roi_px=npx,
                    roi_frac=npx / float(C.BLOCK_PX ** 2), lab_iy=iy, lab_ix=ix,
                    lab_roi_fraction=np.where(iy >= 0, roi_frac[np.maximum(iy, 0), np.maximum(ix, 0)], np.nan),
                    band_zone=band_zone(lon))
        for k in cols:
            cols[k].append(np.asarray(vals[k]).ravel())
        ctx.unit(t, "done", blocks=int(npx.size), inside=int((npx / 256.0 >= C.BLOCK_MIN_ROI_FRACTION).sum()))
    T = {k: np.concatenate(v) for k, v in cols.items()}
    T["valid"] = T["roi_frac"] >= C.BLOCK_MIN_ROI_FRACTION
    T["primary"] = resolve_primary(T, tiles, geos)
    np.savez_compressed(f"{ctx.out}/blocks.npz", tiles=np.array(tiles), **T)
    with open(f"{ctx.out}/tiles.json", "w") as f:
        json.dump({t: geos[t].to_dict() for t in tiles}, f, indent=1)
    summary = dict(tiles=len(tiles), blocks=int(T["valid"].size), valid_blocks=int(T["valid"].sum()),
                   primary_valid=int((T["valid"] & T["primary"]).sum()), frames=agree,
                   area_km2=float(T["roi_px"].sum() * C.PIXEL_M ** 2 / 1e6))
    ctx.complete(summary)
    ctx.log.info("blocks: %s", summary)


def resolve_primary(T, tiles, geos):
    """Primary copy of blocks that exist in two zones (see module docstring)."""
    from pyproj import Transformer
    primary = np.ones(T["zone"].size, bool)
    by_zone_index = {}
    for ti, t in enumerate(tiles):
        z = parse_tile(t)[0]
        sel = (T["tile_idx"] == ti) & T["valid"]
        for br, bc in zip(T["br"][sel], T["bc"][sel]):
            by_zone_index.setdefault(z, set()).add((t, int(br), int(bc)))
    off = np.where(T["valid"] & (T["band_zone"] != T["zone"]))[0]
    for bz in np.unique(T["band_zone"][off]):
        idx = off[T["band_zone"][off] == bz]
        if bz not in by_zone_index:
            continue                          # no tile of that zone selected: this copy stays primary
        x, y = Transformer.from_crs(4326, epsg_of_zone(int(bz)), always_xy=True).transform(T["lon"][idx], T["lat"][idx])
        ti_ = ((x - C.LATTICE_ORIGIN_M) // TILE_M).astype(int)
        tj_ = ((y - C.LATTICE_ORIGIN_M) // TILE_M).astype(int)
        for k, n in enumerate(idx):
            t = tile_key(int(bz), int(ti_[k]), int(tj_[k]))
            _, left, top = tile_corner(t)
            bc = int((x[k] - left) // (C.BLOCK_PX * C.PIXEL_M))
            br = int((top - y[k]) // (C.BLOCK_PX * C.PIXEL_M))
            if (t, br, bc) in by_zone_index[bz]:
                primary[n] = False
    return primary


def load_table(stage_dir):
    with np.load(f"{stage_dir}/blocks.npz", allow_pickle=False) as z:
        T = {k: z[k] for k in z.files}
    with open(f"{stage_dir}/tiles.json") as f:
        G = {t: TileGeom(**d) for t, d in json.load(f).items()}
    return T, G


def tile_slice(T, tile_idx):
    """Rows of the block table that belong to one tile (they are stored contiguously, row-major)."""
    idx = np.where(T["tile_idx"] == tile_idx)[0]
    return slice(int(idx[0]), int(idx[-1]) + 1)
