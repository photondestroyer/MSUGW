#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Map of the ECOSTRESS 70 m MGRS tiles over the Ogallala (High Plains) aquifer polygon.

Shows every tile the catalogue returns for the polygon, which of them the extraction
notebook keeps (>= 5 % of the tile inside the polygon) and how many kept tiles cover each
part of the aquifer, i.e. where the tiles overlap.

Tile footprints are exact, not the catalogue's lat/lon bounding boxes: each tile is a
1568 x 1568 px, 70 m square in its own UTM zone whose upper-left corner sits on the
Sentinel-2 tiling grid. The corner rule below was read from real ECO_L3T_JET GeoTIFFs
(seven tiles, 2026-10-04) and is asserted against them and against the catalogue boxes.

Needs no login: tile names and boxes come from NASA's public CMR search.
Run:  python "plotting codes/plot_ecostress_tiles.py"
Out:  plots/ecostress_tile_plots/ecostress_tiles_over_ogallala.png  (+ per-tile CSV)
"""
import os
import re
import csv
import json
import math
import urllib.parse
import urllib.request

import numpy as np
import geopandas as gpd
import pyproj
import shapely
from shapely.geometry import Polygon, box
from shapely.geometry.polygon import orient
from shapely.ops import unary_union
import rasterio.features
from rasterio.enums import MergeAlg
from rasterio.transform import from_origin
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHP = os.path.join(REPO, "HPA_polygon", "hp_bound2010.shp")
OUT_DIR = os.path.join(REPO, "plots", "ecostress_tile_plots")
OUT_PNG = os.path.join(OUT_DIR, "ecostress_tiles_over_ogallala.png")
OUT_CSV = os.path.join(OUT_DIR, "ecostress_tiles_over_ogallala.csv")

SHORT_NAME, VERSION = "ECO_L3T_JET", "002"
MIN_TILE_ROI_FRAC = 0.05          # the notebook's pre-download rule
TILE_PX, PIX = 1568, 70.0         # 109.76 km square
EA = "EPSG:5070"                  # equal-area CRS for every area figure
# Catalogue windows searched for tile names (tiles are fixed; a few busy months find them all)
WINDOWS = ["2019-07", "2020-07", "2021-01", "2022-06", "2022-07", "2022-08", "2023-07", "2025-07"]

# Upper-left corners read from real GeoTIFFs: the footprint rule must reproduce them exactly.
KNOWN_UL = {"13SFV": (600000.0, 4000020.0), "13TDF": (399960.0, 4600020.0), "13TDG": (399960.0, 4700040.0),
            "14SKF": (199980.0, 4100040.0), "14SME": (399960.0, 4000020.0), "14TLN": (300000.0, 4800000.0),
            "14TNL": (499980.0, 4600020.0)}

# ---- palette (validated: one-hue ordinal ramp for "how many tiles", one accent for the gap) ----
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]      # 1, 2, 3, 4 tiles
GAP = "#eb6834"                                          # aquifer area with no kept tile


def tile_ul(tile):
    """Upper-left corner (x, y) and EPSG code of an MGRS tile on the Sentinel-2 tiling grid."""
    zone, band, col, row = int(tile[:2]), tile[2], tile[3], tile[4]
    e_west = (["ABCDEFGH", "JKLMNPQR", "STUVWXYZ"][(zone - 1) % 3].index(col) + 1) * 100000.0
    idx = "ABCDEFGHJKLMNPQRSTUV".index(row)
    if zone % 2 == 0:
        idx = (idx - 5) % 20
    # the row letter repeats every 2000 km: pick the repeat that falls in the latitude band
    lat_mid = -80 + 8 * "CDEFGHJKLMNPQRSTUVWX".index(band) + 4
    epsg = 32600 + zone if lat_mid >= 0 else 32700 + zone
    _, n_mid = pyproj.Transformer.from_crs(4326, epsg, always_xy=True).transform(-183 + 6 * zone, lat_mid)
    n_south = idx * 100000.0 + 2e6 * round((n_mid - idx * 100000.0 - 50000.0) / 2e6)
    return math.floor(e_west / 60) * 60.0, math.ceil((n_south + 100000.0) / 60) * 60.0, epsg


def tile_footprint(tile):
    x0, y0, epsg = tile_ul(tile)
    sq = shapely.segmentize(box(x0, y0 - TILE_PX * PIX, x0 + TILE_PX * PIX, y0), 2000.0)
    g = gpd.GeoSeries([sq], crs=epsg)
    return g.to_crs(4326).iloc[0], g.to_crs(EA).iloc[0]


def label_xy(tile, fp4):
    """Tile centre, pushed off a UTM zone boundary where two zones' tile columns overlap."""
    c = fp4.centroid
    east = -180 + 6 * int(tile[:2])                      # east edge of the tile's UTM zone
    dx = -0.34 if c.x > east - 0.95 else (0.34 if c.x < east - 6 + 0.95 else 0.0)
    return c.x + dx, c.y


def cmr_tiles(ring):
    """{tile: [W, S, E, N]} for every tile the catalogue returns over the search ring."""
    poly = ",".join(f"{x:.4f},{y:.4f}" for x, y in ring.exterior.coords)
    tiles = {}
    for ym in WINDOWS:
        y, m = int(ym[:4]), int(ym[5:])
        end = f"{y + (m == 12)}-{m % 12 + 1:02d}-01T00:00:00Z"
        page = 1
        while True:
            q = urllib.parse.urlencode(dict(short_name=SHORT_NAME, version=VERSION, polygon=poly,
                                            temporal=f"{ym}-01T00:00:00Z,{end}", page_size=2000, page_num=page))
            with urllib.request.urlopen("https://cmr.earthdata.nasa.gov/search/granules.json?" + q, timeout=180) as r:
                entries = json.load(r)["feed"]["entry"]
            for e in entries:
                m_ = re.search(r"_(\d{2}[A-Z]{3})_\d{8}T\d{6}_", e.get("producer_granule_id") or e.get("title", ""))
                if m_ and e.get("boxes"):
                    s, w, n, ea = [float(v) for v in e["boxes"][0].split()]
                    tiles.setdefault(m_.group(1), [w, s, ea, n])
            if len(entries) < 2000:
                break
            page += 1
        print(f"  catalogue {ym}: {len(tiles)} distinct tiles so far")
    return tiles


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    g4 = gpd.read_file(SHP).to_crs(4326).geometry
    roi = g4.union_all() if hasattr(g4, "union_all") else g4.unary_union
    roi_ea = gpd.GeoSeries([roi], crs=4326).to_crs(EA).iloc[0]
    ring = orient(Polygon(roi.buffer(0.15).simplify(0.10, preserve_topology=True).exterior), 1.0)   # as the notebook

    for t, ul in KNOWN_UL.items():
        assert tile_ul(t)[:2] == ul, f"footprint rule does not reproduce the real corner of {t}: {tile_ul(t)[:2]} vs {ul}"
    boxes = cmr_tiles(ring)

    rows, worst = [], 0.0
    for tile in sorted(boxes):
        fp4, fpea = tile_footprint(tile)
        worst = max(worst, max(abs(a - b) for a, b in zip(fp4.bounds, boxes[tile])))
        rect_ea = gpd.GeoSeries([shapely.segmentize(box(*boxes[tile]), 0.02)], crs=4326).to_crs(EA).iloc[0]
        rule = rect_ea.intersection(roi_ea).area / rect_ea.area          # what the notebook computes
        rows.append(dict(tile=tile, fp4=fp4, fpea=fpea, in_roi=fpea.intersection(roi_ea),
                         frac_rule=rule, frac_exact=fpea.intersection(roi_ea).area / fpea.area,
                         kept=rule >= MIN_TILE_ROI_FRAC))
    assert worst < 0.02, f"footprints disagree with the catalogue boxes by {worst:.4f} deg"
    print(f"{len(rows)} tiles; footprints match the catalogue boxes within {worst:.4f} deg")

    kept = [r for r in rows if r["kept"]]
    dropped = [r for r in rows if not r["kept"]]
    roi_km2 = roi_ea.area / 1e6
    covered = unary_union([r["in_roi"] for r in kept])
    gap_km2 = (roi_ea.area - covered.area) / 1e6
    stored_km2 = sum(r["in_roi"].area for r in kept) / 1e6
    for r in kept:   # share of each kept tile's in-aquifer area that another kept tile also covers
        others = unary_union([o["in_roi"] for o in kept if o is not r and o["fpea"].intersects(r["fpea"])])
        r["shared"] = r["in_roi"].intersection(others).area / r["in_roi"].area

    # ---- how many kept tiles cover each spot: exact areas on a 250 m equal-area grid ----
    minx, miny, maxx, maxy = roi_ea.bounds
    res = 250.0
    tr = from_origin(minx, maxy, res, res)
    shp = (int(math.ceil((maxy - miny) / res)), int(math.ceil((maxx - minx) / res)))
    count = rasterio.features.rasterize([(r["fpea"], 1) for r in kept], out_shape=shp, transform=tr,
                                        merge_alg=MergeAlg.add, dtype="uint8")
    inside = rasterio.features.rasterize([(roi_ea, 1)], out_shape=shp, transform=tr, dtype="uint8").astype(bool)
    n_in = inside.sum()
    share = {k: float(np.count_nonzero(count[inside] == k)) / n_in for k in range(0, 5)}
    share[4] += float(np.count_nonzero(count[inside] > 4)) / n_in
    multi = 1.0 - share[0] - share[1]
    print(f"aquifer {roi_km2:,.0f} km2 | kept {len(kept)} tiles, dropped {len(dropped)}")
    print("share of aquifer area by number of kept tiles covering it:",
          {k: f"{100 * v:.1f} %" for k, v in share.items()})
    print(f"stored in-aquifer tile area {stored_km2:,.0f} km2 = {stored_km2 / (roi_km2 - gap_km2):.3f} x the covered area"
          f" | not covered by any kept tile: {gap_km2:,.0f} km2 ({100 * gap_km2 / roi_km2:.2f} %)")

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tile", "utm_zone", "kept", "frac_in_aquifer_rule", "frac_in_aquifer_exact",
                    "in_aquifer_km2", "share_of_in_aquifer_area_also_in_another_kept_tile"])
        for r in rows:
            w.writerow([r["tile"], r["tile"][:2], int(r["kept"]), f"{r['frac_rule']:.4f}", f"{r['frac_exact']:.4f}",
                        f"{r['in_roi'].area / 1e6:.1f}", f"{r.get('shared', float('nan')):.4f}"])

    # ---- the map (lon/lat, aspect corrected for the mid latitude) ----
    lon0, lat0, lon1, lat1 = roi.bounds
    pad = 0.75
    ext = (lon0 - pad, lon1 + 1.05, lat0 - pad, lat1 + pad)
    step = 0.004
    tr4 = from_origin(ext[0], ext[3], step, step)
    shp4 = (int(round((ext[3] - ext[2]) / step)), int(round((ext[1] - ext[0]) / step)))
    cnt4 = rasterio.features.rasterize([(r["fp4"], 1) for r in kept], out_shape=shp4, transform=tr4,
                                       merge_alg=MergeAlg.add, dtype="uint8")
    in4 = rasterio.features.rasterize([(roi, 1)], out_shape=shp4, transform=tr4, dtype="uint8").astype(bool)
    img = np.ma.masked_where(~in4, np.minimum(cnt4, 4))

    plt.rcParams.update({"font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 9})
    aspect = 1.0 / math.cos(math.radians((lat0 + lat1) / 2))
    # figure sized to the map so that no blank band is left above or below it
    fig_w, left, right, top_in, bottom_in = 9.0, 0.07, 0.985, 1.25, 0.35
    ax_h = fig_w * (right - left) * (ext[3] - ext[2]) * aspect / (ext[1] - ext[0])
    fig_h = ax_h + top_in + bottom_in
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.imshow(img, extent=ext, origin="upper", cmap=ListedColormap([GAP] + RAMP), vmin=-0.5, vmax=4.5,
              interpolation="nearest", aspect=aspect, zorder=1)

    halo = [pe.withStroke(linewidth=1.6, foreground=SURFACE)]
    for zb in (-102.0, -96.0):                                   # UTM zone boundaries
        ax.plot([zb, zb], [ext[2], ext[3]], color=MUTED, lw=0.8, ls=(0, (1, 3)), zorder=2)
    for name, x in (("UTM zone 13", -104.2), ("UTM zone 14", -99.0), ("zone 15", -95.6)):
        ax.text(x, 1.004, name, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=8, color=MUTED)
    for r in dropped:
        x, y = r["fp4"].exterior.xy
        ax.plot(x, y, color=MUTED, lw=0.6, ls=(0, (4, 3)), zorder=3)
    for r in kept:
        x, y = r["fp4"].exterior.xy
        ax.plot(x, y, color=INK2, lw=0.55, zorder=4)
    for geom in getattr(roi, "geoms", [roi]):
        x, y = geom.exterior.xy
        ax.plot(x, y, color=INK, lw=1.1, zorder=5)
    for r in rows:
        lx, ly = label_xy(r["tile"], r["fp4"])
        ax.text(lx, ly, r["tile"], ha="center", va="center", fontsize=5.6, clip_on=True,
                color=INK2 if r["kept"] else MUTED, path_effects=halo, zorder=6)

    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    xt = np.arange(math.ceil(ext[0] / 2) * 2, ext[1], 2)
    yt = np.arange(math.ceil(ext[2] / 2) * 2, ext[3], 2)
    ax.set_xticks(xt, [f"{abs(v):.0f}°W" for v in xt])
    ax.set_yticks(yt, [f"{v:.0f}°N" for v in yt])
    ax.tick_params(colors=MUTED, length=0, labelsize=8)
    ax.grid(color=GRID, lw=0.6, zorder=0)
    for s in ax.spines.values():
        s.set_visible(False)

    handles = [Patch(fc=RAMP[k - 1], ec="none", label=f"{k} tile{'s' if k > 1 else ''}   {100 * share[k]:.1f} % of the aquifer")
               for k in (1, 2, 3, 4) if share[k] > 0]
    handles.append(Patch(fc=GAP, ec="none", label=f"no kept tile   {100 * share[0]:.2f} % ({gap_km2:,.0f} km²)"))
    handles += [Line2D([], [], color=INK2, lw=0.9, label=f"tile kept, {len(kept)} (≥ 5 % inside the aquifer)"),
                Line2D([], [], color=MUTED, lw=0.9, ls=(0, (4, 3)), label=f"tile dropped, {len(dropped)} (< 5 %)"),
                Line2D([], [], color=INK, lw=1.3, label="aquifer boundary"),
                Line2D([], [], color=MUTED, lw=0.9, ls=(0, (1, 3)), label="UTM zone boundary")]
    leg = ax.legend(handles=handles, loc="lower right", frameon=True, framealpha=1.0, facecolor=SURFACE,
                    edgecolor=GRID, fontsize=8, title="Kept tiles covering each spot", title_fontsize=8.5,
                    borderpad=0.9, labelspacing=0.55)
    leg.get_title().set_color(INK)
    for t in leg.get_texts():
        t.set_color(INK2)
    leg._legend_box.align = "left"

    fig.text(0.06, 1 - 0.22 / fig_h, "ECOSTRESS 70 m tiles over the Ogallala aquifer", ha="left", va="top",
             fontsize=14, fontweight="semibold", color=INK)
    fig.text(0.06, 1 - 0.55 / fig_h,
             f"{100 * multi:.0f} % of the aquifer lies in two or more tiles. Tiles are 109.76 km squares on a "
             f"100 km lattice,\nso neighbours share a 9.76 km strip; tiles of two UTM zones overlap more along 102°W.",
             ha="left", va="top", fontsize=9, color=INK2)
    fig.subplots_adjust(left=left, right=right, bottom=bottom_in / fig_h, top=1 - top_in / fig_h)
    fig.savefig(OUT_PNG, dpi=200, facecolor=SURFACE)
    print("saved", OUT_PNG)
    print("saved", OUT_CSV)


if __name__ == "__main__":
    main()
