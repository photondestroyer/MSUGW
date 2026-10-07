"""Land cover, cropland and irrigation with the published legends.

What was wrong before (findings C1-C4):
  * the GFSAD file carries an invented legend; classes 3-5 are rainfed cropland,
    not water / urban / natural vegetation;
  * the irrigation map's class 2 is HIGH irrigation; XGBoost counted only class 1
    as irrigated and the LSTM counted every class above 0;
  * a 4 km cell was represented by one 500 m or 1 km pixel;
  * class codes were fed to the LSTM as if they were quantities.
Here every layer becomes an area fraction of the analysis cell.
"""
import numpy as np

from . import config as C
from . import grid as G
from . import sources as S

# USGS/GFSAD1000_V1, band `landcover` (Earth Engine data catalogue). Nominal year 2010.
GFSAD_LEGEND = {0: "non-cropland", 1: "cropland, irrigation major", 2: "cropland, irrigation minor",
                3: "cropland, rainfed", 4: "cropland, rainfed, minor fragments",
                5: "cropland, rainfed, very minor fragments"}
# Classes 4 and 5 are land with minor cropland fragments. In the aquifer they follow
# MCD12Q1 grassland, not cropland, so they are kept apart from rainfed cropland.
GFSAD_IRRIGATED, GFSAD_RAINFED, GFSAD_FRAGMENTS = (1, 2), (3,), (4, 5)

# Nagaraj et al. (2021), Adv. Water Resour. 152, 103910; band `classification`, ~9 km, 2001-2015.
GIR_LEGEND = {0: "no or very little irrigation",
              1: "low-to-medium irrigation (<= 2000 ha per 86 km2)",
              2: "high irrigation (> 2000 ha per 86 km2)"}

# MCD12Q1 LC_Type1 (IGBP).
IGBP_LEGEND = {1: "evergreen needleleaf forest", 2: "evergreen broadleaf forest",
               3: "deciduous needleleaf forest", 4: "deciduous broadleaf forest", 5: "mixed forest",
               6: "closed shrubland", 7: "open shrubland", 8: "woody savanna", 9: "savanna",
               10: "grassland", 11: "permanent wetland", 12: "cropland", 13: "urban",
               14: "cropland/natural mosaic", 15: "snow and ice", 16: "barren", 17: "water"}
IGBP_NATURAL = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16)

MGMT_IRRIGATED, MGMT_RAINFED_CROP, MGMT_NATURAL, MGMT_MIXED = 0, 1, 2, 3
MGMT_LABELS = {MGMT_IRRIGATED: "irrigated", MGMT_RAINFED_CROP: "rainfed_crop",
               MGMT_NATURAL: "natural", MGMT_MIXED: "mixed_uncertain"}


def gfsad_fractions(grid):
    """Area share of each GFSAD class per cell, plus the pixel count."""
    r = S.Raster(C.SRC["gfsad"], "landcover")
    try:
        ysl, xsl = r.window(grid)
        cls = r.read_classes(0, ysl, xsl)
        cls[(cls < 0) | (cls > 5)] = -1
        cid = G.cell_ids(r.lat[ysl], r.lon[xsl], grid)
        frac, n = G.class_fractions(cls, cid, grid.n, [(k,) for k in range(6)])
    finally:
        r.close()
    return dict(f_noncrop=frac[0], f_irr_major=frac[1], f_irr_minor=frac[2],
                f_irrigated=frac[1] + frac[2], f_rainfed_crop=frac[3],
                f_crop_fragments=frac[4] + frac[5],
                gfsad_n=n)


def igbp_fractions(grid, years=C.LC_YEARS):
    """Mean area share of each PFT group over `years` (pre-window only: no look-ahead),
    and the share of natural cover."""
    r = S.Raster(C.SRC["mcd12q1"], "LC_Type1")
    try:
        ysl, xsl = r.window(grid)
        cid = G.cell_ids(r.lat[ysl], r.lon[xsl], grid)
        yr = r.days.astype("datetime64[Y]").astype(int) + 1970
        sel = np.where((yr >= years[0]) & (yr <= years[1]))[0]
        if sel.size == 0:
            raise S.SourceError(f"MCD12Q1 has no year in {years}")
        groups = list(C.PFT_GROUPS.items()) + [("f_natural", IGBP_NATURAL)]
        acc = np.zeros((len(groups), grid.n), "float64")
        npx = None
        for i in sel:
            frac, npx = G.class_fractions(r.read_classes(int(i), ysl, xsl), cid, grid.n,
                                          [g for _, g in groups])
            acc += frac
        out = {name: (acc[k] / sel.size).astype("float32") for k, (name, _) in enumerate(groups)}
        out["igbp_n"] = npx
        out["igbp_years"] = yr[sel]
    finally:
        r.close()
    return out


def hsg_fractions(grid):
    """Area share of hydrologic soil groups A-D. Dual codes (14, 24, 34: drained /
    undrained complexes) are assigned to their drained group, not dropped."""
    r = S.Raster(C.SRC["hsg"], "b1")
    try:
        ysl, xsl = r.window(grid)
        cls = r.read_classes(0, ysl, xsl)
        cls[cls == 255] = -1
        cid = G.cell_ids(r.lat[ysl], r.lon[xsl], grid)
        frac, n = G.class_fractions(cls, cid, grid.n, [(1, 14), (2, 24), (3, 34), (4,)])
    finally:
        r.close()
    return dict(hsg_A=frac[0], hsg_B=frac[1], hsg_C=frac[2], hsg_D=frac[3], hsg_n=n)


def irrigation_frequency(grid):
    """Share of years 2001-2015 in which the cell's parent ~9 km pixel was class 1
    (low-to-medium) and class 2 (high). The map is coarser than the grid, so each
    cell takes its nearest pixel; `gir_dist` records how far away that is."""
    r = S.Raster(C.SRC["gir"], "classification")
    try:
        lat2, lon2 = np.meshgrid(r.lat, r.lon, indexing="ij")
        idx, dist = G.nearest_source(lat2, lon2, grid.cell_lat, grid.cell_lon)
        low = np.zeros(grid.n, "float64")
        high = np.zeros(grid.n, "float64")
        nyr = np.zeros(grid.n, "float64")
        for i in range(r.days.size):
            cls = r.read_classes(i, slice(None), slice(None)).ravel()[idx]
            ok = np.isin(cls, (0, 1, 2))                   # fill is missing, not "not irrigated"
            low += ok & (cls == 1)
            high += ok & (cls == 2)
            nyr += ok
    finally:
        r.close()
    with np.errstate(invalid="ignore", divide="ignore"):
        return dict(gir_low_freq=np.where(nyr > 0, low / nyr, np.nan).astype("float32"),
                    gir_high_freq=np.where(nyr > 0, high / nyr, np.nan).astype("float32"),
                    gir_years=nyr.astype("int16"), gir_dist=dist)


def classify_management(f_irrigated, f_rainfed_crop, f_natural, gir_high_freq):
    """Management class of each cell from its area fractions.

    Inputs are arrays of shape (N,):
      f_irrigated     GFSAD share of the cell that is irrigated cropland (classes 1-2)
      f_rainfed_crop  GFSAD share that is rainfed cropland (class 3; the fragment classes
                      4-5 are stored separately as f_crop_fragments and are not passed in)
      f_natural       MCD12Q1 share of natural cover (forest, shrub, savanna, grass, wetland, barren)
      gir_high_freq   share of years 2001-2015 the parent ~9 km pixel was "high irrigation"

    Returns an int8 array (N,) of MGMT_IRRIGATED, MGMT_RAINFED_CROP, MGMT_NATURAL
    or MGMT_MIXED. Research Plan 6.2 asks for rainfed, irrigated and
    uncertain-management locations to be analysed separately, so a cell that is
    not clearly one thing belongs in MGMT_MIXED.
    """
    # TODO(human): decide the thresholds and write the classification.
    raise NotImplementedError("classify_management() is waiting for its thresholds")
