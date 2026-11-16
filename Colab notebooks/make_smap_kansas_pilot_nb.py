"""Build SMAP_L3E_Kansas_pilot.ipynb from SMAP_VOD.ipynb (the original is only read).

Why: the existing SMAP L3E files hold data only inside the aquifer polygon, and the corrected one ends in
June 2022. Reference footprints (pilot Task 1.2) lie up to 150 km from Kansas counties, much of it off the
aquifer, so the pilot needs an UNMASKED extraction over a box, through 2024.

Changes (each marked PILOT in the code): the region is a box, the study counties plus 150 km and half a
footprint (lon -104.04 to -95.22, lat 35.49 to 41.50), with no polygon mask; output folder and file prefix.
Everything else (granule choice, flags, values as NASA stores them, resume, checks) is unchanged.
Regenerate after any change here:   python "Colab notebooks/make_smap_kansas_pilot_nb.py"
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "SMAP_VOD.ipynb"
DST = HERE / "SMAP_L3E_Kansas_pilot.ipynb"
BOX = (-104.04, 35.49, -95.22, 41.50)            # lon_min, lat_min, lon_max, lat_max

nb = json.load(open(SRC, encoding="utf-8"))
cells = nb["cells"]


def sub_any(old, new):
    hits = [i for i, c in enumerate(cells) if old in "".join(c["source"])]
    if len(hits) != 1:
        sys.exit(f"expected one cell with {old[:60]!r}, found {len(hits)}")
    i = hits[0]
    s = "".join(cells[i]["source"])
    if s.count(old) != 1:
        sys.exit(f"cell {i}: {old[:60]!r} occurs {s.count(old)} times")
    cells[i]["source"] = s.replace(old, new).splitlines(keepends=True)


sub_any('DRIVE_OUT_DIR = "/content/drive/MyDrive/SMAP_L3E_Ogallala"',
        'DRIVE_OUT_DIR = "/content/drive/MyDrive/MSUGWB/SMAP_L3E_Kansas_pilot"   # PILOT: own folder')
sub_any('OUTPUT_PREFIX = "SPL3SMP_E_Ogallala_reprocessed_"',
        'OUTPUT_PREFIX = "SPL3SMP_E_Kansas_pilot_"   # PILOT')
sub_any('''def load_roi(shp_path, deg=0.1):
    gdf = gpd.read_file(shp_path)
    if gdf.crs is None or gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(epsg=4326)
    roi = gdf.geometry.union_all() if hasattr(gdf.geometry, "union_all") else gdf.geometry.unary_union
    return roi.simplify(deg, preserve_topology=True)''',
        f'''PILOT_BOX = {BOX}   # PILOT: lon_min, lat_min, lon_max, lat_max - Kansas study counties + 150 km + 18 km


def load_roi(shp_path, deg=0.1):
    """PILOT: a box with NO polygon mask: every 9 km cell whose centre lies in the box keeps NASA's values, so
    reference footprints off the aquifer are observed too (pilot Task 1.2). shp_path is not used."""
    return shapely.box(*PILOT_BOX)''')
sub_any('m_var.long_name = "cell centre inside the High Plains aquifer polygon (hp_bound2010, simplified 0.1 deg)"',
        'm_var.long_name = "cell centre inside the pilot box (Kansas study counties + 150 km); no polygon mask"')
sub_any('ds.title = "SMAP Enhanced L3 radiometer soil moisture and VOD (SPL3SMP_E v006), High Plains aquifer subset"',
        'ds.title = "SMAP Enhanced L3 radiometer soil moisture and VOD (SPL3SMP_E v006), Kansas pilot box (unmasked)"')
sub_any('"the aquifer polygon set to _FillValue.")', '"the pilot box set to _FillValue (no polygon mask).")')

HEAD = f"""# PILOT copy: SMAP L3E (SPL3SMP_E) over a box around the Kansas study counties, no polygon mask

A copy of `SMAP_VOD.ipynb` for the groundwater-ET pilot (Tasks 1.2 and 1.3). Built by
`make_smap_kansas_pilot_nb.py`; do not edit by hand.

| Setting | Value | Why |
|---|---|---|
| Region | box lon {BOX[0]} to {BOX[2]}, lat {BOX[1]} to {BOX[3]}; **no polygon mask** | reference footprints lie up to 150 km from the counties, much of it off the aquifer |
| Output | `MyDrive/MSUGWB/SMAP_L3E_Kansas_pilot/SPL3SMP_E_Kansas_pilot_<date>.nc` | never mixed with the aquifer files |
| Period | mission start to today | the pilot needs 2015-2024 (climatology) and 2016-2024 (panel) |

Run: *Run all* (Earthdata login is asked once). Rerun after any stop: it resumes after the last day written.
The pilot reads the newest `SPL3SMP_E_Kansas_pilot_*.nc` in that folder (checkpoint copies are ignored).
"""
cells.insert(0, {"cell_type": "markdown", "metadata": {}, "source": HEAD.splitlines(keepends=True)})
for c in cells:
    if c["cell_type"] == "code":
        c["outputs"], c["execution_count"] = [], None
json.dump(nb, open(DST, "w", encoding="utf-8", newline="\n"), ensure_ascii=False, indent=1)
print("written", DST, len(cells), "cells")
