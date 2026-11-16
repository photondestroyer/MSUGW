"""Build OpenET_v2_1_Kansas_pilot.ipynb from OpenET_Monthly_30m_Ogallala.ipynb (the original is only read).

The pilot notebook is generated, not hand-edited: change the settings here and rerun
    python "Colab notebooks/make_openet_kansas_pilot_nb.py"
Each edit must match the original exactly once, so the script stops if the original notebook changes under it.
Pilot settings: OpenET v2.1, the 11 UTM zone-14 tiles holding the Kansas study area, whole years 2016-2024,
rehearsal June-July 2023 (low deficit), output folder OPENET_v2_1_monthly_30m_tiles_KS_pilot.
"""
import json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "OpenET_Monthly_30m_Ogallala.ipynb"
DST = HERE / "OpenET_v2_1_Kansas_pilot.ipynb"
nb = json.load(open(SRC, encoding="utf-8"))
cells = nb["cells"]


def sub(i, old, new):
    s = "".join(cells[i]["source"])
    n = s.count(old)
    if n != 1:
        sys.exit(f"cell {i}: expected 1 match, found {n}: {old[:70]!r}")
    cells[i]["source"] = s.replace(old, new).splitlines(keepends=True)


TILES = ["z14e01n26", "z14e01n27", "z14e01n28", "z14e02n26", "z14e02n27", "z14e02n28",
         "z14e03n26", "z14e03n27", "z14e03n28", "z14e04n26", "z14e04n27"]

# ---- cell 2: settings ----
sub(2, 'TEST_YEAR     = 2022       # year of the rehearsal',
    'TEST_YEAR     = 2023       # PILOT: a low-deficit year. Kansas GMD May-Sep gridMET ETo - P = 436 mm (2nd lowest of\n'
    '                           # 2016-2024; June 42, July 71 mm). 2022 (783 mm) is the driest season and is not used,\n'
    '                           # as the pilot requires (blinding rule)')
sub(2, 'START_MONTH = "1999-10"    # first month, written "YYYY-MM" (OpenET\'s monthly record starts in October 1999)',
    'START_MONTH = "2016-01"    # PILOT: first month, written "YYYY-MM" (whole years 2016-2024)')
sub(2, 'STOP_MONTH  = None         # last month, written "YYYY-MM" (it is included); None -> the latest month OpenET has.',
    'STOP_MONTH  = "2024-12"    # PILOT: last month, written "YYYY-MM" (it is included); None -> the latest month OpenET has.')
sub(2, '                           # to the yearly files that are already on Drive and never fetches a month twice\n',
    '                           # to the yearly files that are already on Drive and never fetches a month twice\n'
    'ONLY_CALENDAR_MONTHS = None   # PILOT: None -> every month of the range (whole years). A list such as\n'
    '                           # [5, 6, 7, 8, 9] keeps only those calendar months and never fetches the others\n')
sub(2, 'ONLY_TILES = None          # None -> every tile; or e.g. ["z14e02n29", "z14e02n28"] (to split the work over sessions)',
    'ONLY_TILES = ' + json.dumps(TILES) + '\n'
    '                           # PILOT: the 11 UTM zone-14 tiles that hold all of the Kansas study area (GMD 1-5 plus the\n'
    '                           # KGS High Plains aquifer extent in Kansas). Zone-13 tiles overlap Kansas but add no ground,\n'
    "                           # except a ~2 km strip at the Colorado border (west of about -102.027 deg) that OpenET's\n"
    '                           # zone-14 images do not reach; cells there are excluded with that reason. To recover it, add\n'
    '                           # "z13e05n26", "z13e05n27", "z13e05n28" (roughly doubles the work on that ground)')
sub(2, 'OPENET_VERSION = "v2_1"    # which release to fetch.',
    'OPENET_VERSION = "v2_1"    # PILOT: keep v2_1 (one ET release, pilot rule). Which release to fetch.')

# ---- cell 3: own output folder ----
sub(3, 'DRIVE_OUT   = f"{DRIVE_ROOT}/{PRODUCT_TAG}" + ("/_test" if TEST_MODE else "")   # the rehearsal never touches real outputs',
    'DRIVE_OUT   = f"{DRIVE_ROOT}/{PRODUCT_TAG}_KS_pilot" + ("/_test" if TEST_MODE else "")   # PILOT: its own folder, never\n'
    '                                                    # mixed with a full-aquifer run; the rehearsal never touches real outputs')

# ---- cell 5: calendar-month filter ----
sub(5, '''def range_months():
    """Every (year, month) from START_MONTH to STOP_MONTH for which OpenET has data, in order."""
    return sorted(k for k in AVAILABLE if k >= FIRST_MONTH and (LAST_MONTH is None or k <= LAST_MONTH))''',
    '''if ONLY_CALENDAR_MONTHS is not None and (not ONLY_CALENDAR_MONTHS
                                         or any(m not in range(1, 13) for m in ONLY_CALENDAR_MONTHS)):
    raise ValueError(f"ONLY_CALENDAR_MONTHS = {ONLY_CALENDAR_MONTHS!r}: give months 1-12, or None.")


def range_months():
    """Every (year, month) from START_MONTH to STOP_MONTH for which OpenET has data, in order,
    kept to ONLY_CALENDAR_MONTHS when that is set (None in the pilot: whole years)."""
    return sorted(k for k in AVAILABLE if k >= FIRST_MONTH and (LAST_MONTH is None or k <= LAST_MONTH)
                  and (ONLY_CALENDAR_MONTHS is None or k[1] in ONLY_CALENDAR_MONTHS))''')

# ---- cell 9: projection only over the chosen tiles ----
sub(9, '''            n_all = sum(len(v) for P in ALL_TILES for k, v in P["months"].items() if k in wanted)''',
    '''            proj_tiles = [P for P in ALL_TILES if not ONLY_TILES or P["key"] in ONLY_TILES]   # PILOT: chosen tiles only
            n_all = sum(len(v) for P in proj_tiles for k, v in P["months"].items() if k in wanted)''')
sub(9, '''                all_px = sum(v["n_roi"] for v in json.load(f)["frames"].values() if v)''',
    '''                all_px = sum(v["n_roi"] for k, v in json.load(f)["frames"].items()
                             if v and (not ONLY_TILES or k in ONLY_TILES))''')
sub(9, '''            log.info("     %d image pieces on %d tiles", n_all, len(ALL_TILES))''',
    '''            log.info("     %d image pieces on %d tiles", n_all, len(proj_tiles))''')

# ---- cell 10: audit counts the pilot months as complete ----
sub(10, """{'with all 12 months':>18s}""", """{'with all run months':>19s}""")
sub(10, '''            whole += months == list(range(1, 13))''',
    '''            whole += months == (sorted(ONLY_CALENDAR_MONTHS) if ONLY_CALENDAR_MONTHS else list(range(1, 13)))''')
sub(10, '''{whole:18d}''', '''{whole:19d}''')
sub(10, '''#   root  = "/content/drive/MyDrive/MSUGWB/OPENET_v2_1_monthly_30m_tiles"''',
    '''#   root  = "/content/drive/MyDrive/MSUGWB/OPENET_v2_1_monthly_30m_tiles_KS_pilot"''')

# ---- new last cell: 2025 availability (report only) ----
CELL11 = '''# ═══════════════ CELL 11 (PILOT): is 2025 May-September in this release? Report only, nothing is fetched ═══════════════
# Pilot rule: report whether 2025 is available in the same release with complete May-September coverage, and add it
# only after that is agreed and recorded in the decision log. "Images" means OpenET has an image over the tile for
# that month; whether every pixel is covered is seen only after a download.
_want = [(2025, m) for m in range(5, 10)]      # the pilot asks about complete May-September coverage
_tiles = [P for P in ALL_TILES if not ONLY_TILES or P["key"] in ONLY_TILES]
print(f"OpenET {OPENET_VERSION} record over the aquifer: {ym(*min(AVAILABLE))} to {ym(*max(AVAILABLE))}")
_n_ok = 0
for P in _tiles:
    _have = [m for (y, m) in _want if P["months"].get((y, m))]
    _n_ok += len(_have) == len(_want)
    print(f"  {P['key']}: 2025 months with images: {span(_have) or 'none'}"
          + ("" if len(_have) == len(_want) else "   <- incomplete"))
print(f"{_n_ok} of {len(_tiles)} tiles have every month {span([m for _, m in _want])} of 2025")
_all25 = sorted({m for P in _tiles for (y, m) in P["months"] if y == 2025})
print(f"2025 months with images on any of these tiles: {span(_all25) or 'none'}")
'''
cells.append({"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [],
              "source": CELL11.splitlines(keepends=True)})

# ---- new first cell: what this copy is ----
HEAD = '''# PILOT copy: OpenET v2.1, Kansas only, whole years 2016-2024

A copy of `OpenET_Monthly_30m_Ogallala.ipynb` for the pilot (`GW_ET_High_Plains_pilot.docx`, Task 1.4). The
pipeline is unchanged; only these differ (each marked `PILOT` in the code):

| Setting | Value | Why |
|---|---|---|
| `OPENET_VERSION` | `v2_1` | one ET release, as the pilot requires |
| `START_MONTH`, `STOP_MONTH` | `2016-01` to `2024-12` | pilot years, whole years |
| `ONLY_CALENDAR_MONTHS` (new) | `None` | whole years. The pilot's analyses use May-September; the other months give annual and water-year totals (to compare with WIMAS reported use) and an off-season baseline, and are never fetched twice. `[5, 6, 7, 8, 9]` fetches the growing season only (~5/12 of the work) |
| `ONLY_TILES` | 11 zone-14 tiles | they hold the whole Kansas study area (GMD 1-5 and the KGS HPA extent) |
| `TEST_YEAR`, `TEST_MONTHS` | 2023, June-July | the pilot asks for a low-deficit rehearsal year; 2022 was the driest season |
| output folder | `MyDrive/MSUGWB/OPENET_v2_1_monthly_30m_tiles_KS_pilot/` | never mixed with a full-aquifer run |

Pixels kept are still those inside the whole-aquifer polygon, so the tiles also hold the parts of Nebraska,
Colorado and Oklahoma they cover; cut to Kansas cells later.

**Known gap.** OpenET's zone-14 images stop at about -102.027 deg, so a ~2 km strip of Kansas along the Colorado
border has no data in these tiles (the "coverage gap near the Colorado-Kansas border" the pilot mentions).
Exclude those cells with that reason, or add the zone-13 tiles named in `ONLY_TILES`.

**Run order.**
1. *Run all* with `TEST_MODE = True` (2 of the 11 tiles, June-July 2023, into `_test/`). Read the runtime and
   size projection at the end of Cell 9.
2. Set `TEST_MODE = False`, *Run all*; rerun after any stop. `ONLY_TILES` can be split over Colab sessions
   (never two sessions on one tile).
3. Cell 11 reports whether 2025 May-September exists in v2.1 (report only; do not add 2025 until it is logged).
'''
cells.insert(0, {"cell_type": "markdown", "metadata": {}, "source": HEAD.splitlines(keepends=True)})

for c in cells:
    if c["cell_type"] == "code":
        c["outputs"], c["execution_count"] = [], None
json.dump(nb, open(DST, "w", encoding="utf-8", newline="\n"), ensure_ascii=False, indent=1)
print("written", DST, len(cells), "cells")
