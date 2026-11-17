"""Pilot settings: paths, years, tiles, targets. The ONE file to edit for where and how a run happens.

Values that are PRE-REGISTERED (thresholds, windows, exclusions) are NOT here: they live in
pilot/decisions.yaml with their reason and status, and are read with decisions.value("Dxx"). This file only
holds what the pilot document itself fixes (years, targets, months) and technical settings.

Where and how a run happens:
    RUN_ENV   "local" -> data on H: and the Google Drive desktop mount (G:/My Drive)
              "colab" -> data on /content/drive/MyDrive (inputs synced there by `python -m pilot sync-inputs`)
    RUN_MODE  "test"   -> synthetic data root made by pilot/tests/synth.py (seconds; decisions may be proposed)
              "smoke"  -> real data, one county (SMOKE_COUNTY); a code test only. Decisions may be proposed;
                          outputs go to <out>/_smoke and are NOT pilot results (never reported, never used)
              "stage2" -> real data, the processing-test counties of decision D27 (pilot Stage 2)
              "stage3" -> real data, every eligible Kansas cell (pilot Stage 3)
Both can be set without editing this file: environment variables PILOT_RUN_ENV / PILOT_RUN_MODE.
Standard library only.
"""
import os
import tempfile

# ============================================================== where and how
RUN_ENV = os.environ.get("PILOT_RUN_ENV", "local")       # "local" | "colab"
RUN_MODE = os.environ.get("PILOT_RUN_MODE", "test")      # "test" | "smoke" | "stage2" | "stage3"
REAL_MODES = ("smoke", "stage2", "stage3", "main")       # modes that read real data
# "main" -> every eligible Kansas cell, AFTER the freeze (pilot/freeze.json): the frozen workflow with dry cohort ET
#           released by the blinding gate; outputs go to <out>/_main, apart from the blinded pilot outputs

PKG_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(PKG_DIR).replace("\\", "/")  # H:/MSU_GWB/datasets locally

# ---- local roots ----
LOCAL_DRIVE = os.environ.get("PILOT_LOCAL_DRIVE", "G:/My Drive/MSUGWB")   # Google Drive desktop mount
LOCAL_DRIVE_TOP = os.path.dirname(LOCAL_DRIVE)                            # G:/My Drive
LOCAL_OUT = os.environ.get("PILOT_LOCAL_OUT", "H:/MSU_GWB/pilot_out")
LOCAL_SCRATCH = os.environ.get("PILOT_LOCAL_SCRATCH", "H:/MSU_GWB/pilot_scratch")
# Processes that read OpenET tile-years side by side (Tasks 1.4 and et_cells). 1 = one process (the default).
# More workers need about 0.5 GB of memory each on a whole-tile window (Stage 3).
N_WORKERS = int(os.environ.get("PILOT_WORKERS", "1"))
MERGED = f"{REPO_ROOT}/merged_datasets"
PILOT_TEST = f"{REPO_ROOT}/pilot test"
DOWNLOADED = f"{PILOT_TEST}/downloaded_2026-11-13"

# ---- Colab roots ----
COLAB_MYDRIVE = "/content/drive/MyDrive"
COLAB_DRIVE = f"{COLAB_MYDRIVE}/MSUGWB"
COLAB_INPUTS = f"{COLAB_DRIVE}/PILOT_inputs"           # mirror of the local-only inputs (sync-inputs)
COLAB_OUT = f"{COLAB_DRIVE}/PILOT_KS"
COLAB_SCRATCH = "/content/pilot_scratch"

# ---- test mode ----
SYNTH_ROOT = os.path.join(tempfile.gettempdir(), "pilot_synthetic").replace("\\", "/")

SESSION_BUDGET_H = 11.0      # a task stops cleanly (and resumes next session) after this many hours
LOCK_STALE_H = 13.0
CODE_IN_SIGNATURE = True     # a task's output folder changes when its code changes
GATE_OVERRIDE = []           # (kept for the copied runner; the pilot has no gates)

# ============================================================== source datasets
# key -> (local path, Colab path, synthetic name). "{drive}" = the MSUGWB folder of that environment,
# "{inputs}" = local file mirrored to PILOT_inputs on Colab. A path may be a folder (tiles, AIM-HPA years).
_IN = "{inputs}"
DATASETS = {
    # boundaries
    "counties_zip": (f"{LOCAL_DRIVE}/P1_rain_only/_shared/counties/cb_2023_us_county_500k.zip",
                     "{drive}/P1_rain_only/_shared/counties/cb_2023_us_county_500k.zip", "counties.zip"),
    "aquifer_shp": (f"{REPO_ROOT}/HPA_polygon/hp_bound2010.shp", _IN + "/HPA_polygon/hp_bound2010.shp", "aquifer.gpkg"),
    "gmd": (f"{DOWNLOADED}/ks_gmd/ks_gmd.geojson", _IN + "/ks_gmd/ks_gmd.geojson", "gmd.geojson"),
    "plss_sections": (f"{DOWNLOADED}/blm_plss_ks/ks_plss_sections_west_of_97W.geojson",
                      _IN + "/blm_plss_ks/ks_plss_sections_west_of_97W.geojson", "sections.geojson"),
    "kgs_extent": (f"{PILOT_TEST}/kansas geological survey/High_Plains_Aquifer_extent/High_Plains_Aquifer.shp",
                   _IN + "/kgs/High_Plains_Aquifer_extent/High_Plains_Aquifer.shp", "kgs_extent.gpkg"),
    # aquifer
    "kgs_sections": (f"{DOWNLOADED}/kgs_arcgis/kgs_hpa_section_properties.geojson",
                     _IN + "/kgs_arcgis/kgs_hpa_section_properties.geojson", "kgs_sections.geojson"),
    "kgs_satthick_1921": (f"{DOWNLOADED}/kgs_arcgis/kgs_hpa_satthick_2019_2021_by_section.csv",
                          _IN + "/kgs_arcgis/kgs_hpa_satthick_2019_2021_by_section.csv", "kgs_satthick_1921.csv"),
    "kgs_bedrock_points": (f"{PILOT_TEST}/kansas geological survey/High_Plains_Aquifer.csv",
                           _IN + "/kgs/High_Plains_Aquifer.csv", "kgs_bedrock_points.csv"),
    "usgs_bottom": (f"{PILOT_TEST}/geodatabase of high plains aquifer/ds11HGHPLN_non_prop/11HGHPLN_ASCII_files/r11hghpln_bot_a.txt",
                    _IN + "/usgs_hpa_surfaces/r11hghpln_bot_a.txt", "usgs_bottom.tif"),
    "usgs_change": (f"{MERGED}/USGS_HPA_native/USGS_HPA_water_level_change_500m.nc",
                    _IN + "/USGS_HPA_native/USGS_HPA_water_level_change_500m.nc", "usgs_change.nc"),
    "haacker_dir": (f"{DOWNLOADED}/haacker_hydroshare", _IN + "/haacker_hydroshare", "haacker"),
    # wells and water rights
    "kgs_winter": (f"{DOWNLOADED}/kgs_arcgis/kgs_winter_water_levels_1991_2024.csv",
                   _IN + "/kgs_arcgis/kgs_winter_water_levels_1991_2024.csv", "kgs_winter.csv"),
    "wells_ts": (f"{MERGED}/HPA_wells_water_level_timeseries.nc", _IN + "/HPA_wells_water_level_timeseries.nc",
                 "wells_ts.nc"),
    "wimas_pods": (f"{DOWNLOADED}/kgs_arcgis/kgs_wimas_points_of_diversion.csv",
                   _IN + "/kgs_arcgis/kgs_wimas_points_of_diversion.csv", "wimas_pods.csv"),
    "master_wells": (f"{DOWNLOADED}/kgs_arcgis/kgs_master_well_inventory_linked.csv",
                     _IN + "/kgs_arcgis/kgs_master_well_inventory_linked.csv", "master_wells.csv"),
    "wwc5": (f"{PILOT_TEST}/ Water Well Completion Records (WWC5) Database - Kansas/wwc5_wells.txt",
             _IN + "/wwc5/wwc5_wells.txt", "wwc5.txt"),
    # land, irrigation, soils, terrain
    "aim_dir": (f"{DOWNLOADED}/aim_hpa_hydroshare", _IN + "/aim_hpa_hydroshare", "aim"),
    # LANID 1997-2020 on the v3 grid, inside the aquifer polygon (local data processing/build_lanid_1997_2020.py;
    # 1997-2017 are the v3 release value for value, 2018-2020 the extension release). Primary irrigation map (D06).
    "lanid": (f"{MERGED}/LANID_30m_1997-2020_Ogallala.nc", _IN + "/LANID_30m_1997-2020_Ogallala.nc", "lanid.nc"),
    # LANID 2018-2020 as released: whole CONUS (lanid2018.tif ...), so footprints off the aquifer are covered
    "lanid_1820_dir": (f"{PILOT_TEST}/LANID extended", _IN + "/LANID_extended", "lanid_1820"),
    "gssurgo": (f"{MERGED}/USGS data/derived_usgs/GSSURGO_Ogallala_all_attrs_120m.nc",
                "{drive}/USGS data/derived_usgs/GSSURGO_Ogallala_all_attrs_120m.nc", "gssurgo.nc"),
    "dem": (f"{MERGED}/GLO30_DEM_Ogallala_1arcsec.tif",
            f"{COLAB_MYDRIVE}/MGU GWB datasets/merged_datasets (1)/GLO30_DEM_Ogallala_1arcsec.tif", "dem.tif"),
    # weather and soil moisture
    "gridmet": (f"{MERGED}/GRIDMET_Ogallala_lab_latest.nc",
                f"{COLAB_MYDRIVE}/MGU GWB datasets/merged_datasets (1)/GRIDMET_Ogallala_lab_latest.nc", "gridmet.nc"),
    "smap_pilot": (f"{LOCAL_DRIVE}/SMAP_L3E_Kansas_pilot", "{drive}/SMAP_L3E_Kansas_pilot", "smap_pilot"),  # folder
    "smap_reprocessed": (f"{MERGED}/SPL3SMP_E_Ogallala_reprocessed_20261020.nc",      # same bytes as the Drive copy
                         f"{COLAB_MYDRIVE}/SMAP_L3E_Ogallala/SPL3SMP_E_Ogallala_reprocessed_20261020.nc", "smap.nc"),
    "smap_full_old": (f"{MERGED}/SPL3SMP_E_VOD_Ogallala_FULL.nc", _IN + "/SPL3SMP_E_VOD_Ogallala_FULL.nc", "smap.nc"),
    # 30 m tiles (UTM zone 14)
    "cdl_tiles": (f"{LOCAL_DRIVE}/USDA_NASS_CDL_30m_tiles", "{drive}/USDA_NASS_CDL_30m_tiles", "cdl_tiles"),
    "openet_tiles": (f"{LOCAL_DRIVE}/OPENET_v2_1_monthly_30m_tiles_KS_pilot",
                     "{drive}/OPENET_v2_1_monthly_30m_tiles_KS_pilot", "openet_tiles"),
    # benchmark and audit
    "mieno_repo": (f"{DOWNLOADED}/mieno_github/Drought-Production-Risk-Aquifer",
                   _IN + "/mieno_github/Drought-Production-Risk-Aquifer", "mieno"),
}
CDL_TAG = "USDA_NASS_CDL_30m_tiles"
OPENET_TAG = "OPENET_v2_1_monthly_30m_tiles"    # file names inside the _KS_pilot folder keep the release tag
ROI_LABEL = "Ogallala"
# The SMAP file used by Tasks 1.2/1.3: the unmasked Kansas re-extraction when it exists. The fallbacks hold
# data only inside the aquifer polygon (and the reprocessed one ends 2022-06), so results from them are a
# code test, flagged "smap_source=fallback" in every output.
SMAP_ORDER = ["smap_pilot", "smap_reprocessed"]

# Earlier analyses archived by Task 0.1 and audited by Task 0.4 (repository-relative)
EARLIER_OUTPUTS = ["gw_buffering_xgb", "residual_et_store", "tbi_maps_v3_usgs", "XGB model", "temporal buffering",
                   "residual_et", "P1_rain_only"]

# ============================================================== what the pilot document fixes
STATE_FIPS = "20"                               # Kansas
PANEL_YEARS = list(range(2016, 2025))           # growing seasons 2016-2024
LABEL_YEARS = list(range(2011, 2016))           # cohort: irrigated in >= 4 of these
CORN_YEARS = [2013, 2014, 2015]                 # cohort: corn in >= 2 of these
LATER_IRR_YEARS = [2016, 2017, 2018, 2019, 2020]   # rainfed cells: "check whether any were irrigated later"
#                                                    (LANID through 2020; AIM-HPA through 2017)
FOOTPRINT_CHECK_YEARS = [2018, 2019, 2020]      # reference footprints checked against these labels (LANID 2018-20)
SEASON_MONTHS = [5, 6, 7, 8, 9]
PEAK_CANDIDATES = [6, 7, 8]
PRESEASON_MONTHS = [1, 2, 3, 4]
WELL_YEARS = list(range(2011, 2025))            # pre-season water levels kriged (2011-2016 overlap Haacker)
COHORT_MIN_IRR_YEARS = 4
COHORT_MIN_CORN_YEARS = 2
SMAP_GAPS_EXPECTED = [("2019-06-19", "2019-07-23"), ("2022-08-01", "2022-09-30")]   # pilot text; confirmed in 1.3
REP_REGIONS = {                                  # replication candidates (metadata only)
    "texas_high_plains": dict(state="48", lat_max=36.5, lon_min=-103.1, lon_max=-100.0),
    "sw_nebraska": dict(state="31", lat_max=41.4, lon_min=-102.1, lon_max=-99.5),
}
SMOKE_COUNTY = "20081"                          # Haskell, KS (GMD 3): smoke mode reads only this county

# ---- 30 m tiles: the 11 UTM zone-14 tiles that hold the Kansas study area (OpenET pilot notebook) ----
KS_TILES = ["z14e01n26", "z14e01n27", "z14e01n28", "z14e02n26", "z14e02n27", "z14e02n28",
            "z14e03n26", "z14e03n27", "z14e03n28", "z14e04n26", "z14e04n27"]
PIXEL_M = 30.0
TILE_PX = 5120
LATTICE_ORIGIN_M = 15.0                         # Landsat pixel edges sit at 15 m mod 30

# ---- EASE-Grid 2.0 global (EPSG:6933) ----
EASE_X_MIN = -17367530.445161372
EASE_Y_MAX = 7314540.830638504
EASE_M36 = 36032.220840584
EASE_M09 = 9008.055210146

# ---- CDL classes used for crop history ----
CDL_CLASSES = dict(corn=[1, 12, 13], soybean=[5], sorghum=[4], wheat=[22, 23, 24, 26, 225, 236, 238],
                   alfalfa_hay=[36, 37], cotton=[2], other_crop=list(range(3, 4)) + list(range(6, 12))
                   + list(range(14, 22)) + [25, 27, 28, 29, 30] + list(range(31, 36)) + list(range(38, 62))
                   + list(range(66, 81)) + list(range(204, 225)) + list(range(226, 236)) + [237]
                   + list(range(239, 255)))
CROPLAND_CODES = sorted(c for v in CDL_CLASSES.values() for c in v) + [61]     # 61 = fallow/idle cropland

# ---- OpenET ----
ET_MODELS = ["disalexi", "eemetric", "geesebal", "ptjpl", "sims", "ssebop"]
ET_ENSEMBLE = "et_ensemble_mad"

# ---- HPA aquifer codes of the well tables: keep 11x / 12x (Quaternary, Pleistocene, Tertiary Ogallala) ----
HPA_CODE_PREFIXES = ("11", "12")

# ============================================================== targets (pilot Phase 2; reported, never enforced)
TARGETS = dict(
    smap_dry_share=0.70,            # 2.1 dry county-years kept by the soil-moisture mask
    fp_min=3, fp_county_share=0.80,  # 2.2 reference footprints per county
    et_min_models=4, et_cell_share=0.80,   # 2.3
    well_km=8.0, well_share=0.60,   # 2.4 cohort cells within 8 km of a measured well, each year
    both_thirds_counties=15, both_thirds_cells=50,
    thin_thick_se_ratio=3.0,
    ci_halfwidth_mm=15.0,           # 2.6
    min_dry_episodes=3,
)
CAPACITY_N_START = 20
CAPACITY_N_EXTEND = 200
CAPACITY_SEED = 20261114

# ============================================================== figures
COLORS = dict(thin="#c2410c", middle="#a3a3a3", thick="#1d4ed8", rainfed="#15803d", cohort="#7c3aed",
              excluded="#d4d4d4", footprint="#0891b2")


# ============================================================== derived (do not edit)
def check():
    if RUN_ENV not in ("local", "colab"):
        raise ValueError(f"RUN_ENV must be 'local' or 'colab', not {RUN_ENV!r}")
    if RUN_MODE not in ("test",) + REAL_MODES:
        raise ValueError(f"RUN_MODE must be one of test, smoke, stage2, stage3, main, not {RUN_MODE!r}")


def real_data():
    return RUN_MODE in REAL_MODES


def src(key):
    """Absolute path of a source dataset in the current environment and mode."""
    local, colab, synth = DATASETS[key]
    if RUN_MODE == "test":
        return f"{SYNTH_ROOT}/inputs/{synth}"
    if key == "openet_tiles" and RUN_MODE == "smoke" and os.environ.get("PILOT_SMOKE_OPENET"):
        return os.environ["PILOT_SMOKE_OPENET"]      # code test of the ET chain only (never a pilot result)
    if RUN_ENV == "local":
        return local
    return colab.replace("{drive}", COLAB_DRIVE).replace("{inputs}", COLAB_INPUTS)


def source_id(key):
    """What a task signature records about a source it reads: the resolved path, and the size when it is a file.
    A task that reads a different file (the unmasked SMAP extraction instead of the fallback, another OpenET
    folder, a rebuilt LANID file) therefore gets a new signature and runs again instead of being skipped as
    complete. key "smap" = the file io.smap_source() picks. Folders are identified by path only."""
    if key == "smap":
        from .common import io
        try:
            p, label = io.smap_source()
        except Exception:                     # no SMAP file yet: the task itself reports that
            return ["none", None, "none"]
        return [p, os.path.getsize(p), label]
    p = src(key)
    return [p, os.path.getsize(p) if os.path.isfile(p) else None]


def data_root():
    return SYNTH_ROOT if RUN_MODE == "test" else (LOCAL_DRIVE if RUN_ENV == "local" else COLAB_DRIVE)


def scratch_root():
    base = LOCAL_SCRATCH if RUN_ENV == "local" else COLAB_SCRATCH
    return f"{SYNTH_ROOT}/scratch" if RUN_MODE == "test" else base


def out_root():
    if RUN_MODE == "test":
        return f"{SYNTH_ROOT}/out"
    base = LOCAL_OUT if RUN_ENV == "local" else COLAB_OUT
    return f"{base}/_{RUN_MODE}" if RUN_MODE in ("smoke", "stage2", "main") else base


def shared_root():
    """Downloads and caches reused by every mode (cell-id rasters, counties)."""
    return f"{out_root()}/_shared" if RUN_MODE == "test" else (
        f"{LOCAL_OUT if RUN_ENV == 'local' else COLAB_OUT}/_shared")


def snapshot():
    """Every setting (upper-case names) plus every decision value (D01 ...), for task signatures."""
    g = globals()
    snap = {k: g[k] for k in sorted(g) if k.isupper() and not callable(g[k])}
    from . import decisions
    snap.update(decisions.values())
    return snap
