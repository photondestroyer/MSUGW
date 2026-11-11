"""paper1 settings: the ONE file to edit.

Where and how a run happens:
    RUN_ENV   "colab" -> data from Google Drive (/content/drive/MyDrive/MSUGWB)
              "local" -> data from LOCAL_DATA_ROOT (a mirror of Drive's folder layout)
    RUN_MODE  "test"  -> a small window of TEST_TILES for TEST_YEARS (minutes)
              "focus" -> the focus-landscape tiles, all years (the plan's G0 fallback)
              "full"  -> every tile of the aquifer, all years
Both can be overridden without editing this file: environment variables P1_RUN_ENV / P1_RUN_MODE
(the Colab driver notebook sets them).

Every rule threshold and model setting of P1.1 / P1.2 lives here. docs/CONFIG.md is generated from
this file (python -m paper1 docs) and a test fails when it is out of date.
Standard library only: this module is imported by every environment (incl. the torch-only conda env).
"""
import os

# ============================================================== where and how
RUN_ENV = os.environ.get("P1_RUN_ENV", "colab")          # "colab" | "local"
RUN_MODE = os.environ.get("P1_RUN_MODE", "test")         # "test" | "focus" | "full"

COLAB_MYDRIVE = "/content/drive/MyDrive"
COLAB_DATA_ROOT = f"{COLAB_MYDRIVE}/MSUGWB"               # where the dataset notebooks write
COLAB_SCRATCH = "/content/p1_scratch"                     # fast local disk of the Colab VM

LOCAL_DATA_ROOT = os.environ.get("P1_LOCAL_DATA_ROOT", "G:/MSU_GWB/drive_mirror/MSUGWB")
LOCAL_SCRATCH = os.environ.get("P1_LOCAL_SCRATCH", "G:/MSU_GWB/p1_scratch")
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # G:/MSU_GWB/datasets locally
LOCAL_MERGED = os.path.join(REPO_ROOT, "merged_datasets")

OUTPUT_FOLDER = "P1_rain_only"          # outputs: <data root>/P1_rain_only[/_test|/_focus]/<stage>/<signature>/
SHARED_FOLDER = "P1_rain_only/_shared"  # downloads reused by every mode (NHDPlus HR subsets, counties, DEM tiles)

SESSION_BUDGET_H = 11.0      # a stage stops cleanly (and resumes next session) after this many hours
LOCK_STALE_H = 13.0          # a lock older than this is considered abandoned
CODE_IN_SIGNATURE = True     # a stage's output folder changes when its code changes (safer); False freezes it

# ============================================================== which tiles / years
TEST_TILES = ["z14e02n29"]                   # south-west Nebraska (Republican / Platte / Sandhills edge)
TEST_YEARS = [2015, 2016]                    # calendar years of 30 m data used in test mode (Landsat v2 on Drive starts in 2015)
TEST_WINDOW_PX = 1024                         # test mode reads a square of this many 30 m pixels (multiple of 16)
TEST_WINDOW_CENTER = None                     # (row, col) in frame pixels; None -> the frame centre
TEST_SPACE_GROUP_DEG = 0.05                   # spatial CV group size in test mode (the window is ~30 km)
TEST_SPACE_FOLDS = 2
TEST_YEAR_GROUPS = [[2015, 2016]]             # one year group: no year hold-out in test mode

# focus landscapes of the plan (section 3.3): representative points; their tiles are processed in "focus" mode
FOCUS_POINTS = {
    "FL1_central_platte": (-98.9, 40.75), "FL2_sandhills": (-101.0, 42.0), "FL3_republican": (-100.4, 40.2),
    "FL4_arkansas_great_bend": (-98.8, 38.1), "FL5_sheridan6": (-100.6, 39.45),
    "FL6_southern_high_plains": (-101.8, 33.6), "FL7_canadian_breaks": (-101.3, 35.6),
}

YEAR_FIRST, YEAR_LAST = 2008, 2025      # years with CDL + Landsat + OpenET (core period of the plan)
AGG_YEAR_FIRST = 2007                   # blocks are also aggregated for 2007 (pre-season of WY2008)
DRIVER_MONTH_FIRST = "2005-01"          # monthly drivers kept from here (24-month history of WY2008)

# ============================================================== source datasets (registry)
# "drive": path relative to the data root (or absolute); "local": override when RUN_ENV == "local"
# (files already in merged_datasets need no copy). "tiles" datasets follow the 30 m tile layout
# <TAG>/<year>/<tile>/<TAG>_<tile>_<year>_Ogallala.nc (NWI has no year level).
SOURCE_SUBFOLDER = ""        # "" = the finished products; "/_test" would read the notebooks' test outputs
ROI_LABEL = "Ogallala"
LANDSAT_TAG = "LANDSAT_C2L2_indices_30m_monthly_tiles_v2_latest"   # on Drive 2026-11-11: 2015-2026 (2008-2014 to come)
OPENET_TAG = "OPENET_v2_0_monthly_30m_tiles"                       # release 2.0 (1999-2024); the v2_1 folder is empty
CDL_TAG = "USDA_NASS_CDL_30m_tiles"
NWI_TAG = "NWI_wetlands_30m_tiles"

DATASETS = {
    "landsat": dict(kind="tiles", tag=LANDSAT_TAG, per_year=True,
                    vars=["NIRv", "NDMI", "NDVI", "MNDWI", "AWEIsh", "n_obs", "n_obs_oli", "obs_day", "roi_mask"]),
    "openet": dict(kind="tiles", tag=OPENET_TAG, per_year=True,
                   vars=["et_ensemble_mad", "et_ensemble_mad_min", "et_ensemble_mad_max", "et_ensemble_mad_count",
                         "count_ssebop", "roi_mask"]),
    "cdl": dict(kind="tiles", tag=CDL_TAG, per_year=True, vars=["cropland", "roi_mask"]),
    "nwi": dict(kind="tiles", tag=NWI_TAG, per_year=False,
                vars=["wetland_system", "riparian", "nwi_mapped", "roi_mask"]),
    "gridmet": dict(kind="file", drive=f"{COLAB_MYDRIVE}/MGU GWB datasets/merged_datasets (1)/GRIDMET_Ogallala_lab_latest.nc",
                    local=os.path.join(LOCAL_MERGED, "GRIDMET_Ogallala_lab_latest.nc"),
                    vars=["pr", "pet", "vpd", "srad", "tmmn", "tmmx", "roi_fraction"]),
    "dem": dict(kind="file", drive=f"{COLAB_MYDRIVE}/MGU GWB datasets/merged_datasets (1)/GLO30_DEM_Ogallala_1arcsec.tif",
                local=os.path.join(LOCAL_MERGED, "GLO30_DEM_Ogallala_1arcsec.tif")),
    "drought": dict(kind="file", drive=f"{COLAB_MYDRIVE}/MGU GWB datasets/merged_datasets (1)/GRIDMET_DROUGHT_Ogallala_lab_colab.nc",
                    local=os.path.join(LOCAL_MERGED, "GRIDMET_DROUGHT_Ogallala_lab_colab.nc"), vars=["spi90d"]),
    "hres_wtd": dict(kind="file", drive="HRES_WTD/HRES_WTD_Ogallala.nc", vars=["b1", "source_image"],
                     refuse=["HRES-WTD_2015_Ogallala.nc"]),       # the stretched-longitude export
    "lanid": dict(kind="file", drive="LANID/LANID_v3_30m_1997-2017_Ogallala.nc",
                  local=os.path.join(LOCAL_MERGED, "LANID_v3_30m_1997-2017_Ogallala.nc")),
    "gssurgo": dict(kind="file", drive="USGS data/derived_usgs/GSSURGO_Ogallala_all_attrs_120m.nc",
                    local=os.path.join(LOCAL_MERGED, "USGS data", "derived_usgs", "GSSURGO_Ogallala_all_attrs_120m.nc"),
                    vars=["aws0_100", "aws0_150", "rootznaws", "nccpi3all", "soc0_30"]),
    "gde": dict(kind="file", drive="USGS data/derived_usgs/GDE_30arcsec_Ogallala_CF.nc",
                local=os.path.join(LOCAL_MERGED, "USGS data", "derived_usgs", "GDE_30arcsec_Ogallala_CF.nc"),
                vars=["GDE_frac_AA", "AA_frac_GA"]),
    "aquifer": dict(kind="file", drive=f"{COLAB_MYDRIVE}/high_plains_quifer.zip",
                    local=os.path.join(REPO_ROOT, "HPA_polygon", "hp_bound2010.shp"),
                    member="high_plains_quifer/hp_bound2010.shp"),
}
LOCATE_BY_NAME = True        # if a "file" is not at its registry path, search the data root (depth 3) by name

# ============================================================== 30 m tile geometry (from the notebooks)
PIXEL_M = 30.0
TILE_PX = 5120
LATTICE_ORIGIN_M = 15.0      # Landsat pixel edges sit at 15 m mod 30
BLOCK_PX = 16                # 16 x 16 pixels = 480 m blocks, nested in the tile lattice
BLOCK_MIN_ROI_FRACTION = 0.5 # a block is part of the analysis when >= half its pixels are inside the aquifer
LANDSAT_FILL = -32768
NOBS_FILL = 255

# ============================================================== P1.0 block aggregates (plan v3)
HIST_BINS = 16                     # histogram of the 30 m values per block-month (32 is the plan's fallback)
HIST_EDGES = dict(ET=(0.0, 240.0), NIRv=(-0.05, 0.55), NDMI=(-0.45, 0.65))   # outer bins also take values beyond
OPENET_COUNT_VARS = ["count_ssebop"]   # clear-scene count per pixel-month: mean of these OpenET model counts
OBS_REGIMES = [("L5+L7", "2005-01"), ("L7", "2012-01"), ("L7+L8", "2013-04"), ("L8", "2017-01"), ("L8+L9", "2022-01")]
P10_VALID_SHARE_MIN = 0.70         # P1.0 check: valid block-months Apr-Sep 2012 per class x decile

# ============================================================== Step 0
DRIVER_VARS = ["pr", "pet", "vpd", "srad", "tmmn", "tmmx"]   # lab gridMET daily -> monthly
DRIVER_MIN_DAY_FRACTION = 0.9      # a month needs this share of usable days (sums scaled to the full month)
RAIN_WINDOWS_D = [30, 60, 90, 180, 365]   # rain summed over [month start - w, month start - 1]
ARIDITY_BASE_YEARS = (1985, 2014)  # aridity index P / PET climatology (same baseline as residual_et)
QC_CONSTANT, QC_IDENTICAL, QC_NODATA = 1, 2, 4                # lab per-record QC bits
QC_CONSTANT_IS_BAD = ("srad", "tmmn", "tmmx", "vpd", "pet")   # rain keeps constant (all-dry) days
BUCKET_CAPACITY_BINS_MM = [25, 50, 75, 100, 150, 200, 250, 300]
BUCKET_LP = 0.6                    # AET = PET * min(1, S / (LP * capacity))   (HBV soil routine)
BUCKET_SPINUP_FROM = "1979-01-01"
AWS_FOR_BUCKET = "aws0_150"
SPI_VAR = "spi90d"                 # lab drought product variable used for the drought class (stratum only)
DROUGHT_DRY_SPI, DROUGHT_WET_SPI = -0.8, 0.8   # drought class: SPI-90 <= dry | between | >= wet        # gSSURGO available water storage used as bucket capacity (cm -> mm)

STATIC_SOIL_VARS = ["aws0_100", "aws0_150", "rootznaws", "nccpi3all", "soc0_30"]
WTD_SHALLOW_M = [2, 5, 10]         # S0.5 extras: fraction of the block shallower than these depths

# HAND (S0.6)
NHD_GDB_URL = ("https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHDPlusHR/VPU/Current/GDB/"
               "NHDPLUS_H_{hu4}_HU4_GDB.zip")
TNM_API = "https://tnmaccess.nationalmap.gov/api/v1/products"
HU4_CODES = None                   # None -> asked from the TNM Access API for the aquifer bounding box
STREAM_FCODES = [46000, 46003, 46006, 55800, 33400]   # stream/river (unspecified, intermittent, perennial), artificial path, connector
INCLUDE_EPHEMERAL = False          # FCode 46007
MIN_DRAINAGE_KM2 = 10.0            # flowlines draining less than this (VAA TotDASqKm) are not channels
WATERBODY_FTYPES = [390, 436]      # lake/pond, reservoir; playas (361) and swamp/marsh (466) are not drainage
BURN_DEPTH_M = 10.0
HU8_BUFFER_M = 2000.0
CANAL_FTYPES = [336]               # NHD canal / ditch; kept as their own layer for the canal buffer (never channels)
DEM_COG_URL = ("https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM/"
               "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM.tif")
COUNTIES_URL = "https://www2.census.gov/geo/tiger/GENZ2023/shp/cb_2023_us_county_500k.zip"
AQUIFER_STATE_FIPS = ["08", "20", "31", "35", "40", "46", "48", "56"]   # CO KS NE NM OK SD TX WY

# ============================================================== P1.1 rules (block-year)
WTD_MIN_M = 30.0                   # strict depth rung: used wherever DEPTH_ACCEPTED has no entry
WTD_SENSITIVITY_M = [15.0, 60.0]   # rung bits WTD15 / WTD60 (WTD10 and WTD30 come from DEPTH_RUNGS_M)
DEPTH_RUNGS_M = [60.0, 30.0, 15.0, 10.0]   # ladder fixed in advance (plan v3, P1.1)
DEPTH_BY_FUNCTIONAL = dict(winter_cereal=15.0, summer_crop=15.0, other_annual=15.0, perennial_grass=15.0,
                           fallow=15.0, shrub=30.0, hay=30.0)   # class rung; woody blocks always 30 m
WOODY_DEPTH_M = 30.0
WOODY_NOISE_PX = 3                 # a block with >= this many shrub/forest pixels has woody cover
DEPTH_ACCEPTED = {}                # {"<functional>|<decile>": rung m} copied from the ring-test report after review
WTD_MIN_COVERAGE = 0.95            # share of the block's aquifer pixels with a valid HRES-WTD cell
HAND_MIN_M = 10.0
HAND_PASS_FRACTION = 0.9           # share of the block's pixels with HAND > HAND_MIN_M
IRR_BUFFER_M = 800.0               # no LANID ever-irrigated pixel within this distance
LANID_FREQ_MIN_YEARS = 3           # sensitivity set: exclude only pixels irrigated in at least this many years
CANAL_BUFFER_M = 800.0             # no NHD canal / ditch within this distance
FIRE_PERIMETERS = None             # slot: MTBS perimeters (path); None -> FIRE passes and is reported as not assessed
IRRIGATION_REGISTRIES = []         # future: [{"name": "KS_WIMAS", "path": ..., "lon": ..., "lat": ...}]
WETLAND_NOISE_PX = 3               # NWI wetland/riparian pixels tolerated in a block
NWI_MISSING_POLICY = "fail"        # block without an NWI tile: "fail" | "pass" (only for rehearsals)
TEST_NWI_MISSING_POLICY = "pass"   # test mode only (user, 2026-11-11: no NWI tile on Drive yet); reported as not applied
WATER_NOISE_PX = 3                 # open-water pixel-months tolerated in a block-year ...
WATER_MAX_MONTHS = 1               # ... and months with any open water
WATER_MIN_NOBS = 2                 # an open-water pixel-month needs at least this many observations
GDE_SENSITIVITY_MAX = 0.10         # GDE is NOT a main rule (user, 2026-11-04); sensitivity set only
CLASS_PURITY_MIN = 0.8             # dominant CDL group share of the block
DRYLAND_SHARE_MIN = 0.5            # crop is dryland-eligible in a county-year when this share lies outside LANID
DRYLAND_MIN_PIXELS = 1000
MIN_VALID_MONTHS = 6               # a block-year needs this many months with all three targets
TARGET_MIN_PIXEL_FRACTION = 0.5    # a block-month target needs this share of the block's aquifer pixels
MIN_REF_BLOCKS_PER_STRATUM = 500
V2_TOLERANCE = 0.10                # water-year ET <= (P + storage change) * (1 + tolerance)
V2_PASS_SHARE = 0.95
V4_SMD_MAX = 0.25                  # representativeness: |standardised mean difference| per decile
V4_VARS = ["aws0_150", "rootznaws", "nccpi3all", "elevation", "P_annual", "PET_annual"]

# ============================================================== gates (plan v3, section 6)
GATE_OVERRIDE = []                 # gate names (e.g. "G1a") whose failure may not halt focus/full runs

# ============================================================== P1.2 model
TARGETS = ["ET", "NIRv", "NDMI"]
QUANTILES = [0.05, 0.5, 0.95]
ALPHA = 0.10
DYN_INPUTS = ["P", "PET", "VPD", "SRAD", "TMIN", "TMAX", "BUCKET", "P90d", "P365d"]
STATIC_INPUTS = ["aws0_100", "aws0_150", "rootznaws", "nccpi3all", "soc0_30", "elevation"]
MEAS_INPUTS = ["oli_share", "nobs_mean", "obs_day_mean", "ET_scenes_mean"]   # how a value was observed (rule R2 holds)
SPECIES_INPUT = True               # CDL group one-hot as a secondary class input (functional group is primary)
EXPOSURE_VARS = ["WTD", "wtd_min", "wtd_mean", "wtd_p05", "HAND", "hand_frac", "hand_min", "hand_mean",
                 "gde_max", "gde_mean", "nwi_wetland_px", "nwi_riparian_px", "lanid_dist_m", "lanid_assessed",
                 "irrigated", "maxIrrExt", "county", "county_share", "dryland_share", "lat", "lon", "x", "y",
                 "block_id", "M0", "ET_clim", "reference", "flags", "facc_max_km2", "facc_mean_km2", "canal_dist_m",
                 "lanid3_dist_m", "state", "reference_ring", "depth_used_m", "depth_class_m", "SPI90", "spi", "drought"]
PRESEASON_MONTHS = 3
SPACE_GROUP_DEG = 0.5
SPACE_FOLDS = 5
YEAR_GROUPS = [list(range(2008, 2014)), list(range(2014, 2020)), list(range(2020, 2026))]
EMBARGO_YEARS = 1
INNER_SPLIT = (0.70, 0.15, 0.15)   # train / validation (early stopping) / calibration, by spatial group
CAL_YEAR_HOLDOUT = False           # True: also withhold one training year group for calibration
CAL_STRATA = ["target", "functional", "aridity_tercile", "season", "drought", "regime"]
CAL_MERGE_ORDER = ["regime", "aridity_tercile", "drought", "functional", "season"]   # dropped in this order
CAL_MIN_N = 1000                   # calibration block-months a stratum needs before it gets its own margin
TEST_CAL_MIN_N = 30
SAMPLE_BLOCKS_PER_STRATUM = 1500   # blocks sampled per aridity decile x class (all their reference years)
TEST_SAMPLE_BLOCKS_PER_STRATUM = 10**9
M2_WINDOWS_MONTHS = [1, 2, 3, 4, 6, 9, 12]
M2_MIN_SAMPLES = 50
M0_MIN_YEARS = 3
RUN_ABLATION = True                # also train fold (0, 0) without pre-season inputs (leakage check)
RUN_MEAS_ABLATION = True           # ... and without measurement inputs
RUN_RATIO_MODEL = True             # ... and with ET / PET as the ET target (sensitivity; no mass penalty)
CLASS_FLIP = dict(share=0.05, seed=11)   # class-flip test: relabel this share of test blocks
CLASS_CONFUSION = dict(winter_cereal="fallow", fallow="winter_cereal", summer_crop="other_annual",
                       other_annual="summer_crop", perennial_grass="hay", hay="perennial_grass",
                       shrub="perennial_grass")
M1 = dict(hidden=64, tcn_levels=3, kernel=3, static_hidden=64, dropout=0.1, lr=1e-3, weight_decay=1e-4,
          batch=512, max_epochs=40, patience=5, mass_weight=0.1, seed=42, amp=True, num_workers=0)
TEST_M1 = dict(max_epochs=3, patience=2, batch=128)
PLANT = dict(share=0.10, ET_mm=20.0, NIRv=0.02, NDMI=0.04, seed=7)   # synthetic subsidy for verify (3)
PASS = dict(coverage=(0.88, 0.92), plant_bias=0.10, fpr_max=0.07, fpr_nominal=0.05,
            stratum_min_n=1000, mass_violation_max=0.01, drought_bias_sd=0.15, drought_bias_et_mm=3.0,
            regime_step_sd=0.1, ring_mean_sd=0.2, terrain_sd=0.1, crop_layer_sd_ratio=1.15, rotation_sd=0.2)
TEST_PASS_STRATUM_MIN_N = 30
DROUGHT_TEST_YEARS = dict(d2012=[2011, 2012, 2013], d2022=[2022])   # unseen droughts of gate G1b (water years)
DRIFT_YEARS = ((2008, 2017), (2018, 2025))                             # LANID drift test: before / after LANID ends


# ============================================================== derived (do not edit)
def check():
    if RUN_ENV not in ("colab", "local"):
        raise ValueError(f"RUN_ENV must be 'colab' or 'local', not {RUN_ENV!r}")
    if RUN_MODE not in ("test", "focus", "full"):
        raise ValueError(f"RUN_MODE must be 'test', 'focus' or 'full', not {RUN_MODE!r}")
    if TEST_WINDOW_PX % BLOCK_PX:
        raise ValueError("TEST_WINDOW_PX must be a multiple of BLOCK_PX")


def data_root():
    return COLAB_DATA_ROOT if RUN_ENV == "colab" else LOCAL_DATA_ROOT


def scratch_root():
    return COLAB_SCRATCH if RUN_ENV == "colab" else LOCAL_SCRATCH


def out_root():
    sub = {"test": "/_test", "focus": "/_focus", "full": ""}[RUN_MODE]
    return f"{data_root()}/{OUTPUT_FOLDER}{sub}"


def shared_root():
    return f"{data_root()}/{SHARED_FOLDER}"


def years():
    return list(TEST_YEARS) if RUN_MODE == "test" else list(range(YEAR_FIRST, YEAR_LAST + 1))


def agg_years():
    ys = years()
    return sorted(set([ys[0] - 1] + ys)) if RUN_MODE != "test" else ys


def year_groups():
    return [list(g) for g in (TEST_YEAR_GROUPS if RUN_MODE == "test" else YEAR_GROUPS)]


def space_group_deg():
    return TEST_SPACE_GROUP_DEG if RUN_MODE == "test" else SPACE_GROUP_DEG


def space_folds():
    return TEST_SPACE_FOLDS if RUN_MODE == "test" else SPACE_FOLDS


def cal_min_n():
    return TEST_CAL_MIN_N if RUN_MODE == "test" else CAL_MIN_N


def nwi_missing_policy():
    return TEST_NWI_MISSING_POLICY if RUN_MODE == "test" else NWI_MISSING_POLICY


def stratum_min_n():
    return TEST_PASS_STRATUM_MIN_N if RUN_MODE == "test" else PASS["stratum_min_n"]


def m1_settings():
    s = dict(M1)
    if RUN_MODE == "test":
        s.update(TEST_M1)
    return s


def snapshot():
    """Every public setting as plain JSON-able data (written into each output's provenance)."""
    g = globals()
    return {k: g[k] for k in sorted(g) if k.isupper() and not callable(g[k])}
