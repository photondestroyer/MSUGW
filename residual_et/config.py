"""Single source of truth: paths, windows, thresholds (with citations), seeds.

Nothing else in the package hard-codes a path, a threshold or a period.
Imports only the standard library so that both Python environments
(system 3.13 for prepare/XGBoost, conda SSM for the LSTM) can load it.
"""
import os
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parents[1]
MERGED = ROOT / "merged_datasets"
DERIVED = MERGED / "USGS data" / "derived_usgs"
STORE = Path(os.environ.get("RESIDUAL_ET_STORE", ROOT / "residual_et_store"))

SRC = dict(
    gridmet=MERGED / "GRIDMET_Ogallala_lab_latest.nc",
    drought=MERGED / "GRIDMET_DROUGHT_Ogallala_lab_colab.nc",
    ssebop=MERGED / "MODIS_ET_SSEBop_Merged_Ogallala.nc",
    smap=MERGED / "SPL4SMGP_Ogallala_FULL.nc",
    wtd=MERGED / "HRES-WTD_2015_Ogallala.nc",
    mcd12q1=MERGED / "MCD12Q1_Merged_Ogallala.nc",
    gfsad=MERGED / "GFSAD1000_V1_2019_Ogallala.nc",
    gir=MERGED / "Global_irrigation_Area_Merged_Ogallala.nc",
    hsg=MERGED / "Hydrologic_Soil_Group_250m_2019_Ogallala.nc",
    ameriflux=MERGED / "ameriflux.nc",
    gssurgo=DERIVED / "GSSURGO_Ogallala_all_attrs_120m.nc",
    dwl_long=DERIVED / "dWL_predev_to_2019_ft_4km.nc",
    dwl_recent=DERIVED / "dWL_2017_to_2019_ft_4km.nc",
    sites=ROOT / "AmeriFlux_NEON_sites.xlsx",
)

# ---------------------------------------------------------------- profiles
ROI_MIN_FRACTION = 0.5      # a cell is an aquifer cell if >= half of it lies in hp_bound2010


@dataclass(frozen=True)
class Profile:
    """A run configuration. `smoke` proves the code path on a small tile; it
    is not a scientific run (short rainfall baseline, 144 cells)."""
    name: str
    bbox: tuple | None            # (lon_min, lon_max, lat_min, lat_max) or None = whole aquifer
    model_t0: str = "2015-04-01"  # SMAP L4 starts 2015-03-31
    model_t1: str = "2022-06-01"  # exclusive; SSEBop ends with the 2022-05-21 dekad
    forcing_t0: str = "2013-04-01"   # two years of HBV spin-up before the first loss day
    base_met: tuple = (1985, 2014)   # 30-year gridMET baseline, ends before the model window
    base_et: tuple = (2003, 2014)    # SSEBop baseline, outside every split
    n_boot: int = 20                 # bootstrap replicates for the rainfall-anomaly fit

    @property
    def store(self) -> Path:
        return STORE / self.name


PROFILES = {
    "full": Profile("full", None),
    "smoke": Profile("smoke", (-101.25, -100.75, 37.75, 38.25),
                     base_met=(2005, 2014), n_boot=4),
}


def get_profile(name: str) -> Profile:
    if name not in PROFILES:
        raise KeyError(f"unknown profile {name!r}; choose from {sorted(PROFILES)}")
    return PROFILES[name]


# ---------------------------------------------------------------- seasons / windows
BLOCK_GROW_MONTHS = (4, 5, 6, 7, 8, 9)   # split blocks: growing Apr-Sep, dormant Oct-Mar
ANALYSIS_MONTHS = (5, 6, 7, 8, 9)        # buffering statistics (same window as TBI 11.15, A9)
ANTECEDENT_WINDOWS = (30, 90)            # days; trailing windows end the day BEFORE the dekad
MIN_VALID_FRACTION = 0.8                 # of days required in a dekad mean or trailing sum
LC_YEARS = (2010, 2014)                  # land cover averaged over pre-window years only (no look-ahead)
ARIDITY_CAP = 10.0

# ---------------------------------------------------------------- lab QC handling
# Lab flag bits (variable <var>_qc): 1 constant_field, 2 identical_to_previous_record, 4 no_data.
QC_CONSTANT, QC_IDENTICAL, QC_NODATA = 1, 2, 4
# A constant field is natural for rainfall on a dry day, never for these:
QC_CONSTANT_IS_BAD = ("srad", "tmmn", "tmmx", "vpd", "etr", "pet")

# Plausible physical ranges of the decoded values (fail loudly outside them).
RANGES = dict(pr=(0.0, 700.0), tmmn=(200.0, 320.0), tmmx=(210.0, 330.0), etr=(0.0, 40.0),
              vpd=(0.0, 12.0), srad=(0.0, 500.0), et_dekad=(0.0, 150.0), sm=(0.0, 0.7),
              index=(-2.5, 2.5))

# ---------------------------------------------------------------- drought definitions
# Each rule names its index, its accumulation window, and where the threshold
# comes from. `confirmed` records what the threshold was checked against.
DROUGHT_DEFS = {
    "mckee_moderate": dict(
        index="spi90d", window_days=90, op="le", threshold=-1.0,
        cite="McKee, Doesken & Kleist (1993), 8th Conf. Applied Climatology: "
             "SPI <= -1.0 is the 'moderate drought' class.",
        confirmed="recalled from the paper; not re-read here"),
    "usdm_d1": dict(
        index="spi90d", window_days=90, op="le", threshold=-0.8,
        cite="Svoboda et al. (2002), BAMS 83, doi:10.1175/1520-0477-83.8.1181; "
             "US Drought Monitor classification table: D1 = SPI/SPEI -0.8 to -1.29.",
        confirmed="droughtmonitor.unl.edu classification table, fetched 2026-10-07"),
    "spei_moderate": dict(
        index="spei90d", window_days=90, op="le", threshold=-1.0,
        cite="Vicente-Serrano, Begueria & Lopez-Moreno (2010), J. Climate, "
             "doi:10.1175/2009JCLI2909.1; same class scale as SPI.",
        confirmed="recalled from the paper; not re-read here"),
    "eddi_high": dict(
        index="eddi30d", window_days=30, op="ge", threshold=0.84,
        cite="Hobbins et al. (2016), J. Hydrometeor., doi:10.1175/JHM-D-15-0121.1: "
             "EDDI is positive for anomalously high evaporative demand (drier); "
             "0.84 is the 80th percentile of a standard normal.",
        confirmed="sign verified on the lab file (corr. with SPEI -0.80); "
                  "the 80th-percentile category bound is recalled, not re-read"),
}
# Flash drought (Ford & Labosier 2017, Agric. For. Meteorol., doi:10.1016/j.agrformet.2017.08.031;
# Otkin et al. 2018, BAMS, doi:10.1175/BAMS-D-17-0149.1): root-zone soil moisture falling from
# >= 40th to <= 20th percentile within 20 days. Thresholds recalled from the paper.
FLASH = dict(start_pct=40.0, end_pct=20.0, within_dekads=2)

# US Drought Monitor classes for SPI/SPEI (upper bound of each class, inclusive).
USDM_CUTS = (("D0", -0.5), ("D1", -0.8), ("D2", -1.3), ("D3", -1.6), ("D4", -2.0))
# Cut points stated in the lab file's own `category` attribute. They differ from USDM.
LAB_CATEGORY_CUTS = (-0.5, -0.7, -1.2, -1.5, -2.0)
WET_THRESHOLD = 0.8          # placebo periods: SPI-90 >= +0.8 (mirror of D1)
LAB_INDICES = ("spi30d", "spi90d", "spei30d", "spei90d", "eddi30d")
PENTAD_MAX_GAP_DAYS = 6      # lab pentads are 5 days apart, 6 across some year ends

# ---------------------------------------------------------------- model inputs
# Research Plan 6.1: F_hat(P, SMrz, VPD, Rn, T, PFT, season). Nothing else.
XGB_DYNAMIC = ("P_dek_anom", "SPI30_g", "SPI90_g", "VPD_anom", "SRAD_anom", "TMEAN_anom",
               "ETR_anom", "SMrz_anom", "SMrz_lag1_anom")
SEASON_FEATURES = ("doy_sin", "doy_cos")
PFT_GROUPS = dict(f_forest=(1, 2, 3, 4, 5), f_shrub=(6, 7), f_savanna=(8, 9), f_grass=(10,),
                  f_wetland=(11,), f_crop=(12,), f_urban=(13,), f_mosaic=(14,),
                  f_barren=(16,), f_water=(17,))
STATIC_FEATURES = ("PAW", "aridity", "P_annual", "T_annual") + tuple(PFT_GROUPS)
# Exposure and management variables. They stratify results; they never enter a model.
# (ET climatology is excluded too: a high summer normal identifies irrigated cells.)
EXPOSURE_VARS = ("WTD", "WTD_sd", "f_irr_major", "f_irr_minor", "f_irrigated", "f_rainfed_crop",
                 "f_crop_fragments",
                 "f_noncrop", "gir_low_freq", "gir_high_freq", "mgmt", "dWL_long", "dWL_recent",
                 "ET_clim")

# ---------------------------------------------------------------- split
N_FOLDS = 5
SEED = 42
EMBARGO_ACF = 0.2            # embargo = first lag (dekads) where the median ET-anomaly ACF < this
EMBARGO_RANGE = (2, 6)
SPACE_BLOCK_DEG = 0.5        # spatial blocks for the bootstrap
SPLIT_TRAIN, SPLIT_VAL, SPLIT_TEST, SPLIT_EMBARGO = 0, 1, 2, -1

# ---------------------------------------------------------------- uncertainty
ALPHA = 0.10                 # 90 % prediction intervals
N_BOOT_EFFECT = 500
SSEBOP_STEP_MM = 1.0         # SSEBop is stored as whole mm; rounding SD = step / sqrt(12)
SPECIFIC_YIELD = 0.15        # user constant (TBI 11.3)

# ---------------------------------------------------------------- XGBoost
XGB_PARAMS = dict(max_depth=6, learning_rate=0.05, subsample=0.8, colsample_bytree=0.7,
                  reg_lambda=1.0, min_child_weight=10, tree_method="hist", max_bin=256)
XGB_ROUNDS = 800
XGB_EARLY_STOP = 50
XGB_QUANTILES = (0.05, 0.5, 0.95)
XGB_BATCH_DEKADS = 12

# ---------------------------------------------------------------- differentiable LSTM-HBV
# (lo, hi) bounds; raw LSTM outputs pass through a sigmoid into these.
# Literature-typical, not calibrated to the Ogallala (unchanged from delta_gwb/config.py).
HBV_BOUNDS = dict(
    FC=(50.0, 600.0), LP=(0.3, 1.0), BETA=(1.0, 4.0), GAMMA=(0.5, 3.0), PERC=(0.1, 8.0),
    K0=(0.05, 0.9), K1=(0.01, 0.5), K2=(0.001, 0.2), UZL=(0.0, 100.0), TT=(-2.0, 2.0),
    DD=(1.0, 6.0), CWH=(0.0, 0.3), CFR=(0.0, 0.8), ROUT_A=(0.5, 6.0), ROUT_TAU=(0.5, 10.0),
)
DLSTM = dict(hidden=256, layers=2, dropout=0.0, win=730, warm=365, lr=1e-3, epochs=30,
             batch_cells=64, windows_per_step=4, val_cells=512, sm_weight=0.25,
             et_sd_floor=3.0,     # mm/dekad; keeps winter dekads from dominating the loss
             sm_pct=(1.0, 99.0),  # per-cell soil-moisture range used for relative saturation
             grad_clip=1.0)
DLSTM_STATIC = ("PAW", "aridity", "P_annual", "T_annual", "hsg_A", "hsg_B", "hsg_C",
                "hsg_D") + tuple(PFT_GROUPS)


def assert_no_exposure(names):
    """Guard used by both trainers: fail if an exposure variable is a model input."""
    bad = sorted(set(names) & set(EXPOSURE_VARS))
    if bad:
        raise ValueError(f"exposure/management variables used as model inputs: {bad}")
