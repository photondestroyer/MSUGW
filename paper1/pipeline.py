"""The stage graph of paper1 (Step 0 pieces -> P1.1 -> P1.2).

Each stage: name, deps (stages read), keys (config settings that change its result), modules (code that
computes it), func ("module:function"). The runner (core/stages.py) turns this into signatures, output
folders and resume.
"""
CORE = ["core/ncio.py", "core/paths.py"]

STAGES = [
    dict(name="tiles", deps=[], func="step0.tiles:run", modules=CORE + ["step0/tiles.py"],
         keys=["TEST_TILES", "TEST_WINDOW_PX", "TEST_WINDOW_CENTER", "FOCUS_POINTS", "PIXEL_M", "TILE_PX",
               "LATTICE_ORIGIN_M", "BLOCK_PX", "BLOCK_MIN_ROI_FRACTION", "LANDSAT_TAG", "SOURCE_SUBFOLDER"]),
    dict(name="blocks", deps=["tiles"], func="step0.blocks:run", modules=CORE + ["step0/blocks.py", "step0/tiles.py"],
         keys=["TEST_YEARS", "YEAR_FIRST", "YEAR_LAST", "AGG_YEAR_FIRST", "LANDSAT_TAG", "OPENET_TAG",
               "WATER_MIN_NOBS", "LANDSAT_FILL", "NOBS_FILL", "HIST_BINS", "HIST_EDGES", "OPENET_COUNT_VARS"]),
    dict(name="drivers", deps=["tiles"], func="step0.drivers:run", modules=CORE + ["step0/drivers.py"],
         keys=["DRIVER_VARS", "DRIVER_MIN_DAY_FRACTION", "RAIN_WINDOWS_D", "ARIDITY_BASE_YEARS", "QC_CONSTANT_IS_BAD",
               "BUCKET_CAPACITY_BINS_MM", "BUCKET_LP", "BUCKET_SPINUP_FROM", "DRIVER_MONTH_FIRST", "SPI_VAR"]),
    dict(name="counties", deps=[], func="step0.counties:run", modules=CORE + ["step0/counties.py"],
         keys=["COUNTIES_URL", "AQUIFER_STATE_FIPS"]),
    dict(name="hand", deps=["tiles"], func="step0.hand:run", modules=CORE + ["step0/hand.py"],
         keys=["HU4_CODES", "STREAM_FCODES", "INCLUDE_EPHEMERAL", "MIN_DRAINAGE_KM2", "WATERBODY_FTYPES",
               "BURN_DEPTH_M", "HU8_BUFFER_M", "DEM_COG_URL", "NHD_GDB_URL", "CANAL_FTYPES"]),
    dict(name="static", deps=["tiles", "counties", "hand"], func="step0.static:run",
         modules=CORE + ["step0/static.py", "step0/hand.py", "step0/tiles.py"],
         keys=["STATIC_SOIL_VARS", "WTD_SHALLOW_M", "IRR_BUFFER_M", "BUCKET_CAPACITY_BINS_MM", "HAND_MIN_M", "NWI_TAG",
               "LANID_FREQ_MIN_YEARS", "CANAL_BUFFER_M"]),
    dict(name="classes", deps=["tiles", "static"], func="step0.classes:run",
         modules=CORE + ["step0/classes.py", "cdl_groups.py"], keys=["CDL_TAG", "TEST_YEARS", "YEAR_FIRST", "YEAR_LAST"]),
    dict(name="dryland", deps=["tiles", "static", "classes"], func="p1_1.dryland:run",
         modules=["p1_1/dryland.py", "cdl_groups.py"], keys=["DRYLAND_SHARE_MIN", "DRYLAND_MIN_PIXELS"]),
    dict(name="reference", deps=["tiles", "static", "blocks", "classes", "dryland", "drivers"], func="p1_1.build:run",
         modules=CORE + ["p1_1/build.py", "p1_1/rules.py", "cdl_groups.py"],
         keys=["WTD_MIN_M", "WTD_SENSITIVITY_M", "WTD_MIN_COVERAGE", "HAND_MIN_M", "HAND_PASS_FRACTION", "IRR_BUFFER_M",
               "IRRIGATION_REGISTRIES", "WETLAND_NOISE_PX", "NWI_MISSING_POLICY", "WATER_NOISE_PX", "WATER_MAX_MONTHS",
               "GDE_SENSITIVITY_MAX", "CLASS_PURITY_MIN", "MIN_VALID_MONTHS", "TARGET_MIN_PIXEL_FRACTION", "TEST_NWI_MISSING_POLICY",
               "DEPTH_RUNGS_M", "DEPTH_BY_FUNCTIONAL", "WOODY_DEPTH_M", "WOODY_NOISE_PX", "DEPTH_ACCEPTED",
               "LANID_FREQ_MIN_YEARS", "CANAL_BUFFER_M", "FIRE_PERIMETERS"]),
    dict(name="verify_p11", deps=["tiles", "reference", "blocks", "static", "drivers"], func="p1_1.verify:run",
         modules=["p1_1/verify.py", "p1_1/rules.py"], gate="G1a",
         keys=["MIN_REF_BLOCKS_PER_STRATUM", "V2_TOLERANCE", "V2_PASS_SHARE", "V4_SMD_MAX", "V4_VARS",
               "P10_VALID_SHARE_MIN", "DEPTH_ACCEPTED", "DEPTH_RUNGS_M", "IRR_BUFFER_M", "TEST_NWI_MISSING_POLICY",
               "NWI_MISSING_POLICY"]),
    dict(name="samples", deps=["tiles", "static", "reference", "blocks", "drivers", "verify_p11"], gated_by=["verify_p11"],
         func="p1_2.blockset:run",
         modules=["p1_2/blockset.py"], keys=["SAMPLE_BLOCKS_PER_STRATUM", "TEST_SAMPLE_BLOCKS_PER_STRATUM", "TARGETS",
                                             "STATIC_INPUTS", "MEAS_INPUTS"]),
    dict(name="folds", deps=["tiles", "drivers"], func="p1_2.folds:run", modules=["p1_2/folds.py"],
         keys=["SPACE_GROUP_DEG", "SPACE_FOLDS", "YEAR_GROUPS", "EMBARGO_YEARS", "INNER_SPLIT", "TEST_SPACE_GROUP_DEG",
               "TEST_SPACE_FOLDS", "TEST_YEAR_GROUPS"]),
    dict(name="train", deps=["samples", "drivers", "folds"], func="p1_2.train:run",
         modules=["p1_2/train.py", "p1_2/m1.py", "p1_2/features.py", "p1_2/calendars.py", "p1_2/baselines.py",
                  "cdl_groups.py"],
         keys=["DYN_INPUTS", "STATIC_INPUTS", "EXPOSURE_VARS", "PRESEASON_MONTHS", "QUANTILES", "ALPHA", "M1", "TEST_M1",
               "M2_WINDOWS_MONTHS", "M2_MIN_SAMPLES", "M0_MIN_YEARS", "CAL_YEAR_HOLDOUT", "CAL_STRATA", "PLANT",
               "RUN_ABLATION", "V2_TOLERANCE", "MEAS_INPUTS", "SPECIES_INPUT", "CAL_MERGE_ORDER", "CAL_MIN_N",
               "TEST_CAL_MIN_N", "RUN_MEAS_ABLATION", "RUN_RATIO_MODEL", "CLASS_FLIP", "CLASS_CONFUSION", "OBS_REGIMES",
               "DROUGHT_DRY_SPI", "DROUGHT_WET_SPI", "PASS"]),
    dict(name="predict", deps=["tiles", "static", "reference", "blocks", "drivers", "folds", "train"],
         func="p1_2.predict:run", modules=["p1_2/predict.py", "p1_2/features.py", "p1_2/baselines.py"], keys=[]),
    dict(name="verify_p12", deps=["tiles", "folds", "train", "reference", "blocks", "predict", "static", "drivers"],
         func="p1_2.verify:run", gate="G1b", modules=["p1_2/verify.py", "p1_1/rules.py"],
         keys=["PASS", "TEST_PASS_STRATUM_MIN_N", "OBS_REGIMES", "DROUGHT_TEST_YEARS", "DROUGHT_DRY_SPI",
               "DROUGHT_WET_SPI", "DRIFT_YEARS", "DEPTH_BY_FUNCTIONAL", "WTD_MIN_M", "MEAS_INPUTS"]),
]

P11 = ["verify_p11"]
P12 = ["verify_p12"]


def runner():
    from .core.stages import Runner
    return Runner(STAGES)
