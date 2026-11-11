"""USDA NASS Cropland Data Layer codes -> vegetation class groups used by P1.1 and P1.2.

Codes are NASS's published legend (the CDL notebook stores them as class_names on the variable).
Codes 61-65, 81-92 and 171/181 are older legend values that still occur in early years.
FUNCTIONAL (plan v3, P1.2): the class input of the model and the class of the depth ladder. The CDL group stays
a secondary ("species") input. Spring cereals, oilseeds other than sunflower-type summer crops, and the mixed
"other crop" group form "other_annual"; groups never eligible for the reference map to "none".
ELIGIBILITY: "natural" groups are always eligible for the rain-only reference; "crop" groups only where
the crop is mostly grown dryland in that county and year (P1.1 dryland rule); "never" groups never.
"""
import numpy as np

GROUPS = [  # (name, eligibility, codes)
    ("grass_pasture", "natural", [176, 62, 171, 181]),
    ("shrub", "natural", [152, 64]),
    ("winter_wheat", "crop", [24, 26, 225, 236, 238]),
    ("small_grains", "crop", [21, 22, 23, 25, 27, 28, 29, 205]),
    ("corn", "crop", [1, 12, 13]),
    ("sorghum", "crop", [4]),
    ("soybean", "crop", [5]),
    ("cotton", "crop", [2]),
    ("oilseed", "crop", [6, 31, 32, 33]),
    ("hay_alfalfa", "crop", [36, 37, 58, 59, 60]),
    ("fallow", "crop", [61]),
    ("other_crop", "crop", list(range(3, 4)) + list(range(10, 12)) + list(range(14, 21)) + list(range(30, 31))
     + list(range(34, 36)) + list(range(38, 58)) + list(range(66, 81)) + list(range(204, 205))
     + list(range(206, 225)) + list(range(226, 236)) + list(range(237, 238)) + list(range(239, 255))),
    ("forest", "never", [141, 142, 143, 63]),
    ("wetland", "never", [87, 190, 195]),
    ("water", "never", [111, 83]),
    ("developed", "never", [121, 122, 123, 124, 82]),
    ("barren", "never", [131, 65, 112]),
    ("other", "never", []),                      # everything else, incl. 0 background
]
NAMES = [g[0] for g in GROUPS]
ELIGIBILITY = {g[0]: g[1] for g in GROUPS}
N = len(GROUPS)
OTHER = NAMES.index("other")
CROP_GROUPS = [i for i, g in enumerate(GROUPS) if g[1] == "crop"]
NATURAL_GROUPS = [i for i, g in enumerate(GROUPS) if g[1] == "natural"]

LUT = np.full(256, OTHER, "uint8")
for _gi, (_name, _el, _codes) in enumerate(GROUPS):
    for _c in _codes:
        if LUT[_c] != OTHER:
            raise ValueError(f"CDL code {_c} assigned to two groups")
        LUT[_c] = _gi


FUNCTIONAL_OF = dict(grass_pasture="perennial_grass", shrub="shrub", winter_wheat="winter_cereal",
                     small_grains="other_annual", corn="summer_crop", sorghum="summer_crop", soybean="summer_crop",
                     cotton="summer_crop", oilseed="summer_crop", hay_alfalfa="hay", fallow="fallow",
                     other_crop="other_annual", forest="none", wetland="none", water="none", developed="none",
                     barren="none", other="none")
FUNC_NAMES = ["winter_cereal", "summer_crop", "other_annual", "perennial_grass", "shrub", "hay", "fallow", "none"]
NF = len(FUNC_NAMES)
FUNC_NONE = FUNC_NAMES.index("none")
FUNC_LUT = np.full(256, FUNC_NONE, "uint8")          # CDL group index (255 = no class) -> functional index
for _gi, _name in enumerate(NAMES):
    FUNC_LUT[_gi] = FUNC_NAMES.index(FUNCTIONAL_OF[_name])
WOODY_GROUPS = [NAMES.index("shrub"), NAMES.index("forest")]


def functional_of(group):
    """CDL group index array (255 = none) -> functional group index (FUNC_NONE where none)."""
    return FUNC_LUT[np.asarray(group).astype("int64").clip(0, 255)]


def group_of(codes):
    """uint8 CDL codes (255 = fill) -> group index; fill -> 255."""
    codes = np.asarray(codes)
    out = LUT[np.clip(codes, 0, 255).astype("uint8")]
    return np.where(codes == 255, 255, out).astype("uint8")
