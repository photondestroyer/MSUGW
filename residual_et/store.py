"""The prepared feature store: plain .npy arrays plus a JSON manifest.

Written once by features.py (system Python, needs scipy) and read by both
trainers. The LSTM environment only needs numpy to read it.
"""
import json
from pathlib import Path

import numpy as np

# Arrays every complete store must contain. features.py asserts this list on
# write and the trainers on read, so the schema cannot drift silently.
DYNAMIC = ("ET", "ET_sd", "ET_n", "ET_anom", "P_dek", "P_dek_anom", "P30", "P90", "SPI30_g",
           "SPI90_g", "SPI30_se", "SPI90_se", "VPD_anom", "SRAD_anom", "TMEAN_anom", "ETR_anom",
           "SMrz", "SMrz_lag1", "spi30d", "spi90d", "spei30d", "spei90d", "eddi30d")
CLIM = ("ET_clim", "ET_clim_sd", "ET_clim_se")                      # (36, N)
STATIC = ("PAW", "PAW_sd", "WTD", "WTD_sd", "aridity", "P_annual", "T_annual", "f_irrigated",
          "f_rainfed_crop", "f_crop_fragments", "f_noncrop", "f_irr_major", "f_irr_minor", "f_natural",
          "gir_low_freq", "gir_high_freq", "hsg_A", "hsg_B", "hsg_C", "hsg_D", "dWL_long",
          "dWL_recent", "f_forest", "f_shrub", "f_savanna", "f_grass", "f_wetland", "f_crop",
          "f_urban", "f_mosaic", "f_barren", "f_water")
DAILY = ("P_daily", "T_daily", "EP_daily", "SM_daily")              # (ndays, N)
DAILY_CLIM = ("P_dclim", "T_dclim", "EP_dclim")                     # (366, N)
AXES = ("dekad_start", "dek36", "day", "cell_lat", "cell_lon", "cell_iy", "cell_ix",
        "space_block")


class Store:
    def __init__(self, path):
        self.path = Path(path)

    # -- writing
    def create(self, name, shape, dtype="float32", fill=np.nan):
        self.path.mkdir(parents=True, exist_ok=True)
        a = np.lib.format.open_memmap(self.path / f"{name}.npy", mode="w+", dtype=dtype, shape=shape)
        if fill is not None:
            a[:] = fill
        return a

    def save(self, name, arr):
        self.path.mkdir(parents=True, exist_ok=True)
        np.save(self.path / f"{name}.npy", np.asarray(arr))

    def write_json(self, name, obj):
        self.path.mkdir(parents=True, exist_ok=True)
        (self.path / f"{name}.json").write_text(json.dumps(obj, indent=1, default=str))

    # -- reading
    def has(self, name):
        return (self.path / f"{name}.npy").exists()

    def arr(self, name, mmap=True):
        p = self.path / f"{name}.npy"
        if not p.exists():
            raise FileNotFoundError(f"store array missing: {p} (run residual_et.features first)")
        return np.load(p, mmap_mode="r" if mmap else None)

    def read_json(self, name):
        return json.loads((self.path / f"{name}.json").read_text())

    def done(self, stage):
        return (self.path / f"_done_{stage}").exists()

    def mark(self, stage):
        (self.path / f"_done_{stage}").write_text("ok")

    def missing(self):
        want = DYNAMIC + CLIM + STATIC + DAILY + DAILY_CLIM + AXES
        return [n for n in want if not self.has(n)]

    def require_complete(self):
        miss = self.missing()
        if miss:
            raise FileNotFoundError(f"store {self.path} is incomplete; missing {miss[:8]}"
                                    f"{' ...' if len(miss) > 8 else ''}")

    # -- convenience
    @property
    def dekad_start(self):
        return self.arr("dekad_start", mmap=False).astype("datetime64[D]")

    @property
    def days(self):
        return self.arr("day", mmap=False).astype("datetime64[D]")

    @property
    def n_cells(self):
        return int(self.arr("cell_lat").shape[0])
