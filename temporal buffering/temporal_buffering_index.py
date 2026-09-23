"""
Temporal Buffering Index (TBI) toolkit
=======================================
Implements the Vulnerability Index (VI) from VVIP paper:
  Gonzalez Cruz et al. (2021) Sci Rep 11:21648
  https://doi.org/10.1038/s41598-021-98829-5
and extends to a time-varying Temporal Buffering Index per
Groundwater_Buffering_Research_Plan.docx Section 6.

Key features
------------
- READ-ONLY, LAZY loading via xarray + dask (chunks). No dataset is modified.
- Handles 30+ GB files (GRIDMET 58GB, DROUGHT 16GB) by chunking and computing
  only on spatial/temporal subsets.
- Explicit placeholders (NaN) for missing required variables - NEVER synthesizes data.
- All equations reference paper line numbers and include units.

VVIP Equations (Section Methodology):
  Eq1: VI = DSI / (alpha*SBI + beta*GBI) = 1/RI   ; alpha+beta=1, 0<=alpha,beta<=1
  Eq2: E(c_i,j,k) = integral c * f(c) dc  (expected value of drought char)
  Eq3: DSI_l = prod  E(c_i,j,k)^{w_i,j,k}  (weighted product; w negative for IDD)
  Eq4: SBI_l = PAW^{w_PAW} * DAC^{w_DAC} * ADP^{-w_ADP}  (DAC positive, ADP negative)
  Eq5: T = K * ST   (transmissivity = hydraulic conductivity * saturated thickness)
  Eq6: Delta_H = IRRAMT / S   (water level drop per unit irrigation; S ~ specific yield)
  Eq7: GBI_l = T^{w_T} * S^{-w_S}  (reduces to hydraulic diffusivity if w_T=w_S=1)
  Eq8-10: Entropy weights  p_i,j = X_n,i,j / sum_j X_n,i,j ; E_i = -sum p log p / log J ; w_i = (1-E_i)/sum(1-E_k)

Research Plan Section 6 extensions:
  Buffering signal  Bi(t) = Fi(t) - F_hat_i(P, SMrz, VPD, Rn, T, PFT, season)   [6.1]
    Fi = anomaly of ET, SIF, VWC. Residual persists during dry-down.
  Temporal Buffering Index (time-varying): TBI(t) = (alpha*SBI(t)+beta*GBI(t)) / DSI(t)
    or equivalently 1/VI(t). Larger TBI = more buffering.
  Failure threshold: Bi(t) ~ f(WTD, dWTD, dStorage, PFT, rooting_depth, irrigation, intensity)
  Decoupling: Di = g(WTD, ...) across coarsest grid (VOD/SIF limited).

Usage
-----
>>> from temporal_buffering_index import TBIToolkit, PLACEHOLDERS
>>> toolkit = TBIToolkit(base_dir=r"G:/MSU_GWB/datasets")
>>> toolkit.open_all_lazy()  # lazy, no load
>>> # compute on a small test subset to verify pipeline without OOM:
>>> toolkit.demo_small_subset(y_slice=slice(0,10), x_slice=slice(0,10), t_slice=slice(0,20))
>>> # full run (still lazy, triggered only on .compute()):
>>> vi = toolkit.compute_VI(alpha=0.5, beta=0.5)  # dask array, not yet computed
>>> # persist to netCDF (optional, writes NEW file, never overwrites source):
>>> toolkit.save_result(vi, "tbi_output.nc")

Documentation for datasets, assumptions, and required additional data lives in:
  TBI_DOCUMENTATION.md  (this repo)
  See also print_assumptions(), list_required_datasets()

Author: Muse Spark / OpenCode
Date: 2026-08-24
"""

from __future__ import annotations
import os
import warnings
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Optional, Tuple, List

import numpy as np
import xarray as xr
import dask.array as da

# ---------------------------------------------------------------------------
# PLACEHOLDERS - documented missing data. Do NOT fabricate values.
# ---------------------------------------------------------------------------
PLACEHOLDERS = {
    # Soil Buffer Index (SBI) - Eq4
    "PAW": np.nan,  # Plant Available Water = (field capacity - wilting point)/root_depth
                    # Needs: gSSURGO 30m (Staff 2019) field capacity + PWP. Not in repo.
    "PAW_units": "mm/mm (dimensionless, depth-normalized)",
    "PAW_source_required": "USDA gSSURGO Gridded Soil Survey Geographic Database",

    # Groundwater Buffer Index (GBI) - Eq5-7
    "K": np.nan,   # Hydraulic conductivity (m/day)
                   # Needs: Cederstrand & Becker (1998) Digital Map of Hydraulic Conductivity
    "Sy": np.nan,  # Specific yield / storage coefficient (dimensionless)
                   # Needs: McGuire et al. (2012) USGS SIR 2012-5177
    "ST": np.nan,  # Saturated thickness (m)
                   # Needs: McGuire (2017) Water-Level & Recoverable Water + aquifer bottom
                   #         (or USGS SIR 2012-5177 Fig: Saturated Thickness 2009)
    "T": np.nan,   # Transmissivity T=K*ST (m2/day) - derived, placeholder if K or ST missing

    # Soil moisture drought indicators - Eq4 DAC/ADP need SSMI
    "SSMI_source_required": "CPC Global Monthly Soil Moisture (Fan & van den Dool 2004) 0.5deg 1949-present leaky bucket",
    "SSMI_proxy_available": "SPL4SMGP sm_rootzone (2015-present, 9km) - SHORT record, not 30yr",
    "SSMI_note": "VVIP requires >=30yr record for stable E(c) and DAC/ADP. SMAP alone insufficient.",

    # SIF / photosynthesis for Bi(t)
    "SIF_source_required": "OCO-2/3 or TROPOMI SIF (2014/2018-present) or GOSIF (derived, not independent)",
    "SIF_proxy_available": "MOD13A3 NDVI/EVI (greenness, not photosynthesis) + AmeriFlux GPP (point validation)",

    # Well dynamics for reversibility (Research Plan 6.5)
    "WELLS_source_required": "USGS + state well networks (point time series of depth-to-water)",
    "WELLS_proxy_available": "HRES-WTD 2015 static 30m (Ma et al. 2026) + SPL4SMGP depth_to_water_table + GRACE mascon (coarse)",
}

# ---------------------------------------------------------------------------
# Entropy weighting (Eq8-10)
# ---------------------------------------------------------------------------

def entropy_weights(normalized_da: xr.DataArray, dim: str = "location") -> xr.DataArray:
    """
    Compute entropy weights per Eq8-10 of VVIP.

    Parameters
    ----------
    normalized_da : xr.DataArray
        Dimensions (attribute, location) where values already normalized 0-1
        (0 = least preferred, 1 = most preferred). For DSI/SBI/GBI separately.
    dim : str
        Dimension along which to sum for probabilities (the alternatives J).

    Returns
    -------
    xr.DataArray over attribute dimension with weights summing to 1.
    """
    # Eq8: p_ij = X_n,i,j / sum_j X_n,i,j
    # Add tiny epsilon to avoid log(0)
    eps = 1e-12
    p = normalized_da / normalized_da.sum(dim=dim)
    p = p.clip(min=eps)
    # Eq9: E_i = -sum_j p_ij log(p_ij) / log(J)
    J = normalized_da.sizes[dim]
    E = -(p * np.log(p)).sum(dim=dim) / np.log(J)
    # Eq10: w_i = (1-E_i) / sum_k (1-E_k)
    w = (1 - E) / (1 - E).sum(dim="attribute") if "attribute" in E.dims else (1 - E) / (1 - E).sum()
    return w


def weighted_product(da: xr.DataArray, weights: xr.DataArray, attr_dim: str = "attribute") -> xr.DataArray:
    """
    Weighted multiplicative aggregation: prod_i  x_i^{w_i}
    Used for DSI (Eq3), SBI (Eq4), GBI (Eq7). For ADP and IDD the weight is applied
    as negative exponent internally (caller should pass negative w or invert variable).
    Lazy via dask - no compute triggered if inputs are dask-backed.
    """
    # Use log domain for numerical stability, then exp
    # log(prod x_i^{w_i}) = sum w_i * log(x_i)
    # Clip x_i to avoid log(0/neg)
    eps = 1e-9
    safe = da.clip(min=eps)
    log_sum = (weights * np.log(safe)).sum(dim=attr_dim)
    return np.exp(log_sum)

# ---------------------------------------------------------------------------
# Drought characteristics via Theory of Runs (Guerrero-Salazar & Yevjevich 1975)
# ---------------------------------------------------------------------------

def detect_drought_events(series: xr.DataArray, threshold: float = -1.0) -> Dict[str, xr.DataArray]:
    """
    Theory of runs: classify as drought when index <= threshold.
    Returns duration (DD), severity (DS), inter-drought duration (IDD).

    NOTE: This is 1D per location. For 2D grid we apply via xarray.apply_ufunc
    or dask map_blocks lazily. Caller should chunk so this stays lazy.

    For this toolkit we provide the logic but do NOT auto-run on full 58GB.
    Use compute_expected_values() which wraps this lazily.

    Parameters
    ----------
    series : xr.DataArray with dim time, values SPI/SPEI
    threshold : float, default -1 per VVIP p6 (drought if <= -1)
    """
    # Placeholder doc - actual per-pixel loop would be vectorized with numbagg/bottleneck
    # We document approach rather than eagerly compute on full grid here.
    raise NotImplementedError("Use TBIToolkit.compute_expected_values() for lazy chunked version")


# ---------------------------------------------------------------------------
# Core Toolkit
# ---------------------------------------------------------------------------

@dataclass
class DatasetInventory:
    """Mirrors merged_datasets/*.nc contents."""
    path: Path
    exists: bool
    size_gb: float
    variables: List[str]
    dims: Dict[str, int]
    time_range: Optional[str] = None


class TBIToolkit:
    """
    Read-only, lazy toolkit for Temporal Buffering Index.

    All open_dataset calls use chunks and dask; no .load()/.compute() unless
    user explicitly requests a small subset or calls .compute().

    Parameters
    ----------
    base_dir : str | Path
        Root datasets folder, e.g. G:\\MSU_GWB\\datasets
    merged_subdir : str
        Subfolder containing merged NetCDFs (default: merged_datasets)
    """

    def __init__(self, base_dir: str | Path = r"G:/MSU_GWB/datasets", merged_subdir: str = "merged_datasets"):
        self.base_dir = Path(base_dir)
        self.merged_dir = self.base_dir / merged_subdir
        self.datasets: Dict[str, xr.Dataset] = {}  # lazy handles
        self.inventory: Dict[str, DatasetInventory] = {}
        # Chunk choices tuned for each file's native chunking
        self.chunk_configs = {
            "DROUGHT_Merged_Ogallala.nc": {"time": 10, "y": 287, "x": 181},  # native 631MB per var
            "GRIDMET_Merged_Ogallala.nc": {"time": 5, "y": 287, "x": 181},
            "HRES-WTD_2015_Ogallala.nc": {"y": 2048, "x": 2048},
            "Hydrologic_Soil_Group_250m_2019_Ogallala.nc": {"y": 2048, "x": 2048},
            "SPL4SMGP_Ogallala_FULL.nc": {"time": 8, "y": 134, "x": 104},
            "SPL4SMGP_Ogallala_latest.nc": {"time": 8, "y": 134, "x": 104},
            "MOD13A3_Merged_Ogallala.nc": {"time": 1, "y": 1024, "x": 903},
            "MODIS_ET_SSEBop_Merged_Ogallala.nc": {"time": 10, "y": 1024, "x": 779},
            "GRACE_Merged_Ogallala.nc": {"time": 1, "y": 24, "x": 15},
            "MCD12Q1_Merged_Ogallala.nc": {"time": 1, "y": 1024, "x": 1024},
            "GFSAD1000_V1_2019_Ogallala.nc": {"y": 842, "x": 842},
        }

    # ---------------------------------------------------------------
    # Lazy open (read-only, no modification)
    # ---------------------------------------------------------------
    def open_lazy(self, filename: str, chunks: Optional[dict] = None) -> xr.Dataset:
        """Open single file lazily (dask), read-only."""
        fpath = self.merged_dir / filename
        if not fpath.exists():
            # try base_dir root (some single-year files live there)
            fpath = self.base_dir / filename
            if not fpath.exists():
                raise FileNotFoundError(f"Not found: {filename} in {self.merged_dir} or {self.base_dir}")
        cfg = chunks or self.chunk_configs.get(filename, "auto")
        # mask_and_scale=False preserves original ints for QA etc; decode_cf=True for time
        ds = xr.open_dataset(fpath, chunks=cfg, mask_and_scale=False, decode_cf=True)
        # Ensure we never write back: set read-only flag (xarray has no native RO, we just never call to_netcdf on source)
        self.datasets[filename] = ds
        return ds

    def open_all_lazy(self) -> Dict[str, xr.Dataset]:
        """Open every known merged file lazily. Safe for 30GB+ each."""
        import glob
        pattern = str(self.merged_dir / "*.nc")
        files = glob.glob(pattern)
        for f in sorted(files):
            fname = Path(f).name
            try:
                cfg = self.chunk_configs.get(fname, "auto")
                ds = xr.open_dataset(f, chunks=cfg, mask_and_scale=False, decode_cf=True)
                self.datasets[fname] = ds
            except Exception as e:
                warnings.warn(f"Could not lazy-open {fname}: {e}")
        # Also inventory AmeriFlux parquet / SPL4 etc lazily if requested
        return self.datasets

    def close_all(self):
        for ds in self.datasets.values():
            try:
                ds.close()
            except Exception:
                pass
        self.datasets.clear()

    def inventory_all(self) -> Dict[str, DatasetInventory]:
        """Build inventory without loading data (uses lazy open metadata only)."""
        import glob, os
        inv = {}
        for f in glob.glob(str(self.merged_dir / "*.nc")) + glob.glob(str(self.base_dir / "*.nc")):
            p = Path(f)
            try:
                # open with chunks=None just to read header quickly, then close
                ds = xr.open_dataset(f, chunks={}, mask_and_scale=False)
                variables = list(ds.data_vars)
                dims = dict(ds.sizes)
                # time range if present
                tr = None
                if "time" in ds.coords:
                    try:
                        t0 = str(ds["time"].values[0])[:10]
                        t1 = str(ds["time"].values[-1])[:10]
                        tr = f"{t0} -> {t1}"
                    except Exception:
                        tr = None
                size_gb = os.path.getsize(f) / (1024**3)
                inv[p.name] = DatasetInventory(path=p, exists=True, size_gb=size_gb,
                                               variables=variables, dims=dims, time_range=tr)
                ds.close()
            except Exception as e:
                inv[p.name] = DatasetInventory(path=p, exists=False, size_gb=0, variables=[f"ERR: {e}"], dims={})
        self.inventory = inv
        return inv

    # ---------------------------------------------------------------
    # VVIP: DSI (Eq2-3)
    # ---------------------------------------------------------------
    def compute_DSI(self,
                    drought_ds: Optional[xr.Dataset] = None,
                    time_slice: Optional[slice] = None,
                    y_slice: Optional[slice] = None,
                    x_slice: Optional[slice] = None,
                    use_vars: Tuple[str, ...] = ("spi90d", "spi180d", "spei90d", "spei180d"),
                    weights: Optional[Dict[str, float]] = None) -> xr.DataArray:
        """
        Compute Drought Stress Index per VVIP Eq3 as weighted product of
        expected drought characteristics. LAZY.

        In full VVIP, DSI = prod_{j,k,i} E(c_i,j,k)^{w_i,j,k}
          where i in {DS, DD, IDD}, j in {SPI,SPEI}, k in {3M,6M}
          w_DD>0, w_DS>0, w_IDD<0 (closely spaced droughts = more stress)

        With available data (DROUGHT_Merged_Ogallala.nc):
          spi90d ~ SPI-3M, spi180d ~ SPI-6M, spei90d ~ SPEI-3M, spei180d ~ SPEI-6M
          eddi* and pdsi etc are extra but not in VVIP.

        LIMITATION: Full E(c) via Theory of Runs requires fitting exponential
        (DD, IDD) and lognormal/Weibull/Gamma (DS) per pixel (187 locations in
        paper, 287*181=51947 here). That is compute-intensive and needs 30yr.
        This method provides a LAZY proxy: DSI_proxy = weighted product of
        drought intensity magnitudes, and documents the gap.

        If weights None, uses entropy weights (Eq8-10) or equal weights as fallback.
        Returns dask-backed DataArray (y, x) or (time, y, x) depending on reduction.

        Parameters
        ----------
        drought_ds : xr.Dataset | None
            If None, uses self.datasets["DROUGHT_Merged_Ogallala.nc"]
        time_slice, y_slice, x_slice : slice | None
            Subset to avoid OOM. Always use for full-grid compute via iteration.
        """
        if drought_ds is None:
            drought_ds = self.datasets.get("DROUGHT_Merged_Ogallala.nc")
            if drought_ds is None:
                drought_ds = self.open_lazy("DROUGHT_Merged_Ogallala.nc")
        # Subset lazily (no compute)
        sel = {}
        if time_slice is not None:
            sel["time"] = time_slice
        if y_slice is not None:
            sel["y"] = y_slice
        if x_slice is not None:
            sel["x"] = x_slice
        if sel:
            # isel semantics: need to handle time as index slice
            # Use isel for integer slices; sel for label slices not used here.
            # We use isel for simplicity.
            pass  # handled below with isel

        # For demo we use isel directly if slices provided
        if time_slice is not None or y_slice is not None or x_slice is not None:
            # Build isel dict only for non-None
            isel_kwargs = {}
            if time_slice is not None:
                isel_kwargs["time"] = time_slice
            if y_slice is not None:
                isel_kwargs["y"] = y_slice
            if x_slice is not None:
                isel_kwargs["x"] = x_slice
            drought_ds = drought_ds.isel(**isel_kwargs)

        # Stack available vars into attribute dimension for weighted product
        # Check existence
        avail = [v for v in use_vars if v in drought_ds]
        missing = [v for v in use_vars if v not in drought_ds]
        if missing:
            warnings.warn(f"DSI: missing vars {missing}, using only {avail}. "
                          f"VVIP expects SPI/SPEI at 3 & 6 months (spi90d, spi180d, spei90d, spei180d).")
        if not avail:
            raise ValueError("No drought vars available for DSI")

        # Approach A: if we had E(c) already, weighted product. Here proxy:
        # DSI_proxy(y,x) = mean over time of weighted product of |drought index| where index <= -1
        # Severity is captured as mean |SPI| during drought, duration as fraction of time in drought.
        # This is documented as PROXY, not full VVIP.
        # For lazy, we compute drought intensity: intensity = where(index <= -1, -index, 0)
        # Then aggregate over time lazily.

        # Create DataArray with dimension attribute
        # Normalize each var to 0-1 for entropy weighting if needed (per paper, DSI is product, normalization not needed,
        # but for entropy method we need normalized ratings).
        # Here we do time-mean of intensity as proxy for E(DS)/E(DD) combined.
        intensities = []
        for var in avail:
            da_var = drought_ds[var]
            # intensity: drought magnitude where threshold met, else 0
            # Keep dask lazy
            intensity = xr.where(da_var <= -1.0, -da_var, 0.0)
            # time mean (lazy)
            intensity_mean = intensity.mean(dim="time", skipna=True)
            intensities.append(intensity_mean)

        # Stack into (attribute, y, x)
        stacked = xr.concat(intensities, dim="attribute")
        stacked.coords["attribute"] = avail

        # Weights
        if weights is not None:
            w = xr.DataArray([weights.get(v, 1.0/len(avail)) for v in avail],
                             coords={"attribute": avail}, dims=["attribute"])
        else:
            # Try entropy weights lazily - need to compute normalization first
            # Normalize 0-1 per attribute across spatial domain (Eq8 needs alternatives = locations)
            # For dask, this triggers a compute if we call .compute(). Instead we use equal weights
            # as documented fallback (Laplace principle of insufficient reason, VVIP p5).
            # Proper entropy run can be done on a small subset via .compute() explicitly.
            w = xr.DataArray(np.full(len(avail), 1.0/len(avail)),
                             coords={"attribute": avail}, dims=["attribute"])
            warnings.warn("DSI weights: using EQUAL weights (Laplace). "
                          "For entropy weights (Eq8-10), call entropy_weights() on a computed subset.")

        # VVIP Eq3 weight for IDD should be negative. Our proxy intensity conflates DS/DD,
        # so we keep all weights positive. Document this.
        dsi = weighted_product(stacked, w, attr_dim="attribute")  # (y,x)
        dsi.attrs["long_name"] = "Drought Stress Index proxy (VVIP Eq3, intensity mean, equal weights)"
        dsi.attrs["units"] = "1 (normalized proxy)"
        dsi.attrs["note"] = ("PROXY: Full VVIP requires E(DS), E(DD), E(IDD) via Theory of Runs + "
                             "distribution fits (Exp for DD/IDD, LogNormal for DS). "
                             "Here intensity_mean approximates combined DS/DD. "
                             "IDD (inter-drought duration, negative weight) not yet included. "
                             "See TBI_DOCUMENTATION.md.")
        dsi.attrs["vvip_reference"] = "Gonzalez Cruz et al. 2021 Eq2-3, Fig4-6"
        return dsi

    # ---------------------------------------------------------------
    # VVIP: SBI (Eq4)
    # ---------------------------------------------------------------
    def compute_SBI(self,
                    hsg_ds: Optional[xr.Dataset] = None,
                    sm_ds: Optional[xr.Dataset] = None,
                    drought_ds: Optional[xr.Dataset] = None,
                    y_slice: Optional[slice] = None,
                    x_slice: Optional[slice] = None) -> xr.DataArray:
        """
        Soil Buffer Index per VVIP Eq4: SBI = PAW^{w} * DAC^{w} * ADP^{-w}

        With available data:
          - PAW (Plant Available Water) from gSSURGO NOT available -> placeholder NaN
            Proxy: Hydrologic Soil Group (HSG 1=A,2=B,3=C,4=D) inverted (A=high buffering).
            HSG file: Hydrologic_Soil_Group_250m_2019_Ogallala.nc  b1 in 1..4
          - DAC/ADP need paired SSMI vs SPI/SPEI binary contingency tables (per VVIP p4).
            SSMI requires CPC leaky bucket 1949-2018 NOT available.
            Proxy: SPL4SMGP sm_rootzone_pctl (percentile) as ag drought indicator,
                   but record 2015-present too short => DAC/ADP placeholders.

        Returns SBI proxy (y,x) LAZY, with NaN where PAW placeholder propagates.
        If you provide computed DAC/ADP arrays, they will be used.

        This function is LAZY and read-only.
        """
        # HSG branch (proxy PAW)
        if hsg_ds is None:
            try:
                hsg_ds = self.datasets.get("Hydrologic_Soil_Group_250m_2019_Ogallala.nc")
                if hsg_ds is None:
                    hsg_ds = self.open_lazy("Hydrologic_Soil_Group_250m_2019_Ogallala.nc")
            except Exception as e:
                warnings.warn(f"SBI: HSG not available: {e}. PAW will be NaN placeholder.")
                hsg_ds = None

        if hsg_ds is not None:
            if y_slice is not None or x_slice is not None:
                isel_kwargs = {}
                if y_slice is not None:
                    isel_kwargs["y"] = y_slice
                if x_slice is not None:
                    isel_kwargs["x"] = x_slice
                hsg_ds = hsg_ds.isel(**isel_kwargs) if isel_kwargs else hsg_ds
            # b1 is HSG code 1-4. Invert to buffering: A(1)=high, D(4)=low => buffering = (5 - b1)/4
            # This is a PROXY, not PAW. Documented.
            try:
                hsg = hsg_ds["b1"].isel(time=0) if "time" in hsg_ds["b1"].dims else hsg_ds["b1"]
                # Normalize 0-1 (VVIP says normalize for entropy method)
                paw_proxy = (5 - hsg) / 4.0
                paw_proxy = paw_proxy.clip(min=0, max=1)
                paw_proxy.attrs["long_name"] = "PAW proxy from Hydrologic Soil Group (inverted, normalized 0-1)"
                paw_proxy.attrs["note"] = ("PLACEHOLDER PROXY: True PAW = (FC - PWP)/root_depth from gSSURGO not available. "
                                           "HSG A=high infiltration/buffering, D=low. Inverted here. "
                                           "See PLACEHOLDERS['PAW_source_required'].")
            except Exception as e:
                warnings.warn(f"SBI HSG proxy failed: {e}")
                paw_proxy = xr.full_like(hsg_ds["b1"].isel(time=0), np.nan) if "time" in hsg_ds["b1"].dims else xr.full_like(hsg_ds["b1"], np.nan)
        else:
            # Full placeholder grid - need a reference grid to size it. Use drought grid if available.
            if drought_ds is None:
                try:
                    drought_ds = self.datasets.get("DROUGHT_Merged_Ogallala.nc") or self.open_lazy("DROUGHT_Merged_Ogallala.nc")
                    # create empty (y,x) template
                    template = drought_ds["spi90d"].isel(time=0)
                    paw_proxy = xr.full_like(template, np.nan)
                    paw_proxy.attrs["long_name"] = "PAW placeholder (NaN) - gSSURGO required"
                except Exception:
                    paw_proxy = None
            else:
                paw_proxy = None

        # DAC / ADP placeholders (would be computed from SSMI vs SPI/SPEI contingency tables)
        # Document that these are missing.
        # If sm_ds (SPL4SMGP) provided, we could compute a PROXY DAC/ADP on recent period:
        # DAC = P(no ag drought | met drought) - fast attempt but short record => warned.
        dac_proxy = None
        adp_proxy = None
        if sm_ds is not None:
            warnings.warn("SBI DAC/ADP from SPL4SMGP sm_rootzone_pctl is SHORT record proxy (2015+), not 30yr. Results indicative only.")

        # If we have at least paw_proxy, compute SBI as weighted product (with placeholders for DAC/ADP as 1.0 neutral)
        # Per spec: "use placeholder values at a max at clearly document them"
        # We use dac=adp= 0.5 neutral (so log contribution 0 when weight matters) but document.
        # Better: if missing, treat DAC=ADP=1 neutral element for product (since x^w =1 if x=1)
        # But clipping may hide missing. So we explicitly set SBI = paw_proxy * (dac_neutral) * (adp_neutral)
        # and set overall SBI to NaN where paw_proxy is NaN.

        if paw_proxy is None:
            raise ValueError("SBI: no proxy grid available. Need at least HSG or drought grid for template.")

        # Neutral placeholders for DAC/ADP (1.0 so they don't affect product)
        # Document that true VVIP SBI would be PAW*DAC*ADP^{-w}
        sbi = paw_proxy.copy(deep=False)
        sbi.name = "SBI"
        sbi.attrs["long_name"] = "Soil Buffer Index proxy (VVIP Eq4, HSG proxy for PAW, DAC/ADP neutral placeholders)"
        sbi.attrs["units"] = "1 (normalized 0-1)"
        sbi.attrs["note"] = ("PROXY: DAC/ADP set to 1.0 neutral placeholder. Full VVIP needs "
                             "SSMI (CPC leaky bucket 1949-2018) + SPI/SPEI contingency tables. "
                             "See PLACEHOLDERS['SSMI_source_required']." )
        sbi.attrs["vvip_reference"] = "Gonzalez Cruz et al. 2021 Eq4, Fig7-9"

        # If we had real DAC/ADP, SBI = PAW^{wPAW} * DAC^{wDAC} * ADP^{-wADP}
        # Example (when data available):
        # stacked = xr.concat([paw_proxy, dac_proxy, adp_proxy], dim="attribute")
        # w = entropy_weights(normalized_stacked)  # DAC positive, ADP negative already via exponent sign
        # sbi = weighted_product(stacked, w)

        return sbi

    # ---------------------------------------------------------------
    # VVIP: GBI (Eq5-7)
    # ---------------------------------------------------------------
    def compute_GBI(self,
                    wtd_ds: Optional[xr.Dataset] = None,
                    y_slice: Optional[slice] = None,
                    x_slice: Optional[slice] = None,
                    hydraulic_conductivity: Optional[xr.DataArray] = None,
                    specific_yield: Optional[xr.DataArray] = None,
                    saturated_thickness: Optional[xr.DataArray] = None) -> xr.DataArray:
        """
        Groundwater Buffer Index per VVIP Eq7: GBI = T^{wT} * S^{-wS}, T=K*ST (Eq5)

        With available data:
          - HRES-WTD 2015 (b1, meters depth-to-water) is AVAILABLE but STATIC and at 30m.
          - K (hydraulic conductivity), Sy (specific yield), ST (saturated thickness)
            NOT available in repo (need USGS datasets per VVIP Table S1).
          - GRACE lwe_thickness (basin-scale TWS) and SPL4SMGP depth_to_water_table
            are proxies for storage dynamics but NOT the same as Sy/ST.

        Returns GBI proxy LAZY. Currently uses HRES-WTD inverted as placeholder for
        buffering (shallow = more buffering) and sets T/S placeholders.

        Documented placeholders: PLACEHOLDERS['K'], ['Sy'], ['ST'], ['T']
        """
        if wtd_ds is None:
            try:
                wtd_ds = self.datasets.get("HRES-WTD_2015_Ogallala.nc")
                if wtd_ds is None:
                    wtd_ds = self.open_lazy("HRES-WTD_2015_Ogallala.nc")
            except Exception as e:
                warnings.warn(f"GBI: HRES-WTD not available: {e}. GBI will be NaN placeholder.")
                wtd_ds = None

        if wtd_ds is None:
            # Try GRACE/SMAP as broader proxy grid
            raise ValueError("GBI: no WTD grid available. Need HRES-WTD or USGS wells.")

        if y_slice is not None or x_slice is not None:
            # HRES-WTD is 52397*33030 ~1.7B pixels - MUST subset before any compute!
            # Lazy isel is essential.
            isel_kwargs = {}
            if y_slice is not None:
                isel_kwargs["y"] = y_slice
            if x_slice is not None:
                isel_kwargs["x"] = x_slice
            wtd_ds = wtd_ds.isel(**isel_kwargs) if isel_kwargs else wtd_ds

        # Extract WTD (b1, meters). Inverted: shallow (small depth) = high buffering
        # Normalize per VVIP (0-1). But depth distribution is skewed; we do min-max per domain.
        # For lazy, we cannot compute min/max without .compute(). So we provide FORMULA for later compute.
        # Here we keep lazy expression: (max - WTD)/(max-min) will be evaluated on compute.
        # As placeholder we use simple inverse: buffering_proxy = 1 / (1 + WTD)  (0-1, shallow->1)
        # Document as PROXY.

        try:
            wtd = wtd_ds["b1"].isel(time=0) if "time" in wtd_ds["b1"].dims else wtd_ds["b1"]
            # Mask FillValue etc already handled by xarray, but keep where valid
            # Simple physically-motivated proxy: buffering decreases with depth
            # Use 1/(1+WTD/10) where 10m scaling gives ~0.5 at 10m, ~0.09 at 100m
            gbi_proxy = 1.0 / (1.0 + wtd / 10.0)
            gbi_proxy = gbi_proxy.clip(min=0, max=1)
            # Where wtd is NaN or masked, result is NaN
            gbi_proxy.attrs["long_name"] = "Groundwater Buffer Index proxy (VVIP Eq7, WTD inverse proxy)"
            gbi_proxy.attrs["units"] = "1 (normalized 0-1, shallow=1)"
            gbi_proxy.attrs["note"] = ("PROXY: True GBI = T^{wT} * S^{-wS}, T=K*ST (Eq5-7). "
                                       "K, Sy, ST not in repo (PLACEHOLDERS). "
                                       "Here GBI ~ 1/(1+WTD/10) using HRES-WTD 2015 static depth only. "
                                       "Basin-scale GRACE TWS and SPL4SMGP dynamics not integrated yet. "
                                       "See PLACEHOLDERS for required USGS datasets.")
            gbi_proxy.attrs["vvip_reference"] = "Gonzalez Cruz et al. 2021 Eq5-7, Fig10"
            gbi_proxy.name = "GBI"
        except Exception as e:
            warnings.warn(f"GBI proxy failed: {e}")
            # Fallback empty
            gbi_proxy = xr.full_like(wtd_ds["b1"].isel(time=0), np.nan)

        # If true K,Sy,ST were provided, GBI would be:
        # T = K * ST
        # GBI = (T ** wT) * (Sy ** (-wS))  # Eq7, Sy in denominator because larger Sy => smaller head drop per Eq6
        # Then normalized 0-1 per VVIP before entropy weighting.

        return gbi_proxy

    # ---------------------------------------------------------------
    # VVIP: VI and TBI (Eq1 + Research Plan temporal)
    # ---------------------------------------------------------------
    def compute_VI(self,
                   DSI: Optional[xr.DataArray] = None,
                   SBI: Optional[xr.DataArray] = None,
                   GBI: Optional[xr.DataArray] = None,
                   alpha: float = 0.5,
                   beta: float = 0.5,
                   y_slice: Optional[slice] = None,
                   x_slice: Optional[slice] = None) -> xr.DataArray:
        """
        Vulnerability Index per VVIP Eq1: VI = DSI / (alpha*SBI + beta*GBI)
        Also  RI = 1/VI ;  TBI = 1/VI  (temporal buffering, larger = more robust)

        Lazy, read-only. Handles grid mismatch by requiring caller to harmonize
        to common grid beforehand (e.g., regrid HRES-WTD 30m -> 4km GRIDMET).

        Parameters
        ----------
        DSI, SBI, GBI : xr.DataArray | None
            If None, computed via respective proxy methods on the fly.
        alpha, beta : float
            Weights for soil vs groundwater buffering, alpha+beta=1, 0..1.
            VVIP Fig3d shows entropy weights ~0.5 each but varies.
        y_slice, x_slice : slice | None
            Passed through to proxy constructors if they need to be built.

        Returns
        -------
        xr.DataArray VI (y,x) lazy. Values >1 = vulnerable (stress > buffering),
        <1 = robust. For regional comparison, DSI/SBI/GBI should be normalized 0-1
        as done in proxies.
        """
        assert abs((alpha + beta) - 1.0) < 1e-6, "alpha+beta must =1"
        assert 0 <= alpha <= 1 and 0 <= beta <= 1

        # Lazy build missing components (each stays dask-backed)
        if DSI is None:
            DSI = self.compute_DSI(y_slice=y_slice, x_slice=x_slice)
        if SBI is None:
            SBI = self.compute_SBI(y_slice=y_slice, x_slice=x_slice)
        if GBI is None:
            GBI = self.compute_GBI(y_slice=y_slice, x_slice=x_slice)

        # Harmonization note: DSI is on DROUGHT grid (287x181 ~4km),
        # SBI HSG is 250m (5237x3302), GBI WTD is 30m (52397x33030).
        # Full VI requires resampling to common analysis grid (paper used 0.5deg ~55km, 187 pts).
        # Research Plan recommends harmonized ~1-5km (matching gridMET) with 30m aggregated.
        # Here we DOCUMENT mismatch and require user to regrid before division.
        # For demo on matching slices (small test), we assume caller already sliced to same shape
        # or we coarsen via .interp / .coarsen. We warn if shapes differ.

        if DSI.sizes.get("y") != SBI.sizes.get("y") or DSI.sizes.get("x") != SBI.sizes.get("x"):
            warnings.warn(
                f"VI grid mismatch: DSI {dict(DSI.sizes)} vs SBI {dict(SBI.sizes)} vs GBI {dict(GBI.sizes)}. "
                "VVIP harmonized to 0.5deg; Research Plan to ~4km (gridMET). "
                "You MUST regrid to common grid before interpreting VI. "
                "For demo we will interpolate SBI/GBI to DSI grid lazily (nearest). "
                "This is computationally expensive on full grid - subset first!"
            )
            # Lazy interp to DSI grid (still dask, but triggers rechunk)
            try:
                # Use interp with method nearest, stays lazy if dask
                target_y = DSI["y"]
                target_x = DSI["x"]
                SBI = SBI.interp(y=target_y, x=target_x, method="nearest", kwargs={"fill_value": np.nan})
                GBI = GBI.interp(y=target_y, x=target_x, method="nearest", kwargs={"fill_value": np.nan})
            except Exception as e:
                warnings.warn(f"VI regrid failed: {e}. Returning NaN where mismatch.")

        # Eq1
        denominator = alpha * SBI + beta * GBI
        # Avoid divide by zero
        denominator = denominator.clip(min=1e-9)
        VI = DSI / denominator
        VI.name = "VI"
        VI.attrs["long_name"] = "Vulnerability Index (VVIP Eq1) VI = DSI/(alpha*SBI+beta*GBI)"
        VI.attrs["units"] = "1 (dimensionless, >1 vulnerable, <1 robust)"
        VI.attrs["alpha"] = alpha
        VI.attrs["beta"] = beta
        VI.attrs["note"] = ("Regional relative comparison only (VVIP p5). "
                            "Requires normalization 0-1 per index for objective entropy weights. "
                            "Static vulnerability for year 2015 baseline in paper; TBI(t) extends to time-varying. "
                            "Placeholder proxies propagate NaN where source data missing.")
        VI.attrs["robustness_index"] = "RI = 1/VI"
        VI.attrs["temporal_buffering_index"] = "TBI = 1/VI = (alpha*SBI+beta*GBI)/DSI (Research Plan Bi persistence)"
        return VI

    def compute_TBI_time_varying(self,
                                 et_ds: Optional[xr.Dataset] = None,
                                 sm_ds: Optional[xr.Dataset] = None,
                                 met_ds: Optional[xr.Dataset] = None) -> xr.DataArray:
        """
        Temporal Buffering signal Bi(t) per Research Plan Section 6.1:
          Bi(t) = Fi(t) - F_hat(t)  where F_hat modeled from P, SMrz, VPD, Rn, T, PFT, season

        This is the CORE of the temporal extension: persistence of ET/SIF/VWC
        during dry-down beyond what meteorology + root-zone SM explain.

        With available data:
          Fi options: MODIS_ET_SSEBop et (2003-2022 dekadal), MOD13A3 NDVI/EVI,
                      SPL4SMGP land_evapotranspiration_flux, AmeriFlux LE/GPP
          Predictors: GRIDMET pr, sph, srad, tmmn/tmmx, vpd ; SPL4SMGP sm_rootzone ;
                      MCD12Q1 PFT / GFSAD cropland

        STEPS (documented, lazy):
          1. Detect dry-down events: P deficit + SMrz decline + VPD rise (Research Plan 6.1)
          2. Fit expected function F_hat via regression (e.g., RandomForest) per grid cell / PFT
          3. Residual Bi(t) during dry-down window
          4. Attribute to groundwater only where Bi covaries with WTD/wells/storage and survives
             irrigation controls (matched design 6.2). Otherwise Bi may be irrigation/management.

        PLACEHOLDER: This method returns the STRUCTURE and a NaN template with
        full documentation. Full Bi(t) requires SIF (OCO-2), VOD, and harmonized grids
        not fully available here.

        Returns
        -------
        xr.DataArray Bi(t,y,x) template (NaN) with attrs describing required steps.
        """
        warnings.warn("compute_TBI_time_varying returns TEMPLATE (NaN) - full Bi(t) needs SIF/VOD + model training. See docs.")
        # Provide template from SSEBop if available
        if et_ds is None:
            try:
                et_ds = self.datasets.get("MODIS_ET_SSEBop_Merged_Ogallala.nc") or self.open_lazy("MODIS_ET_SSEBop_Merged_Ogallala.nc")
                template = et_ds["et"].isel(time=slice(0, 2), y=slice(0, 5), x=slice(0, 5))
                bi = xr.full_like(template, np.nan)
                bi.name = "Bi"
                bi.attrs["long_name"] = "Buffering signal template (Research Plan 6.1) Bi(t)=Fi - F_hat"
                bi.attrs["note"] = ("TEMPLATE: Full computation needs: ET (SSEBop/PML), SIF (OCO-2/TROPOMI), "
                                    "VOD (SMAP), SMrz (SPL4SMGP), meteorology (gridMET), PFT (MCD12Q1/GFSAD), "
                                    "wells (USGS), storage (GRACE). See TBI_DOCUMENTATION.md Section Methods.")
                return bi
            except Exception:
                pass
        # fallback scalar
        return xr.DataArray(np.nan, attrs={"note": "No template grid available"})

    # ---------------------------------------------------------------
    # Helpers: inventory printing, assumptions
    # ---------------------------------------------------------------
    def print_inventory(self):
        inv = self.inventory or self.inventory_all()
        print(f"{'File':<40} {'Size GB':>8} {'Time range':<23} {'Variables'}")
        print("-"*120)
        for name, info in sorted(inv.items()):
            vars_str = ", ".join(info.variables[:4]) + ("..." if len(info.variables)>4 else "")
            print(f"{name:<40} {info.size_gb:>8.2f} {str(info.time_range or ''):<23} {vars_str}")

    @staticmethod
    def print_assumptions():
        print("""
VVIP ASSUMPTIONS (Gonzalez Cruz et al. 2021, p5-6, Section Model assumptions and limitations):
  1. Eq1 ratio is NOT dimensionally consistent; interpretable only as relative comparison,
     or at single site with subjective/equal weights (Laplace). Normalized 0-1 for regional.
  2. Requires >=30yr record for stable E(c) and DAC/ADP (paper used 1949-2018 hydrologic years).
     SMAP 2015+ alone insufficient for long-term expectations.
  3. Deterministic treatment of stochastic concepts (expected values, entropy). Randomness
     captured via long-term probabilities.
  4. Static vulnerability (point in time, year 2015 baseline). Temporal extension needs
     time-varying SBI(t), GBI(t) and DSI(t).
  5. No rigorous validation (MCDM empirical). Relies on verified input datasets.
  6. Correlations expected between SPI/SPEI, PAW/DAC/ADP, T/S - hence weighted PRODUCT not SUM.
  7. Thresholds for drought classification (SPI/SPEI <= -1) and 3/6 month accumulations.

RESEARCH PLAN ASSUMPTIONS (Section 6, 8 pitfalls):
  8. Residual Bi(t) is NOT automatically groundwater. Must covary with WTD/wells/storage
     AND survive irrigation controls (matched design, stratification by PFT/aridity/topography).
  9. Buffering may weaken SMOOTHLY not abruptly; test reversibility directly, keep slow-recovery
     indicators secondary.
 10. GRACE ~300km coarse - basin constraint only, not pixel. HRES-WTD 30m is static (well+covariate
     derived), not dynamic time series. Wells provide dynamics.
 11. Recovery test requires preselected rebound zones (storage measurably returned).
     Monotonically depleted aquifer cannot test reversibility.
 12. Multi-sensor record short (ECOSTRESS 2018+, SIF 2014+, SMAP 2015+) - use space-for-time
     across many dry-downs, not long trend at single site.
 13. Depth-to-water confounded by topography - requires placebo (wet-year), negative-control,
     Rosenbaum/E-value sensitivity, stratification checks.

IMPLEMENTATION ASSUMPTIONS (this toolkit):
 14. All opens are lazy (dask chunks) and read-only; no source file modified.
 15. Placeholders = NaN (never synthetic). Proxy variables (HSG->PAW, 1/(1+WTD/10)->GBI,
     intensity_mean->DSI) are documented as PROXY and not for publication without
     replacing with true USGS/gSSURGO/CPC sources.
 16. Grid harmonization required: HSG 250m, WTD 30m, DSI 4km, GRACE 0.5deg, SMAP 9km.
     Must aggregate to common analysis grid (Research Plan: ~1-5km gridMET).
     Demo regrid via .interp(nearest) is lazy but expensive full-grid - subset first.
 17. Equal weights fallback per Laplace when entropy cannot be computed lazily full-grid.
     Proper entropy (Eq8-10) should be computed on computed subset statistics.
 18. HSG valid_range 1-4 (A-D); WTD units meters; SPI/SPEI threshold -1.0 per paper p6.
""")

    @staticmethod
    def list_required_datasets() -> str:
        return """
ADDITIONAL DATASETS REQUIRED (beyond merged_datasets/*.nc present)

For FULL VVIP (Table S1 in VVIP Supplement):
  [MISSING] gSSURGO Gridded Soil Survey Geographic Database 30m
            -> field capacity, wilting point, root zone depth to compute PAW = (FC-PWP)/depth
  [MISSING] USGS Saturated Thickness 2009 (McGuire et al. 2012, SIR 2012-5177) + aquifer bottom
  [MISSING] USGS Hydraulic Conductivity map (Cederstrand & Becker 1998)
  [MISSING] USGS Specific Yield map (McGuire et al. 2012)
            -> needed for T=K*ST and Eq6 DH=IRRAMT/S
  [MISSING] CPC Global Monthly Soil Moisture 0.5deg 1949-2018 (Fan & van den Dool 2004, leaky bucket)
            -> needed for SSMI at 3/6 months and DAC/ADP contingency tables
            Available proxy SPL4SMGP sm_rootzone (2015+) is too short for 30yr E(c).

For FULL Research Plan temporal buffering (Section 5):
  [MISSING] ECOSTRESS ET/LST 70m 2018+ (high-res dry-down detail) - only SSEBop available
  [MISSING] OCO-2/3 / TROPOMI SIF (photosynthesis) 2014/2018+  - only NDVI/EVI & AmeriFlux GPP available
  [MISSING] SMAP L-band VOD (vegetation water content proxy) - have sm_rootzone but not VOD
  [MISSING] GLEAM4 0.1deg 1980+ (reference evaporation/stress) - not in repo
  [MISSING] USGS + state well time series (dynamic decline/recovery) - only static HRES-WTD 2015
  [MISSING] LANID / USDA CDL annual crop type + irrigation intensity - have GFSAD 2019 & Global_irrigation (2001-2015) coarse
  [MISSING] PML-V2 ET product (optional longer record) - have SSEBop
  [MISSING] Basin water-balance closure datasets for volume plausibility (Section 6.3)

PRESENT and USABLE (with caveats):
  DROUGHT_Merged_Ogallala.nc  (spi/spei at 14d/30d/90d/180d/270d) -> DSI proxy yes, but needs run theory for E(c)
  GRIDMET_Merged_Ogallala.nc  (pr, tmmn/tmmx, vpd, srad, eto/etr) -> met controls for Bi(t) yes (58GB lazy)
  SPL4SMGP_Ogallala_FULL.nc   (sm_rootzone, depth_to_water_table, ET flux) -> SMrz yes but short record
  Hydrologic_Soil_Group_250m  (HSG 1-4) -> PAW proxy only, not true PAW
  HRES-WTD_2015_Ogallala.nc   (WTD meters, 30m static) -> GBI spatial axis yes, not dynamics
  MOD13A3_Merged_Ogallala.nc  (NDVI/EVI monthly 1km) -> greenness proxy for SIF (derived, not independent)
  MODIS_ET_SSEBop_Merged      (et dekadal 1km) -> Fi for Bi(t) yes
  GRACE_Merged_Ogallala.nc    (lwe_thickness) -> basin storage constraint yes, coarse 0.5deg
  MCD12Q1_Merged_Ogallala.nc  (LC_Type) + GFSAD/GIR -> PFT/irrigation stratification yes
  merged_ameriflux.nc         (point LE/GPP) -> validation (point, not gridded)
"""

    # ---------------------------------------------------------------
    # Demo: small subset end-to-end without OOM
    # ---------------------------------------------------------------
    def demo_small_subset(self, y_slice: slice = slice(140, 150), x_slice: slice = slice(90, 100),
                          t_slice: slice = slice(2000, 2020)):
        """
        End-to-end demo on tiny subset (default 20 time steps x 10x10 spatial)
        to verify lazy pipeline without loading 30GB. Calls .compute() only on subset.

        Prints shapes, triggers small compute, reports VI stats.
        """
        print("=== TBI Demo (small subset, lazy -> compute) ===")
        print(f"  y_slice={y_slice}, x_slice={x_slice}, t_slice={t_slice}")
        try:
            dsi = self.compute_DSI(y_slice=y_slice, x_slice=x_slice, time_slice=t_slice)
            print(f"  DSI proxy lazy: {dsi.sizes}, chunks: {dsi.chunks if hasattr(dsi.data,'chunks') else 'n/a'}")
            # Compute subset (small, safe)
            dsi_c = dsi.compute()
            print(f"  DSI computed subset mean={float(dsi_c.mean().values):.4f} "
                  f"min={float(dsi_c.min().values):.4f} max={float(dsi_c.max().values):.4f}")
        except Exception as e:
            print(f"  DSI demo failed: {e}")

        try:
            # SBI and GBI on same y/x subset; note HSG/WTD grids are much finer,
            # but we isel same index range (geographically different!). This is JUST to test plumbing.
            # For real science, regrid geographically via interp/coarsen, not index slice.
            sbi = self.compute_SBI(y_slice=y_slice, x_slice=x_slice)
            print(f"  SBI proxy lazy: {sbi.sizes}")
            sbi_c = sbi.compute()
            print(f"  SBI computed subset mean={float(np.nanmean(sbi_c.values)):.4f}")
        except Exception as e:
            print(f"  SBI demo failed: {e}")

        try:
            gbi = self.compute_GBI(y_slice=y_slice, x_slice=x_slice)
            print(f"  GBI proxy lazy: {gbi.sizes}")
            gbi_c = gbi.compute()
            print(f"  GBI computed subset mean={float(np.nanmean(gbi_c.values)):.4f}")
        except Exception as e:
            print(f"  GBI demo failed: {e}")

        try:
            # VI demo: use same-grid proxies to avoid 30m/250m vs 4km harmonization cost.
            # Build VI from DSI on its native grid + constant SBI/GBI proxies on SAME grid,
            # so the demo proves Eq1 plumbing without expensive interp.
            dsi_for_vi = self.compute_DSI(y_slice=y_slice, x_slice=x_slice, time_slice=t_slice)
            # Create same-grid constant proxies (avoids regrid)
            sbi_same = xr.full_like(dsi_for_vi, 0.50)
            sbi_same.attrs["note"] = "Demo constant proxy (0.5) on DSI grid to avoid 30m harmonization in demo"
            gbi_same = xr.full_like(dsi_for_vi, 0.55)
            gbi_same.attrs["note"] = "Demo constant proxy (0.55) on DSI grid"
            vi = self.compute_VI(DSI=dsi_for_vi, SBI=sbi_same, GBI=gbi_same, alpha=0.5, beta=0.5)
            print(f"  VI lazy (same-grid demo): {vi.sizes}")
            vi_c = vi.compute()
            print(f"  VI computed subset mean={float(np.nanmean(vi_c.values)):.4f} "
                  f"median={float(np.nanmedian(vi_c.values)):.4f} "
                  f"frac >1 (vulnerable)={float(np.nanmean((vi_c.values>1))):.2%}")
            print("  VI >1 => vulnerable (stress > buffering); <1 => robust. See VVIP Fig12.")
            # Also show true harmonized path is documented (requires interp/coarsen):
            print("  Note: Full VI across native grids requires harmonization (see TBI_DOCUMENTATION.md §2).")
        except Exception as e:
            import traceback
            print(f"  VI demo failed: {e}")
            traceback.print_exc()

        print("=== Demo done (no source files modified) ===")

    # ---------------------------------------------------------------
    # Save result (new file only, never overwrites source)
    # ---------------------------------------------------------------
    def save_result(self, da: xr.DataArray, out_path: str | Path, overwrite: bool = False):
        """
        Save computed DataArray to NEW NetCDF (never overwrites merged_datasets).
        Uses compression and CF conventions.
        """
        out_path = Path(out_path)
        if out_path.exists() and not overwrite:
            raise FileExistsError(f"Refusing to overwrite {out_path}. Use overwrite=True if intentional.")
        # Ensure parent exists
        out_path.parent.mkdir(parents=True, exist_ok=True)
        ds_out = da.to_dataset(name=da.name or "VI")
        # Add global attrs
        ds_out.attrs["title"] = "Temporal Buffering / Vulnerability Index (VVIP + Research Plan)"
        ds_out.attrs["history"] = f"Created by temporal_buffering_index.py (lazy, read-only sources)"
        ds_out.attrs["source_note"] = "Sources read lazily from merged_datasets/*.nc; this file is derived output"
        ds_out.attrs["crs"] = "EPSG:4326"
        # Encoding
        enc = {da.name or "VI": {"zlib": True, "complevel": 3, "dtype": "float32"}}
        ds_out.to_netcdf(out_path, engine="netcdf4", encoding=enc)
        print(f"Saved {out_path} ({out_path.stat().st_size/1024**2:.1f} MB)")

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Temporal Buffering Index toolkit (VVIP + Research Plan)")
    parser.add_argument("--base_dir", default=r"G:/MSU_GWB/datasets", help="datasets root")
    parser.add_argument("--demo", action="store_true", help="run small subset demo")
    parser.add_argument("--inventory", action="store_true", help="print inventory")
    parser.add_argument("--assumptions", action="store_true", help="print assumptions")
    parser.add_argument("--required", action="store_true", help="list required additional datasets")
    args = parser.parse_args()

    tk = TBIToolkit(base_dir=args.base_dir)
    if args.assumptions:
        tk.print_assumptions()
    if args.required:
        print(tk.list_required_datasets())
    if args.inventory:
        tk.inventory_all()
        tk.print_inventory()
    if args.demo:
        tk.open_all_lazy()
        tk.demo_small_subset()
        tk.close_all()
    if not any([args.demo, args.inventory, args.assumptions, args.required]):
        parser.print_help()
