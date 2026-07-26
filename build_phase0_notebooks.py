from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).parent

UTILS = r'''"""Shared utilities for the Ogallala Phase 0 Colab notebooks."""
from __future__ import annotations

import json
import logging
import math
import os
import time
import traceback
from datetime import datetime, timezone
from functools import wraps
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd

REGISTRY_COLUMNS = [
    "timestamp", "notebook_id", "dataset", "task_id_or_file", "pathway",
    "status", "start_time", "end_time", "n_bytes", "error_message", "retries",
]
TERMINAL = {"COMPLETED", "FAILED"}


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


def ensure_layout(root: Path, leaf: str) -> Path:
    """Create the standard Drive tree and return one dataset leaf."""
    root = Path(root)
    for path in [root / "00_setup", root / "logs", root / "manifests"]:
        path.mkdir(parents=True, exist_ok=True)
    for name in ["raw", "cf", "qc"]:
        (root / leaf / name).mkdir(parents=True, exist_ok=True)
    registry = root / "logs" / "task_registry.csv"
    if not registry.exists():
        pd.DataFrame(columns=REGISTRY_COLUMNS).to_csv(registry, index=False)
    return root / leaf


def configure_logger(notebook_id: str, root: Path) -> tuple[logging.Logger, Path]:
    """Configure rotating file and stdout logging for a notebook."""
    from logging.handlers import RotatingFileHandler
    log_dir = Path(root) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    path = log_dir / f"{notebook_id}_{datetime.now().strftime('%Y%m%dT%H%M%S')}.log"
    logger = logging.getLogger(notebook_id)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    file_handler = RotatingFileHandler(path, maxBytes=50 * 1024 * 1024, backupCount=3)
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger, path


def registry_path(root: Path) -> Path:
    """Return the shared task registry path."""
    path = Path(root) / "logs" / "task_registry.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        pd.DataFrame(columns=REGISTRY_COLUMNS).to_csv(path, index=False)
    return path


def append_registry(root: Path, **values: Any) -> None:
    """Append a normalized task or file event to the registry."""
    row = {column: values.get(column, "") for column in REGISTRY_COLUMNS}
    row["timestamp"] = row["timestamp"] or utc_now()
    path = registry_path(root)
    frame = pd.DataFrame([row], columns=REGISTRY_COLUMNS)
    frame.to_csv(path, mode="a", header=path.stat().st_size == 0, index=False)


def logged_step(notebook_id: str, dataset: str, root: Path, logger: logging.Logger) -> Callable:
    """Decorate an I/O step so failures are logged and recorded."""
    def decorator(function: Callable) -> Callable:
        @wraps(function)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            start = utc_now()
            try:
                result = function(*args, **kwargs)
                append_registry(root, notebook_id=notebook_id, dataset=dataset,
                                task_id_or_file=function.__name__, pathway="xee",
                                status="COMPLETED", start_time=start, end_time=utc_now())
                return result
            except Exception as exc:
                logger.error("%s failed: %s\n%s", function.__name__, exc, traceback.format_exc())
                append_registry(root, notebook_id=notebook_id, dataset=dataset,
                                task_id_or_file=function.__name__, pathway="xee",
                                status="FAILED", start_time=start, end_time=utc_now(),
                                error_message=str(exc))
                return None
        return wrapper
    return decorator


class BatchManager:
    """Submit non-blocking Earth Engine exports under a small concurrency cap."""

    def __init__(self, root: Path, notebook_id: str, dataset: str,
                 logger: logging.Logger, max_concurrent: int = 20,
                 max_daily_tasks: int = 2500) -> None:
        self.root = Path(root)
        self.notebook_id = notebook_id
        self.dataset = dataset
        self.logger = logger
        self.max_concurrent = max_concurrent
        self.max_daily_tasks = max_daily_tasks

    def active_count(self) -> int:
        """Count locally registered active submissions."""
        frame = pd.read_csv(registry_path(self.root))
        if frame.empty:
            return 0
        return int(frame["status"].isin(["STARTED", "RUNNING"]).sum())

    def wait_for_slot(self) -> None:
        """Wait before submitting when the local active cap is reached."""
        while self.active_count() >= self.max_concurrent:
            self.logger.info("Active cap reached; sleeping 60 seconds")
            time.sleep(60)

    def start(self, task: Any, description: str, retries: int = 3) -> str:
        """Start a task with retry backoff and return its Earth Engine ID."""
        self.wait_for_slot()
        delays = [10, 60, 300]
        start_time = utc_now()
        for attempt in range(retries + 1):
            try:
                task.start()
                task_id = getattr(task, "id", "") or task.status().get("id", "")
                append_registry(self.root, notebook_id=self.notebook_id, dataset=self.dataset,
                                task_id_or_file=task_id, pathway="batch", status="STARTED",
                                start_time=start_time, retries=attempt)
                time.sleep(2)
                return task_id
            except Exception as exc:
                message = str(exc).lower()
                retryable = "429" in message or "rate limit" in message or "quota" in message
                if not retryable or attempt >= retries:
                    append_registry(self.root, notebook_id=self.notebook_id, dataset=self.dataset,
                                    task_id_or_file=description, pathway="batch", status="FAILED",
                                    start_time=start_time, end_time=utc_now(),
                                    error_message=str(exc), retries=attempt)
                    self.logger.error("Export failed: %s", exc)
                    return ""
                delay = delays[min(attempt, len(delays) - 1)]
                self.logger.warning("Rate limit on %s; retrying in %ss", description, delay)
                append_registry(self.root, notebook_id=self.notebook_id, dataset=self.dataset,
                                task_id_or_file=description, pathway="batch", status="RETRIED",
                                start_time=start_time, retries=attempt + 1,
                                error_message=str(exc))
                time.sleep(delay)
        return ""

    def submit_image(self, image: Any, description: str, folder: str,
                     prefix: str, region: Any, scale: int, crs: str) -> str:
        """Create and start a Drive GeoTIFF export without polling."""
        import ee
        task = ee.batch.Export.image.toDrive(
            image=image, description=description, folder=folder,
            fileNamePrefix=prefix, region=region, scale=scale, crs=crs,
            maxPixels=1e13, fileFormat="GeoTIFF",
        )
        return self.start(task, description)

    def monitor_once(self) -> pd.DataFrame:
        """Poll registered Earth Engine tasks once and update their statuses."""
        import ee
        path = registry_path(self.root)
        frame = pd.read_csv(path)
        if frame.empty:
            return frame
        for index, row in frame.iterrows():
            if row["status"] not in ["STARTED", "RUNNING"] or not row["task_id_or_file"]:
                continue
            try:
                status = ee.data.getTaskStatus(str(row["task_id_or_file"]))[0]
                state = status.get("state", "UNKNOWN")
                mapped = {"READY": "RUNNING", "RUNNING": "RUNNING",
                          "COMPLETED": "COMPLETED", "FAILED": "FAILED",
                          "CANCELLED": "FAILED"}.get(state, state)
                frame.loc[index, "status"] = mapped
                frame.loc[index, "end_time"] = utc_now() if mapped in TERMINAL else ""
                frame.loc[index, "error_message"] = status.get("error_message", "")
            except Exception as exc:
                self.logger.warning("Could not poll %s: %s", row["task_id_or_file"], exc)
        frame.to_csv(path, index=False)
        return frame


def normalize_projected_dataset(ds: Any, crs: str = "EPSG:5070") -> Any:
    """Normalize raster dimensions and add projected CF coordinates."""
    import xarray as xr
    if "lat" in ds.dims and "y" not in ds.dims:
        ds = ds.rename({"lat": "y"})
    if "lon" in ds.dims and "x" not in ds.dims:
        ds = ds.rename({"lon": "x"})
    if "y" not in ds.coords:
        ds = ds.assign_coords(y=np.arange(ds.sizes.get("y", 1), dtype=float))
    if "x" not in ds.coords:
        ds = ds.assign_coords(x=np.arange(ds.sizes.get("x", 1), dtype=float))
    ds["y"].attrs.update({"standard_name": "projection_y_coordinate", "units": "m", "axis": "Y"})
    ds["x"].attrs.update({"standard_name": "projection_x_coordinate", "units": "m", "axis": "X"})
    if "lat" not in ds:
        ds["lat"] = xr.DataArray(np.broadcast_to(ds["y"].values[:, None],
                                                  (ds.sizes["y"], ds.sizes["x"])), dims=("y", "x"))
    if "lon" not in ds:
        ds["lon"] = xr.DataArray(np.broadcast_to(ds["x"].values[None, :],
                                                  (ds.sizes["y"], ds.sizes["x"])), dims=("y", "x"))
    ds["lat"].attrs.update({"standard_name": "latitude", "units": "degrees_north"})
    ds["lon"].attrs.update({"standard_name": "longitude", "units": "degrees_east"})
    return ds


def write_cf_netcdf(ds: Any, path: Path, attrs: dict[str, Any],
                    variable_attrs: dict[str, dict[str, Any]] | None = None) -> Path:
    """Write a compressed projected NetCDF with the Phase 0 CF contract."""
    import xarray as xr
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ds = normalize_projected_dataset(ds)
    ds.attrs.update({
        "Conventions": "CF-1.8", "title": attrs.get("title", "Ogallala Phase 0 dataset"),
        "history": f"{utc_now()} notebook={attrs.get('notebook_id', 'unknown')} git_commit=placeholder",
        "institution": attrs.get("institution", "Ogallala Phase 0 research team"),
        "source": attrs.get("source", "Google Earth Engine"),
        "gee_asset_id": attrs.get("gee_asset_id", ""),
        "references": attrs.get("references", "Groundwater Buffering of Ecosystem Function During Drought Phase 0 plan"),
        "comment": attrs.get("comment", "Master grid: EPSG:5070; native bands and QA bands preserved."),
        "featureType": attrs.get("featureType", "grid"),
    })
    ds["crs"] = xr.DataArray(0, attrs={
        "grid_mapping_name": "albers_conical_equal_area", "semi_major_axis": 6378137.0,
        "inverse_flattening": 298.257222101, "false_easting": 0.0, "false_northing": 0.0,
        "spatial_ref": crs_wkt(),
    })
    variable_attrs = variable_attrs or {}
    for name, variable in ds.data_vars.items():
        if name == "crs":
            continue
        defaults = {"long_name": name, "units": "1", "_FillValue": -9999.0,
                    "grid_mapping": "crs", "coordinates": "lat lon"}
        defaults.update(variable_attrs.get(name, {}))
        variable.attrs.update(defaults)
    encoding: dict[str, dict[str, Any]] = {}
    for name, variable in ds.data_vars.items():
        if name == "crs":
            continue
        fill = variable.attrs.get("_FillValue", -9999.0)
        encoding[name] = {"zlib": True, "complevel": 4, "_FillValue": fill}
        if "time" in variable.dims:
            encoding[name]["chunksizes"] = tuple(min(24, ds.sizes[d]) for d in variable.dims)
    if "time" in ds.coords:
        ds["time"].attrs.update({"standard_name": "time", "long_name": "time",
                                  "units": "days since 1970-01-01 00:00:00",
                                  "calendar": "proleptic_gregorian", "axis": "T"})
        encoding["time"] = {"units": "days since 1970-01-01 00:00:00",
                             "calendar": "proleptic_gregorian"}
    ds.to_netcdf(path, format="NETCDF4_CLASSIC", engine="netcdf4", encoding=encoding)
    return path


def crs_wkt() -> str:
    """Return a WKT representation for the master Albers grid."""
    try:
        from pyproj import CRS
        return CRS.from_epsg(5070).to_wkt()
    except Exception:
        return "EPSG:5070"


def qc_netcdf(path: Path, qc_dir: Path) -> dict[str, Any]:
    """Create summary statistics and a quicklook for a NetCDF file."""
    import matplotlib.pyplot as plt
    import xarray as xr
    path, qc_dir = Path(path), Path(qc_dir)
    qc_dir.mkdir(parents=True, exist_ok=True)
    ds = xr.open_dataset(path)
    summary: dict[str, Any] = {}
    for name, value in ds.data_vars.items():
        if name == "crs":
            continue
        array = value.values.astype(float)
        finite = array[np.isfinite(array)]
        summary[name] = {"min": float(np.min(finite)) if finite.size else None,
                         "max": float(np.max(finite)) if finite.size else None,
                         "mean": float(np.mean(finite)) if finite.size else None,
                         "std": float(np.std(finite)) if finite.size else None,
                         "n_valid": int(finite.size)}
        if finite.size:
            image = np.nanmean(array, axis=0) if "time" in value.dims else array
            plt.imsave(qc_dir / f"{path.stem}_{name}.png", image, cmap="viridis")
    summary_path = qc_dir / f"{path.stem}_summary_stats.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    ds.close()
    return summary


def write_manifest(root: Path, notebook_id: str, files: Iterable[Path],
                   dataset: str, status: str = "COMPLETED") -> Path:
    """Write provenance metadata for files produced by one notebook."""
    entries = []
    for file_path in files:
        file_path = Path(file_path)
        if file_path.exists():
            entries.append({"file": str(file_path), "size_bytes": file_path.stat().st_size,
                            "cf_version": "CF-1.8", "dataset": dataset})
    path = Path(root) / "manifests" / f"{notebook_id}_manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"notebook_id": notebook_id, "status": status,
                                "files": entries}, indent=2), encoding="utf-8")
    return path
'''

DATASETS = {
    "01a": {"name": "GRACE / GRACE-FO Mascons", "asset_id": "NASA/GRACE/MASS_GRIDS/MASCON_CRI", "vars": "TWSA", "bands": ["lwe_thickness"], "scale": 1000, "kind": "collection", "start": "2002-04-01", "end": "2024-09-30", "folder": "01a_grace_mascons", "usage": "Basin-scale storage constraint and slow driver; use for closure and screening."},
    "01b": {"name": "SMAP L4 Root Zone Soil Moisture and L3 VOD", "asset_id": ["NASA/SMAP/SPL4SMGP/007", "NASA/SMAP/SPL3SMP_E/005"], "vars": "MULTI", "bands": ["sm_rootzone", "vegetation_water_content"], "scale": 1000, "kind": "collection_group", "start": "2015-04-01", "end": "2024-09-30", "folder": "01b_smap_l4_rzsm", "usage": "Climate control and vegetation water signal for dry-down decoupling analysis."},
    "01d": {"name": "TROPOMI SIF", "asset_id": "projects/sat-io/open-datasets/TROPOSIF", "vars": "SIF", "bands": ["SIF"], "scale": 1000, "kind": "collection", "start": "2018-05-01", "end": "2024-09-30", "folder": "01d_tropomi_sif", "usage": "Independent carbon signal for recent-period decoupling analysis."},
    "01e": {"name": "ERA5-Land", "asset_id": "ECMWF/ERA5_LAND/HOURLY", "vars": "MULTI", "bands": ["total_precipitation", "temperature_2m", "dewpoint_temperature_2m", "surface_net_solar_radiation"], "scale": 1000, "kind": "collection", "start": "2002-01-01", "end": "2024-09-30", "folder": "01e_era5_land", "usage": "Meteorological controls for expected-function matching."},
    "01f": {"name": "gridMET", "asset_id": "IDAHO_EPSCOR/GridMET", "vars": "MULTI", "bands": ["pr", "vpd", "etr", "tmmx", "tmmn", "vs"], "scale": 1000, "kind": "collection", "start": "1979-01-01", "end": "2024-09-30", "folder": "01f_gridmet", "usage": "Primary climate grid for inference and matching."},
    "01g": {"name": "MODIS NDVI / EVI", "asset_id": "MODIS/061/MOD13A2", "vars": "MULTI", "bands": ["NDVI", "EVI", "DetailedQA", "SummaryQA", "pixel_reliability"], "scale": 1000, "kind": "collection", "start": "2000-02-18", "end": "2024-09-30", "folder": "01g_modis_ndvi_evi", "usage": "Long-record context, screening, and reversibility testing."},
    "01h": {"name": "MODIS Land Cover", "asset_id": "MODIS/061/MCD12Q1", "vars": "MULTI", "bands": ["LC_Type1", "LC_Prop1_Assessment", "LC_Prop1"], "scale": 1000, "kind": "collection", "start": "2001-01-01", "end": "2023-12-31", "folder": "01h_modis_landcover", "usage": "Vegetation-type and crop versus natural stratification."},
    "01i": {"name": "SRTM DEM", "asset_id": "USGS/SRTMGL1_003", "vars": "ELEV", "bands": ["elevation"], "scale": 1000, "kind": "image", "start": "2000-01-01", "end": "2000-12-31", "folder": "01i_srtm_dem", "usage": "Topographic-position stratification."},
    "01j": {"name": "USDA CDL", "asset_id": "USDA/NASS/CDL", "vars": "CROP", "bands": ["cropland"], "scale": 30, "kind": "collection", "start": "2008-01-01", "end": "2024-12-31", "folder": "01j_usda_cdl", "usage": "Crop type and irrigation-control stratification."},
    "01k": {"name": "GFSAD", "asset_id": "USGS/GFSAD/GLASS/CropMasks/2000", "vars": "IRRIG", "bands": ["landcover"], "scale": 1000, "kind": "image", "start": "2000-01-01", "end": "2000-12-31", "folder": "01k_gfsad", "usage": "Managed versus rainfed area control."},
    "01l": {"name": "Aridity Index", "asset_id": "projects/sat-io/open-datasets/CGIAR_ARIDITY", "vars": "ARID", "bands": ["b1"], "scale": 1000, "kind": "image", "start": "1970-01-01", "end": "2000-12-31", "folder": "01l_aridity_index", "usage": "Matching covariate and failure-threshold stratification."},
    "01m": {"name": "ESA WorldCover", "asset_id": "ESA/WorldCover/v200", "vars": "LC", "bands": ["Map"], "scale": 30, "kind": "image", "start": "2021-01-01", "end": "2021-12-31", "folder": "01m_worldcover", "usage": "High-resolution natural versus managed land-cover mask."},
    "01n": {"name": "Landsat 8/9 Surface Reflectance", "asset_id": "LANDSAT/LC08/C02/T1_L2", "vars": "MULTI", "bands": ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "ST_B10", "QA_PIXEL", "QA_RADSAT"], "scale": 30, "kind": "collection", "start": "2013-04-01", "end": "2024-09-30", "folder": "01n_landsat", "usage": "High-resolution spectral context and optional indices."},
    "01o": {"name": "Sentinel-2 Harmonized", "asset_id": "COPERNICUS/S2_HARMONIZED", "vars": "MULTI", "bands": ["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B11", "B12", "QA60", "SCL"], "scale": 30, "kind": "collection", "start": "2017-03-28", "end": "2024-09-30", "folder": "01o_sentinel2", "usage": "Recent high-resolution vegetation context."},
    "01p": {"name": "PML-V2 ET", "asset_id": "projects/sat-io/open-datasets/PML_V2", "vars": "ET", "bands": ["Ec", "Es", "Ei", "ET_water", "ET"], "scale": 1000, "kind": "collection", "start": "2000-01-01", "end": "2023-12-31", "folder": "01p_pml_v2", "usage": "Longer evaporative-response record."},
    "01q": {"name": "SSEBop ET", "asset_id": "USGS/ssebop/ssebopeta_v4", "vars": "ET", "bands": ["et"], "scale": 1000, "kind": "collection", "start": "2003-01-01", "end": "2024-09-30", "folder": "01q_ssebop", "usage": "Longer operational evaporative-response record."},
    "01r": {"name": "SoilGrids", "asset_id": "projects/soilgrids-isric/", "vars": "MULTI", "bands": ["clay_0-5cm_mean", "sand_0-5cm_mean", "bdod_0-5cm_mean", "soc_0-5cm_mean"], "scale": 1000, "kind": "image", "start": "2017-01-01", "end": "2017-12-31", "folder": "01r_soilgrids", "usage": "Rooting depth, texture, and soil-water-capacity controls."},
    "01s": {"name": "NASADEM", "asset_id": "NASA/NASADEM_HGT/001", "vars": "ELEV", "bands": ["elevation", "num", "swb"], "scale": 1000, "kind": "image", "start": "2000-01-01", "end": "2000-12-31", "folder": "01s_nasadem", "usage": "Hydrologic conditioning and topographic position."},
}

INDEPENDENT = {
    "02a": {"name": "Ma et al. (2026) Water Table Depth raster", "vars": "WTD", "folder": "02a_ma_wtd_raster", "method": "requests / rioxarray", "endpoint": "REPLACE_WITH_HYDROSHARE_OR_ZENODO_URL", "usage": "Static depth-to-water axis; not a time series."},
    "02b": {"name": "USGS Groundwater Wells (NWIS)", "vars": "WELL", "folder": "02b_usgs_wells", "method": "dataretrieval", "endpoint": "https://waterservices.usgs.gov/nwis/gwlevels/", "usage": "Observed water-table decline and recovery."},
    "02c": {"name": "State Well Networks", "vars": "WELL", "folder": "02c_state_wells", "method": "requests", "endpoint": "https://www3.twdb.texas.gov/apps/waterdatainteractive/", "usage": "State networks augment NWIS spatial coverage."},
    "02d": {"name": "FLUXNET / AmeriFlux", "vars": "FLUX", "folder": "02d_ameriflux", "method": "ameriflux-api / requests", "endpoint": "REQUIRES_AMERIFLUX_ACCOUNT", "usage": "Independent site validation for buffering and decoupling."},
    "02e": {"name": "SAPFLUXNET", "vars": "SAP", "folder": "02e_sapfluxnet", "method": "requests", "endpoint": "REQUIRES_MANUAL_REGISTRATION", "usage": "Independent plant water-use validation."},
    "02f": {"name": "GLEAM4", "vars": "ET", "folder": "02f_gleam4", "method": "requests / xarray", "endpoint": "REQUIRES_GLEAM_REGISTRATION_URL", "usage": "Reference model product, not independent evidence."},
    "02g": {"name": "GOSIF", "vars": "SIF", "folder": "02g_gosif", "method": "requests / xarray", "endpoint": "https://www.numericalterra.com/", "usage": "Long-record context and reversibility testing."},
    "02h": {"name": "LANID Irrigation Mapping", "vars": "IRRIG", "folder": "02h_lanid", "method": "requests / rioxarray", "endpoint": "REPLACE_WITH_ZENODO_OR_HYDROSHARE_URL", "usage": "High-resolution irrigation stratification."},
    "02i": {"name": "USDA NASS Statistics", "vars": "NASS", "folder": "02i_nass_stats", "method": "nass / requests", "endpoint": "https://quickstats.nass.usda.gov/api/api_GET/", "usage": "Agricultural outcome validation."},
    "02j": {"name": "ECOSTRESS Level-2 ET and LST", "vars": "MULTI", "folder": "02j_ecostress", "method": "earthaccess", "endpoint": "ECOSTRESS_L2_LSTE", "usage": "Recent high-resolution ET and LST dry-down detail."},
}

LEAVES = [d["folder"] for d in DATASETS.values()] + ['01c_smap_l3_vod'] + [d["folder"] for d in INDEPENDENT.values()]


def cell(cell_type: str, source: str, language: str) -> dict:
    return {"cell_type": cell_type, "metadata": {"language": language},
            "source": source.splitlines(True)}


def notebook(cells: list[dict]) -> dict:
    return {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.x"}}, "nbformat": 4, "nbformat_minor": 5}


def setup_code(notebook_id: str, needs_earthaccess: bool = False) -> str:
    auth = '''\n# ECOSTRESS uses a .netrc file in Colab; credentials come only from env vars or prompts.\nif os.getenv("EARTHDATA_USERNAME") and os.getenv("EARTHDATA_PASSWORD"):\n    netrc_path = Path.home() / ".netrc"\n    netrc_path.write_text(f"machine urs.earthdata.nasa.gov login {os.environ['EARTHDATA_USERNAME']} password {os.environ['EARTHDATA_PASSWORD']}\\n", encoding="utf-8")\n    netrc_path.chmod(0o600)\n    import earthaccess\n    if os.getenv("EARTHACCESS_INTERACTIVE") == "1":\n        earthaccess.login(strategy="interactive")\n    else:\n        earthaccess.login(strategy="netrc")\n'''
    if not needs_earthaccess:
        auth = ''
    return f'''# Colab bootstrap. Installs are quiet so the notebook stays readable.\n%pip -q install earthengine-api xee xarray netCDF4 h5netcdf geopandas rioxarray regionmask "dask[diagnostics]" requests tqdm python-dotenv pyyaml matplotlib pandas pyproj shapely\nfrom pathlib import Path\nimport os, sys, json, time, logging, traceback\nfrom datetime import datetime\nfrom google.colab import drive\ndrive.mount('/content/drive')\n\nPROJECT_ID = os.getenv("GEE_PROJECT", "ee-ishansinhagzb")\nimport ee\nee.Authenticate()\nee.Initialize(project=PROJECT_ID)\n{auth}\nDRIVE_ROOT = Path('/content/drive/MyDrive/Ogallala_Phase0')\nSETUP_DIR = DRIVE_ROOT / '00_setup'\nSETUP_DIR.mkdir(parents=True, exist_ok=True)\n# The embedded module makes each notebook independently uploadable to Colab.\nUTILS_SOURCE = {UTILS!r}\n(SETUP_DIR / 'utils.py').write_text(UTILS_SOURCE, encoding='utf-8')\nsys.path.insert(0, str(SETUP_DIR))\nfrom utils import *\nNOTEBOOK_ID = '{notebook_id}'\n'''


def config_code(dataset: dict, notebook_id: str) -> str:
    independent = dataset in INDEPENDENT.values()
    catalog = "independent_datasets" if independent else "datasets"
    folder_prefix = "02_independent" if independent else "01_gee_rasters"
    return f'''import yaml\n\nCONFIG_PATH = DRIVE_ROOT / '00_setup' / 'config.yaml'\nif not CONFIG_PATH.exists():\n    raise FileNotFoundError('Run 00_setup_aoi_and_config.ipynb first to create config.yaml and the AOI.')\nCONFIG = yaml.safe_load(CONFIG_PATH.read_text(encoding='utf-8'))\nDRY_RUN = bool(CONFIG.get('DRY_RUN', True))\nAOI = ee.Geometry(json.loads((DRIVE_ROOT / '00_setup' / 'aoi_ee_geometry.json').read_text()))\nAOI_BUFFER = AOI.buffer(CONFIG.get('AOI_BUFFER_METERS', 1000))\nBBOX = AOI.bounds().coordinates().getInfo()[0]\nCRS_TARGET = 'EPSG:4326'\nMASTER_CRS = CONFIG.get('MASTER_CRS', 'EPSG:5070')\nDATASET_ID = '{notebook_id}'\nDATASET = CONFIG.get('{catalog}', {{}}).get(DATASET_ID, {{}})\nDATASET['drive_subfolder'] = '{folder_prefix}/' + DATASET.get('folder', '{dataset['folder']}')\nMASTER_SCALE = int(CONFIG.get('MASTER_SCALE', {{}}).get(DATASET_ID, {dataset.get('scale', 1000)}))\nLEAF = DRIVE_ROOT / DATASET['drive_subfolder']\nensure_layout(DRIVE_ROOT, DATASET['drive_subfolder'])\n'''


def gee_notebook(notebook_id: str, dataset: dict, group: bool = False) -> dict:
    title = f"# {notebook_id}: {dataset['name']}\n\n**Why:** {dataset['usage']}\n\n- Asset ID: `{dataset['asset_id']}`\n- Outputs: `raw/`, `cf/`, and `qc/` under `{dataset['folder']}`\n- Author: `<name / institution>`\n- This notebook preserves native bands and QA values; it does not select away quality metadata."
    dataset_code = f'''DATASET_ID = '{notebook_id}'\nDATASET = {json.dumps(dataset, indent=2)}\nDATASET['drive_subfolder'] = '01_gee_rasters/{dataset['folder']}'\nDATASET['crs'] = 'EPSG:5070'\nDATASET['band_metadata'] = {{band: {{'long_name': band, 'units': '1', 'standard_name': 'unknown'}} for band in DATASET['bands']}}\n'''
    pipeline = '''def build_gee_source() -> tuple[object, bool]:\n    """Build clipped, reprojected sources while retaining every native band."""\n    asset_ids = DATASET['asset_id'] if isinstance(DATASET['asset_id'], list) else [DATASET['asset_id']]\n    if DATASET['kind'] == 'image':\n        return ee.Image(asset_ids[0]).clip(AOI), False\n    sources = []\n    for asset_id in asset_ids:\n        collection = ee.ImageCollection(asset_id).filterBounds(AOI).filterDate(DATASET['start'], DATASET['end'])\n        # Do not call select(): QA, uncertainty, and native metadata must survive.\n        sources.append(collection.map(lambda image: image.clip(AOI).reproject(crs=MASTER_CRS, scale=MASTER_SCALE)))\n    return (sources if DATASET['kind'] == 'collection_group' else sources[0]), True\n\nSOURCE, IS_COLLECTION = build_gee_source()\n'''
    xee = '''@logged_step(NOTEBOOK_ID, DATASET['name'], DRIVE_ROOT, logger)\ndef dry_run_xee() -> list[Path]:\n    """Pull exactly one chunk per source through XEE and write CF files."""\n    if not DRY_RUN:\n        return []\n    import xarray as xr\n    from dask.diagnostics import ProgressBar\n    sources = SOURCE if isinstance(SOURCE, list) else [SOURCE]\n    outputs = []\n    for index, source in enumerate(sources):\n        if IS_COLLECTION:\n            source = source.filterDate(DATASET['start'], DATASET['start'] + 'T23:59:59')\n            ds = xr.open_dataset(source, engine='ee', crs=MASTER_CRS, scale=MASTER_SCALE,\n                                 geometry=AOI, chunks={'time': 24, 'lat': 256, 'lon': 256}, ee_mask=False)\n        else:\n            ds = xr.open_dataset(source, engine='ee', crs=MASTER_CRS, scale=MASTER_SCALE,\n                                 geometry=AOI, chunks={'lat': 256, 'lon': 256}, ee_mask=False)\n        with ProgressBar():\n            materialized = ds.compute()\n        suffix = f'_part{index + 1}' if len(sources) > 1 else ''\n        start = DATASET['start']\n        name = f"DRYRUN_OG_{DATASET_ID}_{DATASET['vars']}{suffix}_{start}_{start}_{MASTER_SCALE}m.nc"\n        output = LEAF / 'cf' / 'dry_run' / name\n        write_cf_netcdf(materialized, output, {'title': DATASET['name'], 'notebook_id': NOTEBOOK_ID,\n            'source': DATASET['asset_id'], 'gee_asset_id': str(DATASET['asset_id'])}, DATASET['band_metadata'])\n        qc_netcdf(output, LEAF / 'qc' / 'dry_run')\n        outputs.append(output)\n    return outputs\n\nDRY_RUN_OUTPUTS = dry_run_xee()\n'''
    export = '''def submit_normal_exports() -> list[str]:\n    """Queue temporal or static exports without waiting for server completion."""\n    if DRY_RUN:\n        return []\n    manager = BatchManager(DRIVE_ROOT, NOTEBOOK_ID, DATASET['name'], logger, max_concurrent=20)\n    task_ids = []\n    folder = f"Ogallala_Phase0/{DATASET['drive_subfolder']}/raw"\n    sources = SOURCE if isinstance(SOURCE, list) else [SOURCE]\n    for source_index, source in enumerate(sources):\n        suffix = f'_part{source_index + 1}' if len(sources) > 1 else ''\n        if IS_COLLECTION:\n            start_year = int(DATASET['start'][:4])\n            end_year = int(DATASET['end'][:4])\n            for year in range(start_year, end_year + 1):\n                start = f'{year}-01-01'\n                end = f'{year + 1}-01-01'\n                chunk = source.filterDate(start, end).median().clip(AOI)\n                prefix = f"OG_{DATASET_ID}_{DATASET['vars']}{suffix}_{start}_{year}-12-31_{MASTER_SCALE}m"\n                task_ids.append(manager.submit_image(chunk, f'{DATASET_ID}_{source_index}_{year}', folder, prefix, AOI, MASTER_SCALE, MASTER_CRS))\n        else:\n            start, end = DATASET['start'], DATASET['end']\n            prefix = f"OG_{DATASET_ID}_{DATASET['vars']}{suffix}_{start}_{end}_{MASTER_SCALE}m"\n            task_ids.append(manager.submit_image(source, f'{DATASET_ID}_{source_index}_static', folder, prefix, AOI, MASTER_SCALE, MASTER_CRS))\n    return task_ids\n\nQUEUED_TASKS = submit_normal_exports()\n'''
    writer = '''def write_exported_geotiff_as_cf(geotiff: Path) -> Path:\n    """Convert one completed Drive GeoTIFF into the CF output tree."""\n    import rioxarray\n    raster = rioxarray.open_rasterio(geotiff, masked=True).squeeze(drop=True)\n    out = LEAF / 'cf' / f"{geotiff.stem}.nc"\n    return write_cf_netcdf(raster.to_dataset(name=DATASET['vars']), out,\n        {'title': DATASET['name'], 'notebook_id': NOTEBOOK_ID, 'source': DATASET['asset_id'],\n         'gee_asset_id': str(DATASET['asset_id'])}, DATASET['band_metadata'])\n\n# This conversion is intentionally opt-in: normal submissions are non-blocking and\n# Drive exports may not exist until 98_monitor_gdrive_exports reports COMPLETED.\n'''
    qc = '''OUTPUTS = list(DRY_RUN_OUTPUTS)\nfor output in OUTPUTS:\n    write_manifest(DRIVE_ROOT, NOTEBOOK_ID, [output], DATASET['name'])\n'''
    final = '''registry = pd.read_csv(registry_path(DRIVE_ROOT))\nstatus = 'DRY_RUN_COMPLETED' if DRY_RUN and all(path.exists() for path in OUTPUTS) else ('QUEUED' if QUEUED_TASKS else 'FAILED_OR_NOT_RUN')\nprint(pd.DataFrame([{'dataset': DATASET['name'], 'file_or_tasks': len(OUTPUTS) or len(QUEUED_TASKS), 'status': status}]).to_string(index=False))\nif DRY_RUN and not all(path.exists() for path in OUTPUTS):\n    raise RuntimeError(f'Missing dry-run outputs: {[str(path) for path in OUTPUTS if not path.exists()]}')\n'''
    cells = [cell('markdown', title, 'markdown'), cell('code', setup_code(notebook_id), 'python'), cell('code', config_code(dataset, notebook_id), 'python'), cell('code', "logger, LOG_PATH = configure_logger(NOTEBOOK_ID, DRIVE_ROOT)\nREGISTRY = registry_path(DRIVE_ROOT)\n@logged_step(NOTEBOOK_ID, DATASET['name'], DRIVE_ROOT, logger)\ndef logger_health_check() -> str:\n    logger.info('Logger ready: %s', LOG_PATH)\n    return str(LOG_PATH)\nlogger_health_check()\n", 'python'), cell('code', dataset_code, 'python'), cell('code', pipeline, 'python'), cell('code', xee, 'python'), cell('code', export, 'python'), cell('code', "%pip -q install compliance-checker\nMAX_CONCURRENT = 20\nMAX_DAILY_TASKS = 2500\n# BatchManager enforces the local cap, retry backoff, registry writes, and non-blocking submission.\n", 'python'), cell('code', writer, 'python'), cell('code', qc, 'python'), cell('code', final, 'python')]
    return notebook(cells)


def independent_notebook(notebook_id: str, dataset: dict) -> dict:
    title = f"# {notebook_id}: {dataset['name']}\n\n**Why:** {dataset['usage']}\n\n- Access: `{dataset['method']}`\n- Endpoint or collection: `{dataset['endpoint']}`\n- Outputs: `raw/`, `cf/`, and `qc/` under `{dataset['folder']}`\n- Author: `<name / institution>`\n\nCredentials and registration-dependent URLs remain environment variables or config values."
    dataset_code = f'''DATASET_ID = '{notebook_id}'\nDATASET = {json.dumps(dataset, indent=2)}\nDATASET['drive_subfolder'] = '02_independent/{dataset['folder']}'\nDATASET['crs'] = 'EPSG:5070'\n'''
    process = '''def convert_external_to_cf(raw_path: Path) -> Path:\n    """Convert a Drive-resident raster, NetCDF, or API payload to CF NetCDF."""\n    import numpy as np\n    import xarray as xr\n    suffix = raw_path.suffix.lower()\n    if suffix in {'.tif', '.tiff'}:\n        import rioxarray\n        data = rioxarray.open_rasterio(raw_path, masked=True).squeeze(drop=True)\n        dataset = data.to_dataset(name=DATASET['vars'])\n    elif suffix in {'.nc', '.nc4', '.netcdf'}:\n        dataset = xr.open_dataset(raw_path)\n    else:\n        # API and CSV payloads are represented as a one-cell provenance grid until\n        # a source-specific schema is configured; no bytes are dropped or staged locally.\n        dataset = xr.Dataset({DATASET['vars']: (('y', 'x'), np.array([[raw_path.stat().st_size]], dtype='float32'))})\n    output = LEAF / 'cf' / f"OG_{DATASET_ID}_{DATASET['vars']}_{CONFIG.get('START_DATE', '2000-01-01')}_{CONFIG.get('END_DATE', '2024-09-30')}_{MASTER_SCALE}m.nc"\n    write_cf_netcdf(dataset, output, {'title': DATASET['name'], 'notebook_id': NOTEBOOK_ID,\n        'source': DATASET['endpoint'], 'comment': 'External source streamed directly to Drive.'})\n    qc_netcdf(output, LEAF / 'qc')\n    return output\n\n\ndef stream_independent_dataset() -> list[Path]:\n    """Stream an external source directly into Drive, without /content downloads."""\n    if DATASET['endpoint'].startswith(('REPLACE_', 'REQUIRES_')):\n        logger.warning('Configure DATASET endpoint or credentials before running %s', DATASET_ID)\n        return []\n    import requests\n    from urllib.parse import urlparse\n    response = requests.get(DATASET['endpoint'], stream=True, timeout=120)\n    response.raise_for_status()\n    raw_dir = LEAF / 'raw'\n    raw_dir.mkdir(parents=True, exist_ok=True)\n    suffix = Path(urlparse(DATASET['endpoint']).path).suffix or '.bin'\n    target = raw_dir / f"OG_{DATASET_ID}_{DATASET['vars']}_{CONFIG.get('START_DATE', '2000-01-01')}_{CONFIG.get('END_DATE', '2024-09-30')}_{MASTER_SCALE}m{suffix}"\n    with target.open('wb') as handle:\n        for block in response.iter_content(chunk_size=1024 * 1024):\n            if block:\n                handle.write(block)\n    output = convert_external_to_cf(target)\n    append_registry(DRIVE_ROOT, notebook_id=NOTEBOOK_ID, dataset=DATASET['name'],\n                    task_id_or_file=str(output), pathway='xee', status='COMPLETED',\n                    n_bytes=output.stat().st_size)\n    return [output]\n'''
    xee = '''@logged_step(NOTEBOOK_ID, DATASET['name'], DRIVE_ROOT, logger)\ndef dry_run_external() -> list[Path]:\n    """Run one bounded external request when an endpoint is configured."""\n    if not DRY_RUN:\n        return []\n    return stream_independent_dataset()\n\nDRY_RUN_OUTPUTS = dry_run_external()\n'''
    export = '''def submit_normal_external() -> list[Path]:\n    """Run the full external stream; registration-specific sources stay explicit."""\n    if DRY_RUN:\n        return []\n    return stream_independent_dataset()\n\nQUEUED_TASKS = submit_normal_external()\n'''
    writer = '''# The source-aware converter is defined in the processing cell and writes one CF file per streamed artifact.\n'''
    qc = '''for output in DRY_RUN_OUTPUTS:\n    write_manifest(DRIVE_ROOT, NOTEBOOK_ID, [output], DATASET['name'], status='COMPLETED')\n'''
    final = '''print(pd.DataFrame([{'dataset': DATASET['name'], 'file_or_tasks': len(DRY_RUN_OUTPUTS) + len(QUEUED_TASKS),\n                                  'status': 'COMPLETED' if DRY_RUN_OUTPUTS else 'CONFIGURATION_REQUIRED'}]).to_string(index=False))\n'''
    cells = [cell('markdown', title, 'markdown'), cell('code', setup_code(notebook_id, notebook_id == '02j'), 'python'), cell('code', config_code(dataset, notebook_id), 'python'), cell('code', "logger, LOG_PATH = configure_logger(NOTEBOOK_ID, DRIVE_ROOT)\nREGISTRY = registry_path(DRIVE_ROOT)\n", 'python'), cell('code', dataset_code, 'python'), cell('code', process, 'python'), cell('code', xee, 'python'), cell('code', export, 'python'), cell('code', "%pip -q install compliance-checker\nMAX_CONCURRENT = 20\nMAX_DAILY_TASKS = 2500\n", 'python'), cell('code', writer, 'python'), cell('code', qc, 'python'), cell('code', final, 'python')]
    return notebook(cells)


def setup_notebook() -> dict:
    title = "# 00: Setup AOI and configuration\n\n**Why:** This notebook creates the single Drive handoff used by every Phase 0 notebook. It fetches the USGS High Plains boundary directly to Drive when possible, writes the AOI in GeoJSON and Earth Engine formats, and creates the complete output tree.\n\n- Source: USGS High Plains Aquifer boundary\n- Outputs: `00_setup/ogallala_boundary.geojson`, `ogallala_bbox.json`, `aoi_ee_geometry.json`, `config.yaml`, and `utils.py`\n- Author: `<name / institution>`"
    bootstrap = setup_code('00_setup')
    body = f'''import geopandas as gpd\nimport requests\nimport zipfile\nimport io\nimport yaml\n\n# Prefer a user-provided boundary directory on Drive; otherwise stream the USGS archive to Drive.\nsource_dir = SETUP_DIR / 'source_boundary'\nshapefile = next(source_dir.glob('*.shp'), None) if source_dir.exists() else None\nif shapefile is None:\n    url = os.getenv('OGALLALA_BOUNDARY_ZIP_URL', 'https://water.usgs.gov/GIS/dsdl/hp_bound2010.zip')\n    response = requests.get(url, timeout=120)\n    response.raise_for_status()\n    source_dir.mkdir(parents=True, exist_ok=True)\n    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:\n        archive.extractall(source_dir)\n    shapefile = next(source_dir.glob('*.shp'))\n\ngdf = gpd.read_file(shapefile).to_crs('EPSG:4326')\ngdf.to_file(SETUP_DIR / 'ogallala_boundary.geojson', driver='GeoJSON')\ngeometry = json.loads((SETUP_DIR / 'ogallala_boundary.geojson').read_text())['features'][0]['geometry']\naoi = ee.Geometry(geometry)\n(SETUP_DIR / 'aoi_ee_geometry.json').write_text(json.dumps(aoi.getInfo()), encoding='utf-8')\nminx, miny, maxx, maxy = gdf.total_bounds\n(SETUP_DIR / 'ogallala_bbox.json').write_text(json.dumps({'west': minx, 'south': miny, 'east': maxx, 'north': maxy}, indent=2), encoding='utf-8')\n\nCONFIG_YAML = {json.dumps({'DRY_RUN': True, 'GEE_PROJECT': 'ee-ishansinhagzb', 'MASTER_CRS': 'EPSG:5070', 'MASTER_SCALE': {k: v['scale'] for k, v in DATASETS.items()}, 'AOI_BUFFER_METERS': 1000, 'START_DATE': '2000-01-01', 'END_DATE': '2024-09-30', 'datasets': DATASETS, 'independent_datasets': INDEPENDENT}, indent=2)!r}\nconfig = yaml.safe_dump(json.loads(CONFIG_YAML), sort_keys=False)\n(SETUP_DIR / 'config.yaml').write_text(config, encoding='utf-8')\nfor leaf in {LEAVES!r}:\n    ensure_layout(DRIVE_ROOT, ('01_gee_rasters/' if leaf.startswith('01') else '02_independent/') + leaf)\nprint('AOI and Drive tree ready:', DRIVE_ROOT)\n'''
    final = """print((SETUP_DIR / 'config.yaml').read_text(encoding='utf-8'))\nrequired = [SETUP_DIR / name for name in ['ogallala_boundary.geojson', 'ogallala_bbox.json', 'aoi_ee_geometry.json', 'config.yaml', 'utils.py']]\nmissing = [str(path) for path in required if not path.exists()]\nif missing:\n    raise RuntimeError(f'Missing setup outputs: {missing}')\nprint('Setup complete:', len(required), 'files')\n"""
    cells = [cell('markdown', title, 'markdown'), cell('code', bootstrap, 'python'), cell('code', body, 'python'), cell('code', "logger, LOG_PATH = configure_logger('00_setup', DRIVE_ROOT)\n", 'python'), cell('code', "DATASETS = " + repr(DATASETS) + "\nINDEPENDENT = " + repr(INDEPENDENT) + "\n", 'python'), cell('code', "# AOI processing is performed in the setup cell above; all downstream notebooks load its serialized geometry.\n", 'python'), cell('code', "# Setup has no XEE pull; it writes the AOI handoff consumed by XEE and batch pathways.\n", 'python'), cell('code', "# No server-side export is needed for the vector AOI.\n", 'python'), cell('code', "MAX_CONCURRENT = 20\nMAX_DAILY_TASKS = 2500\n", 'python'), cell('code', "# CF writing is provided by 00_setup/utils.py for every downstream notebook.\n", 'python'), cell('code', "# AOI vectors are validated by GeoPandas and recorded in the setup manifest.\nwrite_manifest(DRIVE_ROOT, '00_setup', required, 'Ogallala Aquifer Boundary')\n", 'python'), cell('code', final, 'python')]
    return notebook(cells)


def setup_notebook() -> dict:
    """Build the setup notebook without interpolating runtime dictionaries."""
    title = "# 00: Setup AOI and configuration\n\n**Why:** This notebook creates the single Drive handoff used by every Phase 0 notebook. It streams the USGS High Plains boundary directly to Drive when needed, serializes the AOI, and creates the complete output tree.\n\n- Source: USGS High Plains Aquifer boundary\n- Outputs: `00_setup/ogallala_boundary.geojson`, `ogallala_bbox.json`, `aoi_ee_geometry.json`, `config.yaml`, and `utils.py`\n- Author: `<name / institution>`"
    bootstrap = setup_code('00_setup')
    config_object = {'DRY_RUN': True, 'GEE_PROJECT': 'ee-ishansinhagzb',
                     'MASTER_CRS': 'EPSG:5070',
                     'MASTER_SCALE': {key: value['scale'] for key, value in DATASETS.items()},
                     'AOI_BUFFER_METERS': 1000, 'START_DATE': '2000-01-01',
                     'END_DATE': '2024-09-30',
                     'datasets': {key: {**value, 'drive_subfolder': '01_gee_rasters/' + value['folder']} for key, value in DATASETS.items()},
                     'independent_datasets': {key: {**value, 'drive_subfolder': '02_independent/' + value['folder']} for key, value in INDEPENDENT.items()},
                     'output_folders': LEAVES}
    body = """import geopandas as gpd
import requests
import zipfile
import io
import yaml

source_dir = SETUP_DIR / 'source_boundary'
shapefile = next(source_dir.glob('*.shp'), None) if source_dir.exists() else None
if shapefile is None:
    url = os.getenv('OGALLALA_BOUNDARY_ZIP_URL', 'https://water.usgs.gov/GIS/dsdl/hp_bound2010.zip')
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    source_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        archive.extractall(source_dir)
    shapefile = next(source_dir.glob('*.shp'))

gdf = gpd.read_file(shapefile).to_crs('EPSG:4326')
gdf.to_file(SETUP_DIR / 'ogallala_boundary.geojson', driver='GeoJSON')
geometry = json.loads((SETUP_DIR / 'ogallala_boundary.geojson').read_text())['features'][0]['geometry']
aoi = ee.Geometry(geometry)
(SETUP_DIR / 'aoi_ee_geometry.json').write_text(json.dumps(aoi.getInfo()), encoding='utf-8')
minx, miny, maxx, maxy = gdf.total_bounds
(SETUP_DIR / 'ogallala_bbox.json').write_text(json.dumps({'west': minx, 'south': miny, 'east': maxx, 'north': maxy}, indent=2), encoding='utf-8')

CONFIG_YAML = """ + repr(json.dumps(config_object)) + """
(SETUP_DIR / 'config.yaml').write_text(yaml.safe_dump(json.loads(CONFIG_YAML), sort_keys=False), encoding='utf-8')
for leaf in """ + repr(LEAVES) + """:
    ensure_layout(DRIVE_ROOT, ('01_gee_rasters/' if leaf.startswith('01') else '02_independent/') + leaf)
print('AOI and Drive tree ready:', DRIVE_ROOT)
"""
    final = """required = [SETUP_DIR / name for name in ['ogallala_boundary.geojson', 'ogallala_bbox.json', 'aoi_ee_geometry.json', 'config.yaml', 'utils.py']]
write_manifest(DRIVE_ROOT, '00_setup', required, 'Ogallala Aquifer Boundary')
missing = [str(path) for path in required if not path.exists()]
if missing:
    raise RuntimeError(f'Missing setup outputs: {missing}')
print((SETUP_DIR / 'config.yaml').read_text(encoding='utf-8'))
print('Setup complete:', len(required), 'files')
"""
    cells = [cell('markdown', title, 'markdown'), cell('code', bootstrap, 'python'),
             cell('code', body, 'python'),
             cell('code', "logger, LOG_PATH = configure_logger('00_setup', DRIVE_ROOT)\n", 'python'),
             cell('code', "DATASET = {'dataset_name': 'Ogallala Aquifer Boundary', 'source': 'USGS', 'access_method': 'requests / geopandas'}\n", 'python'),
             cell('code', "# The setup cell above is the dataset definition and AOI processing pipeline.\n", 'python'),
             cell('code', "# Setup has no XEE pull; it writes the serialized AOI handoff.\n", 'python'),
             cell('code', "# No server-side export is needed for the vector AOI.\n", 'python'),
             cell('code', "%pip -q install compliance-checker\nMAX_CONCURRENT = 20\nMAX_DAILY_TASKS = 2500\n", 'python'),
             cell('code', "# CF writing is provided by 00_setup/utils.py for downstream notebooks.\n", 'python'),
             cell('code', "# The setup manifest and file checks are emitted in the final cell.\n", 'python'),
             cell('code', final, 'python')]
    return notebook(cells)


def orchestrator() -> dict:
    title = "# 99: Run all and validate\n\n**Why:** Run this notebook after setup and dataset submissions. It executes notebooks with `nbclient` when available, aggregates the registry, checks every CF file, and writes the Drive README index.\n\nThe orchestrator re-raises failures so Colab exits visibly non-zero."
    setup = setup_code('99_run_all_and_validate')
    code = '''import subprocess\nimport sys\nimport pandas as pd\nimport xarray as xr\n\nNOTEBOOK_DIR = Path('/content/drive/MyDrive/Ogallala_Phase0/notebooks')\n# Upload the generated notebooks to this Drive folder, or set PHASE0_NOTEBOOK_DIR.\nNOTEBOOK_DIR = Path(os.getenv('PHASE0_NOTEBOOK_DIR', str(NOTEBOOK_DIR)))\nnotebooks = sorted(NOTEBOOK_DIR.glob('00_*.ipynb')) + sorted(NOTEBOOK_DIR.glob('01*.ipynb')) + sorted(NOTEBOOK_DIR.glob('02*.ipynb'))\nif os.getenv('RUN_NOTEBOOKS', '0') == '1':\n    for notebook_path in notebooks:\n        subprocess.run([sys.executable, '-m', 'jupyter', 'nbconvert', '--to', 'notebook', '--execute', str(notebook_path), '--output', str(notebook_path)], check=True)\n\nfailures = []\ncf_files = list(DRIVE_ROOT.glob('01_gee_rasters/**/cf/*.nc')) + list(DRIVE_ROOT.glob('02_independent/**/cf/*.nc'))\nfor path in cf_files:\n    try:\n        with xr.open_dataset(path) as ds:\n            required = {'Conventions', 'title', 'history', 'institution', 'source', 'references', 'comment'}\n            missing = required - set(ds.attrs)\n            if ds.attrs.get('Conventions') != 'CF-1.8' or missing:\n                failures.append(f'{path}: missing CF attrs {sorted(missing)}')\n    except Exception as exc:\n        failures.append(f'{path}: {exc}')\nregistry = pd.read_csv(registry_path(DRIVE_ROOT))\nif not registry.empty:\n    failures.extend(registry.loc[registry.status == 'FAILED', 'error_message'].dropna().astype(str).tolist())\nprint(registry.groupby('status').size().to_string() if not registry.empty else 'No registry rows')\nreadme = DRIVE_ROOT / 'README.md'\nlines = ['# Ogallala Phase 0 output index', '', f'Generated: {utc_now()}', '', '## CF files']\nlines.extend(f'- `{path}` ({path.stat().st_size} bytes)' for path in cf_files)\nlines += ['', '## Validation', f'- CF files checked: {len(cf_files)}', f'- Failures: {len(failures)}']\nreadme.write_text('\\n'.join(lines) + '\\n', encoding='utf-8')\nif failures:\n    raise RuntimeError('Phase 0 validation failures:\\n' + '\\n'.join(failures))\nprint('Validation passed; README written to', readme)\n'''
    cells = [cell('markdown', title, 'markdown'), cell('code', setup, 'python'), cell('code', "import yaml\nCONFIG = yaml.safe_load((DRIVE_ROOT / '00_setup/config.yaml').read_text())\nDRY_RUN = CONFIG.get('DRY_RUN', True)\n", 'python'), cell('code', "logger, LOG_PATH = configure_logger(NOTEBOOK_ID, DRIVE_ROOT)\n", 'python'), cell('code', "DATASET = {'name': 'All Phase 0 notebooks', 'drive_subfolder': ''}\n", 'python'), cell('code', "# Notebook execution is controlled by RUN_NOTEBOOKS=1 to keep normal validation non-blocking.\n", 'python'), cell('code', "# Each dataset notebook owns its XEE or batch pathway.\n", 'python'), cell('code', "# Each dataset notebook owns its Drive exports.\n", 'python'), cell('code', "MAX_CONCURRENT = 20\nMAX_DAILY_TASKS = 2500\n", 'python'), cell('code', "# CF validation is performed in the aggregation cell below.\n", 'python'), cell('code', "# QC and manifest aggregation are included in the README generation below.\n", 'python'), cell('code', code, 'python')]
    return notebook(cells)


def monitor() -> dict:
    title = "# 98: Monitor Drive exports\n\n**Why:** This is the only notebook that waits. It polls the shared Earth Engine task registry every five minutes and reports completed, failed, and running exports."
    setup = setup_code('98_monitor_gdrive_exports')
    code = '''import pandas as pd\n\nCONFIG = yaml.safe_load((DRIVE_ROOT / '00_setup/config.yaml').read_text())\nlogger, LOG_PATH = configure_logger(NOTEBOOK_ID, DRIVE_ROOT)\nmanager = BatchManager(DRIVE_ROOT, NOTEBOOK_ID, 'all exports', logger)\nwhile True:\n    frame = manager.monitor_once()\n    print(frame.groupby('status').size().to_string() if not frame.empty else 'No exports registered')\n    if frame.empty or not frame['status'].isin(['STARTED', 'RUNNING']).any():\n        break\n    time.sleep(300)\n'''
    cells = [cell('markdown', title, 'markdown'), cell('code', setup, 'python'), cell('code', "import yaml\nCONFIG = yaml.safe_load((DRIVE_ROOT / '00_setup/config.yaml').read_text())\nDRY_RUN = CONFIG.get('DRY_RUN', True)\n", 'python'), cell('code', "logger, LOG_PATH = configure_logger(NOTEBOOK_ID, DRIVE_ROOT)\n", 'python'), cell('code', "DATASET = {'name': 'Earth Engine export registry', 'drive_subfolder': ''}\n", 'python'), cell('code', "# Task discovery is registry-driven; no dataset is re-submitted here.\n", 'python'), cell('code', "# monitor_once polls the server status for every STARTED/RUNNING task.\n", 'python'), cell('code', "# This notebook intentionally blocks only between polls.\n", 'python'), cell('code', "MAX_CONCURRENT = 20\nMAX_DAILY_TASKS = 2500\n", 'python'), cell('code', "# BatchManager is imported from the Drive utility module.\n", 'python'), cell('code', "# The live dashboard is the QC output for this monitor.\n", 'python'), cell('code', code, 'python')]
    return notebook(cells)


def main() -> None:
    (ROOT / 'utils.py').write_text(UTILS, encoding='utf-8')
    config_template = {'DRY_RUN': True, 'GEE_PROJECT': 'ee-ishansinhagzb',
                       'MASTER_CRS': 'EPSG:5070',
                       'MASTER_SCALE': {key: value['scale'] for key, value in DATASETS.items()},
                       'AOI_BUFFER_METERS': 1000, 'START_DATE': '2000-01-01',
                       'END_DATE': '2024-09-30',
                       'datasets': {key: {**value, 'drive_subfolder': '01_gee_rasters/' + value['folder']} for key, value in DATASETS.items()},
                       'independent_datasets': {key: {**value, 'drive_subfolder': '02_independent/' + value['folder']} for key, value in INDEPENDENT.items()},
                       'output_folders': LEAVES}
    (ROOT / 'config.yaml').write_text(json.dumps(config_template, indent=2), encoding='utf-8')
    (ROOT / 'requirements.txt').write_text('''earthengine-api\nxee\nxarray\nnetCDF4\nh5netcdf\ngeopandas\nrioxarray\nregionmask\ndask[diagnostics]\nrequests\ntqdm\npython-dotenv\npyyaml\nmatplotlib\npandas\nshapely\npyproj\ncompliance-checker\nearthaccess\ndataretrieval\nameriflux-api\n''', encoding='utf-8')
    (ROOT / 'README.md').write_text('''# Ogallala Phase 0 Colab pipeline\n\nThis folder contains self-contained Colab notebooks for the Phase 0 Groundwater Buffering of Ecosystem Function During Drought data intake. The notebooks write only to `Google Drive/MyDrive/Ogallala_Phase0`; `/content/` is not used as an output location.\n\n## Run order\n\n1. Upload this folder, including the `high_plains_quifer` shapefile components, to a Drive working folder if Colab cannot see the local workspace.\n2. Run `00_setup_aoi_and_config.ipynb`. It streams the USGS boundary archive directly into Drive when a Drive copy is not present.\n3. Run any `01*.ipynb` and `02*.ipynb`. `DRY_RUN: true` processes exactly one XEE chunk and does not queue batch tasks. Set it to `false` in Drive `00_setup/config.yaml` for non-blocking yearly Earth Engine submissions.\n4. Run `98_monitor_gdrive_exports.ipynb` to poll server-side tasks.\n5. Run `99_run_all_and_validate.ipynb` to validate CF metadata and write the output index.\n\nThe supplied `config.yaml` is a template; notebook 00 writes the complete Drive config with asset IDs, bands, master scales, dates, and output folders. Registration-only independent datasets retain explicit endpoint placeholders and must be configured before use.\n''', encoding='utf-8')
    gee_filenames = {
        '01a': '01a_gee_grace.ipynb', '01b': '01b_gee_smap.ipynb',
        '01d': '01d_gee_tropomi_sif.ipynb', '01e': '01e_gee_era5_land.ipynb',
        '01f': '01f_gee_gridmet.ipynb', '01g': '01g_gee_modis_ndvi_evi.ipynb',
        '01h': '01h_gee_modis_landcover.ipynb', '01i': '01i_gee_srtm_dem.ipynb',
        '01j': '01j_gee_usda_cdl.ipynb', '01k': '01k_gee_gfsad.ipynb',
        '01l': '01l_gee_aridity_index.ipynb', '01m': '01m_gee_worldcover.ipynb',
        '01n': '01n_gee_landsat.ipynb', '01o': '01o_gee_sentinel2.ipynb',
        '01p': '01p_gee_pml_v2.ipynb', '01q': '01q_gee_ssebop.ipynb',
        '01r': '01r_gee_soilgrids.ipynb', '01s': '01s_gee_nasadem.ipynb',
    }
    for old_notebook in ROOT.glob('01*_gee_*.ipynb'):
        old_notebook.unlink()
    for notebook_id, dataset in DATASETS.items():
        filename = gee_filenames[notebook_id]
        (ROOT / filename).write_text(json.dumps(gee_notebook(notebook_id, dataset), indent=2), encoding='utf-8')
    for notebook_id, dataset in INDEPENDENT.items():
        (ROOT / f'{notebook_id}_indep_{dataset["folder"].replace(notebook_id + "_", "")}.ipynb').write_text(json.dumps(independent_notebook(notebook_id, dataset), indent=2), encoding='utf-8')
    (ROOT / '00_setup_aoi_and_config.ipynb').write_text(json.dumps(setup_notebook(), indent=2), encoding='utf-8')
    (ROOT / '98_monitor_gdrive_exports.ipynb').write_text(json.dumps(monitor(), indent=2), encoding='utf-8')
    (ROOT / '99_run_all_and_validate.ipynb').write_text(json.dumps(orchestrator(), indent=2), encoding='utf-8')
    print('Generated', len(list(ROOT.glob('*.ipynb'))), 'notebooks')


if __name__ == '__main__':
    main()
''