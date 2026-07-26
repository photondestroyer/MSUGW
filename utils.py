"""Shared utilities for the Ogallala Phase 0 Colab notebooks."""
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
