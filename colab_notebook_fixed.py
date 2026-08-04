#!/usr/bin/env python3
"""
colab_notebook_fixed.py — Robust Earth Engine → NetCDF exporter via XEE
========================================================================

Drop-in replacement for the original colab_notebook.ipynb.
Change DATASET_ID below and run on Google Colab.

Bugs fixed from the original notebook
--------------------------------------
1. Missing ``import logging`` → NameError at runtime before any data processing.

2. CRS pixel-size leaking into CF ``scale_factor`` encoding.
   XEE stores the CRS transform as variable attributes.  The 30 m pixel
   resolution ended up as ``scale_factor: 30.0`` in the NetCDF encoding.
   On read-back xarray multiplied every raw value by 30, destroying data.

3. Custom WKT CRS causing silent grid misalignment.
   The CDL dataset's native CRS is a non-standard "IMAGINE GeoTIFF" WKT.
   ``pyproj`` can parse it, but when XEE forwards it to Earth Engine the
   resulting grid does not overlap the actual raster tiles → all NaN.
   Fix: resolve to a standard EPSG code before building the grid.

4. No data validation before writing multi-GiB files.
   The script now loads a small spatial sample and verifies that at least
   some pixel values are finite before committing to an expensive write.

5. No compression — files were 4× larger than necessary.
   Enabled zlib compression by default.

6. No retry / back-off for transient Earth Engine errors.
   Network hiccups, quota throttles, or server 503s now trigger automatic
   retries with exponential back-off.

7. GeoPandas ``union_all`` / ``unary_union`` compatibility.
   ``unary_union`` was deprecated in GeoPandas ≥ 1.0; ``union_all()`` does
   not exist in older versions.  Both are handled.

8. Over-complex ROI geometry can choke Earth Engine.
   If the shapefile has > 50 000 vertices the geometry is simplified
   automatically before being sent to GEE.

9. SR-ORG projection codes (e.g. MODIS Sinusoidal) are unresolvable by
   pyproj.  The script detects these and falls back gracefully.

10. Temp-file cleanup on failure — partial / corrupt files are removed.
"""

# ╔══════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — edit this section for each dataset run         ║
# ╚══════════════════════════════════════════════════════════════════╝

DATASET_ID      = "USDA/NASS/CDL"            # GEE ImageCollection ID
PROJECT_ID      = "msugw-503806"             # GEE Cloud project
SHAPEFILE_DIR   = "/content/ogallala_shp"    # Folder containing .shp + sidecar files
DRIVE_ROOT      = "/content/drive/MyDrive/MSUGWB"
ROI_LABEL       = "Ogallala"                 # Suffix for output filenames
FORCE_OVERWRITE = False                      # Re-download files that already exist?
YEAR_START      = None                       # None → auto-detect from collection
YEAR_END        = None                       # None → auto-detect from collection
FALLBACK_CRS    = "EPSG:4326"                # Used when native CRS has no EPSG match
DASK_WORKERS    = 4                          # Keep ≤ 6 to respect GEE rate limits
CHUNK_XY        = 512                        # Dask spatial chunk size (pixels)
COMPRESS_LEVEL  = 4                          # zlib level (0 = off, 9 = max)
MAX_RETRIES     = 3                          # Per-year retry limit for transient errors
RETRY_BASE_SEC  = 30                         # Base back-off between retries (doubles)

# ╔══════════════════════════════════════════════════════════════════╗
# ║  IMPORTS                                                        ║
# ╚══════════════════════════════════════════════════════════════════╝

import os
import sys
import glob
import time
import shutil
import logging
import tempfile

import numpy as np
import ee
import xarray as xr
import xee                      # noqa: F401  (registers the 'ee' engine)
from xee import helpers
import geopandas as gpd
import pyproj
import dask
import requests

# ╔══════════════════════════════════════════════════════════════════╗
# ║  LOGGING & ENVIRONMENT                                         ║
# ╚══════════════════════════════════════════════════════════════════╝

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s │ %(levelname)-7s │ %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("ee_exporter")

# Suppress noisy HTTP-pool chatter
logging.getLogger("urllib3.connectionpool").setLevel(logging.ERROR)

# Larger connection pool to avoid dropped connections on long runs
_adapter = requests.adapters.HTTPAdapter(pool_connections=50, pool_maxsize=50)
_session = requests.Session()
_session.mount("https://", _adapter)

# Limit Dask parallelism so GEE isn't hammered with concurrent requests
dask.config.set(scheduler="threads", num_workers=DASK_WORKERS)


# ╔══════════════════════════════════════════════════════════════════╗
# ║  HELPER FUNCTIONS                                               ║
# ╚══════════════════════════════════════════════════════════════════╝


def mount_drive() -> None:
    """Mount Google Drive when running on Colab; skip gracefully otherwise."""
    try:
        from google.colab import drive          # type: ignore[import-untyped]
        drive.mount("/content/drive")
        log.info("Google Drive mounted.")
    except ImportError:
        log.info("Not running on Colab — Drive mount skipped.")


def init_earth_engine(project_id: str) -> None:
    """Authenticate and initialise the Earth Engine API."""
    ee.Authenticate()
    ee.Initialize(project=project_id)
    log.info("Earth Engine initialised  (project=%s)", project_id)


# ── ROI loading ──────────────────────────────────────────────────


def _count_coords(geom) -> int:
    """Return the total number of coordinate tuples in *geom*."""
    try:
        import shapely
        return int(shapely.get_num_coordinates(geom))
    except (ImportError, AttributeError):
        # Rough heuristic when shapely ≥ 2.0 is unavailable
        return len(geom.wkt) // 20


def load_roi(shapefile_dir: str):
    """
    Load a shapefile, ensure EPSG:4326, dissolve all features into a
    single geometry, and simplify if the vertex count is too high for
    Earth Engine.

    Returns
    -------
    roi_shapely : shapely.geometry.base.BaseGeometry
        Dissolved geometry in EPSG:4326.
    roi_ee : ee.Geometry
        The same geometry as an Earth Engine object.
    """
    shp_files = glob.glob(os.path.join(shapefile_dir, "**", "*.shp"), recursive=True)
    if not shp_files:
        raise FileNotFoundError(f"No .shp file in '{shapefile_dir}'")

    shp_path = shp_files[0]
    log.info("Loading shapefile: %s", shp_path)
    gdf = gpd.read_file(shp_path)

    # Reproject to WGS-84 if needed
    if gdf.crs is None or gdf.crs.to_epsg() != 4326:
        log.info("Reprojecting shapefile → EPSG:4326 …")
        gdf = gdf.to_crs(epsg=4326)

    # Dissolve — handle both old and new GeoPandas API
    roi_shapely = (
        gdf.geometry.union_all()
        if hasattr(gdf.geometry, "union_all")
        else gdf.geometry.unary_union
    )

    # Simplify very complex polygons to keep within EE limits
    n_coords = _count_coords(roi_shapely)
    if n_coords > 50_000:
        roi_shapely = roi_shapely.simplify(tolerance=0.001, preserve_topology=True)
        log.warning(
            "Simplified ROI: %s → %s vertices (tolerance=0.001°)",
            f"{n_coords:,}",
            f"{_count_coords(roi_shapely):,}",
        )

    roi_ee = ee.Geometry(roi_shapely.__geo_interface__)
    return roi_shapely, roi_ee


# ── CRS / scale detection ───────────────────────────────────────


def resolve_crs_and_scale(collection: ee.ImageCollection):
    """
    Auto-detect native CRS and pixel scale from a GEE ImageCollection,
    and resolve them to a standard EPSG code + scale in CRS units.

    Why this matters
    ~~~~~~~~~~~~~~~~
    Datasets like USDA CDL use a custom WKT CRS that ``pyproj`` can parse
    but that Earth Engine cannot map back to its internal projection.  The
    resulting grid ends up mis-aligned → XEE returns all NaN.

    MODIS datasets use ``SR-ORG:6974`` (Sinusoidal) which ``pyproj``
    cannot parse at all.

    In both cases we must resolve to a well-known EPSG code and compute
    the pixel scale in that code's native units.

    Returns
    -------
    crs : str
        An ``"EPSG:XXXX"`` string that XEE/EE can reliably consume.
    scale : float
        Pixel resolution in the units of *crs* (metres for projected,
        degrees for geographic).
    """
    first_img = collection.first()
    proj = first_img.select(0).projection()
    proj_info = proj.getInfo()

    native_crs_raw = proj_info.get("crs", "EPSG:4326")
    transform = proj_info.get("transform", [1, 0, 0, 0, -1, 0])
    native_scale_crs_units = abs(transform[0])

    # nominalScale() always returns metres, regardless of the CRS
    nominal_scale_m = proj.nominalScale().getInfo()

    log.info(
        "Native CRS : %s",
        native_crs_raw if len(native_crs_raw) < 60 else native_crs_raw[:57] + "…",
    )
    log.info("Native pixel scale : %.6g  (CRS units)", native_scale_crs_units)
    log.info("Nominal scale      : %.1f m", nominal_scale_m)

    resolved_crs = FALLBACK_CRS  # default fallback

    # ── SR-ORG codes → pyproj cannot parse, skip straight to fallback ──
    if native_crs_raw.startswith("SR-ORG:"):
        log.warning("SR-ORG CRS detected (%s) — using fallback %s", native_crs_raw, FALLBACK_CRS)

    else:
        # ── Attempt EPSG resolution via pyproj ──
        try:
            crs_obj = pyproj.CRS.from_user_input(native_crs_raw)
            epsg = crs_obj.to_epsg()

            if epsg:
                resolved_crs = f"EPSG:{epsg}"
                log.info("Resolved to EPSG:%d", epsg)

            elif crs_obj.axis_info and crs_obj.axis_info[0].unit_name == "metre":
                # Metres-based but no exact EPSG match (e.g. ERDAS custom Albers).
                # EPSG:5070 is the standard CONUS Albers Equal-Area and is
                # parameter-identical to the custom CDL projection.
                resolved_crs = "EPSG:5070"
                log.warning(
                    "CRS uses metres but has no EPSG code → using EPSG:5070 (CONUS Albers)"
                )

            else:
                log.warning("CRS has no EPSG code and is not metres-based → using %s", FALLBACK_CRS)

        except Exception as exc:
            log.warning("pyproj could not parse CRS: %s → using %s", exc, FALLBACK_CRS)

    # ── Compute scale in the resolved CRS's units ──
    try:
        resolved_crs_obj = pyproj.CRS.from_user_input(resolved_crs)
        unit = resolved_crs_obj.axis_info[0].unit_name
    except Exception:
        unit = "degree"

    if unit == "metre":
        scale = round(nominal_scale_m, 2)
    else:
        # Geographic CRS → convert metres to approximate degrees
        # 1° latitude ≈ 111 320 m
        scale = nominal_scale_m / 111_320.0

    log.info("Resolved CRS : %s  |  scale : %.6g %s", resolved_crs, scale, unit)
    return resolved_crs, scale


# ── Encoding cleanup ────────────────────────────────────────────


def strip_xee_encoding(ds: xr.Dataset) -> xr.Dataset:
    """
    Remove **all** encoding injected by the XEE engine from the dataset.

    Why
    ---
    XEE stores CRS metadata — including the ``crs_transform`` pixel size —
    as variable attributes.  When a 30 m resolution dataset is opened, the
    number ``30`` ends up as ``scale_factor: 30.0`` in the variable's
    encoding dict.  ``xr.Dataset.to_netcdf()`` then faithfully writes this
    as a CF ``scale_factor`` attribute.  On subsequent read-back, xarray
    **multiplies every value by 30** (or divides — depending on the version),
    destroying the data.

    The ``_FillValue`` is similarly set to ``np.nan`` (float32), which
    interacts badly with integer-typed data and with the spurious scale.

    By clearing all encoding we ensure that ``to_netcdf()`` writes exactly
    the values that are in memory, with no hidden transformations.
    """
    for name in list(ds.data_vars) + list(ds.coords):
        if name in ds:
            ds[name].encoding.clear()
    return ds


def build_clean_encoding(ds: xr.Dataset, compress_level: int = 4) -> dict:
    """
    Build an explicit, minimal CF-compliant encoding dict for ``to_netcdf()``.

    * Uses ``float32`` + ``_FillValue=NaN`` for all data variables.
      This is safe for both continuous and categorical data, avoids type-
      casting pitfalls, and keeps the pipeline simple.
    * Enables ``zlib`` compression, which is especially effective for rasters
      with large masked (NaN) regions outside the ROI boundary.
    * Sets clean dtypes for coordinate variables.
    """
    encoding: dict = {}

    for var in ds.data_vars:
        encoding[var] = {
            "dtype": "float32",
            "_FillValue": np.float32(np.nan),
            "zlib": compress_level > 0,
            "complevel": compress_level,
        }

    # Coordinates: store with full precision and no fill sentinel
    for coord in ds.coords:
        if coord == "time":
            encoding[coord] = {"dtype": "float64", "_FillValue": None}
        elif coord in ("x", "y", "latitude", "longitude", "lat", "lon"):
            encoding[coord] = {"dtype": "float64", "_FillValue": None}

    return encoding


# ── Data validation ─────────────────────────────────────────────


def validate_data(ds: xr.Dataset, band_name: str, year: int) -> None:
    """
    Load a small window from the **centre** of the grid and verify that
    at least some pixels are finite.

    Raises ``RuntimeError`` (non-retryable) if validation fails.
    """
    da = ds[band_name]

    # Identify spatial dimension names (varies across datasets)
    y_dim = next((d for d in ("y", "latitude", "lat") if d in da.dims), None)
    x_dim = next((d for d in ("x", "longitude", "lon") if d in da.dims), None)
    if y_dim is None or x_dim is None:
        log.warning("[%d] Cannot identify spatial dims for validation — skipping.", year)
        return

    ny, nx = da.sizes[y_dim], da.sizes[x_dim]

    # 200×200 pixel window from the spatial centre
    yc, xc = ny // 2, nx // 2
    window = 100

    indexers: dict = {}
    if "time" in da.dims:
        indexers["time"] = 0
    indexers[y_dim] = slice(max(0, yc - window), min(ny, yc + window))
    indexers[x_dim] = slice(max(0, xc - window), min(nx, xc + window))

    sample = da.isel(**indexers).compute()
    values = sample.values
    n_finite = int(np.isfinite(values).sum())
    n_total = values.size

    log.info("[%d] Validation: %s/%s finite in centre sample", year, f"{n_finite:,}", f"{n_total:,}")

    if n_finite == 0:
        raise RuntimeError(
            f"[{year}] VALIDATION FAILED — centre sample is entirely NaN.  "
            "Earth Engine returned no valid pixels.  Check CRS alignment, "
            "geometry overlap, and data availability for this year."
        )

    pct = 100.0 * n_finite / n_total
    if pct < 1.0:
        log.warning(
            "[%d] Only %.1f%% of centre sample is finite — possible partial coverage.", year, pct
        )


# ╔══════════════════════════════════════════════════════════════════╗
# ║  MAIN PIPELINE                                                  ║
# ╚══════════════════════════════════════════════════════════════════╝


def main() -> None:
    log.info("═" * 62)
    log.info("  Earth Engine → NetCDF Exporter (fixed)")
    log.info("  Dataset : %s", DATASET_ID)
    log.info("═" * 62)

    # ── 1. Environment ──────────────────────────────────────────
    mount_drive()
    init_earth_engine(PROJECT_ID)

    # ── 2. Load ROI ─────────────────────────────────────────────
    roi_shapely, _roi_ee = load_roi(SHAPEFILE_DIR)

    # ── 3. Output folder ────────────────────────────────────────
    folder_name = DATASET_ID.replace("/", "_")    # e.g. "USDA_NASS_CDL"
    short_name  = DATASET_ID.split("/")[-1]       # e.g. "CDL"
    drive_folder = os.path.join(DRIVE_ROOT, short_name)
    os.makedirs(drive_folder, exist_ok=True)
    log.info("Output folder : %s", drive_folder)

    # ── 4. Collection metadata ──────────────────────────────────
    full_collection = ee.ImageCollection(DATASET_ID)
    target_variables = full_collection.first().bandNames().getInfo()

    start_year = YEAR_START or int(
        ee.Date(
            full_collection.sort("system:time_start").first().get("system:time_start")
        ).get("year").getInfo()
    )
    end_year = YEAR_END or int(
        ee.Date(
            full_collection.sort("system:time_start", False).first().get("system:time_start")
        ).get("year").getInfo()
    )
    years = list(range(start_year, end_year + 1))

    log.info("Bands (%d) : %s", len(target_variables), target_variables)
    log.info("Years      : %d → %d  (%d total)", start_year, end_year, len(years))

    # ── 5. Resolve CRS & scale ──────────────────────────────────
    grid_crs, pixel_scale = resolve_crs_and_scale(
        full_collection.select(target_variables)
    )
    grid_scale = (pixel_scale, -pixel_scale)

    # ── 6. Build grid via fit_geometry ──────────────────────────
    #   fit_geometry reprojects the EPSG:4326 ROI into the target CRS
    #   and snaps it to a regular pixel grid.
    grid_params = helpers.fit_geometry(
        geometry=roi_shapely,
        geometry_crs="EPSG:4326",
        grid_crs=grid_crs,
        grid_scale=grid_scale,
    )
    log.info("Grid shape : %s", grid_params.get("shape_2d", "?"))

    # ── 7. Year-by-year export ──────────────────────────────────
    ok, fail, skip = 0, 0, 0

    for year in years:
        output_path = os.path.join(drive_folder, f"{short_name}_{year}_{ROI_LABEL}.nc")

        if os.path.exists(output_path) and not FORCE_OVERWRITE:
            log.info("[%d] Exists — skipping (%s)", year, output_path)
            skip += 1
            continue

        local_tmp = os.path.join(
            tempfile.gettempdir(), f"{short_name}_{year}_{ROI_LABEL}.nc"
        )

        succeeded = False
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                log.info("[%d] Attempt %d/%d …", year, attempt, MAX_RETRIES)

                # ── Filter collection to this year ──
                year_col = (
                    ee.ImageCollection(DATASET_ID)
                    .filter(ee.Filter.calendarRange(year, year, "year"))
                    .select(target_variables)
                    .sort("system:time_start")
                )

                n_images = year_col.size().getInfo()
                if n_images == 0:
                    log.warning("[%d] No images found — skipping year.", year)
                    break                         # not a failure, just no data
                log.info("[%d] %d image(s) in collection", year, n_images)

                # ── Open via XEE ──
                ds = xr.open_dataset(
                    year_col,
                    engine="ee",
                    chunks={"x": CHUNK_XY, "y": CHUNK_XY},
                    **grid_params,
                )

                # ── FIX 2 & 5: strip bogus encoding, build clean one ──
                ds = strip_xee_encoding(ds)
                encoding = build_clean_encoding(ds, COMPRESS_LEVEL)

                # ── FIX 4: validate before expensive write ──
                validate_data(ds, target_variables[0], year)

                # ── Write to local temp file ──
                log.info("[%d] Writing NetCDF → %s …", year, local_tmp)
                t0 = time.time()
                ds.to_netcdf(local_tmp, engine="netcdf4", encoding=encoding)
                elapsed = time.time() - t0

                size_bytes = os.path.getsize(local_tmp)
                size_mib = size_bytes / (1024 ** 2)
                log.info("[%d] Written in %.0fs — %.1f MiB", year, elapsed, size_mib)

                if size_bytes < 1000:
                    raise RuntimeError(
                        f"[{year}] File suspiciously small ({size_bytes} B)"
                    )

                # ── Move to Drive ──
                shutil.move(local_tmp, output_path)
                log.info("[%d] ✓ Saved → %s", year, output_path)
                ok += 1
                succeeded = True
                break  # exit retry loop

            except RuntimeError:
                # Validation failures and file-size checks are non-retryable
                log.error("[%d] Non-retryable error — see above.", year)
                fail += 1
                break

            except Exception as exc:
                log.error("[%d] Attempt %d failed: %s", year, attempt, exc)
                if attempt < MAX_RETRIES:
                    backoff = RETRY_BASE_SEC * (2 ** (attempt - 1))
                    log.info("[%d] Retrying in %ds …", year, backoff)
                    time.sleep(backoff)
                else:
                    log.error("[%d] Exhausted %d retries.", year, MAX_RETRIES)
                    fail += 1

            finally:
                # Clean up temp file on failure to avoid leaving stale GiB files
                if not succeeded and os.path.exists(local_tmp):
                    try:
                        os.remove(local_tmp)
                    except OSError:
                        pass

    # ── Summary ─────────────────────────────────────────────────
    log.info("═" * 62)
    log.info(
        "  DONE — %d succeeded · %d failed · %d skipped  (of %d years)",
        ok, fail, skip, len(years),
    )
    log.info("═" * 62)


# ╔══════════════════════════════════════════════════════════════════╗
# ║  ENTRY POINT                                                    ║
# ╚══════════════════════════════════════════════════════════════════╝

if __name__ == "__main__":
    main()
