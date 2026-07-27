# Ogallala Phase 0 Project Guide

This project prepares Phase 0 datasets for the Groundwater Buffering of Ecosystem Function During Drought study. The notebooks run in a local Jupyter environment, authenticate Google Earth Engine locally, and submit server-side GeoTIFF exports to Google Drive.

## What the project does

- Uses the High Plains Aquifer shapefile in `high_plains_quifer/` as the study-area boundary.
- Creates a local setup directory containing the serialized AOI, bounding box, configuration, logs, manifests, and task registry.
- Preserves native GEE bands, QA bands, and uncertainty bands instead of selecting only a headline variable.
- Reprojects GEE exports to the configured master grid, normally `EPSG:5070`.
- Submits non-blocking yearly or static `ee.batch.Export.image.toDrive` tasks.
- Writes exported GeoTIFFs to Google Drive; the local notebooks do not download them or convert them to NetCDF.
- Tracks task IDs, status, retries, timestamps, and errors in the local task registry.
- Provides a monitor notebook for polling Earth Engine task status.

## Prerequisites

Install the local dependencies from [requirements.txt](requirements.txt):

```powershell
python -m pip install -r requirements.txt
```

Use a Python environment with Jupyter or VS Code Notebook support. Set the Earth Engine project before starting:

```powershell
$env:GEE_PROJECT = "ee-ishansinhagzb"
```

The default local output root is `./Ogallala_Phase0`. To choose another local location:

```powershell
$env:OGALLALA_PHASE0_ROOT = "G:\\MSU_GWB\\Ogallala_Phase0"
```

This local root is for configuration and logs. It is not a mounted Google Drive path.

## Run order

1. Open [00_setup_aoi_and_config.ipynb](00_setup_aoi_and_config.ipynb) from the project directory and run its cells in order. Authenticate Earth Engine when the browser opens.
2. Confirm that the local setup root contains `00_setup/aoi_ee_geometry.json`, `00_setup/config.yaml`, `logs/task_registry.csv`, and the dataset folders.
3. Run the required `01*.ipynb` GEE notebooks. Each notebook builds its clipped image or collection and submits GeoTIFF tasks to Google Drive.
4. Run [98_monitor_gdrive_exports.ipynb](98_monitor_gdrive_exports.ipynb) later to poll active Earth Engine tasks. This is the only notebook designed to wait.
5. Run [99_run_all_and_validate.ipynb](99_run_all_and_validate.ipynb) to summarize task status and write the local output index.
6. Use [try.ipynb](try.ipynb) as a small standalone GRACE export example when testing authentication and Drive permissions.

## Google Drive export behavior

The GEE notebooks call the Earth Engine batch API through `BatchManager`. A typical export has this shape:

```python
ee.batch.Export.image.toDrive(
    image=image,
    description=description,
    folder="Ogallala_Phase0/01_gee_rasters/01a_grace_mascons/raw",
    fileNamePrefix=file_name,
    region=AOI,
    scale=MASTER_SCALE,
    crs=MASTER_CRS,
    maxPixels=1e13,
    fileFormat="GeoTIFF",
)
```

The task is started and the notebook continues. It does not wait for Google Drive to receive the file. The monitor notebook can later update the task state to `COMPLETED` or `FAILED`.

The task registry is stored at:

```text
<OGALLALA_PHASE0_ROOT>/logs/task_registry.csv
```

The registry records `STARTED`, `RUNNING`, `COMPLETED`, `FAILED`, and `RETRIED` states. The batch manager limits locally registered active submissions to 20 and retries rate-limit failures with backoff.

## Configuration

[config.yaml](config.yaml) is the editable project template. Notebook 00 writes the active copy under the local setup root. Important settings include:

- `GEE_PROJECT`: Earth Engine Cloud project.
- `DRY_RUN`: controls whether GEE export submission is skipped.
- `MASTER_CRS`: target projection, normally `EPSG:5070`.
- `MASTER_SCALE`: 1000 m for coarse products and 30 m for fine products.
- `datasets`: GEE asset IDs, bands, dates, scales, and output folders.
- `independent_datasets`: external source metadata and registration-dependent endpoints.

## Safe dry run

For a safe suite-wide dry run, temporarily change the active setup configuration to:

```yaml
DRY_RUN: true
```

In local batch mode, this means the GEE notebooks build and inspect their configuration but do not call `task.start()`. No Earth Engine export task is submitted. Restore `DRY_RUN: false` before a production submission run.

The repository preflight command is:

```powershell
python validate_phase0_artifacts.py
```

This checks notebook JSON, cell metadata, Python syntax, required batch-export code, and forbidden Colab paths. It does not authenticate Earth Engine or contact external services.

## Notebook map

- `00_setup_aoi_and_config.ipynb`: AOI, local configuration, output tree, and shared utility setup.
- `01a` through `01s`: GEE raster and image-collection exports. SMAP L4 and L3 VOD are grouped in `01b`.
- `02a` through `02j`: independent datasets. These use their own local streaming or registered-client workflows because their sources are not GEE assets.
- `98_monitor_gdrive_exports.ipynb`: task monitoring and registry updates.
- `99_run_all_and_validate.ipynb`: aggregate status and validation.
- `try.ipynb`: minimal local GRACE-to-Google-Drive export test.

## Important limitations

- Earth Engine authentication, Drive permissions, quotas, and asset availability must be valid for the local user.
- Registration-only independent datasets require their endpoint, account, or API key to be configured before use.
- Exported GeoTIFFs remain in Google Drive. The project intentionally does not run `write_exported_geotiff_as_cf`.
- A successful local submission means that Earth Engine accepted the task, not that the server-side export has completed. Use notebook 98 for completion status.