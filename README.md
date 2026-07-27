# Ogallala Phase 0 local Earth Engine pipeline

These notebooks run in a local Jupyter environment. They authenticate Earth Engine locally and submit non-blocking `ee.batch.Export.image.toDrive` tasks. GeoTIFF outputs are written by Earth Engine to Google Drive; local files contain configuration, logs, manifests, and the task registry. No notebook mounts Google Drive or uses `/content/drive`.

## Run order

1. Install [requirements.txt](requirements.txt) in the local Python environment.
2. Run `00_setup_aoi_and_config.ipynb` from this directory. It uses the supplied `high_plains_quifer` shapefile when available and writes local setup files under `OGALLALA_PHASE0_ROOT` (default `./Ogallala_Phase0`).
3. Run any `01*.ipynb` and `02*.ipynb`. The GEE notebooks submit yearly/static GeoTIFF tasks to the configured Google Drive folder and do not wait for completion.
4. Run `98_monitor_gdrive_exports.ipynb` to poll Earth Engine tasks and update the local registry.
5. Run `99_run_all_and_validate.ipynb` to aggregate task status and write the local README index.

Set `GEE_PROJECT` and optionally `OGALLALA_PHASE0_ROOT` before starting. Set `EE_DRIVE_FOLDER` for the standalone `try.ipynb` example. The notebooks do not run `write_exported_geotiff_as_cf`; exported GeoTIFFs remain in Google Drive. Registration-only independent datasets retain explicit endpoint placeholders and must be configured before use.
