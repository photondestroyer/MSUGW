# Ogallala Phase 0 Colab pipeline

This folder contains self-contained Colab notebooks for the Phase 0 Groundwater Buffering of Ecosystem Function During Drought data intake. The notebooks write only to `Google Drive/MyDrive/Ogallala_Phase0`; `/content/` is not used as an output location.

## Run order

1. Upload this folder, including the `high_plains_quifer` shapefile components, to a Drive working folder if Colab cannot see the local workspace.
2. Run `00_setup_aoi_and_config.ipynb`. It streams the USGS boundary archive directly into Drive when a Drive copy is not present.
3. Run any `01*.ipynb` and `02*.ipynb`. `DRY_RUN: true` processes exactly one XEE chunk and does not queue batch tasks. Set it to `false` in Drive `00_setup/config.yaml` for non-blocking yearly Earth Engine submissions.
4. Run `98_monitor_gdrive_exports.ipynb` to poll server-side tasks.
5. Run `99_run_all_and_validate.ipynb` to validate CF metadata and write the output index.

The supplied `config.yaml` is a template; notebook 00 writes the complete Drive config with asset IDs, bands, master scales, dates, and output folders. Registration-only independent datasets retain explicit endpoint placeholders and must be configured before use.
