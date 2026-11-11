# Ogallala Aquifer Datasets Processing Pipeline

Satellite, reanalysis, model-assimilated and in-situ data over the Ogallala (High Plains) Aquifer for groundwater-buffering and drought-vulnerability analysis. The aquifer boundary is `HPA_polygon/hp_bound2010.*`, used for all clipping and masking. Source data are read-only and every product keeps its native grid and values with minimal change; each file is standalone, CF-compliant and checked for corruption. The corrected datasets are on [Google Drive](https://drive.google.com/drive/folders/1aIxTSmtsB1T7rNAPyNVH-J4UHs4IwEWN?usp=sharing).

## Folder structure

* `Colab notebooks/` holds the data processing pipelines. Run them on Google Colab to regenerate the data.
* `data audit/` holds the audit scripts. `DATA_AUDIT_REPORT_2026-10-20.md` lists what the last audit found.
* `merged_datasets/` holds the data files, with derived USGS products in `USGS data/derived_usgs/` and native USGS files in `USGS_HPA_native/`.
* `MERGED_DATASETS_CATALOG_2026-10-22.md` lists every variable in each dataset and flags redundant, partial and outdated files. It replaces `MERGED_DATASETS_CATALOG.md`, which is kept for history.
* `Ogallala_Project_Summary_2026-10-22.docx` is the one-page summary of the work done.
* `RESIDUAL_ET_MODELS.md` documents the models trained to get the residual ET. They are largely redundant now that more datasets are available. The code is in `residual_et/`.
* `TBI_DOCUMENTATION.md` covers the temporal buffering index analysis; `Groundwater_Buffering_Research_Plan*` is the research plan.
* `local data processing/` holds scripts that run locally (USGS wells and maps, gSSURGO, GOSIF). `landsat indices plot/` and `plots/` hold figures.

## Progress since September

* **Gridded meteorology and drought** now come from the Climatology Lab's own files on its native 1/24° grid (287 × 233), replacing the Earth Engine copies.
* **New datasets:** OCO-2 and OCO-3 SIF soundings, TROPOMI SIF (TROPOSIF and Caltech), ECHO-ET hourly ET, Copernicus GLO-30 DEM (rebuilt), LANID 30 m irrigation, MODIS land cover on its native grid, the USGS well and map files in native form, and a 2.3-million-record well water-level series. AmeriFlux tower data were removed.
* **New notebooks whose output is too large to keep here:** Landsat 30 m monthly indices, Landsat ESPA 4 km indices, USDA CDL, OpenET monthly ET, GEDI biomass and ECOSTRESS 70 m.
* **Audit fixes:** the 4 km USGS products were rebuilt on the gridMET lab grid. The notebooks for SSEBop, GRACE, irrigation, GFSAD, soil group, water-table depth, MODIS and SMAP were rewritten, but several files in `merged_datasets/` still come from before the fixes.

## Open items

* Re-run the rebuilt notebooks in Colab to replace the outdated files, and finish the partial exports (OCO-2, Caltech TROPOMI, Landsat 500 m).
* Several scripts and `residual_et/config.py` still point to removed files.
* Mask the undeclared -4.6e12 value in OCO SIF, and treat `eddi2y` for 1984-86 as missing (a defect in the provider's file).
