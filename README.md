# Ogallala Aquifer Datasets Processing Pipeline

This repository assembles satellite, reanalysis, model-assimilated and
in-situ observations over the Ogallala (High Plains) Aquifer in support of
groundwater-buffering and drought-vulnerability analysis. Gridded holdings
share EPSG:4326; the ~4 km gridMET grid is the common analysis grid. Source
files are read-only; all derived products are written to new files.

## Data processing

Most datasets were processed with Google Colab notebooks using Google Earth
Engine exports, covering: gridMET meteorology and drought indices, MODIS
SSEBop evapotranspiration, MODIS vegetation indices (MOD13A3) and land cover
(MCD12Q1), SMAP Level-4 soil moisture and Level-3 VOD, GRACE/GRACE-FO storage,
OCO-2 solar-induced fluorescence, GFSAD cropland, global irrigation areas,
HiHydroSoil groups, HRES water-table depth, the GDE map, AmeriFlux tower
records and the GLO digital elevation model. The merged outputs reside in
`merged_datasets/`; the export notebooks are retained under
`Colab notebooks/`.

A subset was processed locally (SSM Python environment): the USGS well
networks (water-level change grids, hydraulic conductivity and the
groundwater-wells compilation), the GOSIF 8-day/monthly SIF product and its
extraction mesh, the gSSURGO soil stack and plant-available-water grids, the
clipped GDE layers, and the drought-index distribution fits. These derived
products reside in `merged_datasets/USGS data/derived_usgs/`.

## Documentation

- `MERGED_DATASETS_CATALOG.md` — file-by-file inventory of `merged_datasets/`:
  variables, temporal extents and sampled quality assessments.
- `TBI_DOCUMENTATION.md` — methods, assumptions and results of the temporal
  buffering index analysis.
- `Groundwater_Buffering_Research_Plan.docx` — overarching research plan.
- `AmeriFlux_NEON_data_readme.pdf` — Flux tower data usage policy and definitions. (data is processed)

## Provenance notes

Shapefile composite `HPA_polygon/hp_bound2010.*` defines the aquifer boundary
used for all clipping and masking. Please refer to [[MERGED_DATASETS_CATALOG]] for detailed description of each of the data files. Each of the file are standalone CF 1.8 compliant, and every care has been taken to export data in its native resolution with minimum changes with appropriate checks for data corruption so that there is no problem down the line.  The processed datasets are available at [google drive link](https://drive.google.com/drive/folders/1aIxTSmtsB1T7rNAPyNVH-J4UHs4IwEWN?usp=sharing).
