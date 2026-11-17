# Blinding audit (Task 0.4)

Earlier work in this repository that put ET next to an aquifer variable, found by scanning code, notebooks, documents and figure names. A row is a lead to read, not proof that an analysis was run.

## Known analyses

- **gw_buffering_xgb/**: SSEBop dekadal ET 2015-04 to 2022-05-21 (RESIDUAL_ET_MODELS.md; catalogue item 3); residual ET by water-table depth (buffering_by_wtd.png, effects.csv)
- **residual_et_store/**: feature store of the corrected residual-ET models; dekads 2015-04-01 to 2022-05-21
- **XGB model/**: XGBoost residual-ET script; water-table depth among the summaries
- **temporal buffering/**: temporal buffering index (TBI) scripts over USGS water-level change
- **tbi_maps_v3_usgs/**: TBI maps v3 (USGS water-level data)

## Files (41)

| file | kind | dates in file | year ranges | git first / last | overlaps 2016-2024 |
|---|---|---|---|---|---|
| `temporal buffering/temporal_buffering_index.py` | py | 2026-08-24 .. 2026-08-24 | 2001-2015, 2003-2022 | 2026-09-24 / 2026-09-24 | yes |
| `residual_et_store/full/manifest.json` | json | 2013-04-01 .. 2026-10-07 |  | untracked / untracked |  |
| `residual_et_store/full/prepare.log` | log | 2015-04-01 .. 2022-05-21 |  | untracked / untracked | yes |
| `residual_et_store/smoke/manifest.json` | json | 2013-04-01 .. 2026-10-07 |  | untracked / untracked |  |
| `XGB model/gw_buffering_xgb_et.py` | py | 2016-01-01 .. 2022-01-01 | 2001-2015, 2016-2021 | 2026-09-24 / 2026-09-24 | yes |
| `gw_buffering_xgb/run_log.txt` | txt | 2016-01-01 .. 2021-12-21 |  | untracked / untracked | yes |
| `Colab notebooks/Cropland_&_Ag_Mgmt.ipynb` | ipynb | 2010-01-01 .. 2010-01-01 | 2001-2015, 2007-2012 | 2026-09-24 / 2026-10-22 |  |
| `Colab notebooks/ET_&_Energy_Thermal_Flux.ipynb` | ipynb | 2003-01-01 .. 2022-05-21 |  | 2026-09-24 / 2026-10-23 | yes |
| `Colab notebooks/GridMET_Lab_Ogallala.ipynb` | ipynb |  ..  |  | 2026-10-07 / 2026-10-07 |  |
| `Colab notebooks/SMAP.ipynb` | ipynb | 2000-01-01 .. 2026-10-22 |  | 2026-09-24 / 2026-10-22 |  |
| `Colab notebooks/hydrology_&_soils.ipynb` | ipynb | 2002-04-03 .. 2024-09-30 |  | 2026-09-24 / 2026-10-22 | yes |
| `MERGED_DATASETS_CATALOG.md` | md | 2000-02-01 .. 2025-12-19 | 2001-2015, 2001-2024, 2002-2024, 2003-2022, 2005-2021, 2015-2024, 2017-2019, 2017-2021, 2017-2026 | 2026-09-24 / 2026-09-29 | yes |
| `MERGED_DATASETS_CATALOG_2026-10-22.md` | md | 1980-01-05 .. 2026-10-05 | 2000-2026, 2001-2015, 2001-2024, 2007-2008, 2008-2012, 2008-2016, 2013-2023, 2017-2019, 2017-2023, 2018-2020 | 2026-10-23 / 2026-10-23 | yes |
| `README.md` | md |  ..  |  | 2026-09-24 / 2026-11-11 |  |
| `RESIDUAL_ET_MODELS.md` | md | 1980-01-05 .. 2026-10-07 | 2001-2015, 2001-2024, 2003-2014, 2003-2021, 2010-2014, 2015-2022, 2017-2019 | 2026-10-07 / 2026-10-07 | yes |
| `TBI_DOCUMENTATION.md` | md | 1984-01-05 .. 2026-09-29 | 2001-2015, 2001-2024, 2015-2019, 2015-2022, 2016-2021, 2020-2021, 2021-2022 | 2026-09-24 / 2026-10-07 | yes |
| `docs/superpowers/specs/2026-11-11-paper1-v3-changes-design.md` | md | 2026-11-11 .. 2026-11-11 | 2011-2013, 2014-2016, 2018-2020 | 2026-11-11 / 2026-11-11 | yes |
| `docs/superpowers/specs/2026-11-14-pilot-design.md` | md | 2026-11-14 .. 2026-11-15 | 2011-2015, 2011-2017, 2013-2015, 2015-2024, 2016-2017, 2016-2020, 2016-2024, 2018-2020, 2021-2023 | untracked / untracked | yes |
| `in-this-directory-there-cozy-sonnet.md` | md | 2013-12-30 .. 2026-10-07 | 2001-2024, 2003-2014 | 2026-10-07 / 2026-10-07 | yes |
| `local data processing/smap_l4_extractor.py` | py | 2000-01-01 .. 2020-01-01 |  | 2026-09-24 / 2026-09-24 | yes |
| `mieno nature paper/Drought-Production-Risk-Aquifer-master/Codes/functions.R` | R |  ..  |  | untracked / untracked |  |
| `mieno nature paper/Drought-Production-Risk-Aquifer-master/Data/metadata.Rmd` | Rmd |  ..  |  | untracked / untracked |  |
| `mieno nature paper/metadata.Rmd` | Rmd |  ..  |  | untracked / untracked |  |
| `paper1/P1_colab_driver.ipynb` | ipynb | 2005-01-01 .. 2026-11-11 | 2015-2016, 2015-2026 | 2026-11-11 / 2026-11-17 | yes |
| `paper1/README.md` | md | 2026-11-11 .. 2026-11-11 | 2018-2020 | 2026-11-11 / 2026-11-11 | yes |
| `paper1/checks.py` | py |  ..  |  | 2026-11-11 / 2026-11-11 |  |
| `paper1/config.py` | py | 2026-11-04 .. 2026-11-11 | 2008-2014, 2015-2026 | 2026-11-11 / 2026-11-11 |  |
| `paper1/docs/CONFIG.md` | md | 2026-11-04 .. 2026-11-11 | 2008-2014, 2015-2026 | 2026-11-11 / 2026-11-11 |  |
| `paper1/p1_1/rules.py` | py | 2026-11-11 .. 2026-11-11 |  | 2026-11-11 / 2026-11-11 |  |
| `paper1/tests/synth.py` | py | 2004-01-01 .. 2014-01-01 | 2004-2013 | 2026-11-11 / 2026-11-11 |  |
| `paper1/tests/test_provenance.py` | py |  ..  |  | 2026-11-11 / 2026-11-11 |  |
| `plotting codes/plot_derived_usgs_spatial.ipynb` | ipynb | 2019-06-15 .. 2019-06-15 | 2017-2019 | 2026-09-24 / 2026-09-24 | yes |
| `plotting codes/plot_gw_clusters_yearly.py` | py |  ..  |  | 2026-10-07 / 2026-10-07 |  |
| `plotting codes/plot_gw_shallow_deep_clusters.py` | py |  ..  |  | 2026-10-07 / 2026-10-07 |  |
| `residual_et/buffering.py` | py |  ..  | 2017-2019 | 2026-10-07 / 2026-10-07 | yes |
| `residual_et/config.py` | py | 2013-04-01 .. 2026-10-07 |  | 2026-10-07 / 2026-10-21 |  |
| `residual_et/features.py` | py |  ..  |  | 2026-10-07 / 2026-10-07 |  |
| `residual_et/store.py` | py |  ..  |  | 2026-10-07 / 2026-10-07 |  |
| `residual_et/tests/synth.py` | py | 2013-04-01 .. 2022-06-01 |  | 2026-10-07 / 2026-10-07 | yes |
| `residual_et/tests/test_core.py` | py | 2014-06-01 .. 2022-06-01 |  | 2026-10-07 / 2026-10-07 | yes |
| `residual_et/train_dlstm.py` | py |  ..  |  | 2026-10-07 / 2026-10-21 |  |

Rule for the freeze: results from these analyses for 2016-2024 dry seasons are treated as already seen; the decision log records any choice that could have been influenced by them.