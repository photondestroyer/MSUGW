# Pilot results: main_2026-11-17

Mode `main` (smoke = one county, Haskell 20081: a code test with proposed decisions, not pilot results). Copied from `H:/MSU_GWB/pilot_out/_main`; each folder holds one task's tables (.csv), notes (.md), figures (.png) and sidecars (.meta.json: units, sample mask, commit, decision-log hash).

| task | title | signature | files | summary |
|---|---|---|---|---|
| t0_1 | 0.1 Manifest | 58210a89 | 0 | not run |
| t0_4 | 0.4 Blinding audit | e12d648e | 0 | not run |
| counties | Study counties (D02) | b6041c15 | 4 | {"kansas": 105, "study": 37, "replication": {"texas_high_plains": 49, "sw_nebraska": 13}, "by_gmd": {"1": 5, "2": 4, "3": 12, "4": 9, "5": 7}} |
| t0_2 | 0.2 Cells and masks | e6c19347 | 9 | {"cells": 128156, "cohort": 7840, "rainfed": 17118, "rainfed_matched": 233, "rainfed_corn": 233, "rainfed_corn_scarce": true, "regular_share": 0.9848700021848372, "primary_map": "lanid", "check_map": "aim", "cohort_by_check_map_alone": 7697, "rainfed_irrigated_later": 3396, "crop_classes_cohort": {" ... |
| t0_3 | 0.3 Benchmark reproduction | 7bd97766 | 0 | not run |
| t1_1 | 1.1 Weather | a6d1e0b3 | 9 | {"threshold_mm": 634.2988845022332, "county_years": 333, "dry_county_years": 100, "episodes": 12, "episodes_multiyear": 7, "peak_month_share": {"6": 0.3783783783783784, "7": 0.35135135135135137, "8": 0.2702702702702703}, "missing_months": 0} |
| t1_2 | 1.2 Reference footprints | 815f0eef | 13 | {"radius_km": 150, "share_counties_3_at_100km": 0.6486486486486487, "footprints": 305, "reference": 41, "shared": 39, "smap_source": "pilot", "counties_with_3": 0.8918918918918919, "check_2018_2020": {"status": "done", "map": "LANID 2018-2020 (CONUS release)", "reference_footprints": 41, "fail_2018_ ... |
| t1_3 | 1.3 SMAP dates | 88b423d8 | 9 | {"smap_source": "pilot", "file_ends": "2026-11-14", "footprints": 41, "supported_share": 0.8595121951219512, "county_months_valid": 0.8126126126126126, "gaps_7d": [{"start": "2019-06-19", "end": "2019-07-24", "days": "36"}, {"start": "2022-08-06", "end": "2022-09-23", "days": "49"}]} |
| dry | Dry flags (blinding) | ff522694 | 3 | {"thresholds": {"season_threshold_mm": 634.2988845022332, "s_threshold": -0.39861567598702285}, "county_months": 1665, "dry_share": 0.572972972972973, "dry_season_share": 0.3003003003003003} |
| t1_5 | 1.5 Well inventory | cd03b5eb | 12 | {"well_years": 38447, "usable": 34288, "wells": 3145, "by_flag": {"": 34288, "aquifer code not HPA": 3460, "no DEM": 687, "above land surface": 12}, "usgs": {"usgs_records_ks": 76588, "jan_apr": 40221, "static": 5322, "duplicates_of_kgs": 5293}, "cohort_share_within_8km": {"2016": 1.0, "2017": 1.0,  ... |
| sample | 1.6 Test counties and cells (D27) | 003d28bc | 5 | {"counties": ["20069", "20185", "20055", "20175"], "cells": 300, "by_role": {"cohort": 200, "rainfed": 100}, "seed": 20261114} |
| t1_4 | 1.4 ET coverage | f182fe92 | 8 | {"cell_months": 5661405, "cohort_supported": 0.9529081632653061, "cells_without_openet": 2347, "parts": 81} |
| et_cells | Blinded ET extraction | dd3ad238 | 2 | {"rows": 1123110, "parts": 81, "gate_last": {"rows_in": 1123110, "rows_kept": 1123110, "cohort_dry_dropped": 0, "other_dropped": 0}, "gate_parts": {"rows_in": 1134945, "rows_kept": 1134945, "cohort_dry_dropped": 0, "other_dropped": 0}, "cohort_rows": 352800, "rainfed_rows": 770310, "audit_tables": 8 ... |
| t1_6 | 1.6 Processing test log | 7f5704a4 | 0 | not run |
| t1_inv | Preliminary county-year inventory | a03bf825 | 0 | not run |
| t1_7 | 1.7 Complete support table | bd9e733a | 13 | {"cell_months": 1123110, "masks": {"mask_sm": {"cohort_cells": 7483, "county_years": 244, "dry_episodes": 11}, "mask_peak": {"cohort_cells": 7784, "county_years": 332, "dry_episodes": 12}, "mask_common": {"cohort_cells": 7483, "county_years": 244, "dry_episodes": 11}, "mask_timing": {"cohort_cells": ... |
| rep | Replication region screen (metadata) | ee5fd700 | 0 | not run |
| c2_1 | 2.1 SMAP coverage in dry conditions | dc0a46ac | 0 | not run |
| c2_2 | 2.2 Reference footprints | a5e68abd | 0 | not run |
| c2_3 | 2.3 ET coverage | 720a62b4 | 0 | not run |
| c2_4 | 2.4 Aquifer thickness | 5d859d38 | 21 | {"years": [2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024], "cuts_m": [28.632600418040514, 50.13157604556917], "targets": [{"target": "cohort cells within 8 km of a measured well, worst year", "target_value": 0.6, "observed": 0.9988520408163265}, {"target": "count ... |
| c2_5 | 2.5 Soil-moisture axis | 2d068830 | 0 | not run |
| c2_6 | 2.6 Precision | cc8be7e3 | 0 | not run |
| c2_7 | 2.7 Well capacity | ee249db6 | 0 | not run |
| c2_8 | 2.8 Irrigation status after 2020 | 7ad7b108 | 6 | {"share_known": 0.41597576530612246, "by_year": {"2021": {"irrigated": 0.508, "not irrigated": 0.043, "unknown": 0.449}, "2022": {"irrigated": 0.514, "not irrigated": 0.042, "unknown": 0.443}, "2023": {"irrigated": 0.51, "not irrigated": 0.046, "unknown": 0.444}, "2024": {"irrigated": 0.0, "not irri ... |
| c2_9 | 2.9 Model support | 6339109a | 0 | not run |
| memo | Feasibility memo | 019ca256 | 0 | not run |
| paper | Paper items (pilot drafts) | 10f30c41 | 0 | not run |
| fig_thick | ET by aquifer thickness, simulated (Mieno Figs 3-5) | c2b5ec36 | 0 | not run |
| m_models | Main: models C1, P1, J1, T1 and the two contrasts | a9da36a4 | 8 | {"freeze": "pilot-freeze-2026-11-17", "rows": 65938, "primary": [{"contrast": "soil moisture", "value_mm": -6.083381233283657, "lo_mm": -12.299121150733855, "hi_mm": 0.13235868416654029}, {"contrast": "peak demand", "value_mm": 0.11407340497617408, "lo_mm": -3.3937599689522466, "hi_mm": 3.6219067789 ... |
| m_figs | Main: Figure 1, Figure 3, ET by aquifer thickness (Mieno Figs 3-5) | 6714b40b | 19 | {"freeze": "pilot-freeze-2026-11-17", "fig_thickness": {"rows": 64183, "cells": 7583, "counties": 37, "clusters": 37, "cluster_level": "county", "year_effects": true, "levels_m": [10.0, 40.0, 70.0, 100.0], "A_range_m": [0.011509403688023667, 155.29959667024752], "deficit_range_mm": [-125.53873689994 ... |
