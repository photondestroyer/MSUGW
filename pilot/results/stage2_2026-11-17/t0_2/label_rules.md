# Pilot definitions (written label rules)

Generated from `pilot/definitions.py` and `pilot/decisions.yaml` by `python -m pilot docs`. Each term is one function; every task calls it, none re-implements it.

## Cell

> A quarter-section of about 65 ha.

- Built in: Task 0.2a
- Functions: `definitions.build_cells()`, `definitions.cell_regular()`
- D03 (fixed) Cells: `{"source": "BLM CadNSDI First Division sections", "split": "four quarters by the section's own axes", "min_ha": 55, "max_ha": 75}`
  - `build_cells`: Sections (GeoDataFrame, any projected CRS) -> cells: each section split into 4 quarters along its own axes (common/cells.split_quarters). Adds cell_id = <section id><NE|NW|SW|SE>, area_ha, regular.
  - `cell_regular`: Cell: "A quarter-section of about 65 ha." A quarter whose area lies outside D03's range is kept in the cell table but flagged irregular (short sections, correction lines) and excluded from cohort and rainfed.

## Cohort

> Cells irrigated in at least four of the five years 2011-2015, with corn in at least two of 2013-2015, and not served by surface water. The cohort is fixed.

- Built in: Task 0.2c
- Functions: `definitions.map_share()`, `definitions.map_roles()`, `definitions.is_irrigated_year()`, `definitions.is_corn_year()`, `definitions.maps_agree()`, `definitions.served_by_surface_water()`, `definitions.is_cohort()`
- D04 (fixed) Irrigated cell-year: `{"min_share": 0.5}`
- D05 (fixed) Corn cell-year: `{"min_share": 0.5, "cdl_codes": [1, 12, 13]}`
- D06 (fixed) Primary irrigation map and disagreement: `{"primary": "LANID", "check": "AIM-HPA", "min_agree_years": 3, "min_map_coverage": 0.9}`
- D09 (fixed) Surface-water flag: `{"radius_km": 1.6, "sources": ["S"], "irrigation_use": "any AF_USED_IRR or ACRES > 0 in 1990-2023"}`
  - `map_share`: Irrigated share of a cell in one map and year: irrigated pixels / labelled pixels. Unknown (NaN) when the map labels fewer than D06 min_map_coverage of the cell's pixels, so a cell cut by the edge of a map (LANID holds labels only inside the aquifer polygon) is never labelled from a sliver. Returns (share, coverage).
  - `map_roles`: D06: which irrigation map is primary and which is the check, as the column prefixes of the cell table ('lanid' or 'aim'). Returns (primary, check).
  - `is_irrigated_year`: Irrigated cell-year: irrigated share of the cell's mapped pixels >= D04 min_share. NaN -> False.
  - `is_corn_year`: Corn cell-year: share of the cell's pixels in CDL corn >= D05 min_share. NaN -> False.
  - `maps_agree`: D06: per year, do the primary and the check map give the same irrigated / not-irrigated label (D04 on each)? Returns (n_years_agreeing, agree_ok). A year with either map unknown counts as not agreeing.
  - `served_by_surface_water`: "Not served by surface water" (D09): a WIMAS surface-water point of diversion with irrigation use within D09 radius_km of the cell centre. nearest_sw_km: distance to the nearest such point (inf when none).
  - `is_cohort`: Cohort: "Cells irrigated in at least four of the five years 2011-2015, with corn in at least two of 2013-2015, and not served by surface water. The cohort is fixed: a cell stays in it if it later stops irrigating or changes crop." primary_share, check_share: (n, 5) yearly irrigated shares 2011-2015 of the primary map (the one that decides "irrigated", D06: LANID) and of the check map (AIM-HPA; used only for the agreement test); corn_share: (n, 3) for 2013-2015; surface_water, regular: (n,) bool. Returns (cohort bool (n,), reason str array; '' for cohort cells). The first failing reason is given, in the order of COHORT_REASONS. Later years are never read, which is what keeps the cohort fixed.

## Rainfed cells

> Cropland cells in the same county with no irrigation in 2011-2015, matched to the cohort on crop history and root-zone water storage.

- Built in: Task 0.2c
- Functions: `definitions.is_rainfed()`, `definitions.irrigated_later()`, `definitions.dominant_class()`, `definitions.crop_class()`, `definitions.match_rainfed()`
- D07 (fixed) Rainfed cell: `{"max_irr_share": 0.02, "min_cropland_share": 0.5}`
- D08 (fixed) Crop history and rainfed matching: `{"dominant_min_share": 0.5, "dominant_min_years": 2, "rzws_tol_mm": 25}`
  - `is_rainfed`: Rainfed cells: "Cropland cells in the same county with no irrigation in 2011-2015" (matching is match_rainfed). No irrigation = at most D07 max_irr_share in EVERY year in BOTH maps (both must be known); cropland = at least D07 min_cropland_share in every year. Returns (rainfed bool, reason).
  - `irrigated_later`: "Check whether any were irrigated later": a cell with more than D07 max_irr_share irrigated in any later year of either map. later_shares: (n, years, maps); a map without a year (AIM-HPA ends in 2017) is NaN there and does not count. The flag is reported; it never moves a cell between groups (the groups are fixed from 2011-2015).
  - `dominant_class`: Per cell-year: the crop class covering >= D08 dominant_min_share of the cell, else ''. class_shares: (n, years, k) shares of the k classes in `names`.
  - `crop_class`: Crop history of 2013-2015: the class dominant in at least D08 dominant_min_years of the years, else 'rotation' (or 'none' when no year has a dominant crop).
  - `match_rainfed`: Rainfed matching: a rainfed cell is matched when a cohort cell of the SAME county and crop class has root-zone storage within D08 rzws_tol_mm. Returns (rainfed_matched bool, cohort_n_matches int) per cell.

## Deficit

> County reference ET minus precipitation from gridMET, monthly and May-September sums.

- Built in: Task 1.1
- Functions: `definitions.deficit()`
- D01 (fixed) Reference ET for the deficit: `{"variable": "pet", "name": "grass reference ET (gridMET ETo)"}`
  - `deficit`: Deficit: "County reference ET minus precipitation from gridMET" (D01 picks grass ETo), in mm.

## Peak month

> The month among June, July and August with the largest deficit in that county and year.

- Built in: Task 1.1
- Functions: `definitions.peak_month()`
  - `peak_month`: Peak month: "The month among June, July and August with the largest deficit in that county and year. It is chosen from weather only." Ties go to the earlier month. NaN in any month -> 0 (undefined).

## Dry episode

> A drought that covers neighbouring counties in the same period counts as one episode.

- Built in: Task 1.1
- Functions: `definitions.dry_episodes()`
- D10 (fixed) Dry episode: `{"within": "one season", "link": "counties sharing a border", "multiyear": "report alongside"}`
  - `dry_episodes`: Dry episode: "A drought that covers neighbouring counties in the same period counts as one episode." Dry county-years of one year that are linked through shared borders form one episode (D10). neighbours: {county: set of bordering counties}. Returns (episode id per row, -1 where not dry; multi-year id that also merges consecutive dry years sharing a county, reported alongside).

## Dry conditions

> Seasons above the 70th percentile of seasonal deficit and months below the 30th percentile of S; a month with missing S is dry.

- Built in: Task dry
- Functions: `definitions.season_dry()`, `definitions.month_dry()`
- D20 (fixed) Dry conditions (blinding): `{"season_pct": 70, "s_pct": 30, "missing_s_is_dry": true, "months": [5, 6, 7, 8, 9]}`
  - `season_dry`: Dry season: "seasons above the 70th percentile of seasonal deficit", the percentile taken over all Kansas county-years passed in. Returns (dry bool, threshold).
  - `month_dry`: Dry month: "months below the 30th percentile of the soil-moisture state S ... Treat a month with missing S as dry." Percentile over all Kansas county-months passed in. Returns (dry bool, threshold).

## Aquifer thickness, A

> Pre-season water-table elevation from recovered levels minus the aquifer base elevation, with its uncertainty. Thirds cut once over the cohort, for display.

- Built in: Task 2.4
- Functions: `definitions.thickness()`, `definitions.thickness_thirds()`, `definitions.assign_third()`, `definitions.well_support()`, `definitions.thickness_poor()`
- D22 (fixed) Wells and water levels: `{"kgs": "winter table", "usgs_months": [1, 2, 3, 4], "usgs_status": ["Static"], "keep_code_prefixes": ["11", "12"], "blank_code": "keep inside KGS HPA extent", "land_surface": "GLO30 DEM at the well"}`
- D23 (fixed) Water-table interpolation: `{"method": "regression kriging", "drift": "land-surface elevation", "variogram": "exponential", "neighbours": 32, "max_km": 40, "per": "year"}`
- D24 (fixed) Aquifer base and its uncertainty: `{"base": "KGS section hpa_bedrock", "se": "1.4826 x MAD of (USGS bottom - KGS base) by GMD"}`
- D25 (fixed) Poor thickness estimate: `{"max_rel_se": 0.5, "max_well_km": 8}`
- D26 (fixed) Thickness thirds: `{"true": "cohort cell mean A over 2016-2024", "quantiles": [0.3333333, 0.6666667]}`
  - `thickness`: Aquifer thickness, A: "Pre-season (January-April) water-table elevation, from recovered (static) water levels, minus the aquifer base elevation, with its uncertainty." The two errors are independent.
  - `thickness_thirds`: "Thirds of A (thin, middle, thick) are cut once over the whole cohort and used for display" (D26). Returns the two cut points.
  - `assign_third`: thin (A <= first cut), middle, thick (A > second cut); '' where A is unknown. Cuts from thickness_thirds.
  - `well_support`: A cell-year has a measured well close enough to know its water table: nearest well of that pre-season within D25 max_well_km (the pilot's 8 km support distance).
  - `thickness_poor`: Poor thickness estimate (D25): SE above D25 max_rel_se of A, or nearest measured well beyond max_well_km, or A not positive / unknown.

## Soil-moisture state, S

> SMAP SPL3SMP_E averaged over reference footprints, as an anomaly standardised by calendar month.

- Built in: Task 1.3
- Functions: `definitions.footprint_value()`, `definitions.anomaly()`, `definitions.soil_moisture_state()`
- D11 (fixed) Footprint value: `{"grid": "EASE-Grid 2.0 36 km", "min_nested_9km": 8}`
- D15 (fixed) Footprint weights: `{"weights": "equal"}`
- D16 (fixed) SMAP overpass: `{"overpass": "am"}`
- D17 (fixed) SMAP quality flag: `{"retrieval_qual_flag_bit0": 0}`
- D18 (fixed) SMAP climatology: `{"years": [2015, 2024], "min_years": 5, "anomaly": "z by footprint and calendar month"}`
- D19 (fixed) Valid county month of S: `{"min_supported_footprints": 2}`
  - `footprint_value`: A footprint's soil moisture on a date: mean of its valid nested 9 km retrievals, NaN unless at least D11 min_nested_9km of the 16 are valid. values: (..., 16).
  - `anomaly`: Standardised by calendar month: (value - mean) / SD of that footprint and calendar month.
  - `soil_moisture_state`: Soil-moisture state, S: "SMAP Level-3 soil moisture (SPL3SMP_E) averaged over reference footprints, as an anomaly standardised by calendar month." For one county-month: equal weights (D15) over its SUPPORTED reference footprints; NaN unless at least D19 min_supported_footprints. Returns (S, n_footprints).

## Reference footprint

> A SMAP footprint of about 36 km with less than 5% irrigated area, within 100 km of the county centroid, with similar rainfall and soils.

- Built in: Task 1.2
- Functions: `definitions.is_reference_footprint()`, `definitions.irrigation_free()`, `definitions.footprint_observed()`, `definitions.footprint_radius()`
- D11 (fixed) Footprint value: `{"grid": "EASE-Grid 2.0 36 km", "min_nested_9km": 8}`
- D12 (fixed) Footprint irrigation: `{"map": "AIM-HPA", "years": [2011, 2015], "max_irr_share": 0.05, "min_label_coverage": 0.9}`
- D13 (fixed) Footprint radius and its one widening: `{"radius_km": 100, "fallback_km": 150, "min_footprints": 3, "min_share_counties": 0.8}`
- D14 (fixed) Similar rainfall and soils: `{"normal_years": [1991, 2020], "rain_tol_frac": 0.2, "rzws_tol_frac": 0.3}`
  - `is_reference_footprint`: Reference footprint: "A SMAP footprint of about 36 km with less than 5% irrigated area, within 100 km of the county centroid, with similar rainfall and soils." For one county; arrays over footprints. irr_max: largest yearly irrigated share 2011-2015 (D12); label_cov: share of the footprint the label map covers; rain_ratio / rzws_ratio: footprint / county normal; rzws_known: soils cover enough of it. Returns (bool, reason).
  - `irrigation_free`: The irrigation part of the reference-footprint rule (D12): largest yearly irrigated share under max_irr_share, with the label map covering at least min_label_coverage of the footprint.
  - `footprint_observed`: The SMAP file holds enough of a footprint to ever give it a value: at least D11 min_nested_9km of its 16 nested 9 km cells lie on the file's grid with data.
  - `footprint_radius`: D13: the radius used for every county. n_ref_at_radius: {radius_km: array of reference-footprint counts per county}. 100 km unless fewer than 80% of counties reach 3, then 150 km (the one allowed widening).

## Supported month

> Valid retrievals on at least eight different dates, at least one in each third of the month; footprints are not combined.

- Built in: Task 1.3
- Functions: `definitions.is_supported_month()`, `definitions.month_third()`
- D16 (fixed) SMAP overpass: `{"overpass": "am"}`
- D17 (fixed) SMAP quality flag: `{"retrieval_qual_flag_bit0": 0}`
  - `is_supported_month`: Supported month: "For one footprint: valid retrievals on at least eight different dates, with at least one in each third of the month. Count a date once when it has both a morning and an evening retrieval. Do not combine footprints to reach eight dates." valid_days: day-of-month numbers of ONE footprint (duplicates allowed). Returns (supported, n_dates, thirds_ok).
  - `month_third`: 1 for days 1-10, 2 for 11-20, 3 for 21-end of month.

## Supported ET cell-month

> At least four models cover the cell.

- Built in: Task 1.4
- Functions: `definitions.et_model_ok()`, `definitions.et_supported()`
- D21 (fixed) ET model coverage of a cell-month: `{"min_valid_share": 0.8, "min_models": 4}`
- D29 (fixed) ET value of a cell-month: `{"variable": "et_ensemble_mad", "models": "sensitivity", "cell_value": "mean of valid pixels when D21 holds"}`
  - `et_model_ok`: A model counts for a cell-month when at least D21 min_valid_share of the cell's pixels hold a value with a positive observation count.
  - `et_supported`: Supported ET cell-month: at least D21 min_models models count (pilot target: four).

## Sample masks

> Soil-moisture contrast, peak-demand contrast, common sample, monthly timing.

- Built in: Task 1.7
- Functions: `definitions.sample_masks()`
  - `sample_masks`: The four sample masks of Task 1.7, per cell-month row (May-September): mask_sm      soil-moisture contrast: peak month, ET supported, county S valid, thickness usable (cohort) mask_peak    peak-demand contrast:   peak month, ET supported, thickness usable (cohort) mask_common  in both mask_timing  monthly timing model:   every season month, ET supported, thickness usable (cohort) Rainfed rows follow the same rules without the thickness condition (they give the baseline).

## Contrasts and models

> Thin-minus-thick ET at dry S (or high peak deficit) minus the same at the median; models C1, P1, J1, T1.

- Built in: Task 2.6, 2.9
- Functions: `definitions.contrast_sm()`, `definitions.contrast_peak()`
- D28 (fixed) Precision simulation: `{"n_sim": 200, "dry_variance_factors": [1.5, 2.0], "seed": 20261114}`
  - `contrast_sm`: Soil-moisture contrast: thin-minus-thick peak-month ET at dry S (10th percentile) minus the same at median S. In a model with b * A * S that is b * (A_thin - A_thick) * (S_p10 - S_p50).
  - `contrast_peak`: Peak-demand contrast: thin-minus-thick peak-month ET at a high peak-month deficit (90th percentile) minus the same at the median: b * (A_thin - A_thick) * (D_p90 - D_p50).
