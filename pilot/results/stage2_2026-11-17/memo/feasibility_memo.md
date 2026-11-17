# Feasibility memo

Mode: stage2. Targets are planning targets: a value below target is reported and discussed and does not stop the work (pilot Section 3).

| question | item | observed | target | meets |
|---|---|---|---|---|
| 1 Precision | soil-moisture contrast: median 95% half-width (mm) | 10.2 | 15 | yes |
| 1 Precision | peak-demand contrast: median 95% half-width (mm) | 6.83 | 15 | yes |
| 1 Precision | independent dry episodes in the design | 3 | 3 | yes |
| 2 Thickness | cohort cells within 8 km of a well, worst year | 1 | 0.6 | yes |
| 2 Thickness | counties with >= 50 cells in both thin and thick thirds | 0 | 15 | below target |
| 2 Thickness | thin-to-thick difference / SE | 2.33 | 3 | below target |
| 2 Thickness | water-table RMSE, locations held out (m) | 7.92 |  |  |
| 3 Capacity | cells with a clear link | 0.86 |  |  |
| 3 Capacity | cells with a usable rate (driller yield) | 0.36 |  |  |
| 3 Capacity | cohort cell-years 2021-2024 with known irrigation status | 0.416 |  |  |
| 4 Satellite | dry county-years kept by the soil-moisture mask | 0.74 | 0.7 | yes |
| 4 Satellite | counties with >= 3 usable reference footprints | 0.865 | 0.8 | yes |
| 4 Satellite | peak-month cohort cells with >= 4 ET models | 0.972 | 0.8 | yes |

## Projected precision across scenarios

| model | contrast | dry_variance_factor | base_error | n | cells | counties | cluster_level | clusters | halfwidth_median_mm | halfwidth_p10_mm | halfwidth_p90_mm | target_mm | A_thin_m | A_thick_m | delta_x | reliability_median |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C1 | soil moisture | 1 | cell | 1283 | 186 | 4 | 10 km block | 69 | 10.2 | 7.84 | 13.4 | 15 | 31.2 | 74.4 | -1.18 | 0.758 |
| P1 | peak demand | 1 | cell | 1599 | 186 | 4 | 10 km block | 69 | 6.83 | 5.34 | 8.94 | 15 | 31.2 | 74.4 | 38 | 0.755 |
| C1 | soil moisture | 1 | district | 1283 | 186 | 4 | 10 km block | 69 | 8.19 | 6.45 | 10.9 | 15 | 31.2 | 74.4 | -1.18 | 0.959 |
| P1 | peak demand | 1 | district | 1599 | 186 | 4 | 10 km block | 69 | 5.43 | 4.18 | 7.24 | 15 | 31.2 | 74.4 | 38 | 0.942 |
| C1 | soil moisture | 1.5 | cell | 1283 | 186 | 4 | 10 km block | 69 | 11.9 | 8.95 | 14.8 | 15 | 31.2 | 74.4 | -1.18 | 0.756 |
| P1 | peak demand | 1.5 | cell | 1599 | 186 | 4 | 10 km block | 69 | 7.75 | 6.23 | 9.84 | 15 | 31.2 | 74.4 | 38 | 0.756 |
| C1 | soil moisture | 2 | cell | 1283 | 186 | 4 | 10 km block | 69 | 13 | 10.5 | 17.7 | 15 | 31.2 | 74.4 | -1.18 | 0.755 |
| P1 | peak demand | 2 | cell | 1599 | 186 | 4 | 10 km block | 69 | 8.69 | 6.82 | 11.2 | 15 | 31.2 | 74.4 | 38 | 0.756 |

## Open

- Reference footprints use a 150 km radius (D13); the SMAP file used is named in Task 1.2 (smap_source)
- Authorised rate, pumping hours and verified capacity: wimas_lookup_sheet.csv (Check 2.7)
- Footprint check against 2018-2020 labels (LANID): see check_2018_2020 in the Task 1.2 summary and the irr_max_2018_2020 column of footprints.csv
- Decisions still proposed: []

## Next step (pilot Section 11)

Read the table against Section 11 of the pilot: which row of that table describes these findings decides the next step, with the supervisor.