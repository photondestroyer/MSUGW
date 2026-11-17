# Feasibility memo

Mode: stage3. Targets are planning targets: a value below target is reported and discussed and does not stop the work (pilot Section 3).

| question | item | observed | target | meets |
|---|---|---|---|---|
| 1 Precision | soil-moisture contrast: median 95% half-width (mm) | 2.94 | 15 | yes |
| 1 Precision | peak-demand contrast: median 95% half-width (mm) | 2.16 | 15 | yes |
| 1 Precision | independent dry episodes in the design | 12 | 3 | yes |
| 2 Thickness | cohort cells within 8 km of a well, worst year | 0.999 | 0.6 | yes |
| 2 Thickness | counties with >= 50 cells in both thin and thick thirds | 3 | 15 | below target |
| 2 Thickness | thin-to-thick difference / SE | 3.05 | 3 | yes |
| 2 Thickness | water-table RMSE, locations held out (m) | 7.92 |  |  |
| 3 Capacity | cells with a clear link | 0.667 |  |  |
| 3 Capacity | cells with a usable rate (driller yield) | 0.333 |  |  |
| 3 Capacity | cohort cell-years 2021-2024 with known irrigation status | 0.416 |  |  |
| 4 Satellite | dry county-years kept by the soil-moisture mask | 0.74 | 0.7 | yes |
| 4 Satellite | counties with >= 3 usable reference footprints | 0.865 | 0.8 | yes |
| 4 Satellite | peak-month cohort cells with >= 4 ET models | 0.974 | 0.8 | yes |

## Projected precision across scenarios

| model | contrast | dry_variance_factor | base_error | n | cells | counties | cluster_level | clusters | halfwidth_median_mm | halfwidth_p10_mm | halfwidth_p90_mm | target_mm | A_thin_m | A_thick_m | delta_x | reliability_median |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C1 | soil moisture | 1 | cell | 50873 | 7284 | 32 | county | 32 | 2.94 | 2.17 | 4.25 | 15 | 17.9 | 76.4 | -0.962 | 0.694 |
| P1 | peak demand | 1 | cell | 64183 | 7583 | 37 | county | 37 | 2.16 | 1.64 | 3.16 | 15 | 17.9 | 76.4 | 30.1 | 0.551 |
| C1 | soil moisture | 1 | district | 50873 | 7284 | 32 | county | 32 | 2.65 | 1.83 | 3.83 | 15 | 17.9 | 76.4 | -0.962 | 0.815 |
| P1 | peak demand | 1 | district | 64183 | 7583 | 37 | county | 37 | 1.99 | 1.42 | 2.74 | 15 | 17.9 | 76.4 | 30.1 | 0.6 |
| C1 | soil moisture | 1.5 | cell | 50873 | 7284 | 32 | county | 32 | 3.35 | 2.4 | 4.71 | 15 | 17.9 | 76.4 | -0.962 | 0.693 |
| P1 | peak demand | 1.5 | cell | 64183 | 7583 | 37 | county | 37 | 2.42 | 1.74 | 3.65 | 15 | 17.9 | 76.4 | 30.1 | 0.566 |
| C1 | soil moisture | 2 | cell | 50873 | 7284 | 32 | county | 32 | 3.94 | 2.64 | 5.26 | 15 | 17.9 | 76.4 | -0.962 | 0.692 |
| P1 | peak demand | 2 | cell | 64183 | 7583 | 37 | county | 37 | 2.75 | 1.98 | 3.99 | 15 | 17.9 | 76.4 | 30.1 | 0.57 |

## Open

- Reference footprints use a 150 km radius (D13); the SMAP file used is named in Task 1.2 (smap_source)
- Authorised rate, pumping hours and verified capacity: wimas_lookup_sheet.csv (Check 2.7)
- Footprint check against 2018-2020 labels (LANID): see check_2018_2020 in the Task 1.2 summary and the irr_max_2018_2020 column of footprints.csv
- Decisions still proposed: ['D30', 'D31']

## Next step (pilot Section 11)

Read the table against Section 11 of the pilot: which row of that table describes these findings decides the next step, with the supervisor.