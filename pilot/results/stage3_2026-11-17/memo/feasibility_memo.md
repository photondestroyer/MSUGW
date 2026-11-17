# Feasibility memo

Mode: stage3. Targets are planning targets: a value below target is reported and discussed and does not stop the work (pilot Section 3).

| question | item | observed | target | meets |
|---|---|---|---|---|
| 1 Precision | soil-moisture contrast: median 95% half-width (mm) | 2.77 | 15 | yes |
| 1 Precision | peak-demand contrast: median 95% half-width (mm) | 1.67 | 15 | yes |
| 1 Precision | independent dry episodes in the design | 12 | 3 | yes |
| 2 Thickness | cohort cells within 8 km of a well, worst year | 0.999 | 0.6 | yes |
| 2 Thickness | counties with >= 50 cells in both thin and thick thirds | 1 | 15 | below target |
| 2 Thickness | thin-to-thick difference / SE | 2.69 | 3 | below target |
| 2 Thickness | water-table RMSE, locations held out (m) | 7.92 |  |  |
| 3 Capacity | cells with a clear link | 0.714 |  |  |
| 3 Capacity | cells with a usable rate (driller yield) | 0.238 |  |  |
| 3 Capacity | cohort cell-years 2021-2024 with known irrigation status | 0.416 |  |  |
| 4 Satellite | dry county-years kept by the soil-moisture mask | 0.74 | 0.7 | yes |
| 4 Satellite | counties with >= 3 usable reference footprints | 0.865 | 0.8 | yes |
| 4 Satellite | peak-month cohort cells with >= 4 ET models | 0.974 | 0.8 | yes |

## Projected precision across scenarios

| model | contrast | dry_variance_factor | base_error | n | cells | counties | cluster_level | clusters | halfwidth_median_mm | halfwidth_p10_mm | halfwidth_p90_mm | target_mm | A_thin_m | A_thick_m | delta_x | reliability_median |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C1 | soil moisture | 1 | cell | 42902 | 6131 | 31 | county | 31 | 2.77 | 1.91 | 4.62 | 15 | 23.5 | 75.4 | -0.943 | 0.659 |
| P1 | peak demand | 1 | cell | 53336 | 6299 | 35 | county | 35 | 1.67 | 1.17 | 2.62 | 15 | 23.5 | 75.4 | 29.6 | 0.63 |
| C1 | soil moisture | 1 | district | 42902 | 6131 | 31 | county | 31 | 1.88 | 1.32 | 2.81 | 15 | 23.5 | 75.4 | -0.943 | 0.959 |
| P1 | peak demand | 1 | district | 53336 | 6299 | 35 | county | 35 | 1.18 | 0.817 | 1.95 | 15 | 23.5 | 75.4 | 29.6 | 0.924 |
| C1 | soil moisture | 1.5 | cell | 42902 | 6131 | 31 | county | 31 | 3.23 | 2.07 | 4.56 | 15 | 23.5 | 75.4 | -0.943 | 0.657 |
| P1 | peak demand | 1.5 | cell | 53336 | 6299 | 35 | county | 35 | 1.95 | 1.34 | 3.17 | 15 | 23.5 | 75.4 | 29.6 | 0.634 |
| C1 | soil moisture | 2 | cell | 42902 | 6131 | 31 | county | 31 | 3.3 | 2.33 | 4.92 | 15 | 23.5 | 75.4 | -0.943 | 0.661 |
| P1 | peak demand | 2 | cell | 53336 | 6299 | 35 | county | 35 | 2.26 | 1.46 | 3.69 | 15 | 23.5 | 75.4 | 29.6 | 0.632 |

## Open

- Reference footprints use a 150 km radius (D13); the SMAP file used is named in Task 1.2 (smap_source)
- Authorised rate, pumping hours and verified capacity: wimas_lookup_sheet.csv (Check 2.7)
- Footprint check against 2018-2020 labels (LANID): see check_2018_2020 in the Task 1.2 summary and the irr_max_2018_2020 column of footprints.csv
- Decisions still proposed: []

## Next step (pilot Section 11)

Read the table against Section 11 of the pilot: which row of that table describes these findings decides the next step, with the supervisor.