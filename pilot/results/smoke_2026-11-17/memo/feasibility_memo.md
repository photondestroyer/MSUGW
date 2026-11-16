# Feasibility memo

Mode: smoke. Targets are planning targets: a value below target is reported and discussed and does not stop the work (pilot Section 3).

| question | item | observed | target | meets |
|---|---|---|---|---|
| 1 Precision | soil-moisture contrast: median 95% half-width (mm) | 7.19 | 15 | yes |
| 1 Precision | peak-demand contrast: median 95% half-width (mm) | 17.8 | 15 | below target |
| 1 Precision | independent dry episodes in the design | 3 | 3 | yes |
| 2 Thickness | cohort cells within 8 km of a well, worst year | 1 | 0.6 | yes |
| 2 Thickness | counties with >= 50 cells in both thin and thick thirds | 1 | 15 | below target |
| 2 Thickness | thin-to-thick difference / SE | 1.44 | 3 | below target |
| 2 Thickness | water-table RMSE, locations held out (m) | 7.92 |  |  |
| 3 Capacity | cells with a clear link | 0.857 |  |  |
| 3 Capacity | cells with a usable rate (driller yield) | 0.381 |  |  |
| 3 Capacity | cohort cell-years 2021-2024 with known irrigation status | 0.342 |  |  |
| 4 Satellite | dry county-years kept by the soil-moisture mask | 0.74 | 0.7 | yes |
| 4 Satellite | counties with >= 3 usable reference footprints | 0.865 | 0.8 | yes |
| 4 Satellite | peak-month cohort cells with >= 4 ET models | 1 | 0.8 | yes |

## Projected precision across scenarios

| model | contrast | dry_variance_factor | n | cells | counties | cluster_level | clusters | halfwidth_median_mm | halfwidth_p10_mm | halfwidth_p90_mm | target_mm | A_thin_m | A_thick_m | delta_x | reliability_median |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C1 | soil moisture | 1 | 3694 | 468 | 1 | 10 km block | 24 | 7.19 | 5.39 | 9.95 | 15 | 38 | 70.2 | -0.503 | 0.537 |
| P1 | peak demand | 1 | 4155 | 468 | 1 | 10 km block | 24 | 17.8 | 11.6 | 26.3 | 15 | 38 | 70.2 | 75.6 | 0.391 |
| C1 | soil moisture | 1.5 | 3694 | 468 | 1 | 10 km block | 24 | 7.98 | 6.36 | 11.7 | 15 | 38 | 70.2 | -0.503 | 0.54 |
| P1 | peak demand | 1.5 | 4155 | 468 | 1 | 10 km block | 24 | 18.2 | 14.7 | 22.4 | 15 | 38 | 70.2 | 75.6 | 0.394 |
| C1 | soil moisture | 2 | 3694 | 468 | 1 | 10 km block | 24 | 7.86 | 5.9 | 9.49 | 15 | 38 | 70.2 | -0.503 | 0.548 |
| P1 | peak demand | 2 | 4155 | 468 | 1 | 10 km block | 24 | 22 | 14.4 | 27.7 | 15 | 38 | 70.2 | 75.6 | 0.393 |

## Open

- Reference footprints use a 150 km radius (D13); the SMAP file used is named in Task 1.2 (smap_source)
- Authorised rate, pumping hours and verified capacity: wimas_lookup_sheet.csv (Check 2.7)
- Footprint check against 2018-2020 labels (LANID): see check_2018_2020 in the Task 1.2 summary and the irr_max_2018_2020 column of footprints.csv
- Decisions still proposed: ['D01', 'D02', 'D03', 'D04', 'D05', 'D06', 'D07', 'D08', 'D09', 'D10', 'D11', 'D12', 'D13', 'D14', 'D15', 'D16', 'D17', 'D18', 'D19', 'D20', 'D21', 'D22', 'D23', 'D24', 'D25', 'D26', 'D27', 'D28', 'D29']

## Next step (pilot Section 11)

Read the table against Section 11 of the pilot: which row of that table describes these findings decides the next step, with the supervisor.