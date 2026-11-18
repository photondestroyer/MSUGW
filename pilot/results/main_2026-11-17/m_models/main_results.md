# Main analysis: the two contrasts

Freeze: pilot-freeze-2026-11-17. Expected sign (D31): negative. Evidence: a 95% interval that excludes zero. Ensemble ET, thickness error corrected, primary samples:

| contrast | model | sample | value_mm | lo_mm | hi_mm | n | cells | clusters |
|---|---|---|---|---|---|---|---|---|
| soil moisture | C1 | sm | -6.08 | -12.3 | 0.132 | 50873 | 7284 | 32 |
| peak demand | P1 | peak | 0.114 | -3.39 | 3.62 | 64183 | 7583 | 37 |

Every ET model, the common sample (J1), the district base error, the uncorrected estimates and the known-irrigated subset are in main_contrasts.csv; T1 in main_estimates.csv. The range across ET models is not a confidence interval.

## Support of the evaluation points

| variable | third | eval_point | p05 | p95 | inside |
|---|---|---|---|---|---|
| deficit_peak_mm | thin | 202 | 110 | 210 | True |
| deficit_peak_mm | thick | 202 | 113 | 238 | True |
| S | thin | -1.05 | -1.3 | 0.785 | True |
| S | thick | -1.05 | -1.29 | 0.815 | True |