# Benchmark reproduction: Mieno et al. (2024)

Code: github.com/tmieno2/Drought-Production-Risk-Aquifer (their processed data in Data/data-processed).
R: R version 4.4.3 (2025-02-28 ucrt); fixest 0.14.2; mgcv 1.9.4; bootstrap draws 1000 (theirs: 1000), run serially (theirs: parallel), so bootstrap intervals agree in distribution, not draw for draw. Point estimates do not depend on the bootstrap.

## Sample

| crop | n_obs | n_counties | year_min | year_max |
|---|---|---|---|---|
| corn | 8814 | 305 | 1985 | 2018 |
| soy | 5997 | 233 | 1985 | 2018 |

## Against the paper

| statement | published | balance_used_mm | sat_from_m | sat_to_m | reproduced_loss_t_ha | reproduced_loss_bu_acre | expected_bu_acre |
|---|---|---|---|---|---|---|---|
| corn: yield loss when the aquifer thins from ~70 m to ~40 m (230 to 130 ft) at a 950 mm deficit | about 15 bu/acre | 952 | 70 | 40 | 0.949 | 15.1 | 15 |
| corn: yield loss when the aquifer thins from ~100 m to ~70 m (330 to 230 ft) at a 950 mm deficit | negligible | 952 | 100 | 70 | 0.147 | 2.3 | 0 |

## Figures reproduced (11)

- g_avg_yield.pdf
- g_avg_yield_dif.pdf
- g_ir_yield_dif.pdf
- g_map.pdf
- g_share.pdf
- g_share_dif.pdf
- g_total_impact.pdf
- g_y_wd_all.pdf
- g_yield_intensive.pdf
- g_yield_with_conf_corn.pdf
- g_yield_with_conf_soy.pdf

Differences to look at: any statement whose reproduced loss differs from the published one by more than the rounding of the paper (the paper rounds to whole bushels).