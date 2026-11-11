# paper1: changes for plan v3 (P1.0, P1.1, P1.2)

Source: `Groundwater_Buffering_Papers_1_2_Detailed_Plan_v3_Ogallala.pdf`, cards P1.0-P1.2, sections 6 and 8.
Decided with the user on 2026-11-11:

- Depth rule: class-specific ladder (15 m herbaceous; 30 m shrub, hay/alfalfa and any block with woody cover),
  rungs 60/30/15/10 m stored as bits. The main rule stays at 30 m until a stratum (functional group x aridity
  decile) is accepted in `config.DEPTH_ACCEPTED`.
- Ring test feedback: train on the strict set, the ring test writes a proposed acceptance table, the user
  reviews it and copies it into config, then reruns `reference` and `train`.
- Gates: G1a (verify_p11) and G1b (verify_p12) halt dependent stages in `focus`/`full` mode unless listed in
  `config.GATE_OVERRIDE` (recorded in provenance); in `test` mode they are reported only.

## Slice A - block aggregation (P1.0) and Step 0

| File | Change |
|---|---|
| `step0/blocks.py` | 16-bin histograms of ET, NIRv, NDMI per block-month (edges in `HIST_EDGES`); observation support: share of Landsat 8/9 observations (`n_obs_oli / n_obs`), mean observation day (`obs_day`), OpenET models kept (`et_ensemble_mad_count`), smallest / largest model value (`et_ensemble_mad_min/_max`), clear-scene count (`count_<model>` of `OPENET_COUNT_VARS`). Missing source variables give NaN and are listed in the file attributes. |
| `step0/drivers.py` | Monthly SPI-90 per lab cell from `spi90d` of the lab drought file (value of the last record dated inside the month); a stratum label, never a model input. |
| `step0/hand.py` | Canal / ditch flowlines (`CANAL_FTYPES`, default 336) kept as their own layer; D8 flow accumulation (upstream area, km2) stored next to HAND. |
| `step0/static.py` | `canal_dist_m`; `facc_max_km2`, `facc_mean_km2`; `lanid3_dist_m` (pixels irrigated in >= `LANID_FREQ_MIN_YEARS` years, from `derived/irrFreq`). |

## Slice B - P1.1

- `cdl_groups.py`: functional groups (winter cereal, summer row crop, spring / other annual, perennial grass,
  shrub, hay, fallow) over the CDL groups; woody groups (shrub, forest).
- `p1_1/rules.py`: new bits WTD10, WTD30, CANAL (main), FIRE (main, slot passes when unset), IRR3 (sensitivity),
  CLASSDEPTH (passes the class rung). `WTD` = depth >= max(class rung, accepted rung or 30 m).
  Sets: `reference`, `reference_gde`, `reference_irr3`, `reference_ring` (all main rules except depth pass,
  class rung passes, `WTD` fails). The HAND-only crop set becomes the run-on diagnostic set.
- `p1_1/build.py`: stores functional group, previous-year group, woody pixel count, depth rule used, the sets.
- `p1_1/verify.py`: V1 per functional group x decile; V2 failure diagnosis (by distance to LANID and by
  class); V4 standardised mean differences (soil water, root-zone water, NCCPI, elevation, rain, PET; clay and
  slope listed as not available); P1.0 valid-share table (Apr-Sep 2012 vs 2014-2016); gate G1a.

## Slice C - P1.2

- `p1_2/features.py`: measurement inputs per output month (`MEAS_INPUTS`), functional group one-hot
  (primary), CDL group one-hot (secondary, `SPECIES_INPUT`), previous-year functional group one-hot;
  per item-month drought class (SPI-90 <= -0.8 / between / >= +0.8) and observation regime (`OBS_REGIMES`).
- `p1_2/baselines.py`, `train.py`: conformal strata target x functional x aridity tercile x season x drought x
  regime; strata below `CAL_MIN_N` merge upward (regime, then tercile, then drought, then functional, then
  season). Fold (0,0) variants: no pre-season, no measurement inputs, ET/PET ratio target. Class-flip test.
- `p1_2/verify.py`: per-stratum coverage; mass-bound violations; G1b drought bias on unseen droughts
  (2011-2013, 2022); steps at regime boundaries; residual vs measurement inputs; known-zero sets (held-out
  reference, ring) <= 7 %; run-on set; ring test and the proposed `DEPTH_ACCEPTED`; HAND / flow accumulation /
  canal-distance residual checks; 2008-09 crop-layer test; wheat after fallow; LANID drift after 2017; gate G1b.
- `core/stages.py`, `pipeline.py`: `gated_by` stages; G1a gates `samples`.

Not built (data not on disk): Kansas WIMAS comparison and known-zero set, MTBS fire perimeters, USGS
irrigation 2018-2020. Coupling slope (P1.4) is reported as pending in the ring test.
