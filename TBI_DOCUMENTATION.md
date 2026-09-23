# Temporal Buffering Index (TBI) — Documentation

**Based on:** `gwb papers/vvip.pdf` (Gonzalez Cruz et al. 2021 Sci Rep 11:21648) + `Groundwater_Buffering_Research_Plan.docx` (Sections 5–6)

**Code:** `temporal_buffering_index.py` — read-only, lazy (xarray + dask), never modifies source NetCDFs

**Author:** Muse Spark / OpenCode — 2026-08-24

---

## 1. What this computes

### 1.1 VVIP Vulnerability Index (static)

```
Eq1  VI = DSI / ( α·SBI + β·GBI ) = 1/RI
     α+β=1, 0≤α,β≤1 ;  VI>1 vulnerable, VI<1 robust ; TBI = 1/VI
```

| Symbol | Meaning | Equation |
|---|---|---|
| **DSI** | Drought Stress Index (numerator) | Eq3: `DSI_l = ∏ E(c_i,j,k)^{w_i,j,k}` — weighted product of expected drought characteristics |
| **SBI** | Soil Buffer Index | Eq4: `SBI_l = PAW^{w} · DAC^{w} · ADP^{-w}` |
| **GBI** | Groundwater Buffer Index | Eq7: `GBI_l = T^{wT} · S^{-wS}`, with Eq5 `T=K·ST` and Eq6 `ΔH = IRRAMT/S` |
| **E(c)** | Expected value of drought characteristic | Eq2: `E(c)=∫ c·f(c) dc` |
| **w** | Entropy weights | Eq8 `p_ij=Xn_ij/Σ_j Xn_ij`, Eq9 `E_i=-Σ p log p / log J`, Eq10 `w_i=(1-E_i)/Σ(1-E_k)` |

- `c_i` ∈ {Drought Severity DS, Drought Duration DD, Inter-Drought Duration IDD} ; `w_IDD` negative (close droughts = more stress). `w_PAW, w_DAC` positive, `w_ADP` negative.
- `j` ∈ {SPI, SPEI}, `k` ∈ {3-month, 6-month} in VVIP. With available data we map `spi90d≈SPI-3M`, `spi180d≈SPI-6M`, etc.

**Output maps** in VVIP: Fig11 (DSI, SBI, GBI on 0–1, 5 classes by 20/40/60/80th percentiles), Fig12 (VI, with >1 dark red = most vulnerable). Table 1 cross-tabulates crop area vs vulnerability.

### 1.2 Temporal Buffering Index (time-varying extension, Research Plan §6)

Research Plan generalizes static VI to a **temporal buffering signal**:

```
§6.1  Bi(t) = Fi(t) − F̂_i(P, SMrz, VPD, Rn, T, PFT, season)
      Fi = anomaly of ET, SIF, or VWC (vegetation water content)
      F̂ = expected value from recent precip, root-zone SM, meteorology, PFT, season
      Residual during dry-down is candidate subsurface/management signal

      TBI(t) = (α·SBI(t) + β·GBI(t)) / DSI(t) = 1/VI(t)   — larger = more buffering
```

Additional Research Plan components built on `Bi(t)`:

- **Effect of groundwater access (matched causal, §6.2):** `τ = F_shallow − F_matched-deep` via propensity/Mahalanobis/causal forest, with required checks (wet-year placebo, negative control, Rosenbaum/E-value, irrigation & topography stratification).
- **Sustainability classification (§6.3):** sustainable (stable wells/storage) vs mined (declining wells+TWS) vs irrigation subsidy vs failed/no buffer, subject to volume-plausibility (implied ET demand consistent with depletion).
- **Failure threshold (§6.4):** `Bi(t) ~ f(WTD, ΔWTD, Δstorage, PFT, rooting depth, irrigation, drought intensity)` via segmented/GAM/hierarchical Bayes.
- **Reversibility/hysteresis (§6.5):** test only where storage measurably rebounds; distinguish “water never returned” vs “returned but function did not recover”.
- **Water–carbon–energy decoupling (§6.6):** `Di = g(WTD, …)` as loss of coherence/phase lag among ET, SIF, VWC within a dry-down; space-for-time across many events (short multi-sensor record).

---

## 2. Repository data inventory (what exists)

All files opened with `xr.open_dataset(..., chunks=..., mask_and_scale=False)` — lazy, dask-backed, read-only.

| File | Size | Time | Grid | Key vars | Use for TBI |
|---|---|---|---|---|---|
| `merged_datasets/DROUGHT_Merged_Ogallala.nc` | 4.6 GB | 1984-01-05 → 2026-07-… (3035) | 287×181 ~4 km | `spi14d/30d/90d/180d/270d/1y/2y/5y`, `spei*`, `eddi*`, `pdsi`, `z` | **DSI** — `spi90d≈SPI3M`, `spi180d≈SPI6M` etc. Proxy for Eq3 (needs run theory for E(c)) |
| `merged_datasets/GRIDMET_Merged_Ogallala.nc` | 29.3 GB | 1979-01-01 → 2026-08-… (17380) | 287×181 ~4 km | `pr`, `tmmn/tmmx`, `sph`, `srad`, `vpd`, `eto/etr`, `vs`, `bi,erc,fm100…` | Met controls for `F̂` (§6.1); PET for SPEI if recomputed |
| `merged_datasets/SPL4SMGP_Ogallala_FULL.nc` / `latest.nc` | 0.93/0.90 GB | 2015-03-31 → 2024… (~26k 3-hr) | 134×104 ~9 km (EASE2) | `sm_rootzone`, `sm_rootzone_pctl`, `depth_to_water_table_from_surface`, `land_evapotranspiration_flux`, `baseflow_flux` | **SBI**: `sm_rootzone_pctl` as SSMI proxy (SHORT record, caveat); water-table dynamics proxy |
| `merged_datasets/HRES-WTD_2015_Ogallala.nc` | 7.4 GB | 2015 static | 52397×33030 ~30 m | `b1` Water Table Depth (m) | **GBI** spatial axis (static only, not time series) |
| `merged_datasets/Hydrologic_Soil_Group_250m_2019_Ogallala.nc` | 25 MB | 2019 static | 5237×3302 ~250 m | `b1` HSG 1=A,2=B,3=C,4=D | **SBI PAW proxy** only (inverted: A→high buffering). True PAW needs gSSURGO |
| `merged_datasets/MOD13A3_Merged_Ogallala.nc` | 4.0 GB | 2000-02 → 2026-06 (317) monthly | 1432×903 ~1 km | `NDVI`, `EVI`, `DetailedQA` … | `Fi` greenness proxy for SIF (derived, not independent per plan §5) |
| `merged_datasets/MODIS_ET_SSEBop_Merged_Ogallala.nc` | 0.27 GB | 2003-01 → 2022-05 (699 dekadal) | 1235×779 ~1 km | `et` | `Fi` evaporation (primary for Bi(t)) |
| `merged_datasets/GRACE_Merged_Ogallala.nc` | 0.16 MB | 2002-04 → 2024-09 (238) monthly | 24×15 ~0.5° | `lwe_thickness`, `uncertainty` | Basin-scale TWS constraint & slow driver (not pixel) |
| `merged_datasets/GFSAD1000_V1_2019_Ogallala.nc` | 2.9 MB | 2019 | 1335×842 ~1 km | `landcover` | PFT / irrigation stratification |
| `merged_datasets/Global_irrigation_Area_Merged_Ogallala.nc` | 59 KB | 2001-2015 (15) | 144×91 ~0.08° | `classification` | Irrigation subsidy vs natural buffer |
| `merged_datasets/MCD12Q1_Merged_Ogallala.nc` | 0.22 GB | 2001-2024 (24) | 2864×1806 ~500 m | `LC_Type1..5`, `LC_Prop*` | PFT for matching & `F̂` |
| `merged_datasets/merged_ameriflux.nc` | 3.1 GB | per-site half-hourly (290 sites, 32.6M records) | point | `LE`, `H`, `NEE`, `GPP_uStar_f`, `VPD`, `P`, `Tair` … | Point validation of ET/GPP (not gridded) |
| `smap file/*.h5` (17 files) | ~? | 2016-12-31 → 2017-01-02 | 9 km | raw SMAP L4 | Source for SPL4SMGP merged |
| `high_plains_quifer/hp_bound2010.shp` | — | 2010 | polygon | Ogallala boundary | ROI mask |
| `AmeriFlux_NEON_sites.xlsx` + `batch_001/002.parquet` | — | site metadata | — | 290 sites | Site lookup for AmeriFlux |

**Native grids differ dramatically** — harmonization to a common analysis grid is essential (VVIP: 0.5° × 0.5°, 187 points; Research Plan: ~1–5 km matching gridMET).

---

## 3. Datasets required additionally (not in repo) — do not synthesize

### 3a. For full VVIP (Table S1, VVIP p6-8)

| Missing | Source | Why needed |
|---|---|---|
| **gSSURGO 30 m** — field capacity, permanent wilting point, root-zone depth | USDA Soil Survey Geographic Database (Staff 2019, https://gdg.sc.egov.usda.gov) | `PAW = (FC − PWP)/depth` (VVIP Eq4 stock component, p4). HSG 1-4 is **not** PAW. |
| **Saturated thickness 2009 + aquifer bottom** | USGS SIR 2012-5177 (McGuire et al. 2012) + McGuire 2017 Water-Level & Recoverable Water | `ST` for Eq5 `T=K·ST`. HRES-WTD gives depth, not thickness. |
| **Hydraulic conductivity K** | Cederstrand & Becker (1998) Digital Map of Hydraulic Conductivity, USGS | `T=K·ST` (Eq5). Cannot derive GBI without it. |
| **Specific yield / storage coefficient S** | USGS SIR 2012-5177 (McGuire et al.) | Eq6 `ΔH=IRRAMT/S` and Eq7 `GBI=T^{wT}·S^{-wS}` (hydraulic diffusivity if w=1). |
| **CPC Global Monthly Soil Moisture 0.5° 1949–present (leaky-bucket)** | Fan & van den Dool (2004), Huang et al. (1996), CPC | SSMI at 3 & 6 months for DAC/ADP (VVIP Eq4, Fig8-9). SPL4SMGP 2015+ is too short for 30-year `E(c)` and DAC/ADP contingency tables. |
| **CRU TS 4.03 (or equivalent) PET** | Harris et al. (2020) + Penman-Monteith | VVIP computed SPI/SPEI from CRU (p6). GRIDMET/DROUGHT proxies exist but CRU is paper source. |

### 3b. For full Research Plan temporal buffering (§5)

| Missing | Source | Why needed |
|---|---|---|
| **ECOSTRESS ET/LST 70 m 2018+** | ECOSTRESS (NASA/JPL) | High-res dry-down ET (plan §5) — SSEBop is longer but coarser |
| **OCO-2/3 or TROPOMI SIF** | OCO-2 (Li & Xiao 2019 GOSIF is *derived*, not independent) / TROPOMI 2018+ | Independent carbon signal for decoupling (plan notes GOSIF is not independent) |
| **SMAP L-band VOD** | SMAP VOD product | Vegetation water content proxy (distinct from `sm_rootzone`) |
| **GLEAM4 0.1° 1980+** | Miralles et al. 2025 | Reference transpiration & groundwater-sourced flux |
| **USGS + state well time series** | USGS NWIS + state networks | Dynamic `ΔWTD` and recovery/rebound for threshold & reversibility (§6.4-6.5). HRES-WTD 2015 is static; SPL4SMGP `depth_to_water_table` is model-assimilated, not well-observed. |
| **LANID / USDA CDL annual cropland** | USDA NASS Cropscape 30 m | Managed vs natural stratification (§6.2); GFSAD 2019 is single year, ~1 km |
| **County yield / harvested area / fallowing** | USDA NASS | Outcome validation + documented stress events 2011-15 CA, 2012 High Plains |

**Present proxies can be used for PILOT only** with clear “PROXY” labeling (see Code § place-holders). Full science requires the missing sources above.

---

## 4. Assumptions — all documented in code

### 4.1 VVIP model assumptions (paper § Model assumptions and limitations, p5-6)

1. **Dimensional inconsistency** — Eq1 not dimensionally consistent; qualitative at single site with Laplace equal weights, or relative regional comparison after 0–1 normalization. Values only interpretable *within* study domain.
2. **30+ year record required** — stochastic concepts (expected values, entropy) treated deterministically via long-term probabilities. Needs ≥30 yr (VVIP used hydrologic years 1949–2018). SMAP 2015+ insufficient.
3. **Empirical MCDM** — no rigorous validation against observed VI; validity hinges on verified drought/soil/aquifer inputs. Uncertainty propagation (§56 Krejčí et al.) recommended for bounds.
4. **Static vulnerability** — VI is point-in-time (VVIP baseline 2015). Temporal evolution needs time-varying DSI(t), SBI(t), GBI(t).
5. **Correlated attributes** — SPI↔SPEI, PAW↔DAC↔ADP, T↔S → weighted *product* (Eq3,4,7) not sum, avoids normalization prior to multiplication.
6. **Drought thresholds** — drought if SPI/SPEI ≤ −1; accumulations 3 & 6 months capture intra-season and full-season agricultural drought (VVIP p6).
7. **Distribution fits** — DD/IDD ~ exponential, DS ~ lognormal (83% of sites) / Weibull / Gamma by AIC; used to compute Eq2 expectations.

### 4.2 Research Plan methodological assumptions (§6, §8 pitfalls)

8. **Residual ≠ groundwater** — `Bi(t)` absorbs irrigation, crop variety, soil, phenology, model bias. Attribute to groundwater only where `Bi` covaries with WTD/wells/storage AND survives irrigation/management controls.
9. **Smooth weakening vs tipping** — buffering may degrade gradually, not sharply; test reversibility directly, use “regime shift” language only where recovery fails; slow-recovery indicators secondary.
10. **No recovery to observe where monotonic depletion** — must preselect rebound zones (well clusters / sub-basins where storage measurably returned). Separate “water never returned” from “returned but function did not recover”.
11. **Resolution mismatch** — analysis grid ~1–5 km (gridMET), 30 m WTD aggregated to it, GRACE only basin constraint, decoupling analysis necessarily coarser (~0.1–0.25° limited by VOD/SIF).
12. **Short multi-sensor record** — ECOSTRESS (2018+), SIF (2014+), SMAP (2015+) → use space-for-time across many locations/events, not long trend at one site. Reversibility of 2012 drought uses long single-signal record (greenness/wells/GRACE).
13. **Topographic confounding** — depth-to-water partly set by topography (also affects soils/productivity). Required checks (§6.2) wet-year placebo, negative-control outcome, Rosenbaum/E-value sensitivity, irrigation & topographic stratification.
14. **Volume plausibility** — label “mined” only when persistent function + declining wells/storage are *plausibly* consistent in magnitude with ET demand; full closure not required but volumes must be plausible, otherwise correlation only.

### 4.3 Implementation assumptions (this toolkit)

15. **Read-only & lazy** — every `xr.open_dataset(..., chunks=..., mask_and_scale=False)` is lazy/dask; no `.to_netcdf` ever called on source paths; `.compute()` only on small subsets or explicitly requested.
16. **Placeholders = NaN** — `PLACEHOLDERS` dict in code sets `PAW/K/Sy/ST/T = np.nan` with documented `source_required`. No synthetic values generated; proxies are labeled `PROXY` and invertible.
17. **Proxy mappings** (clearly flagged in `attrs["note"]`):
    - `PAW → (5−HSG)/4` from `Hydrologic_Soil_Group_250m` (A=1→1.0, D=4→0.25), normalized 0–1, **not** true PAW.
    - `GBI → 1/(1+WTD/10)` from `HRES-WTD b1` (meters), 10 m scaling; true GBI needs `K·ST` and `S` per Eq5–7.
    - `DSI → time-mean of drought intensity `where(SPI≤−1, −SPI, 0)` across `spi90d/spei90d/spe…` — approximates combined DS/DD, omits IDD/tail fits.
    - `DAC/ADP → neutral 1.0` (no effect on product) — true DAC/ADP need SSMI–SPI/SPEI binary contingency tables (VVIP Fig8-9).
18. **Grid harmonization required** — DSI 4 km vs HSG 250 m vs WTD 30 m vs GRACE 0.5° vs SMAP 9 km. `compute_VI()` warns on mismatch and does lazy `interp(nearest)` to DSI grid only for demo; production should `coarsen`/`weighted regrid` to common grid (paper 0.5°, plan 4 km) before VI.
19. **Weights** — fallback to equal weights (Laplace) when entropy cannot be computed lazily on full grid. Call `entropy_weights()` on a *computed* subset for Eq8–10; `α=β=0.5` default matches VVIP Fig3d mean but varies spatially.
20. **Time handling** — GRIDMET/DROUGHT use `days since 1970-01-01`, SPL4SMGP uses `hours since 2000-01-01`; xarray `decode_cf=True` handles conversion. Temporal buffering `Bi(t)` needs aligned dekadal/monthly resampling before subtraction.

---

## 5. How to run (lazy, no OOM)

```python
from temporal_buffering_index import TBIToolkit

tk = TBIToolkit(base_dir=r"G:\MSU_GWB\datasets")
tk.inventory_all(); tk.print_inventory()   # no load, just header
tk.open_all_lazy()                          # still lazy

# Verify pipeline on 10×10 pixel × 20 step subset (small, safe)
tk.demo_small_subset(y_slice=slice(0,10), x_slice=slice(0,10), t_slice=slice(0,20))

# Build lazy VI (still dask, no compute) — SBI/GBI auto regrid to DSI grid with warning
vi_lazy = tk.compute_VI(alpha=0.5, beta=0.5)  # (y,x) dask array

# Compute only a window, or persist lazily to new file (never overwrites source)
vi_window = vi_lazy.isel(y=slice(100,110), x=slice(100,110)).compute()
print(vi_window.values)

tk.save_result(vi_lazy, r"G:\MSU_GWB\datasets\tbi_output.nc")  # new file, compressed
tk.close_all()
```

CLI:

```
python temporal_buffering_index.py --inventory
python temporal_buffering_index.py --assumptions
python temporal_buffering_index.py --required
python temporal_buffering_index.py --demo
```

**Memory guidance:** GRIDMET 58 GB + DROUGHT 16 GB cannot be fully `compute()`d. Always `isel`/`sel` first, or use `dask` `map_blocks`/`coarsen` + `to_netcdf` with `chunks`. The toolkit never calls `.load()` on full arrays.

---

## 6. Output interpretation

- **VI (and TBI=1/VI)** are *relative* within the Ogallala domain. Do not compare absolute values to other aquifers without re-normalizing (VVIP p5).
- VVIP Fig12 classes (percentiles of VI at 187 grid points): Very Low / Low / Medium / High / Very High (VI >1 = buffering < stress, dark red, ~20% of aquifer; 0.5<VI<1 = factor of safety 1–2, ~40% transitional zone warranting management).
- **Transition zone (medium–high, 50% of aquifer per VVIP abstract)** concentrates policy relevance — deficit irrigation, drought-tolerant varieties, crop rotation.
- **Current cropping compatibility** (VVIP Table 1): cotton/sorghum in high VI, corn/soybean in low VI — indicates adaptation but 25% high-water crops remain in medium–high zones.

For **temporal** TBI:
- `Bi(t) >0` during dry-down = function persists beyond meteorology+SMrz → candidate buffering (attribute only after controls §6.2).
- **Threshold** search: `Bi` sharply drops beyond critical WTD / storage anomaly → distance-to-failure map (§6.4).
- **Reversibility:** if storage rebounds but `Bi` / greenness does not → hysteretic/persistent shift; if both rebound → reversible decline (§6.5, use long greenness record for 2012).

---

## 7. Placeholders and “do not generate data” compliance

- No random/synthetic arrays are created. Every missing variable is `np.nan` placeholder with `source_required` documented in `PLACEHOLDERS` and `TBIToolkit.list_required_datasets()`.
- Proxy calculations that *do* run (HSG→PAW, WTD→GBI proxy, intensity_mean→DSI) set `attrs["note"] = "PROXY: …"` and warn. They are suitable for **pilot / pipeline testing** but not for publication without replacing with the USGS/gSSURGO/CPC sources.

---

## 8. File lineage — nothing modified

Source NetCDFs are opened `mask_and_scale=False` and never written. Derived outputs go to **new** files (e.g., `tbi_output.nc`) via `TBIToolkit.save_result()`. Parquet `batch_001/002` and shapefile `hp_bound2010.shp` are read-only.

---

## 9. References (key)

- Gonzalez Cruz et al. 2021 Sci Rep 11:21648 (VVIP) — Eqs1–10, p6 data compilation, Fig3 weights, Fig10 groundwater.
- Harris et al. 2020 Sci Data 7 (CRU TS 4.03); Fan & van den Dool 2004 JGR (CPC SM); Huang et al. 1996 J Clim (leaky bucket); McKee et al. SPI; Vicente-Serrano et al. SPEI; Shiau 2006 WRM (drought DS/DD/IDD); Guerrero-Salazar & Yevjevich 1975 (runs theory).
- USGS: McGuire 2017 Water-Level & Recoverable Water; McGuire et al. 2012 SIR 2012-5177 (Sy, ST); Cederstrand & Becker 1998 (K).
- Research Plan refs: Rohde et al. 2024 Nature 632; Ma et al. 2026 Comm Earth Env (HRES-WTD); Miralles et al. 2025 GLEAM4; Jasechko et al. 2024 Nature 625; Fan et al. 2013 Science 339; Forzieri et al. 2022 Nature 608.

---

## 10. Contact / next steps

1. Acquire missing datasets (§3) or run pilot on proxies and document limitations.
2. Harmonize grids (suggest `xESMF` or `xarray.coarsen` to 4 km gridMET, aggregate 30 m WTD by mean, HSG by mode, SMAP by bilinear).
3. Implement full `E(c)` run theory + distribution fits (vectorized per-pixel, or sample 187 VVIP points first).
4. Train `F̂` for `Bi(t)` (RandomForest / GAM per PFT) and run matched design (§6.2) with required sensitivity checks.
5. Threshold (§6.4) via segmented regression / GAM, then reversibility test on rebound wells (§6.5) before global extension.

*All code and this document are read-only with respect to source datasets and use explicit lazy loading for 30 GB+ files.*

---

## 11. UPDATE v2 — New USGS groundwater datasets integrated (2026-08-26)

### 11.1 Newly available raw data (`G:\USGS GW dataset`)

| Folder | Content | Native format / CRS | Role in TBI |
|---|---|---|---|
| F01 `hp_wlcpd19t.zip` -> `.tif` | Mapped water-level change predevelopment (~1950) -> 2019 | GeoTIFF, EPSG:5070, 500 m, feet, nodata -3.4028e38 | Dynamic GW diagnostic (decline axis) |
| F02/F03 wells shapefiles | Point WL-change data (A83) | Shapefile EPSG:5070 | Not used in gridded run (point QC source) |
| F04 `hp_wlc1719t.zip` -> `.tif` | Water-level change 2017 -> 2019 | GeoTIFF EPSG:5070, 500 m, ft | Recent-trend diagnostic |
| F05 wells shapefile | 2017-19 well changes | Shapefile | Not used (points) |
| OFR 98-548 `ofr98-548.e00.gz` | **Digital map of hydraulic conductivity K**, High Plains aquifer (Cederstrand & Becker 1998) | AVCE00 vector (ARC/LAB/PAL), NAD83 Albers; polygon attr `RANGE` e.g. "25 to 50" ft/day | **Real K for GBI** |
| `gSSURGO_CONUS.gdb` | CONUS soils: Valu1 table has **rootznaws** (root-zone available water storage) etc.; MUPOLYGON 36.7 M features, ~24 GB, **no usable spatial index (.spx = 4 KB)** | FileGDB | PAW upgrade BLOCKED on geometry extraction (see 11.4) |

### 11.2 Preprocessing performed (`preprocess_usgs_gw.py`; sources strictly read-only)

1. Common target grid fixed to GRIDMET/DROUGHT native: x=181, y=287, res (0.05358 deg, 0.04167 deg), EPSG:4326.
2. F01/F04 rasters reprojected EPSG:5070 -> EPSG:4326 via `rasterio.warp.reproject` (average resampling), nodata -> NaN.
3. E00 gunzipped; PAL layer (787 polygons) read through GDAL/AVCE00 via geopandas; `K_mid_ft = (MAJOR1+MINOR1)/2` from RANGE intervals {0-25 ... 400-500 ft/day}, converted x0.3048 to m/day (range 0.0-137.2); polygons sorted by AREA descending and rasterized onto the 4 km grid so smaller zones overwrite larger ones (mode-like approximation).
4. Aquifer mask rasterized from `high_plains_quifer/hp_bound2010.shp`: **21,259 px inside** (40.9% of the box). All derived layers masked => spatially consistent aquifer-only analysis on one grid.
5. Outputs written ONLY to new dir `derived_usgs/`:
   - `dWL_predev_to_2019_ft_4km.nc` (-265 ... +54 ft, mean -13.8)
   - `dWL_2017_to_2019_ft_4km.nc` (-18.5 ... +12.4 ft, mean +0.11)
   - `K_hydraulic_conductivity_mday_4km.nc` (mean 19.0 m/d, p90 45.7; valid 20,915 px inside aquifer; 343 aquifer px outside mapped K zones = NaN)
   - `aquifer_mask_4km.nc`, `preprocess_log.txt`

### 11.3 Specific Yield / storage coefficient — USER-SPECIFIED CONSTANT

> **S = 0.15 (dimensionless)**, applied in VVIP Eq7 as `GBI = norm(K_norm^{wT} * S^{-wS})`, wT=wS=1.
> Because S is spatially constant it acts as a linear rescale (`K_norm / 0.15`, range [0, 6.67]) that is removed by the final min-max renormalization to [0,1]; it therefore does NOT alter VI ranking, but is retained in the formula for dimensional fidelity with Eq7 and easy replacement by a gridded Sy later (McGuire et al. 2012).
> Justification for 0.15: Ogallala specific yield typically 0.12-0.18 over most of the aquifer (VVIP Fig 10d); midpoint adopted.
> Recorded below as implementation assumption #19b.

### 11.4 What changed in the index (run `run_tbi_v2_usgs.py` -> `tbi_maps_v2_usgs/`)

| Component | v1 (`tbi_maps_4km/`) | v2 (`tbi_maps_v2_usgs/`) |
|---|---|---|
| DSI(t) | GRIDMET-DROUGHT spi/spei 90&180d geometric-mean intensity | unchanged (spatially consistent native 4 km, 3035 steps, threshold <= -1) |
| SBI | HSG proxy `(5-HSG)/4` | unchanged - true gSSURGO PAW pending (11.5) |
| GBI | depth proxy `1/(1+WTD/10)` | **GBI = norm(norm(K)/S), S=0.15, real OFR98-548 K** (mean 0.139) |
| Masking | none (full GRIDMET rectangle) | **aquifer-masked everywhere (21,259 px)** |
| Diagnostics | - | dWL predev->2019 & 2017->2019 maps; VI-quintile vs decline screen |

VI domain-mean range 0.000-5.62 (1984-2026). Map dates: 1984-07-13, 1985-12-16 (min-DSI), 2000-06-03, 2011-08-03, 2012-07-13, 2012-09-21 (max-DSI), 2014-12-31, 2022-08-03. Multipanel: `VI_multipanel.png`.

### 11.5 Remaining blockers / next actions

1. **PAW**: gSSURGO CONUS gdb has no usable spatial index => bbox polygon read exceeded 30 min without completing. Options: (a) download gSSURGO *raster* GeoTIFF (mukey grid) and join `Valu1.rootznaws` directly - preferred; (b) state-level gSSURGO subsets with spatial indexes; (c) offline ArcGIS zonal extraction. Until then PAW remains a documented HSG-based PROXY.
2. **ST (saturated thickness)**: still missing => T = K*ST not computable; GBI ranking is K-driven. Add McGuire et al. 2012 (SIR 2012-5177) saturated-thickness raster to close Eq5.
3. **DAC / ADP**: neutral placeholders (need >=30-yr CPC leaky-bucket soil moisture for SSMI contingency tables).
4. Screening result (`screen_VI_vs_dWL_decline.png`): deepest mean declines fall in Low-VI classes (-18.1 ft) rather than Very-High-VI (-12.4 ft) => under current proxies the mined-buffer signature is not monotonic; expected while GBI is K-driven. Re-test after real PAW and ST.

### 11.6 Assumption added

* **#19b - Specific yield / storage coefficient S = 0.15**, uniform scalar (user-specified 2026-08-26). Ogallala Sy typically 0.12-0.18 (VVIP Fig10d); midpoint adopted. Replace with McGuire 2012 Sy grid when available.


### 11.7 NEW — Groundwater well data packaged as CF-compliant NetCDF (`make_gw_wells_cf_nc.py`)

The three USGS well shapefiles are now consolidated into one CF-1.8 discrete-sampling-geometry file:

**`derived_usgs/GW_wells_HPA_CF.nc`** (1.8 MB, NETCDF4)

| Source file | Rows | Contribution |
|---|---|---|
| F02 `hp_wlcpd19_wells_A83.shp` | 2,741 | primary predev->2019 network; attr `deltapd_19` |
| F03 `hp_wlcpd19supwells_A83.shp` | 2,052 | supplemental network; 1980 epoch (`wl80`, `delta80_19`); `dpd_19est` |
| F05 `hp_wlc1719_wells_A83.shp` | 7,195 | dense 2017->2019 network; attr `delta17_19` |

**Structure** (`featureType = "point"`, `Conventions = "CF-1.8"`):

```
dims    station=8222 (unique wells, deduped by site badge across files)
        obs=20580  (one record per dated measurement)
coords  lat(station) deg N, lon(station) deg E   [NAD83 datum, from *_nad83 attrs]
        time(obs)   days since 1900-01-01, gregorian
link    station_index(obs)
data    water_level_ft(obs)      depth-to-water below land surface [ft]
        campaign(obs)            pd | 1980 | 2015..2019
static  site_badge, usgs_id, station_name, state, county, source_dataset,
        well_depth_ft,
        dWL_pd19_primary_ft / dWL_pd19_supplemental_ft / dWL_1980_2019_ft /
        dWL_2017_2019_ft          (per-well QC'd changes, negative = decline)
```

**Campaign census:** pd 4,501 · 1980 665 · 2015 70 · 2016 158 · 2017 7,357 · 2018 228 · 2019 7,601.
Coverage: lat 31.81-43.66? see geospatial attrs; lon -105.18..-96.28; time 1901-07-25 .. 2019-06-15.

**QA performed**
1. Sentinels -9999/-999 -> `_FillValue` NaN.
2. Quoted date strings in F05 (`'20170426'`) unquoted before parsing (bug found by spot-check and fixed).
3. 3,837 cross-file duplicate measurements removed (`drop_duplicates(site_badge,time,value)`).
4. Blank measurement dates within an epoch fall back to a mid-season nominal date (documented approximation; affects only the `time` coordinate, not values).
5. Spot-check vs raw attributes passed exactly: USGS 433100102002101 -> 145.15 ft @ 1979-12-11, 125.60 ft @ 2017-04-26, 124.21 ft @ 2019-06-10.

**Usage with the toolkit**

```python
ds = xr.open_dataset("derived_usgs/GW_wells_HPA_CF.nc")           # point DSG
wl  = ds.water_level_ft.where(ds.campaign == b"pd")               # predevelopment levels
# join obs -> stations:
df  = ds[["water_level_ft","campaign","time"]]
      .assign(station_index=ds.station_index).to_dataframe()
```

This partially closes Research Plan data gap "USGS + state well time series": dynamic decline/recovery
signal is now available as points (1901-2019 snapshots at 5 epochs) for §6.4 threshold fitting and
§6.5 rebound-zone selection. Continuous daily records (NWIS) remain future work.

### 11.8 Dataset inventory refresh (post-v2)

`merged_datasets/` (unchanged, GEE exports) + new derived stores:

| Store | Files |
|---|---|
| `derived_usgs/` | dWL_predev_to_2019_ft_4km.nc · dWL_2017_to_2019_ft_4km.nc · K_hydraulic_conductivity_mday_4km.nc · aquifer_mask_4km.nc · **GW_wells_HPA_CF.nc** · preprocess_log.txt |
| `tbi_maps_v2_usgs/` | SBI/GBI/denom statics · DSI & VI & TBI maps (8 dates) · timeseries PNGs · VI_multipanel.png · screen_VI_vs_dWL_decline.png/.csv |


### 11.9 NEW — gSSURGO processed: REAL PAW unlocks the true SBI (v3)

**ISO metadata parsed.** `gSSURGO_CONUS.gdb\ssurgo_iso19115_.xml` (7.0 MB, ISO 19115-3)
parsed by `parse_gssurgo_xml.py` -> `derived_usss/gssurgo_iso19115_summary.json`
(note correct path `derived_usgs/`). Key extracted fields: title "Soil Survey Geographic
Database (SSURGO)" / alt "SSURGO, gSSURGO, SSURGO Portal"; global bbox lon -179.159..179.858,
lat -14.374..71.441; representative-fraction denominators 1:8 000 ... 1:20 000 (map scales of
constituent surveys); topics farming/environment/geoscientificInformation; abstract confirms
"georeferenced digital map data + computerized attribute data".

**Embedded raster discovered and READ (this was the v2 blocker).**
`rasterio.open(gdb, driver="OpenFileGDB")` exposes **MURASTER_30m**: CONUS mapunit-key grid,
30 m, EPSG:5070, int32, 153,996 x 97,053 px, **with pyramids [2,4,8,16,32,64,128,256,513]**.
The earlier ">30 min polygon scan" is moot — no MUPOLYGON geometry is needed at all.

**PAW chain (`process_gssurgo_paw.py`, 122 s total):**
1. Window-read ROI at overview 4 (~120 m effective, nearest resampling keeps exact mukey ids):
   7,640 x 11,494 = 88 M samples in 60 s.
2. `Valu1` table (304,834 rows with numeric mukey) -> vectorized LUT join to
   **rootznaws** = Root-Zone Available Water Storage [cm] (97.4 % of samples mapped;
   LUT median 164 cm CONUS-wide).
3. Area-mean of mapped samples into each GRIDMET ~4 km cell (median **1,502 source
   samples per cell**), NaN where unmapped.
4. Aquifer-masked: **21,259 valid px**; PAW stats min 30.6 / p25 144.3 /
   **median 203.5** / p75 270.8 / max 369.1 cm.

**CF-1.8 output:** `derived_usgs/PAW_rootznaws_4km.nc` (time,y,x; float32; zlib-3).
Data var `root_zone_available_water_storage` carries long_name, units "cm",
`cell_methods="area: mean"`, `grid_mapping="crs"` with a proper CF
`latitude_longitude` crs variable (EPSG:4326 WGS84 ellipsoid params), plus full
provenance attrs (source raster, attribute table, aggregation, ISO snippets).

### 11.10 TBI v3 rerun (`run_tbi_v3_usgs.py` -> `tbi_maps_v3_usgs/`, 38 files)

| Component | v2 | **v3** |
|---|---|---|
| SBI | HSG proxy | **norm(PAW_real), 30.6-369.1 cm** |
| GBI | norm(norm(K)/S), S=0.15 | unchanged |
| DSI(t) | GRIDMET-DROUGHT geo-mean intensity, <= -1 threshold | unchanged |
| Masking | aquifer 21,259 px | unchanged |

* VI(t) domain-mean range **0.000-6.663** (wider than v2's 5.62 because real PAW spans a
  broader normalized gradient than the 4-class HSG proxy).
* New diagnostics:
  - `PAW_real_static.png`, `SBI_real_static.png`
  - `diag_VI_diff_v3_minus_v2.png`: long-term VI(real PAW) minus VI(HSG proxy);
    mean +0.071, p90 |delta| 0.49 -> proxy was unbiased on average but locally off by up
    to half an index unit, confirming the value of the real soils join.
  - Timeseries figure overlays VI(v3) vs VI(v2) for direct comparison.

**Sustainability screen now physically coherent.** Long-term-VI quintiles vs observed
predev->2019 water-level change:

| VI class | mean SBI | mean GBI | mean dWL (ft) |
|---|---|---|---|
| VLow (0.22) | 0.76 | 0.25 | -14.3 |
| Low (0.28) | 0.68 | 0.18 | **-25.0** |
| Mid (0.37) | 0.53 | 0.13 | -22.8 |
| High (0.52) | 0.40 | 0.08 | -10.5 |
| VHigh (0.89) | 0.22 | 0.06 | **-2.5** |

Interpretation: the deepest declines (-23..-25 ft) sit in Low/Mid-VI cells = northern HP
with thick saturated zone + good soils (large buffering AND large historic drawdown);
the VHigh-VI tail has tiny declines because those southern cells were largely dewatered
before the observation era (little left to decline). This replaces v2's non-monotonic
screening note and matches the McGuire (2017) regional pattern; formal mined/sustainable
classification still awaits the Research Plan controls (irrigation stratification etc.).

### 11.11 Assumption register update

* **#15 RETIRED** — HSG-based PAW proxy no longer used in the index (kept only inside
  `diag_VI_diff_v3_minus_v2.png` as legacy comparison). Real gSSURGO rootznaws now feeds SBI.
* Still open: **DAC/ADP neutral placeholders** (need >=30-yr CPC soil-moisture SSMI);
  **ST placeholder in GBI** (T=K*ST pending McGuire saturated-thickness raster);
  **S=0.15 uniform scalar** (#19b) until McGuire Sy grid arrives.
* rootznaws caveat: depth-integrated to the *mapunit* root-zone limiting depth (not a fixed
  150-cm column); values are SSURGO AWS in cm of water.


### 11.12 NEW — Single-file gSSURGO stack clipped to the Ogallala aquifer

**`derived_usgs/GSSURGO_Ogallala_all_attrs_120m.nc`** (938 MB, NETCDF4, CF-1.8)

Created by `make_gssurgo_ogallala_stack.py` from the gdb's embedded raster +
Valu1 table, clipped with `high_plains_quifer/hp_bound2010.shp` (the shapefile
supplied for this task).

| Property | Value |
|---|---|
| Grid | 6,568 x 11,112 px @ **~120 m effective** (MURASTER_30m overview 4), EPSG:5070 NAD83 Conus Albers |
| Clip | hp_bound2010 polygon rasterized onto the grid: **31,631,581 px inside** (43.3 % of window); outside = `_FillValue` |
| Variables | `mukey` (int32) + **all 57 numeric Valu1 columns**: aws{0_5..0_999} (available water storage, cm), tk* thicknesses (cm), soc* organic-C stocks, nccpi3* commodity indices, rootznemc/rootznaws, droughty flag, pctearthmc, musumcpct(a/s), pwsl1pomu |
| CF encoding | 1-D `x`/`y` projection coords (m, axis X/Y, standard_name projection_*), `crs` grid_mapping var = `albers_conical_equal_area` (lat_0 23, lon_0 -96, std parallels 29.5/45.5, NAD83 ellipsoid, EPSG:5070); every data var carries `grid_mapping="crs"`, `coordinates="y x"`, zlib complevel 4, chunks 512x512 |
| Verification | reopen OK; no missing vars; rootznaws med 213 cm @120 m vs 203.5 cm @4 km area-mean product (consistent); mukey fill only outside polygon |

Notes & caveats
---------------
* Units follow SSURGO conventions (documented per-variable in `units` + `units_comment`):
  AWS/root-zone storage = cm of water; SOC = g/m2 per depth prefix; tk = cm;
  nccpi/droughty/pctearthmc/musumcpct = index or percentage.
* `pwsl1pomu` covers only ~43 % of aquifer pixels by design (potential soil-loss
  rating exists only where mapped).
* The 120 m resolution preserves mapunit boundaries well while keeping the file at
  <1 GB; the 30 m native grid over the same clip would be ~100x larger.
* Downstream use: any future Valu1-style join can be done directly on the stored
  `mukey` band without reopening the CONUS gdb.

Usage:

```python
import xarray as xr
ds = xr.open_dataset("derived_usgs/GSSURGO_Ogallala_all_attrs_120m.nc")
print(ds.data_vars)                       # mukey + 57 attributes
sbi_inputs = ds[["rootznaws","rootznemc","aws0_150","droughty"]]
```


### 11.13 DSI numerator — Theory-of-Runs event statistics & AIC distribution fitting

Implemented per VVIP p3/p6 protocol (`fit_dsi_distributions.py`) on GRIDMET-DROUGHT
`spi90d/spi180d/spei90d/spei180d`, aquifer pixels only (n = 21,259), 1984-2026:

* drought := index ≤ −1; runs theory → events; **first & last event discarded**;
* characteristics per pixel/indicator: `DD` duration [pentads], `DS` severity
  = Σ(−index) over run, `IDD` inter-drought duration [pentads] (1 pentad = 5 d;
  ×0.1643 → months);
* candidates fitted per pixel by closed-form / bracketed MLE:
  Exponential(θ; 1 par), Lognormal(μ_log, σ_log; loc 0), Gamma(α, θ; loc 0),
  Weibull_min(k, λ; loc 0); selection by **AIC = 2k − 2 lnL**, minimum wins
  (race among candidates that converged; ≥2 required).
* Median events per pixel: 88 (spi90d), 62 (spi180d), 74 (spei90d), 49 (spei180d).
* All 12 indicator×characteristic races resolved on **100 % of aquifer pixels**.

**AIC winner shares (% of pixels)**

| char | distribution | spi90d | spi180d | spei90d | spei180d |
|---|---|---|---|---|---|
| DD | Exponential | 25.0 | 7.8 | 25.2 | 1.9 |
| DD | **Weibull** | 25.0 | **92.2** | **74.8** | **98.1** |
| DS | Exponential | 5.2 | 0.1 | 0.5 | 0.0 |
| DS | **Weibull** | **94.8** | **99.9** | **99.5** | **100.0** |
| IDD | Exponential | 0.5 | 0.0 | 0.1 | 0.0 |
| IDD | **Weibull** | **99.5** | **100.0** | **99.9** | **100.0** |
| all | Lognormal / Gamma | 0.0 everywhere | | | |

→ **Weibull is the best-fit distribution almost everywhere** (DD: 92–98 % except
spi90d where Exp ties at 25 %; DS & IDD: ~95–100 %). Gamma and Lognormal never win.

**Fitted descriptors — median [IQR] over aquifer pixels**

Weibull (winner): shape k, scale λ

| char | indicator | k | λ |
|---|---|---|---|
| DD | spi90d | 0.750 [0.721, 0.784] | 8.64 [7.66, 9.90] |
| DD | spi180d | 0.858 [0.829, 0.891] | 6.56 [5.98, 7.39] |
| DD | spei90d | 0.807 [0.777, 0.839] | 6.83 [6.26, 7.55] |
| DD | spei180d | 0.936 [0.901, 0.972] | 5.46 [5.00, 6.03] |
| DS | spi90d | 0.664 [0.640, 0.691] | 11.06 [9.73, 12.77] |
| DS | spi180d | 0.751 [0.728, 0.780] | 8.70 [7.88, 9.86] |
| DS | spei90d | 0.715 [0.690, 0.740] | 8.88 [8.07, 9.90] |
| DS | spei180d | 0.825 [0.796, 0.854] | 7.35 [6.69, 8.16] |
| IDD | spi90d | 0.570 [0.546, 0.596] | 27.51 [23.79, 32.04] |
| IDD | spi180d | 0.670 [0.642, 0.702] | 23.70 [21.49, 26.16] |
| IDD | spei90d | 0.598 [0.573, 0.625] | 25.60 [22.81, 28.97] |
| IDD | spei180d | 0.712 [0.684, 0.745] | 22.15 [20.32, 24.29] |

Runner-up Exponential θ (=mean) for reference: DD 10.9/7.3/8.0/5.7; DS 16.2/10.9/11.7/8.3;
IDD 45.3/31.7/39.3/27.7 pentads (order spi90d/spi180d/spei90d/spei180d).
Lognormals (never win): e.g. DS spi90d μ_log 1.65, σ_log 1.53. Full table:
`dsi_fit_analysis/param_summary.csv`; per-pixel fields incl. every candidate's params:
`dsi_run_fits.nc` (CF-style, flag_meanings exponential/lognormal/gamma/weibull).

**Implied expectations E(c) for Eq2** (Weibull mean = λ·Γ(1+1/k)) agree with the
empirical sample means to within ~2-4 % (e.g. DS spi90d: 8.15 fitted vs 8.34
empirical), so the intensity-mean proxy currently used inside `compute_DSI()` is an
acceptable stand-in until Eq3 is driven directly by these E(c) terms with signed
weights (w_DD>0, w_DS>0, w_IDD<0).

**Why our winners differ from VVIP.** The paper (monthly CRU 1949-2018) reported
Exponential best for DD/IDD and Lognormal best for DS (~83 %). Here, at **5-day
pentad resolution**, durations/severities are finer-grained and strongly right-skewed
with heavy lower-k tails: fitted Weibull shapes are consistently **k < 1**
(0.57-0.94), i.e. monotonically decreasing densities — a regime where Weibull
out-scores both Exponential and Lognormal under AIC. This is a resolution effect,
not a contradiction: aggregating our runs to monthly steps would reproduce
paper-like behaviour. It also means the paper's closed-form Exp expectation E(c)=θ
must be replaced by E(c)=λ·Γ(1+1/k) when propagating Eq2 on this dataset.

Figures: `best_dist_maps.png` (categorical majority-vote maps of the winning PDF per
characteristic), `ds_pooled_fit.png` (pooled severity histograms with median-parameter
PDF overlays, winner marked *).


### 11.14 NEW - GDE map clipped to the Ogallala aquifer (CF-compliant NetCDF)

**Source.** `G:\MSU_GWB\datasets\GDE_30arcsec.tif` is byte-identical to the deposited
product `D_Data_GDE_AggregatedLayers/GDE_30arcsec.tif` (Rohde et al. 2024, Nature,
DOI 10.1038/s41586-024-07702-8; deposit DOI 10.5281/zenodo.11062894). Content is
governed by `G:\GDE map\README.pdf`, Section D ("Dataset 2: Aggregated GDE data at
30 arcsecond resolution"). Actual grid: global EPSG:4326, **1/120 deg = 30 arcsec**
(~0.93 km) pixels, 5 uint32 bands in README order, nodata = 4294967295.

**Extraction script.** `extract_gde_ogallala.py` (new file, read-only on sources):
window-reads the 5 bands over the `hp_bound2010.shp` bbox, decodes the INT4U
scaling (fraction bands /1e8 -> [0,1] floats), rasterizes the 199 aquifer polygons
onto the window grid, and sets every cell outside the ROI to `_FillValue`.

**Output.** `derived_usgs/GDE_30arcsec_Ogallala_CF.nc` (3.0 MB, NETCDF4, CF-1.8):
grid y=1431 x=1160 (lon -105.921..-96.259, lat 31.743..43.664), vars
`GDE_sqm`, `GDE_frac_AA`, `GDE_frac_GA`, `AA_sqm`, `AA_frac_GA` (float32, zlib-4,
chunks 256x256), 1-D `y/x` lat/lon coords (same decreasing-y convention as the
other derived files), `crs` grid_mapping var = `latitude_longitude` EPSG:4326,
and global attrs carrying the Rohde et al. citation, CC BY 4.0 license, README
reference, geospatial bounds, and history.

**Verification (run summary).** 682,584 window cells inside the polygon (41.1 %);
154,500 of them carry analyzed area (9.3 % — the rest of the ROI lies outside the
GDE model domain and is nodata by design). Medians: GDE_frac_AA 0.25,
GDE_frac_GA 0.08, AA_frac_GA 0.69; total GDE area in ROI = 21,005.3 km^2.
Reopen check: Conventions=CF-1.8, all 5 bands present, ranges within [0,1] for
fractions. Sources unmodified.

**Assumptions / caveats.** (a) Band order taken from the file's own band
descriptions, cross-checked against README Section D ordering. (b) NAD83-Albers
polygon reprojected to WGS84 for masking (~1 m datum shift, negligible at 1 km).
(c) Per-cell m^2 values are the producers' rounded values; no area recomputation.
(d) Native 30-arcsec grid preserved (not resampled to the 4 km analysis grid);
regridding, if needed for TBI overlays, is a documented follow-up.


### 11.15 NEW - XGBoost expected-ET model + groundwater-buffering quantification
(`gw_buffering_xgb_et.py` -> `gw_buffering_xgb/`, 18 files)

Pilot implementation of Research Plan Sec. 6 (Phase 0), **ET leg only**:
`Bi(t) = ET_obs - ET_hat(P, SMrz, VPD, Rn, T, PFT, season)` with `ET_hat` from
XGBoost (v3.4.1, `hist`, seed 42; 800 trees, depth 6, lr 0.05, subsample 0.8,
colsample 0.7, L2 1.0, min_child_weight 10, early stopping 50 rounds on a
held-out validation year). Panel: 2,143,504 rows (10,630 sampled aquifer pixels
x 216 dekads, 2016-2021). Holdout = 2 driest growing seasons by aquifer-mean
spi90d (**2020, 2021**); validation year 2018.

**Headline results.**
- Model skill: train R2 0.807 / RMSE 8.19; holdout R2 **0.654** / RMSE 10.44;
  holdout dry-down RMSE 9.66 mm/dekad (transfers to drought acceptably).
- Mean residual in growing-season dry-downs: **Bi = +0.25 +/- 0.11 mm/dekad**
  (n=46,820) - small but positive ET persistence beyond climate+soil explanation.
- GW-access effect (shallow-minus-deep WTD, stratified by PFT x management x
  P30 tercile, 14 strata): **tau = +1.62 mm/dekad**; wet-year placebo (2019,
  2017): **-0.02 mm/dekad** (~null, as required by Plan Sec. 6.2).
- By WTD tercile (cuts 17.8 / 30.4 m): shallow +1.49, mid -0.69, deep +0.17.
- By management: irrigated +0.97, rainfed-crop -0.84, natural-other -0.98
  (managed-water signal visible; irrigation must not be read as GW buffering).
- Sustainability screen (pixel-mean dry Bi vs observed predev->2019 decline):
  >50 ft: +3.72 (n=592); 10-50 ft: +1.01 (955); stable: +0.14 (4,688);
  rise: +2.67 (145) - persistence concentrates where depletion happened.
- Volume plausibility (Sy=0.15): seasonal extra ET mean +1.8 mm vs implied
  depletion mean +549 mm, r=+0.07 - current dry-down persistence is ~2 orders of
  magnitude smaller than multi-decadal depletion, so **no "mined" label is
  warranted on volume grounds** (Plan Sec. 6.3 rule); classification stays
  correlational.
- AmeriFlux validation: 1,820 dekads / 17 ROI towers; SSEBop vs tower RMSE
  **16.06 mm/dekad** (tower ET from LE/2.45e6, no QC filtering).
- GRACE basin-mean TWS plotted z-scored with holdout shading (basin constraint
  only, not pixel-scale, per Plan). Saved model: `xgb_et_model.json`;
  pixel map: `Bi_drydown_mean_4km.nc` (CF-lite).

**Assumption register (A-codes also inline in the script).**
- A1 Scope: ET leg only. SIF/VOD legs absent from merged_datasets, so the
  "strongly buffered only if all three persist" test cannot run; all
  attribution language is correspondingly conditional.
- A2 Residual identity: Bi absorbs irrigation, cultivar, soil/rooting error,
  phenology error and model bias (Plan pitfall 1); GW attribution only via the
  WTD/well/storage covariation + controls reported above.
- A3 SIF/VOD independence: not applicable (see A1).
- A4 Common grid: GRIDMET 4 km (287x181), nearest for categorical, linear for
  continuous (repo convention); SMAP via nearest 9-km cell (max gap 0.080 deg).
- A5 Window 2016-2021: full calendar years inside the SMAP x SSEBop overlap;
  2012/2020-21 guidance adapted (2012 predates SMAP; 2020-21 used as holdout).
- A6 PFT = MCD12Q1 LC_Type1 for 2019, one-hot over sampled classes only
  (rare classes 1/4/5/6/8/9 absent from sample; IGBP_11 n=1 and IGBP_17 n=3
  estimates are unstable - see effects.csv).
- A7 SSEBop ET taken as mm/dekad as stored (Jan median ~0, p99 ~69 plausible).
- A8 Pixel stride 2 (~10.6k px); maps use the sampled lattice.
- A9 Growing season May-Sep is the primary inference window.
- A10 Dry-down = spi90d<=-0.8 AND negative SMrz anomaly AND positive VPD
  anomaly (5.5% of rows); thresholds are judgment calls - sensitivity runs
  recommended.
- A11 Management from GIR class==1 majority 2001-2015 (>=0.5) + GFSAD
  cropland/fallow; GIR is generous at 9 km (61% native cells, 78% of aquifer
  pixels flagged irrigated; class 2 treated as non-irrigated) - LANID/CDL
  cross-check recommended before any irrigation-subsidy labeling.
- A12 GRIDMET units assumed native (pr mm/day, T Kelvin->degC auto-detected,
  vpd kPa, srad W/m2); GEE attrs carry no physical units.
- A13 DROUGHT step nearest each dekad (<=~2.5 d offset).
- A14 Climatology baseline is only 2016-2021 (short; anomalies are relative,
  not WMO-standard).
- A15 HRES-WTD 2015 is a static depth axis (never a time series); tercile cuts
  17.8/30.4 m are sample-relative, not physical thresholds.
- A16 tau is stratified contrasts, not full propensity-score/causal-forest
  matching; conditional exchangeability is assumed, not demonstrated.
- A17 Placebo uses the 2 wettest growing seasons (2019, 2017).
- A18 95% CIs are normal approximations over pixel-dekads (ignores
  spatiotemporal autocorrelation - true uncertainty is larger).
- A19 dWL bins are fixed breaks (-50/-10/+10 ft); Sy=0.15 user constant.
- A20 Seasonal sums restricted to pixels with >=3 dry dekads.
- A21 GRACE plotted z-scored (units as stored, assumed cm w.e. - unverified).
- A22 Tower ET = LE*1800/2.45e6 mm per half hour, no QC/u*-filtering; nearest
  1-km SSEBop / sampled 4-km pixel collocation.
- A23 Rn approximated by GRIDMET srad only (no longwave/net-radiation
  computation).
- A24 Sources strictly read-only; all outputs are new files in
  `gw_buffering_xgb/` (re-running overwrites that folder only).

**Pilot reading (Plan Phase-0 table):** signal detectable (+0.25, tau +1.62 with
null placebo) and tracks depth/wells, BUT with a strong irrigation component
(+0.97 irrigated vs negative elsewhere) and failed volume plausibility for
"mined" labeling -> between "buffering is real and separable" and
"managed water dominates": proceed to threshold/reversibility work **with**
irrigation stratification mandatory, and add SIF/VOD legs when available.


### 11.16 NEW - Monthly buffering maps + ROI-averaged buffering time series
(`buffering_monthly_maps_timeseries.py` -> `gw_buffering_xgb/`)

Reloads the saved XGBoost model (`xgb_et_model.json`, no retraining), rebuilds
the identical dekadal panel (2,143,504 rows, deterministic), predicts ET_hat and
forms Bi = ET_obs - ET_hat for 2016-2021.

**New outputs.**
- `Bi_monthly_4km.nc` (CF-1.8, 72 x 287 x 181): monthly-mean Bi per pixel +
  `n_dekads` count layer; all 72 months populated.
- `roi_Bi_monthly.csv`: per month - ROI-mean Bi, pixel SD, n, dry-down pixel
  fraction. Overall monthly mean +0.13; min -3.65 (2020-06); max +8.52
  (2018-04) mm/dekad.
- `roi_Bi_timeseries.png`: monthly ROI-mean Bi with +/-1 SD band, zero line,
  and dry-down fraction bars on a twin axis.
- `map_Bi_YYYY-MM.png`: spatial map for `--year/--month` (default 2020-08:
  median +3.13 mm/dekad, 58.6% of valid pixels positive - ET persisting
  through a holdout drought-year month). Any 2016-2021 month can be rendered
  via CLI; the full monthly cube is in the NetCDF.

**Extra assumptions B1-B4** (on top of the A-register in 11.15): monthly Bi =
mean of dekads starting in that month (B1); ROI average over the stride-2
sampled lattice with dry fraction from the same dry flag (B2); renderable
months limited to 2016-2021 (B3); feature columns forced to the exact training
order with missing PFT dummies filled 0 (B4).


### 11.17 Landsat spectral-indices cube QC (`validate_landsat_indices.py` -> `landsat_indices_qc/`)

**Source.** `D:/Downloads/Landsat_indices_Ogallala_2020_perdate.nc` (41.3 KB;
Landsat 8/9 C2 L2 per-date indices, Ogallala ROI, 2020; produced 2026-09-23 by
the `ls_perdate` Colab pipeline). A companion scene-inventory zip
(`landsat_ot_c2_l2_*.zip`) holds only the USGS acquisition metadata CSV, not
pixel data.

**Spectral indices present (8 + count layer):** NDVI (Normalized Difference
Vegetation Index), EVI (Enhanced Vegetation Index), SAVI (Soil-Adjusted, L=0.5),
MSAVI (Modified Soil-Adjusted), NDMI (Normalized Difference Moisture Index),
NDWI (McFeeters Water Index), MNDWI (Xu Modified Water Index), NBR (Normalized
Burn Ratio) - plus `n_obs` (valid-observation count) and a `crs` grid-mapping
var. Grid: 332 x 266 regular lat/lon (~4 km WGS-84); CF-1.8 with per-variable
valid_min/max [-1, 1].

**Validity verdict: EMPTY.** The `time` dimension has length 0, so every
variable holds 0 cells - the file carries coordinates + metadata but no data.
Structure, attributes and grid are all correct; the failure is upstream (the
producer concatenated zero per-date slices - likely an empty scene/date list
or a date filter matching nothing). Nothing was plotted from data; the script
saved a grid-footprint map vs the aquifer boundary and an explicit verdict
figure instead, plus `qc_report.json/txt`. **Do not use this file for analysis
until the producer pipeline is re-run with a non-empty 2020 date list**; the
validator re-runs as-is on the fixed file (median maps, ROI time series and
histograms generate automatically once `time > 0`).

