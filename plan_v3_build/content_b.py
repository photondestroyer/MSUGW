"""Text of the plan, part B: Paper 2, gates, data to acquire, code changes, open questions, appendices."""
from engine import P, H1, H2, BUL, SP, TABLE, CARD, BOX, FIG, HFIG, ST, M
from content_a import pc


def paper2(D):
    S, SS, MI, WW, G, FIGS = D["S"], D["SS"], D["MI"], D["WW"], D["geom"], D["figs"]
    pan = MI["two_drought_panel"]
    n_both = D["n_both"]
    n3 = pan["<=3 m"]["n"]
    ch3 = pan["<=3 m"]["deeper_gt1m"] + pan["<=3 m"]["shallower_gt1m"]
    sh_n = MI["shallow_panel_n"]
    big2 = sum(pan[k]["abs_gt2m"] for k in ("<=3 m", "3-6 m", "6-10 m"))
    ms = MI["month_share_first200k"]
    rise = MI["dWL_predev_to_2019"]
    big = rise["largest_rise_patches"][0]
    wpy = MI["wells_per_year"]
    c = WW["by_measured_class"]
    s = []
    s += H1("5. Paper 2 in detail: the effect of groundwater access, its limits and its sustainability")
    s.append(P(
        "Paper 2 asks a causal question with observational data. The order of work is therefore: say exactly which effect "
        "is meant (P2.0); measure the exposure and its error (P2.1); show, before any outcome is read, that each design has "
        "the units and the power it needs; then estimate with several designs whose biases point in known directions."))

    s += CARD("P2.0", "Estimands, outcome and inference (new)",
              "Fix what is estimated, on which outcome, and how uncertainty is computed, before any result exists.", [
        dict(d="State the effects for natural and rain-fed land. **E1, total effect:** function during drought on land with "
               "groundwater access against land without, in the same climate, soil group and terrain class, whatever "
               "vegetation grows there. **E2, within-type effect:** the same contrast among blocks of one functional type, "
               "only for types found at both levels (grassland, dryland crops, hay). **E3, dose-response:** E1 as a function "
               "of depth to water, estimated with a doubly robust method for continuous exposures (Kennedy et al. 2017); it "
               "is the input of the threshold tests of Paper 3.",
             why="Wet meadows and riparian woodland exist because the water table is shallow. Matching on vegetation type "
                 "conditions on a consequence of the exposure, which removes part of the effect and can bias the rest "
                 "(Rosenbaum 1984). E1 answers RQ1 for the landscape. E2 answers it for the plant, and holds only under the "
                 "added assumption that nothing unmeasured decides both the type and its drought response.",
             check="Table of functional types by access class and climate decile. E2 is estimated only where both cells hold at "
                   "least 300 blocks.",
             fail="E2 is not reported for that type."),
        dict(d="Primary outcome per block and event: Δ = mean of (observed − expected) over the event months, minus the mean of "
               "the same quantity in the same calendar months of non-drought years. Secondary outcome: the buffering "
               "coefficient of P1.4. Ratios to a pre-event baseline (outcome O2 of revision 2) are not used.",
             why="The difference removes any constant site effect, and for rain-only vegetation it is zero on average if gate "
                 "G1b holds. Ratios are unstable where the baseline is small (winter, arid sites).",
             check="On held-out reference blocks, Δ has a mean of zero within ±0.15 null SD in every stratum.",
             fail="The stratum uses only designs that difference the bias out (P2.3, P2.4)."),
        dict(d="Uncertainty from a block bootstrap over 0.5° cells and events, checked against spatially robust standard errors "
               "(Conley 1999) and against randomisation inference in which the exposure field is shifted by 20 to 100 km. "
               "Before a design is run on outcomes, its minimum detectable effect is computed from the reference variance and "
               "the number of independent clusters.",
             why="Neighbouring blocks share weather and soil, so two million blocks hold far fewer independent observations. A "
                 "design that could not have detected a plausible effect must not be read as showing no effect (R12).",
             check="Gate G2b: minimum detectable effect at 80 % power below half of the effect size fixed in the registered plan "
                   "(proposal: 10 mm per month of ET, 0.5 null SD for the indices).",
             fail="The design is reported as uninformative for that stratum."),
    ], depends="P1.2 (gate G1b), P1.4.")

    s += CARD("P2.1", "Exposure: depth to water at event onset, with its error",
              "A depth to water per block that rests on measurements, uses nothing after onset, and carries its error.", [
        dict(d="Take depth to water from wells: water-table elevation kriged per year from spring measurements, with the DEM "
               "and the water-table map as drift (Desbarats et al. 2002). A block has a measured exposure only within 5 km of "
               "a well measured within two years. Other blocks enter only the designs that need regional depth, not block "
               "depth.",
             why=f"The map (Ma et al. 2026) is a random forest of climate, subsurface and terrain predictors trained on wells. "
                 f"Every design here controls for climate, soil and terrain; what remains of the map after that is mostly model "
                 f"error. Wells carry the information that the predictors do not. //Measured// against {WW['n']:,} aquifer wells "
                 f"(medians of 2000 to 2023; most were probably used to train the map, so this is favourable to it): "
                 f"correlation {WW['pearson']:.2f}, bias {WW['bias']:.1f} m, RMSE {WW['rmse']:.1f} m; RMSE {c['<=5']['rmse']:.1f} m "
                 f"where the measured depth is 5 m or less and {c['>30']['rmse']:.1f} m where it is more than 30 m. A map class "
                 f"of 5 m or less is right in {pc(WW['map_le5_and_measured_le5'], 0)} % of cases, and "
                 f"{pc(WW['measured_le5_found_by_map_le5'], 0)} % of wells with 5 m or less are mapped so.",
             risk="Wells are dense on irrigated uplands and sparse in rangeland and valleys.",
             check="Gate G2a: leave-one-well-out RMSE of at most 3 m for wells with 10 m or less. Table of blocks with a "
                   "measured exposure per state, year and climate decile.",
             fail="Wider depth classes. States without wells contribute only to the terrain-by-regional-depth design."),
        dict(d="Draw 50 realisations of the depth field by conditional simulation and average the effect estimates over them; "
               "the variance between draws is added to the uncertainty (Blackwell et al. 2017).",
             why="Error in a continuous exposure flattens a dose-response curve and moves a threshold. Averaging over plausible "
                 "exposures shows how much of the answer the error allows (R13).",
             check="A dose-response planted on the real well layout is recovered with a slope bias below 10 %.",
             fail="Class contrasts only, no curve."),
        dict(d="Onset depth is the last spring measurement before the event. Classes: 3 m or less, 3 to 6, 6 to 10, 10 to 30, "
               "more than 30 m. Primary contrast: 3 m or less against more than 10 m for herbaceous types, 6 m or less against "
               "more than 15 m for woody types. Rooting-depth maps are not used to define access.",
             why=f"//Measured//: well measurements cluster in January ({pc(ms[0])} %), March and April "
                 f"({pc(ms[2] + ms[3], 0)} %) and October ({pc(ms[9], 0)} %); each summer month holds 3 to 5 %. A spring value "
                 "exists for most wells; a summer value does not. The global rooting-depth product was derived with "
                 "water-table depth as an input (Fan et al. 2017), so defining access with it would put the exposure on both "
                 "sides of the comparison.",
             check="Time-stamp test: no measurement after onset enters. Results for neighbouring class limits (2 m, 4 m) are "
                   "reported."),
        dict(d="Before any outcome is read, publish the count of non-irrigated blocks in each cell of terrain (HAND of 5 m or "
               "less; more than 15 m) by depth (5 m or less; more than 15 m), per climate decile.",
             why=f"//Measured// from the map: never-irrigated land with water within 5 m is {pc(S[0]['nat_le2'] + S[0]['nat_2_5'])} "
                 f"to {pc(S[2]['nat_le2'] + S[2]['nat_2_5'])} % of the four driest deciles and "
                 f"{pc(S[4]['nat_le2'] + S[4]['nat_2_5'], 0)} to {pc(S[8]['nat_le2'] + S[8]['nat_2_5'], 0)} % of the six wetter "
                 "ones (Figure 2b). A design is feasible only where at least three of the four cells are filled.",
             check="At least 300 blocks in at least three cells, per decile.",
             fail="The decile is dropped from the design and listed."),
    ], depends="S0.5, S0.6; state well records (Section 7).")

    s += CARD("P2.2", "Design A: matching inside climate strata (descriptive baseline)",
              "A first estimate of E1 over the whole aquifer, with its expected bias stated in advance.", [
        dict(d="Match inside climate strata only on variables that the water table does not cause: climate normals, soil water "
               "capacity and texture, elevation, slope, HAND class, flow accumulation and event severity. Vegetation type is "
               "matched only for E2. Pre-event function is not a matching variable. Estimation by augmented inverse-probability "
               "weighting with cross-fitting (Chernozhukov et al. 2018).",
             why="Pre-event function is itself raised by groundwater access. Matching on it compares shallow sites with "
                 "unusually productive deep sites and removes part of the effect.",
             risk="With terrain matched, little independent variation in depth remains. The effective sample is small and sits "
                  "where wells depart from what terrain predicts.",
             check="Standardised mean differences below 0.1. Effective sample size after trimming reported. Wet-period placebo "
                   "indistinguishable from zero. Negative-control outcome (dormant-season NIRv). Shifted-exposure randomisation "
                   "p-value. E-value (VanderWeele & Ding 2017).",
             fail="Reported as descriptive only. The expected bias is upward, because soil and run-on favour low ground; this is "
                  "written down before the estimate is made."),
    ], depends="P2.0, P2.1.")

    s += HFIG(2, "Design B: terrain position crossed with depth to water (primary design of the natural thread)",
              FIGS + "/fig4_design.png",
             "Figure 4. The three landscape settings that Design B compares inside one climate decile and soil group. Left and "
             "centre: the difference between valley floor (B) and upland (A) is measured where the valley reaches groundwater "
             "and where it does not; the second difference measures the advantage of low ground alone. Right: on uplands with "
             "a raised water table, depth varies (A against C) while terrain position does not.")
    s += CARD("P2.3", "Design B: valley pairs and raised uplands",
              "Remove the advantage of low ground without having to model it.", [
        dict(d="Compare the valley-minus-upland difference in Δ where the valley floor reaches groundwater with the same "
               "difference where it does not, inside one climate decile and soil group:  τ_{B} = [Δ(B) − Δ(A)] in valleys that "
               "reach groundwater  −  [Δ(B) − Δ(A)] in dry valleys.",
             why="Low ground is wetter, finer in texture and more fertile with or without groundwater (Figure 1). The dry "
                 "valleys measure that advantage alone, and subtracting it leaves the part due to groundwater. This is a "
                 "negative-control exposure used as a difference-in-differences (Lipsitch et al. 2010; Sofer et al. 2016).",
             risk="The assumption is that the advantage of low ground is the same in both kinds of valley. Valleys with "
                  "permanent flow are wider, hold more alluvium and flood more often.",
             check="(1) Valleys paired on contributing area, slope, width and soil texture; standardised differences below 0.1. "
                   "(2) τ_{B} in wet periods indistinguishable from zero. (3) In non-drought months the run-on signature "
                   "(coupling slope at or above 1) is the same in both sets of valleys. (4) Blocks with open water in any month "
                   "are excluded; ET is used only for strips at least 200 m wide, because Landsat thermal pixels are 60 to "
                   "120 m. (5) Pairs are fixed before outcomes are read.",
             fail="If (2) or (3) fails, the valleys differ in more than groundwater. The design is then restricted to pairs "
                  "along one river (connected against disconnected reaches)."),
        dict(d="Classify reaches as reaching groundwater or not from wells within 1 km of the channel, compared with the "
               "channel elevation, as Jasechko et al. (2021) did, supported by the perennial or intermittent code of "
               "NHDPlus HR.",
             why="This uses measured water levels, not vegetation, so the exposure is not circular.",
             check="Agreement between the well-based class and the flow-permanence code of at least 80 %. Reaches where they "
                   "disagree are dropped."),
        dict(d="Add the contrast that needs no valley: uplands with a raised water table against uplands without (A against C "
               "in Figure 4), in the rise area of south-central Nebraska and in Sandhills interdunes.",
             why=f"Here depth varies and terrain position does not. //Measured//: {big['km2']:,} km² with rises of 10 to "
                 f"{big['max_ft']:.0f} ft south of the Platte River, and {rise['rise_ge10ft_km2']:,} km² with rises of 10 ft or "
                 "more in the whole aquifer.",
             risk="A rise of 55 ft from a deep level may still leave the water table out of reach. Most of the area is "
                  "irrigated.",
             check="Before outcomes are read: at least 300 non-irrigated blocks with a measured depth of 6 m or less inside the "
                   "area, and as many matched blocks outside.",
             fail="Used as a case study with exact permutation inference."),
    ], depends="P2.0, P2.1, S0.6; flow-permanence codes and stream elevations (Section 7).")

    s += CARD("P2.4", "Design C: the same block across events",
              "Remove everything about a site that does not change.", [
        dict(d="Keep the within-unit difference-in-differences of revision 2 (drought against non-drought months; access "
               "against no access) with an event-study plot.",
             why="It removes any constant site effect. It does not remove an advantage of low ground that grows in drought; "
                 "P2.3 is for that.",
             check="Leads before onset near zero. Because such pre-tests have low power (Roth 2022), the estimate is also "
                   "reported under violations of parallel trends as large as the leads allow (Rambachan & Roth 2023).",
             fail="The interval is widened by the sensitivity bound, and the headline does not rest on this design alone."),
        dict(d="Changing exposure: relate the change in Δ between the droughts of 2011 to 2013 and 2020 to 2023 to the change in "
               "spring depth at the same block.",
             why="A change of depth at a fixed place removes terrain, soil and vegetation history together. In principle this is "
                 "the cleanest design.",
             risk=f"Little variation. //Measured//: {n_both:,} wells were measured in both periods; {sh_n:,} had water within "
                  f"10 m in the first. Of the {n3:,} within 3 m, only {ch3} changed by more than 1 m, and {big2} of the {sh_n:,} "
                  "by more than 2 m. Nearly all are in Nebraska.",
             check="Minimum detectable effect computed first. Exact permutation inference.",
             fail="Reported as a case series. Contrasts within a year (spring against autumn levels) are added where both are "
                  "measured."),
    ], depends="P2.0, P2.1, S0.4.")

    s += CARD("P2.5", "Design D: the Sheridan-6 LEMA as a management experiment",
              "Use a pumping restriction as variation in water use that weather did not cause (irrigation thread).", [
        dict(d="Keep the difference-in-differences of revision 2 inside against outside the LEMA, before and after 2013, with "
               "outcomes from Paper 1 and a matched comparison ring.",
             why="The restriction changes pumping while climate and soils stay the same.",
             risk=f"//Measured//: the well file holds {wpy['Kansas']['2018']} to {wpy['Kansas']['2024']} Kansas wells a year after "
                  "2016, so water levels after the restriction are missing until the state records arrive. The class-conditioned "
                  "model starts in 2008, which leaves five years before treatment.",
             check="Pre-trends 2008 to 2012 with the sensitivity analysis of P2.4. Direction and size of the pumping change "
                   "agree with published estimates (Deines et al. 2019b; Drysdale & Hendricks 2018).",
             fail="Reported as descriptive."),
    ], depends="P1.4 to P1.7; Kansas records (Section 7).")

    s += CARD("P2.6", "Irrigated land: buffering and its limit (new)",
              "The access question of RQ1 for managed land: how far can pumping hold function through drought?", [
        dict(d="Estimate how the drought increment of irrigation excess depends on saturated thickness, within climate decile "
               "and crop group, across fields and within fields as thickness declines.",
             why="Irrigated fields hold function through drought by pumping more. The limit is well capacity, which falls as the "
                 "saturated thickness shrinks (Foster et al. 2015), and fields then move to dryland farming (Deines et al. "
                 "2020). Saturated thickness is to irrigated land what depth to water is to natural vegetation. It also gives "
                 "Paper 3 its threshold variable for the irrigation thread.",
             risk="Crop choice responds to well capacity (the mediator problem of P2.0 again). Farm management differs.",
             check="(1) First stage: reported pumping in Kansas rises in drought where thickness is large and not where it is "
                   "small. (2) Wet-year placebo. (3) Total effect (crop free to differ) and within-crop effect are reported "
                   "separately.",
             fail="Reported as an association."),
    ], depends="P1.4, P1.7; saturated-thickness grid (Section 7).")

    s += H2("P2.7 Robustness battery (pre-specified; every check is reported, including those that fail)")
    rows = [
        ["Wet-period placebo", "Confounding that does not depend on drought (soil, run-on).",
         "Effect in months with SPI-90 of +0.8 or more indistinguishable from zero."],
        ["Negative-control outcome", "Residual confounding and processing artefacts.",
         "No effect on dormant-season NIRv of grassland or on outcomes before onset."],
        ["Negative-control exposure", "Advantage of low ground.",
         "The valley-against-upland contrast over deep water everywhere is reported; it is the correction term of P2.3."],
        ["Shifted-exposure randomisation", "Spurious significance from spatial autocorrelation.",
         "The estimate lies beyond the 97.5th percentile of 500 estimates with the exposure field shifted by 20 to 100 km."],
        ["E-value, Rosenbaum bounds (Design A)", "Unmeasured confounding.",
         "Reported next to the strength of the measured confounders as a yardstick."],
        ["Leave one landscape or pair out", "One site driving the result.", "Sign unchanged in every run."],
        ["Scale: 30 m, 480 m, 960 m, 1/24°", "Mixing and scale dependence.", "Sign unchanged; size reported by scale."],
        ["Drought definition (SPEI; flash drought; threshold −1.0)", "Dependence on the event catalogue.", "Sign unchanged."],
        ["ET product (each OpenET model; smallest model value)", "Bias of one ET model.", "Sign unchanged; range reported."],
        ["Rain product", "Rain error.", "Estimate changes by less than 20 %."],
        ["Exposure draws (R13)", "Error in depth to water.", "Interval includes the variance between draws."],
        ["Reference rule (15 m against 30 m; GDE sensitivity set)", "Contaminated reference.", "Estimate changes by less than 20 %."],
        ["Estimator (corrected against raw residual; Δ against coupling)", "Artefact of one estimator.", "Sign unchanged."],
        ["Multiple testing", "False discoveries across strata and outcomes.",
         "False-discovery-rate control at 5 % within families fixed in the registered plan."],
    ]
    s += TABLE(["Check", "Bias it would reveal", "Pass criterion"], rows, [52, 52, 76],
               caption="Table 4. Robustness checks of Paper 2.")

    rows = [
        ["A  matching", "Wells that depart from what terrain predicts", "Climate, measured soil and terrain",
         "Run-on, unmeasured soil, exposure error", "Upward", "Small effective sample"],
        ["B  valley pairs", "Valleys that reach groundwater against dry valleys", "Climate; advantage shared by all low ground",
         "Differences between the two kinds of valley", "Upward if wet valleys flood more", "Wells near streams; pairs per decile"],
        ["B  raised uplands", "Canal seepage, interdunes", "Terrain position", "Soil differences, irrigation nearby", "Small",
         "Few non-irrigated blocks"],
        ["C  event study", "Timing of drought", "Constant site effects", "Advantage of low ground that grows in drought",
         "Upward", "Two large events"],
        ["C  changing depth", "Change of depth between droughts", "All constant site factors",
         "Land-use change that goes with the change of depth", "Small", f"{big2} wells with a change above 2 m"],
        ["D  LEMA", "Pumping restriction", "Weather, constant field effects", "Other changes in 2013", "Unknown",
         "Kansas wells after 2016 missing"],
        ["P2.6  thickness", "Well capacity", "Climate, crop group", "Management, farm size", "Unknown",
         "Reported use only in Kansas"],
    ]
    s += TABLE(["Design", "Variation in access it uses", "What it removes", "What can still bias it", "Expected bias",
                "Limit set by the data"], rows, [22, 36, 34, 40, 22, 30],
               caption="Table 5. The designs compared. A result counts as robust where designs with different biases agree "
                       "(gate G2c).")

    s += CARD("P2.8", "Sustainability classes from a water balance (replaces the rule of P2.7 in revision 2)",
              "Classify buffered systems from volumes with units, with a probability for each class.", [
        dict(d="Write a water balance per unit and year. For a unit of area A:  S_{y} · A · Δh = N − U,  where Δh is the yearly "
               "change of water level, U the consumptive use of groundwater from Paper 1 (excess attributed to "
               "groundwater-fed irrigation plus natural groundwater use), and N the net inflow (recharge plus lateral inflow "
               "minus baseflow). Regress Δh on U over the years.",
             why="The slope gives the specific yield and the intercept the net inflow, as Butler et al. (2016) showed with "
                 "reported pumping. The slope can be compared with the USGS specific-yield grid, which is an independent "
                 "check, and the intercept is the use at which water levels would be stable. Classes then follow from numbers "
                 "with units, not from signals that happen to coincide.",
             risk="Units too small for lateral flow to be negligible. Surface-water irrigation enters U and not Δh.",
             check="Gate G2d: (1) specific yield from the slope within the range of the USGS grid in at least 70 % of units. "
                   "(2) Lateral flux across the unit boundary from Darcy's law (conductivity polygons, water-table gradient, "
                   "thickness) below 10 % of U. (3) Units with more than 20 % of excess attributed to surface water are treated "
                   "separately. (4) The aquifer-wide sum against the GRACE storage trend, as a check of order of magnitude.",
             fail="(1) or (2): enlarge the units (management districts, then 0.5° cells)."),
        dict(d="Classes, each with a rule. **Sustainable buffering:** buffering shown (E1 or P2.6 above zero), U not above N, "
               "water level not falling. **Mined buffering:** buffering shown, U above N, water level falling, and the slope "
               "within a factor of 2 of the independent specific yield. **Surface-water subsidy:** excess attributed mainly to "
               "canal or river water. **Failed buffer:** access or irrigation present and Δ at or below zero in drought. "
               "**No buffer:** no excess above the detection limit. Class probabilities by Monte Carlo over specific yield "
               "(grid value and the range 0.10 to 0.20), OpenET model, attribution probability and water-level error; reported "
               "for 2008 to 2016 and for 2017 to 2025.",
             why="A class that flips inside the uncertainty of specific yield is not a result. Two periods show whether units "
                 "move between classes.",
             check="Share of units whose most probable class has a probability of at least 0.8. Share that changes class between "
                   "periods. Comparison with the Kansas water-balance studies (Butler et al. 2016; Whittemore et al. 2016) and "
                   "with published recharge (Scanlon et al. 2012).",
             fail="Units below 0.8 are mapped as undetermined."),
    ], depends="P1.6, P1.7, P2.1; specific-yield and saturated-thickness grids, state well records (Section 7).")
    return s


def back(D):
    S, SS, MI, WW, G = D["S"], D["SS"], D["MI"], D["WW"], D["geom"]
    s = []
    # ------------------------------------------------------------------ section 6
    s += H1("6. Gates and order of work")
    rows = [
        ["G1a", "P1.1", "At least 500 reference blocks per functional group and decile; ring test passed for the chosen rung; "
         "registry comparison in Kansas below 0.1 null SD.", "Stricter rung; registry-based exclusion; strata declared outside "
         "the domain."],
        ["G1b", "P1.2", "Bias in drought months within ±0.15 null SD on unseen droughts; coverage 0.88 to 0.92 per stratum; "
         "steps at regime boundaries below 0.1 null SD; known-zero sets at most 7 %.", "Bias correction; otherwise only "
         "differencing designs in that stratum."],
        ["G1c", "P1.4", "Zero-volume test; false-positive rate 4 to 6 %; planted effects recovered; detection-limit map made.",
         "Probability model refitted; limits only."],
        ["G1d", "P1.7, P1.8", "Excess against pumping (ratio and within-field correlation); pivot contrast within ±25 %; at "
         "least one independent sensor confirms.", "Results per ET model; claim limited to a Landsat-based measure."],
        ["G2a", "P2.1", "Leave-one-well-out RMSE at most 3 m for shallow wells; state records in hand or the state excluded; "
         "support table has three filled cells in each decile used.", "Wider classes; deciles dropped and listed."],
        ["G2b", "P2.0", "Minimum detectable effect below half of the effect of interest, per design.",
         "Design reported as uninformative."],
        ["G2c", "P2.2 to P2.5", "Designs B and C agree in sign for E1; placebos pass; Design A lies on the side its bias "
         "predicts.", "Result reported as unresolved and the failing assumption named (gate G2 of revision 2)."],
        ["G2d", "P2.8", "Specific yield from the regression within the grid range in at least 70 % of units; lateral flux below "
         "10 %.", "Larger units; class undetermined."],
    ]
    s += TABLE(["Gate", "Step", "Passes when", "If it fails"], rows, [12, 20, 92, 56],
               caption="Table 6. Gates. G1a to G1d replace gate G1 of revision 2; G2a to G2d replace gate G2.")
    s.append(P("**Order of work.**"))
    s += BUL([
        "**1.** Add histograms and observation support to the block aggregation (P1.0) before the aquifer-wide run, so the "
        "30 m files are read once.",
        "**2.** Request the state well and water-use records now (Section 7, item 1). Paper 2 cannot pass gate G2a without them "
        "in four states.",
        "**3.** Build the reference with the depth ladder and pass G1a. Train the folds and pass G1b before any target block is "
        "examined.",
        "**4.** Excess, probability and coupling (G1c); validation (G1d). Paper 1 is written from here.",
        "**5.** Paper 2: exposure, support table and power (G2a, G2b), all before outcomes are read; then the designs (G2c) "
        "and the water balance (G2d).",
    ])

    # ------------------------------------------------------------------ section 7
    s += H1("7. Data to acquire, in priority order")
    rows = [
        ["1", "Kansas WIZARD water levels and WIMAS water rights and reported use; Colorado DWR, Oklahoma OWRB and Wyoming "
         "water levels", "P1.1, P1.5, P1.7, P2.1, P2.4, P2.5, P2.8", "Wells after 2016 are absent from the holdings in these states."],
        ["2", "A second rain product: the NCEP Stage IV radar and gauge analysis (4 km, from 2002), or station records of the "
         "state mesonets", "P1.4", "Measures the rain error at block scale."],
        ["3", "AmeriFlux tower records, and image chips (Landsat, OpenET, CDL) at towers outside the polygon", "P1.7",
         "Mead towers are outside the tiles; the tower file is missing."],
        ["4", "USGS irrigation maps for 2018 to 2020 and well completion dates", "P1.1", "LANID ends in 2017."],
        ["5", "USGS saturated-thickness and specific-yield grids (McGuire et al. 2012)", "P2.6, P2.8", "Limit of irrigated "
         "buffering; check of the regression slope."],
        ["6", "Canals and flow-permanence codes from NHDPlus HR (same download as HAND), NWIS discharge, mapped alluvium",
         "P1.1, P2.3", "Canal buffer; reach classification."],
        ["7", "Caltech TROPOMI SIF completed to December 2022; SMAP files rebuilt to 2026", "P1.8", "Independent confirmation "
         "in the 2022 drought."],
        ["8", "MTBS fire perimeters", "P1.1", "Burned blocks leave the reference."],
        ["9", "The new water-table file HRES_WTD_Ogallala.nc copied from Drive", "Appendix B", "The checks here used the "
         "older file; repeat them."],
        ["10", "Optional: USDA Crop Sequence Boundaries (field polygons)", "P1.4", "Local contrast beyond the focus landscapes."],
    ]
    s += TABLE(["#", "Dataset", "Unblocks", "Why"], rows, [6, 86, 36, 52])

    # ------------------------------------------------------------------ section 8
    s += H1("8. What changes in the existing paper1 code")
    s.append(P("The package paper1 already builds P1.1 and P1.2 as revision 2 defined them. The decisions above change it as "
               "follows; nothing has been changed in the code yet."))
    rows = [
        ["step0/blocks.py", "Add the 16-bin histograms per target, the share of Landsat 8 and 9 observations, the mean "
         "observation day, the OpenET scene count and the smallest and largest model value."],
        ["config.py", "Depth rule per class and the ladder; calibration strata; common-mode radius; effect sizes for power."],
        ["p1_1/rules.py", "A 10 m depth bit; distance to canals; fire flag; depth rule by class; the three known-zero sets."],
        ["p1_1/verify.py", "Ring test; V4 representativeness table; registry comparison in Kansas; drift test after 2017."],
        ["p1_2/features.py", "Measurement inputs; class of the previous year; functional group with species as a secondary input."],
        ["p1_2/train.py, baselines.py", "Conformal margins per drought class and observation regime, with the fixed merge order."],
        ["p1_2/verify.py", "Bias inside drought months on unseen droughts; step test at regime boundaries; known-zero sets."],
        ["p1_4/ (new)", "Common-mode removal; two-groups probability; coupling slope; zero-volume test; detection-limit map."],
        ["tests/", "Planted subsidy that fills a share of each deficit; class-flip test; histogram recovery test."],
    ]
    s += TABLE(["File", "Change"], rows, [42, 138])
    s.append(P("The season rule in p1_2/calendars.py (season_window) is still open and is needed before the first real "
               "training run."))

    # ------------------------------------------------------------------ section 9
    s += H1("9. Open questions for you")
    s += BUL([
        "**Depth rule.** Accept 15 m for herbaceous classes, subject to the ring test, or keep 30 m everywhere and accept that "
        "the two wettest deciles lie outside the model domain?",
        "**Second rain product.** Stage IV (radar and gauges, from 2002) or mesonet stations, or both?",
        "**Core period.** Start in 2010 if the 2008 and 2009 crop layers fail their test, at the cost of two years before the "
        "2012 drought?",
        "**Focus landscapes.** Add the rise area of south-central Nebraska as an eighth focus landscape for Design B?",
        "**State records.** Who requests the Kansas, Colorado and Oklahoma well and water-use records, and is there an "
        "existing agreement?",
        "**Effect of interest.** Is 10 mm per month of ET the right smallest effect for the power calculation, or does the "
        "group have a management-relevant number?",
        "**Flux towers.** Which tower records does the group already hold, now that the earlier file is gone from the folder?",
    ])

    # ------------------------------------------------------------------ appendix A
    s += H1("Appendix A. Measured tables")
    rows = []
    for r in S:
        rows.append([str(r["decile"]), f"{r['ai_lo']:.2f} to {r['ai_hi']:.2f}", f"{D['P_dec'][r['decile'] - 1]:.0f}",
                     f"{r['wtd_med']:.1f}", pc(r["never_irr"]), pc(r["ref15"]), pc(r["ref30"]), pc(r["ref60"], 2),
                     pc(r["nat_le2"] + r["nat_2_5"]), pc(r["nat_5_10"]), pc(r["nat_10_30"])])
    s += TABLE(["Decile", "Rain / reference ET", "Rain, mm per year", "Median depth, m", "Never irrigated, %",
                "Reference bound 15 m, %", "30 m, %", "60 m, %", "Never irrigated and depth 5 m or less, %", "5 to 10 m, %",
                "10 to 30 m, %"], rows, [11, 22, 14, 14, 16, 18, 13, 13, 22, 14, 14],
               caption=f"Table A1. Support for the reference population and for exposed natural land by aridity decile. Sample of "
                       f"{SS['n_points']:,} points; depth from the older water-table file; irrigation from LANID 1997 to 2017. "
                       "Reference bound = never irrigated and mapped depth greater than the rule.", font=7.4)
    conf = WW["confusion"]
    labs = ["5 m or less", "5 to 10 m", "10 to 30 m", "more than 30 m"]
    keys = ["<=5", "5-10", "10-30", ">30"]
    rows = []
    for i, lab in enumerate(labs):
        b = WW["by_measured_class"][keys[i]]
        rows.append([lab] + [f"{v:,}" for v in conf[i]] + [f"{b['n']:,}", f"{b['bias']:+.1f}", f"{b['rmse']:.1f}"])
    s += TABLE(["Measured depth (well)", "Map: 5 m or less", "Map: 5 to 10 m", "Map: 10 to 30 m", "Map: more than 30 m", "Wells",
                "Bias, m", "RMSE, m"], rows, [34, 22, 22, 22, 26, 16, 16, 16],
               caption=f"Table A2. Water-table map against {WW['n']:,} wells inside the aquifer (well medians 2000 to 2023). Overall: "
                       f"correlation {WW['pearson']:.2f}, bias {WW['bias']:.1f} m (map shallower), RMSE {WW['rmse']:.1f} m, mean "
                       f"absolute error {WW['mae']:.1f} m. Most of these wells were probably used to train the map.", font=7.6)
    wpy = MI["wells_per_year"]
    yrs = ["2008", "2010", "2012", "2014", "2016", "2018", "2020", "2022", "2024", "2025"]
    rows = [[st] + [f"{wpy[st][y]:,}" for y in yrs] for st in ("Nebraska", "Texas", "Kansas", "Colorado", "New Mexico",
                                                               "Oklahoma", "South Dakota", "Wyoming")]
    s += TABLE(["State"] + yrs, rows, [26] + [15.4] * 10,
               caption="Table A3. Wells inside the aquifer with at least one usable depth measurement in the year "
                       "(HPA_wells_water_level_timeseries.nc).", font=7.6)
    pan = MI["two_drought_panel"]
    rows = []
    for k, lab in (("<=3 m", "3 m or less"), ("3-6 m", "3 to 6 m"), ("6-10 m", "6 to 10 m"), ("10-30 m", "10 to 30 m"),
                   (">30 m", "more than 30 m")):
        p_ = pan[k]
        rows.append([lab, f"{p_['n']:,}", f"{p_['median_change_m']:+.2f}", str(p_["deeper_gt1m"]), str(p_["shallower_gt1m"]),
                     str(p_["abs_gt2m"]), f"{p_['by_state']['Nebraska']:,}", f"{p_['by_state']['Texas']:,}"])
    s += TABLE(["Depth in 2010 to 2013", "Wells in both periods", "Median change to 2020 to 2023, m", "Deeper by more than 1 m",
                "Shallower by more than 1 m", "Change above 2 m", "In Nebraska", "In Texas"], rows,
               [30, 22, 28, 22, 24, 20, 18, 16],
               caption=f"Table A4. Wells measured in both 2010 to 2013 and 2020 to 2023 ({D['n_both']:,} wells), by depth in the "
                       "first period. Positive change = deeper.", font=7.6)
    rows = [[f"{p_['km2']:,}", f"{p_['lat']:.2f}° N, {abs(p_['lon']):.2f}° W", f"{p_['max_ft']:.0f}"]
            for p_ in MI["dWL_predev_to_2019"]["largest_rise_patches"]]
    s += TABLE(["Area, km²", "Centre", "Largest rise, ft"], rows, [40, 80, 40],
               caption="Table A5. The eight largest connected areas where the USGS map shows a rise of 10 ft or more from "
                       "predevelopment to 2019.", font=7.6)

    # ------------------------------------------------------------------ appendix B
    s += H1("Appendix B. How the measurements were made")
    s += BUL([
        f"**Scripts and results.** plan_v3_build/checks: check_support.py, check_wtd_wells.py, check_misc.py and their JSON "
        "outputs. All files were opened read-only; nothing in merged_datasets was changed.",
        f"**Reference support (Figure 2, Table A1).** Every 33rd pixel in both directions of HRES-WTD_2015_Ogallala.nc was read "
        f"(1,589,575 points with a value); {SS['n_points']:,} of them lie in lab cells with at least half their area in the "
        "aquifer and have a LANID value. Aridity is mean yearly rain divided by mean yearly reference ET of the lab gridMET "
        "file, 1985 to 2014; deciles are over the 27,333 aquifer cells. Not applied: HAND, the 800 m buffer, wetlands and "
        "class purity, so the shares are upper bounds.",
        "**Caveat on the water-table file.** The older file was used because the new one is on Drive only. In the older file "
        "the columns are 1.29 times as far apart as the rows; values sit at their stated coordinates, so a sample is not "
        "affected, but the checks should be repeated on HRES_WTD_Ogallala.nc. In either file about 80 % of the aquifer comes "
        "from state images coarser than the finest one, so the map does not resolve narrow valleys everywhere.",
        "**Product description of the map.** Read on 10 November 2026 from the Earth Engine community catalogue entry and the "
        "Zenodo record of the product: random forest, predictors from climate, subsurface properties and terrain, national "
        "test RMSE 14.9 m and correlation 0.79. To be confirmed in Ma et al. (2026) before quoting.",
        "**Map against wells (Table A2).** Wells inside the polygon; measurements of 2000 to 2023 with a depth between 0 and "
        "300 m; duplicates between sources removed; median per well; map value at the nearest pixel.",
        "**Wells per year and the two-drought panel (Tables A3, A4).** Same filters, 2000 to 2026. Months of measurement from "
        "the first 200,000 usable records of 2008 and later.",
        "**Rise areas (Table A5).** Cells of the 500 m USGS map with a change of +10 ft or more, grouped into connected areas.",
        f"**Flux towers.** Coordinates used: US-Ne1 41.1651° N, 96.4766° W; US-Ne2 41.1649° N, 96.4701° W; US-Ne3 41.1797° N, "
        "96.4397° W (to confirm on the AmeriFlux site pages). Distance to the boundary of hp_bound2010 in the Albers "
        "projection.",
        f"**Pivot geometry.** A circle of 400 m radius and a 480 m block of 16 by 16 pixels at a uniformly random offset. A "
        f"block with all pixel centres inside the circle exists for {G['p_inside']:.1f} % of offsets; one with at least 80 % of "
        f"its pixels inside for {G['p80']:.1f} %.",
        "**Synthetic run of the existing code.** The planted-subsidy numbers in P1.2 come from the local test run of the "
        "paper1 package on synthetic data, not from real tiles.",
    ])
    return s
