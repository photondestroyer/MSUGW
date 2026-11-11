"""Text of the plan, part A: front matter, problem statement, data, added rules, Paper 1."""
from engine import P, H1, H2, BUL, SP, TABLE, CARD, BOX, FIG, HFIG, ST, M
from reportlab.platypus import Paragraph, PageBreak


def pc(x, n=1):
    return f"{100 * x:.{n}f}"


def part_a(D):
    S, SS, MI, WW, G, FIGS = D["S"], D["SS"], D["MI"], D["WW"], D["geom"], D["figs"]
    s = []
    tw = MI["towers"]
    kmN1, kmN2, kmN3 = (tw[k]["km_to_boundary"] for k in tw)
    pan = MI["two_drought_panel"]
    n3, ch3 = pan["<=3 m"]["n"], pan["<=3 m"]["deeper_gt1m"] + pan["<=3 m"]["shallower_gt1m"]
    rise = MI["dWL_predev_to_2019"]
    big = rise["largest_rise_patches"][0]

    # ------------------------------------------------------------------ front matter
    s += [Paragraph(M("Groundwater buffering of vegetation function in the High Plains (Ogallala) aquifer"), ST["title"]),
          Paragraph(M("Papers 1 and 2: detailed analysis plan, with a reason and a check for every pipeline decision"), ST["sub"]),
          Paragraph(M("10 November 2026  ·  Detail for Sections 5 and 6 of revision 2 (20 October 2026)  ·  Step identifiers of "
                      "revision 2 are kept; new steps are P1.0, P1.8, P2.0 and P2.8  ·  Repository G:\\MSU_GWB\\datasets"), ST["meta"])]
    s.append(P(
        "**Summary.** The original plan asks whether access to groundwater keeps vegetation transpiring, photosynthesising and "
        "holding canopy water through drought, compared with similar vegetation without access, and whether that buffering is "
        "paid for by renewable water or by depletion. It also names the main difficulty: the part of vegetation function that "
        "rain does not explain is not automatically groundwater. Revision 2 answered with a rain-only expectation model "
        "(Paper 1) and four designs that stay inside one climate band (Paper 2). This document goes through every decision of "
        "those two papers. For each one it states why the decision is made, how it can fail, and which test, with a numeric "
        "pass criterion, would show the failure. Where the files in the repository allow it, the decision was tested on them "
        "on 10 November 2026; those results are marked //measured//."))
    s += BOX("Seven findings that change the plan", [
        "The unexplained part of function has seven components. Three of them (rain error at block scale, persistent site "
        "effects and run-on) can be as large as the groundwater signal. Each gets its own separation method (Section 1), and "
        "every headline number is produced by three estimators that fail in different ways (rule R9).",
        f"//Measured:// a 30 m depth rule leaves at most {pc(S[8]['ref30'])} % and {pc(S[9]['ref30'])} % of the two wettest "
        "climate deciles as reference land, before the terrain, buffer and purity rules remove more; a 60 m rule leaves none. "
        "The reference is therefore built with a depth rule per vegetation class, derived from rooting depth and map error, "
        "and each relaxation is tested before it is accepted (P1.1).",
        "The Landsat record has five observation regimes, and 2012, the anchor drought, was observed by Landsat 7 alone. "
        "Calibration and tests are made per regime, and 2012 is repeated with MODIS (P1.0, P1.2, P1.8).",
        f"//Measured:// the Mead flux towers (US-Ne1, US-Ne2, US-Ne3) lie {kmN3:.0f} to {kmN1:.0f} km outside the aquifer "
        "polygon, so the image tiles hold no pixels there, and the tower file is no longer in the data folder. Validation "
        "needs separate image chips (P1.7).",
        "The water-table map is a statistical model of climate, soil and terrain. Inside one climate and terrain class it "
        "carries little information of its own, so it cannot be the exposure of a design that controls for climate and "
        "terrain. Exposure in Paper 2 comes from wells, with its error carried into the estimates (P2.1), and the main design "
        "uses landscapes where terrain and depth vary separately (P2.3).",
        "Vegetation type is partly caused by groundwater access. Matching on it changes the question that is answered. "
        "Paper 2 therefore states two estimands (P2.0).",
        f"//Measured:// in the well file, Kansas, Colorado, Oklahoma and Wyoming have between 1 and 8 wells a year after 2016. "
        f"Of {n3:,} wells with water within 3 m in 2010 to 2013, only {ch3} changed by more than 1 m by 2020 to 2023. The "
        "state records become a precondition, and the within-site design with changing depth is a small-sample analysis (P2.4).",
    ])
    s.append(P("**How to read this document**", "body"))
    s += BUL([
        "Each step is a card. A card lists decisions. Each decision has the reason (Why), the way it can fail (Risk), the test "
        "with its pass criterion (Check, shaded) and what is done if the test fails.",
        "Pass criteria are proposals. They are frozen in the registered analysis plan (S0.7) before test data are examined. "
        "//Null SD// means the standard deviation of out-of-fold residuals of held-out reference blocks in the same stratum.",
        "//Measured// values were computed from files in merged_datasets with the scripts in plan_v3_build/checks. They are "
        "based on samples and use the older water-table file; Appendix B lists the method and the caveats.",
        "Numbers quoted from papers carry a citation and should be confirmed in the paper before reuse. The 19 references "
        "added here were checked against Crossref by DOI on 10 November 2026; the others were checked for revision 2.",
    ])
    s.append(P("**Contents.** 1 The problem as a measurable quantity · 2 Data re-examined · 3 Rules added · 4 Paper 1 in detail · "
               "5 Paper 2 in detail · 6 Gates and order of work · 7 Data to acquire · 8 Changes to the existing code · "
               "9 Open questions · Appendix A Measured tables · Appendix B How the measurements were made · References", "meta"))

    # ------------------------------------------------------------------ table 1
    s += H2("Table 1. Problems found in revision 2 and how this plan answers them")
    rows = [
        ["1", "The residual mixes seven things; only two are the target.", "Definition (Section 1).",
         "Separate each by its signature in time and space; three estimators.", "P1.4"],
        ["2", "Rain is known per 4 km cell; blocks are 480 m. A storm the grid misses looks like subsidy.",
         "gridMET design (Abatzoglou 2013).", "Remove the part of the residual shared with nearby rain-only blocks; test "
         "against a second rain product.", "P1.4"],
        ["3", "No reference land under the 30 m rule in the wettest fifth of the aquifer.",
         f"Measured: at most {pc(S[8]['ref30'])} % and {pc(S[9]['ref30'])} % of deciles 9 and 10.",
         "Depth rule per class with a tested ladder.", "P1.1"],
        ["4", "LANID is made from greenness, so excluding it trims the green tail of rain-only land.",
         f"Measured: {pc(SS['ever_irr_le5yrs'])} % of its footprint is labelled in 5 or fewer of 21 years.",
         "Compare with a registry-based reference in Kansas; frequency sensitivity; drift test after 2017.", "P1.1"],
        ["5", "A 480 m block mean mixes pivot circle and dry corners.",
         f"Geometry: a pure block exists for about {G['p_inside']:.0f} % of pivots.",
         "Keep within-block histograms.", "P1.0"],
        ["6", "Observation density and sensor change by period; 2012 has one satellite with scan gaps.",
         "Mission record; notebook settings.", "Measurement inputs, calibration per regime, MODIS repeat of 2012.", "P1.0, P1.2, P1.8"],
        ["7", "Flexible models predict too much function in extreme drought.", "Regression to the mean.",
         "Gate on bias inside drought months, tested on unseen droughts.", "P1.2"],
        ["8", "Summing max(0, residual) creates volume out of noise.", "Arithmetic.",
         "Probability-weighted signed excess; zero-volume test.", "P1.4"],
        ["9", "A site that is always better than predicted looks like natural subsidy.", "Unmodelled soil and management.",
         "Coupling slope; differencing within block.", "P1.4, P2.0"],
        ["10", "Blocks failing only the terrain rule were the known-zero set, but low ground gets run-on.", "Hydrology.",
         "New known-zero sets; the old one becomes the run-on diagnostic.", "P1.2"],
        ["11", "Mead towers outside the polygon; tower file missing.", f"Measured: {kmN3:.0f} to {kmN1:.0f} km outside.",
         "Image chips at towers; reacquire tower data.", "P1.7"],
        ["12", "ET and the Landsat indices share inputs; the original plan asks for independent signals.", "OpenET design.",
         "Confirmation with SIF, microwave optical depth, ECOSTRESS and MODIS.", "P1.8"],
        ["13", "The water-table map is a model of the confounders, with large error.",
         f"Measured at wells: RMSE {WW['rmse']:.1f} m.", "Well-anchored exposure with error draws.", "P2.1"],
        ["14", "Vegetation type lies on the path from groundwater to function.", "Ecology of wet meadows and riparian woodland.",
         "Total and within-type estimands.", "P2.0"],
        ["15", "Low ground is wetter with or without groundwater.", "Figure 1.",
         "Valley-against-upland difference, compared between valleys that reach groundwater and dry valleys.", "P2.3"],
        ["16", "Wells missing after 2016 in four states; shallow levels hardly changed between droughts.",
         "Measured (Tables A3 and A4).", "Acquire state records; treat the changing-depth design as small-sample.", "P2.4, Section 7"],
        ["17", "Sustainability classes rested on signals coinciding.", "Revision 2, P2.7.",
         "Water balance per unit: regression gives specific yield and net inflow.", "P2.8"],
    ]
    s += TABLE(["#", "Problem", "Evidence", "Answer in this plan", "Step"], rows, [6, 54, 40, 58, 20])

    # ------------------------------------------------------------------ section 1
    s += H1("1. The problem as a measurable quantity")
    s += H2("1.1 What the residual contains")
    s.append(P(
        "For a block and month, let //Y// be an observed function variable (ET from OpenET, NIRv or NDMI from Landsat) and "
        "//E// the value that the same vegetation, on the same soil, with the same weather history would show if rain were its "
        "only water. The residual //r = Y − E// is what the plan calls excess. It is a sum of seven parts:"))
    s.append(P("r  =  S_{irr}  +  S_{gw}  +  S_{sw}  +  e_{rain}  +  e_{site}  +  e_{model}  +  e_{obs}", "eq"))
    s.append(P(
        "Only the first two are the subject of the project: water applied by irrigation, and groundwater taken up by natural "
        "or rain-fed vegetation. The other five are present in every block-month and are not small. The plan can claim a "
        "groundwater signal only after it has shown how each of the five is removed or bounded. Table 2 gives the property "
        "that separates each part and the step that uses it."))
    rows = [
        ["**S_{irr}**\nirrigation", "Starts and stops abruptly within a season; grows in drought.", "Ends at the field edge; circles and rectangles.",
         "Field-edge contrast, pivot interior against corners, pumping records (P1.4, P1.7)."],
        ["**S_{gw}**\nnatural groundwater use", "Persistent; does not follow recent rain; largest relative to expectation in drought.",
         "Follows valleys and the water table.", "Coupling slope below 1 (P1.4); terrain-by-depth design (P2.3)."],
        ["**S_{sw}**\nsurface water and run-on", "Follows rain and floods; larger after wet spells.", "Low ground whatever the water table; near open water and canals.",
         "Coupling slope at or above 1; water indices; dry valleys as control (P1.4, P2.3)."],
        ["**e_{rain}**\nrain and demand error at block scale", "Changes sign from month to month.", "Shared by all blocks under one grid cell (4 km or more).",
         "Subtract the residual shared with nearby rain-only blocks; second rain product (P1.4)."],
        ["**e_{site}**\npersistent site effect", "Same sign every year; scales with rain.", "Specific to the block.",
         "Difference within the block (drought minus normal months); coupling slope (P1.4, P2.0)."],
        ["**e_{model}**\nmodel error", "Depends on drought severity, class and climate.", "Smooth.",
         "Calibration per drought class; test on unseen droughts (P1.2)."],
        ["**e_{obs}**\nobservation error", "Follows sensor period, number and day of observations.", "Scene and path patterns.",
         "Measurement inputs, calibration per regime, independent sensors (P1.0, P1.2, P1.8)."],
    ]
    s += TABLE(["Part", "Signature in time", "Signature in space", "How it is separated"], rows, [42, 42, 38, 58],
               caption="Table 2. The seven parts of the residual and the property that separates each.")
    s += H2("1.2 What each paper must deliver")
    s += BUL([
        "**Paper 1** delivers, for every block and month: a calibrated probability that function exceeds what rain alone "
        "supports; the excess in millimetres with an interval; the coupling of the block to rain; a source label where one can "
        "be defended; and a map of the smallest excess the method can detect. It answers the enabling question of revision 2.",
        "**Paper 2** delivers the effect of groundwater access on function during drought (RQ1) as two stated estimands, its "
        "change with depth to water, the limit of buffering on irrigated land, and sustainability classes from a water balance "
        "with intervals (RQ5).",
    ])
    s += H2("1.3 Why the comparison of shallow and deep sites is not enough")
    s += FIG(FIGS + "/fig1_dag.png",
             "Figure 1. Causes of function during drought on natural and rain-fed land. Depth to water is caused by climate, "
             "terrain and geology, and the same three also act on function through run-on, soil water and nutrients (orange). "
             "Vegetation type lies on the path from depth to function (green). Variation in depth that terrain does not cause "
             "comes from pumping, canal seepage and whether a stream is connected to the aquifer (blue, top). Climate also acts "
             "on vegetation type and function directly; every contrast in this plan is made inside one climate stratum, so "
             "those arrows are not drawn.")
    s.append(P(
        "Figure 1 has two consequences for Paper 2. First, comparing shallow with deep sites, even after matching on climate "
        "and soil, mixes the groundwater effect with the advantage of low ground. A design is needed in which low ground "
        "exists with and without groundwater (P2.3), or in which depth changes while the site stays the same (P2.4, P2.5). "
        "Second, vegetation type must be handled as a consequence of access, not as a background variable (P2.0)."))

    # ------------------------------------------------------------------ section 2
    s += H1("2. Data re-examined: the properties that decide the design")
    s.append(P("Table 3 lists, for each dataset that Papers 1 and 2 use, the property that can break a step and where the plan "
               "deals with it. Rows marked //measured// rest on the checks of 10 November 2026."))
    wpy = MI["wells_per_year"]
    rows = [
        ["Landsat 30 m monthly medians (2005 on)",
         "The monthly value is the median of one to three observations. Five regimes: 2005 to 2011 Landsat 5 and 7; 2012 to "
         "early 2013 Landsat 7 alone, with scan-line gaps; 2013 to 2016 Landsat 7 and 8; 2017 to 2021 Landsat 8 alone (the "
         "notebook drops Landsat 7 after 2016 because its orbit drifted); from 2022 Landsat 8 and 9 (Wulder et al. 2019).",
         "Noise and bias differ by regime. The change from TM and ETM+ to OLI falls between the 2012 drought and the wet years "
         "after it.", "P1.0, P1.2, P1.8"],
        ["OpenET v2.1 monthly ET",
         "ET is the satellite-derived fraction of reference ET on overpass days, interpolated in time and multiplied by "
         "gridMET reference ET (Melton et al. 2022). The count variables give the clear scenes behind each month; 0 means "
         "interpolated. Thermal pixels are 60 to 120 m. Accuracy is lower outside cropland (Volk et al. 2024).",
         "ET and the Landsat indices share inputs. Drought months without a clear scene are smoothed. Strips narrower than "
         "about 200 m are mixed with their surroundings.", "P1.0, P1.4, P1.7, P2.3"],
        ["CDL crop type (2008 on)",
         "Class is partly a consequence of water (crop choice). The 2008 and 2009 layers were made at 56 m (to confirm in the "
         "NASS metadata).", "A species-level expectation has no rain-fed support in dry deciles; class errors create "
         "residuals.", "P1.1, P1.2"],
        ["gridMET, 1/24° daily", "Rain from a gauge-based climatology with coarse reanalysis timing (Abatzoglou 2013); one "
         "value per cell of about 4 km.", "Local storms and local dry spots are not seen at block scale.", "P1.4"],
        ["Water-table depth map (Ma et al. 2026)",
         f"A random forest of climate, subsurface and terrain predictors trained on wells; a long-term mean. //Measured// "
         f"against {WW['n']:,} aquifer wells: correlation {WW['pearson']:.2f}, bias {WW['bias']:.1f} m, RMSE {WW['rmse']:.1f} m. "
         f"Where the map says deeper than 30 m, the well shows 10 m or less in {pc(WW['map_gt30_but_measured_le10'])} % of cases.",
         "Good enough to exclude shallow water from the reference. Not an independent exposure once climate and terrain are "
         "controlled.", "P1.1, P2.1"],
        ["Wells (HPA_wells_water_level_timeseries.nc)",
         f"//Measured// wells per year, 2012 against 2022: Nebraska {wpy['Nebraska']['2012']:,} and {wpy['Nebraska']['2022']:,}; "
         f"Texas {wpy['Texas']['2012']:,} and {wpy['Texas']['2022']:,}; Kansas {wpy['Kansas']['2012']} and {wpy['Kansas']['2022']}; "
         f"Colorado {wpy['Colorado']['2012']} and {wpy['Colorado']['2022']}; Oklahoma {wpy['Oklahoma']['2012']} and "
         f"{wpy['Oklahoma']['2022']}. Measurements cluster in January, March, April and October.",
         "No dynamic exposure in four states after 2016. Onset depth must be a spring value.", "P2.1, P2.4, P2.5"],
        ["LANID irrigation (1997 to 2017)",
         f"Derived from Landsat greenness (Xie et al. 2021). //Measured//: {pc(SS['ever_irrigated_share'])} % of the aquifer "
         f"sample is ever irrigated; {pc(SS['ever_irr_le5yrs'])} % of that in 5 or fewer years; almost no single-year labels.",
         "Selection on the outcome; nothing after 2017.", "P1.1"],
        ["USGS water-level change, predevelopment to 2019",
         f"A contour map. //Measured//: {rise['rise_ge10ft_km2']:,} km² rose by 10 ft or more; the largest area, "
         f"{big['km2']:,} km² near {big['lat']:.1f}° N, {abs(big['lon']):.1f}° W, rose by up to {big['max_ft']:.0f} ft. "
         f"{rise['decline_ge25ft_km2']:,} km² declined by 25 ft or more.",
         "Uplands with a raised water table exist: a contrast in depth that terrain did not cause.", "P1.1, P2.3"],
        ["GDE map (Rohde et al. 2024b)", "Fractions per 30 arc-second cell; positive in 97 % of natural cells (measured "
         "earlier); derived from satellite vegetation signals.", "Cannot be a selection rule or a strong label.", "P1.1, P1.6"],
        ["Flux towers", f"//Measured//: US-Ne1, US-Ne2 and US-Ne3 are {kmN1}, {kmN2} and {kmN3} km outside the polygon. The "
         "tower file listed in the older catalogue is gone.", "The planned irrigated and rain-fed tower pair is not in the tiles.", "P1.7"],
        ["SIF and microwave optical depth", "Independent of Landsat, but coarse. TROPOSIF in the folder ends in 2021; the "
         "Caltech TROPOMI file ends in October 2019 (provider archive to December 2022).", "Independent confirmation is "
         "possible only at the 1/24° scale; the 2022 drought needs the Caltech file completed.", "P1.8"],
    ]
    s += TABLE(["Dataset", "Property that matters", "Consequence", "Step"], rows, [30, 84, 50, 16],
               caption="Table 3. Dataset properties that decide the design.")
    s += HFIG(2, "2.1 How much rain-only reference land exists", FIGS + "/fig2_support.png",
             f"Figure 2. Measured on a sample of {SS['n_points']:,} points of the water-table map inside the aquifer. (a) Share "
             "of each aridity decile that is never irrigated in LANID and lies over a mapped water table deeper than 15, 30 or "
             "60 m. This is an upper bound: the terrain, buffer, wetland and purity rules remove more. (b) Share that is never "
             "irrigated and lies over a water table of 10 m or less. Deciles are tenths of the aquifer area ordered by the "
             "ratio of rain to reference ET, 1985 to 2014. Full numbers in Appendix A.")
    s.append(P(
        f"Figure 2a shows that the depth rule decides where the method can work. Under the 30 m rule of revision 2 the bound is "
        f"{pc(S[0]['ref30'], 0)} to {pc(S[1]['ref30'], 0)} % in the five driest deciles, {pc(S[7]['ref30'], 0)} % in decile 8, and "
        f"{pc(S[8]['ref30'])} % and {pc(S[9]['ref30'])} % in deciles 9 and 10. A decile holds about 200,000 blocks, so decile 10 "
        f"has at most about {round(200000 * S[9]['ref30'], -2):,.0f} candidate blocks before any other rule, to be shared among ten classes. "
        f"Under 60 m, deciles 7 to 10 have less than 1 %. Under 15 m the bound in deciles 9 and 10 is "
        f"{pc(S[8]['ref15'], 0)} % and {pc(S[9]['ref15'], 0)} %. Figure 2b shows the other side: land where natural groundwater "
        f"use can occur is rare in the dry west (under 2 % of deciles 1 to 4 has water within 5 m) and common in deciles 6 "
        f"to 10. The region with the most exposed land is the region with the least reference land. P1.1 and P1.3 are built "
        "around this fact."))
    return s


def paper1(D):
    S, SS, MI, WW, G, FIGS = D["S"], D["SS"], D["MI"], D["WW"], D["geom"], D["figs"]
    tw = MI["towers"]
    kmN1, kmN2, kmN3 = (tw[k]["km_to_boundary"] for k in tw)
    rise = MI["dWL_predev_to_2019"]
    big = rise["largest_rise_patches"][0]
    s = []
    # ------------------------------------------------------------------ section 3
    s += H1("3. Rules added to R1 to R8 of revision 2")
    s += BUL([
        "**R9 Three estimators.** Every headline quantity is produced by three estimators with different weaknesses: a bound "
        "from the water balance, the rain-only model, and a local contrast that needs no model. A result is reported as robust "
        "only where they agree (Munafò & Davey Smith 2018).",
        "**R10 The null is measured.** Every threshold on a residual is calibrated on held-out reference blocks of the same "
        "class, climate, season, drought class and observation regime. No normal distribution is assumed.",
        "**R11 Every threshold has a derivation, a ladder and a test.** A number such as 30 m is stated with the reasoning "
        "that produced it, the alternatives fixed in advance, and a test on data that can reject it.",
        "**R12 Detection limit before detection, power before test.** No absence of an effect is reported without the smallest "
        "effect the step could have seen.",
        "**R13 Exposure error is carried.** Depth to water enters with an error distribution, and effect estimates are "
        "averaged over draws from it (Blackwell et al. 2017).",
        "**R14 Stop rules.** A failed check halts the step that depends on it. The fallback is written here, before results "
        "exist.",
    ])

    # ------------------------------------------------------------------ section 4
    s += HFIG(1, "4. Paper 1 in detail: when and where vegetation uses water that rain did not supply",
              FIGS + "/fig3_pipeline.png",
              "Figure 3. Steps of Paper 1 with the main check of each (grey) and the four parts of gate G1.")

    s += CARD("P1.0", "Block aggregates that keep what later steps need (new; extends the output of S0.1)",
              "Read the 3.5 TB of 30 m files once, and keep per 480 m block and month everything the later steps need.", [
        dict(d="Store a 16-bin histogram of the 30 m values of ET, NIRv and NDMI per block and month, next to the mean, the "
               "standard deviation and the count.",
             why=f"A quarter-section pivot circle has a radius of about 400 m. A 480 m block lies wholly inside such a circle "
                 f"only if its centre is within about 80 to 105 m of the pivot centre, which is the case for about "
                 f"{G['p_inside']:.0f} % of pivots. For about {G['p80']:.0f} % a block exists that is at least 80 % inside; "
                 "every other block that the pivot touches is a mixture of circle, dry corners and roads (geometry, "
                 "Appendix B). Block means therefore dilute the excess and blur its area. The histogram keeps the "
                 "distribution inside the block at about 40 GB before compression, near 1 % of the 30 m volume, so irrigated "
                 "fractions can be recovered without reading the tiles again.",
             risk="Bin edges fixed in advance may be too coarse where values cluster.",
             check="In three tiles, the irrigated fraction per block recovered from the histograms (share of pixels above the "
                   "pixel-level 95 % reference threshold, minus the 5 % expected by chance) against the LANID fraction of the "
                   "same block, 2008 to 2017: R² at least 0.7 and slope between 0.8 and 1.2.",
             fail="Use 32 bins. If the test still fails, classify at 30 m in the tiles concerned and state the added cost."),
        dict(d="Carry the observation support of every block-month: number of Landsat observations, share from Landsat 8 or 9, "
               "mean day of month of the observations, OpenET clear-scene count, and the smallest and largest OpenET model "
               "value.",
             why="The five observation regimes of Table 3 change the noise and the bias of the targets. OpenET fills months "
                 "without a clear scene from neighbouring dates. The model and the calibration can account for this only if "
                 "the support is stored with the value.",
             check="Aggregates of one tile equal an independent recomputation. Table of median observations per month by regime "
                   "and growing-season month."),
        dict(d="A block-month target is valid when at least half of the block's aquifer pixels are valid. The valid share is "
               "reported per year and month.",
             why="In 2012, scan-line gaps and a single satellite reduce coverage. A hidden drop in coverage would tilt the "
                 "drought year towards the pixels that happened to be seen.",
             check="Share of valid block-months in April to September 2012, per class and aridity decile, at least 70 %; "
                   "compared with 2014 to 2016.",
             fail="For strata below 70 %, results for 2012 come from OpenET, which fills gaps, and from MODIS at the 1/24° scale, "
                  "and are marked."),
    ], depends="S0.1 (tile frames and block grid).")

    s += CARD("P1.1", "Reference population: vegetation that can use only rain",
              "A yearly set of blocks whose only water is rain, selected without the vegetation signal, and large enough in "
              "every climate decile and class.", [
        dict(d="Replace the single 30 m depth rule by a rule per vegetation class: mapped depth to water deeper than 15 m for "
               "herbaceous classes (grassland and pasture, annual crops, fallow) and deeper than 30 m for shrubland, hay and "
               "alfalfa, and any block with woody cover. The ladder 60, 30, 15, 10 m is fixed in advance.",
             why=f"The rule must exceed root reach plus capillary rise plus map error. Reported biome means of maximum rooting "
                 f"depth are about 2.1 m for crops and 2.6 m for temperate grassland, against about 5 m for shrubland (Canadell "
                 f"et al. 1996), with single woody species far deeper (Schenk & Jackson 2002). Below an extinction depth of a "
                 f"few metres, set by soil and cover, groundwater stops contributing to ET (Shah et al. 2007). //Measured// map "
                 f"error at wells: where the map says deeper than 15 m, the well shows 5 m or less in "
                 f"{pc(WW['map_gt15_but_measured_le5'])} % of cases and 10 m or less in {pc(WW['map_gt15_but_measured_le10'])} %; "
                 f"where it says deeper than 30 m, 10 m or less in {pc(WW['map_gt30_but_measured_le10'])} %. //Measured// support "
                 f"(Figure 2a): the 30 m rule leaves at most {pc(S[8]['ref30'])} % and {pc(S[9]['ref30'])} % of deciles 9 and 10, "
                 f"the 15 m rule {pc(S[8]['ref15'], 0)} % and {pc(S[9]['ref15'], 0)} %.",
             risk="Deep-rooted individuals or perched water in the 15 to 30 m band. Wells sit mostly on irrigated uplands, so the "
                  "error rates may be too favourable for valleys.",
             check="Ring test. Train on the stricter rung and predict, out of fold, the blocks that only the looser rung admits. "
                   "Per class and aridity decile the ring must match held-out blocks of the stricter rung: mean residual within "
                   "±0.2 null SD (two one-sided tests); share of block-seasons above the 95 % threshold at most 7 %; coupling "
                   "slope (P1.4) within ±0.1 of 1. The vegetation signal tests the rule as a whole; it never selects single "
                   "blocks, so the selection stays free of the outcome.",
             fail="Keep the stricter rung in that stratum. Where fewer than 500 reference blocks remain, the class and decile "
                  "are declared outside the model domain and handled by P1.3."),
        dict(d="Keep the terrain rule (HAND above 10 m in at least 90 % of the block's pixels) and store flow accumulation as a "
               "second terrain measure.",
             why="Low ground receives run-on and may reach an alluvial water table that the regional map misses. HAND is "
                 "measured against mapped channels whether or not they carry water, so in dry valleys the rule removes land "
                 "that is in fact rain-only. That costs sample, not purity.",
             risk="Run-on on gentle slopes above 10 m.",
             check="On held-out reference blocks, residual against HAND and against flow accumulation: slope not different from "
                   "zero, and the difference between the lowest and highest fifth below 0.1 null SD.",
             fail="Raise the HAND threshold to 15 m or add a flow-accumulation limit; report the change in counts."),
        dict(d="Exclude blocks within 800 m of any pixel that LANID ever labelled irrigated (1997 to 2017), and measure the "
               "bias that this greenness-derived mask introduces, using water-right registries where they exist.",
             why=f"The 800 m buffer equals one pivot diameter and covers misregistration. The difficulty is that LANID labels "
                 f"come from Landsat greenness (Xie et al. 2021). Removing them also removes the greenest dryland fields and "
                 f"lowers the upper tail of the reference distribution. The 95 % threshold then sits too low and rain-only land "
                 f"is flagged as subsidised. //Measured//: {pc(SS['ever_irrigated_share'])} % of the aquifer sample is ever "
                 f"irrigated in LANID, and {pc(SS['ever_irr_le5yrs'])} % of that area is labelled in 5 or fewer of the 21 "
                 "years, which is where occasional mislabelling would sit.",
             risk="Thresholds biased low. New irrigation after 2017 enters the reference unnoticed.",
             check="(1) In Kansas, build the reference twice, from water rights with place of use (WIMAS) and from LANID, and "
                   "compare the 90 % and 95 % quantiles of NIRv and ET in the same class, decile and month: difference below "
                   "0.1 null SD. (2) Repeat with exclusion only where LANID frequency is at least 3 years; report the change in "
                   "thresholds. (3) Drift: the mean residual of reference blocks by year in 2018 to 2025 has no upward trend "
                   "relative to 2008 to 2017 (interval of the slope includes zero), by state.",
             fail="(1) or (2): use registry-based exclusion in Kansas and Nebraska and apply the measured threshold correction "
                  "elsewhere, flagged. (3): add the USGS irrigation maps for 2018 to 2020 and well completion dates; drop the "
                  "affected years of that state from training."),
        dict(d="Exclude blocks with NWI wetland or riparian pixels, with open water in any month of the year, or within 800 m of "
               "a canal or ditch (NHD feature type 336).",
             why=f"Canal seepage raises water tables under uplands. //Measured//: the USGS map shows {rise['rise_ge10ft_km2']:,} "
                 f"km² where water levels rose by 10 ft or more since predevelopment. The largest area, {big['km2']:,} km² near "
                 f"{big['lat']:.1f}° N, {abs(big['lon']):.1f}° W with rises up to {big['max_ft']:.0f} ft, lies south of the Platte "
                 "River in the surface-water irrigation districts of south-central Nebraska (cause to be confirmed in McGuire "
                 "2017). Such land looks like deep-water-table upland in the terrain layers.",
             check="Reference residuals against distance to the nearest canal: flat beyond 800 m.",
             fail="Widen the buffer to the distance at which the residual becomes flat."),
        dict(d="A block-year needs one dominant class group with purity of at least 0.8. A crop qualifies only in counties where "
               "at least half of that crop's pixels lie outside the LANID footprint. The class of the previous year is stored.",
             why="Crop choice follows water, so a crop that is mostly irrigated in a county is unlikely to be rain-fed in a given "
                 "block. The previous class matters because a fallow year stores soil water that the next crop uses (wheat "
                 "after fallow); a model without it would read that water as excess.",
             risk="Class errors create residuals. The 2008 and 2009 layers are coarser.",
             check="(1) Residual SD of reference blocks in 2008 and 2009 against 2010 to 2012: ratio below 1.15. (2) Class-flip "
                   "test: relabel 5 % of reference blocks with the class they are most often confused with; the false-positive "
                   "rate stays below 7 %. (3) Wheat after fallow and wheat after wheat both have a mean residual within ±0.2 "
                   "null SD.",
             fail="(1): start the core period in 2010. (2): raise purity to 0.9 or filter with the NASS confidence layer. (3): add "
                  "months since the last crop as an input."),
        dict(d="The GDE map stays a diagnostic (a reference set with GDE fraction of at most 0.10 is run as a sensitivity). "
               "Burned blocks leave the reference in the fire year and the next (MTBS perimeters, to acquire).",
             why="The GDE map marks a positive fraction in 97 % of the natural cells of the aquifer, so a hard rule would empty "
                 "the natural reference. It is also derived from satellite vegetation signals (Rohde et al. 2024b), which rule "
                 "R2 does not allow in selection. Fire causes a drop and a flush of greenness unrelated to water.",
             check="Thresholds and coverage with and without the GDE sensitivity set differ by less than 0.1 null SD.",
             fail="Both sets are reported side by side in Papers 1 and 2."),
        dict(d="Keep verifications V1 to V3 of revision 2 and add V4, a table of how well the reference represents the land it is "
               "applied to.",
             why="A clean reference is of no use if it does not resemble the target land.",
             check="V1: at least 500 reference blocks per functional class and aridity decile. V2: water-year ET at most rain "
                   "plus soil storage change plus 10 % in at least 95 % of reference block-years. If V2 fails, the failures are "
                   "mapped: clustered near irrigation means contamination; uniform within a class means bias of the ET or rain "
                   "product; the same test at rain-fed flux towers separates the two. V4: standardised mean difference of soil "
                   "water capacity, clay, elevation, slope, rain and reference ET between reference blocks and all "
                   "non-irrigated blocks, per decile, below 0.25.",
             fail="V1 or V4: the stratum goes to P1.3 as outside the domain. V2 through product bias: the mass penalty and the "
                  "bound estimator use the smallest OpenET model value."),
    ], depends="P1.0, S0.5 (depth to water), S0.6 (HAND); registries and fire perimeters (Section 7).")

    s += CARD("P1.2", "Rain-only expectation model",
              "For each block and month, the distribution of ET, NIRv and NDMI that rain-only vegetation would show, unbiased "
              "inside droughts.", [
        dict(d="Targets as in revision 2. Inputs as built (driver history, soil, elevation, class, phenological clock, lagged "
               "pre-season state), plus measurement inputs: share of Landsat 8 or 9 observations, number of observations, mean "
               "observation day, OpenET scene count; and the class of the previous year. Exposure variables stay refused.",
             why="Measurement inputs describe how a value was observed, not whether water is available, so rule R2 holds. "
                 "Without them the sensor change of 2013, which falls between the 2012 drought and the wet years after it, "
                 "would appear as a change in vegetation.",
             risk="A measurement input could act as a proxy for location (zones where satellite paths overlap).",
             check="(1) Step test on held-out reference residuals at the regime boundaries 2011/2012, 2013, 2016/2017 and "
                   "2021/2022: step below 0.1 null SD per target. (2) Residual against each measurement input: flat. (3) "
                   "Ablation without measurement inputs reported.",
             fail="Add an offset per regime fitted on calibration blocks. If a step remains, event contrasts that span that "
                  "boundary are not used."),
        dict(d="The class input is a functional group (winter cereal, summer row crop, perennial grass, shrub, hay, fallow). The "
               "crop species is a secondary input.",
             why="Rain-fed maize hardly exists in the driest deciles, where much irrigated maize grows. A species-specific "
                 "expectation would be pure extrapolation there. A functional group has reference blocks in more strata, and the "
                 "water-balance penalty bounds the answer.",
             check="Reference blocks per functional group and decile (V1). Skill with and without the species input.",
             fail="Strata without support give the water-balance bound only (P1.3)."),
        dict(d="Quantile loss for the 5, 50 and 95 % quantiles with heads that cannot cross, plus the water-year mass penalty "
               "(ET may not exceed rain plus soil storage change, with 10 % tolerance). A sensitivity model predicts ET divided "
               "by reference ET.",
             why="The penalty keeps predictions physically possible where the model extrapolates. OpenET multiplies its "
                 "satellite fraction by gridMET reference ET, which is also a model input, so part of the skill on ET is built "
                 "in. The ratio target removes that part and leaves the limitation by supply.",
             check="Predicted reference block-years violating the bound: below 1 %. Skill of the ratio model reported.",
             fail="Raise the penalty weight within the range fixed in the configuration; report the sensitivity."),
        dict(d="Five spatial folds of 0.5° groups crossed with three six-year groups, one-year embargo, inner split for early "
               "stopping and calibration. One model per fold predicts every block of that fold, reference or not.",
             why="Every residual is out of fold (R3), and calibration is measured on the predictor that is deployed.",
             risk="Residuals correlated over more than one group or more than one year.",
             check="Variogram of out-of-fold reference residuals: range of the correlated part below 50 km. Autocorrelation of "
                   "yearly residuals below 0.2 at the embargo lag.",
             fail="Use 1° groups or a two-year embargo."),
        dict(d="Conformal margins per target, functional group, aridity tercile, season, drought class (SPI-90 at most −0.8; "
               "between; at least +0.8) and observation regime. Strata with fewer than 1,000 calibration block-months are "
               "merged upward in a fixed order (regime first, then tercile).",
             why="A 90 % interval that is right on average can be much too narrow in drought months and too wide in wet months. "
                 "Paper 2 looks only at drought months, so calibration must hold there.",
             check="Coverage of the 90 % interval between 0.88 and 0.92 in every stratum with at least 1,000 test block-months.",
             fail="The stratum is named in the report and its intervals are widened to the measured coverage."),
        dict(d="Gate G1b. Test the model on droughts it has not seen, and require a mean residual of zero on reference blocks in "
               "drought months.",
             why="Flexible models pull predictions towards the average, so in an extreme drought they tend to predict too much "
                 "function. Rain-only blocks then show negative residuals, and every buffered block looks less buffered than it "
                 "is. In Paper 2 this bias would enter every estimate with the same sign. Under the fold design, 2012 is "
                 "predicted by models trained on 2015 to 2025, and 2022 by models trained on 2008 to 2018.",
             risk="The two droughts differ in region and season, so each test covers different climates.",
             check="Mean out-of-fold residual of held-out reference blocks in months with SPI-90 at most −0.8: within ±0.15 null "
                   "SD, and within ±3 mm per month for ET, per functional group and aridity tercile, for 2011 to 2013 and for "
                   "2022 separately. The same for the wettest months.",
             fail="Fit a bias correction as a monotone function of SPI on calibration blocks and test again on test blocks. If "
                  "it still fails, Paper 2 uses for that stratum only designs that difference the bias out (P2.3, P2.4) and the "
                  "local contrast estimator."),
        dict(d="Keep the benchmarks M0 (the block's own climatology) and M2 (rain-use efficiency), the planted-subsidy test in "
               "two variants, and the ablation without pre-season inputs.",
             why="The first local run on synthetic data showed what the second variant is for. A subsidy planted in the season "
                 "months only was recovered within 1 %. One planted in every month was misestimated by 12 to 31 %, because the "
                 "pre-season inputs carried part of it into the expectation.",
             check="M1 beats M0 in every aridity decile with at least 100 test block-months and beats M2 overall. Planted "
                   "season-only subsidy recovered with bias below 10 %.",
             fail="If pre-season inputs absorb more than 10 % of a persistent subsidy, the model without them is used for blocks "
                  "flagged persistent in P1.4, and both versions are reported."),
        dict(d="Replace the known-zero set of revision 2 (dryland crop blocks that fail only the HAND rule) by three sets: "
               "spatially held-out reference blocks; the accepted ring of the depth ladder; and, in Kansas, fields with no "
               "water right and no LANID label. The HAND-only set becomes the run-on diagnostic set.",
             why="Low ground receives run-on. Blocks that fail only the terrain rule are expected to exceed the rain-only "
                 "expectation after wet spells. Used as a zero, they would make a correct model look wrong.",
             check="Share of block-months above the 95 % quantile at most 7 % in each known-zero set. In the run-on set any "
                   "excess must scale with rain (coupling slope of at least 1, P1.4).",
             fail="A known-zero set above 7 % stops gate G1b. The ring test tells contamination from model error."),
    ], depends="P1.1, S0.3 (drivers).")

    s += CARD("P1.3", "Transfer: where the expectation is supported",
              "Know for each block and month whether reference data support the expectation, and keep error rates valid when "
              "target land differs from reference land.", [
        dict(d="Area of applicability per block-month in the importance-weighted input space (Meyer & Pebesma 2021), and a "
               "classifier that tries to tell reference from target blocks in each stratum.",
             why="Reference land is upland over a deep water table. Valleys and irrigated fields differ in soil and terrain. "
                 "Outside the range of the reference data the model output is a guess.",
             check="Share of target block-months inside the area, per functional group and decile, published. A classifier "
                   "accuracy near 0.5 means no shift.",
             fail="Outside the area, statements are limited to the water-balance bound (estimator A of P1.4)."),
        dict(d="Weighted conformal calibration, with weights from the classifier's density ratio, clipped at the 99th "
               "percentile (Barber et al. 2023).",
             why="Conformal coverage assumes that calibration and target data are exchangeable. Weighting the calibration blocks "
                 "towards the target mix restores the guarantee for the covariates the classifier sees.",
             risk="Extreme weights shrink the effective calibration sample.",
             check="Reference blocks resampled to the covariate mix of irrigated land and of valley land keep coverage between "
                   "0.88 and 0.92. Effective sample size at least 500 per stratum.",
             fail="The stratum is treated as outside the area of applicability."),
    ], depends="P1.2.")

    s += CARD("P1.4", "Excess, probability and coupling",
              "A calibrated probability that a block-month exceeds its rain-only expectation, the amount in millimetres with an "
              "interval, and a measure of whether the vegetation still follows rain.", [
        dict(d="Compute excess three ways. **A, bound:** water-year ET minus rain minus the root-zone water capacity; a lower "
               "bound that needs no model and no class. **B, model:** observed minus the rain-only median (P1.2). **C, local "
               "contrast,** at 30 m in the focus landscapes: the pivot interior against the non-irrigated corners of the same "
               "quarter section, and riparian strips against adjacent upland in the same gridMET cell.",
             why="Each fails differently (R9). A is biased by the ET and rain products, not by the model. B depends on the model "
                 "and its class input. C is free of rain error and model error, but compares different crops and exists only "
                 "where a clean neighbour does.",
             check="B at least as large as A in 95 % of irrigated field-years. Seasonal B within ±25 % of C on at least 200 "
                   "pivots in three landscapes, in a dry and in a wet year.",
             fail="Disagreement is reported by landscape and year. Volumes are given as the range from A to B."),
        dict(d="From each block's residual subtract the robust mean of the out-of-fold residuals of reference blocks within "
               "15 km (at least 20 blocks; 30 km if fewer).",
             why="gridMET gives one rain value per cell of about 4 km. A summer storm that the grid misses makes every block "
                 "under it look subsidised. That error is shared by neighbours, while irrigation ends at the field edge and "
                 "groundwater use follows the valley. Neighbouring rain-only blocks measure the shared part directly, and "
                 "because they are predicted by the same fold model, a local bias of that model cancels as well.",
             risk="Few reference blocks in dense irrigation districts. A regional signal, such as cooler and moister air near "
                  "irrigation, is removed with the error; this is accepted and stated.",
             check="(1) Variogram of reference residuals: share of variance with a range above 4 km reported; the removal must "
                   "lower the null SD. (2) The false-positive rate on held-out reference stays between 4 and 6 %. (3) The "
                   "correlation between reference residuals and the difference of gridMET rain from a second rain product "
                   "(Section 7) falls after the removal.",
             fail="Without enough neighbours the uncorrected residual is used and its interval is widened by the correlated "
                  "variance component."),
        dict(d="Turn the corrected residual into a probability with a two-groups model (Efron 2004). The density of all blocks "
               "in a stratum is a mixture of the null density, which is measured on held-out reference blocks, and a subsidised "
               "component. The probability replaces the fixed 95 % threshold.",
             why="With 2.0 million blocks and 216 months, a 5 % threshold flags millions of rain-only block-months. In winter, "
                 "when true subsidy is rare, nearly all flags would be false. The probability accounts for how common subsidy is "
                 "in each stratum and month, and the null is measured, not assumed (R10).",
             check="Mean probability on held-out reference blocks below 0.05. At least 80 % of planted block-months get a "
                   "probability above 0.5 when the planted amount is twice the detection limit.",
             fail="Strata in which the mixture cannot be separated are reported with the detection limit only."),
        dict(d="Joint index over the three targets: a signed Mahalanobis distance using the covariance of reference residuals, "
               "calibrated by its rank among calibration blocks.",
             why="OpenET models use Landsat reflectance and surface temperature, so ET and the indices share inputs and their "
                 "residuals are correlated. The covariance discounts the shared part. Rank calibration avoids assuming normal "
                 "residuals.",
             check="False-positive rate of the joint index on held-out reference between 4 and 6 %.",
             fail="The single targets are reported separately."),
        dict(d="Coupling. For each block and season, regress observed anomalies on expected anomalies over a moving five-year "
               "window. The slope b is the coupling to rain; 1 − b is the buffering coefficient.",
             why="The slope separates what the level of the residual cannot. A persistent site advantage and run-on scale with "
                 "rain (slope at or above 1). Water that rain did not supply fills the deficit, so function stays up when the "
                 "expectation falls (slope well below 1). A block that is always greener than expected and still follows rain "
                 "is a good site, not a groundwater user.",
             risk="Five years with little variation in rain give an unstable slope.",
             check="Reference blocks have slopes centred on 1 with a measured spread. A planted subsidy that fills half of each "
                   "deficit is recovered as a buffering coefficient of 0.5 ± 0.1. Windows in which the expected anomaly varies "
                   "less than a fixed minimum are marked undetermined.",
             fail="Use event contrasts (drought months against normal months) in place of the regression."),
        dict(d="Volume of excess = sum over block-months of probability × corrected residual × area. The form max(0, residual) "
               "of revision 2 is not used in sums.",
             why="Cutting noise at zero gives it a positive mean. Summed over millions of rain-only block-months, that alone "
                 "produces a large volume.",
             check="Zero-volume test (gate G1c): on held-out reference blocks the summed excess per class, decile and year lies "
                   "within ±1 mm per month of zero.",
             fail="The gate stops. The bias is traced to the stratum and the probability model is fitted again."),
        dict(d="Publish a detection-limit map: the smallest excess that reaches a probability of 0.5, per stratum and month.",
             why="No detected subsidy means nothing where the limit is high (few observations, poor transfer). Paper 2 needs to "
                 "know where the method cannot see (R12).",
             check="The limit computed from the null density is confirmed with planted effects at 1, 2 and 3 times the limit."),
    ], depends="P1.2, P1.3; a second rain product (Section 7).")

    s += CARD("P1.5", "Dynamic monthly states (the mask)",
              "A monthly state probability per block that later steps can use as a mask.", [
        dict(d="States: R rain-coupled, S subsidised, D deficit (below the rain-only expectation), B bare or harvested, W open "
               "water or flooded. Revision 2 called D stressed.",
             why="The rain-only expectation already contains the decline that drought causes. A block below it is doing worse "
                 "than rain explains (fallowed, burned, hailed, grazed short). That is a different thing from drought stress, "
                 "and the name should not suggest otherwise.",
             check="Definition; no test."),
        dict(d="Hidden Markov model per functional group with month-specific transitions. Emissions are the likelihood ratios of "
               "P1.4 and the water indices. A season is labelled buffered only if the coupling slope or an abrupt start and end "
               "supports it. Persistent excess with a slope at or above 1 is reported as R+ (favourable site).",
             why="An unmodelled site effect lasts for years. A state model fed with raw residuals would label such a block "
                 "subsidised for the whole record.",
             check="(1) Held-out reference blocks are in state S in at most 2 % of months. (2) Onset months agree within one "
                   "month with the change-point method for at least 70 % of irrigated block-seasons. (3) In Kansas, fields with "
                   "reported use above the median are in S in at least 70 % of peak-season months, and fields with zero reported "
                   "use in fewer than 10 %. (4) Dryland fields are not in S more often in wet years than in dry years.",
             fail="(1) or (4) mean that rain error or site effects leak into S: shorten the common-mode radius and tighten the R+ "
                  "rule. (3) is examined in P1.7."),
    ], depends="P1.4.")

    s += CARD("P1.6", "Source attribution: irrigation, natural groundwater or surface water",
              "Label the source of the excess where the label can be defended, and keep the label out of circular use later.", [
        dict(d="Two stages. First, irrigation against not irrigation, from field geometry, sharp edges, abrupt onset and "
               "registry labels. Then, for the rest, natural groundwater against surface water and run-on, from the coupling "
               "slope, flooding in the water indices, and terrain. Anything else is unresolved.",
             why="Irrigation has strong signatures and good labels. The natural classes are harder and have weaker labels, and "
                 "should not blur the first decision.",
             check="Balanced accuracy by class and by dry, normal and wet years, with spatially blocked folds; irrigation in "
                   "drought years at least 0.85 (H1.1). A run with shuffled labels gives 0.5."),
        dict(d="Labels for natural groundwater come from wells with a measured depth of 3 m or less within 500 m, and from NWI "
               "classes. The GDE map is a weak label only.",
             why="The GDE map is itself derived from satellite vegetation and moisture signals (Rohde et al. 2024b), so "
                 "agreement with it is partly built in.",
             check="Precision of the natural-groundwater class at least 0.7 against well-based labels.",
             fail="The class is reported as unresolved subsidy (gate G1 of revision 2)."),
        dict(d="Depth to water and HAND may be features here. In return, the source label is never used as an outcome against "
               "depth or terrain in Papers 2 and 3. Those papers use excess and coupling, which contain no depth information.",
             why="A label built from depth and then analysed against depth proves nothing.",
             check="An automated test in the Paper 2 code refuses the source label as an outcome."),
    ], depends="P1.5; labels (Section 7).")

    s += CARD("P1.7", "Validation and volumes",
              "Test the excess against measurements that were not used to make it.", [
        dict(d="Validate at flux towers with image chips. For each AmeriFlux site within 100 km of the aquifer, extract a 5 km "
               "chip of Landsat indices, OpenET and CDL with the existing notebooks, and run the model in transfer mode.",
             why=f"//Measured//: US-Ne1, US-Ne2 and US-Ne3 lie {kmN1}, {kmN2} and {kmN3} km outside the aquifer polygon. The "
                 "tiles keep only pixels inside the polygon, so the irrigated and rain-fed tower pair of revision 2 cannot be "
                 "read from them. The tower file named in the older catalogue is no longer in the data folder.",
             check="At the irrigated Mead sites the probability exceeds 0.5 in dry growing-season months and stays below 0.5 at "
                   "the rain-fed site. Excess against tower ET minus rain: slope between 0.7 and 1.3.",
             fail="Reported with the towers' position relative to the area of applicability; the eastern edge is wetter than "
                  "most of the aquifer."),
        dict(d="Compare yearly excess volume with reported pumping by place of use in Kansas, for the OpenET ensemble and for "
               "each single model.",
             why="The ratio of excess to pumping should be below 1, because part of the pumped water returns or is lost. "
                 "Proportionality over years within one field is a stricter test than agreement of totals. One OpenET model "
                 "(SIMS) is a crop-coefficient model meant for well-watered crops (Melton et al. 2022); results must not depend "
                 "on it.",
             check="Within-field correlation of yearly excess and pumping: median at least 0.5. Ratio between 0.5 and 1.0 for "
                   "the ensemble. Sign and ranking unchanged for each single model.",
             fail="Results per model are reported; the model that drives the disagreement is named and left out in a "
                  "sensitivity run."),
        dict(d="Compare dynamic irrigated area-months with LANID and AIM-HPA (Deines et al. 2019a) by year (hypothesis H1.3 of "
               "revision 2).",
             why="This measures what a monthly mask adds over an annual one.",
             check="Share of area-months on which the dynamic and static masks disagree, by year type, with intervals."),
    ], depends="P1.4 to P1.6; Kansas water-use records and tower data (Section 7).")

    s += CARD("P1.8", "Confirmation with independent sensors (new)",
              "Meet the condition of the original plan that persistence is shown independently in ET, photosynthesis and "
              "canopy water.", [
        dict(d="At the 1/24° scale, relate the subsidised area fraction of each cell and month to signals that share no input "
               "with Landsat or OpenET: TROPOMI and OCO SIF, SMAP vegetation optical depth, ECOSTRESS canopy minus air "
               "temperature, and MODIS NIRv.",
             why="OpenET ET and the Landsat indices come from the same scenes, so their agreement is partly built in. SIF and "
                 "microwave optical depth are measured by other instruments with other physics. MODIS has one sensor through "
                 "2012 and 2013, so it also tests the Landsat-7-only year and the sensor change.",
             risk="Coarse footprints. TROPOSIF in the folder ends in 2021 and the Caltech file in October 2019.",
             check="Within aridity terciles, in drought months, the rank correlation between subsidised fraction and each "
                   "independent anomaly is positive with an interval that excludes zero, and it rises from wet to dry months. "
                   "Cell-level anomalies of 2012 from Landsat and from MODIS correlate at 0.8 or more.",
             fail="Paper 1 states that the product is a Landsat-based measure and claims no independent confirmation for the "
                  "signal that failed."),
    ], depends="P1.5; completed Caltech TROPOMI file (Section 7).")
    return s
