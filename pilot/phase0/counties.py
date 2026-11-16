"""Study counties (decision D02) and the replication-region counties.

Physical idea: the pilot's units of weather and drought are county-years ("deficit ... for every county and
year"), and dry episodes link neighbouring counties, so the county table carries borders as well.
Decisions: D02.   Deliverable: counties.gpkg, counties.csv, neighbours.json.
"""
import json

import numpy as np


def run(ctx):
    import geopandas as gpd
    from .. import config as C
    from .. import decisions as D
    from ..common import io, report
    d02 = D.value("D02")

    # ══════ STEP C-1  Census counties of Kansas (2023, 1:500k) ══════
    allc = io.counties().to_crs(5070)
    ks = allc[allc["STATEFP"] == d02["state_fips"]].copy()
    ctx.log.info("STEP C-1  %d Kansas counties", len(ks))

    # ══════ STEP C-2  overlap with GMD 1-5; a county is in when it overlaps by >= D02 min_overlap_km2 ══════
    gmd = io.vector("gmd").to_crs(5070)[["GMD_ID", "NAME", "geometry"]]
    ov = gpd.overlay(ks[["GEOID", "geometry"]], gmd, how="intersection", keep_geom_type=True)
    ov["km2"] = ov.geometry.area / 1e6
    by = ov.groupby("GEOID")["km2"].sum()
    main_gmd = ov.sort_values("km2").groupby("GEOID").last()["GMD_ID"]
    ks["gmd_overlap_km2"] = ks["GEOID"].map(by).fillna(0.0)
    ks["gmd"] = ks["GEOID"].map(main_gmd)
    ks["study"] = ks["gmd_overlap_km2"] >= d02["min_overlap_km2"]
    ctx.log.info("STEP C-2  %d study counties (overlap >= %.1f km2 with a GMD)", int(ks.study.sum()),
                 d02["min_overlap_km2"])

    # ══════ STEP C-3  centroids and areas ══════
    cen = ks.geometry.centroid
    ks["x5070"], ks["y5070"] = cen.x, cen.y
    ll = gpd.GeoSeries(cen, crs=5070).to_crs(4326)
    ks["lon"], ks["lat"] = ll.x, ll.y
    ks["area_km2"] = ks.geometry.area / 1e6
    ks["region"] = "kansas"

    # ══════ STEP C-4  neighbours: counties sharing a border (length > 100 m), among study counties ══════
    st = ks[ks.study].reset_index(drop=True)
    neigh = {g: [] for g in st.GEOID}
    sidx = st.sindex
    for i, (g, geom) in enumerate(zip(st.GEOID, st.geometry)):
        for j in sidx.query(geom.buffer(10), predicate="intersects"):
            if j != i and geom.boundary.intersection(st.geometry.iloc[j].boundary).length > 100:
                neigh[g].append(st.GEOID.iloc[j])
    ctx.log.info("STEP C-4  borders: median %d neighbours per county", int(np.median([len(v) for v in neigh.values()])))

    # ══════ STEP C-5  replication-region counties (metadata only): in the region box and over the aquifer ══════
    aq = io.vector("aquifer_shp").to_crs(5070).union_all()
    reps = []
    for name, r in C.REP_REGIONS.items():
        s = allc[allc.STATEFP == r["state"]].copy()
        c4326 = gpd.GeoSeries(s.geometry.centroid, crs=5070).to_crs(4326)
        inbox = (c4326.y <= r["lat_max"]) & (c4326.x >= r["lon_min"]) & (c4326.x <= r["lon_max"])
        s = s[inbox.to_numpy() & s.geometry.intersects(aq).to_numpy()].copy()
        s["region"], s["study"], s["gmd"], s["gmd_overlap_km2"] = name, False, np.nan, 0.0
        cen = s.geometry.centroid
        s["x5070"], s["y5070"] = cen.x, cen.y
        ll = gpd.GeoSeries(cen, crs=5070).to_crs(4326)
        s["lon"], s["lat"], s["area_km2"] = ll.x, ll.y, s.geometry.area / 1e6
        reps.append(s)
        ctx.log.info("STEP C-5  %s: %d counties", name, len(s))

    # ══════ STEP C-6  write ══════
    out = gpd.GeoDataFrame(__import__("pandas").concat([ks] + reps, ignore_index=True), crs=5070)
    out = out.rename(columns={"GEOID": "county", "NAME": "name"})
    tmp = f"{ctx.scratch}/counties.gpkg"
    import os
    if os.path.exists(tmp):
        os.remove(tmp)
    out.to_file(tmp, layer="counties", driver="GPKG")
    from ..core.drive import deliver
    deliver(tmp, f"{ctx.out}/counties.gpkg")
    report.table(ctx, out.drop(columns="geometry"), "counties", units={"area_km2": "km2", "x5070": "m"})
    report.text(ctx, "neighbours.json", json.dumps(neigh, indent=1))
    ctx.complete(dict(kansas=int(len(ks)), study=int(ks.study.sum()),
                      replication={k: int((out.region == k).sum()) for k in C.REP_REGIONS},
                      by_gmd={str(int(k)): int(v) for k, v in ks[ks.study].groupby("gmd").size().items()}))


def load(stage_dir, study_only=True):
    import geopandas as gpd
    g = gpd.read_file(f"{stage_dir}/counties.gpkg", layer="counties")
    g["county"] = g["county"].astype(str)
    return g[g.study] if study_only else g


def neighbours(stage_dir):
    with open(f"{stage_dir}/neighbours.json") as f:
        return {k: set(v) for k, v in json.load(f).items()}
