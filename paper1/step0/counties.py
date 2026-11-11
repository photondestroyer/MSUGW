"""Census county polygons (cartographic boundary 1:500k, 2023) for the eight aquifer states.

Used for the P1.1 dryland rule (a crop is dryland-eligible where most of its pixels in the county lie
outside LANID's ever-irrigated area). Stored once in <shared>/counties/.
"""
import os

from .. import config as C
from ..core import paths, net


def run(ctx):
    """Stage `counties`."""
    import geopandas as gpd
    d = paths.shared_dir("counties")
    z, meta = net.download(C.COUNTIES_URL, f"{d}/{os.path.basename(C.COUNTIES_URL)}", log=ctx.log)
    gdf = gpd.read_file(f"zip://{z}")
    gdf = gdf[gdf["STATEFP"].isin(C.AQUIFER_STATE_FIPS)][["GEOID", "STATEFP", "NAME", "geometry"]]
    gdf = gdf.sort_values("GEOID").reset_index(drop=True)
    out = f"{ctx.out}/counties.gpkg"
    tmp = f"{ctx.scratch}/counties.gpkg"
    if os.path.exists(tmp):
        os.remove(tmp)
    gdf.to_file(tmp, driver="GPKG", layer="counties")
    from ..core.drive import deliver
    deliver(tmp, out)
    ctx.complete(dict(counties=int(len(gdf)), source=meta))


def load(stage_dir):
    import geopandas as gpd
    return gpd.read_file(f"{stage_dir}/counties.gpkg", layer="counties")
