import os; os.environ["GDAL_CACHEMAX"]="256"
import numpy as np, rasterio, geopandas as gpd, json
from rasterio.features import geometry_mask
from rasterio.windows import Window
from scipy import ndimage
poly=gpd.read_file("G:/MSU_GWB/datasets/HPA_polygon/hp_bound2010.shp").to_crs(4326)
geoms=list(poly.geometry)
f="G:/MSU_GWB/datasets/merged_datasets/GLO_DEM_ogallala.tif"
tot=dict(zero_in=0, zero_in_interior=0, nonzero_out=0, in_px=0, zero_out=0)
ex=[]
with rasterio.open(f) as r:
    step=256
    for r0 in range(0, r.height, step):
        w=Window(0, r0, r.width, min(step, r.height-r0))
        a=r.read(1, window=w)
        inside=~geometry_mask(geoms, out_shape=a.shape, transform=r.window_transform(w), all_touched=False)
        z=a==0
        tot["in_px"]+=int(inside.sum()); tot["zero_in"]+=int((z&inside).sum())
        tot["zero_out"]+=int((z&~inside).sum()); tot["nonzero_out"]+=int((~z&~inside).sum())
        # interior = inside and at least 3 px from the polygon edge (erode the inside mask)
        core=ndimage.binary_erosion(inside, iterations=3)
        zi=z&core; n=int(zi.sum()); tot["zero_in_interior"]+=n
        if n and len(ex)<10:
            rr,cc=np.nonzero(zi); xs,ys=rasterio.transform.xy(r.window_transform(w), rr[:1], cc[:1]); ex.append((float(xs[0]),float(ys[0]),n))
print(json.dumps(tot), ex)
