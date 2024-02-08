# libs

# import geopandas as gp
import xarray as xr
import rioxarray as rio

# lim = gp.read_file("results/limits/limits.shp")
# xmin = min(lim.bounds.minx)
# xmax = max(lim.bounds.maxx)
# ymin = min(lim.bounds.miny)
# ymax = max(lim.bounds.maxx)

ds = xr.open_mfdataset([
    "data/JRC_TMF_AnnualChange_v1_SAM_ID47_N10_W80/forobs/products/tmf_v1/AnnualChange/JRC_TMF_AnnualChange_v1_2001_SAM_ID47_N10_W80.tif",
    "data/JRC_TMF_AnnualChange_v1_SAM_ID48_N10_W70/forobs/products/tmf_v1/AnnualChange/JRC_TMF_AnnualChange_v1_2000_SAM_ID48_N10_W70.tif"]
                       )
ds.isel(band=0).rio.to_raster("test.tif")
ds2 = ds.to_dataarray
ds.to
# mask_lon = (ds.x >= xmin) & (ds.x <= xmax)
# mask_lat = (ds.y >= ymin) & (ds.y <= ymax)
# ds.where(mask_lon & mask_lat, drop = True).load()
    