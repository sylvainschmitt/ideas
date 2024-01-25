# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
filein = "data/guaviare/guaviare.shp"
fileout = "results/data/hansen_guaviare.nc"

# libs
import geopandas as gp
import ee
import xarray as xr
import xesmf as xe
import numpy as np

# grid
area = gp.read_file(filein)
xmin = round(area.bounds.minx[0], 2)
xmax = round(area.bounds.maxx[0], 2)
ymin = round(area.bounds.miny[0], 2)
ymax = round(area.bounds.maxy[0], 2)
ds_out = xr.Dataset(
    {
        "lat": (["lat"], np.arange(ymin, ymax, 0.01), {"units": "degrees_north"}),
        "lon": (["lon"], np.arange(xmin, xmax, 0.01), {"units": "degrees_east"}),
    }
)

# get
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
ee.Image("UMD/hansen/global_forest_change_2022_v1_10")
i = ee.ImageCollection(ee.Image("UMD/hansen/global_forest_change_2022_v1_10"))
ds = xr.open_dataset(i, engine='ee')
ds = ds[["loss"]]

# indices
regridder = xe.Regridder(ds, ds_out, "bilinear")
ds_r = regridder(ds, keep_attrs=True)
ds_r.to_netcdf(fileout)
