# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
filein = "data/guaviare/guaviare.shp"
fileout = "results/data/tmf_guaviare.nc"

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
ds = xr.open_dataset("projects/JRC/TMF/v1_2020/AnnualChanges", engine='ee')
ds_all = ds[["Dec2001"]]

# regid
regridder = xe.Regridder(ds_all, ds_out, "bilinear")
ds_r = regridder(ds_all, keep_attrs=True)
ds_r.to_netcdf(fileout)