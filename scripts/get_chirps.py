# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
filein = "data/guaviare/guaviare.shp"
fileout = "results/data/chirps_guaviare.nc"

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
ds = xr.open_dataset("UCSB-CHG/CHIRPS/DAILY", engine='ee')

ds_old = ds.sel(time=slice("2001-01-01", "2006-12-31")).resample(time="1MS").sum(dim="time").groupby("time.month").mean("time")
ds_new = ds.sel(time=slice("2015-01-01", "2020-12-31")).resample(time="1MS").sum(dim="time").groupby("time.month").mean("time")
ds_diff = ds_new - ds_old

# indices
regridder = xe.Regridder(ds_diff, ds_out, "bilinear")
ds_r = regridder(ds_diff, keep_attrs=True)
ds_r.to_netcdf(fileout)
