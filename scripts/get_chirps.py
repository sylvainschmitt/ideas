# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
filein = "results/limits/limits.shp"
filein2 = "data/capricho/capricho.shp"
# fileout = "results/data/chirps_guaviare.nc"

# libs
import geopandas as gp
import ee
import xarray as xr
import xesmf as xe
import numpy as np

# grid
area = gp.read_file(filein)
xmin = round(min(area.bounds.minx), 2)
xmax = round(max(area.bounds.maxx), 2)
ymin = round(min(area.bounds.miny), 2)
ymax = round(max(area.bounds.maxy), 2)
ds_out = xr.Dataset(
    {
        "lat": (["lat"], np.arange(ymin, ymax, 0.01), {"units": "degrees_north"}),
        "lon": (["lon"], np.arange(xmin, xmax, 0.01), {"units": "degrees_east"}),
    }
)
area = gp.read_file(filein2)
xmin = round(min(area.bounds.minx), 2)
xmax = round(max(area.bounds.maxx), 2)
ymin = round(min(area.bounds.miny), 2)
ymax = round(max(area.bounds.maxy), 2)
ds_out2 = xr.Dataset(
    {
        "lat": (["lat"], np.arange(ymin, ymax, 0.05), {"units": "degrees_north"}),
        "lon": (["lon"], np.arange(xmin, xmax, 0.05), {"units": "degrees_east"}),
    }
)

# get
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
ds = xr.open_dataset("UCSB-CHG/CHIRPS/PENTAD", engine='ee')
ds = ds.sel(time=slice("2000-01-01", "2022-12-31"))
ds = ds.transpose('time', 'lat', 'lon')
ds = ds.chunk(chunks = {'time':100, 'lat': 400, 'lon': 400})
# import matplotlib.pyplot as plt
# ds.sel(time="2001-01-01").precipitation.plot()
# plt.show()

# capricho
regridder = xe.Regridder(ds, ds_out2, "bilinear")
ds_r = regridder(ds, keep_attrs=True)
ds_r.to_netcdf("results/data/chirps_indices_capricho.nc")

# indices
ds = ds.groupby("time.year").sum()
regridder = xe.Regridder(ds, ds_out, "bilinear")
ds_r = regridder(ds, keep_attrs=True)
ds_r.to_netcdf("results/data/chirps_indices_amazon.nc")

# anomalies
ds_anom = ds.sel(year=slice(2020, 2022)).mean("year") -  ds.sel(year=slice(2000, 2002)).mean("year")
regridder = xe.Regridder(ds_anom, ds_out, "bilinear")
ds_r = regridder(ds_anom, keep_attrs=True)
ds_r.to_netcdf("results/data/chirps_anomalies_amazon.nc")
