# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
filein = "results/limits/limits.shp"
# fileout = "results/data/chirps_guaviare.nc"

# libs
import geopandas as gp
import ee
import xarray as xr
import xesmf as xe
import numpy as np
import rioxarray as rio

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

# get
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
ds = xr.open_dataset("UCSB-CHG/CHIRPS/DAILY", engine='ee')
ds = ds.sel(time=slice("2001-01-01", "2020-12-31")).resample(time="1MS").sum(dim="time")
ds = ds.transpose('time', 'lat', 'lon')

# import matplotlib.pyplot as plt
# ds.sel(time="2001-01-01").precipitation.plot()
# plt.show()

# regid
regridder = xe.Regridder(ds, ds_out, "bilinear")
ds_r = regridder(ds, keep_attrs=True)
# ds_r.to_netcdf(fileout2)
ds_r.to_netcdf("results/data/chirps_month.nc")

# indices
ds_jan = ds.sel(time=ds.time.dt.month.isin([1])).groupby("time.year").sum()
ds_jan = ds_jan.rename({"precipitation": "pr_jan"})
ds_jun = ds.sel(time=ds.time.dt.month.isin([6])).groupby("time.year").sum()
ds_jun = ds_jun.rename({"precipitation": "pr_jun"})
ds_year = ds.groupby("time.year").sum()
ds_year = ds_year.rename({"precipitation": "pr"})
ds_all = xr.merge([ds_jan, ds_jun, ds_year])

# regid
regridder = xe.Regridder(ds_all, ds_out, "bilinear")
ds_r = regridder(ds_all, keep_attrs=True)
# ds_r.to_netcdf(fileout2)
ds_r.to_netcdf("results/data/chirps_indices.nc")

# anomalies
ds_anom = ds_all.sel(year=slice(2018, 2020)).mean("year") -  ds_all.sel(year=slice(2001, 2003)).mean("year")
                                                   
# regid
regridder = xe.Regridder(ds_anom, ds_out, "bilinear")
ds_r = regridder(ds_anom, keep_attrs=True)
# ds_r.to_netcdf(fileout2)
ds_r.to_netcdf("results/data/chirps_anomalies.nc")