# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
filein = "results/limits/limits.shp"

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

# get
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
ds = xr.open_dataset("MODIS/061/MOD11A2", engine='ee')

# prep
ds = ds.sel(time=ds.time.dt.month.isin([1])) # janury dryest month
ds = ds.where(ds.QC_Day.isin([int("00000000",2), int("00000100",2)])) 
# we excluded data where the estimated emissivity error was greater than 0.02 and where the LST error was greater than 1 K
    # Bits 0-1: Mandatory QA flags
    #     0: Pixel produced, good quality, not necessary to examine more detailed QA
    # Bits 2-3: Data quality flag
    #     0: Good data quality
    # Bits 4-5: Emissivity error flag
    #     0: Average emissivity error <= 0.01
    #     1: Average emissivity error <= 0.02
    # Bits 6-7: LST error flag
    #     0: Average LST error <= 1K
ds2 = ds[["LST_Day_1km"]]
ds2["LST_Day_1km"].values = ds2["LST_Day_1km"].values * 0.02 - 273.15 # scale factor + degrees

# indices
ds_tas = ds2.groupby("time.year").mean()
ds_tas = ds_tas.rename({"LST_Day_1km": "tas"})
ds_tas = ds_tas.transpose('year', 'lat', 'lon')
# ds_tasmin = ds2.groupby("time.year").min("time")
# ds_tasmin = ds_tasmin.rename({"LST_Day_1km": "tasmin"})
# ds_tasmax = ds2.groupby("time.year").max("time")
# ds_tasmax = ds_tasmax.rename({"LST_Day_1km": "tasmax"})
# ds_all = xr.merge([ds_tas, ds_tasmax, ds_tasmin])
# ds_all = ds_all.transpose('year', 'lat', 'lon')
# import matplotlib.pyplot as plt
# ds_all.sel(year=2001).tas.plot()
# plt.show()
regridder = xe.Regridder(ds_tas, ds_out, "bilinear")
ds_r = regridder(ds_tas, keep_attrs=True)
ds_r.to_netcdf("results/data/modis_indices_amazon.nc")

# anomalies
ds_anom = ds_tas.sel(year=slice(2020, 2022)).mean("year") -  ds_tas.sel(year=slice(2001, 2003)).mean("year")
regridder = xe.Regridder(ds_anom, ds_out, "bilinear")
ds_r = regridder(ds_anom, keep_attrs=True)
ds_r.to_netcdf("results/data/modis_anomalies_amazon.nc")
