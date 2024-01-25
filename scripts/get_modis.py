# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
# filein = "data/guaviare/guaviare.shp"
# fileout = "results/data/modis_guaviare.nc"

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
ds2 = ds2.groupby("time.year").mean() # monthly average because only january
ds2["LST_Day_1km"].values = ds2["LST_Day_1km"].values * 0.02 - 273.15 # scale factor + degrees

# indices
ds_tas = ds2.sel(year=slice(2015, 2020)).mean("year") - ds2.sel(year=slice(2001, 2006)).mean("year")
ds_tas = ds_tas.rename({"LST_Day_1km": "tas"})
ds_tasmin = ds2.sel(year=slice(2018, 2020)).min("year") - ds2.sel(year=slice(2001, 2003)).min("year")
ds_tasmin = ds_tasmin.rename({"LST_Day_1km": "tasmin"})
ds_tasmax = ds2.sel(year=slice(2018, 2020)).max("year") - ds2.sel(year=slice(2001, 2003)).max("year")
ds_tasmax = ds_tasmax.rename({"LST_Day_1km": "tasmax"})
ds_all = xr.merge([ds_tas, ds_tasmax, ds_tasmin])

# regid
regridder = xe.Regridder(ds_all, ds_out, "bilinear")
ds_r = regridder(ds_all, keep_attrs=True)
ds_r.to_netcdf(fileout)
