# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
# filein = "results/data/chelsa_capricho.nc"
# fileout = "results/data/modis_capricho.nc"

# libs
import ee
import xarray as xr
import xesmf as xe
  
# code
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
ds = xr.open_dataset("MODIS/061/MOD11A2", engine='ee')
ds = ds[["LST_Day_1km"]]
ds["LST_Day_1km"].values = ds["LST_Day_1km"].values * 0.02 - 273.15
base = xr.open_dataset(filein)
regridder = xe.Regridder(ds, base, "bilinear")
ds_r = regridder(ds, keep_attrs=True)
ds_r.to_netcdf(fileout)
