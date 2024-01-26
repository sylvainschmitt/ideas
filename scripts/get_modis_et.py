# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
# filein = "data/guaviare/guaviare.shp"
# fileout = "results/data/modis_et_guaviare.nc"

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
ds = xr.open_dataset("MODIS/061/MOD16A2GF", engine='ee')

# prep
ds = ds.sel(time=ds.time.dt.month.isin([1])) # janury dryest month
ds = ds.where(ds.ET_QC.isin([int("0000000",2), int("0100000",2)])) 
# we kept pixels with (1) only good quality (MODLAND), 
# (2) on both Terra and Aqua sensors, 
# (3) detectors fine for up to 50% of channels 1, 2,
# (4) without significant clouds,
# (5) and the main method used for best result possible (no saturation)
    # Bit 0: MODLAND_QC bits
    #     0: Good quality (main algorithm with or without saturation)
    # Bit 1: Sensor
    #     0: Terra
    #     1: Aqua
    # Bit 2: Dead detector
    #     0: Detectors fine for up to 50% of channels 1, 2
    # Bits 3-4: Cloud state
    #     0: Significant clouds NOT present (clear)
    # Bits 5-7: SCF_QC (five level confidence score)
    #     0: Main method used, best result possible (no saturation)

ds2 = ds[["ET"]]
ds2["ET"].values = ds2["ET"].values * 0.1 # scale factor

# indices
ds_et = ds2.sel(time=slice("2015-01-01", "2020-12-31")).groupby("time.year").mean().mean("year") \
    - ds2.sel(time=slice("2001-01-01", "2006-12-31")).groupby("time.year").mean().mean("year")
ds_et = ds_et.rename({"ET": "et"})

ds_etmin = ds2.sel(time=slice("2015-01-01", "2020-12-31")).min("time") - \
    ds2.sel(time=slice("2001-01-01", "2006-12-31")).min("time")
ds_etmin = ds_etmin.rename({"ET": "etmin"})

ds_etmax = ds2.sel(time=slice("2015-01-01", "2020-12-31")).max("time") - \
    ds2.sel(time=slice("2001-01-01", "2006-12-31")).max("time")
ds_etmax = ds_etmax.rename({"ET": "etmax"})

ds_all = xr.merge([ds_et, ds_etmin, ds_etmax])

# regid
regridder = xe.Regridder(ds_all, ds_out, "bilinear")
ds_r = regridder(ds_all, keep_attrs=True)
ds_r.to_netcdf(fileout)
