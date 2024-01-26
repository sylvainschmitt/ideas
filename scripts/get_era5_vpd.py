# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
filein = "data/guaviare/guaviare.shp"
fileout = "results/data/era_vpd_guaviare.nc"

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
ds = xr.open_dataset("ECMWF/ERA5/MONTHLY", engine='ee')

# prep
ds = ds.sel(time=ds.time.dt.month.isin([1])) # janury dryest month
def get_e(tas, sp):
        a = 611.21
        b = 18.678 - (tas / 234.5)
        c = 257.14
        f = 1.00072 + pow(10,-7) * sp * (0.032 + 5.9 * pow(10,-6) * pow(tas,2))
        return(f * a * (np.exp(b * tas / (c + tas))))
def get_vpd(tas, das, sp):
        e = get_e(das, sp)
        esat = get_e(tas, 101325)
        return((esat - e) / 1000)
ds["total_precipitation"].values = get_vpd(ds["mean_2m_air_temperature"].values - 273.15,
                                           ds["dewpoint_2m_temperature"].values - 273.15,
                                           ds["surface_pressure"].values)
ds2 = ds[["total_precipitation"]]
ds2 = ds2.rename({"total_precipitation": "vpd"})
ds_vpd = ds2.sel(time=slice("2015-01-01", "2020-12-31")).groupby("time.year").mean().mean("year") \
    - ds2.sel(time=slice("2001-01-01", "2006-12-31")).groupby("time.year").mean().mean("year")

# regid
regridder = xe.Regridder(ds_vpd, ds_out, "bilinear")
ds_r = regridder(ds_vpd, keep_attrs=True)
ds_r.to_netcdf(fileout)
