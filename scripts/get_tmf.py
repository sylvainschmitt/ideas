# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
filein = "data/capricho/capricho.shp"
fileout = "results/data/tmf_capricho.nc"

# libs
import geopandas as gp
import ee
import xarray as xr
import rioxarray as rio
  
# code with ee
# import ee
# import xarray as xr
# ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
# ds = xr.open_dataset("projects/JRC/TMF/v1_2022/AnnualChanges", engine='ee').rename({"lon" : "x", "lat" : "y"})
# ic = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY").filterDate('1992-10-05', '1992-11-05')
# leg1 = ee.Geometry.Rectangle(113.33, -43.63, 153.56, -10.66)
# ds = xr.open_dataset(
#     ic,
#     engine='ee',
#     projection=ic.first().select(0).projection(),
#     geometry=leg1
# )
# ds.to_netcdf("test.nc")