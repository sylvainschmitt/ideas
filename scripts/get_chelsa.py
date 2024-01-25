# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
area_file = snakemake.input[0]
chelsa_file = snakemake.output[0]

# test
# area_file = "data/capricho/capricho.shp"
# chelsa_file = "results/data/chelsa_capricho.nc"

# libs
import geopandas as gp
import xarray as xr
import rioxarray as rio
import pandas as pd
import datetime
  
# code
area = gp.read_file(area_file)
years = list(range(1980, 2019))
years.remove(2013)
years.remove(2016)
ds_all = []
for year in years:
        a = []
        for month in range(1, 13):
                url = 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/' + "monthly" + '/' + "pr" + "/" \
                        + 'CHELSA_' + "pr" + '_' + '%02d' % (month,) + '_' + str(year) + "_V.2.1.tif"
                ds = rio.open_rasterio(url, decode_coords="all").to_dataset('band').rename_vars({1 : "pr"})
                a.append(ds.rio.clip_box(minx=area.bounds.minx[0], 
                                 miny=area.bounds.miny[0],
                                 maxx=area.bounds.maxx[0], 
                                 maxy=area.bounds.maxy[0]))
        ds_year = xr.concat(a, pd.Index(pd.date_range(datetime.datetime(year, 1, 1), periods=12, freq="M"), name="time"))
        del a
        ds_all.append(ds_year)
        del ds_year
ds = xr.merge(ds_all)
ds = ds[["time", "x", "y", "pr"]]
ds["pr"].values = ds["pr"].values * 0.01
ds.pr.attrs = {'standard_name': 'precipitation', 'long_name': 'Monthly precipitation', 'units': 'mm month-1', 'explanation' : 'Precipitation" in the earth\'s atmosphere means precipitation of water in all phases.'}
ds.chunk({'time':1, 'x':200, 'y':200}).to_netcdf(chelsa_file)
