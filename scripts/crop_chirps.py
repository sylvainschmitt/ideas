# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
area_file = snakemake.input[0]
chirps_raw_file = snakemake.input[1]
chirps_file = snakemake.output[0]

# test
# area_file = "data/capricho/capricho.shp"
# chirps_raw_file = "results/data/chirps.nc"
# chirps_file = "results/data/chirps_capricho.nc"

# libs
import geopandas as gp
import xarray as xr

# code
area = gp.read_file(area_file)
ds = xr.open_dataset(chirps_raw_file)
ds = ds.rio.write_crs("epsg:4326")
ds = ds.rio.clip_box(minx=area.bounds.minx[0],
                     miny=area.bounds.miny[0],
                     maxx=area.bounds.maxx[0],
                     maxy=area.bounds.maxy[0])
ds.to_netcdf(chirps_file)
