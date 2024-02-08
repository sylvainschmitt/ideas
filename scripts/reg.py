# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
filein = snakemake.input[0]
fileout = snakemake.output[0]

# test
filein = "results/data/tmf_deforested_2020.tif"

# libs
import numpy as np
import gstools as gs
import xarray as xr
import rioxarray as rio

pop = rio.open_rasterio(filein)
w_surface_sp = weights.Queen.from_xarray(pop)
pop

# grid
ds = xr.open_dataset(filein)
w_surface_sp = weights.Queen.from_xarray(pop)

data2 = ds["band_data"].to_numpy()
pos = data2.T[:2]  # lat, lon
field = data2.T[2]  # temperature

model = gs.Exponential(dim=2, var=2, len_scale=8)
srf = gs.SRF(model, mean=0, seed=19970221)
bins = np.arange(500)
bin_center, gamma = gs.vario_estimate(pos, field, bins)
fit_model = gs.Stable(dim=2)
fit_model.fit_variogram(bin_center, gamma, nugget=False)
ax = fit_model.plot(x_max=500)
