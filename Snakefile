configfile: "config/config.yml"

# rules #
rule all:
   input:
      expand("results/data/{source}_{area}.nc",
             source=["chirps", "chelsa", "modis"],
             area=config["area"])

## data ##
include: "rules/get_chirps.py"
include: "rules/crop_chirps.py"
include: "rules/get_chelsa.py"
include: "rules/get_tmf.py"
include: "rules/get_modis.py"
