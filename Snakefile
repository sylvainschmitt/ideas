configfile: "config/config.yml"

# rules #
rule all:
   input:
      expand("results/area/{area}",
             area=config["regions"])
      # expand("results/data/{source}_{area}.nc",
      #        source=["hansen", "modis"],
      #        area=config["area"])

## rules ##
include: "rules/get_gadm.py"
include: "rules/get_chirps.py"
include: "rules/crop_chirps.py"
include: "rules/get_chelsa.py"
include: "rules/get_tmf.py"
include: "rules/get_hansen.py"
include: "rules/get_modis.py"
