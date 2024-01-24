configfile: "config/config.yml"

# rules #
rule all:
   input:
      expand("results/data/chirps_{area}.nc",
             area=config["area"])

## data ##
include: "rules/get_chirps.py"
