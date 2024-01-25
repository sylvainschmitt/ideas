rule get_modis:
    input:
        "data/{area}/{area}.shp"
    output:
        "results/data/modis_{area}.nc"
    log:
        "results/logs/get_modis_{area}.log"
    benchmark:
        "results/benchmarks/get_modis_{area}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    script:
      "../scripts/get_modis.py"
