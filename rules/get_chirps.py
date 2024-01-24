rule get_chirps:
    input:
        "data/{area}.shp",
    output:
        "results/data/chirps_{area}.nc"
    log:
        "results/logs/get_chirps_{area}.log"
    benchmark:
        "results/benchmarks/get_chirps_{area}.benchmark.txt"
    threads: 1
    resources:
        mem_mb=1000
    conda:
        "../envs/xarray.yml"
    script:
      "../scripts/get_chirps.py"
