rule get_chirps:
    output:
        # temp("results/data/chirps.nc")
        "results/data/chirps.nc"
    log:
        "results/logs/get_chirps.log"
    benchmark:
        "results/benchmarks/get_chirps.benchmark.txt"
    shell:
      "wget 'https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_monthly/netcdf/chirps-v2.0.monthly.nc' -o {output}"
