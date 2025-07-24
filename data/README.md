# Input data processing notes: 07-24-25

## subset_era_5_to_icechunk.py
This script subsets the ERA5 Zarr store both spatially and temporally, selects a subset a variables and rechunks for a time-series friendly chunking scheme. To reduce the dask task graph size, it writes each variables at a time. This ran for 5 hours on a c8g.16xlarge VM and had a peak memory usage of 87.5 GB. 
This script has some tweaks for potential performance improvements. 
- Uses the Zarr-python v3 obstore backend. # At high concurrency, obstore is more performant then fsspec. 
- Sets Zarr's async concurrency at ~128. # Recs from Earthmover / DevSeed's performance testing of Zarr v3. 
- Use Icechunk to write. # Icechunks has a really performant writing IO + we can checkpoint at each variable.
- The script loops through variables one at a time. You could split each variable up into a separate batch job.



