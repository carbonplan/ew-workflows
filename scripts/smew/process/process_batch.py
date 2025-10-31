# ----------------------------------------
# 
# process batch output from SMEW
# 
# ----------------------------------------
#%% 
import os 

import fsspec
import matplotlib.pyplot as plt
import numpy as np
import pickle
import s3fs
import xarray as xr

import process_batch_fxns as pbf

# ---- set up s3
s3_path = "s3://carbonplan-carbon-removal/SMEW/smew_output"
s3 = s3fs.S3FileSystem(anon=False)  # anon=True if public bucket


batchdir = "/home/tykukla/EWmodel/batch_inputs"
batchname = "site+tstep+mineral_v0.pkl"


# --- define coords
batch_coords = ['site', 'climfn', 'mineral']

# %% 
# --- read in batch dict 
batchdict = pbf.read_batch_dict(batchdir, batchname)

# %% 
# --- read in just one .zarr file (DEBUG)
# thisrun = batchdict[list(batchdict.keys())[1]]['rundir']
# zarrname = "results.zarr"
# store = s3fs.S3Map(
#     root=os.path.join(s3_path, thisrun, zarrname), 
#     s3=s3, check=False
# )
# ds = xr.open_zarr(store, consolidated=True)
# ds['M_rock_in'].values

# --- read in zarr 
zarrname = "results.zarr"

datasets = []    # [1, 5, 25, 45]
# for runkey in [list(batchdict.keys())[i] for i in [1, 5, 25, 45]]:
for runkey in list(batchdict.keys()):
    thisrun = batchdict[runkey]['rundir']
    store = s3fs.S3Map(
        root=os.path.join(s3_path, thisrun, zarrname), 
        s3=s3, check=False
    )
    tmpds = xr.open_zarr(store, consolidated=True)

    tmpds = tmpds.squeeze('runname', drop=True)

    # find the dims we need to add
    add_bool = [thiscoord not in list(tmpds.dims) for thiscoord in batch_coords]
    add_coords = [thiscoord for thiscoord, thisbool in zip(batch_coords, add_bool) if thisbool]
    # broadcast the dims we need to add
    for thiscoord in add_coords:
        tmpds = tmpds.expand_dims(climfile = [tmpds[thiscoord].values])

    # broadcast other dims for merging
    tmpds = tmpds.broadcast_like(tmpds[['site']])
    datasets.append(tmpds)


# --- Step 1: compute the finest spacing ---
def min_time_step(ds):
    times = ds.time.values
    # take differences, drop empty
    diffs = np.diff(np.sort(times))
    return diffs.min() if len(diffs) > 0 else None

steps = [min_time_step(ds) for ds in datasets if ds.sizes["time"] > 1]
finest_step = min([s for s in steps if s is not None])

# --- Step 2: build a regular grid at that spacing ---
tmin = min(ds.time.min().values for ds in datasets)
tmax = max(ds.time.max().values for ds in datasets)
regular_time = np.arange(tmin, tmax + finest_step, finest_step)

# --- Step 3: interpolate each dataset onto that grid ---
datasets_interp = [ds.interp(time=regular_time) for ds in datasets]

# %% 
# --- Step 4: merge or concat ---
merged = xr.merge(datasets_interp)



# %% 
# --- some fun plots
mysites = merged.site.values
myclims = merged.climfile.values
mineraldx = 1
thissite = mysites[0]
thisvar = 'M_rock'


tmpds = merged.sel(site=thissite)

for cf in myclims:
    print(cf)
    plt.plot(tmpds.time, tmpds[thisvar].sel(climfile=cf, app_t_per_ha=20.).isel(mineral=mineraldx, feedstock=mineraldx), label=cf)
    plt.legend()
    plt.title(f"site: {thissite}")





# %%
tmpds['temp_air'].sel(climfile = "yearly.nc").isel(mineral=0, app_t_per_ha=0, feedstock=0).values
# %%
