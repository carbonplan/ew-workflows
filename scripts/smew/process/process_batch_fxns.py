# -------------------------------
# 
# process batch helper fxns
# 
# -------------------------------
# %% 
import os 

import fsspec
import numpy as np
import pandas as pd
import pickle
import s3fs
import xarray as xr

# %% 
def read_batch_dict(
    batchdir: str,
    batchname: str,
) -> dict:
    '''
    Read the batch dictionary and add the run dir to each run

    Parameters 
    ----------
    batchdir : str
        path/to/batch_inputs
    batchname : str
        name of batch .pkl file
    '''
    # --- read in the batch file 
    with open(os.path.join(batchdir, batchname), "rb") as f:
        batchdict = pickle.load(f)

    # --- add rundirs to batch
    batchname_noSuff = batchname.replace(".pkl", "")
    for key in list(batchdict.keys()):
        rundir_tmp = f'{batchname_noSuff}--{key}'
        batchdict[key]['rundir'] = rundir_tmp

    return batchdict 
