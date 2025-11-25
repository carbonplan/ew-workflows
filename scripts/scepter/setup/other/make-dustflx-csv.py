# %%
# -------------------------------------------------------
# 
# Script to make a csv for dust flux data over time 
# 
# Note that as of 10/20/2024 it only works with the 
# SCEPTER/rock_buff_dust-ts_multiyear.py script
# 
# -------------------------------------------------------
import os
import numpy as np
import pandas as pd

from ew_workflows import batch_helperFxns as bhf

# %%
# --- set up the save params
# ***************************************************************
savehere = "s3://carbonplan-carbon-removal/ew-workflows-data/scepter/dust/"
savename = "gbas_100yr_3-yearly_001.csv"
# ***************************************************************

# %% 
# --- set up the function inputs
max_time = 100 # [years] the end of the batch simulation
application_years = list(range(0, max_time, 3))  # start at year zero
# application_years = [0,1,2,3,4,10,13,18,19,20,21,22,23,24,25,33,50]
dustsp = ['gbas'] * len(application_years)
dustrate = ['defer'] * len(application_years)
dustrad = ['defer'] * len(application_years)
dustsp_2nd = [None]
dustrate_2nd = [None]

# %% 
# --- run the function to make the time-dynamics csv
bhf.make_dustflx_csv(
        savehere=savehere,
        savename=savename,
        max_time=max_time,
        application_years=application_years,
        dustsp=dustsp,
        dustrate=dustrate,
        dustrad=dustrad,
        dustsp_2nd=dustsp_2nd,
        dustrate_2nd=dustrate_2nd,
        returndf=True,
        savedf=True,
    )   


# %%
