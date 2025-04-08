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

# %%
# --- set up the save params
# ***************************************************************
savehere = "/home/tykukla/aglime-swap-cdr/scepter/dust-inputs"
savename = "cc_15yr_1app_no2nd_001.csv"
# ***************************************************************

# %% 
# --- DEFINE TIME STEPS
# [1] the total amount of time to split up into sub-runs
max_time = 15 # [years] the end of the batch simulation
# [2] list of the years where a new sub-run starts (values
#     cannot exceed max_time)
start_times = [0, 1]   # [years]

# --- DEFINE DUST APPLICATION
# (note, lists must have same length as start_times)
# [1] define the dust species applied at each timestep 
dustsp = ['cc', 'cc']
# [2] define the dust rates -- note these will override 
#     the default value from the default dict (or from 
#     the batch .csv) unless the entry is non-numeric (
#     (suggest to make those values 'defer' so it's clear
#      we're deferring to the default)
dustrate = ['defer', 0]
# [3] define dust radius 
#     as above, 'defer' or other non-numeric means the default
#     entry will be selected
dustrad = ['defer', 'defer']
# [4] define second dust species
dustsp_2nd = []   # leaving it empty means we use the default (though saying 'defer' should work too)
dustrate_2nd = [] # leaving it empty means we use the defualt (though saying 'defer' should work too)


# %% 
# --- calculate timestep durations 
timestep_dur = []
# Loop through the start_times and calculate the difference
for i in range(len(start_times)):
    if i < len(start_times) - 1:
        # calculate time difference to the next start_time
        interval = start_times[i+1] - start_times[i]
    else:
        # calculate time difference to max_time for the last start_time
        interval = max_time - start_times[i]
    
    # Append the interval to the list
    timestep_dur.append(interval)

# %% 
# --- BRING TOGETHER 
# [1] bring lists into a dictionary
list_dict = {
    "yr_start": start_times,
    "duration": timestep_dur,
    "dustsp": dustsp,
    "dustrate": dustrate,
    "dustrad": dustrad,
    "dustsp_2nd": dustsp_2nd,
    "dustrate_2nd": dustrate_2nd
}
# [2] remove empty lists from the dictionary 
filtered_data = {key: value for key, value in list_dict.items() if value}
# [3] create pd.Dataframe
df = pd.DataFrame(filtered_data)
df
# %%
# --- save result
df.to_csv(os.path.join(savehere, savename), index=False)

# %%
