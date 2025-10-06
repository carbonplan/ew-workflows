# -------------------------------------------
# 
# Script to replace the "run-*.sh" scripts
# Reads in a parameter file, a batch file,
# and generates argo submit commands with 
# some bleed-in delay 
# 
# -------------------------------------------
# 
# Debug note -- this script relies on relative paths...
# if working in an interactive session, make sure the 
# working directory is the file directory (os.getcwd())
# and, if it's not, change it (os.chdir("/path/to/file/dir"))
# 
# %% 
import os
import sys 

# --- read in argo helper functions  
sys.path.append(os.path.abspath('/home/tykukla/ew-workflows/run_scepter'))
import argo_helper_fxns as ahf
# ---


# %% 
# ****************************************************************
run_pars = "batch-base.yaml"
multiyear = False

# run_pars = "batch-base-multiyear.yaml"
# multiyear = True

# run_pars = "tuneup_all.yaml"
# multiyear = False
# ****************************************************************
if multiyear:   # decide whether to skip the model duration check when determining reruns
    skip_duration_check = True
else:
    skip_duration_check = False

# %% 
ahf.run_multiple(parameter_yaml = run_pars)

# %%
ahf.retry_failed_runs_iteratively(
    max_reruns = 5,
    max_delays = 30, 
    rerun_delay = 20,
    parameter_yaml = run_pars,
    multiyear = multiyear,
    maindir = "/home/tykukla/ew-workflows",
    parameter_yaml_subdir = "inputs/scepter/params",
    completed_fn = "completed.res",
    check_results_fn  = "check_results.res",
    check_logs_fn  = "check_logs.res",
    duration_threshold_frac = 0.2,
    workflow_name_runmultiple = "scepter-pyworkflow.yaml",
    bleed_delay_runmultiple = 15,
    stale_threshold_minutes = 80,
    skip_initial_delay=False,
    skip_duration_check=skip_duration_check,
)

# %%

# import os 
# from pathlib import Path
# import re
# import shutil
# import subprocess
# import time 

# import numpy as np
# import pandas as pd
# import yaml

# parameter_yaml = run_pars

# parameter_yaml_subdir="inputs/scepter/params"
# maindir="/home/tykukla/ew-workflows"
# workflow_name = "scepter-pyworkflow.yaml"
# bleed_delay=15
# echo_command=True
# norun_debug=False
# rerun_on=False
# df_reruns=None

# # *****************************************************
# # if it's a rerun, check that we have the df we need
# if rerun_on and (df_reruns is None):
#     raise ValueError("run_multiple fxn expects a rerun dataframe because rerun_on is set to True, but df_reruns is None.")
# # *****************************************************

# # --- read in the parameter file
# # create parameter file path
# parameterfile = os.path.join(maindir, parameter_yaml_subdir, parameter_yaml)
# # check system arguments, or set default
# with open(parameterfile, "r") as file:
#     pars = yaml.safe_load(file)

# # --- read in the batch.csv to get the number of rows
# if rerun_on:
#     df_batch = df_reruns.copy()
# else:
#     df_batch = pd.read_csv(os.path.join(pars['batch-input-dir'], pars['batch-input']))

# # --- loop through rows
# for idx in range(len(df_batch)):
#     if rerun_on:  # then select a custom run index
#         idx_run = int(df_batch.index[idx] + 1)  # add 1 so lowest idx is 1
#     else:
#         # lowest idx is one
#         idx_run = idx + 1
#     # generate the run command
#     command = f'argo submit {workflow_name} --parameter-file {parameterfile} -p batch-index={str(idx_run)}'

#     # run command 
#     if echo_command or norun_debug:
#         print(command)
    
#     if not norun_debug:
#         subprocess.run(command, shell=True)

#     # delay until the next argo submit command
#     time.sleep(bleed_delay)
# # %%
