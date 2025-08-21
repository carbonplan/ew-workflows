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
run_pars = "batch-meanAnnliming_fert_fixedRate-base.yaml"
multiyear = False

# run_pars = "batch-meanAnnliming_fert_dust-multiyear-base.yaml"
# multiyear = True
# ****************************************************************
if multiyear:   # decide whether to skip the model duration check when determining reruns
    skip_duration_check = True
else:
    skip_duration_check = False


# %%
ahf.retry_failed_runs_once(
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
    skip_duration_check=skip_duration_check,
)