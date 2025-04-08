# ---------------------------------------
# 
# Script with helper functions for argo 
# batch simulations 
# 
# ---------------------------------------


# 
# Debug note -- these fxns rely on relative paths...
# if working in an interactive session, make sure the 
# working directory is the file directory (os.getcwd())
# and, if it's not, change it (os.chdir("/path/to/file/dir"))
# 


# %% 
import os 
from pathlib import Path
import re
import shutil
import subprocess
import time 

import numpy as np
import pandas as pd
import yaml


# %% 
def run_multiple(
    parameter_yaml: str,
    parameter_yaml_subdir: str="inputs/scepter/params",
    maindir: str="/home/tykukla/aglime-swap-cdr",
    workflow_name: str="scepter-pyworkflow.yaml",
    bleed_delay: int=15,
    echo_command: bool=True,
    norun_debug: bool=False,
    rerun_on: bool=False,
    df_reruns: pd.DataFrame=None,
):
    '''
    Generate the argo submit command for every row in the batch
    csv file defined by the parameter_yaml. 

    Parameters
    ----------
    parameter_yaml : str
        name of the parameter file for this batch. Something like
        "batch_pars.yaml". 
    parameter_yaml_subdir : str
        location of the parameter_yaml file
    maindir : str
        location of the inputs directory (usually '/my/path/to/aglime-swap-cdr')
    workflow_name : str
        name of the argo workflow .yaml file 
    bleed_delay : int
        [seconds] delay between each argo submit command
    echo_command : bool
        True means the argo submit command is printed before it's submitted.
        False means it's submitted without echoing. 
    norun_debug : bool
        True means we don't run the subprocess, we just print the command
        False means we do run the subprocess. 
    rerun_on : bool
        True means that this is a rerun so we don't use the batch csv file
        from the Parameter input, but instead expect one from the function.
        False means we use the standard batch csv to get df_batch
    df_reruns : pd.DataFrame
        None by default. This is only used (and indeed required) 
        if rerun_on==True. 

    Returns
    -------
    ''' 
    # *****************************************************
    # if it's a rerun, check that we have the df we need
    if rerun_on and (df_reruns is None):
        raise ValueError("run_multiple fxn expects a rerun dataframe because rerun_on is set to True, but df_reruns is None.")
    # *****************************************************

    # --- read in the parameter file
    # create parameter file path
    parameterfile = os.path.join(maindir, parameter_yaml_subdir, parameter_yaml)
    # check system arguments, or set default
    with open(parameterfile, "r") as file:
        pars = yaml.safe_load(file)
    
    # --- read in the batch.csv to get the number of rows
    if rerun_on:
        df_batch = df_reruns.copy()
    else:
        df_batch = pd.read_csv(os.path.join(pars['batch-input-dir'], pars['batch-input']))

    # --- loop through rows
    for idx in range(len(df_batch)):
        if rerun_on:  # then select a custom run index
            idx_run = int(df_batch.index[idx] + 1)  # add 1 so lowest idx is 1
        else:
            # lowest idx is one
            idx_run = idx + 1
        # generate the run command
        command = f'argo submit {workflow_name} --parameter-file {parameterfile} -p batch-index={str(idx_run)}'

        # run command 
        if echo_command or norun_debug:
            print(command)
        
        if not norun_debug:
            subprocess.run(command, shell=True)

        # delay until the next argo submit command
        time.sleep(bleed_delay)


# (check when last file was updated to decide if model stopped running)
def is_last_update_older_than(
    dir_path: str, 
    minutes: int = 20
    ) -> bool:
    """
    Returns True if the most recent update in dir_path is older than `minutes` minutes.
    """
    dir_path = Path(dir_path)  # Convert to a Path object
    files = list(dir_path.glob("*"))  # Get all files

    if not files:  # If directory is empty
        return True

    # Get the most recent modification time of any file in the directory
    latest_mtime = max(file.stat().st_mtime for file in files if file.is_file())

    # Compare with current time
    return (time.time() - latest_mtime) > (minutes * 60)



def check_results_read(
    check_file: str,
    file_on_s3: bool=False,
) -> dict:
    '''
    Read the "check_results.res" file and return the result

    Parameters
    ----------
    check_file : str
        Location of the "check_results.res" file 
        (example: '/path/to/my/scepter_output/runname/check_results.res')
    file_on_s3 : bool
        Whether the file is in an s3 bucket
    
    Returns 
    -------
    dict 
        Each check and the result as a dictionary
    '''
    check_results = {}
    # Read the file and process each line
    if file_on_s3:
        import s3fs
        fs = s3fs.S3FileSystem()

        with fs.open(check_file, "r") as file:
            for line in file:
                # Split the line at the first colon to separate check name and result
                if "--" in line:
                    check_name, result = line.split("--", 1)
                    result = result.strip()  # Remove extra spaces and newline characters
                    result_bool = result.split('\t')[1]
                    check_results[check_name.strip()] = result_bool == "True"  # Store True/False as boolean


    else:
        with open(check_file, "r") as file:
            for line in file:
                # Split the line at the first colon to separate check name and result
                if "--" in line:
                    check_name, result = line.split("--", 1)
                    result = result.strip()  # Remove extra spaces and newline characters
                    result_bool = result.split('\t')[1]
                    check_results[check_name.strip()] = result_bool == "True"  # Store True/False as boolean


    # Print the check results
    # for check, result in check_results.items():
    #     print(f"{check}: {result}")

    return check_results


def check_logs_read(
    check_file: str,
    file_on_s3: bool=False,
) -> dict:
    '''
    Read the "check_logs.res" file and return the model and target durations

    Parameters
    ----------
    check_file : str
        Location of the "check_results.res" file 
        (example: '/path/to/my/scepter_output/runname/check_results.res')
    file_on_s3 : bool
        Whether the file is in an s3 bucket
    
    Returns 
    -------
    dict 
        elements for "target" and "model" with each value being a duration in yrs
    '''
    # Initialize variables for target and model numbers
    target = None
    model = None

    # Read the file and search for target and model numbers
    if file_on_s3:
        import s3fs
        fs = s3fs.S3FileSystem()

        with fs.open(check_file, "r") as file:
            for line in file:
                # Use regular expressions to find the target and model values
                target_match = re.search(r"target:\s*(\d+\.?\d*)", line)
                model_match = re.search(r"model:\s*(\d+\.?\d*)", line)
                
                if target_match:
                    target = float(target_match.group(1))  # Convert target to float
                if model_match:
                    model = float(model_match.group(1))  # Convert model to float
    
    else:
        with open(check_file, "r") as file:
            for line in file:
                # Use regular expressions to find the target and model values
                target_match = re.search(r"target:\s*(\d+\.?\d*)", line)
                model_match = re.search(r"model:\s*(\d+\.?\d*)", line)
                
                if target_match:
                    target = float(target_match.group(1))  # Convert target to float
                if model_match:
                    model = float(model_match.group(1))  # Convert model to float

    # combine as dict
    outdict = {
        'target': target,
        'model': model,
    }

    return outdict



def checkrow_for_rerun(
    row: pd.Series,
    pars: dict,
    duration_threshold_frac: float = 0.2,
    completed_fn: str = "completed.res", 
    check_results_fn: str = "check_results.res",
    check_logs_fn: str = "check_logs.res",
    stale_threshold_minutes: float = 15,
) -> bool:
    '''
    Check the results for a given row in the batch dataFrame to see
    if we need to rerun it. 

    We first check if the output exists on AWS, which only works if 
    the AWS bucket is defined in batch*.csv. If it's not found there,
    we search locally. At each step, we can opt to rerun the case if 
    one of the checks fails: 
        1. The completed.res file is absent from the output directory
        2. One of the checks in the check_results.res file failed (e.g., is False)
        3. The model duration exceeds the allowable threshold difference from the target duration 
           (based on the check_logs.res file)
    
    Parameters
    ----------
    row : pd.Series
        A single row from the batch*.csv file representing a single simulation
    pars : dict
        Dictionary of parameters from the parameters directory (.yaml file) 
        used for finding the model directory
    duration_threshold_frac : float
        The absolute fractional difference allowed between the target duration of 
        the simulation and the actual duration. 
    completed_fn : str
        name of the completed.res file 
    check_results_fn : str
        name of the check_results file 
    check_logs_fn : str
        name of the check_logs file 
    stale_threshold_minutes : float
        [minutes] since the last file update in the localdir/flx directory before 
        we deem the run stale
    
    Returns
    -------
    bool 
        True means we should rerun this case because at least one check failed, 
        False means all checks were passed.
    '''
    # decide where to look for the output 
    output_found = False
    rerun_case = False   # update later if it fails a check
    delay_case = False    # update later if completed file doesn't exist
    delete_case = False   # update later if it's a local dir case that needs a rerun
    localdir = os.path.join(pars['model-dir'], "scepter_output", row['newrun_id_field_full'])
    if "aws_bucket" in row.index:
        if "aws_save" in row.index:
            if row['aws_save'] == "move":
                awsdir = os.path.join(row['aws_bucket'], row['newrun_id_field_full'])
    else:
        awsdir = None

    if awsdir is not None:
        # -------------
        import fsspec
        import s3fs
        # -------------
        fs = s3fs.S3FileSystem()
        if fs.exists(awsdir): # then conduct the three checks
            output_found = True
            # check the completed file
            completed_file_check = fs.exists(os.path.join(awsdir, completed_fn))
            
            # check the results file 
            if fs.exists(os.path.join(awsdir, check_results_fn)):
                check_res_dict = check_results_read(os.path.join(awsdir, check_results_fn), file_on_s3=True)
                check_res_check = all(check_res_dict.values())
            else:
                rerun_case = True
            # check the log file
            # (sometimes the model fails early, but the duration check misses it 
            #  because, somehow, the 'target' and 'model' durations are 10000... 
            #  equivalent to the spinup. So we manually repeat the check here)
            if fs.exists(os.path.join(awsdir, check_logs_fn)):
                check_log_dict = check_logs_read(os.path.join(awsdir, check_logs_fn), file_on_s3=True)
                duration_offset_frac = abs((row['duration'] - check_log_dict['model']) / row['duration'])
                check_dur_check = bool(duration_offset_frac < duration_threshold_frac)
            else:
                rerun_case = True

            # check all the checks if we haven't found a reason to rerun yet
            if not rerun_case: 
                # return True if any are False (e.g., not all are True)
                # otherwise return false
                rerun_case = not all([completed_file_check, check_res_check, check_dur_check])
        
        # ****************************************************************
        # 
        # --- TK some runs are completing correctly, but not being moved to AWS
        #     and I don't know why! this reruns them if those runs can't 
        #     be found on AWS... return to this if we can confirm where 
        #     the issue is...
        # 
        # otherwise it's not on AWS, so we check if it should be ...
        elif (row['aws_save'] == "save") or (row['aws_save'] == "move"):   
            rerun_case = True
        # ****************************************************************

    if not output_found: # then try the localdir 
        if os.path.isdir(localdir): # then conduct the three checks
            output_found = True
            # check the completed file
            completed_file_check = os.path.isfile(os.path.join(localdir, completed_fn))
            
            # check the results file 
            if os.path.isfile(os.path.join(localdir, check_results_fn)):
                check_res_dict = check_results_read(os.path.join(localdir, check_results_fn), file_on_s3=False)
                check_res_check = all(check_res_dict.values())
            else:
                rerun_case = True
            # check the log file
            # (sometimes the model fails early, but the duration check misses it 
            #  because, somehow, the 'target' and 'model' durations are 10000... 
            #  equivalent to the spinup. So we manually repeat the check here)
            if os.path.isfile(os.path.join(localdir, check_logs_fn)):
                check_log_dict = check_logs_read(os.path.join(localdir, check_logs_fn), file_on_s3=False)
                duration_offset_frac = abs((row['duration'] - check_log_dict['model']) / row['duration'])
                check_dur_check = bool(duration_offset_frac < duration_threshold_frac)
            else:
                rerun_case = True

            # check all the checks if we haven't found a reason to rerun yet
            if not rerun_case: 
                # return True if any are False (e.g., not all are True)
                # otherwise return false
                rerun_case = not all([completed_file_check, check_res_check, check_dur_check])
            
            if rerun_case and not delay_case:   # if we need to rerun a local dir case, then make a note to delete it locally
                # check to see if the latest update was recent
                stale = is_last_update_older_than(dir_path = os.path.join(localdir, 'flx'), minutes = stale_threshold_minutes)
                if stale:
                    delete_case = True 
                else: 
                    delay_case = True
        else: # no result was found so we rerun it
            rerun_case = True

    return rerun_case, delay_case, delete_case


def allrows_rerun_check(
    df_batch: pd.DataFrame,
    pars: dict,
    completed_fn: str = "completed.res",
    check_results_fn: str = "check_results.res",
    check_logs_fn: str = "check_logs.res",
    duration_threshold_frac: str = 0.2,
    stale_threshold_minutes: float = 15,
) -> pd.DataFrame:
    '''
    Check all rows in the batch dataframe for whether we need to rerun 
    them. Return a pandas dataframe with columns indicating whether 
    a rerun is needed and how many reruns (including the impending one)
    have been performed already

    Parameters
    ----------
    df_batch : pd.DataFrame
        Pandas dataframe of the batch*.csv file. 
    pars : dict
        dictionary of the parameters in the parameter yaml file
    completed_fn : str
        name of the completed.res file 
    check_results_fn : str
        name of the check_results file 
    check_logs_fn : str
        name of the check_logs file 
    duration_threshold_frac : float
        The absolute fractional difference allowed between the target duration of 
        the simulation and the actual duration. 
    stale_threshold_minutes : float
        [minutes] since the last file update in the localdir/flx directory before 
        we deem the run stale
    
    Returns
    -------
    pd.DataFrame
        Returns the input dataframe with updated columns for rerun_needed and 
        rerun_n 
    '''
    # --- loop through all rows
    rerun_me = []
    delay_me = []  # the "completed.res" file doesn't exist, so we delay until delay_max
    delete_me = []  # local dirs to delete so we don't overwrite them when we re-run them 
    for index, row in df_batch.iterrows():

        rerun_result, delay_result, delete_result = checkrow_for_rerun(
            row = row,
            pars = pars,
            duration_threshold_frac = duration_threshold_frac,
            completed_fn = completed_fn, 
            check_results_fn = check_results_fn,
            check_logs_fn = check_logs_fn,
        )

        rerun_me.append(rerun_result)
        delay_me.append(delay_result)
        delete_me.append(delete_result)

    # add result to df
    df_batch = df_batch.assign(rerun_needed = rerun_me)
    df_batch = df_batch.assign(delay_needed = delay_me)
    df_batch = df_batch.assign(delete_localCase = delete_me)
    
    # track the number of reruns
    if "rerun_n" not in df_batch: 
        # then create the column, setting the number of reruns to zero where
        # rerun_needed is false, and 1 where it's true
        df_batch["rerun_n"] = np.where(df_batch['rerun_needed'] & ~df_batch['delay_needed'], 1, 0)
    else:
        # if the "rerun_n" column already exist, add one everywhere 
        # rerun_needed is True
        df_batch.loc[df_batch['rerun_needed'] & ~df_batch['delay_needed'], "rerun_n"] += 1
    
    # track the number of delays
    if "delay_n" not in df_batch: 
        # then create the column, setting the number of delays to zero where
        # delay_needed is false, and 1 where it's true
        df_batch["delay_n"] = np.where(df_batch['delay_needed'], 1, 0)
    else:
        # if the "delay_n" column already exist, add one everywhere 
        # delay_needed is True
        df_batch.loc[df_batch['delay_needed'], "delay_n"] += 1

    return df_batch



# %% 

# --- TROUBLESHOOT ----
# 
# parameter_yaml = "batch-meanAnnliming_fert_fixedRate-base.yaml"
# multiyear = False
# parameter_yaml_dir = "/home/tykukla/aglime-swap-cdr/scepter/parameters"
# completed_fn = "completed.res"
# check_results_fn = "check_results.res"
# check_logs_fn = "check_logs.res"
# duration_threshold_frac = 0.2
# max_reruns = 5
# max_delays = 5
# rerun_delay = 2
# multiyear = False
# workflow_name_runmultiple = "scepter-pyworkflow.yaml"
# bleed_delay_runmultiple = 15
# 
# ----------------------

# %% 
# .. TK function to 
# [1] read in the df_batch
# [2] run the allrows_rerun_check fxn
# [3] if any rerun_needed where rerun_n < max_reruns, then run `run_multiple` for those cases
# [4] save the batch file with `RERUN_` at the start (overwrite it for every iter )
# [5] FUTURE: if basic rerunning doesn't work, tweak one of the inputs to see if that's a problem...
# 
# (inputs: max_reruns: int,
#          rerun_delay: int,   # (how long to wait after the initial run to check for a rerun...)
#          parameter_yaml: str,
#          multiyear: bool,
#          parameter_yaml_dir: str="parameters",
# )
#

def retry_failed_runs(
    max_reruns: int,
    max_delays: int, 
    rerun_delay: float,
    parameter_yaml: str,
    multiyear: bool,
    maindir: str="/home/tykukla/aglime-swap-cdr",
    parameter_yaml_subdir: str="inputs/scepter/params",
    completed_fn: str = "completed.res",
    check_results_fn: str  = "check_results.res",
    check_logs_fn: str  = "check_logs.res",
    duration_threshold_frac: float = 0.2,
    workflow_name_runmultiple: str = "scepter-pyworkflow.yaml",
    bleed_delay_runmultiple: int = 15,
    stale_threshold_minutes : float = 20,
    skip_initial_delay: bool=False,
):
    '''
    Retry failed simulations. Checks for which simulations failed versus which are still running.
    Retries the failed ones and waits max_delays * rerun_delay for the latest file to be updated
    more than stale_threshold_minutes ago before a retry. 

    Parameters
    ----------
    max_reruns : int
        [] max attempts at rerunning a case before we give up on it
    max_delays : int
        [] max attempts at waiting for a case to finish (defined by last file update more
        than stale_threshold_minutes ago) before we give up on it
    rerun_delay : float
        [minutes] time to wait before we attempt to rerun simulations
    parameter_yaml : str
        name of the parameter file for this batch. Something like
        "batch_pars.yaml". 
    multiyear : bool
        True means it's a multiyear simulation (e.g., multiple iters stitched
        into a composite), necessary to know for crafting the full run IDs. 
        False means it's a standard simulation, not a composite.
    maindir : str
        location of the inputs directory (usually '/my/path/to/aglime-swap-cdr')
    parameter_yaml_subdir : str
        location of the parameter_yaml file
    completed_fn : str
        name of the completed.res file 
    check_results_fn : str
        name of the check_results file 
    check_logs_fn : str
        name of the check_logs file 
    duration_threshold_frac : float
        The absolute fractional difference allowed between the target duration of 
        the simulation and the actual duration. 
    workflow_name_runmultiple : str
        name of the argo workflow .yaml file for the run_multiple function
    bleed_delay_runmultiple : int
        [seconds] delay between each argo submit command for the run_multiple function
    stale_threshold_minutes : float
        [minutes] since the last file update in the localdir/flx directory before 
        we deem the run stale
    skip_initial_delay : bool
        True means we skip the initial rerun_delay. Useful if the runs are done and 
        we are doing this cleanup step after the fact. False (default) keeps the 
        initial delay

    Returns
    -------

    '''
    # *****************************
    # apply the rerun delay
    if not skip_initial_delay:
        time.sleep(rerun_delay * 60)    # convert minutes to seconds
    # *****************************

    # --- read in the parameter file
    # create parameter file path
    parameterfile = os.path.join(parameter_yaml_subdir, parameter_yaml)
    # check system arguments, or set default
    with open(parameterfile, "r") as file:
        pars = yaml.safe_load(file)

    # --- read in the batch.csv file
    df_batch = pd.read_csv(os.path.join(pars['batch-input-dir'], pars['batch-input']))

    # --- add column for full run ID
    if multiyear:
        df_batch["newrun_id_field_full"] = df_batch['newrun_id'] + f"_composite_field"
        df_batch["newrun_id_lab_full"] = df_batch['newrun_id'] + f"_composite_lab"
    else:
        df_batch["newrun_id_field_full"] = df_batch['newrun_id'] + "_" + df_batch['dustsp'] + "_field_tau"+df_batch["duration"].astype(float).astype(str).str.replace(".", "p")  # (duration has to be turned into float first because otherwise we miss the decimal pt)
        df_batch["newrun_id_lab_full"] = df_batch['newrun_id'] + "_" + df_batch['dustsp'] + "_lab_tau"+df_batch["duration"].astype(float).astype(str).str.replace(".", "p")  # (duration has to be turned into float first because otherwise we miss the decimal pt)

    df_batch_rerunCheck = allrows_rerun_check(
        df_batch = df_batch,
        pars = pars,
        completed_fn = completed_fn,
        check_results_fn  = check_results_fn,
        check_logs_fn  = check_logs_fn,
        duration_threshold_frac = duration_threshold_frac,
        stale_threshold_minutes = stale_threshold_minutes,
    )
    # 
    # (Note: don't reset the index for the df_batch_*... the run_multiple function only works
    #  in rerun mode if we retain the original index)
    # 
    # check which simulations to run again
    rerun_eligible = np.where(df_batch_rerunCheck['rerun_needed'] & (df_batch_rerunCheck['rerun_n'] <= max_reruns), True, False)
    df_batch_rerunCheck = df_batch_rerunCheck.assign(rerun_eligible = rerun_eligible)
    df_batch_rerun = df_batch_rerunCheck.loc[df_batch_rerunCheck['rerun_eligible'], :].copy()

    # delete local directories that have problems
    df_delete_dirs = df_batch_rerun.loc[df_batch_rerun['delete_localCase'], :].copy()
    scepter_outdir = os.path.join(pars['model-dir'], 'scepter_output')
    if len(df_delete_dirs) > 0:
        # print("deleting problem dirs ... ")
        delete_dirs = list(df_delete_dirs['newrun_id_field_full'].values) + list(df_delete_dirs['newrun_id_lab_full'].values)
        for ddir in delete_dirs:
            if os.path.exists(os.path.join(scepter_outdir, ddir)):
                print(f"deleting {ddir} ... ")
                shutil.rmtree(os.path.join(scepter_outdir, ddir))

    # --- begin reruns ----------------------------------
    continue_rerun_attempts = True
    rerun_iter = 1

    while continue_rerun_attempts:
        if len(df_batch_rerun) == 0:
            continue_rerun_attempts = False
        else:
            print(f"attempting rerun iteration: {rerun_iter} ...")

            # only run the iterations that don't need to be delayed
            # (no delay needed or delay needed but we already hit the max number of allowed delays)
            df_batch_nodelay = df_batch_rerun.loc[~df_batch_rerun['delay_needed'] & (df_batch_rerun['delay_n'] <= max_delays), :].copy()
            
            # re-run each case 
            run_multiple(
                parameter_yaml = parameter_yaml,
                parameter_yaml_subdir = parameter_yaml_subdir,
                maindir = maindir, 
                workflow_name = workflow_name_runmultiple,
                bleed_delay = bleed_delay_runmultiple,
                echo_command = True,
                norun_debug = False,
                rerun_on = True,     # turn this on because it's a rerun case ! 
                df_reruns = df_batch_nodelay,
            )

            # *****************************
            # apply the rerun delay
            time.sleep(rerun_delay * 60)    # convert minutes to seconds
            # *****************************
            
            # check if the reruns worked 
            df_batch_rerunCheck = allrows_rerun_check(
                df_batch = df_batch_rerun,
                pars = pars,
                completed_fn = "completed.res",
                check_results_fn  = "check_results.res",
                check_logs_fn  = "check_logs.res",
                duration_threshold_frac = 0.2,
                stale_threshold_minutes = stale_threshold_minutes,
            )

            # check which simulations to run again
            rerun_eligible = np.where(df_batch_rerunCheck['rerun_needed'] & (df_batch_rerunCheck['rerun_n'] <= max_reruns), True, False)
            df_batch_rerunCheck = df_batch_rerunCheck.assign(rerun_eligible = rerun_eligible)
            df_batch_rerun = df_batch_rerunCheck.loc[df_batch_rerunCheck['rerun_eligible'], :].copy()

            # save the updated rerun rile 
            fnsave = f"RERUN_{rerun_iter}__{pars['batch-input']}"
            df_batch_rerun.to_csv(os.path.join(pars['batch-input-dir'], fnsave))

            # track progress for printout 
            rerun_iter += 1






# %% 