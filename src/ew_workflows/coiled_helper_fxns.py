# --- helper functions for coiled workflows ---
import os
import json
from pathlib import Path
import re
import time 

import argparse
import numpy as np
import pandas as pd
from urllib.parse import urljoin

import coiled


def _write_json_to_path(path, data):
    """Write JSON to an s3:// path (via fsspec or boto3) or a local path."""
    if str(path).startswith("s3://"):
        try:
            import fsspec

            with fsspec.open(path, "w") as f:
                json.dump(data, f)
            return
        except Exception:
            # try boto3
            try:
                import boto3

                _, rest = str(path).split("s3://", 1)
                parts = rest.split("/", 1)
                bucket = parts[0]
                key = parts[1] if len(parts) > 1 else ""
                boto3.client("s3").put_object(Bucket=bucket, Key=key, Body=json.dumps(data).encode("utf-8"))
                return
            except Exception:
                raise
    else:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w") as f:
            json.dump(data, f)
            

def create_task_dicts(rows):
    # rows is a list of dicts (strings)
    # Prepare per-task env var dicts
    return [{k: ("" if v is None else str(v)) for k, v in r.items()} for r in rows]


def find_repo_root():
    """Walk up directories from this script until we find inputs/scepter/params."""
    # Start from this script's directory
    current = Path(__file__).resolve().parent
    while current.parent != current:  # while we can still go up
        if (current / "inputs" / "scepter" / "params").is_dir():
            return current
        current = current.parent
    raise RuntimeError(
        "Could not find ew-workflows root directory (expected inputs/scepter/params). "
        "Are you running from inside the ew-workflows repository?"
    )


def load_params(param_name):
    # Import the local_coiled.py module
    import importlib.util

    # Find repository root (directory containing inputs/scepter/params)
    repo_root = find_repo_root()
    local_path = repo_root / "inputs" / "scepter" / "params" / "local_coiled.py"

    spec = importlib.util.spec_from_file_location("local_coiled", local_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    params = getattr(mod, param_name, None)
    if params is None:
        raise ValueError(f"param {param_name} not found in {local_path}")
    return params


def read_csv_from_s3_or_local(
        bucket_dir: str,
        name: str,
        return_type: str = "list_of_dicts",
):
    """
    Read CSV either from an S3 prefix (s3://bucket/path/) or a local directory.
    Returns a list of dictionaries (one per row) with string values.
    Requires pandas (and s3fs for s3:// paths).
    """
    full_path = bucket_dir.rstrip("/") + "/" + name
    try:
        # pandas delegates s3:// handling to fsspec/s3fs when installed
        df = pd.read_csv(full_path)
        # add row number as a column
        df.insert(0, "row_number", range(1, len(df) + 1))
    except Exception as e:
        # surface a clearer error for missing s3fs/pandas
        raise RuntimeError(
            f"Failed to read CSV at {full_path}. Ensure pandas and s3fs are installed for s3 paths; original error: {e}"
        ) from e

    # normalize and convert values to strings so they can be set as env vars
    df = df.fillna("")
    for col in df.columns:
        # Ensure all values are strings (env vars must be strings)
        df[col] = df[col].astype(str)

    if return_type == "list_of_dicts":
        rows = df.to_dict(orient="records")
        return rows
    elif return_type == "dataframe":
        return df
    else:
        raise ValueError(f"Unsupported return_type: {return_type}")


def sanitize_job_name(filename):
    """Convert a filename into a valid Coiled job name."""
    # Remove file extension
    name = os.path.splitext(filename)[0]
    # Replace invalid characters with dashes
    import re
    return re.sub(r'[^a-zA-Z0-9-]', '-', name)


def submit_batch(rows, params, worker_script="batch-worker.py"):
    parser = argparse.ArgumentParser()
    parser.add_argument("batchName")
    parser.add_argument("paramName")
    args = parser.parse_args()

    # rows is a list of dicts (strings)
    # Prepare per-task env var dicts
    task_dicts = []
    for r in rows:
        # convert all values to str (Coiled expects strings for env vars)
        task_dicts.append({k: ("" if v is None else str(v)) for k, v in r.items()})

    # The list of keys we passed (for the worker to know which keys to persist)
    keys = list(rows[0].keys()) if rows else []

    # Configure submission
    save_dir = params.get("save-dir")  # try top level first
    if not save_dir:
        # if not at top level, try under singlerun dict
        save_dir = params.get("singlerun", {}).get("save-dir")
    if not save_dir:
        # fallback to default
        save_dir = "/scratch/batch"
    
    job_env = {
        "TASK_INPUT_KEYS": ",".join(keys),
        "SAVE_DIR": save_dir,
        "COILED_SAVE_DIR": save_dir  # set both to ensure it's picked up
    }

    # If DRY_RUN is set, don't call Coiled; just print a summary and return a fake response
    if os.environ.get("DRY_RUN"):
        print("DRY_RUN: would submit with the following parameters:")
        print(f"  tasks: {len(task_dicts)}")
        print(f"  keys: {keys}")
        print(f"  job_env: {job_env}")
        print(f"  vm_type: {params.get('vm_type', 'm8g.xlarge')}")
        return {"dry_run": True, "n_tasks": len(task_dicts)}

    # Submit using map_over_task_var_dicts so each task receives its env vars
    res = coiled.batch.run(
        command=f"python3 {worker_script}",
        map_over_task_var_dicts=task_dicts,
        name=sanitize_job_name(args.batchName),
        vm_type=params.get("vm_type", "m8g.xlarge"),
        scheduler_vm_type=params.get("scheduler_vm_type", params.get("vm_type", "m8g.xlarge")),
        arm=params.get("arm", True),
        region=params.get("region", "us-west-2"),
        spot_policy=params.get("spot_policy", "spot_with_fallback"),
        forward_aws_credentials=True,
        allow_cross_zone=True,
        env=job_env,
    )
    return res


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
    skip_duration_check: bool = False,
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
    skip_duration_check : bool
        True if we want to skip the duration check (used for multi year simulations)
    
    Returns
    -------
    bool 
        True means we should rerun this case because at least one check failed, 
        False means all checks were passed.
    '''
    # decide where to look for the output 
    output_found = False
    rerun_case = False    # update later if it fails a check
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
                if not skip_duration_check: 
                    check_res_check = all(check_res_dict.values()) # account for all checks
                else:
                    check_res_check = all(list(check_res_dict.values())[0:4]) # ignore the duration check
            else:
                rerun_case = True
            # check the log file
            # (sometimes the model fails early, but the duration check misses it 
            #  because, somehow, the 'target' and 'model' durations are 10000... 
            #  equivalent to the spinup. So we manually repeat the check here)
            if fs.exists(os.path.join(awsdir, check_logs_fn)):
                check_log_dict = check_logs_read(os.path.join(awsdir, check_logs_fn), file_on_s3=True)
                if not skip_duration_check:
                    duration_offset_frac = abs((row['duration'] - check_log_dict['model']) / row['duration'])
                    check_dur_check = bool(duration_offset_frac < duration_threshold_frac)
                else:
                    check_dur_check = True

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
                if not skip_duration_check: 
                    check_res_check = all(check_res_dict.values())  # account for all checks
                else:
                    check_res_check = all(list(check_res_dict.values())[0:4]) # ignore the duration check
            else:
                rerun_case = True
            # check the log file
            # (sometimes the model fails early, but the duration check misses it 
            #  because, somehow, the 'target' and 'model' durations are 10000... 
            #  equivalent to the spinup. So we manually repeat the check here)
            if os.path.isfile(os.path.join(localdir, check_logs_fn)):
                check_log_dict = check_logs_read(os.path.join(localdir, check_logs_fn), file_on_s3=False)
                if not skip_duration_check:
                    duration_offset_frac = abs((row['duration'] - check_log_dict['model']) / row['duration'])
                    check_dur_check = bool(duration_offset_frac < duration_threshold_frac)
                else:
                    check_dur_check = True
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
    skip_duration_check: bool=False,
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
        dictionary of the parameters in the parameter dict
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
    skip_duration_check : bool
        True if we skip the model duration check (used for multi_year simulations)
    
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
            stale_threshold_minutes = stale_threshold_minutes,
            skip_duration_check = skip_duration_check,
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



def find_failed_or_stale_runs(
        df_batch: pd.DataFrame,
        multiyear: bool,
        pars: dict,
        completed_fn: str = "completed.res",
        check_results_fn: str  = "check_results.res",
        check_logs_fn: str  = "check_logs.res",
        duration_threshold_frac: float = 0.2,
        skip_duration_check: bool=False
)->pd.DataFrame: 
    """
    Identify failed or incomplete runs from a batch CSV.

    Parameters
    ----------
    df_batch : pd.DataFrame
        DataFrame containing the batch CSV rows.
    multiyear : bool
        Whether the runs are in multiyear mode (affects path building).
    """
    # --- add column for full run ID
    if multiyear:
        df_batch["newrun_id_field_full"] = df_batch['newrun_id'] + f"_composite_field"
        df_batch["newrun_id_lab_full"] = df_batch['newrun_id'] + f"_composite_lab"
    else:
        df_batch["newrun_id_field_full"] = df_batch['newrun_id'] + "_" + df_batch['dustsp'] + "_field_tau"+df_batch["duration"].astype(float).astype(str).str.replace(".", "p")  # (duration has to be turned into float first because otherwise we miss the decimal pt)
        df_batch["newrun_id_lab_full"] = df_batch['newrun_id'] + "_" + df_batch['dustsp'] + "_lab_tau"+df_batch["duration"].astype(float).astype(str).str.replace(".", "p")  # (duration has to be turned into float first because otherwise we miss the decimal pt)

    # --- check each row for whether it needs to be rerun -----
    df_batch_rerunCheck = allrows_rerun_check(
        df_batch = df_batch,
        pars = pars,
        completed_fn = completed_fn,
        check_results_fn  = check_results_fn,
        check_logs_fn  = check_logs_fn,
        duration_threshold_frac = duration_threshold_frac,
        stale_threshold_minutes = 0,
        skip_duration_check = skip_duration_check
    )

    return df_batch_rerunCheck