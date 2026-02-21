#!/usr/bin/env python3
"""
Driver to submit individual dicts as per-task env var dicts to Coiled Batch.
Script is given a name of a batch directory of jsons where each json corresponds
to a single run. The full path to the batch dir comes from combining batchDirName
with `local_coiled.py`.

Usage:
  uv run python3 batch-driver-dicts.py <batchDirName> <paramName>
Example:
  uv run python3 batch-driver-dicts.py cec_simpleTitrate/caOnly/ simplerun

Where <batchDirName> is the name of a directory of .json files to be read in as dicts (1 json per run)
and <paramName> is the key name of a dictionary in inputs/scepter/params/local_coiled.py.

This script will:
 - import local_coiled params dict
 - read filenames from the dir of jsons
 - create commands for the batch of runs
 - set TASK_INPUT_KEYS and SAVE_DIR job-level env vars so the worker knows where to persist
"""
# %% 
import argparse
import pandas as pd
import os
from pathlib import Path
from urllib.parse import urljoin

import coiled

from ew_workflows import coiled_helper_fxns as chf
# %% 

MAX_TASKS_PER_BATCH = 300   # avoid overloading coiled with too many tasks at once

# ====================================================================================================
# # [ DEBUG: fake sys.args ]
# import sys
# sys.argv = ["batch-driver-dicts.py", "cec_simpleTitrate/caOnly", "simplerun"]
# ====================================================================================================

# --- parse inputs
parser = argparse.ArgumentParser()
parser.add_argument("batchName")
parser.add_argument("paramName")
args = parser.parse_args()
# %% 
# --- load params dict
params = chf.load_params(args.paramName)
    
# Allow overriding the bucket/dir via env for local dry-runs
bucket_dir = os.environ.get("BATCH_INPUT_DIR_OVERRIDE") or params.get("batch-input-dir")
if not bucket_dir:
    raise ValueError("batch-input-dir not found in params")

json_list = chf.list_json_files(os.path.join(bucket_dir, args.batchName), return_type='dict')
print(f"Read {len(json_list)} json names from {bucket_dir}/{args.batchName}")

# %% 
# --- submit batch -------------------
worker_script="batch-worker-dicts.py"
container_home = params.get("container-home", "/home/jovyan/")

# Configure submission
save_dir = params.get("save-dir")  # try top level first
if not save_dir:
    # if not at top level, try under singlerun dict
    save_dir = params.get("singlerun", {}).get("save-dir")
if not save_dir:
    # fallback to default
    save_dir = "/scratch/batch"

job_env = {
    "SAVE_DIR": save_dir,
    "COILED_SAVE_DIR": save_dir,  # set both to ensure it's picked up
    # --- scepter-specific
    "DEFAULT_DICT": params.get("default-dict", ""),
    "MODEL_DIR": params.get("model-dir", ""),
    "RUN_SCRIPT": params.get("run-script", ""),
}

# %% 
# If DRY_RUN is set, don't call Coiled; just print a summary and return a fake response
if os.environ.get("DRY_RUN"):
    print("DRY_RUN: would submit with the following parameters:")
    print(f"  tasks: {len(json_list)}")
    print(f"  job_env: {job_env}")
    print(f"  vm_type: {params.get('vm_type', 'm8g.xlarge')}")
    print(f"dry_run: TRUE ; n_tasks: {len(json_list)}")
else:
    all_results = []
    for i, task_chunk in enumerate(chf.tasklist_to_chunks(json_list, MAX_TASKS_PER_BATCH)):
        # Submit using map_over_task_var_dicts so each task receives its env vars
        res = coiled.batch.run(
            # command=f"python3 {worker_script}",
            
            command=(     # [ better for dev, but eventually we should bake this into the container ]
                f"git clone https://github.com/carbonplan/ew-workflows.git {container_home}/ew-workflows && "    # clone in ew-workflows and install
                f"pip install -e {container_home}/ew-workflows && "
                f"python3 {worker_script}"       # run the batch-worker
            ),

            map_over_task_var_dicts=task_chunk,
            name=chf.sanitize_job_name(f"b{i}-{args.batchName}"),
            vm_type=params.get("vm_type", "m8g.xlarge"),
            scheduler_vm_type=params.get("vm_type", "m8g.xlarge"),
            arm=params.get("arm", True),
            container=params.get("container", None),
            region=params.get("region", "us-west-2"),
            tag=params.get("tag", {'Project': 'scepter'}),
            spot_policy=params.get("spot_policy", "spot_with_fallback"),
            forward_aws_credentials=False,   # i have an aws instance profile w/ coiled, so no need to pass creds
            # files=[
            #     str(Path(__file__).resolve().parents[3]),  # sync the whole ew-workflows repo
            # ],
            allow_cross_zone=True,
            # env=job_env,
            env={**job_env, "PYTHONPATH": "ew-workflows/src"},  # ensure ew_workflows importable
        )
        print("Submit result:", res)
        all_results.append(res)

        coiled.batch.status()


