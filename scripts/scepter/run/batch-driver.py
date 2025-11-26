#!/usr/bin/env python3
"""
Driver to submit CSV rows as per-task env var dicts to Coiled Batch.

Usage:
  uv run python3 batch-driver.py <batchName> <paramName>
Example:
  uv run python3 batch-driver.py _test-batch-n6.csv singlerun
  uv run python3 batch-driver.py _test-batch-richards-n6.csv singlerun

Where <batchName> is the CSV filename (in S3 under batch-input-dir configured in local_coiled.py)
and <paramName> is the key name of a dictionary in inputs/scepter/params/local_coiled.py.

This script will:
 - import local_coiled params dict
 - fetch CSV from S3 (or local file when S3 not configured)
 - create one per-row dict, and submit via coiled.batch.run with map_over_task_var_dicts
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
# sys.argv = ["batch-driver.py", "longrun_monthly_ltm_v0_test.csv", "singlerun"]
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

rows = chf.read_csv_from_s3_or_local(bucket_dir, args.batchName)
print(f"Read {len(rows)} rows from {bucket_dir}/{args.batchName}")


# --- submit batch -------------------
worker_script="batch-worker.py"
container_home = params.get("container-home", "/home/jovyan/")

# create task dicts
task_dicts = chf.create_task_dicts(rows)

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
    print(f"  tasks: {len(task_dicts)}")
    print(f"  keys: {keys}")
    print(f"  job_env: {job_env}")
    print(f"  vm_type: {params.get('vm_type', 'm8g.xlarge')}")
    print(f"dry_run: TRUE ; n_tasks: {len(task_dicts)}")

else:
    all_results = []
    for i, task_chunk in enumerate(chf.tasklist_to_chunks(task_dicts, MAX_TASKS_PER_BATCH)):
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
# %%
