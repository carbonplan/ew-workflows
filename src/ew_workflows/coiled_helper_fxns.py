# --- helper functions for coiled workflows ---
import argparse
import json
import pandas as pd
import os
from pathlib import Path
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


def read_csv_from_s3_or_local(bucket_dir, name):
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

    rows = df.to_dict(orient="records")
    return rows


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
