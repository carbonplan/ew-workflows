#!/usr/bin/env python3
"""
Worker script for batch runs.
When submitted via map_over_task_var_dicts, the driver will set per-task env vars for each column.
The driver will also set the TASK_INPUT_KEYS job-level env var listing which keys to persist,
and SAVE_DIR job-level env var specifying where to save per-task outputs.

This script collects those env vars and writes a JSON/text file with the key/values.
"""
import os
import json
import subprocess
import sys
import traceback
import tempfile

from datetime import datetime, timezone


def main():
    # Get keys that the driver indicated
    keys_csv = os.environ.get("TASK_INPUT_KEYS")
    if keys_csv:
        keys = [k for k in keys_csv.split(",") if k]
    else:
        # If not provided, take all env vars except a small blacklist
        blacklist = {"HOME", "PATH", "SHELL", "USER", "JOB_ID", "COILED_BATCH_TASK_INPUT", "TASK_INPUT_KEYS", "SAVE_DIR"}
        keys = [k for k in os.environ.keys() if k not in blacklist]

    # Collect a dict of values (strings)
    result = {k: os.environ.get(k) for k in keys}
    # Add a timestamp
    result["_task_started_at"] = datetime.now(timezone.utc).isoformat()


    # -------------------------------------------------------
    # [ scepter troubleshoot ]
    def build_run_cmd(script_path, results, static_inputs=None):
        # results: dict of key -> value
        # script_path: script name or absolute path
        args = ["python3", script_path]
        args += static_inputs if static_inputs else []
        for k, v in results.items():
            flag = f"--{k.replace('_', '-') }"  # transform underscores to hyphens if you like
            if v is True:
                # boolean True -> flag only
                args.append(flag)
            elif v is False or v is None:
                # skip False/None
                continue
            elif isinstance(v, (list, tuple)):
                # repeat option for each element
                for item in v:
                    args.extend([flag, str(item)])
            else:
                args.extend([flag, str(v)])
        # return the constructed argument list (may be just [python3, script_path])
        return args
    # check if the model-dir path exists and set static inputs
    model_dir = os.environ.get("MODEL_DIR", None)
    run_script = os.environ.get("RUN_SCRIPT", None)
    outdir = os.path.join(model_dir, "scepter_output/")  # (note: the trailing '/' is critical because of lazy coding on the other side)
    default_dict = os.environ.get("DEFAULT_DICT", None)
    # create static inputs list
    static_inputs = [
        f"--modeldir", f"{model_dir}",
        f"--outdir", f"{outdir}",
        f"--default_dict", f"{default_dict}",
    ]
    # log model-dir existence
    if model_dir:
        path_exists = os.path.exists(model_dir)
        result["_model_dir_exists"] = path_exists
        print(f"model-dir: {model_dir}, exists: {path_exists}", flush=True)
    else:
        result["_model_dir_exists"] = "model-dir env var not set" 
        print("model-dir env var not set", flush=True)
    # Only attempt to build and run the command if we have a run_script and model_dir
    if run_script and model_dir:
        runme = os.path.join(model_dir, run_script)
        file_exists = os.path.exists(runme)
        print(f"file exists: {file_exists} -> {runme}")
        if not file_exists:
            print(f"[WARN] run script not found at {runme}; skipping execution", flush=True)
        else:
            mycmd = build_run_cmd(runme, result, static_inputs)
            # mycmd_str = ' '.join(mycmd)
            # mycmd is a list; run without shell to avoid shell-specific parsing
            print(f"Running: {mycmd}", flush=True)
            try:
                subprocess.run(mycmd, check=True)
            except Exception:
                print("[ERROR] subprocess failed", flush=True)
                traceback.print_exc()
    else:
        print("[INFO] RUN_SCRIPT or MODEL_DIR not provided; skipping execution", flush=True)
    # -------------------------------------------------------

    # Determine save dir
    save_dir = os.environ.get("SAVE_DIR") or os.environ.get("COILED_SAVE_DIR") or "/scratch/batch"
    try:
        os.makedirs(save_dir, exist_ok=True)
    except Exception:
        # fallback to temp dir
        import tempfile

        save_dir = os.path.join(tempfile.gettempdir(), "coiled_batch_output")
        os.makedirs(save_dir, exist_ok=True)

    # Determine filename using a compact identifier (timestamp + row_number or first non-empty key)
    id_part = None
    id_part = result.get('row_number', None)
    if not id_part:  # then use timestamp
        id1 = "timestamp"
        id_part = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f")
    else:
        id1 = "rownum"
    # build filename without f-string to avoid accidental tooling filename detection
    out_filename = f"task_{id1}_" + str(id_part) + ".json"

    # Handle S3 paths specially (don't try to open() an s3:// URL)
    if isinstance(save_dir, str) and save_dir.startswith("s3://"):
        # Build bucket/key
        s3_prefix = save_dir.rstrip("/")
        # parse s3://bucket/optional/prefix
        _, rest = s3_prefix.split("s3://", 1)
        parts = rest.split("/", 1)
        bucket = parts[0]
        prefix = parts[1] if len(parts) > 1 else ""
        key = f"{prefix.rstrip('/')}/{out_filename}" if prefix else out_filename
        # First try fsspec (s3fs) since it integrates with fsspec and forwarded creds
        try:
            import fsspec

            s3_path = f"s3://{bucket}/{key}"
            with fsspec.open(s3_path, "w") as f:
                json.dump(result, f, indent=2, ensure_ascii=False)
            print(f"[DONE] wrote {s3_path}", flush=True)
        except Exception:
            # if fsspec isn't installed or fails, try boto3 as a second option
            try:
                import boto3

                s3 = boto3.client("s3")
                body = json.dumps(result, indent=2, ensure_ascii=False).encode("utf-8")
                s3.put_object(Bucket=bucket, Key=key, Body=body)
                print(f"[DONE] wrote s3://{bucket}/{key}", flush=True)
            except Exception:
                print("[ERROR] failed to write output to s3; falling back to local write", flush=True)
                traceback.print_exc()
                # fallback to a local temp write so logs and debugging can continue
                local_out = os.path.join(tempfile.gettempdir(), "coiled_batch_output", out_filename)
                os.makedirs(os.path.dirname(local_out), exist_ok=True)
                try:
                    with open(local_out, "w") as f:
                        json.dump(result, f, indent=2, ensure_ascii=False)
                    print(f"[DONE] wrote fallback {local_out}", flush=True)
                except Exception:
                    print("[ERROR] failed to write fallback output", flush=True)
                    traceback.print_exc()
                    raise
    else:
        out_path = os.path.join(save_dir, out_filename)
        try:
            with open(out_path, "w") as f:
                json.dump(result, f, indent=2, ensure_ascii=False)
            print(f"[DONE] wrote {out_path}", flush=True)
        except Exception:
            print("[ERROR] failed to write output", flush=True)
            traceback.print_exc()
            raise


if __name__ == "__main__":
    main()
