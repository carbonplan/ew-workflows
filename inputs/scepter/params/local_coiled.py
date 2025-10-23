# ---------------------------------
# 
# params for coiled batch runs
# 
# ---------------------------------

singlerun = {
    "control-script": 'singlerun.py',
    "control-script-dir": '/home/tykukla/ew-workflows/scripts/scepter/run/control_scripts/',
    "batch-input-dir": 's3://carbonplan-carbon-removal/ew-workflows-data/scepter/batch/',
    "default-dict": 'singlerun_default',  # see SCEPTER/defaults/dict_singlerun.py
    "model-dir": '/home/jovyan/SCEPTER/',
    "save-dir": 's3://carbonplan-carbon-removal/ew-workflows-data/TMP-TEST/',

    # [ coiled-specific ; optional (the driver has defaults) ]
    "vm_type": "m8g.xlarge",
    "scheduler_vm_type": "m8g.xlarge",
    "arm": True,
    "region": "us-west-2",
    "spot_policy": "spot_with_fallback",
}

