# ---------------------------------
# 
# params for coiled batch runs
# 
# ---------------------------------

singlerun = {
    "control-script": 'singlerun.py',
    "run-script": 'rock_app_singleRun-szn.py',
    "control-script-dir": '/home/tykukla/ew-workflows/scripts/scepter/run/control_scripts/',
    "batch-input-dir": 's3://carbonplan-carbon-removal/ew-workflows-data/scepter/batch/',
    "default-dict": 'singlerun_default',  # see SCEPTER/defaults/dict_singlerun.py
    "model-dir": '/home/jovyan/SCEPTER/',
    "save-dir": 's3://carbonplan-carbon-removal/ew-workflows-data/TMP-TEST/',
    "container": 'quay.io/carbonplan/ew-workflows:latest',

    # [ coiled-specific ; optional (the driver has defaults) ]
    "container-home": "/home/jovyan/",
    "vm_type": "m5.xlarge",      # make sure it's compatible with the container and arm selection
    "arm": False,
    "region": "us-west-2",
    "spot_policy": "spot_with_fallback",
    "tag": {"Project": "scepter"},
}

