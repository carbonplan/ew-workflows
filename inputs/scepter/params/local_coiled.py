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
    "multiyear": False,     # whether we're in multi-year mode (affects path building)

    # [ coiled-specific ; optional (the driver has defaults) ]
    "container-home": "/home/jovyan/",
    "vm_type": "m6i.large",      # make sure it's compatible with the container and arm selection  # (https://instances.vantage.sh)
    "arm": False,
    "region": "us-west-2",
    "spot_policy": "on-demand",    # "spot_with_fallback" [ can be deleted, but much cheaper ]
    "tag": {"Project": "scepter"},
}

multiyear = {
    "run-script": 'rock_buff_dust-ts_multiyear.py',
    "batch-input-dir": 's3://carbonplan-carbon-removal/ew-workflows-data/scepter/batch/',
    "default-dict": 'specifydust_default',  # see SCEPTER/defaults/dict_singlerun.py
    "model-dir": '/home/jovyan/SCEPTER/',
    "save-dir": 's3://carbonplan-carbon-removal/ew-workflows-data/TMP-TEST/',
    "container": 'quay.io/carbonplan/ew-workflows:latest',
    "multiyear": True,     # whether we're in multi-year mode (affects path building)

    # [ coiled-specific ; optional (the driver has defaults) ]
    "container-home": "/home/jovyan/",
    "vm_type": "m6i.large",      # make sure it's compatible with the container and arm selection (https://instances.vantage.sh)
    "arm": False,
    "region": "us-west-2",
    "spot_policy": "on-demand",    # "spot_with_fallback" [ can be deleted, but much cheaper ]
    "tag": {"Project": "scepter"},
}

ph_react = {
    "run-script": 'rock_buff_minPH_multiyear.py',
    "batch-input-dir": 's3://carbonplan-carbon-removal/ew-workflows-data/scepter/batch/',
    "default-dict": 'ph_react_default',  # see SCEPTER/defaults/dict_singlerun.py
    "model-dir": '/home/jovyan/SCEPTER/',
    "save-dir": 's3://carbonplan-carbon-removal/ew-workflows-data/TMP-TEST/',
    "container": 'quay.io/carbonplan/ew-workflows:latest',
    "multiyear": True,     # whether we're in multi-year mode (affects path building)

    # [ coiled-specific ; optional (the driver has defaults) ]
    "container-home": "/home/jovyan/",
    "vm_type": "m6i.large",      # make sure it's compatible with the container and arm selection (https://instances.vantage.sh)
    "arm": False,
    "region": "us-west-2",
    "spot_policy": "on-demand",    # "spot_with_fallback" [ can be deleted, but much cheaper ]
    "tag": {"Project": "scepter"},
}