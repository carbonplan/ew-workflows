# Setting up SCEPTER simulations

## Install dependencies
0. Follow the instructions in the `SCEPTER` repo to set up fortran and other SCEPTER dependencies on your machine. 

1. Install python package dependencies in the `ew-workflows/environment.yaml` file. 

## Update hard-coded paths
1. Update the lines to append your `ew-workflows/run_scepter` path in the following scripts within the `SCEPTER` repository.
    - `rock_buff_minPH_multiyear.py`
    - `tunespin_4_newton_inert_buff_v2_defaultdict.py`
    - `singlerun_postproc_detached.py`
    - `rock_app_singleRun-szn.py`
    - `rock_buff_dust-ts_multiyear.py`

2. Update the lines to append your `ew-workflows/run_scepter` path in the following scripts within the `ew-workflows/scripts/scepter` directory.
    - `*/setup/batch/make_batch_input_grainsize+apprate.py`
    - `*/setup/batch/make_batch_input_grainsize+apprate_multiyear.py`
    - `*/setup/batch/make_batch_input_grainsize+apprate_cdrPot.py`

3. Update the default directories for the following functions (or update the code that runs the functions to input the correct directory).
    - in `ew-workflows/run_scepter/argo_helper_fxns.py`
        - `maindir` in `run_multiple` function should be `/path/to/your/ew-workflows`
        - `maindir` in `retry_failed_runs` function should be `/path/to/your/ew-workflows`
    - in `ew-workflows/scripts/scepter/run/run-multiple.py`
        - `maindir` in `ahf.retry_failed_runs` call should be `/path/to/your/ew-workflows`

## Create a batch file
1. Open `ew-workflows/scripts/scepter/setup/batch/make_batch_input_grainsize+apprate.py`.
    - This file defines values that will override the defaults defined in the `SCEPTER/defaults/dict_singlerun.py`

2. Update system path line if you haven't already.

3. Update the `savepath_batch` variable to your `ew-workflows/inputs/scepter/batch`. 

4. Within `const_dict` (dictionary of values held constant for all batch simulations) update your `aws_save` and `aws_bucket` options. If `aws_save = None` then the results will be saved locally in `SCEPTER/scepter_output/<runname>`.

---
# (Argo only)
The following sections only apply if you are using [Argo workflows](https://argoproj.github.io/workflows/) to run batches of SCEPTER simulations. 

## Update argo parameters
1. Open `ew-workflows/inputs/scepter/params/batch-meanAnnliming_fert_fixedRate-base.yaml` and update the following directories (you should only have to update these directories once):
    - `control-script-dir`
    - `batch-input-dir`
    - `model-dir`

2. Update the `batch-input` variable to the name of the batch file you created in the last step.

3. Update the `default-dict` variable to the name of the default dictionary you want to use in `SCEPTER/defaults/dict_singlerun.py` (or review the current default dictionary to make sure it's in line with what you want).

Note: we had the most success with argo if we imposed a delay (`bleed_delay` in the `run_multiple` function in `ew-workflows/run_scepter/argo_helper_fxns.py`) between each run being submitted. 

## Run with argo
1. In terminal, navigate to `ew-workflows/scripts/scepter/run`. Activate the virtual environment you created with the `ew-workflows/environment.yaml` file. 

2. Open `run-multiple.py` and confirm the `run_pars` variable is set to the `batch-*.yaml` file you updated previously.

2. run `python3 run-multiple.py`

---

# Postprocessing
