# %%
# ---------------------------------------------------
# 
# Generate batch input .csv files for SCEPTER run
# for morris method sensitivity analysis using 
# SALib in python
# 
# T Kukla (CarbonPlan, 2025)
# 
# ---------------------------------------------------

import itertools
import os
import sys

import numpy as np
import pandas as pd
from SALib.sample import morris as morris_sampler
from SALib.analyze import morris as morris_analyzer

# --- read in batch helper functions  
sys.path.append(os.path.abspath('/home/tykukla/ew-workflows/run_scepter'))
import batch_helperFxns as bhf
# ---

SAVE_ON = True    # whether to save the SA settings 
SAsetting_path = "/home/tykukla/ew-workflows/scripts/scepter/setup/batch/SA_inputs"

# %% 
# --- notes on the variables we're sampling
# 
# [ dustrate ] annual mean dust application flux [g m-2 yr-1]
# [ taudust ] duration of the dust spreading event [yr]
# [ dustrad ] radius of the dust particles [micron]
# [ qrun ] water infiltration rate [m yr-1]
# [ mat ] mean annual temperature [degC]
# [ dust_mixdep ] mixing depth of the dust particles [m]
# [ dustrate_2nd ] amount of secondary dusts (amnt) added [g m-2 yr-1]

# %% 
# --- Set up the problem; define model inputs
# 
# note, multiply dustrate by 100 to get from ton ha-1 yr-1 to g m-2 yr-1
# 
problem = {
    "num_vars": 8,
    "names": [ "dustrate",         "taudust",   "dustrad", "qrun",      "mat",   "dust_mixdep", "dustrate_2nd", "soilmoisture_surf"],
    "bounds": [[0.1*100, 20*100],  [0.03, 0.1], [10, 200],  [0.1, 1.5], [3, 30], [0.05, 0.4],   [0.0, 35.0],     [0.1, 0.7] ]  
}
# define control values (a control run is added to the input arrays at the end in the zero index)
param_values_ctrl = {
    "dustrate": 0,
    "taudust": 0.05,
    "dustrad": 100,
    "qrun": 0.3,
    "mat": 8.,
    "dust_mixdep": 0.25,
    "dustrate_2nd": 0.,
    "soilmoisture_surf": 0.63,
}


# --- generate samples using morris sampler
param_values = morris_sampler.sample(
    problem, 
    N=750, 
    num_levels=4, 
    optimal_trajectories=50,
    seed = 1111,
)
len(param_values)


# %% 
# --- save the SA settings --------------------------------------
if SAVE_ON:
    exp_name = f"SAmorris{len(param_values)}_{problem['num_vars']}_v0"
    savesettings_here = os.path.join(SAsetting_path, exp_name)
    if not os.path.exists(savesettings_here):
        os.makedirs(savesettings_here)  
        
    # save array
    np.save(os.path.join(savesettings_here, 'param_values.npy'), param_values)
    # save dict
    pd.to_pickle(problem, os.path.join(savesettings_here, "problem.pkl"))
    # with open(os.path.join(savesettings_here, "problem.pkl"), "wb") as f:
    #     pickle.dump(problem, f)
# ---------------------------------------------------------------



# %% 
# --- USER INPUTS
# [1] vars to update, constant for all runs
fertLevel = "SA"    # name for how much fertilizer is applied
dustsp = "cc"      # name for dust species to apply (must be from accepted list)
extra_tag = f"morris_{problem['num_vars']}_{len(param_values)}"  # another distinguishing tag
pref = f"{fertLevel}Fert_{dustsp}_{extra_tag}"
clim_tag = None   # [string] year-span (e.g., 1950-2020) for clim input if climate files are used
                  # (if clim files are not used, set to None)
# save vars
file_prefix = f"meanAnn_{dustsp}_shortRun_{extra_tag}_{fertLevel}Fert_gs+apprate"  # prefix of output run names
fn = file_prefix + "_v0.csv"
savepath_batch = "/home/tykukla/ew-workflows/inputs/scepter/batch"
multi_run_split = False   # whether to split the csv into multiple files
max_iters_per_set = 20    # [int] the number of runs per csv (only used if multi_run_split is True)

# --- constant values
const_dict = {
    # --- MULTI-YEAR SPECIFIC ---
    # "dust_ts_fn": f"{dustsp}_15yr_1app_no2nd_001.csv",
    # ---

    "duration": 5,  # [yr] duration of run (starts from earliest year)
    "dustsp": dustsp,
    "dustsp_2nd": "amnt",
    "add_secondary": False,
    "imix": 3, # mixing style (1=fickian; 2=homogeneous; 3=tilling)
    "singlerun_seasonality": False,
    "include_psd_full": False,   # False,
    "include_psd_bulk": False,
    'poro_iter_field': False,      # [default='false'] porosity iteration
    'poro_evol': True,            # [default='false'] porosity evolves with solid phase dissolution
    "cec_adsorption_on": False, # True,
    "climatedir": "NA",

    # --- surface area and psd rules
    'sa_rule1': False,       # [True, False, "spinvalue"] SA decreases as porosity increases
    'sa_rule2': True,       # [True, False, "spinvalue"] SA increases as porosity increases
    'psdrain_log10_sd': 0.05, # [] log 10 standard deviation for psd
    'psdrain_wt': 1.0,       # [] weight for the psd
    'use_psdrain_datfile': False,  # False means we construct the PSD based on inputs, rather than copy an existing data file
    'include_roughness_sa': True,  # [True, False] whether to apply roughness calculation from Navarre-Sitchler & Brantley, 2007 
                                   # (roughness factor = (beta / a)^0.33) where beta is particle radius in m and a is BET measurement resolution in m, taken to be 10^-10 
    
    # # --- update CEC?
    # 'cec_update_from_spinup': True,   # [True, False] whether to update CEC and alpha vars relative to the spinup value (False means no change to cec.in is made)
    # 'cec': 15,                   # [cmol kg-1] [only used if cec_update_from_spinup == True] cation exchange capacity
    # 'alpha': 0.1,                   # [only used if cec_update_from_spinup == True]
    
    # --- compute specific
    'aws_save': "move",              # ["move", "copy", None] whether to "move" file to aws, just copy it, or nothing at all
    'aws_bucket': "s3://carbonplan-carbon-removal/SCEPTER/scepter_output_scratch/",  # where to save at AWS (only used if 'aws_save'=True)
    
    # --- which executable to use
    'scepter_exec_name': 'scepter'  # ['scepter', 'scepter_rateA', ...]
}

# %% 
# [2] vars to vary by site
## UNCOMMENT BELOW FOR 2 SITES:
# sites = ['site_311a', 'site_311b']
# by_site = {   # values must have same order as 'sites' var
#     "cec": [21.10329, 6.96125],
#     "spinrun": ["site_311a_pr9_spintuneup4", "site_311b_pr9_spintuneup4"],
#     "climatefiles": ["site_311a", "site_311b"]  # these serve as the site name when there is no cliamte file to use
# }
# -----
## UNCOMMENT BELOW FOR SINGLE SITE:
sites = ['site_311a']
by_site = {   # values must have same order as 'sites' var
    # "cec": [21.10329],
    "spinrun": ["site_311a_pr9_spintuneup4"], # ["spinupLime_cc_v0_site_311a_app_60p0_psize_100_cc_tau100p0"], # ["site_311a_pr9_YoshiCEC_spintuneup4"], # ["spinupLime_cc_v2_site_311a_app_150p0_psize_100_cc_tau100p0"] 
    "climatefiles": ["site_311a"]  # these serve as the site name when there is no cliamte file to use
}

# %% 
# --- BUILD DATAFRAME AND SAVE
df = bhf.build_df_sa(pref, const_dict, sites, by_site, param_values, problem, param_values_ctrl, add_ctrl=True)

# --- add the particle size distribution inputs if relevant 
# 
if const_dict['include_psd_full'] or const_dict['include_psd_bulk']:
    df['psdrain_meanRad'] = df['dustrad'] / 1e6  # (convert um to m)
df


# %% 
# save 
bhf.save_df(df, savepath_batch, fn, multi_run_split, max_iters_per_set)
# %%
fn

# %%
