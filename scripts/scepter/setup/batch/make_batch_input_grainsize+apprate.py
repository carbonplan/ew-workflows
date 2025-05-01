# %%
# ---------------------------------------------------
# 
# Generate batch input .csv files for SCEPTER run
# 
# provide var vectors and assume we want every 
# combination of them, or by site
# 
# T Kukla (CarbonPlan, 2024)
# 
# ---------------------------------------------------

import itertools
import os
import sys

import numpy as np
import pandas as pd

# --- read in batch helper functions  
sys.path.append(os.path.abspath('/home/tykukla/ew-workflows/run_scepter'))
import batch_helperFxns as bhf
# ---

# %% 
# --- USER INPUTS
# [1] vars to update, constant for all runs
fertLevel = "low"    # name for how much fertilizer is applied
dustsp = "cc"      # name for dust species to apply (must be from accepted list)
extra_tag = "base_psdfull"  # another distinguishing tag
pref = f"{fertLevel}Fert_{dustsp}_{extra_tag}"
clim_tag = None   # [string] year-span (e.g., 1950-2020) for clim input if climate files are used
                  # (if clim files are not used, set to None)
# save vars
file_prefix = f"meanAnn_{dustsp}_shortRun_{extra_tag}_{fertLevel}Fert_gs+apprate"  # prefix of output run names
fn = file_prefix + "_v0.csv"
savepath_batch = "/home/tykukla/ew-workflows/inputs/scepter/batch"
multi_run_split = False   # whether to split the csv into multiple files
max_iters_per_set = 20    # [int] the number of runs per csv (only used if multi_run_split is True)
# **************************
# --- do not change
fert_dict = {
    "hi": 35.0,    # 30.0
    "low": 3.0,    # 6.0
    "no": 0.0,
    }
# **************************

const_dict = {
    "duration": 15,  # [yr] duration of run (starts from earliest year)
    "dustsp": dustsp,
    "dustsp_2nd": "amnt",
    "dustrate_2nd": fert_dict[fertLevel],
    "add_secondary": False,
    "imix": 3, # mixing style (1=fickian; 2=homogeneous; 3=tilling)
    "singlerun_seasonality": False,
    "include_psd_full": True,   # False,
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
# [3] vars to vary within site (we'll get every combination of these two)
# dustrate_ton_ha_yr = [0.1, 0.3, 0.6, 1, 2, 5, 7, 10, 15, 25, 35, 45, 60, 100]
dustrate_ton_ha_yr = [0.1, 0.3, 0.6, 1, 1.5, 2, 4, 8, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100]
dustrad_um = [1, 3, 6, 10, 15, 20, 30, 40, 50, 60, 70, 80, 100, 125, 150, 175, 200, 250, 300]
all_combinations = {
    "dustrate": [x * 100 for x in dustrate_ton_ha_yr],  # [ton ha-1 yr-1 * 100 = g m-2]   
    "dustrad": dustrad_um  # [diameter, microns] i think this gets applied to gbas and amnt equally (though amnt is fast-reacting so maybe not a big deal? )
}


# %% 
# --- BUILD DATAFRAME AND SAVE
df = bhf.build_df(pref, const_dict, sites, by_site, all_combinations, add_ctrl=True)

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
