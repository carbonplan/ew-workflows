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
fertLevel = "no"    # name for how much fertilizer is applied
dustsp = "cao"      # name for dust species to apply (must be from accepted list)
n_apps = "5"        # number of annual rock applications
extra_tag = f"cec_depths_{n_apps}yr"  # another distinguishing tag
pref = f"{fertLevel}Fert_{dustsp}_{extra_tag}"
clim_tag = None   # [string] year-span (e.g., 1950-2020) for clim input if climate files are used
                  # (if clim files are not used, set to None)
# save vars
file_prefix = f"clim_{dustsp}_{extra_tag}_{fertLevel}Fert_cec+depth"  # prefix of output run names
fn = file_prefix + "_v0.csv"
savepath_batch = "s3://carbonplan-carbon-removal/ew-workflows-data/scepter/batch"
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
    # --- MULTI-YEAR SPECIFIC ---
    "dust_ts_dir": f"s3://carbonplan-carbon-removal/ew-workflows-data/scepter/dust/",
    "dust_ts_fn": f"{dustsp}_100yr_{n_apps}app_no2nd_001.csv",
    # ---
    # ********************
    # "nz": 30,     # [None, or comment out if we don't want to modify] number of grid cells 
    # ********************
    "duration": 100,  # [yr] duration of run (starts from earliest year)
    "dustsp": dustsp,
    "dustsp_2nd": "amnt",
    "dustrate_2nd": fert_dict[fertLevel],
    "dustrate": 2 * 100,  # [g m-2 yr = t/ha/yr * 100]
    "dustrad_um": 100,
    "add_secondary": False,
    "imix": 3, # mixing style (1=fickian; 2=homogeneous; 3=tilling)
    "singlerun_seasonality": False,
    "include_psd_full": False,   # False,
    "include_psd_bulk": False,
    'poro_iter_field': False,      # [default='false'] porosity iteration
    'poro_evol': False,            # [default='false'] porosity evolves with solid phase dissolution
    "cec_adsorption_on": True, # True,
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
    'cec_update_from_spinup': True,   # [True, False] whether to update CEC and alpha vars relative to the spinup value (False means no change to cec.in is made)
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
# sitenames = ["albany", "atlanta", "central_valley", "minneapolis"]
# sites = [f"{name}_{time_res}" for name in sitenames]

# by_site = {   # values must have same order as 'sites' var
#     "spinrun": ["site_311a_pr9_spintuneup4"]*len(sites),
#     "climatefiles": sites,  # these serve as the site name when there is no cliamte file to use
# }

# -----
## UNCOMMENT BELOW FOR SINGLE SITE:
sites = ['site_311']
by_site = {   # values must have same order as 'sites' var
    # "cec": [21.10329],
    # "spinrun": ["site_311_pr9_cecON_5m_spintuneup4"], # ["spinupLime_cc_v0_site_311a_app_60p0_psize_100_cc_tau100p0"], # ["site_311a_pr9_YoshiCEC_spintuneup4"], # ["spinupLime_cc_v2_site_311a_app_150p0_psize_100_cc_tau100p0"] 
    "climatefiles": ["site_311"]  # these serve as the site name when there is no cliamte file to use
}


# %% 
# [3] vars to vary within site (we'll get every combination of these two)
ztot_field = [0.25, 0.5, 1., 1.5, 2., 3., 4., 5.,]  # must match spinup depths
cec = [5, 10, 15, 20, 25]

all_combinations = {
    "ztot_field": ztot_field,   # [m] 
    'cec': cec  # [diameter, microns] i think this gets applied to gbas and amnt equally (though amnt is fast-reacting so maybe not a big deal? )
}


# %% 
# --- BUILD DATAFRAME AND SAVE
df = bhf.build_df(pref, const_dict, sites, by_site, all_combinations, add_ctrl=True)

# --- add the particle size distribution inputs if relevant 
# 
if const_dict['include_psd_full'] or const_dict['include_psd_bulk']:
    df['psdrain_meanRad'] = df['dustrad'] / 1e6  # (convert um to m)

# --- add the spinruns based on ztot_field
df['spinrun'] = df.apply(lambda row: f"{row['site']}_pr9_cecON_{row['ztot_field']}m", axis=1)

df

# %% 
# save 
bhf.save_df(df, savepath_batch, fn, multi_run_split, max_iters_per_set)
# %%
fn
# %%
