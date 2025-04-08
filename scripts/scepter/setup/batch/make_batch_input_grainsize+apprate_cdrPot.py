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
import os
import numpy as np
import pandas as pd
import itertools
import batch_helperFxns as bhf


# ---------------------------------------------------
# This script is meant to generate application rates for gbas 
# and cc that yield the same CDR potential for the different feedstocks. 
# I include notes on calculating the CDR potential for each here. 
#
# potential constants relative to cations:
# CO2:Ca2+ mass   -->    2.196       # 44.009x2 g CO2/mol over 40.078 g Ca/mol
# CO2:Mg2+ mass   -->    3.621       # 44.009x2 g CO2/mol over 24.305 g Mg/mol
# CO2:Na+ mass    -->    1.914       # 44.009 g CO2/mol over 22.990 g Na/mol
# CO2:K+ mass     -->    1.126       # 44.009 g CO2/mol over 39.098 g K/mol
#
# potential constants relative to cation oxides:
# CO2:CaO mass    -->    1.570       # 44.009x2 g CO2/mol over 56.079 g CaO/mol
# CO2:MgO mass    -->    2.184       # 44.009x2 g CO2/mol over 40.304 g MgO/mol
# CO2:Na2O mass   -->    1.420       # 44.009x2 g CO2/mol over 61.979 g Na2O/mol  # (CO2 mwt is multiplied by 2 to account for 2 moles of Na in the oxide)
# CO2:K2O mass    -->    0.934       # 44.009x2 g CO2/mol over 94.195 g K2O/mol   # (CO2 mwt is multiplied by 2 to account for 2 moles of K in the oxide)
# 
# 
# cdr potential, calcite:
# molar mass of calcite: 100.0869
# CO2 from Ca2+   -->    44.01       # 40.078 * 2.196 / 2 (accounting for half efficiency due to calcite-carbon)
# ----
# Total           -->    44.01
# CDR potential   --> [[ 0.440 ]]    # total g CO2 / g cations over molar mass 
# 
# cdr potential, basalt:
# molar mass of basalt: 120.496
# CO2 from CaO    -->    21.470      # 13.675 g CaO/mol basalt * 1.570 g CO2 / g CaO
# CO2 from MgO    -->    23.956      # 10.969 g MgO/mol basalt * 2.184 g CO2 / g MgO
# CO2 from Na2O   -->     3.560      # 2.507 g Na2O/mol basalt * 1.420 g CO2 / g Na2O
# CO2 from K2O    -->     0.371      # 0.397 g K2O/mol basalt * 0.934 g CO2 / g K2O
# ----
# Total           -->    49.357 
# CDR potential   --> [[ 0.410 ]]    # total g CO2 / g cations over molar mass  
# 
# *******************************************************
# Notes on the basalt cation oxide molar masses (from 
# comments in cflx_proc.py) : 
# 
    # Note: gbas molar mass depends on whether Dmod_bas_cmp is defined in the makefile! 
    # if the makefile has `CPFLAGS       += -Dmod_basalt_cmp`, then scepter.f90
    # uses basalt_defines.h by default (see scepter.f90 lines ~81, 82). Otherwise, 
    # another basalt composition (hard-coded in scepter.f90) is used.  
    # 
    # This flag was set in my makefile, so I use the molar mass we compute from 
    # basalt_defines.h. 
    # 
    # --- molar mass is calculated following line ~358 in scepter.f90: 
    # (note the terms divided by 2 account for the fact that the cation oxide
    #  has 2 moles of the cation, not one)
    # mwtgbas = (fr_si_gbas*mwtamsi + fr_al_gbas/2*mwtal2o3 + fr_na_gbas/2*mwtna2o
    #            + fr_k_gbas/2*mwtk2o + fr_ca_gbas*mwtcao + fr_mg_gbas*mwtmgo
    #            + fr_fe2_gbas*mwtfe2o)
    # where mwt* is molar mass of species, and fr_*_gbas is the fraction from basalt_defines.h
    #
    # si: 1.0 * 60.085                = 60.085
    # al: 0.4683117231/2 * 101.962    = 23.875
    # na: 0.08088553318/2 * 61.979    = 2.5066
    # k:  0.008431137573/2 * 94.195   = 0.3971
    # ca: 0.2438545566 * 56.079       = 13.675
    # mg: 0.2721686927 * 40.304       = 10.969
    # fe2:0.1251095225 * 71.846       = 8.9886
    # 
    # SUM --------------------------> = 120.496

# %% 
# --- SET FEEDSTOCK POTENTIALS
fs_pots = {    # [kg CO2 / kg rock]
    "gbas": 0.41,   
    "cc":   0.44,
}

# %% 
# --- USER INPUTS
# [1] vars to update, constant for all runs
fertLevel = "low"    # name for how much fertilizer is applied
dustsp = "cc"      # name for dust species to apply (must be from accepted list)
pref = f"{fertLevel}Fert_{dustsp}_cdrpot"
clim_tag = None   # [string] year-span (e.g., 1950-2020) for clim input if climate files are used
                  # (if clim files are not used, set to None)
# save vars
file_prefix = f"meanAnn_{dustsp}_shortRun_cdrpot_{fertLevel}Fert_gs+apprate"  # prefix of output run names
fn = file_prefix + "_v0.csv"
savepath_batch = "/home/tykukla/aglime-swap-cdr/scepter/batch-inputs"
multi_run_split = False   # whether to split the csv into multiple files
max_iters_per_set = 20    # [int] the number of runs per csv (only used if multi_run_split is True)

const_dict = {
    "duration": 15,  # duration of run (starts from earliest year)
    "dustsp": dustsp,
    "dustsp_2nd": "amnt",
    "dustrate_2nd": 6.0, # 30.0,
    "add_secondary": False,
    "imix": 1,
    "singlerun_seasonality": False,
    "include_psd_full": False,
    "include_psd_bulk": False,
    "climatedir": "NA"
}

# %% 
# [2] vars to vary by site
sites = ['site_311a', 'site_311b']
by_site = {   # values must have same order as 'sites' var
    "cec": [21.10329, 6.96125],
    "spinrun": ["site_311a_pr9_spintuneup4", "site_311b_pr9_spintuneup4"],
    "climatefiles": ["site_311a", "site_311b"]  # these serve as the site name when there is no cliamte file to use
}

# %% 
# [3] vars to vary within site (we'll get every combination of these two)
dustrate_ton_ha_yr_CO2pot = [0.1, 0.3, 0.5, 1, 2.5, 4, 5.5, 7, 11, 17, 22, 30, 50]
dustrate_ton_ha_yr = [np.round((1 / fs_pots[dustsp]) * x, 4) for x in dustrate_ton_ha_yr_CO2pot]
all_combinations = {
    "dustrate": [x * 100 for x in dustrate_ton_ha_yr],  # [ton ha-1 yr-1 * 100 = g m-2]   
    "dustrad": [1, 10, 30, 50, 75, 100, 125, 150, 200]  # [diameter, microns] i think this gets applied to gbas and amnt equally (though amnt is fast-reacting so maybe not a big deal? )
}


# %% 
# --- BUILD DATAFRAME AND SAVE
df = bhf.build_df(pref, const_dict, sites, by_site, all_combinations, add_ctrl=True)
df
# %% 
# save 
bhf.save_df(df, savepath_batch, fn, multi_run_split, max_iters_per_set)
# %%
