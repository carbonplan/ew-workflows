# ---
# functions to compute cdr
# 
# --- 
# %%
from scipy.integrate import cumulative_trapezoid
from typing import Tuple
from typing import Callable
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import pickle
import math
import os
import re
import glob
import fnmatch
import time

# %% 
# file_path = "/home/tykukla/SCEPTER/scepter_output/hiFert_gbas__site_311a_app_10000p0_psize_200_gbas_field_tau15p0"
# subdir = "postproc_flxs"
# co2_fn = "rockflx_gbas.pkl" # "carbAlk_flxs.pkl"

# tx = pd.read_pickle(os.path.join(file_path, subdir, co2_fn))
# tx

# %%
# SCEPTER/scepter_output/hiFert_gbas__site_311a_app_10000p0_psize_200_gbas_field_tau15p0/postproc_flxs/co2_flxs.pkl

def _read_pickle_retry(
    path: str,
    max_retries: int = 6,
    base_delay: float = 15.0,
) -> pd.DataFrame:
    """
    pd.read_pickle wrapped with exponential backoff for transient S3
    errors (e.g. 503 Service Unavailable), which are usually gone within
    a minute or two. A genuinely missing key (FileNotFoundError) is not
    retried -- it's raised immediately so callers can treat it as "not found".
    """
    for attempt in range(max_retries):
        try:
            return pd.read_pickle(path)
        except FileNotFoundError:
            raise
        except OSError as e:
            if attempt == max_retries - 1:
                raise
            delay = base_delay * (2 ** attempt)
            print(f"  [retry] transient error reading {path} "
                  f"(attempt {attempt + 1}/{max_retries}): {e}. retrying in {delay:.0f}s...")
            time.sleep(delay)


def read_postproc_flux(
    dfin: pd.DataFrame,  
    outdir: str,
    calc_list: list,
    dfin_cols_to_keep: list,
    flx_type: str='int_flx',
    rockdiss_feedstock: str = None,
    subdir: str = "postproc_flxs",
) -> dict:
    """
    Loop through all SCEPTER output dirs in the dfin dataframe
    and read in the .pkl files based on the calc_list. Combine
    all like .pkl files into one and make a dict where each 
    element is a new .pkl file. 

    Parameters
    ----------
    dfin : pd.DataFrame
        batch dataframe used to find and loop through all 
        the casenames
    outdir : str
        location of the output dirs where each casename in 
        dfin is the name of an output dir
    calc_list : list
        list of cdr calculations to generate (which determines
        the .pkl files that need to be read in)
    dfin_cols_to_keep : list
        list of columns from dfin that we want to add to the final dfs
        (generally the groups / dims for later xr datasets)
    flx_type: str
        ['int_flx' | 'flx'] pull out the integrated flux or 
        the transient flux, respectively
    rockdiss_feedstock : str
        the feedstock to use for the rock dissolution files
        (not used if rockdiss is absent from calc_list)
    subdir : str
        name of the directory within file_path that holds the
        .pkl files

    Returns
    -------
    dict
        dictionary where keys are the calcs in calc_list
        and values are the respective collated pd.DataFrames
    """
    # --- create an empty dictionary to hold the results
    outdict = {}

    #
    # iterate through options in the calc_list
    #

    # --- OPTION 1: "co2_flx" ------------------------
    if "co2_flx" in calc_list:
        print("solving co2_flx")
        outdict['co2_flx'] = co2_flx_in(
            dfin, outdir, dfin_cols_to_keep, 
            rockdiss_feedstock, flx_type=flx_type,
        )

    # --- OPTION 2: "camg_flx" -----------------------
    if "camg_flx" in calc_list:
        print("solving camg_flx")
        outdict['camg_flx'] = cation_flx_in(
            dfin, outdir, dfin_cols_to_keep, rockdiss_feedstock,
            cations_to_track = ['ca', 'mg'], flx_type=flx_type,
        )
    
    # --- OPTION 3: "totcat_flx" ---------------------
    if "totcat_flx" in calc_list:
        print("solving totcat_flx")
        outdict['totcat_flx'] = cation_flx_in(
            dfin, outdir, dfin_cols_to_keep, rockdiss_feedstock,
            cations_to_track = ['ca', 'mg', 'na', 'k'], flx_type=flx_type,
        )
    
    # --- OPTION 4: "carbalk_flx" --------------------
    if "carbalk_flx" in calc_list:
        print("solving carbalk_flx")
        outdict['carbalk_flx'] = carbalk_flx_in(dfin, outdir, dfin_cols_to_keep)
    
    # --- OPTION 5: "rockdiss" -----------------------
    if "rockdiss" in calc_list:
        print("solving rockdiss")
        outdict['rockdiss'] = rockdiss_in(dfin, outdir, dfin_cols_to_keep, rockdiss_feedstock)
    
    #
    # update dustrate_ton_ha_yr if there are rockdiss files that exist
    # this is to ensure accurate values for single / variable application
    # runs (rather than assuming the first application applies over all time)
    # 


    # 
    # return the result
    # 
    return outdict


def get_timemean_dustrate(
    outdir: str,
    tdf: pd.DataFrame,
    dustsp: str,
    subdir: str = "postproc_flxs",
    fn_pref: str = "rockflx",
    roundto: int = 3,
    return_df: bool = False,
) -> float:
    '''
    Compute the time mean dustrate from the appropriate rockdiss file
    and return.

    Parameters
    ----------
    outdir : str
        location of the output dirs where each casename in 
        dfin is the name of an output dir
    tdf : str
        single row in the pandas dataframe associated with this run
    dustsp : str
        the feedstock to use for the rock dissolution files
        (not used if rockdiss is absent from calc_list)
    subdir : str
        name of the directory within file_path that holds the
        .pkl files
    fn_pref : str
        prefix for the name of the rockdiss .pkl file
    roundto : int
        number of digits to round the dust value to (helps 
        avoid floating point errors which is important since
        these values must align to a grid)
    return_df : bool
        [default=False] whether to return the rock dataframe along 
        with the calculated timemean dustrate. 
    Returns
    -------
    float
        value for dustrate in ton_ha_yr averaged over run
    '''
    # --- craft the path to the rockdiss file
    this_fn = f"{fn_pref}_{dustsp}.pkl"
    this_path = os.path.join(outdir, tdf['newrun_id_full'], subdir)
    fn_path = os.path.join(this_path, this_fn)

    # read in the pickle if the path exists
    # (for s3 paths, read directly with retry rather than pre-checking
    # existence with a separate HeadObject call)
    rockdf = None
    if fn_path.startswith("s3://"): # then bring it in from s3
        try:
            rockdf = _read_pickle_retry(fn_path)
        except FileNotFoundError:
            rockdf = None
    elif os.path.exists(fn_path):
        rockdf = pd.read_pickle(fn_path)

    if rockdf is not None:
        # compute the time averaged dust (integrated flux divided by time)
        # (note, integrated flux name is 'int_dust_ton_ha_yr' but units should really
        #  be 'ton_ha', which is why we must divide by year)
        dur = np.max(rockdf['time'])
        dusttot = rockdf[rockdf['time'] == dur]['int_dust_ton_ha_yr'].values[0]
        dustrate_mean = round(dusttot / dur, roundto)
    else:
        dustrate_mean = None

    # return result
    if return_df:
        return dustrate_mean, rockdf
    else: 
        return dustrate_mean



def co2_flx_in(dfin: pd.DataFrame,  
               outdir: str,
               dfin_cols_to_keep: list,
               rockdiss_feedstock: str=None,
               flx_type: str = "int_flx",
               subdir: str = "postproc_flxs",
               co2_fn: str = "co2_flxs.pkl",
) -> pd.DataFrame:
    """
    Read in the co2_flxs.pkl file for a all run
    directories in dfin. Separate out just the time-
    integrated fluxe or transient depending on value
    of flx_type. return as a pandas df.

    Parameters
    ----------
    dfin : pd.DataFrame
        batch dataframe used to find and loop through all 
        the casenames
    outdir : str
        location of the output dirs where each casename in 
        dfin is the name of an output dir
    dfin_cols_to_keep : list
        list of columns from dfin that we want to add to the final dfs
        (generally the groups / dims for later xr datasets
    rockdiss_feedstock : str
        the feedstock to use for the rock dissolution files
        (not used if rockdiss is absent from calc_list)
    flx_type: str
        ['int_flx' | 'flx'] pull out the integrated flux or 
        the transient flux, respectively
    subdir : str
        name of the directory within file_path that holds the
        .pkl files
    co2_fn : str
        name of the co2 flx .pkl file

    Returns
    -------
    pd.DataFrame
        single dataframe with cdr fluxes for the given SCEPTER run
    """
    # --- loop through runs
    # track run index
    rundx = 0
    outdf_exist = False
    # loop
    for run in range(len(dfin)):
        tdf = dfin.iloc[run]
        this_path = os.path.join(outdir, tdf['newrun_id_full'], subdir)
        fn_path = os.path.join(this_path, co2_fn)

        # ----------------------------
        # read in the pickle if the path exists
        # (for s3 paths, read directly rather than pre-checking existence
        # with a separate HeadObject call -- pd.read_pickle raises
        # FileNotFoundError for a genuinely missing key)
        tmpdf = None
        if fn_path.startswith("s3://"): # then bring it in from s3
            try:
                tmpdf = _read_pickle_retry(fn_path)
            except FileNotFoundError:
                tmpdf = None
        elif os.path.exists(fn_path):
            tmpdf = pd.read_pickle(fn_path)

        if tmpdf is not None:
            # pull out just the time-integrated data and reset row indices
            tmpdf = tmpdf[tmpdf['flx_type'] == flx_type].reset_index(drop=True)

            # add whether or not it's the control run 
            tmpdf["ctrl"] = tdf["ctrl_run"]

            # add the other columns to keep
            for col in dfin_cols_to_keep:
                tmpdf[col] = tdf[col]

            # update the dustrate if necessary 
            if "dustrate_ton_ha_yr" in dfin_cols_to_keep:
                ds_mean = get_timemean_dustrate(outdir, tdf, rockdiss_feedstock)
                if ds_mean is not None: # then a new dustrate was able to be computed
                    tmpdf["dustrate_ton_ha_yr"] = ds_mean
            
            # bring together
            if (rundx == 0) | (outdf_exist == False):
                outdf = tmpdf.copy()
                outdf_exist = True
            else:
                outdf = pd.concat([outdf.copy(), tmpdf.copy()], ignore_index=True)

        else: 
            print(f"could not find {tdf['newrun_id_full']} -- run {run}")
        # move to the next index
        rundx += 1
     
    # return result
    try:
        return outdf
    except:
        raise ValueError("No co2 flux dataframe was created. All run dirs may be empty, or perhaps the outdir or subdir is incorrect?")



def cation_flx_in(dfin: pd.DataFrame,  
                  outdir: str,
                  dfin_cols_to_keep: list,
                  rockdiss_feedstock: str=None,
                  cations_to_track: list = ['ca', 'mg'],
                  flx_type: str = "int_flx",
                  subdir: str = "postproc_flxs",
                  cation_fn_prefix: str = "cationflx_",
) -> pd.DataFrame:
    """
    Read in the cationflx_*.pkls for calcium and magnesium 
    for all run directories in dfin. Separate out just the time-
    integrated fluxes. return as a pandas df.

    Parameters
    ----------
    dfin : pd.DataFrame
        batch dataframe used to find and loop through all 
        the casenames
    outdir : str
        location of the output dirs where each casename in 
        dfin is the name of an output dir
    dfin_cols_to_keep : list
        list of columns from dfin that we want to add to the final dfs
        (generally the groups / dims for later xr datasets)
    rockdiss_feedstock : str
        the feedstock to use for the rock dissolution files
        (not used if rockdiss is absent from calc_list)
    cations_to_track : list
        list of the cations that we'll read in 
    flx_type: str
        ['int_flx' | 'flx'] pull out the integrated flux or 
        the transient flux, respectively
    subdir : str
        name of the directory within file_path that holds the
        .pkl files
    cation_fn_prefix : str
        prefix for the individual cation .pkl files (e.g., 
        'cationflx_' is the prefix for 'cationflx_ca.pkl')

    Returns
    -------
    pd.DataFrame
        single dataframe with cdr fluxes for the given SCEPTER run
    """
    # decide what the cdrpot columns are (depending on the flx_type)
    # (these get summed across all cations at the end)
    # (note, adv_charge doesn't exist yet, but we make it before cdrpot_cols is called)
    if flx_type == "int_flx":
        cdrpot_cols = ['co2pot_tot_tonHa', 'co2pot_adv_tonHa', 'adv_charge_DICpot']
    elif flx_type == "flx":
        cdrpot_cols = ['co2pot_tot_tonHaYr', 'co2pot_adv_tonHaYr', 'adv_charge_DICpot']

    # decide which columns to keep from the total df
    # since we end up combining the individual cation dfs
    # together
    cols_to_keep = ['time', 'units', 'runname', 'noncarbsld_source', 
                    'carbsld_source', 'flx_type']
    cols_to_keep.extend(cdrpot_cols)
    # columns that we need to append the cation name to 
    append_cation = ['noncarbsld_source', 'carbsld_source']


    # --- loop through runs
    # track run index
    rundx = 0
    outdf_exist = False
    # loop
    for run in range(len(dfin)):
        tdf = dfin.iloc[run]
        this_path = os.path.join(outdir, tdf['newrun_id_full'], subdir)
        # loop through the different cations and read each one directly
        # (for s3 paths, read with retry rather than pre-checking existence
        # with a separate HeadObject call; a missing file is treated as
        # "not found" and skips the whole run, same as before)
        cation_dfs = []
        missing = False
        for cat in cations_to_track:
            fn_thiscat = f"{cation_fn_prefix}{cat}.pkl"
            fn_path = os.path.join(this_path, fn_thiscat)
            try:
                if fn_path.startswith("s3://"): # then bring it in from s3
                    cation_dfs.append(_read_pickle_retry(fn_path))
                elif os.path.exists(fn_path):
                    cation_dfs.append(pd.read_pickle(fn_path))
                else:
                    missing = True
                    break
            except FileNotFoundError:
                missing = True
                break

        # if any cation file is missing, we just move on to the next
        # run (no results are returned for this one)
        if missing:
            continue

        # otherwise loop through each cation and add them together
        cat_dx = 0   # index of the cation to update in loop
        for tmpdf in cation_dfs:
            # pull out just the time-integrated data and reset row indices
            tmpdf = tmpdf[tmpdf['flx_type'] == flx_type].reset_index(drop=True).copy()
            # get the cation
            tmpcat = tmpdf['cation'][0]
            # get the advected flux multiplied by charge and converted to co2 potential
            #  (we need this to sum together at the end for the downstream loss calculation later)
            g_mol = 44.01
            ton_g = 1 / 1e6   # [ton g-1]
            m2_ha = 10e3      # [m2 ha-1]
            conv_factor = g_mol * ton_g * m2_ha 
            tmpdf['adv_charge_DICpot'] = tmpdf['adv'].copy() * tmpdf['charge'].copy() * conv_factor
            
            # extract just the columns we want to keep
            tmpdf = tmpdf[cols_to_keep]

            # multiply by time if flx_type is int (note, for the other
            # .pkl files we multiplied by time in the `cflx_proc.py` processing
            # step. it's just the cations where we need to multiply by time)
            if flx_type == "int_flx":
                dont_multiply_by_time = ['time', 'flx_type', 'units', 'runname', 'cation', 'charge']
                tmpdf = tmpdf.apply(lambda x: x * tmpdf['time'] if x.name not in dont_multiply_by_time else x).copy()

            # append the cation to the relevant columns
            tmpdf.columns = [f'{col}_{tmpcat}' if col in append_cation else col for col in tmpdf.columns]

            # bring together
            if cat_dx == 0:   # use the same tmpdf if it's the first cation
                outdf_site = tmpdf.copy()
            else:
                # add co2 fluxes
                for col in cdrpot_cols:
                    outdf_site[col] = outdf_site.copy()[col] + tmpdf.copy()[col]
                # add the other columns
                for col in append_cation:
                    outdf_site[f'{col}_{tmpcat}'] = tmpdf[f'{col}_{tmpcat}']
            
            cat_dx += 1

        # add whether or not it's the control run 
        outdf_site["ctrl"] = tdf["ctrl_run"]

        # add the other columns to keep
        for col in dfin_cols_to_keep:
            outdf_site[col] = tdf[col]

        # update the dustrate if necessary 
        if "dustrate_ton_ha_yr" in dfin_cols_to_keep:
            ds_mean = get_timemean_dustrate(outdir, tdf, rockdiss_feedstock)
            if ds_mean is not None: # then a new dustrate was able to be computed
                outdf_site["dustrate_ton_ha_yr"] = ds_mean
    
        # bring together
        if (rundx == 0) | (outdf_exist == False):
                outdf = outdf_site.copy()
                outdf_exist = True
        else:
            outdf = pd.concat([outdf.copy(), outdf_site.copy()], ignore_index=True)

        # move to the next index
        rundx += 1
    
    # return result
    try:
        return outdf
    except:
        raise ValueError("No cation flux dataframe was created. All run dirs may be empty, or perhaps the outdir or subdir is incorrect?")


def rockdiss_in(dfin: pd.DataFrame,  
               outdir: str,
               dfin_cols_to_keep: list,
               feedstock: str,
               subdir: str = "postproc_flxs",
               rock_prefix: str = "rockflx_",
               roundto: int=3,
               ) -> pd.DataFrame:
    """
    Read in the roclflx_*.pkl file for a all run
    directories in dfin. Separate out just the time-
    integrated flux or transient depending on value
    of flx_type. return as a pandas df.

    Parameters
    ----------
    dfin : pd.DataFrame
        batch dataframe used to find and loop through all 
        the casenames
    outdir : str
        location of the output dirs where each casename in 
        dfin is the name of an output dir
    dfin_cols_to_keep : list
        list of columns from dfin that we want to add to the final dfs
        (generally the groups / dims for later xr datasets)
    feedstock : str or None
        name of the feedstock to pull from (e.g., 'cc', 'gbas', 'amnt', etc)
        (can usually find this in dfin) if None, then use "dustsp" in dfin. 
    subdir : str
        name of the directory within file_path that holds the
        .pkl files
    rock_prefix : str
        prefix of the .pkl file (e.g., "rockflx_" for "rockflx_gbas.pkl")
    roundto : int
        number of digits to round the dust value to (helps 
        avoid floating point errors which is important since
        these values must align to a grid)

    Returns
    -------
    pd.DataFrame
        single dataframe with cdr fluxes for the given SCEPTER run
    """
    # --- loop through runs
    # track run index
    rundx = 0
    outdf_exist = False
    # loop
    for run in range(len(dfin)):
        tdf = dfin.iloc[run]
        # find feedstock dustsp if not defined
        if feedstock is None:
            if "dustsp" not in dfin.columns:
                raise ValueError("Feedstock not defined and dustsp is not in `dfin` -- can't tell which feedstock to use! ")
            feedstock = dfin['dustsp'].iloc[run]
        this_path = os.path.join(outdir, tdf['newrun_id_full'], subdir)
        tmpfn = f'{rock_prefix}{feedstock}.pkl'
        fn_path = os.path.join(this_path, tmpfn)

        # read in the pickle if the path exists
        # (for s3 paths, read directly with retry rather than pre-checking
        # existence with a separate HeadObject call)
        tmpdf = None
        if fn_path.startswith("s3://"): # then bring it in from s3
            try:
                tmpdf = _read_pickle_retry(fn_path)
            except FileNotFoundError:
                tmpdf = None
        elif os.path.exists(fn_path):
            tmpdf = pd.read_pickle(fn_path)

        if tmpdf is not None:
            # add whether or not it's the control run
            tmpdf["ctrl"] = tdf["ctrl_run"]

            # add the runname since we forgot that earlier
            tmpdf['runname'] = tdf['newrun_id_full']

            # add the other columns to keep
            for col in dfin_cols_to_keep:
                tmpdf[col] = tdf[col]

            # update the dustrate column if necessary 
            if "dustrate_ton_ha_yr" in dfin_cols_to_keep:
                dur = np.max(tmpdf['time'])
                dusttot = tmpdf[tmpdf['time'] == dur]['int_dust_ton_ha_yr'].values[0]
                dustrate_mean = round(dusttot / dur, roundto)
                tmpdf['dustrate_ton_ha_yr'] = dustrate_mean

            # bring together
            if (rundx == 0) | (outdf_exist == False):
                outdf = tmpdf.copy()
                outdf_exist = True
            else:
                outdf = pd.concat([outdf.copy(), tmpdf.copy()], ignore_index=True)

        # move to the next index
        rundx += 1
    
    # return result
    try:
        return outdf
    except:
        raise ValueError("No rockdiss dataframe was created. All run dirs may be empty, or perhaps the outdir or subdir is incorrect?")



# ---- ! 
#      ! WIP because the original cflx code over-wrote 
#      ! this file with the sum of cations file 
#      !
def carbalk_flx_in(dfin: pd.DataFrame,  
                   outdir: str,
                   dfin_cols_to_keep: list,
                   flx_type: str = "int_flx",
                   subdir: str = "postproc_flxs",
                   carbAlk_fn: str = "carbAlk_flxs.pkl",
                   ) -> pd.DataFrame:
    """
    Read in the carbalk.pkl file for a all run
    directories in dfin. Separate out just the time-
    integrated fluxe or transient depending on value
    of flx_type. return as a pandas df.

    Parameters
    ----------
    dfin : pd.DataFrame
        batch dataframe used to find and loop through all 
        the casenames
    outdir : str
        location of the output dirs where each casename in 
        dfin is the name of an output dir
    dfin_cols_to_keep : list
        list of columns from dfin that we want to add to the final dfs
        (generally the groups / dims for later xr datasets)
    flx_type: str
        ['int_flx' | 'flx'] pull out the integrated flux or 
        the transient flux, respectively
    subdir : str
        name of the directory within file_path that holds the
        .pkl files
    carbAlk_fn : str
        name of the carbonate alkalinity flx .pkl file

    Returns
    -------
    pd.DataFrame
        single dataframe with carbalk fluxes for the given SCEPTER run
    """
    return "WIP"  



def cdr_int_per_group(
    flx_dict: dict,
    time_horizon: float,
    calc_list: list,
    dfin_cols_to_keep: list,
    bysite: bool=False,
    ctrl_params: list=None,
) -> Tuple[dict, dict]:
    """
    Compute CDR (or difference from control) for each of the methods
    listed in the calc_list. Return a similarly-structured dictionary.

    flx_dict should be the dict from read_postproc_flux function. 
    Only use this function for integrated fluxes.

    Parameters
    ----------
    flx_dict : dict
        dictionary where each key is an element in calc_list, and 
        values are pd.Dataframes output from the `read_postproc_flux`
        function.
    time_horizon : float
        number of years over which to integrate CDR 
    calc_list : list
        list of cdr calculations to generate (which determines
        the .pkl files that need to be read in)
    dfin_cols_to_keep : list
        list of columns from dfin that we want to keep in the output files
    bysite : bool
        [True] if dfin includes more than one site, otherwise [False]
    
    Returns 
    -------
    dict 
        [1] dictionary with the same structure as the input, but with 
        CDR calcs over time
        [2] dictionary with the same structure as the input, but with 
        CDR calcs integrated to the time horizon
    """
    # --- create an empty dictionary to hold the results
    outdict_full = {}
    outdict_sum = {}

    #
    # iterate through options in the calc_list
    #

    # --- OPTION 1: "co2_flx" ------------------------
    if "co2_flx" in calc_list:
        tc = 'co2_flx'
        print(f"solving {tc}")
        dfin = flx_dict[tc]
        outdict_full[tc], outdict_sum[tc] = co2_flx_cdr(dfin, time_horizon, dfin_cols_to_keep, bysite=bysite, ctrl_params=ctrl_params)

    # --- OPTION 2: "camg_flx" -----------------------
    if "camg_flx" in calc_list:
        tc = 'camg_flx'
        print(f"solving {tc}")
        dfin = flx_dict[tc]
        # (note, we use the same cation_flx_cdr fxn for both camg and totcat
        #  because it just calculates based on whatever cats are provided)
        outdict_full[tc], outdict_sum[tc] = cation_flx_cdr(dfin, time_horizon, dfin_cols_to_keep, cation_tag=tc, bysite=bysite, ctrl_params=ctrl_params)
    
    # --- OPTION 3: "totcat_flx" ---------------------
    if "totcat_flx" in calc_list:
        tc = 'totcat_flx'
        print(f"solving {tc}")
        dfin = flx_dict[tc]
        # (note, we use the same cation_flx_cdr fxn for both camg and totcat
        #  because it just calculates based on whatever cats are provided)
        outdict_full[tc], outdict_sum[tc] = cation_flx_cdr(dfin, time_horizon, dfin_cols_to_keep, cation_tag=tc, bysite=bysite, ctrl_params=ctrl_params)
    
    # --- OPTION 4: "carbalk_flx" --------------------
    if "carbalk_flx" in calc_list:
        print("solving carbalk_flx")
        outdict_full['carbalk_flx'] = carbalk_flx_cdr()
    
    # --- OPTION 5: "rockdiss" -----------------------
    if "rockdiss" in calc_list:
        tc = 'rockdiss'
        print(f"solving {tc}")
        dfin = flx_dict[tc]
        outdict_full[tc], outdict_sum[tc] = rockdiss_synth(dfin, time_horizon, dfin_cols_to_keep, bysite=bysite)

    # --- OPTION 6: "catbudget_flx" ------------------------------------------
    if "catbudget_flx" in calc_list:
        tc = 'catbudget_flx'
        print(f"solving {tc}")
        dfin_cb = flx_dict[tc]
        outdict_full[tc], outdict_sum[tc] = catbudget_cdr(
            dfin_cb, time_horizon, dfin_cols_to_keep,
            bysite=bysite, ctrl_params=ctrl_params,
        )
    # 
    # return the result
    # 
    return outdict_full, outdict_sum


def co2_flx_cdr(
    dfin: pd.DataFrame, 
    time_horizon: float,
    dfin_cols_to_keep: list,
    bysite: bool = False,
    ctrl_params: list = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read in the Pandas DataFrame for the co2flx data 
    generated by the `read_postproc_flux` function and 
    compute the different forms of CDR from it

    Parameters
    ----------
    dfin : pd.DataFrame
        the pd.DataFrame with all flux data from `read_postproc_flux` function
    time_horizon : float
        number of years over which to integrate CDR 
    dfin_cols_to_keep : list
        list of columns from dfin that we want to keep in the output files
    bysite : bool
        [True] if dfin includes more than one site, otherwise [False]
    ctrl_params : list
        OVERWRITES `bysite`! If multiple parameters are needed to find the control case (e.g., if you vary
        across sites and climate resolution, you might want a site-climRes control), 
        list the columns here that we need to match between the case df `tdf` and the
        ctrl df `tdf_ctrl` (the latter is just the runs with no dustflux, so we need 
        to narrow it down furthar). (it extends the bysite functionality to multiple 
        params.)
    
    Returns
    -------
    pd.DataFrame, pd.DataFrame
        [1] timeseries of CDR
        [2] synthesized CDR calculations for each run 
    """
    # columns to carry over to the output df
    cols_to_keep = ['time', 'units', 'flx_type']
    cols_to_keep.extend(dfin_cols_to_keep)

    # get control run
    df_ctrl = dfin.loc[dfin['ctrl'] == True]
    if not bysite: 
        # only the numeric columns for interpolation
        tdf_ctrl = df_ctrl.select_dtypes(include=[np.number])
    # all cases
    df_case = dfin.loc[dfin['ctrl'] == False]

    # get each case name
    case_names = df_case['runname'].unique()

    # track the case index in the loop
    casedx = 0
    # --- loop through each run
    for case in case_names:
        # pull out this case
        tdf = df_case[df_case['runname'] == case].copy()

        # create the temp output df
        tdfout = tdf[cols_to_keep].copy()

        # get control if going by site
        if bysite and not ctrl_params:
            # find the control case (only the numeric columns for integration)
            tdf_ctrl = df_ctrl[df_ctrl['site'] == tdf['site'].values[0]].select_dtypes(include=[np.number])
        elif ctrl_params:
            # find the match across multiple columns
            mask = np.logical_and.reduce(
                [df_ctrl[c] == tdf[c].iloc[0] for c in ctrl_params]
            )
            tdf_ctrl = df_ctrl.loc[mask].select_dtypes(include=[np.number])

        # if case and control are different lengths, we need to interpolate
        # (this happens sometimes due to shifts in how the timesteps are handled in a given run)
        if len(tdf) != len(tdf_ctrl):
            # print(f'{len(tdf)} --- {len(tdf_ctrl)}')
            tdf_ctrl = tdf_ctrl.set_index('time').reindex(tdf['time']).interpolate(method='linear').reset_index().copy()

        # --- DIFFUSIVE CDR: 
        # compute CDR by diffusive flux minus respiration change
        tdfout['cdr_dif_component'] = -1*(tdf['co2flx_dif'].values - tdf_ctrl['co2flx_dif'].values)
        tdfout['cdr_resp_component'] = -1*(tdf['co2flx_resp'].values - tdf_ctrl['co2flx_resp'].values)
        tdfout['cdr_resp_component_noAE'] = np.minimum(tdfout['cdr_resp_component'], 0)
        tdfout['cdr_dif'] = tdfout['cdr_dif_component'] + tdfout['cdr_resp_component_noAE']

        # --- ADVECTIVE CDR (no column inorg formation):
        # compute CDR by the flux of carbon advected out the bottom of the model domain
        # and don't account for any boost from inorganic carbon minerals forming 
        # within the model domain (we do this by making sure the co2flx_inorg is never
        # positive (which would indicate net precipitation)). Note that both advective 
        # cdr equations ensure we aren't counting fossil carbon as CDR
        tdf['co2_adv_no_soilSIC'] = tdf['co2flx_adv'] + np.minimum(tdf['co2flx_inorg'], 0)
        tdf_ctrl['co2_adv_no_soilSIC'] = tdf_ctrl['co2flx_adv'] + np.minimum(tdf_ctrl['co2flx_inorg'], 0)
        tdfout['cdr_adv'] = tdf['co2_adv_no_soilSIC'].values - tdf_ctrl['co2_adv_no_soilSIC'].values

        # --- ADVECTIVE CDR + NEW SOIL SIC:
        # compute CDR by the flux of carbon advected out the bottom of the model domain
        # plus any new soil SIC that forms. we do this with:
        # adv_noinorg = adv + cc (given adv = –dif – resp – cc – tflx)
        # such that positive values of cc (net SIC formation) contribute to CDR. Note that 
        # both advective cdr equations ensure we aren't counting fossil carbon as CDR.
        tdfout['cdr_adv_plus_newSIC'] = (tdf['co2flx_adv_noinorg'].values - tdf_ctrl['co2flx_adv_noinorg'].values)
        
        # --- SOIL SIC CDR
        # compute the CDR that is solely due to new soil SIC forming. If SIC is net
        # dissolving, then this should be zero. 
        tdfout['cdr_SIConly'] = tdfout['cdr_adv_plus_newSIC'] - tdfout['cdr_adv']

        # --- TOTAL ADVECTED CARBON:
        # compute the total amount of C advected out of the system. This is required for 
        # the downstream loss calculations that we'll perform later
        tdfout['tot_adv'] = (tdf['co2flx_adv'].values - tdf_ctrl['co2flx_adv'].values)

        # multiply columns by time to get the time-integrated CDR (don't multiply cols_to_keep
        # since they're not C fluxes)
        # ( !! they're already multiplied by time in cflx_proc.py !! )
        # XX tdfout = tdfout.apply(lambda x: x * tdfout['time'] if x.name not in cols_to_keep else x)

        # throw away data above the time horizon (add a tiny negligible amount to 
        # time horizon to get the exact time horizon value in the output (avoiding 
        # bit rounding errors))
        tdfout = tdfout.loc[tdfout['time'] <= (time_horizon + 1e-6)].reset_index(drop=True).copy()

        # get just the summary row (represents integrated fluxes at the time horizon)
        tdf_summary = tdfout[tdfout['time'] == tdfout['time'].max()]
        tdf_summary = tdf_summary.drop(columns=['time']) # drop time column from time-integrated df

        # --- create outputs
        if casedx == 0:
            dfout_full = tdfout.copy()
            dfout_sum = tdf_summary.copy()
        else:
            dfout_full = pd.concat([dfout_full.copy(), tdfout.copy()], ignore_index=True)
            dfout_sum = pd.concat([dfout_sum.copy(), tdf_summary.copy()], ignore_index=True)
        
        # update the index 
        casedx += 1

    dfout_full['time_horizon'] = time_horizon
    dfout_sum['time_horizon'] = time_horizon
    dfout_full['cdr_fxn'] = "co2_flx"
    dfout_sum['cdr_fxn'] = "co2_flx"

    # return results
    return dfout_full, dfout_sum


def cation_flx_cdr(
    dfin: pd.DataFrame,
    time_horizon: float,
    dfin_cols_to_keep: list,
    cation_tag: str,
    bysite: bool=False,
    ctrl_params: list=None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read in the Pandas DataFrame for the co2flx data 
    generated by the `read_postproc_flux` function and 
    compute the different forms of CDR from it

    Parameters
    ----------
    dfin : pd.DataFrame
        the pd.DataFrame with all flux data from `read_postproc_flux` function
    time_horizon : float
        number of years over which to integrate CDR 
    dfin_cols_to_keep : list
        list of columns from dfin that we want to keep in the output files
    cation_tag : str
        tag for the 'cdr_fxn' var that's saved in the output. Use
        `totcat_flx` for total cations, or `camg_flx` for just calcium and 
        magnesium
    bysite : bool
        [True] if dfin includes more than one site, otherwise [False]
    ctrl_params : list
        OVERWRITES `bysite`! If multiple parameters are needed to find the control case (e.g., if you vary
        across sites and climate resolution, you might want a site-climRes control), 
        list the columns here that we need to match between the case df `tdf` and the
        ctrl df `tdf_ctrl` (the latter is just the runs with no dustflux, so we need 
        to narrow it down furthar). (it extends the bysite functionality to multiple 
        params.)
    
    Returns
    -------
    pd.DataFrame, pd.DataFrame
        [1] timeseries of CDR
        [2] synthesized CDR calculations for each run 
    """
    # columns to carry over to the output df
    cols_to_keep = ['time', 'units', 'flx_type']
    cols_to_keep.extend(dfin_cols_to_keep)

    # get control runs
    df_ctrl = dfin.loc[dfin['ctrl'] == True]
    if not bysite: 
        # only the numeric columns for interpolation
        tdf_ctrl = df_ctrl.select_dtypes(include=[np.number])
    # all cases
    df_case = dfin.loc[dfin['ctrl'] == False]

    # get each case name
    case_names = df_case['runname'].unique()

    # track the case index in the loop
    casedx = 0
    # --- loop through each run
    for case in case_names:
        # pull out this case
        tdf = df_case[df_case['runname'] == case]

        # create the temp output df
        tdfout = tdf[cols_to_keep].copy()

        # get control if going by site
        if bysite and not ctrl_params:
            # find the control case (only the numeric columns for integration)
            tdf_ctrl = df_ctrl[df_ctrl['site'] == tdf['site'].values[0]].select_dtypes(include=[np.number])
        elif ctrl_params:
            # find the match across multiple columns
            mask = np.logical_and.reduce(
                [df_ctrl[c] == tdf[c].iloc[0] for c in ctrl_params]
            )
            tdf_ctrl = df_ctrl.loc[mask].select_dtypes(include=[np.number])

        # if case and control are different lengths, we need to interpolate
        # (this happens sometimes due to shifts in how the timesteps are handled in a given run)
        if len(tdf) != len(tdf_ctrl):
            tdf_ctrl = tdf_ctrl.set_index('time').reindex(tdf['time']).interpolate(method='linear').reset_index().copy()

        # --- TOTAL CDR:
        # advected + change in column storage 
        # first find the column name (it differs depending on the units)
        col_tot_cdr = [col for col in tdf.columns if col.startswith("co2pot_tot_")]
        tdfout[col_tot_cdr] = tdf[col_tot_cdr].values - tdf_ctrl[col_tot_cdr].values

        # --- ADVECTIVE CDR:
        # advected (change in column storage doesn't count)
        # first find the column name (it differs depending on the units)
        col_adv_cdr = [col for col in tdf.columns if col.startswith("co2pot_adv_")]
        tdfout[col_adv_cdr] = tdf[col_adv_cdr].values - tdf_ctrl[col_adv_cdr].values

        # --- get components 
        # first find the column names (should be one per cation)
        col_noncarb = [col for col in tdf.columns if col.startswith("noncarbsld_source")]
        col_carb = [col for col in tdf.columns if col.startswith("carbsld_source")]
        # loop through non-carbonate sources
        for tcol1 in col_noncarb:
            tdfout[tcol1] = tdf[tcol1].values - tdf_ctrl[tcol1].values
        # loop through carbonate sources
        for tcol2 in col_carb:
            tdfout[tcol2] = tdf[tcol2].values - tdf_ctrl[tcol2].values

        # --- GET FLUX OUT OF THE COLUMN
        tdfout['DICpot_adv'] = tdf['adv_charge_DICpot'].values - tdf_ctrl['adv_charge_DICpot'].values


        # throw away data above the time horizon (add a tiny negligible amount to 
        # time horizon to get the exact time horizon value in the output (avoiding 
        # bit rounding errors))
        tdfout = tdfout.loc[tdfout['time'] <= (time_horizon + 1e-6)].reset_index(drop=True).copy()

        # get just the summary row (represents integrated fluxes at the time horizon)
        tdf_summary = tdfout[tdfout['time'] == tdfout['time'].max()]
        tdf_summary = tdf_summary.drop(columns=['time']) # drop time column from time-integrated df

        # --- create outputs
        if casedx == 0:
            dfout_full = tdfout.copy()
            dfout_sum = tdf_summary.copy()
        else:
            dfout_full = pd.concat([dfout_full.copy(), tdfout.copy()], ignore_index=True)
            dfout_sum = pd.concat([dfout_sum.copy(), tdf_summary.copy()], ignore_index=True)

        # update the index 
        casedx += 1

    dfout_full['time_horizon'] = time_horizon
    dfout_sum['time_horizon'] = time_horizon
    dfout_full['cdr_fxn'] = cation_tag
    dfout_sum['cdr_fxn'] = cation_tag
    
    # return results
    return dfout_full, dfout_sum

    
def carbalk_flx_cdr(   # !! WIP !!
    dfin: pd.DataFrame,
    time_horizon: float,
    dfin_cols_to_keep: list,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # WIP
    return "WIP"


def rockdiss_synth(
    dfin: pd.DataFrame,
    time_horizon: float,
    dfin_cols_to_keep: list,
    bysite: bool=False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read in the Pandas DataFrame for the rockdiss data 
    generated by the `read_postproc_flux` function and 
    synthesize into output dfs

    We're not doing any comparison / difference to the 
    control here because we assume there is no application
    in the control. 

    Parameters
    ----------
    dfin : pd.DataFrame
        the pd.DataFrame with all flux data from `read_postproc_flux` function
    time_horizon : float
        number of years over which to integrate CDR 
    dfin_cols_to_keep : list
        list of columns from dfin that we want to keep in the output files
    bysite : bool
        [True] if dfin includes more than one site, otherwise [False]
    
    Returns
    -------
    pd.DataFrame, pd.DataFrame
        [1] timeseries of CDR
        [2] synthesized CDR calculations for each run 
    """
    # columns to carry over to the output df
    # (we keep basically everything, so it's easier to list
    #  the columns we want to get rid of)
    cols_to_discard = ['runname', 'ctrl']
    # if a column to keep is in the discard pile remove it from discard
    overlap = set(dfin_cols_to_keep) & set(cols_to_discard)
    if overlap: # remove the overlap 
        cols_to_discard = [x for x in cols_to_discard if x not in overlap]

    # get control runs
    df_ctrl = dfin.loc[dfin['ctrl'] == True]
    if not bysite: 
        # only the numeric columns for interpolation
        tdf_ctrl = df_ctrl.select_dtypes(include=[np.number])
    # all cases
    df_case = dfin.loc[dfin['ctrl'] == False]

    # get each case name
    case_names = df_case['runname'].unique()

    # track the case index in the loop
    casedx = 0
    # --- loop through each run
    for case in case_names:
        # pull out this case
        tdf = df_case[df_case['runname'] == case]

        # create the temp output df
        tdfout = tdf.drop(columns=cols_to_discard).copy()

        # get control if going by site
        if bysite:
            # find the control case (only the numeric columns for integration)
            tdf_ctrl = df_ctrl[df_ctrl['site'] == tdf['site'].values[0]].select_dtypes(include=[np.number])

        # we're not doing any comparison to the control here 
        # because we assume the control is no application
        
        # throw away data above the time horizon (add a tiny negligible amount to 
        # time horizon to get the exact time horizon value in the output (avoiding 
        # bit rounding errors))
        tdfout = tdfout.loc[tdfout['time'] <= (time_horizon + 1e-6)].reset_index(drop=True).copy()

        # get just the summary row (represents integrated fluxes at the time horizon)
        tdf_summary = tdfout[tdfout['time'] == tdfout['time'].max()]
        tdf_summary = tdf_summary.drop(columns=['time']) # drop time column from time-integrated df

        # --- create outputs
        if casedx == 0:
            dfout_full = tdfout.copy()
            dfout_sum = tdf_summary.copy()
        else:
            dfout_full = pd.concat([dfout_full.copy(), tdfout.copy()], ignore_index=True)
            dfout_sum = pd.concat([dfout_sum.copy(), tdf_summary.copy()], ignore_index=True)

        # update the index 
        casedx += 1

    dfout_full['time_horizon'] = time_horizon
    dfout_sum['time_horizon'] = time_horizon
    
    # return results
    return dfout_full, dfout_sum


def emissions_calculator_df(
    df : pd.DataFrame, 
    dustrate_name: str,
    p80_input: float, 
    truck_km: float, 
    barge_km: float, 
    barge_diesel_km: float, 
    Efactor_org: str, 
    mineral: str,
    bwi_overwrite: float = None,
    efactor_overwrite: float = None,
) -> pd.DataFrame:
    """
    Read in information about amount and grain size of rock, as well as its transport
    pathway and starting grain size (for crushing) to return the CO2 emissions
    per rock, per ha per year, and per ha (over defined time horizon) as a pandas df

    Uses the constants / model of Zhang et al., 2023 
    (https://doi.org/10.1021/acs.est.3c01658)

    Parameters
    ----------
    df : pd.DataFrame
        pandas dataframe output from the cdr_per_group function. Must have 
        columns named "dustrate_ton_ha_yr", "dustrad", and "timehorizon_yr"
        ("dustrate_*" replaces "rockTons_ha_yr"; "dustrad" replaces "p80_output";
        and "timehorizon_yr" replaces "dur" from the regular emissions_calculator fxn)
    dustrate_name : str
        name of the dustrate column in the dataframe (e.g., "int_dust1_ton_ha_yr" or "dustrate_ton_ha_yr")

    # -- rocks -- # 
    p80_input : float
        [microns] p80 diameter before crushing the rock (see Zhang et al., 2023 table S5 for example data)
    
    # -- transport -- # 
    truck_km : float
        [km] distance traveled by truck in transporting material
    barge_km : float
        [km] distance traveled by barge in transporting material
    barge_diesel_km : float
        [km] distance traveled by diesel barge
    
    bwi_overwrite : float
        if we don't want to use the default bond work index we can overwrite it here.
        only overwritten if it's a value > 0
    
    # -- crushing -- # 
    Efactor_org : str
        ["MRO" | "RFC" | "SERC"] which crushing Efactor to use
    mineral : str
        ["cc" | "gbas" | "wls" None] which mineral is being transported / crushed (this sets the bond work index)
    efactor_overwrite : float
        if we don't want to use the default e factor we can overwrite it here.
        only overwritten if it's a value > 0
    """
    # --- CONSTANTS
    # transport -- Zhang et al., 2023, table S1
    truck_factor = 0.0996    # [kg CO2e / ton / km] 
    barge_factor = 0.0282    # [kg CO2e / ton / km] 
    barge_diesel_factor = 0.00534  # [kg CO2e / ton / km] 

    # crushing -- Zhang et al., 2023, table S4 and table S6
    crush_Efactors = {
        "MRO": 0.67,  # [kg CO2e / kWh] emissions per electricity use in crushing per Midwest Reliability Organization
        "RFC": 0.57,  # [kg CO2e / kWh] emissions per electricity use in crushing per ReliabilityFirst Corporation
        "SERC": 0.61  # [kg CO2e / kWh] emissions per electricity use in crushing per SERC Reliability Corporation
    }
    if isinstance(efactor_overwrite, float):
        if efactor_overwrite > 0:
            crush_Efactors[Efactor_org] = efactor_overwrite
    # BOND WORK INDEX NOTES
    # ****************************************************** # 
    # [gbas] 
    # taken from Zhang et al., 2023, which they determine from the median
    # of published estimates ranging from 17.10-20.41 from Kanda and Kotake,
    # 2007; and Bond, 1961 (see Zhang et al's supp. eqn S1)
    # 
    # [cc]
    # Kanda and Kotake, 2007 (Chapter 12 of "Handbook of Powder Technology")
    # list limestone == 12.54; Dolomite == 11.27 (their table 3) 
    # Bond, 1961 part II 
    # list limestone == 11.61; Dolomite == 11.31
    #
    bondwork_indices = {
        "gbas": 18.67,  # Zhang et al., 2023 (see notes above)
        "gbas2": 18.67, # gbas2 differs from gbas only in cation stoichiometry, not grindability
        "baek23": 18.67,# Zhang et al., 2023 (see notes above)
        "baek23nodp": 18.67,# Zhang et al., 2023 (see notes above)
        "bridge": 18.67,# Zhang et al., 2023 (see notes above)
        "fo": 18.67,# Zhang et al., 2023 (see notes above)
        "cc": 12.10,    # mean of limestone estimates in Kanda and Kotake (2007) and Bond, 1961
        "wls": 8.33     # see Marco and Caterina, 2022 (https://www.proquest.com/docview/2781737786?pq-origsite=gscholar&fromopenview=true&sourcetype=Conference%20Papers%20&%20Proceedings)
    }
    if isinstance(bwi_overwrite, float):
        if bwi_overwrite > 0:
            bondwork_indices[mineral] = bwi_overwrite
    # ****************************************************** # 

    # add new empty columns to the dataframe
    new_columns = ["E_transport_tonCO2_tonRock", "E_crush_tonCO2_tonRock", "E_total_tonCO2_tonRock",
                    "E_transport_tonCO2_ha_yr", "E_crush_tonCO2_ha_yr", "E_total_tonCO2_ha_yr",
                    "E_transport_tonCO2_ha", "E_crush_tonCO2_ha", "E_total_tonCO2_ha"
                    ]
    for col in new_columns:
        df[col] = None


    # --- MODEL
    # loop through rows in df
    for index, row in df.iterrows():
        # collect inputs from rows
        p80_output = row['dustrad'] * 2   # convert radius to diameter
        rockTons_ha_yr = row[dustrate_name]
        dur = row['time_horizon']

        # transport per ton of rock
        transport_emissions_perRock = (truck_factor * truck_km) + (barge_factor * barge_km) + (barge_diesel_factor * barge_diesel_km) # [kg CO2e / ton rock]
        # crushing per ton of rock
        crush_Efactor = crush_Efactors[Efactor_org]  # either the "MRO", "RFC", or "SERC" Efactor depending on "Efactor_org"
        crush_energy_perRock = 10 * bondwork_indices[mineral] * ((1 / math.sqrt(p80_output)) - (1 / math.sqrt(p80_input)))   # [kWh / ton rock]
        crush_emissions_perRock = crush_Efactor * crush_energy_perRock  # [kg CO2e / ton rock]
        # total per ton of rock
        total_perRock = transport_emissions_perRock + crush_emissions_perRock
        perRock_dict = {
            "E_transport_tonCO2_tonRock": transport_emissions_perRock/1e3,  # divide by 1e3 to convert kg to tons
            "E_crush_tonCO2_tonRock": crush_emissions_perRock/1e3, 
            "E_total_tonCO2_tonRock": total_perRock/1e3
        }
        
        # transport per hectare per year 
        transport_emissions_perHaYr = transport_emissions_perRock * rockTons_ha_yr # [kg CO2e / ha / yr]
        # crushing per hectare per year
        crush_emissions_perHaYr = crush_emissions_perRock * rockTons_ha_yr  # [kg CO2e / ha / yr]
        # total per hectare per year
        total_perHaYr = transport_emissions_perHaYr + crush_emissions_perHaYr
        # add to dict
        perHaYr_dict = {
            "E_transport_tonCO2_ha_yr": transport_emissions_perHaYr/1e3, # divide by 1e3 to convert kg to tons
            "E_crush_tonCO2_ha_yr": crush_emissions_perHaYr/1e3, 
            "E_total_tonCO2_ha_yr": total_perHaYr/1e3
        }
        
        # transport per hectare (given total duration of rock applicatioin)
        transport_emissions_perHa = transport_emissions_perRock * (rockTons_ha_yr * dur)  # [kg CO2e / ha]
        # crushing per hectare (given total duration of rock application)
        crush_emissions_perHa = crush_emissions_perRock * (rockTons_ha_yr * dur)  # [kg CO2e / ha]
        # total per hectare 
        total_perHa = transport_emissions_perHa + crush_emissions_perHa
        # add to dict
        perHa_dict = {
            "E_transport_tonCO2_ha": transport_emissions_perHa/1e3, # divide by 1e3 to convert kg to tons
            "E_crush_tonCO2_ha": crush_emissions_perHa/1e3, 
            "E_total_tonCO2_ha": total_perHa/1e3
        }

        # --- collate results
        # per rock
        df.loc[index, perRock_dict.keys()] = perRock_dict.values()
        # per ha / yr
        df.loc[index, perHaYr_dict.keys()] = perHaYr_dict.values()
        # per ha
        df.loc[index, perHa_dict.keys()] = perHa_dict.values()
    # add columns for the upstream processes
    df['truck_km'] = truck_km
    df['barge_km'] = barge_km
    df['barge_diesel_km'] = barge_diesel_km
    df['p80_input'] = p80_input
    df['Efactor_org'] = Efactor_org
    df['bondwork_index'] = bondwork_indices[mineral]

    # some columns are being returned as objects when they should be floats
    # convert them if possible 
    object_columns = df.select_dtypes(include=['object']).columns
    # remove columns that are supposed to be strings from this check
    obj_to_remove = "Efactor_org"
    if obj_to_remove in object_columns:
        object_columns = object_columns.drop(obj_to_remove)
    # attempt to convert remaining columns to float
    for col in object_columns:
        try:
            df[col] = pd.to_numeric(df[col], errors="raise")  # Raises error for invalid conversion
        except ValueError as e:
            print(f"Emissions Calculator: Cannot convert column '{col}' to float: {e}. Ignore if expected.")
    
    # return the result
    return df



def emissions_calculator_ds(
    ds: xr.Dataset,
    feedstock: str,
    time_horizon: int | float,
    p80_input: int | float,
    Efactor_org: str,
    truck_km_grid: np.array = np.linspace(0, 400, 20),
    barge_km_grid: np.array = np.linspace(0, 400, 20),
    barge_diesel_km_grid: np.array = np.linspace(0, 400, 20),
    efactor_overwrite: bool=None,
    bwi_overwrite: bool=None,
    expand_over_time: bool=False,
)->xr.Dataset:
    '''
    Given a dataset with coordinates for the amount and grain size of dust,
    compute the emissions scenarios for crushing and transport in a wide range
    of emissions scenarios. 

    Return a dataset with the relevant emissions variables

    Uses the constants / model of Zhang et al., 2023 
    (https://doi.org/10.1021/acs.est.3c01658)

    Parameters
    ----------
    ds : xr.Dataset
        Must have coordinates named "dustrate_ton_ha_yr", "dustrad", and 
        a variable for "timehorizon_yr"
    feedstock : str
        name of the feedstock used (for finding the bond work index)

    # -- rocks -- # 
    p80_input : float
        [microns] p80 diameter before crushing the rock (see Zhang et al., 2023 table S5 for example data)
    
    # -- transport -- # 
    truck_km_grid : np.array
        [km] distance traveled by truck in transporting material
    barge_km_grid : np.array
        [km] distance traveled by barge in transporting material
    barge_diesel_km_grid : np.array
        [km] distance traveled by diesel barge
    
    bwi_overwrite : float
        if we don't want to use the default bond work index we can overwrite it here.
        only overwritten if it's a value > 0
    
    # -- crushing -- # 
    Efactor_org : str
        ["MRO" | "RFC" | "SERC"] which crushing Efactor to use
    efactor_overwrite : float
        if we don't want to use the default e factor we can overwrite it here.
        only overwritten if it's a value > 0

    # -- other -- # 
    expand_over_time : bool
        if True, expand the per-year data over time horizon coordinate
        and compute the cumulative emissions at each time step. Otherwise,
        do nothing
    
    Returns
    -------
    xr.Dataset
        Defined over coordinates of dustrate_ton_ha_yr, dustrad, truck_km, barge_km
        and barge_diesel_km. Variables include emissions for transport, crushing, and 
        total (sum of transport and crushing) per ton of rock; per ha/yr; and per yr. 

    '''
    # --- create a new empty dataset
    dustrate_ton_ha_yr = xr.DataArray(ds.dustrate_ton_ha_yr.values, dims="dustrate_ton_ha_yr")
    dustrad = xr.DataArray(ds.dustrad.values, dims="dustradius")
    truck_km = xr.DataArray(truck_km_grid, dims="truck_km_grid")
    barge_km = xr.DataArray(barge_km_grid, dims="barge_km_grid")
    diesel_km = xr.DataArray(barge_diesel_km_grid, dims="diesel_km_grid")

    dsout = xr.Dataset(
        coords=dict(
            dustrate_ton_ha_yr=ds.dustrate_ton_ha_yr.values,
            dustrad=ds.dustrad.values,
            truck_km=truck_km_grid,
            barge_km=barge_km_grid,
            barge_diesel_km=barge_diesel_km_grid
        )
    )

    # --- CONSTANTS
    # transport -- Zhang et al., 2023, table S1
    truck_factor = 0.0996    # [kg CO2e / ton / km] 
    barge_factor = 0.0282    # [kg CO2e / ton / km] 
    barge_diesel_factor = 0.00534  # [kg CO2e / ton / km] 

    # crushing -- Zhang et al., 2023, table S4 and table S6
    crush_Efactors = {
        "MRO": 0.67,  # [kg CO2e / kWh] emissions per electricity use in crushing per Midwest Reliability Organization
        "RFC": 0.57,  # [kg CO2e / kWh] emissions per electricity use in crushing per ReliabilityFirst Corporation
        "SERC": 0.61  # [kg CO2e / kWh] emissions per electricity use in crushing per SERC Reliability Corporation
    }
    if isinstance(efactor_overwrite, float):
        if efactor_overwrite > 0:
            crush_Efactors[Efactor_org] = efactor_overwrite
    # BOND WORK INDEX NOTES
    # ****************************************************** # 
    # [gbas] 
    # taken from Zhang et al., 2023, which they determine from the median
    # of published estimates ranging from 17.10-20.41 from Kanda and Kotake,
    # 2007; and Bond, 1961 (see Zhang et al's supp. eqn S1)
    # 
    # [cc]
    # Kanda and Kotake, 2007 (Chapter 12 of "Handbook of Powder Technology")
    # list limestone == 12.54; Dolomite == 11.27 (their table 3) 
    # Bond, 1961 part II 
    # list limestone == 11.61; Dolomite == 11.31
    #
    bondwork_indices = {
        "gbas": 18.67,  # Zhang et al., 2023 (see notes above)
        "gbas2": 18.67, # gbas2 differs from gbas only in cation stoichiometry, not grindability
        "baek23": 18.67,# Zhang et al., 2023 (see notes above)
        "baek23nodp": 18.67,# Zhang et al., 2023 (see notes above)
        "bridge": 18.67,# Zhang et al., 2023 (see notes above)
        "fo": 18.67,# Zhang et al., 2023 (see notes above)
        "cc": 12.10,    # mean of limestone estimates in Kanda and Kotake (2007) and Bond, 1961
        "wls": 8.33     # see Marco and Caterina, 2022 (https://www.proquest.com/docview/2781737786?pq-origsite=gscholar&fromopenview=true&sourcetype=Conference%20Papers%20&%20Proceedings)
    }
    if isinstance(bwi_overwrite, float):
        if bwi_overwrite > 0:
            bondwork_indices[mineral] = bwi_overwrite
    # ****************************************************** # 


    # --- COMPUTE TRANSPORT EMISSIONS
    da_transport_perRock = ((truck_factor * dsout['truck_km']) + (barge_factor * dsout['barge_km']) + (barge_diesel_factor * dsout['barge_diesel_km'])) / 1e3   # [tonne CO2e / ton rock]
    da_transport_perYr = da_transport_perRock * dsout['dustrate_ton_ha_yr']  # [tonne CO2e / ha / yr]
    da_transport_perHa = da_transport_perYr * time_horizon  # [tonne CO2e / ha]

    # --- COMPUTE CRUSHING EMISSIONS
    crush_Efactor = crush_Efactors[Efactor_org]  # either the "MRO", "RFC", or "SERC" Efactor depending on "Efactor_org"
    crush_energy_perRock = 10 * bondwork_indices[feedstock] * ((1 / np.sqrt(dsout['dustrad']*2)) - (1 / np.sqrt(p80_input)))   # [kWh / ton rock]
    da_crush_perRock = crush_Efactor * crush_energy_perRock / 1e3  # [tonne CO2e / ton rock]
    da_crush_perYr = da_crush_perRock * dsout['dustrate_ton_ha_yr']  # [tonne CO2e / yr]
    da_crush_perHa = da_crush_perYr * time_horizon  # [tonne CO2e / ha]

    # --- COMPUTE TOTAL EMISSIONS
    da_total_perRock = da_transport_perRock + da_crush_perRock
    da_total_perYr = da_transport_perYr + da_crush_perYr
    da_total_perHa = da_transport_perHa + da_crush_perHa


    # --- add to the dataset
    # [ PER ROCK ]
    dsout['E_transport_perRock'] = da_transport_perRock
    dsout['E_crush_perRock'] = da_crush_perRock
    dsout['E_total_perRock'] = da_total_perRock
    # add attributes
    for var in ['E_transport_perRock', 'E_crush_perRock', 'E_total_perRock']:
        dsout[var].attrs["units"] = "tonne CO2e / tonne rock"
        dsout[var].attrs["long_name"] = "Tonnes of CO2e emitted per tonne of rock"

    # [ PER YEAR ]
    dsout['E_transport_perYr'] = da_transport_perYr
    dsout['E_crush_perYr'] = da_crush_perYr
    dsout['E_total_perYr'] = da_total_perYr
    # add attributes
    for var in ['E_transport_perYr', 'E_crush_perYr', 'E_total_perYr']:
        dsout[var].attrs["units"] = "tonne CO2e / ha / year"
        dsout[var].attrs["long_name"] = "Tonnes of CO2e emitted per hectare of application per year"

    # [ PER HA ]
    dsout['E_transport_perHa'] = da_transport_perHa
    dsout['E_crush_perHa'] = da_crush_perHa
    dsout['E_total_perHa'] = da_total_perHa
    # add attributes
    for var in ['E_transport_perHa', 'E_crush_perHa', 'E_total_perHa']:
        dsout[var].attrs["units"] = "tonne CO2e / ha"
        dsout[var].attrs["long_name"] = "Tonnes of CO2e emitted per hectare, cumulative over the time_horizon"

    # --- expand over time if desired
    if expand_over_time:
        # vars to expand over time 
        vars_to_expand = ["E_transport_perYr", "E_total_perYr", "E_crush_perYr"]
        time = ds["time"]
        # get timestep lengths (yr)
        dt = time.diff("time")
        dt = dt.reindex(time=time, method="bfill")
        # broadcast vars over time
        tmp_dsout = dsout[vars_to_expand].broadcast_like(ds["time"])
        # get cumulative emissions
        for var in vars_to_expand:
            dsout[var] = tmp_dsout[var]
            dsout[f'{var}_int'] = (tmp_dsout[var] * dt).cumsum("time")
            # add attributes
            dsout[f'{var}_int'].attrs["units"] = "cumulative tonne CO2e / ha / year"
            dsout[f'{var}_int'].attrs["long_name"] = "Cumulative tonnes of CO2e emitted per hectare of application per year"


    # ... return result
    return dsout


def df_to_ds_with_time( # for the cdr_ds set of functions
        dims_full: list,
        dfx: list,
        timetol: float=1e-2
)->xr.Dataset:
    """
    convert pandas dataframe to xarray dataset in cases where we have multiple 
    models with the same number of timesteps but slightly different time values
    --> we make timestep index a dimension, then take the mean time per timestep
        then add that as a dimension and return the result
    """
    groupby_list = [x for x in dims_full if x != "time"]
    # replace `time` with `time_step`
    dims_full = ['time_step' if c == 'time' else c for c in dims_full]
    # add time_step to the dataframe
    dfx['time_step'] = (
        dfx.groupby(groupby_list)
        .cumcount()
    )
    # create the cdr dataset
    dfx_full_idx = dfx.set_index(dims_full)
    dsx_full = xr.Dataset.from_dataframe(dfx_full_idx)
    # look at variability in time across experiments
    other_dims = [d for d in dsx_full.dims if d != "time_step"]

    return dsx_full # (keep it in time index terms, convert to time later)
    # # compute std and mean across experiments at each time_step
    # time_std = dsx_full['time'].std(dim=other_dims)
    # time_mean = dsx_full['time'].mean(dim=other_dims)
    # max_abs_diff = time_std.max().item()
    # max_rel_diff = (time_std / time_mean).max().item()
    # # --- [ TROUBLESHOOT ] ---
    # # print("Max absolute difference:", max_abs_diff)
    # # print("Max relative difference:", max_rel_diff)
    # if max_abs_diff < timetol: # add time to cdr dataset
    #     dsx_full = dsx_full.assign_coords(time=time_mean)
    #     dsx_full = dsx_full.swap_dims({"time_step": "time"})
    #     return dsx_full
    # else:
    #     raise ValueError("Simulation timesteps are not similar enough to make an xarray dataset")
    
def convert_time_index_to_time(
        dsx: xr.Dataset,
        timetol: float=1e-2,
)-> xr.Dataset:
    """
    read in a dataset that has the time index as a coordinate
    and time itself as a variable. Convert the variable to a new
    coordinate that replaces the time index so long as the time 
    across all dimensions is similar enough
    """
    # dims across which to compare time data
    other_dims = [d for d in dsx.dims if d != "time_step"]
    # compute std and mean across experiments at each time_step
    time_std = dsx['time'].std(dim=other_dims)
    time_mean = dsx['time'].mean(dim=other_dims)
    max_abs_diff = time_std.max().item()
    max_rel_diff = (time_std / time_mean).max().item()
    # --- [ TROUBLESHOOT ] ---
    # print("Max absolute difference:", max_abs_diff)
    # print("Max relative difference:", max_rel_diff)
    if max_abs_diff < timetol: # add time to cdr dataset
        dsx = dsx.assign_coords(time=time_mean)
        dsx = dsx.swap_dims({"time_step": "time"})
        return dsx
    else:
        raise ValueError("Simulation timesteps are not similar enough to make an xarray dataset")


def cdr_ds(
    cdr_dict: dict,
    dims: list,
    cdr_calc_list: list,
    loss_percents: np.array = np.linspace(100,1,50),
    skip_loss: bool=False,
    convert_time_to_timestep: bool=False,
) -> xr.Dataset:
    """
    Generate an xarray dataset for the removal and emissions fluxes
    for a given case (e.g., just silicate or just calcite runs) 
    and calculate downstream loss effects on CDR.

    Parameters
    ----------
    cdr_dict : dict
        a cdr dictionary, generally the output from `cdr_int_per_group`
    dims : list
        list of dimensions that correspond to column names in all of the 
        dfs within the dictionary (generally equivalent to `dfin_cols_to_keep`)
    cdr_calc_list : list
        list whose elements correspond to the keys in the dictionary above
    loss_percents : np.array
        array of values [0,100] for computing downstream loss
    skip_loss : bool
        whether to skip the loss calculation
    convert_time_to_timestep : bool
        Sometimes time steps slightly differ between runs that have the same number
        of steps and this creates annoying gaps in the dataset. if True, the function 
        converts time to a time step index and uses that for the coordinate, saving 
        time itself as a variable. (only applied if time is in `dims`)
    
    Returns 
    -------
    xr.Dataset
        Dataset with the removals and emissions vars
    """
    # create an empty dict
    dsout_dict = {}

    #
    # iterate through options in the calc_list
    #

    # --- OPTION 1: "co2_flx" ------------------------
    if "co2_flx" in cdr_calc_list:
        tc = 'co2_flx'
        print(f"solving {tc}")
        dfin = cdr_dict[tc]
        dsout_dict[tc] = co2_flx_cdr_ds(
            dfin, dims, loss_percents, skip_loss=skip_loss, 
            convert_time_to_timestep=convert_time_to_timestep,
        )
    
    # --- OPTION 2: "camg_flx" -----------------------
    if "camg_flx" in cdr_calc_list:
        tc = 'camg_flx'
        print(f"solving {tc}")
        dfin = cdr_dict[tc]
        # (note, we use the same cation_flx_cdr fxn for both camg and totcat
        #  because it just calculates based on whatever cats are provided)
        dsout_dict[tc] = cation_flx_cdr_ds(
            dfin, dims, loss_percents, cation_flag='camg', 
            cdr_calc_cols=["co2pot_tot", "co2pot_adv", "DICpot_adv"], skip_loss=skip_loss,
            convert_time_to_timestep=convert_time_to_timestep
        )
    
    # --- OPTION 3: "totcat_flx" ---------------------
    if "totcat_flx" in cdr_calc_list:
        tc = 'totcat_flx'
        print(f"solving {tc}")
        dfin = cdr_dict[tc]
        # (note, we use the same cation_flx_cdr fxn for both camg and totcat
        #  because it just calculates based on whatever cats are provided)
        dsout_dict[tc] = cation_flx_cdr_ds(
            dfin, dims, loss_percents, cation_flag='totcat', 
            cdr_calc_cols=["co2pot_tot", "co2pot_adv", "DICpot_adv"], skip_loss=skip_loss,
            convert_time_to_timestep=convert_time_to_timestep
        )
    
    # --- OPTION 4: "carbalk_flx" --------------------
    if "carbalk_flx" in cdr_calc_list:
        print('carbalk is WIP !!')
        # print("solving carbalk_flx")
        # outdict['carbalk_flx'] = carbalk_flx_cdr()
    
    # --- OPTION 5: "rockdiss" -----------------------
    if "rockdiss" in cdr_calc_list:
        tc = 'rockdiss'
        print(f"solving {tc}")
        dfin = cdr_dict[tc]
        dsout_dict[tc] = rockdiss_ds(dfin, dims, convert_time_to_timestep=convert_time_to_timestep)

    # --- OPTION 6: "catbudget_flx" ------------------------------------------
    if "catbudget_flx" in cdr_calc_list:
        tc = 'catbudget_flx'
        print(f"solving {tc}")
        dfin_cb = cdr_dict[tc]
        dsout_dict[tc] = catbudget_cdr_ds(
            dfin_cb, dims,
            convert_time_to_timestep=convert_time_to_timestep,
        )

    # return dsout_dict
    # merge the dictionary of datasets into a single ds
    dsout = xr.merge(dsout_dict.values())
    if convert_time_to_timestep:
        dsout = convert_time_index_to_time(dsout)
        return dsout
    # 
    # return the result
    #
    return dsout


def co2_flx_cdr_ds(
    dfin: pd.DataFrame, 
    dims: list, 
    loss_percents: np.array,
    cdr_calc_cols: list=["cdr_dif", "cdr_adv", "cdr_adv_plus_newSIC", "cdr_SIConly", "tot_adv"],
    dustrate_name: str = "dustrate_ton_ha_yr",
    skip_loss: bool=False,
    convert_time_to_timestep: bool=False,
) -> xr.Dataset:
    """
    Get the removals from the co2flx calculations in a dataset form
    and calculate downstream loss effects on CDR.
    
    Parameters
    ----------
    dfin : pd.DataFrame
        the respective df within the cdr_dict* from the `cdr_int_per_group` function
    dims : list
        list of dimensions for the xr dataset. must correspond to columns in dfin. generally
        the same as dfin_cols_to_keep 
    loss_percents : np.array
        values indicating the percent downstream loss to calculate
    cdr_calc_cols : list
        columns required for the removal calculation when accounting for 
        downstream loss
    dustrate_name : str
        column name for the dustrate
    skip_loss : bool
        [True] to skip computing the downstream loss factor, False to include it
        (note, set skip_loss to True if it's not time integrated)
    convert_time_to_timestep : bool
        Sometimes time steps slightly differ between runs that have the same number
        of steps and this creates annoying gaps in the dataset. if True, the function 
        converts time to a time step index and uses that for the coordinate, saving 
        time itself as a variable. (only applied if time is in `dims`)

    Returns
    -------
    xr.Dataset
        removals flux dataset
    """
    # get a list of just the cdr vars
    cdr_cols = [col for col in cdr_calc_cols if col.startswith('cdr')]

    # add loss_percents to the dims
    if skip_loss:
        dims_full = dims
    else:
        dims_full = dims + ['loss_percent']


    # determine the columns we want to drop
    cols_to_discard = ['units', 'flx_type', 'cdr_fxn', 'time_horizon']
    dfx = dfin.drop(columns=cols_to_discard)

    # create a cdr df to loop over loss percents
    cdr_cols_andDims = dims + cdr_calc_cols
    dfx_cdr_short = dfx[cdr_cols_andDims]

    if skip_loss:
        dfx_cdr_full = dfx_cdr_short.copy()
    else:
        # track index
        lossdx = 0
        for loss in loss_percents:
            # pull out just the full dat
            tdfx = dfx_cdr_short.copy()

            # add the loss percent
            tdfx['loss_percent'] = loss

            # compute updated loss
            # (calculation is same for sil and cc because we're taking the loss relative
            #  to the TOTAL amount of carbon exported from the domain)
            # (we use np.maximum so that any negative advective fluxes representing
            #  a decrease (unlikely) doesn't increase CDR)
            tdfx = tdfx.apply(lambda x: x - np.maximum(((loss/100) * tdfx['tot_adv']), 0) if x.name in cdr_cols else x)

            # bring together
            if lossdx == 0:
                dfx_cdr_full = tdfx.drop(columns='tot_adv').copy()
            else:
                dfx_cdr_full = pd.concat([dfx_cdr_full.copy(), tdfx.drop(columns='tot_adv').copy()])
            
            lossdx += 1

    # [ add timestep if asked to ]
    if ('time' not in dims_full) or (not convert_time_to_timestep):
        # create the cdr dataset
        dfx_full_idx = dfx_cdr_full.set_index(dims_full)
        dsx_full = xr.Dataset.from_dataframe(dfx_full_idx)
        # create the flux datasets (not defined over loss_percent dimension)
        # using the columns that aren't already in the full ds
        dfx_idx = dfx.drop(columns=list(dsx_full.data_vars)).set_index(dims)
        dsx_x = xr.Dataset.from_dataframe(dfx_idx)
    else:
        dsx_full = df_to_ds_with_time(dims_full, dfx_cdr_full)
        # create the flux datasets (not defined over loss_percent dimension)
        # using the columns that aren't already in the full ds
        dfx_idx = dfx.drop(columns=list(dsx_full.data_vars))
        dsx_x = df_to_ds_with_time(dims, dfx_idx)
    # ----
    # create the flux dataset (not defined over loss_percent dimension)
    # using the columns that aren't already in the full ds
    # dfx_idx = dfx.drop(columns=list(dsx_full.data_vars)).set_index(dims)
    # dsx_x = xr.Dataset.from_dataframe(dfx_idx)

    # bring them together
    dsx = xr.merge([dsx_full, dsx_x])

    # --- add attributes
    dsx['time_horizon'] = dfin['time_horizon'][0]
    dsx.attrs['flx_type'] = dfin['flx_type'][0]

    for var_name in dsx.data_vars:
        dsx[var_name].attrs["units"] = dfin['units'][0]
        dsx[var_name].attrs["cdr_fxn"] = dfin['cdr_fxn'][0]
        if var_name in cdr_cols: # label the CDR columns so we can easily grab them later!
            dsx[var_name].attrs["cdr_var"] = True
        else:
            dsx[var_name].attrs["cdr_var"] = False
    
    # return result
    return dsx


def cation_flx_cdr_ds(
    dfin: pd.DataFrame, 
    dims: list, 
    loss_percents: np.array,
    cation_flag: str,
    cdr_calc_cols: list=["co2pot_tot", "co2pot_adv", "DICpot_adv"],
    dustrate_name: str = "dustrate_ton_ha_yr",
    skip_loss: bool=False,
    convert_time_to_timestep: bool=False,
) -> xr.Dataset:
    """
    Get the removals from the cation calculations in a dataset form

    Parameters
    ----------
    dfin : pd.DataFrame
        the respective df within the cdr_dict* from the `cdr_int_per_group` function
    dims : list
        list of dimensions for the xr dataset. must correspond to columns in dfin. generally
        the same as dfin_cols_to_keep 
    loss_percents : np.array
        values indicating the percent downstream loss to calculate
    cation_flag : str
        flag to add to the cdr column names so we know which calculation we're 
        talking about when we merge the datasets together
    cdr_calc_cols : list
        columns required for the removal calculation when accounting for 
        downstream loss
    dustrate_name : str
        column name for the dustrate
    skip_loss : bool
        [True] to skip computing the downstream loss factor, False to include it
        (note, set skip_loss to True if it's not time integrated)
    convert_time_to_timestep : bool
        Sometimes time steps slightly differ between runs that have the same number
        of steps and this creates annoying gaps in the dataset. if True, the function 
        converts time to a time step index and uses that for the coordinate, saving 
        time itself as a variable. (only applied if time is in `dims`)

    Returns
    -------
    xr.Dataset
        removals flux dataset
    """
    # update the cdr columns based on what's in dfin
    col_tot_cdr = [col for col in dfin.columns if col.startswith("co2pot_tot")]
    col_adv_cdr = [col for col in dfin.columns if col.startswith("co2pot_adv")]
    idx_tot = cdr_calc_cols.index('co2pot_tot')
    idx_adv = cdr_calc_cols.index('co2pot_adv')
    cdr_calc_cols[idx_tot] = col_tot_cdr[0]
    cdr_calc_cols[idx_adv] = col_adv_cdr[0]

    # get a list of just the cdr vars
    cdr_cols = [col for col in cdr_calc_cols if col.startswith('co2pot')]

    # add loss_percents to the dims
    if skip_loss:
        dims_full = dims
    else:
        dims_full = dims + ['loss_percent']

    # determine the columns we want to drop
    cols_to_discard = ['units', 'flx_type', 'cdr_fxn', 'time_horizon']
    dfx = dfin.drop(columns=cols_to_discard)

    # create a cdr df to loop over loss percents
    cdr_cols_andDims = dims + cdr_calc_cols
    dfx_cdr_short = dfx[cdr_cols_andDims]

    if skip_loss:
        dfx_cdr_full = dfx_cdr_short.copy()
    else:
        # track index
        lossdx = 0
        for loss in loss_percents:
            # pull out just the full dat
            tdfx = dfx_cdr_short.copy()

            # add the loss percent
            tdfx['loss_percent'] = loss

            # compute updated loss
            # (calculation is same for sil and cc because we're taking the loss relative
            #  to the TOTAL amount of carbon exported from the domain)
            # (we use np.maximum so that any negative advective fluxes representing
            #  a decrease (unlikely) doesn't increase CDR)
            tdfx = tdfx.apply(lambda x: x - np.maximum(((loss/100) * tdfx['DICpot_adv']), 0) if x.name in cdr_cols else x)

            # bring together
            if lossdx == 0:
                dfx_cdr_full = tdfx.drop(columns='DICpot_adv').copy()
            else:
                dfx_cdr_full = pd.concat([dfx_cdr_full.copy(), tdfx.drop(columns='DICpot_adv').copy()])
            
            lossdx += 1

    # [ add timestep if asked to ]
    if ('time' not in dims_full) or (not convert_time_to_timestep):
        # create the cdr dataset
        dfx_full_idx = dfx_cdr_full.set_index(dims_full)
        dsx_full = xr.Dataset.from_dataframe(dfx_full_idx)
        # create the flux datasets (not defined over loss_percent dimension)
        # using the columns that aren't already in the full ds
        dfx_idx = dfx.drop(columns=list(dsx_full.data_vars)).set_index(dims)
        dsx_x = xr.Dataset.from_dataframe(dfx_idx)
    else:
        dsx_full = df_to_ds_with_time(dims_full, dfx_cdr_full)
        # create the flux datasets (not defined over loss_percent dimension)
        # using the columns that aren't already in the full ds
        dfx_idx = dfx.drop(columns=list(dsx_full.data_vars))
        dsx_x = df_to_ds_with_time(dims, dfx_idx)

    # ----
    # create the flux datasets (not defined over loss_percent dimension)
    # using the columns that aren't already in the full ds
    # dfx_idx = dfx.drop(columns=list(dsx_full.data_vars)).set_index(dims)
    # dsx_x = xr.Dataset.from_dataframe(dfx_idx)

    # bring them together
    dsx = xr.merge([dsx_full, dsx_x])

    # --- add attributes
    dsx['time_horizon'] = dfin['time_horizon'][0]
    dsx.attrs['flx_type'] = dfin['flx_type'][0]

    for var_name in dsx.data_vars:
        dsx[var_name].attrs["cdr_fxn"] = dfin['cdr_fxn'][0]
        if var_name in cdr_cols: # label the CDR columns so we can easily grab them later!
            dsx[var_name].attrs["cdr_var"] = True
        else:
            dsx[var_name].attrs["cdr_var"] = False
            dsx[var_name].attrs["units"] = dfin['units'][0]

        if var_name in cdr_calc_cols:
            newvar_name = f'{var_name}_{cation_flag}'   # update the var name so we don't mix it up with another
                                                        # cation based calculation
            dsx = dsx.rename({var_name: newvar_name}).copy()

    # return result
    return dsx



def rockdiss_ds(
    dfin: pd.DataFrame,
    dims: list,
    convert_time_to_timestep: bool=False,
) -> xr.Dataset:
    """
    Get the rockdiss and emissions data in dataset form

    Parameters
    ----------
    dfin : pd.DataFrame
        the respective df within the cdr_dict* from the `cdr_int_per_group` function
    dims : list
        list of dimensions for the xr dataset. must correspond to columns in dfin. generally
        the same as dfin_cols_to_keep 
    convert_time_to_timestep : bool
        Sometimes time steps slightly differ between runs that have the same number
        of steps and this creates annoying gaps in the dataset. if True, the function 
        converts time to a time step index and uses that for the coordinate, saving 
        time itself as a variable. (only applied if time is in `dims`)
    
    Returns
    -------
    xr.Dataset
        removals flux dataset
    """
    # determine the columns we want to drop
    cols_to_discard = ['truck_km', 'barge_km', 'barge_diesel_km', 'p80_input',
                        'Efactor_org', 'bondwork_index', 'time_horizon']
    cols_to_discard_present = [c for c in cols_to_discard if c in dfin.columns]
    dfx = dfin.drop(columns=cols_to_discard_present)

    # rename the adv column to be clearer that it's rock
    dfx = dfx.rename(columns = {'adv': 'adv_feedstock'})

    # [ add timestep if asked to ]
    if ('time' not in dims) or (not convert_time_to_timestep):
        # create the flux dataset (not defined over loss_percent dimension)
        dfx_idx = dfx.set_index(dims)
        dsx = xr.Dataset.from_dataframe(dfx_idx)
    else:
        # create the flux dataset (not defined over loss_percent dimension)
        dsx = df_to_ds_with_time(dims, dfx)
    # ----

    # --- add attributes back in
    for thiscol in cols_to_discard_present:
        dsx[thiscol] = dfin[thiscol][0]
    
    for var_name in dsx.data_vars:
        dsx[var_name].attrs["cdr_fxn"] = 'rockdiss'

    # return result
    return dsx


def cdr_fs_vs_counterfactual(
    cdr_ds_dict: dict,
    counterfact_fs: str,
    cf_apprate_fixed: float,
    cf_dustrad_fixed: float,
    select_nearest_cfdust: bool,
    E_var: str="E_total_tonCO2_ha",
) -> xr.Dataset:
    """
    Compute netCDR relative to some counterfactual application case. 
    The counterfactual application rate and radius are prescribed. The
    function takes in a dictionary of CDR datasets, one for each 
    feedstock. CDR in these datasets is calculated relative to a zero-
    application baseline. Each dictionary element is named with the 
    feedstock. The counterfactual feedstock name (counterfact_fs) must
    be present in the cdr_ds_dict.keys().

    Parameters
    ----------
    cdr_ds_dict : dict
        dictionary of cdr xr datasets, one for each feedstock (with cdr
        computed relative to zero application baseline). dictionary keys
        should be the feedstock names. 
    counterfact_fs : str
        the name of the feedstock to use as the counterfactual (must be
        a name present in the keys of cdr_ds_dict).
    cf_apprate_fixed : float
        counterfactual rock application rate
    cf_dustrad_fixed : float
        counterfactual radius of rock dust
    select_nearest_cfdust : bool
        [True | False] if True, then select the nearest cf dust flux 
        and radius to the prescribed one. Otherwise, if the prescribed
        values aren't present in the dict, return an error. 
    E_var : str
        The name of the variable for emissions in the xr datasets that are
        within cdr_ds_dict. 

    Returns 
    -------
    xr.Dataset
        Function generates a new xarray dataset with netCDR calculations
        for each feedstock relative to the counterfactual. 
    """
    # pull out the dataset that contains the counterfactual case
    dsall_cf = cdr_ds_dict[counterfact_fs].copy()

    # pull out the counterfactual df
    if select_nearest_cfdust:  # then use method='nearest'
        ds_cf = dsall_cf.sel(dustrad = cf_dustrad_fixed, dustrate_ton_ha_yr = cf_apprate_fixed, method='nearest').copy()
        # update cc pars to what was taken out 
        cf_apprate_fixed = ds_cf['dustrate_ton_ha_yr'].values
        cf_dustrad_fixed = ds_cf['dustrad'].values
    else:
        try:
            ds_cf = dsall_cf.sel(dustrad = cf_dustrad_fixed, dustrate_ton_ha_yr = cf_apprate_fixed).copy()
        except:
            if cf_apprate_fixed not in dsall_cf.dustrate_ton_ha_yr:
                available_dustrates = dsall_cf.dustrate_ton_ha_yr.values
                raise ValueError(f"The counterfactual dust rate is not in the xr.dataset. Please select a value from the following: {available_dustrates}")
            elif cf_dustrad_fixed not in dsall_cf.dustrad:
                available_dustrads = dsall_cf.dustrad.values
                raise ValueError(f"The counterfactual dust radius is not in the xr.dataset. Please select a value from the following: {available_dustrads}")
            else:
                raise ValueError(f"Something went wrong. Consider setting `select_nearest_cfdust` to True.")


    # pull out the cdr vars
    cdr_variables = [var for var in ds_cf.variables if ds_cf[var].attrs.get("cdr_var") == True]

    # empty dataset to hold net CDR results
    ds_anom = xr.Dataset()

    # loop through feedstocks
    fs_dx = 0
    for fs, ds_fs in cdr_ds_dict.items():
        # feedstock-specific dataset that we'll add to the main one at the end
        tmpds_fs_anom = xr.Dataset()   
        # get the data for just this feedstock
        ds_tmp = ds_fs.copy()
        da_emiss_fs = ds_tmp[E_var].copy()  # deployment emissions
        da_emiss_cf = ds_cf[E_var].copy()   # counterfactual emissions

        # loop through CDR vars
        for cvar in cdr_variables:
            # extract just this variable
            tmpda_fs = ds_tmp[cvar].copy()
            tmpda_cf = ds_cf[cvar].copy()

            # --- compute CDR using each approach
            # [1] net R 
            cdrname1 = "netR"
            tmpda1 = tmpda_fs - tmpda_cf
            # [2] net R no negatives
            cdrname2 = "netR_noNeg"
            tmpda2 = tmpda_fs - np.maximum(0, tmpda_cf)
            # [3] simple subtraction
            cdrname3 = "simplesubtract"
            tmpda3 = (tmpda_fs - tmpda_cf) - (da_emiss_fs - da_emiss_cf)
            # [4] simple subtraction, no negatives
            cdrname4 = "simplesubtract_noNeg"
            tmpda4 = (tmpda_fs - np.maximum(0, tmpda_cf)) - (np.maximum(da_emiss_fs - da_emiss_cf, 0))
            # [5] conservative 
            cdrname5 = "conservative"
            tmpda5 = (tmpda_fs - np.maximum(0, tmpda_cf)) - da_emiss_fs

            # all data arrays are the same variable over different `cdr_calc` dimensions.
            # combine them and assign the `cdr_calc` coord
            tmpda = xr.concat([tmpda1, tmpda2, tmpda3, tmpda4, tmpda5], dim='cdr_calc')
            tmpda = tmpda.assign_coords(cdr_calc = [cdrname1, cdrname2, cdrname3, cdrname4, cdrname5]).copy()
            tmpds_cvar_anom = tmpda.to_dataset(name=cvar)
            tmpds_fs_anom[cvar] = tmpds_cvar_anom[cvar]
            # --------------------------------------
        # add feedstock specific variables
        tmpds_fs_anom['time_horizon'] = ds_tmp['time_horizon']
        
        # add the feedstock dimension to the tmpds
        tmpds_fs_anom = tmpds_fs_anom.expand_dims(dim={'feedstock': [fs]}).copy()

        if fs_dx == 0:
            ds_anom = tmpds_fs_anom.copy()
        else:
            ds_anom = xr.merge([ds_anom, tmpds_fs_anom.copy()])
        # update the index
        fs_dx += 1

    # add helpful variables
    ds_anom['cf_apprate'] = cf_apprate_fixed
    ds_anom['cf_dustrad'] = cf_dustrad_fixed

    # return result
    return ds_anom



# %% 
# ---------- POSTPROC PROFILE DATA
def read_profile_nc(
    outdir: str,
    filename_base: str,
    tdf: pd.Series,
    subdir: str="postproc_profs",
) -> xr.Dataset:
    '''
    read in a profile netcdf file to the processing script

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    filename_base : str
        base name of the SCEPTER flx file. Format is "[basename].pkl"
    tdf : pd.Series OR dict
        single row of the batch df containing the info for this specific run. Must have 'newrun_id_full'
    subdir : str
        name of the subdirectory within the run output directory where we'll fine the filename_base file
    '''
    # --- read in postproc profile data
    # get path to data
    fn = f"{filename_base}.nc"
    fn_path = os.path.join(outdir, tdf['newrun_id_full'], subdir, fn)
    # ----------------------------
    # check for S3
    if fn_path.startswith("s3://"): # then bring it in from s3
        import s3fs
        fs = s3fs.S3FileSystem()
        if not fs.exists(fn_path):
            print(f"Warning: batch profile processing {fn_path} could not be found.. returning NA")
            tmpds = None
        else:
            # [ worked on nebari, but not on local ! ]
            # with fs.open(fn_path, mode='rb') as fnx:
            #     tmpds = xr.open_dataset(fnx)
            # [ updated for local (netcdf4 can't read from s3 so we need h5netcdf !) ]
            tmpds = xr.open_dataset(
                fn_path,
                engine="h5netcdf"
            )
    # -----------------------------
    else:
        if not os.path.exists(fn_path):
            print(f"Warning: batch profile processing {fn_path} could not be found.. returning NA")
            tmpds = None
        else:
            tmpds = xr.open_dataset(fn_path)
    
    # return result
    return tmpds


def time_is_close(
    t1: float, 
    t2: float, 
    tol: float=1e-4,
)->bool:
    '''
    Returns True if all time values in t1 and t2 
    are elementwise close within tol.
    '''
    # get to array so we can take len (even if there is only one element)
    t1 = np.atleast_1d(t1)
    t2 = np.atleast_1d(t2)

    if len(t1) != len(t2):
        return False
    return np.all(np.abs(t1 - t2) <= tol)

def group_by_time_values_fuzzy(
    ds_list: list, 
    tol: float=1e-4,
)->list:
    '''
    groups xarray datasets by their time dimension with some 
    "fuzziness" to account for floating point errors

    Parameters
    ----------
    ds_list : list
    '''
    groups = []

    for i, ds in enumerate(ds_list):
        time_i = ds['time'].values
        matched = False
        for group in groups:
            ref_time = group['time']
            if time_is_close(time_i, ref_time, tol):
                group['indices'].append(i)
                matched = True
                break
        if not matched:
            groups.append({'time': time_i, 'indices': [i]})
    
    return groups

def filter_mismatched_time_coords_fuzzy(
    ds_list: list, 
    filename_base: str,
    tol: float=1e-4,
    printwarnings: bool=True,
)->list:
    '''
    Filter out datasets from a list whose time coordinates are mis-matched.
    Groups datasets by time coords and assumes the smallest group must be 
    wrong, removes those. Time coords are compared with "fuzziness" to 
    account for floating point errors. 

    Parameters
    ----------
    ds_list : list
        list of xarray datasets whose time coordinates to check for 
        consistency
    filename_base : str
        base name of the SCEPTER flx file. Format is "[basename].pkl"
    tol : float
        [yr, or time units of ds_list] the tolerance used for comparing
        time coords
    printwarnings: bool
        [default=True] whether to print out warnings about excluded 
        datasets. 
    '''
    groups = group_by_time_values_fuzzy(ds_list, tol=tol)

    # Find largest group
    largest_group = max(groups, key=lambda g: len(g['indices']))
    kept_indices = largest_group['indices']
    ref_time = largest_group['time']

    # identify dropped datasets
    all_indices = set(range(len(ds_list)))
    dropped_indices = sorted(all_indices - set(kept_indices))

    # align time coordinates of kept datasets to reference time
    aligned_ds_list = []
    for i in kept_indices:
        ds = ds_list[i]
        # ensure ref_time is always 1d (even if scalar)
        time_values = np.atleast_1d(ref_time)
        ds.coords['time'] = ('time', time_values)
        aligned_ds_list.append(ds)

    # report results
    if printwarnings:
        if dropped_indices:
            print(f"⚠️ Dropped {len(dropped_indices)} {filename_base} dataset(s) due to mismatched time values, idx: {dropped_indices}")

    # return result
    return aligned_ds_list


def _prof_batchprocess_singlevar_zarr(
    dfin: pd.DataFrame,
    outdir: str,
    batch_axes: list,
    filename_base: str,
    dustsp: str,
    zarr_path: str,
    print_progress: bool = False,
) -> xr.Dataset:
    """
    Two-pass, O(1-run) memory implementation of prof_batchprocess_singlevar.

    Pass 1 (lightweight): read only time coordinates and rockdiss pkl files to
    discover which runs are valid and what the full batch coordinate grid is.
    Pass 2 (data): read one run at a time, write its data to the correct region
    of the zarr store, then immediately discard it.

    Returns xr.open_zarr(zarr_path) (lazy).
    """
    import gc

    # ---- Pass 1: scan metadata ----
    run_meta = []  # (row, dustrate_mean, dustdf, time_array)
    for _, row in dfin.iterrows():
        tmpds = read_profile_nc(outdir, filename_base, row)
        if tmpds is None:
            continue
        t = tmpds["time"].values.copy()
        tmpds.close()
        del tmpds
        dustrate_mean, dustdf = get_timemean_dustrate(
            outdir, row, dustsp, return_df=True
        )
        run_meta.append((row, dustrate_mean, dustdf, t))

    if not run_meta:
        return xr.Dataset()

    # ---- Filter mismatched time coords (replicate filter_mismatched_time_coords_fuzzy) ----
    groups = []
    for i, (*_, t) in enumerate(run_meta):
        placed = False
        for g in groups:
            if time_is_close(t, g["ref_t"]):
                g["indices"].append(i)
                placed = True
                break
        if not placed:
            groups.append({"ref_t": t, "indices": [i]})
    best = max(groups, key=lambda g: len(g["indices"]))
    reference_time = best["ref_t"]
    valid_meta = [run_meta[i] for i in best["indices"]]
    dropped = len(run_meta) - len(valid_meta)
    if dropped:
        print(f"⚠️ {filename_base}: dropped {dropped} run(s) with mismatched time coords")

    # ---- Read first valid run to get inner dimension / variable structure ----
    first_row, first_dr, first_dustdf, _ = valid_meta[0]
    tmpds_first = read_profile_nc(outdir, filename_base, first_row)
    tmpds_first.coords["time"] = ("time", reference_time)
    if filename_base != "soil_ph" and first_dustdf is not None:
        _dd = xr.Dataset.from_dataframe(first_dustdf.set_index("time"))
        _dd = _dd.reindex(time=tmpds_first["time"])
        tmpds_first = xr.merge([tmpds_first, _dd["int_dust_ton_ha_yr"]])
    inner_dims = list(tmpds_first.dims)
    inner_coords = {d: tmpds_first[d].values.copy() for d in inner_dims}
    var_dims_map = {v: list(tmpds_first[v].dims) for v in tmpds_first.data_vars}
    var_dtypes = {v: tmpds_first[v].dtype for v in tmpds_first.data_vars}
    tmpds_first.close()
    del tmpds_first
    gc.collect()

    # ---- Build full batch coordinate grids ----
    batch_coord_vals = {}
    for col in batch_axes:
        if col == "dustrate_ton_ha_yr":
            vals = sorted({
                round(dr, 3) for _, dr, _, _ in valid_meta if dr is not None
            })
            if not vals and col in dfin.columns:
                vals = sorted(set(dfin[col].dropna().tolist()))
        else:
            vals = sorted(set(dfin[col].dropna().tolist())) if col in dfin.columns else []
        batch_coord_vals[col] = vals

    batch_sizes = tuple(len(batch_coord_vals[col]) for col in batch_axes)
    all_coords = {**{col: batch_coord_vals[col] for col in batch_axes}, **inner_coords}

    # ---- Create zarr scaffold (NaN / zero filled) ----
    template_vars = {}
    for var, var_dims in var_dims_map.items():
        full_dims = list(batch_axes) + var_dims
        inner_sizes = tuple(len(inner_coords[d]) for d in var_dims)
        full_shape = batch_sizes + inner_sizes
        coords_for_var = {k: v for k, v in all_coords.items() if k in full_dims}
        fill = np.nan if np.issubdtype(var_dtypes[var], np.floating) else 0
        template_vars[var] = xr.DataArray(
            np.full(full_shape, fill, dtype=var_dtypes[var]),
            dims=full_dims,
            coords=coords_for_var,
        )
    xr.Dataset(template_vars).to_zarr(zarr_path, mode="w")
    del template_vars
    gc.collect()

    # ---- Pass 2: fill zarr one run at a time ----
    canonical_dims = list(batch_axes) + inner_dims
    for run_idx, (row, dustrate_mean, dustdf, _) in enumerate(valid_meta):
        if print_progress and (run_idx % 50 == 0):
            print(f"  {filename_base}: writing run {run_idx + 1}/{len(valid_meta)}")
        tmpds = read_profile_nc(outdir, filename_base, row)
        if tmpds is None:
            continue
        tmpds.coords["time"] = ("time", reference_time)

        if filename_base != "soil_ph" and dustdf is not None:
            _dd = xr.Dataset.from_dataframe(dustdf.set_index("time"))
            _dd = _dd.reindex(time=tmpds["time"])
            tmpds = xr.merge([tmpds, _dd["int_dust_ton_ha_yr"]])

        # find the region slice in the zarr store for each batch axis
        region = {}
        for col in batch_axes:
            val = round(dustrate_mean, 3) if (col == "dustrate_ton_ha_yr" and dustrate_mean is not None) else row[col]
            arr = batch_coord_vals[col]
            try:
                idx = int(np.argmin(np.abs(np.array(arr, dtype=float) - float(val))))
            except (TypeError, ValueError):
                idx = arr.index(val)
            region[col] = slice(idx, idx + 1)

        # expand dims and assign coords to get the right shape for the region write
        for col in batch_axes:
            val = round(dustrate_mean, 3) if (col == "dustrate_ton_ha_yr" and dustrate_mean is not None) else row[col]
            tmpds = tmpds.assign_coords({col: (col, [val])})
            tmpds = tmpds.assign({v: tmpds[v].expand_dims(col) for v in tmpds.data_vars})

        # transpose each variable to match the zarr store's canonical dim order
        tmpds = tmpds.assign({
            v: tmpds[v].transpose(*[d for d in canonical_dims if d in tmpds[v].dims])
            for v in tmpds.data_vars
        })

        tmpds.to_zarr(zarr_path, region=region)
        tmpds.close()
        del tmpds
        gc.collect()

    return xr.open_zarr(zarr_path)


def prof_batchprocess_singlevar(
    dfin: pd.DataFrame,
    outdir: str,
    batch_axes: list,
    filename_base: str,
    dustsp: str=None,
    subdir: str = "postproc_profs",
    zarr_path: str = None,
) -> xr.Dataset:
    '''
    read in individual nc files from the postproc_profs directory
    and combine them into a single file defined over dimensions
    in batch_axes

    Parameters
    ----------
    dfin : pd.DataFrame
        the batch dataframe that defines each run in the batch
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    batch_axes : list
        list of dimensions for the xr dataset. must correspond to columns in dfin.
        the profile equivalent of 'dfin_cols_to_keep' in flux postprocessing land
    filename_base : str
        name of the .nc file excluding '.nc'
    dustsp : str
        name of the dust species (only required if we're reading in the dust
        data to get the integrated dust flux)
    subdir : str
        name of the subdirectory where the .nc files are stored
    zarr_path : str, optional
        If provided, use the two-pass memory-efficient implementation that writes
        one run at a time to this zarr store instead of accumulating all runs in
        a ds_list before merging. Peak RAM is ~1 run instead of all runs combined.
        The store is created fresh; pass a path that does not yet exist.
        Returns xr.open_zarr(zarr_path) (lazy).

    Returns
    -------
    xr.Dataset
        combination of all xarray datasets in the batch
    '''
    if zarr_path is not None:
        return _prof_batchprocess_singlevar_zarr(
            dfin, outdir, batch_axes, filename_base, dustsp, zarr_path
        )

    # loop through all runs
    ds_list = []

    for index, row in dfin.iterrows():
        # read in this nc file
        tmpds = read_profile_nc(outdir, filename_base, row)
        # skip if no profile was returned
        if tmpds is None:
            continue

        # --- read in rockdiss to get true int dustrate 
        dustrate_mean, dustdf = get_timemean_dustrate(outdir, row, dustsp, return_df=True)
        # --- add a rock application variable 
        # (unless the file is "soil_ph" because the field time and 
        #  lab time dims are different (!))
        if filename_base != "soil_ph":
            dustds = xr.Dataset.from_dataframe(dustdf.set_index("time"))
            # get in the same time coordinates
            dustds = dustds.reindex(time=tmpds["time"])
            # merge 
            tmpds = xr.merge([tmpds, dustds['int_dust_ton_ha_yr']])
        # --- assign coords
        for col in batch_axes:
            if (col == "dustrate_ton_ha_yr") and (dustrate_mean is not None):
                tmpds = tmpds.assign_coords({col: (col, [dustrate_mean])})   # add coord
            else:
                tmpds = tmpds.assign_coords({col: (col, [row[col]])})   # add coord
            tmpds = tmpds.assign({var: tmpds[var].expand_dims(col) for var in tmpds.data_vars})    # assign coord to all data vars

        # --- add to the output list
        ds_list.append(tmpds)

    # --- remove any datasets whose time coords are off 
    #     (indicating failed or incomplete run)
    ds_list = filter_mismatched_time_coords_fuzzy(ds_list, filename_base)
    dsout = xr.merge(ds_list)

    # return result
    return dsout


def prof_batchprocess_singlevar_dictpattern(
    rundicts: dict,
    outdir: str,
    batch_axes: list,
    filename_base: str,
    dustsp: str=None,
    subdir: str = "postproc_profs",
    rockAmend_flag: bool=True,
) -> xr.Dataset:
    '''
    read in individual nc files from the postproc_profs directory
    and combine them into a single file defined over dimensions 
    in batch_axes. Works for batches of runs that use a single rundict
    as a json per run rather than a batch.csv and default dict. 

    Parameters
    ----------
    rundicts : dict of dicts
        dict where each element is a rundict for a given run, and the key is the name of the 
        run. 
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    batch_axes : list
        list of dimensions for the xr dataset. must correspond to columns in dfin.
        the profile equivalent of 'dfin_cols_to_keep' in flux postprocessing land
    filename_base : str
        name of the .nc file excluding '.nc'
    dustsp : str
        name of the dust species (only required if we're reading in the dust 
        data to get the integrated dust flux)
    subdir : str
        name of the subdirectory where the .nc files are stored
    rockAmend_flag : bool
        True if you added dust / want to compute dustrates in all profiles. False otherwise. 
    
    Returns
    -------
    xr.Dataset
        combination of all xarray datasets in the batch
    '''

    # loop through all runs
    ds_list = []

    for key, rundict in rundicts.items():
        # read in this nc file
        tmpds = read_profile_nc(outdir, filename_base, rundict)
        # skip if no profile was returned
        if tmpds is None:
            continue

        if rockAmend_flag:
            # --- read in rockdiss to get true int dustrate 
            dustrate_mean, dustdf = get_timemean_dustrate(outdir, rundict, dustsp, return_df=True)
            # --- add a rock application variable 
            # (unless the file is "soil_ph" because the field time and 
            #  lab time dims are different (!))
            if filename_base != "soil_ph":
                dustds = xr.Dataset.from_dataframe(dustdf.set_index("time"))
                # get in the same time coordinates
                dustds = dustds.reindex(time=tmpds["time"])
                # merge 
                tmpds = xr.merge([tmpds, dustds['int_dust_ton_ha_yr']])
        
        # --- assign coords
        for col in batch_axes:
            if (col == "dustrate_ton_ha_yr") and (dustrate_mean is not None):
                tmpds = tmpds.assign_coords({col: (col, [dustrate_mean])})   # add coord
            else:
                tmpds = tmpds.assign_coords({col: (col, [rundict[col]])})   # add coord
            tmpds = tmpds.assign({var: tmpds[var].expand_dims(col) for var in tmpds.data_vars})    # assign coord to all data vars

        # --- add to the output list
        ds_list.append(tmpds)

    # --- remove any datasets whose time coords are off 
    #     (indicating failed or incomplete run)
    ds_list = filter_mismatched_time_coords_fuzzy(ds_list, filename_base)
    dsout = xr.merge(ds_list)

    # return result
    return dsout


# dictionary for which files to process in the batch profile functions
proc_dict_default = {
    "adsorbed": False,
    "adsorbed_percCEC": False,
    "adsorbed_ppm": False,
    "aqueous": True,
    "aqueous_total": False,
    "bulksoil": True,
    "exchange_total": False,
    "gas": True,
    "rate": False,
    "soil_ph": True,
    "solid": False,
    "solid_sp_saturation": False, 
    "solid_volumePercent": False,
    "solid_weightPercent": True,
    "specific_surface_area": False,
    "surface_area": False,
}



def prof_batchprocess_allvars(
    outdir: str,
    dustsp: str,
    dfin: pd.DataFrame,
    batch_axes: list,
    proc_dict: dict=proc_dict_default,
    subdir: str="postproc_profs",
    print_loop_updates: bool=False,
    save_path: str=None,
    filename_suffix: str="batch",
) -> dict:
    '''
    Wrapper around prof_batchprocess_singlevar to repeat that function for all
    variables in the proc_dict. Resulting datasets are stored in a
    dictionary to be saved later, or written to zarr stores immediately
    to keep peak memory to one variable at a time.

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    dustsp : str
        name of the dust species (only required if we're reading in the dust
        data to get the integrated dust flux)
    dfin : pd.DataFrame
        the batch dataframe that defines each run in the batch
    batch_axes : list
        list of dimensions for the xr dataset. must correspond to columns in dfin.
        the profile equivalent of 'dfin_cols_to_keep' in flux postprocessing land
    proc_dict : dict
        keys should be `filename_base` value inputs to prof_batchprocess_singlevar.
        Values should be True | False. True means the prof_batchprocess_singlevar
        function will be run
    subdir : str
        name of the subdirectory where the .nc files are stored (not used for now
        because there's a default set in prof_batchprocess_singlevar)
    print_loop_updates : bool
        True to print out the profile variable currently being compiled.
    save_path : str, optional
        If provided, each variable's dataset is written to a zarr store under
        this directory immediately after it is computed, then freed from memory.
        This keeps peak RAM usage to roughly one variable at a time instead of
        accumulating all variables simultaneously. The zarr store names match
        the convention used by save_batch_postproc_profOnly:
        ``{key}_{dustsp}_{filename_suffix}.zarr``.
        When set, the returned dict maps variable names to zarr store paths
        (strings) rather than in-memory datasets.
    filename_suffix : str
        Suffix appended to zarr store filenames when save_path is set.
        Default ``"batch"`` matches save_batch_postproc_profOnly's default.

    Returns
    -------
    dict
        When save_path is None: dictionary where each key is a filename_base
        from proc_dict elements with values == True, and each value is the
        corresponding xarray Dataset held in memory (original behaviour).
        When save_path is set: dictionary mapping each key to the path of
        the zarr store written to disk/S3, so results can be read back
        individually with xr.open_zarr(path).
    '''
    # empty dict to hold results (datasets or zarr paths)
    outdict = {}
    # update user
    print(f"compiling profile data for {dustsp}")

    # loop through the process dict
    for key, runme in proc_dict.items():
        if print_loop_updates:
            print(key)
        if runme:
            if save_path is not None:
                # Pass the final zarr path directly into prof_batchprocess_singlevar
                # so it writes one run at a time (O(1-run) peak RAM per variable).
                zarr_path = os.path.join(save_path, f"{key}_{dustsp}_{filename_suffix}.zarr")
                ds = prof_batchprocess_singlevar(
                    dfin, outdir, batch_axes, key, dustsp, zarr_path=zarr_path
                )
                if not ds.variables:
                    print(f"Warning: batch profile processing {key} returned no results. Check for typos?")
                    continue
                outdict[key] = zarr_path
                del ds
            else:
                ds = prof_batchprocess_singlevar(dfin, outdir, batch_axes, key, dustsp)
                if not ds.variables:
                    print(f"Warning: batch profile processing {key} returned no results. Check for typos?")
                    continue
                outdict[key] = ds
    # return result
    return outdict


def prof_batchprocess_allvars_dictpattern(
    outdir: str,
    dustsp: str,
    rundicts: dict,
    batch_axes: list,
    proc_dict: dict=proc_dict_default,
    subdir: str="postproc_profs",
    print_loop_updates: bool=False,
    rockAmend_flag: bool = True,

) -> dict:
    '''
    Wrapper around prof_batchprocess_singlevar to repeat that function for all
    variables in the proc_dict. Works for batches that use run-specific dictionaries
    rather than a batch.csv and default dict. 

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    dustsp : str
        name of the dust species (only required if we're reading in the dust 
        data to get the integrated dust flux)
    rundicts : dict of dicts
        dict where each element is a rundict for a given run, and the key is the name of the 
        run. 
    batch_axes : list
        list of dimensions for the xr dataset. must correspond to columns in dfin.
        the profile equivalent of 'dfin_cols_to_keep' in flux postprocessing land
    proc_dict : dict
        keys should be `filename_base` value inputs to prof_batchprocess_singlevar.
        Values should be True | False. True means the prof_batchprocess_singlevar 
        function will be run
    subdir : str
        name of the subdirectory where the .nc files are stored (not used for now 
        because there's a default set in prof_batchprocess_singlevar)
    print_loop_updates : bool
        True to print out the profile variable currently being compiled. 
    rockAmend_flag : bool
        True if you added dust / want to compute dustrates in all profiles. False otherwise. 
    
    Returns
    -------
    dict
        dictionary where each key is a filename_base from proc_dict elements
        with values == True. Each value in the output dict is an xarray 
        dataset with the batch_postprocess result. 
    '''
    # empty dict to hold results
    outdict = {}
    # update user 
    print(f"compiling profile data for {dustsp}")

    # loop through the process dict
    for key, runme in proc_dict.items():
        if print_loop_updates:
            print(key)
        if runme:
            ds = prof_batchprocess_singlevar_dictpattern(rundicts, outdir, batch_axes, key, dustsp, rockAmend_flag=rockAmend_flag)
            # make sure it's not empty
            if not ds.variables: # continue to next iter if it is empty
                print(f"Warning: batch profile processing {key} returned no results. Check for typos?")
                continue
            # ori suggests: print(outdict)
            outdict[key] = ds
            # ds.to_netcdf()
    # return result
    return outdict


def save_batch_postproc_profOnly(
    dsdict: dict,
    dustsp_tag: str,
    filename_suffix: str="batch",
    save_directory: str=None,
    base_path: str=None,
    base_dir_name: str=None,
    use_zarr: bool=True,
):
    '''
    Save the individual files in the batch postprocess profile dictionary.
    Either save_directory or base_path AND base_dir_name must be defined. 
    If save_directory is defined then it is used, otherwise a new 
    directory is created from the base_path and base_dir_name.

    Parameters
    ----------
    dsdict : dict
        dictionary where each key is the filename_base and each element 
        is a dataset returned by the prof_batchprocess_singlevar function
    dustsp_tag : str
        tag to differentiate feedstocks that get saved in the same dir (just dustsp usually)
    filename_suffix : str
        a suffix for the .nc filenames (e.g., to idnicate it's a batch file)
    save_directory : str
        where to save the results (if this is none, then the next two vars
        should be defined)
    base_path : str
        where to create the new directory (if none, save_directory must be defined)
    base_dir_name : 
        name of the new directory to create (if none, save_directory must be defined)
    
    Returns
    -------
    '''
    # create the new directory or set the defined one
    if save_directory is None:
        outdir = create_unique_directory(base_path, base_dir_name)
    else:
        outdir = save_directory

    # --- save results
    if use_zarr: 
        for key, ds in dsdict.items():
            savename = f"{key}_{dustsp_tag}_{filename_suffix}.zarr"
            ds.to_zarr(os.path.join(outdir, savename))
    else:
        for key, ds in dsdict.items():
            savename = f"{key}_{dustsp_tag}_{filename_suffix}.nc"

            # check for AWS
            if outdir.startswith("s3://"):
                # NOTE: I couldn't get the netcdf files to save to aws
                # with the default netcdf4 engine. As a workaround, I'm 
                # saving the file locally in a tempfile, then moving it 
                # to aws and deleting it locally :[]
                import tempfile
                with tempfile.NamedTemporaryFile(delete=False) as tmp_file1:
                    # save the dataset to the temporary file
                    tmp_file1_path = tmp_file1.name
                    ds.to_netcdf(tmp_file1_path)
                    
                    # upload the file to S3 using s3fs
                    import s3fs
                    fs = s3fs.S3FileSystem()
                    fs.put(tmp_file1_path, os.path.join(outdir, savename))
                
            else:
                ds.to_netcdf(os.path.join(outdir, savename))
    
    # return the directory it's saved in
    return outdir




# %% 


# FUNCTION create the dict we'll use to save results
def create_save_dict(
    cdr_calc_list: list,
    group_vars: list,
    csv_fn_sil: str,
    csv_fn_cc: str,
    multiyear_sil: bool, 
    multiyear_cc: bool, 
    time_horizon: float, 
    cf_apprate_fixed: float,
    cf_dustrad_fixed: float, 
    p80_inputsil: float,
    p80_inputcc: float, 
    truck_kmsil: float,
    truck_kmcc: float,
    barge_kmsil: float, 
    barge_kmcc: float, 
    barge_diesel_kmsil: float, 
    barge_diesel_kmcc: float, 
    Efactor_orgsil: float,
    Efactor_orgcc: float,
)->dict:
    """
    Combine all inputs into a pre-formatted dictionary that 
    we'll save as a .res file. 
    """
    savedict = {
        "Setup": {
            # "site": thissite,
            "cdvars": cdr_calc_list,
            "group_vars": group_vars,
            "csv_sil": csv_fn_sil,
            "csv_cc": csv_fn_cc,
            "multiyear_sil": multiyear_sil,
            "multiyear_cc": multiyear_cc
        },
        "Removal accounting": {
            "time_horizon": time_horizon,
            "cf_apprate_fixed": cf_apprate_fixed,
            "cf_dustrad_fixed": cf_dustrad_fixed
        },
        "Emissions (silicate)": {
            "p80_input": p80_inputsil,
            "truck_km": truck_kmsil,
            "barge_km": barge_kmsil,
            "barge_diesel_km": barge_diesel_kmsil,
            "Efactor_org": Efactor_orgsil
        },
        "Emissions (calcite)": {
            "p80_input": p80_inputcc,
            "truck_km": truck_kmcc,
            "barge_km": barge_kmcc,
            "barge_diesel_km": barge_diesel_kmcc,
            "Efactor_org": Efactor_orgcc
        }
    }

    return savedict



def create_save_dict_oneFS_noEmiss(
    EXP_SET: str,
    cdr_calc_list: list,
    group_vars: list,
    collapse_cols: list,
    csv_fn: str,
    multiyear: bool, 
    time_horizon: float, 
    outdir: str,
    batchdir: str,
    EXP_GROUP: str=None,
)->dict:
    """
    Combine all inputs into a pre-formatted dictionary that 
    we'll save as a .res file. 
    """
    savedict = {
        "Setup": {
            # "site": thissite,
            "EXP_GROUP": EXP_GROUP,
            "EXP_SET": EXP_SET,
            "cdvars": cdr_calc_list,
            "group_vars": group_vars,
            "collapse_cols": collapse_cols,
            "csv": csv_fn,
            "multiyear": multiyear,
            "outdir": outdir,
            "batchdir": batchdir,
        },
        "Removal accounting": {
            "time_horizon": time_horizon,
        },
        
    }

    return savedict



def create_save_dict_oneFS_emiss(
    EXP_SET: str,
    cdr_calc_list: list,
    group_vars: list,
    collapse_cols: list,
    csv_fn: str,
    multiyear: bool, 
    time_horizon: float, 
    outdir: str,
    batchdir: str,
    p80_input: float,
    Efactor_org: str,
    EXP_GROUP: str=None,
)->dict:
    """
    Combine all inputs into a pre-formatted dictionary that 
    we'll save as a .res file. 
    """
    savedict = {
        "Setup": {
            # "site": thissite,
            "EXP_GROUP": EXP_GROUP,
            "EXP_SET": EXP_SET,
            "cdvars": cdr_calc_list,
            "group_vars": group_vars,
            "collapse_cols": collapse_cols,
            "p80_input": p80_input,
            "Efactor_org": Efactor_org,
            "csv": csv_fn,
            "multiyear": multiyear,
            "outdir": outdir,
            "batchdir": batchdir,
        },
        "Removal accounting": {
            "time_horizon": time_horizon,
        },
        
    }

    return savedict


def save_batch_postproc(
    savedict: dict, 
    base_path: str, 
    base_dir_name: str,
    fname_res: str,
    cdr_dict_sum_sil: dict,
    cdr_dict_sum_cc: dict,
    ds_sil: xr.Dataset,
    ds_cc: xr.Dataset,
    ds_anom: xr.Dataset,
    save_directory: str=None,
) -> str:
    """
    save the series of files generated when postprocessing the batch 
    results

    Parameters
    ----------
    savedict : dict
        dictionary of variables we will save as a resource file
    base_path : str
        directory where we will create the new directory to hold 
        the results. only used if save_directory is None
    base_dir_name : str
        name of the directory we want to create (the create_unique_directory
        function will make sure nothing is overwritten). only used if 
        save_directory is None
    fname_res : str
        name of the resource file that will be created from savedict
    cdr_dict_sum_sil : dict
        dictionary of pandas dataframes with CDR outputs — these are 
        ultimately input to the emissions calculation step (we save them
        to make it easier to play with other emissions scenarios after
        this step)
    cdr_dict_sum_cc : dict
        as above but for cc, not sil feedstock
    ds_sil : xr.Dataset 
        dataset of the silicate results vs loss_percent and group_vars
    ds_cc : xr.Dataset 
        dataset of the silicate results vs loss_percent and group_vars
    ds_anom : xr.Dataset 
        dataset of the silicate minus cc anomaly for different CDR 
        accounting approaches
    save_directory : str
        path to the directory to save outputs. If this is None, then 
        a unique directory will be created based on base_path and 
        base_dir_name
    
    Returns
    -------
    str
        the path for the save directory
    """
    # create the new directory or set the defined one
    if save_directory is None:
        outdir = create_unique_directory(base_path, base_dir_name)
    else:
        outdir = save_directory

    # [2] save the dict as a resource file
    save_variables_to_file(os.path.join(outdir, fname_res), savedict)

    # [3] save the pandas dataframes as a single pickle file
    savedf_sil_fn = os.path.join(outdir, "cdr_dfs_sil.pkl")
    savedf_cc_fn = os.path.join(outdir, "cdr_dfs_cc.pkl")

    # ----------------------------
    # check for S3
    if savedf_sil_fn.startswith("s3://"): # then bring it in from s3
        import fsspec
        with fsspec.open(savedf_sil_fn, 'wb') as f:
            pickle.dump(cdr_dict_sum_sil, f)
    else:
        with open(savedf_sil_fn, 'wb') as f:
            pickle.dump(cdr_dict_sum_sil, f)
    if savedf_cc_fn.startswith("s3://"): # then bring it in from s3
        import fsspec
        with fsspec.open(savedf_cc_fn, 'wb') as f:
            pickle.dump(cdr_dict_sum_cc, f)
    else:
        with open(savedf_cc_fn, 'wb') as f:
            pickle.dump(cdr_dict_sum_cc, f)
    # ----------------------------
    # NOTE: you can open the dfs like this:
    # with open('dataframes.pkl', 'rb') as f:
        # loaded_dfs = pickle.load(f) 

    # [4] save xr.Datasets (should work w/ aws)
    if outdir.startswith("s3://"):
        # NOTE: I couldn't get the netcdf files to save to aws
        # with the default netcdf4 engine. As a workaround, I'm 
        # saving the file locally in a tempfile, then moving it 
        # to aws and deleting it locally :[]
        import tempfile
        with tempfile.NamedTemporaryFile(delete=False) as tmp_file1, \
            tempfile.NamedTemporaryFile(delete=False) as tmp_file2, \
            tempfile.NamedTemporaryFile(delete=False) as tmp_file3:
            
            # save the dataset to the temporary file
            tmp_file1_path = tmp_file1.name
            tmp_file2_path = tmp_file2.name
            tmp_file3_path = tmp_file3.name

            ds_sil.to_netcdf(tmp_file1_path)
            ds_cc.to_netcdf(tmp_file2_path)
            ds_anom.to_netcdf(tmp_file3_path)

            # upload the file to S3 using s3fs
            import s3fs
            fs = s3fs.S3FileSystem()
            fs.put(tmp_file1_path, os.path.join(outdir, "ds_sil.nc"))
            fs.put(tmp_file2_path, os.path.join(outdir, "ds_cc.nc"))
            fs.put(tmp_file3_path, os.path.join(outdir, "ds_anom.nc"))
        
    else:
        ds_sil.to_netcdf(os.path.join(outdir, "ds_sil.nc"))
        ds_cc.to_netcdf(os.path.join(outdir, "ds_cc.nc"))
        ds_anom.to_netcdf(os.path.join(outdir, "ds_anom.nc"))

    # return the path to the save directory
    return outdir


def save_batch_postproc_oneSet( # repeated for when we only have one set of simulations
    savedict: dict, 
    base_path: str, 
    base_dir_name: str,
    fname_res: str,
    df_dict: dict,
    ds_dict: dict,
    save_directory: str=None,
) -> str:
    """
    save the series of files generated when postprocessing the batch 
    results

    Parameters
    ----------
    savedict : dict
        dictionary of variables we will save as a resource file
    base_path : str
        directory where we will create the new directory to hold 
        the results. only used if save_directory is None
    base_dir_name : str
        name of the directory we want to create (the create_unique_directory
        function will make sure nothing is overwritten). only used if 
        save_directory is None
    fname_res : str
        name of the resource file that will be created from savedict
    df_dict : dict
        dictionary of dictionaries with pd.Dataframes for each run
    ds_dict : xr.Dataset 
        dictionary of xr.Datasets to save where key is the filename and 
        value is the xr.Dataset
    save_directory : str
        path to the directory to save outputs. If this is None, then 
        a unique directory will be created based on base_path and 
        base_dir_name
    
    Returns
    -------
    str
        the path for the save directory
    """
    # create the new directory or set the defined one
    if save_directory is None:
        outdir = create_unique_directory(base_path, base_dir_name)
    else:
        outdir = save_directory

    # [2] save the dict as a resource file
    save_variables_to_file(os.path.join(outdir, fname_res), savedict)

    # [3] save the pandas dataframes as a single pickle file
    for key, df in df_dict.items():
        savedf_fn = os.path.join(outdir, key)
        # ----------------------------
        # check for S3
        if savedf_fn.startswith("s3://"): # then bring it in from s3
            import fsspec
            with fsspec.open(savedf_fn, 'wb') as f:
                pickle.dump(df, f)
        else:
            with open(savedf_fn, 'wb') as f:
                pickle.dump(df, f)
        # ----------------------------
        # NOTE: you can open the dfs like this:
        # with open('dataframes.pkl', 'rb') as f:
            # loaded_dfs = pickle.load(f) 

    # [4] save datasets to zarr
    for key, ds in ds_dict.items():
        ds.to_zarr(os.path.join(outdir, key))

    # return the path to the save directory
    return outdir


def matches_pattern(
    fns: str, 
    mylist: list,
)->bool:
    return any(fnmatch.fnmatch(fns, pattern) for pattern in mylist)



def create_unique_directory(
    base_path: str, 
    base_dir_name: str
)-> str:
    """
    Create a unique directory name to save the results of the batch 
    postprocess experiments

    Parameters
    ----------
    base_path : str
        path to the directory you want to make (for example:
        "/here/is/my/base/path")
    base_dir_name : str
        name of the directory you want to make (without the numbered 
        suffix) (for example: "/here/is/my/base/path/mydir"). The 
        function will then append a suffix to give you something like 
        "mydir_001", ensuring it's unique in the base_path.

    Results
    -------
    str
        returns the full path and dir name of the directory that 
        was created within the function
    """
    # assume we're not using s3 until we check the dir later
    s3_path = False
    
    i = 1
    while True:
        dir_name = os.path.join(base_path, f"{base_dir_name}_{i:03}")
        # ----------------------------
        # check for S3
        if dir_name.startswith("s3://"): # then bring it in from s3
            import s3fs
            fs = s3fs.S3FileSystem()
            s3_path_check = True
        # ----------------------------
            if not fs.exists(dir_name):
                fs.makedirs(dir_name)
                print(f"Directory created: {dir_name}")
                break
        elif not os.path.exists(dir_name):
            os.makedirs(dir_name)
            print(f"Directory created: {dir_name}")
            break
        i += 1
    return dir_name


def save_variables_to_file(
    filename: str, 
    data: dict,
):
    """
    Take in a dictionary and save the results to a .txt file

    Parameters
    ----------
    filename : str
        path and name of the file you want to save
    data : dict
        dictionary of inputs used for the batch processing step

    Returns
    -------

    """
    # assume we're not using aws until we check later
    s3_path = False
    
    # ---
    # check for aws
    if filename.startswith("s3://"): # then bring it in from s3
        import fsspec
        with fsspec.open(filename, 'w') as file:
            for label, variables in data.items():
                file.write(f"{label}:\n")
                for var_name, value in variables.items():
                    file.write(f"  {var_name}: {value}\n")
                file.write("\n")  # Adds a blank line between tests
    # ---
    else:
        with open(filename, 'w') as file:
            for label, variables in data.items():
                file.write(f"{label}:\n")
                for var_name, value in variables.items():
                    file.write(f"  {var_name}: {value}\n")
                file.write("\n")  # Adds a blank line between tests


# quick control run check
def find_control_runs(
        dfin: pd.DataFrame,
        ctrl_conditions: dict,
)->pd.DataFrame:
    """
    Identify control runs based on a set of control 
    conditions in the dictionary

    Parameters
    ----------
    dfin : pd.DataFrame
        the batch.csv file
    ctrl_conditions : dict
        dictionary identifying control conditions
        example: = {"dustrate": 0}
    
    Returns
    -------
    pd.DataFrame
    """
    mask = pd.Series(True, index=dfin.index)
    for col, val in ctrl_conditions.items():
        mask &= (dfin[col] == val)
    
    dfin["ctrl_run"] = mask
    # return result
    return dfin



# ----------------------------------------------------------------------------------
# [ ADD FUNCTIONS FOR GCAM POSTPROCESSING ]
# 
def _int_to_flux(ds_int: xr.Dataset) -> xr.Dataset:
    """
    Differentiate time-integrated SCEPTER outputs to recover instantaneous flux rates.
    Only variables with a 'time' dimension are differentiated; others pass through unchanged.
    Variables whose names start with 'fraction_' are always passed through unchanged —
    they represent dimensionless ratios that should not be differentiated.

    Also renames 'tonHa_' → 'tonHaYr_' (e.g. co2pot_adv_tonHa_camg → co2pot_adv_tonHaYr_camg)
    and appends ' yr-1' to any existing units attribute (for differentiated vars only), so
    the result is a drop-in replacement for a ds_trans slice.
    """
    t = ds_int.time.values

    # Rename cumulative 'tonHa' vars → 'tonHaYr' to match ds_trans naming convention;
    # the derivative d(t/ha)/dt has units t/ha/yr, consistent with the 'Yr' suffix.
    rename_map = {
        v: v.replace('tonHa_', 'tonHaYr_')
        for v in ds_int.data_vars
        if 'tonHa_' in v and 'tonHaYr' not in v
    }
    ds_work = ds_int.rename(rename_map) if rename_map else ds_int

    updates = {}
    for vname in ds_work.data_vars:
        arr = ds_work[vname]
        if 'time' in arr.dims and not vname.startswith('fraction_'):
            new_attrs = dict(arr.attrs)
            if new_attrs.get('units', ''):
                new_attrs['units'] = new_attrs['units'].rstrip() + ' yr-1'
            updates[vname] = xr.DataArray(
                np.gradient(arr.values, t),
                coords=arr.coords, dims=arr.dims, attrs=new_attrs,
            )
    return ds_work.assign(updates)

def gcam_postprocess_region_constApp(
    ds_rep:                  xr.Dataset,
    df_rock:                 pd.DataFrame,
    dsag:                    xr.Dataset,
    region_idx:              int,
    varlist_deploy:          list[str],
    varlist_region:          list[str],
    varlist_region_cumfrac:  list[str],
    varlist_rco2_correction: list[str],
    rco2_correction_factor:  float,
    cdr_var:                 str   = 'cdr_adv',
    t_start:                 float = 0.2,
    basalt_cdrpot:           float = 0.3,
) -> xr.Dataset:
    """
    Scale a representative SCEPTER run to regional CDR totals.

    Parameters
    ----------
    ds_rep : xr.Dataset
        Representative run, already isel'd to one (site, apprate_base, dustrad, cec_factor).
    df_rock : pd.DataFrame
        Rock flux dataframe for this apprate (all regions); filtered internally by region_idx.
    dsag : xr.Dataset
        Cropland dataset for this apprate (all regions); filtered internally by region_idx.
    region_idx : int
        GCAM reg_id used to filter df_rock and dsag. Must match the 'reg_id' column in df_rock.
    varlist_deploy : list of str
        Variables to build as (time, deployment) lag matrices.
    varlist_region : list of str
        Subset of varlist_deploy to aggregate as area-scaled regional totals.
    varlist_region_cumfrac : list of str
        Subset of varlist_deploy to aggregate as cumulative-rock-weighted regional means.
    varlist_rco2_correction : list of str
        Variables to scale by rco2_correction_factor before deployment mapping and regional
        aggregation, to account for GCAM's lower CDR potential per tonne of rock vs. SCEPTER.
    rco2_correction_factor : float
        Multiplicative correction applied to variables in varlist_rco2_correction.
        Defined as RCO2_gcam / RCO2_scepter.
    cdr_var : str
        Primary CDR variable (default 'cdr_adv'); always written as {cdr_var}_region.
    t_start : float
        Time (yr) after which to begin interpolation, skipping the initial numerical blip.
    basalt_cdrpot : float
        g CO2 per g basalt, for converting GCAM CDR to rock mass (default 0.3).

    Returns
    -------
    xr.Dataset
        All deployment-level and region-level variables on a common annual time grid.

    Raises
    ------
    ValueError
        If region_idx is not present in df_rock.
    """

    # ── Step 1: Setup — time grid and base parameters ─────────────────────────────────────────
    time_raw    = ds_rep.time.values
    cdr_raw     = ds_rep[cdr_var].values.copy()
    R0          = float(ds_rep.apprate_base.values)
    region_name = str(ds_rep.site.values)

    n_years   = int(time_raw[-1] - t_start)
    time_vals = np.arange(n_years + 1, dtype=float)
    cdr_vals  = np.interp(time_vals + t_start, time_raw, cdr_raw)
    n_time    = len(time_vals)

    extra_vals = {
        vname: np.interp(time_vals + t_start, time_raw, ds_rep[vname].values.copy())
        for vname in varlist_deploy
    }

    # Apply RCO2 correction before building deployment lag matrices
    for vname in varlist_rco2_correction:
        if vname in extra_vals:
            extra_vals[vname] = extra_vals[vname] * rco2_correction_factor
    if cdr_var in varlist_rco2_correction:
        cdr_vals = cdr_vals * rco2_correction_factor

    # ── Step 2: Build lower-triangular lag matrices ────────────────────────────────────────────
    # mat[j, i] = vals[j - i]  for j >= i, else NaN
    deploy_coords = np.arange(n_time)

    def _lag_matrix(vals):
        mat = np.full((n_time, n_time), np.nan)
        for i in range(n_time):
            mat[i:, i] = vals[:n_time - i]
        return mat

    CDR_flux_2d   = _lag_matrix(cdr_vals)
    extra_flux_2d = {vname: _lag_matrix(v) for vname, v in extra_vals.items()}

    # ── Step 3: R_flux lag matrix + initial ds_out ────────────────────────────────────────────
    df_region = (df_rock[df_rock['reg_id'] == region_idx]
                 .sort_values('time').reset_index(drop=True))

    if df_region.empty:
        raise ValueError(f"reg_id={region_idx} ({region_name}) not found in df_rock")

    gcam_t0    = int(df_region['time'].min())
    gcam_years = df_region['time'].values.astype(float)

    r_gcam = df_region['effective_rate_t_ha'].values
    r_vals = np.interp(time_vals, gcam_years - gcam_t0, r_gcam,
                       left=r_gcam[0], right=r_gcam[-1])
    R_flux_2d = _lag_matrix(r_vals)

    coords_td = {'time': time_vals, 'deployment': deploy_coords}
    coords_t  = {'time': time_vals}

    _rco2_note = ' (scaled by rco2_correction_factor = RCO2_gcam / RCO2_scepter)'

    ds_out = xr.Dataset(
        {
            f'{cdr_var}_deployment': xr.DataArray(
                CDR_flux_2d, coords=coords_td, dims=['time', 'deployment'],
                attrs={
                    'units': 't/ha/yr',
                    'description': ('CDR flux per unit area for each deployment'
                                    + (_rco2_note if cdr_var in varlist_rco2_correction else '')),
                },
            ),
            'rockflx_t_ha_yr': xr.DataArray(
                R_flux_2d, coords=coords_td, dims=['time', 'deployment'],
                attrs={'units': 't/ha/yr', 'description': 'Rock application rate per deployment'},
            ),
            'rco2_correction_factor': xr.DataArray(
                rco2_correction_factor,
                attrs={
                    'units': '-',
                    'description': 'RCO2 correction factor applied to varlist_rco2_correction variables (= RCO2_gcam / RCO2_scepter)',
                },
            ),
            **{
                f'{vname}_deployment': xr.DataArray(
                    extra_flux_2d[vname], coords=coords_td, dims=['time', 'deployment'],
                    attrs={
                        'units': str(ds_rep[vname].attrs.get('units', '')),
                        'description': (f'{vname} per unit area for each deployment'
                                        + (_rco2_note if vname in varlist_rco2_correction else '')),
                    },
                )
                for vname in varlist_deploy
            },
        }
    )
    ds_out.attrs['region'] = region_name

    # ── Step 4: M_GCAM ────────────────────────────────────────────────────────────────────────
    gcam_rock_t = df_region['gcam_gross_cdr_mtCO2'].values / basalt_cdrpot * 1e6
    m_gcam_vals = np.interp(time_vals, gcam_years - gcam_t0, gcam_rock_t,
                            left=0.0, right=gcam_rock_t[-1])

    # ── Step 5: Iterative area calculation ────────────────────────────────────────────────
    A_D_2d    = np.full((n_time, n_time), np.nan)
    M_D_2d    = np.full((n_time, n_time), np.nan)
    mu_arr    = np.zeros(n_time)
    M_new_arr = np.zeros(n_time)

    for j in range(n_time):
        if j > 0:
            A_D_2d[j, :] = A_D_2d[j - 1, :]

        M_D_2d[j, :j]  = R_flux_2d[j, :j] * A_D_2d[j, :j]
        mu_arr[j]       = np.nansum(M_D_2d[j, :j])
        M_new_arr[j]    = m_gcam_vals[j] - mu_arr[j]

        if M_new_arr[j] >= 0:
            A_D_2d[j, j] = M_new_arr[j] / r_vals[0]
            M_D_2d[j, j] = M_new_arr[j]
        else:
            excess = -M_new_arr[j]
            for i_trim in range(j - 1, -1, -1):
                if not (A_D_2d[j, i_trim] > 0):
                    continue
                removable = M_D_2d[j, i_trim]
                if removable <= excess:
                    A_D_2d[j, i_trim] = 0.0
                    M_D_2d[j, i_trim] = 0.0
                    excess -= removable
                else:
                    A_D_2d[j, i_trim] -= excess / R_flux_2d[j, i_trim]
                    M_D_2d[j, i_trim] -= excess
                    excess = 0.0
                if excess <= 0:
                    break

    # ── Steps 6 + 7: Assemble ds_out ───────────────────────────────────────────────────────────
    mask              = dsag['reg_id'].values == region_idx
    total_cropland_ha = float(dsag['cropland_area_by_region'].values[mask][0]) / 1e4

    gcam_cdr_t    = df_region['gcam_gross_cdr_mtCO2'].values * 1e6
    gcam_cdr_vals = np.interp(time_vals, gcam_years - gcam_t0, gcam_cdr_t,
                               left=0.0, right=gcam_cdr_t[-1])

    cumM_D_2d    = np.nancumsum(M_D_2d, axis=0)
    cumM_D_total = np.nansum(cumM_D_2d, axis=1)

    A_total_ERW_vals = np.nansum(A_D_2d, axis=1)

    # precompute area-weighted solid-phase dissolution arrays
    # (total_dissolution and int_dust_ton_ha_yr are cumulative quantities in the SCEPTER
    # output, so the area-weighted sum gives total dissolved / applied mass in the region)
    _total_diss_region = np.nansum(extra_flux_2d['total_dissolution'] * A_D_2d, axis=1)
    _int_dust_region   = np.nansum(extra_flux_2d['int_dust_ton_ha_yr'] * A_D_2d, axis=1)

    # cdr_var is always computed explicitly; skip it in the varlist_region loop to avoid collision
    vlist_region_others = [v for v in varlist_region if v != cdr_var]

    ds_out = ds_out.assign({
        'A_D': xr.DataArray(
            A_D_2d, coords=coords_td, dims=['time', 'deployment'],
            attrs={'units': 'ha', 'description': 'Deployment area (NaN = not yet open)'},
        ),
        'M_D': xr.DataArray(
            M_D_2d, coords=coords_td, dims=['time', 'deployment'],
            attrs={'units': 't/yr', 'description': 'Rock mass per deployment (NaN = not yet open)'},
        ),
        'M_GCAM': xr.DataArray(
            m_gcam_vals, coords=coords_t, dims=['time'],
            attrs={'units': 't/yr', 'description': 'Annual rock mass target from GCAM'},
        ),
        'mu': xr.DataArray(
            mu_arr, coords=coords_t, dims=['time'],
            attrs={'units': 't/yr', 'description': 'Rock from existing deployments'},
        ),
        'M_new': xr.DataArray(
            M_new_arr, coords=coords_t, dims=['time'],
            attrs={'units': 't/yr', 'description': 'Rock allocated to new deployment'},
        ),
        'A_total_ERW': xr.DataArray(
            A_total_ERW_vals, coords=coords_t, dims=['time'],
            attrs={'units': 'ha', 'description': 'Total ERW area'},
        ),
        'cropland_ha_total': xr.DataArray(
            total_cropland_ha,
            attrs={'units': 'ha', 'description': 'Total cropland area in the region (time-invariant)'},
        ),
        'cropfrac_ERW': xr.DataArray(
            A_total_ERW_vals / total_cropland_ha, coords=coords_t, dims=['time'],
            attrs={'units': '-', 'description': 'Fraction of regional cropland under ERW'},
        ),
        'gcam_cdr': xr.DataArray(
            gcam_cdr_vals, coords=coords_t, dims=['time'],
            attrs={'units': 'tCO2/yr', 'description': 'GCAM-derived regional CDR'},
        ),
        f'{cdr_var}_region': xr.DataArray(
            np.nansum(CDR_flux_2d * A_D_2d, axis=1),
            coords=coords_t, dims=['time'],
            attrs={
                'units': 't/yr',
                'description': ('Total CDR for the region'
                                + (_rco2_note if cdr_var in varlist_rco2_correction else '')),
            },
        ),
        **{
            f'{vname}_region': xr.DataArray(
                np.nansum(extra_flux_2d[vname] * A_D_2d, axis=1),
                coords=coords_t, dims=['time'],
                attrs={
                    'units': 't/yr',
                    'description': (f'Area-scaled regional total for {vname}'
                                    + (_rco2_note if vname in varlist_rco2_correction else '')),
                },
            )
            for vname in vlist_region_others
        },
        **{
            f'{vname}_region': xr.DataArray(
                np.where(
                    cumM_D_total > 0,
                    np.nansum(np.where(np.isfinite(extra_flux_2d[vname]), extra_flux_2d[vname], np.nan) * cumM_D_2d, axis=1) / np.where(cumM_D_total > 0, cumM_D_total, 1),
                    np.nan,
                ),
                coords=coords_t, dims=['time'],
                attrs={
                    'units': str(ds_rep[vname].attrs.get('units', '')),
                    'description': f'Cumulative-rock-weighted regional mean for {vname}',
                },
            )
            for vname in varlist_region_cumfrac
        },
        'total_dissolution_region': xr.DataArray(
            _total_diss_region,
            coords=coords_t, dims=['time'],
            attrs={
                'units': 'ton',
                'description': 'Area-weighted regional total of cumulative feedstock dissolution: sum(total_dissolution * A_D) across deployments',
            },
        ),
        'int_dust_ton_ha_yr_region': xr.DataArray(
            _int_dust_region,
            coords=coords_t, dims=['time'],
            attrs={
                'units': 'ton',
                'description': 'Area-weighted regional total of cumulative applied feedstock: sum(int_dust_ton_ha_yr * A_D) across deployments',
            },
        ),
        'fraction_total_dissolved_fromIntFlx_region': xr.DataArray(
            np.where(
                _int_dust_region > 0,
                _total_diss_region / np.where(_int_dust_region > 0, _int_dust_region, 1),
                np.nan,
            ),
            coords=coords_t, dims=['time'],
            attrs={
                'units': '',
                'description': 'Area-weighted solid-phase dissolution fraction: sum(total_dissolution * A_D) / sum(int_dust_ton_ha_yr * A_D) across deployments',
            },
        ),
    })

    return ds_out

def gcam_postprocess_all_dims_constApp(
    ds_trans:                xr.Dataset,
    dfrock_dict:             dict[str, pd.DataFrame],
    dsag_dict:               dict[str, xr.Dataset],
    varlist_deploy:          list[str],
    varlist_region:          list[str],
    varlist_region_cumfrac:  list[str],
    varlist_rco2_correction: list[str],
    rco2_correction_factor:  float,
    cdr_var:                 str                    = 'cdr_adv',
    t_start:                 float                  = 0.2,
    basalt_cdrpot:           float                  = 0.3,
    apprate_key_fn:          Callable[[float], str] = lambda v: f'{float(v)}-tHaYr-appRate',
    site_dim:                str                    = 'site',
    region_col:              str                    = 'region_nospaces',
    verbose:                 bool                   = True,
    use_integrated:          bool                   = True,
    ds:                      xr.Dataset             = None,
) -> xr.Dataset:
    """
    Run postprocess_region for every combination of non-time dims in ds_trans,
    then assemble into a single multi-dimensional Dataset.

    Parameters
    ----------
    ds_trans : xr.Dataset
        Full transient dataset. Used for dimension discovery, iteration, and coordinate/metadata
        lookups. Also serves as the ds_rep source when use_integrated=False.
    dfrock_dict : dict[str, pd.DataFrame]
        Rock flux dataframes keyed by apprate name string.
    dsag_dict : dict[str, xr.Dataset]
        Cropland datasets keyed by apprate name string.
    varlist_deploy, varlist_region, varlist_region_cumfrac : list of str
        Passed through to postprocess_region.
    varlist_rco2_correction : list of str
        Variables to scale by rco2_correction_factor before deployment mapping and regional
        aggregation. Passed through to postprocess_region.
    rco2_correction_factor : float
        Multiplicative correction (RCO2_gcam / RCO2_scepter). Passed through.
    cdr_var, t_start, basalt_cdrpot :
        Passed through to postprocess_region.
    apprate_key_fn : Callable[[float], str]
        Maps an apprate coordinate value to a dfrock_dict / dsag_dict key.
    site_dim : str
        Name of the site/region dimension in ds_trans (default 'site').
    region_col : str
        Column in df_rock that matches ds_trans site coordinate values.
    verbose : bool
        Print a progress line for each combination (default True).
    use_integrated : bool
        If True (default), derive ds_rep by differentiating ds (smooth, no aliasing).
        If False, use ds_trans slices directly (original behavior).
    ds : xr.Dataset or None
        Time-integrated SCEPTER dataset. Required when use_integrated=True.
        Must have the same non-time coordinates as ds_trans.

    Returns
    -------
    xr.Dataset
        Combined dataset with dims (time, deployment, <all loop dims>).
    """
    import itertools

    if use_integrated and ds is None:
        raise ValueError("ds must be provided when use_integrated=True")

    loop_dims  = [d for d in ds_trans.dims if d != 'time']
    dim_ranges = [range(ds_trans.sizes[d]) for d in loop_dims]
    n_total    = 1
    for r in dim_ranges:
        n_total *= len(r)

    results = []
    for k, indices in enumerate(itertools.product(*dim_ranges)):
        isel_dict    = dict(zip(loop_dims, indices))
        ds_rep_meta  = ds_trans.isel(**isel_dict)   # always use ds_trans for coord/meta lookups
        region_name  = str(ds_rep_meta[site_dim].values)

        # Use apprate_base coord value whether it's a dim or was scalar-isel'd out;
        # next(iter(dfrock_dict)) is only a safe fallback when apprate_base doesn't exist at all.
        if 'apprate_base' in ds_rep_meta.coords:
            apprate_key = apprate_key_fn(float(ds_rep_meta['apprate_base'].values))
        else:
            apprate_key = next(iter(dfrock_dict))

        df_this     = dfrock_dict[apprate_key]
        region_rows = df_this[df_this[region_col] == region_name]
        if region_rows.empty:
            if verbose:
                coord_str = ', '.join(f'{d}={ds_rep_meta[d].values}' for d in loop_dims)
                print(f'[{k+1}/{n_total}]  {coord_str}  -> skipped (not in df_rock)')
            continue
        region_idx = int(region_rows['reg_id'].iloc[0])

        if verbose:
            coord_str = ', '.join(f'{d}={ds_rep_meta[d].values}' for d in loop_dims)
            print(f'[{k+1}/{n_total}]  {coord_str}')

        # Select and (optionally) differentiate the representative run
        if use_integrated:
            ds_rep = _int_to_flux(ds.isel(**isel_dict))
        else:
            ds_rep = ds_rep_meta

        ds_out = gcam_postprocess_region_constApp(
            ds_rep                  = ds_rep,
            df_rock                 = dfrock_dict[apprate_key],
            dsag                    = dsag_dict[apprate_key],
            region_idx              = region_idx,
            varlist_deploy          = varlist_deploy,
            varlist_region          = varlist_region,
            varlist_region_cumfrac  = varlist_region_cumfrac,
            varlist_rco2_correction = varlist_rco2_correction,
            rco2_correction_factor  = rco2_correction_factor,
            cdr_var                 = cdr_var,
            t_start                 = t_start,
            basalt_cdrpot           = basalt_cdrpot,
        )

        for dim in loop_dims:
            ds_out = ds_out.assign_coords({dim: ds_rep_meta[dim].values})
        ds_out = ds_out.expand_dims(loop_dims)

        results.append(ds_out)

    ds_combined = xr.combine_by_coords(results, combine_attrs='override')

    ds_combined = ds_combined.assign({
        'gcam_cdr_globe': ds_combined['gcam_cdr'].sum(dim=site_dim).assign_attrs(
            units='tCO2/yr', description='Global GCAM CDR (sum across active regions)'
        ),
        f'{cdr_var}_globe': ds_combined[f'{cdr_var}_region'].sum(dim=site_dim).assign_attrs(
            units='t/yr', description=f'Global modeled CDR ({cdr_var}), sum across active regions'
        ),
        'M_D_region': ds_combined['M_D'].sum(dim='deployment').assign_attrs(
            units='t/yr', description='Total rock spread per region per year (sum of M_D across deployments)'
        ),
    })

    return ds_combined

# --- application -> dissolution -> CDR lag, on a cation-budget basis
def compute_region_lag_curves(ds_r: xr.Dataset, RCO2: float = 0.3) -> dict:
    """
    Build cumulative CDR-equivalent curves for one region, on a cation-budget basis.

    Parameters
    ----------
    ds_r : xr.Dataset
        Single-region slice (already sel'd to one site + apprate_base + roughness_factor +
        cec_factor) of the output of gcam_postprocess_all_dims_constApp. Must contain
        'catbudget_adv_region', 'catbudget_sic_region', 'catbudget_gbas_region', and
        'M_GCAM', plus a 'time' coordinate.
    RCO2 : float
        Stoichiometric CDR potential per unit rock mass (tCO2/tRock). Applied to the raw
        rock-mass curve (M_GCAM) to put it on the same tCO2-equivalent basis as the
        catbudget-derived curves, which are already RCO2-scaled.

    Returns
    -------
    dict
        Keys: 't' (time grid), 'C' (cumulative realized CDR), 'C_sic' (cumulative realized
        CDR incl. secondary carbonate formation), 'theta_A' (cumulative CDR-equivalent
        potential from rock applied), 'theta_D' (cumulative CDR-equivalent potential from
        rock dissolved).

    Notes
    -----
    theta_A uses 'M_GCAM' (the region's annual rock-mass target) rather than summing
    'M_D' (rock mass per deployment cohort) over the 'deployment' dim. The area solver in
    gcam_postprocess_region_constApp enforces sum(M_D, dim='deployment') == M_GCAM exactly
    by construction (mu + M_new == M_GCAM every year, including the area-trimming branch),
    confirmed empirically to ~1e-16 relative difference. Reading M_GCAM directly skips
    loading/summing the much larger (time, deployment) M_D array -- ~2x faster per region.
    """
    t_r = ds_r['time'].values

    cdr_rate = -ds_r['catbudget_adv_region'].compute().values
    C = cumulative_trapezoid(cdr_rate, t_r, initial=0.0)

    cdr_sic_rate = -(ds_r['catbudget_adv_region'].compute().values
                      + ds_r['catbudget_sic_region'].compute().values)
    C_sic = cumulative_trapezoid(cdr_sic_rate, t_r, initial=0.0)

    rockapp_rate = ds_r['M_GCAM'].compute().values
    theta_A = cumulative_trapezoid(rockapp_rate, t_r, initial=0.0) * RCO2

    rockdiss_rate = ds_r['catbudget_gbas_region'].compute().values
    theta_D = cumulative_trapezoid(rockdiss_rate, t_r, initial=0.0)

    return {'t': t_r, 'C': C, 'C_sic': C_sic, 'theta_A': theta_A, 'theta_D': theta_D}


def _lag_record_from_curves(curves: dict, percentiles, suffix: str) -> dict:
    """
    Percentile-crossing lag metrics (years) from one set of cumulative curves.

    For each percentile of the curve's own final value, finds the year each of C / C_sic /
    theta_A / theta_D reaches that same absolute level, and differences the crossing years.
    Shared by compute_lag_table (suffix='region') and compute_global_lag (suffix='global').
    """
    t_r, C, C_sic  = curves['t'], curves['C'], curves['C_sic']
    theta_A, theta_D = curves['theta_A'], curves['theta_D']

    record = {}
    for pct in percentiles:
        p_label = f'P{int(pct * 100)}'

        target = pct * C[-1]
        t_C = np.interp(target, C, t_r)
        t_A = np.interp(target, theta_A, t_r)
        t_D = np.interp(target, theta_D, t_r)
        record[f't_app_cdr_adv_{suffix}_{p_label}'] = t_C - t_A
        record[f't_diss_cdr_adv_{suffix}_{p_label}'] = t_C - t_D
        record[f't_app_diss_{suffix}_{p_label}']     = t_D - t_A

        target_sic = pct * C_sic[-1]
        t_C_sic = np.interp(target_sic, C_sic, t_r)
        t_A_sic = np.interp(target_sic, theta_A, t_r)
        t_D_sic = np.interp(target_sic, theta_D, t_r)
        record[f't_app_cdr_adv+sic_{suffix}_{p_label}']  = t_C_sic - t_A_sic
        record[f't_diss_cdr_adv+sic_{suffix}_{p_label}'] = t_C_sic - t_D_sic

    return record


def compute_lag_table(
    dsall: xr.Dataset,
    RCO2: float = 0.3,
    percentiles: tuple = (0.25, 0.50, 0.75),
    site_dim: str = 'site',
    **sel_kwargs,
) -> Tuple[pd.DataFrame, dict]:
    """
    Per-region application -> dissolution -> CDR lag, for one experimental condition.

    Parameters
    ----------
    dsall : xr.Dataset
        Postprocessed SCEPTER/GCAM output (output of gcam_postprocess_all_dims_constApp),
        with a `site_dim` dimension plus whatever other swept-parameter coords/dims it has
        (e.g. apprate_base, dustrad or roughness_factor, cec_factor).
    RCO2 : float
        Passed through to compute_region_lag_curves.
    percentiles : sequence of float
        Percentiles (0-1) of each region's own final cumulative CDR at which to measure lag.
    site_dim : str
        Name of the region/site dimension in dsall (default 'site').
    **sel_kwargs
        Non-site coordinate values identifying the single experimental condition to select,
        e.g. apprate_base=15.0, dustrad=0.333, cec_factor=1.0 -- passed straight to
        `dsall.sel(**sel_kwargs)`. Left generic since the dustrad/roughness_factor coordinate
        name differs between the raw postprocessed output and downstream analysis notebooks
        that convert dustrad to roughness_factor for readability.

    Returns
    -------
    lag_df : pd.DataFrame
        One row per region (indexed by `site_dim`), columns named '{metric}_region_{Pxx}'.
    region_curves : dict[str, dict]
        Per-region cumulative curves (see compute_region_lag_curves), keyed by region name --
        reusable by compute_global_lag without re-selecting from dsall.
    """
    ds_cond = dsall.sel(**sel_kwargs)

    lag_records = []
    region_curves = {}
    for region in ds_cond[site_dim].values:
        ds_r = ds_cond.sel(**{site_dim: region})
        curves = compute_region_lag_curves(ds_r, RCO2=RCO2)
        region_curves[region] = curves

        record = _lag_record_from_curves(curves, percentiles, suffix='region')
        record[site_dim] = region
        lag_records.append(record)

    lag_df = pd.DataFrame(lag_records).set_index(site_dim)
    return lag_df, region_curves


def compute_global_lag(region_curves: dict, lag_df: pd.DataFrame,
                        percentiles: tuple = (0.25, 0.50, 0.75)) -> pd.DataFrame:
    """
    Roll per-region lag curves up to a global scale, two ways.

    Parameters
    ----------
    region_curves : dict[str, dict]
        Per-region cumulative curves, as returned by compute_lag_table.
    lag_df : pd.DataFrame
        Per-region lag table, as returned by compute_lag_table (same regions/condition as
        region_curves).
    percentiles : sequence of float
        Must match the percentiles used to build lag_df.

    Returns
    -------
    pd.DataFrame
        Indexed by metric name (e.g. 't_app_cdr_adv_global_P50'), with columns:
        - 'global_curve_summed': regional C/theta_A/theta_D curves summed first, then the
          percentile-crossing lag applied to the summed curves -- answers "when does the
          whole global system reach X% of its total CDR?"
        - 'cdr_weighted_regional_avg': each region's own lag (from lag_df), averaged with
          weights equal to that region's share of total CDR. Diverges from the curve-summed
          column when regions ramp up out of phase with each other.
    """
    t_g = next(iter(region_curves.values()))['t']  # all regions share the same time grid
    global_curves = {
        't':       t_g,
        'C':       np.sum([rc['C']       for rc in region_curves.values()], axis=0),
        'C_sic':   np.sum([rc['C_sic']   for rc in region_curves.values()], axis=0),
        'theta_A': np.sum([rc['theta_A'] for rc in region_curves.values()], axis=0),
        'theta_D': np.sum([rc['theta_D'] for rc in region_curves.values()], axis=0),
    }
    global_lag = _lag_record_from_curves(global_curves, percentiles, suffix='global')

    cdr_weights = pd.Series({region: rc['C'][-1] for region, rc in region_curves.items()})
    cdr_weights = cdr_weights / cdr_weights.sum()
    global_lag_weighted_avg = {
        col: float((lag_df[col] * cdr_weights).sum())
        for col in lag_df.columns
    }

    rows = []
    for pct in percentiles:
        p_label = f'P{int(pct * 100)}'
        for kind in ['t_app_cdr_adv', 't_diss_cdr_adv', 't_app_diss',
                     't_app_cdr_adv+sic', 't_diss_cdr_adv+sic']:
            g_key = f'{kind}_global_{p_label}'
            r_key = f'{kind}_region_{p_label}'
            rows.append({
                'metric': g_key,
                'global_curve_summed': global_lag[g_key],
                'cdr_weighted_regional_avg': global_lag_weighted_avg[r_key],
            })

    return pd.DataFrame(rows).set_index('metric')


def compute_lag_table_all_dims(
    dsall: xr.Dataset,
    RCO2: float = 0.3,
    percentiles: tuple = (0.25, 0.50, 0.75),
    site_dim: str = 'site',
    verbose: bool = True,
) -> Tuple[xr.Dataset, xr.Dataset]:
    """
    Run compute_lag_table + compute_global_lag for every combination of dsall's swept
    dims, and assemble the results into two xr.Datasets that share dsall's swept-dim
    structure -- so they merge directly onto ds_all / ds_sum (e.g. via `.assign(...)` or
    `xr.merge([ds_sum, lag_ds, global_lag_ds])`), the same way gcam_postprocess_all_dims_constApp
    assembles per-region results into one Dataset.

    Parameters
    ----------
    dsall : xr.Dataset
        Postprocessed SCEPTER/GCAM output (output of gcam_postprocess_all_dims_constApp),
        with a `site_dim` dimension, a 'time' dimension, and one or more swept-parameter
        dims (e.g. apprate_base/dustrad/cec_factor for a standard sweep, or
        collapsed_inputs for a DGSA run).
    RCO2 : float
        Passed through to compute_lag_table.
    percentiles : sequence of float
        Passed through to compute_lag_table / compute_global_lag.
    site_dim : str
        Name of the region/site dimension in dsall (default 'site').
    verbose : bool
        Print a progress line for each swept-dim combination (default True).

    Returns
    -------
    lag_ds : xr.Dataset
        Per-region lag metrics, dims (site_dim, <swept dims>). One data variable per
        lag_df column (e.g. 't_app_cdr_adv_region_P50').
    global_lag_ds : xr.Dataset
        Global lag metrics, dims ('metric', <swept dims>), where 'metric' indexes the
        15 lag quantities (e.g. 't_app_cdr_adv_global_P50'). Variables
        'global_curve_summed' and 'cdr_weighted_regional_avg' are the two ways of
        rolling regions up to a global number (see compute_global_lag).
    """
    import itertools

    loop_dims = [d for d in dsall.dims if d not in (site_dim, 'time', 'deployment')]
    n_total = 1
    for d in loop_dims:
        n_total *= dsall.sizes[d]

    lag_pieces = []
    global_pieces = []
    for k, combo in enumerate(itertools.product(*[dsall[d].values for d in loop_dims])):
        sel_kwargs = dict(zip(loop_dims, combo))

        if verbose:
            combo_str = ', '.join(f'{d}={v}' for d, v in sel_kwargs.items())
            print(f'[{k + 1}/{n_total}]  {combo_str}')

        lag_df, region_curves = compute_lag_table(
            dsall, RCO2=RCO2, percentiles=percentiles, site_dim=site_dim, **sel_kwargs
        )
        global_lag_df = compute_global_lag(region_curves, lag_df, percentiles=percentiles)

        lag_piece    = xr.Dataset.from_dataframe(lag_df)
        global_piece = xr.Dataset.from_dataframe(global_lag_df)
        for dim, val in sel_kwargs.items():
            lag_piece    = lag_piece.assign_coords({dim: val})
            global_piece = global_piece.assign_coords({dim: val})
        lag_piece    = lag_piece.expand_dims(loop_dims)
        global_piece = global_piece.expand_dims(loop_dims)

        lag_pieces.append(lag_piece)
        global_pieces.append(global_piece)

    lag_ds        = xr.combine_by_coords(lag_pieces, combine_attrs='override')
    global_lag_ds = xr.combine_by_coords(global_pieces, combine_attrs='override')

    return lag_ds, global_lag_ds

# --- cation budget where feedstock is the source and we track all sinks
# Columns `cflx_proc.build_cation_df` adds on top of the raw flx_aq-<cation>.txt
# schema. These MUST be stripped before the budget sees the frame: 'charge' is a
# constant (1 or 2) whose cumulative over a multi-decade run trivially clears
# threshold_frac * gbas_final, and the *_source columns are sums of columns that
# are already counted. Leaving any of them in silently inflates catbudget_other
# by roughly the size of the whole signal.
_CATION_PKL_DERIVED_COLS = (
    'noncarbsld_source', 'carbsld_source',
    'co2pot_adv_tonHaYr', 'co2pot_tot_tonHaYr',
    'co2pot_adv_tonHa', 'co2pot_tot_tonHa',
    'charge', 'units', 'flx_type', 'runname', 'cation',
)


def cation_flx_from_pickle(
    outdir: str,
    runname: str,
    cation: str,
    flx_type: str = 'int_flx',
    subdir: str = 'postproc_flxs',
) -> pd.DataFrame:
    """
    Rebuild the raw ``flx_aq-<cation>.txt`` schema from the postprocessing
    pickle ``<subdir>/cationflx_<cation>.pkl``.

    Used as a fallback when the raw txt file didn't survive the upload to aws.
    The pickle is written by `cflx_proc.build_cation_df` from that very txt
    file and applies no numerical transform to the cation fluxes, so it is a
    faithful stand-in. (The `sumCat_adv` docstring's "int_* fluxes are
    multiplied by time" note describes `co2_flx`/`sld_flx`, not this path.)
    Two differences have to be undone here:

      1. the pickle concatenates the transient ('flx') and time-integrated
         ('int_flx') blocks, so one has to be selected
      2. it carries extra derived and metadata columns that the raw file has
         no equivalent for -- see `_CATION_PKL_DERIVED_COLS`

    `build_cation_df` also drops columns whose values are all below 1e-7. That
    cannot change the budget: such a column's cumulative over an 80 year run is
    at most ~8e-6, four-plus orders below a threshold_frac * gbas_final of order
    0.02-2, so it can never flip the "other" threshold. The one unconditionally
    counted column that can go missing this way is the feedstock column in
    control runs, which `read_cation_budget_flx` re-injects as zero anyway.

    Parameters
    ----------
    outdir : str
        location of the output dirs, as passed to `read_cation_budget_flx`
    runname : str
        name of the run output directory
    cation : str
        cation short name, e.g. 'ca'
    flx_type : str
        which block to return: 'int_flx' (the int_flx_aq-*.txt equivalent) or 'flx'
    subdir : str
        subdirectory of the run directory holding the pickle

    Returns
    -------
    pd.DataFrame or None
        a frame matching the raw txt schema, or None if the pickle is absent
        or holds nothing for `flx_type`
    """
    fn_path = os.path.join(outdir, runname, subdir, f'cationflx_{cation}.pkl')
    try:
        df = _read_pickle_retry(fn_path)
    except Exception:
        return None
    if df is None:
        return None

    if 'flx_type' in df.columns:
        df = df[df['flx_type'] == flx_type]
    if df.empty:
        return None

    keep = [col for col in df.columns if col not in _CATION_PKL_DERIVED_COLS]
    return df[keep].reset_index(drop=True)


def read_cation_budget_flx(
    dfin: pd.DataFrame,
    outdir: str,
    dfin_cols_to_keep: list,
    cation_list: list = None,
    flx_subdir: str = 'flx',
    flx_fn_prefix: str = 'int_flx_aq-',
    threshold_frac: float = 0.02,
    sic_minerals: list = None,
    verbose: bool = False,
    pkl_fallback: bool = True,
) -> tuple:
    """
    For each run in dfin, read raw SCEPTER cation flux txt files and compute
    the CDR-potential budget summed over all cations in cation_list.

    Always includes: dustsp (primary source), adv (CDR), tflx (storage), and any
    sic_minerals column present (cc, dlm, arg by default) -- these are named,
    always-tracked pathways grouped into 'catbudget_sic'. Conditionally includes
    any other column (except 'res') in 'catbudget_other' if its absolute
    cumulative value at t_end >= threshold_frac * |dustsp cumulative at t_end|.

    Control runs: the dustsp column (e.g. 'gbas') is absent when no rock is
    applied. If detected, a zero column is injected so controls are processed
    consistently with case runs (catbudget_gbas = 0 for controls by construction).
    Because dustsp is forced to 0, controls have no reference scale for the
    'other' threshold check -- rather than trivially including every column
    (0 always clears a zero threshold), 'other' columns are excluded entirely
    for runs with zero dustsp. SIC minerals are unconditional precisely so this
    asymmetry can't affect them: both case and control always include the same
    SIC columns, which is what lets control subtraction cleanly cancel the large
    shared background (no-rock) SIC/exchange signal instead of leaving it stuck
    on only one side of the subtraction.

    Sign convention — uniform negation of raw SCEPTER file values:
      Positive = source (adds cations to aqueous)
      Negative = sink   (removes cations from aqueous / column)

      catbudget_gbas     = -dustsp_rate  (+ : dissolution adds to aqueous)
      catbudget_adv      = -adv_rate     (− : advection removes from column)
      catbudget_tflx     = -tflx_rate    (+ when releasing (adds to aqueous), − when accumulating (removes from aqueous))
      catbudget_sic      = -sic_rate     (+ when dissolving (adds to aqueous), − when forming (removes from aqueous))
      catbudget_other    = -other_rate   (sign tracks aqueous gain/loss)
      catbudget_residual = sum of all above  (≈ 0 if mass balance is closed)

    Output columns mirror cfp.read_postproc_flux conventions:
      flx_type : 'int_flx' (integrated df) | 'flx' (transient df)
      units    : 't CO2-equiv ha-1' | 't CO2-equiv ha-1 yr-1'
      ctrl     : bool, from dfin['ctrl_run']
      + all columns in dfin_cols_to_keep

    Parameters
    ----------
    dfin : pd.DataFrame
        Batch DataFrame; must have columns 'newrun_id_full', 'dustsp', 'ctrl_run'.
    outdir : str
        Base directory (local or s3://) for run outputs.
    dfin_cols_to_keep : list
        dfin columns to carry over to output DataFrames.
    cation_list : list, optional
        Cations to read. Default: ['ca', 'mg', 'na', 'k'].
    flx_subdir : str
        Subdirectory within each run folder holding flux txt files.
    flx_fn_prefix : str
        Filename prefix: '{flx_fn_prefix}{cat}.txt' → e.g. 'int_flx_aq-ca.txt'.
    pkl_fallback : bool
        If True, when a raw cation txt file can't be read, fall back to
        rebuilding it from postproc_flxs/cationflx_<cat>.pkl (see
        `cation_flx_from_pickle`). Recoveries are always printed, so a run
        sourced from the pickle is visible in the output rather than silent.
    threshold_frac : float
        Column is included only if |cum[-1]| >= threshold_frac * |gbas_cum[-1]|.
    sic_minerals : list, optional
        Column names treated as secondary inorganic carbonate minerals.
    verbose : bool
        If True, print which 'other' columns are included per run.

    Returns
    -------
    df_int : pd.DataFrame
        Time-integrated (cumulative) budget in t CO₂-equiv ha⁻¹.
    df_transient : pd.DataFrame
        Instantaneous rate budget in t CO₂-equiv ha⁻¹ yr⁻¹.
    """
    if cation_list is None:
        cation_list = ['ca', 'mg', 'na', 'k']
    if sic_minerals is None:
        sic_minerals = ['cc', 'dlm', 'arg']

    charge_dict = {'ca': 2, 'mg': 2, 'na': 1, 'k': 1}
    M_CO2 = 44.0
    conv  = 10000.0 / 1e6   # mol m⁻² → t CO₂-equiv ha⁻¹  (× charge × M_CO2)

    int_dfs   = []
    trans_dfs = []

    for run in range(len(dfin)):
        tdf     = dfin.iloc[run]
        dustsp  = tdf['dustsp']
        is_ctrl = bool(tdf['ctrl_run'])
        runname = tdf['newrun_id_full']
        flx_dir = os.path.join(outdir, runname, flx_subdir)

        # --- load cation txt files ------------------------------------------
        # the raw txt is the primary source; a run that lost it in transfer to
        # aws can still be recovered from its postprocessing pickle, which was
        # built from the same file before the upload
        pkl_flx_type = 'int_flx' if flx_fn_prefix.startswith('int_') else 'flx'
        dfs_flx = {}
        for cat in cation_list:
            fn_path = os.path.join(flx_dir, f'{flx_fn_prefix}{cat}.txt')
            try:
                dfs_flx[cat] = pd.read_csv(fn_path, sep=r'\s+', engine='python')
            except Exception as e:
                df_pkl = cation_flx_from_pickle(
                    outdir, runname, cat, flx_type=pkl_flx_type,
                ) if pkl_fallback else None
                if df_pkl is not None:
                    dfs_flx[cat] = df_pkl
                    print(f"  [run {run}] {cat}: raw txt missing, recovered from "
                          f"cationflx_{cat}.pkl")
                else:
                    print(f"  [run {run}] missing {cat}: {e}")

        if not dfs_flx:
            print(f"  [run {run}] no cation files found in {flx_dir}, skipping")
            continue

        # --- control runs: dustsp column absent (no rock applied) -----------
        # inject a zero column so controls are processed identically to cases
        # (catbudget_gbas = 0 for controls by construction)
        if is_ctrl:
            for cat, df in dfs_flx.items():
                if dustsp not in df.columns:
                    dfs_flx[cat] = df.copy()
                    dfs_flx[cat][dustsp] = 0.0

        first_cat = next(iter(dfs_flx))
        t_col = dfs_flx[first_cat].columns[0]
        _t    = dfs_flx[first_cat][t_col].values
        n_t   = len(_t)

        # --- accumulate raw-file rates across all cations (pre-negation) ----
        gbas_rate  = np.zeros(n_t)
        adv_rate   = np.zeros(n_t)
        tflx_rate  = np.zeros(n_t)
        sic_rate   = np.zeros(n_t)
        other_rate = np.zeros(n_t)
        other_cols_used = set()

        for cat in cation_list:
            if cat not in dfs_flx:
                continue
            df  = dfs_flx[cat]
            tc  = df.columns[0]
            ch  = charge_dict.get(cat, 1)
            fac = ch * M_CO2 * conv   # mol m⁻² yr⁻¹ → t CO₂-equiv ha⁻¹ yr⁻¹

            # cumulative (rate_avg × time) at each timestep for threshold check
            cum = {c: df[c].values * _t for c in df.columns if c != tc}
            gbas_final = abs(cum[dustsp][-1]) if dustsp in cum else 1.0
            # control runs apply no rock, so dustsp should carry no real signal --
            # but the raw file's own 'gbas' column isn't always an exact 0.0 (it can
            # sit at floating-point noise, e.g. ~1e-22, rather than being injected
            # fresh), so gbas_final ends up as a tiny nonzero number rather than
            # exactly 0. "abs(x) >= threshold_frac * gbas_final" is then still
            # trivially true for virtually any column, pulling extra species into
            # catbudget_other for controls that a matched case run (with a real,
            # meaningful gbas_final) would exclude. That asymmetry corrupts
            # catbudget_other/residual after control subtraction (case_other −
            # ctrl_other mixes different column sets). Gate on is_ctrl directly
            # instead of gbas_final's magnitude -- with no rock applied there's no
            # reference scale to judge relative significance, so default to
            # excluding "other" columns rather than trivially including all of them.
            has_ref_scale = (not is_ctrl) and (gbas_final > 0)

            for col in df.columns:
                if col == tc or col == 'res':
                    continue
                vals = df[col].values * fac
                if col == dustsp:
                    gbas_rate  += vals
                elif col == 'adv':
                    adv_rate   += vals
                elif col == 'tflx':
                    tflx_rate  += vals
                elif col in sic_minerals:
                    # SIC formation/dissolution is a named, always-tracked pathway
                    # (like adv/tflx above), not a "maybe relevant" one -- included
                    # unconditionally so case and control treat it identically
                    # regardless of gbas_final. This matters because SIC has a large
                    # shared background (no-rock) signal in both case and control;
                    # unconditional inclusion on both sides is what lets control
                    # subtraction cancel that background rather than leaving it
                    # asymmetrically in just one side.
                    sic_rate += vals
                else:
                    if has_ref_scale and col in cum and abs(cum[col][-1]) >= threshold_frac * gbas_final:
                        other_rate += vals
                        other_cols_used.add(col)

        if verbose and other_cols_used:
            print(f"  [run {run} | {runname[-40:]}] other cols: {sorted(other_cols_used)}")

        # --- uniform negation: positive = source, negative = sink -----------
        # Raw SCEPTER aqueous-flux files: positive = out of aqueous/column,
        # negative = into aqueous (dissolution) or within-column phase change.
        # Negating gives the intuitive sign: sources positive, sinks negative.
        cdp_gbas  = -gbas_rate   # + : dissolution source
        cdp_adv   = -adv_rate    # − : advection removes cations from column
        cdp_tflx  = -tflx_rate   # + when exchange releases, − when fills (captures)
        cdp_sic   = -sic_rate    # + when SIC dissolves, − when forms (captures)
        cdp_other = -other_rate  # sign follows aqueous gain/loss

        # residual = sum of all included terms; ≈ 0 if mass balance is closed
        # (any non-zero residual equals the contribution of 'res' + excluded cols)
        cdp_resid = cdp_gbas + cdp_adv + cdp_tflx + cdp_sic + cdp_other

        # --- integrated (rate_avg × time = cumulative) -----------------------
        run_int = pd.DataFrame({
            'time':               _t,
            'flx_type':           'int_flx',
            'units':              't CO2-equiv ha-1',
            'ctrl':               is_ctrl,
            'catbudget_gbas':     cdp_gbas  * _t,
            'catbudget_adv':      cdp_adv   * _t,
            'catbudget_tflx':     cdp_tflx  * _t,
            'catbudget_sic':      cdp_sic   * _t,
            'catbudget_other':    cdp_other * _t,
            'catbudget_residual': cdp_resid * _t,
            'runname':            runname,
        })
        for col in dfin_cols_to_keep:
            run_int[col] = tdf[col]
        int_dfs.append(run_int)

        # --- transient (instantaneous rate = d(cumulative)/dt) ---------------
        run_tr = pd.DataFrame({'time': _t, 'flx_type': 'flx',
                               'units': 't CO2-equiv ha-1 yr-1',
                               'ctrl': is_ctrl, 'runname': runname})
        for col_name, arr in [
            ('catbudget_gbas',     cdp_gbas  * _t),
            ('catbudget_adv',      cdp_adv   * _t),
            ('catbudget_tflx',     cdp_tflx  * _t),
            ('catbudget_sic',      cdp_sic   * _t),
            ('catbudget_other',    cdp_other * _t),
            ('catbudget_residual', cdp_resid * _t),
        ]:
            run_tr[col_name] = np.gradient(arr, _t)
        for col in dfin_cols_to_keep:
            run_tr[col] = tdf[col]
        trans_dfs.append(run_tr)

    if not int_dfs:
        raise ValueError("No output produced — check paths and cation file locations.")

    df_int       = pd.concat(int_dfs,   ignore_index=True)
    df_transient = pd.concat(trans_dfs, ignore_index=True)
    return df_int, df_transient


def catbudget_cdr(
    dfin: pd.DataFrame,
    time_horizon: float,
    dfin_cols_to_keep: list,
    bysite: bool = False,
    ctrl_params: list = None,
) -> tuple:
    """
    Post-process catbudget_flx DataFrame from read_cation_budget_flx:
    subtract control from each case run (if control data is present) and
    clip to time_horizon.

    Mirrors co2_flx_cdr / cation_flx_cdr in cdr_fxns_postproc.py.
    Control subtraction is attempted but skipped gracefully when no matching
    control is found — analogous to rockdiss_synth, since control runs may not
    have raw aqueous-flux txt files (no rock applied → no dissolution budget).

    Parameters
    ----------
    dfin : pd.DataFrame
        Output from read_cation_budget_flx (flx_type='int_flx' or 'flx').
        May contain both case (ctrl=False) and control (ctrl=True) rows.
    time_horizon : float
        Clip to this time (yr) and take summary row at the horizon.
    dfin_cols_to_keep : list
        Coordinate columns to carry forward.
    bysite : bool
        Match control by site when ctrl_params is None.
    ctrl_params : list
        Columns used to uniquely match a case to its control run (overrides bysite).

    Returns
    -------
    dfout_full : pd.DataFrame
        Budget timeseries for each case run (control-subtracted if available).
    dfout_sum : pd.DataFrame
        Budget values at t = time_horizon (one row per case run).
    """
    catbudget_cols = ['catbudget_gbas', 'catbudget_adv', 'catbudget_tflx',
                      'catbudget_sic', 'catbudget_other', 'catbudget_residual']
    cols_to_keep = ['time', 'flx_type', 'runname'] + dfin_cols_to_keep

    df_ctrl = dfin.loc[dfin['ctrl'] == True]
    df_case = dfin.loc[dfin['ctrl'] == False]
    ctrl_available = len(df_ctrl) > 0

    if not ctrl_available:
        print("  NOTE: no control rows in dfin — using raw case values (no subtraction).")
    if not bysite and ctrl_params is None and ctrl_available:
        _ctrl_num = df_ctrl.select_dtypes(include=[np.number])

    case_names = df_case['runname'].unique()
    casedx = 0

    for case in case_names:
        tdf    = df_case[df_case['runname'] == case].copy()
        tdfout = tdf[cols_to_keep].copy()

        # copy budget columns from case (default — no subtraction)
        for col in catbudget_cols:
            if col in tdf.columns:
                tdfout[col] = tdf[col].values

        # attempt control subtraction if control data exists
        if ctrl_available:
            if ctrl_params:
                mask = np.logical_and.reduce(
                    [df_ctrl[c] == tdf[c].iloc[0] for c in ctrl_params]
                )
                tdf_ctrl = df_ctrl.loc[mask].select_dtypes(include=[np.number])
            elif bysite:
                tdf_ctrl = df_ctrl[
                    df_ctrl['site'] == tdf['site'].values[0]
                ].select_dtypes(include=[np.number])
            else:
                tdf_ctrl = _ctrl_num

            if len(tdf_ctrl) == 0:
                print(f"  WARNING: no control match for {case.split('/')[-1][-50:]} — skipping subtraction")
            else:
                # handle time-grid mismatch by interpolation (mirrors co2_flx_cdr)
                if len(tdf) != len(tdf_ctrl):
                    tdf_ctrl = (
                        tdf_ctrl.set_index('time')
                        .reindex(tdf['time'].values)
                        .interpolate(method='linear')
                        .reset_index()
                    )
                for col in catbudget_cols:
                    if col in tdf.columns and col in tdf_ctrl.columns:
                        tdfout[col] = tdf[col].values - tdf_ctrl[col].values

        # clip to time horizon
        tdfout = tdfout.loc[
            tdfout['time'] <= (time_horizon + 1e-6)
        ].reset_index(drop=True)

        # add pipeline-standard metadata columns
        tdfout['time_horizon'] = time_horizon
        tdfout['cdr_fxn']      = 'catbudget_flx'
        tdfout['units']        = dfin['units'].iloc[0]

        # summary row: budget values at the time horizon
        tdf_summary = (
            tdfout[tdfout['time'] == tdfout['time'].max()]
            .drop(columns=['time'])
            .copy()
        )

        if casedx == 0:
            dfout_full = tdfout.copy()
            dfout_sum  = tdf_summary.copy()
        else:
            dfout_full = pd.concat([dfout_full, tdfout],      ignore_index=True)
            dfout_sum  = pd.concat([dfout_sum,  tdf_summary], ignore_index=True)

        casedx += 1

    dfout_full['cdr_fxn'] = 'catbudget_flx'
    dfout_sum['cdr_fxn']  = 'catbudget_flx'
    return dfout_full, dfout_sum

def catbudget_cdr_ds(
    dfin: pd.DataFrame,
    dims: list,
    convert_time_to_timestep: bool = False,
) -> xr.Dataset:
    """
    Convert catbudget_cdr output DataFrame to an xarray Dataset.
    Analogous to co2_flx_cdr_ds, without the loss_percent dimension.

    Parameters
    ----------
    dfin : pd.DataFrame
        Output from catbudget_cdr (timeseries or summary form).
    dims : list
        Index columns for the Dataset (e.g. dfin_cols_to_keep or
        dfin_cols_to_keep_time).  Must be columns in dfin.
    convert_time_to_timestep : bool
        When True, replaces time coordinate with a step index via
        df_to_ds_with_time (for runs with slightly different time grids).

    Returns
    -------
    xr.Dataset
        Dataset with catbudget_* variables and coordinate metadata attributes.
    """
    catbudget_cols = ['catbudget_gbas', 'catbudget_adv', 'catbudget_tflx',
                      'catbudget_sic', 'catbudget_other', 'catbudget_residual']

    # drop metadata columns that shouldn't become data variables
    cols_to_discard = ['units', 'flx_type', 'cdr_fxn', 'time_horizon', 'runname']
    dfx = dfin.drop(columns=[c for c in cols_to_discard if c in dfin.columns])

    if not convert_time_to_timestep:
        dsx = xr.Dataset.from_dataframe(dfx.set_index(dims))
    else:
        # Drop 'time' before df_to_ds_with_time: catbudget reads raw txt files
        # whose time grids can differ slightly from pkl-based datasets (co2_flx,
        # etc.).  Keeping 'time' here causes a MergeError on the shared 'time'
        # variable.  The merged output already gets 'time' from the other datasets.
        dfx_ts = dfx.drop(columns=['time'], errors='ignore')
        dsx = df_to_ds_with_time(dims, dfx_ts)

    # attach attributes (mirrors co2_flx_cdr_ds)
    dsx.attrs['flx_type'] = dfin['flx_type'].iloc[0]
    units_str = dfin['units'].iloc[0]
    for var_name in dsx.data_vars:
        if var_name in catbudget_cols:
            dsx[var_name].attrs['units']   = units_str
            dsx[var_name].attrs['cdr_fxn'] = 'catbudget_flx'

    if 'time_horizon' in dfin.columns:
        dsx['time_horizon'] = dfin['time_horizon'].iloc[0]

    return dsx