# ---
# functions to compute cdr
# 
# --- 
# %%
from scipy.integrate import cumulative_trapezoid
from typing import Tuple
import numpy as np
import pandas as pd
import xarray as xr
import pickle
import math
import os
import re
import glob
import fnmatch

# %% 
# file_path = "/home/tykukla/SCEPTER/scepter_output/hiFert_gbas__site_311a_app_10000p0_psize_200_gbas_field_tau15p0"
# subdir = "postproc_flxs"
# co2_fn = "rockflx_gbas.pkl" # "carbAlk_flxs.pkl"

# tx = pd.read_pickle(os.path.join(file_path, subdir, co2_fn))
# tx

# %% 
# SCEPTER/scepter_output/hiFert_gbas__site_311a_app_10000p0_psize_200_gbas_field_tau15p0/postproc_flxs/co2_flxs.pkl

def read_postproc_flux(
    dfin: pd.DataFrame,  
    outdir: str,
    calc_list: list,
    dfin_cols_to_keep: list,
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
        outdict['co2_flx'] = co2_flx_in(dfin, outdir, dfin_cols_to_keep, rockdiss_feedstock)

    # --- OPTION 2: "camg_flx" -----------------------
    if "camg_flx" in calc_list:
        print("solving camg_flx")
        outdict['camg_flx'] = cation_flx_in(dfin, outdir, dfin_cols_to_keep, rockdiss_feedstock,
                                            cations_to_track = ['ca', 'mg'])
    
    # --- OPTION 3: "totcat_flx" ---------------------
    if "totcat_flx" in calc_list:
        print("solving totcat_flx")
        outdict['totcat_flx'] = cation_flx_in(dfin, outdir, dfin_cols_to_keep, rockdiss_feedstock,
                                              cations_to_track = ['ca', 'mg', 'na', 'k'])
    
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
    # --- assume we're not working in s3 (aws) until we check the outdir
    s3_path = False

    # --- craft the path to the rockdiss file 
    this_fn = f"{fn_pref}_{dustsp}.pkl"
    this_path = os.path.join(outdir, tdf['newrun_id_full'], subdir)
    fn_path = os.path.join(this_path, this_fn)

    # check if it exists
    # ----------------------------
    # check for S3
    if fn_path.startswith("s3://"): # then bring it in from s3
        import s3fs
        fs = s3fs.S3FileSystem()
        s3_path = fs.exists(fn_path)
    # ----------------------------

    # read in the pickle if the path exists
    if (os.path.exists(fn_path) or s3_path):
        rockdf = pd.read_pickle(fn_path) # pd can read in aws if fsspec is imported
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
    # --- assume we're not working in s3 (aws) until we check the outdir
    s3_path = False
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
        # check for S3
        if fn_path.startswith("s3://"): # then bring it in from s3
            import s3fs
            fs = s3fs.S3FileSystem()
            s3_path = fs.exists(fn_path)
        # ----------------------------

        # read in the pickle if the path exists
        if (os.path.exists(fn_path) or s3_path):
            tmpdf = pd.read_pickle(fn_path) # pd can read in aws if fsspec is imported

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
    # --- assume we're not using AWS until we check the outdir later
    s3_path = False

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
        # loop through the different cations to check that they all exist
        existing_files = []
        for cat in cations_to_track:
            fn_thiscat = f"{cation_fn_prefix}{cat}.pkl"
            fn_path = os.path.join(this_path, fn_thiscat)
            # ----------------------------
            # check for S3
            if fn_path.startswith("s3://"): # then bring it in from s3
                import s3fs
                fs = s3fs.S3FileSystem()
                s3_path = fs.exists(fn_path)
            # ----------------------------
            if (os.path.exists(fn_path) or s3_path):
                existing_files.append(fn_path)
        
        # if the length of existing_files is != cations_to_track
        # then we assume some files are missing and we just move on 
        # to the next run (no results are returned for this one)
        if len(existing_files) != len(cations_to_track):
            continue
        
        # otherwise loop through each cation and add them together
        cat_dx = 0   # index of the cation to update in loop
        for fpath in existing_files:
            tmpdf = pd.read_pickle(fpath)  # works with aws if fsspec is imported

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
    feedstock : str
        name of the feedstock to pull from (e.g., 'cc', 'gbas', 'amnt', etc)
        (can usually find this in dfin)
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
    # --- assume we're not using aws until we check the outdir later 
    s3_path = False

    # --- loop through runs
    # track run index
    rundx = 0
    outdf_exist = False
    # loop
    for run in range(len(dfin)):
        tdf = dfin.iloc[run]
        this_path = os.path.join(outdir, tdf['newrun_id_full'], subdir)
        tmpfn = f'{rock_prefix}{feedstock}.pkl'
        fn_path = os.path.join(this_path, tmpfn)
        # ----------------------------
        # check for S3
        if fn_path.startswith("s3://"): # then bring it in from s3
            import s3fs
            fs = s3fs.S3FileSystem()
            s3_path = fs.exists(fn_path)
        # ----------------------------
        # read in the pickle if the path exists
        if (os.path.exists(fn_path) or s3_path):
            tmpdf = pd.read_pickle(fn_path)

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

    Returns 
    -------
    dict 
        [1] dictionary with the same structure as the input, but with 
        CDR calcs over time
        [2] dictionary with the same structure as the input, but with 
        CDR calcs integrated to the time horizon
    """
    # --- create an empty dictionary to hold the results
    outdict_full = outdict_sum = {}

    #
    # iterate through options in the calc_list
    #

    # --- OPTION 1: "co2_flx" ------------------------
    if "co2_flx" in calc_list:
        tc = 'co2_flx'
        print(f"solving {tc}")
        dfin = flx_dict[tc]
        outdict_full[tc], outdict_sum[tc] = co2_flx_cdr(dfin, time_horizon, dfin_cols_to_keep)

    # --- OPTION 2: "camg_flx" -----------------------
    if "camg_flx" in calc_list:
        tc = 'camg_flx'
        print(f"solving {tc}")
        dfin = flx_dict[tc]
        # (note, we use the same cation_flx_cdr fxn for both camg and totcat
        #  because it just calculates based on whatever cats are provided)
        outdict_full[tc], outdict_sum[tc] = cation_flx_cdr(dfin, time_horizon, dfin_cols_to_keep, cation_tag=tc)
    
    # --- OPTION 3: "totcat_flx" ---------------------
    if "totcat_flx" in calc_list:
        tc = 'totcat_flx'
        print(f"solving {tc}")
        dfin = flx_dict[tc]
        # (note, we use the same cation_flx_cdr fxn for both camg and totcat
        #  because it just calculates based on whatever cats are provided)
        outdict_full[tc], outdict_sum[tc] = cation_flx_cdr(dfin, time_horizon, dfin_cols_to_keep, cation_tag=tc)
    
    # --- OPTION 4: "carbalk_flx" --------------------
    if "carbalk_flx" in calc_list:
        print("solving carbalk_flx")
        outdict['carbalk_flx'] = carbalk_flx_cdr()
    
    # --- OPTION 5: "rockdiss" -----------------------
    if "rockdiss" in calc_list:
        tc = 'rockdiss'
        print(f"solving {tc}")
        dfin = flx_dict[tc]
        outdict_full[tc], outdict_sum[tc] = rockdiss_synth(dfin, time_horizon, dfin_cols_to_keep)

    # 
    # return the result
    # 
    return outdict_full, outdict_sum


def co2_flx_cdr(
    dfin: pd.DataFrame, 
    time_horizon: float,
    dfin_cols_to_keep: list,
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

        # if case and control are different lengths, we need to interpolate
        # (this happens sometimes due to shifts in how the timesteps are handled in a given run)
        if len(tdf) != len(tdf_ctrl):
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

    # get control run
    df_ctrl = dfin.loc[dfin['ctrl'] == True]
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
    mineral: str
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
    
    # -- crushing -- # 
    Efactor_org : str
        ["MRO" | "RFC" | "SERC"] which crushing Efactor to use
    mineral : str
        ["cc" | "gbas" | "wls" None] which mineral is being transported / crushed (this sets the bond work index)
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
        "cc": 12.10,    # mean of limestone estimates in Kanda and Kotake (2007) and Bond, 1961
        "wls": 8.33     # see Marco and Caterina, 2022 (https://www.proquest.com/docview/2781737786?pq-origsite=gscholar&fromopenview=true&sourcetype=Conference%20Papers%20&%20Proceedings)
    }
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
            print(f"Emissions Calculator: Cannot convert column '{col}' to float: {e}")
    
    # return the result
    return df



def cdr_ds(
    cdr_dict: dict,
    dims: list,
    cdr_calc_list: list,
    loss_percents: np.array = np.linspace(100,1,50),
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
        dsout_dict[tc] = co2_flx_cdr_ds(dfin, dims, loss_percents)

    # --- OPTION 2: "camg_flx" -----------------------
    if "camg_flx" in cdr_calc_list:
        tc = 'camg_flx'
        print(f"solving {tc}")
        dfin = cdr_dict[tc]
        # (note, we use the same cation_flx_cdr fxn for both camg and totcat
        #  because it just calculates based on whatever cats are provided)
        dsout_dict[tc] = cation_flx_cdr_ds(dfin, dims, loss_percents, cation_flag='camg', cdr_calc_cols=["co2pot_tot", "co2pot_adv", "DICpot_adv"])
    
    # --- OPTION 3: "totcat_flx" ---------------------
    if "totcat_flx" in cdr_calc_list:
        tc = 'totcat_flx'
        print(f"solving {tc}")
        dfin = cdr_dict[tc]
        # (note, we use the same cation_flx_cdr fxn for both camg and totcat
        #  because it just calculates based on whatever cats are provided)
        dsout_dict[tc] = cation_flx_cdr_ds(dfin, dims, loss_percents, cation_flag='totcat', cdr_calc_cols=["co2pot_tot", "co2pot_adv", "DICpot_adv"])
    
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
        dsout_dict[tc] = rockdiss_ds(dfin, dims)

    # merge the dictionary of datasets into a single ds
    dsout = xr.merge(dsout_dict.values())
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

    Returns
    -------
    xr.Dataset
        removals flux dataset
    """
    # get a list of just the cdr vars
    cdr_cols = [col for col in cdr_calc_cols if col.startswith('cdr')]

    # add loss_percents to the dims
    dims_full = dims + ['loss_percent']


    # determine the columns we want to drop
    cols_to_discard = ['units', 'flx_type', 'cdr_fxn', 'time_horizon']
    dfx = dfin.drop(columns=cols_to_discard)

    # create a cdr df to loop over loss percents
    cdr_cols_andDims = dims + cdr_calc_cols
    dfx_cdr_short = dfx[cdr_cols_andDims]

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

    # create the cdr dataset
    dfx_full_idx = dfx_cdr_full.set_index(dims_full)
    dsx_full = xr.Dataset.from_dataframe(dfx_full_idx)

    # create the flux dataset (not defined over loss_percent dimension)
    # using the columns that aren't already in the full ds
    dfx_idx = dfx.drop(columns=list(dsx_full.data_vars)).set_index(dims)
    dsx_x = xr.Dataset.from_dataframe(dfx_idx)

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
    dims_full = dims + ['loss_percent']

    # determine the columns we want to drop
    cols_to_discard = ['units', 'flx_type', 'cdr_fxn', 'time_horizon']
    dfx = dfin.drop(columns=cols_to_discard)

    # create a cdr df to loop over loss percents
    cdr_cols_andDims = dims + cdr_calc_cols
    dfx_cdr_short = dfx[cdr_cols_andDims]

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

    # create the cdr datasets
    dfx_full_idx = dfx_cdr_full.set_index(dims_full)
    dsx_full = xr.Dataset.from_dataframe(dfx_full_idx)

    # create the flux datasets (not defined over loss_percent dimension)
    # using the columns that aren't already in the full ds
    dfx_idx = dfx.drop(columns=list(dsx_full.data_vars)).set_index(dims)
    dsx_x = xr.Dataset.from_dataframe(dfx_idx)

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
    
    Returns
    -------
    xr.Dataset
        removals flux dataset
    """
    # determine the columns we want to drop
    cols_to_discard = ['truck_km', 'barge_km', 'barge_diesel_km', 'p80_input',
                        'Efactor_org', 'bondwork_index', 'time_horizon']
    dfx = dfin.drop(columns=cols_to_discard)

    # rename the adv column to be clearer that it's rock
    dfx = dfx.rename(columns = {'adv': 'adv_feedstock'})

    # create the flux datasets (not defined over loss_percent dimension)
    dfx_idx = dfx.set_index(dims)
    dsx = xr.Dataset.from_dataframe(dfx_idx)

    # --- add attributes
    dsx['time_horizon'] = dfin['time_horizon'][0]
    dsx['truck_km'] = dfin['truck_km'][0]
    dsx['barge_km'] = dfin['barge_km'][0]
    dsx['barge_diesel_km'] = dfin['barge_diesel_km'][0]
    dsx['p80_input'] = dfin['p80_input'][0]
    dsx['Efactor_org'] = dfin['Efactor_org'][0]
    dsx['bondwork_index'] = dfin['bondwork_index'][0]

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
    tdf : pd.Series
        single row of the batch df containing the info for this specific run
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
            with fs.open(fn_path, mode='rb') as fnx:
                tmpds = xr.open_dataset(fnx)
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
    

def prof_batchprocess_singlevar(
    dfin: pd.DataFrame,
    outdir: str,
    batch_axes: list,
    filename_base: str,
    dustsp: str=None,
    subdir: str = "postproc_profs",
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
    
    Returns
    -------
    xr.Dataset
        combination of all xarray datasets in the batch
    '''

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
) -> dict:
    '''
    Wrapper around prof_batchprocess_singlevar to repeat that function for all
    variables in the proc_dict. Resulting datasets are stored in a 
    dictionary to be saved later

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
            ds = prof_batchprocess_singlevar(dfin, outdir, batch_axes, key, dustsp)
            # make sure it's not empty
            if not ds.variables: # continue to next iter if it is empty
                print(f"Warning: batch profile processing {key} returned no results. Check for typos?")
                continue
            outdict[key] = ds
    
    # return result
    return outdict



def save_batch_postproc_profOnly(
    dsdict: dict,
    dustsp_tag: str,
    filename_suffix: str="batch",
    save_directory: str=None,
    base_path: str=None,
    base_dir_name: str=None,
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

