# ---
# functions to compute cdr
# 
# --- 
# %%
from scipy.integrate import cumulative_trapezoid
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

# --- FUNCTION to preprocess .txt files for consistent delimiters
def preprocess_txt(file_path):
    data = []  # Initialize a list to store the processed data

    # Initialize a flag to determine if we are reading the header
    is_header = True

    # Read the file line by line and process the data
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()  # Remove leading/trailing whitespace
            if is_header:
                # Split the first line into column names
                column_names = re.split(r'\s+', line)
                is_header = False
            else:
                # Split the other lines into data values
                values = re.split(r'\s+', line)
                data.append(values)

    # Create a DataFrame with the processed data and set column names
    df = pd.DataFrame(data, columns=column_names)
    # return
    return df


# --- FUNCTION to read in the flux data for CDR calculation
def read_flux(dfin,      # dataframe with batch inputs
              var_fn,    # filename to read in 
              outdir,    # base dir for model output
              fn_varInclude = [],   # will read all species if empty, otherwise just the ones defined
             ):
    # --- read in flux data 
    # define file name pattern
    df = pd.DataFrame()  # initialize empty df to store dat
    fn_pref = [("int_"+var_fn), var_fn] # ["int_flx_co2sp", "flx_co2sp"]   # will read in all species unless fn_varInclude is defined
    varCheck = True if len(fn_varInclude) > 0 else False
    fn_ext = ".txt"
    
    
    # loop through runs
    for run in range(len(dfin)):
        tdf = dfin.iloc[run]
        this_path = os.path.join(outdir, tdf['newrun_id_full'])
        flx_path, prof_path = os.path.join(this_path, "flx"), os.path.join(this_path, "prof")
    
        print("now reading in run " + str(run + 1) + " of " + str(len(dfin)) + "...")
        for fset in fn_pref:
            # set pattern
            fn_pattern = f"{fset}-*{fn_ext}"
            # get list of filenames
            file_paths = glob.glob(f"{flx_path}/{fn_pattern}")
        
            # read in data and concatenate
            for file_path in file_paths:
                # get the variable 
                varpattern = re.escape(fset) + r'-(.*?).txt'
                varmatch = re.search(varpattern, file_path)
                var = varmatch.group(1)
                # skip this step if it's not in the include arr
                if varCheck:
                    if var not in fn_varInclude:
                        continue
                # read in
                dfi = preprocess_txt(file_path)
                # apply pd.to_numeric to all columns using the "map" method
                dfi = dfi.map(pd.to_numeric)
                # add set, var, spinrun, ctrl
                dfi["set"] = fset
                dfi["var"] = var
                dfi["spinrun"] = tdf["spinrun"]
                dfi["runname"] = tdf["newrun_id_full"]
                if 'site' in dfin.columns:
                    dfi['site'] = tdf['site']
                else:
                    try:  # in case "site" isn't in newrun_id we try it first and otherwise pull from climatefiles
                        dfi["site"] = tdf['climatefiles']
                    except:
                        dfi["site"] = re.search(r'site_(\d{3})', tdf['newrun_id_full'])[0]
                dfi["ctrl"] = tdf["ctrl_run"]
                # add dustrate vars
                if "dustrate_ton_ha_yr" in tdf.index:
                    dfi["dustrate_ton_ha_yr_constant"] = tdf["dustrate_ton_ha_yr"]
                if "dustrad" in tdf.index:
                    dfi["dustrad"] = tdf["dustrad"]
                if "dustsp" in tdf.index:
                    dfi["dustsp"] = tdf["dustsp"]
                # combine
                df = pd.concat([df, dfi], ignore_index=True)
    
    # drop all time slices dangerously close to zero (these produce astronomical (like 10^10 or higher) residuals)
    df = df.loc[df['time'] > 1e-3]
    
    # sort by time and depth
    df = df.sort_values(by=["runname", "site", "var", "time"])
    # --- return result
    return df


# --- FUNCTION to add dust flux over time from dust.txt output file 
def add_dust_from_file(outdir: str, 
                       df: pd.DataFrame,
                       dust_fn: str="dust.txt")->pd.DataFrame:
    """
    take a dataframe produced by the read_flux function and append the dust flux data
    to it using the dust.txt files in the output dirs. We don't append the transient 
    dust applications because they can happen in short time steps that don't show up 
    in the output flux data, so interpolating to the flux timesteps can erase entire
    applications. Instead, we also output the dust data itself as its own df with the 
    correct timesteps

    Parameters
    ----------
    outdir : str
        path to output dir (ex: /home/name/SCEPTER/scepter_output)
    df : pd.DataFrame
        pandas dataframe output by read_flux function
    dust_fn : str
        the name of the dust output file we want ("dust.txt" by default)
    
    Returns 
    -------
    pd.DataFrame 
        the input dataframe with the dust flux data appended
    pd.DataFrame
        the compiled dust.txt data for all runs
    """
    # create empty new columns
    # df['dust1_ton_ha_yr'], df['dust2_ton_ha_yr'] = np.nan, np.nan
    df['int_dust1_ton_ha_yr'], df['int_dust2_ton_ha_yr'] = np.nan, np.nan
    df['dust1_avg_ton_ha_yr'], df['dust2_avg_ton_ha_yr'] = np.nan, np.nan

    # group by vars to get individual runs
    grouped = df.groupby(['site', 'runname', 'set', 'var'])

    # empty list of compiled dust dfs
    dfdust_all = []

    # ... loop through runs
    for (site, runname, set_, var), group in grouped:    
        # read in dust data
        dust_path = os.path.join(outdir, runname, "flx", "dust.txt")
        dfdust = preprocess_txt(dust_path)
        # apply pd.to_numeric to all columns using the "map" method
        dfdust = dfdust.map(pd.to_numeric)
        # re-calculate integral because right now it's integrated by timestep, not
        # by the entire run itself
        dfdust['int_dust1_g_m2_yr'] = cumulative_trapezoid(dfdust['dust1_g_m2_yr'], dfdust['time'], initial=0)
        dfdust['int_dust2_g_m2_yr'] = cumulative_trapezoid(dfdust['dust2_g_m2_yr'], dfdust['time'], initial=0)
        # get average dust flux over the whole run
        tmp_maxTime = dfdust['time'].max()
        tmp_maxTime_dx = dfdust['time'].idxmax()   # max timestep (for getting integrated dust)
        tmp_dust1_int = dfdust.loc[tmp_maxTime_dx, 'int_dust1_g_m2_yr']
        tmp_dust2_int = dfdust.loc[tmp_maxTime_dx, 'int_dust2_g_m2_yr']
        tmp_dust1_rate = tmp_dust1_int / tmp_maxTime  # mean dust rate over entire run
        tmp_dust2_rate = tmp_dust2_int / tmp_maxTime  # mean dust rate over entire run

        # add metadata
        dfdust['dust1_avg_ton_ha_yr'] = tmp_dust1_rate / 100  # divide by 100 to convert g/m2/yr to ton/ha/yr
        dfdust['dust2_avg_ton_ha_yr'] = tmp_dust2_rate / 100  # divide by 100 to convert g/m2/yr to ton/ha/yr
        dfdust['site'] = site
        dfdust['runname'] = runname
        dfdust['set'] = set_
        dfdust['var'] = var

        # add to list
        dfdust_all.append(dfdust)

        # add dust data 
        if len(dfdust) != group.shape[0]: # group.shape[0] is number of rows in group
            # drop duplicates in the 'time' column, keeping first occurrence
            dfdust_nodup = dfdust.drop_duplicates(subset='time', keep='first')
            # then interpolate the data 
            # merge the DataFrames on 'time' using outer join to keep all time points
            df_merged = pd.merge(group, dfdust_nodup, on='time', how='outer', suffixes=('', '_orig'))
            # interpolate the columns of interest 
            # df_merged['dust1_g_m2_yr'] = df_merged['dust1_g_m2_yr'].interpolate()
            # df_merged['dust2_g_m2_yr'] = df_merged['dust2_g_m2_yr'].interpolate()
            df_merged['int_dust1_g_m2_yr'] = df_merged['int_dust1_g_m2_yr'].interpolate()
            df_merged['int_dust2_g_m2_yr'] = df_merged['int_dust2_g_m2_yr'].interpolate()
            # keep only the points in the group time steps
            # dust1 = df_merged[df_merged['time'].isin(group['time'])]['dust1_g_m2_yr'].values / 100
            # dust2 = df_merged[df_merged['time'].isin(group['time'])]['dust2_g_m2_yr'].values / 100
            intdust1 = df_merged[df_merged['time'].isin(group['time'])]['int_dust1_g_m2_yr'].values / 100 # divide by 100 to convert g/m2/yr to ton/ha/yr
            intdust2 = df_merged[df_merged['time'].isin(group['time'])]['int_dust2_g_m2_yr'].values / 100 # divide by 100 to convert g/m2/yr to ton/ha/yr
        else:
            # dust1 = dfdust['dust1_g_m2_yr'].values / 100 # divide by 100 to convert g/m2/yr to ton/ha/yr
            # dust2 = dfdust['dust2_g_m2_yr'].values / 100
            intdust1 = dfdust['int_dust1_g_m2_yr'].values / 100 # divide by 100 to convert g/m2/yr to ton/ha/yr
            intdust2 = dfdust['int_dust2_g_m2_yr'].values / 100 # divide by 100 to convert g/m2/yr to ton/ha/yr

        # add back to df
        # df.loc[group.index, 'dust1_ton_ha_yr'] = dust1
        # df.loc[group.index, 'dust2_ton_ha_yr'] = dust2
        df.loc[group.index, 'int_dust1_ton_ha_yr'] = intdust1
        df.loc[group.index, 'int_dust2_ton_ha_yr'] = intdust2
        df.loc[group.index, 'dust1_avg_ton_ha_yr'] = tmp_dust1_rate / 100  # divide by 100 to convert g/m2/yr to ton/ha/yr
        df.loc[group.index, 'dust2_avg_ton_ha_yr'] = tmp_dust2_rate / 100  # divide by 100 to convert g/m2/yr to ton/ha/yr

    # compile output dust df
    dfdust_out = pd.concat(dfdust_all, ignore_index=True)
    # return result
    return df, dfdust_out




# --- FUNCTION to compute diffusive and advective CDR 
def cdr_dif_adv(df,   # df produced by the read_flux function
               cdr_var,  # variable to use for CDR calculation (usually pco2) - # if using *flx_co2sp: [DIC, co2g] ; if using *flx_gas: [pco2]
                feedstock,   # name of feedstock (only matters if cc or not cc for now)
               ):
    # --- calculate CDR 
    allsites = df['site'].unique() 
    
    # initialize columns to fill
    # define the variable names for saving 
    dif_c = "cdr_dif_component"
    resp_c = "cdr_resp_component"
    adv_c = "cdr_adv_component"
    dif = "cdr_dif"
    adv = "cdr_adv"
    
    df[dif_c], df[resp_c], df[adv_c], df[dif], df[adv] = -9999., -9999., -9999., -9999., -9999.
    
    # loop through each site
    for site in allsites:
        dfsite = df.loc[df['site'] == site]
        # get control run
        dfsite_ctrl = dfsite.loc[dfsite['ctrl'] == True]
        dfsite_notctrl = dfsite.loc[dfsite['ctrl'] == False]
    
        # loop through non-control runs
        nonctrl_runs = dfsite_notctrl['runname'].unique()
        for trun in nonctrl_runs:
            print(trun)
            # split into case and control
            tdf_case = dfsite_notctrl.loc[(dfsite_notctrl['var'] == cdr_var) & (dfsite_notctrl['runname'] == trun)]
            tdf_ctrl = dfsite_ctrl.loc[(dfsite_ctrl['var'] == cdr_var)]
            # compute cdr 
            # loop through sets
            sets = tdf_case['set'].unique()
            for thisset in sets:
                ttdf_case = tdf_case.loc[tdf_case['set'] == thisset]
                ttdf_ctrl = tdf_ctrl.loc[tdf_ctrl['set'] == thisset]
                
                # if case and control are different lengths, we need to interpolate
                # (this happens sometimes due to shifts in how the timesteps are handled in a given run)
                if len(ttdf_case) != len(tdf_ctrl):
                    ctrl_dif = ttdf_ctrl.set_index('time')['dif'].reindex(ttdf_case['time']).interpolate(method='linear').values
                    ctrl_resp = ttdf_ctrl.set_index('time')['g2'].reindex(ttdf_case['time']).interpolate(method='linear').values
                    ctrl_adv = ttdf_ctrl.set_index('time')['adv'].reindex(ttdf_case['time']).interpolate(method='linear').values
                else:
                    ctrl_dif = ttdf_ctrl['dif'].values
                    ctrl_resp = ttdf_ctrl['g2'].values
                    ctrl_adv = ttdf_ctrl['adv'].values
    
                # --------
                # [TROUBLESHOOT] -- what happens if we set the ctrl equal to some mean initial state??
                # ctrl_dif = ttdf_ctrl['dif'][3:10].values.mean()
                # ctrl_resp = ttdf_ctrl['g2'][3:10].values.mean()
                # ctrl_adv = ttdf_ctrl['adv'][3:10].values.mean()
                # [TROUBLESHOOT] -- what happens ctrl fluxes are zero
                # ctrl_dif = 0
                # ctrl_resp = 0
                # ctrl_adv = 0
                # --------
                
                # get each component
                dif_component = -1*(ttdf_case['dif'].values - ctrl_dif)
                resp_component = -1*(ttdf_case['g2'].values - ctrl_resp)
                adv_component = ttdf_case['adv'].values - ctrl_adv
                
                # compute dif and adv versions
                cdr_dif = dif_component + np.minimum(resp_component, 0)
                cdr_adv = adv_component + np.minimum(resp_component, 0)

                # apply stoichiometric factor if cc feedstock (TK update this to something more accurate...)
                if feedstock == "cc":
                    cdr_adv = cdr_adv / 2   # only half the advected carbon comes from the atmosphere
    
                # add back to df
                cond_case = (df["site"] == site) & (df['ctrl'] == False) & (df['runname'] == trun) & (df['set'] == thisset) & (df['var'] == cdr_var)
                df.loc[cond_case, dif_c] = dif_component
                df.loc[cond_case, resp_c] = resp_component
                df.loc[cond_case, adv_c] = adv_component
                df.loc[cond_case, dif] = cdr_dif
                df.loc[cond_case, adv] = cdr_adv
    
    # return result     
    return df


# --- FUNCTION to compute cc-silicate cdr difference at different loss rates 
def cdr_compare_lossrange(df_cc,    # cdr dataframe for cc feedstock runs
                          df_sil,   # cdr dataframe for sil feedstock runs
                          cc_app_fixed,   # counterfactual cc application rate (must be in df_cc)
                         time_horizon,    # time duration over which to compare results
                         cdvar,     # cdr variable to compare
                         thissite,  # site to compare
                          var_fn,   # variable filename
                          cdr_var,  # species used for tracking cdr (e.g., DIC; pco2)
                         ):
    # --- troubleshoot; check that application rate is in df_cc
    if not cc_app_fixed in df_cc['dustrate_ton_ha_yr'].unique():
        available_dustrates = df_cc['dustrate_ton_ha_yr'].unique()
        raise ValueError(f"The CC application rate is not in the cc dataframe. Please select a value from the following: {available_dustrates}")
    # --- get CDR for set timestep and control case
    # ... get counterfactual cc dat
    df_cc_cf = df_cc.loc[(df_cc["ctrl"] == False) & (df_cc["set"] == ("int_"+var_fn)) & (df_cc["var"] == cdr_var) & 
                        (df_cc["site"] == thissite) & (df_cc['dustrate_ton_ha_yr'] == cc_app_fixed) &
                        (df_cc['time'] <= time_horizon)]
    # multiply by time to get integrated cdr at the last timestep
    ts_max_idx = df_cc_cf['time'].idxmax()
    cdr_cc_cf = (df_cc_cf[cdvar]*df_cc_cf['time'])[ts_max_idx]
    
    # ... loop through silicate application rates
    appRates = np.sort(df_sil['dustrate_ton_ha_yr'].unique())
    appRates = appRates[appRates > 0] # remove control
    loss_percents = np.linspace(100, 1, 50)
    outarr = np.zeros((len(loss_percents), len(appRates)))
    # indices for tracking
    lpdx, appdx = 0,0
    
    for loss_percent in loss_percents:
        for apprate in appRates:
            df_sil_case = df_sil.loc[(df_sil["ctrl"] == False) & (df_sil["set"] == ("int_"+var_fn)) & (df_sil["var"] == cdr_var) & 
                            (df_sil["site"] == thissite) & (df_sil['dustrate_ton_ha_yr'] == apprate) &
                            (df_sil['time'] <= time_horizon)]
            # get integrated CDR at last timestep (multiply by time to get integrated value)
            ts_max_idx = df_sil_case['time'].idxmax()
            cdr_sil = (df_sil_case[cdvar]*df_sil_case['time'])[ts_max_idx]
    
            # apply loss
            sil_cdr = cdr_sil - ((loss_percent/100)*cdr_sil)
            cc_cdr = cdr_cc_cf - (2*(loss_percent/100)*cdr_cc_cf)
            # disallow avoided emissions accounting
            if cc_cdr < 0:
                cc_cdr = 0
    
            # add result to outarr
            outarr[lpdx, appdx] = sil_cdr - cc_cdr
            
            # update index
            appdx += 1
            if appdx >= len(appRates):
                appdx = 0
        # update index
        lpdx += 1
    
    # get into data array
    D_apprate = appRates - cc_app_fixed
    mycoords = {'loss_percent': loss_percents, 'D_appRate': D_apprate}
    mydims=['loss_percent', 'D_appRate']
    da = xr.DataArray(outarr, coords=mycoords, dims=mydims)

    # return result
    return da


# --- FUNCTION get cc and silicate cdr for different app rates at a given time horizon
def cdr_per_apprate(df_cc,    # cdr dataframe for cc feedstock runs
                    df_sil,   # cdr dataframe for sil feedstock runs
                    time_horizon,    # time duration over which to compare results
                    cdvar,     # cdr variable to compare
                    thissite,  # site to compare
                    var_fn,   # variable filename
                    cdr_var,  # species used for tracking cdr (e.g., DIC; pco2)
                   ):
    # get cdr vs apprate df
    df_cc_ar = df_cc.loc[(df_cc["ctrl"] == False) & (df_cc["set"] == ("int_"+var_fn)) & (df_cc["var"] == cdr_var) & 
                        (df_cc["site"] == thissite) & (df_cc['time'] <= time_horizon)]
    df_sil_ar = df_sil.loc[(df_sil["ctrl"] == False) & (df_sil["set"] == ("int_"+var_fn)) & (df_sil["var"] == cdr_var) & 
                        (df_sil["site"] == thissite) & (df_sil['time'] <= time_horizon)]
    
    grouped_cc = df_cc_ar.groupby('dustrate_ton_ha_yr')
    grouped_sil = df_sil_ar.groupby('dustrate_ton_ha_yr')
    
    # loop through each application rate case
    cdr_cc_ar, cdr_sil_ar, cc_ar, sil_ar = [], [], [], []
    # [CALCITE]
    for name, group in grouped_cc:
        # get integrated CDR at last timestep (multiply by time to get integrated value)
        ts_max_idx = group['time'].idxmax()
        cdr_tmp = (group[cdvar]*group['time'])[ts_max_idx]
        ar_tmp = group['dustrate_ton_ha_yr'][ts_max_idx]
        cdr_cc_ar.append(cdr_tmp)
        cc_ar.append(ar_tmp)
    
    # [SILICATE]
    for name, group in grouped_sil:
        # get integrated CDR at last timestep (multiply by time to get integrated value)
        ts_max_idx = group['time'].idxmax()
        cdr_tmp = (group[cdvar]*group['time'])[ts_max_idx]
        ar_tmp = group['dustrate_ton_ha_yr'][ts_max_idx]
        cdr_sil_ar.append(cdr_tmp)
        sil_ar.append(ar_tmp)
    
    # create pandas dfs
    dfcc_ar = pd.DataFrame({
        'apprate': cc_ar,
        'cdr': cdr_cc_ar,
        'fs': 'cc'
    })
    
    dfsil_ar = pd.DataFrame({
        'apprate': sil_ar,
        'cdr': cdr_sil_ar,
        'fs': 'sil'
    })

    # return results
    return dfcc_ar, dfsil_ar



# --- FUNCTION same as above but only takes in one df at a time and solves for 
#     an arbitrary number of user-defined groups
def cdr_per_group(df_cdr: pd.DataFrame,    # cdr dataframe for cc feedstock runs
                  dust_dict: dict,   # dictionary describing how to treat the integrated mean dust flux calculation
                    time_horizon: float,    # time duration over which to compare results
                    cdvar: str,     # cdr variable to compare
                    thissite: str,  # site to compare
                    var_fn: str,   # variable filename
                    cdr_var: str,  # species used for tracking cdr (e.g., DIC; pco2)
                    group_vars: list,   # list of variable names to group the CDR output by
                    dfdust: pd.DataFrame=None, # dust df, none by default because you can run this without it 
                   ) -> pd.DataFrame:
    """
    Given a dataframe produced by cdr_dif_adv function, compute time integrated CDR
    for a prescribed time horizon and do this for every group in the dataframe. 
    Groups are defined by column names listed in group_vars. 

    Parameters
    ----------
    df_cdr : pd.DataFrame
        dataframe with CDR values from the cdr_dif_adv function
    dust_dict : dict
        dictionary with inputs regarding how to treat the dust calculation. keys include:
        'dust_from_file' : [True | False] whether to use the dfdust file to constrain mean dust flux
        'dustname_dustfile' : [str] name of the dust var in the dust file (only matters if dust_from_file==True)
        'dustname_cdrfile_idx' : [int] index of group_vars for the dust name (only matters if dust_from_file==True)
        'runname_idx' : [int] idex of group_vars for the run name (only matters if dust_from_file==True)
        'convert_to_tons' : [bool] whether the dfdust dust flux needs to be converted to tons (only matters if dust_from_file==True)
        'round_dust_digits' : [int] how many digits to round dust flux to (after converting, if convert_to_tons==True) (only matters if dust_from_file==True)
    time_horizon : float
        [years] time horizon over which to integrate total CDR 
    cdvar : str
        cdr variable to use (e.g., "cdr_dif")
    thissite : str
        site name for filtering df_cdr (e.g., "site_311a")
    var_fn : str
        filename that has the variable of interest (e.g., "flx_gas" or "flx_co2sp")
    cdr_var : str
        variable within the var_fn filename to use. if using *flx_co2sp: [DIC, co2g] ; 
        if using *flx_gas: [pco2] (Note we filter by the filename and variable because the
        same variable name can exist in multiple files where it means different things)
    group_vars : list
        list of column names in df_cdr with which to group the CDR calculations and output
        (e.g., ["dustrate_ton_ha_yr", "dustrad"] for grouping by dustrate and radius)
    dfdust : pd.DataFrame
        result of the add_dust_from_file function. only used if dust_dict['dust_from_file']==True
    """
    # get cdr vs apprate df
    df_cdr_ar = df_cdr.loc[(df_cdr["ctrl"] == False) & (df_cdr["set"] == ("int_"+var_fn)) & (df_cdr["var"] == cdr_var) & 
                        (df_cdr["site"] == thissite) & (df_cdr['time'] <= time_horizon)]
    
    # confirm that runname is one of the groups
    if "runname" not in group_vars:
        print("Warning: 'runname' needs to be one of the group_vars to align with dfdust!")

    # # check if we need to round 
    # for dx in range(len(group_vars)):
    #     if round_on[dx]: # apply round if this index is true
    #         thisvar = group_vars[dx]
    #         round_to = round_digits[dx]
    #         var_rounded = df_cdr_ar[thisvar].round(round_to)
    #         # using same index for group_vars, round_digits, and round_on
    #         df_cdr_ar.loc[:, thisvar] = var_rounded
    

    # group by the variables listed
    df_grouped = df_cdr_ar.groupby(group_vars)
    if dust_dict['dust_from_file']:
        # select arbitrary set / var in dfdust for filtering 
        # (values are repeated for each instance, so doesn't matter which we choose, except for maybe site)
        thisset_dust = dfdust['set'].unique()[0]
        thisvar_dust = dfdust['var'].unique()[0]
    # collect summary data by group
    summary_data = []
    # loop through groups
    for group_values, group_df in df_grouped:

        # ------------------------------------------------------------------
        # get the dust flux from file if needed
        if dust_dict['dust_from_file']: 
            dustname = dust_dict['dustname_dustfile'] # dust name for the dfdust file
            dustname_group = group_vars[dust_dict['dustname_cdrfile_idx']] # dustname for the cdr file
            runname_name = group_vars[dust_dict['runname_idx']]
            # filter out runname and times that are > time_horizon
            dfdust_filt = dfdust.loc[(dfdust[runname_name] == group_df[runname_name].values[0]) & (dfdust['time'] <= time_horizon) & 
                                    (dfdust['set'] == thisset_dust) & (dfdust['var'] == thisvar_dust) & 
                                    (dfdust['site'] == thissite)]
            # integrate dust across time slice
            var_full = cumulative_trapezoid(dfdust_filt[dustname], dfdust_filt['time'], initial=0)
            if dust_dict['convert_to_tons']:
                #group_df[dustname_group] = (var_full/100).round(dust_dict['round_dust_digits'])[-1]
                tmp_dust_int = (var_full/100).round(dust_dict['round_dust_digits'])[-1] 
                tmp_dust_int_mean = tmp_dust_int / time_horizon
            else:
                #group_df[dustname_group] = var_full.round(dust_dict['round_dust_digits'])[-1]
                tmp_dust_int = var_full.round(dust_dict['round_dust_digits'])[-1]
                tmp_dust_int_mean = tmp_dust_int / time_horizon
        # ------------------------------------------------------------------

        # get integrated CDR at last timestep (multiply by time to get integrated value)
        ts_max_idx = group_df['time'].idxmax()  # get the index of max time step
        cdr_tmp = (group_df[cdvar]*group_df['time'])[ts_max_idx]  # the CDR value (time integrated)

        # combine the group_values with the calculated CDR
        summary_row = list(group_values) + [cdr_tmp]
        # update the dust flux
        if dust_dict['dust_from_file']:
            summary_row[dust_dict['dustname_cdrfile_idx']] = tmp_dust_int_mean

        # append the summary row to the data list
        summary_data.append(summary_row)

    # bring together in pd.DataFrame
    outdf = pd.DataFrame(summary_data, columns=group_vars + ['cdr'])
    # add duration and cdvar 
    outdf['timehorizon_yr'] = time_horizon
    outdf['cdvar'] = cdvar
    # add dust species for reference
    if "dustsp" in df_cdr.columns:
        outdf['dustsp'] = df_cdr['dustsp'].values[0]
    # return
    return outdf





# --- FUNCTION: CO2 emissions from crushing and transport (ignoring spreading) 
# based off of values / model in Zhang et al., 2023
def emissions_calculator(
    rockTons_ha_yr: float, dur: float, p80_output: float,
    p80_input: float, truck_km: float, barge_km: float, 
    barge_diesel_km: float, Efactor_org: str, mineral: str
    ) -> pd.DataFrame:
    """
    Read in information about amount and grain size of rock, as well as its transport
    pathway and starting grain size (for crushing) to return the CO2 emissions
    per rock, per ha per year, and per ha (over defined time horizon) as a pandas df

    Uses the constants / model of Zhang et al., 2023 
    (https://doi.org/10.1021/acs.est.3c01658)

    Parameters
    ----------
    # -- rocks -- # 
    rockTons_ha_yr : float
        [tons / ha / yr] rock application rate in tons / ha / yr
    dur : float
        [yr] duration of application (used with rockTons_ha_yr to derive total rock use)
    p80_output : float
        [microns] p80 diameter after crushing the rock (that used in the simulation)
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
        ["cc" | "gbas"] which mineral is being transported / crushed (this sets the bond work index)
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

    # --- MODEL
    # transport per ton of rock
    transport_emissions_perRock = (truck_factor * truck_km) + (barge_factor * barge_km) + (barge_diesel_factor * barge_diesel_km) # [kg CO2e / ton rock]
    # crushing per ton of rock
    crush_Efactor = crush_Efactors[Efactor_org]  # either the "MRO", "RFC", or "SERC" Efactor depending on "Efactor_org"
    crush_energy_perRock = 10 * bondwork_indices[mineral] * ((1 / math.sqrt(p80_output)) - (1 / math.sqrt(p80_input)))   # [kWh / ton rock]
    crush_emissions_perRock = crush_Efactor * crush_energy_perRock  # [kg CO2e / ton rock]
    # total per ton of rock
    total_perRock = transport_emissions_perRock + crush_emissions_perRock
    perRock_dict = {
        'transport': transport_emissions_perRock,
        'crushing': crush_emissions_perRock,
        'total': total_perRock,
        'units': "kgCO2e kgRock-1"
    }

    # transport per hectare per year 
    transport_emissions_perHaYr = transport_emissions_perRock * rockTons_ha_yr # [kg CO2e / ha / yr]
    # crushing per hectare per year
    crush_emissions_perHaYr = crush_emissions_perRock * rockTons_ha_yr  # [kg CO2e / ha / yr]
    # total per hectare per year
    total_perHaYr = transport_emissions_perHaYr + crush_emissions_perHaYr
    # add to dict
    perHaYr_dict = {
        'transport': transport_emissions_perHaYr,
        'crushing': crush_emissions_perHaYr,
        'total': total_perHaYr,
        'units': "kgCO2e ha-1 yr-1"
    }

    # transport per hectare (given total duration of rock applicatioin)
    transport_emissions_perHa = transport_emissions_perRock * (rockTons_ha_yr * dur)  # [kg CO2e / ha]
    # crushing per hectare (given total duration of rock application)
    crush_emissions_perHa = crush_emissions_perRock * (rockTons_ha_yr * dur)  # [kg CO2e / ha]
    # total per hectare 
    total_perHa = transport_emissions_perHa + crush_emissions_perHa
    # add to dict
    perHa_dict = {
        'transport': transport_emissions_perHa,
        'crushing': crush_emissions_perHa,
        'total': total_perHa,
        'units': "kgCO2e ha-1"
    }

    # --- collate results
    outdf = pd.DataFrame([perRock_dict, perHaYr_dict, perHa_dict])
    # return
    return outdf


# --- FUNCTION: Same as `emissions_calculator` but built to handle 
#     the pandas.DataFrame output from the cdr_per_group fxn
# based off of values / model in Zhang et al., 2023
def emissions_calculator_df( df : pd.DataFrame, dustrate_name: str,
    p80_input: float, truck_km: float, barge_km: float, 
    barge_diesel_km: float, Efactor_org: str, mineral: str
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
        ["cc" | "gbas" | None] which mineral is being transported / crushed (this sets the bond work index)
        (NOTE we take the "dustsp" from df by default if it's offered)
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

    # define the mineral (overwrite provided value with dataframe if we can)
    if "dustsp" in df.columns:
        mineral = df['dustsp'].values[0]

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
        dur = row['timehorizon_yr']

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

    # return
    return df



# FUNCTION to match control filename pattern to the full filename
def matches_pattern(fns, mylist):
    return any(fnmatch.fnmatch(fns, pattern) for pattern in mylist)



# FUNCTION to get net CDR from calcite and silicate comparison
# for a bunch of different silicate deployments and different
# CDR accounting approaches 
def cdr_batchCalc_sil_vs_cc(dfcc: pd.DataFrame, dfsil: pd.DataFrame,
                            dustrate_name: str,
                            cc_apprate_fixed: float, cc_dustrad_fixed: float,
                            include_cc_2d: bool=True,
                            select_nearest_ccdust: bool=True,
                            loss_percents:np.array = np.linspace(100,1,50))->xr.Dataset:
    """
    Reads in synthesized calcite and silicate weathering results from 
    SCEPTER plus the emissions data from emissions_calculator_df fxn
    and returns two xarray datasets with CDR and net CDR results.

    The input pd.DataFrames must have columns for dustrate_name (whatever we end up calling it), 
    'dustrad', 'cdr', 'E_total_tonCO2_ha'

    The output anomaly dataset (see below) includes four accounting options
    defined below, where R=removals, E=emissions, sil=silicate, cc=calcite:
    [1] net_R: Rsil - Rcc
    [2] simple_subtract: (Rsil - Rcc) - (Esil - Ecc)
    [3] no_negatives: (Rsil - MAX(0,Rcc)) - MAX(0, Esil - Ecc)
    [4] conservative: (Rsil - MAX(0,Rcc)) - Esil

    Parameters
    ----------
    dfcc : pd.DataFrame
        calcite dataframe with results from the emissions_calculator_df fxn
    dfsil : pd.DataFrame
        silicate dataframe with results from the emissions_calculator_df fxn
    dustrate_name : str
        name of the dustrate column in dfcc and dfsil
    cc_apprate_fixed : float
        [ton / ha / yr] counterfactual calcite application rate. Must be a value
        that exists in dfcc[dustrate_name]
    cc_dustrad_fixed : float
        [micron] counterfactual calcite dust radius. Must be a value that exists 
        in dfcc['dustrad']
    include_cc_2d : bool
        [True | False] whether we should make a 2d variable for cc removals in the `ds` dataset
        (note, i think this only works if cc is *also* defined over particle size *and* 
         application rate; and we assume cc and sil have the same apprates and dust rads)
    select_nearest_ccdust : bool
        [True | False] if cc_apprate_fixed is not one of the options, setting this to true
        will force the code to select the nearest cc apprate instaed. 
    loss_percents : np.arr
        1D numpy array of values for the loss_percent axis

    Returns
    -------
    xr.Dataset 1
        Variables are CDR due to silicate application (sil_R) and CDR due to 
        silicate application minus silicate-sourcing emissions (sil_R_E). 
        Dimensions are (1) the sil application rate; (2) the sil particle size;
        (3) the percent of CDR that's lost downstream
    
    xr.Dataset 2
        netCDR silicate minus aglime. Variables are netCDR due to accouting 
        option [1], [2], [3], and [4]. Dimensions are (1) the change in 
        application rate (sil-cc), (2) the change in particle size (sil-cc) 
        and (3) the percent of cdr that's lost downstream
    """
    # --- troubleshoot; check that application rate is in df_cc
    if not cc_apprate_fixed in dfcc[dustrate_name].unique():
        if select_nearest_ccdust:
            # get index where dustrate_name is closest to prescribed rate
            closest_index = (dfcc[dustrate_name] - cc_apprate_fixed).abs().idxmin()
            # select the corresponding value
            cc_apprate_fixed = dfcc.loc[closest_index, dustrate_name]
        else:
            available_dustrates = dfcc[dustrate_name].unique()
            raise ValueError(f"The CC application rate is not in the cc dataframe. Please select a value from the following: {available_dustrates}")
    if not cc_dustrad_fixed in dfcc['dustrad'].unique():
        available_dustrads = dfcc['dustrad'].unique()
        raise ValueError(f"The CC dust radius is not in the cc dataframe. Please select a value from the following: {available_dustrads}")
    # --- make sure that the two dfs have the same grid
    if set(dfcc[dustrate_name]) != set(dfsil[dustrate_name]):
        print("Warning: the dust rate columns in the cc and sil dfs do not match!!")
    if set(dfcc['dustrad']) != set(dfsil['dustrad']):
        print("Warning: the dust radius columns in the cc and sil dfs do not match!!")
    # ---

    # get counterfactual cdr and emissions (time integrated)
    cc_cdr_cf = dfcc[(dfcc[dustrate_name]==cc_apprate_fixed) & (dfcc['dustrad']==cc_dustrad_fixed)]['cdr'].values[0]
    cc_emiss_cf = dfcc[(dfcc[dustrate_name]==cc_apprate_fixed) & (dfcc['dustrad']==cc_dustrad_fixed)]['E_total_tonCO2_ha'].values[0]

    # create a 3D array with dimensions for loss percent ; app rate ; dust radius
    appRates = np.sort(dfsil[dustrate_name].unique())
    dustrads = np.sort(dfsil['dustrad'].unique())
    outarr_sil = np.zeros((len(loss_percents), len(appRates), len(dustrads)))
    outarr_sil_E = np.zeros((len(loss_percents), len(appRates), len(dustrads))) 
    outarr_sil_net = np.zeros((len(loss_percents), len(appRates), len(dustrads))) 
    if include_cc_2d:
        outarr_cc = np.zeros((len(loss_percents), len(appRates), len(dustrads)))
        outarr_cc_E = np.zeros((len(loss_percents), len(appRates), len(dustrads))) 
        outarr_cc_net = np.zeros((len(loss_percents), len(appRates), len(dustrads))) 
    # net removal accounting options
    # [1] net R: Rsil - Rcc
    outarr_dif1 = np.zeros((len(loss_percents), len(appRates), len(dustrads)))
    # [1a] net R no negative RCC: Rsil - max(0, Rcc)
    outarr_dif1a = np.zeros((len(loss_percents), len(appRates), len(dustrads)))
    # [2] simple subtraction: (Rsil - Rcc) - (Esil - Ecc)
    outarr_dif2 = np.zeros((len(loss_percents), len(appRates), len(dustrads)))
    # [3] No negative Rcc or E: (Rsil - MAX(0,Rcc)) - MAX(0, Esil - Ecc)
    outarr_dif3 = np.zeros((len(loss_percents), len(appRates), len(dustrads)))
    # [4] Conservative: (Rsil - MAX(0,Rcc)) - Esil
    outarr_dif4 = np.zeros((len(loss_percents), len(appRates), len(dustrads)))


    # indices for tracking
    lpdx, appdx, radx = 0,0,0

    # -- TROUBLESHOOT -- 
    print("n_iters: " + str(len(loss_percents) * len(appRates) * len(dustrads)))
    # --- LOOP THROUGH
    for loss_percent in loss_percents:
        for apprate in appRates:
            for dustrad in dustrads:
                # -- TROUBLESHOOT -- 
                # print(str(loss_percent) + " -- " + str(apprate) + " -- " + str(dustrad))
                # --

                # --- pull out integrated CDR and emissions
                # check if row exists in cc
                if include_cc_2d:
                    if len(dfcc[(dfcc[dustrate_name]==apprate) & (dfcc['dustrad']==dustrad)]['E_total_tonCO2_ha']) == 0:
                        outarr_cc[lpdx, appdx, radx], outarr_cc_E[lpdx, appdx, radx], outarr_cc_net[lpdx, appdx, radx] = np.nan, np.nan, np.nan
                    else:
                        # before loss
                        cc_cdr_tmp = dfcc[(dfcc[dustrate_name]==apprate) & (dfcc['dustrad']==dustrad)]['cdr'].values[0]
                        # apply loss
                        if cc_cdr_tmp > 0:
                            cc_cdr_tmp = cc_cdr_tmp - (2*(loss_percent/100)*cc_cdr_tmp)
                        else:
                            cc_cdr_tmp = cc_cdr_tmp 
                        # fill in arrays
                        cc_emiss = dfcc[(dfcc[dustrate_name]==apprate) & (dfcc['dustrad']==dustrad)]['E_total_tonCO2_ha'].values[0]
                        outarr_cc[lpdx, appdx, radx] = cc_cdr_tmp
                        outarr_cc_E[lpdx, appdx, radx] = cc_emiss
                        outarr_cc_net[lpdx, appdx, radx] = cc_cdr_tmp - cc_emiss

                # see if the row exists (sometimes it's missing if a run failed)
                if len(dfsil[(dfsil[dustrate_name]==apprate) & (dfsil['dustrad']==dustrad)]['E_total_tonCO2_ha']) == 0:
                    # then all arrays are nan here
                    outarr_sil[lpdx, appdx, radx], outarr_sil_E[lpdx, appdx, radx], outarr_sil_net[lpdx, appdx, radx] = np.nan, np.nan, np.nan
                    outarr_dif1[lpdx, appdx, radx], outarr_dif2[lpdx, appdx, radx] = np.nan, np.nan
                    outarr_dif3[lpdx, appdx, radx], outarr_dif4[lpdx, appdx, radx] = np.nan, np.nan
                                    
                else:
                    # otherwise finish computation
                    sil_emiss = dfsil[(dfsil[dustrate_name]==apprate) & (dfsil['dustrad']==dustrad)]['E_total_tonCO2_ha'].values[0]
                    # before loss
                    sil_cdr = dfsil[(dfsil[dustrate_name]==apprate) & (dfsil['dustrad']==dustrad)]['cdr'].values[0]
                    # apply loss to positive CDR 
                    # NOTE: (WIP) not totally sure this is correct... we're assuming that if CDR is negative in the field we 
                    # can ignore downstream losses... but i think we need a more explicit way of calculating this 
                    # (perhaps based on the advected DIC flux)
                    if sil_cdr > 0:
                        sil_cdr = sil_cdr - ((loss_percent/100)*sil_cdr)
                    if cc_cdr_cf > 0:
                        cc_cdr = cc_cdr_cf - (2*(loss_percent/100)*cc_cdr_cf)
                    else:
                        cc_cdr = cc_cdr_cf 
                    
                    # add to arrays (note we don't always have a cc array because it can have a constant appdx and radx...)
                    outarr_sil[lpdx, appdx, radx] = sil_cdr
                    outarr_sil_E[lpdx, appdx, radx] = sil_emiss
                    outarr_sil_net[lpdx, appdx, radx] = sil_cdr - sil_emiss
                    
                    # fill in the accounting arrays
                    # [1] net R
                    outarr_dif1[lpdx, appdx, radx] = sil_cdr - cc_cdr  
                    # [1a] net R no negatives
                    outarr_dif1a[lpdx, appdx, radx] = sil_cdr - max(0, cc_cdr)
                    # [2] simple subtraction
                    outarr_dif2[lpdx, appdx, radx] = (sil_cdr - cc_cdr) - (sil_emiss - cc_emiss_cf)
                    # [3] no negative Rcc or E
                    outarr_dif3[lpdx, appdx, radx] = (sil_cdr - max(0, cc_cdr)) - max(0, (sil_emiss - cc_emiss_cf))
                    # [4] conservative
                    outarr_dif4[lpdx, appdx, radx] = (sil_cdr - max(0, cc_cdr)) - sil_emiss

                # update indices
                radx += 1
                if radx >= len(dustrads):  # reset
                    radx = 0
            appdx += 1
            if appdx >= len(appRates):  # reset
                    appdx = 0
        lpdx += 1
        if lpdx >= len(loss_percents): # reset
            lpdx = 0
            
    # --- build the anomaly dataset
    # define coordinates
    D_appRate = appRates - cc_apprate_fixed
    D_dustrad = dustrads - cc_dustrad_fixed
    mycoords = {'loss_percent': loss_percents, 'D_apprate': D_appRate, 'D_dustrad': D_dustrad}
    ds_anom = xr.Dataset({
        'net_R': (['loss_percent', 'D_apprate', 'D_dustrad'], outarr_dif1),
        'net_R_no_negatives': (['loss_percent', 'D_apprate', 'D_dustrad'], outarr_dif1a),
        'simple_subtract': (['loss_percent', 'D_apprate', 'D_dustrad'], outarr_dif2),
        'no_negatives': (['loss_percent', 'D_apprate', 'D_dustrad'], outarr_dif3),
        'conservative': (['loss_percent', 'D_apprate', 'D_dustrad'], outarr_dif4),
        # add the cc counterfactual data 
        'cc_apprate_fixed': cc_apprate_fixed,
        'cc_dustrad_fixed': cc_dustrad_fixed,
        # and contextual dat
        'timehorizon_yr': dfcc['timehorizon_yr'].values[0],
        'cdvar': dfcc['cdvar'].values[0],
        # add emissions dat
        'p80_input_cc': dfcc['p80_input'].values[0],
        'truck_cc': dfcc['truck_km'].values[0],
        'barge_cc': dfcc['barge_km'].values[0],
        'barge_diesel_cc': dfcc['barge_diesel_km'].values[0],
        'Efactor_org_cc': dfcc['Efactor_org'].values[0],
        'bondwork_index_cc': dfcc['bondwork_index'].values[0],
        'p80_input_sil': dfsil['p80_input'].values[0],
        'truck_sil': dfsil['truck_km'].values[0],
        'barge_sil': dfsil['barge_km'].values[0],
        'barge_diesel_sil': dfsil['barge_diesel_km'].values[0],
        'Efactor_org_sil': dfsil['Efactor_org'].values[0],
        'bondwork_index_sil': dfsil['bondwork_index'].values[0]
    }, coords = mycoords)

    # --- build the sil dataset
    mycoords2 = {'loss_percent': loss_percents, 'apprate': appRates, 'dustrad': dustrads}
    ds = xr.Dataset({
        'sil_R': (['loss_percent', 'apprate', 'dustrad'], outarr_sil),
        'sil_E': (['loss_percent', 'apprate', 'dustrad'], outarr_sil_E),
        'sil_R_E': (['loss_percent', 'apprate', 'dustrad'], outarr_sil_net),
        # and contextual dat
        'timehorizon_yr': dfsil['timehorizon_yr'].values[0],
        'cdvar': dfsil['cdvar'].values[0],
        # add emissions dat
        'p80_input_sil': dfsil['p80_input'].values[0],
        'truck_sil': dfsil['truck_km'].values[0],
        'barge_sil': dfsil['barge_km'].values[0],
        'barge_diesel_sil': dfsil['barge_diesel_km'].values[0],
        'Efactor_org_sil': dfsil['Efactor_org'].values[0],
        'bondwork_index_sil': dfsil['bondwork_index'].values[0]
    }, coords = mycoords2)
    # --- build the cc dataset
    if include_cc_2d:
        dscc = xr.Dataset({
            'cc_R': (['loss_percent', 'apprate', 'dustrad'], outarr_cc),
            'cc_E': (['loss_percent', 'apprate', 'dustrad'], outarr_cc_E),
            'cc_R_E': (['loss_percent', 'apprate', 'dustrad'], outarr_cc_net),
            # and contextual dat
            'timehorizon_yr': dfcc['timehorizon_yr'].values[0],
            'cdvar': dfcc['cdvar'].values[0],
            # add emissions dat
            'p80_input_cc': dfcc['p80_input'].values[0],
            'truck_cc': dfcc['truck_km'].values[0],
            'barge_cc': dfcc['barge_km'].values[0],
            'barge_diesel_cc': dfcc['barge_diesel_km'].values[0],
            'Efactor_org_cc': dfcc['Efactor_org'].values[0],
            'bondwork_index_cc': dfcc['bondwork_index'].values[0]
        }, coords = mycoords2)
        # combine
        ds = xr.merge([ds, dscc])

    # return results
    return ds, ds_anom



# FUNCTION to save results of batch synthesis calculations
def save_batch_postproc(savedict: dict, base_path: str, base_dir_name: str,
                        fname_res: str, df_sil: pd.DataFrame, dfsil_sum: pd.DataFrame, 
                        df_cc: pd.DataFrame, dfcc_sum: pd.DataFrame, ds: xr.Dataset, 
                        ds_anom: xr.Dataset
                        ):
    """
    save the series of files generated when postprocessing the batch 
    results

    Parameters
    ----------
    savedict : dict
        dictionary of variables we will save as a resource file
    base_path : str
        directory where we will create the new directory to hold 
        the results
    base_dir_name : str
        name of the directory we want to create (the create_unique_directory
        function will make sure nothing is overwritten)
    fname_res : str
        name of the resource file that will be created from savedict
    df_sil : pd.DataFrame
        pandas dataframe of the full silicate results
    df_cc : pd.DataFrame
        pandas dataframe of the full calcite results
    dfsil_sum : pd.DataFrame
        pandas dataframe of the summary silicate results
    dfcc_sum : pd.DataFrame
        pandas dataframe of the summary calcite results
    ds : xr.Dataset 
        dataset of the silicate results vs loss_percent and group_vars
    ds_anom : xr.Dataset 
        dataset of the silicate minus cc anomaly for different CDR 
        accounting approaches
    """
    # [1] create the new directory
    outdir = create_unique_directory(base_path, base_dir_name)

    # [2] save the dict as a resource file
    save_variables_to_file(os.path.join(outdir, fname_res), savedict)

    # [3] save the pandas dataframes as a single pickle file
    dfs = {'df_sil': df_sil, 'dfsil_sum': dfsil_sum, 'df_cc': df_cc,
            'dfcc_sum': dfcc_sum}
    savedf_fn = os.path.join(outdir, "dataframes.pkl")
    # save results
    with open(savedf_fn, 'wb') as f:
        pickle.dump(dfs, f)
    # NOTE: you can open the dfs like this:
    # with open('dataframes.pkl', 'rb') as f:
        # loaded_dfs = pickle.load(f) 

    # [4] save xr.Datasets
    ds.to_netcdf(os.path.join(outdir, "ds.nc"))
    ds_anom.to_netcdf(os.path.join(outdir, "ds_anom.nc"))



# FUNCTION to create a save directory for batch results
def create_unique_directory(base_path: str, base_dir_name: str)-> str:
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
    i = 1
    while True:
        dir_name = os.path.join(base_path, f"{base_dir_name}_{i:03}")
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
            print(f"Directory created: {dir_name}")
            break
        i += 1
    return dir_name


def save_variables_to_file(filename: str, data: dict):
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
    with open(filename, 'w') as file:
        for label, variables in data.items():
            file.write(f"{label}:\n")
            for var_name, value in variables.items():
                file.write(f"  {var_name}: {value}\n")
            file.write("\n")  # Adds a blank line between tests


# FUNCTION create the dict we'll use to save results
def create_save_dict(thissite, cdvar, group_vars, csv_fn_sil,
                    csv_fn_cc, multiyear_sil, multiyear_cc, 
                    time_horizon, cc_apprate_fixed,
                    cc_dustrad_fixed, p80_input, truck_km,
                    barge_km, barge_diesel_km, Efactor_org)->dict:
    savedict = {
        "Setup": {
            "site": thissite,
            "cdvar": cdvar,
            "group_vars": group_vars,
            "csv_sil": csv_fn_sil,
            "csv_cc": csv_fn_cc,
            "multiyear_sil": multiyear_sil,
            "multiyear_cc": multiyear_cc
        },
        "Removal accounting": {
            "time_horizon": time_horizon,
            "cc_apprate_fixed": cc_apprate_fixed,
            "cc_dustrad_fixed": cc_dustrad_fixed
        },
        "Emissions": {
            "p80_input": p80_input,
            "truck_km": truck_km,
            "barge_km": barge_km,
            "barge_diesel_km": barge_diesel_km,
            "Efactor_org": Efactor_org
        }
    }

    return savedict
# %%
