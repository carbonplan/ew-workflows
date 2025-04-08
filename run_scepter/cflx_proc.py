# ------------------------------------------------
# 
# Functions to calculate / synthesize CDR-relevant
# fluxes for single SCEPTER simulations. 
# 
# In order to make the functions usable for a wide 
# range of SCEPTER run files, fluxes are NOT 
# compared to a control scenario. Only absolute 
# fluxes for the single run are considered.
# 
# These functions are called within the SCEPTER
# run .py scripts, creating a `postproc_*` dir 
# in the run output directory.
# 
# ------------------------------------------------
# %%
import math
import os
import re
import glob
import fnmatch
# from typing import Union, Literal

from scipy.integrate import cumulative_trapezoid
import numpy as np
import pandas as pd
import xarray as xr
import pickle

# %% 
# --- Dict for molar masses of solid species (from Kanzaki et al., 2022; table 1)
molar_mass_dict = {
    'amsi': 60.085,
    'qtz': 60.085,
    'gb': 78.004,
    'gt': 88.854,
    'hm': 159.692,
    'gps': 172.168,
    'arg': 100.089,
    'cc': 100.089,
    'dlm': 184.403,
    'ab': 262.225,
    'kfs': 278.33,
    'an': 278.311,
    'fo': 140.694,
    'fa': 203.778,
    'en': 100.389,
    'fer': 131.931,
    'dp': 216.553,
    'hb': 248.09,
    'tm': 812.374,
    'antp': 780.976,
    'mscv': 398.311,
    'plgp': 417.262,
    'ct': 277.113,
    'ka': 258.162,
    'anl': 220.155,
    'nph': 142.055,
    'nabd': 367.609,
    'kbd': 372.978,
    'cabd': 366.625,
    'mgbd': 363.996,
    'ill': 383.90,
    'g1': 30,
    'g2': 30,
    'g3': 30,
    'amnt': 80.043,  # Kanzaki et al., 2023, table 3
    # --- see treatment of gbas below (not in table 1 of Kanzaki et al., 2022)
    'gbas': 120.496
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
}


# %%
# --- FUNCTION to preprocess .txt files for consistent delimiters
def preprocess_txt(
    outdir: str,
    runname: str,
    fn: str,
    run_subdir: str = "flx",
    map_numeric: bool = True
) -> pd.DataFrame:
    """
    Convert SCEPTER output .txt files to pandas DataFrame for further
    analysis

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <RUNNAME>_field.
    fn : str
        name of the file in the flux directory to input
    run_subdir : str
        name of the subdirectory in the run output directory to find the file (either "prof" or "flux")
        default is "flx". 
    map_numeric : bool
        if True, all columns are mapped to numeric values. This should only be false for dust.txt which has 
        columns for the dust species

    Returns
    -------
    pd.DataFrame
        same format as the SCEPTER .txt file
    """
                    
    data = []  # Initialize a list to store the processed data

    # Initialize a flag to determine if we are reading the header
    is_header = True

    # get file path
    file_path = os.path.join(outdir, runname, run_subdir, fn)

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
    # apply pd.to_numeric to all columns using the "map" method
    if map_numeric:
        df = df.map(pd.to_numeric)
    # return
    return df

def get_data(
    outdir: str,
    runname: str,
    var_fn: str,
    cdvar: str,
    run_subdir: str = "flx",
    get_int: bool = True
) -> pd.DataFrame:
    """
    Get the SCEPTER data from defined text files. both timeseries
    and integrated flux files (e.g., <var_fn>-<cdvar> and int_<var_fn>-<cdvar>)
    are brought in. Both pandas dataframes are returned.
    

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <RUNNAME>_field.
    var_fn : str
        base name of the SCEPTER flx file. Format is "[basename]-[cdvar].txt"
    cdvar : str
        name of the variable used for CDR metric. Should align with var_fn such that format is "[basename]-[cdvar].txt"
    run_subdir : str
        name of the subdirectory in the run output directory to find the file (either "prof" or "flux")
        default is "flx".
    get_int : bool
        if True, then get the time-integrated version of the file and return that too

    Returns
    -------
    pd.DataFrame
        timeseries of the carbon flux metric
    pd.DataFrame
        integrated timeseries of the carbon flux metric
    """
    # get the txt file as a pandas dataframe
    fn = f"{var_fn}-{cdvar}.txt"   # flux timeseries
    df = preprocess_txt(outdir, runname, fn, run_subdir = run_subdir)
    if get_int:
        fn_int = f"int_{fn}"       # integrated flux timeseries
        dfint = preprocess_txt(outdir, runname, fn_int, run_subdir = run_subdir)
    # return result
    if get_int:
        return df, dfint
    else:
        return df, _

def get_data_prof(
    outdir: str,
    runname: str,
    var_prefix: str,
    time_index: int,
    run_subdir: str = "prof",
    multi_iter: bool = False,
    multi_iter_idx: int = 0,
) -> pd.DataFrame:
    '''
    Return data from a SCEPTER/scepter_output/prof/*.txt file 
    in pandas dataframe format.

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <NAME>_field.
    var_prefix : str
        prefix of the file to read in the 'prof' subdirectory. Includes all characters up to `*-xxx.txt`
        where xxx is a 3 digit numeric from 001 to 020. 
    time_index : int 
        [1:20] single integer between 1 and 20 denoting a timestep that was output as a profile. There is
        one profile file per timestep.
    run_subdir : str
        name of the subdirectory that holds the profile files ("prof" by default)
    multi_iter : bool 
        True if the simulation is multiple runs stitched together, False if it is 
        a single run (this matters for the file naming convention)
    multi_iter_idx : int
        If multi_iter is True, this integer denotes the index to append to the time_index 
        for the file name convention
    
    Returns 
    -------
    pd.DataFrame
        profile values defined over depth and time.
    '''
    # --- handle the multi-iter case 
    if multi_iter:
        fn = f"{var_prefix}-{time_index:03d}-{str(multi_iter_idx)}.txt"    # create filename
        infile = os.path.join(outdir, runname, run_subdir, fn)   # paste filename to path
        try:
            # read in data
            returnme = preprocess_txt(outdir, runname, fn, run_subdir = run_subdir)
            # rename z to depth
            returnme = returnme.rename(columns={"z": "depth"})

            # return result
            return returnme
        except:
            # print(f"Warning: preprocess_txt failed for profile data {runname} / {fn}")    
            return "nodata"

    # --- handle the standard case 
    else:
        # generate file path
        fn = f"{var_prefix}-{time_index:03d}.txt"    # create filename
        infile = os.path.join(outdir, runname, run_subdir, fn)   # paste filename to path
        # read in data
        returnme = preprocess_txt(outdir, runname, fn, run_subdir = run_subdir)
        # rename z to depth
        returnme = returnme.rename(columns={"z": "depth"})

        # return result
        return returnme


# %%
# ******************************************************************************
# ----------------------- CDR METRIC FUNCTIONS ---------------------------------
def co2_flx(
    outdir: str,
    runname: str,
    var_fn: str = "flx_gas",
    cdvar: str = "pco2",
    organic_sp_list: list = ["g1", "g2", "g3"],
    inorganic_sp_list: list = ["arg", "cc", "dlm"],
    convert_units: bool = True,
    co2_g_mol: float=44.01, # molar_mass_dict not used bc the input file should be in mol co2 / m2 / yr
) -> pd.DataFrame:
    """
    Get the CO2 diffusive and advective flux over time. Uses the *flx_gas* 
    files and pco2 variable by default. files are in mol/m2/yr and output 
    is in g/m2/yr if convert_units=False; ton/ha/yr if convert_units=True.
    
    You can use *flx_co2sp* files and DIC variable as well but it does 
    not include aqueous complexation between CO2 species and some cations, 
    so there can be some difference from the flx_gas results. 
    
    int_* fluxes are multiplied by time so the output is the time-integral.

    **************************************************************
    --- derivation notes ---
    the flux balance can be written as: adv = -dif - [sources] - tflx
    where positive values indicate net advection out. 
    
    if our only sources are respiration and cc dissolution, we get:
    adv = -dif - resp - cc - tflx

    so the advection flux without the inorganic contribution would simply be
    adv_noinorg = -dif - resp - tflx  
    [or] adv_noinorg = adv + cc
    **************************************************************

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <RUNNAME>_field.
    var_fn : str
        base name of the SCEPTER flx file. Format is "[basename]-[cdvar].txt"
    cdvar : str
        name of the variable used for CDR metric. Should align with var_fn such that format is "[basename]-[cdvar].txt"
    organic_sp_list : dict
        dict where keys are the IDs for all possible organic species (required for 
        resp calculation) and values are molar masses in g/mol 
        (required to go from mol/m2/yr to mass/m2/yr) (for keys and values, see Kanzaki 
        et al., 2022; table 1)
    inorganic_sp_list : dict
        dict where keys are the IDs for all carbon-bearing minerals 
        and values are their molar masses in g/mol (required to go 
        from mol/m2/yr to mass/m2/yr) (for keys and values, see Kanzaki 
        et al., 2022; table 1)
    convert_units : bool
        if True, convert mol/m2/yr to ton/ha/yr, if False convert mol/m2/yr to g/m2/yr
    co2_g_mol : float
        molar mass of CO2 [g mol-1] used to convert mol/m2/yr to g/m2/yr

    Returns
    -------
    pd.DataFrame
        timeseries of the carbon flux metric
    """
    # *****************************************
    # define unit conversion constants
    # for g/m2/yr to ton/ha/yr
    ton_g = 1 / 1e6   # [ton g-1]
    m2_ha = 10e3      # [m2 ha-1]
    conv_factor = ton_g * m2_ha 
    # *****************************************

    # get the txt file as a pandas dataframe 
    df, dfint = get_data(outdir, runname, var_fn, cdvar)

    # find how many organic species are present
    org_sp = list(set(organic_sp_list) & set(df.columns))
    org_sp_int = list(set(organic_sp_list) & set(dfint.columns))
    # find out how many carbon-bearing minerals are present
    inorg_sp = list(set(inorganic_sp_list) & set(df.columns))
    inorg_sp_int = list(set(inorganic_sp_list) & set(dfint.columns))

    # pull out the relevant components
    tdf = df.loc[:, ['time', 'dif', 'tflx', 'adv']].rename(columns={'time': 'time', 
                                                        'dif': 'co2flx_dif', 'tflx': 'co2flx_tflx',
                                                        'adv': 'co2flx_adv'})
    tdfint = dfint.loc[:, ['time', 'dif', 'tflx', 'adv']].rename(columns={'time': 'time', 
                                                        'dif': 'co2flx_dif', 'tflx': 'co2flx_tflx',
                                                        'adv': 'co2flx_adv'})

    
    # --- organic and inorganic contributions
    # compute resp component and add to tdf(int)
    if len(org_sp) > 0:  # if no org species are present, skip
        for sp in org_sp:
            tdf[sp] = df[sp]
        tdf['co2flx_resp'] = df[org_sp].sum(axis=1)  # add inorg fluxes together
    if len(org_sp_int) > 0:  # if no org species are present, skip
        for sp in org_sp_int:
            tdfint[sp] = dfint[sp]
        tdfint['co2flx_resp'] = dfint[org_sp_int].sum(axis=1)  # add inorg fluxes together
        
    # compute inorg component and add to tdf(int)
    if len(inorg_sp) > 0:  # if no org species are present, skip
        for sp in inorg_sp:
            tdf[sp] = df[sp]
        tdf['co2flx_inorg'] = df[inorg_sp].sum(axis=1)  # add inorg fluxes together
        tdf['co2flx_adv_noinorg'] = tdf['co2flx_adv'] + tdf['co2flx_inorg']
    if len(inorg_sp_int) > 0:  # if no org species are present, skip
        for sp in inorg_sp_int:
            tdfint[sp] = dfint[sp]
        tdfint['co2flx_inorg'] = dfint[inorg_sp_int].sum(axis=1)  # add inorg fluxes together
        tdfint['co2flx_adv_noinorg'] = tdfint['co2flx_adv'] + tdfint['co2flx_inorg']
    # --- convert units 
    if convert_units:  # convert mol/m2/yr to ton/ha/yr
        tdfint = tdfint.apply(lambda x: x * co2_g_mol * conv_factor if x.name not in ['time', 'units'] else x)
        tdf = tdf.apply(lambda x: x * co2_g_mol * conv_factor if x.name not in ['time', 'units'] else x)
        # add units columns
        tdfint['units'] = 'ton ha-1'
        tdf['units'] = 'ton ha-1 yr-1'
    else: # convert mol/m2/yr to g/m2/yr
        tdfint = tdfint.apply(lambda x: x * co2_g_mol if x.name not in ['time', 'units'] else x)
        tdf = tdf.apply(lambda x: x * co2_g_mol if x.name not in ['time', 'units'] else x)
        # add units columns
        tdfint['units'] = 'g m-2'
        tdf['units'] = 'g m-2 yr-1'
    
    # multiply tdfint columns by time (required to output time-integrated fluxes)
    tdfint = tdfint.apply(lambda x: x * tdfint['time'] if x.name not in ['time', 'units'] else x)

    
    # --- tidy up output
    tdfint['flx_type'] = "int_flx"
    tdf['flx_type'] = "flx"
    tdfint['runname'] = tdf['runname'] = runname
    tdfint['var'] = tdf['var'] = cdvar
    
    # combine 
    outdf = pd.concat([tdf, tdfint], axis=0, ignore_index=True)
    return outdf


def sld_flx(
    outdir: str,
    runname: str,
    feedstock: str,
    var_fn: str = "flx_sld",
    dust_from_file: bool = True,
    molar_mass_dict: dict = molar_mass_dict
) -> pd.DataFrame:
    """
    Get the feedstock fluxes. Uses the *flx_sld* files for the feedstock. Only
    time-integrated fluxes are returned. integrated dust application can be 
    computed from the dust.txt file if "dust_from_file"=True (this is 
    required for re-application runs which are composites of multiple 
    individual runs)


    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <RUNNAME>_field.
    feedstock : str
        name of the feedstock of rock applied to the system. Format is "[basename]-[feedstock].txt"
    var_fn : str
        base name of the SCEPTER flx file. Format is "[basename]-[cdvar].txt"
    dust_from_file : bool
        if True, then compute integrated dust fluxes from the dust.txt file. Should be 
        set to True for re-application runs since the int_sld* files are patched 
        together from multiple 1-year runs, so they don't reflect the true time integral
        (just the timeSTEP integral). This ensures that the integral is correct
    molar_mass_dict : dict
        dictionary where keys are species IDs and values are molar masses in [g mol-1]
        taken from Kanzaki et al., 2022 GMD; Table 1. 

    Returns
    -------
    pd.DataFrame
        timeseries of the carbon flux metric
    """
    # *****************************************
    # define unit conversion constants
    # for g/m2/yr to ton/ha/yr
    ton_g = 1 / 1e6   # [ton g-1]
    m2_ha = 10e3      # [m2 ha-1]
    conv_factor = ton_g * m2_ha 
    # *****************************************

    # get the txt file as a pandas dataframe 
    df, dfint = get_data(outdir, runname, var_fn, cdvar=feedstock)

    # multiply tdfint columns by time (required to output time-integrated fluxes)
    dfint = dfint.apply(lambda x: x * dfint['time'] if x.name != 'time' else x)
    
    # add dust if needed
    # [note] no unit conversion from mol/m2/yr to g/m2/yr is needed if 
    # dust_from_file because that is calculated based on application in g/m2/yr 
    if dust_from_file: 
        dfdust0 = preprocess_txt(outdir, runname, fn='dust.txt',
                                run_subdir = "flx", map_numeric = False)
        # map numeric columns to numerics
        dfnum = dfdust0.drop(columns=['dustsp1', 'dustsp2']).map(pd.to_numeric)
        # add other columns back
        dfdust = dfnum.copy()
        dfdust['dustsp1'] = dfdust0['dustsp1'].copy() 
        dfdust['dustsp2'] = dfdust0['dustsp2'].copy() 
        
        # --- identify the feedstock index
        columns_with_fs = [col for col in dfdust.columns if (dfdust[col] == feedstock).any()]
        if len(columns_with_fs) > 0:
            fscol = columns_with_fs[0] # should only return one column, so we take the 0 index
        else:
            fscol = "not found"
        if fscol == "dustsp1":
            dust_dx = '1'
        elif fscol == "dustsp2":
            dust_dx = '2'
        else:    # assume it's the first index if it doesn't exist
            dust_dx = '1'
        # re-calculate integral because right now it's integrated by timestep, not
        # by the entire run itself
        dfdust['int_dust_g_m2_yr'] = cumulative_trapezoid(dfdust[f'dust{dust_dx}_g_m2_yr'], dfdust['time'], initial=0)
        # add dust data 
        if len(dfdust) != len(df): # if mis-matched, then match them
            # drop duplicates in the 'time' column, keeping first occurrence
            dfdust_nodup = dfdust.drop_duplicates(subset='time', keep='first')
            # then interpolate the data 
            # merge the DataFrames on 'time' using outer join to keep all time points
            df_merged = pd.merge(df, dfdust_nodup, on='time', how='outer', suffixes=('', '_orig'))
            # interpolate the columns of interest 
            df_merged['int_dust_g_m2_yr'] = df_merged['int_dust_g_m2_yr'].interpolate()
            # keep only the points in the df time steps
            intdust = df_merged[df_merged['time'].isin(df['time'])]['int_dust_g_m2_yr'].values / 100 # divide by 100 to convert g/m2/yr to ton/ha/yr
        else:
            intdust = dfdust['int_dust_g_m2_yr'].values / 100 # divide by 100 to convert g/m2/yr to ton/ha/yr
        # add integrated dust
        dfint['int_dust_ton_ha_yr'] = intdust
    else:
        # use rain dust (already multiplied by time, so it's time-integrated)
        # multiply by -1 to get positive values into the soil column
        # multiply by molar_mass_dict[feedstock] to get mol/m2/yr in g/m2/yr
        # multiply by conv_factor go get g/m2/yr to ton/ha/yr
        dfint['int_dust_ton_ha_yr'] = -1 * dfint['rain'] * molar_mass_dict[feedstock] * conv_factor # (note, we multiplied rain by time earlier so we don't have to do it here)

    # convert other variables to ton/ha/yr 
    dfint['adv'] = dfint['adv'] * molar_mass_dict[feedstock] * conv_factor
    dfint[feedstock] = dfint[feedstock] * molar_mass_dict[feedstock] * conv_factor

    # pull out just the columns we want
    tdfint = dfint.loc[:, ['time', 'int_dust_ton_ha_yr', 'adv', feedstock]]
    # compute dissolution
    tdfint['dust_minus_adv'] = tdfint['int_dust_ton_ha_yr'] - tdfint['adv'] # dust that's left after solid advection
    tdfint['total_dissolution'] = tdfint[feedstock]  # net dissolution
    tdfint['fraction_sld_advected'] = tdfint['adv'] / tdfint['int_dust_ton_ha_yr']
    tdfint['fraction_sld_remaining'] = (tdfint['int_dust_ton_ha_yr'] - tdfint['adv'] - tdfint[feedstock]) / tdfint['int_dust_ton_ha_yr']
    # fraction of non-advected rock that gets dissolved
    tdfint['fraction_remaining_dissolved'] = tdfint['total_dissolution'] / tdfint['dust_minus_adv']
    # fraction of total applied rock that gets dissolved
    tdfint['fraction_total_dissolved'] = tdfint['total_dissolution'] / tdfint['int_dust_ton_ha_yr']

    # --- return result
    return tdfint
    
    
def carbAlk_adv(
    outdir: str,
    runname: str,
    var_fn: str = "flx_co2sp",
    cdvar: str = "ALK",
    convert_units: bool = True,
    co2potential_g_mol_sil: float=88.02, # potential grams of CO2 per mole of alkalinity assuming 2:1 DIC_CO2:ALK
    co2potential_g_mol_cc: float=44.01, # potential grams of CO2 per mole of alkalinity assuming 1:1 DIC_CO2:ALK
) -> pd.DataFrame:
    """
    Get the advective and storage fluxes of carbonate alkalinity over time. 
    Uses the *flx_co2sp* files and ALK variable by default. files are in 
    mol/m2/yr and output is same or in [mass]/m2/yr for co2potential.
    
    carbonate alkalinity is defined as: ALK = [HCO3]- + 2[CO3]--. 
    
    int_* fluxes are multiplied by time so the output is the time-integral.

    Note, diffusive flux is generally negligibly small and source / sink 
    contributions are not computed in the input files so they're ignored here. 

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <RUNNAME>_field.
    var_fn : str
        base name of the SCEPTER flx file. Format is "[basename]-[cdvar].txt"
    cdvar : str
        name of the variable used for CDR metric. Should align with var_fn such that format is "[basename]-[cdvar].txt"
    convert_units : bool
        if True, convert mol/m2/yr to ton/ha/yr, if False convert mol/m2/yr to g/m2/yr
    co2potential_g_mol_sil : float
        molar mass of CO2 [g mol-1] used to convert mol/m2/yr to g/m2/yr. This refers to the potential 
        CO2 sequestration by a silicate feedstock where you could get two moles of DIC from CO2 per mole 
        alkalinity
    co2potential_g_mol_cc : float
        molar mass of CO2 [g mol-1] used to convert mol/m2/yr to g/m2/yr. This refers to the potential 
        CO2 sequestration by a carbonate feedstock where you could get one mole of DIC from CO2 per mole 
        alkalinity

    Returns
    -------
    pd.DataFrame
        timeseries of the carbonate alkalinity flux metrics
    """
    # *****************************************
    # define unit conversion constants
    # for g/m2/yr to ton/ha/yr
    ton_g = 1 / 1e6   # [ton g-1]
    m2_ha = 10e3      # [m2 ha-1]
    conv_factor = ton_g * m2_ha 
    # *****************************************

    # get the txt file as a pandas dataframe 
    df, dfint = get_data(outdir, runname, var_fn, cdvar)

    # pull out the relevant components
    tdf = df.loc[:, ['time', 'tflx', 'adv']].rename(columns={'time': 'time', 
                                                        'tflx': 'calkflx_tflx',
                                                        'adv': 'calkflx_adv'})
    tdfint = dfint.loc[:, ['time', 'tflx', 'adv']].rename(columns={'time': 'time', 
                                                        'tflx': 'calkflx_tflx',
                                                        'adv': 'calkflx_adv'})
    # get total carbonate alkalinity flux (add the storage back in)
    tdf['calkflx_tot'] = tdf['calkflx_adv'] + tdf['calkflx_tflx']
    tdfint['calkflx_tot'] = tdfint['calkflx_adv'] + tdfint['calkflx_tflx']

    # compute potential flux
    if convert_units:  # get co2 potentials in ton/ha/yr
        # advective flx
        tdf['co2pot_adv_tonHaYr_sil'] = tdf['calkflx_adv'] * co2potential_g_mol_sil * conv_factor
        tdf['co2pot_adv_tonHaYr_cc'] = tdf['calkflx_adv'] * co2potential_g_mol_cc * conv_factor
        tdfint['co2pot_adv_tonHa_sil'] = tdfint['calkflx_adv'] * co2potential_g_mol_sil * conv_factor
        tdfint['co2pot_adv_tonHa_cc'] = tdfint['calkflx_adv'] * co2potential_g_mol_cc * conv_factor
        # advective + storage flx
        tdf['co2pot_tot_tonHaYr_sil'] = tdf['calkflx_tot'] * co2potential_g_mol_sil * conv_factor
        tdf['co2pot_tot_tonHaYr_cc'] = tdf['calkflx_tot'] * co2potential_g_mol_cc * conv_factor
        tdfint['co2pot_tot_tonHa_sil'] = tdfint['calkflx_tot'] * co2potential_g_mol_sil * conv_factor
        tdfint['co2pot_tot_tonHa_cc'] = tdfint['calkflx_tot'] * co2potential_g_mol_cc * conv_factor

    else: # get co2 potentials in g/m2/yr
        # advective flx
        tdf['co2pot_adv_gm2Yr_sil'] = tdf['calkflx_adv'] * co2potential_g_mol_sil 
        tdf['co2pot_adv_gm2Yr_cc'] = tdf['calkflx_adv'] * co2potential_g_mol_cc 
        tdfint['co2pot_adv_gm2_sil'] = tdfint['calkflx_adv'] * co2potential_g_mol_sil 
        tdfint['co2pot_adv_gm2_cc'] = tdfint['calkflx_adv'] * co2potential_g_mol_cc 
        # advective + storage flx
        tdf['co2pot_tot_gm2Yr_sil'] = tdf['calkflx_tot'] * co2potential_g_mol_sil 
        tdf['co2pot_tot_gm2Yr_cc'] = tdf['calkflx_tot'] * co2potential_g_mol_cc 
        tdfint['co2pot_tot_gm2_sil'] = tdfint['calkflx_tot'] * co2potential_g_mol_sil 
        tdfint['co2pot_tot_gm2_cc'] = tdfint['calkflx_tot'] * co2potential_g_mol_cc 

    # multiply tdfint columns by time (required to output time-integrated fluxes)
    tdfint = tdfint.apply(lambda x: x * tdfint['time'] if x.name not in ['time'] else x)

    
    # --- tidy up output
    tdf['units'] = 'mol m-2 yr'
    tdfint['units'] = 'mol m-2'
    tdfint['flx_type'] = "int_flx"
    tdf['flx_type'] = "flx"
    tdfint['runname'] = tdf['runname'] = runname
    tdfint['var'] = tdf['var'] = cdvar
    
    # combine 
    outdf = pd.concat([tdf, tdfint], axis=0, ignore_index=True)
    return outdf


def sumCat_adv(
    outdir: str,
    runname: str,
    var_fn: str = "flx_aq",
    catvars_charge: dict = {"ca": 2, "mg": 2, "k": 1, "na": 1},
    convert_units: bool = True
) -> pd.DataFrame:
    """
    Get the advective and storage fluxes of the sum of cations over time. 
    Uses the *flx_aq* files and cation variable by default. files are in 
    mol/m2/yr and output is same or in [mass]/m2/yr for co2potential. 
    
    int_* fluxes are multiplied by time so the output is the time-integral.

    Note, diffusive flux is generally negligibly small and source / sink 
    contributions are not computed in the input files so they're ignored here. 

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <RUNNAME>_field.
    var_fn : str
        base name of the SCEPTER flx file. Format is "[basename]-[cdvar].txt"
    catvars_charge : dict
        dict of variable used for cation metric. key is the variable, value is the charge. 
        Should align with var_fn such that format is "[basename]-[cdvar].txt"
    convert_units : bool
        if True, convert mol/m2/yr to ton/ha/yr, if False convert mol/m2/yr to g/m2/yr

    Returns
    -------
    dict
        dict of pd.DataFrames, one for each species
    pd.DataFrame
        summary timeseries cation fluxes
    """
    # *****************************************
    # define unit conversion constants
    # for g/m2/yr to ton/ha/yr
    ton_g = 1 / 1e6   # [ton g-1]
    m2_ha = 10e3      # [m2 ha-1]
    conv_factor = ton_g * m2_ha 
    # *****************************************

    # find which cations are tracked
    catvars_exist = find_cations(outdir, runname, var_fn=var_fn, catvars_charge=catvars_charge)
    
    # loop through cation files and compile
    out_dict = {}
    for cat in catvars_exist:
        # get cation charge
        ccharge = catvars_charge[cat]
        # get the txt file as a pandas dataframe 
        tmpdf, tmpdfint = get_data(outdir, runname, var_fn, cat)
        # get the summary file for this cation that we'll save 
        outdf_cat, nonintdf, intdf = build_cation_df(runname, cat, ccharge, tmpdf, tmpdfint)
        # add to the output dictionary
        out_dict[cat] = outdf_cat


        # get the columns we want and append with charge because we'll multiply them by their charge shortly
        tdf = nonintdf.loc[:, ['time', 'tflx', 'adv', 'carbsld_source', 'noncarbsld_source']].rename(
            columns={'time': 'time', 'tflx': 'tflx_charge', 'adv': 'adv_charge',
                    'carbsld_source': 'carbsld_source_charge', 'noncarbsld_source': 'noncarbsld_source_charge'}
        )
        tdfint = intdf.loc[:, ['time', 'tflx', 'adv', 'carbsld_source', 'noncarbsld_source']].rename(
            columns={'time': 'time', 'tflx': 'tflx_charge', 'adv': 'adv_charge',
                    'carbsld_source': 'carbsld_source_charge', 'noncarbsld_source': 'noncarbsld_source_charge'}
        )

        # get their charge flux
        tdf = tdf.apply(lambda x: x * ccharge if x.name not in ['time'] else x)
        tdfint = tdfint.apply(lambda x: x * ccharge if x.name not in ['time'] else x)

        # --- compute charge adjusted cation flux
        # advective flx (sum of -tflx -noncarbsld_source -carbsld_source/2)
        # compute potential flux
        if convert_units:  # get co2 potentials in ton/ha/yr
            # advective flx (sum of -tflx -noncarbsld_source -carbsld_source/2)
            tdf['co2pot_adv_tonHaYr'] = nonintdf['co2pot_adv_tonHaYr']
            tdfint['co2pot_adv_tonHa'] = intdf['co2pot_adv_tonHa']
            # advective plus storage flx (sum of -noncarbsld_source -carbsld_source/2)
            tdf['co2pot_tot_tonHaYr'] = nonintdf['co2pot_tot_tonHaYr']
            tdfint['co2pot_tot_tonHa'] = intdf['co2pot_tot_tonHa']

        else: # get co2 potentials in g/m2/yr
            # advective flx (sum of -tflx -noncarbsld_source -carbsld_source/2)
            tdf['co2pot_adv_gm2Yr'] = nonintdf['co2pot_adv_gm2Yr']
            tdfint['co2pot_adv_gm2'] = intdf['co2pot_adv_gm2']
            # advective plus storage flx (sum of -noncarbsld_source -carbsld_source/2)
            tdf['co2pot_tot_gm2Yr'] = nonintdf['co2pot_tot_gm2Yr']
            tdfint['co2pot_tot_gm2'] = intdf['co2pot_tot_gm2']

        # if it's the first loop, create the output dataframe
        if cat == catvars_exist[0]:
            outdf = tdf.copy()
            outdfint = tdfint.copy()
        else:
            # add all columns except time (effectively adding the new cation into the mix)
            outdf = outdf.drop(columns='time').copy() + tdf.drop(columns='time').copy()
            outdfint = outdfint.drop(columns='time').copy() + tdfint.drop(columns='time').copy()
            # add time back in
            outdf.insert(0, 'time', tdf['time'])
            outdfint.insert(0, 'time', tdfint['time'])
            
    # multiply tdfint columns by time (required to output time-integrated fluxes)
    outdfint = outdfint.apply(lambda x: x * outdfint['time'] if x.name not in ['time'] else x)

    # --- tidy up output
    outdf['units'] = 'mol m-2 yr x charge'
    outdfint['units'] = 'mol m-2 x charge'
    outdfint['flx_type'] = "int_flx"
    outdf['flx_type'] = "flx"
    outdfint['runname'] = outdf['runname'] = runname
    outdfint['var'] = outdf['var'] = '+'.join(catvars_exist)
    
    # combine 
    outdfx = pd.concat([outdf, outdfint], axis=0, ignore_index=True)
    return outdfx, out_dict
        

def find_cations(
    outdir: str,
    runname: str,
    var_fn: str = "flx_aq",
    catvars_charge: dict = {"ca": 2, "mg": 2, "k": 1, "na": 1}
) -> list:
    """
    Get a list of cation files present in the flux output.

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <RUNNAME>_field.
    var_fn : str
        base name of the SCEPTER flx file. Format is "[basename]-[cdvar].txt"
    catvars_charge : dict
        dict of variable used for cation metric. key is the variable, value is the charge. 
        Should align with var_fn such that format is "[basename]-[cdvar].txt"
    
    Returns
    -------
    list 
        list of catvars that are presenet in the output. 
    """
    # create empty list to hold existing catvars
    outlist = []
    for key, val in catvars_charge.items():
        # build the path / fn
        fn = f"{var_fn}-{key}.txt"   # flux timeseries
        file_path = os.path.join(outdir, runname, "flx", fn)
        if os.path.exists(file_path):
            outlist.append(key)
    # return result
    return outlist


def build_cation_df(
    runname: str, 
    cation: str,
    ccharge: float,
    tmpdf: pd.DataFrame, 
    tmpdfint: pd.DataFrame,
    inorganic_sp_list: list = ['cc', 'dlm', 'arg'],
    negligible_val: float = 1e-7,
    co2potential_g_mol: float=44.01,
    convert_units: bool=True
)-> pd.DataFrame:
    """
    read in the two cation dataframes, clean them up, return them

    Parameters
    ----------
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <RUNNAME>_field.
    cation : str
        name of the cation in the cation list
    ccharge : float
        value of cation charge
    tmpdf : pd.DataFrame
        dataframe generated in the for loop of sumCat_adv function
    tmpdfint : pd.DataFrame
        time integrated dataframe generated in the for loop of sumCat_adv function
    inorganic_sp_list : list
        list of all the inorganic species that we deduct half c potential removal from
    negligible_val : float
        if any column has all its values below this value, we remove it from the df
    co2potential_g_mol : float
        molar mass of co2
    convert_units : bool
        if True, convert mol/m2/yr to ton/ha/yr, if False convert mol/m2/yr to g/m2/yr

    Returns
    -------
    pd.DataFrame
        one for tmpdf, one for tmpdfint
    """
    # *****************************************
    # define unit conversion constants
    # for g/m2/yr to ton/ha/yr
    ton_g = 1 / 1e6   # [ton g-1]
    m2_ha = 10e3      # [m2 ha-1]
    conv_factor = ton_g * m2_ha 
    # *****************************************

    # list of columns we cannot trim
    columns_to_keep = ["time", "adv", "tflx", "res"]
    # remove empty or negligible columns (to simplify) (but keep "time" and "res" no matter what)
    tmpdf = tmpdf.loc[:, (np.abs(tmpdf) > 1e-7).any(axis=0) | tmpdf.columns.isin(columns_to_keep)].copy()
    tmpdfint = tmpdfint.loc[:, (np.abs(tmpdfint) > 1e-7).any(axis=0) | tmpdfint.columns.isin(columns_to_keep)].copy()
    
    # find out how many carbon-bearing minerals are present
    inorg_sp = list(set(inorganic_sp_list) & set(tmpdf.columns))
    inorg_sp_int = list(set(inorganic_sp_list) & set(tmpdfint.columns))
    # get other sources
    exclude_cols = inorg_sp + ['time', 'tflx', 'adv', 'res']
    exclude_cols_int = inorg_sp_int + ['time', 'tflx', 'adv', 'res']
    noninorg_sources = tmpdf.columns.difference(exclude_cols, sort=False)
    noninorg_sources_int = tmpdfint.columns.difference(exclude_cols, sort=False)

    # add columns for individual sources times the charge
    # -- noninorganic sources
    if len(noninorg_sources) > 0:
        tmpdf.loc[:, 'noncarbsld_source'] = tmpdf[noninorg_sources].sum(axis=1) 
    else:
        tmpdf['noncarbsld_source'] = 0.
    if len(noninorg_sources_int) > 0:
        tmpdfint.loc[:, 'noncarbsld_source'] = tmpdfint[noninorg_sources_int].sum(axis=1) 
    else:
        tmpdfint['noncarbsld_source'] = 0.
    
    # -- inorganic sources
    if len(inorg_sp) > 0:
        tmpdf.loc[:, 'carbsld_source'] = tmpdf[inorg_sp].sum(axis=1)
    else:
        tmpdf['carbsld_source'] = 0.
    if len(inorg_sp_int) > 0:
        tmpdfint.loc[:, 'carbsld_source'] = tmpdfint[inorg_sp_int].sum(axis=1)
    else:
        tmpdfint['carbsld_source'] = 0.
    # advective flx (sum of -tflx -noncarbsld_source -carbsld_source/2)
    adv_cat = -tmpdf['tflx'] - tmpdf['noncarbsld_source'] - tmpdf['carbsld_source']/2
    adv_cat_int = -tmpdfint['tflx'] - tmpdfint['noncarbsld_source'] - tmpdfint['carbsld_source']/2
    # advective plus storage flx (sum of -noncarbsld_source -carbsld_source/2)
    tot_cat = -tmpdf['noncarbsld_source'] - tmpdf['carbsld_source']/2
    tot_cat_int = -tmpdfint['noncarbsld_source'] - tmpdfint['carbsld_source']/2

    # compute potential flux
    if convert_units:  # get co2 potentials in ton/ha/yr
        # advective flx (sum of -tflx -noncarbsld_source -carbsld_source/2)
        tmpdf['co2pot_adv_tonHaYr'] = adv_cat * ccharge * co2potential_g_mol * conv_factor
        tmpdfint['co2pot_adv_tonHa'] = adv_cat_int * ccharge * co2potential_g_mol * conv_factor
        # advective plus storage flx (sum of -noncarbsld_source -carbsld_source/2)
        tmpdf['co2pot_tot_tonHaYr'] = tot_cat * ccharge * co2potential_g_mol * conv_factor
        tmpdfint['co2pot_tot_tonHa'] = tot_cat_int * ccharge * co2potential_g_mol * conv_factor

    else: # get co2 potentials in g/m2/yr
        # advective flx (sum of -tflx -noncarbsld_source -carbsld_source/2)
        tmpdf['co2pot_adv_gm2Yr'] = adv_cat * ccharge * co2potential_g_mol 
        tmpdfint['co2pot_adv_gm2'] = adv_cat_int * ccharge * co2potential_g_mol 
        # advective plus storage flx (sum of -noncarbsld_source -carbsld_source/2)
        tmpdf['co2pot_tot_gm2Yr'] = tot_cat * ccharge * co2potential_g_mol 
        tmpdfint['co2pot_tot_gm2'] = tot_cat_int * ccharge * co2potential_g_mol 

    # add qualitative columns
    tmpdf['units'] = 'mol m-2 yr'
    tmpdfint['units'] = 'mol m-2'
    tmpdfint['flx_type'] = "int_flx"
    tmpdf['flx_type'] = "flx"
    tmpdfint['runname'] = tmpdf['runname'] = runname
    tmpdfint['cation'] = tmpdf['cation'] = cation
    tmpdfint['charge'] = tmpdf['charge'] = ccharge

    # --- return
    # combine 
    outdf = pd.concat([tmpdf, tmpdfint], axis=0, ignore_index=True)
    return outdf, tmpdf, tmpdfint




# %% 
# ******************************************************************************
# -------------------------- CDR METRIC MAIN -----------------------------------
def cflx_calc(
    outdir: str,
    runname: str,
    feedstock: str,
    dust_from_file: bool = True,
    convert_units: bool = True,
    save_dir: str = 'postproc_flxs',
    calc_list: list = ['co2_flx', 'carbAlk_adv', 'sumCat_adv', 'sld_flx']
):
    """
    compute various cdr-relevant c flux metrics in the postprocessing
    of a SCEPTER run. Save the results

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <RUNNAME>_field.
    feedstock : str
        ID of the feedstock used (for rock dissolution calculation)
    dust_from_file : bool
        if True, then compute integrated dust fluxes from the dust.txt file. Should be 
        set to True for re-application runs since the int_sld* files are patched 
        together from multiple 1-year runs, so they don't reflect the true time integral
        (just the timeSTEP integral). This ensures that the integral is correct
    convert_units : bool
        if True, convert mol/m2/yr to ton/ha/yr, if False convert mol/m2/yr to g/m2/yr
    save_dir : str
        the name of the subdirectory in the run dir where the cflux files get saved
    calc_list : list
        list of all the functions we're going to run

    Results
    -------
    [ no output - fxn saves the results ]
    """
    # where to save the results
    savehere = os.path.join(outdir, runname, save_dir)
    # make dir if it doesn't exist
    if not os.path.exists(savehere):
        os.makedirs(savehere)


    # [1] --- CO2 FLUXES ------------------------------------------------
    if "co2_flx" in calc_list:
        co2df = co2_flx(outdir, runname, var_fn = "flx_gas", cdvar= "pco2", 
                        organic_sp_list = ["g1", "g2", "g3"],
                        inorganic_sp_list = ["arg", "cc", "dlm"],
                        convert_units = convert_units,
                        co2_g_mol = 44.01
                        ) 
        # save
        savename = "co2_flxs.pkl"
        co2df.to_pickle(os.path.join(savehere, savename))

    
    # [2] --- CARBONATE ALKALINITY FLUXES ---------------------------------
    if "carbAlk_adv" in calc_list:
        cAlkdf = carbAlk_adv(outdir, runname, var_fn = "flx_co2sp", cdvar = "ALK",
                            convert_units = convert_units,
                            co2potential_g_mol_sil = 88.02, 
                            co2potential_g_mol_cc = 44.01
                            )
        # save
        savename = "carbAlk_flxs.pkl"
        cAlkdf.to_pickle(os.path.join(savehere, savename))


    # [3] --- SUM OF CATION FLUXES -----------------------------------------
    if "sumCat_adv" in calc_list: 
        sumCatdf, sumCatdict = sumCat_adv(outdir, runname, var_fn = "flx_aq",
                                          catvars_charge = {"ca": 2, "mg": 2, "k": 1, "na": 1},
                                          convert_units = True)
        # save the main file
        savename_main = "cationflx_sum.pkl"
        sumCatdf.to_pickle(os.path.join(savehere, savename_main))
        
        # save the output dictionary
        for name, df in sumCatdict.items():
            # create file name using the dictionary key
            savename = f'cationflx_{name}.pkl'
            # save
            df.to_pickle(os.path.join(savehere, savename))
    

    # [4] --- ROCK DISSOLUTION ----------------------------------------------
    if "sld_flx" in calc_list: 
        # feedstock must be in list to loop through
        if not isinstance(feedstock, list):
            feedstock = [feedstock]
        for fs in feedstock:
            rockdf = sld_flx(outdir, runname, fs, var_fn = "flx_sld",
                            dust_from_file = True,
                            molar_mass_dict = molar_mass_dict)
            # save
            savename = f"rockflx_{fs}.pkl"
            rockdf.to_pickle(os.path.join(savehere, savename))





# %% 
# ******************************************************************************
# ----------------------- PROFILE METRIC FUNCTIONS ---------------------------------

# dictionary of postprocess function inputs
postproc_prof_dict = {

    "adsorbed": {
        "var_prefix": "prof_aq(ads)",
        "var_units": "mol/L",
    },
    
    # --- requires special function (due to base saturation calc)
    "adsorbed_percCEC": {
        "var_prefix": "prof_aq(ads%cec)",
        "var_units": "%cec",
    },
    # ----------------------------------------------------------

    "adsorbed_ppm": {
        "var_prefix": "prof_aq(adsppm)",
        "var_units": "ppm",
    },

    "aqueous": {
        "var_prefix": "prof_aq",
        "var_units": "mol/L",
    },

    "aqueous_total": {
        "var_prefix": "prof_aq(tot)",
        "var_units": "mol/L",
    },

    # --- requires special fxn (due to variable units)
    "bulksoil": {
        "var_prefix": "bsd"
    },
    # -------------------------------------------------

    "exchange_total": {
        "var_prefix": "prof_ex(tot)",
        "var_units": "mol/L",
    },
    
    "gas": {
        "var_prefix": "prof_gas",
        "var_units": "atm",
        "calculate_mean": False,
    },
    
    "rate": {
        "var_prefix": "rate",
        "var_units": "mol/m2/yr",
    },

    # --- requires special treatment because it's from the lab run
    "soil_ph": {
        "var_prefix": "prof_aq",
        "var_units": "mol/L",
    },
    # -------------------------------------------------------------

    "solid": {
        "var_prefix": "prof_sld",
        "var_units": "mol/m3",
    },

    "solid_sp_saturation": {
        "var_prefix": "sat_sld",
        "var_units": "X",
    }, 

    "solid_volumePercent": {
        "var_prefix": "prof_sld(v%)",
        "var_units": "%vol",
    },

    "solid_weightPercent": {
        "var_prefix": "prof_sld(wt%)",
        "var_units": "wt%",
    },

    "specific_surface_area": {
        "var_prefix": "ssa",
        "var_units": "m2/g(?)",
    },

    "surface_area": {
        "var_prefix": "sa",
        "var_units": "X",
    },
}

def profile_to_ds_optMean(
    outdir: str,
    runname: str,
    var_prefix: str,
    var_units: str,
    time_indices: np.array=np.arange(1,21,1),
    calculate_mean: bool=True,
    depth_mean_suffix: str="coredep",
    run_subdir: str="prof",
    multi_iter: bool = False,
) -> xr.Dataset:
    '''
    Collects the profile data from a given var_prefix for the defined
    time_indices. Converts from pandas dataframes into xr datasets 
    defined over depth and time. If calculate_mean is true, then 
    the mean of each variable is calculated from depth 0 to depth i. 

    Note, the averaging calculation only works for equally spaced depth 
    grids (which SCEPTER is by default)

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <NAME>_field.
    var_prefix : str
        prefix of the file to read in the 'prof' subdirectory. Includes all characters up to `*-xxx.txt`
        where xxx is a 3 digit numeric from 001 to 020. 
    var_units : str
        units for the variables you're reading in (solid profiles, for example, are often mol/m3)
    time_indices : np.array
        SCEPTER outputs 20 time indices as profiles. This array selects which to include in the 
        output dataset. Default is all 20.
    calculate_mean : bool
        [True | False] true means the mean of all variables is calculated at each depth (as 
        though someone took a sample to that depth and calculated the mean)
    depth_mean_suffix : str
        Only used if calculate_mean == True. This is the variable name suffix to append if 
        we calculate the depth mean
    run_subdir : str
        name of the subdirectory that holds the profile files ("prof" by default)
    multi_iter : bool 
        True if the simulation is multiple runs stitched together, False if it is 
        a single run (this matters for the file naming convention)
    
    Returns 
    -------
    xr.Dataset
        all columns of a given variable defined over time and depth. 
    '''
    # loop through time indices to extract data
    for ts in time_indices:
        # --- handle the case where we're reading a composite of multiple runs
        if multi_iter:
            all_iters_found = False   # whether to keep looking for more iters
            idx = 0      # iteration index
            while not all_iters_found:
                if bool((idx == 0) & (ts == time_indices[0])):
                    outdf = get_data_prof(outdir, runname, var_prefix, ts, run_subdir, multi_iter, multi_iter_idx = idx)
                    if type(outdf) == str: # then no data was returned and we assume we have no more indices to check
                        all_iters_found = True
                else:    
                    tmpdf = get_data_prof(outdir, runname, var_prefix, ts, run_subdir, multi_iter, multi_iter_idx = idx)
                    if type(tmpdf) == str:  # then no data was returned and we assume we have no more indices to check
                        all_iters_found = True
                    else:
                        outdf = pd.concat([outdf, tmpdf])

                idx += 1 
            # outdf is not sorted by time, so fix that
            outdf = outdf.sort_values(by=["time", "depth"]).copy()

        # --- handle the standard case of a single run
        else:
            if ts == time_indices[0]:
                outdf = get_data_prof(outdir, runname, var_prefix, ts, run_subdir)
            else:
                tmpdf = get_data_prof(outdir, runname, var_prefix, ts, run_subdir)
                outdf = pd.concat([outdf, tmpdf])


    # convert to an xarray dataset
    ds = outdf.set_index(['depth', 'time']).to_xarray()
    
    if calculate_mean:
        # get the mean wt% from depth 0 to depth i for every depth
        # calculate the cumulative sum and count over the 'depth' dimension
        cumsum_ds = ds.cumsum(dim='depth')
        cumcount_ds = xr.ones_like(ds).cumsum(dim='depth') # get the cumulative count (cumsum / cumcount = mean; since depth is equally spaced)

        # calculate the average by dividing the cumulative sum by the cumulative count
        average_ds = cumsum_ds / cumcount_ds
        average_ds_renamed = average_ds.rename({var: f"{var}_{depth_mean_suffix}" for var in average_ds.data_vars})

        # add new variables to the original dataset
        ds = xr.merge([ds, average_ds_renamed])

    # add units
    for var in ds.data_vars:
        ds[var].attrs["units"] = var_units  # set the same units for all variables
    
    # add output file type
    ds.attrs["outfile"] = var_prefix

    # return result
    return ds


def get_bsd_prof_units(
    var_name
)-> str:
    '''
    Read in a variable name and extract the units from it. 
    Return just the units as a string. Units are denoted by
    [brackets]. If no brackets, return 'NA'

    Parameters
    ----------
    var_name : str
        name of the data_var
    
    Returns
    -------
    str
        the value within the brackets of var_name (or "NA" if no brackets)
    '''
    # check if the variable name contains units in brackets (e.g., var1[m/yr])
    match = re.search(r'\[(.*?)\]', var_name)
    if match:
        # If units are found, return them
        return match.group(1)
    else:
        # If no units are found, return 'NA'
        return 'NA'


def bsd_profile_to_ds(
    outdir: str,
    runname: str,
    var_prefix: str,
    calculate_mean: bool=True,
    time_indices: np.array=np.arange(1,21,1),
    depth_mean_suffix: str="coredep",
    run_subdir: str="prof",
    multi_iter: bool=False,
) -> xr.Dataset:
    '''
    Collects the bsd profile data and turns it into an xarray defined
    over time and depth. Units that are in the brackets of column names
    are also defined at each attribute using the get_bsd_prof_units function

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <NAME>_field.
    var_prefix : str
        prefix of the file to read in the 'prof' subdirectory. Includes all characters up to `*-xxx.txt`
        where xxx is a 3 digit numeric from 001 to 020. 
    calculate_mean : bool
        [True | False] true means the mean of all variables is calculated at each depth (as 
        though someone took a sample to that depth and calculated the mean)
    time_indices : np.array
        SCEPTER outputs 20 time indices as profiles. This array selects which to include in the 
        output dataset. Default is all 20.
    depth_mean_suffix : str
        Only used if calculate_mean == True. This is the variable name suffix to append if 
        we calculate the depth mean
    run_subdir : str
        name of the subdirectory that holds the profile files ("prof" by default)
    multi_iter : bool 
        True if the simulation is multiple runs stitched together, False if it is 
        a single run (this matters for the file naming convention)

    Returns
    -------
    xr.Dataset
        all bsd variables defined over time and depth. 
    '''

    # loop through time indices to extract data
    for ts in time_indices:
        # --- handle the case where we're reading a composite of multiple runs
        if multi_iter:
            all_iters_found = False   # whether to keep looking for more iters
            idx = 0      # iteration index
            while not all_iters_found:
                if bool((idx == 0) & (ts == time_indices[0])):
                    outdf = get_data_prof(outdir, runname, var_prefix, ts, run_subdir, multi_iter, multi_iter_idx = idx)
                    if type(outdf) == str: # then no data was returned and we assume we have no more indices to check
                        all_iters_found = True
                else:    
                    tmpdf = get_data_prof(outdir, runname, var_prefix, ts, run_subdir, multi_iter, multi_iter_idx = idx)
                    if type(tmpdf) == str:  # then no data was returned and we assume we have no more indices to check
                        all_iters_found = True
                    else:
                        outdf = pd.concat([outdf, tmpdf])

                idx += 1 
            # outdf is not sorted by time, so fix that
            outdf = outdf.sort_values(by=["time", "depth"]).copy()

        # --- handle the standard case of a single run
        else:
            if ts == time_indices[0]:
                outdf = get_data_prof(outdir, runname, var_prefix, ts, run_subdir)
            else:
                tmpdf = get_data_prof(outdir, runname, var_prefix, ts, run_subdir)
                outdf = pd.concat([outdf, tmpdf])

    # convert to an xarray dataset
    ds = outdf.set_index(['depth', 'time']).to_xarray()

    # create dictionary with units for each variable
    units_dict = {var: get_bsd_prof_units(var) for var in ds.data_vars}

    # add the units as an attribute to each variable in the dataset
    for var in ds.data_vars:
        ds[var].attrs['units'] = units_dict[var]
    
    # add depth mean base saturation
    if calculate_mean:
        # get the mean wt% from depth 0 to depth i for every depth
        # calculate the cumulative sum and count over the 'depth' dimension
        cumsum_ds = ds.cumsum(dim='depth')
        cumcount_ds = xr.ones_like(ds).cumsum(dim='depth') # get the cumulative count (cumsum / cumcount = mean; since depth is equally spaced)

        # calculate the average by dividing the cumulative sum by the cumulative count
        average_ds = cumsum_ds / cumcount_ds
        average_ds_renamed = average_ds.rename({var: f"{var}_{depth_mean_suffix}" for var in average_ds.data_vars})

        # add new variables to the original dataset
        ds = xr.merge([ds, average_ds_renamed])

    # return result
    return ds


def ads_percCec_prof_baseSat(
    outdir: str,
    runname: str,
    var_prefix: str,
    var_units: str,
    time_indices: np.array=np.arange(1,21,1),
    calculate_mean: bool=True,
    depth_mean_suffix: str="coredep",
    run_subdir: str="prof",
    multi_iter: bool=False,
) -> xr.Dataset:
    '''
    Get adsorbed species as a percent of CEC in an xarray dataset. This is just a wrapper 
    around the profile_to_ds_optMean function which computes base saturation after getting 
    the ads%cec profile as a dataset.
    
    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <NAME>_field.
    var_prefix : str
        prefix of the file to read in the 'prof' subdirectory. Includes all characters up to `*-xxx.txt`
        where xxx is a 3 digit numeric from 001 to 020. 
    var_units : str
        name of the units for the variables to add as an attribute to the dataset
    time_indices : np.array
        SCEPTER outputs 20 time indices as profiles. This array selects which to include in the 
        output dataset. Default is all 20.
    calculate_mean : bool
        [True | False] true means the mean of all variables is calculated at each depth (as 
        though someone took a sample to that depth and calculated the mean)
    depth_mean_suffix : str
        Only used if calculate_mean == True. This is the variable name suffix to append if 
        we calculate the depth mean
    run_subdir : str
        name of the subdirectory that holds the profile files ("prof" by default)
    multi_iter : bool 
        True if the simulation is multiple runs stitched together, False if it is 
        a single run (this matters for the file naming convention)

    Returns
    -------
    xr.Dataset
        all bsd variables defined over time and depth. 
    '''
    # get the profile as a dataset
    ds = profile_to_ds_optMean(outdir, runname, var_prefix, var_units, time_indices, calculate_mean, depth_mean_suffix, run_subdir, multi_iter = multi_iter)

    # compute base saturation
    bs_sp_list = ['ca','mg','k','na']     # bases that we want to sum
    valid_vars = [var for var in bs_sp_list if var in ds]

    if valid_vars:
        ds["base_saturation"] = ds[valid_vars].to_array(dim="vars").sum(dim="vars")
        ds["base_saturation"].attrs["units"] = var_units  # add units
        # add depth mean base saturation
        if calculate_mean:
            # get the mean wt% from depth 0 to depth i for every depth
            # calculate the cumulative sum and count over the 'depth' dimension
            cumsum_ds = ds.cumsum(dim='depth')
            cumcount_ds = xr.ones_like(ds).cumsum(dim='depth') # get the cumulative count (cumsum / cumcount = mean; since depth is equally spaced)

            # calculate the average by dividing the cumulative sum by the cumulative count
            average_ds = cumsum_ds / cumcount_ds
            average_ds_renamed = average_ds.rename({var: f"{var}_{depth_mean_suffix}" for var in average_ds.data_vars})

            # add new variables to the original dataset
            ds = xr.merge([ds, average_ds_renamed[f'base_saturation_{depth_mean_suffix}']])
    else:
        ds['base_saturation'] = np.nan
        ds["base_saturation"].attrs["units"] = var_units  # add units
    
    # return result
    return ds


def prof_postproc_save(
    outdir: str, 
    runname_field: str, 
    runname_lab: str,
    postproc_prof_list: list=["all"],
    save_dir: str="postproc_profs",
    multi_iter: bool=False,
):
    '''
    Convert SCEPTER/scepter_output/myrun/prof/* files to .nc files. Relies on other profile
    postproc embedded functions. Only the profile files listed in postproc_prof_list will 
    be processed. 

    Parameters
    ----------
    outdir : str
        path to the SCEPTER output directory. Usually something like 'my/path/SCEPTER/scepter_output'
    runname_field : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <NAME>_field.
    runname_lab : str
        name of the SCEPTER run (equivalent to the directory within outdir). Generally <NAME>_lab.
    postproc_prof_list : list 
        list of postprocess files that are also keywords in postproc_prof_dict (see `cflx_proc.py`)
    save_dir : str
        name of the subdirectory where the .nc files are saved
    multi_iter : bool 
        True if the simulation is multiple runs stitched together, False if it is 
        a single run (this matters for the file naming convention)

    Returns
    -------

    '''
    # where to save the results
    savehere = os.path.join(outdir, runname_field, save_dir)
    # make dir if it doesn't exist
    if not os.path.exists(savehere):
        os.makedirs(savehere)

    # create the list if user asks for all postproc files
    if (postproc_prof_list == "all") or (postproc_prof_list == ["all"]):
        postproc_prof_list = list(postproc_prof_dict.keys())

    # check that all listed prof names are in the dictionary
    missing_keys = [key for key in postproc_prof_list if key not in postproc_prof_dict]
    if missing_keys:
        print(f"Warning: The following postprocess profile names are not compatible: {missing_keys}; check for typos!")
        # remove the missing ones
        pplist_new = [key for key in postproc_prof_list if key in postproc_prof_dict]
        postproc_prof_list = pplist_new

    # --- loop through postproc list, save result
    for pp in postproc_prof_list:
        
        # --- check for special cases first
        if pp == "bulksoil":
            ds = bsd_profile_to_ds(outdir, runname_field, **postproc_prof_dict[pp], multi_iter=multi_iter)
        elif pp == "adsorbed_percCEC":
            ds = ads_percCec_prof_baseSat(outdir, runname_field, **postproc_prof_dict[pp], multi_iter=multi_iter)
        elif pp == "soil_ph":
            ds = profile_to_ds_optMean(outdir, runname_lab, **postproc_prof_dict[pp], multi_iter=multi_iter)
            ds = ds.sel(time = np.max(ds.time.values)) # keep only the max time (lab run time is not aligned with field!)
        # ---------------------------------
        else:
            ds = profile_to_ds_optMean(outdir, runname_field, **postproc_prof_dict[pp], multi_iter=multi_iter)
        
        # --- save 
        savename = f"{pp}.nc"
        ds.to_netcdf(os.path.join(savehere, savename))
        



# ------- SCRATCH ----------------------------------------------------- # 
# %%
# outdir = "/home/tykukla/SCEPTER/scepter_output"
# # runname = "noFert_cc_multiyear_site_311b_app_10000p0_psize_75_composite_field"
# runname = "midfert_cc_311a_0p0_cc_field_tau15p0"
# feedstock = ['cc', 'amnt']

# cflx_calc(outdir, runname, feedstock)







# %% 

# outdir = "/home/tykukla/SCEPTER/scepter_output"
# runname = "noFert_cc_multiyear_site_311b_app_10000p0_psize_75_composite_field"
# # runname = "lowfert_cc_311b_30p0_cc_field_tau15p0"
# var_fn = 'flx_co2sp'
# cdvar = 'ALK'

# df, dfint = get_data(outdir, runname, var_fn, cdvar)

# df = co2_flx(outdir, runname)
# df
# %% 
# import matplotlib.pyplot as plt

# plt.plot(dfint['time'], dfint['adv']*dfint['time'], label='adv')
# plt.plot(dfint['time'], dfint['tflx']*dfint['time'], label='tflx')
# plt.legend()
# %% 
# plt.plot(df['time'], df['adv'], label='adv')
# plt.plot(df['time'], df['tflx'], label='tflx')
# plt.plot(df['time'], df['tflx'] + df['adv'], label='tflx+adv')
# plt.legend()

# %% 

# outdir = "/home/tykukla/SCEPTER/scepter_output"
# runname = "noFert_cc_multiyear_site_311b_app_10000p0_psize_75_composite_field"
# var_fn = "flx_sld"
# cdvar = "cc"
# organic_list = ["g1", "g2", "g3"]


# # get the txt file as a pandas dataframe
# fn = f"{var_fn}-{cdvar}.txt"   # flux timeseries
# fn_int = f"int_{fn}"       # integrated flux timeseries
# df = preprocess_txt(outdir, runname, fn, run_subdir = "flx")
# dfint = preprocess_txt(outdir, runname, fn_int, run_subdir = "flx")

# # %% 
# plt.plot(dfint['time'], dfint['tflx']*dfint['time'])
# plt.plot(dfint['time'], dfint['rain']*dfint['time'])

# # %% 
# plt.plot(dfint['time'], dfint['time']*(-1*(dfint['rain']) - dfint['tflx']))
# plt.plot(dfint['time'], dfint['cc']*dfint['time'])

# # %% 
# plt.plot(dfint['time'], dfint['time']*(-1*(dfint['rain']) - dfint['cc']))
# plt.plot(dfint['time'], (dfint['tflx']+dfint['adv'])*dfint['time'])

# # %% 
# plt.plot(dfint['time'], dfint['cc']*dfint['time'], label='diss')
# plt.plot(dfint['time'], -1*dfint['rain']*dfint['time'], label='rain')
# plt.plot(dfint['time'], dfint['tflx']*dfint['time'], label='tflx')
# plt.plot(dfint['time'], dfint['adv']*dfint['time'], label='adv')
# plt.legend()

# # %% 
# plt.plot(dfint['time'], dfint['time']*(-1*dfint['rain']-dfint['adv']), label='rain-adv')
# plt.plot(dfint['time'], dfint['cc']*dfint['time'], label='diss')
# plt.plot(dfint['time'], dfint['tflx']*dfint['time'], label='tflx')
# plt.legend()

# # %% 
# plt.plot(df['time'], df['cc'], label='diss')
# plt.plot(df['time'], -1*df['rain'], label='rain')
# plt.plot(df['time'], df['tflx'], label='tflx')
# plt.plot(df['time'], df['adv'], label='adv')
# plt.ylim(-10, 10)
# plt.legend()



# # %%
# plt.plot(dfint['time'], dfint['cc']*dfint['time'])
# plt.plot(dfint['time'], -1*dfint['rain']*dfint['time'])

# # %%
# outdir = "/home/tykukla/SCEPTER/scepter_output"
# runname = "noFert_cc_multiyear_site_311b_app_10000p0_psize_75_composite_field"
# fn = "flx_aq-ca.txt"
# df = preprocess_txt(outdir, runname, fn)
# # %%
# df
# # %%
# import matplotlib.pyplot as plt
# plt.plot(df['time'], df['adv'])
# plt.plot(df['time'], df['tflx'])
# # %%



# %% 

# # .. inputs minus outputs
# inputs_x = (dfint['cc'] + dfint['g2'])*dfint['time']
# outputs_x = (dfint['dif'] + dfint['adv'])*dfint['time']


## -- CONFIRMED
## 
## cc + g2 - (dif + adv) = tflx
## --


# plt.plot(dfint['time'],  -1*dfint['cc']*dfint['time'], label='cc_in')
# plt.plot(dfint['time'],  -1*dfint['g2']*dfint['time'], label='g2_in')
# plt.plot(dfint['time'],  -1*dfint['dif']*dfint['time'], label='dif_out')
# plt.plot(dfint['time'],  -1*dfint['adv']*dfint['time'], label='adv_out')
# plt.legend()
# # %%


# plt.plot(dfint['time'], -1*dfint['adv']*dfint['time'], label='adv_out')
# plt.plot(dfint['time'], -1*dfint['cc']*dfint['time'], label='cc_in')
# plt.plot(dfint['time'], -1*dfint['g2']*dfint['time'], label='g2_in')
# plt.plot(dfint['time'], -1*dfint['tflx']*dfint['time'], label='tflx_')
# plt.plot(dfint['time'], -1*(dfint['adv'] + dfint['cc'])*dfint['time'], label='adv - cc')
# plt.plot(dfint['time'], -1*(dfint['adv'] + (dfint['cc']-dfint['tflx']))*dfint['time'], label='adv - (cc-tflx)')
# plt.legend()
# # %%
# plt.plot(dfint['time'], dfint['tflx']*dfint['time'], label='tflx_')
# plt.plot(dfint['time'], dfint['adv']*dfint['time'], label='adv_out')
# # plt.plot(dfint['time'], dfint['dif']*dfint['time'], label='dif_out')
# plt.plot(dfint['time'], dfint['cc']*dfint['time'], label='cc_in')
# plt.legend()


# # %% 
# plt.plot(dfint['time'], (dfint['cc']*dfint['time'] / ((dfint['g2'] + dfint['cc'])*dfint['time'])))
# # %%

# plt.plot(dfint['time'], (dfint['cc'])*dfint['time'], label='cc')
# plt.plot(dfint['time'], (dfint['g2'])*dfint['time'], label='g2')
# plt.legend()
# # %%
# plt.plot(dfint['time'], (dfint['cc'])*dfint['time'] / ((dfint['cc'])*dfint['time'] + (dfint['g2'])*dfint['time']))


# # %%
# plt.plot(df['time'], df['adv'] - df['adv'][0], label='adv_out')
# plt.plot(df['time'], df['cc'] - df['cc'][0], label='cc_in')
# plt.plot(df['time'], df['g2'] - df['g2'][0], label='g2_in')
# plt.plot(df['time'], df['tflx'] - df['tflx'][0], label='tflx_')
# # plt.plot(df['time'], (df['adv'] + dfint['cc']), label='adv - cc')
# # plt.plot(df['time'], (df['adv'] + (dfint['cc']-dfint['tflx'])), label='adv - (cc-tflx)')
# plt.ylim(-4,4)
# plt.legend()

# # %%
# plt.plot(df['time'], df['adv'] - df['adv'][0], label='adv_out')
# plt.plot(df['time'], df['cc'] , label='cc_in')
# # plt.plot(df['time'], df['g2'] , label='g2_in')
# plt.plot(df['time'], df['g2'] - df['g2'][0], label='g2_in')
# plt.plot(df['time'], df['tflx'] , label='tflx_')
# # plt.plot(df['time'], (df['adv'] + dfint['cc']), label='adv - cc')
# # plt.plot(df['time'], (df['adv'] + (dfint['cc']-dfint['tflx'])), label='adv - (cc-tflx)')
# plt.ylim(-4,4)
# plt.legend()
# # %%

# plt.plot(df['time'], df['adv'] - df['adv'][0], label='adv_out')
# plt.plot(df['time'], -1*df['cc'] , label='cc_in')
# # plt.plot(df['time'], df['g2'] , label='g2_in')
# plt.plot(df['time'], -1*df['g2'] - -1*df['g2'][0], label='g2_in')
# plt.plot(df['time'], -1*df['dif'] - -24.2, label='dif_in')
# # plt.plot(df['time'], df['tflx'] , label='tflx_')
# # plt.plot(df['time'], (df['adv'] + dfint['cc']), label='adv - cc')
# # plt.plot(df['time'], (df['adv'] + (dfint['cc']-dfint['tflx'])), label='adv - (cc-tflx)')
# plt.ylim(-4,4)
# plt.legend()

# # %% 
# plt.plot(dfint['time'], -1*dfint['dif']*dfint['time'] - dfint['g2']*dfint['time'], label='-dif-resp')
# plt.plot(dfint['time'], dfint['adv']*dfint['time'], label='adv')
# plt.plot(dfint['time'], -1*dfint['cc']*dfint['time'], label='cc')
# plt.plot(dfint['time'], -1*dfint['dif']*dfint['time'] - dfint['g2']*dfint['time'] - dfint['cc']*dfint['time'], label='-dif-resp+cc')
# plt.plot(dfint['time'], -1*dfint['dif']*dfint['time'] - dfint['g2']*dfint['time'] - dfint['cc']*dfint['time'] - dfint['tflx']*dfint['time'], label='-dif-resp+cc+tflx',linestyle='dashed')
# #plt.plot(dfint['time'], -1*dfint['dif']*dfint['time'], label='-dif')
# #plt.plot(dfint['time'], dfint['g2']*dfint['time'], label='resp')
# plt.legend()

# # --- CONFIRMED -- 
# # adv = (-dif-resp) - cc - tflx
# # where:
# #    -dif-resp = respired carbon that's stored (it's positive when resp > diffusion, indicating respired carbon staying in the system)
# #    cc = c flux due to cc (cc is negative directionally, so subtracting it is the same as adding its absolute value (e.g., cc added to the system))
# #    tflx = change in storage (tflx is negative directionally, so subtracting it is the same as adding its absolute value (e.g., added storage to the system))


# # %%

# # plt.plot(df['time'], -1*df['dif'] - df['g2'], label='-dif-resp')
# plt.plot(df['time'], df['adv'], label='adv')
# # plt.plot(df['time'], -1*df['cc'], label='cc')
# # plt.plot(df['time'], -1*df['dif'] - df['g2'] - df['cc'], label='-dif-resp+cc')
# plt.plot(df['time'], -1*df['dif'] - df['g2'] - df['cc'] - df['tflx'], label='-dif-resp-cc-tflx',linestyle='dashed')
# #plt.plot(dfint['time'], -1*dfint['dif']*dfint['time'], label='-dif')
# #plt.plot(dfint['time'], dfint['g2']*dfint['time'], label='resp')
# plt.legend()

# # %%
# plt.plot(df['time'], df['tflx'], label='tflx')
# plt.plot(df['time'], -1*df['dif'] - df['g2'] - df['cc'] - df['adv'], label='-dif-resp-cc-adv',linestyle='dashed')
# plt.legend()
# # %%
