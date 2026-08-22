# --------------------------------------------
#
# helper functions for building the batch
# input .csv files
#
# --------------------------------------------
import io
import itertools
import os
import warnings

import fsspec
import numpy as np
import pandas as pd
import xarray as xr


# %% 
## [1] BUILD DATAFRAME
def build_df(pref: str, const_dict: dict, sites: list, 
            by_site: dict, all_combinations: dict, add_ctrl: bool) -> pd.DataFrame:
    """
    Turn dictionaries into a .csv file that serves as a batch input for 
    running SCEPTER. Vars are either constant for all runs, only vary by 
    site, or vary from one run to the next (within a site). 
    All dicts should be structured such that the key
    is the column name, the value is a cell value

    Parameters
    ----------
    pref : str
        prefix for each run name
    const_dict : dict
        dictionary for the values that are held constant 
        across all simulations
    sites : list
        list of the site names across which to run the simulations
    by_site : dict
        dictionary for all the values that only vary by site
        (each value should be a list of len(sites) where the 
        first value corresponds to the first site indexed, and 
        so on)
    all_combinations : dict
        dictionary for all the values to vary such that every
        unique combination of these values is tested for each 
        site
    add_ctrl : bool
        [True | False] whether to add a control simulation with zero 
        dust application for each site (NOTE, you may not want to add 
        0 to your dust app rate list because it will needlessly repeat
        for every other var in all_combinations)

    Returns
    -------
    pd.DataFrame
        This is the file that will become the .csv batch input. Each column
        should be a variable that the SCEPTER python scripts can recognize 
        (no typos!)
    """
    # [1] generate all combinations from the dict's vectors
    all_combos_list = list(itertools.product(*all_combinations.values()))
    # make dataframe and repeat values for the number of sites
    df = pd.DataFrame(all_combos_list * len(sites), columns=all_combinations.keys())

    # [2] add site-specific vars
    # add site labels
    df['site'] = [site for site in sites for _ in range(len(all_combos_list))]
    # for each key in by_site, alternate the values for the corresponding site
    for key, values in by_site.items():
        df[key] = [values[site_idx] for site_idx in range(len(sites)) for _ in range(len(all_combos_list))]

    # [3] add constant vars
    for key, value in const_dict.items():
        df[key] = value
    
    # [4] add control cases (no dust application) if add_ctrl is True
    if add_ctrl:
        # loop through sites
        for thissite in reversed(sites):
            tmp_row = df[df['site'] == thissite].iloc[0]
            # set dust to zero
            tmp_row['dustrate'] = 0.0
            # concat to top row
            df = pd.concat([pd.DataFrame([tmp_row]), df], ignore_index=True)

    # [5] add the newrun ID
    df = newrun_id_fxn(df, pref, None)
    # check that all the run ids are unique and return a warning if not
    if not df['newrun_id'].is_unique:
        warnings.warn("Column newrun_id contains duplicate entries. The latter run ID will likely overwrite the former")

    return df


# [1a] DATAFRAME FOR SENSITIVITY ANALYSIS 
def build_df_sa(
    pref: str,
    const_dict: dict,
    sites: list,
    by_site: dict,
    param_values: list,
    problem: dict,
    param_values_ctrl: dict,
    add_ctrl: bool=True,
) -> pd.DataFrame: 
    '''
    Turn SA problem (array of input parameters) into a batch dataframe.

    Parameters
    ----------
    pref : str
        prefix for each run name
    const_dict : dict
        dictionary for the values that are held constant 
        across all simulations
    sites : list
        list of the site names across which to run the simulations
    by_site : dict
        dictionary for all the values that only vary by site
        (each value should be a list of len(sites) where the 
        first value corresponds to the first site indexed, and 
        so on)
    param_values : list
        array of input parameters produced by the SALib sampler
    problem : dict
        dictionary of parameter inputs for the SALib sampler
    param_values_ctrl : dict
        dictionary of control parameter inputs for the control run 
        (note, control run is placed at the zero index in the df)
    add_ctrl : bool
        [True | False] whether to add a control simulation with zero 
        dust application for each site (NOTE, you may not want to add 
        0 to your dust app rate list because it will needlessly repeat
        for every other var in all_combinations)
    '''
    # [1] generate the param_values df
    df_init = pd.DataFrame(param_values, columns = problem['names'])
    # repeat values for the number of sites
    df = pd.concat(
        [df_init.assign(site=site) for site in sites],
        ignore_index=True
    )

    # [2] add site-specific vars
    # for each key in by_site, alternate the values for the corresponding site
    for key, values in by_site.items():
        df[key] = [values[site_idx] for site_idx in range(len(sites)) for _ in range(len(df))]

    # [3] add constant vars
    for key, value in const_dict.items():
        df[key] = value

    # [4] add control cases (no dust application) if add_ctrl is True
    if add_ctrl:
        # loop through sites
        for thissite in reversed(sites):
            tmp_row = df[df['site'] == thissite].iloc[0]
            # set dust to zero
            tmp_row['dustrate'] = 0.0
            # set other vars
            for key, value in param_values_ctrl.items():
                tmp_row[key] = value
            # concat to top row
            df = pd.concat([pd.DataFrame([tmp_row]), df], ignore_index=True)

    # [5] add the newrun ID
    df = newrun_id_fxn(df, pref, None, add_index = True)
    # check that all the run ids are unique and return a warning if not
    if not df['newrun_id'].is_unique:
        warnings.warn("Column newrun_id contains duplicate entries. The latter run ID will likely overwrite the former")

    # return result
    return df

# [1b] BUILD DF FOR ALL_COMBINATIONS (+allowing a flexible control)
# (key difference from `build_df` is that "const_dict" is used for 
#  control values, and does not overwrite the other dict)
def build_df_all_combinations(
        pref: str,
        const_dict: dict,
        sites: list,
        by_site: dict,
        all_combinations: dict,
        add_ctrl: bool,
        add_index: bool=True,
) -> pd.DataFrame:
    '''
    Turn dictionaries into a .csv file that serves as a batch input for 
    running SCEPTER. Vars are either constant for all runs, only vary by 
    site, or vary from one run to the next (within a site). 
    All dicts should be structured such that the key
    is the column name, the value is a cell value

    Parameters
    ----------
    pref : str
        prefix for each run name
    const_dict : dict
        dictionary for the values that are held constant 
        across all simulations
    sites : list
        list of the site names across which to run the simulations
    by_site : dict
        dictionary for all the values that only vary by site
        (each value should be a list of len(sites) where the 
        first value corresponds to the first site indexed, and 
        so on)
    all_combinations : dict
        dictionary for all the values to vary such that every
        unique combination of these values is tested for each 
        site
    add_ctrl : bool
        [True | False] whether to add a control simulation with zero 
        dust application for each site (NOTE, you may not want to add 
        0 to your dust app rate list because it will needlessly repeat
        for every other var in all_combinations)

    Returns
    -------
    pd.DataFrame
        This is the file that will become the .csv batch input. Each column
        should be a variable that the SCEPTER python scripts can recognize 
        (no typos!)
    '''
    # [1] generate all combinations from the dict's vectors
    all_combos_list = list(itertools.product(*all_combinations.values()))
    # make dataframe and repeat values for the number of sites
    df = pd.DataFrame(all_combos_list * len(sites), columns=all_combinations.keys())

    # [2] add site-specific vars
    # add site labels
    df['site'] = [site for site in sites for _ in range(len(all_combos_list))]
    # for each key in by_site, alternate the values for the corresponding site
    for key, values in by_site.items():
        df[key] = [values[site_idx] for site_idx in range(len(sites)) for _ in range(len(all_combos_list))]

    # [3] add constant vars that do not appear in "all_combinations"
    unique_keys = {k: v for k, v in const_dict.items() if k not in df.columns}
    control_keys = {k: v for k, v in const_dict.items() if k in df.columns}
    for key, value in unique_keys.items():
        df[key] = value
    
    # [4] add control cases (no dust application) if add_ctrl is True
    if add_ctrl:
        # loop through sites
        for thissite in reversed(sites):
            tmp_row = df[df['site'] == thissite].iloc[0]
            # set dust to zero
            tmp_row['dustrate'] = 0.0
            # add control keys 
            for key, value in control_keys.items():
                tmp_row[key] = value
            # concat to top row
            df = pd.concat([pd.DataFrame([tmp_row]), df], ignore_index=True)

    # [5] add the newrun ID
    df = newrun_id_fxn(df, pref, None, add_index=add_index)
    # check that all the run ids are unique and return a warning if not
    if not df['newrun_id'].is_unique:
        warnings.warn("Column newrun_id contains duplicate entries. The latter run ID will likely overwrite the former")

    return df


def build_df_one_at_a_time(
        pref: str,
        const_dict: dict,
        sites: list,
        by_site: dict,
        all_combinations: dict,
        individual_cases: dict,
        add_ctrl: bool=True,
        add_index: bool = True,
)-> pd.DataFrame:
    """
    build dataframe where we perturb one param at a time
    from the "individual_cases". Each individual case is 
    repeated for each site and each "all_combinations" 
    combo. 

    Parameters
    ----------
    pref : str
        prefix for each run name
    const_dict : dict
        dictionary for the values that are held constant 
        across all simulations
    sites : list
        list of the site names across which to run the simulations
    by_site : dict
        dictionary for all the values that only vary by site
        (each value should be a list of len(sites) where the 
        first value corresponds to the first site indexed, and 
        so on)
    all_combinations : dict
        dictionary for all the values to vary such that every
        unique combination of these values is tested for each 
        site
    individual_cases : dict
        dictionary for each one-at-a-time perturbation
    add_ctrl : bool
        [True | False] whether to add a control simulation with zero 
        dust application for each site (NOTE, you may not want to add 
        0 to your dust app rate list because it will needlessly repeat
        for every other var in all_combinations)

    Returns
    -------
    pd.DataFrame
        This is the file that will become the .csv batch input. Each column
        should be a variable that the SCEPTER python scripts can recognize 
        (no typos!)
    """
    # [1] generate all combinations from the dict's vectors
    all_combos_list = list(itertools.product(*all_combinations.values()))
    # make dataframe and repeat values for the number of sites
    df = pd.DataFrame(all_combos_list * len(sites), columns=all_combinations.keys())

    # [2] add site-specific vars
    # add site labels
    df['site'] = [site for site in sites for _ in range(len(all_combos_list))]
    # for each key in by_site, alternate the values for the corresponding site
    for key, values in by_site.items():
        df[key] = [values[site_idx] for site_idx in range(len(sites)) for _ in range(len(all_combos_list))]

    # [3] add constant vars that do not appear in "all_combinations"
    unique_keys = {k: v for k, v in const_dict.items() if k not in df.columns}
    control_keys = {k: v for k, v in const_dict.items() if k in df.columns}
    for key, value in unique_keys.items():
        df[key] = value

    # [4] repeat for all individual cases
    # each key in `individual_cases` creates a *set* of additional simulations
    individual_override_list = []
    for key, values in individual_cases.items():
        for v in values:
            individual_override_list.append({key: v})
    dfs = [df]
    for override in individual_override_list:
        df_new = df.copy()
        # apply override to all rows
        for k, v in override.items():
            df_new[k] = v
        dfs.append(df_new)
    df = pd.concat(dfs, ignore_index=True)

    # [5] add control cases (no dust application) if add_ctrl is True
    if add_ctrl:
        # loop through sites
        for thissite in reversed(sites):
            tmp_row = df[df['site'] == thissite].iloc[0]
            # set dust to zero
            tmp_row['dustrate'] = 0.0
            # add control keys 
            for key, value in control_keys.items():
                tmp_row[key] = value
            # concat to top row
            df = pd.concat([pd.DataFrame([tmp_row]), df], ignore_index=True)

    # [6] add the newrun ID
    df = newrun_id_fxn(df, pref, None, add_index=add_index)
    # check that all the run ids are unique and return a warning if not
    if not df['newrun_id'].is_unique:
        print("YIKES!!!")
        # warnings.warn("Column newrun_id contains duplicate entries. The latter run ID will likely overwrite the former")

    # return result
    return df



def build_df_multi_spinup(
        pref: str,
        const_dict: dict,
        sites: list,
        by_site: dict,
        all_combinations: dict,
        individual_cases: dict,
        add_ctrl: bool=True,
        add_index: bool = True,
)-> pd.DataFrame:
    """
    build dataframe where we have multiple "by_site" entries for each
    site. These are repeated across "all_combinations", which generally
    includes zero dust flux to capture the control. 

    Parameters
    ----------
    pref : str
        prefix for each run name
    const_dict : dict
        dictionary for the values that are held constant 
        across all simulations
    sites : list
        list of the site names across which to run the simulations
    by_site : dict
        dictionary for all the values that only vary by site
        (each value should be a list of len(sites) where the 
        first value corresponds to the first site indexed, and 
        so on)
    all_combinations : dict
        dictionary for all the values to vary such that every
        unique combination of these values is tested for each 
        site
    individual_cases : dict
        dictionary for each one-at-a-time perturbation
    add_ctrl : bool
        [True | False] whether to add a control simulation with zero 
        dust application for each site (NOTE, you may not want to add 
        0 to your dust app rate list because it will needlessly repeat
        for every other var in all_combinations)

    Returns
    -------
    pd.DataFrame
        This is the file that will become the .csv batch input. Each column
        should be a variable that the SCEPTER python scripts can recognize 
        (no typos!)
    """
        
    # [1] generate all combinations from the dict's vectors
    all_combos_list = list(itertools.product(*all_combinations.values()))
    # make dataframe and repeat values for the number of sites
    df = pd.DataFrame(all_combos_list, columns=all_combinations.keys())

    # [2] add site-specific vars
    df_site = pd.DataFrame(by_site)
    # merge with all_combinations
    df = df_site.merge(df, how='cross')

    # [3] add constant vars that do not appear in "all_combinations"
    unique_keys = {k: v for k, v in const_dict.items() if k not in df.columns}
    control_keys = {k: v for k, v in const_dict.items() if k in df.columns}
    for key, value in unique_keys.items():
        df[key] = value

    # [4] repeat for all individual cases
    # each key in `individual_cases` creates a *set* of additional simulations
    individual_override_list = []
    for key, values in individual_cases.items():
        for v in values:
            individual_override_list.append({key: v})
    dfs = [df]
    for override in individual_override_list:
        df_new = df.copy()
        # apply override to all rows
        for k, v in override.items():
            df_new[k] = v
        dfs.append(df_new)
    df = pd.concat(dfs, ignore_index=True)

    # [5] add control cases (no dust application) if add_ctrl is True
    if add_ctrl:
        # loop through sites
        for thissite in reversed(sites):
            tmp_row = df[df['site'] == thissite].iloc[0]
            # set dust to zero
            tmp_row['dustrate'] = 0.0
            # add control keys 
            for key, value in control_keys.items():
                tmp_row[key] = value
            # concat to top row
            df = pd.concat([pd.DataFrame([tmp_row]), df], ignore_index=True)

    # [6] add the newrun ID
    df = newrun_id_fxn(df, pref, None, add_index=add_index)
    # check that all the run ids are unique and return a warning if not
    if not df['newrun_id'].is_unique:
        print("YIKES!!!")
        # warnings.warn("Column newrun_id 
    # return result
    return df



# [2] MAKE NEWRUN ID
def newrun_id_fxn(df: pd.DataFrame, pref: str, clim_tag: str, add_index: bool=False) -> pd.DataFrame:
    """
    Generate an ID string for a given run based on the inputs (especially
    the application rate, duration, and dust radius -- put the dust type in 
    the prefix, though note it can get added in SCEPTER py script anyway.)

    Parameters
    ----------
    df : pd.DataFrame
        dataframe generated within build_df function which has all inputs
        represented for each row.
    pref : str
        prefix for the run (this will often include the dust type; e.g., 'gbas')
        it can be None, but that's not recommended
    clim_tag : str
        the tag for the years of climate model output used (e.g., '1950-2020'). 
        If no climate model data is used, set this to None
    add_index : bool
        True to add the row index to the runname (ensures no duplicate names, esp for
        sensitivity analysis applications). Defaults to false 
    
    Returns
    -------
    pd.DataFrame
        This is the file that will become the .csv batch input. Each column
        should be a variable that the SCEPTER python scripts can recognize 
        (no typos!)
    """
    # create empty column
    df['newrun_id'] = None
    # loop through each row
    for index, row in df.iterrows():
        # assign the dust flux if it exists
        dstflx = "app_" + str(row['dustrate']).replace('.', 'p') if 'dustrate' in row else None
        # assign the particle size if it exists
        psize = "psize_" + str(row['dustrad']).replace('.', 'p') if 'dustrad' in row else None

        # pull out other vars
        try:
            site = row['climatefiles']
        except:
            site = row['site_TMP']
        dur = "tau_" + str(row['duration']).replace('.', 'p') if 'duration' in row else None
        # (not using dur for now because the python script in SCEPTER repo adds it in itself)

        # create the new ID
        if add_index: 
            this_id = '_'.join([s for s in [pref, site, clim_tag, dstflx, psize, str(index)] if s is not None])
        else:
            this_id = '_'.join([s for s in [pref, site, clim_tag, dstflx, psize] if s is not None])
        
        # add it to the pandas df
        df.at[index, 'newrun_id'] = this_id
    
    # return result
    return df

# [3] SAVE DATAFRAME AS CSV
def save_df(df: pd.DataFrame, savepath_batch: str, fn: str,
            multi_run_split: bool, max_iters_per_set: int=None):
    """
    Save the pandas DataFrame as a .csv file to use for batch inputs
    for SCEPTER. If you're filename already exists, the function will
    warn you and ask if you want to proceed with overwriting it. 

    Parameters
    ----------
    df : pd.DataFrame
        dataframe generated within build_df function which has all inputs
        represented for each row.
    savepath_batch : str
        path to the directory where you want to save the output csv files
    fn : str
        filename of the .csv file we want to save (if we split into multiple,
        this is the core name that holds it all together)
    multi_run_split : bool
        T/F; whether or not to split the df into multiple CSV files
    max_iters_per_set : int
        if multi_run_split is true, then how many iters should we allow per 
        csv file (default is None, assuming multi_run_split is False)
    
    Returns
    -------

    """
    if multi_run_split:
        # calculate the number of chunks needed
        num_chunks = -(-len(df) // max_iters_per_set)  # round up the division result
        # splitting df into chunks of approximately 20 rows each
        dfs = [
            df.iloc[i * max_iters_per_set : (i + 1) * max_iters_per_set]
            for i in range(num_chunks)
        ]
        # saving each df to a different file
        for i, df_chunk in enumerate(dfs):
            tot_dfs = len(dfs)
            savefn_nosuff = fn.rstrip(".csv")
            savefn_out = f"{savefn_nosuff}_set{i+1}of{tot_dfs}.csv"
            # print(os.path.join(savepath_batch, savefn_out))
            savepath = os.path.join(savepath_batch, savefn_out)
            # save and prevent accidentally overwriting another file
            prevent_accidental_overwrite(df_chunk, savepath, savefn_out)
        # also save total
        savepath_all = os.path.join(savepath_batch, fn)
        prevent_accidental_overwrite(df, savepath_all, fn)
    else:
        savepath_all = os.path.join(savepath_batch, fn)
        prevent_accidental_overwrite(df, savepath_all, fn)
        

# FUNCTION check if a file already exists
def prevent_accidental_overwrite(df: pd.DataFrame, path_fn: str, fn: str):
    """
    This function is used within the save_df function to check that we're 
    not going to accidentally overwrite another file upon save. If we are,
    we ask the user to give permission to proceed. 

    Parameters
    ----------
    df : pd.DataFrame
        dataframe generated within build_df function which has all inputs
        represented for each row.
    path_fn : str
        full path with filename of the file we want to save (but haven't yet)
        (ex: "/this/is/my/path/to/file.csv")
    fn : str
        just the name of the file
        (ex: "file.csv")

    Returns 
    -------

    """
    # check if the file already exists at this path
    if os.path.exists(path_fn):
        # ask the user whether they want to continue / overwrite the file
        response = input(f"A file with this name already exists ({fn}). Do you want to overwrite it? (Y/N): ").strip().upper()
        if response == "N":
            warnings.warn("File save canceled to prevent overwrite")
            return
        elif response != "Y":
            warnings.warn("Invalid response. File save canceled.")
            return
        else:
            df.to_csv(path_fn, index=False)
            print("Thanks for your response, the file was saved.")
    # if it doesn't exist, allow the save to happen
    else:        
        df.to_csv(path_fn, index=False)
        print("File successfully saved")


# [ make_dustflx_csv helper ]
def expand_to_years(lst, n, lstname="list"):
    if len(lst) == 1 and lst:
        return lst * n
    elif len(lst) != n:
        print(f"Warning: {lstname} length {len(lst)} does not match expected {n}. Will likely cause problems.")
    return lst

# [ make_dustflx_csv helper ]
def add_dust_row(
        y, duration, sp, rate, sp_2nd, rate_2nd,
        dustsp=None, dustrate=None,
        dustsp_2nd=None, dustrate_2nd=None
):
    tmprow = {
        "yr_start": y,
        "duration": duration,
    }
    # add keys if lists are not all none
    if dustsp and not all(v is None for v in dustsp):
        tmprow["dustsp"] = sp
    if dustrate and not all(v is None for v in dustrate):
        tmprow["dustrate"] = rate
    if dustsp_2nd and not all(v is None for v in dustsp_2nd):
        tmprow["dustsp_2nd"] = sp_2nd
    if dustrate_2nd and not all(v is None for v in dustrate_2nd):
        tmprow["dustrate_2nd"] = rate_2nd
    return tmprow


def make_dustflx_csv(
        savehere: str,
        savename: str,
        max_time: int,
        application_years: list,
        dustsp: list,
        dustrate: list,
        dustrad: list,
        dustsp_2nd: list,
        dustrate_2nd: list,
        returndf: bool=False,
        savedf: bool=True,
    ):
    '''
    Create a dust flux input .csv file for SCEPTER multi-year runs.
    The application_years list defines where a new dust application starts.
    All other lists must either be length 1 (for constant values) or
    length equal to the number of application years (for varying values).

    Parameters
    ----------
    savehere : str
        path to save the output .csv file
    savename : str
        name of the output .csv file
    max_time : int
        maximum time (years) for the total simulation
    application_years : list
        list of years where a new dust application starts 
        (all other years are assumed to have no dust application)
    dustsp : list
        list of dust species applied at each application year
        (use length 1 for constant dustsp or 
        len(application_years) for varying dustsp)
    dustrate : list
        list of dust application rates at each application year
        (use length 1 for constant dustrate or 
        len(application_years) for varying dustrate)
    dustrad : list
        list of dust particle sizes at each application year
        (use length 1 for constant dustrad or 
        len(application_years) for varying dustrad)
    dustsp_2nd : list
        list of second dust species applied at each application year
        (use length 0 for default, length 1 for constant dustsp_2nd or 
        len(application_years) for varying dustsp_2nd)
    dustrate_2nd : list
        list of second dust application rates at each application year
        (use length 0 for default, length 1 for constant dustrate_2nd or 
        len(application_years) for varying dustrate_2nd)
    '''
    # sort application_years
    yrs = sorted(set(int(y) for y in application_years))
    # --- handle lists ( determine if singleton or per-app values )
    dustsp       = expand_to_years(dustsp, len(yrs), lstname="dustsp")
    dustrate     = expand_to_years(dustrate, len(yrs), lstname="dustrate")
    dustrad      = expand_to_years(dustrad, len(yrs), lstname="dustrad")
    dustsp_2nd   = expand_to_years(dustsp_2nd, len(yrs), lstname="dustsp_2nd")
    dustrate_2nd = expand_to_years(dustrate_2nd, len(yrs), lstname="dustrate_2nd")
    # --- set up loop
    rows = []
    i=0
    n=len(yrs)
    # [ loop ]
    while i < n:
        y = yrs[i]
        sp = dustsp[i]
        rate = dustrate[i]
        sp_2nd = dustsp_2nd[i]
        rate_2nd = dustrate_2nd[i]

        if i == n - 1:
            # last application year (subtract 1 from max_time bc duration starts at 0)
            duration = max_time - 1 - y
            if duration < 1: # no zero-duration rows
                i += 1
                continue
            else:
                rows.append(add_dust_row(y, duration, sp, rate, sp_2nd, rate_2nd,
                                         dustsp=dustsp, dustrate=dustrate, dustsp_2nd=dustsp_2nd, dustrate_2nd=dustrate_2nd))
                i += 1
                continue

        gap = yrs[i + 1] - y

        if gap == 1:
            # consecutive block: walk forward until the block ends
            j = i
            while j + 1 < n and yrs[j + 1] - yrs[j] == 1:
                j += 1
            last_in_block = yrs[j]
            duration = last_in_block - y  # could be 0 if only single year, otherwise >=1
            rows.append(add_dust_row(y, duration, sp, rate, sp_2nd, rate_2nd,
                                         dustsp=dustsp, dustrate=dustrate, dustsp_2nd=dustsp_2nd, dustrate_2nd=dustrate_2nd))
            i = j + 1
            # add the next row for zero application
            duration = yrs[i] - last_in_block
            if duration > 0:
                rows.append(add_dust_row(last_in_block, duration, sp, rate, sp_2nd, rate_2nd,
                                         dustsp=dustsp, dustrate=dustrate, dustsp_2nd=dustsp_2nd, dustrate_2nd=dustrate_2nd))
        elif gap > 1:
            # big gap: create two rows:
            # 1) application in first year (duration = 1)
            rows.append(add_dust_row(y, 1, sp, rate, sp_2nd, rate_2nd,
                                         dustsp=dustsp, dustrate=dustrate, dustsp_2nd=dustsp_2nd, dustrate_2nd=dustrate_2nd))
            # 2) post-application period with zero dustrate covering the rest of the gap
            rows.append(add_dust_row(y + 1, gap - 1, sp, 0, sp_2nd, 0,
                                         dustsp=dustsp, dustrate=dustrate, dustsp_2nd=dustsp_2nd, dustrate_2nd=dustrate_2nd))
            i += 1
        else:
            # gap == 0 should not happen because yrs are unique, but handle defensively
            rows.append(add_dust_row(y, 0, sp, rate, sp_2nd, rate_2nd,
                                         dustsp=dustsp, dustrate=dustrate, dustsp_2nd=dustsp_2nd, dustrate_2nd=dustrate_2nd))
            i += 1

    # create dataframe
    df_out = pd.DataFrame(rows)
    # save to csv
    if savedf:
        savepath = os.path.join(savehere, savename)
        df_out.to_csv(savepath, index=False)
        print(f"Dust flux .csv saved to: {savepath}")
    if returndf:
        return df_out

# %% 

# --- [ FUNCTIONS TO COMPUTE ROCK AMOUNT AND GRAIN SIZE GIVEN A FIXED SA ]
def collect_data_vars(**kwargs):
    data_vars={
            "lambda": (
                ["radius"], 
                kwargs.get('lambda_arr'),
                {
                    "long_name": "particle roughness factor (function of radius except for smooth or geometric surface area)",
                    "units": "[]",
                }
            ),
            "volume_sphere": (
                ["radius"], 
                kwargs.get('v_sphere'),
                {
                    "long_name": "volume of a single sphere",
                    "units": "m3",
                }
            ),
            "mass_sphere": (
                ["radius"], 
                kwargs.get('m_sphere'),
                {
                    "long_name": "mass of one sphere",
                    "units": "kg",
                }
            ),
            "sa_sphere": (
                ["radius"], 
                kwargs.get('sa_sphere'),
                {
                    "long_name": "surface area of one sphere",
                    "units": "m2",
                }
            ),
            "n_spheres": (
                ["radius"], 
                kwargs.get('n_spheres'),
                {
                    "long_name": "number of spheres required to hit target surface area (accounting for roughness)",
                    "units": "[]",
                }
            ),
            "mass_rock": (
                ["radius"], 
                kwargs.get('m_rock'),
                {
                    "long_name": "mass of rock required to hit target surface area",
                    "units": "kg",
                }
            ),
        }
    if 'sa_fixed' not in kwargs.get('other_output_dims'):
        data_vars['sa_fixed'] = (
                [],
                kwargs.get('SA_FIXED'),
                {
                    "long_name": "constant surface area for all radius - app amount combinations",
                    "units": "m2",
                }
            )
    if 'lambda_model' not in kwargs.get('other_output_dims'):
        data_vars['lambda_model'] = (
            [],
            kwargs.get('lambda_model'),
            {
                "long_name": "Name of the model used to calculate lambda",
                "units": "[]"
            }
        )
    if 'feedstock' not in kwargs.get('other_output_dims'):
        data_vars['feedstock'] = (
            [],
            kwargs.get('feedstock'),
            {
                "long_name": "Name of rock feedstock",
                "units": "[]"
            }
        )
    return data_vars

def collect_coords(**kwargs):
    coords={
        "radius": (
            ["radius"],
            kwargs.get('radii'),
            {
                "long_name": "grain radius",
                "units": "m"
            }
        )
    }
    return coords


# [ prepare dicts for calculation ]
# set densities
scepter_densities = { # [g cm-3]
    'amnt': 1.725,   
    'cc': 2.71,
    'gbas': 3.0,   # confirmed in scepter_richards.f90 "! assuming 3.0 g/cm3 particle density"
    'gbas2': 3.0,  # gbas2 holds mwt fixed at mwtgbas, so mv = mwt/3.0 gives the same density as gbas
}

# [ roughness factor equations ]
def nsb07(r):
    return 10**3.3 * r**0.33
def bm00(r):
    return 10**0.7 * r**-0.1
def l21a(r):
    return np.maximum(2.02*np.log10(r) + 10.126, 1)
def l21b(r): # ridiculously high numbers
    return 10**(154.25 * np.exp(1.0219 * np.log10(r)))
def smooth(r):
    return np.repeat(1., len(r))



# pull lambda functions together
lambda_fxns = {
    'nsb07': nsb07,
    'bm00': bm00,
    'l21a': l21a,
    'l21b': l21b,
    'smooth': smooth,
}

def rockmass_given_sa(
        SA_FIXED: float,   # [m2]
        FEEDSTOCK: str,    # []
        radii: np.array,     # [m]
        lambda_model: str, # [] name of the roughness factor eq'n
        scepter_densities: dict = scepter_densities, # [g cm-3]
        lambda_fxns: dict = lambda_fxns,
        other_output_dims: list = ['lambda_model', 'feedstock', 'sa_fixed'],
):
    """
    Given a target surface area, radius, and roughness factor, calculate the
    amount of rock we need to add. Account for rock density and roughness
    factor equation (defined in external fxns).

    Parameters
    ----------
    SA_FIXED : float
        [m2] target total surface area of applied rock in m2
    FEEDSTOCK : str
        [] name of the feedstock to use (must be in "scepter_densities" dict)
    radii : np.array
        [m] 1d array of particle radii input
    lambda_model : str
        [] name of the roughness factor (lambda) parameterization to use
    scepter_densities : dict
        [g cm-3] feedstock densities; ex: {'gbas': 3.0}; where values 
        are in g/cm3 and converted to kg/m3 in the function
    lambda_fxns : dict
        [] dictionary of functions that read in radius array and return lambda array
    other_output_dims : list
        [] list of strings naming the output dimensions to add to the xr.dataset. 
        Radius is a guaranteed dim. OPTIONS: ["lambda_model", "feedstock", "sa_fixed"]
    """
    # [ roughness factor ]
    lambda_arr = lambda_fxns[lambda_model](radii)
    # [ volume of one sphere ]
    v_sphere = (4/3) * 3.14 * radii**3    # [m3]
    # [ mass of one sphere ]
    m_sphere = (scepter_densities[FEEDSTOCK] * 1e3) * v_sphere # [kg]
    # [ surface area of one sphere ]
    sa_sphere = 4 * 3.14 * radii**2    # [m2]
    # [ number of spheres to hit target SA ]
    n_spheres = SA_FIXED / (sa_sphere * lambda_arr)
    # [ required mass of rock ]
    m_rock = n_spheres * m_sphere

    # put together the output data
    data_vars = collect_data_vars(**locals())
    coords = collect_coords(**locals())
    ds = xr.Dataset(data_vars=data_vars, coords=coords)
    # add extra dims
    if 'lambda_model' in other_output_dims:
        ds = ds.assign_coords(lambda_model=lambda_model).expand_dims("lambda_model")
    if 'feedstock' in other_output_dims:
        ds = ds.assign_coords(feedstock=FEEDSTOCK).expand_dims("feedstock")
    if 'sa_fixed' in other_output_dims:
        ds = ds.assign_coords(sa_fixed=SA_FIXED).expand_dims("sa_fixed")

    return ds


# %% ──────────────────────────────────────────────────────────────────────────
# [ DGSA CLIMATE DIRECTORY BUILDER ]
# ─────────────────────────────────────────────────────────────────────────────

def build_dgsa_climatedirs(
    df_treat: pd.DataFrame,
    base_clim_s3: str,
    clim_selection_label: str,
    exp_set: str,
    dgsa_version: str,
    ref_apprate: float,
    apprate_suffix: str,
    runtype: str,
    deg_resolution: int = 1,
    dryrun: bool = False,
    save_validation_plots: bool = False,
    verbose: bool = True,
) -> pd.Series:
    """
    For each row in df_treat, copy the reference fixed-apprate climate dir to a
    DGSA-specific destination named by a unique integer identifier (paramSet_N)
    and write scaled versions of Dust_temp.in and q_temp.in.

    df_treat must contain a 'paramSet_id' column (1-indexed integers) that
    provides the stable identifier for each treatment row.

    Source dir:
        {base_clim_s3}/{clim_selection_label}/{exp_set}/{ref_apprate}{apprate_suffix}/
            {region}_{deg}_{runtype}/

    Destination dir:
        {base_clim_s3}/_dgsa/{dgsa_version}/{clim_selection_label}/{exp_set}/
            paramSet_{paramSet_id}/{region}_{deg}_{runtype}/

    A summary CSV is also saved to the parent climatedir:
        {base_clim_s3}/_dgsa/{dgsa_version}/{clim_selection_label}/{exp_set}/
            _dgsa-treatment-cases.csv

    Modifications applied to the destination:
        Dust_temp.in  — dust column × (apprate_base / ref_apprate)
        q_temp.in     — runoff column × qrun_factor
        All other files (T_temp.in, Wet_temp.in, clim-timeseries.png) are
        copied unchanged.

    Parameters
    ----------
    df_treat : pd.DataFrame
        Treatment rows from the DGSA batch dataframe. Must contain columns
        'region_nospaces', 'apprate_base', 'qrun_factor', and 'paramSet_id'.
    base_clim_s3 : str
        Base S3 path for climate dirs (no trailing slash).
        e.g. "s3://carbonplan-carbon-removal/ew-workflows-data/scepter/clim/one-cell-per-region_1deg"
    clim_selection_label : str
        Label for the climate selection (e.g. "repcells1").
    exp_set : str
        Experiment set name (e.g. "exe_central1p5").
    dgsa_version : str
        Version tag for this DGSA run so repeated runs don't overwrite each
        other (e.g. "v1.0").
    ref_apprate : float
        Application rate of the existing source climate dir. Must match a dir
        already present at {base_clim_s3}/{clim_selection_label}/{exp_set}/.
    apprate_suffix : str
        Suffix appended to the apprate in directory names
        (e.g. "-tHaYr-appRate-FIXED").
    runtype : str
        Climate run type (e.g. "yearly_ltm").
    deg_resolution : int
        Spatial resolution label used in region sub-directory names (default 1).
    dryrun : bool
        If True, compute paths and print what would be done, but skip all S3
        reads and writes. Useful for verifying path logic before committing.
    save_validation_plots : bool
        If True, generate and save two PNG validation plots to the destination
        directory alongside the climate files:
          scaling-validation_dusttemp.png — instantaneous rates, cumulative
              rock application, and per-year integrals for Dust_temp.in.
          scaling-validation_qtemp.png   — runoff time series and ratio check
              for q_temp.in.
        Ignored when dryrun=True.
    verbose : bool
        Print progress for each unique combo processed.

    Returns
    -------
    pd.Series
        Indexed like df_treat; each entry is the full S3 path of the
        destination climatedir for that row (trailing slash included).
    """
    import s3fs

    def _load_tabfile(path):
        with fsspec.open(path, "r") as fh:
            raw = fh.read()
        header = "\n".join(ln for ln in raw.splitlines() if ln.startswith("#"))
        data = np.loadtxt(io.StringIO(raw), comments="#")
        return header, data[:, 0], data[:, 1]

    def _format_tabfile(header, times, values):
        lines = [header] if header else []
        for t, v in zip(times, values):
            lines.append(f"{t:18.8e}\t{v:18.8e}")
        return "\n".join(lines) + "\n"

    def _save_fig_to_s3(fig, dest_path, fs):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        buf.seek(0)
        with fs.open(dest_path.replace("s3://", ""), "wb") as fh:
            fh.write(buf.read())

    def _plot_dust_validation(time_yr, dust_src, dust_scaled, ref_apprate, apprate, region):
        import matplotlib.pyplot as plt

        dt = np.empty_like(time_yr)
        dt[:-1] = np.diff(time_yr)
        dt[-1]  = dt[-2]
        cum_src    = np.cumsum(dust_src    * dt) / 100
        cum_scaled = np.cumsum(dust_scaled * dt) / 100

        n_years = int(np.floor(time_yr[-1])) + 1
        yrs, ann_src, ann_scaled = [], [], []
        for yr in range(n_years):
            mask = (time_yr >= yr) & (time_yr < yr + 1)
            if not mask.any():
                continue
            yrs.append(yr)
            ann_src.append(   np.sum(dust_src[mask]    * dt[mask]) / 100)
            ann_scaled.append(np.sum(dust_scaled[mask] * dt[mask]) / 100)
        yrs        = np.array(yrs)
        ann_src    = np.array(ann_src)
        ann_scaled = np.array(ann_scaled)
        actual_ratio   = ann_scaled.mean() / ann_src.mean() if ann_src.mean() > 0 else np.nan
        expected_ratio = apprate / ref_apprate
        pct_err        = 100 * abs(actual_ratio - expected_ratio) / expected_ratio

        fig, axes = plt.subplots(1, 3, figsize=(16, 4))

        n_plot = min(3, n_years)
        mask   = time_yr < n_plot
        ax = axes[0]
        ax.plot(time_yr[mask], dust_src[mask]    / 100, lw=1,       label=f"source  ({ref_apprate} t/ha/yr)")
        ax.plot(time_yr[mask], dust_scaled[mask] / 100, lw=1, ls="--", label=f"scaled  ({apprate!r} t/ha/yr)")
        ax.set_xlabel("Time (yr)")
        ax.set_ylabel("Rate (t/ha/yr)")
        ax.set_title(f"Dust_temp.in — rates (first {n_plot} yr)")
        ax.legend(fontsize=8)

        t_end = time_yr[-1]
        ax = axes[1]
        ax.plot(time_yr, cum_src,    lw=1.5,          label="source")
        ax.plot(time_yr, cum_scaled, lw=1.5, ls="--", label="scaled")
        ax.plot([0, t_end], [0, ref_apprate * t_end], color="gray",   ls=":", lw=1, label="expected ref slope")
        ax.plot([0, t_end], [0, apprate     * t_end], color="orange", ls=":", lw=1, label="expected target slope")
        ax.set_xlabel("Time (yr)")
        ax.set_ylabel("Cumulative (t/ha)")
        ax.set_title("Time-integrated rock application")
        ax.legend(fontsize=8)

        ax = axes[2]
        ax.plot(yrs, ann_src,    "o", ms=2, label=f"source  (expect {ref_apprate})")
        ax.plot(yrs, ann_scaled, "s", ms=2, label=f"scaled  (expect {apprate!r})")
        ax.axhline(ref_apprate, color="tab:blue",   ls="--", lw=0.8)
        ax.axhline(apprate,     color="tab:orange", ls="--", lw=0.8)
        ax.set_xlabel("Year")
        ax.set_ylabel("Annual rock (t/ha)")
        ax.set_title(
            f"Per-year integrals\n"
            f"actual×{actual_ratio:.5f}  expected×{expected_ratio:.5f}  err={pct_err:.3f}%"
        )
        ax.legend(fontsize=8)

        fig.suptitle(
            f"Dust_temp.in scaling — {region}  "
            f"({ref_apprate} → {apprate!r} t/ha/yr)",
            y=1.02,
        )
        plt.tight_layout()
        return fig

    def _plot_qtemp_validation(q_time, q_src, q_scaled, qrun_factor, region):
        import matplotlib.pyplot as plt

        n_plot_yr = min(5, int(q_time[-1]))
        mask      = q_time < n_plot_yr
        ratio     = q_scaled / q_src

        fig, axes = plt.subplots(1, 2, figsize=(11, 3.5))

        ax = axes[0]
        ax.plot(q_time[mask], q_src[mask],    lw=1.5,          label="source  (×1.0)")
        ax.plot(q_time[mask], q_scaled[mask], lw=1.5, ls="--", label=f"scaled  (×{qrun_factor!r})")
        ax.set_xlabel("Time (yr)")
        ax.set_ylabel("Runoff (mm/month)")
        ax.set_title(f"q_temp.in (first {n_plot_yr} yr)")
        ax.legend(fontsize=8)

        ax = axes[1]
        ax.plot(q_time[mask], ratio[mask], lw=1)
        ax.axhline(qrun_factor, color="red", ls="--", lw=1, label=f"expected {qrun_factor!r}")
        ax.set_xlabel("Time (yr)")
        ax.set_ylabel("scaled / source")
        ax.set_title(f"Ratio check — mean={ratio.mean():.6f}  (expect {qrun_factor!r})")
        ax.legend(fontsize=8)

        fig.suptitle(
            f"q_temp.in scaling — {region}  (qrun_factor={qrun_factor!r})",
            y=1.02,
        )
        plt.tight_layout()
        return fig

    if "paramSet_id" not in df_treat.columns:
        raise ValueError("df_treat must have a 'paramSet_id' column before calling build_dgsa_climatedirs")

    fs = None if dryrun else s3fs.S3FileSystem(anon=False)

    parent_dir = (
        f"{base_clim_s3}/_dgsa/{dgsa_version}/{clim_selection_label}/{exp_set}/"
    )
    id_to_dest: dict = {}

    for _, row in df_treat.iterrows():
        param_set_id = row["paramSet_id"]
        region       = row["region_nospaces"]
        apprate      = row["apprate_base"]
        qrun_factor  = row["qrun_factor"]
        dust_scale   = apprate / ref_apprate

        region_subdir = f"{region}_{deg_resolution}_{runtype}"
        source_dir = (
            f"{base_clim_s3}/{clim_selection_label}/{exp_set}/"
            f"{ref_apprate}{apprate_suffix}/{region_subdir}/"
        )
        paramset_dir = f"{parent_dir}paramSet_{param_set_id}/"
        dest_dir = f"{paramset_dir}{region_subdir}/"
        id_to_dest[(param_set_id, region)] = paramset_dir

        if verbose:
            print(
                f"  paramSet_{param_set_id} ({region}): apprate={apprate!r}, qrun_factor={qrun_factor!r} "
                f"→ dust×{dust_scale:.4f}"
            )
        if dryrun:
            if verbose:
                print(f"    [dryrun] dest: {dest_dir}")
            continue

        # Load + scale Dust_temp.in
        dust_header, time_yr, dust_src = _load_tabfile(source_dir + "Dust_temp.in")
        dust_scaled = dust_src * dust_scale

        # Load + scale q_temp.in
        q_header, q_time, q_src = _load_tabfile(source_dir + "q_temp.in")
        q_scaled = q_src * qrun_factor

        # Write scaled files
        for dest_path, header, times, values in [
            (dest_dir + "Dust_temp.in", dust_header, time_yr, dust_scaled),
            (dest_dir + "q_temp.in",    q_header,    q_time,  q_scaled),
        ]:
            content = _format_tabfile(header, times, values)
            with fs.open(dest_path.replace("s3://", ""), "w") as fh:
                fh.write(content)
            if verbose:
                print(f"    wrote {dest_path.split('/')[-1]}")

        # Copy unchanged files
        for fn in ["T_temp.in", "Wet_temp.in", "clim-timeseries.png"]:
            src_key = (source_dir + fn).replace("s3://", "")
            dst_key = (dest_dir   + fn).replace("s3://", "")
            if fs.exists(src_key):
                fs.copy(src_key, dst_key)
                if verbose:
                    print(f"    copied {fn}")
            elif verbose:
                print(f"    skipped {fn} (not found at source)")

        # Validation plots
        if save_validation_plots:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig_dust = _plot_dust_validation(time_yr, dust_src, dust_scaled, ref_apprate, apprate, region)
            _save_fig_to_s3(fig_dust, dest_dir + "scaling-validation_dusttemp.png", fs)
            plt.close(fig_dust)

            fig_q = _plot_qtemp_validation(q_time, q_src, q_scaled, qrun_factor, region)
            _save_fig_to_s3(fig_q, dest_dir + "scaling-validation_qtemp.png", fs)
            plt.close(fig_q)

            if verbose:
                print("    saved scaling-validation_dusttemp.png")
                print("    saved scaling-validation_qtemp.png")

    # Save treatment cases CSV parallel to the paramSet_* dirs.
    csv_path = parent_dir + "_dgsa-treatment-cases.csv"
    df_save = df_treat.copy()
    df_save["climatedir"] = df_treat.apply(
        lambda r: id_to_dest[(r["paramSet_id"], r["region_nospaces"])], axis=1
    )
    if dryrun:
        if verbose:
            print(f"[dryrun] would save treatment cases CSV → {csv_path}")
    else:
        csv_bytes = df_save.to_csv(index=False).encode()
        with fs.open(csv_path.replace("s3://", ""), "wb") as fh:
            fh.write(csv_bytes)
        if verbose:
            print(f"saved treatment cases CSV → {csv_path}")

    if dryrun and verbose:
        print("build_dgsa_climatedirs: dryrun=True — no files written")

    # Map every row in df_treat back to its destination path.
    return df_treat.apply(
        lambda r: id_to_dest[(r["paramSet_id"], r["region_nospaces"])], axis=1
    )
