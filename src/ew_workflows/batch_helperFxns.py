# --------------------------------------------
# 
# helper functions for building the batch
# input .csv files
# 
# --------------------------------------------
import os
import pandas as pd
import itertools
import warnings

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
        site = row['climatefiles']
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