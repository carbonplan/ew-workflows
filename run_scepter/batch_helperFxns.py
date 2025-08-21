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