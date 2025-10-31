# ------------------------------------
# functions to support making 
# batch files for SMEW
# ------------------------------------
import os
import itertools
import numpy as np
import pandas as pd
from scipy.stats import qmc
import warnings

# --- function to sample parameter space with latin hypercube
def latin_hypercube_sampler(
        parameter_ranges: dict,
        n_samples: int,
        nonnum_repeat_type: str,
        round_to: int=2,
)-> pd.DataFrame:
    '''
    Read in a dictionary of parameter ranges where each key is the parameter
    name (as listed in defaults/default_dicts.py) and each value is a 
    list (e.g., `[0,1]`). Given a number of samples, return parameter values 
    for each sample. 

    Parameters
    ----------
    parameter_ranges : dict
        dictionary with keys equal to default_dicts.py parameters 
        (must be numeric!) and values indicating a min / max range.
    n_samples : int
        number of samples to generate (each sample amounts to a given
        simulation of the PSM)
    nonnum_repeat_type : str
        ["prescribed_cases" | "all_combinations"] how to repeat the numeric 
        latin hypercube samples across the non-numeric parameters. "prescribed_cases"
        means we repeat the lhs once for each index in the non-numeric params. 
        "all_combinations" means we repeat the lhs once for every possible 
        combination of the non-numeric params. 
    round_to : int
        number of decimal places to round to. 

    Returns
    -------
    pd.DataFrame 
        dataframe where each column is a parameter and each row is a sample
    '''
    # separate numeric and non-numeric elements
    numeric_params = {key: value for key, value in parameter_ranges.items() if all(isinstance(x, (int, float)) for x in value)}
    non_numeric_params = {key: value for key, value in parameter_ranges.items() if not all(isinstance(x, (int, float)) for x in value)}

    # -----
    # make sure each value has only two elements
    for key, values in numeric_params.items():
        if len(values) != 2:
            raise ValueError(f"The key '{key}' has a length of {len(values)}, but it must be 2.")
    # -----
    # get the number of parameters
    num_parameters = len(numeric_params)

    # create a latin hypercube sampler
    sampler = qmc.LatinHypercube(d=num_parameters)

    # generate the samples in the unit hypercube
    sample = sampler.random(n=n_samples)
    
    # extract the ranges from the dictionary
    keys = list(numeric_params.keys())
    min_vals = [numeric_params[key][0] for key in keys]
    max_vals = [numeric_params[key][1] for key in keys]

    # scale the samples to the desired parameter ranges
    scaled_sample = qmc.scale(sample, min_vals, max_vals)

    # convert back to a dict
    parameter_values = {keys[i]: np.round(scaled_sample[:, i], round_to) for i in range(num_parameters)}
    # and now to a pandas dataframe
    df_num = pd.DataFrame(parameter_values)

    # --- handle the non-numeric values
    if len(non_numeric_params) > 0:
        if nonnum_repeat_type == "prescribed_cases": 
            df_nonnum = pd.DataFrame(non_numeric_params)
        elif nonnum_repeat_type == "all_combinations":
            # generate all combinations of parameter values
            combinations = list(itertools.product(*non_numeric_params.values()))
            df_nonnum = pd.DataFrame(combinations, columns=non_numeric_params.keys())
        # merge dfs
        # create a dummy key for cross join
        df_num['key'] = 1
        df_nonnum['key'] = 1
        # merge the DataFrames on the dummy key
        dfout = pd.merge(df_num, df_nonnum, on='key').drop('key', axis=1)
    else:
        dfout = df_num

    # return result
    return dfout


# --- function to sample all possible values
def all_combinations_sampler(
        parameter_values: dict,
)-> pd.DataFrame:
    '''
    Read in a dictionary of parameter values where each key is the parameter
    name (as listed in defaults/default_dicts.py) and each value is a 
    list of values to test. Output a pd.DataFrame that includes all 
    possible combinations of the listed parameter values.

    Parameters
    ----------
    parameter_ranges : dict
        dictionary with keys equal to default_dicts.py parameters 
        (must be numeric!) and values are values to test.
    
    Returns
    -------
    pd.DataFrame 
        dataframe where each column is a parameter and each row is a sample
    '''
    # check if all parameter value lengths are 2 (this might
    # indicate that the dict is for latin hypercube sampling)
    all_len_2 = all(len(value) == 2 for value in parameter_values.values())
    if all_len_2:
        warnings.warn("All lists in the dict have a length of 2. If these are min / max values, you may have meant to use the Latin hypercube sampler!", UserWarning)

    # generate all combinations of parameter values
    combinations = list(itertools.product(*parameter_values.values()))

    # create a DataFrame with parameter names as columns
    df = pd.DataFrame(combinations, columns=parameter_values.keys())
    # return result
    return df

# --- prescribed cases
def prescribed_cases(
        parameter_cases: dict,
)->pd.DataFrame:
    '''
    Create pandas dataframe where each row is an individual experimental
    setup that aligns 1:1 with the structure of the parameter_cases
    dictionary

    Parameters
    ----------
    parameter_cases : dict
        dictionary with keys equal to default_dicts.py parameters 
        (must be numeric!) and values are values to test. Each 
        parameter value array must be the same length. 
    
    Returns
    -------
    pd.DataFrame 
        dataframe where each column is a parameter and each row is a sample
    '''
    # confirm that all lists are the same length
    all_same_length = all(len(value) == len(next(iter(parameter_cases.values()))) for value in parameter_cases.values())
    if not all_same_length:
        raise ValueError(f"All parameter value lists must be the same length for prescribed_cases!")
    
    # convert to pandas dataframe
    dfout = pd.DataFrame(parameter_cases)

    # return result
    return dfout

# --- one at a time
def one_at_a_time_sensitivity(
        parameter_values: dict,
        default_str: str = "**default**",
)-> pd.DataFrame:
    '''
    Create pandas dataframe where each row is a test of a single value from a 
    single parameter, setting all else to default.

    Parameters
    ----------
    parameter_cases : dict
        dictionary with keys equal to default_dicts.py parameters 
        (must be numeric!) and values are values to test. Each 
        parameter value array must be the same length. 
    default_str : str
        name to use if the value for the default dictionary should be used.
        CAUTION: this might be hard-coded in the helper_functions.py to 
        expect a certain value! Only change if you know what you're doing.
    
    Returns
    -------
    pd.DataFrame 
        dataframe where each column is a parameter and each row is a sample
    '''
    # initialize an empty row
    rows = []
    # loop through each key and value in the parameter dict
    for key, values in parameter_values.items():
        for val in values:
            # set all rows to default
            row = {col: default_str for col in parameter_values.keys()}
            row[key] = val # over-write the default value with the given parameter value
            rows.append(row) 
    # return the dataframe
    return pd.DataFrame(rows)


# --- function to add constant values to the dict
def add_constant_parameters(
        df_batch: pd.DataFrame,
        constant_dict: dict,
)->pd.DataFrame:
    '''
    Take in the existing batch dataframe and add the constant parameter 
    values. 

    Parameters
    ----------
    df_batch : pd.DataFrame
        the pandas dataframe that is output from one of the other sample 
        functions (latin_hypercube_sampler, all_combinations_sampler,
        prescribed_cases).
    constant_dict : dict
        dictionary where keys are parameter names and each has a single 
        value that is held constant for all rows. 

    Returns
    -------
    pd.DataFrame
        the final batch .csv that gets saved
    '''
    # check that all dicts have only one value
    all_len_1 = all((isinstance(value, (list, tuple)) and len(value) == 1) or not isinstance(value, (list, tuple)) for value in constant_dict.values())
    if not all_len_1:
        warnings.warn("Expected all dict parameters to have one value but at least one parameter has more. This may lead to unintended results.", UserWarning)

    # constant DataFrames
    df2 = pd.DataFrame([constant_dict])

    # Create a dummy key for cross join
    df_batch['key'] = 1
    df2['key'] = 1

    # Merge the DataFrames on the dummy key
    combined_df = pd.merge(df_batch, df2, on='key').drop('key', axis=1)
    return combined_df
