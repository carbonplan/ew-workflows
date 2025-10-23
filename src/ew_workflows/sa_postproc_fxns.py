from scipy.integrate import cumulative_trapezoid
from typing import Tuple
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
from SALib.sample import morris as morris_sampler
from SALib.analyze import morris as morris_analyzer



# [ save dict for sensitivity analysis (or single feedstock) case ]
def create_save_dict_sa(
    cdr_calc_list: list,
    group_vars: list,
    csv_fn: str,
    multiyear: bool, 
    time_horizon: float, 
    cf_apprate_fixed: float,
    cf_dustrad_fixed: float, 
    p80_input: float,
    truck_km: float,
    barge_km: float, 
    barge_diesel_km: float, 
    Efactor_org: float,
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
            "csv_fn": csv_fn,
            "multiyear": multiyear,
        },
        "Removal accounting": {
            "time_horizon": time_horizon,
            "cf_apprate_fixed": cf_apprate_fixed,
            "cf_dustrad_fixed": cf_dustrad_fixed
        },
        "Emissions": {
            "p80_input": p80_input,
            "truck_km": truck_km,
            "barge_km": barge_km,
            "barge_diesel_km": barge_diesel_km,
            "Efactor_org": Efactor_org
        },
        
    }

    return savedict




def save_batch_postproc_sa(
    savedict: dict, 
    base_path: str, 
    base_dir_name: str,
    fname_res: str,
    cdr_dict_full: dict,
    cdr_dict_sum: dict,
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
    cdr_dict_full : dict
        dictionary of pandas dataframes with full CDR outputs (over time)
    cdr_dict_sum : dict
        dictionary of pandas dataframes with CDR outputs — these are 
        ultimately input to the emissions calculation step (we save them
        to make it easier to play with other emissions scenarios after
        this step)
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
    savedf_fn_sum = os.path.join(outdir, "cdr_dfs_sum.pkl")
    savedf_fn_full = os.path.join(outdir, "cdr_dfs_full.pkl")

    # ----------------------------
    # check for S3
    if savedf_fn_sum.startswith("s3://"): # then bring it in from s3
        import fsspec
        with fsspec.open(savedf_fn_sum, 'wb') as f:
            pickle.dump(cdr_dict_sum, f)
    else:
        with open(savedf_fn_sum, 'wb') as f:
            pickle.dump(cdr_dict_sum, f)

    # repeat for the full dataset 
    if savedf_fn_full.startswith("s3://"): # then bring it in from s3
        import fsspec
        with fsspec.open(savedf_fn_full, 'wb') as f:
            pickle.dump(cdr_dict_sum, f)
    else:
        with open(savedf_fn_full, 'wb') as f:
            pickle.dump(cdr_dict_sum, f)
    # ----------------------------
    # NOTE: you can open the dfs like this:
    # with open('dataframes.pkl', 'rb') as f:
        # loaded_dfs = pickle.load(f) 

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



# --- sensitivity analysis functions -----------------

varname_dict_in = {
    "dustrate" : "Application flux",
    "taudust" : "Application duration",
    "dustrad" : "Dust radius",
    "qrun": "Runoff",
    "mat" : "Mean annual temperature",
    "dust_mixdep" : "Dust mixing depth",
    "dustrate_2nd": "Fertilizer flux"
}



def reorder_rows(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
)->pd.DataFrame:
    '''
    Reorder rows of df2 based on the order of rows in the matching
    columns of df1. Return updated df2. 
    '''
    # --- get column in the right order 
    common_cols = df1.columns.intersection(df2.columns)
    # check that they're in the same order 
    # (np.allclose should account for floating point errors)
    same_order = np.allclose(df1[common_cols].values, df2[common_cols].values, equal_nan=True)
    print("Same order:", same_order)
    if not same_order:
        # set the index temporarily to those matching columns
        df1_indexed = df1.set_index(list(common_cols))
        df2_indexed = df2.set_index(list(common_cols))

        # reorder df2 using df1's index
        df2 = df2_indexed.loc[df1_indexed.index].reset_index()
    
    # return result
    return df2

# -----------------------------------------------------------
def sa_quickplot_scatter(
        SA_result,
        problem: dict,
        ptitle: str= "Morris Sensitivity Analysis",
        saveme: bool=False,
        savedir: str=None,
        savename: str=None,
):
    '''
    quickly visualize SA results
    '''
    # plot
    fig, ax = plt.subplots()
    ax.errorbar(SA_result["mu_star"], SA_result["sigma"], xerr=SA_result["mu_star_conf"], fmt="o")
    ax.set_xlabel("μ* (mean absolute effect)")
    ax.set_ylabel("σ (standard deviation of effects)")
    if ptitle is not None:
        ax.set_title(ptitle)
    ax.grid(True)
    for i, name in enumerate(problem["names"]):
        ax.annotate(name, (SA_result["mu_star"][i], SA_result["sigma"][i]), textcoords="offset points", xytext=(5,5))
    plt.tight_layout()
    
    if saveme:
        os.makedirs(savedir, exist_ok=True)
        plt.savefig(os.path.join(savedir, savename), dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
        plt.close()


def sa_quickplot_scatter_compareMu(
        SA_result,
        problem: dict,
        ptitle: str= "Morris Sensitivity Analysis",
        saveme: bool=False,
        savedir: str=None,
        savename: str=None,
):
    '''
    quickly visualize SA results
    '''
    # 1:1 and -1:1 inputs
    lims = [0, np.max([SA_result["mu"], SA_result["mu_star"]])] 

    # plot
    fig, ax = plt.subplots()
    # 1:1 line
    ax.plot(lims, lims, color='gray', linestyle='--', label='1:1')
    # -1:1 line (negative slope)
    ax.plot([-x for x in lims], lims, color='gray', linestyle='--', label='1:1')
    # results
    ax.errorbar(SA_result["mu"], SA_result["mu_star"], xerr=SA_result["mu_star_conf"], fmt="o")
    ax.set_xlabel("μ (mean effect)")
    ax.set_ylabel("μ* (mean absolute effect)")
    if ptitle is not None:
        ax.set_title(ptitle)
    ax.grid(True)
    for i, name in enumerate(problem["names"]):
        ax.annotate(name, (SA_result["mu"][i], SA_result["mu_star"][i]), textcoords="offset points", xytext=(5,5))
    plt.tight_layout()
    
    if saveme:
        os.makedirs(savedir, exist_ok=True)
        plt.savefig(os.path.join(savedir, savename), dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
        plt.close()


def sa_quickplot_bar(
        SA_result,
        varname_dict: dict=varname_dict_in,
        barcolor: str= 'steelblue',
        linecolor: str='black',
        bar_order: list=None,
        ptitle: str= "Morris Sensitivity Analysis",
        fs_ylab: int | float = 12,
        fs_xlab: int | float = 12,
        fs_xticks: int | float = 10,
        saveme: bool=False,
        savedir: str=None,
        savename: str=None,
): 
    '''
    quick bar plot with vars ranked by magnitude
    '''
    # Extract values
    mu_star = np.array(SA_result['mu_star'])
    mu_star_conf = np.array(SA_result['mu_star_conf'])
    param_names = np.array(SA_result['names'])

    if bar_order is not None:
        param_names_list = list(param_names)
        order_idx = [param_names_list.index(p) for p in bar_order]
        mu_star_sorted = mu_star[order_idx]
        mu_star_conf_sorted = mu_star_conf[order_idx]
        param_names_sorted = param_names[order_idx]
    else:
        # Sort by descending mu_star
        sorted_idx = np.argsort(mu_star)[::-1]
        mu_star_sorted = mu_star[sorted_idx]
        mu_star_conf_sorted = mu_star_conf[sorted_idx]
        param_names_sorted = param_names[sorted_idx]

    # set y parameter names 
    if varname_dict is not None:
        # use mapped names for y-axis labels
        param_names_sorted_pretty = [varname_dict.get(name, name) for name in param_names_sorted]
    # generic names for mapping
    y_pos = np.arange(len(param_names))


    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.barh(y_pos, mu_star_sorted, xerr=mu_star_conf_sorted, align='center', color=barcolor, ecolor=linecolor)
    ax.set_yticks(y_pos)
    if varname_dict is not None:
        ax.set_yticklabels(param_names_sorted_pretty, fontsize = fs_ylab)
    else:
        ax.set_yticklabels(param_names_sorted, fontsize = fs_ylab)
    ax.invert_yaxis()  # Largest at the top
    ax.set_xlabel('Mean absolute effect (μ*)', fontsize = fs_xlab)
    ax.tick_params(axis='x', labelsize=fs_xticks)
    if ptitle is not None:
        ax.set_title(ptitle)

    plt.tight_layout()

    if saveme:
        os.makedirs(savedir, exist_ok=True)
        plt.savefig(os.path.join(savedir, savename), dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
        plt.close()


def Si_and_plot(
    problem: dict,
    param_values: list,
    thismetric: str,
    thiscase: str,
    Yout: np.array,
    saveplots: bool,
    plot_dict: dict=None,
    conf_level: float=0.95,
    print_si_to_console: bool=False,
):
    '''
    Calculate Sensitivity results and save the plot if asked
    '''
    # --- conduct Si and plot 
    Si = morris_analyzer.analyze(problem, param_values, Yout, conf_level=conf_level, print_to_console=print_si_to_console)
    tmpdf = pd.DataFrame(Si, index = problem['names'])
    tmpds = xr.Dataset.from_dataframe(tmpdf).drop_vars('names')
    tmpds = tmpds.expand_dims(metric = [thismetric])

    # --- check if we should plot this result 
    if saveplots:
        saveme = thismetric in plot_dict['plot_metrics']

    # --- plot it 
    if saveme:
        # make title and save location
        ptitle = f"{thiscase} | {thismetric}"
        savedir = os.path.join(plot_dict['savepath'], thismetric)
        
        # make and save plots 
        # [ scatter plot ]
        if "scatter" in plot_dict['plot_types']:
            savename = f"scatter---{thiscase}+{thismetric}.png"
            sa_quickplot_scatter(
                Si,
                problem,
                ptitle,
                saveme = True,
                savedir = savedir,
                savename = savename
            )
        
        # [ scatter plot compare mu ]
        if "scatter_mu" in plot_dict['plot_types']:
            savename = f"scatterMu---{thiscase}+{thismetric}.png"
            sa_quickplot_scatter_compareMu(
                Si,
                problem,
                ptitle,
                saveme = True,
                savedir = savedir,
                savename = savename
            )

        # [ scatter plot compare mu ]
        if "bar" in plot_dict['plot_types']:
            savename = f"bar---{thiscase}+{thismetric}.png"
            sa_quickplot_bar(
                Si,
                ptitle = ptitle,
                saveme = True,
                savedir = savedir,
                savename = savename
            )
        
    # return result
    return tmpds    



def Si_plot_batch_anoms(
    anom_outputs_dict: dict,
    outputs_dict: dict,
    outdir: str,
    cdr_dfs: str,
    problem: dict,
    param_values: list,
    saveplots: bool,
    plot_dict: dict=None,
    conf_level: float=0.95,
    print_si_to_console: bool=False,
)-> xr.Dataset:
    '''
    Batch sensitivity analysis and produce / save plots.

    Parameters
    ----------
    anom_outputs_dict : dict
        dictionary where the key is the name of two keys in outputs_dict,
        indicating which is subtracted from the other (e.g., `gbas_1yr - cc_1yr`)
        and Values are bools indicating whether to run that analysis
    outputs_dict : dict
        dictionary of outputs to consider. Keys are names of the output case 
        (e.g., `cc_1yr`) and values are the directory where outputs are 
        saved (e.g., `meanAnn_shortRun_SAFert_cc_morris_7_400_5yrInt_001`)
    outdir : str
        location of the output data (already processed)
    cdr_dfs : str
        name of the file with the dict of cdr data
    problem : dict
        SA problem dictionary
    param_values : list
        array of the parameter input values 
    saveplots : bool
        [True] to save figures (using info from the plot_dict)
    plot_dict : dict
        plotting information including `plot_types`; `plot_metrics`
        and `savepath`. 
    conf_level : float
        value to use for the confidence level in the SA 
    print_si_to_console : bool
        whether to print the SA results to console 
    
    Returns
    -------
    xr.Dataset
        Dataset with SA outputs for all cases / metrics etc. 
    '''
    # --- create param dataframe
    df_params = pd.DataFrame(param_values, columns = problem['names'])

    # --- loop through the extra SA steps 
    no_ds = True 

    for case_compare, runme in anom_outputs_dict.items(): 
        if runme: 
            case1, case2 = case_compare.split(" - ")
            # --- single case 1
            df_dict1 = pd.read_pickle(os.path.join(outdir, outputs_dict[case1], cdr_dfs))
            df_co21 = reorder_rows(df_params, df_dict1['co2_flx'])
            df_rockdiss1 = reorder_rows(df_params, df_dict1['rockdiss'])
            # --- single case 2
            df_dict2 = pd.read_pickle(os.path.join(outdir, outputs_dict[case2], cdr_dfs))
            df_co22 = reorder_rows(df_params, df_dict2['co2_flx'])
            df_rockdiss2 = reorder_rows(df_params, df_dict2['rockdiss'])
            # --- metrics
            Y_cdr_dif = df_co21['cdr_dif'].values - df_co22['cdr_dif'].values
            Y_cdr_adv = df_co21['cdr_adv'].values - df_co22['cdr_adv'].values
            Y_cdr_dif_component = df_co21['cdr_dif_component'].values - df_co22['cdr_dif_component'].values
            Y_cdr_resp_component = df_co21['cdr_resp_component'].values - df_co22['cdr_resp_component'].values
            Y_fraction_remaining_dissolved = df_rockdiss1['fraction_remaining_dissolved'].values - df_rockdiss2['fraction_remaining_dissolved'].values

            # --- SI results
            # [ cdr dif ]
            tmpds1 = Si_and_plot(
                problem,
                param_values,
                thismetric='cdr_dif',
                thiscase=case_compare,
                Yout=Y_cdr_dif,
                saveplots=True,
                plot_dict=plot_dict,
                conf_level=0.95,
                print_si_to_console = False,
            )

            # [ cdr adv ]
            tmpds2 = Si_and_plot(
                problem,
                param_values,
                thismetric='cdr_adv',
                thiscase=case_compare,
                Yout=Y_cdr_adv,
                saveplots=True,
                plot_dict=plot_dict,
                conf_level=0.95,
                print_si_to_console = False,
            )

            # [ cdr Y_cdr_dif_component ]
            tmpds3 = Si_and_plot(
                problem,
                param_values,
                thismetric='cdr_dif_component',
                thiscase=case_compare,
                Yout=Y_cdr_dif_component,
                saveplots=True,
                plot_dict=plot_dict,
                conf_level=0.95,
                print_si_to_console = False,
            )
            
            # [ cdr cdr_resp_component ]
            tmpds4 = Si_and_plot(
                problem,
                param_values,
                thismetric='cdr_resp_component',
                thiscase=case_compare,
                Yout=Y_cdr_resp_component,
                saveplots=True,
                plot_dict=plot_dict,
                conf_level=0.95,
                print_si_to_console = False,
            )

            # [ fraction_remaining_dissolved ]
            tmpds5 = Si_and_plot(
                problem,
                param_values,
                thismetric='fraction_remaining_dissolved',
                thiscase=case_compare,
                Yout=Y_fraction_remaining_dissolved,
                saveplots=True,
                plot_dict=plot_dict,
                conf_level=0.95,
                print_si_to_console = False,
            )

            # bring together 
            if no_ds:
                ds = xr.concat([tmpds1, tmpds2, tmpds3, tmpds4, tmpds5], dim="metric")
                ds = ds.expand_dims(case = [case_compare])
                no_ds = False
            else:
                tmpds = xr.concat([tmpds1, tmpds2, tmpds3, tmpds4, tmpds5], dim="metric")
                tmpds = tmpds.expand_dims(case = [case_compare])
                ds = xr.concat([ds, tmpds], dim='case')

    return ds


def Si_plot_batch(
    outputs_dict: dict,
    outdir: str,
    cdr_dfs: str,
    problem: dict,
    param_values: list,
    saveplots: bool,
    plot_dict: dict=None,
    conf_level: float=0.95,
    print_si_to_console: bool=False,
)-> xr.Dataset:
    '''
    Batch sensitivity analysis and produce / save plots.

    Parameters
    ----------
    outputs_dict : dict
        dictionary of outputs to consider. Keys are names of the output case 
        (e.g., `cc_1yr`) and values are the directory where outputs are 
        saved (e.g., `meanAnn_shortRun_SAFert_cc_morris_7_400_5yrInt_001`)
    outdir : str
        location of the output data (already processed)
    cdr_dfs : str
        name of the file with the dict of cdr data
    problem : dict
        SA problem dictionary
    param_values : list
        array of the parameter input values 
    saveplots : bool
        [True] to save figures (using info from the plot_dict)
    plot_dict : dict
        plotting information including `plot_types`; `plot_metrics`
        and `savepath`. 
    conf_level : float
        value to use for the confidence level in the SA 
    print_si_to_console : bool
        whether to print the SA results to console 
    
    Returns
    -------
    xr.Dataset
        Dataset with SA outputs for all cases / metrics etc. 
    '''
    # --- create param dataframe
    df_params = pd.DataFrame(param_values, columns = problem['names'])

    # --- get metrics for a single case 
    no_ds = True   # whether we've created the output dataset yet 

    for thiscase, thispath in outputs_dict.items():
        # --- single case 
        df_dict = pd.read_pickle(os.path.join(outdir, outputs_dict[thiscase], cdr_dfs))
        df_co2 = reorder_rows(df_params, df_dict['co2_flx'])
        df_rockdiss = reorder_rows(df_params, df_dict['rockdiss'])
        # --- metrics
        Y_cdr_dif = df_co2['cdr_dif'].values
        Y_cdr_adv = df_co2['cdr_adv'].values
        Y_cdr_dif_component = df_co2['cdr_dif_component'].values
        Y_cdr_resp_component = df_co2['cdr_resp_component'].values
        Y_fraction_remaining_dissolved = df_rockdiss['fraction_remaining_dissolved'].values

        # --- SI results
        # [ cdr dif ]
        tmpds1 = Si_and_plot(
            problem,
            param_values,
            thismetric='cdr_dif',
            thiscase=thiscase,
            Yout=Y_cdr_dif,
            saveplots=True,
            plot_dict=plot_dict,
            conf_level=0.95,
            print_si_to_console = False,
        )

        # [ cdr adv ]
        tmpds2 = Si_and_plot(
            problem,
            param_values,
            thismetric='cdr_adv',
            thiscase=thiscase,
            Yout=Y_cdr_adv,
            saveplots=True,
            plot_dict=plot_dict,
            conf_level=0.95,
            print_si_to_console = False,
        )

        # [ cdr Y_cdr_dif_component ]
        tmpds3 = Si_and_plot(
            problem,
            param_values,
            thismetric='cdr_dif_component',
            thiscase=thiscase,
            Yout=Y_cdr_dif_component,
            saveplots=True,
            plot_dict=plot_dict,
            conf_level=0.95,
            print_si_to_console = False,
        )
        
        # [ cdr cdr_resp_component ]
        tmpds4 = Si_and_plot(
            problem,
            param_values,
            thismetric='cdr_resp_component',
            thiscase=thiscase,
            Yout=Y_cdr_resp_component,
            saveplots=True,
            plot_dict=plot_dict,
            conf_level=0.95,
            print_si_to_console = False,
        )

        # [ fraction_remaining_dissolved ]
        tmpds5 = Si_and_plot(
            problem,
            param_values,
            thismetric='fraction_remaining_dissolved',
            thiscase=thiscase,
            Yout=Y_fraction_remaining_dissolved,
            saveplots=True,
            plot_dict=plot_dict,
            conf_level=0.95,
            print_si_to_console = False,
        )

        # bring together 
        if no_ds:
            ds = xr.concat([tmpds1, tmpds2, tmpds3, tmpds4, tmpds5], dim="metric")
            ds = ds.expand_dims(case = [thiscase])
            no_ds = False
        else:
            tmpds = xr.concat([tmpds1, tmpds2, tmpds3, tmpds4, tmpds5], dim="metric")
            tmpds = tmpds.expand_dims(case = [thiscase])
            ds = xr.concat([ds, tmpds], dim='case')

    return ds
