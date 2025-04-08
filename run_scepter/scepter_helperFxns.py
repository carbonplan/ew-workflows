# ----------------------------------------------------
#
# Functions used to integrate aglime-swap-cdr repo
# with SCEPTER simulations and storage in aws buckets
#
# T Kukla (CarbonPlan, 2024)
#
# ----------------------------------------------------
# %%
import os
import re
import shutil
import subprocess
from scipy.integrate import cumulative_trapezoid

import pandas as pd

# %%
# --------------------------------------------------------------------------
# --- RUN SETUP FUNCTIONS:
#     these are used to help build the model directory and set
#     the model inputs


def parse_arguments(args: list) -> dict:
    """
    Turn system arguments that over-write the default values
    into a dictionary that can be compared to the default dictionary

    Parameters
    ----------
    args : list
        result of `sys.argv`
        system args should have format --dictKey1 dictValue1 --dictKey2 dictValue2 ...


    Returns
    -------
    dict
        keywords are variable names, values are variable values
    """

    parsed_args = {}  # create empty dict
    i = 1  # start from index 1 to skip the script name (which lies at sys.argv[0])
    while i < len(args):
        if args[i].startswith("--"):  # then it's a keyword
            key = args[i][2:]  # Remove '--' prefix
            if i + 1 < len(args) and not args[i + 1].startswith("--"):
                value = args[i + 1]
                parsed_args[key] = value
                i += 1  # Skip the next item as it's the value for the current key
            else:
                parsed_args[key] = None  # If no value provided, set to None
        i += 1
    return parsed_args


# Function to set global variables from defaults / system args
def set_vars(default_args: dict, system_args: dict) -> dict:
    """
    Set global variables using defaults first, then over-writing
    with the system args.
    (NOTE: setting global vars doesn't work when they're needed
    outside of this helper function module.. add in main script:
    for key, value in combined_dict.items():
        globals()[key] = value
    to address this limitation.

    Parameters
    ----------
    default_args : dict
        dictionary of default arguments, generally from
        `getattr(defaults.dict_singlerun, import_dict)`
    system_args : dict
        dictionary of values to override the defaults,
        generally the output of parse_arguments()


    Returns
    -------
    dict
        the final variable values used in the run. The
        function also sets the values as global variables
        (no return necessary for that step)
    """

    # define pattern for identifying floats in sys.args
    float_pattern = r"^[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?$"
    # dict to save
    save_vars = {}
    for key, value in default_args.items():
        if key in system_args:
            float_test1 = (
                re.match(float_pattern, system_args[key]) is not None
            )  # (captures all cases but "2.")
            float_test2 = (
                system_args[key].replace(".", "", 1).isdigit()
            )  # (misses negatives but gets others including "2.")
            bool_test_true = (
                system_args[key] == "True"
            )  # check if we should turn string True into boolean
            bool_test_false = (
                system_args[key] == "False"
            )  # check if we should turn string False into boolean
            if bool_test_true:
                save_vars[key] = True
                globals()[key] = True
            elif bool_test_false:
                save_vars[key] = False
                globals()[key] = False
            elif float_test1 or float_test2:  # check is sys arg is a float
                save_vars[key] = float(system_args[key])
                globals()[key] = float(system_args[key])
            else:
                save_vars[key] = system_args[key]
                globals()[key] = system_args[key]
        else:
            save_vars[key] = value
            globals()[key] = value
    return save_vars


def save_dict_to_text_file(dictionary: dict, filename: str, delimiter: str = "\t"):
    """
    Write the final dictionary to a text file in a given destination,
    generally in the model run directory

    Parameters
    ----------
    dictionary : dict
        dictionary where keyword is the variable and value is
        the value. the dictionary is generally the output of `set_vars()`
    filename : str
        name of the variable file to save, appended to the path
        ('/my/file/goes/here/filename.res')


    Returns
    -------
    """
    with open(filename, "w") as file:
        file.write("*** variables set by dictionary and system args\n")
        file.write("    (note not all vars are used!!)\n")
        for key, value in dictionary.items():
            file.write(f"{key}{delimiter}{value}\n")


def copy_files(src_dir: str, dst_dir: str):
    """
    Copy all files from a source directory to a destination directory
    (assumes neither are aws)

    Parameters
    ----------
    src_dir : str
        location of the source directory (often the spinup run)
        ('/my/source_dir')
    dst_dir : str
        location of the destination directory (often the case run)
        the destination directory must already exist
        ('/my/dst_dir')


    Returns
    -------
    """
    # directory must exist
    if not os.path.exists(dst_dir):
        print("Cannot find " + dst_dir)

    # iterate over source files
    for filename in os.listdir(src_dir):
        src_file = os.path.join(src_dir, filename)
        dst_file = os.path.join(dst_dir, filename)

        # copy them over (skipping any directories)
        if os.path.isfile(src_file):
            shutil.copy2(src_file, dst_file)
            # print(f"Copied {src_file} to {dst_file}")


def generate_timesteps(total_duration: float, timestep: float) -> list:
    """
    Create a list of time steps given a total duration for
    a scepter run and a timestep. (for use in multi-run only).

    Parameters
    ----------
    total_duration : float
        [years] value denoting the total duration of the
        SCEPTER simulation
    timestep : float


    Returns
    -------
    list
        list of start times for SCEPTER simulations
        for ex: [0,2,4,6,...]
    """
    timesteparr = []
    current_time = 0
    while current_time <= (total_duration - timestep):
        timesteparr.append(current_time)
        current_time += timestep
    return timesteparr


def write_iter_file_with_marker(
    values: list,
    current_index: int,
    filename: str,
):
    """
    Write a file showing all timesteps and an arrow
    toward the specific timestep handled by the given
    scepter run (for use in multi-run only)

    Parameters
    ----------
    values : list
        list of timesteps, generally the output of the
        generate_timesteps function
    current_index : int
        the run number in the sequence of simulations
        (set by the 'counter' in the multi-run loop)
    filename : str
        the path to the file that notes which timestep is
        used. default in the scepter multi-run script is
        "multiyear-iter.res"


    Returns
    -------
    """
    with open(filename, "w") as file:
        file.write("year \n")
        for i, value in enumerate(values):
            if i == current_index:
                file.write(f"{value} <--- \n")
            else:
                file.write(f"{value}\n")


def update_clim(
    inputfile: str,
    outputfile: str,
    timezero: float,
):
    """
    SCEPTER climatology always starts at year zero. When you run a multi-year
    run you need to modify t=0 depending on the timestep (e.g., for the run
    spanning years 2-4, climatology year 2 must become year zero). this function
    updates the climatology to match the multi-run index and saves the result
    in the new run's directory. (for use in multi-run only)

    Parameters
    ----------
    inputfile : str
        path to the source climatology file (in the old timestep dir)
    outputfile : src
        location where we'll write the output climatology file (the new timestep dir)
    timezero : float
        the time value that gets set to zero in this iteration


    Returns
    -------
    """
    # open the input and output files
    with open(inputfile, "r") as f_in, open(outputfile, "w") as f_out:
        # skip the header line and copy it to new file
        header = next(f_in)
        f_out.write(header)  # write to output
        # read the file line by line
        for line in f_in:
            # split each line into year and climate value
            year, clim = line.strip().split("\t")  # Adjust the delimiter as needed

            # convert year to an integer and subtract 5
            year = float(year) - timezero

            # check if the adjusted year is greater than or equal to tstep
            if year >= 0:
                formatted_yr = "{:.7f}".format(float(year))
                formatted_clim = "{:.7f}".format(float(clim))
                # write the adjusted year and temperature to the output file
                # print(f"{formatted_yr}\t{formatted_clim}\n")
                f_out.write(
                    f"{formatted_yr}\t{formatted_clim}\n"
                )  # Adjust the delimiter as needed


def remove_duplicates(input_file: str):
    """
    Read the input file, identify and delete duplicate lines.
    Written for use in the solid input files, where adding
    solid phases can lead to duplicating the defaults.

    Parameters
    ----------
    input_file : str
        location and name of the input file, like
        `/my/input/file/is/here.in`


    Returns
    -------
    """
    with open(input_file, "r") as f:
        lines = f.readlines()
        header = lines[0]  # Save the header
        lines_without_tabs = [
            line.replace("\t", "") for line in lines
        ]  # remove tabs before comparing
        # add '\n' if a line doesn't have it (just for comparison sake)
        lines_without_tabs2 = [
            element if element.endswith("\n") else element + "\n"
            for element in lines_without_tabs
        ]
        unique_lines = set(
            lines_without_tabs2[1:]
        )  # remove duplicates from mineral names

    # check if the last entry ends with "\n" and remove it if needed
    last_entry = list(unique_lines)[-1]

    # write lines to output file, keeping the header at the top
    with open(input_file, "w") as f:
        f.write(header)  # write the header first
        # write unique mineral names
        for line in unique_lines:
            if line != last_entry:
                if not line.endswith(
                    "\n"
                ):  # if the last item got moved around, it won't have a newline indicator so we must add it..
                    line += "\n"
                f.write(line)
            else:
                f.write(line.rstrip("\n"))  # don't add newline for the last entry


# --------------------------------------------------------------------------
# --- DIAGNOSTIC FUNCTIONS
#     these are run after a SCEPTER run to check that the model completed
#     and outputs were saved


def preprocess_txt(
    file_path: str,
) -> pd.DataFrame:
    """
    Preprocess a SCEPTER output .txt file to get into
    a pandas dataframe (used in duration_check fxn)

    Parameters
    ----------
    file_path : str
        path to the file we want to read in and convert


    Returns
    -------
    pd.DataFrame
        pandas DataFrame with data in the .txt file
    """
    data = []  # Initialize a list to store the processed data

    # Initialize a flag to determine if we are reading the header
    is_header = True

    # Read the file line by line and process the data
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()  # Remove leading/trailing whitespace
            if is_header:
                # Split the first line into column names
                column_names = re.split(r"\s+", line)
                is_header = False
            else:
                # Split the other lines into data values
                values = re.split(r"\s+", line)
                data.append(values)

    # Create a DataFrame with the processed data and set column names
    df = pd.DataFrame(data, columns=column_names)
    # get df to numeric if it's not already
    df = df.map(pd.to_numeric)
    # return
    return df


def duration_check(
    directory: str,
    target_duration: float,
    file_to_check: str = "flx_aq-ca.txt",
    duration_tol: float = 0.1,
) -> tuple[bool, list]:
    """
    Compare the target duration of the run (in years) to the
    duration of saved output to see if the run made it
    to completion

    Parameters
    ----------
    directory : str
        path to the run directory that we want to check
    empty_dirs : str
        output of empty_dir_check to make sure the flx dir exists
    target_duration : float
        [years] duration we set the model to run
    file_to_check : str
        name of the file we want to compare the duration to
    duration_tol : float
        [percent] tolerance for the difference between target
        and actual duration (target is reached if percent difference
        is less than duration_tol


    Returns
    -------
    tuple(bool, list)
        bool: [True, False] denoting whether the model run duration is
               within tolerance of the target duration
        list: list of format [target_duration, actual_duration]
    """

    # read in data
    dat_path = os.path.join(directory, "flx", file_to_check)
    if os.path.isfile(dat_path):
        df = preprocess_txt(dat_path)

        # compare to target
        model_duration = df["time"].values.max()
        diff_percent = abs(((model_duration - target_duration) / target_duration) * 100)
        diff_result = True if diff_percent < duration_tol else False

        # make output list
        dur_list = [target_duration, model_duration]

    else:
        diff_result = f"could not complete check; file {file_to_check} not found"
        dur_list = ["NA", "NA"]

    return diff_result, dur_list


def empty_file_check(
    directory: str,
    omit_saveSuff: bool = True,
) -> tuple[bool, list]:
    """
    Check that all files in all subdirectories (e.g., "prof", "flx")
    are not empty. Empty files are appended to a list which is returned.
    Does not check whether files in 'directory' are empty (these are
    used to run scepter, we only want to check the outputs)

    Parameters
    ----------
    directory : str
        path to the run directory that we want to check
    omit_saveSuff : bool
        True or False - whether to omit files with suffix .save
        (we might be okay / expect some of these to be empty)


    Returns
    -------
    bool
        boolean (True, False) denoting whether no empty files
        were found. (True means no empty files were found)
    list
        list of form ['empty_filename1', 'empty_filename2', ...]
    """
    # check that dir exists and is a dir
    if not os.path.exists(directory):
        raise FileNotFoundError(f"The directory '{directory}' does not exist.")
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"The path '{directory}' is not a directory.")

    empty_files = []  # create empty list to hold results

    # walk through the directory and its subdirectories
    for root, dirs, files in os.walk(directory):
        # -- troubleshoot
        # print(f"Current root: {root}")
        # print(f"Subdirectories: {dirs}")
        # print(f"Files: {files}")

        # skip files in the main directory, only check subdirectories
        if root == directory:
            continue

        for file in files:
            if omit_saveSuff:  # then check if the suffix is ".save"
                # get the file extension
                _, fn_ext = os.path.splitext(file)
                # skip this file if suffix is .save
                if fn_ext.lower() == ".save":
                    # print('found a .save')
                    continue
            file_path = os.path.join(root, file)
            # check if the file is empty (ignoring spaces)
            with open(file_path, "r") as f:
                if (file_path.endswith(".pkl") or file_path.endswith(".nc")):
                    print("ignoring .pkl and .nc files")
                else:
                    # print(file_path)
                    content = f.read().strip()
                    if not content:
                        empty_files.append(file_path)

    # set up return
    if len(empty_files) == 0:
        no_empties = True
    else:
        no_empties = False

    return no_empties, empty_files


def empty_dir_check(
    directory: str,
    omit_ipynb: bool = True,
) -> tuple[bool, list, bool, bool]:
    """
    Check that all subdirectories (e.g., "prof", "flx")
    are not empty. Empty dirs are appended to a list which is returned.


    Parameters
    ----------
    directory : str
        path to the run directory that we want to check
    omit_ipynb : bool
        True will omit .ipynb_checkpoints file when searching for directories


    Returns
    -------
    bool1
        boolean (True, False) denoting whether any empty dirs
        were found. (True means empty files were found)
    list
        list of form ['empty_dirname1', 'empty_dirname2', ...]
    bool2
        boolean (True, False) denoting whether the "prof" subdir exists
    bool3
        boolean (True, False) denoting whether the "flx" subdir exists
    """
    # check that dir exists and is a dir
    if not os.path.exists(directory):
        raise FileNotFoundError(f"The directory '{directory}' does not exist.")
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"The path '{directory}' is not a directory.")

    # check for prof and flx to exist
    prof_path = os.path.join(directory, "prof")
    flx_path = os.path.join(directory, "flx")
    prof_exists = True if os.path.exists(prof_path) else False
    flx_exists = True if os.path.exists(flx_path) else False

    empty_subdirs = []  # create empty list to hold results

    # walk through the directory and its subdirectories
    for root, dirs, files in os.walk(directory):
        # -- troubleshoot
        # print(f"Current root: {root}")
        # print(f"Subdirectories: {dirs}")
        # print(f"Files: {files}")

        # skip files in the main directory, only check subdirectories
        if root == directory:
            continue

        # omit ipynb if asked
        if omit_ipynb and (".ipynb" in root):
            # print('found ipynb')
            continue

        # check if the current directory is empty
        # list all subdirectories in the main directory
        for entry in os.listdir(directory):
            entry_path = os.path.join(directory, entry)
            if os.path.isdir(entry_path):  # check if it is a subdirectory
                contents = [
                    item
                    for item in os.listdir(entry_path)
                    if item != ".ipynb_checkpoints"
                ]
                if not contents:  # check if the subdirectory is empty
                    if (
                        entry_path not in empty_subdirs
                    ):  # only append if we don't have it already (we had some duplication issues)
                        empty_subdirs.append(entry_path)

        # if not dirs and not files:
        #     empty_subdirs.append(root)

    # set up return
    if len(empty_subdirs) == 0:
        no_empties = True
    else:
        no_empties = False

    return no_empties, empty_subdirs, prof_exists, flx_exists


def run_complete_check(
    runname_field: str,
    runname_lab: str,
    outdir: str,
    target_duration: float = 0.00001,
    include_duration_check: bool = True,
    omit_saveSuff: bool = True,
    omit_ipynb: bool = True,
):
    """
    Run basic diagnostics and save the results in a
    resource file in the run directory. Diagnostics are
    (1) ensuring output dirs are populated;
    (2) ensuring output files are not empty;
    (3) ensuring model reached target duration;
    (4) TK..

    Parameters
    ----------
    runname_lab : str
        name of the lab run, also the name of the run directory which stores
        the lab results
    runname_field : str
        name of the field run, also the name of the run directory which stores
        the field results
    outdir : str
        path to local directory that holds the model run directories to be moved
    target_duration : float
        [years] the duration we wanted the model to run
    include_duration_check : bool
        True or False - whether to check that model duration equals target duration
        (may want to be False for spinup diagnostics if target duration isn't easy to access)
    omit_saveSuff : bool
        True or False - whether to omit files with suffix .save
        (we might be okay / expect some of these to be empty)
        (this is used in the empty_file_check fxn)
    omit_ipynb : bool
        True or False - whether to omit ipynb_checkpoints when
        searching for empty dirs
        (this is used in the empty_dir_check fxn)


    Returns
    -------
    """

    # loop through dirs
    for runname in [runname_field, runname_lab]:
        if runname == runname_lab:
            include_duration_check_internal = (
                False  # don't do duration checks on lab runs
            )
        else:
            include_duration_check_internal = include_duration_check

        # set destination
        dst = os.path.join(outdir, runname)
        # check that prof / flx dirs exist
        none_emptyDir_tf, empty_dirs, prof_exists, flx_exists = empty_dir_check(
            dst, omit_ipynb
        )
        # check that files aren't empty
        none_emptyFn_tf, empty_files = empty_file_check(dst, omit_saveSuff)
        # check that model reached target duration
        if include_duration_check_internal:
            dur_tf, dur_list = duration_check(dst, target_duration)
        else:
            dur_tf = "NA; check deliberately not run"
            dur_list = ["NA", "NA"]
        # write the check results file
        fn_checkres = os.path.join(dst, "check_results.res")
        with open(fn_checkres, "w") as file:
            file.write("*** results of postprocess diagnostic checks\n")
            file.write(f"check1 -- Profile subdir exists: \t{prof_exists}\n")
            file.write(f"check2 -- Flx subdir exists: \t{flx_exists}\n")
            file.write(
                f"check3 -- Prof and flx dirs are not empty: \t{none_emptyDir_tf}\n"
            )
            file.write(f"check4 -- No output files are empty: \t{none_emptyFn_tf}\n")
            file.write(f"check5 -- Model reached target duration: \t{dur_tf}\n")

        # write the check log file
        fn_checkres = os.path.join(dst, "check_logs.res")
        with open(fn_checkres, "w") as file:
            file.write("*** log of postprocess diagnostic checks\n")
            file.write("\n")
            file.write("check1 -- no output (see check_results.res)\n")
            file.write("\n\n")
            file.write("check2 -- no output (see check_results.res)\n")
            file.write("\n\n")
            file.write("check3 -- list of empty output directories\n")
            for empty_dir in empty_dirs:
                file.write(f"{empty_dir}\n")
            file.write("\n\n")
            file.write("check4 -- list of empty output files\n")
            for empty_file in empty_files:
                file.write(f"{empty_file}\n")
            file.write("\n\n")
            file.write("check5 -- target versus model duration (years) \n")
            file.write(f"target: {str(dur_list[0])}\t model: {str(dur_list[1])}")

        # if we got here, the run completed, so mint that empty file
        fn_completed = os.path.join(dst, "completed.res")
        with open(fn_completed, "w") as f:
            pass  # pass does nothing, creating an empty file


# modify the flx/dust.txt file to get total dust flux over time 
def dustflx_calc(
    outdir: str,
    runname_field: str, 
    fdust1: float,
    fdust2: float,
    dustsp: str,
    dustsp_2nd: str=None):
    """
    SCEPTER's default dust file (*/flx/dust.txt) shows time and the relative amount
    of dust that gets applied. We use the dust flux to compute the timeseries
    of dust application and the integrated dust application. Nothing gets returned,
    but the updated dust.txt file gets saved. 

    Parameters
    ----------
    outdir : str
        output directory for SCEPTER results (ex: "home/name/SCEPTER/scepter_output")
    runname_field : str
        name of the directory for the SCEPTER field run (for now, we're not applying 
        any changes to the lab dust fluxes)
    fdust1 : float
        [g m2 yr] amount of primary dust applied each year 
    fdust2 : float
        [g m2 yr] amount of secondary dust applied each year 
    dustsp : str
        name of the primary dust species (e.g., "cc" or "gbas")
    dustsp_2nd : str
        name of the secondary dust species (e.g., "amnt")

    Results
    -------
    """
    # read in the dust flux dataframe
    file_path = os.path.join(outdir, runname_field, 'flx', 'dust.txt')
    df = preprocess_txt(file_path)

    # make sure dust fluxes are numbers
    if not isinstance(fdust1, (int, float)):
        try:
            fdust1 = float(fdust1)
        except:
            fdust1 = 0
    if not isinstance(fdust2, (int, float)):
        try:
            fdust2 = float(fdust2)
        except:
            fdust2 = 0
    
    # add a column for the dust flux
    df['dust1_STATIC'] = fdust1
    df['dust2_STATIC'] = fdust2

    # compute dust timeseries 
    df['dust1_g_m2_yr'] = df['dust(relative_to_average)'] * df['dust1_STATIC']
    df['dust2_g_m2_yr'] = df['dust(relative_to_average)'] * df['dust2_STATIC']

    # compute integrated dust timeseries
    df['int_dust1_g_m2_yr'] = cumulative_trapezoid(df['dust1_g_m2_yr'], df['time'], initial=0)
    df['int_dust2_g_m2_yr'] = cumulative_trapezoid(df['dust2_g_m2_yr'], df['time'], initial=0)

    # add dust species info
    df['dustsp1'] = dustsp
    df['dustsp2'] = dustsp_2nd
    
    # save the result
    df.to_csv(  # default is mode='w' which will overwrite the existing file (that's fine because we've merged it with the new data)
        file_path, index=None, sep="\t"
        )  


# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# --- AWS FUNCTIONS:
#     these are used to move/copy files to or from an aws bucket


def to_aws(
    aws_save: str,
    aws_bucket: str,
    outdir: str,
    runname_lab: str | list[str],
    runname_field: str | list[str],
):
    """
    Copy or move files in runname_lab and runname_field
    directories from local to an aws bucket or,
    depending on value of aws_save, do nothing.

    Parameters
    ----------
    aws_save : str
        ['move', 'copy', 'None']. 'move' will empty the local
        directories and move them to aws. 'copy' will retain the
        local copy and add to aws. 'None', or any other input,
        results in no action being taken.
    aws_bucket : str
        path to the aws bucket where content will be moved ('s3://my/bucket/is/here').
        Only used if aws_save == 'move' or 'copy'
    outdir : str
        path to local directory that holds the model run directories to be moved
    runname_lab : str | list[str]
        name of the lab run, also the name of the run directory which stores
        the lab results. If it's a string, then it gets converted to a one-element
        list of strings in the function and we loop through that list. 
    runname_field : str | list[str]
        name of the field run, also the name of the run directory which stores
        the field results. If it's a string, then it gets converted to a one-element
        list of strings in the function and we loop through that list. 

    Returns
    -------
    """
    # --- check if runnames are strings and, if so, convert to lists
    if isinstance(runname_field, str):
        runname_field = [runname_field]
    if isinstance(runname_lab, str):
        runname_lab = [runname_lab]

    # --- check if runname_field and runname_lab are not the same length
    if len(runname_field) != len(runname_lab):
        raise ValueError("Number of runname_field inputs != number of runname_lab inputs in `shf.to_aws`. Can't loop through each run properly since each run must have a field and lab component")
    else:
        nruns = len(runname_field)  # get the total number of runs to loop through

    
    # loop through runs
    for rundx in range(nruns):
        tmp_runname_field = runname_field[rundx]
        tmp_runname_lab = runname_lab[rundx]    
        
        # loop through dirs
        for runname in [tmp_runname_field, tmp_runname_lab]:
            src = os.path.join(outdir, runname)
            dst_aws = os.path.join(aws_bucket)
            cmd_activate = "conda run -n cdr-scepter1p0_env"  # activate the virtual environment that has s5cmd
            # check if the destination directory exists and remove it if it does
            # if os.path.exists(dst_aws):
            #     shutil.rmtree(dst_aws)

            if aws_save == "move":
                cmd_run = "s5cmd mv " + src + " " + dst_aws
                result1 = subprocess.run(cmd_activate, shell=True, check=True)
                result2 = subprocess.run(cmd_run, shell=True, check=True)
                outdir_final = dst_aws
                # ... remove empty sub directories
                for dirpath, dirnames, filenames in os.walk(src, topdown=False):
                    for dirname in dirnames:
                        dir_to_check = os.path.join(dirpath, dirname)
                        if not os.listdir(dir_to_check):  # check if the directory is empty
                            os.rmdir(dir_to_check)
                            print(f"Removed empty directory: {dir_to_check}")
                # ... remove the src directory if it's empty
                if not os.listdir(src):  # Check if src is empty
                    os.rmdir(src)

            elif aws_save == "copy":
                cmd_run = "s5cmd cp " + src + " " + dst_aws
                result1 = subprocess.run(cmd_activate, shell=True, check=True)
                result2 = subprocess.run(cmd_run, shell=True, check=True)
                outdir_final = dst_aws
            else:
                outdir_final = outdir
        
    # return the new directory for postproc purposes
    return outdir_final


def aws_copy_to_local(
    aws_bucket: str,
    aws_src: tuple[str, str],
    outdir: str,
):
    """
    Copy files from aws to local. Each directory name in
    aws_src tuple will be its own directory in the local
    outdir. If outdir already has this directory, similar file
    names are over-written.

    Parameters
    ----------
    aws_bucket : str
        path to the aws bucket where content will be moved ('s3://my/bucket/is/here').
    aws_src : tuple[str, str]
        tuple including the names of the directories within aws_bucket to
        be copied back to local (e.g., ('runname_lab', 'runname_field'))
    outdir : str
        path to the model output directory where directories with aws_src
        names will be stored.


    Returns
    -------
    """
    # command to activate virtuan env with s5cmd
    cmd_activate = "conda run -n cdr-scepter1p0_env"

    # loop through source dirs
    for src in aws_src:
        src_tmp = os.path.join(aws_bucket, src, "*")
        dst_tmp = os.path.join(outdir, src)
        cmd_run = "s5cmd cp " + src_tmp + " " + dst_tmp
        # run the commands
        result1 = subprocess.run(cmd_activate, shell=True, check=True)
        result2 = subprocess.run(cmd_run, shell=True, check=True)


# --------------------------------------------------------------------------
