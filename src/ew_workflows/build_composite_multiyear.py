# --------------------------------
# functions to build composite
# output from multiple model
# iterations
#
# --------------------------------
# %%

# script to compile output into one dataset
import os
import re
import shutil

import pandas as pd


# --- function to find relevant subdirectories
def find_subdirs_1level(root_dir, target_string, field_only=False):
    subdirectories = [
        d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))
    ]

    matching_subdirectories = []
    for dirname in subdirectories:
        if target_string in dirname:
            matching_subdirectories.append(os.path.join(root_dir, dirname))
    if field_only:
        # remove elements with "lab" in the directory name
        filtered_directories = [
            directory for directory in matching_subdirectories if "lab" not in directory
        ]
        # print(filtered_directories)
    else:
        filtered_directories = matching_subdirectories
    return filtered_directories


# --- function to preprocess .txt files for consistent delimiters
def preprocess_txt(file_path):
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
    # return
    return df


# --- function to build composite output directory
def build_composite(basename, outdir):
    # ... first get all the outdir paths
    # (these will be passed to the function in the main script!)
    # outdir = "/home/tykukla/SCEPTER/scepter_output"
    # -- we're having an issue where a basename like 
    #    mydir_10 returns "mydir_10_startyear*" but it also 
    #    returns "mydir_100_startyear*". So we add an under-
    #    score to prevent this but only for finding alldirs
    runname_base_underscore = basename + "_"
    runname_base = basename
    alldirs = sorted(find_subdirs_1level(outdir, runname_base_underscore))

    # remove paths with the word "composite"
    alldirs = [path for path in alldirs if "composite" not in os.path.basename(path)]

    # split into field and lab
    field_paths = [path for path in alldirs if path.endswith("_field")]
    lab_paths = [path for path in alldirs if path.endswith("_lab")]

    # ... create the new output directory (for field and lab)
    newname_field = runname_base + "_composite_field"
    newname_lab = runname_base + "_composite_lab"
    dst_main_field = os.path.join(outdir, newname_field)
    dst_main_lab = os.path.join(outdir, newname_lab)
    os.makedirs(dst_main_field, exist_ok=True)
    os.makedirs(dst_main_lab, exist_ok=True)

    # counters for the loop
    lab_counter = field_counter = -1

    # --- loop through each file
    for src_tmp in alldirs:
        # ... get start year and iteration index
        # expression pattern to extract startyear and iter values
        # pattern = r"startyear-(\d+)_iter-(\d+)"
        # pattern = r"startyear-(\d+)p\d+_iter-(\d+)"
        pattern = r"startyear-(\d+)(p\d+)?_iter-(\d+)"

        # Search for the pattern in the directory name
        match = re.search(pattern, src_tmp)
        if (
            match.group(2) is None
        ):  # if there's no decimal value, just use the first value
            startyear_tmp = match.group(1)
        else:  # otherwise add decimal and second value
            if "p" in match.group(2):
                startyear_tmp = match.group(1) + match.group(2).replace("p", ".")
            else:
                startyear_tmp = match.group(1) + "." + match.group(2)
        iter_tmp = match.group(3)

        # ... move all model files into a subdirectory within the composite outdir
        #     (note, results files are moved / composited later)
        # check if lab or field run
        if src_tmp in field_paths:
            runtype = "field"
            field_counter += 1
            dst_main = dst_main_field
        elif src_tmp in lab_paths:
            runtype = "lab"
            lab_counter += 1
            dst_main = dst_main_lab

        # make the subdir
        dst_subdir_tmp = os.path.join(dst_main, os.path.basename(src_tmp))
        os.makedirs(dst_subdir_tmp, exist_ok=True)

        # loop through all items in the source directory
        for item in os.listdir(src_tmp):
            source_item = os.path.join(src_tmp, item)
            # check if the item is a file (not a directory)
            if os.path.isfile(source_item):
                # move to dest dir
                shutil.copy(source_item, dst_subdir_tmp)

        # ----------------------------------------------------------
        # ... move the profile data
        # make 'prof' dir if it doesn't exist
        dst_main_prof = os.path.join(dst_main, "prof")
        os.makedirs(dst_main_prof, exist_ok=True)

        # get all the filenames
        src_tmp_prof = os.path.join(src_tmp, "prof")
        prof_fns_tmp = os.listdir(src_tmp_prof)

        # exclude filenames with "save" or "restart" in them
        prof_fns_tmp = [f for f in prof_fns_tmp if "save" not in f]
        prof_fns_tmp = [f for f in prof_fns_tmp if "restart" not in f]

        # exclude any .ipynb files
        prof_fns_tmp = [f for f in prof_fns_tmp if ".ipynb" not in f]

        # update time for relevant files and save
        for fn in prof_fns_tmp:
            # read in
            thispath = os.path.join(src_tmp_prof, fn)
            dfsrc = preprocess_txt(thispath)

            # update time if its a column
            if "time" in dfsrc.columns:
                dfsrc["time"] = pd.to_numeric(dfsrc["time"]) + float(startyear_tmp)

            # add the iter index to each filename
            fn_save = fn.replace(".txt", "-" + str(iter_tmp) + ".txt")
            # save file in destination
            savedst_tmp = os.path.join(dst_main_prof, fn_save)
            dfsrc.to_csv(savedst_tmp, index=None, sep="\t")

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # ... move the flx data
        # make 'flx' dir if it doesn't exist
        dst_main_flx = os.path.join(dst_main, "flx")
        os.makedirs(dst_main_flx, exist_ok=True)

        # get all the filenames
        src_tmp_flx = os.path.join(src_tmp, "flx")
        flx_fns_tmp = os.listdir(src_tmp_flx)

        # exclude filenames with "save" or "restart" in them
        flx_fns_tmp = [f for f in flx_fns_tmp if "save" not in f]
        flx_fns_tmp = [f for f in flx_fns_tmp if "restart" not in f]

        # exclude any .ipynb files
        flx_fns_tmp = [f for f in flx_fns_tmp if ".ipynb" not in f]

        # update time for relevant files and save
        for fn in flx_fns_tmp:
            # read in
            thispath = os.path.join(src_tmp_flx, fn)
            dfsrc = preprocess_txt(thispath)

            # update time if its a column
            if "time" in dfsrc.columns:
                dfsrc["time"] = pd.to_numeric(dfsrc["time"]) + float(startyear_tmp)
                # order values by time
                dfsrc = dfsrc.sort_values(by="time")

            # ... save file in destination
            # check if file exists
            savedst_tmp = os.path.join(dst_main_flx, fn)
            filecheck = os.path.isfile(savedst_tmp)
            if filecheck:  # then read in the existing df and append the new one 
                # read the existing file into a dataframe
                existing_df = preprocess_txt(savedst_tmp)
                # join the existing and source dfs together by column
                # (join='outer' means a column that only exists in df2 will be kept, with nans in df1)
                #  and axis=0 means it's concatenated row-wise, not column-wise)
                new_df = pd.concat([existing_df, dfsrc], axis=0, join='outer')
                new_df.to_csv(  # default is mode='w' which will overwrite the existing file (that's fine because we've merged it with the new data)
                    savedst_tmp, index=None, sep="\t"
                )  
            else:  # if the file doesn't exist, save and include the header
                dfsrc.to_csv(
                    savedst_tmp, index=None, sep="\t", mode="a"
                )  # mode = "a" will append to end of file if exists
        # ----------------------------------------------------------
    # return the new output dirs
    return newname_field, newname_lab # dst_main_field, dst_main_lab

# # save result or append if one exists
#     if filecheck:
#         dfsrc.to_csv(dst, header=None, index=None, sep='\t', mode='a')
#     else:
#         # write the first line
#         with open(src) as f:
#             lines = f.readlines() #read
#         with open(dst, "w") as f:
#             f.writelines(lines[0]) #write back
#         # append result
#         dfsrc.to_csv(dst, index=None, sep='\t', mode='a')   # mode = "a" will append to end of file if exists


# %%
