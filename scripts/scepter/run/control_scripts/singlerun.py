# -------------------------------------------
# 
# Script to replace the "singlerun.sh" script
# Reads in system args from the argo submit
# command and runs the appropriate python 
# run script  
# 
# -------------------------------------------
import os
import subprocess
import sys

import pandas as pd


# --- set the variables from system args
# name of run script
runscript = os.path.join(sys.argv[1], "rock_app_singleRun-szn.py")
# model directory
modeldir=f"{sys.argv[1]}"
# output directory
outdir = os.path.join(modeldir, "scepter_output/")  # (note: the trailing '/' is critical because of lazy coding on the other side)
# default dictionary
default_dict=f"{sys.argv[5]}"  # see SCEPTER/defaults/dict_singlerun.py


# --- COLLECT INPUTS ---
# input csv
input_dir=f'{sys.argv[2]}'
input_name=f'{sys.argv[3]}'
input_index=int(sys.argv[4])
# subtract 1 from the index so first element is zero
input_index -= 1 

# read in the batch csv
df_batch = pd.read_csv(os.path.join(input_dir, input_name))
# pull out just the row of interest 
row = df_batch.iloc[input_index]

# construct the first part of the command 
python_cmd_start = f'python3 {runscript} --modeldir {modeldir} --outdir {outdir} --default_dict {default_dict}'
# add all the other variables whose default values we want to overwrite
python_cmd_end = " ".join([f'--{col} {row[col]}' for col in df_batch.columns])
python_cmd = f'{python_cmd_start} {python_cmd_end}'

# --- run the command
subprocess.run(python_cmd, shell=True)


# --- TROUBLESHOOT 
# print(python_cmd)