# %% 
# -------------------------------------
# 
# Run a batch of simulations based
# on a batch input file
# 
# -------------------------------------
import os
import subprocess
import sys

import pickle
# %% 

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# USER INPUTS
# 
batch_filename = "site+tstep+mineral_v0.pkl"
batchdir = "./batch_inputs"
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# %% 
# --- read in the batch file 
with open(os.path.join(batchdir, batch_filename), "rb") as f:  # "rb" = read binary
    batch_dict = pickle.load(f)
    n_runs = len(batch_dict)

# %%
# --- loop through the runs
for thisrun in range(n_runs):
    subprocess.run([sys.executable, 'singlerun.py', batch_filename, str(thisrun)])

# %%
