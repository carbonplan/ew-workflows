# -----------------------------------------------
# 
# create output variable dictionaries
# 
# -----------------------------------------------
# %% 
import os

import numpy as np
import pandas as pd

# %% 
# --- read in the mapping spreadsheet
mappath = "/home/tykukla/EWmodel/other"
mapfn = "smew-output-mapping.csv"

df = pd.read_csv(os.path.join(mappath, mapfn), skiprows=5)

# %% 
# --- split outdict and omit unnecessary columns
df_outdictsave = df[(df['dict'] == 'outdict') & (df['save'] == "Y")].reset_index(drop=True)
df_outdictsave = df_outdictsave.drop(columns=['dict', 'Unnamed: 3'])

# %%
# --- split rundict and omit unnecessary columns
df_rundictsave = df[(df['dict'] == 'rundict') & (df['save'] == "Y")].reset_index(drop=True)
df_rundictsave = df_rundictsave.drop(columns=['dict', 'Unnamed: 3'])

# %%
# -------------------------------------------------
# [ save the dataframes ]
df_outdictsave.to_csv(os.path.join(mappath, "outdict_map.csv"), index=False)
df_rundictsave.to_csv(os.path.join(mappath, "rundict_map.csv"), index=False)
# %%
