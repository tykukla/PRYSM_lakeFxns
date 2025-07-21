# --------------------------------------------------
# 
# run batch simulations of the lake driver psm
# 
# T Kukla, 2025
# --------------------------------------------------
# %% 
import re
import os
import gc
import shutil
import subprocess  # for running fortran model
import numpy as np
import pandas as pd
import xarray as xr
from math import pi, sqrt, exp
import matplotlib
import matplotlib.pyplot as plt
import warnings

# helper functions (T Kukla)
import helper_functions as lhf

# sensor models
# import sensor_gdgt as gdgt
# import sensor_carbonate as carb


# %% 
# [CHANGE THESE] ------
batch_dir = "/Users/tylerkukla/Documents/GitHub/PRYSM/psm/lake_v2/batch_inputs"
batch_file = "batch_CP_sensitivity_morris-40iter_v1.csv"
# ---------------------


# %% 
# ******************************************************************
# LOOP THROUGH EACH BATCH.CSV ROW
# 
# batch csv
batch_loc = os.path.join(batch_dir, batch_file)
df = pd.read_csv(batch_loc)
# loop
for idx in range(len(df)):
    # track which case we're on in the df
    df['this_case'] = "N"
    df.loc[idx, 'this_case'] = "Y"
    # extract the current row
    row = df.iloc[[idx]]

    # --- read in the batch and default files
    rundict, inc_list = lhf.read_parameter_files(row)
    # --- overwrite the default values
    rundict = lhf.override_defaults(row, rundict)
    # --- update climate if necessary
    if rundict.get('update_clim'):
        rundict = lhf.update_clim(rundict)

    # --- set up the case directory
    lhf.setup_case(rundict, inc_list, df=df)

    # --- run the lake ebm
    lhf.run_ebm(rundict, idx=(idx+1), idx_tot=len(df))
    
    # --- process the lake ebm results and save in results dir
    casedir = lhf.process_ebm_results(rundict)

    # --- run the clumped sensor
    lhf.process_lake_carbonate(casedir, rundict)

    # --- run garbage collectsor before next iter
    gc.collect()




# %%
# ------------------------------------------------------------------
# 
# +++ SCRATCH +++
# 
# ------------------------------------------------------------------



# extract column names and values
keys = df.columns.tolist()  # columns are the keys
values = df.iloc[0].tolist()  # values are the row

default_tag = "**default**"

# %% 
# Update default values
# for key, value in zip(keys, values):
# %% 

key = "b_area"
value = values[3]
value
#%% 

if (key in rundict) & (value != default_tag):
    # attempt to infer type from default dictionary
    default_type = type(rundict[key])
    try:
        if default_type is list:
            # Convert comma-separated string to a list
            rundict[key] = value.split(',') if isinstance(value, str) else list(value)
        else:
            rundict[key] = default_type(value)
    except ValueError:
        raise ValueError(f"Could not convert value '{value}' to type {default_type} for key '{key}'")
elif value == default_tag:
    pass # then do nothing
else:
    if key in ["default_dict_path", "dict_name", "this_case"]:
        # these keys are intentionally omitted from the default dict because
        # they define the default dict
        rundict[key] = value
    else:
        # if a column is not in the default_dict then
        # there's probably a typo
        warnings.warn(f"batch.csv column {key} is not in the default dictionary. Is there a typo?", UserWarning)