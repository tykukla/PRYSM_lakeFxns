# --------------------------------------------------
# 
# run batch simulations of the lake driver psm
# 
# T Kukla, 2025
# --------------------------------------------------
# %% 
import re
import os
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
batch_file = "batch_AMJJASO.csv"
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

    # --- set up the case directory
    lhf.setup_case(rundict, inc_list, df=df)

    # --- run the lake ebm
    lhf.run_ebm(rundict, idx=(idx+1), idx_tot=len(df))
    
    # --- process the lake ebm results and save in results dir
    casedir = lhf.process_ebm_results(rundict)

    # --- run the clumped sensor
    lhf.process_lake_carbonate(casedir, rundict)





# %%
