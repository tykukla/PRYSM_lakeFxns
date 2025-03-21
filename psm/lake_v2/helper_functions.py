# --------------------------------------------------------
# 
# helper functions for the batch and singlerun versions of
# the lake PSM 
# 
# T Kukla; 2025
# --------------------------------------------------------

# %% 
import os
import re
import time
import shutil
import pickle
import numpy as np
import pandas as pd
import xarray as xr
import subprocess  # for running fortran model
import warnings
from typing import Tuple

# %%
# --------------------------------------------------------------
# 
# RUN THE MODEL
# 
# --------------------------------------------------------------
# --- function to read in the batch.csv and default dictionary
def read_parameter_files(
        df: pd.DataFrame,
)->Tuple[pd.DataFrame, dict]: 
    '''
    Given a row of a batch.csv file, use the columns in that file
    to read in the default dictionary. 

    Parameters
    ----------
    df : pd.DataFrame
        for of the batch.csv dataframe for this single iter
    
    Returns
    -------
    pd.DataFrame
        pandas dataframe of the batch .csv file
    dict 
        dictionary of the default values (not yet overwritten by pd.DataFrame)
    '''
    # --- check that the necessary columns exist
    necessary_df_cols = ["default_dict_path", "dict_name"]
    for col in necessary_df_cols:
        if col not in df.columns:
            raise KeyError(f'Could not find required column "{col}" in the batch.csv file')
    # ---

    # --- read in the defaults
    default_loc = os.path.join(df['default_dict_path'].values[0], 'default_dicts.py')
    # execute the file and access the dict
    default_namespace = {}
    with open(default_loc, "r") as file:
        exec(file.read(), default_namespace)
    # access the specific dictionary
    rundict = default_namespace[df['dict_name'].values[0]]

    # --- read in the .inc list
    inc_loc = os.path.join(df['default_dict_path'].values[0], 'varlists_fixed.py')
    # execute the file and access the dict
    varlist_namespace = {}
    with open(inc_loc, "r") as file:
        exec(file.read(), varlist_namespace)
    # access the specific dictionary
    inc_list = varlist_namespace['inc_vars']


    # --- return results
    return rundict, inc_list


# --- function to overwrite the defaults based on the batch.csv
def override_defaults(
        df: pd.DataFrame,
        rundict: dict,
) -> dict:
    '''
    Take values from the batch dataframe and use them to overwrite 
    the default dictionary. Return the new dictionary.

    Parameters
    ----------
    df : pd.DataFrame
        dataframe of inputs from the batch.csv
    rundict : dict
        dictionary of default values

    Returns
    -------
    dict 
        dictionary after overrides

    '''
    # extract column names and values
    keys = df.columns.tolist()  # columns are the keys
    values = df.iloc[0].tolist()  # values are the row

    # Update default values
    for key, value in zip(keys, values):
        if key in rundict:
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
        else:
            if key in ["default_dict_path", "dict_name", "this_case"]:
                # these keys are intentionally omitted from the default dict because
                # they define the default dict
                rundict[key] = value
            else:
                # if a column is not in the default_dict then
                # there's probably a typo
                warnings.warn(f"batch.csv column {key} is not in the default dictionary. Is there a typo?", UserWarning)
    # return results
    return rundict


# --- function to setup the case directory
def setup_case(
        rundict: dict,
        inc_list: list,
        # f90_list: list,
        df: pd.DataFrame=None,
        casetype: str='batch',
):
    '''
    Create the case directory, copy over and modify the relevant files.
    Note, rundict, inc_list, f90_list all come from inputs in the 
    /outdir/defaults directory.

    Function creates and populates the case directory, it doesn't return anything.

    Parameters
    ----------
    rundict : dict
        dictionary of defaults (if batch mode, this should be already 
        modified to include any overwrites from a batch.csv file)
    inc_list : list
        list of variable names that might be in the .inc file to 
        modify based on rundict
    f90_list : list
        list of variable names that might be in the .inc file to 
        modify based on rundict
    df : pd.DataFrame
        the batch dataframe to save (only if casetype=="batch"). this
        dataframe should have a column indicating which row the current 
        case iteration belongs to
    casetype : str
        ["batch" | "singlerun"]

    Returns
    -------
    '''
    # --- assign variables 
    outdir = rundict['outdir']
    casename = rundict['casename']
    clim_datafile = rundict['datafile']

    # --- create case dir 
    casedir = os.path.join(outdir, 'output', casetype, casename)
    # create the new directory
    os.makedirs(casedir, exist_ok=True)

    # --- modify the template.inc file
    modify_inc_file(rundict, inc_list)

    # --- modify the env_*.f90 file
    modify_f90_file(rundict)# , f90_list)
    
    # --- copy the climate input .txt file
    fn_clim = os.path.join(outdir, 'clim_inputs', clim_datafile)
    shutil.copy(fn_clim, casedir)

    # --- save the dictionary values for reference
    pickel_path = os.path.join(casedir, 'rundict.pkl')
    with open(pickel_path, 'wb') as f:
        pickle.dump(rundict, f)
    
    # --- save the batch.csv values for reference
    if casetype == "batch":
        df.to_csv(os.path.join(casedir, 'batch.csv'), index=False)
    # ---------------------------------------------

    
# --- function to modify a `.inc` file template for a given run
def modify_inc_file(
        rundict: dict,
        inc_list: list,
        casetype: str='batch',
):
    '''
    Read in the .inc file from the main directory and 
    save the modified version in the case directory.

    Function does not return anything.

    Parameters
    ----------
    rundict : dict
        dictionary of defaults (if batch mode, this should be already 
        modified to include any overwrites from a batch.csv file)
    inc_list : list
        list of variable names that might be in the .inc file to 
        modify based on rundict
    casetype : str
        ["batch" | "singlerun"]
        
    Returns
    -------

    '''
    # --- assign vars
    outdir = rundict['outdir']
    template_name = rundict['lake_inc_template']

    # --- read in the template
    file_path = os.path.join(outdir, template_name)
    # read the file into a list of lines
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # process the `parameter` values in each line
    updated_lines = []
    for line in lines:
        for key in inc_list:
            if re.search(rf"\b{key}\b", line):
                if key in rundict['special_format_1']:
                    value = f"'{rundict[key]}'\n"
                    line = re.sub(rf"(\b{key}\s*=\s*)([^,)]+)", rf"\g<1>{value}", line)
                else:
                    value = rundict[key]
                    # replace the numeric value for the matching key
                    line = re.sub(rf"(\b{key}\s*=\s*)([^,)]+)", rf"\g<1>{value}", line)
                    
        updated_lines.append(line)

    # write updated lines to a new file
    output_file = os.path.join(rundict['outdir'], 'output', casetype, rundict['casename'], rundict['lake_inc_filename'])
    with open(output_file, "w") as file:
        file.writelines(updated_lines)


# --- function to modify a `.f90` file template for a given run
def modify_f90_file(
        rundict: dict,
        # f90_list: list,
        casetype: str='batch',
):
    '''
    Modify the f90 template file for a given run. For now, 
    f90_list is commented out. But if we ever update more than one 
    item in f90 file we will create a loop through the f90_list.

    Function does not return anything.

    Parameters
    ----------
    rundict : dict
        dictionary of defaults (if batch mode, this should be already 
        modified to include any overwrites from a batch.csv file)
    casetype : str
        ["batch" | "singlerun"]

    Returns
    -------
    '''
    # --- create the old and new .inc input filenames
    old_incfile = f"'{rundict['env_f90_default_inc']}'"
    new_incfile = f"'{rundict['lake_inc_filename']}'"
    
    # --- read file into a list of lines
    file_path = os.path.join(rundict['outdir'], rundict['env_f90_template'])
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # replace the target string in each line
    updated_lines = [line.replace(old_incfile, new_incfile) for line in lines]

    # write the updated lines back to a new file
    output_file = os.path.join(rundict['outdir'], 'output', casetype, rundict['casename'], rundict['env_f90_filename'])
    with open(output_file, "w") as file:
        file.writelines(updated_lines)


# --- function to run the ebm in fortran 
def run_ebm(
        rundict: dict,
        casetype: str='batch',
        idx: int=0,
        idx_tot: int=-9999,
): 
    '''
    Run the lake ebm in fortran. Can take 10-30+ mins depending on 
    duration and settings of run. 

    Nothing returned. 

    Parameters
    ----------
    rundict : dict
        dictionary of defaults (if batch mode, this should be already 
        modified to include any overwrites from a batch.csv file)
    casetype : str
        ["batch" | "singlerun"]
    idx : int
        the index of the current simulation in the batch (for printout purposes)
    idx_tot : int
        the total number of simulations in the batch (for printout purposes)

    Returns
    -------
    '''
    # --- assign vars
    outdir = rundict['outdir']
    casename = rundict['casename']
    # get the case directory
    casedir = os.path.join(outdir, 'output', casetype, casename)
    # change directory to the case we created in the setups tep
    os.chdir(casedir)

    # --- run the lake model
    # 
    # Note, this can take a while to run depending on the number of timesteps in the 
    # lake model. The example run took my machine ~13 minutes
    # 
    # Specify the fortran source file and the output binary
    # 
    fortran_file = rundict['env_f90_filename']
    output_binary = "lakepsm"       # user-selected name of the compiled executable

    # Construct the gfortran command
    command = ["gfortran", fortran_file, "-o", output_binary]

    # --- Compile
    try:
        # Run the command
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        # print("Compilation successful!")
        # print("Output:", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Compilation failed!")
        print("Error output:", e.stderr)

    # -----------------------------------
    # --- Run (can take ~10-30 mins ish)
    try:
        print(f"Now running case {idx} of {idx_tot} ... ")
        # run the compiled program
        start_time = time.time()  # record the start time
        # ---
        result = subprocess.run(["./" + output_binary], check=True, text=True, capture_output=True)
        # ---
        end_time = time.time()    # record the end time
        elapsed_time = end_time - start_time
        # convert to minutes and seconds
        minutes, seconds = divmod(elapsed_time, 60)
        print(f"--> Completed in {int(minutes)}:{seconds:.2f}")
    except subprocess.CalledProcessError as e:
        print("Execution failed!")
        print("Error output:", e.stderr)





# --------------------------------------------------------------
# 
# PROCESS DATA 
# 
# --------------------------------------------------------------
# %%
# --- function to process the profile data and return as xyz dataframe
#
# TODO: Check my day-of-year averaging at the bottom... Sylvia seems to be 
#       assuming 360 days per year by re-shaping the data to have 12 months
#       when each timestep is separated by 30 days. Is that what the EBM is
#       doing as well? If not, are there implications for the long-term 
#       climatology as each "timestep" drifts by ~5 days per year?
#
def profile_read_process(
        prof_path: str,
        prof_fn: str,
        omit_spinup: bool=True,
        min_threshold: float=0.1,
        return_dataset: bool=True,
        avg_year: bool=True,
        days_per_year: int=360,
) -> pd.DataFrame:
    '''
    Read in the lake energy balance model profile data and shape it into 
    a pandas dataframe with columns for the depth index, temperature,
    and day of year. 

    One tricky problem is the first value is "day of year" and then rows
    fill up with temperature values over depth until the next "day of year" 
    begins. That means we need to figure out if a row is a new day of year, 
    or a continuation of temperature data from the previous row. 

    We do this by comparing the last value in the previous row with the first 
    value of the next one. If they're within min_threshold (or some larger 
    number based on derivative of the previous row) then we assume it's 
    a continuation. Otherwise we assume the first value of the next row is 
    a new day of year. 

    Parameters
    ----------
    prof_path : str
        path to the profile data that you want to read in
    prof_fn : str
        name of the profile data (generally "*_prof.dat")
    omit_spinup : bool
        use True to omit spinup data from the returned pd.DataFrame
    min_threshold : float
        the minimum threshold for a difference in values to detect where
        the next set of day-of-year temperatures begins.
    return_dataset : bool
        whether to return the result as an xr.Dataset (True) or xyz 
        pd.DataFrame (False)
    avg_year : bool
        True to return the average year (avg by timestep)
    days_per_year : int
        days per year to convert day-of-year to month. day-of-year timeslices
        are offset by 30 days each, so we assume a 360 day year for simplicity
        (each month is 30 days)

    Returns
    -------
    pd.DataFrame
        results over depth, temperature, day of year
    '''
    # --- read in the data
    prof_loc = os.path.join(prof_path, prof_fn)
    df = pd.read_csv(prof_loc, delim_whitespace=True, header=None)

    # --- loop through data to construct output dataframe
    tmpdx = 0   # index for rows in a given doy
    doy_row = True  # assume first row includes day of year information
    dfout = pd.DataFrame(columns=["depth_index", "temp_c", "day"]) # empty to hold the final output
    temps_doy = pd.DataFrame()  # empty to hold the temperature output for a given doy

    for rowdx in range((len(df)-1)):      
        if doy_row:
            doy = df.iloc[rowdx,0].copy()
            temps_row = df.iloc[rowdx,1:]
        else:
            temps_row = df.iloc[rowdx,:]

        # create the full row if it doesn't exist
        if tmpdx == 0:
            temps_doy = temps_row
        else:  # add these temps to the last
            temps_doy = pd.concat([temps_doy, temps_row], ignore_index=True)

        # get the last temperature of the current row (avoiding the trailing nans)
        tlast_i = df.iloc[rowdx,df.iloc[rowdx,].last_valid_index()]
        # get the first temperature of the next row
        tfirst_ip1 = df.iloc[(rowdx+1),0]
        # get the difference
        tdiff = tlast_i - tfirst_ip1

        # to get an upper bound on the difference we'd expect between
        # tfirst_ip1 and tlast_i, we take the abs mean diff of row i and 
        # multiply by 2
        diff_threshold = max(np.abs(temps_row.diff()[1:]).mean() * 2, min_threshold)
        temps_row_continue = (np.abs(tdiff) < diff_threshold) & (tdiff >= 0)

        # --- troubleshoot
        # print(f'diff: {tdiff} ; continue: {temps_row_continue}')
        # --- 
        # set up the next row
        if temps_row_continue:
            # then move to the next row knowing it's a continuation of the previous
            doy_row = False
            tmpdx += 1 
        else: # add this temperature round to the full df
            tmpdfout = pd.DataFrame({
                                        "depth_index": np.arange(len(temps_doy)) + 1,
                                        "temp_c": temps_doy,
                                        "day": doy
                                    })
            dfout = pd.concat([dfout, tmpdfout], ignore_index=True)
            doy_row = True
            tmpdx = 0

    # --- check if we should omit spinup data
    if omit_spinup:
        # find the end of the spinup based on the last doy 
        # that is less than the one before it 
        for rowdx in range(len(dfout)-1):
            doy_i = dfout.iloc[rowdx, :]['day']
            doy_ip1 = dfout.iloc[(rowdx+1), :]['day']
            if doy_ip1 < doy_i: 
                end_spinup_dx = rowdx
                
        dfout = dfout[end_spinup_dx+1:].reset_index(drop=True).copy()

    # --- check if we should average years over timesteps
    if avg_year: # take annual average
        dfout['doy'] = dfout['day'] % days_per_year
        dfout = dfout.groupby(['doy', 'depth_index']).mean().reset_index() # take mean
        time_col = "doy"
    else: 
        time_col = "day"

    # --- check if we should convert to xr dataset
    if return_dataset:
        dfout = profile_df_to_ds(dfout, time_col=time_col)

    return dfout


# %% 
# --- function to turn profile dataframe (from profile_read_process) to an xarray dataset
def profile_df_to_ds(
        df: pd.DataFrame,
        depth_col: str="depth_index",
        time_col: str="day",
        value_col: str="temp_c",
) -> xr.Dataset:
    """
    Read in lake model profile data as a pandas dataframe (generally the
    output of profile_read_process) and return the data as an xarray dataset
    with coordinates depth, time, and variable for value

    Parameters
    ----------
    df : pd.DataFrame
        dataframe of lake temperatures over depth and time
    depth_col : str
        name of the depth column in the dataframe
    time_col : str
        name of the time column in the dataframe
    value_col : str
        name of the value column in the dataframe (generally temperature, 'temp_c')

    Returns
    -------
    xr.Dataset
        coordinates depth, time, and variable for value
    """
    # make sure depth index is a numeric
    df[depth_col] = pd.to_numeric(df[depth_col], errors='coerce')
    # pivot the dataframe to create a grid (depth vs. day)
    pivot_df = df.pivot(index=depth_col, columns=time_col, values=value_col)

    # convert the pivoted DataFrame to an xarray Dataset
    dataset = xr.Dataset(
        {
            'temp_c': ([depth_col, time_col], pivot_df.values)  # Define the variable with coordinates
        },
        coords={
            depth_col: pivot_df.index.values,  # depth as coordinate
            time_col: pivot_df.columns.values   # day as coordinate
        }
    )

    return dataset


# %% 
# --- function to process the ebm results
def process_ebm_results(
        rundict: dict,
        casetype: str="batch",
):
    '''
    Process the raw lake ebm output files (*prof.dat and *surf.dat)
    and save the results.

    Parameters
    ----------
    rundict : dict
        dictionary of defaults (if batch mode, this should be already 
        modified to include any overwrites from a batch.csv file)
    casetype : str
        ["batch" | "singlerun"]

    Returns
    -------
    '''
    # --- assign vars
    outdir = rundict['outdir']
    casename = rundict['casename']
    fn_output_surf = rundict['output_surf']
    fn_output_prof = rundict['output_prof']

    # --- create case dir 
    casedir = os.path.join(outdir, 'output', casetype, casename)
    
    # create the results directory
    results_dir = os.path.join(casedir, 'results')
    os.makedirs(results_dir, exist_ok=True)

    # move output .dat files to the results directory
    move_files = [fn_output_surf, fn_output_prof]
    for fn in move_files:
        src = os.path.join(casedir, fn)
        if not os.path.exists(os.path.join(results_dir, fn)):
            shutil.move(src, results_dir)

    # --- PROFILE PROCESS
    # process the profile data 
    dat_prof = profile_read_process(
        prof_path=results_dir,
        prof_fn=fn_output_prof,
        omit_spinup=True,
        return_dataset=True,
        avg_year=True,
    ) 

    # save the result in the results_dir
    nc_filename = fn_output_prof.replace(".dat", ".nc")
    dat_prof.to_netcdf(os.path.join(results_dir, nc_filename))
    # ---

    # --- SURFACE PROCESS
    # process the surface data
    # TODO ....
    # ---


# --- [hailey carbonate sensor functions]    
# --- function to calculate cap47 from temperature
def clumped_sensor(
        temp_c: pd.DataFrame,
        model: str,
) -> np.array:
    
    if model == 'I-CDES90':
        cap47 = 0.0004 * 10**6 / (temp_c + 273.15)**2 + 0.154 # Anderson et al 2021 I-CDES90 ref frame (w/ conversion from Celcius to Kelvin) # not including +/- uncertainties, can add later
    elif model == 'CDES90': 
        cap47 = 1  # will add other functions
    else:
        raise ValueError('Model not recognized')
        
    return cap47.to_numpy()  
 
    # TODO add more models
    
# --- function to calculate temperature from cap47
def clumped_temperature(
     cap47: pd.DataFrame,
     model: str,
) -> np.array:
    
    if model == 'I-CDES90':
        temp_c = ((0.0004 * (10**6) / (cap47 - 0.154))**0.5) - 273.15
    elif model == 'CDES90': 
        temp_c = 1  # will add other functions
    else:
        raise ValueError('Model not recognized')
        
    return temp_c.to_numpy()
 
    # TODO add more models
    
# --- function to generate depth weights  
def generate_depth_weights(
    depth: pd.DataFrame, 
    weight_type: str,
) -> np.ndarray:
    """
    Generate a depth-weighting function.

    Parameters
    ----------
    depth: pd.DataFrame
        DataFrame column of depth values.
    weight_type: str
        Equation type for computing the weighting with depth.

    Returns:
        np.ndarray: Array of weights aligned with the depth values.
    """
    if weight_type == 'uniform':
            weights = [1 if not pd.isna(val) else np.nan for val in depth] # Initialize all weights to 1
            weights = pd.DataFrame(weights, columns=['Weight'])
        
    elif weight_type == 'surface':
        weights = [0 if not pd.isna(val) else np.nan for val in depth] # Initialize all weights to 0
        weights[0] = 1  # Assign weight of 1 to the first element
        weights = pd.DataFrame(weights)

    elif weight_type == 'step':
        depth_min = X  # Set your desired minimum depth
        depth_max = X  # Set your desired maximum depth

        conditions = [
            depth.isna(),
            (depth >= depth_min) & (depth <= depth_max)
        ]
        choices = [np.nan, 1]
        weights = np.select(conditions, choices, default=0)
        weights = pd.DataFrame(weights, index=depth.index)
        
    elif weight_type == 'normal_dist':
        depth_index_mean = round(depth.values.mean(), 2)
        depth_index_std = round(depth.values.std(), 2)
        weights = (1 / (depth_index_std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((depth - depth_index_mean) / depth_index_std) ** 2)
        # Normalize the weights so they sum to 1
        weights /= weights.sum()
    else:
        raise ValueError('Weight type not recognized')
    
    return weights

# --- function to generate time weights
def generate_time_weights(
    timesteps:np.ndarray, 
    weight_type: str, 
    selected_months= list,
) -> np.ndarray: 
    """
    Generate a time-weighting function.
    
    Parameters
    ----------
    depth: pd.DataFrame
        DataFrame column of depth values.
    weight_type: str
        Equation type for computing the weighting with doy.
    selected_months: list 
        List of months from dict to apply the 'range' time weighting function to (e.g. selected_months=['January', 'February']), otherwise 'None'.

    Returns:
       np.ndarray: Array of weights aligned with doy values.
    """
    months = {
        'MAMJ': range(60, 182),  # March to June
        'AMJ': range(91, 182),  # April to June  
        'AMJJ': range(91, 213),  # April to July
        'AMJJASO': range(91, 305),  # April to October
        'JJA': range(152, 244),  # July to August 
        'ASO': range(213, 305),  # August to October
        'January': range(0, 32),
        'February': range(32, 60),
        'March': range(60, 91),
        'April': range(91, 121),
        'May': range(121, 152),
        'June': range(152, 182),
        'July': range(182, 213),
        'August': range(213, 244),
        'September': range(244, 274),
        'October': range(274, 305),
        'November': range(305, 335),
        'December': range(335, 366)
        }
    
    if weight_type == 'uniform':
        weights = [1] * len(timesteps)  # Initialize all weights to 1
    
    elif weight_type == 'range':
        # ... check to see if any selected months don't occur in the months dict
        #     return a message to let the user know (maybe they mis-spelled something)
        missing_keys = [key for key in selected_months if key not in months]
        if missing_keys:
            print(f"Warning: selected month option(s) {missing_keys} is not valid. Please select from {list(months.keys())}")

        # ... filter the months dict for just the selected months
        selected_months_dict = {key: months[key] for key in selected_months if key in months}

        # ... return 1 if the timestep is within one of the ranges, return 0 otherwise
        weights = [1 if any(start <= num <= end for start, end in selected_months_dict.values()) else 0 for num in timesteps]
        
    elif weight_type == 'normal_dist':
        time_index_mean = round(timesteps.values.mean(), 2)
        time_index_std = round(timesteps.values.std(), 2)
        weights = (1 / (time_index_std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((timesteps - time_index_mean) / time_index_std) ** 2)
        # Normalize the weights so they sum to 1
        weights /= weights.sum()
    else:
        raise ValueError('Weight type not recognized')
    
    return weights

# --- function to process lake profile data for clumped carbonate sensor
def process_lake_data(
    lakedata: xr.Dataset,
    ) -> xr.Dataset:

    timesteps = lakedata['doy'].values.tolist()
    cap47 = []

    for ts in (timesteps) :
        tmpds = lakedata.sel(doy = ts) # pull out this timestep
        tmpdf = tmpds.to_dataframe().reset_index()  # convert it to pandas df

        temp_c = pd.DataFrame(tmpdf['temp_c'])['temp_c']
        depth = pd.DataFrame(tmpdf['depth_index'])['depth_index']
        
        ts_cap47 = clumped_sensor(temp_c, 'I-CDES90')
        ts_cap47 = pd.DataFrame(ts_cap47).copy()
        ts_cap47 = ts_cap47.rename(columns={'temp_c': 'cap47'})
        cap47.append(ts_cap47)
        nan_count = temp_c[::-1].isna().cumprod().sum()
        depth.iloc[-nan_count:] = np.nan 

        depth_weights = generate_depth_weights(depth, 'uniform')  
        depth_weights = pd.DataFrame(depth_weights.values.flatten())
       
        if ts == timesteps[0]:
            df = pd.DataFrame(tmpdf['doy'].values, columns=['doy'])
            cap47_depth_wtd_mean = (ts_cap47.values * depth_weights).sum() / depth_weights.sum()
            temp_c_depth_wtd_mean = (temp_c * depth_weights).sum() / depth_weights.sum()
            df = df.assign(temp_c_depth_wtd_mean = temp_c_depth_wtd_mean).dropna()
            df = df.assign(cap47_depth_wtd_mean = cap47_depth_wtd_mean).dropna()
        else:
            tdf = pd.DataFrame(tmpdf['doy'].values, columns=['doy'])
            cap47_depth_wtd_mean = (ts_cap47.values * depth_weights).sum() / depth_weights.sum()
            temp_c_depth_wtd_mean = (temp_c * depth_weights).sum() / depth_weights.sum()
            tdf = tdf.assign(temp_c_depth_wtd_mean = temp_c_depth_wtd_mean).dropna()
            tdf = tdf.assign(cap47_depth_wtd_mean = cap47_depth_wtd_mean).dropna()
            df = pd.concat([tdf, df], ignore_index=True)
            df = df.dropna()      

    T47_c_depth_wtd_mean = clumped_temperature(df['cap47_depth_wtd_mean'], 'I-CDES90')

    df = df.set_index('doy', drop=True)
    df = df.assign(T47_c_depth_wtd_mean = T47_c_depth_wtd_mean)
    df = xr.Dataset.from_dataframe(df)
    time_weights = generate_time_weights(timesteps, 'uniform') 
    time_weights = time_weights[::-1]
    df = df.assign(time_weights=(['doy'], time_weights))
    df['cap47_time_wtd_mean'] = (df['cap47_depth_wtd_mean']* df['time_weights']).sum() / df['time_weights'].sum()
    df['T47_c_MAT'] = (df['T47_c_depth_wtd_mean'] * df['time_weights']).sum() / df['time_weights'].sum()    

    final_output_init = xr.combine_by_coords([df, lakedata])

    final_output_init = final_output_init.assign(depth_weights=(['depth_index'], depth_weights.to_numpy().flatten()))
    
    ds_cap47 = xr.Dataset(
    {
        "cap47": (["doy", "depth_index"], np.stack([df[0].values for df in cap47]))
    },
    coords={"doy": final_output_init['doy'].values, "depth_index": final_output_init['depth_index'].values}
    )

    final_output = xr.merge([ds_cap47, final_output_init])
    depth_weights_expanded = np.tile(final_output['depth_weights'], (len(final_output['doy']), 1))
    final_output["depth_weights"] = (['doy', 'depth_index'], depth_weights_expanded)
    final_output = final_output.dropna(dim='depth_index', how='all')
    return final_output