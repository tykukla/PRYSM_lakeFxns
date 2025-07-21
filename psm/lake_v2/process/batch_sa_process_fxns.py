# --------------------------------------------------
# 
# functions to process sensitivity analysis 
# simulations (and other batches)
# 
# --------------------------------------------------
# %% 
import glob
import os 
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import xarray as xr


# %%

# [ create dictionary of all output paths in the batch ]
def create_output_dict(
        resultdir: str,
        saname: str,
)->dict:
    '''
    Create a dictionary of output directories where each
    key represents the run index integer. (this gets passed
    to other steps in the postprocessing)

    Parameters
    ----------
    resultdir : str
        path to all of the directories we want to grab
    saname : str
        name of the sensitivity analysis, or the pattern used
        for identifying the directories of interest
    
    Returns
    -------
    dict
        dictionary of paths to output directories, labeled with
        run index integers
    '''

    # --- create output dict (key is the run integer, value is the output directory)
    saname_pattern = f'{saname}*'
    # get all matches to the pattern
    sa_dir_matches = glob.glob(os.path.join(resultdir, saname_pattern))
    # filter for dirs only 
    sa_dirs = [d for d in sa_dir_matches if os.path.isdir(d)]
    # build a dictionary: {int_suffix: dir_name}
    dir_dict = {}
    pattern = re.compile(re.escape(saname) + r"(\d+)$")

    for d in sa_dirs:
        base = os.path.basename(d)  # get just the final directory name
        match = pattern.search(base)
        if match:
            key = int(match.group(1))
            dir_dict[key] = d  # keep full path as value

    # sort the dictionary by run integer
    dir_dict = dict(sorted(dir_dict.items()))

    # return result
    return dir_dict


# [ identify the variables we want to add to the xr dataset ]
def is_numeric_or_default(col):
    col = col.dropna()
    for val in col:
        if val == "**default**":
            continue
        # Exclude boolean values explicitly
        if isinstance(val, bool):
            return False
        try:
            float(val)
        except (ValueError, TypeError):
            return False
    return True


# [ read in and format nc result from a single case ]
def read_and_format_nc_result(
        this_int: int,
        outdir: str,
        ncname: str = "clumpedSens_CP_control_prof.nc",
        resultsdir: str = "results",
        batchname: str = "batch.csv",
        rundictname: str = "rundict.pkl",
)-> xr.Dataset:
    '''
    Read in the model results for a single run. Add variables 
    for the values that were modified by the batch.csv file. 
    Add a dimension for the run index integer (this_int)

    Parameters
    ----------
    this_int : int
        run index integer for this run
    outdir : str
        path to the output directory for this run
    ncname : str
        name of the .nc file to read in. should be located 
        in `outdir/resultsdir/*`
    resultsdir : str
        name of the results directory within the output dir
        (usually just "results")
    batchname : str
        name of the batch file located in the output dir
    rundictname : str
        name of the run dictionary file located in the output dir

    Returns
    -------
    xr.Dataset
        The .nc dataset that we read in formatted with the relevant
        variables and the run index integer set as a new coordinate 
    '''
    # --- read in the batch .csv
    dfbatch = pd.read_csv(os.path.join(outdir, batchname))
    # pull out this row
    thisrow = dfbatch[dfbatch['this_case'] == "Y"]
    # find the columns that we modified with numeric values
    numeric_cols = [col for col in thisrow.columns if is_numeric_or_default(thisrow[col])]

    # --- read in the rundict 
    rundict = pd.read_pickle(os.path.join(outdir, rundictname))

    # --- read in the .nc file
    ds = xr.open_dataset(os.path.join(outdir, resultsdir, ncname))
    # add coordinate
    ds = ds.expand_dims(case_int=[this_int])
    # add vars
    for key in numeric_cols:
        values = rundict[key]  # could be scalar or array
        ds[key] = xr.DataArray([values], dims=["case_int"])

    # return result
    return ds


# [ read in all .ncs and merge into single dataset ]
def read_and_merge_all_results(
        dir_dict: dict,
        ncname: str = "clumpedSens_CP_control_prof.nc",
        resultsdir: str = "results",
        batchname: str = "batch.csv",
        rundictname: str = "rundict.pkl",
        append_by: str = "case_int",
        troubleshoot_dims: bool = False,
        expected_doy: list = [15.,  45.,  75., 105., 135., 165., 195., 225., 255., 285., 315., 345.],
)-> xr.Dataset:
    '''
    Read in all datasets in the dir_dict and merge them into one.

    Parameters
    ----------
    dir_dict : dict
        dictionary of dirs with keys set to the run index integer
        (output of the create_output_dict fxn)
    ncname : str
        name of the .nc file to read in. should be located 
        in `outdir/resultsdir/*`
    resultsdir : str
        name of the results directory within the output dir
        (usually just "results")
    batchname : str
        name of the batch file located in the output dir
    rundictname : str
        name of the run dictionary file located in the output dir
    append_by : str
        name of the coordinate that we want to append by when
        we merge all vars (should be the unique valued coord, e.g.,
        the run index integer)
    troubleshoot_dims : bool
        [True | False] whether to check all individual datasets
        for off dimensions that will create NAs upon merge
    expected_doy : list
        expected values for the day of year (doy) dim used for 
        troubleshooting
    
    Returns
    -------
    xr.Dataset
        dataset of all results from the N runs merged together
    '''
    # --- loop through each run
    ds_list = []
    for this_int, outdir in dir_dict.items():
        # get the single .nc file
        tmpds = read_and_format_nc_result(
            this_int,
            outdir, 
            ncname,
            resultsdir,
            batchname,
            rundictname,
        )
        # add it to the output list
        ds_list.append(tmpds)

    if troubleshoot_dims:
        messed_up = []
        for i, ds in enumerate(ds_list):
            doy_vals = ds.coords["doy"].values
            if not np.array_equal(np.sort(doy_vals), np.sort(expected_doy)):
                messed_up.append(ds.case_int.values[0])
                print(f"⚠️ Mismatch in experiment {i}:")
                print(f"  Found: {np.sort(doy_vals)}")
                print(f"  Missing: {set(expected_doy) - set(doy_vals)}")
                print(f"  Extra: {set(doy_vals) - set(expected_doy)}")
        
        return messed_up
    
    else:
        # merge results together
        dsout = xr.concat(ds_list, dim = append_by)
        # return result
        return dsout


# [ get lake surface and bottom temperatures ]
def get_lakesurface_mean_temp(da):
    # da: temp_c for a single case_int, dims: [doy, depth]
    # 1. Mask out depths where all values are NA
    valid_depths = da.notnull().any(dim='doy')
    
    # 2. Find shallowest valid depth
    shallowest_depth = da['depth_index'].where(valid_depths).min(dim='depth_index')
    
    # 3. Select temp at shallowest depth and average over doy
    temp_at_shallow = da.sel(depth_index=shallowest_depth)
    return temp_at_shallow.mean(dim='doy')

def get_lakebottom_mean_temp(da):
    # da: temp_c for a single case_int, dims: [doy, depth]
    # 1. Mask out depths where all values are NA
    valid_depths = da.notnull().any(dim='doy')
    
    # 2. Find shallowest valid depth
    shallowest_depth = da['depth_index'].where(valid_depths).max(dim='depth_index')
    
    # 3. Select temp at shallowest depth and average over doy
    temp_at_shallow = da.sel(depth_index=shallowest_depth)
    return temp_at_shallow.mean(dim='doy')

def get_unwtd_lake_topbottom_temps(
        ds: xr.Dataset,
        coord_group: str='case_int'
)->xr.Dataset:
    '''
    add variables defined solely over `case_int` to dataset
    for lake surface and bottom temperature (and difference)

    Parameters
    ----------
    ds : xr.Dataset
        input dataset including all batch simulations
        must have a `case_int` coordinate (or coord_group=X coord)
    coord_group : str
        name of the coordinate that we group data by before 
        applying the function
    
    Returns
    -------
    xr.Dataset
        input dataset with added temperature variables
    '''
    # get surface and bottom temperature arrays
    lst_unwtd_mean = ds['temp_c'].groupby('case_int').map(get_lakesurface_mean_temp)
    lakebottom_temp_unwtd_mean = ds['temp_c'].groupby('case_int').map(get_lakebottom_mean_temp)
    # add to ds
    ds['LST_unwtd_mean'] = lst_unwtd_mean
    ds['lakebottom_temp_unwtd_mean'] = lakebottom_temp_unwtd_mean
    ds['LST_minus_lakebottom_temp'] = ds['LST_unwtd_mean'] - ds['lakebottom_temp_unwtd_mean']
    # return result
    return ds


# [ get mean air temperature (weighted and unweighted) for a given case ]
def get_mean_air_temperature(
        tmpds: xr.Dataset,
        outdir: str,
        climfile_prefix: str="CP_SLIM_modernTopo_280ppm_input+",
        tempcolumn: int=2,
        daycolumn: int=1,
        days_per_year: int=360,
        convert_K_to_C: bool=True
):
    '''
    Calculate the mean air temperature, unweighted and 
    using the calcite weighting.
    '''

    # --- get the runname
    runname = os.path.basename(outdir)
    clime_fn = f"{climfile_prefix}{runname}.txt"

    # --- read in the climate file
    dfclim = pd.read_csv(os.path.join(outdir, clime_fn), sep='\s+', header=None)
    # add a doy column
    dfclim['doy'] = ((dfclim[daycolumn] - 1) % days_per_year) + 1
    # get mean by doy
    dfclim = dfclim.groupby('doy').mean(numeric_only=True).reset_index()

    # --- weight temperature by time
    meantemp_wtd = np.sum(dfclim[tempcolumn] * tmpds['time_weights'].values) / np.sum(tmpds['time_weights'].values)
    meantemp = dfclim[tempcolumn].mean()

    # --- convert to C if asked
    if convert_K_to_C:
        meantemp_wtd -= 273.15
        meantemp -= 273.15

    # return results
    return meantemp_wtd, meantemp


# [ add mean air temperature (weighted and unweighted) for each case to main dataset ]
def add_air_temperature_to_ds(
        ds: xr.Dataset,
        dir_dict: dict,
        climfile_prefix: str="CP_SLIM_modernTopo_280ppm_input+",
        tempcolumn: int=2,
        daycolumn: int=1,
        days_per_year: int=360,
        convert_K_to_C: bool=True
)-> xr.Dataset:
    '''
    Loop through dir_dict and add the air temperature to the dataset 
    '''
    # --- loop through
    wtd_templist, templist = [], []
    for this_int, outdir in dir_dict.items():
        # get the temporary dataset
        tmpds = ds.sel(case_int=this_int).copy()
        meantemp_wtd, meantemp = get_mean_air_temperature(
                tmpds,
                outdir,
                climfile_prefix=climfile_prefix,
                tempcolumn=tempcolumn,
                daycolumn=daycolumn,
                days_per_year=days_per_year,
                convert_K_to_C=convert_K_to_C
        ) 
        wtd_templist.append(meantemp_wtd)
        templist.append(meantemp)
        
    # add to xarray dataset
    ds['airtemp_mean_wtd'] = xr.DataArray(wtd_templist, dims='case_int', coords={'case_int': ds['case_int']})
    ds['airtemp_mean'] = xr.DataArray(templist, dims='case_int', coords={'case_int': ds['case_int']})

    # return result
    return ds


def add_extra_airtemp_vars(
        ds: xr.Dataset,
):
    '''
    Add extra variables based on the air temperature data
    to the main dataset. 
    
    Assumes the "control" scenario 
    is the last in the case_int list. 
    '''
    # --- get the delta in air and lake temperature
    ctrl_idx = len(ds.case_int.values)
    # get control data
    # [ air ]
    ds['ctrl_airteamp_mean_wtd'] = ds.sel(case_int = ctrl_idx)['airtemp_mean_wtd']
    ds['ctrl_airteamp_mean'] = ds.sel(case_int = ctrl_idx)['airtemp_mean']
    # [ lake surface ]
    ds['ctrl_lst_mean_wtd'] = ds.sel(case_int = ctrl_idx)['Actual_MALST_PRYSM']
    ds['ctrl_lst_mean'] = ds.sel(case_int = ctrl_idx)['LST_unwtd_mean']
    # [ lake column ]
    ds['ctrl_lakeCol_mean'] = ds.sel(case_int = ctrl_idx)['T47C_MALT']

    # get difference
    # [ air ]
    ds['Delta_airtemp_mean_wtd'] = ds['airtemp_mean_wtd'] - ds['ctrl_airteamp_mean_wtd']
    ds['Delta_airtemp_mean'] = ds['airtemp_mean'] - ds['ctrl_airteamp_mean']
    # [ lake surface ]
    ds['Delta_lst_mean_wtd'] = ds['Actual_MALST_PRYSM'] - ds['ctrl_lst_mean_wtd']
    ds['Delta_lst_mean'] = ds['LST_unwtd_mean'] - ds['ctrl_lst_mean']
    # [ lake column ]
    ds['Delta_lakeCol_mean_wtd'] = ds['T47C_MALT'] - ds['ctrl_lakeCol_mean']

    # ratio between delta air and delta lake surface
    ds['Delta_lst_per_Delta_air_wtd'] = ds['Delta_lst_mean_wtd'] / ds['Delta_airtemp_mean_wtd']
    ds['Delta_lst_per_Delta_air'] = ds['Delta_lst_mean'] / ds['Delta_airtemp_mean']
    ds['Delta_lakeCol_per_Delta_air'] = ds['Delta_lakeCol_mean_wtd'] / ds['Delta_airtemp_mean_wtd']

    # lake temp minus air temp
    ds['lst_min_air_temp'] = ds['LST_unwtd_mean'] - ds['airtemp_mean']

    # return result
    return ds


# -----------------------------------------------------
def batch_sa_workflow(
        resultpath: str,
        batchname: str,
        save_ds: bool,
        savehere: str='/Users/tylerkukla/Documents/GitHub/PRYSM/psm/lake_v2/output/batch/_SA-proc'
)->xr.Dataset:
    '''
    Complete full processing workflow to prepare 
    batch simulations for SA. Note, each sub-function can have
    some default params that we haven't touched here. Add 
    to the function input if necessary. 

    Parameters
    ----------
    resultpath : str
        path to dir with all result output dirs (e.g., /lake_v2/output/batch/)
    batchname : str
        name of the set of batch runs (if a single run is called "MyName_1" then
        the batchname is likely "MyName")
    save_ds : bool
        if True, save the dataset produced by the workflow
    savehere : str
        name of the path where dataset would be saved (filename is batchname)

    Returns
    -------
    xr.Dataset
        compiled and processed data from batch
    '''
    # [ create dictionary of output directories ] 
    dir_dict = create_output_dict(resultpath, f'{batchname}_')
    
    # [ read and merge all .nc results ]
    ds = read_and_merge_all_results(dir_dict)

    # [ get lake surface and bottom temps ]
    ds = get_unwtd_lake_topbottom_temps(ds)

    # [ add air temp ]
    ds = add_air_temperature_to_ds(
                ds,
                dir_dict,
        )
    
    # [ add extra temperature calculations ]
    ds = add_extra_airtemp_vars(ds)

    # [ save if asked ]
    if save_ds:    
        # save dataset
        ds.to_netcdf(os.path.join(savehere, f"{batchname}.nc"))

    # --- return result
    return ds 


# %%
# -----------------------------------------------------------
def sa_quickplot_scatter(
        SA_result,
        problem: dict,
        ptitle: str= "Morris Sensitivity Analysis",
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
    plt.show()


def sa_quickplot_scatter_compareMu(
        SA_result,
        problem: dict,
        ptitle: str= "Morris Sensitivity Analysis",
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
    plt.show()



varname_dict_in = {
    "depth_begin": "Water depth",
    "salty_begin": "Salinity",
    "oblq": "Orbital obliquity",
    "xlat": "Site latitude",
    "xlon": "Site longitude",
    "cdrn": "Neutral drag coeff",
    "eta": "SW extinction coeff",
    "f": "Fraction advected air", 
    "TSplus": "Air temperature",
    "RHplus": "Relative humidity",
    "FSDSplus": "Incoming SW",
    "FLDSplus": "Incoming LW",
    "Ufac": "Wind speed",
    "PSplus": "Surface pressure",
    "PRECTplus": "Precipitation",
    "runoffplus": "Runoff"
}

def sa_quickplot_bar(
        SA_result,
        varname_dict: dict=varname_dict_in,
        barcolor: str= 'steelblue',
        linecolor: str='black',
        bar_order: list=None,
        ptitle: str= "'Morris Sensitivity Analysis",
        fs_ylab: int | float = 12,
        fs_xlab: int | float = 12,
        fs_xticks: int | float = 10,
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
    plt.show()
