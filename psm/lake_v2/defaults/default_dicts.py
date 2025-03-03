# ---------------------------------------------------
# 
# default dictionaries for the lake PSM 
# 
# ---------------------------------------------------

defaults_main = {
    # ---------------------------------------------------------------------
    # --- general model setup
    # 
    # [change]
    'casename': '_test_case',  # name of the run case
    'datafile': 'CP_SLIM_lowTopo_500ppm_test_input.txt',  # name of the climate input data file 
    #
    # --- model directory structure
    # [change]
    'outdir': '/Users/tylerkukla/Documents/GitHub/PRYSM/psm/lake_v2',  # location to save output data
    # [do not change]
    'clim_subdir': 'clim_inputs',    # where climate model .txt input files are stored
    'batch_subdir': 'batch_inputs',  # where batch .csv files are stored
    #
    # --- model files [change as needed..]
    # templates (must be located in 'outdir')
    'lake_inc_template': 'template.inc',     # name of the '.inc' template to modify
    'env_f90_template': 'env_heatflux_locked_template.f90',  # name of the .f90 template to read in
    'env_f90_default_inc': 'Tanganyika.inc',  # name of the .inc file included in the env_f90 template (that we need to replace)
    # files to generate (will be created in outdir/output/*/casename)
    'lake_inc_filename': 'lake_setup.inc',   # name of the '.inc' file to generate
    'env_f90_filename': 'env_heatflux_locked.f90',  # name of the .f90 file to save 
    # output file names
    'output_prof': '_test_prof.dat',  # name of the output profile file
    'output_surf': '_test_surf.dat',  # name of the output surface file
    # [do not change]
    'special_format_1': ['datafile', 'output_prof', 'output_surf'],  # keys that require single quotes w/in double and \n at end of single quote

    # ---------------------------------------------------------------------
    # --- template.inc inputs
    # 
    # --- lake-specific parameters
    "oblq": 23.4,      # obliquity
    "xlat": -6.30,     # latitude (negative for south)
    "xlon": 29.5,      # longitude (negative for west)
    "gmt": "+3",       # [hrs] local time relative to gmt in hours
    "max_dep": 999,    # [m] depth of lake at sill 
    "basedep": 733.,   # [m] elevation of basin at bottom 
    "b_area": 23100000, # [ha] area of catchment + lake
    "cdrn": 2.0e-3,    # neutral drag coefficient (1.8 HAD 1.7GISS 1.2CCSM)
    "eta": 0.065,      # [m-1] shortwave extinction coefficient
    "f": 0.3,          # fraction of advected air over lake
    "alb_slush": 0.4,  # [] albedo of melting snow
    "alb_snow": 0.7,   # [] albedo of non-melting snow
    "depth_begin": 570, # [m] prescribed depth
    "salty_begin": 0.0, # [ppt] prescribed salinity
    "o18air": -14.0,   # [per mil VSMOW] d18O of air above lake
    "deutair": -96.,   # [per mil VSMOW] dD of air above lake
    "tempinit": 18.,   # [degC] temperature to initialize lake at in INIT_LAKE subroutine
    "deutinit": -100., # [per mil VSMOW] dD to initialize lake at in INIT_LAKE subroutine
    "o18init": -10.,   # [per mil VSMOW] d18O to initialize lake at in INIT_LAKE subroutine
    # 
    # --- simulation-specific parameters
    "nspin": 10,          # [yrs] number of years for spinup
    "bndry_flag": ".false.", # [".true." | ".false."] true for explicit boundary layer computation (presently only for sigma coord climate models)
    "sigma": 0.9925561,   # sigma level for boundary flag
    "wb_flag": ".false.", # [".true." | ".false."] true for variable lake depth
    "iceflag": ".false.", # [".true." | ".false."] true for variable ice cover
    "s_flag": ".false.",  # [".true." | ".false."] true for variable salinity
    "o18flag": ".false.", # [".true." | ".false."] true for variable d18O
    "deutflag": ".false.", # [".true." | ".false."] true for variable dD
    "z_screen": 5.0,      # [m] height of meteorological inputs

    # ---------------------------------------------------------------------
    # --- helper function inputs
    #
    # --- `profile_read_process` (postprocess profile output)
    # only the vars that aren't generated by the run-script are handled here
    'prof_omit_spinup': True,   # [bool] whether to omit the spinup data from the output file
    'prof_min_threshold': 0.1,  # [float] minimum numerical threshold for determining if a new row in the prof.dat file starts with a temperature value or day-of-year
    'prof_return_dataset': True, # [bool] whether to return xarray dataset (true) or pandas dataframe in xyz (false)
    'prof_avg_year': True,      # [bool] whether to return data for each year (false) or the long-term averaged monthly data (true)
    'prof_days_per_year': 360,  # [int] days per year used for the long term averaging (not used if `prof_avg_year` is false) (360 is consistent with Dee et al. assuming 12 months/yr from 30-day timesteps)
}
    
