# --------------------------------------------------------
# 
# lists of names of relevant vars for updating each of 
# the two lake ebm files (the `.inc` and `.f90` files)
# 
# --------------------------------------------------------
#
# [DO NOT CHANGE without good reason]
#

# -- variable names to loop through when modifying the .f90 file
f90_vars = [
    "lake_inc_filename", # name of the .inc file to update in .f90
]

# -- variable names to loop through when modifying the .inc file
inc_vars = [
    # --- files
    "datafile",  # name of the climate input .txt file
    "output_prof",  # name to use for the profile output file
    "output_surf",  # name to use for the surface output file
    # --- lake-specific parameters
    "oblq",      # obliquity
    "xlat",      # latitude (negative for south)
    "xlon",      # longitude (negative for west)
    "gmt",       # [hrs] local time relative to gmt in hours
    "max_dep",   # [m] depth of lake at sill 
    "basedep",   # [m] elevation of basin at bottom 
    "b_area",    # [ha] area of catchment + lake
    "cdrn",      # neutral drag coefficient (1.8 HAD 1.7GISS 1.2CCSM)
    "eta",       # [m-1] shortwave extinction coefficient
    "f",         # fraction of advected air over lake
    "alb_slush", # [] albedo of melting snow
    "alb_snow",  # [] albedo of non-melting snow
    "depth_begin", # [m] prescribed depth
    "salty_begin", # [ppt] prescribed salinity
    "o18air",    # [per mil VSMOW] d18O of air above lake
    "deutair",   # [per mil VSMOW] dD of air above lake
    "tempinit",  # [degC] temperature to initialize lake at in INIT_LAKE subroutine
    "deutinit",  # [per mil VSMOW] dD to initialize lake at in INIT_LAKE subroutine
    "o18init",   # [per mil VSMOW] d18O to initialize lake at in INIT_LAKE subroutine
    # 
    # --- simulation-specific parameters
    "nspin",        # [yrs] number of years for spinup
    "bndry_flag",   # [".true." | ".false."] true for explicit boundary layer computation (presently only for sigma coord climate models)
    "sigma",        # sigma level for boundary flag
    "wb_flag",      # [".true." | ".false."] true for variable lake depth
    "iceflag",      # [".true." | ".false."] true for variable ice cover
    "s_flag",       # [".true." | ".false."] true for variable salinity
    "o18flag",      # [".true." | ".false."] true for variable d18O
    "deutflag",     # [".true." | ".false."] true for variable dD
    "z_screen",     # [m] height of meteorological inputs
]
