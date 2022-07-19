# ---hcustom build command on linux---
# 'cd' to houdini's install directory and source the following environment with 
# the following command

source houdini_setup

# (optional) set some optimization flags for gcc build (assuming bash shell)
export HCUSTOM_CFLAGS="-ffast-math -O3 -fopenmp"

# now 'cd' to the folder containing aaOceanSOP.c and run

hcustom -I../../../externals/aaOcean/src -I../../../externals/helpers aaOceanSOP.c