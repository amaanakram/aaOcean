# hcustom build command on windows
# please make sure that you are using the compiler required for your version of houdini
# also, you need VC directories in the PATH variable for windows. See HDK docs
# Launch Houdini Command Line Tools
# run the following


# now 'cd' to aaOceanSOP.c dir
set HCUSTOM_CFLAGS=-fp:fast -Ox -Oy -GL -D_OPENMP -openmp
hcustom -I ../../../externals/aaOcean/src -I ../../../externals/helpers aaOceanSOP.c


# plugin installed to %UserProfile%\Documents\houdini<version>\dso