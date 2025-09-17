
# Compiler settings
# diag-disable=10448 flag is needed to disable warning about how ifort is deprecated
#FFLAGS=-diag-disable=10448 -traceback
FFLAGS=-diag-disable=10448 -traceback -g -warn all -check all

# Components of program
SRC=kinds.f90 \
    constants.f90 \
    mpasio.f90 \
    map_proj_helper_mod.f90
OBJ=${SRC:.f90=.o}

# Program dependencies
map_util_dir=./wps_map_utils
INCLD=-I${netcdf_fortran_ROOT}/include \
      -I$(map_util_dir)

# Compile individual f90 files
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(INCLD)

main: $(OBJ)

clean:
	rm -f *.o *.mod
