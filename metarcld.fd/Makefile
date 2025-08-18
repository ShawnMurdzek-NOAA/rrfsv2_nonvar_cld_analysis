
# Compiler settings
# diag-disable=10448 flag is needed to disable warning about how ifort is deprecated
FC=mpiifort
FFLAGS=-diag-disable=10448 -traceback

# Components of program
SRC=kinds.f90 \
    constants.f90 \
    mpasio.f90 \
    cld_parm_array_module.f90 \
    reorg_metar_cloud_regular.f90 \
    read_prepbufr_metarcld.f90 \
    process_metar_cloud.f90
OBJ=${SRC:.f90=.o}

# Program dependencies
map_util_dir=../wps_map_utils
INCLD=-I${netcdf_fortran_ROOT}/include \
      -I$(map_util_dir)
LIBS=-L${netcdf_fortran_ROOT}/lib -lnetcdff \
     -L${bufr_ROOT}/lib64 -lbufr_4 \
     -L${w3nco_ROOT}/lib -lw3nco_4
MAP=$(map_util_dir)/misc_definitions_module.o \
    $(map_util_dir)/module_map_utils.o

# Compile individual f90 files
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(INCLD)

# Linking to create final exe
process_metarcld.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(MAP) $(LIBS)

clean:
	rm -f *.o *.mod process_metarcld.exe 
