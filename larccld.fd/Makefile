
# Compiler settings
# diag-disable=10448 flag is needed to disable warning about how ifort is deprecated
FFLAGS=-diag-disable=10448 -traceback
#FFLAGS=-diag-disable=10448 -traceback -g -warn all -check all

# Components of program
SRC=larccld_utils_mod.f90 \
    larccld_bufr_io_mod.f90 \
    process_NASALaRC_cloud.f90
OBJ=${SRC:.f90=.o}

# Program dependencies
map_util_dir=../share/wps_map_utils
share_dir=../share
INCLD=-I${netcdf_fortran_ROOT}/include \
      -I$(map_util_dir) \
      -I$(share_dir)
LIBS=-L${netcdf_fortran_ROOT}/lib -lnetcdff \
     -L${bufr_ROOT}/lib64 -lbufr_4
MAP=$(map_util_dir)/misc_definitions_module.o \
    $(map_util_dir)/module_map_utils.o
SHARE=$(share_dir)/kinds.o \
      $(share_dir)/constants.o \
      $(share_dir)/mpasio.o

# Compile individual f90 files
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(INCLD)

# Linking to create final exe
process_larccld.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(MAP) $(SHARE) $(LIBS)

clean:
	rm -f *.o *.mod process_larccld.exe 
