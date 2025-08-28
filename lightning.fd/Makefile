
# Compiler settings
# diag-disable=10448 flag is needed to disable warning about how ifort is deprecated
FFLAGS=-diag-disable=10448 -traceback

# Components of program
SRC=read_lightning_bufr.f90 \
    Check_Lightning_QC.f90 \
    netCDFsub_lightning.f90 \
    Check_NLDN.f90 \
    write_bufr_lght.f90 \
    process_Lightning.f90
OBJ=${SRC:.f90=.o}

# Program dependencies
map_util_dir=../share/wps_map_utils
share_dir=../share
INCLD=-I${netcdf_fortran_ROOT}/include \
      -I$(map_util_dir) \
      -I$(share_dir)
LIBS=-L${netcdf_fortran_ROOT}/lib -lnetcdff \
     -L${bufr_ROOT}/lib64 -lbufr_4 \
     -L${w3nco_ROOT}/lib -lw3nco_4
MAP=$(map_util_dir)/misc_definitions_module.o \
    $(map_util_dir)/module_map_utils.o
SHARE=$(share_dir)/kinds.o \
      $(share_dir)/mpasio.o

# Compile individual f90 files
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(INCLD)

# Linking to create final exe
process_Lightning.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(MAP) $(SHARE) $(LIBS)

clean:
	rm -f *.o *.mod process_Lightning.exe 
