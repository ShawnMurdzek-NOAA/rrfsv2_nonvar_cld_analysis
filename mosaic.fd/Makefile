
# Compiler settings
# diag-disable=10448 flag is needed to disable warning about how ifort is deprecated
FFLAGS=-diag-disable=10448 -traceback
#FFLAGS=-diag-disable=10448 -traceback -g -warn all -check all

# Components of program
SRC=netCDFsub_mod.f90 \
    read_ncep_binary_mod.f90 \
    read_nssl_binary_mod.f90 \
    read_grib2_mod.f90 \
    module_read_NSSL_mosaic.f90 \
    write_nsslref_mod.f90 \
    mosaic_interp_mod.f90 \
    process_NSSL_mosaic.f90
OBJ=${SRC:.f90=.o}

# Program dependencies
map_util_dir=../share/wps_map_utils
share_dir=../share
INCLD=-I${netcdf_fortran_ROOT}/include \
      -I${g2_ROOT}/include_4 \
      -I$(map_util_dir) \
      -I$(share_dir)
LIBS=-L${netcdf_fortran_ROOT}/lib -lnetcdff \
     -L${jasper_ROOT}/lib64 -ljasper \
     -L${libpng_ROOT}/lib -lpng \
     -L${zlib_ROOT}/lib -lz \
     -L${bacio_ROOT}/lib -lbacio_4 \
     -L${bufr_ROOT}/lib64 -lbufr_4 \
     -L${g2_ROOT}/lib64 -lg2_4
MAP=$(map_util_dir)/misc_definitions_module.o \
    $(map_util_dir)/module_map_utils.o
SHARE=$(share_dir)/kinds.o \
      $(share_dir)/constants.o \
      $(share_dir)/mpasio.o

# Compile individual f90 files
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(INCLD)

# Linking to create final exe
process_NSSL_mosaic.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(MAP) $(SHARE) $(LIBS)

clean:
	rm -f *.o *.mod process_NSSL_mosaic.exe 
