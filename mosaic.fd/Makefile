
# Compiler settings
# diag-disable=10448 flag is needed to disable warning about how ifort is deprecated
FC=mpiifort
FFLAGS=-diag-disable=10448 -traceback

# Components of program
SRC=kinds.f90 \
    constants.f90 \
    mpasio.f90 \
    netCDFsub.f90 \
    read_ncep_binary.f90 \
    read_nssl_binary.f90 \
    read_grib2_mod.f90 \
    module_read_NSSL_mosaic.f90 \
    write_netcdf_ref.f90 \
    write_bufr_ref.f90 \
    process_NSSL_mosaic.f90
OBJ=${SRC:.f90=.o}

# Program dependencies
map_util_dir=../wps_map_utils
INCLD=-I${netcdf_fortran_ROOT}/include \
      -I${g2_ROOT}/include_4 \
      -I$(map_util_dir)
LIBS=-L${netcdf_fortran_ROOT}/lib -lnetcdff \
     -L${jasper_ROOT}/lib64 -ljasper \
     -L${libpng_ROOT}/lib -lpng \
     -L${zlib_ROOT}/lib -lz \
     -L${bacio_ROOT}/lib -lbacio_4 \
     -L${bufr_ROOT}/lib64 -lbufr_4 \
     -L${g2_ROOT}/lib64 -lg2_4
MAP=$(map_util_dir)/misc_definitions_module.o \
    $(map_util_dir)/module_map_utils.o

# Compile individual f90 files
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(INCLD)

# Linking to create final exe
process_NSSL_mosaic.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(MAP) $(LIBS)

clean:
	rm -f *.o *.mod process_NSSL_mosaic.exe 
