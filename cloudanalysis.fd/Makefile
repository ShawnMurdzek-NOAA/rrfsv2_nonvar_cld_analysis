
# Compiler settings
# diag-disable=10448 flag is needed to disable warning about how ifort is deprecated
FFLAGS=-diag-disable=10448 -traceback

# Components of program
SRC=pbl_height.f90 \
    rapidrefresh_cldsurf_mod.f90 \
    namelist_mod.f90 \
    get_mpas_bk_mod.f90 \
    cloudanalysis_tmp_driver.f90
OBJ=${SRC:.f90=.o}

# Program dependencies
share_dir=../share
INCLD=-I${netcdf_fortran_ROOT}/include \
      -I$(share_dir)
LIBS=-L${netcdf_fortran_ROOT}/lib -lnetcdff
SHARE=$(share_dir)/kinds.o \
      $(share_dir)/constants.o \
      $(share_dir)/mpasio.o

# Compile individual f90 files
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(INCLD)

# Linking to create final exe
mpas_nonvarcldana.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(SHARE) $(LIBS)

clean:
	rm -f *.o *.mod mpas_nonvarcldana.exe 
