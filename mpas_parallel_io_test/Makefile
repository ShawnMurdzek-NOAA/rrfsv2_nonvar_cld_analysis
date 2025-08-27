
# Compiler settings
# diag-disable=10448 flag is needed to disable warning about how ifort is deprecated
FFLAGS=-diag-disable=10448 -traceback

# Components of program
SRC=read_mpas_test.f90
OBJ=${SRC:.f90=.o}

# Program dependencies
INCLD=-I${netcdf_fortran_ROOT}/include
LIBS=-L${netcdf_fortran_ROOT}/lib -lnetcdff

# Compile individual f90 files
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(INCLD)

# Linking to create final exe
read_mpas.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIBS)

clean:
	rm -f *.o *.mod read_mpas.exe 
