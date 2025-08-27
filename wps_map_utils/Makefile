
# Compiler settings
# diag-disable=10448 flag is needed to disable warning about how ifort is deprecated
FFLAGS=-diag-disable=10448 -traceback

# Components of program
SRC=constants_module_wps.f90 \
    misc_definitions_module.f90 \
    module_map_utils.f90 \
    test_map_utils.f90
OBJ=${SRC:.f90=.o}

# Compile individual f90 files
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(INCLD)

# Linking to create final exe
test_map_utils.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(MAP) $(LIBS)

clean:
	rm -f *.o *.mod test_map_utils.exe
