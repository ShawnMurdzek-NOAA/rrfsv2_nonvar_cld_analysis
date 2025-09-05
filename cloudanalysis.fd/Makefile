
# Compiler settings
# diag-disable=10448 flag is needed to disable warning about how ifort is deprecated
FFLAGS=-diag-disable=10448 -traceback

# Components of program
SRC=NonVarCldLib/map_ctp.f90 \
    NonVarCldLib/map_ctp_lar.f90 \
    NonVarCldLib/read_Surface.f90 \
    NonVarCldLib/read_radar_ref.f90 \
    NonVarCldLib/read_Lightning_cld.f90 \
    NonVarCldLib/read_nasalarc_cld.f90 \
    NonVarCldLib/vinterp_radar_ref.f90 \
    NonVarCldLib/build_missing_REFcone.f90 \
    NonVarCldLib/convert_lghtn2ref.f90 \
    NonVarCldLib/adaslib.f90 \
    NonVarCldLib/BckgrndCC.f90 \
    NonVarCldLib/cloudCover_Surface.f90 \
    NonVarCldLib/cloudCover_NESDIS.f90 \
    NonVarCldLib/cloudLayers.f90 \
    NonVarCldLib/mthermo.f90 \
    NonVarCldLib/ARPS_cldLib.f90 \
    NonVarCldLib/cloudType.f90 \
    NonVarCldLib/cloudLWC.f90 \
    NonVarCldLib/PrecipType.f90 \
    NonVarCldLib/hydro_mxr_thompson.f90 \
    NonVarCldLib/cloud_saturation.f90 \
    pbl_height.f90 \
    rapidrefresh_cldsurf_mod.f90 \
    namelist_mod.f90 \
    get_mpas_bk_mod.f90 \
    cloudanalysis_tmp_driver.f90
OBJ=${SRC:.f90=.o}

# Program dependencies
share_dir=../share
nonvar_cld_lib_dir=./NonVarCldLib
INCLD=-I${netcdf_fortran_ROOT}/include \
      -I$(share_dir) \
      -I$(nonvar_cld_lib_dir)
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
	rm -f *.o *.mod ./NonVarCldLib/*.o NonVarCLdLib/*.mod mpas_nonvarcldana.exe 
