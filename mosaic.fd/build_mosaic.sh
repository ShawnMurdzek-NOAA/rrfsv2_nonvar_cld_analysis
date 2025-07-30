
env_dir='../env/'
map_util_dir='../wps_map_utils/'
machine='ursa'

echo 'Cleaning...'
rm *.o
rm *.mod

echo 'Configuring Environment...'
source ${env_dir}/${machine}.env
module list
echo

echo 'Building Program...'
compiler=mpiifort

# diag-disable flag disables a warning about how ifort is deprecated
flags='-diag-disable=10448 -traceback'
${compiler} -c ${flags} kinds.f90
${compiler} -c ${flags} constants.f90
${compiler} -c ${flags} mpasio.f90 -I ${netcdf_fortran_ROOT}/include
${compiler} -c ${flags} netCDFsub.f90 -I ${netcdf_fortran_ROOT}/include
${compiler} -c ${flags} read_ncep_binary.f90
${compiler} -c ${flags} read_nssl_binary.f90
${compiler} -c ${flags} read_grib2_mod.f90 -I ${g2_ROOT}/include_4  # Also need to add to linking step below
${compiler} -c ${flags} module_read_NSSL_mosaic.f90
${compiler} -c ${flags} write_netcdf_ref.f90 -I ${netcdf_fortran_ROOT}/include
${compiler} -c ${flags} write_bufr_ref.f90

# NCEPLIBS-g2 (required for GRIB2 I/O) requires the following: jasper, libpng, zlib, NCEPLIBS-bacio
if [[ ${machine} == 'ursa' ]]; then
  ${compiler} ${flags} process_NSSL_mosaic.f90 -o process_NSSL_mosaic.exe \
  	  ./kinds.o ./constants.o ./mpasio.o ./netCDFsub.o \
	  ./read_ncep_binary.o ./read_nssl_binary.o ./read_grib2_mod.o \
	  ./module_read_NSSL_mosaic.o ./write_netcdf_ref.o ./write_bufr_ref.f90 \
  	  ${map_util_dir}/misc_definitions_module.o ${map_util_dir}/module_map_utils.o \
 	  -I${map_util_dir} \
	  -L${netcdf_fortran_ROOT}/lib -lnetcdff \
	  -L${jasper_ROOT}/lib64 -ljasper \
	  -L${libpng_ROOT}/lib -lpng \
	  -L${zlib_ROOT}/lib -lz \
	  -L${bacio_ROOT}/lib -lbacio_4 \
	  -L${bufr_ROOT}/lib64 -lbufr_4 \
	  -L${g2_ROOT}/lib64 -lg2_4
else
  echo "${machine} is not supported yet. Skipping linking step"
fi
