
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
${compiler} -c ${flags} read_NASALaRC_cloud.f90
${compiler} -c ${flags} write_bufr_NASALaRC.f90
if [[ ${machine} == 'ursa' ]]; then
  echo 'success'
  ${compiler} ${flags} process_NASALaRC_cloud.f90 -o process_larccld.exe \
  	  ./kinds.o ./constants.o ./mpasio.o ./read_NASALaRC_cloud.o ./write_bufr_NASALaRC.o \
  	  ${map_util_dir}/misc_definitions_module.o ${map_util_dir}/module_map_utils.o \
 	  -I${map_util_dir} \
	  -L${netcdf_fortran_ROOT}/lib -lnetcdff \
	  -L${bufr_ROOT}/lib64 -lbufr_4
else
  ${compiler} ${flags} process_NASALaRC_cloud.f90 -o process_larccld.exe \
  	  ./kinds.o ./constants.o ./mpasio.o ./read_NASALaRC_cloud.o ./write_bufr_NASALaRC.o \
  	  ${map_util_dir}/misc_definitions_module.o ${map_util_dir}/module_map_utils.o \
 	  -I${map_util_dir} \
	  -L${netcdf_fortran_ROOT}/lib -lnetcdff \
	  -L${bufr_ROOT}/lib64 -lbufr_d
fi
