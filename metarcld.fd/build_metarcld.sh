
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
${compiler} -c ${flags} cld_parm_array_module.f90
${compiler} -c ${flags} reorg_metar_cloud_regular.f90
${compiler} -c ${flags} read_prepbufr_metarcld.f90 -I${map_util_dir}
if [[ ${machine} == 'ursa' ]]; then
  #${compiler} -c ${flags} read_prepbufr_metarcld.f90 -I${map_util_dir} -L${bufr_ROOT}/lib64 -lbufr_4
  ${compiler} ${flags} process_metar_cloud.f90 -o process_metarcld.exe \
  	  ./kinds.o ./constants.o ./mpasio.o ./cld_parm_array_module.o \
	  ./reorg_metar_cloud_regular.o ./read_prepbufr_metarcld.o \
  	  ${map_util_dir}/misc_definitions_module.o ${map_util_dir}/module_map_utils.o \
 	  -I${map_util_dir} \
	  -L${netcdf_fortran_ROOT}/lib -lnetcdff \
	  -L${bufr_ROOT}/lib64 -lbufr_4 \
	  -L${w3nco_ROOT}/lib -lw3nco_4
else
  #${compiler} -c ${flags} read_prepbufr_metarcld.f90 -I${map_util_dir} -L${bufr_ROOT}/lib64 -lbufr_d
  ${compiler} ${flags} process_metar_cloud.f90 -o process_metarcld.exe \
  	  ./kinds.o ./constants.o ./mpasio.o ./cld_parm_array_module.o \
	  ./reorg_metar_cloud_regular.o ./read_prepbufr_metarcld.o \
  	  ${map_util_dir}/misc_definitions_module.o ${map_util_dir}/module_map_utils.o \
 	  -I${map_util_dir} \
	  -L${netcdf_fortran_ROOT}/lib -lnetcdff \
	  -L${bufr_ROOT}/lib64 -lbufr_d \
	  -L${w3nco_ROOT}/lib -lw3nco_d
fi
