
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
flags='-diag-disable=10448 -traceback -g -check all'
${compiler} -c ${flags} kinds.f90
${compiler} -c ${flags} mpasio.f90 -I ${netcdf_fortran_ROOT}/include
${compiler} -c ${flags} read_lightning_bufr.f90
${compiler} -c ${flags} Check_Lightning_QC.f90
${compiler} -c ${flags} netCDFsub_lightning.f90 -I ${netcdf_fortran_ROOT}/include
${compiler} -c ${flags} Check_NLDN.f90
${compiler} -c ${flags} write_bufr_lght.f90

if [[ ${machine} == 'ursa' ]]; then
  ${compiler} ${flags} process_Lightning.f90 -o process_Lightning.exe \
  	  ./kinds.o ./mpasio.o \
	  ./read_lightning_bufr.o ./Check_Lightning_QC.o \
	  ./netCDFsub_lightning.o ./Check_NLDN.o ./write_bufr_lght.f90 \
  	  ${map_util_dir}/misc_definitions_module.o ${map_util_dir}/module_map_utils.o \
 	  -I${map_util_dir} \
	  -L${netcdf_fortran_ROOT}/lib -lnetcdff \
	  -L${bufr_ROOT}/lib64 -lbufr_4 \
          -L${w3nco_ROOT}/lib -lw3nco_4
else
  echo "${machine} is not supported yet. Skipping linking step"
fi
