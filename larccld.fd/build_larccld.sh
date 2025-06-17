
rrfsv1_code_dir="/scratch4/BMC/wrfruc/murdzek/nonvar_cld_analysis_testing/RRFSv1/rrfs-workflow"
map_util_dir='../wps_map_utils/'

echo 'Cleaning...'
rm *.o
rm *.mod

echo 'Configuring Environment...'
module purge
module use ${rrfsv1_code_dir}/modulefiles
module load build_hera_intel
module list
echo

echo 'Building Program...'
compiler=mpiifort
${compiler} -c kinds.f90
${compiler} -c constants.f90
${compiler} -c mpasio.f90 -I ${netcdf_fortran_ROOT}/include
${compiler} -c read_NASALaRC_cloud.f90
${compiler} -c write_bufr_NASALaRC.f90
${compiler} process_NASALaRC_cloud_NEW.f90 -o process_larccld.exe \
	./kinds.o ./constants.o ./mpasio.o ./read_NASALaRC_cloud.o ./write_bufr_NASALaRC.o \
	${map_util_dir}/misc_definitions_module.o ${map_util_dir}/module_map_utils.o \
	-I${map_util_dir} \
	-L${netcdf_fortran_ROOT}/lib -lnetcdff \
	-L${bufr_ROOT}/lib64 -lbufr_d
