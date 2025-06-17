
rrfsv1_code_dir="/scratch4/BMC/wrfruc/murdzek/nonvar_cld_analysis_testing/RRFSv1/rrfs-workflow"

echo 'Cleaning...'
rm *.o
rm *.mod
rm *.exe

echo 'Configuring Environment...'
module purge
module use ${rrfsv1_code_dir}/modulefiles
module load build_hera_intel
module list
echo

echo 'Building Program...'
compiler=mpiifort
${compiler} -c constants_module_wps.f90
${compiler} -c misc_definitions_module.f90
${compiler} -c module_map_utils.f90
${compiler} test_map_utils.f90 -o test_map_utils.exe \
	./misc_definitions_module.o ./module_map_utils.o
