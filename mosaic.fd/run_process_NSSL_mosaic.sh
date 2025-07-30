#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 00:30:00
#SBATCH --ntasks=48
#SBATCH --nodes=2

# Script to run process_NSSL_mosaic.exe for MPAS for a single time

# Input Parameters
# ================
machine='ursa'
env_dir='../env/'
fix_dir='../fix/'
obs_dir="/scratch4/BMC/wrfruc/murdzek/nonvar_cld_analysis_testing/RRFSv2/test_data/obs"
mpas_mesh_file="/scratch4/BMC/wrfruc/murdzek/nonvar_cld_analysis_testing/RRFSv2/test_data/conus3km.grid.nc"
valid=2024050813


# Run Program
# =========================

# Load environment
echo "Configuring environment"
source ${env_dir}/${machine}.env
module list
echo

#ulimit -s unlimited
#ulimit -a

# Get input data
cp ${fix_dir}/prepobs_prep_RAP.bufrtable prepobs_prep.bufrtable
cp ${obs_dir}/reflectivity/MergedReflectivityQC*${valid::8}-${valid:8:2}0* .
cp ${mpas_mesh_file} mesh.nc
ls MergedReflectivityQC* > filelist_mrms

# Create namelists
cat << EOF > namelist.mosaic
   &setup
    tversion=1,
    analysis_time = ${valid},
    dataPath = './',
   /
EOF

cat << EOF > namelist.mosaic_netcdf
   &setup_netcdf
    output_netcdf = .true.,
    max_height = 11001.0,
    use_clear_air_type = .true.,
    precip_dbz_thresh = 10.0,
    clear_air_dbz_thresh = 5.0,
    clear_air_dbz_value = 0.0,
    precip_dbz_horiz_skip = 0,
    precip_dbz_vert_skip = 0,
    clear_air_dbz_horiz_skip = 0,
    clear_air_dbz_vert_skip = 0,
    remove_bdy = .false.,
   /
EOF

# Clean old text files
clean_files=( 'Gridded_ref.nc' 'RefInGSI3D.dat' )
for f in ${clean_files[@]}; do
  if [[ -f ${f} ]]; then
    rm ${f}
  fi
done

# Run program
date
echo
srun --export=ALL process_NSSL_mosaic.exe
echo
date
