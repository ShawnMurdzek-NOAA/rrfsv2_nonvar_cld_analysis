#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 00:30:00
#SBATCH --ntasks=1
#SBATCH --mem=20G

# Script to run process_larccld.exe for FV3 LAM for a single time

# Input Parameters
# ================
obs_dir="/scratch2/BMC/zrtrr/RRFS_RETRO_DATA/obs_rap"
rrfsv1_home="/scratch2/BMC/wrfruc/murdzek/nonvar_cld_analysis_testing/RRFSv1"
rrfsv1_code_dir="${rrfsv1_home}/rrfs-workflow"
mpas_mesh_file="/scratch2/BMC/wrfruc/murdzek/nonvar_cld_analysis_testing/RRFSv2/test_data/conus3km.grid.nc"
valid=2024050813


# Run Program
# =========================

# Load environment
echo "Configuring environment"
module purge
module use ${rrfsv1_code_dir}/modulefiles
module load build_hera_intel
module list
echo

#ulimit -s unlimited
#ulimit -a

# Get input data
cp ${rrfsv1_code_dir}/fix/gsi/prepobs_prep_RAP.bufrtable prepobs_prep.bufrtable
cp ${obs_dir}/${valid}.rap.t${valid:8:2}z.lgycld.tm00.bufr_d lgycld.bufr_d
cp ${mpas_mesh_file} mesh.nc

# Create namelist
cat << EOF > namelist.nasalarc
 &setup
  analysis_time = ${valid},
  bufrfile='NASALaRCCloudInGSI_bufr.bufr',
  rad_km=9,
  ioption = 2,
 /
EOF

# Run program
date
echo
srun --export=ALL process_larccld.exe
echo
date
