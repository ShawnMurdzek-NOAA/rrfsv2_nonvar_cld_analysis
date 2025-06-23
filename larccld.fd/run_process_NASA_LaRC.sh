#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 00:30:00
#SBATCH --ntasks=1
#SBATCH --mem=20G

# Script to run process_larccld.exe for FV3 LAM for a single time

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
cp ${obs_dir}/${valid}.rap.t${valid:8:2}z.lgycld.tm00.bufr_d lgycld.bufr_d
cp ${mpas_mesh_file} mesh.nc

# Create namelist
cat << EOF > namelist.nasalarc
 &setup
  analysis_time = ${valid},
  bufrfile='NASALaRCCloudInGSI_bufr.bufr',
  npts_rad=3,
  ioption = 2,
 /
EOF

# Clean old text files
clean_files=( 'NASALaRC_cloud4fv3.bin' 'nasa_larc_obs_interp.txt' 'nasa_larc_obs_raw.txt' )
for f in ${clean_files[@]}; do
  if [[ -f ${f} ]]; then
    rm ${f}
  fi
done

# Run program
date
echo
srun --export=ALL process_larccld.exe
echo
date
