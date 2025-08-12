#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 00:30:00
#SBATCH --ntasks=1
#SBATCH --mem=20G

# Script to run process_Lightning.exe for MPAS for a single time

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
cp ${obs_dir}/${valid}.rap.t${valid:8:2}z.lghtng.tm00.bufr_d lghtngbufr
cp ${mpas_mesh_file} mesh.nc

# Create namelist
cat << EOF > namelist.lightning
 &setup
  analysis_time = ${valid},
  minute=00,
  trange_start=-10,
  trange_end=10,
  obs_type = "bufr",
  debug=1
 /

EOF

# Clean old text files
clean_files=( 'LightningInGSI_bufr.bufr' 'LightningInMPAS.dat' 
              'lightning_interp.txt' 'lightning_raw_1.txt' )
for f in ${clean_files[@]}; do
  if [[ -f ${f} ]]; then
    rm ${f}
  fi
done

# Run program
date
echo
srun --export=ALL process_Lightning.exe
echo
date
