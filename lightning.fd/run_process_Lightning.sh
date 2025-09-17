#!/bin/sh

#SBATCH -M c6
#SBATCH -A bil-pmp
#SBATCH -t 00:10:00
#SBATCH --ntasks=1
#SBATCH --qos=debug

# Script to run process_Lightning.exe for MPAS for a single time

# Input Parameters
# ================
machine='gaeaC6'
env_dir='../env/'
fix_dir='../fix/'
obs_dir="../../test_data/obs"
mpas_mesh_file="../../test_data/mpas_south_3.5km.grid.nc"
valid=2024050813
proj_name='CONUS'

# Build Program
# =========================

# Load environment
echo "Configuring environment"
source ${env_dir}/${machine}.env
module list
echo

# Build program
make clean
make


# Run Program
# =========================

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
  proj_name = '${proj_name}',
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
