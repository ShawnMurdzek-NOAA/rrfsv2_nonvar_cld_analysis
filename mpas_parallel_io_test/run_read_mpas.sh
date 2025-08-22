#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 00:10:00
#SBATCH --ntasks=8

# Script to run read_mpas.exe

# Input Parameters
# ================
machine='ursa'
env_dir='../env/'
mpas_file="/scratch4/BMC/wrfruc/murdzek/nonvar_cld_analysis_testing/RRFSv2/test_data/mpas_south3.5km/mpasout.2024-05-08_13.00.00.nc"


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

# Get input data
cp ${mpas_file} ./mpas_atm.nc

# Run program a few different ways
diff=()
for i in $(seq 0 2); do

echo
echo "============================================================"
echo "use_mpi = ${i}"
echo

cat << EOF > namelist.input
 &setup
  use_mpi=${i},
  fname=mpas_atm.nc
 /

EOF

s=$(date +%s%N)
srun --export=ALL read_mpas.exe
e=$(date +%s%N)
diff+=($((e-s)))

done

# Results
echo
echo "============================================================"
for i in "${!diff[@]}"; do
  echo "Elapsed time for use_mpi ${i} = ${diff[i]} ns"
done
echo
