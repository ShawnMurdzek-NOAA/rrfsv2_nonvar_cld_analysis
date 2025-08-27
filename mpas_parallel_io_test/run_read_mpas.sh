#!/bin/sh

#SBATCH -A wrfruc
#SBATCH -t 00:30:00
#SBATCH --ntasks=8

# Script to run read_mpas.exe

# Input Parameters
# ================
machine='ursa'
env_dir='../env/'
#mpas_file="../../test_data/mpas_south3.5km/mpasout.2024-05-08_13.00.00.nc"
mpas_file="../../test_data/mpas_conus3km/hrrrv5.history.2025-08-24_13.00.00.nc"


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
for i in $(seq 0 5); do

echo "////////////////////////////////////////////////////////////"
echo "Iteration ${i}"
echo

diff=()
exp=()
for j in $(seq -1 2); do

echo
echo "============================================================"
echo "use_mpi = ${j}"

cat << EOF > namelist.input
 &setup
  use_mpi=${j},
  fname=mpas_atm.nc
 /

EOF

s=$(date +%s%N)
srun --export=ALL read_mpas.exe
e=$(date +%s%N)
diff+=($((e-s)))
exp+=(${j})

done

# Results
echo
echo "============================================================"
for j in "${!diff[@]}"; do
  echo "Elapsed time for use_mpi ${exp[j]} = ${diff[j]} ns"
done
echo

done
