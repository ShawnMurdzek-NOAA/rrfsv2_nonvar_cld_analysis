
module load gnu

# Compile programs
gfortran compute_haversine.f90 -o haversine.exe
gfortran compute_euclidean_dist.f90 -o euclidean.exe

# Run each program
echo 'Haversine'
time ./haversine.exe

echo
echo 'Euclidean'
time ./euclidean.exe
