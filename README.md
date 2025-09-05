# Nonvariational Cloud Analysis for MPAS-Based RRFSv2

The nonvariational cloud analysis is an algorithm for updating the model state based on cloud and precipitation observations from ceilometers, satellite, radar, and lightning detection networks. The ultimate goal of this analysis is to better represent clouds in the model initial conditions and improve forecasts of cloud ceilings. More detailed information about the design and value of the nonvariational cloud analysis can be found in [Benjamin et al. (2021)](https://doi.org/10.1175/MWR-D-20-0319.1)

This code is based off of [rrfs-utl](https://github.com/NOAA-GSL/rrfs_utl), which contains the nonvariational cloud analysis for the FV3-based RRFSv1. Notable differences between the FV3 and MPAS implementations can be found in the "Comparison to FV3" section below. All files that have been altered to accomodate MPAS have a corresponding unaltered file that has "FV3" or "fv3" in the name. Although these files are still included here, they are excluded from the Makefiles, and are, therefore, not compiled when building the programs.

## Organization

The nonvariational cloud analysis consists of 5 primary programs:

- `larccld.fd`: Interpolate NASA LaRC cloud tops to MPAS mesh
- `lightning.fd`: Interpolate lightning observations to MPAS mesh
- `metarcld.fd`: Interpolate ceilometer observations to MPAS mesh
- `mosaic.fd`: Interpolate NSSL MRMS reflectivity to MPAS mesh
- `cloudanalysis.fd`: Main program for performing the nonvariational cloud analysis

Some code is shared by all or some of the 5 aforementioned programs. That code is found in `share`. Code in `share` must be compiled first before attempting to compile any of the nonvariational cloud analysis programs. This is automatically handled by `build_all.sh`.

There is also some code for visualization and sandbox testing:

- `test_plots`: Includes various Python scripts for visualizing output from the various nonvariational cloud analysis programs
- `mpas_parallel_io_test`: Sandbox for testing parallel I/O with MPAS

The remaining directories are as follows:

- `env`: Files containing the build environment for various NOAA RDHPCS machines
- `fix`: Various static files used by the nonvariational cloud analysis

## Building

0. Add an `.env` file to `./env/` if there is not one for the machine you are working on
1. Edit `build_all.sh` to include the correct machine
2. Run `bash build_all.sh`

## Running

0. Follow the above steps to build all 5 programs
1. `cd run`
2. Edit the top of `run_mpas_nonvar_cld_analysis.sh` to have the proper machine and slurm specs (e.g., allocation)
3. `sbatch run_mpas_nonvar_cld_analysis.sh`

## Comparison to FV3

The chief difference between MPAS and FV3 that is relevant for the nonvariational cloud analysis is that FV3 uses a 3D dimensional grid (x, y, z) for 3D fields whereas MPAS uses a 2D dimensional grid (cell, z). In an effort to minimize the number of files that need to be changed when adapting the nonvariational cloud analysis for MPAS, many of the arrays in `cloudanalysis.fd` are still 3D with dimensions `(lon2, lat2, nsig)`. For the MPAS implementation here, `lon2 = 1, lat2 = nCell, nsig = nz`. Thus, these 3D arrays are effectively 2D.

Another notable difference is that there are several instances in the nonvar cloud analysis for FV3 where it is assumed that there is a ring of buffer points around the 2D spatial domain (e.g., loops over the x and y dimensions go from 2 to lon2-1). There are no buffer points in the MPAS domain, so these loops are changed so that they go from 1 to lon2 or lat2.
