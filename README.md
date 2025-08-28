# Nonvariational Cloud Analysis for MPAS-Based RRFSv2

The nonvariational cloud analysis is an algorithm for updating the model state based on cloud and precipitation observations from ceilometers, satellite, radar, and lightning detection networks. The ultimate goal of this analysis is to better represent clouds in the model initial conditions and improve forecasts of cloud ceilings. More detailed information about the design and value of the nonvariational cloud analysis can be found in [Benjamin et al. (2021)](https://doi.org/10.1175/MWR-D-20-0319.1)

This code is based off of [rrfs-utl](https://github.com/NOAA-GSL/rrfs_utl), which contains the nonvariational cloud analysis for the FV3-based RRFSv1.

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
