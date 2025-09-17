# NASA LaRC Preprocessing Program

This program takes satellite-based cloud top observations from NASA and interpolates them to the MPAS mesh.

## Compiling

1. Load the proper environment from `../env/`
2. `make clean && make`

## Running

To run larccld.fd by itself, one can use the sample slurm submission script provided in this directory (`run_process_NASA_LaRC.sh`). The program is designed to run on a single processor, but there is some MPI code included so that if multiple processors are used to run the program, it will only run on the root processor.

### Inputs

- `lgycld.bufr_d`: NASA LaRC obs in BUFR format.
- `prepobs_prep.bufrtable`: BUFR table (found in `../fix/prepobs_prep_RAP.bufrtable`).
- `mesh.nc`: An MPAS mesh file that includes the (lat, lon) coordinates of each MPAS cell.
- `namelist.nasalarc`: Namelist file (see options below).

### Outputs

- `NASALaRC_cloud4mpas.bin`: Binary file with NASA LaRC data mapped to the MPAS mesh. This is the file used by cloudanalysis.fd.
- `NASALaRCCloudInGSI_bufr.bufr`: Data as the previous bullet, but in BUFR format.
- `nasa_larc_obs_raw.txt`: Text file with the raw NASA LaRC observations (if `debug = 1`).
- `nasa_larc_obs_interp.txt`: Text file with the NASA LaRC observations interpolated to MPAS mesh (if `debug = 1`).

### Namelist Options

All namelist options are in a single section titles `&setup`

#### Scalar Options

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| `analysis_time` | 2018051718 | Analysis valid time in YYYYMMDDHH format. |
| `bufrfile` | NASALaRCCloudInGSI\_bufr.bufr | Name of output BUFR file. |
| `ioption` | 2 | Interpolation option. 1 = nearest neighbor, 2 = median. Only ioption = 2 has been tested. |
| `npts_rad` | 1 | Half length of the square box used to interpolate NASA LaRC observations (see "code overview" section). Units are grid boxes of the map projection. |
| `userDX` | 3000. | Model mesh spacing in meters |
| `proj_name` | CONUS | Map projection to use. Must be defined in `../share/map_proj_helper_mod.f90`. |
| `debug` | 0 | Option to print additional output for debugging. Set to 0 to not print any additional output |

#### Vector Options

The following parameters are vectors with length boxMAX = 10. By default, all values in these vectors are the same. The three parameters listed here control the configuration of the variable box used for interpolation, which is turned off by default. Additional details are found in the "code overview" section.

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| `boxlat0` | 999. | Edges of the latitude bins (in deg N) used for the variable interpolation box. Set to a large number (e.g., 999.) to not use a variable box. |
| `boxhalfx` | -1 | Number of map projection grid points in the x direction for the variable interpolation box for each latitude bin defined by boxlat0. |
| `boxhalfy` | -1 | Number of map projection grid points in the y direction for the variable interpolation box for each latitude bin defined by boxlat0. |

## Code Overview

Below is a brief explanation of each of the fortran files used in this program.

### larccld\_bufr\_io\_mod.f90

Module with various subroutines for reading NASA LaRC observations in BUFR format and writing processed NASA LaRC observations to a BUFR file.

### larccld\_utils\_mod.f90

Module with various utilities used by larccld.fd.

### process\_NASALaRC\_cloud.f90

This is the main driver for this program. The steps followed by the driver are as follows:

1. Define the map projection used as an intermediary between the NASA LaRC observations and MPAS mesh.
2. Read in the MPAS mesh information and NASA LaRC observations
3. Interpolate the NASA LaRC observations to the MPAS mesh using the map projection
    1. Define a box around each point in the map projection. The box size is defined by either `npts_rad` for a constant-sized box or `boxlat0`, `boxhalfx`, and `boxhalfy` for a variable-sized box.
    2. Determine all the observations that fall into each box. The observations in each box are saved in arrays that contain "xx" in the name.
    3. Match each MPAS cell with a map projection grid point and corresponding box. The NASA LaRC observations interpolated to that MPAS cell is a reduction of all observations within the corresponding box, with the reduction determined by `ioption`.
4. Write out results to a binary file and BUFR file.

### process\_NASALaRC\_cloud\_FV3.f90

Older version of process\_NASALaRC\_cloud.f90 for the FV3 dynamical core. This file is included here for reference and is not compiled as part of the program.

### read\_NASALaRC\_cloud.f90

Includes a single subroutine for reading NASA LaRC observations in netCDF format. This file is included here for reference (as it was part of the FV3 version of larccld.fd), but is not compiled as part of the program. 
