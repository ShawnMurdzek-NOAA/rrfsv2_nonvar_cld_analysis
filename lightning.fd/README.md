# Lightning Observation Preprocessing Program

This program takes lightning observations and interpolates them to the MPAS mesh.

## Compiling

1. Load the proper environment from `../env/`
2. `make clean && make`

## Running

To run lightning.fd by itself, one can use the sample slurm submission script provided in this directory (`run_process_Lightning.sh`). The program is designed to run on a single processor, but there is some MPI code included so that if multiple processors are used to run the program, it will only run on the root processor.

### Inputs

- `lghtngbufr`: Lightning observations in BUFR format.
- `prepobs_prep.bufrtable`: BUFR table (found in `../fix/prepobs_prep_RAP.bufrtable`).
- `mesh.nc`: An MPAS mesh file that includes the (lat, lon) coordinates of each MPAS cell.
- `namelist.lightning`: Namelist file (see options below).

### Outputs

- `LightningInMPAS.dat`: Binary file with lightning observations mapped to the MPAS mesh. This is the file used by cloudanalysis.fd.
- `LightningInGSI_bufr.bufr`: Same data as the previous bullet, but in BUFR format. Can also use lightning observations in netCDF format, but that has not been tested.
- `lightning_raw_N.txt`: Text file with the raw lightning observations from file number N (if `debug = 1`).
- `lightning_interp.txt`: Text file with the lightning observations interpolated to MPAS mesh (if `debug = 1`).

### Namelist Options

All namelist options are in a single section titles `&setup`

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| `analysis_time` | N/A | Analysis valid time in YYYYMMDDHH format. |
| `minute` | 0 | Minute of the analysis time. |
| `trange_start` | 0 | Lower bound of observation temporal processing window (minutes). Only observations with a timestamp > (obs time + trange\_start) are used. |
| `trange_end` | 0 | Upper bound of observation temporal processing window (minutes). Only observations with a timestamp < (obs time + trange\_end) are used. |
| `obs_type` | "none" | Observation file format. Options: "bufr" or "nldn\_nc". Only "bufr" has been tested. |
| `proj_name` | CONUS | Map projection to use. Must be defined in `../share/map_proj_helper_mod.f90`. |
| `debug` | 0 | Option to print additional output for debugging. Set to 0 to not print any additional output. |

## Code Overview

Below is a brief explanation of each of the fortran files used in this program.

### Check\_Lght\_QC\_mod.f90

Subroutines for lightning observation quality control checks.

### lightning\_bufr\_io\_mod.f90

Subroutines for reading lightning observations in BUFR format and writing lightning observations processed to the MPAS mesh in BUFR format.

### netCDFsub\_lightning\_mod.f90

Subroutunes for reading lightning observations in netCDF format.

### process\_Lightning.f90

This is the main driver for this program. The steps followed by the driver are as follows:

1. Define the map projection used as an intermediary between the lightning observations and MPAS mesh.
2. Read in MPAS mesh information and map each MPAS cell to the closest map projection integer coordinate.
3. Read lightning observations and perform QC check.
4. Map lightning observations to MPAS mesh:
    1. Match each lightning observation with the closest map projection integer coordinate.
    2. Loop over each map projection integer coordinate. If there are lightning observations associated with that coordinate, compute the distance between the observation and all MPAS cells associated with a 3x3 box of map projection integer coordinates centered on the coordinate of interest. Assign the observation to the closest MAPS cell.
5. Write out results.  

### process\_Lightning\_FV3.f90

Older version of process\_Lightning.f90 for the FV3 dynamical core. This file is included here for reference and is not compiled as part of the program.
