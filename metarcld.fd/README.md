# Ceilometer Preprocessing Program

This program takes METAR ceilometer observations, adds a radius of impact, thens interpolates them to the MPAS mesh.

## Compiling

1. Load the proper environment from `../env/`
2. `make clean && make`

## Running

To run metarcld.fd by itself, one can use the sample slurm submission script provided in this directory (`run_process_metarcld.sh`). The program is designed to run on a single processor, but there is some MPI code included so that if multiple processors are used to run the program, it will only run on the root processor.

### Inputs

- `prepbufr`: Ceilometer obs in BUFR format.
- `prepobs_prep.bufrtable`: BUFR table (found in `../fix/prepobs_prep_RAP.bufrtable`).
- `mesh.nc`: An MPAS mesh file that includes the (lat, lon) coordinates of each MPAS cell.
- `namelist.metarcld`: Namelist file (see options below).

### Outputs

- `mpas_metarcloud.bin`: Binary file with NASA LaRC data mapped to the MPAS mesh. This is the file used by cloudanalysis.fd.
- `processed_metar_obs_mpas.txt`: Text file with the ceilometer observations processed to the MPAS mesh (if `debug = 1`).

### Namelist Options

All namelist options are in a single section titles `&setup`

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| `analysis_time` | 2018051718 | Analysis valid time in YYYYMMDDHH format. |
| `analysis_minute` | 0 | Valid minute of analysis time. |
| `prepbufrfile` | prepbufr | BUFR file containing ceilometer observations. |
| `twindin` | 0.5 | Only use ceilometer observations that are within `twindin` hours of the analysis time. |
| `metar_impact_radius` | 16. | Radius of impact for each ceilometer observation, in grid boxes of the map projection (so radius in meters = `(metar_impact_radius) X (grid spacing of map projection)`) |
| `l_metar_impact_radius_change` | .false. | Option to use a value for `metar_impact_radius` that increases with height. |
| `metar_impact_radius_max` | 50000. | Maximum ceilometer radius of impact when `l_metar_impact_radius = .true.` (m). |
| `metar_impact_radius_min` | 20000. | Minimum ceilometer radius of impact when `l_metar_impact_radius = .true.` (m). |
| `metar_impact_radius_max_height` | 3000. | Height above which the ceilometer radius of impact is `metar_impact_radius_max` when `l_metar_impact_radius = .true.` (m). |
| `metar_impact_radius_min_height` | 200. | Height below which the ceilometer radius of impact is `metar_impact_radius_min` when `l_metar_impact_radius = .true.` (m). |
| `region_dx` | 3000. | Model mesh spacing in meters |
| `debug` | 0 | Option to print additional output for debugging. Set to 0 to not print any additional output |

## Code Overview

Below is a brief explanation of each of the fortran files used in this program.

### cld\_parm\_array\_module.f90  

Module that defines a subset of the parameters found in the program namelist.

### process\_metar\_cloud.f90    

This is the main driver for this program. It reads the MPAS mesh information, calls the subroutine in read\_prepbufr\_metar\_mod.f90, then writes the results to a binary file.

### process\_metar\_cloud\_FV3.f90	   

Older version of process\_metar\_cloud.f90 for the FV3 dynamical core. This file is included here for reference and is not compiled as part of the program.

### read\_prepbufr\_metar\_mod.f90	

This module contains a single subroutine that:

1. Defines the map projection used as an intermediary between the ceilometer observations and MPAS mesh. 
2. Reads ceilometer observations from the BUFR file.
3. Determines ceilometer observation locations in the desired map projection.
4. Writes all desired ceilometer fields to the array `cdata_out`.
5. Determines the MPAS cell locations in the desired map projection.
6. Calls the subroutine in reorg\_metar\_cloud\_mod.f90.

### read\_prepbufr\_metarcld\_FV3.f90  

Older version of read\_prepbufr\_metar\_mod.f90 for the FV3 dynamical core. This file is included here for reference and is not compiled as part of the program.

### reorg\_metar\_cloud\_mod.f90

This module contains a single subroutine that:

1. Performs various quality control checks of the ceilometer observations (QC flags are saved to `cdata(22,:)`).
2. Determines the closest map projection grid point for each ceilometer observation and MPAS cell.
3. Creates a halo around each ceilometer observation (with radius `metar_impact_radius`) and interpolates to the MPAS mesh:
    1. Loop over a grid coarser than the map projection with a grid size of `metar_impact_radius`.
    2. For each grid box of the grid specified in (i), find all ceilometers with an encompassing box centered on the grid box with length `3 X metar_impact_radius`.
    3. For each MPAS cell within a given grid box of the grid specified in (i), find the closest ceilometer within the encompassing box defined in (ii). If the distance is less than `metar_impact_radius`, then assign the observations from the closest ceilometer to the MPAS cell.


### reorg\_metar\_cloud\_regular\_FV3.f90

Older version of reorg\_metar\_cloud\_mod.f90 for the FV3 dynamical core. This file is included here for reference and is not compiled as part of the program.
