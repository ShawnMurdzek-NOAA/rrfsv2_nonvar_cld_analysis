# Radar Reflectivity Mosaic Preprocessing Program

This program takes radar reflectivity mosaic observations (e.g., MRMS) and interpolates them to the MPAS mesh.

## Compiling

1. Load the proper environment from `../env/`
2. `make clean && make`

## Running

To run mosaic.fd by itself, one can use the sample slurm submission script provided in this directory (`run_process_NSSL_mosaic.sh`). The program is designed to run on multiple processors, and depending on the input data type, there are a minimum number of processors that must be used:

| tversion | Description | Processors |
| -------- | ----------- | ---------- |
| 1 | NSSL 1 tile GRIB2 | 33+ |
| 4 | NSSL 4 tiles binary | 4+ |
| 81 | NCEP 8 tiles binary | 8+ |
| 8 | NSSL 8 tiles netCDF | 8+ |

### Inputs

- `MergedReflectivityQC_XX.XX_YYYYMMDD-HHMMSS.grib2`: MRMS radar reflectivity observations. Each file contains 1 vertical level at `XX.XX` km AGL. `YYYYMMDD-HHMMSS` is the timestamp for the file.
  - This is the `tversion = 1` option referenced above. Using other observation formats (e.g., NSSL 4 tiles binary) has not been tested.
- `prepobs_prep.bufrtable`: BUFR table (found in `../fix/prepobs_prep_RAP.bufrtable`).
- `mesh.nc`: An MPAS mesh file that includes the (lat, lon) coordinates of each MPAS cell.
- `filelist_mrms`: Text file containing the names of all MRMS radar reflectivity files.
- `namelist.mosaic`: Namelist file (see options below).
- `namelist.mosaic_netcdf`: Namelist file (see options below).

### Outputs

- `RefInGSI3D.dat`: Binary file with radar reflectivity observations mapped to the MPAS mesh. This is the file used by cloudanalysis.fd.
- `Gridded_ref.nc`: NetCDF file with radar reflectivity observations mapped to the MPAS mesh. Only produced if `output_netcdf = .true.` in `namelist.mosaic_netcdf`.

### Namelist Options

#### namelist.mosaic

All namelist options are in a single section titled `&setup`.

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| `tversion` | N/A | Date format for the radar reflectivity observations (see table above). |
| `analysis_time` | N/A | Analysis valid time in YYYYMMDDHH format. |
| `dataPath` | N/A | Path to radar reflectivity observations (only used if `tversion = 4, 8, or 81`).  |

#### namelist.mosaic\_netcdf

All namelist options are in a single section titled `&setup_netcdf`. These options only apply to `Gridded_ref.nc` and do not impact the output dumped to `RefInGSI3D.dat`.

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| `output_netcdf` | .false. | Option to output a netCDF file. |
| `max_height` | 20000. | Maximum height (m MSL) to retain observations. |
| `use_clear_air_type` | .false. | Option to output clear-air reflectivity observations. |
| `precip_dbz_thresh` | 15. | Minimum reflectivity threshold that is considered precipitation (dBZ). |
| `clear_air_dbz_thresh` | 0. | Maximum reflectivity threshold that is considered clear air (dBZ). |
| `clear_air_dbz_value` | 0. | Value to assign to clear-air reflectivity observations (dBZ). |
| `precip_dbz_horiz_skip` | 0 | Horizontal thinning factor for precipitation observations. Hardcoded to always be 0. Thinning is commented out in the program. |
| `precip_dbz_vert_skip` | 0 | Vertical thinning factor for precipitation observations. Hardcoded to always be 0. Thinning is commented out in the program. |
| `clear_air_dbz_horiz_skip` | 0 | Horizontal thinning factor for clear-air observations. Hardcoded to always be 0. Thinning is commented out in the program. |
| `clear_air_dbz_vert_skip` | 0 | Vertical thinning factor for clear-air observations. Hardcoded to always be 0. Thinning is commented out in the program. |
| `remove_bdy` | .false. | Option to set reflectivity values along the boundary to missing (-999). |

## Code Overview

Below is a brief explanation of each of the fortran files used in this program.

### module\_read\_NSSL\_mosaic.f90  

Module that defines a derived type and associated subroutines for reading NSSL MRMS reflectivity observations.

### mosaic\_interp\_mod.f90

Module that defines a single subroutine for interpolating gridded radar reflectivity observations to an array of model (latitude, longitude) coordinates. The interpolation is designed so that the model coordinates can be irregularly spaced (i.e., not on a grid).

### netCDFsub\_mod.f90

Module that defines various subroutines for reading radar reflectivity observations in netCDF format. Used by `module_read_NSSl_mosaic.f90`.

### process\_NSSL\_mosaic.f90

This is the main driver for this program. The steps followed by the driver include:

1. Read in radar reflectivity observations and MPAS mesh information.
2. Interpolate reflectivity observations to the model using `mosaic2grid` from `mosaic_interp_mod.f90`.
3. Dump output to `RefInGSI.dat`.
4. If `output_netcdf = .true.`, perform additional processing and dump results to a netCDF file.

### process\_NSSL\_mosaic\_FV3.f90

Older version of `process_NSSL_mosaic.f90` for the FV3 dynamical core. This file is included here for reference and is not compiled as part of the program.

### read\_grib2\_mod.f90

Module that defines various subroutines for reading radar reflectivity observations in GRIB2 format. Used by `module_read_NSSl_mosaic.f90`.

### read\_ncep\_binary\_mod.f90

Module that defines various subroutines for reading radar reflectivity observations in NCEP binary format. Used by `module_read_NSSl_mosaic.f90`.

### read\_nssl\_binary\_mod.f90

Module that defines various subroutines for reading radar reflectivity observations in NSSL binary format. Used by `module_read_NSSl_mosaic.f90`.

### write\_bufr\_ref\_FV3.f90

Defines a single subroutine for writing output to a BUFR file. This subroutine has since been moved to `write_nsslref_mod.f90`. This file is included here for reference and is not compiled as part of the program.

### write\_nsslref\_mod.f90

Module that defines various subroutines for writing output to either a BUFR ot netCDF file. Writing output to a BUFR file is manually disabled in `process_NSSL_mosaic.f90`.
