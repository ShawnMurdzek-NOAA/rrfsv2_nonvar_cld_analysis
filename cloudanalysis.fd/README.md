# Cloud Analysis Program

This program performs the nonvariational cloud analysis on the MPAS mesh using preprocessed observations from larccld.fd, lightning.fd, metarcld.fd, and mosaic.fd.

## Compiling

1. Load the proper environment from `../env/`
2. `make clean && make`

## Running

To run cloudanalysis.fd by itself, one can use the sample slurm submission script provided in this directory (`run_mpas_nonvarcldana.sh`). The program is designed to run on an arbitrary number of processors. The parallelization strategy is as follows:

1. Read in all MPAS fields using the root processor, then equally distribute to all processors. Each field is distributed right after it is read so that the root processor does not need to hold all full 3D model fields in memory at the same time.
2. Each processor performs the cloud analysis using a subset of all MPAS columns. The cloud analysis operates on vertical columns, so the columns distributed to a particular processor do not need to be adjacent.
3. Send output fields back to the root processor to be written to the output netCDF file. Similar to step (1), each output field is written right after it is sent to the root processor in order to reduce the memory load on the root processor.

### Inputs

- `NASALaRC_cloud4mpas.bin`: NASA LaRC observations interpolated to the MPAS mesh by larccld.fd.
- `LightningInMPAS.dat`: Lightning observations interpolated to the MPAS mesh by lightning.fd.
- `mpas_metarcloud.bin`: Ceilometer observations interpolated to the MPAS mesh by metarcld.fd.
- `RefInGSI3D.dat`: Radar reflectivity observations interpolated to the MPAS mesh by mosaic.fd.
- `invariant.nc`: An MPAS invariant file that includes the (lat, lon) coordinates, land mask, terrain height, and height MSL of each MPAS cell.
- `mpasout.nc`: An MPAS background netCDF file that contains the following fields (variable name in MPAS given in parenthesis, fields that are updated have an asterisk):
  - Potential temperature (`theta`)*
  - Surface pressure (`surface_pressure`)
  - Pressure (`pressure_p` and `pressure_base`)
  - Water vapor mass mixing ratio (`qv`)*
  - Cloud water mass mixing ratio (`qc`)*
  - Cloud ice mass mixing ratio (`qi`)*
  - Rain mass mixing ratio (`qr`)*
  - Snow mass mixing ratio (`qs`)*
  - Graupel mass mixing ratio (`qg`)*
  - Rain number concentration (`nr`)*
  - Cloud ice number concentration (`ni`)*
  - Cloud water number concentration (`nc`)*
  - Cloud fraction (`cldfrac`)
  - Skin temperature (`skintemp`)
- `gsiparm.anl`: Namelist file (see options below).

### Outputs

- `mpasout.nc`: Same file as the `mpasout.nc` input file, but with the fields listed with an asterisk above updated.

### Namelist Options

Namelist options are included in either `&SETUP` or `&RAPIDREFRESH_CLDSURF`

#### SETUP

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| `iyear` | 2020 | Year of analysis time. |
| `imonth` | 07 | Month of analysis time. |
| `iday` | 29 | Day of analysis time. |
| `ihour` | 10 | Hour of analysis time. |
| `iminute` | 0 | Minute of analysis time. |
| `isecond` | 0 | Second of analysis time. |
| `dump_cld_cover_3d` | 0 | Add cld\_cover\_3d field to output netCDF file. Useful for debugging. |

#### RAPIDREFRESH\_CLDSURF

A full list of all parameters can be found in `namelist_mod.f90` and `rapidrefresh_cldsurf_mod.f90`. Many of these parameters are no longer in use, so they are not included here.

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| `l_conserve_thetaV` | .false. | Option to conserve virtual potential temperature during moisture adjustment. |
| `r_cleanSnow_WarmTs_threshold` | 8 | Surface temperature threshold for cleaning snow. Only keep snow if surface temperature is less than this threshold. |
| `i_conserve_thetaV_iternum` | 3 | Iteration number for conserving virtual potential temperature during moisture adjustment. |
| `l_cld_bld` | .false. | Option to turn on cloud building using NASA LaRC observations.  |
| `cld_bld_hgt` | 1200 | Only build clouds below this height (m). |
| `build_cloud_frac_p` | 0.95 | NASA LaRC cloud fraction threshold for building clouds. |
| `clear_cloud_frac_p` | 0.1 | NASA LaRC cloud fraction threshold for clearing clouds. |
| `iclean_hydro_withRef` | 1 | Option to clean hydrometeors if a grid point has no ehco and maximum reflectivity = 0 dBZ. |
| `iclean_hydro_withRef_allcol` | 0 | Option to clean the whole column of hydrometeors if the maximum reflectivity = 0 dBZ and satellite observations show no clouds. |
| `l_use_hydroretrieval_all` | .false. | Alternative option for the hydrometeor analysis that is typically used for 3DRTMA. Overrides `l_precip_clear_only`. |
| `i_lightpcp` | 0 | Option to add light precipitation in the warm section. |
| `l_numconc` | .false. | Option to update cloud water and cloud ice number concentration. |
| `qv_max_inc` | 0.005 | Maximum allowed water vapor increment (kg/kg). |
| `l_precip_clear_only` | .false. | Option to only clear precipitating hydrometeors.  |
| `l_fog_off` | .false. | Option to turn off the impact of METAR visibility on model cloud fields. |
| `cld_bld_coverage` | 0.6 | Minimum cloud coverage required for building cloud ice and cloud water. |
| `cld_clr_coverage` | 0.6 | Maximum cloud coverage required for clearing cloud ice and cloud water. |
| `i_T_Q_adjust` | 1 | Option for adjusting temperature and humidity. 0 = No adjustment. 1 = Adjustment. 2 = Only adjust for cloud clearing. |
| `i_precip_vertical_check` | 0 | Various options for adjusting rain, snow, and graupel fields after analysis. See `rapidrefresh_cldsurf_mod.f90` for details. |
| `l_rtma3d` | .false. | Option to turn on the configuration for 3D RTMA. |
| `l_qnr_from_qr` | .false. | Option to compute rain number concentration from rain mixing ratio. |
| `n0_rain` | 1e8 | Intercept parameter for raindrop size distribution. |
| `r_cloudfrac_threshold` | 0.45 | Model cloud fraction threshold for cloud building. |

## Code Overview

Refer to the `Makefile` for the files actually compiled as part of this program. The main driver is `cloudanalysis_mpas_driver.f90` and `get_mpas_bk_mod.f90` contains all the subroutines for reading and writing to the MPAS netCDF files. The `NonVarCldLib/` directory contains various subroutines used the cloud analysis that are model agnostic.
