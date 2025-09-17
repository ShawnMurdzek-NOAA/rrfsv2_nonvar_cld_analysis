# Nonvariational Cloud Analysis Utilities

This directory contains a variety of modules used by the five main nonvariational cloud analysis programs. These modules include:

- `constants.f90`: Defines various physical constants.
- `kinds.f90`: Defines various variable types.
- `map_proj_helper_mod.f90`: Helper subroutines for map projections used in larccld.fd, lightning.fd, and metarcld.fd.
- `mpasio.f90`: Defines various subroutines for reading and writing MPAS netCDF output.
- `wps_map_utils`: Define various subroutines related to map projections. This code comes from the WRF Preprocessing System (WPS) and can be found in that program at WPS/geogrid/src

The modules here must be compiled prior to compiling any of the five main nonvariational cloud analysis programs.
