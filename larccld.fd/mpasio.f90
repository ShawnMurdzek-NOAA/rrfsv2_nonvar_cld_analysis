module mpasio

contains

subroutine read_MPAS_nCell(mpasfile, nCell)
!---------------------------------------------------------------------------------------------------
!
! Determine the number of MPAS cells ina file
!
! Inputs
!   mpasfile : string
!     MPAS netCDF file name
!
! Outputs
!   nCell : integer
!     Number of cells
! 
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

  use netcdf

  implicit none

  character(len=100), intent(in) :: mpasfile
  integer, intent(out) :: nCell

  integer :: ncid, stat, cellid
  character(len=50) :: cellname

  stat = nf90_open(mpasfile, nf90_nowrite, ncid)

  stat = nf90_inq_dimid(ncid, 'nCells', cellid)
  stat = nf90_inquire_dimension(ncid, cellid, cellname, nCell)

  stat = nf90_close(ncid)

end subroutine read_MPAS_nCell

subroutine read_MPAS_lat_lon(mpasfile, nCell, lat, lon)
!---------------------------------------------------------------------------------------------------
!
! Read in (lat, lon) coordinates for MPAS mesh grid cells
!
! Inputs
!   mpasfile : string
!     MPAS netCDF file name
!   nCell : integer
!     Number of mesh cells
!
! Outputs
!   lat : real
!     MPAS cell latitudes (deg N)
!   lon : real
!     MPAS cell longitudes (deg E)
! 
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

  use netcdf

  implicit none

  character(len=100), intent(in) :: mpasfile
  integer, intent(in) :: nCell
  real, intent(out) :: lat(nCell), lon(nCell)

  integer :: ncid, stat, latid, lonid
  real, parameter :: pi=3.1415927

  stat = nf90_open(mpasfile, nf90_nowrite, ncid)

  stat = nf90_inq_varid(ncid, 'latCell', latid)
  stat = nf90_inq_varid(ncid, 'lonCell', lonid)

  stat = nf90_get_var(ncid, latid, lat)
  stat = nf90_get_var(ncid, lonid, lon)

  lat = lat * 180 / pi
  lon = lon * 180 / pi

  stat = nf90_close(ncid)

end subroutine read_MPAS_lat_lon

end module mpasio
