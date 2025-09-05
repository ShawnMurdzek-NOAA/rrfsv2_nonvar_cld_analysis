module mpasio

use netcdf

contains

subroutine handle_nc_error(stat, name)
!---------------------------------------------------------------------------------------------------
!
! Handle errors from netCDF-fortran library by taking the error code and printing a human-readable
! error message
!
! Inputs
!   stat : integer
!     Error code from netCDF-fortran function
!   name : string 
!     Name of subroutine where error occurred
!
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, intent(in) :: stat
  character(len=50), intent(in) :: name

  if (stat /= nf90_noerr) then
    write(*,*) 'ERROR in ', trim(name)
    write(*,*) trim(nf90_strerror(stat))
    stop stat
  endif

end subroutine handle_nc_error

subroutine read_MPAS_dim(mpasfile, namedim, ndim)
!---------------------------------------------------------------------------------------------------
!
! Determine the length of a dimension from an MPAS netCDF file
!
! Inputs
!   mpasfile : string
!     MPAS netCDF file name
!   namedim : string
!     Dimension name
!
! Outputs
!   ndim : integer
!     Dimension length
! 
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=100), intent(in) :: mpasfile
  character(len=50), intent(in) :: namedim
  integer, intent(out) :: ndim

  integer :: ncid, stat, dimid
  character(len=50) :: tmpname,routine_name

  routine_name = 'read_MPAS_dim'

  stat = nf90_open(mpasfile, nf90_nowrite, ncid)
  call handle_nc_error(stat, routine_name)

  stat = nf90_inq_dimid(ncid, trim(namedim), dimid)
  call handle_nc_error(stat, routine_name)
  stat = nf90_inquire_dimension(ncid, dimid, tmpname, ndim)
  call handle_nc_error(stat, routine_name)

  stat = nf90_close(ncid)
  call handle_nc_error(stat, routine_name)

end subroutine read_MPAS_dim

subroutine read_MPAS_lat_lon(mpasfile, nCell, lat, lon)
!---------------------------------------------------------------------------------------------------
!
! Read in (lat, lon) coordinates for MPAS mesh grid cells from a grid.nc file
!
! Inputs
!   mpasfile : string
!     MPAS netCDF file name
!   nCell : integer
!     Number of mesh cells
!
! Outputs
!   lat : real
!     MPAS cell latitudes (deg N, range -90 to 90)
!   lon : real
!     MPAS cell longitudes (deg E, range -180 to 180)
! 
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=100), intent(in) :: mpasfile
  integer, intent(in) :: nCell
  real, intent(out) :: lat(nCell), lon(nCell)

  integer :: ncid, stat, latid, lonid, i
  real, parameter :: pi=3.1415927
  character(len=50)  :: routine_name

  routine_name = 'read_MPAS_lat_lon'

  stat = nf90_open(mpasfile, nf90_nowrite, ncid)
  call handle_nc_error(stat, routine_name)

  stat = nf90_inq_varid(ncid, 'latCell', latid)
  call handle_nc_error(stat, routine_name)
  stat = nf90_inq_varid(ncid, 'lonCell', lonid)
  call handle_nc_error(stat, routine_name)

  stat = nf90_get_var(ncid, latid, lat)
  call handle_nc_error(stat, routine_name)
  stat = nf90_get_var(ncid, lonid, lon)
  call handle_nc_error(stat, routine_name)

  lat = lat * 180 / pi
  lon = lon * 180 / pi

  ! Ensure that longitudes are in the range (-180, 180)
  do i=1,nCell
    if (lon(i) > 180) then
      lon(i) = lon(i) - 360.
    endif
  enddo

  stat = nf90_close(ncid)
  call handle_nc_error(stat, routine_name)

end subroutine read_MPAS_lat_lon

subroutine read_MPAS_1D_int(mpasfile, ndim, field, data)
!---------------------------------------------------------------------------------------------------
!
! Read in 1D integer field from an MPAS netCDF file
!
! Inputs
!   mpasfile : string
!     MPAS netCDF file name
!   ndim : integer
!     Length of the one dimension
!   field : string
!     1D integer field to read
!
! Outputs
!   data : integer
!     Integer field
! 
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=100), intent(in) :: mpasfile
  integer, intent(in) :: ndim
  character(len=50), intent(in) :: field
  integer, intent(out) :: data(ndim)

  integer :: ncid, stat, fieldid
  character(len=50) :: routine_name

  routine_name = 'read_MPAS_1D_int'

  stat = nf90_open(mpasfile, nf90_nowrite, ncid)
  call handle_nc_error(stat, routine_name)
  stat = nf90_inq_varid(ncid, field, fieldid)
  call handle_nc_error(stat, routine_name)
  stat = nf90_get_var(ncid, fieldid, data)
  call handle_nc_error(stat, routine_name)
  stat = nf90_close(ncid)
  call handle_nc_error(stat, routine_name)

end subroutine read_MPAS_1D_int

subroutine read_MPAS_1D_real(mpasfile, ndim, field, data, rem_dim)
!---------------------------------------------------------------------------------------------------
!
! Read in 1D real field from an MPAS netCDF file
!
! Inputs
!   mpasfile : string
!     MPAS netCDF file name
!   ndim : integer
!     Length of the one dimension
!   field : string
!     1D real field to read
!   rem_dim : boolean, optional
!     Option to explicitly read only the first element in the last dimension (usually time)
!
! Outputs
!   data : real
!     Real field
! 
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=100), intent(in) :: mpasfile
  integer, intent(in) :: ndim
  character(len=50), intent(in) :: field
  logical, intent(in), optional :: rem_dim
  real, intent(out) :: data(ndim)

  integer :: ncid, stat, fieldid
  logical :: remove_dim
  character(len=50) :: routine_name

  remove_dim=.false.
  if (present(rem_dim)) remove_dim=rem_dim

  routine_name = 'read_MPAS_1D_real'

  stat = nf90_open(mpasfile, nf90_nowrite, ncid)
  call handle_nc_error(stat, routine_name)
  stat = nf90_inq_varid(ncid, field, fieldid)
  call handle_nc_error(stat, routine_name)
  if ( remove_dim ) then
    stat = nf90_get_var(ncid, fieldid, data, start=(/ 1, 1 /), count=(/ ndim, 1 /))
  else
    stat = nf90_get_var(ncid, fieldid, data)
  endif
  call handle_nc_error(stat, routine_name)
  stat = nf90_close(ncid)
  call handle_nc_error(stat, routine_name)

end subroutine read_MPAS_1D_real

subroutine read_MPAS_2D_real(mpasfile, ndim1, ndim2, field, data, rem_dim)
!---------------------------------------------------------------------------------------------------
!
! Read in 2D real field from an MPAS netCDF file
!
! Inputs
!   mpasfile : string
!     MPAS netCDF file name
!   ndim1 : integer
!     Length of the first dimension for field
!   ndim2 : integer
!     Length of the second dimension for field
!   field : string
!     2D real field to read
!   rem_dim : boolean, optional
!     Option to explicitly read only the first element in the last dimension (usually time)
!
! Outputs
!   data : real
!     Real field
! 
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=100), intent(in) :: mpasfile
  integer, intent(in) :: ndim1,ndim2
  character(len=50), intent(in) :: field
  logical, intent(in), optional :: rem_dim
  real, intent(out) :: data(ndim1,ndim2)

  integer :: ncid, stat, fieldid
  logical :: remove_dim
  character(len=50) :: routine_name

  remove_dim=.false.
  if (present(rem_dim)) remove_dim=rem_dim

  routine_name = 'read_MPAS_2D_real'

  stat = nf90_open(mpasfile, nf90_nowrite, ncid)
  call handle_nc_error(stat, routine_name)
  stat = nf90_inq_varid(ncid, field, fieldid)
  call handle_nc_error(stat, routine_name)
  if ( remove_dim ) then
    stat = nf90_get_var(ncid, fieldid, data, start=(/ 1, 1, 1 /), count=(/ ndim1, ndim2, 1 /))
  else
    stat = nf90_get_var(ncid, fieldid, data)
  endif
  call handle_nc_error(stat, routine_name)
  stat = nf90_close(ncid)
  call handle_nc_error(stat, routine_name)

end subroutine read_MPAS_2D_real

subroutine update_MPAS_2D_real(mpasfile, ndim1, ndim2, field, data, add_dim)
!---------------------------------------------------------------------------------------------------
!
! Write a 2D real field to an MPAS netCDF file
!
! Note: The field must already exist in the netCDF file
!
! Inputs
!   mpasfile : string
!     MPAS netCDF file name
!   ndim1 : integer
!     Length of the first dimension for field
!   ndim2 : integer
!     Length of the second dimension for field
!   field : string
!     2D real field to read
!   add_dim : boolean, optional
!     Option to explicitly add a single element in the last dimension (usually time)
!
! Outputs
!   data : real
!     Real field
! 
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=100), intent(in) :: mpasfile
  integer, intent(in) :: ndim1,ndim2
  character(len=50), intent(in) :: field
  logical, intent(in), optional :: add_dim
  real, intent(in) :: data(ndim1,ndim2)

  integer :: ncid, stat, fieldid
  logical :: add_last_dim
  character(len=50) :: routine_name

  add_last_dim=.false.
  if (present(add_dim)) add_last_dim=add_dim

  routine_name = 'update_MPAS_2D_real'

  stat = nf90_open(mpasfile, nf90_write, ncid)
  call handle_nc_error(stat, routine_name)
  stat = nf90_inq_varid(ncid, field, fieldid)
  call handle_nc_error(stat, routine_name)
  if ( add_last_dim ) then
    stat = nf90_put_var(ncid, fieldid, data, start=(/ 1, 1, 1 /), count=(/ ndim1, ndim2, 1 /))
  else
    stat = nf90_put_var(ncid, fieldid, data)
  endif
  call handle_nc_error(stat, routine_name)
  stat = nf90_close(ncid)
  call handle_nc_error(stat, routine_name)

end subroutine update_MPAS_2D_real

end module mpasio
