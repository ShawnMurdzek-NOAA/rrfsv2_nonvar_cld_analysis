program test_map_utils

  use misc_definitions_module
  use map_utils

  implicit none

  type(proj_info) :: proj

  ! Map projection parameters
  real :: lat1, lon1, truelat1, truelat2, stdlon, dx, knowni, knownj
  !real, parameter :: r_earth=6370.

  ! Sample data
  integer :: n, i
  real, allocatable :: lat(:), lon(:), x(:), y(:)

  ! Map projection parameters (come from the HRRR)
  lat1 = 38.5
  lon1 = -97.5
  truelat1 = 38.5
  truelat2 = 38.5
  stdlon = -97.5
  dx = 3000.
  knowni = 899
  knownj = 529

  ! Initialize map_proj structure, with all variables initialized as missing
  call map_init(proj)

  ! Set desired mpa projection specifications
  call map_set(PROJ_LC, &
               proj, &
               lat1=lat1, &
               lon1=lon1, &
               truelat1=truelat1, &
               truelat2=truelat2, &
               stdlon=stdlon, &
               dx=dx, &
               knowni=knowni, &
               knownj=knownj)

  ! Test our map projection
  n = 3
  allocate (lat(n))
  allocate (lon(n))
  allocate (x(n))
  allocate (y(n))
  lat = [40.0190,   38.9072,  25.7617]
  lon = [-105.2747, -77.0369, -80.1918]
  do i=1,n
    call latlon_to_ij(proj, lat(i), lon(i), x(i), y(i))
  enddo

  ! Print results
  print*, 'lat =', lat
  print*, 'lon =', lon
  print*, 'x =', x
  print*, 'y =', y

end program 
