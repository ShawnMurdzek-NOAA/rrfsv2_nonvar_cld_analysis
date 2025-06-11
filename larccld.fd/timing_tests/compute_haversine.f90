program compute_havsersine

  implicit none

  integer :: i
  real :: lat1, lon1, lat2, lon2
  real :: lat1r, lon1r, lat2r, lon2r
  real :: dist, dlat, dlon, asin_arg
  real :: deg2rad, rearth

  ! Define constants
  deg2rad = acos(-1.) / 180.
  rearth = 6371200.

  ! Define (lat, lon) values
  lat1 = 29.9969
  lon1 = -110.9375
  lat2 = 71.78191
  lon2 = -76.21359

  do i=1,1000000000
    ! Convert inputs to radians
    lat1r = lat1 * deg2rad
    lon1r = lon1 * deg2rad
    lat2r = lat2 * deg2rad
    lon2r = lon2 * deg2rad
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r

    asin_arg = sin(dlat/2.)**2 + cos(lat1r) * cos(lat2r) * sin(dlon/2.)**2
    dist = rearth * 2 * asin(sqrt(asin_arg))
  enddo

  print*, dist

end program
