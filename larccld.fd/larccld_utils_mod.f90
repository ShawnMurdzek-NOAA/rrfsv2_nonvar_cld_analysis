module larccld_utils

contains

subroutine sortmed(p,n,is)

      implicit none

      integer n
      real p(n)
      integer is(n)
! * count cloudy fov
      real    f
      integer cfov,temp
      integer i,j,nm1,ip1,iold

      cfov = 0
      do i=1,n
! - changed for NASA LaRC, p set = -9 for clear FOVs
         if(p(i) .gt. 0.) cfov = cfov + 1
      enddo
      f = float(cfov)/(max(1,n))
! cloud-top pressure is sorted high cld to clear
      nm1 = n-1
      do 10 i=1,nm1
      ip1 = i+1
        do 10 j=ip1,n
        if(p(i).le.p(j)) cycle
          temp = p(i)
          p(i) = p(j)
          p(j) = temp
          iold  = is(i)
          is(i) = is(j)
          is(j) = iold
   10 continue
      return
end subroutine sortmed

subroutine compute_haversine_dist(lat1, lon1, lat2, lon2, dist)
!
! Compute the distance between two (lat, lon) coordinates in m using the haversine distance formula
!
! https://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise.haversine_distances.html
!
! Cross-check using this web tool: https://www.nhc.noaa.gov/gccalc.shtml

  implicit none

  real, intent(in) :: lat1, lon1, lat2, lon2
  real, intent(out) :: dist

  real :: lat1r, lon1r, lat2r, lon2r, dlat, dlon, asin_arg
  real :: deg2rad, rearth

  ! Define constants
  deg2rad = acos(-1.) / 180.
  rearth = 6371200.

  ! Convert inputs to radians
  lat1r = lat1 * deg2rad
  lon1r = lon1 * deg2rad
  lat2r = lat2 * deg2rad
  lon2r = lon2 * deg2rad
  dlat = lat2r - lat1r
  dlon = lon2r - lon1r

  asin_arg = sin(dlat/2.)**2 + cos(lat1r) * cos(lat2r) * sin(dlon/2.)**2
  dist = rearth * 2 * asin(sqrt(asin_arg))

end subroutine compute_haversine_dist

end module larccld_utils
