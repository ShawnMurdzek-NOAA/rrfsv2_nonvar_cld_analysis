program compute_euclidean_dist

  implicit none

  integer :: i
  real :: x1, y1, x2, y2, dist

  x1 = 34.2
  y1 = 12.2
  x2 = 10.
  y2 = 17.4

  do i=1,1000000000
    dist = sqrt((x2 - x1)**2 + (y2 - y1)**2)
  enddo

  print*, dist

end program
