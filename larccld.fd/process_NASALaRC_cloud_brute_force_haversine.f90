program  process_NASALaRC_cloud
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-01-18
!
! ABSTRACT: 
!     This routine read in NASA LaRC cloud products and 
!     interpolate them into ESG grid
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  
!
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90 + EXTENSIONS
!   MACHINE:  wJET
!
!$$$
!
!_____________________________________________________________________
!
  use mpi
  use kinds, only: r_kind,i_kind,r_single
  use mpasio, only: read_MPAS_nCell,read_MPAS_lat_lon

  implicit none
!
  INCLUDE 'netcdf.inc'
!
! MPI variables
  integer :: npe, mype,ierror
!SATID
!  integer, parameter :: satidgoeswest=259  ! GOES 15  Stopped after March 2nd, 2020
!  integer, parameter :: satidgoeswest=271  ! GOES 17 stopped after January 4th,2022
  integer, parameter :: satidgoeswest=272  ! GOES 17 stopped after January 4th,2022
  integer, parameter :: satidgoeseast=270  ! GOES 16
  integer,parameter  :: boxMAX=10
!
  character*256 output_file
!
! MPAS mesh
  integer(i_kind) :: nCell
  REAL(r_single), allocatable :: lat_m(:)
  REAL(r_single), allocatable :: lon_m(:)
  CHARACTER*180   meshfile
!
!  For NASA LaRC 
!
  CHARACTER*180   workPath
  CHARACTER*80   satfile
  INTEGER ::   nxp, nyp  ! dimension
  
!     ****VARIABLES FOR THIS NETCDF FILE****
!
  CHARACTER*24 :: cbase_time
  INTEGER(i_kind) ::  base_time
  INTEGER(i_kind) ::  ibase_year,ibase_month,ibase_day,ibase_hour,ihour
  INTEGER(i_kind) ::  icycle_year,icycle_month,icycle_day,icycle_hour
  REAL*8      time_offset
  REAL(r_single), allocatable ::   lat_l(:)
  REAL(r_single), allocatable ::   lon_l(:)
  REAL(r_single), allocatable ::   lwp_l(:)
  REAL(r_single), allocatable ::   teff_l(:)
  REAL(r_single), allocatable ::   ptop_l(:)
  integer(i_kind), allocatable ::   phase_l(:)
!
!  array for RR
!
  REAL(r_single), allocatable ::   w_pcld(:,:)
  REAL(r_single), allocatable ::   w_tcld(:,:)
  REAL(r_single), allocatable ::   w_frac(:,:)
  REAL(r_single), allocatable ::   w_lwp (:,:)
  integer(i_kind),allocatable ::   nlev_cld(:,:)
!
! Working
!
  integer  nfov
  parameter (nfov=160)
  real, allocatable ::     Pxx(:,:),Txx(:,:),WPxx(:,:)
  real,allocatable  ::     xdist(:,:), xxxdist(:)
  real     fr,sqrt, qc, type
  integer,allocatable  ::  PHxx(:,:),index(:), jndex(:)
  integer  ixx,ii,jj,med_pt,igrid,jgrid  &
               ,ncount,ncount1,ncount2,ii1,jj1,nobs,n

!
! namelist
!
  integer            :: analysis_time
  integer            :: ioption
  character(len=100) :: bufrfile
  integer(i_kind)    :: rad_km, nptsx, nptsy
  integer(i_kind)    :: boxhalfx(boxMAX), boxhalfy(boxMAX)
  real (r_kind)      :: boxlat0(boxMAX)
  character(len=20)  :: grid_type
  real (r_kind)      :: userDX
  namelist/setup/ grid_type,analysis_time, ioption, rad_km,bufrfile, &
                  boxhalfx, boxhalfy, boxlat0,userDX
!
!
!  ** misc
      
!  real(edp)     :: rlat,rlon
!  real(edp)     :: xc,yc

  integer i,j,k,ipt,jpt,cfov,ibox,imesh
  Integer nf_status,nf_fid,nf_vid

  integer :: NCID
  integer(8) :: east_time, west_time
  integer :: isat, maxobs,numobs

  integer :: status
  character*10  atime
  logical :: ifexist

! For testing haversine subroutine
  real, allocatable :: dist(:)


!**********************************************************************
!
!            END OF DECLARATIONS....start of program
!
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  if(mype==0) then
!
!  get namelist
!
     analysis_time=2018051718
     bufrfile='NASALaRCCloudInGSI_bufr.bufr'
     rad_km=3
     boxhalfx=-1
     boxhalfy=-1
     boxlat0= 999.0 !don't use variable box by default
      ! * ioption = 1 is nearest neighrhood
      ! * ioption = 2 is median of cloudy fov
     ioption = 2
     grid_type="none"
     userDX=3000.0
 
     inquire(file='namelist.nasalarc', EXIST=ifexist )
     if(ifexist) then
       open(10,file='namelist.nasalarc',status='old')
          read(10,setup)
       close(10)
       write(*,*) 'Namelist setup are:'
       write(*,setup)
     else
       write(*,*) 'No namelist file exist, use default values'
       write(*,*) "analysis_time,bufrfile,rad_km,ioption"
       write(*,*) analysis_time, trim(bufrfile),rad_km,ioption
       write(*,*) "boxhalfx,boxhalfy,boxlat0"
       write(*,*) boxhalfx,boxhalfy,boxlat0
     endif

!
! read in MPAS mesh information
!

     meshfile='mesh.nc'
     call read_MPAS_nCell(meshfile, nCell)
     write(6,*) 'model nCell   =', nCell
     allocate(lat_m(nCell))
     allocate(lon_m(nCell))
     call read_MPAS_lat_lon(meshfile, nCell, lat_m, lon_m)
     write(6,*) 'max model lat =', maxval(lat_m)
     write(6,*) 'min model lat =', minval(lat_m)
     write(6,*) 'max model lon =', maxval(lon_m)
     write(6,*) 'min model lon =', minval(lon_m)

!
!  read in the NASA LaRC cloud data
!
     satfile='lgycld.bufr_d'
     call read_NASALaRC_cloud_bufr_survey(satfile,satidgoeseast,satidgoeswest,east_time, west_time,maxobs)
     if(maxobs==0) then
        write(*,*) "WARNING: no observation available"
        stop 0
     endif
     allocate(lat_l(maxobs))
     allocate(lon_l(maxobs))
     allocate(ptop_l(maxobs))
     allocate(teff_l(maxobs))
     allocate(phase_l(maxobs))
     allocate(lwp_l(maxobs))
     ptop_l=-9.
     teff_l=-9.
     lat_l =-9.
     lon_l =-9.
     lwp_l =-9.
     phase_l=-9
     call read_NASALaRC_cloud_bufr(satfile,atime,satidgoeseast,satidgoeswest,east_time, west_time,   &   
            maxobs,numobs, ptop_l, teff_l, phase_l, lwp_l,lat_l, lon_l)

     write(6,*)
     write(6,*) 'number of obs = ', numobs
     write(6,'(6a12)')'ptop', 'teff', 'lat', 'lon', 'lwp', 'phase'
     do j=1,numobs,numobs/50
        write(6,'(4f12.3,f12.4,I12)') ptop_l(j),teff_l(j),lat_l(j),lon_l(j),lwp_l(j),phase_l(j)
     enddo

     do i=1,numobs
        if (phase_l(i).eq.4) phase_l(i) = 0   ! clear
        if (phase_l(i).eq.5) phase_l(i) = -9  ! bad data
        if (phase_l(i).eq.0) ptop_l(i) = -20.
     enddo

! -----------------------------------------------------------
! -----------------------------------------------------------
!     Map each FOV onto model mesh
! -----------------------------------------------------------
! -----------------------------------------------------------
!
     allocate (Pxx(nCell,nfov),Txx(nCell,nfov),WPxx(nCell,nfov))
     allocate (xdist(nCell,nfov), xxxdist(nfov))
     allocate (PHxx(nCell,nfov),index(nCell), jndex(nfov))
     allocate (dist(nCell))
     index=0

     do ipt=1,numobs
       if ( mod(ipt, 10).eq.0 ) then
         write(6,*) 'In ipt loop, ipt =', ipt
       endif
       if (phase_l(ipt).ge.0) then
!  Indicates there is some data (not missing)

! Compute distances between ob and all mesh cells
         call compute_haversine_dist(nCell, lat_l(ipt), lon_l(ipt), lat_m, lon_m, dist)

         if ( minval(dist).le.rad_km ) then
! We have at least one cell within rad_km of the ob
           do imesh=1,nCell
             if ( dist(imesh).le.rad_km ) then
               if (index(imesh).lt.nfov) then
                 index(imesh) = index(imesh) + 1
                 Pxx(imesh,index(imesh))   = Ptop_l(ipt)
                 Txx(imesh,index(imesh))   = Teff_l(ipt)
                 PHxx(imesh,index(imesh))  = phase_l(ipt)
                 WPxx(imesh,index(imesh))  = lwp_l(ipt)
                 xdist(imesh,index(imesh)) = dist(nCell)
               else
                 write(6,*) 'ALERT: too many data in one grid, increase nfov'
                 write(6,*) nfov, imesh
               endif
             endif ! mesh cell is within rad_km of ob
           enddo ! imesh
         endif ! at least one mesh cell within rad_km

       endif   ! phase_l >= 0
     enddo   ! ipt

     deallocate(lat_l,lon_l,ptop_l,teff_l,phase_l,lwp_l)
     write(6,*) 'The max index number is: ', maxval(index)

!=========================================================================
! There should be additional code here that has not been added yet. See 
! process_NASALaRC_cloud.f90
!
! Development on this approach stopped b/c it is WAY too slow
!
!=========================================================================

     write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="
  endif ! if mype==0 

  call MPI_FINALIZE(ierror)

end program process_NASALaRC_cloud

subroutine read_NASALaRC_cloud_bufr(satfile,atime,satidgoeseast,satidgoeswest,east_time, west_time, &
             maxobs,numobs,ptop, teff, phase, lwp_iwp,lat, lon)
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2010-07-09
!
! ABSTRACT: 
!     This routine read in NASA LaRC cloud products 
!     from a bufr file                      
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  
!
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90 + EXTENSIONS
!   MACHINE:  wJET
!
!$$$
!
!_____________________________________________________________________
!
  use kinds, only: r_kind,i_kind

  implicit none
!
!
!
  character(80):: hdstr='YEAR  MNTH  DAYS HOUR  MINU  SECO SAID'
  character(80):: obstr='CLATH  CLONH CLDP HOCT CDTP EBBTH VILWC'
! CLDP     |  CLOUD PHASE
! HOCB     | HEIGHT OF BASE OF CLOUD
! HOCT     | HEIGHT OF TOP OF CLOUD             (METERS)
! CDBP     | PRESSURE AT BASE OF CLOUD
! CDTP     | PRESSURE AT TOP OF CLOUD           (PA)
! EBBTH    | EQUIVALENT BLACK BODY TEMPERATURE  (KELVIN)
! VILWC    | VERTICALLY-INTEGRATED LIQUID WATER CONTENT

  real(8) :: hdr(7),obs(7,1)

  INTEGER :: ireadmg,ireadsb

  character(8) subset
  integer :: unit_in=10,idate,iret,nmsg,ntb
  integer :: satid

!
!  For NASA LaRC 
!
  CHARACTER*40,intent(in) ::   satfile
!SATID
  integer,intent(in) :: satidgoeswest
  integer,intent(in) :: satidgoeseast  
  integer(8),intent(in) :: east_time, west_time

  INTEGER,intent(in)  ::   maxobs! dimension
  INTEGER,intent(out) ::   numobs  ! dimension
  INTEGER(8) ::  obs_time
  REAL*8      time_offset
  REAL*4      lat                            (  maxobs)
  REAL*4      lon                            (  maxobs)
  integer     phase                          (  maxobs)
  REAL*4      lwp_iwp                        (  maxobs)
  REAL*4      teff                           (  maxobs)
  REAL*4      ptop                           (  maxobs)
!
!
!  ** misc

!  integer i,j,k
!  Integer nf_status,nf_fid,nx,ny,nf_vid
!
!  integer :: status
  character*10  atime
!
!**********************************************************************
!
 open(24,file='NASA.bufrtable')
 open(unit_in,file=trim(satfile),form='unformatted')
 call openbf(unit_in,'IN',unit_in)
 call dxdump(unit_in,24)
 call datelen(10)
   nmsg=0
   ntb = 0
   msg_report: do while (ireadmg(unit_in,subset,idate) == 0)
     nmsg=nmsg+1
     sb_report: do while (ireadsb(unit_in) == 0)
       call ufbint(unit_in,hdr,7,1,iret,hdstr)
       obs_time=(hdr(1)-2000.0_8)*100000000.0_8+hdr(2)*1000000.0_8+hdr(3)*10000.0_8+hdr(4)*100.0_8+hdr(5)
       satid=int(hdr(7))
       if( (obs_time == east_time .and. satid==satidgoeseast ) .or.  &
           (obs_time == west_time .and. (satid==satidgoeswest .or. satid==259 .or. satid==271) ) ) then
         call ufbint(unit_in,obs,7,1,iret,obstr)
         if(abs(obs(3,1) -4.0) < 1.e-4) then
           obs(7,1)=99999. ! clear
           obs(6,1)=99999. ! clear
           obs(5,1)=101300.0  ! clear (hpa)
         endif
         if(obs(5,1) < 1.e7 .and. obs(5,1) > 100.0 ) then
         if(obs(6,1) < 1.e7 .and. obs(6,1) > 10.0) then
           ntb = ntb+1
           if(ntb > maxobs) then
             write(*,*) 'ALERT: need to increase maxobs',maxobs, ntb
             ntb = maxobs
           endif

           lwp_iwp(ntb)=99999.0
           lat(ntb)=99999.0
           lon(ntb)=99999.0
           phase(ntb)=99999
           teff(ntb)=99999.0
           ptop(ntb)=99999.0
           if(obs(1,1) < 1.e9) lat(ntb)=real(obs(1,1))
           if(obs(2,1) < 1.e9) lon(ntb)=real(obs(2,1))
           if(obs(3,1) < 1.e9) phase(ntb)=int(obs(3,1))
           if(obs(7,1) < 1.e9) lwp_iwp(ntb)=real(obs(7,1))
           if(obs(6,1) < 1.e9) teff(ntb)=real(obs(6,1))
           if(obs(5,1) < 1.e9) ptop(ntb)=real(obs(5,1))/100.0 ! pa to hpa

         endif
         endif
       endif   ! east_time, west_time
     enddo sb_report
   enddo msg_report
   write(*,*) 'message/reports num=',nmsg,ntb
 call closbf(unit_in)
 numobs=ntb
 write(atime,'(I10)') idate

end subroutine read_NASALaRC_cloud_bufr

subroutine sortmed(p,n,is)
      real p(n)
      integer is(n)
! * count cloudy fov
      real    f
      integer cfov
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

subroutine read_NASALaRC_cloud_bufr_survey(satfile,satidgoeseast,satidgoeswest,east_time, west_time,maxobs)
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2010-07-09
!
! ABSTRACT: 
!     This routine read in NASA LaRC cloud products 
!     from a bufr file                      
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  
!
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90 + EXTENSIONS
!   MACHINE:  wJET
!
!$$$
!
!_____________________________________________________________________
!
  use kinds, only: r_kind,i_kind

  implicit none
!
  character(80):: hdstr='YEAR  MNTH  DAYS HOUR  MINU  SECO  SAID'
  real(8) :: hdr(7)

  INTEGER :: ireadmg,ireadsb

  character(8) subset
  integer :: unit_in=10,idate,iret,nmsg,ntb

!
!  For NASA LaRC 
!
  CHARACTER*40, intent(in)    ::   satfile
!SATID
  integer,intent(in) :: satidgoeswest
  integer,intent(in) :: satidgoeseast  
  integer(8),intent(out) :: east_time, west_time
  integer,intent(out) :: maxobs

  INTEGER ::   numobs  ! dimension
  INTEGER(8) ::  obs_time
  REAL*8      time_offset

  INTEGER(i_kind),parameter :: max_obstime=30
  integer(8) :: num_obstime_all(max_obstime)
  integer(i_kind) :: num_subset_all(max_obstime)
  integer(i_kind) :: num_obstime_hh(max_obstime)
  integer(i_kind) :: num_satid(max_obstime)
  integer(i_kind) :: num_obstime

!
  character*10  atime
  integer :: i,ii,hhh
  integer :: numobs_east, numobs_west
  integer :: satid
!
!**********************************************************************
!
 num_obstime=0
 num_satid=0
 num_obstime_all=0
 num_subset_all=0
 hhh=99
 open(24,file='NASA.bufrtable')
 open(unit_in,file=trim(satfile),form='unformatted')
 call openbf(unit_in,'IN',unit_in)
 call dxdump(unit_in,24)
 call datelen(10)
   nmsg=0
   msg_report: do while (ireadmg(unit_in,subset,idate) == 0)
     ntb = 0
     nmsg=nmsg+1
     sb_report: do while (ireadsb(unit_in) == 0)
       call ufbint(unit_in,hdr,7,1,iret,hdstr)
       obs_time=(hdr(1)-2000.0_8)*100000000.0_8+hdr(2)*1000000.0_8+hdr(3)*10000.0_8+hdr(4)*100.0_8+hdr(5)
       hhh=int(hdr(5))
       ntb=ntb+1
       satid=int(hdr(7))
     enddo sb_report
! message inventory
     if(num_obstime == 0 ) then
       num_obstime=1
       num_obstime_all(num_obstime)=obs_time
       num_obstime_hh(num_obstime)=hhh
       num_subset_all(num_obstime)= ntb
       num_satid(num_obstime)=satid
     else
       ii=0
       DO i=1,num_obstime
          if(num_obstime_all(i) == obs_time .and. num_satid(i)== satid) ii=i
       ENDDO
       if( ii > 0 .and. ii <=num_obstime) then
          num_subset_all(ii)=num_subset_all(ii) + ntb
       else
          num_obstime=num_obstime+1
          if(num_obstime> max_obstime) then
             write(*,*) 'Error: too many message types'
             write(*,*) 'Need to increase :max_obstime'
             stop 1234
          endif
          num_obstime_all(num_obstime)=obs_time
          num_obstime_hh(num_obstime)=hhh
          num_subset_all(num_obstime)=num_subset_all(num_obstime)+ntb
          num_satid(num_obstime)=satid
       endif
     endif
   enddo msg_report
   write(*,*) 'message/reports num=',nmsg,ntb
 call closbf(unit_in)

 write(*,'(2x,a15,a15,a15,a15)') 'time_level','satid','subset_num','hour'
 DO i=1,num_obstime
   write(*,'(i2,i15,i15,i15,i15)') i,num_obstime_all(i),num_satid(i),num_subset_all(i),num_obstime_hh(i)
 ENDDO
!  GOES EAST  : 1815, 1845, 1915, 2045, no anymore, changed to 1830, 1900 just like WEST
!  GOES WEST  : 1830, 1900, 2030
 east_time=0
 west_time=0
 numobs_east=0
 numobs_west=0
 DO i=1,num_obstime
   if(num_subset_all(i) > 10) then
      if(num_satid(i) == satidgoeseast ) then  
         if(east_time < num_obstime_all(i)) then
              east_time=num_obstime_all(i)
              numobs_east=num_subset_all(i)
         endif
      endif
      if(num_satid(i) == satidgoeswest .or. num_satid(i)==259 .or. num_satid(i)==271 ) then
         if(west_time < num_obstime_all(i)) then
             west_time=num_obstime_all(i)
             numobs_west=num_subset_all(i)
         endif
      endif
   endif
 ENDDO
 write(*,*) 'east_time and number=',east_time,numobs_east
 write(*,*) 'west_time and number=',west_time,numobs_west
 
 maxobs=numobs_west+numobs_east
 maxobs=maxobs+int(maxobs*0.2)
 write(*,*) 'maxobs=',maxobs

end subroutine read_NASALaRC_cloud_bufr_survey

subroutine compute_haversine_dist(n, latpt, lonpt, lat, lon, dist)
!
! Compute the distance between two (lat, lon) coordinates in m using the haversine distance formula
!
! https://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise.haversine_distances.html
!
! Cross-check using this web tool: https://www.nhc.noaa.gov/gccalc.shtml

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: latpt, lonpt
  real, intent(in) :: lat(n), lon(n)  
  real, intent(out) :: dist(n)

  real :: latptr, lonptr
  real :: latr(n), lonr(n), dlat(n), dlon(n), asin_arg(n)
  real :: deg2rad, rearth

  ! Define constants
  deg2rad = acos(-1.) / 180.
  rearth = 6371200.

  ! Convert inputs to radians
  latptr = latpt * deg2rad
  lonptr = lonpt * deg2rad
  latr = lat * deg2rad
  lonr = lon * deg2rad
  dlat = latr - latptr
  dlon = lonr - lonptr

  asin_arg = sin(dlat/2.)**2 + cos(latptr) * cos(latr) * sin(dlon/2.)**2
  dist = rearth * 2 * asin(sqrt(asin_arg))

end subroutine compute_haversine_dist
