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

  ! Modules from WPS
  use map_utils
  use misc_definitions_module

  implicit none
!
  INCLUDE 'netcdf.inc'
!
! Map projection
  type(proj_info) :: proj
  integer :: nlat, nlon
  real :: lat1, lon1, truelat1, truelat2, stdlon, dx, knowni, knownj
  real :: ll_lat, ur_lat, left_lat, bot_lat, right_lat, top_lat
  real :: ll_lon, ur_lon, left_lon, bot_lon, right_lon, top_lon
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
  REAL(r_single), allocatable ::   w_pcld(:)
  REAL(r_single), allocatable ::   w_tcld(:)
  REAL(r_single), allocatable ::   w_frac(:)
  REAL(r_single), allocatable ::   w_lwp (:)
  integer(i_kind),allocatable ::   nlev_cld(:)
!
! Working
!
  integer  nfov
  parameter (nfov=160)
  real, allocatable ::     Pxx(:,:,:),Txx(:,:,:),WPxx(:,:,:)
  real,allocatable  ::     xdist(:,:,:), xxxdist(:)
  real, allocatable ::     latxx(:,:,:),lonxx(:,:,:)
  real     fr,sqrt, qc, type
  integer,allocatable  ::  PHxx(:,:,:),index(:,:), jndex(:)
  integer  ixx,ii,jj,med_pt,igrid,jgrid  &
               ,ncount,ncount1,ncount2,ii1,jj1,nobs,n,im

!
! namelist
!
  integer            :: analysis_time
  integer            :: ioption
  character(len=100) :: bufrfile
  integer(i_kind)    :: npts_rad, nptsx, nptsy
  integer(i_kind)    :: boxhalfx(boxMAX), boxhalfy(boxMAX)
  real (r_kind)      :: boxlat0(boxMAX)
  character(len=20)  :: grid_type
  real (r_kind)      :: userDX
  integer            :: debug
  namelist/setup/ grid_type,analysis_time, ioption, npts_rad,bufrfile, &
                  boxhalfx, boxhalfy, boxlat0,userDX,debug
!
!
!  ** misc
      
  real :: rlat,rlon
  real :: xc,yc
  real :: dist

  integer i,j,k,ipt,jpt,cfov,ibox,imesh
  Integer nf_status,nf_fid,nf_vid

  integer :: NCID
  integer(8) :: east_time, west_time
  integer :: isat, maxobs,numobs

  integer :: status
  character*10  atime
  logical :: ifexist

  integer :: noutside

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
     npts_rad=1
     boxhalfx=-1
     boxhalfy=-1
     boxlat0= 999.0 !don't use variable box by default
      ! * ioption = 1 is nearest neighrhood
      ! * ioption = 2 is median of cloudy fov
     ioption = 2
     grid_type="none"
     userDX=3000.0
     debug=0
 
     inquire(file='namelist.nasalarc', EXIST=ifexist )
     if(ifexist) then
       open(10,file='namelist.nasalarc',status='old')
          read(10,setup)
       close(10)
       write(*,*) 'Namelist setup are:'
       write(*,setup)
     else
       write(*,*) 'No namelist file exist, use default values'
       write(*,*) "analysis_time,bufrfile,npts_rad,ioption"
       write(*,*) analysis_time, trim(bufrfile),npts_rad,ioption
       write(*,*) "boxhalfx,boxhalfy,boxlat0,debug"
       write(*,*) boxhalfx,boxhalfy,boxlat0,debug
     endif

!
! define map projection (Lambert Conformal from HRRR
!

     lat1 = 38.5
     lon1 = -97.5
     truelat1 = 38.5
     truelat2 = 38.5
     stdlon = -97.5
     dx = 3000.
! HRRR grid
!     nlat = 1060
!     nlon = 1800
! Slightly larger grid for MPAS
     nlat = 1200
     nlon = 2000
     knowni = 0.5 * nlon - 1.
     knownj = 0.5 * nlat - 1.

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

! get corners of map projection

     call ij_to_latlon(proj, 0., 0., ll_lat, ll_lon)
     call ij_to_latlon(proj, real(nlon), real(nlat), ur_lat, ur_lon)
     call ij_to_latlon(proj, 0., real(knowni), left_lat, left_lon)
     call ij_to_latlon(proj, real(knownj), 0., bot_lat, bot_lon)
     call ij_to_latlon(proj, real(nlon), knowni, right_lat, right_lon)
     call ij_to_latlon(proj, real(knownj), real(nlat), top_lat, top_lon)

     write(6,*)
     write(6,*) 'map projection:'
     write(6,*) 'llcrnr    =', ll_lat, ll_lon
     write(6,*) 'urcrnr    =', ur_lat, ur_lon
     write(6,*) 'left ctr  =', left_lat, left_lon
     write(6,*) 'bot ctr   =', bot_lat, bot_lon
     write(6,*) 'right ctr =', right_lat, right_lon
     write(6,*) 'top ctr   =', top_lat, top_lon

!
! read in MPAS mesh information
!

     meshfile='mesh.nc'
     call read_MPAS_nCell(meshfile, nCell)
     write(6,*)
     write(6,*) 'model nCell   =', nCell
     allocate(lat_m(nCell))
     allocate(lon_m(nCell))
     call read_MPAS_lat_lon(meshfile, nCell, lat_m, lon_m)
     write(6,*) 'min model lat =', minval(lat_m)
     write(6,*) 'min model lon =', minval(lon_m)
     write(6,*) 'max model lat =', maxval(lat_m)
     write(6,*) 'max model lon =', maxval(lon_m)
     write(6,*)

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

!
! DEBUGGING: Write out raw NASA LaRC data to text file
!
     if (debug > 0) then
       open(12, file="nasa_larc_obs_raw.txt", status="new", action="write")
       write(12,'(6a12)') 'lat', 'lon', 'ptop', 'teff', 'lwp', 'phase'
       do j=1,numobs
          write(12,'(4f12.3,f12.4,I12)') lat_l(j),lon_l(j),ptop_l(j),teff_l(j),lwp_l(j),phase_l(j)
       enddo
       close(12)
     endif

! -----------------------------------------------------------
! -----------------------------------------------------------
!     Map each FOV onto model mesh
! -----------------------------------------------------------
! -----------------------------------------------------------
!
     allocate (Pxx(nlon,nlat,nfov),Txx(nlon,nlat,nfov),WPxx(nlon,nlat,nfov))
     allocate (xdist(nlon,nlat,nfov), xxxdist(nfov))
     allocate (latxx(nlon,nlat,nfov), lonxx(nlon,nlat,nfov))
     allocate (PHxx(nlon,nlat,nfov),index(nlon,nlat), jndex(nfov))
     index=0

     do ipt=1,numobs
       if (phase_l(ipt).ge.0) then
!  Indicates there is some data (not missing)

!  to determine npts
         nptsx=npts_rad !by default
         nptsy=npts_rad !by default
         if (lat_l(ipt) > boxlat0(1) ) then
           do ibox=1,boxMAX !to get the largest possible npts
             if (lat_l(ipt) > boxlat0(ibox)) then
               if (boxhalfx(ibox)>0) nptsx=boxhalfx(ibox)
               if (boxhalfy(ibox)>0) nptsy=boxhalfy(ibox)
             endif
           enddo
         endif

! * Compute RR grid x/y at lat/lon of cloud data
         rlon=lon_l(ipt)
         rlat=lat_l(ipt)
         call latlon_to_ij(proj,rlat,rlon,xc,yc)

         ii1 = int(xc+0.5)
         jj1 = int(yc+0.5)
         if ( (jj1-1.ge.1 .and. jj1+1.le.nlat) .and.  &
              (ii1-1.ge.1 .and. ii1+1.le.nlon) )    then
            do jj = max(1,jj1-nptsy), min(nlat,jj1+nptsy)
               do ii = max(1,ii1-nptsx), min(nlon,ii1+nptsx)
! * We check multiple data within gridbox
                  if (index(ii,jj).lt.nfov) then
                      index(ii,jj) = index(ii,jj) + 1

                     Pxx(ii,jj,index(ii,jj))   = Ptop_l(ipt)
                     Txx(ii,jj,index(ii,jj))   = Teff_l(ipt)
                     PHxx(ii,jj,index(ii,jj))  = phase_l(ipt)
                     WPxx(ii,jj,index(ii,jj))  = lwp_l(ipt)
                     latxx(ii,jj,index(ii,jj)) = rlat
                     lonxx(ii,jj,index(ii,jj)) = rlon
                     !xdist(ii,jj,index(ii,jj)) =       &
                     !    sqrt( (XC+1-ii)**2 + (YC+1-jj)**2 )
                  else
                     write(6,*) 'ALERT: too many data in one grid, increase nfov'
                     write(6,*) nfov, ii,jj
                  endif
               enddo ! ii
            enddo  ! jj
         endif ! obs is inside the analysis domain
!
       endif   ! phase_l >= 0
     enddo   ! ipt

     deallocate(lat_l,lon_l,ptop_l,teff_l,phase_l,lwp_l)
     write(6,*) 'The max index number is: ', maxval(index)

!
!  Now, map the observations to the MPAS mesh
!
     allocate(w_pcld(nCell))
     allocate(w_tcld(nCell))
     allocate(w_frac(nCell))
     allocate(w_lwp(nCell))
     allocate(nlev_cld(nCell))
     w_pcld=99999.
     w_tcld=99999.
     w_frac=99999.
     w_lwp=99999.
     nlev_cld = 99999

     noutside = 0

     do im=1,nCell

! Compute map projection x/y at lat/lon of MPAS cell
       rlon=lon_m(im)
       rlat=lat_m(im)
       call latlon_to_ij(proj,rlat,rlon,xc,yc)
      
! Determine closest x/y integer coordinate to MPAS cell 
       ii1 = int(xc+0.5)
       jj1 = int(yc+0.5)

       if (ii1 < 1 .or. jj1 < 1 .or. ii1 > nlon .or. jj1 > nlat) then
         noutside = noutside + 1
       else

         if ((index(ii1,jj1) >= 1 .and. index(ii1,jj1) < 3) .and. userDX < 7000.0) then
            w_pcld(im) = Pxx(ii1,jj1,1) ! hPa
            w_tcld(im) = Txx(ii1,jj1,1) ! K
            w_lwp(im) = WPxx(ii1,jj1,1) ! g/m^2
            w_frac(im) = 1
            nlev_cld(im) = 1
            if (w_pcld(im).eq.-20) then
                 w_pcld(im) = 1013. ! hPa - no cloud
                 w_frac(im)=0.0
                 nlev_cld(im) = 0
            end if
         elseif(index(ii1,jj1) .ge. 3) then

! * We decided to use nearest neighborhood for ECA values,
! *     a kind of convective signal from GOES platform...
!
! * Sort to find closest distance if more than one sample
            if(ioption == 1) then    !nearest neighborhood
              do i=1,index(ii1,jj1)
                jndex(i) = i
                call compute_haversine_dist(rlat, rlon, latxx(ii1,jj1,i), lonxx(ii1,jj1,i), dist) 
                xxxdist(i) = dist
                !xxxdist(i) = xdist(ii1,jj1,i)
              enddo
              call sortmed(xxxdist,index(ii1,jj1),jndex,fr)
              w_pcld(im) = Pxx(ii1,jj1,jndex(1))
              w_tcld(im) = Txx(ii1,jj1,jndex(1))
              w_lwp(im) = WPxx(ii1,jj1,jndex(1))
            endif
! * Sort to find median value 
            if(ioption .eq. 2) then    !pick median 
              do i=1,index(ii1,jj1)
                jndex(i) = i
                xxxdist(i) = Pxx(ii1,jj1,i)
              enddo
              call sortmed(xxxdist,index(ii1,jj1),jndex,fr)
              med_pt = index(ii1,jj1)/2  + 1
              w_pcld(im) = Pxx(ii1,jj1,jndex(med_pt)) ! hPa
              w_tcld(im) = Txx(ii1,jj1,jndex(med_pt)) ! K
              w_lwp(im) = WPxx(ii1,jj1,jndex(med_pt)) !  g/m^2
            endif   ! pick median
!
! missing pcld
            if (w_pcld(im).eq.-20) then
               w_pcld(im) = 1013. ! hPa - no cloud
               w_frac(im)=0.0
               nlev_cld(im) = 0
! cloud fraction based on phase (0 are clear), what about -9 ????
            elseif( w_pcld(im) < 1012.99) then
               cfov = 0
               do i=1,index(ii1,jj1)
                 if(PHxx(ii1,jj1,i) .gt. 0.1) cfov = cfov + 1
               enddo
               w_frac(im) = float(cfov)/(max(1,index(ii1,jj1)))     !  fraction
               if( w_frac(im) > 0.01 ) nlev_cld(im) = 1
            endif
         endif   ! index > 3
       endif   ! check to see if index in map projection domain
     enddo  !im

     deallocate (Pxx,Txx,WPxx)
     deallocate (xdist, xxxdist)
     deallocate (PHxx, jndex)
!
!  write results
!

     write(6,*)
     write(6,*) 'number of MPAS cells outside of map projection =', noutside
     write(6,*) 'percentage of MPAS cells outside of map projection =', 100. * real(noutside) / real(nCell)
     write(6,*)
     write(6,'(8a12)') 'im', 'lat', 'lon', 'w_pcld', 'w_tcld', 'w_frac', 'w_lwp ', 'nlev_cld'
     do im=1,nCell,nCell/50
        write(6,'(I12,5f12.3,f12.4,I12)') im,lat_m(im),lon_m(im),w_pcld(im),&
                w_tcld(im),w_frac(im),w_lwp(im),nlev_cld(im)
     enddo
!
     open(15, file='NASALaRC_cloud4fv3.bin',form='unformatted')
        write(15)  nlon,nlat,nCell
        write(15)  index
        write(15)  w_pcld
        write(15)  w_tcld
        write(15)  w_frac
        write(15)  w_lwp
        write(15)  nlev_cld
     close(15)
!
! DEBUGGING: Write out interpolated NASA LaRC data to text file
!
     if (debug > 0) then
       open(13, file="nasa_larc_obs_interp.txt", status="new", action="write")
       write(13,'(7a12)') 'lat', 'lon', 'ptop', 'teff', 'lwp', 'frac', 'nlev_cld'
       do j=1,nCell
          write(13,'(4f12.3,f12.4,f12.4,I12)') lat_m(j),lon_m(j),w_pcld(j),w_tcld(j),w_lwp(j),w_frac(j),nlev_cld(j)
       enddo
       close(13)
     endif

!=========================================================================
! There should be additional code here that has not been added yet. See 
! process_NASALaRC_cloud.f90
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
