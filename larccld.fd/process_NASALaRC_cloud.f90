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
  use mpasio, only: read_MPAS_dim,read_MPAS_lat_lon
  use map_utils
  use misc_definitions_module
  use larccld_utils, only: sortmed,compute_haversine_dist
  use larccld_bufr_io, only: read_NASALaRC_cloud_bufr,read_NASALaRC_cloud_bufr_survey
  use larccld_bufr_io, only: write_bufr_NASALaRC

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
! MPAS mesh
  integer(i_kind) :: nCell
  REAL(r_single), allocatable :: lat_m(:)
  REAL(r_single), allocatable :: lon_m(:)
  CHARACTER*180   meshfile
  CHARACTER*50    mpasfield
!
!  For NASA LaRC 
!
  CHARACTER*80   satfile
  
!     ****VARIABLES FOR THIS NETCDF FILE****
!
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
  real, allocatable ::     Pxx(:,:,:),Txx(:,:,:),WPxx(:,:,:)
  real,allocatable  ::     xdist(:,:,:), xxxdist(:)
  real, allocatable ::     latxx(:,:,:),lonxx(:,:,:)
  integer,allocatable  ::  PHxx(:,:,:),index(:,:), jndex(:)
  integer  ii,jj,med_pt,ii1,jj1,im

!
! namelist
!
  integer            :: analysis_time
  integer            :: ioption
  character(len=100) :: bufrfile
  integer(i_kind)    :: npts_rad, nptsx, nptsy
  integer(i_kind)    :: boxhalfx(boxMAX), boxhalfy(boxMAX)
  real (r_kind)      :: boxlat0(boxMAX)
  real (r_kind)      :: userDX
  integer            :: debug
  namelist/setup/ analysis_time, ioption, npts_rad,bufrfile, &
                  boxhalfx, boxhalfy, boxlat0,userDX,debug
!
!
!  ** misc
      
  real :: rlat,rlon
  real :: xc,yc
  real :: dist

  integer i,j,ipt,cfov,ibox

  integer(8) :: east_time, west_time
  integer :: maxobs,numobs

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
! define map projection (Lambert Conformal, slightly larger than the HRRR domain)
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
     mpasfield='nCells'
     call read_MPAS_dim(meshfile, mpasfield, nCell)
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
     allocate(w_pcld(1,nCell))
     allocate(w_tcld(1,nCell))
     allocate(w_frac(1,nCell))
     allocate(w_lwp(1,nCell))
     allocate(nlev_cld(1,nCell))
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
            w_pcld(1,im) = Pxx(ii1,jj1,1) ! hPa
            w_tcld(1,im) = Txx(ii1,jj1,1) ! K
            w_lwp(1,im) = WPxx(ii1,jj1,1) ! g/m^2
            w_frac(1,im) = 1
            nlev_cld(1,im) = 1
            if (w_pcld(1,im).eq.-20) then
                 w_pcld(1,im) = 1013. ! hPa - no cloud
                 w_frac(1,im)=0.0
                 nlev_cld(1,im) = 0
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
              call sortmed(xxxdist,index(ii1,jj1),jndex)
              w_pcld(1,im) = Pxx(ii1,jj1,jndex(1))
              w_tcld(1,im) = Txx(ii1,jj1,jndex(1))
              w_lwp(1,im) = WPxx(ii1,jj1,jndex(1))
            endif
! * Sort to find median value 
            if(ioption .eq. 2) then    !pick median 
              do i=1,index(ii1,jj1)
                jndex(i) = i
                xxxdist(i) = Pxx(ii1,jj1,i)
              enddo
              call sortmed(xxxdist,index(ii1,jj1),jndex)
              med_pt = index(ii1,jj1)/2  + 1
              w_pcld(1,im) = Pxx(ii1,jj1,jndex(med_pt)) ! hPa
              w_tcld(1,im) = Txx(ii1,jj1,jndex(med_pt)) ! K
              w_lwp(1,im) = WPxx(ii1,jj1,jndex(med_pt)) !  g/m^2
            endif   ! pick median
!
! missing pcld
            if (w_pcld(1,im).eq.-20) then
               w_pcld(1,im) = 1013. ! hPa - no cloud
               w_frac(1,im)=0.0
               nlev_cld(1,im) = 0
! cloud fraction based on phase (0 are clear), what about -9 ????
            elseif( w_pcld(1,im) < 1012.99) then
               cfov = 0
               do i=1,index(ii1,jj1)
                 if(PHxx(ii1,jj1,i) .gt. 0.1) cfov = cfov + 1
               enddo
               w_frac(1,im) = float(cfov)/(max(1,index(ii1,jj1)))     !  fraction
               if( w_frac(1,im) > 0.01 ) nlev_cld(1,im) = 1
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
        write(6,'(I12,5f12.3,f12.4,I12)') im,lat_m(im),lon_m(im),w_pcld(1,im),&
                w_tcld(1,im),w_frac(1,im),w_lwp(1,im),nlev_cld(1,im)
     enddo
!
     open(15, file='NASALaRC_cloud4mpas.bin',form='unformatted')
        write(15)  1,nCell
        write(15)  index
        write(15)  w_pcld
        write(15)  w_tcld
        write(15)  w_frac
        write(15)  w_lwp
        write(15)  nlev_cld
     close(15)
!
!  write out results
!
     call write_bufr_NASALaRC(bufrfile,analysis_time,1,nCell,userDX,index,w_pcld,w_tcld,w_frac,w_lwp,nlev_cld)
!
! DEBUGGING: Write out interpolated NASA LaRC data to text file
!
     if (debug > 0) then
       open(13, file="nasa_larc_obs_interp.txt", status="new", action="write")
       write(13,'(7a12)') 'lat', 'lon', 'ptop', 'teff', 'lwp', 'frac', 'nlev_cld'
       do j=1,nCell
          write(13,'(4f12.3,f12.4,f12.4,I12)') lat_m(j),lon_m(j),w_pcld(1,j),w_tcld(1,j),w_lwp(1,j),w_frac(1,j),nlev_cld(1,j)
       enddo
       close(13)
     endif

!
     deallocate(w_pcld)
     deallocate(w_tcld)
     deallocate(w_frac)
     deallocate(w_lwp)
     deallocate(nlev_cld)
     deallocate(index)
!
     write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="
  endif ! if mype==0 

  call MPI_FINALIZE(ierror)

end program process_NASALaRC_cloud
