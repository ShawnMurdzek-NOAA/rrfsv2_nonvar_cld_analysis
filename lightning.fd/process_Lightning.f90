program process_Lightning
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2022-01-18
!
! ABSTRACT: 
!     This routine read in lightning data and 
!     map them into FV3LAM ESG grid
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  NLDN lightning data
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
  use kinds, only: r_kind,i_kind
  use mpasio, only: read_MPAS_nCell,read_MPAS_lat_lon

! Modules from WPS
  use map_utils
  use misc_definitions_module

  implicit none
  INCLUDE 'netcdf.inc'
!
!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror

!
  character*256 output_file

! Map projection
  type(proj_info) :: proj
  integer :: nlat, nlon
  real :: lat1, lon1, truelat1, truelat2, stdlon, dx, knowni, knownj
  real :: ll_lat, ur_lat, left_lat, bot_lat, right_lat, top_lat
  real :: ll_lon, ur_lon, left_lon, bot_lon, right_lon, top_lon

!  For MPAS mesh
  integer(i_kind) :: nCell
  real, allocatable :: lat_m(:),lon_m(:)
  CHARACTER*180   meshfile

!
!  For lightning data
!
  integer,parameter ::    max_numStrike=1000000
  integer ::    numStrike
  character*180   lightsngle
  real,allocatable:: llon(:)    !
  real,allocatable:: llat(:)    !
  integer,allocatable:: ltime(:)   !
  integer,allocatable:: lstrike(:) !
  character*21,allocatable:: ctime(:)   !
  real :: rtmp
  integer,allocatable:: lquality(:) !

  real, allocatable :: lightning(:)   ! lightning  strakes
  real(r_kind), allocatable :: lightning_out(:,:)   ! lightning  strakes

  integer :: numNLDN_all, numNLDN_used
!
!! Declare namelists 
!
! SETUP (general control namelist) :
!
  character*10 :: analysis_time
  real :: trange_start,trange_end
  integer :: minute,debug
  character(len=20) :: obs_type

  integer      :: NLDN_filenum
  logical      :: IfAlaska
  character(len=20) :: grid_type
  namelist/setup/analysis_time,minute,trange_start,trange_end,&
                 obs_type,grid_type,NLDN_filenum,IfAlaska,debug
!
!  ** misc
      
  CHARACTER*180   workpath

  integer, allocatable :: cell_id(:,:,:),lght_id(:,:,:),index_m(:,:),index_l(:,:)
  real, allocatable :: x_m(:),y_m(:),x_l(:),y_l(:)
  integer :: nCell_mp,nlght_mp
  real :: rlon,rlat,xc,yc
  integer i,j,ilght,nt,icell,c_id,l_id
  integer :: ii,jj
  real :: d,d2
  integer :: NCID, istatus
  integer :: numlightning,idate,filenum
  logical :: ifexist


!**********************************************************************
!
!            END OF DECLARATIONS....start of program
! MPI setup
  call MPI_INIT(ierror) 
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  if(mype==0) then
     numNLDN_all=0
     numNLDN_used=0

     NLDN_filenum=1  
     trange_start=0
     trange_end=0
     minute=0
     obs_type="none"
     grid_type="none"
     debug=0
     inquire(file='namelist.lightning', EXIST=ifexist )
     if(ifexist) then
        open(15, file='namelist.lightning')
          read(15,setup)
        close(15)
     else
        write(*,*) "ERROR: cannot find namelist file"
        stop 123
     endif

!
! define map projection (Lambert Conformal, similar to HRRR)
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

     write(*,*)
     write(*,*) 'map projection:'
     write(*,*) 'llcrnr    =', ll_lat, ll_lon
     write(*,*) 'urcrnr    =', ur_lat, ur_lon
     write(*,*) 'left ctr  =', left_lat, left_lon
     write(*,*) 'bot ctr   =', bot_lat, bot_lon
     write(*,*) 'right ctr =', right_lat, right_lon
     write(*,*) 'top ctr   =', top_lat, top_lon

!
! get model domain dimension
!
     meshfile='mesh.nc'
     call read_MPAS_nCell(meshfile, nCell)
     allocate(lat_m(nCell))
     allocate(lon_m(nCell))
     call read_MPAS_lat_lon(meshfile, nCell, lat_m, lon_m)
     write(*,*)
     write(*,*) 'model nCell   =', nCell
     write(*,*) 'min model lat =', minval(lat_m)
     write(*,*) 'min model lon =', minval(lon_m)
     write(*,*) 'max model lat =', maxval(lat_m)
     write(*,*) 'max model lon =', maxval(lon_m)
     write(*,*)

!
! Map each MPAS cell to the closest map projection integer coordinate
!
    nCell_mp=10
    allocate(cell_id(nlon,nlat,nCell_mp))
    allocate(index_m(nlon,nlat))
    allocate(x_m(nCell_mp))
    allocate(y_m(nCell_mp))
    cell_id = -99
    index_m = 0

    do i=1,nCell
      rlon=lon_m(i)
      rlat=lat_m(i)
      call latlon_to_ij(proj,rlat,rlon,xc,yc)
      x_m(i) = xc
      y_m(i) = yc
      ic = int(xc+0.5)
      jc = int(yc+0.5)
      if ( (ic.ge.1 .and. ic.le.nlon) .and. &
           (jc.ge.1 .and. jc.le.nlat) ) then
        if ( index_m(ic,jc).lt.nCell_mp ) then
          index_m(ic,jc) = index_m(ic,jc) + 1
          cell_id(ic,jc,index_m(ic,jc)) = i
        else
          write(*,*)
          write(*,*) 'ERROR: nCell_mp exceeded when mapping MPAS to the map projection. Please increase nCell_mp'
          stop 
        endif
      endif
    enddo 

    if (debug.gt.0) then
      write(*,*)
      write(*,*) 'Max number of MPAS cells mapped to a single map projection point =', maxval(index_m)
      write(*,*)
    endif

!
! Read lightning observations
!

     allocate(lightning(nCell))
     lightning=0

     if(obs_type=="bufr") then
        filenum=1
     elseif(obs_type=="nldn_nc") then
        filenum=NLDN_filenum
     else
        write(*,*) 'ERROR: unknown data type:',obs_type
     endif
!
     read(analysis_time,'(I10)') idate
     do nt=1,filenum
        if(obs_type=="bufr") then
   
           allocate(llon(max_numStrike))
           allocate(llat(max_numStrike))
           allocate(ltime(max_numStrike))
           allocate(lStrike(max_numStrike))

           lightsngle='lghtngbufr'
           call read_lightning_bufr(lightsngle,max_numStrike,analysis_time,minute,trange_start,trange_end,&
                             numStrike,llon,llat,ltime,lStrike)
           numNLDN_all=numStrike

           allocate(lquality(numStrike))
           lquality = 0    ! 0 good data,  > 0 bad data
           call Check_Lightning_QC(numStrike,llon,llat,ltime,lstrike,lquality)

!
!  process NLDN data
!
        elseif(obs_type=="nldn_nc") then
           write(lightsngle,'(a,I1)') 'NLDN_lightning_',nt
           write(*,*) trim(lightsngle)
           call ifexist_file(trim(lightsngle),istatus)
           if (ISTATUS .NE. NF_NOERR) CYCLE

           call GET_DIM_ATT_NLDN(lightsngle,numStrike)
           write(*,*) 'number of strikes=', nt, numStrike

           allocate(llon(numStrike))
           allocate(llat(numStrike))
           allocate(ltime(numStrike))
           allocate(lStrike(numStrike))

           call GET_lightning_NLDN(lightsngle,numStrike,llon,llat,ltime,lStrike)
! check quality
           allocate(lquality(numStrike))
           lquality = 0    ! 0 good data,  > 0 bad data
           call Check_NLDN(numStrike,llon,llat,ltime,lstrike,lquality)
        endif

! -----------------------------------------------------------
! -----------------------------------------------------------
!     Map each lightning strike to the closest MPAS cell
! -----------------------------------------------------------
! -----------------------------------------------------------
!
        nlght_mp=2000
        allocate(lght_id(nlon,nlat,nlght_mp))
        allocate(index_l(nlon,nlat))
        allocate(x_l(nlght_mp))
        allocate(y_l(nlght_mp))
        lght_id = -99
        index_l = 0

        do i=1,numStrike
          if(lquality(i) == 0) then
            rlon=llon(i)
            rlat=llat(i)
            call latlon_to_ij(proj,rlat,rlon,xc,yc)
            x_l(i) = xc
            y_l(i) = yc
            ic = int(xc+0.5)
            jc = int(yc+0.5)
            if ( (ic.ge.1 .and. ic.le.nlon) .and. &
                 (jc.ge.1 .and. jc.le.nlat) ) then
              if ( index_l(ic,jc).lt.nlght_mp ) then
                index_l(ic,jc) = index_l(ic,jc) + 1
                lght_id(ic,jc,index_l(ic,jc)) = i
              else
                write(*,*)
                write(*,*) 'WARNING: nlght_mp exceeded when mapping lightning to the map projection'
                write(*,*) 'nlght_mp =', nlght_mp
                write(*,*) 'map projection location =', ic, jc
              endif
            endif
          endif
        enddo

        if (debug.gt.0) then
          write(*,*)
          write(*,*) 'Max number of lightning strikes mapped to a single map projection point =', maxval(index_l)
          write(*,*)
        endif

        deallocate(ltime)
        deallocate(lStrike)
        deallocate(lquality)

!
!  Use nearest-neighbor to interpolate lightning to MPAS mesh
!
        do i=1,nlon
          do j=1,nlat
            if (index_l(i,j).gt.0) then
              do ilght=1,index_l(i,j)
                d = 1.e9
                nearest_id = -99
                do ii=max(1,i-1), min(nlon,i+1)
                  do jj=max(1,j-1), min(nlon,j+1)
                    do icell=1,index_m(ii,jj)
                      c_id = cell_id(ii,jj,icell)
                      l_id = lght_id(ii,jj,ilght)
                      d2 = ((x_m(c_id) - x_l(l_id))**2 + &
                            (y_m(c_id) - y_l(l_id))**2)
                      if (d2.lt.d) then
                        nearest_id = c_id
                      endif 
                    enddo ! icell
                  enddo ! jj
                enddo ! ii
                if (nearest_id.gt.0) then
                  lightning(nearest_id) = lightning(nearest_id) + 1
                  numNLDN_used=numNLDN_used+1
                endif
              enddo ! ilght
            endif ! index_l > 0
          enddo ! nlat
        enddo ! nlon

        deallocate(llon)
        deallocate(llat)
        deallocate(lght_id)
        deallocate(index_l)
        deallocate(x_l)
        deallocate(y_l)
     enddo ! nt
!
     deallocate(cell_id)
     deallocate(index_m)
     deallocate(x_m)
     deallocate(y_m)
!
!  report statistic
!
     write(*,*) ' The total number of lightning obs is:', numNLDN_all
     write(*,*) ' The number of obs used is:', numNLDN_used

!  Write out results

     allocate(lightning_out(3,nCell))
     numlightning=0
     do j=1,nCell
       if(lightning(i) > 0 ) then
         numlightning=numlightning+1
         lightning_out(1,numlightning)=float(i)
         lightning_out(2,numlightning)=float(minute)/60.0
         lightning_out(3,numlightning)=lightning(i)
         if(lightning_out(4,numlightning) > 1000.0 ) then
            lightning_out(4,numlightning)=1000.0
            write(6,*) 'high lightning strokes=',lightning(i),i
         endif
!        write(*,*) numlightning,i,j,lightning(i,j)
       endif
     enddo

     write(*,*) 'Write out results for MPAS:',numlightning
     OPEN(10,file='LightningInMPAS.dat',form='unformatted')
      write(10) 3,nCell,numlightning,1,2
      write(10) ((real(lightning_out(i,j)),i=1,3),j=1,numlightning)
      write(10) lightning
     close(10)

     write(6,*) ' write lightning in BUFR for cycle time ',idate
     call write_bufr_lightning(1,nCell,numlightning,lightning_out,idate)
     deallocate(lightning_out)

  endif ! mype
  call MPI_FINALIZE(ierror)
!
end program process_Lightning 
