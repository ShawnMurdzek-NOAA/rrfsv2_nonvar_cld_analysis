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
  use mpasio, only: read_MPAS_dim,read_MPAS_lat_lon
  use map_utils
  use map_proj_helper, only: init_proj,write_corners
  use lightning_bufr_io, only: read_lightning_bufr,write_bufr_lightning
  use check_lght_qc, only: Check_Lightning_QC,Check_NLDN
  use netCDFsub_lightning, only: get_dim_att_nldn,ifexist_file,get_lightning_nldn

  implicit none
  INCLUDE 'netcdf.inc'
!
!
! MPI variables
  integer :: npe, mype, ierror

! Map projection
  type(proj_info) :: proj

!  For MPAS mesh
  integer(i_kind) :: nCell
  real, allocatable :: lat_m(:),lon_m(:)
  CHARACTER*180   meshfile
  CHARACTER*50    mpasfield

!
!  For lightning data
!
  integer,parameter ::    max_numStrike=1000000
  integer ::    numStrike
  character*180   lightsngle,lght_out_file
  real,allocatable:: llon(:)    !
  real,allocatable:: llat(:)    !
  integer,allocatable:: ltime(:)   !
  integer,allocatable:: lstrike(:) !
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
  character(len=25) :: proj_name
  namelist/setup/analysis_time,minute,trange_start,trange_end,&
                 obs_type,NLDN_filenum,IfAlaska,proj_name,debug
!
!  ** misc
  integer, allocatable :: cell_id(:,:,:),lght_id(:,:,:),index_m(:,:),index_l(:,:)
  real, allocatable :: x_m(:),y_m(:),x_l(:),y_l(:)
  integer :: nCell_mp,nlght_mp
  real :: rlon,rlat,xc,yc,d,d2
  integer i,j,ii,jj,ilght,nt,icell,c_id,l_id,ic,jc,nearest_id
  integer :: istatus
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
     proj_name='CONUS'
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

     call init_proj(proj, proj_name)
     call write_corners(proj)

!
! get model domain dimension
!
     meshfile='mesh.nc'
     mpasfield='nCells'
     call read_MPAS_dim(meshfile, mpasfield, nCell)
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
    nCell_mp=8
    allocate(cell_id(proj%nlon,proj%nlat,nCell_mp))
    allocate(index_m(proj%nlon,proj%nlat))
    allocate(x_m(nCell))
    allocate(y_m(nCell))
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
      if ( (ic.ge.1 .and. ic.le.proj%nlon) .and. &
           (jc.ge.1 .and. jc.le.proj%nlat) ) then
        if ( index_m(ic,jc).lt.nCell_mp ) then
          index_m(ic,jc) = index_m(ic,jc) + 1
          cell_id(ic,jc,index_m(ic,jc)) = i
        else
          write(*,*)
          write(*,*) 'ERROR: nCell_mp exceeded when mapping MPAS to the map projection.'
          write(*,*) 'Please increase nCell_mp.'
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

!
! DEBUGGING: Write out raw lightning obs
!
        if (debug > 0) then
          write(lght_out_file,'(a,I1,a)') 'lightning_raw_',nt,'.txt'
          open(12, file=lght_out_file, status="new", action="write")
          write(12,'(5a12)') 'lat', 'lon', 'time', 'lstrike', 'quality'
          do j=1,numStrike
            write(12,'(2f12.3,3I12)') llat(j),llon(j),ltime(j),lstrike(j),lquality(j)
          enddo
          close(12)
        endif


! -----------------------------------------------------------
! -----------------------------------------------------------
!     Map each lightning strike to the closest MPAS cell
! -----------------------------------------------------------
! -----------------------------------------------------------
!
        nlght_mp=2000
        allocate(lght_id(proj%nlon,proj%nlat,nlght_mp))
        allocate(index_l(proj%nlon,proj%nlat))
        allocate(x_l(numStrike))
        allocate(y_l(numStrike))
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
            if ( (ic.ge.1 .and. ic.le.proj%nlon) .and. &
                 (jc.ge.1 .and. jc.le.proj%nlat) ) then
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
        do i=1,proj%nlon
          do j=1,proj%nlat
            if (index_l(i,j).gt.0) then
              do ilght=1,index_l(i,j)
                l_id = lght_id(i,j,ilght)
                d = 1.e9
                nearest_id = -99
                do ii=max(1,i-1), min(proj%nlon,i+1)
                  do jj=max(1,j-1), min(proj%nlon,j+1)
                    if (index_m(ii,jj).gt.0) then
                      do icell=1,index_m(ii,jj)
                        c_id = cell_id(ii,jj,icell)
                        d2 = ((x_m(c_id) - x_l(l_id))**2 + &
                              (y_m(c_id) - y_l(l_id))**2)
                        if (d2.lt.d) then
                          nearest_id = c_id
                          d = d2
                        endif 
                      enddo ! icell
                    endif
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
!  Entries in lightning_out include:
!    1: Always 1 (formerly longitude)
!    2: Cell index (formerly latitude)
!    3: Time
!    4: Lightning number

     allocate(lightning_out(4,nCell))
     numlightning=0
     do i=1,nCell
       if(lightning(i) > 0 ) then
         numlightning=numlightning+1
         lightning_out(1,numlightning)=1
         lightning_out(2,numlightning)=float(i)
         lightning_out(3,numlightning)=float(minute)/60.0
         lightning_out(4,numlightning)=lightning(i)
         if(lightning_out(4,numlightning) > 1000.0 ) then
            lightning_out(4,numlightning)=1000.0
            write(6,*) 'high lightning strokes=',lightning(i),i
         endif
!        write(*,*) numlightning,i,j,lightning(i,j)
       endif
     enddo

     write(*,*) 'Write out results for MPAS:',numlightning
     OPEN(10,file='LightningInMPAS.dat',form='unformatted')
      write(10) 4,1,nCell,numlightning,1,2
      write(10) ((real(lightning_out(i,j)),i=1,4),j=1,numlightning)
      write(10) lightning
     close(10)

     write(6,*) ' write lightning in BUFR for cycle time ',idate
     call write_bufr_lightning(1,nCell,numlightning,lightning_out,idate)

!
! DEBUGGING: Write out interpolated lightning data to text file
!
     if (debug > 0) then
       open(13, file="lightning_interp.txt", status="new", action="write")
       write(13,'(4a12)') 'cell_id', 'lat', 'lon', 'nstrikes'
       do j=1,numlightning
          c_id = int(lightning_out(2,j))
          write(13,'(I12,2f12.3,I12)') c_id,lat_m(c_id),lon_m(c_id),int(lightning_out(4,j))
       enddo
       close(13)
     endif

     deallocate(lightning_out)

  endif ! mype
  call MPI_FINALIZE(ierror)
!
end program process_Lightning 
