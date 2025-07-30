program process_NSSL_mosaic
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2007-12-17
!
! ABSTRACT: 
!     This routine read in NSSL reflectiivty mosaic fiels and 
!     interpolate them into GSI mass grid
!
!     tversion=8  : NSSL 8 tiles netcdf
!     tversion=81 : NCEP 8 tiles binary
!     tversion=4  : NSSL 4 tiles binary
!     tversion=1  : NSSL 1 tile grib2
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  mosaic_files
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
  use module_read_NSSL_refmosaic, only: read_nsslref
  use mpasio, only: read_MPAS_nCell,read_MPAS_lat_lon,read_MPAS_bdyMaskCell

  implicit none
!
  INCLUDE 'netcdf.inc'
!
  type(read_nsslref) :: readref 

!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!
  character*256 output_file
!
!  gridded reflectivity
  REAL, allocatable :: ref3d(:,:)   ! 3D reflectivity
  REAL, allocatable :: ref0(:,:)   ! 3D reflectivity
  REAL(r_kind), allocatable :: ref3d_column(:,:)   ! 3D reflectivity in column
!
!  MPAS mesh
  integer(i_kind) :: nCell
  real, allocatable :: lat_m(:),lon_m(:)
  integer, allocatable :: bdyMask(:)
  CHARACTER*180   meshfile
!
!  namelist files
!
  integer      ::  tversion
  character*10 :: analysis_time
  CHARACTER*180   dataPath
  namelist/setup/ tversion,analysis_time,dataPath
  integer(i_kind)  ::  idate
!
!  namelist and other variables for netcdf output
!
!  output_netcdf             logical controlling whether netcdf output file should be created
!  max_height                maximum height (m MSL) for data to be retained
!  use_clear_air_type        logical controlling whether to output clear-air (non-precipitation) reflectivity obs
!  precip_dbz_thresh         threshold (dBZ) for minimum reflectivity that is considered precipitation
!  clear_air_dbz_thresh      threshold (dBZ) for maximum reflectivity that is considered clear air
!  clear_air_dbz_value       value (dBZ) assigned to clear-air reflectivity obs
!  precip_dbz_horiz_skip     horizontal thinning factor for reflectivity data in precipitation
!  precip_dbz_vert_skip      vertical thinning factor for reflectivity data in precipitation
!  clear_air_dbz_horiz_skip  horizontal thinning factor for clear air reflectivity data
!  clear_air_dbz_vert_skip   vertical thinning factor for clear air reflectivity data
!
  logical :: output_netcdf = .false.
  real :: max_height = 20000.0
  logical :: use_clear_air_type = .false.
  real :: precip_dbz_thresh = 15.0
  real :: clear_air_dbz_thresh = 0.0
  real :: clear_air_dbz_value = 0.0
  integer :: precip_dbz_horiz_skip = 0
  integer :: precip_dbz_vert_skip = 0
  integer :: clear_air_dbz_horiz_skip = 0
  integer :: clear_air_dbz_vert_skip = 0
  namelist/setup_netcdf/ output_netcdf, max_height,                       &
                         use_clear_air_type, precip_dbz_thresh,           &
                         clear_air_dbz_thresh, clear_air_dbz_value,       &
                         precip_dbz_horiz_skip, precip_dbz_vert_skip,     &
                         clear_air_dbz_horiz_skip, clear_air_dbz_vert_skip
  logical, allocatable :: precip_ob(:,:)
  logical, allocatable :: clear_air_ob(:,:)
  integer, parameter :: maxMosaiclvl=33
  real :: height_real(maxMosaiclvl)
  integer :: levelheight(maxMosaiclvl)
  data levelheight /500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500,         &
                    2750, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, &
                    7500, 8000, 8500, 9000, 10000, 11000, 12000, 13000, 14000,  &
                    15000, 16000, 17000, 18000, 19000/
  integer num_precip_obs, num_clear_air_obs, num_obs
  logical :: fileexist
!
!  ** misc
!      
  character*80 outfile
  character*256 outfile_netcdf
  integer i,ii,j,jj,k,kk
  integer :: id
  INTEGER(i_kind)  ::  maxlvl
  INTEGER(i_kind)  ::  numlvl,numref
  integer :: maxcores

!**********************************************************************
!
!            END OF DECLARATIONS....start of program
! MPI setup
  call MPI_INIT(ierror) 
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  if(mype==0) write(*,*) mype, 'deal with mosaic'

  datapath="./"
  open(15, file='namelist.mosaic')
    read(15,setup)
  close(15)

  inquire(file='namelist.mosaic_netcdf', exist=fileexist)
  if(fileexist) then
    if(mype==0) write(*,*) 'reading namelist.mosaic_netcdf'
    open(15, file='namelist.mosaic_netcdf')
    read(15, setup_netcdf)
    close(15)
  endif

  if(mype==0) then
    write(6,*)
    write(6,*) 'tversion = ', tversion
    write(6,*) 'analysis_time = ', analysis_time
    write(6,*) 'dataPath = ', dataPath
    write(6,*) 'output_netcdf = ', output_netcdf
    write(6,*) 'max_height = ', max_height
    write(6,*) 'use_clear_air_type = ', use_clear_air_type
    write(6,*) 'precip_dbz_thresh = ', precip_dbz_thresh
    write(6,*) 'clear_air_dbz_thresh = ', clear_air_dbz_thresh
    write(6,*) 'clear_air_dbz_value = ', clear_air_dbz_value
    write(6,*) 'precip_dbz_horiz_skip = ', precip_dbz_horiz_skip
    write(6,*) 'precip_dbz_vert_skip = ', precip_dbz_vert_skip
    write(6,*) 'clear_air_dbz_horiz_skip = ', clear_air_dbz_horiz_skip
    write(6,*) 'clear_air_dbz_vert_skip = ', clear_air_dbz_vert_skip
    write(6,*)
  endif

!
!  safty check for cores used in this run
!
  read(analysis_time,'(I10)') idate
  if(mype==0) write(6,*) 'cycle time is :', idate

  if( tversion == 8 .or. tversion == 14) then
     maxcores=8
  elseif( tversion == 81 ) then
     maxcores=8
  elseif( tversion == 4 ) then
     maxcores=4
  elseif( tversion == 1 ) then
     maxcores=33
  else
     write(*,*) 'unknow tversion !'
     stop 1234
  endif

  if(mype==0) write(6,*) 'total cores for this run is ',npe
  if(npe < maxcores) then
     write(6,*) 'ERROR, this run must use ',maxcores,' or more cores !!!'
     call MPI_FINALIZE(ierror)
     stop 1234
  endif
!
! read NSSL mosaic 
!
  mypeLocal=mype+1
  call readref%init(tversion,mypeLocal,datapath)
!
! deal with certain tile
!
  call readref%readtile(mypeLocal)
  call mpi_barrier(MPI_COMM_WORLD,ierror)
!
  maxlvl=readref%maxlvl
!
! get model domain dimension
!
  meshfile='mesh.nc'
  call read_MPAS_nCell(meshfile, nCell)
  allocate(lat_m(nCell))
  allocate(lon_m(nCell))
  allocate(bdyMask(nCell))
  call read_MPAS_lat_lon(meshfile, nCell, lat_m, lon_m)
  call read_MPAS_bdyMaskCell(meshfile, nCell, bdyMask)
  if(mype==0) then
    write(6,*)
    write(6,*) 'model nCell   =', nCell
    write(6,*) 'min model lat =', minval(lat_m)
    write(6,*) 'min model lon =', minval(lon_m)
    write(6,*) 'max model lat =', maxval(lat_m)
    write(6,*) 'max model lon =', maxval(lon_m)
    write(6,*)
  endif

  allocate(ref3d(nCell,maxlvl))
  ref3d=-999.0
  call mosaic2grid(readref,nCell,maxlvl,lon_m,lat_m,ref3d)
  call mpi_barrier(MPI_COMM_WORLD,ierror)
!
!  collect data from all processes to root (0)
!
  if(mype==0) then
     allocate( ref0(nCell,maxlvl) )
  endif
  call MPI_REDUCE(ref3d, ref0, nCell*maxlvl, MPI_REAL, MPI_MAX, 0, &
                  MPI_COMM_WORLD, ierror)
  deallocate(ref3d)
!
  if(mype==0) then
     write(outfile,'(a,a)') './', 'RefInGSI3D.dat'
     write(outfile_netcdf,'(a,a)') './', 'Gridded_ref.nc'
     OPEN(10,file=trim(outfile),form='unformatted')
        write(10) maxlvl,nCell
        write(10) ref0
     close(10)
     DO k=1,maxlvl
        write(*,*) k,maxval(ref0(:,k)),minval(ref0(:,k))
     ENDDO

! turn this part off to speed up the process for RRFS.
! Writing to a BUFR file has been adapted for MPAS, but has not been tested. Use at your own risk!!
     if(1==2) then
!
        allocate(ref3d_column(maxlvl+1,nCell))
        ref3d_column=-999.0
        numref=0
        DO i=1,nCell
          numlvl=0
          DO k=1,maxlvl
            if(abs(ref0(i,k)) < 888.0 ) numlvl=numlvl+1
          ENDDO
          if(numlvl > 0 ) then
            numref=numref+1
            ref3d_column(1,numref)=float(i)
            DO k=1,maxlvl
               ref3d_column(1+k,numref)=ref0(i,k)
            ENDDO
          endif
        ENDDO

        write(*,*) 'Dump out results', numref, 'out of', nCell
        OPEN(10,file='./'//'RefInGSI.dat',form='unformatted')
          write(10) maxlvl,nCell,numref,1,2
          write(10) ((ref3d_column(k,i),k=1,maxlvl+2),i=1,numref)
        close(10)
  
        write(*,*) 'Start write_bufr_nsslref'
        call write_bufr_nsslref(maxlvl,nCell,numref,ref3d_column,idate)
        deallocate(ref3d_column)
     endif

     if ( output_netcdf .and. (maxlvl.eq.maxMosaiclvl) ) then

        allocate( precip_ob(nCell,maxlvl) )
        allocate( clear_air_ob(nCell,maxlvl) )

        ! Don't produce any netcdf radar observations along the lateral boundaries
        ! bdyMaskCell = 7 values indicate cells along the edge of the domain
        do i=1,nCell
          if (bdyMask(i) .gt. 6) then
            ref0(i,:) = -999.0
          endif
        enddo

        ! Identify precip and clear-air reflectivity observations
        precip_ob(:,:) = .false.
        clear_air_ob(:,:) = .false.
        num_precip_obs = 0
        num_clear_air_obs = 0
        do i=1,nCell
          do k=1,maxlvl
            if ( (levelheight(k) .le. max_height) .and. (ref0(i,k) .ge. precip_dbz_thresh) ) then
              precip_ob(i,k) = .true.
              num_precip_obs = num_precip_obs + 1
            else if ( use_clear_air_type .and. (levelheight(k) .le. max_height) .and. &
                      (ref0(i,k) .gt. -900.0) .and. (ref0(i,k) .le. clear_air_dbz_thresh) ) then
              clear_air_ob(i,k) = .true.
              ref0(i,k) = clear_air_dbz_value
              num_clear_air_obs = num_clear_air_obs + 1
            endif
          enddo
        enddo
        write(*,*) 'number of precip obs found, before thinning, = ', num_precip_obs
        write(*,*) 'number of clear air obs found, before thinning, = ', num_clear_air_obs

        ! Thin precip reflectivity observations
        if (precip_dbz_vert_skip .gt. 0) then
          do k=1,maxlvl
            if (mod(k-1, precip_dbz_vert_skip+1) .ne. 0) then
              write(*,*) 'Thinning:  removing precip obs at level ', k
              precip_ob(:,k) = .false.
            endif
          enddo
        endif
        !if (precip_dbz_horiz_skip .gt. 0) then
        !  write(*,*) 'Horizontal thinning of precip obs'
        !  do i=1,nCell
        !    if ( ref0(i,1) .gt. -900.0 )  then
        !      do k=1,maxlvl
        !        if (precip_ob(i,k)) then
        !          do jj=max(2, j-precip_dbz_horiz_skip), min(nlat-1, j+precip_dbz_horiz_skip)
        !            do ii=max(2, i-precip_dbz_horiz_skip), min(nlon-1, i+precip_dbz_horiz_skip)
        !              precip_ob(ii,k) = .false.
        !            enddo
        !          enddo
        !          precip_ob(i,k) = .true.
        !        endif
        !      enddo
        !    endif
        !  enddo
        !endif

        ! Thin clear-air reflectivity observations
        if (use_clear_air_type .and. (clear_air_dbz_vert_skip .gt. 0) ) then
          do k=1,maxlvl
            if (mod(k-1, clear_air_dbz_vert_skip+1) .ne. 0) then
              write(*,*) 'Thinning:  removing clear air obs at level ', k
              clear_air_ob(:,k) = .false.
            endif
          enddo
        endif
        !if (use_clear_air_type .and. (clear_air_dbz_vert_skip .lt. 0) ) then
        !  write(*,*) 'Thinning:  removing clear air obs at all but two levels'
        !  clear_air_ob(:, 1:12) = .false.
        !  clear_air_ob(:, 14:21) = .false.
        !  clear_air_ob(:, 23:maxlvl) = .false.
        !endif
        !if (use_clear_air_type .and. (clear_air_dbz_horiz_skip .gt. 0) ) then
        !  do i=1,nCell
        !    if ( ref0(i,1) .gt. -900.0 )  then
        !      do k=1,maxlvl
        !        if (clear_air_ob(i,j,k)) then
        !          do jj=max(2, j-clear_air_dbz_horiz_skip), min(nlat-1, j+clear_air_dbz_horiz_skip)
        !            do ii=max(2, i-clear_air_dbz_horiz_skip), min(nlon-1, i+clear_air_dbz_horiz_skip)
        !              clear_air_ob(ii,k) = .false.
        !            enddo
        !          enddo
        !          clear_air_ob(i,k) = .true.
        !        endif
        !      enddo
        !    endif
        !  enddo
        !endif

        ! Count number of valid obs
        num_obs = 0
        do i=1,nCell
          if ( ref0(i,1) .gt. -900.0 )  then
            do k=1,maxlvl
              if ( precip_ob(i,k) .or. clear_air_ob(i,k) ) then
                num_obs = num_obs + 1
              endif
            enddo
          endif
        enddo
        write(*,*) 'num_obs = ', num_obs

        ! Write obs to netcdf file
        do k=1,maxlvl
          height_real(k) = levelheight(k)
          do i=1,nCell
            if ( .not. precip_ob(i,k) .and. .not. clear_air_ob(i,k) ) then
              ref0(i,k) = -999.0
            endif
          enddo
        enddo
        call write_netcdf_nsslref( outfile_netcdf,maxlvl,nCell,ref0,idate,lon_m,lat_m,height_real )

        write(*,*) 'Finish netcdf output'

        deallocate(precip_ob)
        deallocate(clear_air_ob)

     else if (output_netcdf) then

        write(*,*) 'unknown vertical levels'
        write(*,*) 'maxlvl = ', maxlvl
        write(*,*) 'maxMosaiclvl = ', maxMosaiclvl
        write(*,*) 'no netcdf output'

     endif ! output_netcdf

     deallocate(ref0)

     endif ! mype==0

  deallocate(lon_m)
  deallocate(lat_m)

  call readref%close()

  if(mype==0)  write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="

  call MPI_FINALIZE(ierror)
!
end program process_NSSL_mosaic

subroutine mosaic2grid(readref,nCell,maxlvl,lon_m,lat_m,ref3d)
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-01-20
!
! ABSTRACT: 
!     This routine interpolate NSSL reflectiivty mosaic fields to grid
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:  mosaic_files
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
  use module_read_NSSL_refmosaic, only: read_nsslref

  implicit none
!
  type(read_nsslref),intent(in) :: readref 
  integer,intent(in) :: nCell,maxlvl
  real,intent(in)    :: lon_m(nCell)    !
  real,intent(in)    :: lat_m(nCell)    !
  real,intent(inout) :: ref3d(nCell,maxlvl)   ! 3D reflectivity
!
!  For reflectiivty mosaic
!
  REAL   :: lonMin,latMin,lonMax,latMax
  REAL*8 :: dlon,dlat
!
!  ** misc
!      
  integer i,k,kk
  integer :: tversion

  REAL ::  rlat,rlon
  INTEGER  :: ip,jp,ipp1,jpp1
  REAL ::  rip,rjp
  REAL ::  dip,djp
  REAL ::  w1,w2,w3,w4
  REAL ::  ref1,ref2,ref3,ref4,refl_ltng
  real, allocatable :: mscValue(:,:)

!**********************************************************************
!
!
  tversion=readref%tversion
  dlon=readref%dlon
  dlat=readref%dlat
  latMax=readref%latMax
  lonMax=readref%lonMax
  lonMin=readref%lonMin
  latMin=readref%latMin
!
  if(readref%if_fileexist) then
!
     allocate(mscValue(readref%mscNlon,readref%mscNlat))
     do k=1, readref%mscNlev
          if(tversion == 1) then
             kk=readref%ilevel
             mscValue(:,:) = readref%mscValue3d(:,:,1)
          else
             kk=k
             mscValue(:,:) = readref%mscValue3d(:,:,k)
          endif
!
          DO i=1,nCell
             rlat=lat_m(i)
             rlon=lon_m(i)

             if(tversion == 14 ) then
               rip=(rlon-lonMin)/dlon+1
               rjp=(rlat-latMin)/dlat+1
               ip=int(rip)
               jp=int(rjp)
               dip=rip-ip
               djp=rjp-jp
             elseif(tversion == 8 .or. tversion == 81) then
               rip=(rlon-lonMin)/dlon+1
               rjp=(latMax-rlat)/dlat+1
               ip=int(rip)
               jp=int(rjp)
               dip=rip-ip
               djp=rjp-jp 
             elseif(tversion == 4 ) then
               rip=(rlon-lonMin)/dlon+1
               rjp=(rlat-latMin)/dlat+1
               ip=int(rip)
               jp=int(rjp)
               dip=rip-ip
               djp=rjp-jp
             elseif(tversion == 1 ) then
               if(rlon<0.0) rlon=360.0+rlon
               rip=(rlon-lonMin)/dlon+1
               rjp=(latMax-rlat)/dlat+1
               ip=int(rip)
               jp=int(rjp)
               dip=rip-ip
               djp=rjp-jp
             else
               write(*,*) ' Unknown Mosaic format !!'
               stop 123
             endif
             if( ip >= 1 .and. ip <= readref%mscNlon ) then
             if( jp >= 1 .and. jp <= readref%mscNlat ) then
! inside mosaic domain
               ipp1=min(ip+1,readref%mscNlon)
               jpp1=min(jp+1,readref%mscNlat)
               w1=(1.0-dip)*(1.0-djp)
               w2=dip*(1.0-djp)
               w3=dip*djp
               w4=(1.0-dip)*djp
               ref1=mscValue(ip,jp)
               ref2=mscValue(ipp1,jp)
               ref3=mscValue(ipp1,jpp1)
               ref4=mscValue(ip,jpp1)
               if(ref1 > readref%rthresh_ref .and. ref2 > readref%rthresh_ref .and.  &
                  ref3 > readref%rthresh_ref .and. ref4 > readref%rthresh_ref ) then
                  ref3d(i,kk)=(ref1*w1+ref2*w2+ref3*w3+ref4*w4)/float(readref%var_scale)
               elseif(ref1 > readref%rthresh_miss .and. ref2 > readref%rthresh_miss .and.  &
                  ref3 > readref%rthresh_miss .and. ref4 > readref%rthresh_miss ) then
                  ref3d(i,kk)=-99.0   ! clear
               else
                  ref3d(i,kk)=-999.0  ! no observation
               endif
             endif
             endif
          ENDDO
      ENDDO  ! mscNlev

      deallocate(mscValue)
  else
     ref3d=-999.0
!     write(*,*) trim(readref%mosaicfile), '   does not exist!!!'
  endif

end subroutine mosaic2grid
