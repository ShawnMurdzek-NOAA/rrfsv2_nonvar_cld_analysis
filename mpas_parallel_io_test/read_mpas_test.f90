program read_mpas_test
!
! Read multiple MPAS 3D fields with varying levels of parallelization
!

  use mpi
  use netcdf

  implicit none

  ! netCDF variables
  integer :: ncid,stat
  integer :: timeid,cellid,vertid,thid,rhid,cldid,rhoid,varid
  integer :: nTime,nCell,nVert
  character(len=50) :: timename,cellname,vertname
  real(8), allocatable :: theta(:,:),relhum(:,:),cldfrac(:,:),rho(:,:)
  real(8), allocatable :: var(:,:)

  ! MPI variables
  integer :: mype,npe,ierr

  ! Namelist variables
  integer :: use_mpi
  character(len=50) :: fname
  logical :: ifexist
  namelist/setup/use_mpi,fname

  ! Misc
  real(8) :: start_read, end_read


!===============================================================================

  ! Setup MPI
  call MPI_INIT(ierr) 
  call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)

  ! Read namelist
  if (mype == 0) then
    use_mpi = 0
    fname = 'mpas_atm.nc'
    inquire(file='namelist.input', EXIST=ifexist)
    if (ifexist) then
      open(15, file='namelist.input')
        read(15,setup)
      close(15)
    else
      write(*,*) 'WARNING: namelist.input missing. Using default values'
    endif
    write(*,*)
    write(*,'(a12,I5)') 'use_mpi = ', use_mpi
    write(*,'(a10,a50)') 'fname = ', fname
    write(*,*)

    write(*,'(a30,I5)') 'Number of processors =', npe
  endif

  call MPI_Bcast(use_mpi, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fname, 50, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  ! Open netCDF file
  if ((use_mpi == 0) .or. (use_mpi == -1)) then
    ! No parallelization
    if (mype == 0) then
      stat = nf90_open(fname, nf90_nowrite, ncid)
    endif
    call MPI_Bcast(stat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  elseif (use_mpi == 1) then
    ! Parallelization
    stat = nf90_open(fname, nf90_nowrite, ncid, comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
  elseif (use_mpi == 2) then
    ! Non-explicit parallelization
    stat = nf90_open(fname, nf90_nowrite, ncid)
  endif

  ! Check for errors
  if (stat /= nf90_noerr) then
    write(*,*) 'Error opening netCDF file. Stopping'
    write(*,*) stat
    stop
  endif

  ! Determine field dimensions and allocate fields
  if (mype == 0) then
    stat = nf90_inq_dimid(ncid, 'Time', timeid)
    stat = nf90_inquire_dimension(ncid, timeid, timename, nTime)
    stat = nf90_inq_dimid(ncid, 'nCells', cellid)
    stat = nf90_inquire_dimension(ncid, cellid, cellname, nCell)
    stat = nf90_inq_dimid(ncid, 'nVertLevels', vertid)
    stat = nf90_inquire_dimension(ncid, vertid, vertname, nVert)
    allocate(theta(nVert,nCell))
    allocate(relhum(nVert,nCell))
    allocate(cldfrac(nVert,nCell))
    allocate(rho(nVert,nCell))
    theta = 0
  endif

  call MPI_Bcast(nTime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nCell, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nVert, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  allocate(var(nCell,nVert))

  ! Read fields from MPAS netCDF
  if (use_mpi == -1) then

    ! Read 1 field in serial
    if (mype == 0) then
      start_read = MPI_Wtime()
      stat = nf90_inq_varid(ncid, 'theta', thid)
      stat = nf90_get_var(ncid, thid, theta, start=(/ 1, 1, 1 /), count=(/ nVert, nCell, 1 /))
      end_read = MPI_Wtime()
      write(*,*) 'Reading time for 1 field in serial (s) =', end_read - start_read
      stat = nf90_close(ncid)
    endif

  elseif (use_mpi == 0) then

    ! Read 4 fields in serial
    if (mype == 0) then
      start_read = MPI_Wtime()
      stat = nf90_inq_varid(ncid, 'theta', thid)
      stat = nf90_get_var(ncid, thid, theta, start=(/ 1, 1, 1 /), count=(/ nVert, nCell, 1 /))
      end_read = MPI_Wtime()
      write(*,*) 'Reading time for 1 field in serial (s) =', end_read - start_read
      stat = nf90_inq_varid(ncid, 'relhum', rhid)
      stat = nf90_get_var(ncid, rhid, relhum, start=(/ 1, 1, 1 /), count=(/ nVert, nCell, 1 /))
      stat = nf90_inq_varid(ncid, 'cldfrac', cldid)
      stat = nf90_get_var(ncid, cldid, cldfrac, start=(/ 1, 1, 1 /), count=(/ nVert, nCell, 1 /))
      stat = nf90_inq_varid(ncid, 'rho', rhoid)
      stat = nf90_get_var(ncid, rhoid, rho, start=(/ 1, 1, 1 /), count=(/ nVert, nCell, 1 /))
      end_read = MPI_Wtime()
      write(*,*) 'Reading time for 4 fields in serial (s) =', end_read - start_read
      stat = nf90_close(ncid)
    endif

  elseif (use_mpi > 0) then

    ! Use a separate processor to read each variable
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    start_read = MPI_Wtime()
    if (mype < 4) then
      if (mype == 0) then
        stat = nf90_inq_varid(ncid, 'theta', varid)
      elseif (mype == 1) then
        stat = nf90_inq_varid(ncid, 'relhum', varid)
      elseif (mype == 2) then
        stat = nf90_inq_varid(ncid, 'cldfrac', varid)
      elseif (mype == 3) then
        stat = nf90_inq_varid(ncid, 'rho', varid)
      endif
      stat = nf90_get_var(ncid, varid, var, start=(/ 1, 1, 1 /), count=(/ nVert, nCell, 1 /))
      end_read = MPI_Wtime()
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end_read = MPI_Wtime()
    if (mype == 0) then
      write(*,*) 'Reading time for 4 fields in parallel (s) =', end_read - start_read
      theta = var
    endif
    stat = nf90_close(ncid)

    ! Send data to root processor
    if (mype == 0) then
      theta = var
      call MPI_RECV(var,nCell*nVert,MPI_DOUBLE_PRECISION,1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      relhum = var
      call MPI_RECV(var,nCell*nVert,MPI_DOUBLE_PRECISION,2,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      cldfrac = var
      call MPI_RECV(var,nCell*nVert,MPI_DOUBLE_PRECISION,3,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      rho = var
    elseif (mype == 1) then
      call MPI_SEND(var,nCell*nVert,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
    elseif (mype == 2) then
      call MPI_SEND(var,nCell*nVert,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
    elseif (mype == 3) then
      call MPI_SEND(var,nCell*nVert,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
    endif  

  endif

  ! Compute the total of theta, relhum, cldfrac, and rho
  if (mype == 0) then
    write(*,*)
    write(*,'(a20,E20.8)') 'Sum theta =', sum(theta)
    write(*,'(a20,E20.8)') 'Sum relhum =', sum(relhum)
    write(*,'(a20,E20.8)') 'Sum cldfrac =', sum(cldfrac)
    write(*,'(a20,E20.8)') 'Sum rho =', sum(rho)
  endif

  call MPI_FINALIZE(ierr)

end program read_mpas_test 
