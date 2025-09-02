module get_mpas_bk_mod
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! program:  get_mpas_bk_mod driver to get MPAS background
! calculation
!
!   PRGMMR: Ming Hu          ORG: GSD/AMB        DATE: 2020-8-10
!
! ABSTRACT: 
!
! PROGRAM HISTORY LOG:
!    2020-08-10  Hu  Add NCO document block
!
! USAGE:
!   INPUT FILES: 
!
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90 
!   MACHINE:  Linux cluster (JET)
!
!$$$
!
!_____________________________________________________________________
!

  use mpi
  use kinds,   only: r_single,i_kind, r_kind
  use constants, only: init_constants,init_constants_derived
  use constants, only: rd,h1000,rd_over_cp,grav
  use mpasio, only: read_MPAS_dim,read_MPAS_2D_real,read_MPAS_1D_real

  implicit none
  private

  public :: t_bk,h_bk,p_bk,ps_bk,zh,q_bk,pblh
  public :: ges_ql,ges_qi,ges_qr,ges_qs,ges_qg,ges_qnr,ges_qni,ges_qnc,ges_qcf
  public :: nCell_full,nCell,nz
  public :: displ_1d,displ_2d,counts_send_1d,counts_send_2d
!
  public :: read_mpas_init
  public :: read_mpas_bk
!
!
! background
!
  real(r_single),allocatable:: t_bk(:,:,:)     ! potential temperature
  real(r_single),allocatable:: h_bk(:,:,:)     ! height above ground
  real(r_single),allocatable:: p_bk(:,:,:)     ! pressure
  real(r_single),allocatable:: ps_bk(:,:)      ! surface pressure
  real(r_single),allocatable:: zh(:,:)         ! terrain height
  real(r_single),allocatable:: q_bk(:,:,:)     ! water vapor mass mixing ratio

  real(r_single),allocatable:: pblh(:,:)       ! PBL height (grid coordinate)
!
! hydrometeors
  real(r_single),allocatable :: ges_ql(:,:,:)  ! cloud water
  real(r_single),allocatable :: ges_qi(:,:,:)  ! could ice
  real(r_single),allocatable :: ges_qr(:,:,:)  ! rain
  real(r_single),allocatable :: ges_qs(:,:,:)  ! snow
  real(r_single),allocatable :: ges_qg(:,:,:)  ! graupel
  real(r_single),allocatable :: ges_qnr(:,:,:) ! rain number concentration
  real(r_single),allocatable :: ges_qni(:,:,:) ! cloud ice number concentration
  real(r_single),allocatable :: ges_qnc(:,:,:) ! cloud water number concentration
  real(r_single),allocatable :: ges_qcf(:,:,:) ! cloud fraction
!
! MPAS mesh information
! Note that both the full background array and the smaller arrays on each processor
! have a vertical dimension of nz
  integer(i_kind) :: nCell_full,nz   ! Dimensions of full background array
  integer(i_kind) :: nCell           ! Dimensions of background arrays on each processor
  character(len=50) :: dim_name
!
! background files
  character(len=100) :: mpasout_fname          ! MPAS netCDF output file
  character(len=100) :: mpas_invariant_fname   ! MPAS invariant file (contains heights MSL)
!
! MPI variables for parallelization
  integer(i_kind), allocatable :: counts_send_1d(:),displ_1d(:)
  integer(i_kind), allocatable :: counts_send_2d(:),displ_2d(:)
!

contains
!

  subroutine read_mpas_init(mype, npe)
!---------------------------------------------------------------------------------------------------
!
! Read in MPAS mesh information on root processor, then determine size of background arrays on each
! processor
!
! Inputs
!   mype : integer
!     Processor number
!   npe : integer
!     Total number of processors
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: mype,npe

    integer :: ierror
    integer :: min_cell,remainder
    integer :: i

    mpasout_fname = 'mpasout.nc'
    mpas_invariant_fname = 'invariant.nc'

    if (mype == 0) then

! Determine dimensions of full MPAS background
      dim_name = 'nCells'
      call read_MPAS_dim(mpasout_fname, dim_name, nCell_full)
      dim_name = 'nVertLevels'
      call read_MPAS_dim(mpasout_fname, dim_name, nz)

! Determine dimensions of background arrays on each processor
      allocate(counts_send_1d(npe))
      allocate(counts_send_2d(npe))
      allocate(displ_1d(npe))
      allocate(displ_2d(npe))
      min_cell = nCell_full / npe
      remainder = mod(nCell_full, npe)
      do i=1,npe
        if ((i-1) < remainder) then
          counts_send_1d(i) = min_cell + 1
          counts_send_2d(i) = (min_cell + 1) * nz
        else
          counts_send_1d(i) = min_cell
          counts_send_2d(i) = min_cell * nz
        endif
        if (i > 1) then
          displ_1d(i) = displ_1d(i-1) + counts_send_1d(i-1)
          displ_2d(i) = displ_2d(i-1) + counts_send_2d(i-1)
        else
          displ_1d(i) = 0
          displ_2d(i) = 0
        endif
      enddo 

    endif

! Determine size of the background arrays on each processor
    call MPI_BCAST(nCell_full,1,mpi_integer,0,mpi_comm_world,ierror)
    call MPI_BCAST(min_cell,1,mpi_integer,0,mpi_comm_world,ierror)
    call MPI_BCAST(remainder,1,mpi_integer,0,mpi_comm_world,ierror)
    call MPI_BCAST(nz,1,mpi_integer,0,mpi_comm_world,ierror)
    if (mype < remainder) then 
      nCell = min_cell + 1
    else
      nCell = min_cell
    endif

! Write info (useful for debugging)
    write(6,*)
    write(6,*) '------------------------------------------------------------'
    write(6,*) 'Initializing information for reading MPAS background'
    if (mype == 0) then
      write(6,*)
      write(6,*) 'Dimensions of full MPAS background array:'
      write(6,'(A8,I10,A8,I5)') 'nCell =', nCell_full, 'nz =', nz
      write(6,*)
      write(6,*) 'Dimensions of MPAS background on each processor:'
      write(6,'(A5,4A18)') 'mype','displ_1d','displ_2d','counts_send_1d','counts_send_2d'
      do i=1,npe
        write(6,'(I5,5I18)') i-1,displ_1d(i),displ_2d(i),counts_send_1d(i),counts_send_2d(i)
      enddo
    endif
    write(6,*)
    write(6,*) 'Dimensions of MPAS background array on this processor:'
    write(6,'(A8,I10,A8,I5)') 'nCell =', nCell, 'nz =', nz

  end subroutine read_mpas_init

  subroutine read_mpas_bk(mype)
!---------------------------------------------------------------------------------------------------
!
! Read in MPAS background fields, scatter to all processors, and convert to proper units
!
! Parallel I/O is not possible with MPAS netCDF files, so the strategy here is to read in each field
! using the root process and immediately scatter the field to all processors before reading the
! next field. This makes it so the root processor will have at most 1 full 3D MPAS field stored in
! memory at a time (except when reading height AGL). After running this subroutine, each processor 
! will have possess all necessary background fields for a subset of the MPAS horizontal domain. 
! read_mpas_init must be run prior to running this subroutine.
!
! Inputs
!   mype : integer
!     Processor number
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: mype

    character(len=50) :: varname
    integer :: i,ierror
    real(r_single), allocatable :: tmp1_full(:,:),tmp2_full(:,:),tmp(:,:)
    real(r_single) :: p_base(1,nCell,nz)

    write(6,*)
    write(6,*) '------------------------------------------------------------'
    write(6,*) 'Reading MPAS background'
    write(6,*)
    write(6,'(4A12)') 'variable','lvl','max','min'

! Potential temperature (K)
    allocate(t_bk(1,nCell,nz))
    varname = 'theta'
    call read_scatter_2d_field(mype,mpasout_fname,varname,t_bk)
    do i=1,nz
      write(6,'(A12,I12,2F12.4)') 'theta', i, maxval(t_bk(1,:,i)), minval(t_bk(1,:,i))
    enddo 

! Surface pressure (hPa)
    allocate(ps_bk(1,nCell))
    varname = 'surface_pressure'
    call read_scatter_1d_field(mype,mpasout_fname,varname,ps_bk)
    ps_bk = ps_bk * 0.01
    write(6,'(A12,I12,2F12.4)') 'ps', -1, maxval(ps_bk(1,:)), minval(ps_bk(1,:))

! Terrain height (m)
    allocate(zh(1,nCell))
    varname = 'ter'
    call read_scatter_1d_field(mype,mpas_invariant_fname,varname,zh,rem_dim=.false.)
    write(6,'(A12,I12,2F12.4)') 'zh', -1, maxval(zh(1,:)), minval(zh(1,:))

! Height above ground (m AGL)
! Note that zgrid is the MSL height of the interface between two layers. We want the AGL height of 
! the center of the layer
    if (mype == 0) then
      allocate(tmp1_full(nz+1, nCell_full))
      allocate(tmp2_full(nz, nCell_full))
    endif
    allocate(tmp(nz, nCell))
    allocate(h_bk(1,nCell,nz))

    call MPI_BARRIER(mpi_comm_world,ierror)
    if (mype == 0) then
      varname = 'zgrid'
      call read_MPAS_2D_real(mpas_invariant_fname, nz+1, nCell_full, varname, tmp1_full, rem_dim=.true.)
      do i=1,nz
        tmp2_full(i,:) = 0.5*(tmp1_full(i,:) + tmp1_full(i+1,:))
      enddo
      deallocate(tmp1_full)
    endif
    call MPI_BARRIER(mpi_comm_world,ierror)
    call MPI_SCATTERV(tmp2_full,counts_send_2d,displ_2d,MPI_REAL, &
                      tmp,nCell*nz,MPI_REAL,0,mpi_comm_world,ierror)

    do i=1,nz
      h_bk(1,:,i) = tmp(i,:) - zh(1,:)
    enddo

    call MPI_BARRIER(mpi_comm_world,ierror)
    deallocate(tmp)
    if (mype == 0) then
      deallocate(tmp2_full)
    endif

    do i=1,nz
      write(6,'(A12,I12,2E12.4)') 'hgt', i, maxval(h_bk(1,:,i)), minval(h_bk(1,:,i))
    enddo 

! Pressure (hPa)
    allocate(p_bk(1,nCell,nz))
    varname = 'pressure_p'
    call read_scatter_2d_field(mype,mpasout_fname,varname,p_bk)
    varname = 'pressure_base'
    call read_scatter_2d_field(mype,mpasout_fname,varname,p_base)
    p_bk = (p_bk + p_base) * 0.01
    do i=1,nz
      write(6,'(A12,I12,2E12.4)') 'p', i, maxval(p_bk(1,:,i)), minval(p_bk(1,:,i))
    enddo 

! Water vapor mass mixing ratio (kg/kg)
    allocate(q_bk(1,nCell,nz))
    varname = 'qv'
    call read_scatter_2d_field(mype,mpasout_fname,varname,q_bk)
    do i=1,nz
      write(6,'(A12,I12,2E12.4)') 'qv', i, maxval(q_bk(1,:,i)), minval(q_bk(1,:,i))
    enddo 

! PBL height (level number)
    allocate(pblh(1,nCell))
    call calc_pbl_height(nCell,1,nz,q_bk,t_bk,p_bk,pblh)
    write(6,'(A12,I12,2F12.4)') 'pblh', -1, maxval(pblh(1,:)), minval(pblh(1,:))

! Cloud water mass mixing ratio (kg/kg)
    allocate(ges_ql(1,nCell,nz))
    varname = 'qc'
    call read_scatter_2d_field(mype,mpasout_fname,varname,ges_ql)
    do i=1,nz
      write(6,'(A12,I12,2E12.4)') 'ql', i, maxval(ges_ql(1,:,i)), minval(ges_ql(1,:,i))
    enddo 

! Cloud ice mass mixing ratio (kg/kg)
    allocate(ges_qi(1,nCell,nz))
    varname = 'qi'
    call read_scatter_2d_field(mype,mpasout_fname,varname,ges_qi)
    do i=1,nz
      write(6,'(A12,I12,2E12.4)') 'qi', i, maxval(ges_qi(1,:,i)), minval(ges_qi(1,:,i))
    enddo 

! Rain mass mixing ratio (kg/kg)
    allocate(ges_qr(1,nCell,nz))
    varname = 'qr'
    call read_scatter_2d_field(mype,mpasout_fname,varname,ges_qr)

! Snow mass mixing ratio (kg/kg)
    allocate(ges_qs(1,nCell,nz))
    varname = 'qs'
    call read_scatter_2d_field(mype,mpasout_fname,varname,ges_qs)

! Graupel mass mixing ratio (kg/kg)
    allocate(ges_qg(1,nCell,nz))
    varname = 'qg'
    call read_scatter_2d_field(mype,mpasout_fname,varname,ges_qg)

! Rain number concentration (number/kg)
    allocate(ges_qnr(1,nCell,nz))
    varname = 'nr'
    call read_scatter_2d_field(mype,mpasout_fname,varname,ges_qnr)

! Cloud ice number concentration (number/kg)
    allocate(ges_qni(1,nCell,nz))
    varname = 'ni'
    call read_scatter_2d_field(mype,mpasout_fname,varname,ges_qni)

! Cloud water number concentration (number/kg)
    allocate(ges_qnc(1,nCell,nz))
    varname = 'nc'
    call read_scatter_2d_field(mype,mpasout_fname,varname,ges_qnc)

! Cloud fraction (unitless)
    allocate(ges_qcf(1,nCell,nz))
    varname = 'cldfrac_bl'
    call read_scatter_2d_field(mype,mpasout_fname,varname,ges_qcf)
    do i=1,nz
      write(6,'(A12,I12,2E12.4)') 'cldfrac', i, maxval(ges_qcf(1,:,i)), minval(ges_qcf(1,:,i))
    enddo 

  end subroutine read_mpas_bk

  subroutine read_scatter_1d_field(mype, fname, varname, sub, rem_dim)
!---------------------------------------------------------------------------------------------------
!
! Read in a single 1D MPAS background field and scatter to all processors
!
! Inputs
!   mype : integer
!     Processor number
!   fname : string
!     netCDF file to read
!   varname : string
!     MPAS variable name to read
!   sub : 2D array
!     Subset of the full MPAS 1D field on each processor. First dimension should always be 1
!   rem_dim : boolean, optional
!     Option to explicitly read only the first element in the last dimension (usually time)
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: mype
    character(len=100), intent(in) :: fname
    character(len=50), intent(in) :: varname
    logical, intent(in), optional :: rem_dim
    real(r_single), intent(out) :: sub(1,nCell)

    real(r_single), allocatable :: tmp_full(:)
    real(r_single), allocatable :: tmp(:)
    logical :: remove_dim
    integer :: ierror

    remove_dim=.true.
    if (present(rem_dim)) remove_dim=rem_dim

! Allocate temporary arrays
    if (mype == 0) then
      allocate(tmp_full(nCell_full))
    endif
    allocate(tmp(nCell))

! Read full 1D field, then scatter to all processors
    call MPI_BARRIER(mpi_comm_world,ierror)
    if (mype == 0) then
      call read_MPAS_1D_real(fname, nCell_full, varname, tmp_full, rem_dim=remove_dim)
    endif
    call MPI_BARRIER(mpi_comm_world,ierror)
    call MPI_SCATTERV(tmp_full,counts_send_1d,displ_1d,MPI_REAL, &
                      tmp,nCell,MPI_REAL,0,mpi_comm_world,ierror)

! Reshape array  
    sub(1,:) = tmp(:)

! Clean up
    call MPI_BARRIER(mpi_comm_world,ierror)
    deallocate(tmp)
    if (mype == 0) then
      deallocate(tmp_full)
    endif

  end subroutine read_scatter_1d_field

  subroutine read_scatter_2d_field(mype, fname, varname, sub)
!---------------------------------------------------------------------------------------------------
!
! Read in a single 2D MPAS background field and scatter to all processors
!
! Inputs
!   mype : integer
!     Processor number
!   fname : string
!     netCDF file to read
!   varname : string
!     MPAS variable name to read
!   sub : 3D array
!     Subset of the full MPAS 2D field on each processor. First dimension should always be 1
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: mype
    character(len=100), intent(in) :: fname
    character(len=50), intent(in) :: varname
    real(r_single), intent(out) :: sub(1,nCell,nz)

    real(r_single), allocatable :: tmp_full(:,:)
    real(r_single), allocatable :: tmp(:,:)

    integer :: i,ierror

! Allocate temporary arrays
    if (mype == 0) then
      allocate(tmp_full(nz, nCell_full))
    endif
    allocate(tmp(nz, nCell))

! Read full 2D field, then scatter to all processors
    call MPI_BARRIER(mpi_comm_world,ierror)
    if (mype == 0) then
      call read_MPAS_2D_real(fname, nz, nCell_full, varname, tmp_full, rem_dim=.true.)
    endif
    call MPI_BARRIER(mpi_comm_world,ierror)
    call MPI_SCATTERV(tmp_full,counts_send_2d,displ_2d,MPI_REAL, &
                      tmp,nCell*nz,MPI_REAL,0,mpi_comm_world,ierror)

! Reshape array
    do i=1,nz
      sub(1,:,i) = tmp(i,:)
    enddo

! Clean up
    call MPI_BARRIER(mpi_comm_world,ierror)
    deallocate(tmp)
    if (mype == 0) then
      deallocate(tmp_full)
    endif

  end subroutine read_scatter_2d_field

end module get_mpas_bk_mod
