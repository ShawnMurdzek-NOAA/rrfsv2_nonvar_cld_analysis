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
  character(len=100) :: mpasout_fname   ! MPAS netCDF output file
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

  end subroutine read_mpas_init

  subroutine read_mpas_bk(mype)
!---------------------------------------------------------------------------------------------------
!
! Read in MPAS background fields, scatter to all processors, and convert to proper units
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
    integer :: ierror
    integer :: i
    real(r_single), allocatable :: tmp2d_full(:,:),tmp1d_full(:)
    real(r_single), allocatable :: tmp2d(:,:),tmp1d(:)

    if (mype == 0) then
      allocate(tmp2d_full(nz, nCell_full))
      allocate(tmp1d_full(nCell_full))
    endif
    allocate(tmp2d(nz, nCell))
    allocate(tmp1d(nCell))
    
    allocate(t_bk(1,nCell,nz))
    allocate(ps_bk(1,nCell))

! Potential temperature (K)
    call MPI_BARRIER(mpi_comm_world,ierror)
    if (mype == 0) then
      varname = 'theta'
      call read_MPAS_2D_real(mpasout_fname, nz, nCell_full, varname, tmp2d_full, rem_dim=.true.)
      write(*,*) sum(tmp2d_full)
    endif
    call MPI_BARRIER(mpi_comm_world,ierror)
    call MPI_SCATTERV(tmp2d_full,counts_send_2d,displ_2d,MPI_REAL, &
                      tmp2d,nCell*nz,MPI_REAL,0,mpi_comm_world,ierror)
    do i=1,nz
      t_bk(1,:,i) = tmp2d(i,:)
    enddo

! Surface pressure (hPa)
    call MPI_BARRIER(mpi_comm_world,ierror)
    if (mype == 0) then
      varname = 'surface_pressure'
      call read_MPAS_1D_real(mpasout_fname, nCell_full, varname, tmp1d_full, rem_dim=.true.)
    endif
    call MPI_BARRIER(mpi_comm_world,ierror)
    call MPI_SCATTERV(tmp1d_full,counts_send_1d,displ_1d,MPI_REAL, &
                      tmp1d,nCell,MPI_REAL,0,mpi_comm_world,ierror)
    ps_bk(1,:) = tmp1d * 0.01

! Clean up
    call MPI_BARRIER(mpi_comm_world,ierror)
    deallocate(tmp2d)
    deallocate(tmp1d)
    if (mype == 0) then
      deallocate(tmp2d_full)
      deallocate(tmp1d_full)
    endif

  end subroutine read_mpas_bk

end module get_mpas_bk_mod
