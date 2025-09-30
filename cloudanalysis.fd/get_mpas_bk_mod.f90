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
  use constants, only: rad2deg
  use mpasio, only: read_MPAS_dim,read_MPAS_2D_real,read_MPAS_1D_real,read_MPAS_1D_int
  use mpasio, only: update_MPAS_2D_real,create_MPAS_2D_real

  implicit none
  private

  public :: t_bk,h_bk,p_bk,ps_bk,zh,q_bk,pblh
  public :: ges_ql,ges_qi,ges_qr,ges_qs,ges_qg,ges_qnr,ges_qni,ges_qnc,ges_qcf
  public :: nCell_full,nCell,nz,lon2,lat2,nsig
  public :: displ_1d,displ_2d,counts_send_1d,counts_send_2d
  public :: xlat,xlon,xland,soiltbk
!
  public :: read_mpas_init
  public :: read_mpas_bk
  public :: update_mpas
  public :: release_mem_mpas
  public :: write_new_2d_field_mpas
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
! information needed for cloudCover_NESDIS subroutine
  real(r_single),allocatable :: xlon(:,:)      ! Longitude
  real(r_single),allocatable :: xlat(:,:)      ! Latitude
  real(r_single),allocatable :: xland(:,:)     ! Surface type (water or land)
  real(r_single),allocatable :: soiltbk(:,:)   ! Background soil temperature
!
! MPAS mesh information
! Note that both the full background array and the smaller arrays on each processor
! have a vertical dimension of nz
  integer(i_kind) :: nCell_full,nz   ! Dimensions of full background array
  integer(i_kind) :: nCell           ! Dimensions of background arrays on each processor
  integer(i_kind) :: lon2,lat2,nsig  ! Dimensions used by cloudanalysis.fd
                                     ! For MPAS, lon2=1, lat2=nCell, nsig=nz
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

    allocate(displ_1d(npe))

    if (mype == 0) then

! Determine dimensions of full MPAS background
      dim_name = 'nCells'
      call read_MPAS_dim(mpasout_fname, dim_name, nCell_full)
      dim_name = 'nVertLevels'
      call read_MPAS_dim(mpasout_fname, dim_name, nz)

! Determine dimensions of background arrays on each processor
      allocate(counts_send_1d(npe))
      allocate(counts_send_2d(npe))
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
    call MPI_BCAST(displ_1d,npe,mpi_integer,0,mpi_comm_world,ierror)
    if (mype < remainder) then 
      nCell = min_cell + 1
    else
      nCell = min_cell
    endif

! Define lon2, lat2, and nsig (needed in cloudanalysis.fd)
    lon2 = 1
    lat2 = nCell
    nsig = nz

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
    write(6,'(A8,I5,A8,I10,A8,I5)') 'lon2 =', lon2, 'lat2 =', lat2, 'nsig =', nsig

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
    integer, allocatable :: tmp_int_full(:),tmp_int(:)
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
! Original nonvar cloud analysis used cloud fraction from MYNN
! Here, we use the generic cloud fraction that includes cloud fraction from MYNN
    allocate(ges_qcf(1,nCell,nz))
    varname = 'cldfrac'
    call read_scatter_2d_field(mype,mpasout_fname,varname,ges_qcf)
    do i=1,nz
      write(6,'(A12,I12,2E12.4)') 'cldfrac', i, maxval(ges_qcf(1,:,i)), minval(ges_qcf(1,:,i))
    enddo 

! Latitude (deg N)
    allocate(xlat(1,nCell))
    varname = 'latCell'
    call read_scatter_1d_field(mype,mpas_invariant_fname,varname,xlat,rem_dim=.false.)
    xlat = xlat * rad2deg
    write(6,'(A12,I12,2F12.4)') 'xlat', -1, maxval(xlat(1,:)), minval(xlat(1,:))

! Longitude (deg E, range 0 to 360)
    allocate(xlon(1,nCell))
    varname = 'lonCell'
    call read_scatter_1d_field(mype,mpas_invariant_fname,varname,xlon,rem_dim=.false.)
    xlon = xlon * rad2deg
    do i=1,nCell
      if (xlon(1,i) < 0) then
        xlon(1,i) = 360. + xlon(1,i)
      endif
    enddo
    write(6,'(A12,I12,2F12.4)') 'xlon', -1, maxval(xlon(1,:)), minval(xlon(1,:))

! Land/water mask (1 = land, 0 = water)
! Cannot use read_scatter_1d_field b/c landmask is an integer
    if (mype == 0) then
      allocate(tmp_int_full(nCell_full))
    endif
    allocate(tmp_int(nCell))
    allocate(xland(1,nCell))

    call MPI_BARRIER(mpi_comm_world,ierror)
    if (mype == 0) then
      varname = 'landmask'
      call read_MPAS_1D_int(mpas_invariant_fname, nCell_full, varname, tmp_int_full)
    endif
    call MPI_BARRIER(mpi_comm_world,ierror)
    call MPI_SCATTERV(tmp_int_full,counts_send_1d,displ_1d,MPI_INTEGER, &
                      tmp_int,nCell,MPI_INTEGER,0,mpi_comm_world,ierror)

    xland(1,:) = real(tmp_int(:))

    call MPI_BARRIER(mpi_comm_world,ierror)
    deallocate(tmp_int)
    if (mype == 0) then
      deallocate(tmp_int_full)
    endif

    write(6,'(A12,I12,2F12.4)') 'xland', -1, maxval(xland(1,:)), minval(xland(1,:))

! Skin temperature (K)
    allocate(soiltbk(1,nCell))
    varname = 'skintemp'
    call read_scatter_1d_field(mype,mpasout_fname,varname,soiltbk)
    write(6,'(A12,I12,2F12.4)') 'soiltbk', -1, maxval(soiltbk(1,:)), minval(soiltbk(1,:))

  end subroutine read_mpas_bk

  subroutine update_mpas(mype)
!---------------------------------------------------------------------------------------------------
!
! Update MPAS background using output from cloudanalysis.fd
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

    write(6,*)
    write(6,*) '------------------------------------------------------------'
    write(6,*) 'Updating MPAS background'
    write(6,*)

! Potential temperature
    varname = 'theta'
    call gather_write_2d_field(mype,mpasout_fname,varname,t_bk,write_minmax=.true.)

! Water vapor mass mixing ratio
    varname = 'qv'
    call gather_write_2d_field(mype,mpasout_fname,varname,q_bk,write_minmax=.true.)

! Cloud water mass mixing ratio
    varname = 'qc'
    call gather_write_2d_field(mype,mpasout_fname,varname,ges_ql)

! Cloud ice mass mixing ratio
    varname = 'qi'
    call gather_write_2d_field(mype,mpasout_fname,varname,ges_qi)

! Rain mass mixing ratio
    varname = 'qr'
    call gather_write_2d_field(mype,mpasout_fname,varname,ges_qr)

! Snow mass mixing ratio
    varname = 'qs'
    call gather_write_2d_field(mype,mpasout_fname,varname,ges_qs)

! Graupel mass mixing ratio
    varname = 'qg'
    call gather_write_2d_field(mype,mpasout_fname,varname,ges_qg)

! Rain number concentration
    varname = 'nr'
    call gather_write_2d_field(mype,mpasout_fname,varname,ges_qnr)

! Cloud ice number concentration
    varname = 'ni'
    call gather_write_2d_field(mype,mpasout_fname,varname,ges_qni)

! Cloud water number concentration
    varname = 'nc'
    call gather_write_2d_field(mype,mpasout_fname,varname,ges_qnc)

  end subroutine update_mpas

  subroutine write_new_2d_field_mpas(mype, field, varname)
!---------------------------------------------------------------------------------------------------
!
! Write a new 2D field to an existing MPAS netCDF file
!
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: mype
    real(r_single), intent(in) :: field(1,nCell,nz)
    character(len=50), intent(in) :: varname

    character(len=50) :: namedim1,namedim2

! Create field
    namedim1='nCells'
    namedim2='nVertLevels'
    if (mype == 0) then
      call create_MPAS_2D_real(mpasout_fname,namedim1,namedim2,varname)    
    endif

! Fill field
    call gather_write_2d_field(mype,mpasout_fname,varname,field)

  end subroutine write_new_2d_field_mpas

  subroutine release_mem_mpas
!---------------------------------------------------------------------------------------------------
!
! Deallocate MPAS background arrays used by cloudanalysis.fd
!
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

    implicit none

! Thermodynamic background fields
    if(allocated(t_bk)) deallocate(t_bk)
    if(allocated(h_bk)) deallocate(h_bk)
    if(allocated(p_bk)) deallocate(p_bk)
    if(allocated(ps_bk))deallocate(ps_bk)
    if(allocated(zh))   deallocate(zh)
    if(allocated(q_bk)) deallocate(q_bk)
    if(allocated(pblh)) deallocate(pblh)

! Hydrometeor fields
    if(allocated(ges_ql))  deallocate(ges_ql)
    if(allocated(ges_qi))  deallocate(ges_qi)
    if(allocated(ges_qr))  deallocate(ges_qr)
    if(allocated(ges_qs))  deallocate(ges_qs)
    if(allocated(ges_qg))  deallocate(ges_qg)
    if(allocated(ges_qnr)) deallocate(ges_qnr)
    if(allocated(ges_qni)) deallocate(ges_qni)
    if(allocated(ges_qnc)) deallocate(ges_qnc)
    if(allocated(ges_qcf)) deallocate(ges_qcf)

! Hydrometeor uncertainties
!    if(l_cld_uncertainty) then
!      if(allocated(unc_ql))  deallocate(unc_ql)
!      if(allocated(unc_qi))  deallocate(unc_qi)
!      if(allocated(unc_qr))  deallocate(unc_qr)
!      if(allocated(unc_qs))  deallocate(unc_qs)
!      if(allocated(unc_qg))  deallocate(unc_qg)
!    endif

! Extra fields for cloudCover_NESDIS subroutine
    if(allocated(xlon))     deallocate(xlon)
    if(allocated(xlat))     deallocate(xlat)
    if(allocated(xland))    deallocate(xland)
    if(allocated(soiltbk))  deallocate(soiltbk)

! Arrays used for scattering fields using MPI
    if(allocated(displ_1d))       deallocate(displ_1d)
    if(allocated(displ_2d))       deallocate(displ_2d)
    if(allocated(counts_send_1d)) deallocate(counts_send_1d)
    if(allocated(counts_send_2d)) deallocate(counts_send_2d)

  end subroutine release_mem_mpas

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

  subroutine gather_write_2d_field(mype, fname, varname, sub, add_dim, write_minmax)
!---------------------------------------------------------------------------------------------------
!
! Gather a 2D MPAS field from all processors, then write to a netCDF file
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
!   add_dim : boolean, optional
!     Option to explicitly create a third and final dimension with length 1 (usually time)
!   write_minmax : boolean, optional
!     Option to print min/max values from the root processor
!
! shawn.s.murdzek@noaa.gov
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: mype
    character(len=100), intent(in) :: fname
    character(len=50), intent(in) :: varname
    logical, intent(in), optional :: add_dim, write_minmax
    real(r_single), intent(in) :: sub(1,nCell,nz)

    real(r_single), allocatable :: tmp_full(:,:)
    real(r_single), allocatable :: tmp(:,:)
    integer :: i,ierror
    logical :: add_last_dim,do_write

    add_last_dim=.true.
    if (present(add_dim)) add_last_dim=add_dim

    do_write=.false.
    if (present(write_minmax)) do_write=write_minmax

! Allocate temporary arrays
    if (mype == 0) then
      allocate(tmp_full(nz,nCell_full))
    endif
    allocate(tmp(nz,nCell))

! Send data from all processors to the root processor
    do i=1,nz
      tmp(i,:) = sub(1,:,i)
    enddo
    call MPI_BARRIER(mpi_comm_world,ierror)
    call MPI_GATHERV(tmp,nCell*nz,MPI_REAL, &
                     tmp_full,counts_send_2d,displ_2d,MPI_REAL,0,mpi_comm_world,ierror)
    call MPI_BARRIER(mpi_comm_world,ierror)


! Write field from root processor to netCDF file
    if (mype == 0) then
      call update_MPAS_2D_real(fname, nz, nCell_full, varname, tmp_full, add_dim=add_last_dim)
      if (do_write) then
        write(6,'(4A12)') 'variable','lvl','max','min'
        do i=1,nz
          write(6,'(A12,I12,2E12.4)') trim(varname), i, maxval(tmp_full(i,:)), minval(tmp_full(i,:))
        enddo
      endif
    endif

! Clean up
    call MPI_BARRIER(mpi_comm_world,ierror)
    if (mype == 0) then
      deallocate(tmp_full)
    endif
    deallocate(tmp)

  end subroutine gather_write_2d_field

end module get_mpas_bk_mod
