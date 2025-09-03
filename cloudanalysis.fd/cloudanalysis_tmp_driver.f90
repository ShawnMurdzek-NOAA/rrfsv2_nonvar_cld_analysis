program cloudanalysis
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  cloudanalysis      driver for generalized No Varational cloud/hydrometeor analysis
!
!   PRGMMR: Ming Hu          ORG: GSL/AVID        DATE: 2020-08-10
!
! ABSTRACT: 
!  This subroutine serves as a driver for generalized No Varational cloud/hydrometeor analysis
!
! PROGRAM HISTORY LOG:
!    2008-12-20  Hu  Add NCO document block
!
!   input argument list:
!     mype     - processor ID that does this IO
!
!   output argument list:
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
!   MACHINE:  Linux cluster (WJET) at NOAA/ESRL - Boulder, CO
!
!$$$
!
!_____________________________________________________________________
!
! 
  use mpi
  use kinds,   only: r_single,i_kind, r_kind
!  use wrf_mass_guess_mod, only: soil_temp_cld,isli_cld,ges_xlon,ges_xlat,ges_tten

  use constants, only: init_constants,init_constants_derived
  use constants, only: rd_over_cp, h1000
  use constants, only: zero,one,rad2deg,fv

  use rapidrefresh_cldsurf_mod, only: init_rapidrefresh_cldsurf
  use rapidrefresh_cldsurf_mod, only: dfi_radar_latent_heat_time_period,   &
                                      metar_impact_radius,                 &
                                      l_cleanSnow_WarmTs,l_conserve_thetaV,&
                                      r_cleanSnow_WarmTs_threshold,        &
                                      i_conserve_thetaV_iternum,           &
                                      l_cld_bld, cld_bld_hgt,              &
                                      build_cloud_frac_p, clear_cloud_frac_p, &
                                      nesdis_npts_rad, &
                                      iclean_hydro_withRef, iclean_hydro_withRef_allcol, &
                                      l_use_hydroretrieval_all, &
                                      i_lightpcp, l_numconc, qv_max_inc,ioption, &
                                      l_precip_clear_only,l_fog_off,cld_bld_coverage,cld_clr_coverage,&
                                      i_T_Q_adjust,l_saturate_bkCloud,i_precip_vertical_check,l_rtma3d, &
                                      l_qnr_from_qr, n0_rain, &
                                      r_cloudfrac_threshold,l_cld_uncertainty

  use namelist_mod, only: load_namelist
  use namelist_mod, only: iyear,imonth,iday,ihour,iminute,isecond

  use get_mpas_bk_mod, only: nCell_full,nCell,nz,lon2,lat2,nsig
  use get_mpas_bk_mod, only: displ_1d,displ_2d,counts_send_1d,counts_send_2d
  use get_mpas_bk_mod, only: t_bk,h_bk,p_bk,ps_bk,zh,q_bk,pblh
  use get_mpas_bk_mod, only: ges_ql,ges_qi,ges_qr,ges_qs,ges_qg,ges_qnr,ges_qni,ges_qnc,ges_qcf
  use get_mpas_bk_mod, only: read_mpas_init,read_mpas_bk
!
  implicit none

! MPI variables
  integer :: npe, mype, ierror

! Declare passed variables
!
  integer :: regional_time(6)
! background
!
!  real(r_single),allocatable:: xlon(:,:)        ! 2D longitude in each grid
!  real(r_single),allocatable:: xlat(:,:)        ! 2D latitude in each grid
!  real(r_single),  allocatable:: xland(:,:)
!  real(r_single),allocatable:: soiltbk(:,:)
!
!  surface observation
!
  character(len=7) :: obstype
  character(len=20):: isis
  integer :: nreal,nchanl,ilat,ilon,ndata_model
  integer :: istart,jstart
!
  integer(i_kind) :: nvarcld_p
  parameter (nvarcld_p=13)

  integer(i_kind)              :: numsao
  real(r_single), allocatable  :: oi(:)
  real(r_single), allocatable  :: oj(:)
  integer(i_kind),allocatable  :: ocld(:,:)
  character*10,   allocatable  :: owx(:)
  real(r_single), allocatable  :: oelvtn(:)
  real(r_single), allocatable  :: odist(:)
  character(8),   allocatable  :: cstation(:)
  real(r_single), allocatable  :: oistation(:)
  real(r_single), allocatable  :: ojstation(:)
  real(r_single), allocatable  :: wimaxstation(:)
!
  integer(i_kind),allocatable  :: osfc_station_map(:,:)
!
!  lightning observation: 2D field in RR grid
!
  real(r_single),allocatable  :: lightning(:,:)
!
!  GOES - NASA LaRc cloud products: several 2D fields in RR grid
!
  real(r_single),allocatable  :: nasalarc_cld(:,:,:)

!
!  radar observation : 3D reflectvity in RR grid
!
  real(r_single),allocatable :: ref_mos_3d(:,:,:)
  real(r_single),allocatable :: ref_mos_3d_tten(:,:,:)
  real(r_single),allocatable :: ref_mosaic31(:,:,:)
  integer(i_kind)          :: nmsclvl_radar
!
!  GOES - NESDIS cloud products : 2d fields
!
  real(r_single), allocatable :: sat_ctp(:,:)
  real(r_single), allocatable :: sat_tem(:,:)
  real(r_single), allocatable :: w_frac(:,:)
  integer(i_kind),allocatable :: nlev_cld(:,:)
!
! cloud/hydrometeor analysis variables
!
!=========================================================================
!  cld_cover_3d in the Generalized Cloud/Hydrometeor Analysis code
!   Definition:  3-d gridded observation-based information
!      including 0-1 cloud-fraction (with zero value for clear)
!      and negative values indicating "unknown" status.
!   cld_cover_3d is initialized with negative values.
!   cld_type_3d, pcp_type_3d, wthr_type_2d - similar to cld_cover_3d
!=========================================================================

  real(r_single), allocatable :: cld_cover_3d(:,:,:)  ! cloud cover
  integer(i_kind),allocatable :: cld_type_3d(:,:,:)   ! cloud type
  integer(i_kind),allocatable :: pcp_type_3d(:,:,:)   ! precipitation type
  integer(i_kind),allocatable :: wthr_type_2d(:,:)    ! weather type type
  integer(i_kind),allocatable :: cloudlayers_i(:,:,:) ! 5 different layers 
                                                      ! 1= the number of layers
                                                      ! 2,4,... bottom
                                                      ! 3,5,... top
!
  real(r_single),allocatable :: cldwater_3d(:,:,:)    ! cloud water
  real(r_single),allocatable :: nwater_3d(:,:,:)      ! cloud water number concentration
  real(r_single),allocatable :: cldice_3d(:,:,:)      ! cloud ice
  real(r_single),allocatable :: nice_3d(:,:,:)        ! cloud ice number concentration
  real(r_single),allocatable :: rain_3d(:,:,:)        ! rain
  real(r_single),allocatable :: nrain_3d(:,:,:)       ! rain number concentration
  real(r_single),allocatable :: snow_3d(:,:,:)        ! snow
  real(r_single),allocatable :: graupel_3d(:,:,:)     ! graupel
  real(r_single),allocatable :: cldtmp_3d(:,:,:)      ! cloud temperature

  real(r_single),allocatable :: rain_1d_save(:)       ! rain
  real(r_single),allocatable :: nrain_1d_save(:)      ! rain number concentration    
  real(r_single),allocatable :: snow_1d_save(:)       ! snow
  real(r_single),allocatable :: vis2qc(:,:)           ! fog

  real(r_kind)    ::  thunderRadius=2.5_r_kind
  integer(i_kind) :: miss_obs_int
  real(r_kind)    :: miss_obs_real
  parameter ( miss_obs_int = -99999  )
  parameter ( miss_obs_real = -99999.0_r_kind )
  real(r_single)  ::  krad_bot          ! radar bottom level

!
! collect cloud
  real(r_kind)    :: cloud_def_p
  data  cloud_def_p       / 0.000001_r_kind/
  real(r_kind),allocatable :: sumqci(:,:,:)  ! total liquid water
  real(r_kind),allocatable :: watericemax(:,:)  ! max of total liquid water
  integer(i_kind),allocatable :: kwatericemax(:,:)  ! lowest level of total liquid water
  real(r_single),allocatable::temp1(:,:),tempa(:)
  real(r_single),allocatable::all_loc(:,:)
  real(r_single),allocatable::strp(:)
  integer(i_kind) :: im,jm
!
! option in namelist
!
  integer(i_kind) :: opt_cloudwaterice_retri  ! method for cloud water retrieval
  integer(i_kind) :: opt_hydrometeor_retri    ! method for precipitation retrieval
  integer(i_kind) :: opt_cloudtemperature     ! if open temperature adjustment scheme
  integer(i_kind) :: istat_surface,istat_nesdis,istat_radar    ! 1 has observation
  integer(i_kind) :: istat_nasalarc,istat_lightning            ! 0 no observation
  integer(i_kind) :: imerge_nesdis_nasalarc  !  =1 merge NASA LaRC with NESDIS
                                             !  =2 use NASA LaRC only
                                             !  = other, use NESDIS only
!
!
  real(r_kind), pointer :: ges_z (:,:  )=>NULL()  ! geopotential height
  real(r_kind), pointer :: ges_ps(:,:  )=>NULL()  ! surface pressure
  real(r_single), pointer :: ges_tv(:,:,:)=>NULL()  ! virtual temperature
  real(r_single), pointer :: ges_q (:,:,:)=>NULL()  ! specifici humidity
!
!  misc.
!
  integer(i_kind) :: ytotal,ybegin,yend
  integer(i_kind) :: i,j,k,kk
  integer(i_kind) :: iglobal,jglobal,ilocal,jlocal
  logical :: ifindomain
  integer(i_kind) :: imaxlvl_ref
  real(r_kind)    :: max_retrieved_qrqs,max_bk_qrqs,ratio_hyd_bk2obs
  real(r_kind)    :: qrqs_retrieved
  real(r_kind)    :: qrlimit,qrlimit_lightpcp
  real(r_kind)    :: qnr_limit
  real(r_kind)    :: dbz_clean_graupel
  integer(i_kind) :: ilat1s,ilon1s
  integer(i_kind) :: clean_count,build_count,part_count,miss_count
  integer :: sss,rrr

  real(r_kind)    :: refmax,snowtemp,raintemp,nraintemp,graupeltemp
  real(r_kind)    :: snowadd,ratio2
  integer(i_kind) :: imax, jmax, ista, iob, job
  real(r_kind)    :: dfi_lhtp, qmixr, tsfc
  real(r_kind)    :: Temp, watwgt
  real(r_kind)    :: cloudwater, cloudice

  real(r_kind),parameter    :: pi = 4._r_kind*atan(1._r_kind)
  real(r_kind),parameter    :: rho_w = 999.97_r_kind, rho_a = 1.2_r_kind
  real(r_kind),parameter    :: cldDiameter = 10.0E3_r_kind

  real(r_kind),parameter :: am_r = pi * 1000.0_r_kind / 6.0_r_kind
  real(r_kind)           :: lambda

! local variables used for adjustment of qr/qs for RTMA_3D to alleviate ghost reflectivity
  logical         :: print_verbose
  logical         :: verbose
  integer(i_kind) :: k_cap            ! highest level when adjument is done (used for adjust qr/qs for RTMA_3D)
  logical         :: fileexist
  character(len=80) :: obsfile
  integer         :: lunin
  integer         :: nsat1
  integer         :: istatus
!
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  write(obsfile,'(a,I4.4)') 'stdout_cloudanalysis.d',mype
  open(6, file=trim(obsfile),form='formatted',status='unknown')
  write(6,*) '===> cloud analysis over subdomain = ', mype
!
  clean_count=0
  build_count=0
  part_count=0
  miss_count=0
!!
  write(6,*) '========================================'
  write(6,*) 'gsdcloudanalysis: Start generalized cloud analysis', mype
  write(6,*) '========================================'
!
!
  call init_constants(.true.)
  call init_constants_derived

  call init_rapidrefresh_cldsurf
  call load_namelist(mype)
!
! get background ready
! Note that MPAS still uses 3D arrays even though the data is 2D, and the arrays
! still have dimensions (lon2, lat2, nsig)
! For MPAS, lon2 = 1, lat2 = nCell, nsig = nz
!
  call read_mpas_init(mype,npe)
  call MPI_BARRIER(mpi_comm_world,ierror)
  call read_mpas_bk(mype)
!
!
  krad_bot=7.0_r_single

  opt_hydrometeor_retri=3       ! 1=Kessler 2=Lin 3=Thompson
  opt_cloudtemperature=3        ! 3=latent heat, 4,5,6 = adiabat profile
  opt_cloudwaterice_retri=1     ! 1 = RUC layer saturation and autoconvert
                                ! 2 = convective 
  imerge_nesdis_nasalarc=2      !  =1 merge NASA LaRC with NESDIS
                                !  =2 use NASA LaRC only
                                !  =3 No Satellite cloud top used
                                !  = other, use NESDIS only

!
! initialize the observation flag  
!
  istat_surface=0
  istat_nesdis=0
  istat_radar=0
  istat_lightning=0
  istat_nasalarc=0

  print_verbose=.false.
  if (verbose) print_verbose=.true.

  write(6,*)
  write(6,*) "analysis time is=",iyear,imonth,iday,ihour,iminute,isecond
  regional_time(1)=iyear
  regional_time(2)=imonth
  regional_time(3)=iday
  regional_time(4)=ihour
  regional_time(5)=iminute
  regional_time(6)=isecond

!  call load_gsdpbl_hgt(mype)
!
!  check consistency of the options
!

! Now either stratiform or cumulus cloud is considered in the cloud
!  water calculation. This leads to a limitation for the temperature
!  adjustment when stratiform cloud is chosen because adiabat profile
!  scheme based on the convection. This limitation may change when 
!  stratiform and cumulus cloud are both considered at the same time in the future.

  if(opt_cloudwaterice_retri == 1 .and. opt_cloudtemperature >= 4) then
     write(6,*) 'gsdcloudanalysis: ',&
       'inconsistent option for opt_cloudwaterice_retri and opt_cloudtemperature'
     write(6,*) 'gsdcloudanalysis: ',&
       'opt_cloudtemperature must be set to 3 when opt_cloudwaterice_retri =1'
     call stop2(113)
  endif
!!
!!----------------------------------------------
!! 2. read observations                  
!!----------------------------------------------
!!
!! 1.1   allocate observation fields
!!
!
  allocate(ref_mos_3d(lon2,lat2,nsig))
  allocate(ref_mos_3d_tten(lon2,lat2,nsig))
  ref_mos_3d=miss_obs_real
  ref_mos_3d_tten=miss_obs_real

  allocate(lightning(lon2,lat2))
  lightning=-9999.0_r_kind

  allocate(sat_ctp(lon2,lat2))
  allocate(sat_tem(lon2,lat2))
  allocate(w_frac(lon2,lat2))
  allocate(nlev_cld(lon2,lat2))
  sat_ctp=miss_obs_real
  sat_tem=miss_obs_real
  w_frac=miss_obs_real
  nlev_cld=miss_obs_int

  allocate(osfc_station_map(lon2,lat2))
  osfc_station_map=miss_obs_int
!
!!
!! 1.2 start to read observations                 
!!
  nmsclvl_radar = -999
  lunin=55

!!  1.2.2 read in surface observations
!!
!!  Note that OI is always 1 and OJ is the MPAS cell index
!!
  istart=2
  jstart=displ_1d(mype+1)+2
  numsao=0
!
  fileexist=.false.
  obsfile='mpas_metarcloud.bin'
  inquire(file=trim(obsfile),exist=fileexist)
  if(fileexist) then
     open(lunin, file=trim(obsfile),form='unformatted')
     read(lunin) obstype,isis,nreal,nchanl,ilat,ilon,ndata_model
     write(6,*) 'metar cloud=',obstype,isis,nreal,nchanl,ilat,ilon,ndata_model
     numsao=ndata_model
     allocate(oi(numsao))
     allocate(oj(numsao))
     allocate(ocld(nvarcld_p,numsao))
     allocate(owx(numsao))
     allocate(oelvtn(numsao))
     allocate(odist(numsao))
     allocate(cstation(numsao))
     allocate(oistation(numsao))
     allocate(ojstation(numsao))
     allocate(wimaxstation(numsao))
     call read_Surface(mype,lunin,istart,jstart,lon2,lat2, &
                       numsao,nvarcld_p,oi,oj,ocld,owx,oelvtn,&
                       odist,cstation,oistation,ojstation)

! Immediately deallocate oistation and ojstation b/c they contain incorrect information
     deallocate(oistation)
     deallocate(ojstation)

     if(mype == 0) write(6,*) 'gsdcloudanalysis: ',                                  &
                      'Surface cloud observations are read in successfully'
     istat_surface=1
     close(lunin)
  endif


  write(6,*)
  write(6,*) '========================================'
  write(6,*) 'Nonvar cloud analysis ran without error!'

  call MPI_FINALIZE(ierror)

end program cloudanalysis
