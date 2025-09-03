subroutine read_prepbufr_metarcld(infile,analysis_time,analysis_minute,&
                                  twindin,nCell,lat_m,lon_m)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_prepbuf_emtarcld        read metar cld obs from prepbufr file
!   prgmmr: Ming Hu          org: GSL                date: 2020-08-02
!
! abstract:  This routine reads METAR cloud  data found in the prepbufr
!            file. 
!
! program history log:
!   2020-08-02  Hu

!   input argument list:
!     infile   - unit from which to read BUFR data
!     twindin  - input group time window (hours)
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_single,r_kind,i_kind,r_double
  use cld_parm_array_mod, only : obstype, sis, nchanl,nreal,ilat,ilon,ndata
  use cld_parm_array_mod, only : cdata_regular

! Modules from WPS
  use map_utils
  use misc_definitions_module

  implicit none

! Declare passed variables
  character(len=*)                      ,intent(in   ) :: infile
  real(r_kind)                          ,intent(in   ) :: twindin
  integer                               ,intent(in   ) :: analysis_time
  integer                               ,intent(in   ) :: analysis_minute
  integer                               ,intent(in   ) :: nCell
  real                                  ,intent(in   ),dimension(nCell) :: lat_m,lon_m

! Declare local parameters
  real(r_kind),parameter :: bmiss= 10.e10_r_kind
  real(r_kind),parameter :: r0_1_bmiss=0.1_r_kind*bmiss

! Declare local variables

  character(40) hdstr,qcstr
  character(40) metarcldstr,goescldstr,metarvisstr,metarwthstr
  character(80) obstr
  character(10) date
  character(8) subset
  character(8) c_station_id
  character(8) stnid

  integer(i_kind) ireadmg,ireadsb,iiout
  integer(i_kind) lunin,i,maxobs,j,nmsgmax,mxtb
  integer(i_kind) metarcldlevs,metarwthlevs
  integer(i_kind) nmsg                ! message index
  integer(i_kind) idx                 ! order index of aircraft temperature bias
  integer(i_kind) iyyyymm
  integer(i_kind),dimension(5):: idate5
  real(r_kind),allocatable,dimension(:,:):: cdata_all,cdata_out

  integer(i_kind) :: nread,ntb

  real(r_kind) rstation_id
  real(r_double) qcmark_huge
  real(r_double),dimension(8):: hdr,hdrtsb
  real(r_double),dimension(13,255):: obsdat
  real(r_double),dimension(8,255):: qcmark 
  real(r_double),dimension(2,10):: metarcld
  real(r_double),dimension(1,10):: metarwth
  real(r_double),dimension(2,1) :: metarvis

!  equivalence to handle character names
  equivalence(rstation_id,c_station_id)

!  data statements
  data hdstr  /'SID XOB YOB DHR TYP ELV SAID T29'/
  data obstr  /'POB QOB TOB ZOB UOB VOB PWO MXGS HOVI CAT PRSS TDO PMO' /
  data qcstr  /'PQM QQM TQM ZQM WQM NUL PWQ PMQ'/
  data metarcldstr /'CLAM HOCB'/      ! cloud amount and cloud base height
  data metarwthstr /'PRWE'/           ! present weather
  data metarvisstr /'HOVI TDO'/       ! visibility and dew point
                                                !   cloud top temperature, cloud top temp. qc mark

  data lunin / 13 /

  logical :: metarcldobs
  integer(i_kind) iret,irec
  integer(i_kind) kk,nc,kx
  integer(i_kind) iy,im,idd,ihh
  integer(i_kind) idate
  integer(i_kind) levs,k,qm,iout

  real(r_kind),parameter:: r90  = 90.0_r_kind
  real(r_kind),parameter:: r360 = 360.0_r_kind
  real(r_kind),parameter:: r180 = 180.0_r_kind

  real(r_kind) :: dlat_earth_deg,dlon_earth_deg
  real(r_kind) :: dlat_earth,dlon_earth
  real(r_kind) dlat,dlon,stnelev
  real(r_kind) timeobs
  real(r_kind) usage
  real(r_kind) deg2rad
  real            :: rlon,rlat,xc,yc
  integer(i_kind) :: ivalid
  integer(i_kind) :: minobs, minan

! Map projection
  type(proj_info) :: proj
  integer :: nlat, nlon
  real :: lat1, lon1, truelat1, truelat2, stdlon, dx, knowni, knownj
  real :: ll_lat, ur_lat, left_lat, bot_lat, right_lat, top_lat
  real :: ll_lon, ur_lon, left_lon, bot_lon, right_lon, top_lon

! MPAS coordinates in map projection space
  real,dimension(nCell) :: x_mp_m, y_mp_m

! Initialize variables

  nchanl=0
  nreal=27
  maxobs=10000
  obstype='mta_cld'
  metarcldobs = obstype == 'mta_cld'
  sis ="metarcld"
  deg2rad=3.1415926_r_kind/180.0_r_kind

!------------------------------------------------------------------------

! define map projection (similar to Lambert Conformal from HRRR)

  lat1 = 38.5
  lon1 = -97.5
  truelat1 = 38.5
  truelat2 = 38.5
  stdlon = -97.5
  dx = 3000.
! HRRR grid
!  nlat = 1060
!  nlon = 1800
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

!------------------------------------------------------------------------

! loop over convinfo file entries; operate on matches
  
  allocate(cdata_all(nreal,maxobs))
  cdata_all=0.0_r_kind
  nread=0
  ilon=2
  ilat=3

! read cloud observations from METAR
  open(lunin,file=infile,form='unformatted')
     call openbf(lunin,'IN',lunin)
     call datelen(10)

!    Big loop over prepbufr file	

     ntb  = 0
     nmsg = 0
     irec = 0
     iout = 0
     loop_msg: do while (ireadmg(lunin,subset,idate)== 0)
        irec = irec + 1
        loop_readsb: do while(ireadsb(lunin) == 0)
!          use msg lookup table to decide which messages to skip
!          use report id lookup table to only process matching reports
           ntb = ntb+1
                 
!          Extract type, date, and location information
           call ufbint(lunin,hdr,8,1,iret,hdstr)
           kx=hdr(5)
           if(kx < 180 .or. kx > 188) cycle loop_readsb 

           if(abs(hdr(3))>r90 .or. abs(hdr(2))>r360) cycle loop_readsb
           if(hdr(2) == r360)hdr(2)=hdr(2)-r360
           if(hdr(2) < 0.0_r_kind)hdr(2)=hdr(2)+r360
           dlon_earth_deg=hdr(2)
           dlat_earth_deg=hdr(3)
           dlon_earth=hdr(2)*deg2rad
           dlat_earth=hdr(3)*deg2rad
           rlon=dlon_earth_deg
           rlat=dlat_earth_deg

! convert to grid coordinate and check if it is inside the MPAS domain
           call latlon_to_ij(proj,rlat,rlon,xc,yc)
           dlon=xc
           dlat=yc
           if( (int(dlon) <= 0 .or. int(dlon) >= nlon) .or. &
               (int(dlat) <= 0 .or. int(dlat) >= nlat) ) cycle loop_readsb


!------------------------------------------------------------------------
!             in time correction for observations to account for analysis
!                      time being different from obs file time.
           write(date,'( i10)') idate
           read (date,'(i4,3i2)') iy,im,idd,ihh
           idate5(1)=iy
           idate5(2)=im
           idate5(3)=idd
           idate5(4)=ihh
           idate5(5)=0
           call w3fs21(idate5,minobs)    !  obs ref time in minutes relative to historic date
           minobs=minobs+hdr(4)*60.0_r_kind
!
           idate5(1)=analysis_time/1000000
           idate5(2)=(analysis_time-idate5(1)*1000000)/10000
           idate5(3)=(analysis_time-idate5(1)*1000000-idate5(2)*10000)/100
           idate5(4)=analysis_time-idate5(1)*1000000-idate5(2)*10000-idate5(3)*100
           idate5(5)=analysis_minute
           call w3fs21(idate5,minan)    !  analysis ref time in minutes relative to historic date
!           write(*,*) idate,hdr(4)
!           write(*,*) analysis_time,analysis_minute
           timeobs=minobs-minan
!           write(*,*) timeobs,twindin
           timeobs=timeobs/60.0_r_kind 
           if(abs(timeobs) > twindin) cycle loop_readsb

!          Extract data information on levels
           call ufbint(lunin,obsdat,13,255,levs,obstr)
           call ufbint(lunin,qcmark,8,255,levs,qcstr)
           if(metarcldobs) then
              metarcld=bmiss
              metarwth=bmiss
              metarvis=bmiss
              call ufbint(lunin,metarcld,2,10,metarcldlevs,metarcldstr)
              call ufbint(lunin,metarwth,1,10,metarwthlevs,metarwthstr)
              call ufbint(lunin,metarvis,2,1,iret,metarvisstr)
              if(levs /= 1 ) then
                 write(6,*) 'READ_PREPBUFR: error in Metar observations, levs sould be 1 !!!'
                 !call stop2(110)
                 stop 110
              endif
           endif
           ivalid=0
           do k=1,10
             if(metarcld(1,k) < 1.0e9) ivalid=ivalid+1 
             if(metarcld(2,k) < 1.0e9) ivalid=ivalid+1 
             if(metarwth(1,k) < 1.0e9) ivalid=ivalid+1 
           enddo
           if(metarvis(1,1) < 1.0e9) ivalid=ivalid+1
           if(ivalid==0)  cycle loop_readsb
           iout=iout+1
!
!          Set station ID
           rstation_id=hdr(1)
           stnelev=hdr(6)
           usage = 0.0_r_kind
           nc=kx
           qm=0      
!             extract aircraft profile information
           if(iout > maxobs) then
              write(6,*)'READ_PREPBUFR:  ***WARNING*** iout > maxobs',iout,maxobs
              iout=maxobs
           end if

! METAR cloud observation
           if(metarcldobs) then
                 cdata_all(1,iout)=rstation_id    !  station ID
                 cdata_all(2,iout)=dlon           !  map projection x coordinate (later changed to 1)
                 cdata_all(3,iout)=dlat           !  map projection y coordinate (later changed to closest MPAS cell index)
                 cdata_all(4,iout)=stnelev        !  station  elevation
                 if(metarvis(1,1) < r0_1_bmiss) then
                    cdata_all(5,iout)=metarvis(1,1)  !  visibility (m)
                 else
                    cdata_all(5,iout) = -99999.0_r_kind
                 endif
                 do kk=1, 6
                    if(metarcld(1,kk) < r0_1_bmiss) then
                       cdata_all(5+kk,iout) =metarcld(1,kk)  !  cloud amount
                    else
                       cdata_all(5+kk,iout) = -99999.0_r_kind
                    endif
                    if(metarcld(2,kk) < r0_1_bmiss) then
                       cdata_all(11+kk,iout)=metarcld(2,kk)  !  cloud bottom height (m)
                    else
                       cdata_all(11+kk,iout)= -99999.0_r_kind
                    endif
                 enddo
                 do kk=1, 3
                    if(metarwth(1,kk) < r0_1_bmiss) then
                       cdata_all(17+kk,iout)=metarwth(1,kk)  !  weather
                    else
                       cdata_all(17+kk,iout)= -99999.0_r_kind
                    endif
                 enddo
                 cdata_all(21,iout)=timeobs  !  time observation
                 cdata_all(22,iout)=usage
                 cdata_all(23,iout)=0.0_r_kind  ! reserved for distance between obs and grid
!     Calculate dewpoint depression from surface obs, to be used later
!         with haze and ceiling logic to exclude dust-caused ceiling obs
!         from cloud analysis
                 if(metarvis(2,1)  < 1.e10_r_kind) then
                    cdata_all(24,iout)=obsdat(3,1)-metarvis(2,1)  ! temperature - dew point
                 else
                    cdata_all(24,iout)=-99999.0_r_kind  ! temperature - dew point
                 endif
                 cdata_all(25,iout)=nc                     ! type
                 cdata_all(26,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
                 cdata_all(27,iout)=dlat_earth_deg         ! earth relative latitude (degrees)
           end if

!
!    End k loop over levs
        end do loop_readsb

!
!   End of bufr read loop
     enddo loop_msg
!    Close unit to bufr file
  call closbf(lunin)

! Normal exit

  write(*,*) 'number of valid cloud obs =',iout
  ndata=iout
  allocate(cdata_out(nreal,ndata))
  do i=1,ndata
     do k=1,nreal
        cdata_out(k,i)=cdata_all(k,i)
     end do
  end do
  deallocate(cdata_all)

! apply map projection to MPAS mesh

  do i=1,nCell
     rlon=lon_m(i)
     rlat=lat_m(i)
     call latlon_to_ij(proj,rlat,rlon,xc,yc)
     x_mp_m(i)=xc
     y_mp_m(i)=yc
  enddo

! define a closest METAR cloud observation for each grid point

  if(metarcldobs .and. ndata > 0) then
     maxobs=2000000
     allocate(cdata_all(nreal,maxobs))
     call reorg_metar_cloud_regular(cdata_out,nreal,ndata,nlat,nlon,nCell,lat_m,lon_m,&
                                    x_mp_m,y_mp_m,cdata_all,maxobs,iout)
     ndata=iout
     deallocate(cdata_out)
     allocate(cdata_regular(nreal,ndata))
     do i=1,nreal
        do j=1,ndata
          cdata_regular(i,j)=cdata_all(i,j)
        end do
     end do
     deallocate(cdata_all)
  endif

  call closbf(lunin)

  close(lunin)

! End of routine
  return

end subroutine read_prepbufr_metarcld

