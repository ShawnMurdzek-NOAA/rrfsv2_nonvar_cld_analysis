module larccld_bufr_io

contains

subroutine read_NASALaRC_cloud_bufr(satfile,atime,satidgoeseast,satidgoeswest,east_time, west_time, &
             maxobs,numobs,ptop, teff, phase, lwp_iwp,lat, lon)
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2010-07-09
!
! ABSTRACT: 
!     This routine read in NASA LaRC cloud products 
!     from a bufr file                      
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
  use kinds, only: r_kind,i_kind

  implicit none

! Explicitly declare external functions and subroutines
  external :: openbf,dxdump,datelen,ufbint,closbf
  integer, external :: ireadmg,ireadsb
!
!
!
  character(80):: hdstr='YEAR  MNTH  DAYS HOUR  MINU  SECO SAID'
  character(80):: obstr='CLATH  CLONH CLDP HOCT CDTP EBBTH VILWC'
! CLDP     |  CLOUD PHASE
! HOCB     | HEIGHT OF BASE OF CLOUD
! HOCT     | HEIGHT OF TOP OF CLOUD             (METERS)
! CDBP     | PRESSURE AT BASE OF CLOUD
! CDTP     | PRESSURE AT TOP OF CLOUD           (PA)
! EBBTH    | EQUIVALENT BLACK BODY TEMPERATURE  (KELVIN)
! VILWC    | VERTICALLY-INTEGRATED LIQUID WATER CONTENT

  real(8) :: hdr(7),obs(7,1)

  character(8) subset
  integer :: unit_in=10,idate,iret,nmsg,ntb
  integer :: satid

!
!  For NASA LaRC 
!
  CHARACTER*40,intent(in) ::   satfile
!SATID
  integer,intent(in) :: satidgoeswest
  integer,intent(in) :: satidgoeseast  
  integer(8),intent(in) :: east_time, west_time

  INTEGER,intent(in)  ::   maxobs! dimension
  INTEGER,intent(out) ::   numobs  ! dimension
  INTEGER(8) ::  obs_time
  REAL*4      lat                            (  maxobs)
  REAL*4      lon                            (  maxobs)
  integer     phase                          (  maxobs)
  REAL*4      lwp_iwp                        (  maxobs)
  REAL*4      teff                           (  maxobs)
  REAL*4      ptop                           (  maxobs)
!
!
!  ** misc

!  integer i,j,k
!  Integer nf_status,nf_fid,nx,ny,nf_vid
!
!  integer :: status
  character*10  atime
!
!**********************************************************************
!
 open(24,file='NASA.bufrtable')
 open(unit_in,file=trim(satfile),form='unformatted')
 call openbf(unit_in,'IN',unit_in)
 call dxdump(unit_in,24)
 call datelen(10)
   nmsg=0
   ntb = 0
   msg_report: do while (ireadmg(unit_in,subset,idate) == 0)
     nmsg=nmsg+1
     sb_report: do while (ireadsb(unit_in) == 0)
       call ufbint(unit_in,hdr,7,1,iret,hdstr)
       obs_time=(hdr(1)-2000.0_8)*100000000.0_8+hdr(2)*1000000.0_8+hdr(3)*10000.0_8+hdr(4)*100.0_8+hdr(5)
       satid=int(hdr(7))
       if( (obs_time == east_time .and. satid==satidgoeseast ) .or.  &
           (obs_time == west_time .and. (satid==satidgoeswest .or. satid==259 .or. satid==271) ) ) then
         call ufbint(unit_in,obs,7,1,iret,obstr)
         if(abs(obs(3,1) -4.0) < 1.e-4) then
           obs(7,1)=99999. ! clear
           obs(6,1)=99999. ! clear
           obs(5,1)=101300.0  ! clear (hpa)
         endif
         if(obs(5,1) < 1.e7 .and. obs(5,1) > 100.0 ) then
         if(obs(6,1) < 1.e7 .and. obs(6,1) > 10.0) then
           ntb = ntb+1
           if(ntb > maxobs) then
             write(*,*) 'ALERT: need to increase maxobs',maxobs, ntb
             ntb = maxobs
           endif

           lwp_iwp(ntb)=99999.0
           lat(ntb)=99999.0
           lon(ntb)=99999.0
           phase(ntb)=99999
           teff(ntb)=99999.0
           ptop(ntb)=99999.0
           if(obs(1,1) < 1.e9) lat(ntb)=real(obs(1,1))
           if(obs(2,1) < 1.e9) lon(ntb)=real(obs(2,1))
           if(obs(3,1) < 1.e9) phase(ntb)=int(obs(3,1))
           if(obs(7,1) < 1.e9) lwp_iwp(ntb)=real(obs(7,1))
           if(obs(6,1) < 1.e9) teff(ntb)=real(obs(6,1))
           if(obs(5,1) < 1.e9) ptop(ntb)=real(obs(5,1))/100.0 ! pa to hpa

         endif
         endif
       endif   ! east_time, west_time
     enddo sb_report
   enddo msg_report
   write(*,*) 'message/reports num=',nmsg,ntb
 call closbf(unit_in)
 numobs=ntb
 write(atime,'(I10)') idate

end subroutine read_NASALaRC_cloud_bufr

subroutine read_NASALaRC_cloud_bufr_survey(satfile,satidgoeseast,satidgoeswest,east_time, west_time,maxobs)
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2010-07-09
!
! ABSTRACT: 
!     This routine read in NASA LaRC cloud products 
!     from a bufr file                      
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
  use kinds, only: r_kind,i_kind

  implicit none

! Explicitly declare external functions and subroutines
  external :: openbf,dxdump,datelen,ufbint,closbf
  integer, external :: ireadmg,ireadsb
!
  character(80):: hdstr='YEAR  MNTH  DAYS HOUR  MINU  SECO  SAID'
  real(8) :: hdr(7)

  character(8) subset
  integer :: unit_in=10,idate,iret,nmsg,ntb

!
!  For NASA LaRC 
!
  CHARACTER*40, intent(in)    ::   satfile
!SATID
  integer,intent(in) :: satidgoeswest
  integer,intent(in) :: satidgoeseast  
  integer(8),intent(out) :: east_time, west_time
  integer,intent(out) :: maxobs

  INTEGER(8) ::  obs_time

  INTEGER(i_kind),parameter :: max_obstime=30
  integer(8) :: num_obstime_all(max_obstime)
  integer(i_kind) :: num_subset_all(max_obstime)
  integer(i_kind) :: num_obstime_hh(max_obstime)
  integer(i_kind) :: num_satid(max_obstime)
  integer(i_kind) :: num_obstime

!
  integer :: i,ii,hhh
  integer :: numobs_east, numobs_west
  integer :: satid
!
!**********************************************************************
!
 num_obstime=0
 num_satid=0
 num_obstime_all=0
 num_subset_all=0
 obs_time=0
 hhh=99
 satid=0
 open(24,file='NASA.bufrtable')
 open(unit_in,file=trim(satfile),form='unformatted')
 call openbf(unit_in,'IN',unit_in)
 call dxdump(unit_in,24)
 call datelen(10)
   nmsg=0
   msg_report: do while (ireadmg(unit_in,subset,idate) == 0)
     ntb = 0
     nmsg=nmsg+1
     sb_report: do while (ireadsb(unit_in) == 0)
       call ufbint(unit_in,hdr,7,1,iret,hdstr)
       obs_time=(hdr(1)-2000.0_8)*100000000.0_8+hdr(2)*1000000.0_8+hdr(3)*10000.0_8+hdr(4)*100.0_8+hdr(5)
       hhh=int(hdr(5))
       ntb=ntb+1
       satid=int(hdr(7))
     enddo sb_report
! message inventory
     if(num_obstime == 0 ) then
       num_obstime=1
       num_obstime_all(num_obstime)=obs_time
       num_obstime_hh(num_obstime)=hhh
       num_subset_all(num_obstime)= ntb
       num_satid(num_obstime)=satid
     else
       ii=0
       DO i=1,num_obstime
          if(num_obstime_all(i) == obs_time .and. num_satid(i)== satid) ii=i
       ENDDO
       if( ii > 0 .and. ii <=num_obstime) then
          num_subset_all(ii)=num_subset_all(ii) + ntb
       else
          num_obstime=num_obstime+1
          if(num_obstime> max_obstime) then
             write(*,*) 'Error: too many message types'
             write(*,*) 'Need to increase :max_obstime'
             stop 1234
          endif
          num_obstime_all(num_obstime)=obs_time
          num_obstime_hh(num_obstime)=hhh
          num_subset_all(num_obstime)=num_subset_all(num_obstime)+ntb
          num_satid(num_obstime)=satid
       endif
     endif
   enddo msg_report
   write(*,*) 'message/reports num=',nmsg,ntb
 call closbf(unit_in)

 write(*,'(2x,a15,a15,a15,a15)') 'time_level','satid','subset_num','hour'
 DO i=1,num_obstime
   write(*,'(i2,i15,i15,i15,i15)') i,num_obstime_all(i),num_satid(i),num_subset_all(i),num_obstime_hh(i)
 ENDDO
!  GOES EAST  : 1815, 1845, 1915, 2045, no anymore, changed to 1830, 1900 just like WEST
!  GOES WEST  : 1830, 1900, 2030
 east_time=0
 west_time=0
 numobs_east=0
 numobs_west=0
 DO i=1,num_obstime
   if(num_subset_all(i) > 10) then
      if(num_satid(i) == satidgoeseast ) then  
         if(east_time < num_obstime_all(i)) then
              east_time=num_obstime_all(i)
              numobs_east=num_subset_all(i)
         endif
      endif
      if(num_satid(i) == satidgoeswest .or. num_satid(i)==259 .or. num_satid(i)==271 ) then
         if(west_time < num_obstime_all(i)) then
             west_time=num_obstime_all(i)
             numobs_west=num_subset_all(i)
         endif
      endif
   endif
 ENDDO
 write(*,*) 'east_time and number=',east_time,numobs_east
 write(*,*) 'west_time and number=',west_time,numobs_west
 
 maxobs=numobs_west+numobs_east
 maxobs=maxobs+int(maxobs*0.2)
 write(*,*) 'maxobs=',maxobs

end subroutine read_NASALaRC_cloud_bufr_survey

subroutine write_bufr_NASALaRC(bufrfile,idate,nlon,nlat,dx,index,w_pcld,w_tcld,w_frac,w_lwp,nlev_cld)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:   write_bufr_NASALaRC 
!   prgmmr: hu           org: essl/gsd                date: 2008-12-01
!   
! abstract: write NASA LaRC in RR grid into bufr
!   
! program history log:
!   2009-09-18  Hu
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  linux 
!
!$$$
    use constants, only: zero, one
    use kinds, only: r_kind,i_kind,r_single

    implicit none

! Explicitly declare external functions and subroutines
    external :: datelen,openbf,openmb,ufbint,writsb,closbf
!
    character(len=*), intent(in) :: bufrfile
    integer(i_kind),intent(in)   :: idate
    integer(i_kind), intent(in) :: nlon,nlat
    integer, intent(in)  ::   index(nlon,nlat)
    REAL(r_single), intent(in) ::   w_pcld(nlon,nlat)
    REAL(r_single), intent(in) ::   w_tcld(nlon,nlat)
    REAL(r_single), intent(in) ::   w_frac(nlon,nlat)
    REAL(r_single), intent(in) ::   w_lwp (nlon,nlat)
    REAL(r_kind),   intent(in) ::   dx
    INTEGER(i_kind),   intent(in) ::   nlev_cld(nlon,nlat)

    real(r_kind) :: hdr(5),obs(1,5)
    character(80):: hdrstr='SID XOB YOB DHR TYP'
    character(80):: obsstr='POB'

    REAL(i_kind),PARAMETER ::  MXBF = 160000_i_kind
    INTEGER(i_kind) :: ibfmsg = MXBF/4_i_kind

    character(8) subset,sid
    integer(i_kind) :: ludx,lendian_in

    INTEGER(i_kind)  ::  maxlvl, numref
    INTEGER(i_kind)  ::  i,j,iret


    write(6,*) 'cycle time is :', idate
    subset='ADPUPA'
    sid='NASALaRC'
    ludx=22
    lendian_in=10

    open(ludx,file='prepobs_prep.bufrtable',action='read')
    open(lendian_in,file=trim(bufrfile),action='write',form='unformatted')

    call datelen(10)
    call openbf(lendian_in,'OUT',ludx)
    maxlvl=5
    numref=0
!mhu    do j=1,nlat
!mhu    do i=1,nlon
! GSI has peroidic boundary which can crash the cloud analysis if there are data along the boundary.
! wait for GSI to fix this issue.
    do j=2,nlat-1
    do i=2,nlon-1
      if((index(i,j) .ge. 3) .or.  &
         (index(i,j) .ge. 1 .and. dx < 4000.0) ) then
        numref = numref + 1
        hdr(1)=transfer(sid,hdr(1))
        hdr(2)=float(i)/10.0_r_kind
        hdr(3)=float(j)/10.0_r_kind
        hdr(4)=0
        hdr(5)=500

        if( w_pcld(i,j) > 88888.0 .or. w_pcld(i,j) < -0.001) then
          obs(1,1)=9999.0
        else
          obs(1,1)=w_pcld(i,j)
        endif
        if( w_tcld(i,j) > 88888.0 .or. w_tcld(i,j) < -0.001) then
          obs(1,2)=9999.0
        else
          obs(1,2)=w_tcld(i,j)
        endif
        if( w_frac(i,j) > 88888.0 .or. w_frac(i,j) < -0.001) then
          obs(1,3)=9999.0
        else
          obs(1,3)=w_frac(i,j)*100.0
        endif
        if( w_lwp(i,j) > 88888.0 .or. w_lwp(i,j) < -0.001) then
          obs(1,4)=9999.0
        else
          obs(1,4)=w_lwp(i,j)*1000.0
        endif
        if( nlev_cld(i,j) > 88888.0 .or. nlev_cld(i,j) < -0.001) then
          obs(1,5)=9999.0
        else
          obs(1,5)=float(nlev_cld(i,j))
        endif


        call openmb(lendian_in,subset,idate)
        call ufbint(lendian_in,hdr,5,   1,iret,hdrstr)
        call ufbint(lendian_in,obs,1,maxlvl,iret,obsstr)
        call writsb(lendian_in,ibfmsg,iret)
      endif
    enddo   !i
    enddo   !j
    call closbf(lendian_in)
    write(6,*) 'write_to file',trim(bufrfile)
    write(6,*) 'write_bufr_nasaLaRC, DONE: write columns:',numref

end subroutine  write_bufr_NASALaRC

end module larccld_bufr_io
