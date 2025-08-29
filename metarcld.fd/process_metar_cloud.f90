program  process_metar_cloud
!
!   PRGMMR: Ming Hu          ORG: GSD        DATE: 2009-09-04
!
! ABSTRACT: 
!     This routine read in NASA LaRC cloud products and 
!     interpolate them into GSI mass grid
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
  use mpi
  use kinds, only: r_kind,i_kind,r_single
  use mpasio, only: read_MPAS_dim,read_MPAS_lat_lon

  use cld_parm_array_mod, only: region_dy,region_dx
  use cld_parm_array_mod, only: metar_impact_radius
  use cld_parm_array_mod, only: l_metar_impact_radius_change, &
                            metar_impact_radius_max,metar_impact_radius_min,&
                            metar_impact_radius_max_height,metar_impact_radius_min_height
  use cld_parm_array_mod, only: init_cld_parm

  use cld_parm_array_mod, only : obstype, sis, nchanl,nreal,ilat,ilon,ndata
  use cld_parm_array_mod, only : cdata_regular

  implicit none
!
  INCLUDE 'netcdf.inc'
!
! MPI variables
  integer :: npe, mype,ierror
!
! MPAS mesh
  integer(i_kind) :: nCell
  real, allocatable :: lat_m(:),lon_m(:)
  CHARACTER*180   meshfile
  CHARACTER*50    mpasfield
!
!     ****VARIABLES FOR THIS NETCDF FILE****
!
! namelist
!
  integer :: analysis_time
  integer :: analysis_minute
  integer :: debug
  character(len=100) :: prepbufrfile
  real(r_kind)       :: twindin
  namelist/setup/ analysis_time,analysis_minute,prepbufrfile,twindin,&
                  metar_impact_radius,l_metar_impact_radius_change, &
                  metar_impact_radius_max,metar_impact_radius_min, &
                  metar_impact_radius_max_height,metar_impact_radius_min_height,&
                  debug
!
!  ** misc
      
  logical :: ifexist
  integer :: lunout

! For debugging
  integer :: i,icell
  real(r_kind) :: sta_lon,rstation_id
  character(8) :: c_station_id
  equivalence(rstation_id,c_station_id)

!**********************************************************************
!
!            END OF DECLARATIONS....start of program
!
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

  if(mype==0) then
     lunout=15
     call init_cld_parm
!
!  get namelist
!
     analysis_time=2018051718
     prepbufrfile='prepbufr'
     analysis_minute=0
     twindin=0.5
     debug=0
 
     inquire(file='namelist.metarcld', EXIST=ifexist )
     if(ifexist) then
       open(10,file='namelist.metarcld',status='old')
          read(10,setup)
       close(10)
       write(*,*) 'Namelist setup are:'
       write(*,setup)
     else
       write(*,*) 'No namelist file exist, use default values'
       write(*,*) analysis_time
     endif
!
! read MPAS mesh information
!
     meshfile='mesh.nc'
     mpasfield='nCells'
     call read_MPAS_dim(meshfile, mpasfield, nCell)
     write(6,*)
     write(6,*) 'model nCell   =', nCell
     allocate(lat_m(nCell))
     allocate(lon_m(nCell))
     call read_MPAS_lat_lon(meshfile, nCell, lat_m, lon_m)
     write(6,*) 'min model lat =', minval(lat_m)
     write(6,*) 'min model lon =', minval(lon_m)
     write(6,*) 'max model lat =', maxval(lat_m)
     write(6,*) 'max model lon =', maxval(lon_m)
     write(6,*)
!
     call read_prepbufr_metarcld(prepbufrfile,analysis_time,analysis_minute,&
                                 twindin,nCell,lat_m,lon_m)

     write(*,*) 'number of cloud on MPAS mesh=',ndata
     write(*,*) obstype,nreal,nchanl,ilat,ilon,sis
     open(lunout,file='mpas_metarcloud.bin',form='unformatted')
        write(lunout) obstype,sis,nreal,nchanl,ilat,ilon,ndata
        write(lunout) cdata_regular
     close(lunout)
!
! DEBUGGING: Write out processed METAR obs to text file
!
     if (debug > 0) then
       open(13, file="processed_metar_obs_mpas.txt", status="new", action="write")
       write(13,'(13a12)') 'icell', 'MPAS_lat', 'MPAS_lon', 'sta_lat', 'sta_lon', 'min_dist', 'ID', &
              'cover1', 'cover2', 'cover3', 'base1', 'base2', 'base3' 
       do i=1,ndata
          icell = cdata_regular(24,i)
          sta_lon = cdata_regular(26,i) - 360.
          rstation_id = cdata_regular(1,i)
          write(13,'(I12,5f12.3,a12,6f12.3)') icell, lat_m(icell), lon_m(icell), &
                 cdata_regular(27,i), sta_lon, cdata_regular(23,i), trim(c_station_id), &
                 cdata_regular(6,i), cdata_regular(7,i), cdata_regular(8,i), &
                 cdata_regular(12,i), cdata_regular(13,i), cdata_regular(14,i)
                  
       enddo
       close(13)
     endif

     deallocate(cdata_regular)
!
     write(6,*) "=== RAPHRRR PREPROCCESS SUCCESS ==="

  endif ! if mype==0 

  call MPI_FINALIZE(ierror)

end program process_metar_cloud
