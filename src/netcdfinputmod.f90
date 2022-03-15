module netcdfinputmod

use parametersmod, only : i2,i4,sp,dp
use statevarsmod,  only : monvars
use errormod,      only : ncstat,netcdf_err
use getdatamod,    only : readdata
use netcdf

implicit none

! character(100), parameter :: soilfile = './soil30m_list-formatted.nc'
! character(100), parameter :: topofile = './world_slope_1800arcsec_list-formatted.nc'
! character(100), parameter :: LAIfile  = './LAI_mean_monthly_1981-2015_halfdeg_list-formatted.nc'
! character(100), parameter :: topofile = '/home/leolaio/grid2list_hktopo/topofile_integrated_list-formatted.nc'

contains

!---------------------------------------------------------------------

subroutine metdatainput(infile,calcyrs)

! Read in the full array (all months for all gridcells) of monthly variables from netcdf file

character(100), intent(in) :: infile
integer(i4),    intent(in) :: calcyrs

integer :: ifid
integer :: dimid
integer :: varid

! monthly input driver variables
real(sp), allocatable, dimension(:) :: tmp        ! mean monthly temperature (degC)
real(sp), allocatable, dimension(:) :: dtr        ! mean monthly diurnal temperature range (degC)
real(sp), allocatable, dimension(:) :: pre        ! total monthly precipitation (mm)
real(sp), allocatable, dimension(:) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
real(sp), allocatable, dimension(:) :: cld        ! mean monthly cloud cover (fraction)
real(sp), allocatable, dimension(:) :: wnd        ! mean monthly 10m windspeed (m s-1)

integer(i4), allocatable, dimension(:) :: var_in

integer(i4), dimension(3) :: srt
integer(i4), dimension(3) :: cnt
integer(i4) :: cntt
integer(i4) :: p0
integer(i4) :: p1

real(sp)    :: scale_factor
real(sp)    :: add_offset
integer(i2) :: missing_value

integer :: i


!----------------------------------------------------
! Pointers to mpi derived type variables
! gridcount = 1

! -----------------------------------------------------
! INPUT: Read dimension IDs and lengths of dimensions

! ncstat = nf90_open(infile,nf90_nowrite,ifid,comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)           ! Open netCDF-file (inpput file name, no writing rights, assigned file number)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)      ! Check for errors (after every step)
!
! ncstat = nf90_inq_dimid(ifid,'lon',dimid)              ! get dimension ID from dimension 'lon' in the input file
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inquire_dimension(ifid,dimid,len=xlen)   ! get dimension name and length from input file for dimension previously inquired
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inq_dimid(ifid,'lat',dimid)              ! get dimension ID from dimension 'lon' in the input file
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inquire_dimension(ifid,dimid,len=ylen)   ! get dimension name and length from input file for dimension previously inquired
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inq_dimid(ifid,'index',dimid)             ! get dimension ID from dimension 'index' in the input file
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inquire_dimension(ifid,dimid,len=ilen)   ! get dimension name and length from input file for dimension previously inquired
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inq_dimid(ifid,'time',dimid)             ! Get dimension ID for time
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inquire_dimension(ifid,dimid,len=tlen)   ! Get length of dimension 'time' and assign it to variable tlen
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----------------------------------------------------
! Read variable IDs and values

! allocate(lon(xlen))         ! Allocate length to longitude array
! allocate(lat(ylen))         ! Allocate length to latitude array
! allocate(time(tlen))        ! Allocate length to time array
! allocate(indx(xlen,ylen))   ! Allocate length to index array

! ncstat = nf90_inq_varid(ifid,"lon",varid)                ! Get variable ID for longitude
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_var(ifid,varid,lon)                    ! Get variable values for longitude
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_inq_varid(ifid,"lat",varid)                ! Get variable ID for latitude
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_var(ifid,varid,lat)                    ! Get variable values for latitude
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_inq_varid(ifid,"index",varid)               ! Get variable ID for index
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_var(ifid,varid,indx)                    ! Get variable values for index
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ncstat = nf90_inq_varid(ifid,"time",varid)               ! Get variable ID for time
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_var(ifid,varid,time)                   ! Get variable values for time
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----------------------------------------------------
! Read in elevation data from input file

! allocate(var_in(ilen))
! allocate(elev(ilen))      ! Allocate length to index array
!
! ncstat = nf90_inq_varid(ifid,'elv',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_var(ifid,varid,var_in,start=[srt(1)],count=[cnt(1)])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_att(ifid,varid,"missing_value",missing_value)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_att(ifid,varid,"scale_factor",scale_factor)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_att(ifid,varid,"add_offset",add_offset)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! where (var_in /= missing_value)
!   elev = real(var_in) * scale_factor + add_offset
! ! elsewhere
! !   elev = -9999.
! end where

! ----------------------------------------------------
! Read in variables from netcdf file by each gridcell
srt = [1,1,1]
cnt = [1,1,(calcyrs*12)]
cntt = (calcyrs * 12) + 24
p0 = 13
p1 = cntt - 12

ncstat = nf90_open(infile,nf90_nowrite,ifid)           ! Open netCDF-file (inpput file name, no writing rights, assigned file number)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(tmp(cntt))
allocate(dtr(cntt))
allocate(pre(cntt))
allocate(wet(cntt))
allocate(cld(cntt))
allocate(wnd(cntt))

call readdata(ifid,'tmp',srt,cnt,tmp(p0:p1))
call readdata(ifid,'dtr',srt,cnt,dtr(p0:p1))
call readdata(ifid,'pre',srt,cnt,pre(p0:p1))
call readdata(ifid,'wet',srt,cnt,wet(p0:p1))
call readdata(ifid,'cld',srt,cnt,cld(p0:p1))
call readdata(ifid,'wnd',srt,cnt,wnd(p0:p1))

! write(0,*) 'ok read data', xlen, ylen, tlen, ilen, cnt

!---------------------------------------------------------------------
! copy first and last year into buffers at each end
if (p0 == 13) then
  tmp(1:12) = tmp(13:24)
  dtr(1:12) = dtr(13:24)
  pre(1:12) = pre(13:24)
  wet(1:12) = wet(13:24)
  cld(1:12) = cld(13:24)
  wnd(1:12) = wnd(13:24)
end if

if (cntt > p1) then
  tmp(p1+1:cntt) = tmp(p1-11:p1)
  dtr(p1+1:cntt) = dtr(p1-11:p1)
  pre(p1+1:cntt) = pre(p1-11:p1)
  wet(p1+1:cntt) = wet(p1-11:p1)
  cld(p1+1:cntt) = cld(p1-11:p1)
  wnd(p1+1:cntt) = wnd(p1-11:p1)
end if

!---------------------------------------------------------------------
! enforce reasonable values of prec and wetdays

where (pre > 0.1)
  wet = max(wet,1.)
  wet = min(wet,10. * pre)
elsewhere
  pre = 0.
  wet = 0.
end where

!---------------------------------------------------------------------
! Copy monthly time series into the mon_metvar derived type structure
monvars%tmp(:) = tmp(:)
monvars%dtr(:) = dtr(:)
monvars%pre(:) = pre(:)
monvars%wet(:) = wet(:)
monvars%cld(:) = cld(:)
monvars%wnd(:) = wnd(:)

ncstat = nf90_close(ifid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)


end subroutine metdatainput

!---------------------------------------------------------------------

end module netcdfinputmod
