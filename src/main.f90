program fireHK

! Program to calculate fire danger index time series in Hong Kong (Leo Lai, Mar 2022 HKU)
! Code adapted from original weathergen
! Currently only read time series input (single grid netCDF file)

! To compile
! ~ make clean && make

! To run
! ~ ./fireHK input.nc

! To print output to text file
! ~ ./fireHK input.nc > output.txt

use parametersmod,  only : sp,dp,i4,ndaymonth
use statevarsmod,   only : monvars,dayvars
use netcdfinputmod, only : metdatainput
use gwgenmod,       only : gwgen_new
use orbitmod,       only : orbit,calcorbitpars
use radiationmod,   only : elev_Ratm,calcPjj,radpet
use diurnaltempmod, only : diurnaltemp,humidity
use ffdimod,        only : calcFFDI
use gfdimod,        only : calcGFDI

implicit none

character(200) :: infile

real(dp) :: lat
real(sp) :: elev

integer(i4) :: yr
integer(i4) :: calcyrs
integer(i4) :: startyr
integer(i4) :: endyr
integer(i4) :: s
integer(i4) :: e
integer(i4) :: ndyear

integer :: mon
integer :: day
integer :: m
integer :: d

!---------------------------------------------------------------------

! Get netCDF infile name
call getarg(1,infile)

!---------------
! Currently hard coded with start and end year of input file
startyr = 1971
endyr   = 2100
calcyrs = endyr - startyr + 1

lat     = 22.3        ! Prescribe latitude of Hong Kong
elev    = 0.          ! Prescribe elevation of Hong Kong

!---------------
! Initialize monthly input variables by allocating monvars to number of months
! Monthly and daily variables saved in statevarmod
call initmonvars(monvars,startyr,calcyrs)

!---------------
! Read and store input data
call metdatainput(infile,calcyrs)

!---------------
! Begin year loop for calculation of daily meterology
! Decision to run on a yearly time-step because Grassland Fire Danger Index (GFDI) calculation
! requires daily precipitation accumulation within the year
yearloop: do yr = 1, calcyrs

  ! Find start and end array index of the first and last month of current year
  s = (yr * 12) + 1
  e = s + 11

  ! Find the number of days in current year (365 or 366)
  ndyear = sum(monvars%nd(s:e))

  ! Allocate / deallocate dayvars within each yearloop iteration
  allocate(dayvars(ndyear))

  ! Generate daily meterology
  call gwgen_new(yr,s,e)

  ! Calculate orbital parameters for year for radpet routine
  call calcorbitpars(startyr,yr,orbit)

  !---------------

  dayloop : do d = 1, ndyear

    ! Calculate diurnal temperature and daylength
    ! This routine requires tmin of the next day, in ALVAR, gwgen interpolate Jan of next year as well to obtain tmin of 1 Jan
    ! here I am just copying the tmin of 31 Dec (for now)
    call diurnaltemp(d,ndyear,lat)

    ! Calculate relative atmospheric pressure for radpet routine (constant at the moment)
    call elev_Ratm(d,elev)

    ! Calculate precipitation equitability index for radpet routine
    call calcPjj(d,s,e)

    ! Calculate dpet and hence tdew for humidity routine
    call radpet(d,s,e,lat)

    ! Calculate relative humidity
    call humidity(d)
    !
    ! Calculate daily Forest Fire Danger Index
    call calcFFDI(yr,d,s,e)
    !
    ! Calculate daily Grassland Fire Danger Index
    call calcGFDI(d,ndyear)

  end do dayloop

  !---------------
  ! Print daily output to terminal
  mon = 1
  day = 1

  do m = s, e

    do d = 1, monvars%nd(m)

      print *, (startyr+yr-1), &
               mon,&
               d,&
               day,&
               dayvars(day)%tmin,&
               dayvars(day)%tmax,&
               dayvars(day)%tmean,&
               dayvars(day)%tdew,&
               dayvars(day)%prec,&
               dayvars(day)%cldf,&
               dayvars(day)%wind,&
               dayvars(day)%rhum,&
               dayvars(day)%prec_accum,&
               dayvars(day)%curing,&
               dayvars(day)%fuelmc,&
               dayvars(day)%KBDI,&
               dayvars(day)%FFDI,&
               dayvars(day)%FFDI_cat,&
               count(dayvars%FFDI_cat >= 3),&
               dayvars(day)%GFDI,&
               dayvars(day)%GFDI_cat,&
               count(dayvars%GFDI_cat >= 3)

      day = day + 1

    end do

    mon = mon + 1

  end do

  !---------------

  deallocate(dayvars)

  !---------------

end do yearloop


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

contains

subroutine initmonvars(monvars,startyr,calcyrs)

use parametersmod, only : ndaymonth
use statevarsmod,  only : mon_metvars

implicit none

type(mon_metvars), intent(inout) :: monvars
integer(i4)      , intent(in)    :: startyr
integer(i4)      , intent(in)    :: calcyrs

integer :: cntt     ! Number of months including buffer years
integer :: year
integer :: yr
integer :: i

cntt = (calcyrs * 12) + 24

allocate(monvars%tmp(cntt))
allocate(monvars%dtr(cntt))
allocate(monvars%pre(cntt))
allocate(monvars%wet(cntt))
allocate(monvars%cld(cntt))
allocate(monvars%wnd(cntt))
allocate(monvars%nd(cntt))

i = 1

do yr = 1, calcyrs+2    ! Account for start and end buffer years

  year = startyr + yr - 2

  do m = 1, 12

    monvars%nd(i) = ndaymonth(year,m)

    i = i + 1

  end do

end do

end subroutine initmonvars

!---------------------------------------------------------------------


end program fireHK
