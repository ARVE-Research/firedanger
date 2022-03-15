module gwgenmod

! Use the Makefile to compile this program

!! Program to run gwgen with gridded input, provide global lon/lat indexed list output of
!! absolute minimum and maximum temperature
! JO Kaplan, HKU, 2019; Leo O Lai, HKU, 2021

! Terminal command line: mpirun -np 18 ./src/gwgen ~/path/to/input.file startyr/calcyr ~/path/to/outfile

! List of modules that will be used and the variables within these modules that are used in this program:

use parametersmod, only : sp,dp,i4,i2
use statevarsmod,  only : monvars,dayvars
use errormod,      only : ncstat,netcdf_err
use geohashmod,    only : geohash
use randomdistmod,
use newsplinemod,  only : newspline_all
use weathergenmod, only : metvars_in,metvars_out,weathergen,roundto
use netcdf

implicit none

public :: gwgen_new

contains

! Inquire about the dimensions of the input file

!-------------------------------------------------------

subroutine gwgen_new(yr,smon,emon)

implicit none

integer(i4), intent(in) :: yr
integer(i4), intent(in) :: smon
integer(i4), intent(in) :: emon

! Pointer variables to genvars

real(sp), pointer, dimension(:) :: tmp        ! mean monthly temperature (degC)
real(sp), pointer, dimension(:) :: dtr        ! mean monthly diurnal temperature range (degC)
real(sp), pointer, dimension(:) :: pre        ! total monthly precipitation (mm)
real(sp), pointer, dimension(:) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
real(sp), pointer, dimension(:) :: cld        ! mean monthly cloud cover (fraction)
real(sp), pointer, dimension(:) :: wnd        ! mean monthly 10m windspeed (m s-1)
integer,  pointer, dimension(:)  :: nd

! monthly derived driver variables
! Process 12 months monthly series (with +/- 4 months buffer) in each iteration

real(sp), dimension(20) :: mtmin      ! maximum monthly temperature (degC)
real(sp), dimension(20) :: mtmax      ! monthly minimum temperature (degC)
real(sp), dimension(20) :: wetf       ! fraction of wet days in a month

! Elements to calculate current year and amount of days in current month

integer :: i_count,outd
integer :: i,k,t,d,m
integer :: d0
integer :: d1
integer :: ndm


! Variables for the smoothing process

integer, parameter :: w = 3              ! filter half-width for smoothing of monthly mean climate variables to pseudo-daily values (months)
integer, parameter :: wbuf = 31*(1+2*w)  ! length of the buffer in which to hold the smoothed pseudo-daily  meteorological values (days)

integer,  dimension(-w:w) :: ndbuf       ! number of days in the month
real(sp), dimension(-w:w) :: mtminbuf    ! monthly minimum temperature
real(sp), dimension(-w:w) :: mtmaxbuf    ! monthly maximum temperature
real(sp), dimension(-w:w) :: cldbuf      ! monthly cloud fractions
real(sp), dimension(-w:w) :: wndbuf      ! monthly wind speed

real(sp), dimension(wbuf) :: tmin_sm     ! smoothed pseudo-daily values of min temperature
real(sp), dimension(wbuf) :: tmax_sm     ! smoothed pseudo-daily values of max temperature
real(sp), dimension(wbuf) :: cld_sm      ! smoothed pseudo-daily values of cloudiness
real(sp), dimension(wbuf) :: wnd_sm      ! smoothed pseudo-daily values of wind speed

! quality control variables

integer  :: mwetd_sim    ! simulated number of wet days
real(sp) :: mprec_sim    ! simulated total monthly precipitation (mm)

integer  :: pdaydiff     ! difference between input and simulated wet days
real(sp) :: precdiff     ! difference between input and simulated total monthly precipitation (mm)

real(sp) :: prec_t       ! tolerance for difference between input and simulated total monthly precipitation (mm)
integer, parameter  :: wetd_t = 10  ! tolerance for difference between input and simulated wetdays (days)

integer  :: pdaydiff1 = huge(i4)       ! stored value of the best match difference between input and simulated wet days
real(sp) :: precdiff1 = huge(0.0_sp)   ! stored value of the difference between input and simulated total monthly precipitation (mm)

! data structures for meteorology

type(metvars_in)  :: met_in   ! structure containing one day of meteorology input to weathergen
type(metvars_out) :: met_out  ! structure containing one day of meteorology output from weathergen

type(metvars_out), dimension(31) :: month_met  ! buffer containing one month of simulated daily meteorology

real(sp) :: mtmin_sim
real(sp) :: mtmax_sim
real(sp) :: mcldf_sim
real(sp) :: mwind_sim

real(sp) :: prec_corr
real(sp) :: tmin_corr
real(sp) :: tmax_corr
real(sp) :: cldf_corr
real(sp) :: wind_corr

integer :: start
integer :: end

integer :: baddata_check

!---------------------------------------------------------------------
! point to variables

tmp => monvars%tmp(smon-4:emon+4)
dtr => monvars%dtr(smon-4:emon+4)
pre => monvars%pre(smon-4:emon+4)
wet => monvars%wet(smon-4:emon+4)
cld => monvars%cld(smon-4:emon+4)
wnd => monvars%wnd(smon-4:emon+4)
nd  => monvars%nd(smon-4:emon+4)

! Assign the module random state to the current met_in
! The randomstate generated in last month of the year (Dec) will be saved back into the module variable
! This will carry onto the next year iteration, until all calcyears are finished
! at which the random state will be re-initiated by the genrndstate() in modelmod after moving on to next gridcell
if (yr == 1) met_in%rndst = georndst

!---------------------------------------------------------------------
! calculate derived climate variables

! Correct for negative diurnal temperature range values (bad data in transient file???) (Leo Lai Apr 2021)
dtr(:) = abs(dtr(:))

mtmin = tmp(:) - 0.5 * dtr(:)
mtmax = tmp(:) + 0.5 * dtr(:)
wetf  = wet(:) / nd

!--- Checking bad data (Leo)
!--- Cycling bad data (Leo)

baddata_check = 0
!
do k = 1, 20

  if(mtmin(k) < -273.15) then
    mtmin(k) = (mtmin(k-1) + mtmin(k+1)) / 2
    mtmax(k) = (mtmax(k-1) + mtmax(k+1)) / 2

    baddata_check = baddata_check + 1
  end if
  !---
  if (pre(k) < 0.) then
    pre(k) = 0.
  end if
  !---
  if (wnd(k) < 0.) then
    wnd(k) = 0.
  end if
  !---
  if (cld(k) < 0. .or. cld(k) > 100.) then
    cld(k) = 0.
  end if

end do
!
! if(baddata_check /= 0) then
!   cycle
! end if


!---------------------------------------------------------------------
! prepare pseudo-daily smoothed meteorological variables
! initialize the smoothing buffer variables with all values from the first month of the input
! since we will always start in year 2, we can do:

t = 5     ! First month (i.e. Jan) located at index 5 due to 20 months buffered 'genvars'

ndbuf    = nd(t-w:t+w)
mtminbuf = mtmin(t-w:t+w)
mtmaxbuf = mtmax(t-w:t+w)
cldbuf   = cld(t-w:t+w)
wndbuf   = wnd(t-w:t+w)

met_out%pday(1) = .false.
met_out%pday(2) = .false.
met_out%resid = 0.

start = 1     ! initiate index value for storing into dayvars

! start month loop
! 13 iteration of month loop from Jan to next year Jan
! Extra month (Jan) required for simulation of diurnal temperature, where 1st Jan tmin is needed for 31st Dec

monthloop : do m = 1, 12

  t = m + 4   ! first month (Jan) at index 5 due to 4 months of buffer at the beginning

  ! generate pseudo-daily smoothed meteorological variables (using means-preserving algorithm)

  call newspline_all(mtminbuf,ndbuf,tmin_sm(1:sum(ndbuf)))
  call newspline_all(mtmaxbuf,ndbuf,tmax_sm(1:sum(ndbuf)))
  call newspline_all(cldbuf,ndbuf,cld_sm(1:sum(ndbuf)))
  call newspline_all(wndbuf,ndbuf,wnd_sm(1:sum(ndbuf)))

  where (tmin_sm > tmax_sm)
    tmax_sm = tmin_sm + 0.1
  end where

  do i = 1,sum(ndbuf)
    if (tmin_sm(i) > tmax_sm(i)) print *, tmin_sm(i), tmax_sm(i)
  end do

  ! calculcate start and end positons of the current month pseudo-daily buffer

  d0 = sum(ndbuf(-w:-1)) + 1
  d1 = d0 + ndbuf(0) - 1

  ndm = d1 - d0 + 1

  ! restrict simulated total monthly precip to +/-10% or 1 mm of observed value

  prec_t = max(1.,0.1 * pre(t))

  i_count = 0

      !---------------------------------------------------------------------------------
      ! quality control loop calling the weathergen - this loop principally checks that
      ! the number of wet days and total precip stayed close to the input data

  qualityloop : do

    i_count = i_count + 1    ! increment iteration number

    mwetd_sim = 0
    mprec_sim = 0.

    outd = 1

    dayloop : do d = d0,d1  ! day loop

      met_in%prec  = pre(t)
      met_in%wetd  = wet(t)
      met_in%wetf  = wetf(t)
      met_in%tmin  = tmin_sm(d)
      met_in%tmax  = tmax_sm(d)
      met_in%cldf  = real(cld_sm(d))
      met_in%wind  = real(wnd_sm(d))
      met_in%pday  = met_out%pday
      met_in%resid = met_out%resid

      call weathergen(met_in,met_out)

      met_in%rndst = met_out%rndst
      month_met(outd) = met_out    ! save this day into a month holder

      if (met_out%prec > 0.) then

        mwetd_sim = mwetd_sim + 1
        mprec_sim = mprec_sim + met_out%prec

      end if

      outd = outd + 1

    end do dayloop ! day loop

    ! quality control checks

    if (pre(t) == 0.) then ! if there is no precip in this month a single iteration is ok

      pdaydiff = 0
      precdiff = 0.

      exit

    else if (i_count < 2) then

      cycle  !enforce at least two times over the month to get initial values for residuals ok

    else if (pre(t) > 0. .and. mprec_sim == 0.) then

      cycle  ! need to get at least some precip if there is some in the input data

    end if

    pdaydiff = abs(mwetd_sim - wet(t))

    precdiff = (mprec_sim - pre(t)) / pre(t)

    ! if (i_count > 990) print *, mwetd_sim, wet(t), pdaydiff1, pdaydiff, pre(t), precdiff, t, m      ! Debugging

    if (pdaydiff <= wetd_t .and. precdiff <= prec_t) then

      exit

    else if (pdaydiff < pdaydiff1 .and. precdiff < precdiff1) then

      ! save the values you have in a buffer in case you have to leave the loop
      ! should save the entire monthly state so that the "closest" acceptable value
      ! could be used in the event of needing a very large number of iteration cycles

      pdaydiff1 = pdaydiff
      precdiff1 = precdiff

    else if (i_count > 1000) then

      write (*,*) "No good solution found after 1000 iterations."

      stop

    end if

  end do qualityloop

  ! end of quality control loop
  !---------------------------------------------------------------------------------

  ! adjust meteorological values to match the input means following Richardson & Wright 1984

  mtmin_sim = sum(month_met(1:ndm)%tmin) / ndm
  mtmax_sim = sum(month_met(1:ndm)%tmax) / ndm
  mcldf_sim = sum(month_met(1:ndm)%cldf) / ndm
  mwind_sim = sum(month_met(1:ndm)%wind) / ndm

  if (mprec_sim == 0.) then
    if (pre(t) > 0.) stop 'simulated monthly prec = 0 but input prec > 0'
    prec_corr = 1.
  else
    prec_corr = pre(t) / mprec_sim
  end if

  tmin_corr = mtmin(t) - mtmin_sim
  tmax_corr = mtmax(t) - mtmax_sim

  mcldf_sim = 100. * mcldf_sim    ! Convert fraction back to percentage (0 to 100) (Leo)

  if (mcldf_sim == 0.) then
    if (cld(t) > 0.) stop 'simulated monthly cloud = 0 but input cloud > 0'
    cldf_corr = 1.
  else
    cldf_corr = cld(t) / mcldf_sim
  end if

  ! if (mwind_sim == 0.) then
  !   if (wnd(t) > 0.) stop 'simulated monthly wind = 0 but input wind > 0'
  !   wind_corr = 1.
  ! else
  !   wind_corr = wnd(t) / mwind_sim
  ! end if

  !--- Replaced "stop" command for now since only tmin and tmax is wanted (Leo)
  if (mwind_sim == 0.) then
    if (wnd(t) > 0.) then
      ! write(0,*) 'Warning: simulated monthly wind = 0 but input wind > 0'
      wind_corr = 1.
    end if
  else
    wind_corr = wnd(t) / mwind_sim
  end if

  !---

  month_met(1:ndm)%prec = month_met(1:ndm)%prec * prec_corr
  month_met(1:ndm)%tmin = month_met(1:ndm)%tmin + tmin_corr
  month_met(1:ndm)%tmax = month_met(1:ndm)%tmax + tmax_corr
  month_met(1:ndm)%cldf = month_met(1:ndm)%cldf * cldf_corr
  month_met(1:ndm)%wind = month_met(1:ndm)%wind * wind_corr
  month_met(1:ndm)%cldf = min(max(month_met(1:ndm)%cldf,0.),1.)
  month_met(1:ndm)%wind = max(month_met(1:ndm)%wind,0.)

  month_met(1:ndm)%prec = roundto(month_met(1:ndm)%prec,1)
  month_met(1:ndm)%tmin = roundto(month_met(1:ndm)%tmin,1)
  month_met(1:ndm)%tmax = roundto(month_met(1:ndm)%tmax,1)
  month_met(1:ndm)%cldf = roundto(month_met(1:ndm)%cldf,3)
  month_met(1:ndm)%wind = roundto(month_met(1:ndm)%wind,2)

  ! months with very small DTR can lead to tmin > tmax at daily timestep (Leo Lai Apr 2021)
  where (month_met%tmin > month_met%tmax)
    month_met%tmax = month_met%tmin + 0.1
  end where

  !-----------------------------------------------------------
  ! add the current monthly values on to the smoothing buffer
  ! write(0,*)cntt,t,w,t+w+1

  if (m < 12) then

    mtminbuf = eoshift(mtminbuf,1,mtmin(t+w+1))
    mtmaxbuf = eoshift(mtmaxbuf,1,mtmax(t+w+1))
    ndbuf    = eoshift(ndbuf,1,nd(t+w+1))
    cldbuf   = eoshift(cldbuf,1,cld(t+w+1))
    wndbuf   = eoshift(wndbuf,1,wnd(t+w+1))

  end if

  !-----------------------------------------------------------
  ! save month_met into dayvars module variable

  end = start + ndm - 1

  dayvars(start:end)%prec = month_met(1:ndm)%prec
  dayvars(start:end)%tmin = month_met(1:ndm)%tmin
  dayvars(start:end)%tmax = month_met(1:ndm)%tmax
  dayvars(start:end)%cldf = month_met(1:ndm)%cldf
  dayvars(start:end)%wind = month_met(1:ndm)%wind

  start = end + 1

end do monthloop ! month loop

!-----------------------------------------------------------
! Save the randomstate from last output for next year iteration
! georndst(grid) = met_out%rndst


end subroutine gwgen_new


end module gwgenmod
