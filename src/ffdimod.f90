module ffdimod

use parametersmod, only : i4,sp

implicit none


! Module variables that are saved here for calculations in the next time-step
integer(i4) :: dryd       ! Number of dry days preceding current day / day since rain
real(sp)    :: KBDI_0     ! Previous day KBDI
real(sp)    :: prec_acc   ! Rainfall accumulation of the current wet period (mm)


contains

!---------------------------------------------------------------------

subroutine calcFFDI(year,day,s,e)

! Subroutine to calculate McArthur's fire-danger meters
! Equations from Noble et al. (1980) McArthur's fire danger meters expressed as equations https://doi.org/10.1111/j.1442-9993.1980.tb01243.x
! Coded by Leo Lai as part of ALVAR model extension (Jul 2021)

use statevarsmod, only : monvars,dayvars

integer(i4), intent(in) :: year    ! Year
integer(i4), intent(in) :: day     ! Day of year
integer(i4), intent(in) :: s       ! Index of first month (Jan) of current year
integer(i4), intent(in) :: e       ! Index of last month (Dec) of current year

! Pointers
real(sp),    pointer :: DF          ! Drought factor
real(sp),    pointer :: KBDI        ! Keetch-Byram Drought Index (mm equivalent)
real(sp),    pointer :: FFDI        ! Forest fire danger index Mark 5 meter
integer(i4), pointer :: FFDI_cat    ! Forest fire danger index category

! Local variables
real(sp) :: tmean         ! 24 hour mean temperature (degC)
real(sp) :: tmax          ! Daily maximum temperature (degC)
real(sp) :: rhum          ! Relative humidity (%)
real(sp) :: wind          ! 10m windspeed (m s-1)
real(sp) :: prec          ! Daily precipitation (mm)
real(sp) :: aprec         ! Annual precipitation (mm)
real(sp) :: sloperad
real(sp) :: Qsi

!----------
! Point to state variables

DF       => dayvars(day)%DF
KBDI     => dayvars(day)%KBDI
FFDI     => dayvars(day)%FFDI
FFDI_cat => dayvars(day)%FFDI_cat


! sloperad = topovars%sloperad   ! Disabled for now because no topography data

tmean = dayvars(day)%tmean !* sloperad
tmax  = dayvars(day)%tmax
rhum  = dayvars(day)%rhum
wind  = dayvars(day)%wind
prec  = dayvars(day)%prec

aprec = sum(monvars%pre(s:e))

!----------
! Initialize the module variables to avoid numerical errors
if (year == 1 .and. day == 1) then

  KBDI_0   = 0.
  dryd     = 0.
  prec_acc = 0.

end if

!----------
! Calculation of the Keetch_Bryam Drought Index (KBDI)
! Equations taken from Janis et al. (2002 Intl J. of Wildland Fire, 2002, 11, 281-289 http://climate.geog.udel.edu/~climate/publication_html/Pdf/JJF_IJWF_02.pdf

if (prec > 0) then      ! if wet day

  dryd     = 0                 ! reset consecutive dry day
  prec_acc = prec_acc + prec   ! accumulate rainfall for wet day (mm)

else

  dryd     = dryd + 1          ! add one to consecutive dry day count
  prec_acc = 0.                ! reset precipitation accumulation (mm)

end if

!----------

DF = (800. - KBDI_0) * (0.968 * exp(0.0875 * tmax + 1.5552) - 8.30) * 1.e-3  &
      / (1. + 10.88 * exp(-1.74e-3 * aprec))

!----------

if (prec == 0. .and. tmax <= 6.78) then

  KBDI = KBDI_0

else if  (prec == 0. .and. tmax > 6.78) then

  KBDI = KBDI_0 + DF

else if (prec > 0 .and. prec_acc < 5.1) then

  KBDI = KBDI_0 + DF

else if (prec > 0 .and. prec_acc >= 5.1) then

  DF = (800. - KBDI_0 + 3.937 * prec_acc) * (0.968 * exp(0.0875 * tmax + 1.5552) - 8.30) * 1.e-3  &
        / (1. + 10.88 * exp(-1.74e-3 * aprec))

  KBDI = KBDI_0 - 3.937 * prec_acc + DF

  prec_acc = 0.                  ! Reset rainfall accumulation after subtracting once >> Eq. 2d text explanation

end if

!----------

KBDI = max(0., KBDI)

KBDI_0 = KBDI        ! Save current KBDI for calculation of next time-step

!----------
! Here I implemented the KBDI equations from Yui Wang's code (equation source?)
! Seems less complicated and similar to what I have above

Qsi = KBDI_0 - prec

KBDI = Qsi + (((203.2 - Qsi) * (0.968 * exp((0.0875 * tmean) + 1.5552) - 8.3) * 1.* 0.01) / (1. + (10.88 * exp(-0.001736 * aprec))))

!----------
! Calculation of Forest Fire Danger Index
! Equations from Noble et al. (1980) McArthur's fire danger meters expressed as equations https://doi.org/10.1111/j.1442-9993.1980.tb01243.x

wind = wind * 3.6         ! Convert to km h-1

DF = 0.191 * (KBDI + 104.) * (dryd + 1) ** 1.5 / (3.52 * (dryd + 1) ** 1.5 + prec - 1.)

FFDI = 2.0 * exp(-4.50 + 0.987 * log(DF) - 0.0345 * rhum + 0.0338 * tmean + 0.0234 * wind)

! FFDI = 1.25 * DF * exp((tmean - rhum) / 30.0 + 0.0234 * wind)       ! Simplified euqation from Noble et al.

! if (ForFireMk5 > 120.) ForFireMk5 = 120.

!----------
! Calculate FFDI categories
if (FFDI <= 5.) then

  FFDI_cat = 1

else if (FFDI <= 12.) then

  FFDI_cat = 2

else if (FFDI <= 25.) then

  FFDI_cat = 3

else if (FFDI <= 50.) then

  FFDI_cat = 4

else

  FFDI_cat = 5

end if


end subroutine calcFFDI

!---------------------------------------------------------------------

end module ffdimod
