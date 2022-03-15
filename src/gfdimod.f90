module gfdimod

use parametersmod, only : i4,sp

implicit none

contains

!-------------------------------------------------------------

subroutine calcGFDI(day,ndyear)

! Subroutine to calculate daily Grassland Fire Danger Index

use statevarsmod,  only : dayvars

!arguments
integer,  intent(in) :: day       ! Julian day (1 to 366)
integer,  intent(in) :: ndyear

!pointers
real(sp),    pointer :: prec_accum
real(sp),    pointer :: curing
real(sp),    pointer :: fuelmc
real(sp),    pointer :: GFDI
integer(i4), pointer :: GFDI_cat

! Local variables
real(sp) :: tmean
real(sp) :: rhum
real(sp) :: wind
real(sp) :: fuelw
integer  :: dsFeb
integer  :: Mar1_jd
integer  :: Jul31_jd

fuelw = 5.0

!-------------------------------------------------------------
! Get variables from dayvars
prec_accum => dayvars(day)%prec_accum
curing     => dayvars(day)%curing
fuelmc     => dayvars(day)%fuelmc
GFDI       => dayvars(day)%GFDI
GFDI_cat   => dayvars(day)%GFDI_cat

tmean = dayvars(day)%tmean
rhum  = dayvars(day)%rhum
wind  = dayvars(day)%wind * 3.6   ! Convert ms-1 to kmh-1

!------------------
! Find Julian day for Mar 1 and Jul 31
if (ndyear == 365) then

  Mar1_jd  = 60
  Jul31_jd = 212

else                ! Leap year

  Mar1_jd  = 61
  Jul31_jd = 213

end if

!------------------
! Calculate degree of curing for grass (based on equations given by Yui Wang)
if (day < Mar1_jd) then   ! Before 1 Mar or after 31 Jul

  prec_accum = 0.

  curing = 0.246 * day + 76.58     ! 76.58 is curing at 31 Dec of previous year

else if (day > Jul31_jd) then

  prec_accum = 0.

  curing = 0.246 * (day - Jul31_jd) + 38.7

else

  prec_accum = sum(dayvars(Mar1_jd:day)%prec)    ! Prec accumulation since 1 Mar

  curing = -0.0324 * prec_accum + 82.5407

end if

!------
! Calclate fuel moisture content
fuelmc = ((97.7 + (4.06 * rhum))) / (tmean + 6.0) - (0.00854 * rhum) + (3000.0 / curing) - 30.0

!------
! Calculate Grassland fire danger index
if (fuelmc < 18.8) then

  GFDI = 3.35 * fuelw * exp((-0.0897 * fuelmc) + (0.0403 * wind))

else if (fuelmc < 30) then

  GFDI = 0.299 * fuelw * exp (-1.686 + (0.0403 * wind)) * (30 - fuelmc)

else

  GFDI = 0.

end if

!------
! Calculate GFDI category
if (GFDI <= 3) then

  GFDI_cat = 1

else if (GFDI > 3 .and. GFDI <= 8) then

  GFDI_cat = 2

else if (GFDI > 8 .and. GFDI <= 21) then

  GFDI_cat = 3

else if (GFDI > 21 .and. GFDI <= 51) then

  GFDI_cat = 4

else if (GFDI > 51) then

  GFDI_cat = 5

end if


end subroutine calcGFDI

! -------------------------------------------------------------

end module gfdimod
