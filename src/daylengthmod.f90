module daylengthmod

! Copied from arve-dgvm https://github.com/jedokaplan/arve-dgvm_original.git
! Edited by Leo O Lai, HKU, 2021

implicit none

public :: daylength

!---------------------------------

contains

!---------------------------------
subroutine daylength(day,lat,dayl_n)

!for any given latitude and day of year, calculate the daylength in seconds
!for use by ARVEpoint only! NOTE: this is set up to calculate the
!NEXT days values (i.e. d+1)
!NOTE: the solar declination value herein is not the Berger one so is not valid for
!paleosimulations. ARVE-Point is likely not to be used for that anyway. JM Oct 30 08

use parametersmod, only : dp,pi,d2r,a2s,solarc,dayspy
! use statevars, only : sv,dayl

implicit none

!arguments
integer,  intent(in)  :: day      ! Day in julian date (1 to 366)
real(dp), intent(in)  :: lat      ! Lat in degrees
real(dp), intent(out) :: dayl_n   ! Daylength of the next day (definition adapted from original arve-gdvm code)

!pointers
! real(dp), pointer :: lat        !lat in degrees
! real(dp), pointer :: delta      !Next day's solar declination (deg)

!local variables
real(dp) :: sinlat
real(dp) :: coslat
real(dp) :: rdelta      ! Solar declination (rad)
real(dp) :: u
real(dp) :: v
real(dp) :: hh
! real(dp) :: sinhh
! real(dp) :: qo
! real(dp) :: w
integer  :: d,e
real(dp) :: theta               ! Earth orbit seasonal angle (radians)

!parameters
real(dp), parameter :: beta   =    0.17   ! Global average shortwave albedo fraction
real(dp), parameter :: ac     =    0.25   ! Angstrom parameters for sunlight transmission through clouds
real(dp), parameter :: ad     =    0.50

!point pointers
! lat => sv(1)%lat
! delta => sv(1)%sdecl_n

!------
!calculations begin

!correct for the next timestep adjustment
if (day == 366) then
 d = 1
 e = 365
else if (day == 0) then
 d = 1
 e = 1
else
 d = day
 e = d - 1
end if

sinlat = sin(d2r * lat)
coslat = cos(d2r * lat)

!find solar declination in rad (old way)
!rdelta = -23.4 * d2r * cos(d2r * 360.0 * (real(d) + 10.0) / 365.0)

!find the Earth orbit seasonal angle in radians
 theta = 2._dp * pi * real(d) / dayspy

!Solar declination in radians (from CLM 3.0 and many other sources):
 rdelta = 0.006918_dp - 0.399912_dp * cos(theta) + 0.070257_dp * sin(theta) - &
          0.006758_dp * cos(2._dp * theta) + 0.000907_dp * sin(2._dp * theta) - &
          0.002697_dp * cos(3._dp * theta) + 0.001480_dp * sin(3._dp * theta)

u = sinlat * sin(rdelta) !Eqn 9
v = coslat * cos(rdelta) !Eqn 10

!Calculate half-day length in angular units, hh
!In Eqn (11), hh defined for u in range -v to v
!For u >= v, hh = pi (12 hours, i.e. polar day)
!For u <= -v, hh = 0 (i.e. polar night)

if (u >= v) then         !polar day
  hh = pi
else if (u <= -v) then   !polar night
  hh = 0.0
else                     !normal day and night
  hh = acos(-u / v)
end if

!calculate daylength in sec from hh
dayl_n = 86400. * hh / pi

!solar rad calculations
! qo = solarc * (1.0 + 2.0 * 0.01675 * cos(2.0 * pi * real(e) / 365.0))
!
! w = (ac + ad * sunf) * (1.0 - beta) * qo
!
! sinhh = sin(hh)
!
! !Calculate total net downward shortwave radiation for this day
! calcrad = 2.0 * w * (u * hh + v * sinhh) * a2s  !(kJ cm-2 d-1)
!
! !delta is expected elsewhere in degrees so convert.
! delta = rdelta / d2r

end subroutine daylength


end module daylengthmod
