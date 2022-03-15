module parametersmod

! Simple module defining some types and parameters

use iso_fortran_env, only : int8,int16,int32,real32,real64,output_unit

implicit none

public :: ndaymonth

integer, parameter :: i1 = int8    ! 1 byte integer
integer, parameter :: i2 = int16   ! 2 byte integer
integer, parameter :: i4 = int32   ! 4 byte integer
integer, parameter :: sp = real32  ! 4 byte real
integer, parameter :: dp = real64  ! 8 byte real

integer, parameter :: so = output_unit  ! unit number for standard output

integer, parameter :: baseyr = 1971

real(sp), parameter :: Tfreeze = 273.15 ! freezing temperature of freshwater (K)

real(dp), parameter :: dayspy = 365._dp     !number of days in one year
real(dp), parameter :: daysec = 86400._dp   !number of seconds in 24hrs
real(dp), parameter :: hr1    = 3600._dp    !one hour in seconds
real(dp), parameter :: hr23   = 82800._dp   !23 hours in seconds

real(dp), parameter :: pi     = 3.14159265358979323846_dp
real(dp), parameter :: solarc = 1360.            !Solar constant (1360 W/m2)
real(dp), parameter :: d2r    = pi / 180.        !conversion factor from degrees to radians
real(dp), parameter :: r2d    = 180. / pi        !conversion factor from degrees to radians
real(dp), parameter :: a2s    = 1. / (300. * pi) !conversion factor from solar angular to seconds

real(sp), parameter :: hsp = huge(sp)    ! largest positive 4-byte real

integer, parameter, dimension(12) :: midday    = [ 16,44,75,105,136,166,197,228,258,289,319,350 ]  !day number of mid-month day

real(sp),    parameter :: missing_sp = -9999.
integer(i2), parameter :: missing_i2 = -32768

contains

!-----------------------------------------------------------------------------

integer function ndaymonth(yr,mon)

! Function to find out the number of days in a month, considering leap years and the year given as AD

! Input: Current year and month

integer, intent(in) :: yr
integer, intent(in) :: mon

! Arrays defining standard and leap years

integer, parameter, dimension(12) :: std_year = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
integer, parameter, dimension(12) :: leapyear = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

!---------------------------------

if (mod(yr,400) == 0) then        ! If year can be divided by 400, then it's a leap year

  ndaymonth = leapyear(mon)       ! Choose amount of days from leap year array

else if (mod(yr,100) == 0) then   ! If year can't be divided by 400 but can by 100, then it's a standard year

  ndaymonth = std_year(mon)

else if (mod(yr,4) == 0) then     ! If year can't be divided by 400 or 100, but can be divided by 4, it's a leap year

  ndaymonth = leapyear(mon)

else                              ! Any other case, it's a standard year

  ndaymonth = std_year(mon)

end if

end function ndaymonth

!-----------------------------------------------------------------------------

end module parametersmod
