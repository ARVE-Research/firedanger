module diurnaltempmod

! Copied from arve-dgvm https://github.com/jedokaplan/arve-dgvm_original.git
! Edited by Leo O Lai, HKU, 2021

use parametersmod, only : i2,i4,sp,dp,pi,d2r,daysec
use daylengthmod,  only : daylength

implicit none

public :: diurnaltemp

real(dp) :: sunrise     !time of (relative to solar noon (rad))
real(dp) :: sundown        !time of (relative to solar noon (rad))
real(dp) :: peakt       !time of peak temperature
real(dp) :: t0                !temperature at the sunset hour
real(dp) :: r                !difference between max temp and temp at sunset (C)
real(dp) :: a           !amplitude of temp change (max - min) in that day (C)
real(dp) :: b                !eqn 10 from Cesaraccio paper (ref below)
real(dp) :: hni         !length of the night (hrs)
real(dp) :: tfollow     !the next days minimum temperature (assigned based on grid or point mode)

!-------
contains

! Formulas used herein are adapted from Cesaraccio, C. Spano, D., Pierpaolo, D., Snyder R. Int. J. Biometerol.
! 2001 (45) 161-169
!NOTE there is an error in the paper,  7a should have a sin after the alpha.
!--

!-----------------------------------------------------------

subroutine diurnaltemp(day,ndyear,lat)

use statevarsmod,  only : dayvars

!arguments
integer,  intent(in) :: day       ! Julian day (1 to 366)
integer,  intent(in) :: ndyear
real(dp), intent(in) :: lat

!pointers
real(sp), pointer :: tmin                       !min temperature (C)
real(sp), pointer :: tmax                       !max temperature (C)
real(sp), pointer :: tday                       !daytime temperature (C)
real(sp), pointer :: tnight                     !nighttime temperature (C)
real(sp), pointer :: tmean                      !mean temperature (C)

!parameter
real(dp), parameter :: tfpk = 1./6.                 !delay between solar noon and peak temperature (fraction)

!local variables
real(dp), dimension(2) :: dayl
real(dp) :: hdl,hdlnext             !half day length (sec) (this day and next)
real(dp) :: ti1                     !midnight till sunup
real(dp) :: ti                      !sundown till midnight
real(dp) :: tam                     !daytime till noon
real(dp) :: tpm                     !daytime post noon till sundown
real(dp) :: morn                    !sunrise - peakt
real(dp) :: sunrise_next            !relative to solar noon (rad)

integer(i4) :: dayhour
integer(i4) :: nighthour
integer(i4) :: dayhour_n

!point pointers to global data
tmin     => dayvars(day)%tmin
tmax     => dayvars(day)%tmax
tmean    => dayvars(day)%tmean
tday     => dayvars(day)%tday
tnight   => dayvars(day)%tnight

if (day == ndyear) then
  tfollow = dayvars(ndyear)%tmin
else
  tfollow = dayvars(day+1)%tmin
end if

!------
!initial assignments



! Calculate current day and next day length
! Subroutine calc dayl of next day, so -1 will give current day
call daylength(day-1, lat, dayl(1))
call daylength(day, lat, dayl(2))

! For seperating daytime and nighttime data points: (Leo Lai May 2019)
! Define daylength of polar night to be 1 hour
! Define daylength of polar day to be 23 hours
if (dayl(1) >= daysec - 3600.) dayl(1) = daysec - 3600.
if (dayl(2) >= daysec - 3600.) dayl(2) = daysec - 3600.

if (dayl(1) < 3600.) dayl(1) = 3600.
if (dayl(2) < 3600.) dayl(2) = 3600.

!these are prep for other calculations below
!tfollow = minimum temperature of the following day
t0 = tmax - 0.39 * (tmax - tfollow)
a = tmax - tmin
r = tmax - t0

!find sunrise and sundown
hdl = 0.5 * dayl(1) / 3600.  !hrs
sunrise = 12. - hdl  !hrs, fixes time from noon till sun up by using half the day length
sundown = 12. + hdl

!find next days sunrise
hdlnext = 0.5 * dayl(2) / 3600.
sunrise_next = 12 - hdlnext

!this gets the time of peak temperature
peakt = 12. + 2. * hdl * tfpk

! if (dayl(1) <= 3600.) then
!   !FLAG this could be a bit of an ugly way of doing this... it could bias the weather to be mostly
!   !hot since the day (which is the longest) will be set to the max temp of the day. The main reason for
!   !this problem is that for gridded data we do not have the mean daily temp calced out in weathergen
!   !like we do for the max and min. Will think of better ways to do this -JM Oct 29 08
!       !FLAG check on this again after we fix weathergen JM Jan 10 09
!   tday = tmax
!   tnight = tmin

! if (dayl(1) == dayl(2) .AND. dayl(1) /= 0.) print *, dayl

if (dayl(1) < daysec) then !daylength is less than 24 hours

        !has a night and a day, calculate night first

          !find the length of the night
          hni = (43200. - 0.5 * dayl(1) + 43200. - 0.5 * dayl(2))  / 3600. !hrs

          b = (tfollow - t0) / sqrt(hni) !eqn 10 from Cesaraccio paper

          ti  = t0 * sundown                                            !sundown
          ti1 = t0 * (sundown + hni) + 2./3. * b * (hni**(3./2.))   !sunrise (next morn)

          tnight = (ti1 - ti) / hni

          if (dayl(1) > 0.) then  !regular night and day

                    !morning integral (ti is at sunrise, ti1 is at temperature peak time)
                    morn = sunrise - peakt

                    ti  = (tmin * sunrise) + (1. / pi) * (2. * a * morn)
                    ti1 = (tmin * peakt)   + (1. / pi) * (2. * a * morn * cos(pi/2. * (sunrise - peakt) / morn))

                    tam = (ti1 - ti) / (-morn)

                    !afternoon integral (ti is at temperature peak time, ti1 is at sundown)
                    ti  = t0 * peakt   - (1. / pi) * 8. * r * cos(pi / 8. * (-4.))
                    ti1 = t0 * sundown - (1. / pi) * 8. * r * cos(pi / 8. * (peakt - sundown - 4.))

                    tpm = (ti1 - ti) / (sundown - peakt)

                    tday = (tam + tpm) / 2.

          ! else
          !
          !       tday = tnight      !only night, day = night

          end if

end if

! COMMENTED OUT BY Leo Lai (May 2019) because polar day and polar night is redefined above
! else !no night, only day
!
!           !morning integral (ti is at sunrise, ti1 is at temperature peak time)
!           morn = sunrise - peakt
!
!           ti  = tmin * sunrise + 1. / pi * (2. * a * morn)
!           ti1 = tmin * peakt   + 1. / pi * (2. * a * morn * cos(pi / 2. * (sunrise - peakt) / morn))
!
!           tam = (ti1 - ti) / (-morn)
!
!           !afternoon integral (t10 is at temperature peak time, ti1 is at sundown)
!           ti  = t0 * peakt   - 1. / pi * 8. * r * cos(pi / 8. * (-4.))
!           ti1 = t0 * sundown - 1. / pi * 8. * r * cos(pi / 8. * (peakt - sundown - 4.))
!
!
!           tpm = (ti1 - ti) / (sundown - peakt)
!
!           tday = (tam + tpm) / 2.
!
!           tnight = tday
!
! end if

! Calculate mean temperature from tday and tnight to stored into 'dayvars'
tmean = tday * (dayl(1) / daysec) + tnight * ((daysec - dayl(1)) / daysec)

! Claculate and save daylength variables to dayvars
dayvars(day)%dayl = dayl(1) / 3600.      ! second to hour
dayvars(day)%dayl_n = dayl(2) / 3600.    ! second to hour

! Find number of hours with sunlight on current and next day (closest integer for indexing in sub-daily processes)
dayhour = ceiling(dayl(1) / 3600.)
dayhour_n = ceiling(dayl(2) / 3600.)

! Find current day sunrise and sunset hour
! Variable adapted for indexing on a 24 hourly array such that
!     sunrise = the first hour at which sunlight is observed (inclusive) --> same for sunrise_n
!     sunset = the first hour at which there is NO sunlight (i.e. if sunset = 18:00, meaning light is last observed in 17:00 hour interval)
! NOTE: for indexing, day == (sunrise:sunset-1) and night == (sunset:24) + (1:sunrise_n-1)

if (mod(dayhour,2) /= 0) then

  dayvars(day)%sunrise = 12 - (dayhour-1) / 2
  dayvars(day)%sunset = 12 + (dayhour-1) / 2 + 1

else

  dayvars(day)%sunrise = 12 - (dayhour/2) + 1
  dayvars(day)%sunset = 12 + (dayhour/2) + 1

end if

!---

if (mod(dayhour_n,2) /= 0) then

  dayvars(day)%sunrise_n = 12 - (dayhour_n-1) / 2

else

  dayvars(day)%sunrise_n = 12 - (dayhour_n/2) + 1

end if

!---

nighthour = (24 - dayvars(day)%sunset + 1) + (dayvars(day)%sunrise_n - 1)

!---

dayvars(day)%dayhour = dayhour
dayvars(day)%nighthour = nighthour

! if (lprint .and. grid==gprint) &
!   print*, dayhour, nighthour, dayvars(day)%sunset-dayvars(day)%sunrise

! if (dayvars(day)%sunset-dayvars(day)%sunrise /= dayhour) &
!   print *, dayhour, nighthour, dayvars(day)%sunrise, dayvars(day)%sunrise_n

end subroutine diurnaltemp

!-----------------------------------------------------------

subroutine humidity(day)

! Calculate relative humidity from dew point temperature
! From Lawrence (2005) The relationship between relative humidity and the dewpoint temperature. A. MetSoc (https://doi.org/10.1175/BAMS-86-2-225)
! Equation (11)

use statevarsmod, only : dayvars

implicit none

integer, intent(in) :: day

real(sp), pointer :: tmean
real(sp), pointer :: tdew
real(sp), pointer :: rhum

real(sp) :: tmean_K
real(sp) :: tdew_K
real(sp) :: tmean_F
real(sp) :: tdew_F

real(dp), parameter :: L  = 2257000
real(dp), parameter :: Rw = 461.5

!------
! point pointers to global data
tmean => dayvars(day)%tmean
tdew  => dayvars(day)%tdew
rhum  => dayvars(day)%rhum

tmean_K = tmean + 273.15
tdew_K  = tdew + 273.15

tmean_F = tmean * 1.8 + 32.
tdew_F  = tdew * 1.8 + 32.

!------

rhum = 100. * exp(((1. - (tmean_K / tdew_K)) * (L / Rw)) / tmean_K)

! rhum = 100. - (25./9.) * (tmean_F - tdew_F)

if (rhum > 100.) rhum = 100.    ! When temp is lower than dew point


end subroutine humidity


end module diurnaltempmod
