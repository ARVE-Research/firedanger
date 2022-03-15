module radiationmod

! Copied from LPJ-LMFire https://github.com/ARVE-Research/LPJ-LMfire.git
! Edited by Leo O Lai, HKU, 2021

!consolidated module for calculation airmass, surface downwelling shortwave, net longwave, net radiation and PET.

use parametersmod, only : i2,i4,sp,dp,pi,d2r

implicit none

!module subroutines and functions

public  :: initairmass
public  :: calcPjj
public  :: elev_corr
public  :: radpet

private :: airmass
private :: surf_sw
private :: surf_lw
private :: netrad_pet

private :: m
private :: F
private :: esat
private :: desdT

!module parameters

real(sp), parameter :: w   = 15.  !solar angular velocity (degrees hr-1)
real(sp), parameter :: rw  = d2r * w !solar angular velocity (radians hr-1)

real(sp), parameter :: m0  =  1.   !air mass at 0 degree solar zenith angle
real(sp), parameter :: m80 =  5.6  !air mass at 80 degree solar zenith angle
real(sp), parameter :: m90 = 39.7  !air mass at 90 degree solar zenith angle

real(sp), parameter :: cos80 = cos(80. * d2r)  !(degrees)

real(sp), parameter :: albedo = 0.17  !surface shortwave albedo (fraction)

!module shared variables

real(sp), dimension(3) :: c00  !air mass coefficients for solar zenith angle <=80 degrees
real(sp), dimension(3) :: c80  !air mass coefficients for solar zenith angle  >80 degrees

type airmasspars
  real(sp) :: Ratm    !relative atmospheric pressure 1=sea level
  real(sp) :: mbar    !daytime mean optical air mass (unitless, 1 at equatorial noon)
  real(sp) :: mo      !air mass at cosine zenith angle maximum
  real(sp) :: mc      !air mass at cosine zenith angle medium
  real(sp) :: ml      !air mass at cosine zenith angle bottom quarter range point
end type airmasspars

real(sp) :: lvap    !Latent heat of vaporization of water (temperature dependent) (kJ kg-1)
real(sp) :: gamma   !psychrometer constant (Pa K-1)
real(sp) :: ss      !rate of increase of saturated vapor pressure with temperature (desdT) (Pa K-1)

contains   !the following subroutines and functions

!----------------------------------------------------------------------------------------------------------------

subroutine initairmass()

!calculate parameters used in the airmass calculations

implicit none

c00(1) = 0.008307
c00(2) = (m0 - m80) * (c00(1) + 1.) * (c00(1) + cos80) / (cos80 - 1.)
c00(3) = m0 - c00(2) / (c00(1) + 1.)

c80(1) = 0.037160
c80(2) = (m90 - m80) * c80(1) * (c80(1) + cos80) / cos80
c80(3) = m90 - c80(2) / c80(1)

end subroutine initairmass

!----------------------------------------------------------------------------------------------------------------

subroutine calcPjj(day,s,e)

use statevarsmod, only : monvars,dayvars    ! modified for weathergen (Leo O Lai Apr 2021)

implicit none

!arguments
integer, intent(in) :: day
integer, intent(in) :: s
integer, intent(in) :: e

real(sp), dimension(12) :: temp  !annual time series of monthly temperature
real(sp), dimension(12) :: prec  !annual time series of monthly precipitation

real(sp), pointer :: Pjj

!local variables

integer, dimension(1) :: wm  !index position of the warmest month
integer, dimension(1) :: cm  !index position of the coldest month

real(sp) :: p_wm
real(sp) :: p_cm

!-----------------------------------------------------
! Define pointer variables to genvar and dayvars calculated by 'gwgen'
temp = monvars%tmp(s:e)
prec = monvars%pre(s:e)

Pjj  => dayvars(day)%Pjj

!-----------------------------------------------------
!calculation of the "precipitation equitability index"

wm   = maxloc(temp)  !temperature of the warmest month (should be carried as the last warmest month, or running average)
cm   = minloc(temp)  !temperature of the coldest month (should be carried as the last coldest month, or running average)
p_wm = prec(wm(1))   !total precipitation in the warmest month
p_cm = prec(cm(1))   !total precipitation in the coldest month

if (p_wm + p_cm > 0.) then
  Pjj = 2. * (p_wm - p_cm) / (p_wm + p_cm)
  Pjj = max(Pjj,0.)
else
  Pjj = 0.
end if

end subroutine calcPjj

!----------------------------------------------------------------------------------------------------------------

real(sp) function elev_corr(elevation)

implicit none

real(sp), intent(in)  :: elevation

real(sp), parameter :: z0 = 1. / 8000.

!----

elev_corr = exp(-elevation * z0)

end function elev_corr

!----------------------------------------------------------------------------------------------------------------

subroutine elev_Ratm(day,elev)

! Calculate variables for atmospheric pressure relative to sea-level
! Equations copied from ARVE-DGVM (Leo Lai, Jul 2021)

use statevarsmod, only : dayvars

implicit none

real(sp), parameter :: Pstd = 101325.     ! Standard atmospheric pressure (Pa)

integer,  intent(in) :: day
real(sp), intent(in) :: elev

real(sp), pointer :: Patm           ! Atmospheric pressure at elevation (Pa)
real(sp), pointer :: Patm30         ! Atmospheric pressure at elevation adjusted by reference 30m height (Pa)
real(sp), pointer :: Ratm           ! Relative atmospheric pressure to sea-level at elevation (fraction)
real(sp), pointer :: Ratm30         ! Relative atmospheric pressure to 30m ground at elevation (fraction)

Patm   => dayvars(day)%Patm
Patm30 => dayvars(day)%Patm30
Ratm   => dayvars(day)%Ratm
Ratm30 => dayvars(day)%Ratm30

!----------

! Elevation is used to get the atmospheric pressure (Diehl 1925 (from Lloyd and Farquhar 1994))
! Code copied from ARVE-DGVM by Leo Lai (Jul 2021)
Ratm = (1.0 - elev / 44308.) ** 5.2568         ! Multiply fraction by standard atm for pressure level in Pa

! Also calculate the reference height pressure (30m above ground) (same ref as above)
Ratm30 = (1.0 - (elev + 30._dp) / 44308.) ** 5.2568

Patm = Pstd * Ratm

Patm30 = Pstd * Ratm30

end subroutine elev_Ratm

!----------------------------------------------------------------------------------------------------------------

subroutine radpet(day,s,e,lat)

! Calculate total downwelling shortwave radiation and PET (potential evapotranspiration)

use parametersmod, only : midday
use orbitmod,      only : orbitpars,toa_insolation,orbit
use statevarsmod,  only : monvars,dayvars

implicit none

! Arguments
integer, intent(in) :: day
integer, intent(in) :: s
integer, intent(in) :: e
real(dp),intent(in) :: lat       ! Latitude (degrees)

!local variables

real(sp) :: temp      ! Daytime mean temperature (degC)
real(sp) :: prec      ! Total 24h precipitation (mm)
real(sp) :: cldf      ! Cloud cover fraction (percent)
real(sp) :: tdew      ! Dewpoint temperature (degC)
real(sp) :: tcm       ! Temperature of the coldest month (of the year) (degC)
real(sp) :: Ratm      ! Relative atmospheric pressure based on elevation (fraction)
real(sp) :: Pjj       ! Water equitability index

type(airmasspars) :: air

real(sp) :: toa_sw    ! Top of the atmosphere downwelling shortwave rad (kJ m-2 d-1)
real(sp) :: delta     ! Solar declination (degrees)
real(sp) :: pet0      ! Previous value for PET (mm d-1)
real(sp) :: direct    ! Direct beam surface downwelling shortwave (kJ m-2 d-1)
real(sp) :: diffuse   ! Diffuse surface downwelling shortwave (kJ m-2 d-1)
real(sp) :: lw_rad    ! Net longwave radiation (kJ m-2 d-1)

real(sp) :: dayl      ! Daylength (h)
real(sp) :: sw_rad    ! Total surface downwelling shortwave (kJ m-2 d-1)
real(sp) :: pet       ! Daily potential evapotranspiraton (mm) (per day)

real(sp) :: sloperad  ! Ratio of direct raditaion on sloped and flat surface

!counters

integer :: i

!----------------------------------------------------------------------------------
! Get met variables for the day from module variable 'dayvars'
temp = dayvars(day)%tmean
prec = dayvars(day)%prec
cldf = dayvars(day)%cldf
Ratm = dayvars(day)%Ratm
Pjj  = dayvars(day)%Pjj
tcm  = minval(monvars%tmp(s:e))    ! genvar structure with +/- 4 months buffer

call toa_insolation(orbit,day,lat,toa_sw,dayl,delta)

call airmass(lat,delta,dayl,Ratm,air)

call surf_lw(day,temp,dayvars(day)%tmin,cldf,dayl,lw_rad,tdew)

!------

pet  = 0.
pet0 = 0.

i = 1

do ! Because of the weak dependence of surface shortwave on PET, we equilibrate PET and surf_sw

  call surf_sw(Pjj,Ratm,toa_sw,cldf,dayl,air,prec,tcm,pet,direct,diffuse)

  sw_rad = direct + diffuse

  call netrad_pet(sw_rad,lw_rad,pet)

  if(abs(pet - pet0) < 0.01) exit

  if (i > 100) then
    write(0,*) 'No good stable solution for sw_rad and pet after 100 iterations'
    exit
  end if

  pet0 = pet

  i = i + 1

end do

!write(stdout,'(7f10.3)')temp,prec,cldf,toa_sw,sw_rad,lw_rad,pet

! Calculate the slope and aspect effect on direct radiation budget (sloperad = slope/flat radiation ratio)
! call radslope(grid,day,dayl,delta,lat,toa_sw,direct,diffuse,sloperad)
!
! topovars(grid)%sloperad = sloperad

dayvars(day)%tdew     = tdew
dayvars(day)%dsol     = delta
dayvars(day)%dayl     = dayl
dayvars(day)%srad     = sw_rad
dayvars(day)%srad_dir = direct
dayvars(day)%srad_dif = diffuse
dayvars(day)%lrad     = lw_rad
dayvars(day)%dpet     = pet

! Calculate tdew after getting final stable solution for pet and sw_rad
call calctdew(day,s,e)

end subroutine radpet

!----------------------------------------------------------------------------------------------------------------

subroutine airmass(lat,delta,dayl,Ratm,air)

!This code is based on the paper:
!X. Yin (1997) Optical air mass: Daily integration and its applications, Meteorol. Atmos. Phys. 63, 227-233
!Jed Kaplan, EPFL, 2008

implicit none

!arguments

real(dp), intent(in) :: lat    !latitude (degrees)
real(sp), intent(in) :: delta  !solar declination (degrees)
real(sp), intent(in) :: dayl   !day length (hours)
real(sp), intent(in) :: Ratm   !relative atmospheric pressure

type(airmasspars), intent(out) :: air

!local variables

real(sp) :: mbar    !daytime mean optical air mass (unitless, 1 at equatorial noon)
real(sp) :: mo      !air mass at cosine zenith angle maximum
real(sp) :: mc      !air mass at cosine zenith angle medium
real(sp) :: ml      !air mass at cosine zenith angle bottom quarter range point

real(sp) :: rlat      !latitude (radians)
real(sp) :: rdelta    !solar declination (radians)

real(sp) :: t1        !number of hours between sunrise/sunset and solar noon (hr)
real(sp) :: t80       !solar hour corresponding to the 80 degree zenith angle

! real(sp) :: t    !solar hour (hr)       ! Wunused variable

real(sp) :: Z    !solar zenith angle (degrees)
real(sp) :: Zn   !lesser of solar zenith angle at sunset or at midnight (degrees)
real(sp) :: Z0   !zenith angle at solar noon (degrees)
real(sp) :: cosZ !cosine solar zenith angle (fraction), used in calculation of instantaneous air mass

! real(sp) :: l     ! Wunused variable

! integer :: steps  !integer number of time steps   ! Wunused variable
! integer :: i      !counter                        ! Wunused variable

real(sp) :: sinlat
real(sp) :: coslat
real(sp) :: sindel
real(sp) :: cosdel

real(sp)               :: a   !values in equation 2.6b
real(sp)               :: b
real(sp), dimension(3) :: c

real(sp) :: tmp1
real(sp) :: tmp2
real(sp) :: tmp3

real(sp) :: tinv

real(sp) :: rZ0
real(sp) :: rZn

real(sp), parameter :: mindayl = 2. * tiny(0._sp)

!-------------------------------------
!calculate daily mean air mass (mbar)

call initairmass()

if (dayl == 0.) then

  mbar = m90
  mc   = m90
  ml   = m90
  mo   = m90

else

  !basic setup

  rlat   = d2r * lat
  rdelta = d2r * delta

  sinlat = sin(rlat)
  sindel = sin(rdelta)
  coslat = cos(rlat)
  cosdel = cos(rdelta)

  !------

  ! Eqn. 2.5 -- commented out because of floating point invalid in rare cases,
  ! plus we already have the day length calculated elsewhere (JOK 06.2016)

  ! if (abs(lat - delta) < 90. .and. abs(lat + delta) >= 90.) then
  !  t1 = 12.
  ! else
  !  t1 = (12. / pi) * acos(-tan(rlat) * tan(rdelta))
  ! end if

  if (dayl > mindayl) then
    t1 = 0.5 * dayl
  else
    t1 = dayl
  end if

  !Eqn. 2.9
  if (abs(lat + delta) >= 90.) then
    Zn = acos(sinlat * sindel - coslat * cosdel) / d2r
  else
    Zn = 90.
  end if

  !Eqn. 2.10
  if (abs(lat - delta) >= 90.) then
    Z0 = 90.
  else
    Z0 = lat - delta
  end if

  rZ0 = Z0 * d2r  !convert to radians
  rZn = Zn * d2r

  !--------------------------

  b = coslat * cosdel

  if (t1 == 0.) then

    mbar = m90

  else if (abs(Zn) <= 80.) then

    tinv = 1. / t1

    c = c00
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else if (abs(Z0) >= 80.) then

    tinv = 1. / t1

    c = c80
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else

    t80 = 1. / w * acos((cos80 - sinlat * sindel) / (coslat * cosdel)) / d2r  !Eqn. 2.8

    c = c00
    a = c(1) + sinlat * sindel

    !write(*,*)'crash',t80,a,b,c

    tmp1 = F(t80,a,b,c)

    c = c80
    a = c(1) + sinlat * sindel
    tmp2 = F(t1,a,b,c)

    c = c80
    a = c(1) + sinlat * sindel
    tmp3 = F(t80,a,b,c)

    tinv = 1. / t1


    mbar = tinv * (tmp1 + tmp2 - tmp3)

  end if

  !--------------------------
  !calculate instantaneous air mass at max, mid, and bottom quarter solar zenith angle (m0, mc, ml)

  Z = Z0

  cosZ = cos(Z * d2r)

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if

  mo = m(cosZ,c)

  !--

  Z = (Z0 + Zn) / 2.

  cosz = (cos(rZ0) + cos(rZn)) / 2.

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if

  mc = m(cosZ,c)

  !--

  Z = (Z0 + 3. * Zn) / 4.

  cosz = (cos(rZ0) + 3. * cos(rZn)) / 4.

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if

  ml = m(cosZ,c)

end if

!-------------------------------------
!correct calculated air mass for elevation

air%mbar = Ratm * mbar
air%mo   = Ratm * mo
air%mc   = Ratm * mc
air%ml   = Ratm * ml


end subroutine airmass

!----------------------------------------------------------------------------------------------------------------

subroutine surf_sw(Pjj,Ratm,r0,cldf,dayl,air,prec,tcm,pet,direct,diffuse)

!This code is based on the paper:
!X. Yin (1998) Temporally-aggregated atmospheric optical properties as a function of common climatic information:
!Systems development and application, Meteorol. Atmos. Phys. 68, 99-113
!Jed Kaplan, EPFL, 2008

implicit none

!arguments

real(sp), intent(in)  :: Pjj
real(sp), intent(in)  :: Ratm
real(sp), intent(in)  :: r0    !top-of-atmospere insolation (kJ m-2 d-1)
real(sp), intent(in)  :: cldf  !bright sunshine duration fraction, n/N (percent)
real(sp), intent(in)  :: dayl  !daylength (hr)

type(airmasspars), intent(in) :: air

real(sp), intent(in)  :: prec  !precipitation mm/day
real(sp), intent(in)  :: tcm   !temperature of the coldest month (used as tropics indicator)
real(sp), intent(in)  :: pet   !potential evapotranspiration mm/day

real(sp), intent(out) :: direct  !direct-beam downwelling shortwave (kJ m-2 d-1)
real(sp), intent(out) :: diffuse !diffuse downwelling shortwave (kJ m-2 d-1)

!local variables

real(sp) :: mbar    !daytime mean optical air mass (unitless, 1 at equatorial noon)
real(sp) :: mo      !air mass at cosine zenith angle maximum
real(sp) :: mc      !air mass at cosine zenith angle medium
real(sp) :: ml      !air mass at cosine zenith angle bottom quarter range point

real(sp) :: tau   !direct insolation atmospheric turbidity factor
real(sp) :: zeta0 !diffuse insolation atmospheric turbidity factor
real(sp) :: x     !tropics indicator (tropical = 1, else 0)
real(sp) :: fm    !atmospheric transmittance function

! real(sp) :: j2w
! real(sp) :: fdif
! real(sp) :: stmp
real(sp) :: sunf   !bright sunshine duration fraction, n/N (fraction)

!-----------------------------------
!parameters

real(sp), parameter :: kp  = 0.500 !links absorption coeff. to trans. coeff.
real(sp), parameter :: kag = 3.300
real(sp), parameter :: kan = 2.320
real(sp), parameter :: kn  = 0.686 !cloud parameter

!----------------------------------------------------------------------------

mbar = air%mbar
mo   = air%mo
mc   = air%mc
ml   = air%ml

!------

sunf = 1. - cldf

if (tcm < 10.) then
  x = 0.
else if (tcm > 20.) then
  x = 1.
else
  x = sin(pi / 2. * (tcm / 10. - 1.))
end if

!Yin Eqn. 4.1
tau = exp(-0.115 * Ratm * ((2.15 - 0.713 * x + exp(-6.74 / (prec + 1.))) * exp(0.0971 * pet) - 0.650 * (1. - x) * Pjj))

fm = 0.01452 * (mbar + ml) * exp(1.403 * tau) - 0.1528 * mo + mc + 0.48700 * (mc - ml) + 0.2323   !Eqn. 2.4 2nd term

direct = sunf * tau**kp * r0 * tau**fm   !Eqn. 2.4

!Yin Eqn. 4.2
zeta0 = 0.503 * exp(-1.20 * Ratm * exp(-0.633 / (prec + 1.) - 0.226 * pet)) * kag**albedo * kan**(1. - sunf) * &
        (1. - kn * (1. - sunf))

diffuse = zeta0 * kag**albedo * kan**(1. - sunf) * (1 - kn * (1. - sunf)) * (tau**kp * r0 - direct)   !Eqn. 2.5

!write(stdout,'(a,6f12.4)')'shortwave',r0,Ratm,prec,Pjj,sunf,direct+diffuse

end subroutine surf_sw

!----------------------------------------------------------------------------------------------------------------

subroutine surf_lw(day,temp,tmin,cldf,dayl,lw_rad,tdew)

!This code is based on the paper:
!A. Haxeltine and Prentice, I.C., BIOME3..., Glob. Biogeochem. Cycles, 10, 693-709
!With a new calculations of:
!downwelling longwave   (Josey et al., 2003. J. Geophys. Res., 108(C4), 3108, doi:10.1029/2002JC001418)
!dEsat/dT               (Oleson et al., 2004, CLM 3.0 technical note)
!lvap                   (Henderson-Sellers, 1984. Quart. J. R. Met. Soc., 110, 1186-1190)
!Other references:
!Linacre (1968) Agr. Meteorol., 5, 49-63
!Prentice et al. (1993) Ecological Modelling, 65, 51-70.
!Jed Kaplan, EPFL, 2008, 2011

use parametersmod, only : Tfreeze
use statevarsmod,    only : dayvars

implicit none

!arguments

integer,  intent(in)  :: day     !day of year
real(sp), intent(in)  :: temp    !surface air (2m) temperature (C)
real(sp), intent(in)  :: tmin    !surface air (2m) temperature (C)
real(sp), intent(in)  :: cldf    !cloud cover fraction
real(sp), intent(in)  :: dayl    !daylength (h)

real(sp), intent(out) :: lw_rad  !daytime net longwave radiation (kJ m-2 d-1)
real(sp), intent(out) :: tdew    !dew point temperature (based on input temperature) (C)

!parameters

real(sp), parameter :: sb = 5.6704e-8  !Stefan-Bolzmann constant (W m-2 K-4)
real(sp), parameter :: e  = 0.98     !emissivity ()
real(sp), parameter :: al = 0.045    !longwave reflectivity (lw albedo), Josey et al., pg 5-9

real(sp), parameter :: a  =  10.77   !parameters in Josey et al.
real(sp), parameter :: b  =   2.34
real(sp), parameter :: c  = -18.44

real(sp), parameter :: cs = 1.5 !shape parameter for the curve relating fractional cloud cover to fractional sunshine duration

!local variables

real(sp) :: Tk     !surface air temperature (K)
real(sp) :: Ts     !ground surface temperature (K)
real(sp) :: TdewK  !dewpoint temperature (K)
real(sp) :: D      !dew point depression (K)
real(sp) :: es     !saturation vapor pressure

real(sp) :: f      !Linacre parameter (function of sunshine fraction)

real(sp) :: Ql     !net longwave radiation (W m-2)
real(sp) :: Ql_up  !upwelling longwave radiation (W m-2)
real(sp) :: Ql_dn  !downwelling longwave radiation (W m-2)

real(sp) :: sunf   !bright sunshine duration fraction, n/N (fraction)

integer :: iter

!-------------------------------------------------

sunf = 1. - cldf

Tk = temp + Tfreeze

!calculate gamma, lvap

gamma = 65.05 + temp * 0.064  !psychrometer constant

lvap = 0.001 * 1.91846e6 * (Tk / (Tk - 33.91))**2  !(kJ kg-1) Eqn. from Henderson-Sellers (1984)

ss = desdT(Tk)

f = 0.2 + 0.8 * sunf  !Linacre Eqn. 7

!-------------------------------------------------
!calculate longwave radiation

Ts = Tk !approximation that mean daily surface temperature equals air temp.

!black body upwelling longwave (W m-2)  !various sources e.g., Oleson et al.

Ql_up = e * sb * Ts**4

!--
!Josey formulation for downwelling longwave

!To estimate dewpoint temperature we use the day's minimum temperature
!this makes the asumption that there is a close correlation between Tmin and dewpoint
!see, e.g., Glassy & Running, Ecological Applications, 1994
es = 0.01 * esat(tmin+Tfreeze) !saturation vapor pressure (mbar)
! es = 0.01 * esat(Tk) !saturation vapor pressure (mbar)

TdewK = 34.07 + 4157. / log(2.1718e8 / es)  !Josey et al., Eqn. 10

!----

D = TdewK - Tk

Ql_dn = sb * (Tk + a*cldf**2 + b*cldf + c + 0.84 * (D + 4.01))**4  !downwelling longwave (W m-2) Josey et al. Eqn. 14,J2

Ql = Ql_up - (1. - al) * Ql_dn   !Josey et al., Eqn 1

!----

lw_rad = 0.001 * 3600. * dayl * Ql  !daytime net longwave (kJ m-2 d-1)

tdew = TdewK - Tfreeze

!write(stdout,*)'longwave',cldf,Tdewk

end subroutine surf_lw

!----------------------------------------------------------------------------------------------------------------

subroutine netrad_pet(sw_rad,lw_rad,pet)

implicit none

real(sp), intent(in)  :: sw_rad   !downwelling shortwave radiation (kJ m-2 d-1)
real(sp), intent(in)  :: lw_rad   !net longwave radiation (kJ m-2 d-1)
real(sp), intent(out) :: pet      !potential evapotranspiration (mm d-1)

!local variable

real(sp) :: netrad !net radiation (kJ m-2 d-1)

!----

netrad = (1. - albedo) * sw_rad - lw_rad             !(kJ m-2 d-1)

pet = max((ss / (ss + gamma)) * netrad / lvap, 0.)   !(mm d-1)

end subroutine netrad_pet

!----------------------------------------------------------------------------------------------------------------

subroutine calcVPD(day)

! Subroutine to calculate vapor pressure deficit (Leo Lai Apr 2021)
! Reference: Thornton et al. (2000) doi:10.1016/S0168-1923(00)00170-2

use parametersmod, only : Tfreeze
use statevarsmod,    only : dayvars

implicit none

integer,  intent(in)  :: day     !day of year

real(sp), pointer :: vpd

real(sp) :: tmin
real(sp) :: tday

tmin = dayvars(day)%tmin
tday = dayvars(day)%tday
vpd  => dayvars(day)%vpd

vpd = esat(tday+Tfreeze) - esat(tmin+Tfreeze)

! Adjust for small numerical difference between tday and tmin which leads to negative VPD
if (vpd < 0.) vpd = 0.

end subroutine calcVPD

!----------------------------------------------------------------------------------------------------------------

function esat(temp)

  !Function to calculate saturation vapor pressure in water and ice
  !From CLM formulation, table 5.2, after Flatau et al. 1992

  use parametersmod, only : dp,Tfreeze

  implicit none

  real(sp) :: esat  !saturation vapor pressure (Pa)
  real(sp), intent(in) :: temp !temperature in K

  real(sp), dimension(9) :: al !coefficients for liquid water
  real(sp), dimension(9) :: ai !coefficients for ice

  real(sp), dimension(0:8) :: a !coefficients

  real(sp) :: T

  integer :: i

  al(1) = 6.11213476
  al(2) = 4.44007856e-1
  al(3) = 1.43064234e-2
  al(4) = 2.64461437e-4
  al(5) = 3.05903558e-6
  al(6) = 1.96237241e-8
  al(7) = 8.92344772e-11
  al(8) =-3.73208410e-13
  al(9) = 2.09339997e-16

  ai(1) = 6.11123516
  ai(2) = 5.03109514e-1
  ai(3) = 1.88369801e-2
  ai(4) = 4.20547422e-4
  ai(5) = 6.14396778e-6
  ai(6) = 6.02780717e-8
  ai(7) = 3.87940929e-10
  ai(8) = 1.49436277e-12
  ai(9) = 2.62655803e-15

  if (temp <= Tfreeze) then   !these coefficients are for temperature values in Celcius
    a(0:8) = ai
  else
    a(0:8) = al
  end if

  T = temp - Tfreeze

  if (abs(T) < 1e-2) T = 0.   ! Avoid underflow at exponential operation (Leo O Lai, Jun 2021)

  esat = a(0)

  do i = 1,8
    esat = esat + a(i) * T**i
  end do

  esat = 100. * esat

end function esat

!----------------------------------------------------------------------------------------------------------------

function desdT(temp)

  !Function to calculate the first derivative of saturation vapor pressure in water and ice vs. temperature
  !From CLM formulation, table 5.3, after Flatau et al. 1992

  use parametersmod, only : dp,Tfreeze

  implicit none

  real(sp) :: desdT    !derivative of saturation vapor pressure
  real(sp), intent(in) :: temp !temperature in K

  real(sp), dimension(9) :: bl !coefficients for liquid water
  real(sp), dimension(9) :: bi !coefficients for ice

  real(sp), dimension(0:8) :: b !coefficients

  ! real(sp) :: tmp

  real(sp) :: T

  integer :: i

  bl(1) = 4.44017302e-1
  bl(2) = 2.86064092e-2
  bl(3) = 7.94683137e-4
  bl(4) = 1.21211669e-5
  bl(5) = 1.03354611e-7
  bl(6) = 4.04125005e-10
  bl(7) =-7.88037859e-13
  bl(8) =-1.14596802e-14
  bl(9) = 3.81294516e-17

  bi(1) = 5.03277922e-1
  bi(2) = 3.77289173e-2
  bi(3) = 1.26801703e-3
  bi(4) = 2.49468427e-5
  bi(5) = 3.13703411e-7
  bi(6) = 2.57180651e-9
  bi(7) = 1.32268878e-11
  bi(8) = 3.94116744e-14
  bi(9) = 4.98070196e-17

  if (temp <= Tfreeze) then
    b(0:8) = bi
  else
    b(0:8) = bl
  end if

  T = temp - Tfreeze  !these coefficients are for temperature values in Celcius

  if (abs(T) < 1e-2) T = 0.   ! Avoid underflow at exponential operation (Leo O Lai, Jun 2021)

  desdT = b(0)

  do i = 1,8
    desdT = desdT + b(i) * T**i
  end do

  desdT = 100. * desdT

end function desdT

!----------------------------------------------------------------------------------------------------------------

real(sp) function m(cosZ,c)

!Instantaneous air mass m, equation 2.1 in Yin, 1997

implicit none

real(sp),               intent(in) :: cosZ
real(sp), dimension(:), intent(in) :: c

m = c(2) / (c(1) + cosZ) + c(3)

end function m

!----------------------------------------------------------------------------------------------------------------

real(sp) function F(t1,a,b,c)

!integral air mass function F, equation 2.6b in Yin, 1997
!section inside curly braces only - multiply result by 1/t1 to get mbar

implicit none

real(sp),               intent(in) :: t1
real(sp),               intent(in) :: a
real(sp),               intent(in) :: b
real(sp), dimension(:), intent(in) :: c

real(sp) :: wt1
real(sp) :: wpi

real(sp) :: e1
real(sp) :: e2

wpi  = 180. / (pi * w)
wt1  = rw * t1

if (a > b) then

  F = wpi * c(2) / sqrt(a**2 - b**2) * acos((b + a * cos(wt1)) / (a + b * cos(wt1))) + c(3) * t1

else if (a < b) then

  e1 = sqrt((b + a) * (1. + cos(wt1))) + sqrt((b - a) * (1. - cos(wt1)))
  e2 = sqrt((b + a) * (1. + cos(wt1))) - sqrt((b - a) * (1. - cos(wt1)))

  F = wpi * c(2) / sqrt(b**2 - a**2) * log(e1 / e2) + c(3) * t1

else

  F = wpi * c(2) / a * tan(wt1 / 2.) + c(3) * t1

end if

!write(stdout,*)'F ab ',a,b
!write(stdout,*)'F X  ',wpi * c(2) / sqrt(b**2 - a**2) * log(e1 / e2)
!write(stdout,*)'Fc3t1',c(3) * t1

end function F

!----------------------------------------------------------------------------------------------------------------

! subroutine radslope(grid,day,dayl,delta,lat,toa_sw,direct,diffuse,sloperad)
!
! ! Subroutine to calculate slope and aspect effect on direct beam shortwave radiation
! ! Equations taken from MT-CLIM model (mtclim43.c), coded for ALVAR by Leo O Lai (Aug, 2021)
!
! use parametersmod, only : d2r,r2d,pi
! use statevarsmod,    only : topovars
!
! integer(i4), intent(in)  :: grid          ! Grid number
! integer(i4), intent(in)  :: day           ! Julian day (1 to 366)
! real(sp),    intent(in)  :: dayl          ! Daylength (hours)
! real(sp),    intent(in)  :: delta         ! Solar declination angle (degrees)
! real(dp),    intent(in)  :: lat           ! Local latitude (degrees)
! real(sp),    intent(in)  :: toa_sw        ! Top of the atmosphere downwelling shortwave rad (kJ m-2 d-1)
! real(sp),    intent(in)  :: direct        ! Direct beam surface downwelling shortwave (kJ m-2 d-1)
! real(sp),    intent(in)  :: diffuse       ! Diffuse surface downwelling shortwave (kJ m-2 d-1)
! real(sp),    intent(out) :: sloperad      ! Ratio of shortwave radiation on slope to flat surface (sloperad/flatrad)
!
! real(sp), parameter :: secperrad = 13750.9871         ! seconds per raidna of hour angle
! real(sp), parameter :: mindecl   = -0.4092729         ! Minimum declination angle (radians)
! real(sp), parameter :: radperday = 0.017214           ! Radians per day
! real(sp), parameter :: radperdeg = 0.01745329         ! Radians per degree
! real(sp), parameter :: daysoff   = 11.25              ! Julian day offset of winter solstice
!
! real(sp), parameter, dimension(21) :: optam = [2.90,3.05,3.21,3.39,3.69,3.82,4.07,4.37,4.72,5.12,5.60, &
!                                               6.18,6.88,7.77,8.90,10.39,12.44,15.36,19.79,26.96,30.00]
!
! ! Local variables
! real(sp) :: dslope            ! Topographic slope (degrees)
! real(sp) :: daspect           ! Topographic aspect (degrees)
! real(sp) :: rlat              ! Local latitude (radians)
! real(sp) :: rdsol             ! Solar declination angle (radians)
! real(sp) :: rslope            ! Topographic slope (radians)
! real(sp) :: raspect           ! Topographic aspect (radians)
!
! real(sp) :: coslat            ! Cosine of latitude
! real(sp) :: sinlat            ! Sine of latitude
! real(sp) :: cosslp            ! Cosine of slope
! real(sp) :: sinslp            ! Sine of slope
! real(sp) :: cosasp            ! Cosine of aspect
! real(sp) :: sinasp            ! Sine of aspect
! real(sp) :: cosdecl           ! Cosine of solar declination angle
! real(sp) :: sindecl           ! Sine of solar declination angle
! real(sp) :: coszeh            ! Cosine of zenith angle for east horizon
! real(sp) :: coszwh            ! Cosine of zenith angle for west horizon
! real(sp) :: cosegeom
! real(sp) :: sinegeom
!
! real(sp) :: cosh              ! Cosine of hour angle
! real(sp) :: sinh              ! Sine of the hour angle
! real(sp) :: dt                ! Duration of hour angle time-step (seconds)
! real(sp) :: dh                ! Angle of hour angle time-step (radians)
! real(sp) :: h                 ! Hour angle (radians)
!
! real(sp) :: hss               ! Hour angle at sunset (radians) >> -hss is the hour angle at sunrise
! real(sp) :: bsg1              ! Beam-slope geometry component 1
! real(sp) :: bsg2              ! Beam-slope geometry component 2
! real(sp) :: bsg3              ! Beam-slope geometry component 3
!
! real(sp) :: cza               ! Cosine of solar zenith angle
! real(sp) :: cbsa              ! Cosine of beam-slope angle
!
! real(sp) :: dir_flat_topa     ! Direct top-of-atmosphere shortwave radiation on flat surface (kJ m-2 d-1)
! real(sp) :: dir_beam_topa     ! Direct top-of-atmosphere shortwave radiation on sloped surface (kJ m-2 d-1)
!
! real(sp) :: sum_flatrad       ! Sum of direct shortwave radiation on flat surface (kJ m-2 d-1)
! real(sp) :: sum_sloperad      ! Sum of direct shortwave radiation on sloped surface during hours with direct beam (kJ m-2 d-1)
!
! real(sp) :: am
! integer(i4) :: ami
!
! !------
! ! Point to variables
! dslope  = topovars(grid)%slope
! daspect = topovars(grid)%aspect
!
! ! Convert degrees to radians
! rlat    = lat * d2r
! rslope  = dslope * d2r
! raspect = daspect * d2r
!
! rdsol = mindecl * cos((day + daysoff) * radperday)
!
! !------
! ! Set up the trigonometric variables for lat, slope and aspect
! if (rlat > 1.5707)  rlat = 1.5707
! if (rlat < -1.5707) rlat = -1.5707
!
! coslat = cos(rlat)
! sinlat = sin(rlat)
! cosslp = cos(rslope)
! sinslp = sin(rslope)
! cosasp = cos(raspect)
! sinasp = sin(raspect)
!
! cosdecl = cos(rdsol)
! sindecl = sin(rdsol)
!
! cosegeom = coslat * cosdecl
! sinegeom = sinlat * sindecl
!
! !------
! ! Calculate hour angle of sunset
! hss = (dayl * 3600 / 2.0) / secperrad
!
! ! Calculate bean-slope geometry (bsg)
! bsg1 = -sinslp * sinasp * cosdecl;
! bsg2 = (-cosasp * sinslp * sinlat + cosslp * coslat) * cosdecl
! bsg3 = (cosasp * sinslp * coslat + cosslp * sinlat) * sindecl
!
! !------
! ! Set time step (3600 seconds = 1 hour)
! dt = 3600.
! dh = dt / secperrad
!
! h = -hss
!
! sum_sloperad = 0.0
! sum_flatrad = 0.0
!
! coszeh = cos(1.570796 - (rlat - rdsol))
! coszwh = cos(-1.570796 - (rlat - rdsol))
!
! ! coszeh = cos(0.0*d2r + rslope)
! ! coszwh = cos(180.0*d2r - rslope)
!
! do
!
!   cosh = cos(h)
!   sinh = sin(h)
!
!   cza = cosegeom * cosh + sinegeom
!
!   cbsa = sinh * bsg1 + cosh * bsg2 + bsg3
!
!   dir_flat_topa = 0.0
!
!   if (cza > 0.0) then
!
!     dir_flat_topa = toa_sw * cza
!
!     ! am = 1.0 / (cza + 0.0000001)
!     !
!     ! if (am > 2.9) then
!     !   ami = int(acos(cza) * d2r) - 69
!     !   if (ami < 0) ami = 0
!     !   if (ami > 20) ami = 20
!     !   am = optam(ami)
!     ! end if
!
!     sum_flatrad = sum_flatrad + dir_flat_topa
!
!   end if
!
!   !---
!
!   ! if (cbsa < 0) print *, h, acos(cza)*r2d,acos(coszeh)*r2d, acos(coszwh)*r2d,dslope,daspect
!
!   if ((h < 0.0 .and. acos(cza) > acos(coszeh) .and. cbsa > 0.0) .or. &
!       (h >= 0.0 .and. acos(cza) < acos(coszwh) .and. cbsa > 0.0)) then
!
!     sum_sloperad = sum_sloperad + toa_sw * cbsa
!
!   end if
!
!   ! print *, acos(cza)*r2d, acos(coszeh)*r2d, acos(coszwh)*r2d
!
!   !---
!
!   h = h + dh
!
!   if (h > hss) exit
!
! end do
!
! sloperad = sum_sloperad / sum_flatrad
!
! ! if (topovars(grid)%elev > 920.) print *, topovars(grid)%elev, dayl, sloperad, sum_flatrad, sloperad/sum_flatrad, dslope, daspect, &
! !             acos(coszeh)*r2d, acos(coszwh)*r2d
!
! ! print *, lat, dayl,delta,rdsol,dslope,rslope,daspect,raspect,coszeh,coszwh,hss
!
!
! end subroutine radslope

!----------------------------------------------------------------------------------------------------------------

subroutine calctdew(day,s,e)

use parametersmod, only : pi,d2r,r2d
use statevarsmod,  only : monvars,dayvars

implicit none

integer, intent(in) :: day
integer, intent(in) :: s
integer, intent(in) :: e

real(sp), parameter :: a = 1.26   ! Constant for Ep formula in Kimball et al. (1997)
real(sp), parameter :: c = 0.66   ! Constant for Ep formula in Kimball et al. (1997)

real(sp) :: tmin
real(sp) :: tmax
real(sp) :: tmean
real(sp) :: aprec    ! Annual rainfall in m

real(sp) :: tmin_K
real(sp) :: tmax_K
real(sp) :: tmean_K

real(sp) :: Lv    ! Latent heat of vaporization in J K-1
real(sp) :: Pw    ! Water density in kg m-3
real(sp) :: dSVP  ! Rate of change of saturation vapour pressure in Pa K-1
real(dp) :: Ep    ! Potential evapotranspiration in kg m-2 s-1
real(dp) :: EF    ! Ratio of Ep to annual precipitation (m)

! real(sp) :: ampl  ! Seasonal variation in daylength in hour
! real(sp) :: sunrise  ! Sunrise time
! real(sp) :: sunset  ! Sunset time
real(sp) :: dayl  ! Daylength in seconds
! real(sp) :: dsol  ! Solar declination angle in degrees
! real(sp) :: w  ! Elevation angle in degrees
! real(sp) :: Ma    ! Mean anomaly of orbit in rad
! real(sp) :: va    ! True anomaly of orbit in rad
! real(sp) :: Rd    ! Actual distance between sun and Earth at yearday in Gm

real(sp) :: Rn    ! Daily average insolation in W m-2
real(sp) :: Gn    ! Daily average surface conductive energy flux in W m-2

real(sp) :: tdew_K    ! Dew point temperature in degree Celcius
real(sp) :: tdew
! real(sp) :: RH
! real(sp) :: Tw

!---

tmin = dayvars(day)%tmin
tmax = dayvars(day)%tmax
tmean = dayvars(day)%tmean
aprec = sum(monvars%pre(s:e))           ! Annual rainfall in mm
dayl = dayvars(day)%dayl * 3600.    ! Daylength in seconds

!------

tmean_K = tmean + 273.15
tmin_K  = tmin + 273.15
tmax_K  = tmax + 273.15

!------

! Density of water Graf (2009)
Pw = 1000. * (1 - (((tmean - 3.9863) ** 2) / 508929.2) * ((tmean + 288.9414) / (tmean + 68.12963)))


! dSVP = 6.1078 * exp((17.269 * tavg) / (237.3 + tavg)) * ((237.3 + tavg) * 17.269 - 17.269 * tavg) / ((237.3 + tavg) ** 2)   ! Running & Coughlan (1988)

! dSVP = 100. * dSVP     ! Convert mbar to Pa


Rn    = dayvars(day)%srad * 1000. / (3600. * 24)
Gn    = 0.1 * Rn                                          ! Kimball et al. (1997) estimation of surface conductive energy flux
dSVP  = desdT(tmean_K)
gamma = 65.05 + tmean_K * 0.064                           ! Psychrometer constant (Pa K-1)
Lv    = 1.91846e6 * (tmean_K / (tmean_K - 33.91)) ** 2    ! Henderson-Sellers (1984) in Davis et al. (2017)

!------

Ep = a * (dSVP / (dSVP + gamma)) * (Rn - Gn) / Lv   ! Kimball et al. (1997)

Ep = dayvars(day)%dpet / 86400.      ! Convert dpet from mm d-1 to mm s-1

! if (dayvars(day)%dpet > 0.) print *, Ep * 86400, dayvars(day)%dpet

!------

aprec = max(aprec, 30.0)
aprec = aprec / 1000.     ! Convert mm to m of rainfall

EF = ((Ep / Pw) * dayl) / aprec     ! Kimball et al. (1997)

!------

tdew_K = tmin_K * (-0.127 + 1.121 * (1.003 - 1.444 * EF + 12.312 * (EF**2) &      ! Kimball et al. (1997)
         - 32.766 * (EF**3)) + 0.0006 * (tmax_K - tmin_K))

tdew = tdew_K - 273.15

! RH = 100 * (EXP((17.625*Td)/(243.04+Td))/EXP((17.625*tavg)/(243.04+tavg)))
!
! Tw = tavg * atan(0.151977 * sqrt(RH + 8.313659)) + atan(tavg + RH)            &
!      - atan(RH - 1.676331) + 0.00391838 * (RH ** 1.5) * atan(0.023101 * RH)   &
!      - 4.686035                 ! From Stull (2011) Wet-Bulb Temperature from Relative Humidity and Air Temperature


dayvars(day)%tdew = tdew

end subroutine calctdew


!----------------------------------------------------------------------------------------------------------------


end module radiationmod
