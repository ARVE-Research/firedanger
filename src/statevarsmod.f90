module statevarsmod

use parametersmod, only : i2,i4,sp,dp

! This module contains the data structures that are read, calculated and stored by other subroutines

!---------------------------------------------------------------------

type mon_metvars
  ! Derived datatype for monthly met variables input
  ! Dimension allocated from number of months

  real(sp), allocatable, dimension(:) :: tmp        ! mean monthly temperature (degC)
  real(sp), allocatable, dimension(:) :: dtr        ! mean monthly diurnal temperature range (degC)
  real(sp), allocatable, dimension(:) :: pre        ! total monthly precipitation (mm)
  real(sp), allocatable, dimension(:) :: wet        ! number of days in the month with precipitation > 0.1 mm (days)
  real(sp), allocatable, dimension(:) :: cld        ! mean monthly cloud cover (fraction)
  real(sp), allocatable, dimension(:) :: wnd        ! mean monthly 10m windspeed (m s-1)
  integer(i4), allocatable, dimension(:) :: nd

end type mon_metvars

type(mon_metvars), target :: monvars     ! Dimension (inside type) allocate to number of months needed for calculation

!---------------------------------------------------------------------

type day_metvars
  ! Derived datatype for the daily met variables output
  real(sp) :: Ratm        ! Relative atmospheric pressure to sea-level (fraction)
  real(sp) :: Ratm30      ! Relative atmospheric pressure to 30m above sea-level (fraction)
  real(sp) :: Patm        ! Atmospheric pressure to sea-level (Pa)
  real(sp) :: Patm30      ! Atmospheric pressure to 30m above sea-level (Pa)

  real(sp) :: prec        ! 24 hour total precipitation (mm)
  real(sp) :: tmin        ! 24 hour mean minimum temperature (degC)
  real(sp) :: tmax        ! 24 hour mean maximum temperature (degC)
  real(sp) :: cldf        ! 24 hour mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
  real(sp) :: wind        ! wind speed (m s-1)

  ! New added variables (Leo O Lai 16 Apr 2021)
  real(sp) :: dayl           ! daylength (h)
  real(sp) :: dayl_n         ! daylength of the next day (h)
  integer(i4) :: sunrise     ! Current day sunrise hour (index of 24h hourly array)
  integer(i4) :: sunset      ! Current day sunset hour (index of 24h hourly array)
  integer(i4) :: sunrise_n   ! next day sunrise hour (index of 24h hourly array)
  integer(i4) :: dayhour     ! day time hours (h) --> from current day sunrise to sunset
  integer(i4) :: nighthour   ! night time hours (h) --> from current day sunset to next day sunrise

  real(sp) :: tmean       ! 24 hour mean temperature (degC)
  real(sp) :: tday        ! mean daytime temperature (degC)
  real(sp) :: tnight      ! mean nighttime temperature (degC)
  real(sp) :: Pjj         ! precipitation equitability index for calculating PET
  real(sp) :: tdew        ! dew point temperature (degC)
  real(sp) :: rhum        ! relative humidity (%)
  real(sp) :: dsol        ! solar declination angle (degree)
  real(sp) :: srad        ! downwelling surface shortwave radiation (kJ m-2 d-1)
  real(sp) :: srad_dir    ! direct beam downwelling shortwave raditaion (kJ m-2 d-1)
  real(sp) :: srad_dif    ! diffuse downwelling shortwave raditaion (kJ m-2 d-1)
  real(sp) :: lrad        ! upswelling surface longwave radiation (kJ m-2 d-1)
  real(sp) :: dpet        ! total potential evapotranspiration (mm)
  real(sp) :: daet        ! total actual evapotranspiration (mm)
  real(sp) :: alpha       ! ratio of daily AET/PET (fraction)
  real(sp) :: vpd         ! average daytime saturation vapor pressure deficit (Pa)

  ! Fire index variables
  real(sp) :: DF          ! Drought factor
  real(sp) :: KBDI        ! Keetch-Byram Drought Index (mm equivalent)
  real(sp) :: FFDI        ! Forest fire danger index Mark 5 meter

  real(sp) :: prec_accum  ! Accumulation of precipitation since 1 Mar (mm)
  real(sp) :: curing      ! Degree of grass curing (%)
  real(sp) :: fuelmc      ! Fuel moisture content
  real(sp) :: GFDI        ! Grassland fire danger index

  integer(i4) :: FFDI_cat ! Forest fire danger index category
  integer(i4) :: GFDI_cat ! Grassland fire danger index category

  real(dp), dimension(24) :: hprec    ! hourly precipitation (mm)

end type day_metvars

type(day_metvars), target, allocatable, dimension(:) :: dayvars ! Dimension allocate from number of days in year (365 or 366)

!---------------------------------------------------------------------



end module statevarsmod
