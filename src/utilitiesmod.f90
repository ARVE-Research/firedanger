module utilitiesmod

use parametersmod, only : i4,sp,dp

implicit none

public :: matsol
public :: findloc

contains

! ------------------------------------------------------------------------------------------------------------------

subroutine matsol(mat,sol)

! Provides matrix solution to X matrix in a A * X = B system using LU decomposition
! Code adapted from Press et al. (1996) Numerical Recipes in Fortran 90 : The Art of Parallel Scientific Computing
! 2nd Edition (P.1016-1017)

implicit none

real(sp), dimension(:,:), intent(inout) :: mat
real(sp), dimension(:)  , intent(inout) :: sol

integer(i4), dimension(size(sol)) :: indx
real(sp)   , dimension(size(sol)) :: mv
real(sp)   , dimension(:,:), allocatable   :: prod
integer(i4), dimension(1)   :: maxl

real(sp)   , parameter      :: tiny_sp = 1.0e-38_sp

integer(i4)                 :: i, n, k, ll
integer(i4)                 :: max
real(sp)                    :: summ

! ----------------------------------

n = size(sol)

mv = 1. / maxval(abs(mat), dim=2)

!---

do i = 1, n

  maxl = maxloc(mv(i:n) * abs(mat(i:n,i)))

  max = (i - 1) + maxl(1)

  indx(i) = max

  !---

  if (mat(i,i) == 0.) mat(i,i) = tiny_sp

  mat(i+1:n, i) = mat(i+1:n, i) / mat(i,i)

  !---

  allocate(prod(i+1:n, i+1:n))

  prod = spread(mat(i+1:n, i), dim=2, ncopies=size(mat(i, i+1:n)))

  prod = prod * spread(mat(i, i+1:n), dim=1, ncopies=size(mat(i+1:n, i)))

  !---

  mat(i+1:n, i+1:n) = mat(i+1:n, i+1:n) - prod

  deallocate(prod)

end do

!---

k = 0

do i = 1, n

  ll = indx(i)
  summ = sol(ll)
  sol(ll) = sol(i)

  if (k /= 0) then

    summ = summ - dot_product(mat(i, k:i-1), sol(k:i-1))

  else if (summ /= 0.) then

    k = i

  end if

  sol(i) = summ

end do

!---

do i = n, 1, -1

  sol(i) = (sol(i) - dot_product(mat(i, i+1:n), sol(i+1:n))) / mat(i,i)

end do

end subroutine matsol

! ------------------------------------------------------------------------------------------------------------------

subroutine findloc(mat,x,loc)

! returns the index (loc) of a value (x) in a 1D array (mat)
! adapted from Press et al. (1996) Numerical Recipes in Fortran 90 : The Art of Parallel Scientific Computing
! 2nd Edition

implicit none

real(sp), dimension(:), intent(in)  :: mat
real(sp),               intent(in)  :: x
integer(i4),            intent(out) :: loc

real(sp), dimension(:), allocatable :: diff
integer(i4) :: len

!----

len = size(mat)

allocate(diff(len))

diff = abs(mat - x)

loc = minloc(diff,dim=1)

end subroutine findloc

! ------------------------------------------------------------------------------------------------------------------

subroutine tridiag(a,b,c,r,u)

! Subroutine to solve triadiagonal system of equations

! Solves for a vector u of size N the tridiagonal linear set using given by equation (2.4.1) using a
! a serial algorithm. Input vectors b (diagonal elements) and r (right-hand side) have size N,
! while a and c (off-diagonal elements) are not defined in the first and last elements, respectively.
! Based on Numerical Recipes in F77/F90

! Copied from ARVE-DGVM (Leo Lai, Jul 2021)

implicit none


real(dp), dimension(:), intent(in) :: a
real(dp), dimension(:), intent(in) :: b
real(dp), dimension(:), intent(in) :: c
real(dp), dimension(:), intent(in) :: r
real(dp), dimension(:), intent(out) :: u

integer :: n
integer :: k
real(dp) :: bet
real(dp), dimension(size(b)) :: gam

!----

n = size(b)

bet = b(1)

u(1) = r(1) / bet

!decomposition and forward substitution

do k = 2,n

  gam(k) = c(k-1) / bet

  bet = b(k) - a(k) * gam(k)

  u(k) = (r(k) - a(k) * u(k-1)) / bet

end do

!backsubstitution

do k = n-1,1,-1
  u(k) = u(k) - gam(k+1) * u(k+1)
end do

end subroutine tridiag

! ------------------------------------------------------------------------------------------------------------------

subroutine getmonth(day,ndyear,month,startday,endday)

! Get the month (1-12) from input Julian day and get the start and end Julian
! date of the month

implicit none

integer(i4), intent(in)  :: day               ! Julian day of the year (1-366)
integer(i4), intent(in)  :: ndyear            ! Mumber of days in the year (365 or 366)
integer(i4), intent(out) :: month             ! Month of the input day
integer(i4), intent(out) :: startday          ! Startday of the month
integer(i4), intent(out) :: endday            ! End day of the month

logical :: leap

!------

! Check if current year is leap year or not
if (ndyear == 365) leap = .FALSE.
if (ndyear == 366) leap = .TRUE.

! Assign month based on input Julian date
if (.not.leap) then

  !------

  if (day <= 31) then

    month    = 1
    startday = 1
    endday   = 31

  else if (day >= 32 .AND. day <= 59) then

    month    = 2
    startday = 32
    endday   = 59

  else if (day >= 60 .AND. day <= 90) then

    month    = 3
    startday = 60
    endday   = 90

  else if (day >= 91 .AND. day <= 120) then

    month    = 4
    startday = 91
    endday   = 120

  else if (day >= 121 .AND. day <= 151) then

    month    = 5
    startday = 121
    endday   = 151

  else if (day >= 152 .AND. day <= 181) then

    month    = 6
    startday = 152
    endday   = 181

  else if (day >= 182 .AND. day <= 212) then

    month    = 7
    startday = 182
    endday   = 212

  else if (day >= 213 .AND. day <= 243) then

    month    = 8
    startday = 213
    endday   = 243

  else if (day >= 244 .AND. day <= 273) then

    month    = 9
    startday = 244
    endday   = 273

  else if (day >= 274 .AND. day <= 304) then

    month    = 10
    startday = 274
    endday   = 304

  else if (day >= 305 .AND. day <= 334) then

    month    = 11
    startday = 305
    endday   = 334

  else if (day >= 335) then

    month    = 12
    startday = 335
    endday   = 365

  end if

  !------

else if (leap) then

  !------

  if (day <= 31) then

    month    = 1
    startday = 1
    endday   = 31

  else if (day >= 32 .AND. day <= 60) then

    month    = 2
    startday = 32
    endday   = 60

  else if (day >= 61 .AND. day <= 91) then

    month    = 3
    startday = 61
    endday   = 91

  else if (day >= 92 .AND. day <= 121) then

    month    = 4
    startday = 92
    endday   = 121

  else if (day >= 122 .AND. day <= 152) then

    month    = 5
    startday = 122
    endday   = 152

  else if (day >= 153 .AND. day <= 182) then

    month    = 6
    startday = 153
    endday   = 182

  else if (day >= 183 .AND. day <= 213) then

    month    = 7
    startday = 183
    endday   = 213

  else if (day >= 214 .AND. day <= 244) then

    month    = 8
    startday = 214
    endday   = 244

  else if (day >= 245 .AND. day <= 274) then

    month    = 9
    startday = 245
    endday   = 274

  else if (day >= 275 .AND. day <= 305) then

    month    = 10
    startday = 275
    endday   = 305

  else if (day >= 306 .AND. day <= 335) then

    month    = 11
    startday = 306
    endday   = 335

  else if (day >= 336) then

    month    = 12
    startday = 336
    endday   = 366

  end if

  !------

end if

end subroutine getmonth

! ------------------------------------------------------------------------------------------------------------------

end module utilitiesmod
