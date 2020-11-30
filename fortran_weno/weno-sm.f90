module weno_smoothness_functions

    implicit none

    ! module variable definitions

contains

!*****************************************************************************
!* A collection of smoothness functions,
!*
!*        written by Mohamad Ibrahim Cheikh 
!*
!*------------------------
!* List of subroutines:
!
!   smoothness_3
!
!*------------------------
!*
!*****************************************************************************

subroutine smoothness_3 (f, n, sigma)

  integer, intent(in) :: n! number of faces in the array
  real*8, intent(in) , dimension(n) :: f
  real*8, intent(out) , dimension(n,3) :: sigma

  integer :: i
  real*8  :: sigma0, sigma1, sigma2

  sigma = 0.0d0

  do i = 4, n -3

      sigma0 = (+3.33333333333333) * f(i) * f(i) + &
               (-10.3333333333333) * f(i) * f(i + 1) + &
               (+3.66666666666667) * f(i) * f(i + 2) + &
               (+8.33333333333333) * f(i + 1) * f(i + 1) + &
               (-6.33333333333333) * f(i + 1) * f(i + 2) + &
               (+1.33333333333333) * f(i + 2) * f(i + 2)

      sigma1 = (+1.33333333333333) *  f(i - 1) * f(i - 1) + &
               (-4.33333333333333) *  f(i - 1) * f(i) +     &
               (+1.66666666666667) *  f(i - 1) * f(i + 1) + &
               (+4.33333333333333) *  f(i) * f(i) +         &        
               (-4.33333333333333) *  f(i) * f(i + 1) +     &
               (+1.33333333333333) *  f(i + 1) * f(i + 1)

      sigma2 = (+1.33333333333333) *  f(i - 2) * f(i - 2) + &
               (-6.33333333333333) *  f(i - 2) * f(i - 1) + &
               (+3.66666666666667) *  f(i - 2) * f(i)      + &
               (+8.33333333333333) *  f(i - 1) * f(i - 1) + &
               (-10.3333333333333) *  f(i - 1) * f(i)  +     &
               (+3.33333333333333) *  f(i) * f(i) 

      sigma(i , 1)  = sigma0
      sigma(i , 2)  = sigma1
      sigma(i , 3)  = sigma2

  end do

end subroutine smoothness_3

end module weno_smoothness_functions