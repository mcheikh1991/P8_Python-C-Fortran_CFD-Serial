module weno_weighting_functions

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
!   weights_left_3
!   weights_right_3
!
!*------------------------
!*
!*****************************************************************************

 subroutine weights_left_3 (sigma, n, omega)

  integer, intent(in) :: n 
  real*8, intent(in) , dimension(n,3) :: sigma
  real*8, intent(out), dimension(n,3) :: omega

  integer :: i
  real*8  :: sigma0, sigma1, sigma2, acc, omega0 , omega1, omega2

  omega = 0.0d0

  do i = 4, n -3

      sigma0 = sigma(i , 1)
      sigma1 = sigma(i , 2)
      sigma2 = sigma(i , 3)
   
      acc = 0.0
      omega0 = (+0.1) / ((sigma0 + 1.0e-6) * (sigma0 + 1.0e-6))
      acc = acc + omega0
      omega1 = (+0.6) / ((sigma1 + 1.0e-6) * (sigma1 + 1.0e-6))
      acc = acc + omega1
      omega2 = (+0.3) / ((sigma2 + 1.0e-6) * (sigma2 + 1.0e-6))
      acc = acc + omega2

      omega0 = omega0 / acc
      omega1 = omega1 / acc
      omega2 = omega2 / acc

      omega(i , 1)  = omega0
      omega(i , 2)  = omega1
      omega(i , 3)  = omega2

  end do
 end subroutine weights_left_3

subroutine weights_right_3 (sigma, n, omega)

  integer, intent(in) :: n
  real*8, intent(in) , dimension(n,3) :: sigma
  real*8, intent(out), dimension(n,3) :: omega

  integer :: i
  real*8  :: sigma0, sigma1, sigma2, acc, omega0 , omega1, omega2

  omega = 0.0d0

  do i = 4, n -3

      sigma0 = sigma(i , 1)
      sigma1 = sigma(i , 2)
      sigma2 = sigma(i , 3)
   
      acc = 0.0
      omega0 = (+0.3) / ((sigma0 + 1.0e-6) * (sigma0 + 1.0e-6))
      acc = acc + omega0
      omega1 = (+0.6) / ((sigma1 + 1.0e-6) * (sigma1 + 1.0e-6))
      acc = acc + omega1
      omega2 = (+0.1) / ((sigma2 + 1.0e-6) * (sigma2 + 1.0e-6))
      acc = acc + omega2

      omega0 = omega0 / acc
      omega1 = omega1 / acc
      omega2 = omega2 / acc
      
      omega(i , 1)  = omega0
      omega(i , 2)  = omega1
      omega(i , 3)  = omega2

  end do
end subroutine weights_right_3

end module weno_weighting_functions