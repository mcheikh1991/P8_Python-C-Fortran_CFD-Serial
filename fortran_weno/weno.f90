module weno

    use weno_smoothness_functions
    use weno_weighting_functions
    use weno_reconstruction_functions
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


 subroutine weno_left_3(q,n,qr)

  integer, intent(in) :: n
  real*8, intent(in) , dimension(n) :: q
  real*8, intent(out), dimension(n) :: qr

  real*8, dimension(n,3) :: sigma , omega

  call smoothness_3 (q, n, sigma)
  call weights_left_3 (sigma, n, omega)
  call reconstruct_left_3(q, omega, n, qr)

  end subroutine weno_left_3


 subroutine weno_right_3(q,n,qr)

  integer, intent(in) :: n 
  real*8, intent(in) , dimension(n) :: q
  real*8, intent(out), dimension(n) :: qr

  real*8, dimension(n,3) :: sigma , omega

  call smoothness_3 (q, n, sigma)
  call weights_right_3 (sigma, n, omega)
  call reconstruct_right_3(q, omega, n, qr)

  end subroutine weno_right_3


end module weno