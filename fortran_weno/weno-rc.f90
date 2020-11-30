module weno_reconstruction_functions

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
!   smoothness003
!
!*------------------------
!*
!*****************************************************************************

subroutine reconstruct_left_3 (f, omega, n, fr)

  integer, intent(in) :: n! number of faces in the array
  real*8, intent(in) , dimension(n) :: f
  real*8, intent(in) , dimension(n,3) :: omega
  real*8, intent(out), dimension(n) :: fr

  integer :: i
  real*8  :: acc, omega0 , omega1, omega2, fr0, fr1, fr2

  fr(1) = 0.0d0
  fr(2) = 0.0d0
  fr(3) = 0.0d0
  fr(n-2) = 0.0d0
  fr(n-1) = 0.0d0
  fr(n) = 0.0d0
  
  do i = 4, n -3

      omega0 = omega(i , 1)
      omega1 = omega(i , 2)
      omega2 = omega(i , 3)
   
      fr0 = (+1.83333333333333)  * f(i + 0)  + (-1.16666666666667) * f(i + 1)  + &
        (+0.333333333333333) * f(i + 2) 
      fr1 = (+0.333333333333333) * f(i - 1)  + (+0.833333333333333) * f(i + 0)  + &
        (-0.166666666666667) * f(i + 1)  
      fr2 = (-0.166666666666667) * f(i - 2)  + (+0.833333333333333) * f(i - 1)  + &
        (+0.333333333333333) * f(i + 0) 

      fr(i) = (omega0 * fr0) + (omega1 * fr1) + (omega2 * fr2)

  end do
end subroutine reconstruct_left_3

subroutine reconstruct_right_3 (f, omega, n, fr)
  
  integer, intent(in) :: n! number of faces in the array
  real*8, intent(in) , dimension(n) :: f
  real*8, intent(in) , dimension(n,3) :: omega
  real*8, intent(out), dimension(n) :: fr

  integer :: i
  real*8  :: acc, omega0 , omega1, omega2, fr0, fr1, fr2

  fr(1) = 0.0d0
  fr(2) = 0.0d0
  fr(3) = 0.0d0
  fr(n-2) = 0.0d0
  fr(n-1) = 0.0d0
  fr(n) = 0.0d0
  
  do i = 4, n -3

      omega0 = omega(i , 1)
      omega1 = omega(i , 2)
      omega2 = omega(i , 3)
   
      fr0 = (+0.333333333333333) * f(i + 0) + (+0.833333333333333)*f(i + 1)   + (-0.166666666666667) * f(i + 2)
      fr1 = (-0.166666666666667) * f(i - 1) + (+0.833333333333333) * f(i + 0) + (+0.333333333333333) * f(i + 1)
      fr2 = (+0.333333333333333) * f(i - 2) + (-1.16666666666667) * f(i - 1)  + (+1.83333333333333) * f(i + 0)

      fr(i) = (omega0 * fr0) + (omega1 * fr1) + (omega2 * fr2)

  end do
end subroutine reconstruct_right_3

end module weno_reconstruction_functions