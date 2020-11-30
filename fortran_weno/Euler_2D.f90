module Euler

    use weno
    implicit none

contains

 subroutine Euler_2D(t,rho,rhoU,rhoV,rhoE,n,m)

  integer, intent(in) :: n, m, t
  real*8, intent(in) , dimension(n,m) :: rho,rhoU,rhoV,rhoE

  call smoothness_3 (q, n, sigma)
  call weights_left_3 (sigma, n, omega)
  call reconstruct_left_3(q, omega, n, qr)

 end subroutine Euler_2D


end module Euler