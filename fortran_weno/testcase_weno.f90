program Validation_Weno 

  use compute_faces
  use weno
  use weno_smoothness_functions
  use weno_weighting_functions
  implicit none

  real, parameter        :: Pi = 3.1415927
  integer, parameter     :: n  = 21

  real*8, dimension(n)   :: x
  real*8, dimension(n-1) :: phi, phi_l, phi_r
  real*8, dimension(n-1,3) :: sigma, omega_l , omega_r

  integer :: i
  real :: dx

  dx = (2*Pi)/(n-1)

  do i=1,n
    x(i) = 0 + (i-1)*dx
  end do

  phi = (sin(x(2:n))- sin(x(1:n-1))) / (x(2:n) - x(1:n-1))

  ! Validation of Smoothness Function
  call smoothness_3 (phi, n-1, sigma)
  open (unit = 1, file = "output/fortran-smooth-smoothness.txt")
 
  do i=1,n-1
    write(1,*) sigma(i,1) , sigma(i,2) ,  sigma(i,3) 
  end do

  ! Validation of Weighting Function
  call weights_left_3 (sigma, n-1, omega_l)
  call weights_right_3 (sigma, n-1, omega_r)
  open (unit = 2, file = "output/fortran-smooth-weights-left.txt")
  open (unit = 3, file = "output/fortran-smooth-weights-right.txt")
 
  do i=1,n-1
    write(2,*) omega_l(i,1) , omega_l(i,2) ,  omega_l(i,3) 
    write(3,*) omega_r(i,1) , omega_r(i,2) ,  omega_r(i,3) 
  end do

  ! Validation of Weighting Function
  call weno_left_3 (phi,n-1,phi_l)
  call weno_right_3(phi,n-1,phi_r)

  open (unit = 7, file = "output/fortran-smooth-left.txt")
  open (unit = 8, file = "output/fortran-smooth-right.txt")


  do i=1,n-1
    write(7,*) phi_l(i)
  end do

  do i=1,n-1
    write(8,*) phi_r(i)
  end do




end program Validation_Weno