module compute_faces

    use weno
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


 subroutine compute_lr_faces(phi,phi_b,boundary,nx,ny,phi_l,phi_r)

  integer, intent(in) :: nx, ny
  real*8,  intent(in) , dimension(nx,ny)   :: phi
  integer, intent(in) , dimension(2,ny)    :: phi_b_type ! Type of boundary
  real*8,  intent(in) , dimension(2,ny)    :: phi_b ! 1 : left boundary value, 2 : right boundary value
  real*8, intent(out) , dimension(nx+1,ny) :: phi_l , phi_r

  integer :: j, i
  real*8, dimension(nx+6) :: phi_1D, phi_1D_l_test , phi_1D_r_test
  real*8, dimension(nx+1) :: phi_1D_l , phi_1D_r

  ! Loop over the x direction (y is fixed):
  do j = 1,ny

    phi_1D(4:nx+3) =  phi(:,j)   ! Loop over the x direction (y is fixed):

    ! Apply Boundary Conditions Left Boundary:
    if (phi_b_type(1,j) .eq. 1) then ! Periodic
      phi_1D(3) = phi_1D(nx+3)
      phi_1D(2) = phi_1D(nx+2)
      phi_1D(1) = phi_1D(nx+1)
    else if (phi_b_type(1,j) .eq. 2) then ! Dirichlet
      phi_1D(3) = phi_b(1,j)
      phi_1D(2) = phi_b(1,j)
      phi_1D(1) = phi_b(1,j)
    else if (phi_b_type(1,j) .eq. 3) then ! ZeroGradient
      phi_1D(3) = phi_1D(4)
      phi_1D(2) = phi_1D(4)
      phi_1D(1) = phi_1D(4)
    else
      print *, 'Wrong boundary type: ', phi_b_type(1,j),' it should be either (1=Periodic, 2=Dirichlet, 3=ZeroGradient)'
    end if

    if (phi_b_type(1,j) .eq. 1) then ! Periodic
      phi_1D(nx+4) = phi_1D(4)
      phi_1D(nx+5) = phi_1D(5)
      phi_1D(nx+6) = phi_1D(6)
    else if (phi_b_type(1,j) .eq. 2) then ! Dirichlet
      phi_1D(nx+4) = phi_b(2,j)
      phi_1D(nx+5) = phi_b(2,j)
      phi_1D(nx+6) = phi_b(2,j)
    else if (phi_b_type(1,j) .eq. 3) then ! ZeroGradient
      phi_1D(nx+4) = phi_1D(nx+3)
      phi_1D(nx+5) = phi_1D(nx+3)
      phi_1D(nx+6) = phi_1D(nx+3)
    else
      print *, 'Wrong boundary type: ', phi_b_type(2,j),' it should be either (1=Periodic, 2=Dirichlet, 3=ZeroGradient)'
    end if

    ! Finding the face value
    call weno_left_3(phi_1D,nx,phi_1D_l_test)
    call weno_right_3(phi_1D,nx,phi_1D_r_test)

    ! Fixing the face value at the boundaries

    phi_1D_l    = phi_1D_l_test(4:nx+4)
    phi_1D_r    = phi_1D_r_test(3:nx+3)

    if (phi_b_type(1,j) .eq. 1) then ! Periodic
      phi_1D_r(1) = phi_1D_r(nx+3)
    else if (phi_b_type(1,j) .eq. 2) then ! Dirichlet
      phi_1D_r(1) = phi_b(1,j) ! Boundary at left
    else if (phi_b_type(1,j) .eq. 3) then ! ZeroGradient
      phi_1D_r(1) = phi_1D(4) ! Boundary at left
      phi_1D_l(1) = phi_1D(4)
    end if

    if (phi_b_type(2,j) .eq. 1) then ! Periodic
      phi_1D_l(nx+1) = phi_1D_l(1)
    else if (phi_b_type(2,j) .eq. 2) then ! Dirichlet
      phi_1D_l(nx+1) = phi_b(2,j) ! Boundary at right
    else if (phi_b_type(2,j) .eq. 3) then ! ZeroGradient
      phi_1D_l(nx+1) = phi_1D(nx+3) ! Boundary at right
      phi_1D_r(nx+1) = phi_1D(nx+3)
    end if

    ! Saving the Face Value
      phi_l(:,j) = phi_1D_l
      phi_r(:,j) = phi_1D_r
  
  end do
end subroutine compute_lr_faces

end module compute_faces