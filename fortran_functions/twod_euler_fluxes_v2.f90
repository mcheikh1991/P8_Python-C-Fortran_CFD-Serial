module Convective_Schemes_2D

    implicit none

    ! module variable definitions

contains

!*****************************************************************************
!* A collection of two-dimensional Euler numerical fluxes, Version 2 (2010),
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!*------------------------
!* List of Flux Functions:
!
!   Roe
!   Rotated-RHLL
!
!*------------------------
!*
!* These F90 routines were written and made available for download
!* for an educational purpose. Detailed descripstion of each numerical
!* flux can be found in the original paper, or in popular textbooks, or
!* in the (will-be-available) second volume of "I do like CFD". 
!*
!* Note that all routines are not efficietly implemented for clarity; you can 
!* improve the efficiecy and also you can covert it to double precision 
!* version if you wish.
!*
!* NOTES: There were bugs in the Rotated-RHLL solver in version 1.
!*        It has been fixed and tested for a shock-diffraction problem.
!*
!* This file may be updated in future. (to add more flux functions.)
!*****************************************************************************

subroutine Test(rho_rl, rhoU_rl, rhoV_rl, rhoE_rl, rho_ud, rhoU_ud, rhoV_ud, rhoE_ud, Fn_rl,Fn_ud, dx, dy, nx, ny,Answer)

  integer, intent(in) :: nx , ny! number of faces in the array
  real*8, intent(in) :: dx , dy
  real*8, intent(in) , dimension(2,nx,ny+1) :: rho_rl, rhoU_rl, rhoV_rl, rhoE_rl , Fn_rl
  real*8, intent(in) , dimension(2,nx+1,ny) :: rho_ud, rhoU_ud, rhoV_ud, rhoE_ud , Fn_ud
  real*8, intent(out), dimension(4,nx,ny+1) :: Answer

  Answer(1,:,:) = rho_rl(1,:,:)
  Answer(2,:,:) = rho_rl(2,:,:)
  Answer(3,:,:) = rho_rl(1,:,:)
  Answer(4,:,:) = rho_rl(2,:,:)

end subroutine Test

subroutine Roe_Matrix(rho_rl,rhoU_rl,rhoV_rl,rhoE_rl,rho_ud, rhoU_ud, rhoV_ud, rhoE_ud, Fn_rl,Fn_ud,dxy,nx,ny,Answer)

  integer, intent(in) :: nx , ny! number of faces in the array
  real*8, intent(in) , dimension(2) :: dxy
  real*8, intent(in) , dimension(2,nx,ny+1) :: rho_rl, rhoU_rl, rhoV_rl, rhoE_rl , Fn_rl
  real*8, intent(in) , dimension(2,nx+1,ny) :: rho_ud, rhoU_ud, rhoV_ud, rhoE_ud , Fn_ud
  real*8, intent(out), dimension(4,nx,ny):: Answer
  
  integer :: i , j
  real*8 :: dx, dy       
  real*8 :: uL(4), uR(4) !  conservative variables rho*[1, u, v, E]
  real*8 :: Flux(4)      !  flux
  real*8 :: fx, fy       !  face normal vector, 
  real*8, dimension(nx,ny+1) :: x_phi , x_phiUp , x_phiVp , x_phiEp ! x-flux
  real*8, dimension(nx+1,ny) :: y_phi , y_phiUp , y_phiVp , y_phiEp ! y-flux
  real*8, dimension(nx,ny):: L_rho, L_rho_u, L_rho_v, L_E
  real*8:: one

  one = 1.0
  dx = dxy(1)
  dy = dxy(2)

  ! Loop Over the faces in x-direction (y-fixed)

  do j = 1, nx
    do i   = 1, ny + 1

      uL(1) = rho_rl(1,j,i)
      uL(2) = rhoU_rl(1,j,i)
      uL(3) = rhoV_rl(1,j,i)
      uL(4) = rhoE_rl(1,j,i)

      uR(1) = rho_rl(2,j,i)
      uR(2) = rhoU_rl(2,j,i)
      uR(3) = rhoV_rl(2,j,i)
      uR(4) = rhoE_rl(2,j,i)

      fx    = Fn_rl(1,j,i)
      fy    = Fn_rl(2,j,i)

      Flux = Roe(uL, uR, fx, fy)

      x_phi(j,i)   = Flux(1)
      x_phiUp(j,i) = Flux(2)
      x_phiVp(j,i) = Flux(3)
      x_phiEp(j,i) = Flux(4)

   end do
  end do

  ! Loop Over the faces in y-direction (x-fixed)

  do j = 1, nx + 1 
    do i   = 1, ny 

      uL(1) = rho_ud(1,j,i)
      uL(2) = rhoU_ud(1,j,i)
      uL(3) = rhoV_ud(1,j,i)
      uL(4) = rhoE_ud(1,j,i)

      uR(1) = rho_ud(2,j,i)
      uR(2) = rhoU_ud(2,j,i)
      uR(3) = rhoV_ud(2,j,i)
      uR(4) = rhoE_ud(2,j,i)

      fx    = Fn_ud(1,j,i)
      fy    = Fn_ud(2,j,i)

      Flux = Roe(uL, uR, fx, fy)

      y_phi(j,i)   = Flux(1)
      y_phiUp(j,i) = Flux(2)
      y_phiVp(j,i) = Flux(3)
      y_phiEp(j,i) = Flux(4)

   end do
  end do

  ! FIX THESE
  L_rho   = (-one/(dx))*(x_phi(:,2:ny+1)   - x_phi(:,1:ny))   + (-one/(dy))*(y_phi(2:ny+1,:)   - y_phi(1:ny,:))
  L_rho_u = (-one/(dx))*(x_phiUp(:,2:ny+1) - x_phiUp(:,1:ny)) + (-one/(dy))*(y_phiUp(2:ny+1,:) - y_phiUp(1:ny,:))
  L_rho_v = (-one/(dx))*(x_phiVp(:,2:ny+1) - x_phiVp(:,1:ny)) + (-one/(dy))*(y_phiVp(2:ny+1,:) - y_phiVp(1:ny,:))
  L_E     = (-one/(dx))*(x_phiEp(:,2:ny+1) - x_phiEp(:,1:ny)) + (-one/(dy))*(y_phiEp(2:ny+1,:) - y_phiEp(1:ny,:))

  Answer(1,:,:) = L_rho
  Answer(2,:,:) = L_rho_u
  Answer(3,:,:) = L_rho_v
  Answer(4,:,:) = L_E

end subroutine Roe_Matrix







subroutine Roe_Array(rho,rhoU,rhoV,rhoE,Fn,n,Answer)

  integer, intent(in) :: n! number of faces in the array
  real*8, intent(in) , dimension(2,n) :: rho
  real*8, intent(in) , dimension(2,n) :: rhoU
  real*8, intent(in) , dimension(2,n) :: rhoV
  real*8, intent(in) , dimension(2,n) :: rhoE
  real*8, intent(in) , dimension(2,n) :: Fn
  real*8, intent(out), dimension(4,n) :: Answer

  integer :: i 
  real*8:: uL(4), uR(4) !  conservative variables rho*[1, u, v, E]
  real*8:: fx, fy       !  face normal vector, [nx, ny] (Left-to-Right)

   do i   = 1, n

    uL(1) = rho (1,i)
    uL(2) = rhoU(1,i)
    uL(3) = rhoV(1,i)
    uL(4) = rhoE(1,i)

    uR(1) = rho (2,i)
    uR(2) = rhoU(2,i)
    uR(3) = rhoV(2,i)
    uR(4) = rhoE(2,i)

    fx    = Fn(1,i)
    fy    = Fn(2,i)

    Answer(:,i) = Roe(uL, uR, fx, fy)
   end do

end subroutine Roe_Array

!*****************************************************************************
!* -- Roe's Flux Function ---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!* 
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
function Roe(uL, uR, nx, ny)
 real*8:: uL(4), uR(4) !  Input: conservative variables rho*[1, u, v, E]
 real*8:: nx, ny       !  Input: face normal vector, [nx, ny] (Left-to-Right)
 real*8:: Roe(4)       ! Output: Roe flux function (upwind)
!Local constants
 real*8:: gamma                          ! Ratio of specific heat.
 real*8:: zero, fifth, half, one, two    ! Numbers
!Local variables
 real*8:: tx, ty       ! Tangent vector (perpendicular to the face normal)
 real*8:: vxL, vxR, vyL, vyR             ! Velocity components.
 real*8:: rhoL, rhoR, pL, pR             ! Primitive variables.
 real*8:: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
 real*8:: aL, aR, HL, HR                 ! Speeds of sound.
 real*8:: RT,rho,vx,vy,H,a,vn, vt        ! Roe-averages
 real*8:: drho,dvx,dvy,dvn,dvt,dp,dV(4)  ! Wave strenghs
 real*8:: ws(4),dws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
 real*8:: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
 integer :: i, j

!Constants.
     gamma = 1.4
      zero = 0.0
     fifth = 0.2
      half = 0.5
       one = 1.0
       two = 2.0

!Tangent vector (Do you like it? Actually, Roe flux can be implemented 
! without any tangent vector. See "I do like CFD, VOL.1" for details.)
  tx = -ny
  ty = nx

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
     vxL = uL(2)/uL(1)
     vyL = uL(3)/uL(1)
     vnL = vxL*nx+vyL*ny
     vtL = vxL*tx+vyL*ty
      pL = (gamma-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(4) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
     vxR = uR(2)/uR(1)
     vyR = uR(3)/uR(1)
     vnR = vxR*nx+vyR*ny
     vtR = vxR*tx+vyR*ty
      pR = (gamma-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(4) + pR ) / rhoR

!First compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
    vx = (vxL+RT*vxR)/(one+RT)
    vy = (vyL+RT*vyR)/(one+RT)
     H = ( HL+RT* HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*(vx*vx+vy*vy)) )
    vn = vx*nx+vy*ny
    vt = vx*tx+vy*ty

!Wave Strengths
   drho = rhoR - rhoL 
     dp =   pR - pL
    dvn =  vnR - vnL
    dvt =  vtR - vtL

  dV(1) = (dp - rho*a*dvn )/(two*a*a)
  dV(2) = rho*dvt
  dV(3) =  drho - dp/(a*a)
  dV(4) = (dp + rho*a*dvn )/(two*a*a)

!Wave Speed
  ws(1) = abs(vn-a)
  ws(2) = abs(vn)
  ws(3) = abs(vn)
  ws(4) = abs(vn+a)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
! only for the nonlinear fields.
  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(4) = fifth
   if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

!Right Eigenvectors
  Rv(1,1) = one    
  Rv(2,1) = vx - a*nx
  Rv(3,1) = vy - a*ny
  Rv(4,1) =  H - vn*a

  Rv(1,2) = zero
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = vt*a

  Rv(1,3) = one
  Rv(2,3) = vx
  Rv(3,3) = vy 
  Rv(4,3) = half*(vx*vx+vy*vy)

  Rv(1,4) = one
  Rv(2,4) = vx + a*nx
  Rv(3,4) = vy + a*ny
  Rv(4,4) =  H + vn*a

!Dissipation Term
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
   end do
  end do

!Compute the flux.
  fL(1) = rhoL*vnL
  fL(2) = rhoL*vnL * vxL + pL*nx
  fL(3) = rhoL*vnL * vyL + pL*ny
  fL(4) = rhoL*vnL *  HL

  fR(1) = rhoR*vnR
  fR(2) = rhoR*vnR * vxR + pR*nx
  fR(3) = rhoR*vnR * vyR + pR*ny
  fR(4) = rhoR*vnR *  HR

  Roe = half * (fL + fR - diss)

end function Roe

!*****************************************************************************
!* -- Rotated-Roe-HLL Flux Function ---
!*
!* H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free, Boundary-Layer
!* Resolving, Rotated-Hybrid Riemann Solvers,
!* Journal of Computational Physics, 227, pp. 2560-2581, 2008.
!*
!* The most robust Riemann solver known to the author (in terms of nonlinear
!* instability such as carbuncle).
!*
!* NB: At a boundary face, need to switch to a geometric normal vector:
!*               (nx2,ny2)=(nx, ny), (nx1,ny1)=(-ny,nx).
!*     This is not implemented here. It requires information on whether
!*     the geometric normal, (nx,ny), is on a boundary face or not.
!*     It shouldn't be difficult for you to implement it.
!*
!* Katate Masatsuka, February 2010. http://www.cfdbooks.com
!*****************************************************************************
function Rotated_RHLL(uL, uR, nx, ny)
 real*8:: uL(4), uR(4)    !  Input: conservative variables rho*[1, u, v, E]
 real*8:: nx, ny          !  Input: face normal vector, [nx, ny] (Left-to-Right)
 real*8:: Rotated_RHLL(4) ! Output: Rotated_RHLL flux function.
!Local constants
 real*8:: gamma                          ! Ratio of specific heat.
 real*8:: zero, fifth, half, one, two    ! Numbers
 real*8:: eps                            ! 
!Local variables
 real*8:: nx1, ny1, nx2, ny2             ! Rotated normals, n1 and n2
 real*8:: tx, ty                         ! Tangent vector (taken as n1)
 real*8:: alpha1, alpha2                 ! Projections of the new normals
 real*8:: vxL, vxR, vyL, vyR             ! Velocity components.
 real*8:: rhoL, rhoR, pL, pR             ! Primitive variables.
 real*8:: vnL, vnR, vtL, vtR             ! Normal and tagent velocities
 real*8:: aL, aR, HL, HR                 ! Speeds of sound and total enthalpy
 real*8:: RT,rho,vx,vy,H,a               ! Roe-averages
 real*8:: vn, vt                         ! Normal and tagent velocities(Roe-average)
 real*8:: drho,dvx,dvy,dvn,dvt,dp,dV(4)  ! Wave strenghs
 real*8:: abs_dq                         ! Magnitude of the velocity difference
 real*8:: abs_ws(4),ws(4),dws(4), Rv(4,4)! Wave speeds and right-eigevectors
 real*8:: SRp,SLm                        ! Wave speeds for the HLL part
 real*8:: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
 real*8:: temp
 integer :: i, j

!Constants.
     gamma = 1.4
      zero = 0.0
     fifth = 0.2
      half = 0.5
       one = 1.0
       two = 2.0
       eps = 1.0e-5 ! 1.0e-12 in the original paper (double precision)

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
     vxL = uL(2)/uL(1)
     vyL = uL(3)/uL(1)
      pL = (gamma-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(4) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
     vxR = uR(2)/uR(1)
     vyR = uR(3)/uR(1)
      pR = (gamma-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(4) + pR ) / rhoR

     vnL = vxL*nx + vyL*ny
     vnR = vxR*nx + vyR*ny

!Compute the flux.
   fL(1) = rhoL*vnL
   fL(2) = rhoL*vnL * vxL + pL*nx
   fL(3) = rhoL*vnL * vyL + pL*ny
   fL(4) = rhoL*vnL *  HL

   fR(1) = rhoR*vnR
   fR(2) = rhoR*vnR * vxR + pR*nx
   fR(3) = rhoR*vnR * vyR + pR*ny
   fR(4) = rhoR*vnR *  HR

!Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
!(NB: n1 and n2 may need to be frozen at some point during 
!     a steady calculation to fully make it converge. For time-accurate 
!     calculation, this is fine.)
! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    abs_dq = sqrt( (vxR-vxL)**2+(vyR-vyL)**2 )
  if ( abs_dq > eps) then
       nx1 = (vxR-vxL)/abs_dq
       ny1 = (vyR-vyL)/abs_dq
  else
    nx1 = -ny 
    ny1 =  nx
  endif
    alpha1 = nx * nx1 + ny * ny1 
!   To make alpha1 always positive.
      temp = sign(one,alpha1)
       nx1 = temp * nx1
       ny1 = temp * ny1
    alpha1 = temp * alpha1

! Take n2 as perpendicular to n1.
       nx2 = -ny1
       ny2 =  nx1
    alpha2 = nx * nx2 + ny * ny2
!   To make alpha2 always positive.
      temp = sign(one,alpha2)
       nx2 = temp * nx2
       ny2 = temp * ny2
    alpha2 = temp * alpha2

!Now we are going to compute the Roe flux with n2 as the normal
!and n1 as the tagent vector, with modified wave speeds (5.12)

!Compute the Roe Averages
     RT = sqrt(rhoR/rhoL)
    rho = RT*rhoL
     vx = (vxL+RT*vxR)/(one+RT)
     vy = (vyL+RT*vyR)/(one+RT)
      H = ( HL+RT* HR)/(one+RT)
      a = sqrt( (gamma-one)*(H-half*(vx*vx+vy*vy)) )
     vn = vx*nx2+vy*ny2
     vt = vx*nx1+vy*ny1

!Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    vnL = vxL*nx2 + vyL*ny2
    vnR = vxR*nx2 + vyR*ny2
    vtL = vxL*nx1 + vyL*ny1
    vtR = vxR*nx1 + vyR*ny1

   drho = rhoR - rhoL 
     dp =   pR - pL
    dvn =  vnR - vnL
    dvt =  vtR - vtL

  dV(1) = (dp - rho*a*dvn )/(two*a*a)
  dV(2) =  rho*dvt/a
  dV(3) =  drho - dp/(a*a)
  dV(4) = (dp + rho*a*dvn )/(two*a*a)

!Wave Speeds for Roe flux part.
    ws(1) = vn-a
    ws(2) = vn
    ws(3) = vn
    ws(4) = vn+a
  abs_ws  = abs(ws)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
!only for the nonlinear fields.
  dws(1) = fifth
   if (abs_ws(1)<dws(1)) abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1))
  dws(4) = fifth
   if (abs_ws(4)<dws(4)) abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4)+dws(4))

!HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
   SRp = max( zero, vtR + aR, vt + a)
   SLm = min( zero, vtL - aL, vt - a)

!Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
   ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + two*alpha1*SRp*SLm )/ (SRp-SLm)

!Right Eigenvectors: with n2 as normal and n1 as tangent.
  tx = nx1
  ty = ny1

  Rv(1,1) = one    
  Rv(2,1) = vx - a*nx2
  Rv(3,1) = vy - a*ny2
  Rv(4,1) =  H - vn*a

  Rv(1,2) = zero
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = a*vt

  Rv(1,3) = one
  Rv(2,3) = vx
  Rv(3,3) = vy 
  Rv(4,3) = half*(vx*vx+vy*vy)

  Rv(1,4) = one
  Rv(2,4) = vx + a*nx2
  Rv(3,4) = vy + a*ny2
  Rv(4,4) =  H + vn*a

!Dissipation Term: Roe dissipation with the modified wave speeds.
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
   end do
  end do

!Compute the Rotated-RHLL flux.
  Rotated_RHLL = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss

end function Rotated_RHLL


end module Convective_Schemes_2D