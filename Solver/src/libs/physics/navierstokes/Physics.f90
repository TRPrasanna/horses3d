!      Physics.f90
!      Created: 2011-07-20 09:17:26 -0400 
!      By: David Kopriva
!      From DSEM Code
!
!!     The variable mappings for the Navier-Stokes Equations are
!!
!!              Q(1) = rho
!!              Q(2) = rhou
!!              Q(3) = rhov
!!              Q(4) = rhow
!!              Q(5) = rhoe
!!     Whereas the gradients are:
!!              grad(1) = grad(u)
!!              grad(2) = grad(v)
!!              grad(3) = grad(w)
!!              grad(4) = grad(T)
!
!////////////////////////////////////////////////////////////////////////
!    
!@mark -
!
#include "Includes.h"
!  **************
   Module Physics 
!  **************
!
      USE SMConstants
      USE PhysicsStorage
      use VariableConversion
      IMPLICIT NONE

      private
      public  InviscidFlux, ViscousFlux
      public  ViscousFlux0D, ViscousFlux2D, ViscousFlux3D
      public  InviscidFlux0D, InviscidFlux3D
      public  InviscidJacobian
      public  getStressTensor, SutherlandsLaw
!
!     ---------
!     Constants
!     ---------
!
      INTEGER, PARAMETER   :: WALL_BC = 1, RADIATION_BC = 2
      REAL(KIND=RP)        :: waveSpeed
      INTEGER              :: boundaryCondition(4), bcType

     interface InviscidFlux
         module procedure InviscidFlux0D, InviscidFlux3D
     end interface InviscidFlux

     interface ViscousFlux
         module procedure ViscousFlux0D, ViscousFlux2D, ViscousFlux3D
         module procedure ViscousFlux0D_withSGS, ViscousFlux3D_withSGS
     end interface ViscousFlux
!
!     ========
      CONTAINS 
!     ========
!
!     
!
!//////////////////////////////////////////////////////////////////////////////
!
!           INVISCID FLUXES
!           ---------------   
!
!//////////////////////////////////////////////////////////////////////////////
!
      pure subroutine InviscidFlux0D(Q, F)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)           :: u , v , w , p

         associate ( gammaMinus1 => thermodynamics % gammaMinus1 ) 

         u = Q(IRHOU) / Q(IRHO)
         v = Q(IRHOV) / Q(IRHO)
         w = Q(IRHOW) / Q(IRHO)
         p = gammaMinus1 * (Q(IRHOE) - 0.5_RP * ( Q(IRHOU) * u + Q(IRHOV) * v + Q(IRHOW) * w ) )
!
!        X-Flux
!        ------         
         F(IRHO , IX ) = Q(IRHOU)
         F(IRHOU, IX ) = Q(IRHOU) * u + p
         F(IRHOV, IX ) = Q(IRHOU) * v
         F(IRHOW, IX ) = Q(IRHOU) * w
         F(IRHOE, IX ) = ( Q(IRHOE) + p ) * u
!
!        Y-Flux
!        ------
         F(IRHO , IY ) = Q(IRHOV)
         F(IRHOU ,IY ) = F(IRHOV,IX)
         F(IRHOV ,IY ) = Q(IRHOV) * v + p
         F(IRHOW ,IY ) = Q(IRHOV) * w
         F(IRHOE ,IY ) = ( Q(IRHOE) + p ) * v
!
!        Z-Flux
!        ------
         F(IRHO ,IZ) = Q(IRHOW)
         F(IRHOU,IZ) = F(IRHOW,IX)
         F(IRHOV,IZ) = F(IRHOW,IY)
         F(IRHOW,IZ) = Q(IRHOW) * w + P
         F(IRHOE,IZ) = ( Q(IRHOE) + p ) * w
      
         end associate

      end subroutine InviscidFlux0D

      pure subroutine InviscidFlux3D(N, Q, F)
         implicit none
         integer,       intent(in)  :: N(3)
         real(kind=RP), intent(in)  :: Q(1:NCONS,0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(out) :: F(1:NCONS,0:N(1),0:N(2),0:N(3),1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                 :: i, j, k
         real(kind=RP)           :: u(0:N(1),0:N(2),0:N(3)) , v(0:N(1),0:N(2),0:N(3)) , w(0:N(1),0:N(2),0:N(3)) , p(0:N(1),0:N(2),0:N(3))

         associate ( gammaMinus1 => thermodynamics % gammaMinus1 ) 

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            u(i,j,k) = Q(IRHOU,i,j,k) / Q(IRHO,i,j,k)
            v(i,j,k) = Q(IRHOV,i,j,k) / Q(IRHO,i,j,k)
            w(i,j,k) = Q(IRHOW,i,j,k) / Q(IRHO,i,j,k)
            p(i,j,k) = gammaMinus1 * (Q(IRHOE,i,j,k) - 0.5_RP * ( Q(IRHOU,i,j,k) * u(i,j,k) + Q(IRHOV,i,j,k) * v(i,j,k) + Q(IRHOW,i,j,k) * w(i,j,k) ) )
            
            F(IRHO,i,j,k , IX ) = Q(IRHOU,i,j,k)
            F(IRHOU,i,j,k, IX ) = Q(IRHOU,i,j,k) * u(i,j,k) + p(i,j,k)
            F(IRHOV,i,j,k, IX ) = Q(IRHOU,i,j,k) * v(i,j,k)
            F(IRHOW,i,j,k, IX ) = Q(IRHOU,i,j,k) * w(i,j,k)
            F(IRHOE,i,j,k, IX ) = ( Q(IRHOE,i,j,k) + p(i,j,k) ) * u(i,j,k)

         end do   ; end do          ; end do
   
         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            F(IRHO,i,j,k , IY ) = Q(IRHOV,i,j,k)
            F(IRHOU,i,j,k ,IY ) = Q(IRHOU,i,j,k) * v(i,j,k)
            F(IRHOV,i,j,k ,IY ) = Q(IRHOV,i,j,k) * v(i,j,k) + p(i,j,k)
            F(IRHOW,i,j,k ,IY ) = Q(IRHOV,i,j,k) * w(i,j,k)
            F(IRHOE,i,j,k ,IY ) = ( Q(IRHOE,i,j,k) + p(i,j,k) ) * v(i,j,k)
         end do   ; end do          ; end do
   
         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            F(IRHO,i,j,k ,IZ) = Q(IRHOW,i,j,k)
            F(IRHOU,i,j,k,IZ) = Q(IRHOW,i,j,k) * u(i,j,k)
            F(IRHOV,i,j,k,IZ) = Q(IRHOW,i,j,k) * v(i,j,k)
            F(IRHOW,i,j,k,IZ) = Q(IRHOW,i,j,k) * w(i,j,k) + p(i,j,k)
            F(IRHOE,i,j,k,IZ) = ( Q(IRHOE,i,j,k) + p(i,j,k) ) * w(i,j,k)
         end do   ; end do          ; end do

         end associate

      end subroutine InviscidFlux3D
!
!     -------------------------------------------------------------------------------
!     Subroutine for computing the Jacobian of the inviscid flux when it has the form 
!
!        F = f*iHat + g*jHat + h*kHat
!
!     First index indicates the flux term and second index indicates the conserved 
!     variable term. 
!     ***** This routine is necessary for computing the analytical Jacobian. *****
!     -------------------------------------------------------------------------------
      pure subroutine InviscidJacobian(q,dfdq,dgdq,dhdq)
         implicit none
         !-------------------------------------------------
         real(kind=RP), intent (in)  :: q(NCONS)
         real(kind=RP), intent (out) :: dfdq(NCONS,NCONS)
         real(kind=RP), intent (out) :: dgdq(NCONS,NCONS)
         real(kind=RP), intent (out) :: dhdq(NCONS,NCONS)
         !-------------------------------------------------
         real(kind=RP)  :: u,v,w ! Velocity components
         real(kind=RP)  :: V2    ! Total velocity squared
         real(kind=RP)  :: p     ! Pressure
         real(kind=RP)  :: H     ! Total enthalpy
         !-------------------------------------------------
         
         associate( gammaMinus1 => thermodynamics % gammaMinus1, & 
                    gamma => thermodynamics % gamma )
         
         u  = q(IRHOU) / q(IRHO)
         v  = q(IRHOV) / q(IRHO)
         w  = q(IRHOW) / q(IRHO)
         V2 = u*u + v*v + w*w
         p  = Pressure(q)
         H  = (q(IRHOE) + p) / q(IRHO)
!
!        Flux in the x direction (f)
!        ---------------------------

         dfdq(1,1) = 0._RP
         dfdq(1,2) = 1._RP
         dfdq(1,3) = 0._RP
         dfdq(1,4) = 0._RP
         dfdq(1,5) = 0._RP
         
         dfdq(2,1) = -u*u + 0.5_RP*gammaMinus1*V2
         dfdq(2,2) = (3._RP - gamma) * u
         dfdq(2,3) = -gammaMinus1 * v
         dfdq(2,4) = -gammaMinus1 * w
         dfdq(2,5) = gammaMinus1
         
         dfdq(3,1) = -u*v
         dfdq(3,2) = v
         dfdq(3,3) = u
         dfdq(3,4) = 0._RP
         dfdq(3,5) = 0._RP
         
         dfdq(4,1) = -u*w
         dfdq(4,2) = w
         dfdq(4,3) = 0._RP
         dfdq(4,4) = u
         dfdq(4,5) = 0._RP
         
         dfdq(5,1) = u * (0.5_RP*gammaMinus1*V2 - H)
         dfdq(5,2) = H - gammaMinus1 * u*u
         dfdq(5,3) = -gammaMinus1 * u*v
         dfdq(5,4) = -gammaMinus1 * u*w
         dfdq(5,5) = gamma * u
         
!
!        Flux in the y direction (g)
!        ---------------------------
         
         dgdq(1,1) = 0._RP
         dgdq(1,2) = 0._RP
         dgdq(1,3) = 1._RP
         dgdq(1,4) = 0._RP
         dgdq(1,5) = 0._RP
         
         dgdq(2,1) = -u*v
         dgdq(2,2) = v
         dgdq(2,3) = u
         dgdq(2,4) = 0._RP
         dgdq(2,5) = 0._RP
         
         dgdq(3,1) = -v*v + 0.5_RP*gammaMinus1*V2
         dgdq(3,2) = -gammaMinus1 * u
         dgdq(3,3) = (3._RP - gamma) * v
         dgdq(3,4) = -gammaMinus1 * w
         dgdq(3,5) = gammaMinus1
         
         dgdq(4,1) = -v*w
         dgdq(4,2) = 0._RP
         dgdq(4,3) = w
         dgdq(4,4) = v
         dgdq(4,5) = 0._RP
         
         dgdq(5,1) = v * (0.5_RP*gammaMinus1*V2 - H)
         dgdq(5,2) = -gammaMinus1 * u*v
         dgdq(5,3) = H - gammaMinus1 * v*v
         dgdq(5,4) = -gammaMinus1 * v*w
         dgdq(5,5) = gamma * v
!
!        Flux in the z direction (h)
!        ---------------------------
         
         dhdq(1,1) = 0._RP
         dhdq(1,2) = 0._RP
         dhdq(1,3) = 0._RP
         dhdq(1,4) = 1._RP
         dhdq(1,5) = 0._RP
         
         dhdq(2,1) = -u*w
         dhdq(2,2) = w
         dhdq(2,3) = 0._RP
         dhdq(2,4) = u
         dhdq(2,5) = 0._RP
         
         dhdq(3,1) = -v*w
         dhdq(3,2) = 0._RP
         dhdq(3,3) = w
         dhdq(3,4) = v
         dhdq(3,5) = 0._RP
         
         dhdq(4,1) = -w*w + 0.5_RP*gammaMinus1*V2
         dhdq(4,2) = -gammaMinus1 * u
         dhdq(4,3) = -gammaMinus1 * v
         dhdq(4,4) = (3._RP - gamma) * w
         dhdq(4,5) = gammaMinus1
         
         dhdq(5,1) = w * (0.5_RP*gammaMinus1*V2 - H)
         dhdq(5,2) = -gammaMinus1 * u*w
         dhdq(5,3) = -gammaMinus1 * v*w
         dhdq(5,4) = H - gammaMinus1 * w*w
         dhdq(5,5) = gamma * w
         
         end associate
         
      end subroutine InviscidJacobian
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!>        VISCOUS FLUXES
!         --------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine ViscousFlux0D(Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         real(kind=RP), intent(in)  :: Q   (1:NCONS     )
         real(kind=RP), intent(in)  :: U_x (1:N_GRAD_EQN)
         real(kind=RP), intent(in)  :: U_y (1:N_GRAD_EQN)
         real(kind=RP), intent(in)  :: U_z (1:N_GRAD_EQN)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:NCONS, 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                    :: T , sutherLaw
         real(kind=RP)                    :: divV
         real(kind=RP)                    :: u , v , w

         u = Q(IRHOU) / Q(IRHO)
         v = Q(IRHOV) / Q(IRHO)
         w = Q(IRHOW) / Q(IRHO)

         T     = Temperature(Q)
         sutherLaw = SutherlandsLaw(T)

         divV = U_x(IGU) + U_y(IGV) + U_z(IGW)

         F(IRHO,IX)  = 0.0_RP
         F(IRHOU,IX) = mu * sutherLaw * (2.0_RP * U_x(IGU) - 2.0_RP/3.0_RP * divV ) 
         F(IRHOV,IX) = mu * sutherLaw * ( U_x(IGV) + U_y(IGU) ) 
         F(IRHOW,IX) = mu * sutherLaw * ( U_x(IGW) + U_z(IGU) ) 
         F(IRHOE,IX) = F(IRHOU,IX) * u + F(IRHOV,IX) * v + F(IRHOW,IX) * w + kappa * sutherLaw * U_x(IGT) 

         F(IRHO,IY) = 0.0_RP
         F(IRHOU,IY) = F(IRHOV,IX) 
         F(IRHOV,IY) = mu * sutherLaw * (2.0_RP * U_y(IGV) - 2.0_RP / 3.0_RP * divV ) 
         F(IRHOW,IY) = mu * sutherLaw * ( U_y(IGW) + U_z(IGV) ) 
         F(IRHOE,IY) = F(IRHOU,IY) * u + F(IRHOV,IY) * v + F(IRHOW,IY) * w + kappa * sutherLaw * U_y(IGT) 

         F(IRHO,IZ) = 0.0_RP
         F(IRHOU,IZ) = F(IRHOW,IX) 
         F(IRHOV,IZ) = F(IRHOW,IY) 
         F(IRHOW,IZ) = mu * sutherLaw * ( 2.0_RP * U_z(IGW) - 2.0_RP / 3.0_RP * divV ) 
         F(IRHOE,IZ) = F(IRHOU,IZ) * u + F(IRHOV,IZ) * v + F(IRHOW,IZ) * w + kappa * sutherLaw *U_z(IGT)

         ! with Pr = constant, dmudx = dkappadx
      end subroutine ViscousFlux0D

      pure subroutine ViscousFlux2D( N, Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         integer         , intent(in)  :: N(2)
         real(kind=RP),    intent(in)  :: Q  (1:NCONS, 0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: U_x(1:N_GRAD_EQN, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_y(1:N_GRAD_EQN, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_z(1:N_GRAD_EQN, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2))
         real(kind=RP),    intent(out) :: F   (1:NCONS, 1:NDIM, 0:N(1), 0:N(2))
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: T(0:N(1),0:N(2)) , sutherLaw(0:N(1),0:N(2))
         real(kind=RP) :: divV(0:N(1),0:N(2))
         real(kind=RP) :: u(0:N(1),0:N(2)) , v(0:N(1),0:N(2)) , w(0:N(1),0:N(2))
         integer       :: i , j , k

         associate( gammaM2 => dimensionless % gammaM2, &
                    gammaMinus1 => thermodynamics % gammaMinus1 ) 

         do j = 0, N(2) ; do i = 0, N(1)
            u(i,j) = Q(IRHOU,i,j) / Q(IRHO,i,j)
            v(i,j) = Q(IRHOV,i,j) / Q(IRHO,i,j)
            w(i,j) = Q(IRHOW,i,j) / Q(IRHO,i,j)
   
   
            T(i,j) = gammaM2 * gammaMinus1 * (Q(IRHOE,i,j)  & 
                   - 0.5_RP * ( Q(IRHOU,i,j) * u(i,j) + Q(IRHOV,i,j) * v(i,j) + Q(IRHOW,i,j) * w(i,j) ) ) / Q(IRHO,i,j)
   

            sutherLaw(i,j) = SutherlandsLaw(T(i,j))

            divV(i,j) = U_x(IGU,i,j) + U_y(IGV,i,j) + U_z(IGW,i,j)
   
            F(IRHO ,IX,i,j) = 0.0_RP
            F(IRHOU,IX,i,j) = mu(i,j) * sutherLaw(i,j) * (2.0_RP * U_x(IGU,i,j) - 2.0_RP/3.0_RP * divV(i,j) ) 
            F(IRHOV,IX,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_x(IGV,i,j) + U_y(IGU,i,j) ) 
            F(IRHOW,IX,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_x(IGW,i,j) + U_z(IGU,i,j) ) 
            F(IRHOE,IX,i,j) = F(IRHOU,IX,i,j) * u(i,j) + F(IRHOV,IX,i,j) * v(i,j) + F(IRHOW,IX,i,j) * w(i,j) &
                  + sutherLaw(i,j) * kappa(i,j) * U_x(IGT,i,j) 
   
            F(IRHO, IY,i,j) = 0.0_RP
            F(IRHOU,IY,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_x(IGV,i,j) + U_y(IGU,i,j) )  
            F(IRHOV,IY,i,j) = mu(i,j) * sutherLaw(i,j) * (2.0_RP * U_y(IGV,i,j) - 2.0_RP / 3.0_RP * divV(i,j) ) 
            F(IRHOW,IY,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_y(IGW,i,j) + U_z(IGV,i,j) ) 
            F(IRHOE,IY,i,j) = F(IRHOU,IY,i,j) * u(i,j) + F(IRHOV,IY,i,j) * v(i,j) + F(IRHOW,IY,i,j) * w(i,j) &
                  + sutherLaw(i,j) * kappa(i,j) * U_y(IGT,i,j) 
   
            F(IRHO, IZ,i,j ) = 0.0_RP
            F(IRHOU,IZ,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_x(IGW,i,j) + U_z(IGU,i,j) ) 
            F(IRHOV,IZ,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_y(IGW,i,j) + U_z(IGV,i,j) ) 
            F(IRHOW,IZ,i,j) = mu(i,j) * sutherLaw(i,j) * ( 2.0_RP * U_z(IGW,i,j) - 2.0_RP / 3.0_RP * divV(i,j) ) 
            F(IRHOE,IZ,i,j) = F(IRHOU,IZ,i,j) * u(i,j) + F(IRHOV,IZ,i,j) * v(i,j) + F(IRHOW,IZ,i,j) * w(i,j) &
                  + sutherLaw(i,j) * kappa(i,j) * U_z(IGT,i,j) 
   
         end do    ; end do

         end associate

      end subroutine ViscousFlux2D

      pure subroutine ViscousFlux3D( N, Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         integer         , intent(in)  :: N(3)
         real(kind=RP),    intent(in)  :: Q  (1:NCONS, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: U_x(1:N_GRAD_EQN, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_y(1:N_GRAD_EQN, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_z(1:N_GRAD_EQN, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(out) :: F   (1:NCONS, 0:N(1), 0:N(2), 0:N(3), 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: T(0:N(1),0:N(2),0:N(3)) , sutherLaw(0:N(1),0:N(2),0:N(3))
         real(kind=RP) :: divV(0:N(1),0:N(2),0:N(3))
         real(kind=RP) :: u(0:N(1),0:N(2),0:N(3)) , v(0:N(1),0:N(2),0:N(3)) , w(0:N(1),0:N(2),0:N(3))
         integer       :: i , j , k

         associate( gammaM2 => dimensionless % gammaM2, &
                    gammaMinus1 => thermodynamics % gammaMinus1 ) 

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            u(i,j,k) = Q(IRHOU,i,j,k) / Q(IRHO,i,j,k)
            v(i,j,k) = Q(IRHOV,i,j,k) / Q(IRHO,i,j,k)
            w(i,j,k) = Q(IRHOW,i,j,k) / Q(IRHO,i,j,k)
   
   
            T(i,j,k) = gammaM2 * gammaMinus1 * (Q(IRHOE,i,j,k)  & 
                   - 0.5_RP * ( Q(IRHOU,i,j,k) * u(i,j,k) + Q(IRHOV,i,j,k) * v(i,j,k) + Q(IRHOW,i,j,k) * w(i,j,k) ) ) / Q(IRHO,i,j,k)
   

            sutherLaw(i,j,k) = SutherlandsLaw(T(i,j,k))

            divV(i,j,k) = U_x(IGU,i,j,k) + U_y(IGV,i,j,k) + U_z(IGW,i,j,k)
   
            F(IRHO,i,j,k ,IX) = 0.0_RP
            F(IRHOU,i,j,k,IX) = mu(i,j,k) * sutherLaw(i,j,k) * (2.0_RP * U_x(IGU,i,j,k) - 2.0_RP/3.0_RP * divV(i,j,k) ) 
            F(IRHOV,i,j,k,IX) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_x(IGV,i,j,k) + U_y(IGU,i,j,k) ) 
            F(IRHOW,i,j,k,IX) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_x(IGW,i,j,k) + U_z(IGU,i,j,k) ) 
            F(IRHOE,i,j,k,IX) = F(IRHOU,i,j,k,IX) * u(i,j,k) + F(IRHOV,i,j,k,IX) * v(i,j,k) + F(IRHOW,i,j,k,IX) * w(i,j,k) &
                  + sutherLaw(i,j,k) * kappa(i,j,k) * U_x(IGT,i,j,k) 
         end do      ; end do    ; end do

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            F(IRHO,i,j,k ,IY) = 0.0_RP
            F(IRHOU,i,j,k,IY) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_x(IGV,i,j,k) + U_y(IGU,i,j,k) )  
            F(IRHOV,i,j,k,IY) = mu(i,j,k) * sutherLaw(i,j,k) * (2.0_RP * U_y(IGV,i,j,k) - 2.0_RP / 3.0_RP * divV(i,j,k) ) 
            F(IRHOW,i,j,k,IY) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_y(IGW,i,j,k) + U_z(IGV,i,j,k) ) 
            F(IRHOE,i,j,k,IY) = F(IRHOU,i,j,k,IY) * u(i,j,k) + F(IRHOV,i,j,k,IY) * v(i,j,k) + F(IRHOW,i,j,k,IY) * w(i,j,k) &
                  + sutherLaw(i,j,k) * kappa(i,j,k) * U_y(IGT,i,j,k) 
         end do      ; end do    ; end do

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            F(IRHO,i,j,k,IZ ) = 0.0_RP
            F(IRHOU,i,j,k,IZ) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_x(IGW,i,j,k) + U_z(IGU,i,j,k) ) 
            F(IRHOV,i,j,k,IZ) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_y(IGW,i,j,k) + U_z(IGV,i,j,k) ) 
            F(IRHOW,i,j,k,IZ) = mu(i,j,k) * sutherLaw(i,j,k) * ( 2.0_RP * U_z(IGW,i,j,k) - 2.0_RP / 3.0_RP * divV(i,j,k) ) 
            F(IRHOE,i,j,k,IZ) = F(IRHOU,i,j,k,IZ) * u(i,j,k) + F(IRHOV,i,j,k,IZ) * v(i,j,k) + F(IRHOW,i,j,k,IZ) * w(i,j,k) &
                  + sutherLaw(i,j,k) * kappa(i,j,k) * U_z(IGT,i,j,k) 
         end do      ; end do    ; end do

         end associate

      end subroutine ViscousFlux3D

      pure subroutine ViscousFlux0D_withSGS(Q, U_x, U_y, U_z, mu, kappa, tauSGS, qSGS, F)
         implicit none
         real(kind=RP), intent(in)  :: Q   (1:NCONS     )
         real(kind=RP), intent(in)  :: U_x (1:N_GRAD_EQN)
         real(kind=RP), intent(in)  :: U_y (1:N_GRAD_EQN)
         real(kind=RP), intent(in)  :: U_z (1:N_GRAD_EQN)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(in)  :: tauSGS(NDIM, NDIM)
         real(kind=RP), intent(in)  :: qSGS(NDIM)
         real(kind=RP), intent(out) :: F(1:NCONS, 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                    :: T , sutherLaw
         real(kind=RP)                    :: divV
         real(kind=RP)                    :: u , v , w

         u = Q(IRHOU) / Q(IRHO)
         v = Q(IRHOV) / Q(IRHO)
         w = Q(IRHOW) / Q(IRHO)

         T     = Temperature(Q)
         sutherLaw = SutherlandsLaw(T)

         divV = U_x(IGU) + U_y(IGV) + U_z(IGW)

         F(IRHO,IX)  = 0.0_RP
         F(IRHOU,IX) = mu * sutherLaw * (2.0_RP * U_x(IGU) - 2.0_RP/3.0_RP * divV ) - tauSGS(1,1)
         F(IRHOV,IX) = mu * sutherLaw * ( U_x(IGV) + U_y(IGU) ) - tauSGS(2,1)
         F(IRHOW,IX) = mu * sutherLaw * ( U_x(IGW) + U_z(IGU) ) - tauSGS(3,1)
         F(IRHOE,IX) = F(IRHOU,IX) * u + F(IRHOV,IX) * v + F(IRHOW,IX) * w + kappa * sutherLaw * U_x(IGT) - qSGS(1)

         F(IRHO,IY) = 0.0_RP
         F(IRHOU,IY) = F(IRHOV,IX) 
         F(IRHOV,IY) = mu * sutherLaw * (2.0_RP * U_y(IGV) - 2.0_RP / 3.0_RP * divV ) - tauSGS(2,2)
         F(IRHOW,IY) = mu * sutherLaw * ( U_y(IGW) + U_z(IGV) ) - tauSGS(3,2)
         F(IRHOE,IY) = F(IRHOU,IY) * u + F(IRHOV,IY) * v + F(IRHOW,IY) * w + kappa * sutherLaw * U_y(IGT) - qSGS(2)

         F(IRHO,IZ) = 0.0_RP
         F(IRHOU,IZ) = F(IRHOW,IX) 
         F(IRHOV,IZ) = F(IRHOW,IY) 
         F(IRHOW,IZ) = mu * sutherLaw * ( 2.0_RP * U_z(IGW) - 2.0_RP / 3.0_RP * divV ) - tauSGS(3,3)
         F(IRHOE,IZ) = F(IRHOU,IZ) * u + F(IRHOV,IZ) * v + F(IRHOW,IZ) * w + kappa * sutherLaw *U_z(IGT) - qSGS(3)

      end subroutine ViscousFlux0D_withSGS

      pure subroutine ViscousFlux2D_withSGS( N, Q, U_x, U_y, U_z, mu, kappa, tauSGS, qSGS, F)
         implicit none
         integer         , intent(in)  :: N(2)
         real(kind=RP),    intent(in)  :: Q  (1:NCONS, 0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: U_x(1:N_GRAD_EQN, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_y(1:N_GRAD_EQN, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_z(1:N_GRAD_EQN, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: tauSGS(1:NDIM, 1:NDIM, 0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: qSGS(1:NDIM, 0:N(1), 0:N(2))
         real(kind=RP),    intent(out) :: F   (1:NCONS, 1:NDIM, 0:N(1), 0:N(2))
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: T(0:N(1),0:N(2)) , sutherLaw(0:N(1),0:N(2))
         real(kind=RP) :: divV(0:N(1),0:N(2))
         real(kind=RP) :: u(0:N(1),0:N(2)) , v(0:N(1),0:N(2)) , w(0:N(1),0:N(2))
         integer       :: i , j , k

         associate( gammaM2 => dimensionless % gammaM2, &
                    gammaMinus1 => thermodynamics % gammaMinus1 ) 

         do j = 0, N(2) ; do i = 0, N(1)
            u(i,j) = Q(IRHOU,i,j) / Q(IRHO,i,j)
            v(i,j) = Q(IRHOV,i,j) / Q(IRHO,i,j)
            w(i,j) = Q(IRHOW,i,j) / Q(IRHO,i,j)
   
   
            T(i,j) = gammaM2 * gammaMinus1 * (Q(IRHOE,i,j)  & 
                   - 0.5_RP * ( Q(IRHOU,i,j) * u(i,j) + Q(IRHOV,i,j) * v(i,j) + Q(IRHOW,i,j) * w(i,j) ) ) / Q(IRHO,i,j)
   

            sutherLaw(i,j) = SutherlandsLaw(T(i,j))

            divV(i,j) = U_x(IGU,i,j) + U_y(IGV,i,j) + U_z(IGW,i,j)
   
            F(IRHO ,IX,i,j) = 0.0_RP
            F(IRHOU,IX,i,j) = mu(i,j) * sutherLaw(i,j) * (2.0_RP * U_x(IGU,i,j) - 2.0_RP/3.0_RP * divV(i,j) ) - tauSGS(1,1,i,j)
            F(IRHOV,IX,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_x(IGV,i,j) + U_y(IGU,i,j) ) - tauSGS(2,1,i,j)
            F(IRHOW,IX,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_x(IGW,i,j) + U_z(IGU,i,j) ) - tauSGS(3,1,i,j)
            F(IRHOE,IX,i,j) = F(IRHOU,IX,i,j) * u(i,j) + F(IRHOV,IX,i,j) * v(i,j) + F(IRHOW,IX,i,j) * w(i,j) &
                  + sutherLaw(i,j) * kappa(i,j) * U_x(IGT,i,j) - qSGS(1,i,j)
   
            F(IRHO, IY,i,j) = 0.0_RP
            F(IRHOU,IY,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_x(IGV,i,j) + U_y(IGU,i,j) )  - tauSGS(1,2,i,j)
            F(IRHOV,IY,i,j) = mu(i,j) * sutherLaw(i,j) * (2.0_RP * U_y(IGV,i,j) - 2.0_RP / 3.0_RP * divV(i,j) ) - tauSGS(2,2,i,j)
            F(IRHOW,IY,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_y(IGW,i,j) + U_z(IGV,i,j) ) - tauSGS(3,2,i,j)
            F(IRHOE,IY,i,j) = F(IRHOU,IY,i,j) * u(i,j) + F(IRHOV,IY,i,j) * v(i,j) + F(IRHOW,IY,i,j) * w(i,j) &
                  + sutherLaw(i,j) * kappa(i,j) * U_y(IGT,i,j) - qSGS(2,i,j)
   
            F(IRHO, IZ,i,j ) = 0.0_RP
            F(IRHOU,IZ,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_x(IGW,i,j) + U_z(IGU,i,j) ) - tauSGS(1,3,i,j)
            F(IRHOV,IZ,i,j) = mu(i,j) * sutherLaw(i,j) * ( U_y(IGW,i,j) + U_z(IGV,i,j) ) - tauSGS(2,3,i,j)
            F(IRHOW,IZ,i,j) = mu(i,j) * sutherLaw(i,j) * ( 2.0_RP * U_z(IGW,i,j) - 2.0_RP / 3.0_RP * divV(i,j) ) - tauSGS(3,3,i,j)
            F(IRHOE,IZ,i,j) = F(IRHOU,IZ,i,j) * u(i,j) + F(IRHOV,IZ,i,j) * v(i,j) + F(IRHOW,IZ,i,j) * w(i,j) &
                  + sutherLaw(i,j) * kappa(i,j) * U_z(IGT,i,j) - qSGS(3,i,j)
   
         end do    ; end do

         end associate

      end subroutine ViscousFlux2D_withSGS

      pure subroutine ViscousFlux3D_withSGS( N, Q, U_x, U_y, U_z, mu, kappa, tauSGS, qSGS, F)
         implicit none
         integer         , intent(in)  :: N(3)
         real(kind=RP),    intent(in)  :: Q  (1:NCONS, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: U_x(1:N_GRAD_EQN, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_y(1:N_GRAD_EQN, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_z(1:N_GRAD_EQN, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: tauSGS(1:NDIM, 1:NDIM, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: qSGS(1:NDIM, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(out) :: F   (1:NCONS, 0:N(1), 0:N(2), 0:N(3), 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: T(0:N(1),0:N(2),0:N(3)) , sutherLaw(0:N(1),0:N(2),0:N(3))
         real(kind=RP) :: divV(0:N(1),0:N(2),0:N(3))
         real(kind=RP) :: u(0:N(1),0:N(2),0:N(3)) , v(0:N(1),0:N(2),0:N(3)) , w(0:N(1),0:N(2),0:N(3))
         integer       :: i , j , k

         associate( gammaM2 => dimensionless % gammaM2, &
                    gammaMinus1 => thermodynamics % gammaMinus1 ) 

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            u(i,j,k) = Q(IRHOU,i,j,k) / Q(IRHO,i,j,k)
            v(i,j,k) = Q(IRHOV,i,j,k) / Q(IRHO,i,j,k)
            w(i,j,k) = Q(IRHOW,i,j,k) / Q(IRHO,i,j,k)
   
   
            T(i,j,k) = gammaM2 * gammaMinus1 * (Q(IRHOE,i,j,k)  & 
                   - 0.5_RP * ( Q(IRHOU,i,j,k) * u(i,j,k) + Q(IRHOV,i,j,k) * v(i,j,k) + Q(IRHOW,i,j,k) * w(i,j,k) ) ) / Q(IRHO,i,j,k)
   

            sutherLaw(i,j,k) = SutherlandsLaw(T(i,j,k))

            divV(i,j,k) = U_x(IGU,i,j,k) + U_y(IGV,i,j,k) + U_z(IGW,i,j,k)
   
            F(IRHO,i,j,k ,IX) = 0.0_RP
            F(IRHOU,i,j,k,IX) = mu(i,j,k) * sutherLaw(i,j,k) * (2.0_RP * U_x(IGU,i,j,k) - 2.0_RP/3.0_RP * divV(i,j,k) ) - tauSGS(1,1,i,j,k)
            F(IRHOV,i,j,k,IX) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_x(IGV,i,j,k) + U_y(IGU,i,j,k) ) - tauSGS(2,1,i,j,k)
            F(IRHOW,i,j,k,IX) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_x(IGW,i,j,k) + U_z(IGU,i,j,k) ) - tauSGS(3,1,i,j,k)
            F(IRHOE,i,j,k,IX) = F(IRHOU,i,j,k,IX) * u(i,j,k) + F(IRHOV,i,j,k,IX) * v(i,j,k) + F(IRHOW,i,j,k,IX) * w(i,j,k) &
                  + sutherLaw(i,j,k) * kappa(i,j,k) * U_x(IGT,i,j,k) - qSGS(1,i,j,k)
   
         end do      ; end do    ; end do

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            F(IRHO,i,j,k ,IY) = 0.0_RP
            F(IRHOU,i,j,k,IY) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_x(IGV,i,j,k) + U_y(IGU,i,j,k) )  - tauSGS(1,2,i,j,k)
            F(IRHOV,i,j,k,IY) = mu(i,j,k) * sutherLaw(i,j,k) * (2.0_RP * U_y(IGV,i,j,k) - 2.0_RP / 3.0_RP * divV(i,j,k) ) - tauSGS(2,2,i,j,k)
            F(IRHOW,i,j,k,IY) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_y(IGW,i,j,k) + U_z(IGV,i,j,k) ) - tauSGS(3,2,i,j,k)
            F(IRHOE,i,j,k,IY) = F(IRHOU,i,j,k,IY) * u(i,j,k) + F(IRHOV,i,j,k,IY) * v(i,j,k) + F(IRHOW,i,j,k,IY) * w(i,j,k) &
                  + sutherLaw(i,j,k) * kappa(i,j,k) * U_y(IGT,i,j,k) - qSGS(2,i,j,k)
   
         end do      ; end do    ; end do

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            F(IRHO,i,j,k,IZ ) = 0.0_RP
            F(IRHOU,i,j,k,IZ) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_x(IGW,i,j,k) + U_z(IGU,i,j,k) ) - tauSGS(1,3,i,j,k)
            F(IRHOV,i,j,k,IZ) = mu(i,j,k) * sutherLaw(i,j,k) * ( U_y(IGW,i,j,k) + U_z(IGV,i,j,k) ) - tauSGS(2,3,i,j,k)
            F(IRHOW,i,j,k,IZ) = mu(i,j,k) * sutherLaw(i,j,k) * ( 2.0_RP * U_z(IGW,i,j,k) - 2.0_RP / 3.0_RP * divV(i,j,k) ) - tauSGS(3,3,i,j,k)
            F(IRHOE,i,j,k,IZ) = F(IRHOU,i,j,k,IZ) * u(i,j,k) + F(IRHOV,i,j,k,IZ) * v(i,j,k) + F(IRHOW,i,j,k,IZ) * w(i,j,k) &
                  + sutherLaw(i,j,k) * kappa(i,j,k) * U_z(IGT,i,j,k) - qSGS(3,i,j,k)
   
         end do      ; end do    ; end do

         end associate

      end subroutine ViscousFlux3D_withSGS

!
!---------------------------------------------------------------------
!! Compute the molecular diffusivity by way of Sutherland's law
!---------------------------------------------------------------------
!
      PURE FUNCTION SutherlandsLaw(T) RESULT(mu)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), INTENT(IN) :: T !! The temperature
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: mu !! The diffusivity
!      
      mu = (1._RP + tRatio)/(T + tRatio)*T*SQRT(T)


      END FUNCTION SutherlandsLaw

      pure subroutine getStressTensor(Q,U_x,U_y,U_z,tau)
         implicit none
         real(kind=RP), intent(in)      :: Q   (1:NCONS         )
         real(kind=RP), intent(in)      :: U_x (1:N_GRAD_EQN    )
         real(kind=RP), intent(in)      :: U_y (1:N_GRAD_EQN    )
         real(kind=RP), intent(in)      :: U_z (1:N_GRAD_EQN    )
         real(kind=RP), intent(out)     :: tau (1:NDIM, 1:NDIM   )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: T , muOfT
         real(kind=RP) :: divV

         associate ( mu0 => dimensionless % mu )

         T     = Temperature(Q)
         muOfT = SutherlandsLaw(T)

         divV = U_x(IGU) + U_y(IGV) + U_z(IGW)

         tau(IX,IX) = mu0 * muOfT * (2.0_RP * U_x(IGU) - 2.0_RP/3.0_RP * divV )
         tau(IY,IX) = mu0 * muOfT * ( U_x(IGV) + U_y(IGU) ) 
         tau(IZ,IX) = mu0 * muOfT * ( U_x(IGW) + U_z(IGU) ) 
         tau(IX,IY) = tau(IY,IX)
         tau(IY,IY) = mu0 * muOfT * (2.0_RP * U_y(IGV) - 2.0_RP/3.0_RP * divV )
         tau(IZ,IY) = mu0 * muOfT * ( U_y(IGW) + U_z(IGV) ) 
         tau(IX,IZ) = tau(IZ,IX)
         tau(IY,IZ) = tau(IZ,IY)
         tau(IZ,IZ) = mu0 * muOfT * (2.0_RP * U_z(IGW) - 2.0_RP/3.0_RP * divV )

         end associate

      end subroutine getStressTensor
   END Module Physics
!@mark -
!
! /////////////////////////////////////////////////////////////////////
!
!----------------------------------------------------------------------
!! This routine returns the maximum eigenvalues for the Euler equations 
!! for the given solution value in each spatial direction. 
!! These are to be used to compute the local time step.
!----------------------------------------------------------------------
!
      SUBROUTINE ComputeEigenvaluesForState( Q, eigen )
      
      USE SMConstants
      USE PhysicsStorage
      USE VariableConversion, ONLY:Pressure
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=Rp), DIMENSION(N_EQN) :: Q
      REAL(KIND=Rp), DIMENSION(3)     :: eigen
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=Rp) :: u, v, w, p, a
!      
      associate ( gamma => thermodynamics % gamma ) 

      u = ABS( Q(2)/Q(1) )
      v = ABS( Q(3)/Q(1) )
      w = ABS( Q(4)/Q(1) )
      p = Pressure(Q)
      a = SQRT(gamma*p/Q(1))
      
      eigen(1) = u + a
      eigen(2) = v + a
      eigen(3) = w + a

      end associate
      
      END SUBROUTINE ComputeEigenvaluesForState
