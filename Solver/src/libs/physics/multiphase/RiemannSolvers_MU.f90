!
#include "Includes.h"
module RiemannSolvers_MU
   use SMConstants
   use Physics_MU
   use PhysicsStorage_MU
   use VariableConversion_MU
   use FluidData_MU

   private 
   public RiemannSolver, SetRiemannSolver, ExactRiemannSolver

   abstract interface
      subroutine RiemannSolverFCN(QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL,fR)
         use SMConstants
         use PhysicsStorage_MU
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: rhoL
         real(kind=RP), intent(in)       :: rhoR
         real(kind=RP), intent(in)       :: muL
         real(kind=RP), intent(in)       :: muR
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: fL(1:NCONS)
         real(kind=RP), intent(out)      :: fR(1:NCONS)
      end subroutine RiemannSolverFCN
   end interface

   procedure(RiemannSolverFCN)     , pointer  :: RiemannSolver      => NULL()

   contains
      SUBROUTINE SetRiemannSolver(which, splitType)
!
!        **************************************************************
!              This subroutine is to set which Riemann solver is used.
!           the user cannot decide amongst the averaging function.
!           It is automatically selected depending on which split
!           form is enabled.
!           The user can choose the dissipation type:
!              None (central), Roe, Lax-Friedrichs, Rusanov
!
!           And the dissipation intensity, with the lambda stabilization
!           parameter. By default it is set to 1 (whole dissipation),
!           instead for central fluxes, which is 0 (no dissipation).
!        **************************************************************
!
         IMPLICIT NONE
         integer, intent(in) :: which
         integer, intent(in) :: splitType
         
         
         select case ( which )
         case(RIEMANN_CENTRAL)
            RiemannSolver => CentralRiemannSolver

         case(RIEMANN_EXACT)
            RiemannSolver => ExactRiemannSolver

         case default
            print*, "Undefined choice of Riemann Solver."
            print*, "Options available are:"
            print*, "   * Central"
            print*, "   * Exact"
            errorMessage(STD_OUT)
            STOP
         end select
      
      END SUBROUTINE SetRiemannSolver
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!        Riemann solvers
!        ---------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine CentralRiemannSolver(QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL, fR)
         implicit none 
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: rhoL
         real(kind=RP), intent(in)       :: rhoR
         real(kind=RP), intent(in)       :: muL
         real(kind=RP), intent(in)       :: muR
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: fL(1:NCONS)
         real(kind=RP), intent(out)      :: fR(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: cL, uL, vL, wL, pL, invSqrtRhoL
         real(kind=RP) :: cR, uR, vR, wR, pR, invSqrtRhoR
!
!        Left state variables and fluxes
!        -------------------------------
         invSqrtRhoL = 1.0_RP / sqrt(rhoL)
         cL = QLeft(IMC)
         uL = invSqrtRhoL * (QLeft(IMSQRHOU) * nHat(1) + QLeft(IMSQRHOV) * nHat(2) + QLeft(IMSQRHOW) * nHat(3))
         vL = invSqrtRhoL * (QLeft(IMSQRHOU) * t1(1)   + QLeft(IMSQRHOV) * t1(2)   + QLeft(IMSQRHOW) * t1(3))
         wL = invSqrtRhoL * (QLeft(IMSQRHOU) * t2(1)   + QLeft(IMSQRHOV) * t2(2)   + QLeft(IMSQRHOW) * t2(3))
         pL = QLeft(IMP)

         fL(IMC)      = uL*cL
         fL(IMSQRHOU) = 0.5_RP*rhoL*uL*uL + pL
         fL(IMSQRHOV) = 0.5_RP*rhoL*uL*vL
         fL(IMSQRHOW) = 0.5_RP*rhoL*uL*wL
         fL(IMP)      = 0.0_RP

!
!        Right state variables and fluxes
!        --------------------------------
         invSqrtRhoR = 1.0_RP / sqrt(rhoR)
         cR = QRight(IMC)
         uR = invSqrtRhoR * (QRight(IMSQRHOU) * nHat(1) + QRight(IMSQRHOV) * nHat(2) + QRight(IMSQRHOW) * nHat(3))
         vR = invSqrtRhoR * (QRight(IMSQRHOU) * t1(1)   + QRight(IMSQRHOV) * t1(2)   + QRight(IMSQRHOW) * t1(3))
         wR = invSqrtRhoR * (QRight(IMSQRHOU) * t2(1)   + QRight(IMSQRHOV) * t2(2)   + QRight(IMSQRHOW) * t2(3))
         pR = QRight(IMP)

         fR(IMC)      = uR*cR
         fR(IMSQRHOU) = 0.5_RP*rhoR*uR*uR + pR
         fR(IMSQRHOV) = 0.5_RP*rhoR*uR*vR
         fR(IMSQRHOW) = 0.5_RP*rhoR*uR*wR
         fR(IMP)      = 0.0_RP
!
!        Perform the average and rotation
!        --------------------------------
         fL = 0.5_RP*(fL + fR)
         fR = fL
!
!        Add the non-conservative term
!        -----------------------------          
         fL(IMSQRHOU) = fL(IMSQRHOU) + 0.5_RP*cL*(muR-muL) + 0.25_RP*rhoL*uL*(uR-uL)
         fL(IMSQRHOV) = fL(IMSQRHOV) + 0.25_RP*rhoL*uL*(vR-vL)
         fL(IMSQRHOW) = fL(IMSQRHOW) + 0.25_RP*rhoL*uL*(wR-wL)
         fL(IMP)      = fL(IMP)      + 0.5_RP*dimensionless % invMa2*(uR-uL)

         fR(IMSQRHOU) = fR(IMSQRHOU) + 0.5_RP*cR*(muL-muR) + 0.25_RP*rhoR*uR*(uL-uR)
         fR(IMSQRHOV) = fR(IMSQRHOV) + 0.25_RP*rhoR*uR*(vL-vR)
         fR(IMSQRHOW) = fR(IMSQRHOW) + 0.25_RP*rhoR*uR*(wL-wR)
         fR(IMP)      = fR(IMP)      + 0.5_RP*dimensionless % invMa2*(uL-uR)

         fL(IMSQRHOU:IMSQRHOW) = nHat*fL(IMSQRHOU) + t1*fL(IMSQRHOV) + t2*fL(IMSQRHOW)
         fR(IMSQRHOU:IMSQRHOW) = nHat*fR(IMSQRHOU) + t1*fR(IMSQRHOV) + t2*fR(IMSQRHOW)

      end subroutine CentralRiemannSolver

      subroutine ExactRiemannSolver(QLeft, QRight, rhoL, rhoR, muL, muR, nHat, t1, t2, fL, fR)
         implicit none 
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: rhoL
         real(kind=RP), intent(in)       :: rhoR
         real(kind=RP), intent(in)       :: muL
         real(kind=RP), intent(in)       :: muR
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: fL(1:NCONS)
         real(kind=RP), intent(out)      :: fR(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: cL,uL, vL, wL, pL, invRhoL, invSqrtRhoL, lambdaMinusL, lambdaPlusL
         real(kind=RP)  :: cR,uR, vR, wR, pR, invRhoR, invSqrtRhoR, lambdaMinusR, lambdaPlusR
         real(kind=RP)  :: rhoStarL, rhoStarR, uStar, pStar, rhoStar, vStar, wStar, cuStar, halfRhouStar
         real(kind=RP)  :: QLRot(NCONS), QRRot(NCONS) 
         real(kind=RP)  :: lambda_mu = 0.0_RP
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         invRhoL     = 1.0_RP / rhoL
         invSqrtRhoL = sqrt(invRhoL)
         cL = QLeft(IMC)
         uL = invSqrtRhoL * (QLeft(IMSQRHOU) * nHat(1) + QLeft(IMSQRHOV) * nHat(2) + QLeft(IMSQRHOW) * nHat(3))
         vL = invSqrtRhoL * (QLeft(IMSQRHOU) * t1(1)   + QLeft(IMSQRHOV) * t1(2)   + QLeft(IMSQRHOW) * t1(3))
         wL = invSqrtRhoL * (QLeft(IMSQRHOU) * t2(1)   + QLeft(IMSQRHOV) * t2(2)   + QLeft(IMSQRHOW) * t2(3))
         pL = QLeft(IMP)

         invRhoR     = 1.0_RP / rhoR
         invSqrtRhoR = sqrt(invRhoR)
         cR = QRight(IMC)
         uR = invSqrtRhoR * (QRight(IMSQRHOU) * nHat(1) + QRight(IMSQRHOV) * nHat(2) + QRight(IMSQRHOW) * nHat(3))
         vR = invSqrtRhoR * (QRight(IMSQRHOU) * t1(1)   + QRight(IMSQRHOV) * t1(2)   + QRight(IMSQRHOW) * t1(3))
         wR = invSqrtRhoR * (QRight(IMSQRHOU) * t2(1)   + QRight(IMSQRHOV) * t2(2)   + QRight(IMSQRHOW) * t2(3))
         pR = QRight(IMP)
!
!        Compute the Star Region
!        -----------------------
         lambdaMinusR = 0.5_RP * (uR - sqrt(uR*uR + 4.0_RP*dimensionless % invMa2/rhoR))
         lambdaPlusR  = 0.5_RP * (uR + sqrt(uR*uR + 4.0_RP*dimensionless % invMa2/rhoR))

         lambdaMinusL = 0.5_RP * (uL - sqrt(uL*uL + 4.0_RP*dimensionless % invMa2/rhoL))
         lambdaPlusL  = 0.5_RP * (uL + sqrt(uL*uL + 4.0_RP*dimensionless % invMa2/rhoL))

         uStar = (pR-pL+rhoR*uR*lambdaMinusR-rhoL*uL*lambdaPlusL)/(rhoR*lambdaMinusR - rhoL*lambdaPlusL)
         pStar = pR + rhoR*lambdaMinusR*(uR-uStar)
         rhoStarL = (rhoL*lambdaPlusL)/(uStar-lambdaMinusL)
         rhoStarR = (rhoR*lambdaMinusR)/(uStar - lambdaPlusR)

         if ( uStar .ge. 0.0_RP ) then
            rhoStar = rhoStarL
            vStar   = vL
            wStar   = wL

         else
            rhoStar = rhoStarR
            vStar   = vR
            wStar   = wR

         end if

         cuStar = 0.5_RP*(cL*uL + cR*uR)
         halfRhouStar = 0.5_RP*rhoStar*uStar
!
!      - Add first the common (conservative) part
         fL = [cuStar+lambda_mu*(muL-muR), rhoStar*uStar*uStar + pStar, rhoStar*uStar*vStar, rhoStar*uStar*wStar, dimensionless % invMa2 * uStar]
         fR = fL
!
!      - Add the non--conservative part
         fL = fL + [0.0_RP, cL*0.5_RP*(muR-muL)-halfRhouStar*uL,-halfRhouStar*vL, -halfRhouStar*wL, -dimensionless % invMa2*uL]
         fR = fR + [0.0_RP, cR*0.5_RP*(muL-muR)-halfRhouStar*uR,-halfRhouStar*vR, -halfRhouStar*wR, -dimensionless % invMa2*uR]
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         fL(2:4) = nHat*fL(2) + t1*fL(3) + t2*fL(4)
         fR(2:4) = nHat*fR(2) + t1*fR(3) + t2*fR(4)

      end subroutine ExactRiemannSolver

end module RiemannSolvers_MU
