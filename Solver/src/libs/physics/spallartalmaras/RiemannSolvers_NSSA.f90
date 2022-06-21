!
#include "Includes.h"
module RiemannSolvers_NSSA
   use SMConstants
   use Physics_NSSA
   use PhysicsStorage_NSSA
   use VariableConversion_NSSA
   use FluidData_NSSA

   private 
   public RiemannSolver, SetRiemannSolver, RiemannSolver_dFdQ

   abstract interface
      subroutine RiemannSolverFCN(QLeft, QRight, nHat, t1, t2, flux)
         use SMConstants
         use PhysicsStorage_NSSA
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
      end subroutine RiemannSolverFCN
      subroutine AveragedStatesFCN(QLeft, QRight, pL, pR, invRhoL, invRhoR, flux) 
         use SMConstants
         use PhysicsStorage_NSSA
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
      end subroutine AveragedStatesFCN
      subroutine RiemannSolver_dFdQFCN(ql,qr,nHat,dfdq_num,side)
         use SMConstants
         use PhysicsStorage_NSSA
         implicit none
         real(kind=RP), intent (in)  :: ql(NCONS)                 !<  Current solution on the left
         real(kind=RP), intent (in)  :: qr(NCONS)                 !<  Current solution on the right
         real(kind=RP), intent (in)  :: nHat(NDIM)                !<  Normal vector
         real(kind=RP), intent(out)  :: dfdq_num(NCONS,NCONS)     !>  Numerical flux Jacobian 
         integer      , intent (in)  :: side   
      end subroutine RiemannSolver_dFdQFCN
   end interface

   procedure(RiemannSolverFCN)     , pointer  :: RiemannSolver      => NULL()
   procedure(AveragedStatesFCN)    , pointer  :: AveragedStates     => NULL()
   procedure(RiemannSolver_dFdQFCN), pointer  :: RiemannSolver_dFdQ => NULL()
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
         
         RiemannSolver_dFdQ => BaseClass_RiemannSolver_dFdQ
         
         select case ( which )
#if defined(SPALARTALMARAS)

         case ( RIEMANN_ROE )
            RiemannSolver => RoeRiemannSolver

         case ( RIEMANN_LXF )
            RiemannSolver => LxFRiemannSolver

         case ( RIEMANN_RUSANOV)
            RiemannSolver => RusanovRiemannSolver

         CASE ( RIEMANN_STDROE)
            RiemannSolver => StdRoeRiemannSolver

         case ( RIEMANN_CENTRAL )
            RiemannSolver => CentralRiemannSolver

         case ( RIEMANN_ROEPIKE )
            RiemannSolver => RoePikeRiemannSolver

         case ( RIEMANN_LOWDISSROE )
            RiemannSolver => LowDissipationRoeRiemannSolver

         case ( RIEMANN_MATRIXDISS )
            RiemannSolver => MatrixDissipationRiemannSolver

         case default
            print*, "Undefined choice of Riemann Solver."
            print*, "Options available are:"
            print*, "   * Central"
            print*, "   * Roe"
            print*, "   * Standard Roe"
            print*, "   * Roe-Pike"
            print*, "   * Low dissipation Roe"
            print*, "   * Lax-Friedrichs"
            print*, "   * Rusanov"
            print*, "   * Matrix dissipation"
            errorMessage(STD_OUT)
            STOP

#endif
         end select

!        Set up an averaging function
!        ----------------------------
         select case ( splitType )
         case (STANDARD_SPLIT)
            AveragedStates => StandardAverage
            whichAverage = STANDARD_SPLIT

         case (DUCROS_SPLIT)
            AveragedStates => DucrosAverage
            whichAverage = DUCROS_SPLIT

         case (MORINISHI_SPLIT)
            AveragedStates => MorinishiAverage
            whichAverage = MORINISHI_SPLIT
   
         case (PIROZZOLI_SPLIT)
            AveragedStates => PirozzoliAverage
            whichAverage = PIROZZOLI_SPLIT

         case (KENNEDYGRUBER_SPLIT)
            AveragedStates => KennedyGruberAverage
            whichAverage = KENNEDYGRUBER_SPLIT

         case (ENTROPYCONS_SPLIT)
            AveragedStates => EntropyConservingAverage
            whichAverage = ENTROPYCONS_SPLIT

         case (CHANDRASEKAR_SPLIT)
            AveragedStates => ChandrasekarAverage
            whichAverage = CHANDRASEKAR_SPLIT

         case default
            print*, "Split form not recognized"
            errorMessage(STD_OUT)
            stop
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
      subroutine BaseClass_RiemannSolver_dFdQ(ql,qr,nHat,dfdq_num,side)
         implicit none
         !--------------------------------------------
         real(kind=RP), intent (in)  :: ql(NCONS)                 !<  Current solution on the left
         real(kind=RP), intent (in)  :: qr(NCONS)                 !<  Current solution on the right
         real(kind=RP), intent (in)  :: nHat(NDIM)                !<  Normal vector
         real(kind=RP), intent(out)  :: dfdq_num(NCONS,NCONS)     !>  Numerical flux Jacobian 
         integer      , intent (in)  :: side                      !<  Either LEFT or RIGHT
         !--------------------------------------------
         
         ERROR stop 'Requested Riemann solver not implemented for implicit time-integration'
      end subroutine BaseClass_RiemannSolver_dFdQ

      subroutine CentralRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         implicit none 
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: rhoL, rhouL, rhovL, rhowL, rhoeL, pL, rhoV2L, thetaL
         real(kind=RP) :: rhoR, rhouR, rhovR, rhowR, rhoeR, pR, rhoV2R, thetaR
         real(kind=RP) :: invRhoL, invRhoR
         real(kind=RP) :: QLRot(NCONS), QRRot(NCONS), oflux(NCONS)

         associate(gm1 => thermodynamics % GammaMinus1)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(1)                  ; rhoR = QRight(1)
         invRhoL = 1.0_RP / rhoL          ; invRhoR = 1.0_RP / rhoR
         
         thetaL = QLeft(IRHOTHETA)/rhoL       ; thetaR = QRight(IRHOTHETA)/rhoR 

         rhouL = QLeft(2) * nHat(1) + QLeft(3) * nHat(2) + QLeft(4) * nHat(3)
         rhovL = QLeft(2) * t1(1)   + QLeft(3) * t1(2)   + QLeft(4) * t1(3)
         rhowL = QLeft(2) * t2(1)   + QLeft(3) * t2(2)   + QLeft(4) * t2(3)

         rhouR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)
         rhovR = QRight(2) * t1(1)   + QRight(3) * t1(2)   + QRight(4) * t1(3)
         rhowR = QRight(2) * t2(1)   + QRight(3) * t2(2)   + QRight(4) * t2(3)

         rhoeL = QLeft(5) ; rhoeR = QRight(5)                

         rhoV2L = (POW2(rhouL) + POW2(rhovL) + POW2(rhowL)) * invRhoL
         rhoV2R = (POW2(rhouR) + POW2(rhovR) + POW2(rhowR)) * invRhoR

         pL = gm1 * (rhoeL - 0.5_RP * rhoV2L)
         pR = gm1 * (rhoeR - 0.5_RP * rhoV2R)
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         QLRot = (/ rhoL, rhouL, rhovL, rhowL, rhoeL, 0.0_RP /)
         QRRot = (/ rhoR, rhouR, rhovR, rhowR, rhoeR, 0.0_RP /)

         call AveragedStates(QLRot, QRRot, pL, pR, invRhoL, invRhoR, flux)

         flux(6) = (thetaR + thetaL)*0.5_RP
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate

      end subroutine CentralRiemannSolver

      subroutine StdRoeRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         implicit none 
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i
         real(kind=RP)  :: QLRot(NCONS), QRRot(NCONS), VL(NPRIM), VR(NPRIM), aL, aR
         real(kind=RP)  :: dQ(5), lambda(5), K(5,5), V2abs, alpha(5), dLambda
         real(kind=RP)  :: rho, u, v, w, V2, H, a
         real(kind=RP)  :: stab(5)     ! Careful with this variable
         real(kind=RP)  :: thetaL, thetaR
         associate(gm1 => thermodynamics % gammaMinus1)
!
!        ********************
!        Perform the rotation
!        ********************
!
         QLRot(1) = QLeft(1)  ; QRRot(1) = QRight(1)

         QLRot(2) = QLeft (2) * nHat(1) + QLeft (3) * nHat(2) + QLeft (4) * nHat(3)
         QRRot(2) = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         QLRot(3) = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         QRRot(3) = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         QLRot(4) = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         QRRot(4) = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         QLRot(5) = QLeft(5) ; QRRot(5) = QRight(5)    

         thetaL = QLeft(IRHOTHETA)/QLeft(1)  ;  thetaR = QRight(IRHOTHETA)/QRight(1)            
!
!        ***************************
!        Compute primitive variables
!        ***************************
!
         call getPrimitiveVariables(QLRot, VL)
         call getPrimitiveVariables(QRRot, VR)

         aL = sqrt(VL(IPA2))  ; aR = sqrt(VR(IPA2))
!
!        *********************
!        Compute Roe variables: [rho, u, v, w, H, a]
!        *********************
!
         call getRoeVariables(QLRot, QRRot, VL, VR, rho, u, v, w, V2, H, a)
!
!        Eigenvalues
!        -----------
         lambda(1)   = u-a
         lambda(2:4) = u
         lambda(5)   = u+a
!
!        Eigenvectors
!        ------------
         K(:,1) = (/ 1.0_RP, u-a, v, w, H-u*a /)
         K(:,2) = (/ 1.0_RP, u, v, w, 0.5_RP*V2 /)
         K(:,3) = (/ 0.0_RP, 0.0_RP, 1.0_RP, 0.0_RP, v /)
         K(:,4) = (/ 0.0_RP, 0.0_RP, 0.0_RP, 1.0_RP, w /)
         K(:,5) = (/ 1.0_RP, u+a, v, w, H+u*a /)
!
!        Projections
!        -----------
         dQ = QRRot(1:5) - QLRot(1:5)

         alpha(3) = dQ(3) - v * dQ(1)  ; alpha(4) = dQ(4) - w * dQ(1)

         dQ(5) = dQ(5) - alpha(3) * v - alpha(4) * w

         alpha(2) = gm1 * ( dQ(1)*(H - u*u) + u * dQ(2) - dQ(5) ) / (POW2(a))
         alpha(1) = 0.5_RP * (dQ(1)*lambda(5) - dQ(2) - a*alpha(2)) / a
         alpha(5) = dQ(1) - alpha(1) - alpha(2)
!
!        **********************
!        Perform an entropy fix. Here we use Van Leer's modification of Harten's entropy fix, derived
!        in: A. Harten, "High resolution schemes for hyperbolic conservation laws". To recover the
!        Harten entropy fix, set dLambda to 0.5
!        **********************
!
!        Wave #1
!        -------
         dLambda = max((VR(IPU)-aR) - (VL(IPU)-aL), 0.0_RP)
         if ( abs(lambda(1)) .ge. 2.0_RP * dLambda ) then
            lambda(1) = abs(lambda(1))
         
         else
            lambda(1) = POW2(lambda(1)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        Wave #5
!        -------
         dLambda = max((VR(IPU)+aR) - (VL(IPU)+aL), 0.0_RP)
         if ( abs(lambda(5)) .ge. 2.0_RP * dLambda ) then
            lambda(5) = abs(lambda(5))
         
         else
            lambda(5) = POW2(lambda(5)) / (4.0_RP * dLambda) + dLambda

         end if

!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         call AveragedStates(QLRot, QRRot, VL(IPP), VR(IPP), VL(IPIRHO), VR(IPIRHO), flux)
!
!        Compute the Roe stabilization
!        -----------------------------
         select case (whichAverage)
         case(PIROZZOLI_SPLIT, KENNEDYGRUBER_SPLIT)
!
!           ***************************************************************************
!           Eigenvalue matrix is corrected for PI and KG variants, see Winters et. al. 
!           "A comparative study on polynomial dealiasing and split form discontinuous 
!           Galerkin schemes for under-resolved turbulence computations"
!           ***************************************************************************
!
            lambda(1) = lambda(5)
         end select
         
         stab = 0.0_RP
         do i = 1, 5
            stab = stab + 0.5_RP * alpha(i) * abs(lambda(i)) * K(:,i)
         end do
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux(1:5) = flux(1:5) - lambdaStab * stab

         flux(6) = flux(1)*(thetaL+thetaR)*0.5_RP - abs(flux(1))*(thetaR-thetaL)*0.5_RP
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate
         
      end subroutine StdRoeRiemannSolver

      subroutine MatrixDissipationRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         use Utilities, only: logarithmicMean
         implicit none 
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k
         real(kind=RP)  :: QLRot(NCONS), QRRot(NCONS), betaL, betaR
         real(kind=RP)  :: EVL(5), EVR(5)
         real(kind=RP)  :: dQ(5)
         real(kind=RP)  :: a_bar, h_bar
         real(kind=RP)  :: R1(5,5), T(5), Lambda(5)
         real(kind=RP)  :: stab(5)
         real(kind=RP)  :: rhoLogMean, betaLogMean, pMean, uMean, vMean, wMean, V2abs
         real(kind=RP)  :: uL, vL, wL, uR, vR, wR, vtotL, vtotR, pL, pR
         real(kind=RP)  :: invRhoL, invRhoR
         real(kind=RP)  :: thetaL, thetaR

         associate(gm1 => thermodynamics % gammaMinus1, gamma => thermodynamics % gamma, &
                   invGamma => thermodynamics % invGamma, cp => thermodynamics % GammaDivGammaMinus1, &
                   gammaMinus1Div2g => thermodynamics % gammaMinus1Div2g)
!
!        ********************
!        Perform the rotation
!        ********************
!
         QLRot(1) = QLeft(1)  ; QRRot(1) = QRight(1)

         QLRot(2) = QLeft (2) * nHat(1) + QLeft (3) * nHat(2) + QLeft (4) * nHat(3)
         QRRot(2) = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         QLRot(3) = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         QRRot(3) = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         QLRot(4) = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         QRRot(4) = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         QLRot(5) = QLeft(5) ; QRRot(5) = QRight(5)     

         thetaL = QLeft(IRHOTHETA)/QLeft(1)  ;  thetaR = QRight(IRHOTHETA)/QRight(1)            

!
!        *************************
!        Compute Entropy variables
!        *************************
!
         call NSGradientVariables_ENTROPY(NCONS,NGRAD,QLRot,EVL)
         call NSGradientVariables_ENTROPY(NCONS,NGRAD,QRRot,EVR)

         invRhoL = 1.0_RP / QLRot(IRHO) ; invrhoR = 1.0_RP / QRRot(IRHO)
         uL      = QLRot(IRHOU)*invRhoL ; uR      = QRRot(IRHOU)*invRhoR
         vL      = QLRot(IRHOV)*invRhoL ; vR      = QRRot(IRHOV)*invRhoR
         wL      = QLRot(IRHOW)*invRhoL ; wR      = QRRot(IRHOW)*invRhoR
         vtotL   = uL*uL + vL*vL + wL*wL; vtotR   = uR*uR + vR*vR + wR*wR

         pL      = gm1*(QLRot(IRHOE)-0.5_RP*QLRot(IRHO)*vtotL)
         pR      = gm1*(QRRot(IRHOE)-0.5_RP*QRRot(IRHO)*vtotR)

         betaL   = -0.5_RP*EVL(IRHOE) ; betaR = -0.5_RP*EVR(IRHOE)

         call logarithmicMean(betaL      , betaR      , betaLogMean)
         call logarithmicMean(QLRot(IRHO), QRRot(IRHO), rhoLogMean)

         pMean = 0.5_RP*(QLRot(IRHO)+QRRot(IRHO))/(betaL + betaR)
         a_bar = sqrt(gamma * pMean / rhoLogMean)
         uMean = AVERAGE(uL, uR)
         vMean = AVERAGE(vL, vR)
         wMean = AVERAGE(wL, wR)
         v2Abs =   2.0_RP * ( POW2(uMean) + POW2(vMean) + POW2(wMean) )      &
                 - 0.5_RP * ( vtotL + vtotR )

         h_bar = 0.5_RP * ( cp / betaLogMean + v2Abs ) 
!
!        ***********************
!        Compute the eigenvalues
!        ***********************
!
         lambda(1)   = abs(uMean - a_bar)
         lambda(2:4) = abs(uMean        )
         lambda(5)   = abs(uMean + a_bar)
!
!        Eigenvectors
!        ------------
         R1(:, 1) = (/ 1.0_RP, uMean-a_bar, vMean , wMean , h_bar-uMean*a_bar/) 
         R1(:, 2) = (/ 1.0_RP, uMean      , vMean , wMean , 0.5_RP*V2abs     /) 
         R1(:, 3) = (/ 0.0_RP, 0.0_RP     , 1.0_RP, 0.0_RP, vMean            /) 
         R1(:, 4) = (/ 0.0_RP, 0.0_RP     , 0.0_RP, 1.0_RP, wMean            /) 
         R1(:, 5) = (/ 1.0_RP, uMean+a_bar, vMean , wMean , h_bar+uMean*a_bar/) 
!
!        Intensities
!        -----------
         T(1) = 0.5_RP * rhoLogMean * invGamma
         T(2) = 2.0_RP * gammaMinus1Div2g * rhoLogMean
         T(3) = pMean
         T(4) = pMean
         T(5) = T(1)
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         call AveragedStates(QLRot, QRRot, pL, pR, invRhoL, invRhoR, flux)
!
!        Compute the Roe stabilization
!        -----------------------------
         select case (whichAverage)
         case(PIROZZOLI_SPLIT, KENNEDYGRUBER_SPLIT)
!
!           ***************************************************************************
!           Eigenvalue matrix is corrected for PI and KG variants, see Winters et. al. 
!           "A comparative study on polynomial dealiasing and split form discontinuous 
!           Galerkin schemes for under-resolved turbulence computations"
!           ***************************************************************************
!
            lambda(1) = lambda(5)
         end select
         
         stab = 0.0_RP
         do i = 1, 5
            do j = 1, 5 ; do k = 1, 5
               stab(i) = stab(i) + 0.5_RP * R1(i,j) * lambda(j) * T(j) * R1(k,j) * (EVR(k) - EVL(k))
            end do      ; end do
         end do
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux(1:5) = flux(1:5) - lambdaStab * stab
!
         flux(6) = flux(1)*(thetaL+thetaR)*0.5_RP - abs(flux(1))*(thetaR-thetaL)*0.5_RP

!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate
         
      end subroutine MatrixDissipationRiemannSolver

      subroutine RoePikeRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
         implicit none 
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i
         real(kind=RP)  :: QLRot(NCONS), QRRot(NCONS), VL(NPRIM), VR(NPRIM), aL, aR
         real(kind=RP)  :: dQ(5), lambda(5), K(5,5), V2abs, alpha(5), dLambda
         real(kind=RP)  :: rho, u, v, w, V2, H, a
         real(kind=RP)  :: stab(5)    
         real(kind=RP)  :: thetaL, thetaR

         associate(gm1 => thermodynamics % gammaMinus1)
!
!        ********************
!        Perform the rotation
!        ********************
!
         QLRot(1) = QLeft(1)  ; QRRot(1) = QRight(1)

         QLRot(2) = QLeft (2) * nHat(1) + QLeft (3) * nHat(2) + QLeft (4) * nHat(3)
         QRRot(2) = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         QLRot(3) = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         QRRot(3) = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         QLRot(4) = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         QRRot(4) = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         QLRot(5) = QLeft(5) ; QRRot(5) = QRight(5)                

         thetaL = QLeft(IRHOTHETA)/QLeft(1)  ;  thetaR = QRight(IRHOTHETA)/QRight(1)            
!
!        ***************************
!        Compute primitive variables
!        ***************************
!
         call getPrimitiveVariables(QLRot, VL)
         call getPrimitiveVariables(QRRot, VR)

         aL = sqrt(VL(IPA2))  ; aR = sqrt(VR(IPA2))
!
!        *********************
!        Compute Roe variables: [rho, u, v, w, H, a]
!        *********************
!
         call getRoeVariables(QLRot, QRRot, VL, VR, rho, u, v, w, V2, H, a)
!
!        Eigenvalues
!        -----------
         lambda(1)   = u-a
         lambda(2:4) = u
         lambda(5)   = u+a
!
!        Eigenvectors
!        ------------
         K(:,1) = (/ 1.0_RP, u-a, v, w, H-u*a /)
         K(:,2) = (/ 1.0_RP, u, v, w, 0.5_RP*V2 /)
         K(:,3) = (/ 0.0_RP, 0.0_RP, 1.0_RP, 0.0_RP, v /)
         K(:,4) = (/ 0.0_RP, 0.0_RP, 0.0_RP, 1.0_RP, w /)
         K(:,5) = (/ 1.0_RP, u+a, v, w, H+u*a /)
!
!        Projections
!        -----------
         alpha(1) = ((VR(IPP)-VL(IPP)) - rho * a * (VR(IPU)-VL(IPU)))/(2.0_RP * a * a)
         alpha(2) = (QRight(IRHO)-QLeft(IRHO)) - (VR(IPP)-VL(IPP))/(a*a)
         alpha(3) = rho * (VR(IPV)-VL(IPV))
         alpha(4) = rho * (VR(IPW)-VL(IPW))
         alpha(5) = ((VR(IPP)-VL(IPP)) + rho * a * (VR(IPU)-VL(IPU)))/(2.0_RP * a * a)
!
!        **********************
!        Perform an entropy fix. Here we use Van Leer's modification of Harten's entropy fix, derived
!        in: A. Harten, "High resolution schemes for hyperbolic conservation laws". To recover the
!        Harten entropy fix, set dLambda to 0.5
!        **********************
!
!        Wave #1
!        -------
         dLambda = max((VR(IPU)-aR) - (VL(IPU)-aL), 0.0_RP)
         if ( abs(lambda(1)) .ge. 2.0_RP * dLambda ) then
            lambda(1) = abs(lambda(1))
         
         else
            lambda(1) = POW2(lambda(1)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        Wave #5
!        -------
         dLambda = max((VR(IPU)+aR) - (VL(IPU)+aL), 0.0_RP)
         if ( abs(lambda(5)) .ge. 2.0_RP * dLambda ) then
            lambda(5) = abs(lambda(5))
         
         else
            lambda(5) = POW2(lambda(5)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         call AveragedStates(QLRot, QRRot, VL(IPP), VR(IPP), VL(IPIRHO), VR(IPIRHO), flux)
!
!        Compute the Roe stabilization
!        -----------------------------
         select case (whichAverage)
         case(PIROZZOLI_SPLIT, KENNEDYGRUBER_SPLIT)
!
!           ***************************************************************************
!           Eigenvalue matrix is corrected for PI and KG variants, see Winters et. al. 
!           "A comparative study on polynomial dealiasing and split form discontinuous 
!           Galerkin schemes for under-resolved turbulence computations"
!           ***************************************************************************
!
            lambda(1) = lambda(5)
         end select
         
         stab = 0.0_RP
         do i = 1, 5
            stab = stab + 0.5_RP * alpha(i) * abs(lambda(i)) * K(:,i)
         end do
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux(1:5) = flux(1:5) - lambdaStab * stab
!
         flux(6) = flux(1)*(thetaL+thetaR)*0.5_RP - abs(flux(1))*(thetaR-thetaL)*0.5_RP
!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate
      end subroutine RoePikeRiemannSolver
 
      subroutine LowDissipationRoeRiemannSolver(QLeft, QRight, nHat, t1, t2, flux)
!
!        ***********************************************************************
!           This version, presented by Oßwald et. al. [*], is a modification of 
!           Roe-Pike solver, decreasing the velocity jumps intensity for low
!           Mach numbers. Normal velocities are scaled such that
!                         du <- z*du
!           where z tends to zero as the Mach number tends to zero. Precisely:
!                         z = min(1, max(ML, MR))
!
!           These normal velocities are just scaled to compute Roe dissipation's
!           intensities (alpha coefficients), not the fluxes (as in the original
!           Low dissipation method by Thornber et. al.)
!
!           * Oßwald et. al. L2Roe: a low dissipation version of Roe’s a
!              pproximate Riemann solver for low Mach numbers
!        ***********************************************************************
!
         implicit none 
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i
         real(kind=RP)  :: rhoL, rhouL, rhovL, rhowL, rhoeL, pL, rhoHL, rhoV2L, ML
         real(kind=RP)  :: rhoR, rhouR, rhovR, rhowR, rhoeR, pR, rhoHR, rhoV2R, MR
         real(kind=RP)  :: uL, vL, wL, uR, vR, wR, aL, aR, dLambda, z, du, dv, dw
         real(kind=RP)  :: dp
         real(kind=RP)  :: QLRot(NCONS), QRRot(NCONS)
         real(kind=RP)  :: sqrtRhoL, sqrtRhoR, invSumSqrtRhoLR
         real(kind=RP)  :: invSqrtRhoL, invSqrtRhoR, invRhoL, invRhoR
         real(kind=RP)  :: rho, u, v, w, H, a, dQ(5), lambda(5), K(5,5), V2abs, alpha(5)
         real(kind=RP)  :: stab(5)
         real(kind=RP)  :: thetaL, thetaR

         associate(gamma => thermodynamics % gamma, gm1 => thermodynamics % gammaMinus1)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(1)                  ; rhoR = QRight(1)
         invRhoL = 1.0_RP/ rhoL           ; invRhoR = 1.0_RP / rhoR
         sqrtRhoL = sqrt(rhoL)            ; sqrtRhoR = sqrt(rhoR)
         invSqrtRhoL = 1.0_RP / sqrtRhoL  ; invSqrtRhoR = 1.0_RP / sqrtRhoR
         invSumSqrtRhoLR = 1.0_RP / (sqrtRhoL + sqrtRhoR)

         rhouL = QLeft (2) * nHat(1) + QLeft (3) * nHat(2) + QLeft (4) * nHat(3)
         rhouR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         rhovL = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         rhovR = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         rhowL = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         rhowR = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         rhoeL = QLeft(5) ; rhoeR = QRight(5)                

         thetaL = QLeft(IRHOTHETA)/QLeft(1)  ;  thetaR = QRight(IRHOTHETA)/QRight(1)            


         uL = rhouL * invRhoL    ; uR = rhouR * invRhoR
         vL = rhovL * invRhoL    ; vR = rhovR * invRhoR
         wL = rhowL * invRhoL    ; wR = rhowR * invRhoR

         rhoV2L = (POW2(uL) + POW2(vL) + POW2(wL)) * rhoL
         rhoV2R = (POW2(uR) + POW2(vR) + POW2(wR)) * rhoR
!
!        Compute the enthalpy: here defined as rhoH = gogm1 p + 0.5 rho V^2
!        --------------------
         rhoHL = gamma*rhoeL - 0.5_RP*gm1*rhoV2L
         rhoHR = gamma*rhoeR - 0.5_RP*gm1*rhoV2R

         pL = gm1 * (rhoeL - 0.5_RP * rhoV2L)
         pR = gm1 * (rhoeR - 0.5_RP * rhoV2R)
         
         aL = sqrt(gamma * pL * invRhoL)
         aR = sqrt(gamma * pR * invRhoR)
!
!        Compute Roe - Pike variables
!        ----------------------------
         rho = sqrtRhoL * sqrtRhoR
         u = (invSqrtRhoL * rhouL + invSqrtRhoR * rhouR) * invSumSqrtRhoLR
         v = (invSqrtRhoL * rhovL + invSqrtRhoR * rhovR) * invSumSqrtRhoLR
         w = (invSqrtRhoL * rhowL + invSqrtRhoR * rhowR) * invSumSqrtRhoLR
         H = (invSqrtRhoL * rhoHL + invSqrtRhoR * rhoHR) * invSumSqrtRhoLR
         V2abs = POW2(u) + POW2(v) + POW2(w)
         a = sqrt(gm1*(H - 0.5_RP*V2abs))
!
!        Eigenvalues
!        -----------
         lambda(1)   = u-a
         lambda(2:4) = u
         lambda(5)   = u+a
!
!        Eigenvectors
!        ------------
         K(:,1) = (/ 1.0_RP, u-a, v, w, H-u*a /)
         K(:,2) = (/ 1.0_RP, u, v, w, 0.5_RP*V2abs /)
         K(:,3) = (/ 0.0_RP, 0.0_RP, 1.0_RP, 0.0_RP, v /)
         K(:,4) = (/ 0.0_RP, 0.0_RP, 0.0_RP, 1.0_RP, w /)
         K(:,5) = (/ 1.0_RP, u+a, v, w, H+u*a /)
!
!        Projections
!        -----------
!
!        ----------------------------------------------------------------------------
!        Low dissipation Roe-Pike Riemann solver: Reduce the dissipation associated
!        to the jump in normal velocity. See Obwald et. al. L2Roe: a low dissipation 
!        version of Roe’s approximate Riemann solver for low Mach numbers
!        ----------------------------------------------------------------------------

         ML = abs(uL) / aL    ; MR = abs(uR) / aR
         z  = min(1.0_RP, max(ML,MR))

         du = z * (uR - uL)
         dv = z * (vR - vL)
         dw = z * (wR - wL)
         dp = pR - pL

         alpha(1) = (dp - rho * a * du)/(2.0_RP * a * a)
         alpha(2) = (rhoR-rhoL) - dp/(a*a)
         alpha(3) = rho * dv
         alpha(4) = rho * dw
         alpha(5) = (dp + rho * a * du)/(2.0_RP * a * a)
!
!        **********************
!        Perform an entropy fix. Here we use Van Leer's modification of Harten's entropy fix, derived
!        in: A. Harten, "High resolution schemes for hyperbolic conservation laws". To recover the
!        Harten entropy fix, set dLambda to 0.5
!        **********************
!
!        Wave #1
!        -------
         dLambda = max((uR-aR) - (uL-aL), 0.0_RP)
         if ( abs(lambda(1)) .ge. 2.0_RP * dLambda ) then
            lambda(1) = abs(lambda(1))
         
         else
            lambda(1) = POW2(lambda(1)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        Wave #5
!        -------
         dLambda = max((uR+aR) - (uL+aL), 0.0_RP)
         if ( abs(lambda(5)) .ge. 2.0_RP * dLambda ) then
            lambda(5) = abs(lambda(5))
         
         else
            lambda(5) = POW2(lambda(5)) / (4.0_RP * dLambda) + dLambda

         end if
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         QLRot = (/ rhoL, rhouL, rhovL, rhowL, rhoeL, 0.0_RP /)
         QRRot = (/ rhoR, rhouR, rhovR, rhowR, rhoeR, 0.0_RP /)
         call AveragedStates(QLRot, QRRot, pL, pR, invRhoL, invRhoR, flux)
!
!        Compute the Roe stabilization
!        -----------------------------
         select case (whichAverage)
         case(PIROZZOLI_SPLIT, KENNEDYGRUBER_SPLIT)
!
!           ***************************************************************************
!           Eigenvalue matrix is corrected for PI and KG variants, see Winters et. al. 
!           "A comparative study on polynomial dealiasing and split form discontinuous 
!           Galerkin schemes for under-resolved turbulence computations"
!           ***************************************************************************
!
            lambda(1) = lambda(5)
         end select
         
         stab = 0.0_RP
         do i = 1, 5
            stab = stab + 0.5_RP * alpha(i) * abs(lambda(i)) * K(:,i)
         end do
!
!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux(1:5) = flux(1:5) - lambdaStab * stab
!
         flux(6) = flux(1)*(thetaL+thetaR)*0.5_RP - abs(flux(1))*(thetaR-thetaL)*0.5_RP!
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!

         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate

      end subroutine LowDissipationRoeRiemannSolver

!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LxFRiemannSolver( QLeft, QRight, nHat, t1, t2, flux ) 
         implicit none 
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         real(kind=RP)  :: rhoL, rhouL, rhovL, rhowL, rhoeL, pL , aL, rhoV2L, thetaL
         real(kind=RP)  :: rhoR, rhouR, rhovR, rhowR, rhoeR, pR, aR, rhoV2R, thetaR
         real(kind=RP)  :: QLRot(NCONS), QRRot(NCONS)
         real(kind=RP)  :: invRhoL, invRhoR
         real(kind=RP)  :: lambda, stab(NCONS)

         associate(gamma => thermodynamics % gamma, gm1 => thermodynamics % gammaMinus1)
!
!        Rotate the variables to the face local frame using normal and tangent vectors
!        -----------------------------------------------------------------------------
         rhoL = QLeft(1)                  ; rhoR = QRight(1)
         invRhoL = 1.0_RP/ rhoL           ; invRhoR = 1.0_RP / rhoR

         rhouL = QLeft(2)  * nHat(1) + QLeft(3)  * nHat(2) + QLeft(4)  * nHat(3)
         rhouR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

         rhovL = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
         rhovR = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

         rhowL = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
         rhowR = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

         rhoV2L = (POW2(rhouL) + POW2(rhovL) + POW2(rhowL)) * invRhoL
         rhoV2R = (POW2(rhouR) + POW2(rhovR) + POW2(rhowR)) * invRhoR

         rhoeL = QLeft(5) ; rhoeR = QRight(5)                

         thetaL = QLeft(6)  * invRhoL
         thetaR = QRight(6) * invRhoR

         pL = gm1 * (rhoeL - 0.5_RP * rhoV2L)
         pR = gm1 * (rhoeR - 0.5_RP * rhoV2R)

         aL = sqrt(gamma * pL * invRhoL)
         aR = sqrt(gamma * pR * invRhoR)
!
!        Eigenvalues: lambda = max(|uL|,|uR|) + max(aL,aR)
!        -----------
!         lambda = max(abs(rhouL*invRhoL),abs(rhouR*invRhoR)) + max(aL, aR)   ! This is a more dissipative version (not consistent with the Jacobian below)

         lambda = max(abs(rhouL*invRhoL) + aL,abs(rhouR*invRhoR) + aR)
!
!        ****************
!        Compute the flux
!        ****************
!
!        Perform the average using the averaging function
!        ------------------------------------------------
         QLRot = (/ rhoL, rhouL, rhovL, rhowL, rhoeL, 0.0_RP /)
         QRRot = (/ rhoR, rhouR, rhovR, rhowR, rhoeR, 0.0_RP /)

         call AveragedStates(QLRot, QRRot, pL, pR, invRhoL, invRhoR, flux)
!
!        Compute the Lax-Friedrichs stabilization
!        ----------------------------------------
         stab(1:5) = 0.5_RP * lambda * (QRRot(1:5) - QLRot(1:5))

!        Compute the flux: apply the lambda stabilization here.
!        ----------------
         flux = flux - lambdaStab * stab
         flux(6) = flux(1)*(thetaL+thetaR)*0.5_RP - abs(flux(1))*(thetaR-thetaL)*0.5_RP

!        
!        ************************************************
!        Return momentum equations to the cartesian frame
!        ************************************************
!
         flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)

         end associate
         
      END SUBROUTINE LxFRiemannSolver

!
      SUBROUTINE RoeRiemannSolver( QLeft, QRight, nHat, t1, t2, flux )
!
!        **************************************************************
!           This Roe Riemann solver implementation does not uses the
!           averaging function, nor the lambda stabilization
!        **************************************************************
!
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         REAL(KIND=RP) :: rho , rhou , rhov , rhow  , rhoe
         REAL(KIND=RP) :: rhon, rhoun, rhovn, rhown , rhoen
         REAL(KIND=RP) :: ul  , vl   , wl   , pleft , ql  , hl  , betal
         REAL(KIND=RP) :: ur  , vr   , wr   , pright, qr  , hr  , betar
         REAL(KIND=RP) :: rtd , utd  , vtd  , wtd   , htd , atd2, atd, qtd
         REAL(KIND=RP) :: dw1 , sp1  , sp1m , hd1m  , eta1, udw1, rql
         REAL(KIND=RP) :: dw4 , sp4  , sp4p , hd4   , eta4, udw4, rqr
         REAL(KIND=RP)                   :: ds = 1.0_RP
         REAL(KIND=RP) :: rhotheta, rhothetan, thetal, thetar, ttd

         associate ( gamma => thermodynamics % gamma )

         rho  = QLeft(1)
         rhou = QLeft(2)
         rhov = QLeft(3)
         rhow = QLeft(4)
         rhoe = QLeft(5)
         rhotheta= QLeft(6)

         rhon  = QRight(1)
         rhoun = QRight(2)
         rhovn = QRight(3)
         rhown = QRight(4)
         rhoen = QRight(5)
         rhothetan= QRight(6)
         
         ul = rhou/rho 
         vl = rhov/rho 
         wl = rhow/rho 
         pleft = (gamma-1._RP)*(rhoe - 0.5_RP/rho*                        &
        &                           (rhou**2 + rhov**2 + rhow**2 )) 
         thetal = rhotheta/rho

!
         ur = rhoun/rhon 
         vr = rhovn/rhon 
         wr = rhown/rhon 
         pright = (gamma-1._RP)*(rhoen - 0.5_RP/rhon*                    &
        &                           (rhoun**2 + rhovn**2+ rhown**2)) 
         thetar = rhothetan/rhon

!
         ql = nHat(1)*ul + nHat(2)*vl + nHat(3)*wl
         qr = nHat(1)*ur + nHat(2)*vr + nHat(3)*wr
         hl = 0.5_RP*(ul*ul + vl*vl + wl*wl) +                               &
        &                 gamma/(gamma-1._RP)*pleft/rho 
         hr = 0.5_RP*(ur*ur + vr*vr + wr*wr) +                               &
        &                  gamma/(gamma-1._RP)*pright/rhon 
!
!        ---------------------
!        Square root averaging  
!        ---------------------
!
         rtd = sqrt(rho*rhon) 
         betal = rho/(rho + rtd) 
         betar = 1._RP - betal 
         utd = betal*ul + betar*ur 
         vtd = betal*vl + betar*vr 
         wtd = betal*wl + betar*wr 
         htd = betal*hl + betar*hr 
         atd2 = (gamma-1._RP)*(htd - 0.5_RP*(utd*utd + vtd*vtd + wtd*wtd)) 
         atd = sqrt(atd2) 
         qtd = utd*nHat(1) + vtd*nHat(2)  + wtd*nHat(3)
!
         IF(qtd >= 0.0_RP)     THEN
   
            dw1 = 0.5_RP*((pright - pleft)/atd2 - (qr - ql)*rtd/atd) 
            sp1 = qtd - atd 
            sp1m = min(sp1,0.0_RP) 
            hd1m = ((gamma+1._RP)/4._RP*atd/rtd)*dw1 
            eta1 = max(-abs(sp1) - hd1m,0.0_RP) 
            udw1 = dw1*(sp1m - 0.5_RP*eta1) 
            rql = rho*ql 
            flux(1) = ds*(rql + udw1) 
            flux(2) = ds*(rql*ul + pleft*nHat(1) + udw1*(utd - atd*nHat(1))) 
            flux(3) = ds*(rql*vl + pleft*nHat(2) + udw1*(vtd - atd*nHat(2))) 
            flux(4) = ds*(rql*wl + pleft*nHat(3) + udw1*(wtd - atd*nHat(3))) 
            flux(5) = ds*(rql*hl + udw1*(htd - qtd*atd)) 
            flux(6) = flux(1)*(thetar+thetal)*0.5_RP - abs(flux(1))*(thetar-thetal)*0.5_RP

         ELSE 
   
            dw4 = 0.5_RP*((pright - pleft)/atd2 + (qr - ql)*rtd/atd) 
            sp4 = qtd + atd 
            sp4p = max(sp4,0.0_RP) 
            hd4 = ((gamma+1._RP)/4._RP*atd/rtd)*dw4 
            eta4 = max(-abs(sp4) + hd4,0.0_RP) 
            udw4 = dw4*(sp4p + 0.5_RP*eta4) 
            rqr = rhon*qr 
            flux(1) = ds*(rqr - udw4) 
            flux(2) = ds*(rqr*ur + pright*nHat(1) - udw4*(utd + atd*nHat(1))) 
            flux(3) = ds*(rqr*vr + pright*nHat(2) - udw4*(vtd + atd*nHat(2))) 
            flux(4) = ds*(rqr*wr + pright*nHat(3) - udw4*(wtd + atd*nHat(3))) 
            flux(5) = ds*(rqr*hr - udw4*(htd + qtd*atd)) 
            flux(6) = flux(1)*(thetar+thetal)*0.5_RP - abs(flux(1))*(thetar-thetal)*0.5_RP

         ENDIF

         end associate
         
      END SUBROUTINE RoeRiemannSolver

!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RusanovRiemannSolver( QLeft, QRight, nHat, t1, t2, flux )
      
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: nHat(1:NDIM)
         real(kind=RP), intent(in)       :: t1(1:NDIM)
         real(kind=RP), intent(in)       :: t2(1:NDIM)
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         REAL(KIND=RP) :: rho , rhou , rhov , rhow  , rhoe
         REAL(KIND=RP) :: rhon, rhoun, rhovn, rhown , rhoen
         REAL(KIND=RP) :: ul  , vl   , wl   , pleft , ql  , hl  , betal, al, al2
         REAL(KIND=RP) :: ur  , vr   , wr   , pright, qr  , hr  , betar, ar, ar2
         REAL(KIND=RP) :: rtd , utd  , vtd  , wtd   , htd , atd2, atd, qtd
         REAL(KIND=RP) :: dw1 , sp1  , sp1m , hd1m  , eta1, udw1, rql
         REAL(KIND=RP) :: dw4 , sp4  , sp4p , hd4   , eta4, udw4, rqr
         REAL(KIND=RP)                   :: ds = 1.0_RP
         REAL(KIND=RP) :: rhotheta, rhothetan, thetal, thetar
         REAL(KIND=RP) :: smax, smaxL, smaxR
         REAL(KIND=RP) :: Leigen(2), Reigen(2)
      
         associate ( gamma => thermodynamics % gamma ) 

         rho  = QLeft(1)
         rhou = QLeft(2)
         rhov = QLeft(3)
         rhow = QLeft(4)
         rhoe = QLeft(5)
         rhotheta= QLeft(6)

         rhon  = QRight(1)
         rhoun = QRight(2)
         rhovn = QRight(3)
         rhown = QRight(4)
         rhoen = QRight(5)
         rhothetan= QRight(6)

         ul = rhou/rho 
         vl = rhov/rho 
         wl = rhow/rho 
         pleft = (gamma-1._RP)*(rhoe - 0.5_RP/rho*                        &
        &                           (rhou**2 + rhov**2 + rhow**2 )) 
         thetal = rhotheta/rho
!        
         ur = rhoun/rhon 
         vr = rhovn/rhon 
         wr = rhown/rhon 
         pright = (gamma-1._RP)*(rhoen - 0.5_RP/rhon*                    &
        &                           (rhoun**2 + rhovn**2+ rhown**2)) 
         thetar = rhothetan/rhon

!
         ql = nHat(1)*ul + nHat(2)*vl + nHat(3)*wl
         qr = nHat(1)*ur + nHat(2)*vr + nHat(3)*wr
         hl = 0.5_RP*(ul*ul + vl*vl + wl*wl) +                               &
        &                 gamma/(gamma-1._RP)*pleft/rho 
         hr = 0.5_RP*(ur*ur + vr*vr + wr*wr) +                               &
        &                  gamma/(gamma-1._RP)*pright/rhon 
!
!        ---------------------
!        Square root averaging  
!        ---------------------
!
         rtd = sqrt(rho*rhon) 
         betal = rho/(rho + rtd) 
         betar = 1._RP - betal 
         utd = betal*ul + betar*ur 
         vtd = betal*vl + betar*vr 
         wtd = betal*wl + betar*wr 
         htd = betal*hl + betar*hr 
         atd2 = (gamma-1._RP)*(htd - 0.5_RP*(utd*utd + vtd*vtd + wtd*wtd)) 
         atd = sqrt(atd2) 
         qtd = utd*nHat(1) + vtd*nHat(2)  + wtd*nHat(3)
         !Rusanov
         ar2 = (gamma-1.d0)*(hr - 0.5d0*(ur*ur + vr*vr + wr*wr)) 
         al2 = (gamma-1.d0)*(hl - 0.5d0*(ul*ul + vl*vl + wl*wl)) 
         ar = sqrt(ar2)
         al = sqrt(al2)
!           
         rql = rho*ql 
         rqr = rhon*qr             
         flux(1) = ds*(rql + rqr) 
         flux(2) = ds*(rql*ul + pleft*nHat(1) + rqr*ur + pright*nHat(1)) 
         flux(3) = ds*(rql*vl + pleft*nHat(2) + rqr*vr + pright*nHat(2))
         flux(4) = ds*(rql*wl + pleft*nHat(3) + rqr*wr + pright*nHat(3)) 
         flux(5) = ds*(rql*hl + rqr*hr) 
         flux(6) = flux(1)*(thetar+thetal)*0.5_RP - abs(flux(1))*(thetar-thetal)*0.5_RP

         smax = MAX(ar+ABS(qr),al+ABS(ql))

         flux = (flux - ds*smax*(QRight-QLeft))/2.d0

         RETURN 

         end associate
         
      END SUBROUTINE RusanovRiemannSolver           
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
!        Averaged states functions
!        -------------------------
!
!     To this averages, the states QLeft and QRight velocities have already been rotated,
!  so that only the first flux component has to be computed (Rotational invariance, see
!  E.F. Toro - Riemann solvers and numerical methods for fluid dynamics. Page 105).
!
!  Implemented two-point averages are:
!     -> Standard: Central fluxes
!     -> Ducros
!     -> Morinishi
!     -> Kennedy and Gruber
!     -> Pirozzoli
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine StandardAverage(QLeft, QRight, pL, pR, invRhoL, invRhoR, flux) 
!
!        *********************************************************************
!           Computes the standard average of the two states:
!              F* = {{F}} = 0.5 * (FL + FR)
!
!           State vectors are rotated.
!        *********************************************************************
!
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: uL, vL, wL 
         real(kind=RP)     :: uR, vR, wR

         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = 0.5_RP * ( QLeft(IRHOU) + QRight(IRHOU) )
         flux(IRHOU) = 0.5_RP * ( QLeft(IRHOU) * uL + QRight(IRHOU) * uR + pL + pR )
         flux(IRHOV) = 0.5_RP * ( QLeft(IRHOU) * vL + QRight(IRHOU) * vR )
         flux(IRHOW) = 0.5_RP * ( QLeft(IRHOU) * wL + QRight(IRHOU) * wR )
         flux(IRHOE) = 0.5_RP * ( uL*(QLeft(IRHOE) + pL) + uR*(QRight(IRHOE) + pR) )
         flux(IRHOTHETA) = 0.5_RP * (uL * QLeft(IRHOTHETA) + uR * QRight(IRHOTHETA)  )
      end subroutine StandardAverage

      subroutine DucrosAverage(QLeft, QRight, pL, pR, invRhoL, invRhoR, flux) 
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: uL, vL, wL
         real(kind=RP)     :: uR, vR, wR

         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = 0.25_RP * ( QLeft(IRHO) + QRight(IRHO) ) * (uL + uR)
         flux(IRHOU) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * (uL + uR) + 0.5_RP * (pL + pR)
         flux(IRHOV) = 0.25_RP * ( QLeft(IRHOV) + QRight(IRHOV) ) * (uL + uR)
         flux(IRHOW) = 0.25_RP * ( QLeft(IRHOW) + QRight(IRHOW) ) * (uL + uR)
         flux(IRHOE) = 0.25_RP * ( QLeft(IRHOE) + pL + QRight(IRHOE) + pR ) * (uL + uR)

      end subroutine DucrosAverage

      subroutine MorinishiAverage(QLeft,QRight, pL, pR, invRhoL, invRhoR, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: uL, vL, wL, hL
         real(kind=RP)     :: uR, vR, wR, hR

         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Here the enthalpy does not contain the kinetic energy
!        -----------------------------------------------------
         hL = dimensionless % cp * pL  ; hR = dimensionless % cp * pR 
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = 0.5_RP * ( QLeft(IRHOU) + QRight(IRHOU) )
         flux(IRHOU) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * ( uL + uR ) + 0.5_RP * ( pL + pR )
         flux(IRHOV) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * ( vL + vR )
         flux(IRHOW) = 0.25_RP * ( QLeft(IRHOU) + QRight(IRHOU) ) * ( wL + wR )
         flux(IRHOE) = 0.5_RP * ( uL*hL + uR*hR) + 0.25_RP * ( QLeft(IRHOU)*uL + QRight(IRHOU)*uR ) * ( uL + uR ) &
                                                 + 0.25_RP * ( QLeft(IRHOU)*vL + QRight(IRHOU)*vR ) * ( vL + vR ) &
                                                 + 0.25_RP * ( QLeft(IRHOU)*wL + QRight(IRHOU)*wR ) * ( wL + wR ) &
                                                 - 0.25_RP * ( QLeft(IRHOU)*POW2(uL) + QRight(IRHOU)*POW2(uR)   ) &
                                                 - 0.25_RP * ( QLeft(IRHOU)*POW2(vL) + QRight(IRHOU)*POW2(vR)   ) &
                                                 - 0.25_RP * ( QLeft(IRHOU)*POW2(wL) + QRight(IRHOU)*POW2(wR)   ) 

      end subroutine MorinishiAverage

      subroutine KennedyGruberAverage(QLeft,QRight, pL, pR, invRhoL, invRhoR, flux) 
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: uL, vL, wL
         real(kind=RP)     :: uR, vR, wR
         real(kind=RP)     :: rho, u, v, w, e, p

         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Compute KG averages
!        -------------------   
         rho = 0.5_RP * (QLeft(IRHO) + QRight(IRHO))
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)
         e   = 0.5_RP * (QLeft(IRHOE)*invRhoL + QRight(IRHOE)*invRhoR)
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = rho * u
         flux(IRHOU) = rho * u * u + p
         flux(IRHOV) = rho * u * v
         flux(IRHOW) = rho * u * w
         flux(IRHOE) = rho * u * e + p * u
         
      end subroutine KennedyGruberAverage

      subroutine PirozzoliAverage(QLeft,QRight, pL, pR, invRhoL, invRhoR, flux)
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: uL, vL, wL
         real(kind=RP)     :: uR, vR, wR
         real(kind=RP)     :: rho, u, v, w, h, p

         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
   
         rho = 0.5_RP * (QLeft(IRHO) + QRight(IRHO))
         u   = 0.5_RP * (uL + uR)
         v   = 0.5_RP * (vL + vR)
         w   = 0.5_RP * (wL + wR)
         p   = 0.5_RP * (pL + pR)
         h   = 0.5_RP * ((QLeft(IRHOE)+pL)*invRhoL + (QRight(IRHOE)+pR)*invRhoR)
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = rho * u
         flux(IRHOU) = rho * u * u + p
         flux(IRHOV) = rho * u * v
         flux(IRHOW) = rho * u * w
         flux(IRHOE) = rho * u * h
         
      end subroutine PirozzoliAverage

      subroutine EntropyConservingAverage(QLeft,QRight, pL, pR, invRhoL, invRhoR, flux)
         use SMConstants
         use Utilities, only: logarithmicMean
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: rhoL, uL, vL, wL
         real(kind=RP)     :: rhoR, uR, vR, wR
         real(kind=RP)     :: rho, u, v, w, h, p, p2
         real(kind=RP)     :: zL(NCONS), zR(NCONS), zSum(NCONS), invZ1Sum
         real(kind=RP)     :: z5Log, z1Log

         associate ( gammaPlus1Div2      => thermodynamics % gammaPlus1Div2, &
                     gammaMinus1Div2     => thermodynamics % gammaMinus1Div2, &
                     gammaDivGammaMinus1 => thermodynamics % gammaDivGammaMinus1, &
                     invGamma            => thermodynamics % invGamma ) 

         rhoL = QLeft(IRHO)               ; rhoR = QRight(IRHO)
         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Compute Ismail and Roe parameter vector
!        ---------------------------------------
         zL(5) = sqrt(rhoL*pL)      ; zR(5) = sqrt(rhoR*pR)
         zL(1) = rhoL / zL(5)       ; zR(1) = rhoR / zR(5)
         zL(2) = zL(1) * uL         ; zR(2) = zR(1) * uR
         zL(3) = zL(1) * vL         ; zR(3) = zR(1) * vR
         zL(4) = zL(1) * wL         ; zR(4) = zR(1) * wR

         zSum = zL + zR
         invZ1Sum = 1.0_RP / zSum(1)

         call logarithmicMean(zL(1),zR(1), z1Log)
         call logarithmicMean(zL(5),zR(5), z5Log)

         rho = 0.5_RP * zSum(1) * z5Log
         u   = zSum(2) * invZ1Sum
         v   = zSum(3) * invZ1Sum
         w   = zSum(4) * invZ1Sum
         p   = zSum(5) * invZ1Sum
         p2  = (gammaPlus1Div2 * z5Log / z1Log + gammaMinus1Div2 * p) * invGamma
         h   = gammaDivGammaMinus1 * p2 / rho + 0.5_RP*(POW2(u) + POW2(v) + POW2(w))
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = rho * u
         flux(IRHOU) = rho * u * u + p
         flux(IRHOV) = rho * u * v
         flux(IRHOW) = rho * u * w
         flux(IRHOE) = rho * u * h
         
         end associate

      end subroutine EntropyConservingAverage

      subroutine ChandrasekarAverage(QLeft,QRight, pL, pR, invRhoL, invRhoR, flux)
         use SMConstants
         use Utilities, only: logarithmicMean
         implicit none
         real(kind=RP), intent(in)       :: QLeft(1:NCONS)
         real(kind=RP), intent(in)       :: QRight(1:NCONS)
         real(kind=RP), intent(in)       :: pL, pR
         real(kind=RP), intent(in)       :: invRhoL, invRhoR
         real(kind=RP), intent(out)      :: flux(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: rhoL, uL, vL, wL, betaL
         real(kind=RP)     :: rhoR, uR, vR, wR, betaR
         real(kind=RP)     :: rho, u, v, w, h, p, betaLog

         associate ( gammaMinus1 => thermodynamics % gammaMinus1 ) 

         rhoL = QLeft(IRHO)               ; rhoR = QRight(IRHO)
         uL = invRhoL * QLeft(IRHOU)      ; uR = invRhoR * QRight(IRHOU)
         vL = invRhoL * QLeft(IRHOV)      ; vR = invRhoR * QRight(IRHOV)
         wL = invRhoL * QLeft(IRHOW)      ; wR = invRhoR * QRight(IRHOW)
!
!        Compute Chandrasekar's variables
!        --------------------------------
         betaL = 0.5_RP * rhoL / pL    ; betaR = 0.5_RP * rhoR / pR
         call logarithmicMean(betaL, betaR, betaLog)

         call logarithmicMean(rhoL,rhoR,rho)
         u   = AVERAGE(uL, uR)
         v   = AVERAGE(vL, vR)
         w   = AVERAGE(wL, wR)
         p   = 0.5_RP * (rhoL + rhoR) / (betaL + betaR)
         h   =   0.5_RP/(betaLog*(gammaMinus1)) &
               - 0.5_RP*AVERAGE(POW2(uL)+POW2(vL)+POW2(wL), POW2(uR)+POW2(vR)+POW2(wR)) &
               + p/rho + POW2(u) + POW2(v) + POW2(w)
!
!        Compute the flux
!        ----------------
         flux(IRHO)  = rho * u
         flux(IRHOU) = rho * u * u + p
         flux(IRHOV) = rho * u * v
         flux(IRHOW) = rho * u * w
         flux(IRHOE) = rho * u * h

         end associate

      end subroutine ChandrasekarAverage

end module RiemannSolvers_NSSA
