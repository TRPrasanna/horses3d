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
#include "Includes.h"
      Module PhysicsKeywordsModule
         IMPLICIT NONE 
         INTEGER, PARAMETER :: KEYWORD_LENGTH = 132
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: MACH_NUMBER_KEY           = "mach number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: REYNOLDS_NUMBER_KEY       = "reynolds number"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_THETA_KEY             = "aoa theta"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: AOA_PHI_KEY               = "aoa phi"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: FLOW_EQUATIONS_KEY        = "flow equations"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: RIEMANN_SOLVER_NAME_KEY   = "riemann solver"
         
         CHARACTER(LEN=KEYWORD_LENGTH), DIMENSION(2) :: physicsKeywords = [MACH_NUMBER_KEY, FLOW_EQUATIONS_KEY]
         
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: ROE_SOLVER_NAME           = "roe"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: RUSANOV_SOLVER_NAME       = "rusanov"
         CHARACTER(LEN=KEYWORD_LENGTH), PARAMETER :: LAXFRIEDRICHS_SOLVER_NAME = "lax friedrichs"
         
      END MODULE PhysicsKeywordsModule
!
!////////////////////////////////////////////////////////////////////////
!    
!    ******
     MODULE PhysicsStorage
!    ******
!
     USE SMConstants
     use FluidData
     
     IMPLICIT NONE

     SAVE
!
!    ----------------------------
!    Either NavierStokes or Euler
!    ----------------------------
!
     LOGICAL :: flowIsNavierStokes = .true.
!
!    --------------------------
!!   The sizes of the NS system
!    --------------------------
!
     INTEGER, PARAMETER :: N_EQN = 5, N_GRAD_EQN = 4
!
!    -----------------------------
!    Number of physical dimensions
!    -----------------------------
!
     INTEGER, PARAMETER       :: NDIM = 3
     INTEGER, PARAMETER       :: IX = 1 , IY = 2 , IZ = 3
!
!    -------------------------------------------
!!   The positions of the conservative variables
!    -------------------------------------------
!
     INTEGER, PARAMETER       :: NCONS = 5
     INTEGER, PARAMETER       :: IRHO = 1 , IRHOU = 2 , IRHOV = 3 , IRHOW = 4 , IRHOE = 5
!
!    ---------------------------------------
!!   The positions of the gradient variables
!    ---------------------------------------
!
     INTEGER, PARAMETER  :: IGU = 1 , IGV = 2 , IGW = 3 , IGT = 4
!
!    --------------------------------------------
!!   The temperature scale in the Sutherland law:
!!   198.6 for temperatures in R, 110.3 for
!!   temperatures in K.
!    --------------------------------------------
!
     REAL( KIND=RP ) :: TScale
!
!    ------------------------------------------------
!!   The ratio of the scale and reference tempartures
!    ------------------------------------------------
!
     REAL( KIND=RP ) :: TRatio 
!    ----------------------------------
!
!    ------------------------------------
!    Riemann solver associated quantities
!    ------------------------------------
!
     INTEGER, PARAMETER :: ROE = 0, LXF = 1, RUSANOV = 2
     INTEGER            :: riemannSolverChoice = ROE

     type(Thermodynamics_t), target, private :: ThermodynamicsAir = Thermodynamics_t( &
                                                              "Air", & ! Name
                                    287.15_RP * 5.0_RP / 9.0_RP, & ! R
                                                         1.4_RP, & ! gamma
                                                   sqrt(1.4_RP), & ! sqrtGamma
                                                1.4_RP - 1.0_RP, & ! gammaMinus1         
                                     (1.4_RP - 1.0_RP) / 2.0_RP, & ! gammaMinus1Div2
                                     (1.4_RP + 1.0_RP) / 2.0_RP, & ! gammaPlus1Div2
                    (1.4_RP - 1.0_RP) / (2.0_RP * sqrt(1.4_RP)), & ! gammaMinus1Div2sg 
                          (1.4_RP - 1.0_RP) / (2.0_RP * 1.4_RP), & ! gammaMinus1Div2g 
                                     2.0_RP / (1.4_RP + 1.0_RP), & ! InvGammaPlus1Div2 
                                     1.0_RP / (1.4_RP - 1.0_RP), & ! InvGammaMinus1
                                                1.0_RP / 1.4_RP, & ! InvGamma
                                   1.4_RP / ( 1.4_RP - 1.0_RP ), & ! gammaDivGammaMinus1
     287.15_RP * 5.0_RP / 9.0_RP * 1.4_RP / ( 1.4_RP - 1.0_RP ), & ! cp
              287.15_RP * 5.0_RP / 9.0_RP / ( 1.4_RP - 1.0_RP ), & ! cp
                                                         0.0_RP  & ! Bulk viscosity ratio
)
!
!    ========
     CONTAINS
!    ========
!
!     ///////////////////////////////////////////////////////
!
!     --------------------------------------------------
!!    Constructor: Define default values for the physics
!!    variables.
!     --------------------------------------------------
!
      SUBROUTINE ConstructPhysicsStorage( controlVariables, success )
      USE FTValueDictionaryClass
      USE PhysicsKeywordsModule
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FTValueDictionary) :: controlVariables
      LOGICAL                 :: success
!
!     ---------------
!     Local variables
!     ---------------
!
      CHARACTER(LEN=KEYWORD_LENGTH) :: keyword
      type(Thermodynamics_t), pointer  :: thermodynamics_
      type(RefValues_t)                :: refValues_
      type(Dimensionless_t)            :: dimensionless_
!
!     --------------------
!     Collect input values
!     --------------------
!
      success = .TRUE.
      CALL CheckPhysicsInputIntegrity(controlVariables,success)
      IF(.NOT. success) RETURN 
!
!
!     ---------------------
!     Set the gas to be air
!     ---------------------
!
      thermodynamics_ => thermodynamicsAir
!
!     ------------------------
!     Dimensionless quantities
!     ------------------------
!
      dimensionless_ % cp = thermodynamics_ % gamma * thermodynamics_ % InvGammaMinus1
      dimensionless_ % cv = thermodynamics_ % InvGammaMinus1
      dimensionless_ % Mach = controlVariables % doublePrecisionValueForKey(MACH_NUMBER_KEY)
      dimensionless_ % Pr   = 0.72_RP
      dimensionless_ % Fr   = 0.0_RP

      keyword = controlVariables % stringValueForKey(FLOW_EQUATIONS_KEY,KEYWORD_LENGTH)
      CALL toLower(keyword)
      IF ( keyword == "euler" )     THEN
         flowIsNavierStokes = .FALSE.
         dimensionless_ % Re = 0.0_RP     ! Just to avoid non-initialized variables
         dimensionless_ % mu = 0.0_RP
         dimensionless_ % kappa = 0.0_RP
      ELSE 
         flowIsNavierStokes = .TRUE.
         IF ( controlVariables % containsKey(REYNOLDS_NUMBER_KEY) )     THEN
            dimensionless_ % Re = controlVariables % doublePrecisionValueForKey(REYNOLDS_NUMBER_KEY) 
         ELSE 
            PRINT *, "Input file is missing entry for keyword: ",REYNOLDS_NUMBER_KEY
            success = .FALSE.
            RETURN 
         END IF 

         dimensionless_ % mu   = 1.0_RP / dimensionless_ % Re
         dimensionless_ % kappa = 1.0_RP / ( thermodynamics_ % gammaMinus1 * &
                                              POW2( dimensionless_ % Mach) * &
                                      dimensionless_ % Re * dimensionless_ % Pr )
      END IF 

      dimensionless_ % gammaM2 = thermodynamics_ % gamma * POW2( dimensionless_ % Mach )
      dimensionless_ % invFroudeSquare = 0.0_RP
!
!     --------------------
!     Set reference values: TODO read from parameter file
!     --------------------
!
      refValues_ % L = 1.0_RP       ! m
      refValues_ % T = 520.0_RP     ! Rankine
      refValues_ % rho = 101325.0_RP / (thermodynamics_ % R * refValues_ % T)
      refValues_ % V =   dimensionless_ % Mach &
                       * sqrt( thermodynamics_ % gamma * thermodynamics_ % R * refValues_ % T )
      refValues_ % p = refValues_ % rho * POW2( refValues_ % V )
      if ( flowIsNavierStokes ) then
         refValues_ % mu = refValues_ % rho * refValues_ % V * refValues_ % L / dimensionless_ % Re
         refValues_ % kappa = refValues_ % mu * thermodynamics_ % cp / dimensionless_ % Pr

      else
         refValues_ % mu = 0.0_RP
         refValues_ % kappa = 0.0_RP
      
      end if

      refValues_ % time = refValues_ % L / refValues_ % V

!
!     --------------------------------------------------------------------
!     The riemann solver is also optional. Set it to Roe if not requested.
!     --------------------------------------------------------------------
!
      IF ( controlVariables % containsKey(RIEMANN_SOLVER_NAME_KEY) )     THEN
         keyword = controlVariables % stringValueForKey(key             = RIEMANN_SOLVER_NAME_KEY,&
                                                        requestedLength = KEYWORD_LENGTH)
         CALL toLower(keyword)
         SELECT CASE ( keyword )
            CASE( ROE_SOLVER_NAME ) 
               riemannSolverChoice = ROE
            CASE( LAXFRIEDRICHS_SOLVER_NAME )
               riemannSolverChoice = LXF 
            CASE( RUSANOV_SOLVER_NAME )
               riemannSolverChoice = RUSANOV
            CASE DEFAULT 
               PRINT *, "Unknown Riemann solver choice: ", TRIM(keyword), ". Defaulting to Roe"
               riemannSolverChoice = ROE
         END SELECT 
      ELSE 
         PRINT *, "Input file is missing keyword 'riemann solver'. Using Roe by default"
         riemannSolverChoice = ROE 
      END IF 
!
!     ------------------------------------------------------------------------------
!     The angle of attack parameters are optional. If not present, set them to zero.
!     ------------------------------------------------------------------------------
!
      IF ( controlVariables % containsKey(AOA_PHI_KEY) )     THEN
         refValues_ % AOAPhi = controlVariables % doublePrecisionValueForKey(AOA_PHI_KEY) 
      ELSE
         refValues_ % AOAPhi = 0.0_RP
      END IF 
      IF ( controlVariables % containsKey(AOA_THETA_KEY) )     THEN
         refValues_ % AOATheta = controlVariables % doublePrecisionValueForKey(AOA_THETA_KEY) 
      ELSE
         refValues_ % AOATheta = 0.0_RP
      END IF 
!
!     --------------------------
!     Sutherland's law constants
!     --------------------------
!
      TScale          = 198.6_RP
      TRatio          = TScale/ refValues_ % T

      call setThermodynamics( thermodynamics_ )
      call setDimensionless( dimensionless_ )
      call setRefValues( refValues_ )


      CALL DescribePhysicsStorage()
!
      END SUBROUTINE ConstructPhysicsStorage
!
!     ///////////////////////////////////////////////////////
!
!     -------------------------------------------------
!!    Destructor: Does nothing for this storage
!     -------------------------------------------------
!
      SUBROUTINE DestructPhysicsStorage
      
      END SUBROUTINE DestructPhysicsStorage
!
!     //////////////////////////////////////////////////////
!
!     -----------------------------------------
!!    Descriptor: Shows the gathered data
!     -----------------------------------------
!
      SUBROUTINE DescribePhysicsStorage()
         USE Headers
         IMPLICIT NONE
         real(kind=RP)  :: pRef

         pRef = thermodynamics % R * refValues % rho * refValues % T

         write(STD_OUT,'(/,/)')
         if (flowIsNavierStokes) then
            call Section_Header("Loading Navier-Stokes physics")
         else
            call Section_Header("Loading Euler physics")
         end if

         write(STD_OUT,'(/)')
         call SubSection_Header("Fluid data")
         write(STD_OUT,'(30X,A,A22,A10)') "->" , "Gas: " , "Air"
         write(STD_OUT,'(30X,A,A22,F10.3,A)') "->" , "State constant: " , thermodynamics % R, " I.S."
         write(STD_OUT,'(30X,A,A22,F10.3)') "->" , "Specific heat ratio: " , thermodynamics % gamma

         write(STD_OUT,'(/)')
         call SubSection_Header("Reference quantities")
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference Temperature: " , refValues % T, " K."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference pressure: " , pRef, " Pa."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference density: " , refValues % rho , " kg/m^3."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference velocity: " , refValues % V , " m/s."
         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reynolds length: " , refValues % L , " m."
         
         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference viscosity: ",refValues % mu , " Pa·s."
            write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference conductivity: ", refValues % kappa, " W/(m·K)."
         end if

         write(STD_OUT,'(30X,A,A30,F10.3,A)') "->" , "Reference time: ", refValues % time, " s."

         write(STD_OUT,'(/)')
         call SubSection_Header("Dimensionless quantities")
         write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Mach number: " , dimensionless % Mach
         if ( flowIsNavierStokes ) then
            write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Reynolds number: " , dimensionless % Re
            write(STD_OUT,'(30X,A,A20,F10.3)') "->" , "Prandtl number: " , dimensionless % Pr
         end if

      END SUBROUTINE DescribePhysicsStorage
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckPhysicsInputIntegrity( controlVariables, success )  
         USE FTValueDictionaryClass
         USE PhysicsKeywordsModule
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(FTValueDictionary) :: controlVariables
         LOGICAL                 :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(FTObject), POINTER :: obj
         INTEGER                  :: i
         success = .TRUE.
         
         DO i = 1, SIZE(physicsKeywords)
            obj => controlVariables % objectForKey(physicsKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",physicsKeywords(i)
               success = .FALSE. 
            END IF  
         END DO  
         
      END SUBROUTINE CheckPhysicsInputIntegrity
!
!    **********       
     END MODULE PhysicsStorage
!    **********
!@mark -
!
!  **************
   Module Physics 
!  **************
!
      USE SMConstants
      USE PhysicsStorage
      IMPLICIT NONE
!
!     ---------
!     Constants
!     ---------
!
      INTEGER, PARAMETER   :: WALL_BC = 1, RADIATION_BC = 2
      REAL(KIND=RP)        :: waveSpeed
      INTEGER              :: boundaryCondition(4), bcType


!
!    ---------------
!    Interface block
!    ---------------
!
     interface GradientValuesForQ
         module procedure GradientValuesForQ_0D , GradientValuesForQ_3D
     end interface GradientValuesForQ

     interface InviscidFlux
         module procedure InviscidFlux0D , InviscidFlux1D , InviscidFlux2D , InviscidFlux3D
     end interface InviscidFlux

     interface ViscousFlux
         module procedure ViscousFlux0D , ViscousFlux1D , ViscousFlux2D , ViscousFlux3D
     end interface ViscousFlux
!
!     ========
      CONTAINS 
!     ========
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RiemannSolver( QLeft, QRight, nHat, flux )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN)  :: Qleft, Qright, flux
         REAL(KIND=RP), DIMENSION(3)      :: nHat
         
         SELECT CASE ( riemannSolverChoice )
            CASE ( ROE )
               CALL RoeSolver( QLeft, QRight, nHat, flux )
            CASE (LXF)
               PRINT *, "3D LXF to be implemented."
               STOP !DEBUG
               CALL LxFSolver( QLeft, QRight, nHat, flux )
            CASE (RUSANOV)
               CALL RusanovSolver( QLeft, QRight, nHat, flux )               
            CASE DEFAULT
               PRINT *, "Undefined choice of Riemann Solver. Abort"
               STOP
         END SELECT

      
      END SUBROUTINE RiemannSolver
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RoeSolver( QLeft, QRight, nHat, flux )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN) :: Qleft, Qright, flux
         REAL(KIND=RP), DIMENSION(3)     :: nHat
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
      
         associate ( gamma => thermodynamics % gamma )
            
         rho  = Qleft(1)
         rhou = Qleft(2)
         rhov = Qleft(3)
         rhow = Qleft(4)
         rhoe = Qleft(5)
   
         rhon  = Qright(1)
         rhoun = Qright(2)
         rhovn = Qright(3)
         rhown = Qright(4)
         rhoen = Qright(5)
         
         ul = rhou/rho 
         vl = rhov/rho 
         wl = rhow/rho 
         pleft = (gamma-1._RP)*(rhoe - 0.5_RP/rho*                        &
        &                           (rhou**2 + rhov**2 + rhow**2 )) 
!
         ur = rhoun/rhon 
         vr = rhovn/rhon 
         wr = rhown/rhon 
         pright = (gamma-1._RP)*(rhoen - 0.5_RP/rhon*                    &
        &                           (rhoun**2 + rhovn**2+ rhown**2)) 
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
         ENDIF

         end associate
         
      END SUBROUTINE RoeSolver
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LxFSolver( QLeft, QRight, nHat, flux ) 
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN) :: Qleft, Qright, flux
         REAL(KIND=RP), DIMENSION(3)     :: nHat
         REAL(KIND=RP)                   :: ds = 1.0_RP
!
!        ---------------
!        Local Variables
!        ---------------
!
!
      REAL(KIND=RP) :: rho , rhou , rhov  , rhoe
      REAL(KIND=RP) :: rhon, rhoun, rhovn , rhoen
      REAL(KIND=RP) :: ul  , vl   , pleft , ql, cl
      REAL(KIND=RP) :: ur  , vr   , pright, qr, cr
      REAL(KIND=RP) :: sM
      REAL(KIND=RP), DIMENSION(N_EQN) :: FL, FR
      
      associate ( gamma => thermodynamics % gamma ) 

      rho  = Qleft(1)
      rhou = Qleft(2)
      rhov = Qleft(3)
      rhoe = Qleft(4)

      rhon  = Qright(1)
      rhoun = Qright(2)
      rhovn = Qright(3)
      rhoen = Qright(4)

      ul = rhou/rho 
      vl = rhov/rho 
      pleft = (gamma-1.d0)*(rhoe - 0.5d0/rho*(rhou**2 + rhov**2)) 
!
      ur = rhoun/rhon 
      vr = rhovn/rhon 
      pright = (gamma-1.d0)*(rhoen - 0.5d0/rhon*(rhoun**2 + rhovn**2)) 
!
      ql = nHat(1)*ul + nHat(2)*vl 
      qr = nHat(1)*ur + nHat(2)*vr 
      cl = SQRT( gamma*pleft/rho )
      cr = SQRT( gamma*pright/rhon )
!
      FL(1) = rho*ql
      FL(2) = rhou*ql + pleft*nHat(1)
      FL(3) = rhov*ql + pleft*nHat(2)
      FL(4) = (rhoe + pleft)*ql
!
      FR(1) = rhon*qr
      FR(2) = rhoun*qr + pright*nHat(1)
      FR(3) = rhovn*qr + pright*nHat(2)
      FR(4) = (rhoen + pright)*qr
!
      sM = MAX( ABS(ql) + cl, ABS(qr) + cr )
!
      flux = ds * 0.5_RP * ( FL + FR - sM*(Qright - Qleft) )      

      end associate
         
      END SUBROUTINE LxFSolver
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RusanovSolver( QLeft, QRight, nHat, flux )
      
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN) :: Qleft, Qright, flux
         REAL(KIND=RP), DIMENSION(3)     :: nHat
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
         
         REAL(KIND=RP) :: smax, smaxL, smaxR
         REAL(KIND=RP) :: Leigen(2), Reigen(2)
         !REAL(KIND=RP) :: gamma = 1.4_RP
      
         associate ( gamma => thermodynamics % gamma ) 

         rho  = Qleft(1)
         rhou = Qleft(2)
         rhov = Qleft(3)
         rhow = Qleft(4)
         rhoe = Qleft(5)
   
         rhon  = Qright(1)
         rhoun = Qright(2)
         rhovn = Qright(3)
         rhown = Qright(4)
         rhoen = Qright(5)
   
         ul = rhou/rho 
         vl = rhov/rho 
         wl = rhow/rho 
         pleft = (gamma-1._RP)*(rhoe - 0.5_RP/rho*                        &
        &                           (rhou**2 + rhov**2 + rhow**2 )) 
!
         ur = rhoun/rhon 
         vr = rhovn/rhon 
         wr = rhown/rhon 
         pright = (gamma-1._RP)*(rhoen - 0.5_RP/rhon*                    &
        &                           (rhoun**2 + rhovn**2+ rhown**2)) 
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
         ar = SQRT(ar2)
         al = SQRT(al2)
!           
         rql = rho*ql 
         rqr = rhon*qr             
         flux(1) = ds*(rql + rqr) 
         flux(2) = ds*(rql*ul + pleft*nHat(1) + rqr*ur + pright*nHat(1)) 
         flux(3) = ds*(rql*vl + pleft*nHat(2) + rqr*vr + pright*nHat(2))
         flux(4) = ds*(rql*wl + pleft*nHat(3) + rqr*wr + pright*nHat(3)) 
         flux(5) = ds*(rql*hl + rqr*hr) 

         smax = MAX(ar+ABS(qr),al+ABS(ql))

         flux = (flux - ds*smax*(Qright-Qleft))/2.d0

         RETURN 

         end associate
         
      END SUBROUTINE RusanovSolver           
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE xFlux( Q, f )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN) :: Q
         REAL(KIND=RP), DIMENSION(N_EQN) :: f
!
!        ---------------
!        Local Variables
!        ---------------
!
         REAL(KIND=RP) :: u, v, w, rho, rhou, rhov, rhoe, rhow, p
!      
         associate ( gammaMinus1 => thermodynamics % gammaMinus1 ) 

         rho  = Q(1)
         rhou = Q(2)
         rhov = Q(3)
         rhow = Q(4)
         rhoe = Q(5)
!
         u = rhou/rho 
         v = rhov/rho
         w = rhow/rho
         p = gammaMinus1*(rhoe - 0.5_RP*rho*(u**2 + v**2 + w**2)) 
!
         f(1) = rhou 
         f(2) = p + rhou*u 
         f(3) = rhou*v 
         f(4) = rhou*w 
         f(5) = u*(rhoe + p) 

         end associate
         
      END SUBROUTINE xFlux
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE yFlux( Q, g )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN) :: Q
         REAL(KIND=RP), DIMENSION(N_EQN) :: g
!
!        ---------------
!        Local Variables
!        ---------------
!
         REAL(KIND=RP) :: u, v, w, rho, rhou, rhov, rhoe, rhow, p
!      
         associate ( gammaMinus1 => thermodynamics % gammaMinus1 ) 

         rho  = Q(1)
         rhou = Q(2)
         rhov = Q(3)
         rhow = Q(4)
         rhoe = Q(5)
!
         u = rhou/rho 
         v = rhov/rho 
         w = rhow/rho
         p = gammaMinus1*(rhoe - 0.5_RP*rho*(u**2 + v**2 + w**2)) 
!
         g(1) = rhov 
         g(2) = rhou*v 
         g(3) = p + rhov*v 
         g(4) = rhow*v 
         g(5) = v*(rhoe + p) 

         end associate
         
      END SUBROUTINE yFlux
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE zFlux( Q, h )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN) :: Q
         REAL(KIND=RP), DIMENSION(N_EQN) :: h
!
!        ---------------
!        Local Variables
!        ---------------
!
         REAL(KIND=RP) :: u, v, w, rho, rhou, rhov, rhoe, rhow, p
!      
         associate ( gammaMinus1 => thermodynamics % gammaMinus1 ) 
  
         rho  = Q(1)
         rhou = Q(2)
         rhov = Q(3)
         rhow = Q(4)
         rhoe = Q(5)
!
         u = rhou/rho 
         v = rhov/rho 
         w = rhow/rho
         p = gammaMinus1*(rhoe - 0.5_RP*rho*(u**2 + v**2 + w**2)) 
!
         h(1) = rhow 
         h(2) = rhou*w 
         h(3) = rhov*w
         h(4) = p + rhow*w 
         h(5) = w*(rhoe + p) 

         end associate
         
      END SUBROUTINE zFlux
   
      pure function InviscidFlux0D( Q ) result ( F )
         implicit none
         real(kind=RP), intent(in)           :: Q(1:NCONS)
         real(kind=RP)           :: F(1:NCONS , 1:NDIM)
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

      end function InviscidFlux0D

      pure function InviscidFlux1D( N , Q ) result ( F )
         implicit none
         integer,       intent (in) :: N
         real(kind=RP), intent (in) :: Q(0:N , 1:NCONS)
         real(kind=RP)              :: F(0:N , 1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)           :: u(0:N) , v(0:N) , w(0:N) , p(0:N)

         associate ( gammaMinus1 => thermodynamics % gammaMinus1 ) 

         u = Q(:,IRHOU) / Q(:,IRHO)
         v = Q(:,IRHOV) / Q(:,IRHO)
         w = Q(:,IRHOW) / Q(:,IRHO)
         p = gammaMinus1 * (Q(:,IRHOE) - 0.5_RP * ( Q(:,IRHOU) * u + Q(:,IRHOV) * v + Q(:,IRHOW) * w ) )
         
         F(:,IRHO , IX ) = Q(:,IRHOU)
         F(:,IRHOU, IX ) = Q(:,IRHOU) * u + p
         F(:,IRHOV, IX ) = Q(:,IRHOU) * v
         F(:,IRHOW, IX ) = Q(:,IRHOU) * w
         F(:,IRHOE, IX ) = ( Q(:,IRHOE) + p ) * u

         F(:,IRHO , IY ) = Q(:,IRHOV)
         F(:,IRHOU ,IY ) = F(:,IRHOV,IX)
         F(:,IRHOV ,IY ) = Q(:,IRHOV) * v + p
         F(:,IRHOW ,IY ) = Q(:,IRHOV) * w
         F(:,IRHOE ,IY ) = ( Q(:,IRHOE) + p ) * v

         F(:,IRHO ,IZ) = Q(:,IRHOW)
         F(:,IRHOU,IZ) = F(:,IRHOW,IX)
         F(:,IRHOV,IZ) = F(:,IRHOW,IY)
         F(:,IRHOW,IZ) = Q(:,IRHOW) * w + P
         F(:,IRHOE,IZ) = ( Q(:,IRHOE) + p ) * w
   
         end associate

      end function InviscidFlux1D

      pure function InviscidFlux2D( N , Q ) result ( F )
         implicit none
         integer,       intent (in) :: N
         real(kind=RP), intent (in) :: Q(0:N , 0:N , 1:NCONS)
         real(kind=RP)              :: F(0:N , 0:N , 1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)           :: u(0:N,0:N) , v(0:N,0:N) , w(0:N,0:N) , p(0:N,0:N)

         associate ( gammaMinus1 => thermodynamics % gammaMinus1 ) 

         u = Q(:,:,IRHOU) / Q(:,:,IRHO)
         v = Q(:,:,IRHOV) / Q(:,:,IRHO)
         w = Q(:,:,IRHOW) / Q(:,:,IRHO)
         p = gammaMinus1 * (Q(:,:,IRHOE) - 0.5_RP * ( Q(:,:,IRHOU) * u + Q(:,:,IRHOV) * v + Q(:,:,IRHOW) * w ) )
         
         F(:,:,IRHO , IX ) = Q(:,:,IRHOU)
         F(:,:,IRHOU, IX ) = Q(:,:,IRHOU) * u + p
         F(:,:,IRHOV, IX ) = Q(:,:,IRHOU) * v
         F(:,:,IRHOW, IX ) = Q(:,:,IRHOU) * w
         F(:,:,IRHOE, IX ) = ( Q(:,:,IRHOE) + p ) * u

         F(:,:,IRHO , IY ) = Q(:,:,IRHOV)
         F(:,:,IRHOU ,IY ) = F(:,:,IRHOV,IX)
         F(:,:,IRHOV ,IY ) = Q(:,:,IRHOV) * v + p
         F(:,:,IRHOW ,IY ) = Q(:,:,IRHOV) * w
         F(:,:,IRHOE ,IY ) = ( Q(:,:,IRHOE) + p ) * v

         F(:,:,IRHO ,IZ) = Q(:,:,IRHOW)
         F(:,:,IRHOU,IZ) = F(:,:,IRHOW,IX)
         F(:,:,IRHOV,IZ) = F(:,:,IRHOW,IY)
         F(:,:,IRHOW,IZ) = Q(:,:,IRHOW) * w + P
         F(:,:,IRHOE,IZ) = ( Q(:,:,IRHOE) + p ) * w

         end associate

      end function InviscidFlux2D

      pure function InviscidFlux3D( Nx, Ny, Nz, Q ) result ( F )
         implicit none
         integer,       intent (in) :: Nx, Ny, Nz
         real(kind=RP), intent (in) :: Q(0:Nx,0:Ny,0:Nz , 1:NCONS)
         real(kind=RP)              :: F(0:Nx,0:Ny,0:Nz , 1:NCONS , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)           :: u(0:Nx,0:Ny,0:Nz) , v(0:Nx,0:Ny,0:Nz) , w(0:Nx,0:Ny,0:Nz) , p(0:Nx,0:Ny,0:Nz)

         associate ( gammaMinus1 => thermodynamics % gammaMinus1 ) 

         u = Q(:,:,:,IRHOU) / Q(:,:,:,IRHO)
         v = Q(:,:,:,IRHOV) / Q(:,:,:,IRHO)
         w = Q(:,:,:,IRHOW) / Q(:,:,:,IRHO)
         p = gammaMinus1 * (Q(:,:,:,IRHOE) - 0.5_RP * ( Q(:,:,:,IRHOU) * u + Q(:,:,:,IRHOV) * v + Q(:,:,:,IRHOW) * w ) )
         
         F(:,:,:,IRHO , IX ) = Q(:,:,:,IRHOU)
         F(:,:,:,IRHOU, IX ) = Q(:,:,:,IRHOU) * u + p
         F(:,:,:,IRHOV, IX ) = Q(:,:,:,IRHOU) * v
         F(:,:,:,IRHOW, IX ) = Q(:,:,:,IRHOU) * w
         F(:,:,:,IRHOE, IX ) = ( Q(:,:,:,IRHOE) + p ) * u

         F(:,:,:,IRHO , IY ) = Q(:,:,:,IRHOV)
         F(:,:,:,IRHOU ,IY ) = F(:,:,:,IRHOV,IX)
         F(:,:,:,IRHOV ,IY ) = Q(:,:,:,IRHOV) * v + p
         F(:,:,:,IRHOW ,IY ) = Q(:,:,:,IRHOV) * w
         F(:,:,:,IRHOE ,IY ) = ( Q(:,:,:,IRHOE) + p ) * v

         F(:,:,:,IRHO ,IZ) = Q(:,:,:,IRHOW)
         F(:,:,:,IRHOU,IZ) = F(:,:,:,IRHOW,IX)
         F(:,:,:,IRHOV,IZ) = F(:,:,:,IRHOW,IY)
         F(:,:,:,IRHOW,IZ) = Q(:,:,:,IRHOW) * w + P
         F(:,:,:,IRHOE,IZ) = ( Q(:,:,:,IRHOE) + p ) * w

         end associate

      end function InviscidFlux3D
!
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! DiffusionRiemannSolution computes the coupling on the solution for
!! the calculation of the gradient terms.
!---------------------------------------------------------------------
!
      SUBROUTINE DiffusionRiemannSolution( nHat, QLeft, QRight, Q )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(N_EQN) :: Qleft, Qright, Q
         REAL(KIND=RP), DIMENSION(3)     :: nHat
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER :: j
!
!        -----------------------------------------------
!        For now, this is simply the Bassi/Rebay average
!        -----------------------------------------------
!
         DO j = 1, N_EQN
            Q(j) = 0.5_RP*(Qleft(j) + Qright(j))
         END DO

      END SUBROUTINE DiffusionRiemannSolution
!
! /////////////////////////////////////////////////////////////////////////////
!
!-----------------------------------------------------------------------------
!! DiffusionRiemannSolution computes the coupling on the gradients for
!! the calculation of the contravariant diffusive flux.
!-----------------------------------------------------------------------------
!
      SUBROUTINE DiffusionRiemannFlux(nHat, ds, Q, gradLeft, gradRight, flux)
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(N_EQN)   :: Q,flux
      REAL(KIND=RP), DIMENSION(3,N_EQN) :: gradLeft, gradRight
      REAL(KIND=RP), DIMENSION(3)       :: nHat
      REAL(KIND=RP)                     :: ds
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                           :: j,k
      REAL(KIND=RP), DIMENSION(3,N_EQN) :: grad
      REAL(KIND=RP), DIMENSION(N_EQN)   :: fx, fy, fz
!
!     -------------------------------------------------
!     For now, this simply uses the Bassi/Rebay average
!     -------------------------------------------------
!
      DO j = 1, N_GRAD_EQN
         DO k = 1,3
            grad(k,j) = 0.5_RP*(gradLeft(k,j) + gradRight(k,j))
         END DO
      END DO
!
!     ----------------------------
!     Compute the component fluxes
!     ----------------------------
!
      CALL xDiffusiveFlux( Q, grad, fx )
      CALL yDiffusiveFlux( Q, grad, fy )
      CALL zDiffusiveFlux( Q, grad, fz )
!
!     ------------------------------
!     Compute the contravariant flux
!     ------------------------------
!
      DO j = 1, N_EQN
         flux(j) = ds*(nHat(1)*fx(j) + nHat(2)*fy(j) + nHat(3)*fz(j))
      END DO
      
      END SUBROUTINE DiffusionRiemannFlux
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! xDiffusiveFlux computes the x viscous flux component.
!---------------------------------------------------------------------
!
      SUBROUTINE xDiffusiveFlux( Q, grad, f )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
!!    Q contains the solution values
!
      REAL(KIND=RP), DIMENSION(N_EQN)      :: Q
!
!!    grad contains the (physical) gradients needed for the
!!    equations. For the Navier-Stokes equations these are
!!    grad(u), grad(v), grad(w), grad(T).
!
      REAL(KIND=RP), DIMENSION(3,N_GRAD_EQN) :: grad
!
!!     f is the viscous flux in the physical x direction returned by
!!     this routine.
! 
      REAL(KIND=RP), DIMENSION(N_EQN)      :: f
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: tauXX, tauXY, tauXZ
      REAL(KIND=RP) :: T, muOfT, kappaOfT, divVelocity
      REAL(KIND=RP) :: u, v, w
!      
      associate ( Re => dimensionless % Re , &
                  Pr => dimensionless % Pr , &
                  gammaM2 => dimensionless % gammaM2, &
                  gammaDivGammaMinus1 => thermodynamics % gammaDivGammaMinus1 ) 

      T        = Temperature(Q)
      muOfT    = MolecularDiffusivity(T)
      kappaOfT = ThermalDiffusivity(T)
      u        = Q(2)/Q(1)
      v        = Q(3)/Q(1)
      w        = Q(4)/Q(1)
      
      divVelocity = grad(1,1) + grad(2,2) + grad(3,3)
      tauXX       = 2.0_RP*muOfT*(grad(1,1) - divVelocity/3._RP)
      tauXY       = muOfT*(grad(1,2) + grad(2,1))
      tauXZ       = muOfT*(grad(1,3) + grad(3,1))
      
      f(1) = 0.0_RP
      f(2) = tauXX/RE
      f(3) = tauXY/RE
      f(4) = tauXZ/RE
      f(5) = (u*tauXX + v*tauXY + w*tauXZ + &
     &        gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*grad(1,4))/RE

      end associate

      END SUBROUTINE xDiffusiveFlux
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! yDiffusiveFlux computes the y viscous flux component.
!---------------------------------------------------------------------
!
      SUBROUTINE yDiffusiveFlux( Q, grad, f )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
!!    Q contains the solution values
!
      REAL(KIND=RP), DIMENSION(N_EQN)      :: Q
!
!!    grad contains the (physical) gradients needed for the
!!    equations. For the Navier-Stokes equations these are
!!    grad(u), grad(v), grad(w), grad(T).
!
      REAL(KIND=RP), DIMENSION(3,N_GRAD_EQN) :: grad
!
!!     f is the viscous flux in the physical x direction returned by
!!     this routine.
! 
      REAL(KIND=RP), DIMENSION(N_EQN)      :: f
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: tauYX, tauYY, tauYZ
      REAL(KIND=RP) :: T, muOfT, kappaOfT, divVelocity
      REAL(KIND=RP) :: u, v, w
!      

      associate ( Re => dimensionless % Re , &
                  Pr => dimensionless % Pr , &
                  gammaM2 => dimensionless % gammaM2, &
                  gammaDivGammaMinus1 => thermodynamics % gammaDivGammaMinus1 ) 

      T        = Temperature(Q)
      muOfT    = MolecularDiffusivity(T)
      kappaOfT = ThermalDiffusivity(T)
      u        = Q(2)/Q(1)
      v        = Q(3)/Q(1)
      w        = Q(4)/Q(1)
      
      divVelocity = grad(1,1) + grad(2,2) + grad(3,3)
      tauYX       = muOfT*(grad(1,2) + grad(2,1))
      tauYY       = 2.0_RP*muOfT*(grad(2,2) - divVelocity/3._RP)
      tauYZ       = muOfT*(grad(2,3) + grad(3,2))
      
      f(1) = 0.0_RP
      f(2) = tauYX/RE
      f(3) = tauYY/RE
      f(4) = tauYZ/RE
      f(5) = (u*tauYX + v*tauYY + w*tauYZ + &
     &        gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*grad(2,4))/RE

      end associate

      END SUBROUTINE yDiffusiveFlux
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! yDiffusiveFlux computes the y viscous flux component.
!---------------------------------------------------------------------
!
      SUBROUTINE zDiffusiveFlux( Q, grad, f )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
!!    Q contains the solution values
!
      REAL(KIND=RP), DIMENSION(N_EQN)      :: Q
!
!!    grad contains the (physical) gradients needed for the
!!    equations. For the Navier-Stokes equations these are
!!    grad(u), grad(v), grad(w), grad(T).
!
      REAL(KIND=RP), DIMENSION(3,N_GRAD_EQN) :: grad
!
!!     f is the viscous flux in the physical x direction returned by
!!     this routine.
! 
      REAL(KIND=RP), DIMENSION(N_EQN)      :: f
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP)           :: tauZX, tauZY, tauZZ
      REAL(KIND=RP)           :: T, muOfT, kappaOfT, divVelocity
      REAL(KIND=RP)           :: u, v, w
!      

      associate ( Re => dimensionless % Re , &
                  Pr => dimensionless % Pr , &
                  gammaM2 => dimensionless % gammaM2, &
                  gammaDivGammaMinus1 => thermodynamics % gammaDivGammaMinus1 ) 

      T        = Temperature(Q)
      muOfT    = MolecularDiffusivity(T)
      kappaOfT = ThermalDiffusivity(T)
      u        = Q(2)/Q(1)
      v        = Q(3)/Q(1)
      w        = Q(4)/Q(1)
      
      divVelocity = grad(1,1) + grad(2,2) + grad(3,3)
      tauZX       = muOfT*(grad(1,3) + grad(3,1))
      tauZY       = muOfT*(grad(2,3) + grad(3,2))
      tauZZ       = 2.0_RP*muOfT*(grad(3,3) - divVelocity/3._RP)
      
      f(1) = 0.0_RP
      f(2) = tauZX/RE
      f(3) = tauZY/RE
      f(4) = tauZZ/RE
      f(5) = (u*tauZX + v*tauZY + w*tauZZ + &
     &        gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*grad(3,4))/RE

      end associate

      END SUBROUTINE zDiffusiveFlux

      pure function ViscousFlux0D( Q , U_x , U_y , U_z ) result (F)
         implicit none
         real ( kind=RP ) , intent ( in ) :: Q    ( 1:NCONS          ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 1:N_GRAD_EQN     ) 
         real(kind=RP)                    :: F    ( 1:NCONS , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                    :: T , muOfT , kappaOfT
         real(kind=RP)                    :: divV
         real(kind=RP)                    :: u , v , w

         associate ( Re => dimensionless % Re , &
                     Pr => dimensionless % Pr , &
                     gammaM2 => dimensionless % gammaM2, &
                     gammaDivGammaMinus1 => thermodynamics % gammaDivGammaMinus1 ) 

         u = Q(IRHOU) / Q(IRHO)
         v = Q(IRHOV) / Q(IRHO)
         w = Q(IRHOW) / Q(IRHO)

         T     = Temperature(Q)
         muOfT = MolecularDiffusivity(T)
         kappaOfT = ThermalDiffusivity(T)

         divV = U_x(IGU) + U_y(IGV) + U_z(IGW)

         F(IRHO,IX)  = 0.0_RP
         F(IRHOU,IX) = muOfT * (2.0_RP * U_x(IGU) - 2.0_RP/3.0_RP * divV ) / RE
         F(IRHOV,IX) = muOfT * ( U_x(IGV) + U_y(IGU) ) / RE
         F(IRHOW,IX) = muOfT * ( U_x(IGW) + U_z(IGU) ) / RE
         F(IRHOE,IX) = F(IRHOU,IX) * u + F(IRHOV,IX) * v + F(IRHOW,IX) * w + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_x(IGT) / RE

         F(IRHO,IY) = 0.0_RP
         F(IRHOU,IY) = F(IRHOV,IX)
         F(IRHOV,IY) = muOfT * (2.0_RP * U_y(IGV) - 2.0_RP / 3.0_RP * divV ) / RE
         F(IRHOW,IY) = muOfT * ( U_y(IGW) + U_z(IGV) ) / RE
         F(IRHOE,IY) = F(IRHOU,IY) * u + F(IRHOV,IY) * v + F(IRHOW,IY) * w + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_y(IGT) / RE

         F(IRHO,IZ) = 0.0_RP
         F(IRHOU,IZ) = F(IRHOW,IX)
         F(IRHOV,IZ) = F(IRHOW,IY)
         F(IRHOW,IZ) = muOfT * ( 2.0_RP * U_z(IGW) - 2.0_RP / 3.0_RP * divV ) / RE
         F(IRHOE,IZ) = F(IRHOU,IZ) * u + F(IRHOV,IZ) * v + F(IRHOW,IZ) * w + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_z(IGT) / RE

         end associate

      end function ViscousFlux0D

      pure function ViscousFlux1D( N , Q , U_x , U_y , U_z ) result (F)
         implicit none
         integer          , intent ( in ) :: N
         real ( kind=RP ) , intent ( in ) :: Q    ( 0:N , 1:NCONS          ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 0:N , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 0:N , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 0:N , 1:N_GRAD_EQN     ) 
         real(kind=RP)                    :: F    ( 0:N , 1:NCONS , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                    :: T(0:N) , muOfT(0:N) , kappaOfT(0:N)
         real(kind=RP)                    :: divV(0:N)
         real(kind=RP)                    :: u(0:N) , v(0:N) , w(0:N)
         integer                          :: i

         associate ( Re => dimensionless % Re , &
                     Pr => dimensionless % Pr , &
                     gammaM2 => dimensionless % gammaM2, &
                     gammaMinus1 => thermodynamics % gammaMinus1, &
                     gammaDivGammaMinus1 => thermodynamics % gammaDivGammaMinus1 ) 

         u = Q(:,IRHOU) / Q(:,IRHO)
         v = Q(:,IRHOV) / Q(:,IRHO)
         w = Q(:,IRHOW) / Q(:,IRHO)


         T = gammaM2 * gammaMinus1 * (Q(:,IRHOE)  & 
               - 0.5_RP * ( Q(:,IRHOU) * u + Q(:,IRHOV) * v + Q(:,IRHOW) * w ) ) / Q(:,IRHO)


         do i = 0 , N
            muOfT(i) = MolecularDiffusivity(T(i))
            kappaOfT(i) = ThermalDiffusivity(T(i))
         end do 

         divV = U_x(:,IGU) + U_y(:,IGV) + U_z(:,IGW)

         F(:,IRHO ,IX) = 0.0_RP
         F(:,IRHOU,IX) = muOfT * (2.0_RP * U_x(:,IGU) - 2.0_RP/3.0_RP * divV ) / RE
         F(:,IRHOV,IX) = muOfT * ( U_x(:,IGV) + U_y(:,IGU) ) / RE
         F(:,IRHOW,IX) = muOfT * ( U_x(:,IGW) + U_z(:,IGU) ) / RE
         F(:,IRHOE,IX) = F(:,IRHOU,IX) * u + F(:,IRHOV,IX) * v + F(:,IRHOW,IX) * w &
               + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_x(:,IGT) / RE

         F(:,IRHO ,IY) = 0.0_RP
         F(:,IRHOU,IY) = F(:,IRHOV,IX)
         F(:,IRHOV,IY) = muOfT * (2.0_RP * U_y(:,IGV) - 2.0_RP / 3.0_RP * divV ) / RE
         F(:,IRHOW,IY) = muOfT * ( U_y(:,IGW) + U_z(:,IGV) ) / RE
         F(:,IRHOE,IY) = F(:,IRHOU,IY) * u + F(:,IRHOV,IY) * v + F(:,IRHOW,IY) * w &
               + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_y(:,IGT) / RE

         F(:,IRHO,IZ ) = 0.0_RP
         F(:,IRHOU,IZ) = F(:,IRHOW,IX)
         F(:,IRHOV,IZ) = F(:,IRHOW,IY)
         F(:,IRHOW,IZ) = muOfT * ( 2.0_RP * U_z(:,IGW) - 2.0_RP / 3.0_RP * divV ) / RE
         F(:,IRHOE,IZ) = F(:,IRHOU,IZ) * u + F(:,IRHOV,IZ) * v + F(:,IRHOW,IZ) * w &
               + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_z(:,IGT) / RE

         end associate

      end function ViscousFlux1D

      pure function ViscousFlux2D( N , Q , U_x , U_y , U_z ) result (F)
         implicit none
         integer          , intent ( in ) :: N
         real ( kind=RP ) , intent ( in ) :: Q    ( 0:N , 0:N , 1:NCONS          ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 0:N , 0:N , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 0:N , 0:N , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 0:N , 0:N , 1:N_GRAD_EQN     ) 
         real(kind=RP)                    :: F    ( 0:N , 0:N , 1:NCONS , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)                    :: T(0:N,0:N) , muOfT(0:N,0:N) , kappaOfT(0:N,0:N)
         real(kind=RP)                    :: divV(0:N,0:N)
         real(kind=RP)                    :: u(0:N,0:N) , v(0:N,0:N) , w(0:N,0:N)
         integer                          :: i , j 

         associate ( Re => dimensionless % Re , &
                     Pr => dimensionless % Pr , &
                     gammaM2 => dimensionless % gammaM2, &
                     gammaDivGammaMinus1 => thermodynamics % gammaDivGammaMinus1 ) 

         u = Q(:,:,IRHOU) / Q(:,:,IRHO)
         v = Q(:,:,IRHOV) / Q(:,:,IRHO)
         w = Q(:,:,IRHOW) / Q(:,:,IRHO)


         T = gammaM2 * (Q(:,:,IRHOE)  & 
               - 0.5_RP * ( Q(:,:,IRHOU) * u + Q(:,:,IRHOV) * v + Q(:,:,IRHOW) * w ) ) / Q(:,:,IRHO)

         do i = 0 , N ;    do j = 0 , N 
            muOfT    ( i,j )  = MolecularDiffusivity ( T ( i,j )  ) 
            kappaOfT ( i,j )  = ThermalDiffusivity   ( T ( i,j )  ) 
         end do       ;    end do 

         divV = U_x(:,:,IGU) + U_y(:,:,IGV) + U_z(:,:,IGW)

         F(:,:,IRHO ,IX) = 0.0_RP
         F(:,:,IRHOU,IX) = muOfT * (2.0_RP * U_x(:,:,IGU) - 2.0_RP/3.0_RP * divV ) / RE
         F(:,:,IRHOV,IX) = muOfT * ( U_x(:,:,IGV) + U_y(:,:,IGU) ) / RE
         F(:,:,IRHOW,IX) = muOfT * ( U_x(:,:,IGW) + U_z(:,:,IGU) ) / RE
         F(:,:,IRHOE,IX) = F(:,:,IRHOU,IX) * u + F(:,:,IRHOV,IX) * v + F(:,:,IRHOW,IX) * w &
               + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_x(:,:,IGT) / RE

         F(:,:,IRHO ,IY) = 0.0_RP
         F(:,:,IRHOU,IY) = F(:,:,IRHOV,IX)
         F(:,:,IRHOV,IY) = muOfT * (2.0_RP * U_y(:,:,IGV) - 2.0_RP / 3.0_RP * divV ) / RE
         F(:,:,IRHOW,IY) = muOfT * ( U_y(:,:,IGW) + U_z(:,:,IGV) ) / RE
         F(:,:,IRHOE,IY) = F(:,:,IRHOU,IY) * u + F(:,:,IRHOV,IY) * v + F(:,:,IRHOW,IY) * w &
               + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_y(:,:,IGT) / RE

         F(:,:,IRHO,IZ ) = 0.0_RP
         F(:,:,IRHOU,IZ) = F(:,:,IRHOW,IX)
         F(:,:,IRHOV,IZ) = F(:,:,IRHOW,IY)
         F(:,:,IRHOW,IZ) = muOfT * ( 2.0_RP * U_z(:,:,IGW) - 2.0_RP / 3.0_RP * divV ) / RE
         F(:,:,IRHOE,IZ) = F(:,:,IRHOU,IZ) * u + F(:,:,IRHOV,IZ) * v + F(:,:,IRHOW,IZ) * w &
               + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_z(:,:,IGT) / RE

         end associate

      end function ViscousFlux2D

      pure function ViscousFlux3D( Nx, Ny, Nz , Q , U_x , U_y , U_z ) result (F)
         implicit none
         integer          , intent ( in ) :: Nx
         integer          , intent ( in ) :: Ny
         integer          , intent ( in ) :: Nz
         real ( kind=RP ) , intent ( in ) :: Q    ( 0:Nx , 0:Ny , 0:Nz , 1:NCONS          ) 
         real ( kind=RP ) , intent ( in ) :: U_x  ( 0:Nx , 0:Ny , 0:Nz , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_y  ( 0:Nx , 0:Ny , 0:Nz , 1:N_GRAD_EQN     ) 
         real ( kind=RP ) , intent ( in ) :: U_z  ( 0:Nx , 0:Ny , 0:Nz , 1:N_GRAD_EQN     ) 
         real ( kind=RP )                 :: F    ( 0:Nx , 0:Ny , 0:Nz , 1:NCONS , 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: T(0:Nx,0:Ny,0:Nz) , muOfT(0:Nx,0:Ny,0:Nz) , kappaOfT(0:Nx,0:Ny,0:Nz)
         real(kind=RP) :: divV(0:Nx,0:Ny,0:Nz)
         real(kind=RP) :: u(0:Nx,0:Ny,0:Nz) , v(0:Nx,0:Ny,0:Nz) , w(0:Nx,0:Ny,0:Nz)
         integer       :: i , j , k

         associate ( Re => dimensionless % Re , &
                     Pr => dimensionless % Pr , &
                     gammaM2 => dimensionless % gammaM2, &
                     gammaMinus1 => thermodynamics % gammaMinus1, &
                     gammaDivGammaMinus1 => thermodynamics % gammaDivGammaMinus1 ) 

         u = Q(:,:,:,IRHOU) / Q(:,:,:,IRHO)
         v = Q(:,:,:,IRHOV) / Q(:,:,:,IRHO)
         w = Q(:,:,:,IRHOW) / Q(:,:,:,IRHO)


         T = gammaM2 * gammaMinus1 * (Q(:,:,:,IRHOE)  & 
               - 0.5_RP * ( Q(:,:,:,IRHOU) * u + Q(:,:,:,IRHOV) * v + Q(:,:,:,IRHOW) * w)) / Q(:,:,:,IRHO)

         do i = 0 , Nx ;    do j = 0 , Ny ;    do k = 0 , Nz
            muOfT    ( i,j,k )  = MolecularDiffusivity ( T ( i,j,k )  ) 
            kappaOfT ( i,j,k )  = ThermalDiffusivity   ( T ( i,j,k )  ) 
         end do       ;    end do       ;    end do 

         divV = U_x(:,:,:,IGU) + U_y(:,:,:,IGV) + U_z(:,:,:,IGW)

         F(:,:,:,IRHO ,IX) = 0.0_RP
         F(:,:,:,IRHOU,IX) = muOfT * (2.0_RP * U_x(:,:,:,IGU) - 2.0_RP/3.0_RP * divV ) / RE
         F(:,:,:,IRHOV,IX) = muOfT * ( U_x(:,:,:,IGV) + U_y(:,:,:,IGU) ) / RE
         F(:,:,:,IRHOW,IX) = muOfT * ( U_x(:,:,:,IGW) + U_z(:,:,:,IGU) ) / RE
         F(:,:,:,IRHOE,IX) = F(:,:,:,IRHOU,IX) * u + F(:,:,:,IRHOV,IX) * v + F(:,:,:,IRHOW,IX) * w &
               + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_x(:,:,:,IGT) / RE

         F(:,:,:,IRHO ,IY) = 0.0_RP
         F(:,:,:,IRHOU,IY) = F(:,:,:,IRHOV,IX)
         F(:,:,:,IRHOV,IY) = muOfT * (2.0_RP * U_y(:,:,:,IGV) - 2.0_RP / 3.0_RP * divV ) / RE
         F(:,:,:,IRHOW,IY) = muOfT * ( U_y(:,:,:,IGW) + U_z(:,:,:,IGV) ) / RE
         F(:,:,:,IRHOE,IY) = F(:,:,:,IRHOU,IY) * u + F(:,:,:,IRHOV,IY) * v + F(:,:,:,IRHOW,IY) * w &
               + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_y(:,:,:,IGT) / RE

         F(:,:,:,IRHO,IZ ) = 0.0_RP
         F(:,:,:,IRHOU,IZ) = F(:,:,:,IRHOW,IX)
         F(:,:,:,IRHOV,IZ) = F(:,:,:,IRHOW,IY)
         F(:,:,:,IRHOW,IZ) = muOfT * ( 2.0_RP * U_z(:,:,:,IGW) - 2.0_RP / 3.0_RP * divV ) / RE
         F(:,:,:,IRHOE,IZ) = F(:,:,:,IRHOU,IZ) * u + F(:,:,:,IRHOV,IZ) * v + F(:,:,:,IRHOW,IZ) * w &
               + gammaDivGammaMinus1*kappaOfT/(PR*gammaM2)*U_z(:,:,:,IGT) / RE

         end associate

      end function ViscousFlux3D
!
!
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      SUBROUTINE GradientValuesForQ_0D( Q, U )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(N_EQN)     , INTENT(IN)  :: Q
      REAL(KIND=RP), DIMENSION(N_GRAD_EQN), INTENT(OUT) :: U
!
!     ---------------
!     Local Variables
!     ---------------
!     
      U(1) = Q(2)/Q(1)
      U(2) = Q(3)/Q(1)
      U(3) = Q(4)/Q(1)
      U(4) = Temperature(Q)

      END SUBROUTINE GradientValuesForQ_0D

      SUBROUTINE GradientValuesForQ_3D( Q, U )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), INTENT(IN)  :: Q(0:,0:,0:,:)
      REAL(KIND=RP), INTENT(OUT) :: U(0:,0:,0:,:)
      integer                    :: Nx, Ny, Nz
!
!     ---------------
!     Local Variables
!     ---------------
!     
      associate ( gammaM2 => dimensionless % gammaM2, &
                  gammaMinus1 => thermodynamics % gammaMinus1 ) 
      Nx = size(Q , 1) - 1
      Ny = size(Q , 2) - 1
      Nz = size(Q , 3) - 1
      
      U(0:Nx,0:Ny,0:Nz,IGU) = Q(0:Nx,0:Ny,0:Nz,IRHOU) / Q(0:Nx,0:Ny,0:Nz,IRHO) 
      U(0:Nx,0:Ny,0:Nz,IGV) = Q(0:Nx,0:Ny,0:Nz,IRHOV) / Q(0:Nx,0:Ny,0:Nz,IRHO) 
      U(0:Nx,0:Ny,0:Nz,IGW) = Q(0:Nx,0:Ny,0:Nz,IRHOW) / Q(0:Nx,0:Ny,0:Nz,IRHO) 
      U(0:Nx,0:Ny,0:Nz,IGT) = gammaM2 * gammaMinus1 * ( Q(0:Nx,0:Ny,0:Nz,IRHOE) / Q(0:Nx,0:Ny,0:Nz,IRHO) &
                  - 0.5_RP * ( U(0:Nx,0:Ny,0:Nz,IGU) * U(0:Nx,0:Ny,0:Nz,IGU) &
                             + U(0:Nx,0:Ny,0:Nz,IGV) * U(0:Nx,0:Ny,0:Nz,IGV) &
                             + U(0:Nx,0:Ny,0:Nz,IGW) * U(0:Nx,0:Ny,0:Nz,IGW) ) )

      end associate

      END SUBROUTINE GradientValuesForQ_3D
!
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! Compute the pressure from the state variables
!---------------------------------------------------------------------
!
      PURE FUNCTION Pressure(Q) RESULT(P)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(N_EQN), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: P
      
      P = thermodynamics % gammaMinus1*(Q(5) - 0.5_RP*(Q(2)**2 + Q(3)**2 + Q(4)**2)/Q(1))

      END FUNCTION Pressure
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the molecular diffusivity by way of Sutherland's law
!---------------------------------------------------------------------
!
      PURE FUNCTION MolecularDiffusivity(T) RESULT(mu)
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


      END FUNCTION MolecularDiffusivity
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the thermal diffusivity by way of Sutherland's law
!---------------------------------------------------------------------
!
      PURE FUNCTION ThermalDiffusivity(T) RESULT(kappa)
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
      REAL(KIND=RP) :: kappa !! The diffusivity
!      
      kappa = (1._RP + tRatio)/(T + tRatio)*T*SQRT(T)


      END FUNCTION ThermalDiffusivity
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the temperature from the state variables
!---------------------------------------------------------------------
!
      PURE FUNCTION Temperature(Q) RESULT(T)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(N_EQN), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: T
!
      T = dimensionless % gammaM2*Pressure(Q)/Q(1)

      END FUNCTION Temperature
      
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
      USE Physics, ONLY:Pressure
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
