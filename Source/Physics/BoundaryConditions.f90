!
!////////////////////////////////////////////////////////////////////////
!
!      BoundaryConditions.f90
!      Created: 2011-07-21 10:30:14 -0400 
!      By: David Kopriva
!      
!      This file has the BCs for Navier-Stokes computations.
!      The two entry points are:
!
!         SUBROUTINE ExternalStateForBoundaryName( x, y, t, nHat, Q, boundaryName )
!
!            Returns the state Q for the boundary with name boundaryName
!
!         ExternalGradientForBoundaryName( x, y, t, nHat, gradU, boundaryName )
!
!            Returns the gradient gradU for the boundary named boundaryName
!
!      The two entry points call boundary routines depending on the type 
!      of boundary conditions. Currently supported BCs are:
!
!      (1) Free-Stream Inflow-Outflow (constant external state)
!          Dirichlet: Free stream value
!          Neumann  : All gradients = 0
!
!          Set up in:
!
!             SUBROUTINE UniformFlowState( x, t, Q )
!             SUBROUTINE UniformFlowNeumann( x, t, nHat, gradU )
!
!      (2) Free-Slip Wall (symmetry boundary, normal velocity and normal gradients vanish)
!
!             Dirichlet: \vec q^ext \cdot \hat n = - \vec q^int \cdot \hat n
!             Neumann  : \vec \nabla q^ext \cdot \hat n = - \vec \nabla q^int \cdot \hat n
!
!          Set up in:
!
!             SUBROUTINE FreeSlipWallState ( x, t, nHat, Q )
!             SUBROUTINE FreeSlipNeumann ( x, t, nHat, grad )
!
!      (3) Adiabataic Wall (velocity = 0, normal temperature gradient = 0)
!
!             Dirichlet: \vec q^ext = - \vec q^int
!             Neumann  : \vec \nabla q(3)^ext \cdot \hat n = - \vec \nabla q(3)^int \cdot \hat n
!
!          Set up in:
!
!             SUBROUTINE NoSlipAdiabaticWallState ( x, t, Q )
!             SUBROUTINE NoSlipAdiabaticNeumann( x, t, nHat, gradU )
!
!      (4) Isothermal Wall (velocity = 0, temperature set to wall value, gradient from interior)
!
!          Set up in:
!
!             SUBROUTINE NoSlipIsothermalWallState( x, t, Q )
!             SUBROUTINE NoSlipIsothermalWallNeumann( x, t, gradU )
!
!
!////////////////////////////////////////////////////////////////////////
!
      MODULE BoundaryConditionFunctions
         USE SMConstants
         USE Physics
         USE SharedBCModule
      
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      FUNCTION WallTemperatureForBoundaryNamed(boundaryName) RESULT(T)
         IMPLICIT NONE
         CHARACTER(LEN=32), INTENT(IN)    :: boundaryName
         REAL(KIND=RP)                    :: T
         
         ! Choose different temperatures according to boundary name, if desired...
         ! For now, just use an input value.
         
         T = bcValueDictionary % realValueForKey(boundaryName)
      
      END FUNCTION WallTemperatureForBoundaryNamed
!
!////////////////////////////////////////////////////////////////////////
!
      FUNCTION ExternalPressureForBoundaryNamed(boundaryName) RESULT(p)
         IMPLICIT NONE
         CHARACTER(LEN=32), INTENT(IN)    :: boundaryName
         REAL(KIND=RP)                    :: p
         
         ! Choose different temperatures according to boundary name, if desired...
         ! For now, just use the UninformFlow value.
         
         p    = 1.0_RP/(gammaM2)
      
      END FUNCTION ExternalPressureForBoundaryNamed
!
!////////////////////////////////////////////////////////////////////////
!
      FUNCTION WallTemperature(x) RESULT(T)
         IMPLICIT NONE
         REAL(KIND=RP) :: x(3)
         REAL(KIND=RP) :: T
         
         ! just a simple set the temperature value functionn
         T = 1.0_RP
      
      END FUNCTION WallTemperature
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE FreeSlipWallState( x, t, nHat, Q )
!
!     ----------------------------------------------------
!     Set up the conditions for a wall boundary condition.
!     This is the no-slip condition with normal velocity 
!     of the exterior fictitious cell set to the negative
!     of the interior value. Used for Euler and Symmetry
!     BCs for Navier-Stokes.
!     ----------------------------------------------------
!
      IMPLICIT NONE 
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), INTENT(IN)    :: x(3), t
      REAL(KIND=RP), INTENT(IN)    :: nHat(3)
      REAL(KIND=RP), INTENT(INOUT) :: Q(N_EQN)
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: qNorm, qTanx, qTany, qTanz
      REAL(KIND=RP) :: QInt(N_EQN)
      
      QInt = Q
      
      qNorm =  nHat(1)*QInt(2) + nHat(2)*QInt(3) + nHat(3)*QInt(4)
      qTanx = QInt(2) - qNorm*nHat(1)
      qTany = QInt(3) - qNorm*nHat(2)
      qTanz = QInt(4) - qNorm*nHat(3)

      Q(1) = QInt(1)
      Q(2) = qTanx - qNorm*nHat(1)
      Q(3) = qTany - qNorm*nHat(2)
      Q(4) = qTanz - qNorm*nHat(3)
      Q(5) = QInt(5)

      END SUBROUTINE FreeSlipWallState
!
!     /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!
!!    This routine sets the Neumann conditions for a free slip wall.
!!    In this case it means setting the normal gradients on the exterior
!!    (right) to be the negative of the interior (Left) gradients.
!!    so that they average out to zero.
!     ----------------------------------------------------------------
!
      SUBROUTINE FreeSlipNeumann(x, t, nHat, U_x, U_y, U_z )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), INTENT(IN)    :: x(3), t
      REAL(KIND=RP), INTENT(IN)    :: nHat(3)
      REAL(KIND=RP), INTENT(INOUT) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: gradUNorm, UTanx, UTany, UTanz
      INTEGER       :: k
!
      DO k = 1, N_GRAD_EQN 
         gradUNorm =  nHat(1)*U_x(k) + nHat(2)*U_y(k) + nHat(3)*U_z(k)
         UTanx = U_x(k) - gradUNorm*nHat(1)
         UTany = U_y(k) - gradUNorm*nHat(2)
         UTanz = U_z(k) - gradUNorm*nHat(3)
   
         U_x(k) = UTanx - gradUNorm*nHat(1)
         U_y(k) = UTany - gradUNorm*nHat(2)
         U_z(k) = UTanz - gradUNorm*nHat(3)
      END DO

      END SUBROUTINE FreeSlipNeumann
!
!---------------------------------------------------------------------
!!   SUBROUTINE NoSlipWall( thisMortar, Q, time ): Enforce
!!   the no velocity condition through the Riemann solver.
!---------------------------------------------------------------------
!
      SUBROUTINE NoSlipAdiabaticWallState( x, t, Q )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), INTENT(IN)    :: x(3), t
         REAL(KIND=RP), INTENT(INOUT) :: Q(N_EQN)
!
!        -----------------------------------------------
!        Generate the external flow along the face, that
!        represents a solid wall.
!        -----------------------------------------------
!
         Q(1) =  Q(1)
         Q(2) = -Q(2)
         Q(3) = -Q(3)
         Q(4) = -Q(4)
         Q(5) =  Q(5)

      END SUBROUTINE NoSlipAdiabaticWallState
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE NoSlipAdiabaticWallNeumann( x, t, nHat, U_x, U_y, u_z )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), INTENT(IN)    :: x(3), t
         REAL(KIND=RP), INTENT(IN)    :: nHat(3)
         REAL(KIND=RP), INTENT(INOUT) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER :: k = 3 ! = temperature
!
         REAL(KIND=RP) :: gradUNorm, UTanx, UTany, UTanz
!
         gradUNorm =  nHat(1)*U_x(k) + nHat(2)*U_y(k)+ nHat(3)*U_z(k)
         UTanx = U_x(k) - gradUNorm*nHat(1)
         UTany = U_y(k) - gradUNorm*nHat(2)
         UTanz = U_z(k) - gradUNorm*nHat(3)
   
         U_x(k) = UTanx - gradUNorm*nHat(1)
         U_y(k) = UTany - gradUNorm*nHat(2)
         U_z(k) = UTanz - gradUNorm*nHat(3)
      
      END SUBROUTINE NoSlipAdiabaticWallNeumann
!
!---------------------------------------------------------------------
!!   SUBROUTINE NoSlipWall( thisMortar, Q, time ): Enforce
!!   the no velocity condition through the Riemann solver.
!---------------------------------------------------------------------
!
      SUBROUTINE NoSlipIsothermalWallState( x, t, Q )
      
      USE SMConstants
      USE PhysicsStorage
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), INTENT(IN)    :: x(3), t
      REAL(KIND=RP), INTENT(INOUT) :: Q(N_EQN)
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP)           :: wallTemp
!
!     -----------------------------------------------
!     Generate the external flow along the face, that
!     represents a solid wall at a fixed temperature.
!     Note that this requires an external routine 
!     called WallTemperature
!     -----------------------------------------------
!
         wallTemp = wallTemperature(x) !Swap this out later with BC version
         Q(1) =  Q(1)
         Q(2) = -Q(2)
         Q(3) = -Q(3)
         Q(4) = -Q(4)
         Q(5) =  Q(1)*wallTemp/gammaMinus1/gammaM2

      END SUBROUTINE NoSlipIsothermalWallState
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE NoSlipIsothermalWallNeumann( x, t, nHat, U_x, U_y, U_z  )
      USE SMConstants
      USE PhysicsStorage
      IMPLICIT NONE
      
      REAL(KIND=RP)                :: x(3), t, nHat(3)
      REAL(KIND=RP), INTENT(INOUT) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)
      
      !do nothing - specifying Temperature is sufficient.
      
      END SUBROUTINE NoSlipIsothermalWallNeumann
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE UniformFlowState( x, t, Q )
      USE SMConstants
      USE PhysicsStorage
      IMPLICIT NONE
      
      REAL(KIND=RP) :: x(3), t
      REAL(KIND=RP) :: Q(N_EQN)
      
      REAL(KIND=RP) :: theta, phi, qq
      REAL(KIND=RP) :: u, v, w, p
      
      theta = AOATheta*(PI/180.0_RP)
      phi   = AOAPhi*(PI/180.0_RP)
      
      qq = 1.0_RP
      u  = qq*cos(theta)*COS(phi)
      v  = qq*sin(theta)*COS(phi)
      w  = qq*SIN(phi)
      
      Q(1) = 1.0_RP
      p    = 1.0_RP/(gammaM2)
      Q(2) = Q(1)*u
      Q(3) = Q(1)*v
      Q(4) = Q(1)*w
      Q(5) = p/(gamma - 1._RP) + 0.5_RP*Q(1)*(u**2 + v**2 + w**2)
      
      END SUBROUTINE UniformFlowState
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE UniformFlowNeumann( x, t, nHat, U_x, U_y, U_z  )
!
!     ----------------------------
!     Set normal gradients to zero
!     ----------------------------
!
      USE SMConstants
      USE PhysicsStorage
      IMPLICIT NONE
      
      REAL(KIND=RP) :: x(3), t
      REAL(KIND=RP) :: nHat(3)
      REAL(KIND=RP), INTENT(INOUT) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)
      
      INTEGER :: k
!
      REAL(KIND=RP) :: gradUNorm, UTanx, UTany, UTanz
!
      DO k = 1, N_GRAD_EQN
         gradUNorm =  nHat(1)*U_x(k) + nHat(2)*U_y(k) + nHat(3)*U_z(k)
         UTanx = U_x(k) - gradUNorm*nHat(1)
         UTany = U_y(k) - gradUNorm*nHat(2)
         UTanz = U_z(k) - gradUNorm*nHat(3)
   
         U_x(k) = UTanx - gradUNorm*nHat(1)
         U_y(k) = UTany - gradUNorm*nHat(2)
         U_z(k) = UTanz - gradUNorm*nHat(3)
      END DO
      
      END SUBROUTINE UniformFlowNeumann
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExternalPressureState( x, t, nHat, Q, pExt )
!
!     -------------------------------------------------------
!     Compute the external state given the internal state and
!     the exterior pressure
!     -------------------------------------------------------
!
      IMPLICIT NONE 
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP) :: x(3), t, nHat(3)
      REAL(KIND=RP) :: Q(N_EQN), pExt
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: qDotN, qTanx, qTany, qTanz, p, a, a2
      REAL(KIND=RP) :: rPlus, entropyConstant, u, v, w, rho, normalMachNo
!      
      qDotN = (nHat(1)*Q(2) + nHat(2)*Q(3) + nHat(3)*Q(4))/Q(1)
      qTanx = Q(2)/Q(1) - qDotN*nHat(1)
      qTany = Q(3)/Q(1) - qDotN*nHat(2)
      qTanz = Q(4)/Q(1) - qDotN*nHat(3)
      
      p            = gammaMinus1*( Q(5) - 0.5_RP*(Q(2)**2 + Q(3)**2 + Q(4)**2)/Q(1) )
      a2           = gamma*p/Q(1)
      a            = SQRT(a2)
      normalMachNo = ABS(qDotN/a)
      
      IF ( normalMachNo <= 1.0_RP )     THEN
!
!        -------------------------------
!        Quantities coming from upstream
!        -------------------------------
!
         rPlus           = qDotN + 2.0_RP*a/gammaMinus1
         entropyConstant = p - a2*Q(1)
!
!        ----------------
!        Resolve solution
!        ----------------
!
         rho   = -(entropyConstant - pExt)/a2
         a     = SQRT(gamma*pExt/rho)
         qDotN = rPlus - 2.0_RP*a/gammaMinus1
         u     = qTanx + qDotN*nHat(1)
         v     = qTany + qDotN*nHat(2)
         w     = qTanz + qDotN*nHat(3)
         
         Q(1) = rho
         Q(2) = rho*u
         Q(3) = rho*v
         Q(4) = rho*w
         Q(5) = pExt/gammaMinus1 + 0.5_RP*rho*(u*u + v*v + w*w)
        
      END IF

      END SUBROUTINE ExternalPressureState
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE computeBoundaryFlux(elementOnLeft, faceID, time, externalStateProcedure)
      USE ElementClass
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(Element)           :: elementOnLeft
      INTEGER                 :: faceID
      REAL(KIND=RP)           :: time
      EXTERNAL                :: externalStateProcedure
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER                         :: N
      REAL(KIND=RP)                   :: bvExt(N_EQN), flux(N_EQN)
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
            
      N            = elementOnLeft % N
      boundaryType = bcTypeDictionary % stringValueForKey(key = elementOnLeft % boundaryName(faceID),&
                                                          requestedLength = BC_STRING_LENGTH)
      
      IF ( boundaryType == "OutflowSpecifyP" )     THEN
         DO j = 0, N
            DO i = 0, N
               bvExt = elementOnLeft % Qb(:,i,j,faceID)
               CALL PressureRiemannFlux(Q    = bvExt, &
                                        nHat = elementOnLeft % geom % normal(:,i,j,faceID),&
                                        flux = flux)
               elementOnLeft % FStarb(:,i,j,faceID) = flux*elementOnLeft % geom % scal(i,j,faceID)
            END DO   
         END DO   
      ELSE 
         DO j = 0, N
            DO i = 0, N
               bvExt = elementOnLeft % Qb(:,i,j,faceID)
               CALL externalStateProcedure( elementOnLeft % geom % xb(:,i,j,faceID), &
                                            time, &
                                            elementOnLeft % geom % normal(:,i,j,faceID), &
                                            bvExt,&
                                            boundaryType )
               CALL RiemannSolver(QLeft  = elementOnLeft % Qb(:,i,j,faceID), &
                                  QRight = bvExt, &
                                  nHat   = elementOnLeft % geom % normal(:,i,j,faceID), &
                                  flux   = flux)
               elementOnLeft % FStarb(:,i,j,faceID) = flux*elementOnLeft % geom % scal(i,j,faceID)
            END DO   
         END DO   
      END IF 
      
      
      

      END SUBROUTINE computeBoundaryFlux
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PressureRiemannFlux( Q, nHat, flux )
!
!     -------------------------------------------------------------------
!     Compute the normal flux given the resolved *state* from the pressure 
!     solver.
!     -------------------------------------------------------------------
!
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP) :: Q(N_EQN), nHat(3), flux(N_EQN)
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP) :: fx(N_EQN), fy(N_EQN), fz(N_EQN)
      
      CALL xFlux( Q, fx)
      CALL yFlux( Q, fy)
      CALL zFlux( Q, fz)
      
      flux = fx*nHat(1) + fy*nHat(2) + fz*nHat(3)

      END SUBROUTINE PressureRiemannFlux

!     ===========
      END MODULE BoundaryConditionFunctions
!     ===========
!
!=====================================================================================================
!=====================================================================================================
!
!
      SUBROUTINE externalStateForBoundaryName( x, t, nHat, Q, boundaryName )
!
!     ----------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external state for each boundary.
!     ----------------------------------------------
!
      USE BoundaryConditionFunctions
      USE MeshTypes
      
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: Q(N_EQN)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
!
!     ---------------
!     Local variables
!     ---------------
!
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
      REAL(KIND=RP)                   :: pExt
      
      boundaryType = bcTypeDictionary % stringValueForKey(key             = boundaryName, &
                                                          requestedLength = BC_STRING_LENGTH)
      
      IF ( boundarytype == "FreeSlipWall" )             THEN
         CALL FreeSlipWallState( x, t, nHat, Q )
      ELSE IF ( boundaryType == "NoSlipAdiabaticWall" ) THEN 
         CALL  NoSlipAdiabaticWallState( x, t, Q)
      ELSE IF ( boundarytype == "NoSlipIsothermalWall") THEN 
         CALL NoSlipIsothermalWallState( x, t, Q )
      ELSE IF ( boundaryType == "OutflowSpecifyP" )     THEN 
         pExt =  ExternalPressureForBoundaryNamed(boundaryName)
         CALL ExternalPressureState ( x, t, nHat, Q, pExt )
      ELSE
         CALL UniformFlowState( x, t, Q )
      END IF

      END SUBROUTINE externalStateForBoundaryName
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExternalGradientForBoundaryName( x, t, nHat, U_x, U_y, U_z, boundaryName )
!
!     ------------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external gradients on each boundary.
!     ------------------------------------------------
!
      USE BoundaryConditionFunctions
      USE MeshTypes
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
!
!     ---------------
!     Local variables
!     ---------------
!
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
      
      boundaryType = bcTypeDictionary % stringValueForKey(key             = boundaryName, &
                                                          requestedLength = BC_STRING_LENGTH)
      
      IF ( boundarytype == "FreeSlipWall" )                   THEN
         CALL FreeSlipNeumann( x, t, nHat, U_x, U_y, U_z )
      ELSE IF ( boundaryType == "NoSlipAdiabaticWall" )       THEN 
         CALL  NoSlipAdiabaticWallNeumann( x, t, nHat, U_x, U_y, U_z)
      ELSE IF ( boundarytype == "NoSlipIsothermalWall")       THEN 
         CALL NoSlipIsothermalWallNeumann( x, t, nHat, U_x, U_y, U_z )
      ELSE
         CALL UniformFlowNeumann( x, t, nHat, U_x, U_y, U_z )
      END IF

      END SUBROUTINE ExternalGradientForBoundaryName
