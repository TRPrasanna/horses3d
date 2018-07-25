!
!//////////////////////////////////////////////////////
!
!   @File:    BoundaryConditions_iNS.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jun 19 17:39:25 2018
!   @Last revision date: Mon Jul 23 10:59:35 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: b0edd55b642212b62cae102b966c37b726378791
!
!//////////////////////////////////////////////////////
!
MODULE BoundaryConditionFunctions_iNS
   USE SMConstants
   USE Physics_iNS
   USE SharedBCModule
   use PhysicsStorage_iNS
   use FluidData_iNS

   private

   public implementediNSBCNames
   public NoSlipWallState, NoSlipWallNeumann, UserDefinedNeumann, UserDefinedState
   public InflowState, InflowNeumann
   public OutflowState, OutflowNeumann
   public FreeSlipWallState, FreeSlipWallNeumann

   CHARACTER(LEN=BC_STRING_LENGTH), DIMENSION(7) :: implementediNSBCNames = &
         ["periodic-           ", &
          "periodic+           ", &
          "user-defined        ", &
          "noslipwall          ", &
          "inflow              ", &
          "outflow             ", &
          "freeslipwall        "]

   enum, bind(C)
      enumerator :: PERIODIC_PLUS_INDEX=1, PERIODIC_MINUS_INDEX
      enumerator :: USERDEFINED_INDEX, NOSLIPWALL_INDEX
      enumerator :: INFLOW_INDEX, OUTFLOW_INDEX
      enumerator :: FREESLIPWALL_INDEX
   end enum
!
!  ========         
   contains
!  ========
! 
      SUBROUTINE NoSlipWallState( x, t, nHat, Q )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), INTENT(IN)    :: x(3), t
         real(kind=RP), intent(in)    :: nHat(3)
         REAL(KIND=RP), INTENT(INOUT) :: Q(NINC)
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

      END SUBROUTINE NoSlipWallState

      subroutine InflowState(x, t, nHat, Q)
         implicit none
         REAL(KIND=RP), INTENT(IN)    :: x(3), t
         real(kind=RP), intent(in)    :: nHat(3)
         REAL(KIND=RP), INTENT(INOUT) :: Q(NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: u, v, w, theta, phi

         theta = refValues % AoATheta * PI / 180.0_RP
         phi   = refValues % AoAPhi   * PI / 180.0_RP

         u = cos(theta) * cos(phi)
         v = sin(theta) * cos(phi)
         w = sin(phi)

         Q(INSRHO) = 1.0_RP
         Q(INSRHOU) = Q(INSRHO)*u
         Q(INSRHOV) = Q(INSRHO)*v
         Q(INSRHOW) = Q(INSRHO)*w

      end subroutine InflowState

      subroutine OutflowState(x, t, nHat, Q)
         implicit none
         REAL(KIND=RP), INTENT(IN)    :: x(3), t
         real(kind=RP), intent(in)    :: nHat(3)
         REAL(KIND=RP), INTENT(INOUT) :: Q(NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: u, v, w, theta, phi, un

         un = Q(INSRHOU)*nHat(IX) + Q(INSRHOV)*nHat(IY) + Q(INSRHOW)*nHat(IZ)
         
         if ( un .ge. -1.0e-4_RP ) then
         
            Q(INSP) = 0.0_RP

         else

            theta = refValues % AoATheta * PI / 180.0_RP
            phi   = refValues % AoAPhi   * PI / 180.0_RP
   
            u = cos(theta) * cos(phi)
            v = sin(theta) * cos(phi)
            w = sin(phi)

            Q(INSRHOU) = Q(INSRHO)*u
            Q(INSRHOV) = Q(INSRHO)*v
            Q(INSRHOW) = Q(INSRHO)*w

         end if

      end subroutine OutflowState

      SUBROUTINE FreeSlipWallState( x, t, nHat, Q )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), INTENT(IN)    :: x(3), t
         real(kind=RP), intent(in)    :: nHat(3)
         REAL(KIND=RP), INTENT(INOUT) :: Q(NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: vn
!
!        -----------------------------------------------
!        Generate the external flow along the face, that
!        represents a solid wall.
!        -----------------------------------------------
!
         vn = sum(Q(INSRHOU:INSRHOW)*nHat)

         Q(INSRHO)  = Q(INSRHO)
         Q(INSRHOU) = Q(INSRHOU) - 2.0_RP * vn * nHat(IX)
         Q(INSRHOV) = Q(INSRHOV) - 2.0_RP * vn * nHat(IY)
         Q(INSRHOW) = Q(INSRHOW) - 2.0_RP * vn * nHat(IZ)
         Q(INSP)    = Q(INSP)

      END SUBROUTINE FreeSlipWallState
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE NoSlipWallNeumann( x, t, nHat, U_x, U_y, u_z )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), INTENT(IN)    :: x(3), t
         REAL(KIND=RP), INTENT(IN)    :: nHat(3)
         REAL(KIND=RP), INTENT(INOUT) :: U_x(NINC), U_y(NINC), U_z(NINC)
!
!        ---------------
!        Local Variables
!        ---------------
!
!
         REAL(KIND=RP) :: gradUNorm, UTanx, UTany, UTanz
!
!
!        Remove the normal component of the density gradient
!        ---------------------------------------------------
         gradUNorm =  nHat(1)*U_x(INSRHO) + nHat(2)*U_y(INSRHO)+ nHat(3)*U_z(INSRHO)
         UTanx = U_x(INSRHO) - gradUNorm*nHat(1)
         UTany = U_y(INSRHO) - gradUNorm*nHat(2)
         UTanz = U_z(INSRHO) - gradUNorm*nHat(3)
   
         U_x(INSRHO) = UTanx - gradUNorm*nHat(1)
         U_y(INSRHO) = UTany - gradUNorm*nHat(2)
         U_z(INSRHO) = UTanz - gradUNorm*nHat(3)
!
!        Remove the normal component of the pressure gradient
!        ----------------------------------------------------
         gradUNorm =  nHat(1)*U_x(INSP) + nHat(2)*U_y(INSP)+ nHat(3)*U_z(INSP)
         UTanx = U_x(INSP) - gradUNorm*nHat(1)
         UTany = U_y(INSP) - gradUNorm*nHat(2)
         UTanz = U_z(INSP) - gradUNorm*nHat(3)
   
         U_x(INSP) = UTanx - gradUNorm*nHat(1)
         U_y(INSP) = UTany - gradUNorm*nHat(2)
         U_z(INSP) = UTanz - gradUNorm*nHat(3)

      
      END SUBROUTINE NoSlipWallNeumann

      subroutine InflowNeumann(x, t, nHat, U_x, U_y, U_z)
!
!        **********************************
!           Set drho/dn = 0, and dp/dn = 0
!        **********************************
!
         implicit none
         REAL(KIND=RP), INTENT(IN)    :: x(3), t
         REAL(KIND=RP), INTENT(IN)    :: nHat(3)
         REAL(KIND=RP), INTENT(INOUT) :: U_x(NINC), U_y(NINC), U_z(NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: dpdn, drhodn

         drhodn = (U_x(INSRHO)*nHat(IX) + U_y(INSRHO)*nHat(IY) + U_z(INSRHO)*nHat(IZ))
         U_x(INSRHO) = U_x(INSRHO) - 2.0_RP * drhodn * nHat(IX)
         U_y(INSRHO) = U_y(INSRHO) - 2.0_RP * drhodn * nHat(IY) 
         U_z(INSRHO) = U_z(INSRHO) - 2.0_RP * drhodn * nHat(IZ) 

         dpdn = (U_x(INSP)*nHat(IX) + U_y(INSP)*nHat(IY) + U_z(INSP)*nHat(IZ))
         U_x(INSP) = U_x(INSP) - 2.0_RP * dpdn * nHat(IX)
         U_y(INSP) = U_y(INSP) - 2.0_RP * dpdn * nHat(IY) 
         U_z(INSP) = U_z(INSP) - 2.0_RP * dpdn * nHat(IZ) 

      end subroutine InflowNeumann

      subroutine OutflowNeumann(x, t, nHat, U_x, U_y, U_z)
!
!        ********************************************************
!           Set drho/dn = 0, and du/dn = 0
!        ********************************************************
!
         implicit none
         REAL(KIND=RP), INTENT(IN)    :: x(3), t
         REAL(KIND=RP), INTENT(IN)    :: nHat(3)
         REAL(KIND=RP), INTENT(INOUT) :: U_x(NINC), U_y(NINC), U_z(NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: drhodn, dudn, dvdn, dwdn

         drhodn = (U_x(INSRHO)*nHat(IX) + U_y(INSRHO)*nHat(IY) + U_z(INSRHO)*nHat(IZ))
         U_x(INSRHO) = U_x(INSRHO) - 2.0_RP * drhodn * nHat(IX)
         U_y(INSRHO) = U_y(INSRHO) - 2.0_RP * drhodn * nHat(IY) 
         U_z(INSRHO) = U_z(INSRHO) - 2.0_RP * drhodn * nHat(IZ) 

         dudn = (U_x(INSRHOU)*nHat(IX) + U_y(INSRHOU)*nHat(IY) + U_z(INSRHOU)*nHat(IZ))
         U_x(INSRHOU) = U_x(INSRHOU) - 2.0_RP * dudn * nHat(IX)
         U_y(INSRHOU) = U_y(INSRHOU) - 2.0_RP * dudn * nHat(IY) 
         U_z(INSRHOU) = U_z(INSRHOU) - 2.0_RP * dudn * nHat(IZ) 

         dvdn = (U_x(INSRHOV)*nHat(IX) + U_y(INSRHOV)*nHat(IY) + U_z(INSRHOV)*nHat(IZ))
         U_x(INSRHOV) = U_x(INSRHOV) - 2.0_RP * dvdn * nHat(IX)
         U_y(INSRHOV) = U_y(INSRHOV) - 2.0_RP * dvdn * nHat(IY) 
         U_z(INSRHOV) = U_z(INSRHOV) - 2.0_RP * dvdn * nHat(IZ) 

         dwdn = (U_x(INSRHOW)*nHat(IX) + U_y(INSRHOW)*nHat(IY) + U_z(INSRHOW)*nHat(IZ))
         U_x(INSRHOW) = U_x(INSRHOW) - 2.0_RP * dwdn * nHat(IX)
         U_y(INSRHOW) = U_y(INSRHOW) - 2.0_RP * dwdn * nHat(IY) 
         U_z(INSRHOW) = U_z(INSRHOW) - 2.0_RP * dwdn * nHat(IZ) 

      end subroutine OutflowNeumann

      subroutine UserDefinedState(x, t, nHat, Q)
         implicit none
         real(kind=RP)  :: x(NDIM)
         real(kind=RP)  :: t
         real(kind=RP)  :: nHat(NDIM)
         real(kind=RP)  :: Q(NINC)
         interface
            subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
               use SMConstants
               use PhysicsStorage_iNS
               use FluidData_iNS
               implicit none
               real(kind=RP)  :: x(NDIM)
               real(kind=RP)  :: t
               real(kind=RP)  :: nHat(NDIM)
               real(kind=RP)  :: Q(NINC)
               type(Thermodynamics_t), intent(in)  :: thermodynamics_
               type(Dimensionless_t),  intent(in)  :: dimensionless_
               type(RefValues_t),      intent(in)  :: refValues_
            end subroutine UserDefinedState1
         end interface

         call UserDefinedState1(x, t, nHat, Q, thermodynamics, dimensionless, refValues)

      end subroutine UserDefinedState

      subroutine UserDefinedNeumann(x, t, nHat, U_x, U_y, U_z)
         implicit none
         real(kind=RP)  :: x(NDIM)
         real(kind=RP)  :: t
         real(kind=RP)  :: nHat(NDIM)
         real(kind=RP)  :: U_x(NINC), U_y(NINC), U_z(NINC)
         interface
            subroutine UserDefinedNeumann1(x, t, nHat, U_x, U_y, U_z, thermodynamics_, dimensionless_, refValues_)
               use SMConstants
               use PhysicsStorage_iNS
               use FluidData_iNS
               implicit none
               real(kind=RP)  :: x(NDIM)
               real(kind=RP)  :: t
               real(kind=RP)  :: nHat(NDIM)
               real(kind=RP)  :: U_x(NINC), U_y(NINC), U_z(NINC)
               type(Thermodynamics_t), intent(in)  :: thermodynamics_
               type(Dimensionless_t),  intent(in)  :: dimensionless_
               type(RefValues_t),      intent(in)  :: refValues_
            end subroutine UserDefinedNeumann1
         end interface

         call UserDefinedNeumann1(x, t, nHat, U_x, U_y, U_z, thermodynamics, dimensionless, refValues)

      end subroutine UserDefinedNeumann
!
!////////////////////////////////////////////////////////////////////////
!
!     ===========
      END MODULE BoundaryConditionFunctions_iNS
!     ===========
!
!=====================================================================================================
!=====================================================================================================
!
!
      SUBROUTINE externalStateForBoundaryName_iNS( nEqn, x, t, nHat, Q, boundaryType, boundaryName )
!
!     ----------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external state for each boundary.
!     ----------------------------------------------
!
      use SMConstants
      USE BoundaryConditionFunctions_iNS
      use PhysicsStorage_iNS
      
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      integer         , intent(in)    :: nEqn
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: Q(nEqn)
      CHARACTER(LEN=BC_STRING_LENGTH), INTENT(IN)    :: boundaryType
      CHARACTER(LEN=BC_STRING_LENGTH), INTENT(IN)    :: boundaryName
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP)   :: pExt
      LOGICAL         :: success

      select case (trim(boundarytype))
      case("noslipwall")
         call NoSlipWallState(x, t, nHat, Q)
      case("freeslipwall")
         call FreeSlipWallState(x, t, nHat, Q)
      case("inflow")
         call InflowState(x, t, nHat, Q)
      case("outflow")
         call OutflowState(x, t, nHat, Q)
      case("user-defined")
         call UserDefinedState(x, t, nHat, Q)
      end select

      END SUBROUTINE externalStateForBoundaryName_iNS
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExternalGradientForBoundaryName_iNS( nGradEqn, x, t, nHat, GradU, boundaryType, boundaryName )
!
!     ------------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external gradients on each boundary.
!     ------------------------------------------------
!
      use SMConstants
      USE BoundaryConditionFunctions_iNS
      use PhysicsStorage_iNS
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      integer,          intent(in)    :: nGradEqn
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: GradU(3,nGradEqn)
      CHARACTER(LEN=BC_STRING_LENGTH), INTENT(IN)    :: boundaryType
      CHARACTER(LEN=BC_STRING_LENGTH), INTENT(IN)    :: boundaryName
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP) :: U_x(NINC), U_y(NINC), U_z(NINC)

      U_x(:) = GradU(1,:)
      U_y(:) = GradU(2,:)
      U_z(:) = GradU(3,:)

      select case (trim(boundarytype))
      case("noslipwall")
         call NoSlipWallNeumann(x, t, nHat, U_x, U_y, U_z)
      case("freeslipwall")
         call NoSlipWallNeumann(x, t, nHat, U_x, U_y, U_z)
      case("inflow")
         call InflowNeumann(x, t, nHat, U_x, U_y, U_z)
      case("outflow")
         call OutflowNeumann(x, t, nHat, U_x, U_y, U_z)
      case("user-defined")
         call UserDefinedNeumann1(x, t, nHat, U_x, U_y, U_z)
      end select

      GradU(1,:) = U_x(:)
      GradU(2,:) = U_y(:)
      GradU(3,:) = U_z(:)

      END SUBROUTINE ExternalGradientForBoundaryName_iNS
