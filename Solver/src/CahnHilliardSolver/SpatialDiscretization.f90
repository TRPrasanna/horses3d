!
!//////////////////////////////////////////////////////
!
!   @File:    SpatialDiscretization.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 17:14:44 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
!      DGTimeDerivativeRoutines.f95
!      Created: 2008-07-13 16:13:12 -0400 
!      By: David Kopriva  
!
!      3D version by D.A. Kopriva 6/17/15, 12:35 PM
!
!
!////////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use InviscidMethods
      use ViscousMethods
      use LESModels
      use SpectralVanishingViscosity
      use DGWeakIntegrals
      use MeshTypes
      use HexMeshClass
      use ElementClass
      use PhysicsStorage
      use MPI_Face_Class
      use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif


      abstract interface
         SUBROUTINE computeElementInterfaceFluxF(f)
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f   
         end subroutine computeElementInterfaceFluxF

         SUBROUTINE computeMPIFaceFluxF(f)
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f   
         end subroutine computeMPIFaceFluxF

         SUBROUTINE computeBoundaryFluxF(f, time, externalStateProcedure , externalGradientsProcedure)
            use SMConstants
            use FaceClass
            IMPLICIT NONE
            type(Face),    intent(inout) :: f
            REAL(KIND=RP)                :: time
            EXTERNAL                     :: externalStateProcedure
            EXTERNAL                     :: externalGradientsProcedure
         end subroutine computeBoundaryFluxF
      end interface

      procedure(computeElementInterfaceFluxF), pointer :: computeElementInterfaceFlux => computeElementInterfaceFlux_NS
      procedure(computeMPIFaceFluxF),          pointer :: computeMPIFaceFlux          => computeMPIFaceFlux_NS
      procedure(computeBoundaryFluxF),         pointer :: computeBoundaryFlux         => computeBoundaryFlux_NS

!
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods(controlVariables, mesh)
         use PhysicsStorage
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         implicit none
         class(FTValueDictionary),  intent(in)  :: controlVariables
         class(HexMesh)                         :: mesh
         character(len=LINE_LENGTH)       :: inviscidDiscretization
         character(len=LINE_LENGTH)       :: viscousDiscretization

         if ( MPI_Process % isRoot ) then
            write(STD_OUT,'(/)')
            call Section_Header("Spatial discretization scheme")
            write(STD_OUT,'(/)')
         end if
!
!        Initialize inviscid discretization
!        ----------------------------------
         inviscidDiscretization = controlVariables % stringValueForKey(inviscidDiscretizationKey,requestedLength = LINE_LENGTH)

         call toLower(inviscidDiscretization)
      
         select case ( trim(inviscidDiscretization) )

         case ( "standard" )
            if (.not. allocated(InviscidMethod)) allocate( StandardDG_t  :: InviscidMethod )

         case ( "split-form")
            if (.not. allocated(InviscidMethod)) allocate(SplitDG_t   :: InviscidMethod)

         case default
            write(STD_OUT,'(A,A,A)') 'Requested inviscid discretization "',trim(inviscidDiscretization),'" is not implemented.'
            write(STD_OUT,'(A)') "Implemented discretizations are:"
            write(STD_OUT,'(A)') "  * Standard"
            write(STD_OUT,'(A)') "  * Split-Form"
            errorMessage(STD_OUT)
            stop 

         end select
            
         call InviscidMethod % Initialize(controlVariables)
!
!        Initialize viscous discretization
!        ---------------------------------         
         if ( flowIsNavierStokes ) then
            viscousDiscretization = controlVariables % stringValueForKey(viscousDiscretizationKey, requestedLength = LINE_LENGTH)
            call toLower(viscousDiscretization)
            
            select case ( trim(viscousDiscretization) )
            case("br1")
               if (.not. allocated(ViscousMethod)) allocate( BassiRebay1_t :: ViscousMethod  ) 

            case("br2")
               if (.not. allocated(ViscousMethod)) allocate( BassiRebay2_t :: ViscousMethod  ) 

            case("ip")
               if (.not. allocated(ViscousMethod)) allocate( InteriorPenalty_t :: ViscousMethod  ) 

            case default
               write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(viscousDiscretization),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented discretizations are:"
               write(STD_OUT,'(A)') "  * BR1"
               write(STD_OUT,'(A)') "  * BR2"
               write(STD_OUT,'(A)') "  * IP"
               errorMessage(STD_OUT)
               stop 

            end select
   
         else
            if (.not. allocated(ViscousMethod)) allocate( ViscousMethod_t  :: ViscousMethod )
            
         end if

         call ViscousMethod % Initialize(controlVariables)
!
!        Initialize models
!        -----------------
         call InitializeLESModel(LESModel, controlVariables)
!
!        Compute wall distances
!        ----------------------
         call mesh % ComputeWallDistances
!
!        Initialize SVV
!        --------------
         call InitializeSVV(SVV, controlVariables, mesh)

         if ( SVV % enabled ) then
            computeElementInterfaceFlux => computeElementInterfaceFlux_SVV
            computeMPIFaceFlux          => computeMPIFaceFlux_SVV
            computeBoundaryFlux         => computeBoundaryFlux_SVV

         end if
            
         
      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_ComputeQDot( mesh , t, externalState, externalGradients )
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
         external                   :: externalState, externalGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID
         interface
            subroutine UserDefinedSourceTerm(mesh, time, thermodynamics_, dimensionless_, refValues_)
               USE HexMeshClass
               use PhysicsStorage
               IMPLICIT NONE
               CLASS(HexMesh)                        :: mesh
               REAL(KIND=RP)                         :: time
               type(Thermodynamics_t),    intent(in) :: thermodynamics_
               type(Dimensionless_t),     intent(in) :: dimensionless_
               type(RefValues_t),         intent(in) :: refValues_
            end subroutine UserDefinedSourceTerm
         end interface
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_VolumetricContribution( mesh % elements(eID) , t)
         end do
!$omp end do nowait
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) 
         do fID = 1, size(mesh % faces) 
            associate( f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_INTERIOR) 
               CALL computeElementInterfaceFlux( f ) 
 
            case (HMESH_BOUNDARY) 
               CALL computeBoundaryFlux(f, t, externalState, externalGradients) 
 
            end select 
            end associate 
         end do 
!$omp end do 
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
! 
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( e % hasSharedFaces ) cycle
            call TimeDerivative_FacesContribution(e, t, mesh) 
 
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
            end do         ; end do          ; end do 
            end associate 
         end do
!$omp end do
!
!        ****************************
!        Wait until messages are sent
!        ****************************
!
!$omp single
         if ( MPI_Process % doMPIAction ) then
            if ( flowIsNavierStokes ) then 
               call mesh % GatherMPIFacesGradients
            else  
               call mesh % GatherMPIFacesSolution
            end if          
         end if
!$omp end single
!
!        **************************************
!        Compute Riemann solver of shared faces
!        **************************************
!
!$omp do schedule(runtime) 
         do fID = 1, size(mesh % faces) 
            associate( f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_MPI) 
               CALL computeMPIFaceFlux( f ) 
            end select 
            end associate 
         end do 
!$omp end do 
!
!        ***********************************************************
!        Surface integrals and scaling of elements with shared faces
!        ***********************************************************
! 
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( .not. e % hasSharedFaces ) cycle
            call TimeDerivative_FacesContribution(e, t, mesh) 
 
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
            end do         ; end do          ; end do 
            end associate 
         end do
!$omp end do
!
!        Add a source term
!        -----------------
         call UserDefinedSourceTerm(mesh, t, thermodynamics, dimensionless, refValues)
!
!        Add a MPI Barrier
!        -----------------
#ifdef _HAS_MPI_
!$omp single
         if ( MPI_Process % doMPIAction ) call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
#endif

      end subroutine TimeDerivative_ComputeQDot
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_VolumetricContribution( e , t )
         use HexMeshClass
         use ElementClass
         use PhysicsStorage
         implicit none
         type(Element)      :: e
         real(kind=RP)      :: t

!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: inviscidContravariantFlux ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: viscousContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: SVVContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call InviscidMethod % ComputeInnerFluxes ( e , inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         if ( .not. LESModel % active ) then
!
!           Without LES model
!           -----------------
            call ViscousMethod  % ComputeInnerFluxes ( e , viscousContravariantFlux  ) 

         else
!
!           With LES model
!           --------------
            call ViscousMethod  % ComputeInnerFluxesWithSGS ( e , viscousContravariantFlux  ) 

         end if
!
!        Compute the SVV dissipation
!        ---------------------------
         if ( .not. SVV % enabled ) then
            SVVcontravariantFlux = 0.0_RP
         else
            call SVV % ComputeInnerFluxes(e, SVVContravariantFlux)
         end if
!
!        ************************
!        Perform volume integrals
!        ************************
!
         select type ( InviscidMethod )
         type is (StandardDG_t)
!
!           Compute the total Navier-Stokes flux
!           ------------------------------------
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux - SVVContravariantFlux
!
!           Perform the Weak Volume Green integral
!           --------------------------------------
            e % storage % QDot = ScalarWeakIntegrals % StdVolumeGreen ( e , contravariantFlux ) 

         type is (SplitDG_t)
!
!           Compute sharp fluxes for skew-symmetric approximations
!           ------------------------------------------------------
            call InviscidMethod % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!
!           Peform the Weak volume green integral
!           -------------------------------------
            viscousContravariantFlux = viscousContravariantFlux + SVVContravariantFlux

            e % storage % QDot = -ScalarWeakIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

         end select

      end subroutine TimeDerivative_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_FacesContribution( e , t , mesh)
         use HexMeshClass
         use PhysicsStorage
         implicit none
         type(Element)           :: e
         real(kind=RP)           :: t
         type(HexMesh)           :: mesh

         e % storage % QDot = e % storage % QDot - ScalarWeakIntegrals % StdFace( e, &
                      mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % fStar, &
                      mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % fStar, &
                      mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % fStar, &
                      mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % fStar, &
                      mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % fStar, &
                      mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % fStar )

      end subroutine TimeDerivative_FacesContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
!        Riemann solver drivers 
!        ---------------------- 
! 
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE computeElementInterfaceFlux_NS(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: inv_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))

         if ( .not. LESModel % active ) then
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousMethod % RiemannSolver(f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do

         else
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousMethod % RiemannSolverWithSGS(f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do
         end if

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Invscid fluxes
!              --------------
!      
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )

               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         call f % ProjectFluxToElements(flux, (/1,2/))

      END SUBROUTINE computeElementInterfaceFlux_NS

      SUBROUTINE computeMPIFaceFlux_NS(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
!
!        --------------
!        Invscid fluxes
!        --------------
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousMethod % RiemannSolver(f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

!
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(flux, (/thisSide, HMESH_NONE/))

      end subroutine ComputeMPIFaceFlux_NS

      SUBROUTINE computeBoundaryFlux_NS(f, time, externalStateProcedure , externalGradientsProcedure)
      USE ElementClass
      use FaceClass
      USE ViscousMethods
      USE Physics
      use PhysicsStorage
      USE BoundaryConditionFunctions
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
      EXTERNAL                     :: externalStateProcedure
      EXTERNAL                     :: externalGradientsProcedure
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(N_EQN)
      REAL(KIND=RP)                   :: UGradExt(NDIM , N_GRAD_EQN) 
      real(kind=RP)                   :: visc_flux(N_EQN, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(N_EQN, 0:f % Nf(1), 0: f % Nf(2))
      
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
            
      boundaryType = f % boundaryType
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL externalStateProcedure( f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j),&
                                      boundaryType )

      end do               ; end do

      if ( flowIsNavierStokes ) then
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            UGradExt(IX,:) = f % storage(1) % U_x(:,i,j)
            UGradExt(IY,:) = f % storage(1) % U_y(:,i,j)
            UGradExt(IZ,:) = f % storage(1) % U_z(:,i,j)

            CALL externalGradientsProcedure(  f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              UGradExt,&
                                              boundaryType )

            f % storage(2) % U_x(:,i,j) = UGradExt(IX,:)
            f % storage(2) % U_y(:,i,j) = UGradExt(IY,:)
            f % storage(2) % U_z(:,i,j) = UGradExt(IZ,:)

!   
!           --------------
!           Viscous fluxes
!           --------------
!   
            CALL ViscousMethod % RiemannSolver(f = f, &
                                               QLeft = f % storage(1) % Q(:,i,j), &
                                               QRight = f % storage(2) % Q(:,i,j), &
                                               U_xLeft = f % storage(1) % U_x(:,i,j), &
                                               U_yLeft = f % storage(1) % U_y(:,i,j), &
                                               U_zLeft = f % storage(1) % U_z(:,i,j), &
                                               U_xRight = f % storage(2) % U_x(:,i,j), &
                                               U_yRight = f % storage(2) % U_y(:,i,j), &
                                               U_zRight = f % storage(2) % U_z(:,i,j), &
                                               nHat = f % geom % normal(:,i,j) , &
                                               dWall = f % geom % dWall(i,j), &
                                               flux  = visc_flux(:,i,j) )

         end do               ; end do
      else
         visc_flux = 0.0_RP

      end if

      DO j = 0, f % Nf(2)
         DO i = 0, f % Nf(1)
!
!           Inviscid part
!           -------------
            CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &   
                               nHat   = f % geom % normal(:,i,j), &
                               t1     = f % geom % t1(:,i,j), &
                               t2     = f % geom % t2(:,i,j), &
                               flux   = inv_flux)

            fStar(:,i,j) = (inv_flux - visc_flux(:,i,j) ) * f % geom % jacobian(i,j)
         END DO   
      END DO   

      call f % ProjectFluxToElements(fStar, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux_NS

      SUBROUTINE computeElementInterfaceFlux_SVV(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: inv_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: SVV_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
!
!        ----------
!        SVV fluxes
!        ----------
!
         if ( SVV % enabled ) then 
         CALL SVV % RiemannSolver(f = f, &
                              QLeft = f % storage(1) % Q, &
                             QRight = f % storage(2) % Q, &
                            U_xLeft = f % storage(1) % U_x, &
                            U_yLeft = f % storage(1) % U_y, &
                            U_zLeft = f % storage(1) % U_z, &
                           U_xRight = f % storage(2) % U_x, &
                           U_yRight = f % storage(2) % U_y, &
                           U_zRight = f % storage(2) % U_z, &
                              flux  = SVV_flux               )
         else
            SVV_Flux = 0.0_RP

         end if
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousMethod % RiemannSolver(f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

!      
!              --------------
!              Invscid fluxes
!              --------------
!      
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )

               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j) - SVV_flux(:,i,j) ) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         call f % ProjectFluxToElements(flux, (/1,2/))

      END SUBROUTINE computeElementInterfaceFlux_SVV

      SUBROUTINE computeMPIFaceFlux_SVV(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: SVV_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
!
!        ----------
!        SVV fluxes
!        ----------
!
         CALL SVV % RiemannSolver(f = f, &
                              QLeft = f % storage(1) % Q, &
                             QRight = f % storage(2) % Q, &
                            U_xLeft = f % storage(1) % U_x, &
                            U_yLeft = f % storage(1) % U_y, &
                            U_zLeft = f % storage(1) % U_z, &
                           U_xRight = f % storage(2) % U_x, &
                           U_yRight = f % storage(2) % U_y, &
                           U_zRight = f % storage(2) % U_z, &
                              flux  = SVV_flux               )

!
!        --------------
!        Invscid fluxes
!        --------------
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousMethod % RiemannSolver(f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

!
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j) - SVV_flux(:,i,j) ) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(flux, (/thisSide, HMESH_NONE/))


      end subroutine ComputeMPIFaceFlux_SVV

      SUBROUTINE computeBoundaryFlux_SVV(f, time, externalStateProcedure , externalGradientsProcedure)
      USE ElementClass
      use FaceClass
      USE ViscousMethods
      USE Physics
      use PhysicsStorage
      USE BoundaryConditionFunctions
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
      EXTERNAL                     :: externalStateProcedure
      EXTERNAL                     :: externalGradientsProcedure
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(N_EQN)
      REAL(KIND=RP)                   :: UGradExt(NDIM , N_GRAD_EQN) 
      real(kind=RP)                   :: visc_flux(N_EQN, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: SVV_flux(N_EQN, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(N_EQN, 0:f % Nf(1), 0: f % Nf(2))
      
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
            
      boundaryType = f % boundaryType
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL externalStateProcedure( f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j),&
                                      boundaryType )

      end do               ; end do

      if ( flowIsNavierStokes ) then
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            UGradExt(IX,:) = f % storage(1) % U_x(:,i,j)
            UGradExt(IY,:) = f % storage(1) % U_y(:,i,j)
            UGradExt(IZ,:) = f % storage(1) % U_z(:,i,j)

            CALL externalGradientsProcedure(  f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              UGradExt,&
                                              boundaryType )

            f % storage(2) % U_x(:,i,j) = UGradExt(IX,:)
            f % storage(2) % U_y(:,i,j) = UGradExt(IY,:)
            f % storage(2) % U_z(:,i,j) = UGradExt(IZ,:)

!   
!           --------------
!           Viscous fluxes
!           --------------
!   
            CALL ViscousMethod % RiemannSolver(f = f, &
                                               QLeft = f % storage(1) % Q(:,i,j), &
                                               QRight = f % storage(2) % Q(:,i,j), &
                                               U_xLeft = f % storage(1) % U_x(:,i,j), &
                                               U_yLeft = f % storage(1) % U_y(:,i,j), &
                                               U_zLeft = f % storage(1) % U_z(:,i,j), &
                                               U_xRight = f % storage(2) % U_x(:,i,j), &
                                               U_yRight = f % storage(2) % U_y(:,i,j), &
                                               U_zRight = f % storage(2) % U_z(:,i,j), &
                                               nHat = f % geom % normal(:,i,j) , &
                                               dWall = f % geom % dWall(i,j), &
                                               flux  = visc_flux(:,i,j) )

         end do               ; end do
      else
         visc_flux = 0.0_RP

      end if
!
!     ----------
!     SVV fluxes
!     ----------
!
      CALL SVV % RiemannSolver(f = f, &
                           QLeft = f % storage(1) % Q, &
                          QRight = f % storage(2) % Q, &
                         U_xLeft = f % storage(1) % U_x, &
                         U_yLeft = f % storage(1) % U_y, &
                         U_zLeft = f % storage(1) % U_z, &
                        U_xRight = f % storage(2) % U_x, &
                        U_yRight = f % storage(2) % U_y, &
                        U_zRight = f % storage(2) % U_z, &
                           flux  = SVV_flux               )

      DO j = 0, f % Nf(2)
         DO i = 0, f % Nf(1)
!
!           Inviscid part
!           -------------
            CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &   
                               nHat   = f % geom % normal(:,i,j), &
                               t1     = f % geom % t1(:,i,j), &
                               t2     = f % geom % t2(:,i,j), &
                               flux   = inv_flux)

            fStar(:,i,j) = (inv_flux - visc_flux(:,i,j) - SVV_flux(:,i,j) ) * f % geom % jacobian(i,j)
         END DO   
      END DO   

      call f % ProjectFluxToElements(fStar, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux_SVV
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_ComputeGradient( mesh , time , externalStateProcedure , externalGradientsProcedure )
         use HexMeshClass
         use PhysicsStorage
         implicit none
         type(HexMesh)                  :: mesh
         real(kind=RP),      intent(in) :: time
         interface
            subroutine externalStateProcedure(x,t,nHat,Q,boundaryName)
               use SMConstants
               real(kind=RP)   , intent(in)    :: x(3), t, nHat(3)
               real(kind=RP)   , intent(inout) :: Q(:)
               character(len=*), intent(in)    :: boundaryName
            end subroutine externalStateProcedure
            
            subroutine externalGradientsProcedure(x,t,nHat,gradU,boundaryName)
               use SMConstants
               real(kind=RP)   , intent(in)    :: x(3), t, nHat(3)
               real(kind=RP)   , intent(inout) :: gradU(:,:)
               character(len=*), intent(in)    :: boundaryName
            end subroutine externalGradientsProcedure
         end interface

         call ViscousMethod % ComputeGradient( mesh , time , externalStateProcedure , externalGradientsProcedure )

      end subroutine DGSpatial_ComputeGradient
!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization