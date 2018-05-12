!
!//////////////////////////////////////////////////////
!
!   @File:    SpatialDiscretization.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 17:14:44 2018
!   @Last revision date: Sat May 12 20:53:24 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: ece79010cbff566c377be7e7026f86a2889a191e
!
!//////////////////////////////////////////////////////
!
#ifdef HAS_PETSC
#include "petsc/finclude/petsc.h"
#endif
#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use Physics
      use EllipticDiscretizations
      use DGWeakIntegrals
      use MeshTypes
      use FaceClass
      use ElementClass
      use HexMeshClass
      use PhysicsStorage
      use MPI_Face_Class
      use MPI_Process_Info
      use DGSEMClass
      use BoundaryConditionFunctions, only: C_BC, MU_BC
      use GradientsStabilization
      use FluidData
      use VariableConversion
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public  ComputeLaplacian, DGSpatial_ComputeGradient
      public  Initialize_SpaceAndTimeMethods, ComputeTimeDerivative, ComputeTimeDerivativeIsolated
      public  ComputeTimeDerivative_onlyLinear, ComputetimeDerivative_onlyNonLinear

      interface GetPoiseuilleFlow
         module procedure GetPoiseuilleFlow_Element, GetPoiseuilleFlow_Face
      end interface 

      logical, parameter   :: enable_speed = .false. 
!
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods(controlVariables, sem)
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         implicit none
         class(FTValueDictionary),  intent(in)  :: controlVariables
         class(DGSEM)                           :: sem
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH) :: inviscidDiscretization
         character(len=LINE_LENGTH) :: viscousDiscretization
         integer                    :: eID, fID

         if ( MPI_Process % isRoot ) then
            write(STD_OUT,'(/)')
            call Section_Header("Spatial discretization scheme")
            write(STD_OUT,'(/)')
         end if
!
!        Initialize viscous discretization
!        ---------------------------------         
         call BassiRebay1     % Initialize(controlVariables)
         call BassiRebay2     % Initialize(controlVariables)
         call InteriorPenalty % Initialize(controlVariables)

         viscousDiscretization = controlVariables % stringValueForKey(viscousDiscretizationKey, requestedLength = LINE_LENGTH)
         call toLower(viscousDiscretization)
         
         select case ( trim(viscousDiscretization) )
         case("br1")
            EllipticDiscretization => BassiRebay1

         case("br2")
            EllipticDiscretization => BassiRebay2

         case("ip")
            EllipticDiscretization => InteriorPenalty

         case default
            write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(viscousDiscretization),'" is not implemented.'
            write(STD_OUT,'(A)') "Implemented discretizations are:"
            write(STD_OUT,'(A)') "  * BR1"
            write(STD_OUT,'(A)') "  * BR2"
            write(STD_OUT,'(A)') "  * IP"
            errorMessage(STD_OUT)
            stop 

         end select

         call EllipticDiscretization % Describe
!
!        Compute wall distances
!        ----------------------
         call sem % mesh % ComputeWallDistances
!
!        Get Poiseuille flow
!        -------------------
         do eID = 1, sem % mesh % no_of_elements
            call GetPoiseuilleFlow(sem % mesh % elements(eID))
         end do

         do fID = 1, size(sem % mesh % faces)   
            call GetPoiseuilleFlow(sem % mesh % faces(fID))
         end do
      
      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeTimeDerivative( mesh, particles, time, BCFunctions)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target      :: mesh
         logical                    :: particles
         REAL(KIND=RP)              :: time
         type(BCFunctions_t), intent(in)  :: BCFunctions(no_of_BCsets)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: i, j, k, eID, fID
         class(Element), pointer  :: e
         class(Face),    pointer  :: f
         interface
            subroutine UserDefinedSourceTermCH(x, time, S, multiphase_)
               use SMConstants
               USE HexMeshClass
               use FluidData
               IMPLICIT NONE
               real(kind=RP),             intent(in)  :: x(NDIM)
               real(kind=RP),             intent(in)  :: time
               real(kind=RP),             intent(out) :: S(NCONS)
               type(Multiphase_t)                     :: multiphase_               
            end subroutine UserDefinedSourceTermCH
         end interface

!
!        **************************************
!        Compute chemical potential: Q stores c
!        **************************************
!
!$omp parallel shared(mesh, time) private(e, i, j, k, eID, fID)
!
!        ------------------------------
!        Change memory to concentration
!        ------------------------------
!
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % storage % SetStorageToCH_c
         end do

         do fID = 1, size(mesh % faces)
            call mesh % faces(fID) % storage(1) % SetStorageToCH_c
            call mesh % faces(fID) % storage(2) % SetStorageToCH_c
         end do
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         CALL DGSpatial_ComputeGradient( mesh , time , BCFunctions(C_BC) % externalState)

#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
!        Linear part
!        -----------
         call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(C_BC) % externalState, &
                  externalGradients = BCFunctions(C_BC) % externalGradients )

!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % mu = - POW2(multiphase % eps) * e % storage % QDot
            call AddQuarticDWPDerivative(e % storage % c, e % storage % mu)
!
!           Move storage to chemical potential
!           ----------------------------------
            call e % storage % SetStorageToCH_mu
         end do
!$omp end do

!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            call mesh % faces(fID) % storage(1) % SetStorageToCH_mu
            call mesh % faces(fID) % storage(2) % SetStorageToCH_mu
         end do
!$omp end do
!
!        *************************
!        Compute cDot: Q stores mu
!        *************************
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         CALL DGSpatial_ComputeGradient( mesh , time , BCFunctions(MU_BC) % externalState)

#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
         call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(MU_BC) % externalState, &
                  externalGradients = BCFunctions(MU_BC) % externalGradients )
!
!        Scale QDot with the Peclet number
!        ---------------------------------
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % QDot = (1.0_RP / multiphase % Pe) * e % storage % QDot
         end do
!$omp end do

!
!        Add a source term
!        -----------------
!$omp do schedule(runtime) private(i,j,k)
      do eID = 1, mesh % no_of_elements
         associate ( e => mesh % elements(eID) )
         do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call UserDefinedSourceTermCH(e % geom % x(:,i,j,k), time, e % storage % S(:,i,j,k), multiphase)
         end do                  ; end do                ; end do
         end associate
      end do
!$omp end do
!
!        *****************************
!        Return the concentration to Q
!        *****************************
!
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % storage % SetStorageToCH_c
         end do
!$omp end do

!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            call mesh % faces(fID) % storage(1) % SetStorageToCH_c
            call mesh % faces(fID) % storage(2) % SetStorageToCH_c
         end do
!$omp end do

!
!        ***********************************
!        Compute the concentration advection
!        ***********************************
!
         if ( enable_speed ) then
!
!        Perform the stabilization
!        -------------------------
         call StabilizeGradients(mesh, time, BCFunctions(C_BC) % externalState)
!
!        Add the velocity field
!        ----------------------
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
         
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)   ; do i = 0, e % Nxyz(1)
               e % storage % cDot(1,i,j,k) = e % storage % cDot(1,i,j,k) &
                                 - e % storage % v(IX,i,j,k) * e % storage % c_x(1,i,j,k) &
                                 - e % storage % v(IY,i,j,k) * e % storage % c_y(1,i,j,k) &
                                 - e % storage % v(IZ,i,j,k) * e % storage % c_z(1,i,j,k) 
            end do                ; end do                  ; end do
         end do
!$omp end do

         end if
!$omp end parallel

      end subroutine ComputeTimeDerivative

      subroutine ComputeTimeDerivative_OnlyLinear( mesh, particles, time, BCFunctions)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target      :: mesh
         logical                    :: particles
         REAL(KIND=RP)              :: time
         type(BCFunctions_t), intent(in)  :: BCFunctions(no_of_BCsets)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: i, j, k, eID, fID
         class(Element), pointer  :: e
         class(Face),    pointer  :: f
!
!        **************************************
!        Compute chemical potential: Q stores c
!        **************************************
!
!$omp parallel shared(mesh, time) private(e, i, j, k, eID, fID)
!
!        ------------------------------
!        Change memory to concentration
!        ------------------------------
!
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % storage % SetStorageToCH_c
         end do

         do fID = 1, size(mesh % faces)
            call mesh % faces(fID) % storage(1) % SetStorageToCH_c
            call mesh % faces(fID) % storage(2) % SetStorageToCH_c
         end do
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         CALL DGSpatial_ComputeGradient( mesh , time , BCFunctions(C_BC) % externalState)

#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
!        Linear part
!        -----------
         call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(C_BC) % externalState, &
                  externalGradients = BCFunctions(C_BC) % externalGradients )

!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % mu = - POW2(multiphase % eps) * e % storage % QDot
!
!           Move storage to chemical potential
!           ----------------------------------
            call e % storage % SetStorageToCH_mu
         end do
!$omp end do

!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            call mesh % faces(fID) % storage(1) % SetStorageToCH_mu
            call mesh % faces(fID) % storage(2) % SetStorageToCH_mu
         end do
!$omp end do
!
!        *************************
!        Compute cDot: Q stores mu
!        *************************
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         CALL DGSpatial_ComputeGradient( mesh , time , BCFunctions(MU_BC) % externalState)

#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
         call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(MU_BC) % externalState, &
                  externalGradients = BCFunctions(MU_BC) % externalGradients )
!
!        Scale QDot with the Peclet number
!        ---------------------------------
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % QDot = (1.0_RP / multiphase % Pe) * e % storage % QDot
         end do
!$omp end do

!
!        *****************************
!        Return the concentration to Q
!        *****************************
!
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % storage % SetStorageToCH_c
         end do
!$omp end do

!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            call mesh % faces(fID) % storage(1) % SetStorageToCH_c
            call mesh % faces(fID) % storage(2) % SetStorageToCH_c
         end do
!$omp end do
!$omp end parallel

      end subroutine ComputeTimeDerivative_OnlyLinear

      subroutine ComputeTimeDerivative_OnlyNonLinear( mesh, particles, time, BCFunctions)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         logical                         :: particles
         REAL(KIND=RP)                   :: time
         type(BCFunctions_t), intent(in) :: BCFunctions(no_of_BCsets)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                 :: i, j, k, eID, fID
         class(Element), pointer :: e
         class(Face),    pointer :: f
         interface
            subroutine UserDefinedSourceTermCH(x, time, S, multiphase_)
               use SMConstants
               USE HexMeshClass
               use FluidData
               IMPLICIT NONE
               real(kind=RP),             intent(in)  :: x(NDIM)
               real(kind=RP),             intent(in)  :: time
               real(kind=RP),             intent(out) :: S(NCONS)
               type(Multiphase_t)                     :: multiphase_               
            end subroutine UserDefinedSourceTermCH
         end interface

!
!        **************************************
!        Compute chemical potential: Q stores c
!        **************************************
!
!$omp parallel shared(mesh, time) private(e, i, j, k, eID, fID)
!
!        ------------------------------
!        Change memory to concentration
!        ------------------------------
!
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % storage % SetStorageToCH_c
         end do

         do fID = 1, size(mesh % faces)
            call mesh % faces(fID) % storage(1) % SetStorageToCH_c
            call mesh % faces(fID) % storage(2) % SetStorageToCH_c
         end do
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         CALL DGSpatial_ComputeGradient( mesh , time , BCFunctions(C_BC) % externalState)

#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
!        Linear part: only Neumann boundary conditions contribution
!        -----------
         call ComputeLaplacianNeumannBCs(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(C_BC) % externalState, &
                  externalGradients = BCFunctions(C_BC) % externalGradients )

!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % mu = - POW2(multiphase % eps) * e % storage % QDot
            call AddQuarticDWPDerivative(e % storage % c, e % storage % mu)
!
!           Move storage to chemical potential
!           ----------------------------------
            call e % storage % SetStorageToCH_mu
         end do
!$omp end do

!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            call mesh % faces(fID) % storage(1) % SetStorageToCH_mu
            call mesh % faces(fID) % storage(2) % SetStorageToCH_mu
         end do
!$omp end do
!
!        *************************
!        Compute cDot: Q stores mu
!        *************************
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         CALL DGSpatial_ComputeGradient( mesh , time , BCFunctions(MU_BC) % externalState)

#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
         call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(MU_BC) % externalState, &
                  externalGradients = BCFunctions(MU_BC) % externalGradients )
!
!        Scale QDot with the Peclet number
!        ---------------------------------
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % QDot = (1.0_RP / multiphase % Pe) * e % storage % QDot
         end do
!$omp end do

!
!        Add a source term
!        -----------------
!$omp do schedule(runtime) private(i,j,k)
      do eID = 1, mesh % no_of_elements
      associate ( e => mesh % elements(eID) )
      do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
         call UserDefinedSourceTermCH(e % geom % x(:,i,j,k), time, e % storage % S(:,i,j,k), multiphase)
      end do                  ; end do                ; end do
      end associate
   end do
!$omp end do
!
!        *****************************
!        Return the concentration to Q
!        *****************************
!
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % storage % SetStorageToCH_c
         end do
!$omp end do

!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            call mesh % faces(fID) % storage(1) % SetStorageToCH_c
            call mesh % faces(fID) % storage(2) % SetStorageToCH_c
         end do
!$omp end do

!
!        ***********************************
!        Compute the concentration advection
!        ***********************************
!
         if ( enable_speed ) then
!
!        Perform the stabilization
!        -------------------------
         call StabilizeGradients(mesh, time, BCFunctions(C_BC) % externalState)
!
!        Add the velocity field
!        ----------------------
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
         
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)   ; do i = 0, e % Nxyz(1)
               e % storage % cDot(1,i,j,k) = e % storage % cDot(1,i,j,k) &
                                 - e % storage % v(IX,i,j,k) * e % storage % c_x(1,i,j,k) &
                                 - e % storage % v(IY,i,j,k) * e % storage % c_y(1,i,j,k) &
                                 - e % storage % v(IZ,i,j,k) * e % storage % c_z(1,i,j,k) 
            end do                ; end do                  ; end do
         end do
!$omp end do

         end if
         
!$omp end parallel

      end subroutine ComputeTimeDerivative_OnlyNonLinear
      
      subroutine ComputeTimeDerivativeIsolated( mesh, particles, time, BCFunctions)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target      :: mesh
         logical                    :: particles
         REAL(KIND=RP)              :: time
         type(BCFunctions_t), intent(in)  :: BCFunctions(no_of_BCsets)
         
         ERROR stop 'ComputeTimeDerivativeIsolated not implemented for Cahn-Hilliard'
      end subroutine ComputeTimeDerivativeIsolated
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!           Procedures to compute the state variables Laplacian
!           ---------------------------------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeLaplacian( mesh , t, externalState, externalGradients )
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
            
            case (HMESH_MPI)

            case default
               print*, "Unrecognized face type"
               errorMessage(STD_OUT)
               stop
                
            end select 
            end associate 
         end do 
!$omp end do 
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
! 
!$omp do schedule(runtime) private(i, j, k)
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
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
            call mesh % GatherMPIFacesGradients
!$omp end single
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
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
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
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
!           Add a MPI Barrier
!           -----------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif

      end subroutine ComputeLaplacian

      subroutine ComputeLaplacianNeumannBCs( mesh , t, externalState, externalGradients )
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
!
!        **************************
!        Reset QDot and face fluxes
!        **************************
!
         do eID = 1, mesh % no_of_elements
            mesh % elements(eID) % storage % QDot = 0.0_RP
         end do
   
         do fID = 1, size(mesh % faces)
            mesh % faces(fID) % storage(1) % genericInterfaceFluxMemory = 0.0_RP
            mesh % faces(fID) % storage(2) % genericInterfaceFluxMemory = 0.0_RP
         end do
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) 
         do fID = 1, size(mesh % faces) 
            associate( f => mesh % faces(fID)) 
            select case (f % faceType) 
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
!$omp do schedule(runtime) private(i, j, k)
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
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
            call mesh % GatherMPIFacesGradients
!$omp end single
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
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
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
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
!           Add a MPI Barrier
!           -----------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif

      end subroutine ComputeLaplacianNeumannBCs
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_VolumetricContribution( e , t )
         use HexMeshClass
         use ElementClass
         implicit none
         type(Element)      :: e
         real(kind=RP)      :: t

!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: contravariantFlux  ( 1:NCOMP, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute contravariant flux
!        --------------------------
         call EllipticDiscretization  % ComputeInnerFluxes (NCOMP, NCOMP, e , CHDivergenceFlux3D, contravariantFlux  ) 
!
!        ************************
!        Perform volume integrals
!        ************************
!
         e % storage % QDot = - ScalarWeakIntegrals % StdVolumeGreen ( e , NCOMP, contravariantFlux ) 

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

         e % storage % QDot = e % storage % QDot + ScalarWeakIntegrals % StdFace(e, NCOMP, &
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
      subroutine computeElementInterfaceFlux(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: flux(1:NCOMP,0:f % Nf(1),0:f % Nf(2))

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL EllipticDiscretization % RiemannSolver(nEqn = NCOMP, nGradEqn = NCOMP, &
                                                  f = f, &
                                                  EllipticFlux = CHDivergenceFlux0D, &
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
                                                  flux  = flux(:,i,j) )

               flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

            end do
         end do
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         call f % ProjectFluxToElements(NCOMP, flux, (/1,2/))

      end subroutine computeElementInterfaceFlux

      subroutine computeMPIFaceFlux(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: flux(1:NCOMP,0:f % Nf(1),0:f % Nf(2))

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL EllipticDiscretization % RiemannSolver(nEqn = NCOMP, nGradEqn = NCOMP, &
                                                  f = f, &
                                                  EllipticFlux = CHDivergenceFlux0D, &
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
                                                  flux  = flux(:,i,j) )

               flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(NCOMP, flux, (/thisSide, HMESH_NONE/))

      end subroutine ComputeMPIFaceFlux

      subroutine computeBoundaryFlux(f, time, externalStateProcedure , externalGradientsProcedure)
      USE ElementClass
      use FaceClass
      USE EllipticDiscretizations
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
      procedure(BCState_FCN)     :: externalState
      procedure(BCGradients_FCN) :: externalGradients
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: UGradExt(NDIM , NCOMP) 
      real(kind=RP)                   :: flux(NCOMP, 0:f % Nf(1), 0:f % Nf(2))
      
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

      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         UGradExt(IX,:) = f % storage(1) % U_x(:,i,j)
         UGradExt(IY,:) = f % storage(1) % U_y(:,i,j)
         UGradExt(IZ,:) = f % storage(1) % U_z(:,i,j)

         CALL externalGradientsProcedure(  f % geom % x(:,i,j), &
                                           time, &
                                           f % geom % normal(:,i,j), &
                                           UGradExt,&
                                           boundaryType )

         f % storage(1) % U_x(:,i,j) = UGradExt(IX,:)
         f % storage(1) % U_y(:,i,j) = UGradExt(IY,:)
         f % storage(1) % U_z(:,i,j) = UGradExt(IZ,:)

         f % storage(2) % U_x(:,i,j) = UGradExt(IX,:)
         f % storage(2) % U_y(:,i,j) = UGradExt(IY,:)
         f % storage(2) % U_z(:,i,j) = UGradExt(IZ,:)

!   
!           --------------
!           Viscous fluxes
!           --------------
!   
         CALL EllipticDiscretization % RiemannSolver(nEqn = NCOMP, nGradEqn = NCOMP, &
                                            f = f, &
                                            EllipticFlux = CHDivergenceFlux0D, &
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
                                            flux  = flux(:,i,j) )

         flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

      end do               ; end do

      call f % ProjectFluxToElements(NCOMP, flux, (/1, HMESH_NONE/))

      end subroutine computeBoundaryFlux
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_ComputeGradient( mesh , time , externalStateProcedure)
         implicit none
         type(HexMesh)                  :: mesh
         real(kind=RP),      intent(in) :: time
         procedure(BCState_FCN)         :: externalStateProcedure

         call EllipticDiscretization % ComputeGradient( NCOMP, NCOMP, mesh , time , externalStateProcedure, CHGradientValuesForQ_0D, CHGradientValuesForQ_3D)

      end subroutine DGSpatial_ComputeGradient
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           GET POISEUILLE FLOW
!           -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine GetPoiseuilleFlow_Element(e)
!
!        ****************************************************************************
!              To ensure a continuous flow (OK, impossible in NS, but at least
!           in the advection eqn) we first interpolate the poiseuille flow in
!           Chebyshev-Gauss-Lobatto points, and then traduce the result
!           to Gauss-(Lobatto) points.
!        ****************************************************************************
!
         implicit none
         type(Element)  :: e
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k, n, m, l
         real(kind=RP)  :: xCGL(1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)  :: vCGL(1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
!!
!!        Recover CGL coordinates
!!        -----------------------
!         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!            xCGL(:,i,j,k) = e % hexMap % transfiniteMapAt([e % spAxi % xCGL(i), e % spAeta % xCGL(j), &
!                                                           e % spAzeta % xCGL(k)])
!         end do         ; end do          ; end do
!!
!!        Compute the velocity in CGL points
!!        ----------------------------------
!         do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!            call PoiseuilleFlow(xCGL(:,i,j,k), vCGL(:,i,j,k))
!         end do                  ; end do                ; end do
!!
!!        Return to Gauss points
!!        ----------------------
!         do k = 0, e % Nxyz(3)  ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
!            do n = 0, e % Nxyz(3) ; do m = 0, e % Nxyz(2) ; do l = 0, e % Nxyz(1)
!               e % storage % v(:,i,j,k) = e % storage % v(:,i,j,k) + vCGL(:,l,m,n) & 
!                                          * e % spAxi   % TCheb2Gauss(i,l) &
!                                          * e % spAeta  % TCheb2Gauss(j,m) &
!                                          * e % spAzeta % TCheb2Gauss(k,n)
!            end do              ; end do              ; end do
!         end do               ; end do               ; end do

         do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call PoiseuilleFlow(e % geom % x(:,i,j,k), e % storage % v(:,i,j,k))
         end do                  ; end do                ; end do
         
      end subroutine GetPoiseuilleFlow_Element

      subroutine GetPoiseuilleFlow_Face(f)
         implicit none
         type(Face)  :: f
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j
         logical, save :: shown = .false.

         if ( .not. shown ) then
            write(STD_OUT,'(A)' )"Poiseuille flow computed in face with Gauss points directly."
            shown = .true.
         end if
!
!        Compute the velocity    TODO: compute it in Chebyshev points first 
!        --------------------
         do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
            call PoiseuilleFlow(f % geom % x(:,i,j), f % storage(1) % v(:,i,j))
            call PoiseuilleFlow(f % geom % x(:,i,j), f % storage(2) % v(:,i,j))
         end do                ; end do

      end subroutine GetPoiseuilleFlow_Face

end module SpatialDiscretization
