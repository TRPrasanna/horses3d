!
!////////////////////////////////////////////////////////////////////////
!
!      ProblemFile.f90
!      Created: June 26, 2015 at 8:47 AM 
!      By: David Kopriva  
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedSetUp
!      UserDefinedInitialCondition(mesh)
!      UserDefinedPeriodicOperation(mesh)
!      UserDefinedFinalize(mesh)
!      UserDefinedTermination
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(mesh , thermodynamics_, &
                                                 dimensionless_, &
                                                     refValues_ )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            IMPLICIT NONE
            CLASS(HexMesh)             :: mesh
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedInitialCondition(mesh, thermodynamics_, &
                                                      dimensionless_, &
                                                          refValues_  )
!
!           ------------------------------------------------
!           Called to set the initial condition for the flow
!              - By default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            USE SMConstants
            use PhysicsStorage
            use HexMeshClass
            implicit none
            class(HexMesh)                        :: mesh
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
!
!           ---------------
!           Local variables
!           ---------------
!
            integer        :: eID, i, j, k
            real(kind=RP)  :: qq, u, v, w, p
            real(kind=RP)  :: Q(N_EQN), phi, theta

            associate ( gammaM2 => dimensionless_ % gammaM2, &
                        gamma => thermodynamics_ % gamma )
            theta = refValues_ % AOATheta*(PI/180.0_RP)
            phi   = refValues_ % AOAPhi*(PI/180.0_RP)
      
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          Ny => mesh % elements(eID) % Nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
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

                  mesh % elements(eID) % storage % Q(:,i,j,k) = Q 
               end do;        end do;        end do
               end associate
            end do

            end associate
            
         END SUBROUTINE UserDefinedInitialCondition

         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(N_EQN)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
         end subroutine UserDefinedState1

         subroutine UserDefinedNeumann(x, t, nHat, U_x, U_y, U_z)
!
!           --------------------------------------------------------
!           Used to define a Neumann user defined boundary condition
!           --------------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: U_x(N_GRAD_EQN)
            real(kind=RP), intent(inout)  :: U_y(N_GRAD_EQN)
            real(kind=RP), intent(inout)  :: U_z(N_GRAD_EQN)
         end subroutine UserDefinedNeumann

!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedSourceTerm(mesh, time, thermodynamics_, dimensionless_, refValues_)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
!
!           Usage example (by default no source terms are added)
!           ----------------------------------------------------
!           do eID = 1, mesh % no_of_elements
!              associate ( e => mesh % elements(eID) )
!              do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!                 e % QDot(:,i,j,k) = e % QDot(:,i,j,k) + Source(:)
!              end do                  ; end do                ; end do
!           end do
   
         end subroutine UserDefinedSourceTerm
!
!//////////////////////////////////////////////////////////////////////// 
! 

         SUBROUTINE UserDefinedPeriodicOperation(mesh, time)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            USE HexMeshClass
            IMPLICIT NONE
            CLASS(HexMesh)  :: mesh
            REAL(KIND=RP) :: time
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual, thermodynamics_, &
                                                    dimensionless_, &
                                                        refValues_   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            use FTAssertions
            USE HexMeshClass
            use PhysicsStorage
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
!
!           ---------------
!           Local variables
!           ---------------
!
            CHARACTER(LEN=29)                  :: testName           = "Different orders cylinder"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
!
!           -----------------------------------------------------------------------------------------
!           Expected solutions. 
!           InnerCylinder 0.0 NoSlipAdiabaticWall
!           Front 0.0 Inflow
!           bottom 0.0 FreeSlipWall
!           top 0.0 FreeSlipWall
!           Back 0.0 Inflow
!           Left 0.0 Inflow
!           Right 0.0 OutflowSpecifyP 
!           -----------------------------------------------------------------------------------------
!
!
!           ------------------------------------------------
!           Expected Solutions: Wall conditions on the sides
!           Number of iterations are for CFL of 0.2, for
!           the roe solver and mach = 0.3
!           ------------------------------------------------
!
            INTEGER                            :: iterations = 100
            REAL(KIND=RP)                      :: residuals = 264.64042924505213
!            
            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = iterations, &
                               actualValue   = iter, &
                               msg           = "Number of time steps to tolerance")
            CALL FTAssertEqual(expectedValue = residuals, &
                               actualValue   = maxResidual, &
                               tol           = 1.d-3, &
                               msg           = "Final maximum residual")
            
            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
               WRITE(6,*) "This test case has no expected solution yet, only checks the residual after 100 iterations."
            ELSE
               WRITE(6,*) testName, " ... Failed"
               WRITE(6,*) "NOTE: Failure is expected when the max eigenvalue procedure is changed."
               WRITE(6,*) "      If that is done, re-compute the expected values and modify this procedure"
                STOP 99
            END IF 
            WRITE(6,*)
            
            CALL finalizeSharedAssertionsManager
            CALL detachSharedAssertionsManager
            

         END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      
