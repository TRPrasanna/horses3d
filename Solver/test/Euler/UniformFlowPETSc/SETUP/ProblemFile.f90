!
!//////////////////////////////////////////////////////
!
!   @File:    ProblemFile.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Sat May 12 20:54:09 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
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
         SUBROUTINE UserDefinedFinalSetup(mesh, thermodynamics_, &
                                                 dimensionless_, &
                                                     refValues_ )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            use SMConstants
            use HexMeshClass
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
!              - By default it sets a uniform initial
!                 condition.
!           ------------------------------------------------
!
            USE SMConstants
            use PhysicsStorage
            use HexMeshClass
            implicit none
            class(HexMesh)                      :: mesh
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

#if defined(NAVIERSTOKES)
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
!
!              -------------------------------------------------
!              Perturb mean flow in the expectation that it will
!              relax back to the mean flow
!              -------------------------------------------------
!
               mesh % elements(eID) % storage % Q(1,3,3,3) = 1.05_RP*mesh % elements(eID) % storage % Q(1,3,3,3)

            end do

            end associate
#endif
            
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
         subroutine UserDefinedSourceTerm(x, time, S, thermodynamics_, dimensionless_, refValues_)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            IMPLICIT NONE
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(out) :: S(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
!
!           Usage example
!           -------------
!           S(:) = x(1) + x(2) + x(3) + time
   
         end subroutine UserDefinedSourceTerm
!
!//////////////////////////////////////////////////////////////////////// 
! 

         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, monitors)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)  :: mesh
            REAL(KIND=RP) :: time
            type(Monitor_t),  intent(in)  :: monitors
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual, thermodynamics_, &
                                                    dimensionless_, &
                                                        refValues_, &
                                                          monitors, &
                                                         elapsedTime, &
                                                            CPUTime   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            use SMConstants
            USE FTAssertions
            use PhysicsStorage
            use HexMeshClass
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer,                   intent(in) :: iter
            real(kind=RP),             intent(in) :: maxResidual
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
            type(Monitor_t),           intent(in) :: monitors
            real(kind=RP),          intent(in) :: elapsedTime
            real(kind=RP),          intent(in) :: CPUTime
!
!           ---------------
!           Local variables
!           ---------------
!
            INTEGER                            :: numberOfFailures
            CHARACTER(LEN=29)                  :: testName           = "27 element uniform flow tests"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            real(kind=RP)                      :: qq, u, v, w, p, Q(N_EQN), theta, phi
            TYPE(FTAssertionsManager), POINTER :: sharedManager
!
!           -----------------------------------------------------------------------
!           Expected Values. Note they will change if the run parameters change and
!           when the eigenvalue computation for the time step is fixed. These 
!           results are for the Mach 0.5 and rusanov solvers.
!           -----------------------------------------------------------------------
!
#if defined(NAVIERSTOKES)
            INTEGER                            :: expectedIterations = 8
            REAL(KIND=RP)                      :: expectedResidual   = 1.0700773600547120E-011_RP
            
            N = mesh % elements(1) % Nxyz(1) ! This works only because this is an isotropic case.
            
            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue= expectedIterations, &
                               actualValue   =  iter, &
                               msg           = "Number of time steps to tolerance")
            CALL FTAssertEqual(expectedValue = expectedResidual, &
                               actualValue   = maxResidual, &
                               tol           = 1.d-16, &
                               msg           = "Final maximum residual")
            
            ALLOCATE(QExpected(N_EQN,0:N,0:N,0:N))
            
            maxError = 0.0_RP
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

                  QExpected(:,i,j,k) = Q 
               end do;        end do;        end do
               end associate
               maxError = MAXVAL(ABS(QExpected - mesh % elements(eID) % storage % Q))
            end do
            end associate

            CALL FTAssertEqual(expectedValue = 0.0_RP, &
                               actualValue   = maxError, &
                               tol           = 1.d-10, &
                               msg           = "Maximum error")
            
            
            CALL sharedManager % summarizeAssertions(title = testName,iUnit = 6)
   
            IF ( sharedManager % numberOfAssertionFailures() == 0 )     THEN
               WRITE(6,*) testName, " ... Passed"
            ELSE
               WRITE(6,*) testName, " ... Failed"
               WRITE(6,*) "NOTE: Failure is expected when the max eigenvalue procedure is fixed."
               WRITE(6,*) "      When that is done, re-compute the expected values and modify this procedure"
               STOP 99
            END IF 
            WRITE(6,*)
            
            CALL finalizeSharedAssertionsManager
            CALL detachSharedAssertionsManager
#endif

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
      