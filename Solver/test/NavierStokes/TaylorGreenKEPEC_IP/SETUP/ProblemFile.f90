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
         SUBROUTINE UserDefinedFinalSetup(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            CLASS(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedInitialCondition(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        )
!
!           ------------------------------------------------
!           called to set the initial condition for the flow
!              - by default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            use fluiddata
            implicit none
            class(hexmesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            REAL(KIND=RP) :: x(3)        
            INTEGER       :: i, j, k, eID
            REAL(KIND=RP) :: rho , u , v , w , p
            REAL(KIND=RP) :: L, u_0, rho_0, p_0
            integer       :: Nx, Ny, Nz

#if defined(NAVIERSTOKES)
            
            L     = 1.0_RP
            u_0   = 1.0_RP
            rho_0 = 1.0_RP 
            p_0   = 100.0_RP

            associate( gamma => thermodynamics_ % gamma ) 
            DO eID = 1, SIZE(mesh % elements)
               Nx = mesh % elements(eID) % Nxyz(1)
               Ny = mesh % elements(eID) % Nxyz(2)
               Nz = mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx 

                         x = mesh % elements(eID) % geom % x(:,i,j,k)
                       
                         rho = rho_0
                         u   =  u_0 * sin(x(1)/L) * cos(x(2)/L) * cos(x(3)/L) 
                         v   = -u_0 * cos(x(1)/L) * sin(x(2)/L) * cos(x(3)/L)
                         w   =  0.0_RP
                         p   =   p_0 + rho_0 / 16.0_RP * (                          &
                               cos(2.0_RP*x(1)/L)*cos(2.0_RP*x(3)/L) +                  &
                               2.0_RP*cos(2.0_RP*x(2)/L) + 2.0_RP*cos(2.0_RP*x(1)/L) +  &
                               cos(2.0_RP*x(2)/L)*cos(2.0_RP*x(3)/L)                    &
                               )

                         mesh % elements(eID) % storage % Q(1,i,j,k) = rho
                         mesh % elements(eID) % storage % Q(2,i,j,k) = rho*u
                         mesh % elements(eID) % storage % Q(3,i,j,k) = rho*v
                         mesh % elements(eID) % storage % Q(4,i,j,k) = rho*w
                         mesh % elements(eID) % storage % Q(5,i,j,k) = p / (gamma - 1.0_RP) + 0.5_RP * rho * (u*u + v*v + w*w)

                     END DO
                  END DO
               END DO 
               
            END DO 
            end associate
#endif            
         end subroutine UserDefinedInitialCondition
#if defined(NAVIERSTOKES)
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
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
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: U_x(NGRAD)
            real(kind=RP), intent(inout)  :: U_y(NGRAD)
            real(kind=RP), intent(inout)  :: U_z(NGRAD)
         end subroutine UserDefinedNeumann
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
<<<<<<< HEAD
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, Monitors)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
#if defined(NAVIERSTOKES)
            use MonitorsClass
#endif
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
#if defined(NAVIERSTOKES)
            type(Monitor_t), intent(in) :: monitors
#else
            logical, intent(in) :: monitors
#endif
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedSourceTerm(mesh, time      &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_ &
                                        , dimensionless_  &
                                        , refValues_ & 
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
)
=======
         subroutine UserDefinedSourceTerm(x, time, S, thermodynamics_, dimensionless_, refValues_)
>>>>>>> master
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
<<<<<<< HEAD
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
=======
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(out) :: S(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
>>>>>>> master
!
!           Usage example
!           -------------
!           S(:) = x(1) + x(2) + x(3) + time
   
         end subroutine UserDefinedSourceTerm
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual &
#if defined(NAVIERSTOKES)
                                                    , thermodynamics_ &
                                                    , dimensionless_  &
                                                    , refValues_ & 
#endif   
#if defined(CAHNHILLIARD)
                                                    , multiphase_ &
#endif
                                                    , monitors, &
                                                      elapsedTime, &
                                                      CPUTime   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            use SMConstants
            use FTAssertions
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)    :: thermodynamics_
            type(Dimensionless_t),  intent(in)    :: dimensionless_
            type(RefValues_t),      intent(in)    :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)    :: multiphase_
#endif
            type(Monitor_t),        intent(in)    :: monitors
            real(kind=RP),             intent(in) :: elapsedTime
            real(kind=RP),             intent(in) :: CPUTime
!
!           ---------------
!           Local variables
!           ---------------
!
#if defined(NAVIERSTOKES)
            CHARACTER(LEN=29)                  :: testName           = "Taylor-Green vortex with Kinetic Energy preserving and entrophy conserving split form + IP"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success
            integer                            :: rank
            real(kind=RP), parameter           :: kinEn = 0.12499744094576197_RP
            real(kind=RP), parameter           :: kinEnRate = -4.2793777270020350E-004_RP
            real(kind=RP), parameter           :: enstrophy = 0.37499382750807664_RP 
            real(kind=RP), parameter           :: res(5) = [7.8469872929594112E-005_RP, &
                                                            0.12839595215106728_RP, & 
                                                            0.12839595254623606_RP, & 
                                                            0.24998296113284613_RP, &
                                                            0.61597739698335341_RP ]

            CALL initializeSharedAssertionsManager
            sharedManager => sharedAssertionsManager()
            
            CALL FTAssertEqual(expectedValue = monitors % residuals % values(1,1) + 1.0_RP, &
                               actualValue   = res(1) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "continuity residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(2,1) + 1.0_RP, &
                               actualValue   = res(2) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "x-momentum residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(3,1) + 1.0_RP, &
                               actualValue   = res(3) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "y-momentum residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(4,1) + 1.0_RP, &
                               actualValue   = res(4) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "z-momentum residual")

            CALL FTAssertEqual(expectedValue = monitors % residuals % values(5,1) + 1.0_RP, &
                               actualValue   = res(5) + 1.0_RP, &
                               tol           = 1.0e-7_RP, &
                               msg           = "energy residual")

            CALL FTAssertEqual(expectedValue = monitors % volumeMonitors(1) % values(1), &
                               actualValue   = kinEn, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Kinetic Energy")

            CALL FTAssertEqual(expectedValue = monitors % volumeMonitors(2) % values(1) + 1.0_RP, &
                               actualValue   = kinEnRate + 1.0_RP, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Kinetic Energy Rate")

            CALL FTAssertEqual(expectedValue = monitors % volumeMonitors(3) % values(1), &
                               actualValue   = enstrophy, &
                               tol           = 1.0e-11_RP, &
                               msg           = "Enstrophy")

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
      
