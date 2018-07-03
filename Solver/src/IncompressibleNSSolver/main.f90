!
!//////////////////////////////////////////////////////
!
!   @File:    main.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Jun 20 18:14:45 2018
!   @Last revision date: Tue Jul  3 17:26:17 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 96905b05f7c99a4dc1a38da8202804d6dfef8cb3
!
!//////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
!      HORSES3DMain.f90
!      Created: May 21, 2015 at 12:56 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
      PROGRAM HORSES3DMainiNS
      
      USE SMConstants
      use FTValueDictionaryClass
      USE PhysicsStorage
      USE SharedBCModule
      USE zoneClass
      USE DGSEMClass
      USE BoundaryConditionFunctions
      USE TimeIntegratorClass
      USE mainKeywordsModule
      USE Headers
      USE pAdaptationClass
      use StopwatchClass
      use MPI_Process_Info
      use SpatialDiscretization
      use pAdaptationClass          , only: GetMeshPolynomialOrders
      use NodalStorageClass
      use FluidData
      use FileReaders               , only: ReadControlFile 
      use FileReadingUtilities      , only: getFileName
      use ProblemFileFunctions
#ifdef _HAS_MPI_
      use mpi
#endif
      
      IMPLICIT NONE
      TYPE( FTValueDictionary)            :: controlVariables
      TYPE( DGSem )                       :: sem
      TYPE( TimeIntegrator_t )            :: timeIntegrator
      LOGICAL                             :: success, saveGradients
      integer                             :: initial_iteration
      INTEGER                             :: ierr
      real(kind=RP)                       :: initial_time
      type(BCFunctions_t)                 :: BCFunctions(1)
      procedure(BCState_FCN)              :: externalStateForBoundaryName_iNS
      procedure(BCGradients_FCN)          :: ExternalGradientForBoundaryName_iNS
      character(len=LINE_LENGTH)          :: solutionFileName
      integer, allocatable                :: Nx(:), Ny(:), Nz(:)
      integer                             :: Nmax
      type(pAdaptation_t)                 :: pAdaptator
      procedure(UserDefinedStartup_f)     :: UserDefinedStartup
      procedure(UserDefinedFinalSetup_f)  :: UserDefinedFinalSetup
      procedure(UserDefinedFinalize_f)    :: UserDefinedFinalize
      procedure(UserDefinedTermination_f) :: UserDefinedTermination

      solver = "incompressible navier-stokes"
!
!     ---------------
!     Initializations
!     ---------------
!
      call MPI_Process % Init
      call CheckIfTheVersionIsRequested
!
!     ----------------------------------------------------------------------------------
!     The main is always compiled, so that __DATE__ and __TIME__ are updated accordingly
!     ----------------------------------------------------------------------------------
!
      if ( MPI_Process % doMPIAction ) then
         CALL Main_Header("HORSES3D High-Order (DG) Spectral Element Parallel Incompressible Navier-Stokes Solver",__DATE__,__TIME__)

      else
         CALL Main_Header("HORSES3D High-Order (DG) Spectral Element Sequential Incompressible Navier-Stokes Solver",__DATE__,__TIME__)

      end if

      CALL controlVariables % initWithSize(16)
      CALL UserDefinedStartup
      CALL ConstructSharedBCModule
      
      CALL ReadControlFile( controlVariables )
      CALL CheckInputIntegrity(controlVariables, success)
      IF(.NOT. success)   ERROR STOP "Control file reading error"
      
!
!     ----------------
!     Set up the DGSEM
!     ----------------
!      
      CALL ConstructPhysicsStorage( controlVariables, success )
      IF(.NOT. success)   ERROR STOP "Physics parameters input error"
      
      ! Initialize manufactured solutions if necessary
      
      call GetMeshPolynomialOrders(controlVariables,Nx,Ny,Nz,Nmax)
      call InitializeNodalStorage(Nmax)
      call pAdaptator % construct (Nx,Ny,Nz,controlVariables)      ! If not requested, the constructor returns doing nothing
      
      BCFunctions(1) % externalState => externalStateForBoundaryName_iNS
      BCFunctions(1) % externalGradients => externalGradientForBoundaryName_iNS

      call sem % construct (  controlVariables  = controlVariables,                                         &
                              BCFunctions = BCFunctions, &
                                 Nx_ = Nx,     Ny_ = Ny,     Nz_ = Nz,                                                 &
                                 success           = success)

      call Initialize_SpaceAndTimeMethods(controlVariables, sem % mesh)
                           
      IF(.NOT. success)   ERROR STOP "Mesh reading error"
      CALL checkBCIntegrity(sem % mesh, success)
      IF(.NOT. success)   ERROR STOP "Boundary condition specification error"
      CALL UserDefinedFinalSetup(sem % mesh, thermodynamics, dimensionless, refValues)
!
!     -------------------------
!     Set the initial condition
!     -------------------------
!
      call sem % SetInitialCondition(controlVariables, initial_iteration, initial_time)
!
!     -----------------------------
!     Construct the time integrator
!     -----------------------------
!
      CALL timeIntegrator % construct (controlVariables, initial_iteration, initial_time)
!
!     -----------------
!     Integrate in time
!     -----------------
!
      CALL timeIntegrator % integrate(sem, controlVariables, sem % monitors, pAdaptator, ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
!
!     ----------------------------------
!     Export particles to VTK (temporal)
!     ----------------------------------
!TODO
!      call sem % particles % ExportToVTK()
!
!     --------------------------
!     Show simulation statistics
!     --------------------------
!
      call DisplaySimulationStatistics(sem % numberOftimeSteps, sem % mesh)
!
!     -----------------------------------------------------
!     Let the user perform actions on the computed solution
!     -----------------------------------------------------
!
      CALL UserDefinedFinalize(sem % mesh, timeIntegrator % time, sem % numberOfTimeSteps, &
                              sem % maxResidual, thermodynamics, dimensionless, refValues, &
                              sem % monitors, Stopwatch % ElapsedTime("Solver"), &
                              Stopwatch % CPUTime("Solver"))
#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         call mpi_barrier(MPI_COMM_WORLD, ierr)
      end if
#endif
!
!     -------------------------------------
!     Save the results to the solution file
!     -------------------------------------
!
      IF(controlVariables % stringValueForKey(solutionFileNameKey,LINE_LENGTH) /= "none")     THEN 
         solutionFileName = trim(getFileName(controlVariables % stringValueForKey(solutionFileNameKey,LINE_LENGTH))) // ".hsol"
         saveGradients    = controlVariables % logicalValueForKey(saveGradientsToSolutionKey)
         CALL sem % mesh % SaveSolution(sem % numberOfTimeSteps, timeIntegrator % time, solutionFileName, saveGradients)
      END IF
!
!     ---------
!     Finish up
!     ---------
!
      if (pAdaptator % Constructed) call pAdaptator % destruct()
      CALL timeIntegrator % destruct()
      CALL sem % destruct()
      call Finalize_SpaceAndTimeMethods
      call DestructGlobalNodalStorage()
      CALL destructSharedBCModule
      
      CALL UserDefinedTermination

      call MPI_Process % Close
      
      END PROGRAM HORSES3DMainiNS
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckBCIntegrity(mesh, success)
!
         use SMConstants
         use MeshTypes
         USE HexMeshClass
         use FTValueDictionaryClass
         USE SharedBCModule
         USE BoundaryConditionFunctions_iNS, ONLY:implementediNSBCNames
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh) :: mesh
         LOGICAL       :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                              :: i, j
         INTEGER                              :: faceID, eId
         CHARACTER(LEN=BC_STRING_LENGTH)      :: bcName, namedBC
         CHARACTER(LEN=BC_STRING_LENGTH)      :: bcType
         real(kind=RP)                        :: bcValue
         TYPE(FTMutableObjectArray), POINTER :: bcObjects
         CLASS(FTValue)             , POINTER :: v
         CLASS(FTObject), POINTER             :: obj
         
         success = .TRUE.
!
!        ----------------------------------------------------------
!        Check to make sure that the boundaries defined in the mesh
!        have an associated name in the control file.
!        ----------------------------------------------------------
         
         DO eID = 1, SIZE( mesh % elements )
            DO faceID = 1, 6
               namedBC = mesh % elements(eId) % boundaryName(faceID)
               IF( namedBC == emptyBCName ) CYCLE
               
               bcName = bcTypeDictionary % stringValueForKey(key             = namedBC,         &
                                                             requestedLength = BC_STRING_LENGTH)
               IF ( LEN_TRIM(bcName) == 0 )     THEN
                  PRINT *, "Control file does not define a boundary condition for boundary name = ", &
                            mesh % elements(eId) % boundaryName(faceID)
                  success = .FALSE.
                  return 
               END IF 
            END DO   
         END DO
!
!        --------------------------------------------------------------------------
!        Check that the boundary conditions to be applied are implemented
!        in the code. Keep those updated in the boundary condition functions module
!        --------------------------------------------------------------------------
!
         bcObjects => bcTypeDictionary % allObjects()
         DO j = 1, bcObjects % COUNT()
            obj => bcObjects % objectAtIndex(j)
            CALL castToValue(obj,v)
            bcType = v % stringValue(requestedLength = BC_STRING_LENGTH)
            DO i = 1, SIZE(implementediNSBCNames)
               IF ( bcType == implementediNSBCNames(i) )     THEN
                  success = .TRUE. 
                  EXIT 
               ELSE 
                  success = .FALSE. 
               END IF 
            END DO
            
            IF ( .NOT. success )     THEN
               PRINT *, "Boundary condition ", TRIM(bcType)," not implemented in this code"
               CALL release(bcObjects)
               RETURN 
            END IF  
            
         END DO
         
         CALL release(bcObjects)
         
      END SUBROUTINE checkBCIntegrity
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckInputIntegrity( controlVariables, success )  
         use SMConstants
         use Utilities, only: toLower
         USE FTValueDictionaryClass
         USE mainKeywordsModule
         use FTValueClass
         use MPI_Process_Info
         use SpatialDiscretization, only: viscousDiscretizationKey
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
         character(len=LINE_LENGTH)    :: inviscidDiscretization, discretizationNodes
         
         success = .TRUE.
!
!        Control variables with default value
!        ------------------------------------
         obj => controlVariables % objectForKey(saveGradientsToSolutionKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey(".false.",saveGradientsToSolutionKey)
         end if

         obj => controlVariables % objectForKey(discretizationNodesKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey("Gauss",discretizationNodesKey)
         end if

         obj => controlVariables % objectForKey(inviscidDiscretizationKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey("Standard",inviscidDiscretizationKey)
         end if

         obj => controlVariables % objectForKey(viscousDiscretizationKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey("BR1",viscousDiscretizationKey)
         end if

         obj => controlVariables % objectForKey(splitFormKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey("Skew-symmetric",splitFormKey)
         end if
!
!        Check for inconsistencies in the input variables
!        ------------------------------------------------
         inviscidDiscretization = trim(controlVariables % stringValueForKey(inviscidDiscretizationKey, LINE_LENGTH))
         discretizationNodes = trim(controlVariables % stringValueForKey(discretizationNodesKey, LINE_LENGTH))

         call toLower(inviscidDiscretization)
         call toLower(discretizationNodes)

         if ( (trim(inviscidDiscretization) .eq. "split-form") .and. (trim(discretizationNodes) .eq. "gauss") ) then
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(A)') "*** WARNING:    Only Gauss-Lobatto nodes are available for Split-Form discretizations"
               write(STD_OUT,'(A)') "*** WARNING:    Automatically switched to Gauss-Lobatto points"
            end if
            call controlVariables % removeObjectForKey(discretizationNodesKey)
            call controlVariables % addValueForKey("Gauss-Lobatto",discretizationNodesKey)
         end if
!
!        Check the controlVariables created
!        ----------------------------------        
         DO i = 1, SIZE(mainKeywords)
            obj => controlVariables % objectForKey(mainKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",mainKeywords(i)
               success = .FALSE. 
            END IF  
         END DO  
         
      END SUBROUTINE checkInputIntegrity

      subroutine DisplaySimulationStatistics(iter,mesh)
         use SMConstants
         use HexMeshClass
         use StopwatchClass
         use Headers
         use MPI_Process_Info
#ifdef _HAS_MPI_
         use mpi
#endif
         implicit none
         integer,    intent(in)      :: iter
         type(HexMesh),   intent(in) :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: eID
         integer                    :: NDOF, localNDOF, ierr
         real(kind=RP)              :: Naverage, localNaverage
         real(kind=RP)              :: t_elaps, t_cpu
   
         if ( MPI_Process % isRoot ) write(STD_OUT,'(/)')
         call Section_Header("Simulation statistics")
         if ( MPI_Process % isRoot ) write(STD_OUT,'(/)')
!
!        Get mesh-related quantities
!        ---------------------------
         NDOF = 0
         Naverage = 0
   
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            NDOF = NDOF + (e % Nxyz(1) + 1)*(e % Nxyz(2) + 1)*(e % Nxyz(3) + 1)      
            Naverage = Naverage + e % Nxyz(1) + e % Nxyz(2) + e % Nxyz(3)
            end associate
         end do

         Naverage = Naverage / (3.0_RP * mesh % no_of_elements)
!
!        Perform a broadcast for the MPI solver
!        --------------------------------------
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
            localNDOF = NDOF
            localNaverage = Naverage * 3.0_RP * mesh % no_of_elements
            call mpi_allreduce(localNDOF, NDOF, 1, MPI_INT, MPI_SUM, &
                               MPI_COMM_WORLD, ierr)

            call mpi_allreduce(localNaverage, Naverage, 1, MPI_DOUBLE, MPI_SUM, &
                               MPI_COMM_WORLD, ierr)

            Naverage = Naverage / (3.0_RP * mesh % no_of_allElements)

         end if
#endif

         if ( .not. MPI_Process % isRoot ) return
!
!        Show preprocessing time
!        -----------------------
         t_elaps = Stopwatch % Elapsedtime("Preprocessing")
         t_cpu   = Stopwatch % CPUTime("Preprocessing")

         call Subsection_Header("Preprocessing")

         write(STD_OUT,'(30X,A,I0,A,F5.2,A,I0,A)')      "->   ", mesh % no_of_elements, &
                                                      " elements with polynomial order ",Naverage," (NDOF = ",NDOF,")."
         write(STD_OUT,'(30X,A,A30,ES10.3,A,ES10.3,A)') "->", "Preprocessing time: ",t_elaps," seconds (total CPU time: ",t_cpu,")."

!
!        Show simulation time
!        --------------------
         write(STD_OUT,'(/)')
         call Subsection_Header("Solver")
         if ( iter .le. 0 ) return

         t_elaps = Stopwatch % ElapsedTime("Solver")
         t_cpu   = Stopwatch % CPUTime("Solver")

         write(STD_OUT,'(30X,A,A30,ES10.3,A)') "->", "Simulation elapsed time: ",t_elaps," seconds."
         write(STD_OUT,'(30X,A,A30,ES10.3,A,ES10.3,A)') "->", "Simulation CPU time: ",t_cpu," seconds (ratio is ",t_cpu/t_elaps ,")."
         write(STD_OUT,'(30X,A,A30,ES10.3,A)') "->", "Solver efficiency: " , t_elaps/(NDOF * iter)*1.0e6_RP, " seconds/(1 Million DOF·iter)."

      end subroutine DisplaySimulationStatistics

      subroutine CheckIfTheVersionIsRequested
         use SMConstants
         use MPI_Process_Info
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: nArgs, i
         character(len=128)    :: arg

         if ( .not. MPI_Process % isRoot ) return

         nArgs = command_argument_count()

         do i = 1, nArgs
            call get_command_argument(i, arg)

            if ( trim(arg) .eq. "--version" ) then
               print*, "Current HORSES version: ", trim(VERSION)
               stop
            end if
         end do
            
      end subroutine CheckIfTheVersionIsRequested
