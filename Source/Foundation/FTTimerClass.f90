!
!////////////////////////////////////////////////////////////////////////
!
!      FTTimerClass.f90
!      Created: September 2, 2013 1:58 PM 
!      By: NocturnalAviationSoftware  
!
!      Defines a class for timing Fortran program
!      execution
!      
!      Usage
!
!         * Initialization *
!
!               CALL timer % init()
!
!         * Starting the timer *
!
!               CALL timer % start()
!
!         * Stopping the timer *
!
!               CALL timer % stop()
!
!         * Reading the time *
!
!               time = timer % elapsedTime(units) ! Wall clock time
!               time = timer % totalTime(units)   ! Total processor time
!
!           units (optional) = TC_SECONDS or TC_MINUTES or TC_HOURS
!
!////////////////////////////////////////////////////////////////////////
!
      MODULE FTTimerClass 
      IMPLICIT NONE
      PRIVATE
!
!     -----------------
!     Private constants
!     -----------------
!
      INTEGER, PARAMETER, PRIVATE :: d = 15
!
!     ----------------
!     Public constants
!     ----------------
!
      INTEGER, PARAMETER, PUBLIC  :: TP = SELECTED_REAL_KIND(d)
      INTEGER, PARAMETER, PUBLIC  :: TC_SECONDS = 0, TC_MINUTES = 1, TC_HOURS = 2
!
!     ---------------------
!     Class type definition
!     ---------------------
!
      TYPE, PUBLIC :: FTTimer
         LOGICAL      , PRIVATE :: started    = .FALSE., stopped = .FALSE.
         REAL(KIND=TP), PRIVATE :: startTime  = 0.0_TP
         REAL(KIND=TP), PRIVATE :: finishTime = 0.0_TP
         INTEGER :: &
           initialTicks = 0, &
           finalTicks   = 0, & ! final value of the clock tick counter
           countMax     = 0, & ! maximum value of the clock counter
           countRate    = 0, & ! number of clock ticks per second
           numTicks     = 0    ! number of clock ticks of the code
!
!        ========
         CONTAINS
!        ========
!
         PROCEDURE, PASS :: init  => initTimer
         PROCEDURE, PASS :: start => startTimer
         PROCEDURE, PASS :: stop  => stopTimer
         PROCEDURE, PASS :: elapsedTime
         PROCEDURE, PASS :: totalTime
         
      END TYPE FTTimer
!
!     ========
      CONTAINS
!     ========
! 
!
!////////////////////////////////////////////////////////////////////////  
!  
      SUBROUTINE initTimer(self)  
         IMPLICIT NONE
         CLASS(FTTimer) :: self
         self % started = .FALSE.
         CALL SYSTEM_CLOCK(COUNT_RATE = self % countRate, &
                           COUNT_MAX  = self % countMax)
         
      END SUBROUTINE initTimer
!
!////////////////////////////////////////////////////////////////////////  
!  
      SUBROUTINE startTimer(self)  
         IMPLICIT NONE
         CLASS(FTTimer) :: self
         self % started = .TRUE.
         CALL CPU_TIME(self % startTime)         
         CALL SYSTEM_CLOCK(COUNT = self % initialTicks)
      END SUBROUTINE startTimer
!
!////////////////////////////////////////////////////////////////////////  
!  
      SUBROUTINE stopTimer(self)  
         IMPLICIT NONE
         CLASS(FTTimer) :: self
         self % stopped = .TRUE.
         CALL CPU_TIME(self % finishTime)
         CALL SYSTEM_CLOCK(COUNT = self % finalTicks)
      END SUBROUTINE stopTimer
!
!//////////////////////////////////////////////////////////////////////// 
! 
       FUNCTION totalTime(self,units)  
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FTTimer)    :: self
         INTEGER, OPTIONAL :: units
         REAL(KIND=TP)     :: totalTime
!
!        ------------------------------------------
!        Return zero if the timer was never started
!        ------------------------------------------
!
         IF ( .NOT.self % started )     THEN
            totalTime = 0.0_TP
            RETURN
         END IF 
!
!        ----------------------------------------------
!        If the timer was not stopped, then return the 
!        current time elapsed
!        ----------------------------------------------
!
         IF ( .NOT.self % stopped )     THEN
            CALL self % stop() 
         END IF 

         totalTime =  self % finishTime - self % startTime
!
!        -------------------------------------
!        Convert to requested units if present
!        -------------------------------------
!
         IF ( PRESENT(units) )     THEN
         
            SELECT CASE ( units )
               CASE( TC_MINUTES ) 
                  totalTime = totalTime/60.0_TP
               CASE( TC_HOURS )
                  totalTime = totalTime/3600.0_TP
               CASE DEFAULT 
               
            END SELECT 
         END IF 
      
      END FUNCTION totalTime
!
!//////////////////////////////////////////////////////////////////////// 
! 
       FUNCTION elapsedTime(self,units)  
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FTTimer)    :: self
         INTEGER, OPTIONAL :: units
         REAL(KIND=TP)     :: elapsedTime
!
!        ------------------------------------------
!        Return zero if the timer was never started
!        ------------------------------------------
!
         IF ( .NOT.self % started )     THEN
            elapsedTime = 0.0_TP
            RETURN
         END IF 
!
!        ----------------------------------------------
!        If the timer was not stopped, then return the 
!        current time elapsed
!        ----------------------------------------------
!
         IF ( .NOT.self % stopped )     THEN
            CALL self % stop() 
         END IF 

         self % numTicks = self % finalTicks - self % initialTicks
         IF (self % finalTicks < self % initialTicks)  self % numTicks = self % numTicks + self % countMax
         elapsedTime   = REAL(self % numTicks, KIND=TP) / self % countRate
!
!        -------------------------------------
!        Convert to requested units if present
!        -------------------------------------
!
         IF ( PRESENT(units) )     THEN
         
            SELECT CASE ( units )
               CASE( TC_MINUTES ) 
                  elapsedTime = elapsedTime/60.0_TP
               CASE( TC_HOURS )
                  elapsedTime = elapsedTime/3600.0_TP
               CASE DEFAULT 
               
            END SELECT 
         END IF 
      
      END FUNCTION elapsedTime
      
      END MODULE FTTimerClass
