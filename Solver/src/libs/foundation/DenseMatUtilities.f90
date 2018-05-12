!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      DenseMatUtilities.f90
!      Created: 2017-05-23 10:06:00 +0100 
!      By: Andrés Rueda
!
!      Module containing useful routines for dense matrices
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE DenseMatUtilities
   USE SMConstants
   IMPLICIT NONE

   private
   public inverse, Mat2File, ComputeLU, SolveLU
!========
 CONTAINS
!========

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeLU(A,ALU,LUpivots)
      implicit none
      !------------------------------------------------------------
      real(kind=RP), dimension(:,:)                , intent(in)  :: A          !<  Matrix to factorize
      real(kind=RP), dimension(size(A,1),size(A,1)), intent(out) :: ALU        !>  Factorized LU matrix
      integer      , dimension(size(A,1))          , intent(out) :: LUpivots   !>  LU pivots for factorization
      !------------------------------------------------------------
      integer :: n      ! Matrix size
      integer :: info   ! Lapack error code
      !------------------------------------------------------------


      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      ALU = A
      n = size(A,1)
      
#ifdef HAS_LAPACK
      call DGETRF(n, n, ALU, n, LUpivots, info)
#else
      call LuDecomp(n, ALU, LUpivots, info)
#endif

      if (info /= 0) then
         print*, 'Matrix is numerically singular. ERROR:', info
         stop
      end if

   end subroutine
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SolveLU(ALU,LUpivots,x,b)
      implicit none
      !------------------------------------------------------------
      real(kind=RP), dimension(:,:)        , intent(in)  :: ALU              !<  Factorized LU matrix
      integer      , dimension(size(ALU,1)), intent(in)  :: LUpivots   !<  LU pivots for factorization
      real(kind=RP), dimension(size(ALU,1)), intent(in)  :: b
      real(kind=RP), dimension(size(ALU,1)), intent(out) :: x
      !------------------------------------------------------------
      integer :: n      ! Matrix size
      integer :: info   ! Lapack error code
      !------------------------------------------------------------
      
      n = size(ALU,1)
      x = b
      
#ifdef HAS_LAPACK
      call dgetrs ( 'N' , n, 1   , ALU, n  , LUpivots, x ,n  , info )
#else
      call LUsolve(n, ALU, LUpivots , x, info)
#endif
      
      if (info /= 0) then
         print*,  '*** System could not be solved. ERROR:', info
         stop
      end if

   end subroutine
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   !------------------------------------------------------------------------------
   ! Returns the inverse of a matrix calculated by LU decomposition.  Depends on LAPACK.
   !   Modified from: http://fortranwiki.org/fortran/show/Matrix+inversion
   function inverse(A) result(Ainv)
   !------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------
      real(KIND=RP), dimension(:,:), intent(in) :: A
      real(KIND=RP), dimension(size(A,1),size(A,2)) :: Ainv
      
      real(KIND=RP), dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info
      !------------------------------------------------------------
      
#ifdef HAS_LAPACK
      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI
      
      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)
      
      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)
      
      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if
      
      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)
      
      if (info /= 0) then
         stop 'Lapack matrix inversion failed!'
      end if
#else
      STOP ':: Matrix inversion routine needs LAPACK.'
#endif
   !------------------------------------------------------------------------------
   end function
   !------------------------------------------------------------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   !---------------------------------------------------------------
   SUBROUTINE Mat2File(Mat,FileName,Formato)
   !       Writes a dense matrix to a File
   !---------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------
      REAL(KIND=RP)               :: Mat(:,:)   !< Matrix
      CHARACTER(len=*)            :: FileName
      CHARACTER(len=*),OPTIONAL   :: Formato
      !------------------------------------------
      INTEGER                     :: i, j    ! Counters
      CHARACTER(len=50)           :: Forma   ! Final format
      INTEGER                     :: fd      ! File unit for writing
      !------------------------------------------
      
      IF (PRESENT(Formato)) THEN
         Forma = Formato
      ELSE
         Forma = '(300000F14.8)' !Case default
      END IF
      
      OPEN(newunit=fd,file=FileName)
      DO i=1, SIZE(Mat,1)
         WRITE(fd,Forma) (Mat(i,j),j=1,SIZE(Mat,2))
      END DO
      ClOSE(fd)
      
   !---------------------------------------------------------------
   END SUBROUTINE Mat2File
   !---------------------------------------------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------------
!  Perform LU decomposition (taken from the DSEM code)
!  ---------------------------------------------------
   SUBROUTINE LuDecomp(dim, A, pivots, info)
!  ----------
!
   IMPLICIT NONE
!  
!  ---------
!  Arguments
!  ---------
!
   INTEGER                              , INTENT(IN)    :: dim
   REAL(KIND=RP), DIMENSION(dim,dim), INTENT(INOUT) :: A
   INTEGER          , DIMENSION(dim)    , INTENT(OUT)   :: pivots
   INTEGER                              , INTENT(OUT)   :: info
!
!  ---------------
!  Local Variables
!  ---------------
!
   INTEGER                 ::  k,i,j,n,kmax,temp
   INTEGER, DIMENSION(dim) :: p
!
   n = dim
   if ( n /= size(A,1) ) then
      print *, "LU Decomposition failed, incorrect dimensions."
      stop
   end if 
!
   do i=1,n
      p(i) = i
   end do
   do i=1,n
      kmax = i
      do k=i+1,n
         if ( abs(A(p(k),i)) > abs(A(p(kmax),i)) ) then
            kmax = k
         end if
      end do
      temp = p(i)
      p(i) = p(kmax)
      p(kmax) = temp
      do j=2,i
         if ( abs(A(p(j-1),j-1)) < 10*tiny(1._RP) ) then
            print *, "LUdecomp failed, stopping."
            stop
         end if
         A(p(i),j-1) = A(p(i),j-1) / A(p(j-1),j-1)
         do k=1,j-1 
            A(p(i),j) = A(p(i),j) - A(p(i),k)*A(p(k),j)
         end do
      end do
      do j=i+1,n 
         do k=1,i-1 
            A(p(i),j) = A(p(i),j) - A(p(i),k)*A(p(k),j)
         end do
      end do
   end do
   pivots = p
   
   info = 0
!  
!  --------------
   END SUBROUTINE LuDecomp
!  --------------
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------------
!  Solve a linear system wuth the LU factorized matrix (taken from the DSEM code)
!  ------------------------------------------------------------------------------
   SUBROUTINE LUsolve(dim, A, pivots , b, info)
!  ----------
!
!  ---------
!  Arguments
!  ---------
!
   INTEGER                          , INTENT(IN)     ::  dim
   REAL(KIND=RP), DIMENSION(:,:), INTENT(IN)     ::  A
   REAL(KIND=RP), DIMENSION(:)  , INTENT(INOUT)  ::  b
   INTEGER          , DIMENSION(:)  , INTENT(IN)     ::  pivots
   INTEGER                          , INTENT(OUT)    ::  info
!
!  ---------------
!  Local Variables
!  ---------------
!
   REAL(KIND=RP), DIMENSION(SIZE(b)) ::  bTemp
   INTEGER                               ::  i
!
   bTemp = b
!
   call LUFsolve(A,pivots,bTemp)
   call LUBsolve(A,pivots,bTemp)
!
   do i=1,size(b)
      b(i) = bTemp(pivots(i))
   end do
   
   info = 0
!
!  --------------
   END SUBROUTINE LUsolve
!  --------------
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------------
!  Solve for B (taken from the DSEM code)
!  ---------------------------------------------------
   subroutine LUBsolve(A,p,b)
!  ----------
!
   real(KIND=RP), dimension(:,:), intent(in)     ::  A
   integer , dimension(:),   intent(in)              ::  p
   real(KIND=RP), dimension(:),   intent(inout)  ::  b
!
   integer  ::  j,k,n
!
   n = size(A(1,:))
   if ( n /= size(A(:,1)) .or. n /= size(b) .or. n /= size(p) ) then
      print *, "LUBsolve failed, incorrect dimensions."
      stop
   end if
!
   do k=1,n
      if ( abs(A(p(k),k)) < 10*tiny(1._RP) ) then
         print *, "LUBsolve failed, check LU decomposition; stopping"
         stop
      end if
   end do
!
   do k=n,1,-1
      b(p(k)) = b(p(k))/A(p(k),k)
      do j=1,k-1
         b(p(j)) = b(p(j)) - A(p(j),k)*b(p(k))
      end do
   end do
!
!  --------------
   end subroutine LUBsolve
!  --------------
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------
!  Solve for F (taken from the DSEM code)
!  --------------------------------------
   subroutine LUFSolve(A,p,b)
!  ----------
!
   real(KIND=RP), dimension(:,:), intent(in)     ::  A
   integer , dimension(:),   intent(in)              ::  p
   real(KIND=RP), dimension(:),   intent(inout)  ::  b
!
   integer  ::  j,k,n
!
   n = size(A(1,:))
   if ( n /= size(A(:,1)) .or. n /= size(b) .or. n /= size(p) ) then
      print *, "LUFsolve failed, incorrect dimensions."
      stop
   end if
!
   do k=1,n
      do j=k+1,n
         b(p(j)) = b(p(j)) - A(p(j),k)*b(p(k))
      end do
   end do
!
!  --------------
   end subroutine LUFsolve
!  --------------
END MODULE DenseMatUtilities
