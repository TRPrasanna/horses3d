!
! //////////////////////////////////////////////////////////////////////////////
!
!
!     FacePatchClass.F
!
!!
!!     Modification History:
!!        version 0.0 July 30, 2002
!!        Updated for NSLite 5/19/15, 10:51 AM
!!
!!     self class stores data needed to define a 2D interpolant. In
!!     self context, self means an iterpolant that defines a surface.
!
!      TYPE FacePatch
!
!      PUBLIC METHODS:
!         SUBROUTINE ConstructFacePatch         ( self, points, uKnots, vKnots )
!         SUBROUTINE DestructFacePatch          ( self )
!         SUBROUTINE ComputeFacePoint           ( self, u, p )
!         SUBROUTINE ComputeFaceDerivative      ( self, u, grad )
!         SUBROUTINE PrintFacePatch             ( self )
!         LOGICAL FUNCTION FaceIs4CorneredQuad  ( self )
!
!!    @author David A. Kopriva
!!    
! //////////////////////////////////////////////////////////////////////////////
!
!  ******
   MODULE FacePatchClass
!  ******
!
     USE SMConstants
     IMPLICIT NONE
     PRIVATE
!
!    ---------------
!    Type definition
!    ---------------
!
!!   Stores the data needed to specify a surface through a 2D
!!   interpolant. The surface is a function of (u,v) in [-1,1]x[-1,1].
!    Arrays are dimensioned 1:nKnots.
!
     TYPE FacePatch
         REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: points
         REAL(KIND=RP), DIMENSION(:)    , ALLOCATABLE :: uKnots,vKnots
         REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: divTable
         INTEGER      , DIMENSION(2)                  :: noOfKnots
!
!        ========         
         CONTAINS 
!        ========         
!
         PROCEDURE :: construct => ConstructFacePatch
         PROCEDURE :: destruct  => DestructFacePatch
         PROCEDURE :: setFacePoints
     END TYPE FacePatch
     
     PUBLIC:: FacePatch
     PUBLIC:: ConstructFacePatch, DestructFacePatch
     PUBLIC:: ComputeFacePoint  , ComputeFaceDerivative
     PUBLIC:: PrintFacePatch    , FaceIs4CorneredQuad
!
!    ========
     CONTAINS
!    ========
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    The constructor for a 2D interpolant takes an array of points and
!!    an array of knots to define a sursurface.
!     -----------------------------------------------------------------
!
      SUBROUTINE ConstructFacePatch( self, uKnots, vKnots, points )
!
!     ---------
!     Arguments
!     ---------
!
      CLASS(FacePatch)                          :: self
      REAL(KIND=RP), DIMENSION(:,:,:), OPTIONAL :: points
      REAL(KIND=RP), DIMENSION(:)               :: uKnots,vKnots
!
!     ---------------
!     Local Variables
!     ---------------
!

!
!     ---------------------------------------------------------
!     The dimensions of the interpolant are given by the number 
!     of points in each  direction. Allocate memory to store
!     self data in.
!     ---------------------------------------------------------
!
      self % noOfKnots(1) = SIZE(uKnots)
      self % noOfKnots(2) = SIZE(vKnots)
!      
      ALLOCATE( self % points(3, self % noOfKnots(1), self % noOfKnots(2)) )
      ALLOCATE( self % uKnots(self % noOfKnots(1)) )
      ALLOCATE( self % vKnots(self % noOfKnots(2)) )
!
!     -------------------------
!     Save the points and knots
!     -------------------------
!
      self % uKnots  = uKnots
      self % vKnots  = vKnots

      IF(PRESENT(points)) CALL self % setFacePoints( points )
!
      END SUBROUTINE ConstructFacePatch
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    Deletes allocated memory
!     -----------------------------------------------------------------
!
      SUBROUTINE DestructFacePatch(self)
!
!     ---------
!     Arguments
!     ---------
!
      CLASS(FacePatch) :: self
      
      IF ( ALLOCATED( self % points ) )   DEALLOCATE( self % points )
      IF ( ALLOCATED( self % uKnots ) )   DEALLOCATE( self % uKnots )
      IF ( ALLOCATED( self % vKnots ) )   DEALLOCATE( self % vKnots )
      IF ( ALLOCATED( self % divTable ) ) DEALLOCATE( self % divTable )
      self % noOfKnots = 0
!
      END SUBROUTINE DestructFacePatch
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE setFacePoints(self,points)  
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(FacePatch)                :: self
         REAL(KIND=RP), DIMENSION(:,:,:) :: points
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: n, k, j
!
         self % points  = points
         
         IF ( .NOT. FaceIs4CorneredQuad(self) )     THEN
            IF(ALLOCATED(self % divTable)) DEALLOCATE( self % divTable)
            ALLOCATE( self % divTable(self % noOfKnots(1), self % noOfKnots(2), 3) )
!
!           -------------------------------------
!           Compute the divided difference tables
!           -------------------------------------
!
            DO j = 1, self % noOfKnots(2)
               DO k = 1,3
                  DO n = 1, self % noOfKnots(1)
                      self % divTable(n,j,k) = points(k,n,j)
                  END DO
                  CALL divdif(self % noOfKnots(1), self % uKnots, self % divTable(:,j,k))
               END DO
            END DO
         END IF 
         
      END SUBROUTINE setFacePoints
!
!     //////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!
!!     ComputeFacePoint: Compute the nodes on a surface
!!     by interpolating
!!
!!>
!!     u       = computational space variable (u, Zeta)
!!     p        = resultant physical space variable (x, y, z)
!!     self     = Interpolation data that defines a surface
!!<
!!
!     --------------------------------------------------------------------------
!
      SUBROUTINE ComputeFacePoint(self, u, p)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FacePatch), INTENT(IN)                 :: self
      REAL(KIND=RP)  , DIMENSION(2), INTENT(IN)   :: u
      REAL(KIND=RP)  , DIMENSION(3), INTENT(OUT)  :: p
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: j
!
!     ------------------------------------------------------------------
!     Use bi-linear mapping if self is a flat side (for speed) otherwise
!     do general Lagrange interpolation
!     ------------------------------------------------------------------
!
      IF( FaceIs4CorneredQuad(self))     THEN
         DO j = 1,3
            p(j)  = self % points(j,1,1)*(1._RP - u(1))*(1._RP - u(2)) &
                  + self % points(j,2,1)*(1._RP + u(1))*(1._RP - u(2)) &
                  + self % points(j,2,2)*(1._RP + u(1))*(1._RP + u(2))  &
                  + self % points(j,1,2)*(1._RP - u(1))*(1._RP + u(2))
         END DO
         p = 0.25_RP*p
      ELSE
         CALL ComputePoly2D(self, u, p)
      END IF
      
      RETURN
      END SUBROUTINE ComputeFacePoint
!
!///////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!
!!     ComputePoly2D: Compute the lagrange interpolant through the points 
!!     p(i, j) at the location (x, y)
!!
!!>
!!     u   = computational space variable (u, Zeta)
!!     p    = resultant physical space variable (x, y, z)
!!     self = Interpolation data that defines a surface
!!<
!!
!     --------------------------------------------------------------------------
!
      SUBROUTINE ComputePoly2D(self, u, p)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FacePatch)             :: self
      REAL(KIND=RP), DIMENSION(2) :: u
      REAL(KIND=RP), DIMENSION(3) :: p
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: w
      REAL(KIND=RP)                              :: l_j
      INTEGER                                    :: j, k
!
      ALLOCATE( w(3,self % noOfKnots(2)) )
!
      DO j = 1, self % noOfKnots(2)
         DO k = 1,3
            w(k,j) = EvaluateNewtonPolynomial( u(1), self % noOfKnots(1), &
                                               self % uKnots, self % divTable(:,j,k) )
         END DO
      END DO
!
      p = 0.0_RP
      DO j = 1, self % noOfKnots(2)
         l_j = EvaluateLagrangePoly(j,u(2), self % noOfKnots(2),self % vKnots)
         DO k = 1,3
            p(k) = p(k) + w(k,j)*l_j
         END DO
      END DO

      DEALLOCATE (w)
      
      RETURN
      END SUBROUTINE ComputePoly2D
!
!     ///////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!
!!     ComputeSurfaceDerivative: Compute the derivative in the two local  
!!     coordinate directions on a subdomain surface by interpolating
!!     the surface data
!!
!     --------------------------------------------------------------------------
!
      SUBROUTINE ComputeFaceDerivative(self, u, grad)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FacePatch)       :: self
      REAL(KIND=RP) , DIMENSION(2)   :: u
      REAL(KIND=RP) , DIMENSION(3,2) :: grad
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: j
!
!     ------------------------------------------------------------------
!     Use bi-linear mapping of self is a flat side (for speed) otherwise
!     do general lagrange interpolation
!     ------------------------------------------------------------------
!
      IF( FaceIs4CorneredQuad(self))     THEN
         DO j = 1,3
            grad(j,1)  =   -self % points(j,1,1)*(1._RP - u(2)) &
                          + self % points(j,2,1)*(1._RP - u(2)) &
                          + self % points(j,2,2)*(1._RP + u(2)) &
                          - self % points(j,1,2)*(1._RP + u(2))

            grad(j,2)  =  - self % points(j,1,1)*(1._RP - u(1)) &
                          - self % points(j,2,1)*(1._RP + u(1))      &
                          + self % points(j,2,2)*(1._RP + u(1))      &
                          + self % points(j,1,2)*(1._RP - u(1))
         END DO
         grad = 0.25*grad
      ELSE
         CALL Compute2DPolyDeriv(self, u, grad)
      END IF

      RETURN
      END SUBROUTINE ComputeFaceDerivative
!
!     ///////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!
!!     Compute2DPolyDeriv: Compute the gradient of the 
!!                         lagrange interpolant through the points 
!!                         p(i, j) at the location (u, v). The result, grad,
!!                         is (dX/du,dX/dv) where X = (x,y,z)
!!
!     --------------------------------------------------------------------------
!
      SUBROUTINE Compute2DPolyDeriv(self, u, grad)
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FacePatch)      , INTENT(IN)  :: self
      REAL(KIND=RP), DIMENSION(2)   , INTENT(IN)  :: u
      REAL(KIND=RP), DIMENSION(3,2) , INTENT(OUT) :: grad
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE :: l_i
      REAL(KIND=RP)                            :: l_j
      INTEGER                                  :: i, j, k
!
      grad = 0.0_RP
      ALLOCATE( l_i(self % noOfKnots(1)) )
!
!     --------------------------
!     First direction derivative
!     --------------------------
!
      DO i = 1, self % noOfKnots(1)
         l_i(i) = EvaluateLagrangePolyDeriv(i, u(1), self % noOfKnots(1),&
     &                                    self % uKnots)
      END DO
!
      DO j = 1, self % noOfKnots(2)
         l_j = EvaluateLagrangePoly(j, u(2), self % noOfKnots(2),&
     &                              self % vKnots)
         DO i = 1, self % noOfKnots(1)
            DO k = 1,3
               grad(k,1) = grad(k,1) + self % points(k, i, j)*l_i(i)*l_j
            END DO
         END DO
      END DO
!
!     ---------------------------
!     Second direction derivative
!     ---------------------------
!
      DO i = 1, self % noOfKnots(1)
         l_i(i) = EvaluateLagrangePoly(i, u(1), self % noOfKnots(1),&
     &                                 self % uKnots)
      END DO
!
      DO j = 1, self % noOfKnots(2)
         l_j = EvaluateLagrangePolyDeriv(j, u(2), self % noOfKnots(2), &
     &                                 self % vKnots)
         DO i = 1, self % noOfKnots(1)
            DO k = 1,3
               grad(k,2) = grad(k,2) + self % points(k,i,j)*l_i(i)*l_j
            END DO
         END DO
      END DO

      DEALLOCATE (l_i )
      
      RETURN
      END SUBROUTINE Compute2DPolyDeriv
!
!     /////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!!    For debugging purposes, print the information in the surface
!!    interpolant.
!     -----------------------------------------------------------------
!
      SUBROUTINE PrintFacePatch( self )
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(FacePatch)        :: self
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: i,j
!
      PRINT *, "-------------Surface Interpolant--------------"
      DO j = 1, self % noOfKnots(2)
         DO i = 1, self % noOfKnots(1)
            PRINT *, i,j, self % uKnots(i), self % vKnots(j), self % points(:,i,j)
         END DO
      END DO
!
      END SUBROUTINE PrintFacePatch
!
!////////////////////////////////////////////////////////////////////////
!
!
!-------------------------------------------------------------------------
! FaceIs4CorneredQuad returns .TRUE. if the face is represented only by
! the four corners.
!-------------------------------------------------------------------------
!
      LOGICAL FUNCTION FaceIs4CorneredQuad( self )
      IMPLICIT NONE
!
!      ---------
!      Arguments
!      ---------
!
       TYPE(FacePatch)        :: self
!
!      ---------------
!      Local Variables
!      ---------------
!
      IF ( (self % noOfKnots(1) == 2) .AND. (self % noOfKnots(2) == 2) )     THEN
         FaceIs4CorneredQuad = .TRUE.
      ELSE
         FaceIs4CorneredQuad = .FALSE.
      END IF 
      
      END FUNCTION FaceIs4CorneredQuad
!@mark -
!
! /////////////////////////////////////////////////////////////////////
!
!!    Compute the newton divided difference table.
!!    The node values are destroyed while creating the table.
!
! /////////////////////////////////////////////////////////////////////
!
      subroutine divdif( nc, x, table )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                          :: nc
      REAL(KIND=RP), DIMENSION(nc) :: x, table
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                          :: i, j
!
      do 100 i = 2,nc
         do 100 j = nc,i,-1
            table(j) = (table(j) - table(j-1))/(x(j) - x(j-i+1))
 100  continue
!
      return
      END subroutine divdif
!
! /////////////////////////////////////////////////////////////////////
!
!!    Compute the newton form interpolant at the point x
!
! /////////////////////////////////////////////////////////////////////
!
      FUNCTION EvaluateNewtonPolynomial( x, nc, xvals, table )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                          :: nc
      REAL(KIND=RP), DIMENSION(nc) :: xVals, table
      REAL(KIND=RP)                :: x
      REAL(KIND=RP)                :: EvaluateNewtonPolynomial
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP)                :: w, p
      INTEGER                          :: k
      
      p = table(1)
      w = 1.0_RP
      do 100 k = 2,nc
         w = (x - xVals(k-1))*w
         p = p + table(k)*w
 100  continue
      EvaluateNewtonPolynomial = p
!
      END FUNCTION EvaluateNewtonPolynomial
!
! /////////////////////////////////////////////////////////////////////
!
!!    Compute the derivative of the newton form interpolant at 
!!    the point x
!
! /////////////////////////////////////////////////////////////////////
!
      FUNCTION EvaluateNewtonPolynomialDeriv( x, nc, xvals, table )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                      :: nc
      REAL(KIND=RP), DIMENSION(nc) :: xVals, table
      REAL(KIND=RP)                :: x
      REAL(KIND=RP)                :: EvaluateNewtonPolynomialDeriv
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(nc) :: s
      REAL(KIND=RP)                :: w, pp, z
      INTEGER                          :: i,j
      
      pp = table(2)
      EvaluateNewtonPolynomialDeriv = pp
      IF( nc == 2 )     RETURN
      s(1) = x - xVals(1)
      s(2) = x - xVals(2)
      pp = pp + (s(1) + s(2))*table(3)
      EvaluateNewtonPolynomialDeriv = pp
      IF( nc == 3 )     RETURN
      w = s(1)*s(2)
      
      DO i = 4,nc
         z = 0.0_RP
         DO j = 1, i-1
            s(j) = s(j)*(x - xVals(i-1))
            z = z + s(j)
         END DO
         s(i) = w
         z    = z + w
         w    = w*(x - xVals(i-1))
         pp   = pp + z*table(i)
      END DO
      
      EvaluateNewtonPolynomialDeriv = pp
!
      END FUNCTION EvaluateNewtonPolynomialDeriv
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute at x the lagrange polynomial L_k of degree n-1    
!!    whose zeros are given by the z(i) 
!---------------------------------------------------------------------
!
      FUNCTION EvaluateLagrangePoly(k,x,n,z)
!
      IMPLICIT NONE
      REAL(KIND=RP) :: EvaluateLagrangePoly
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER, INTENT(IN)                         :: k !! Which poly
      INTEGER, INTENT(IN)                         :: n !! Order+1
      REAL(KIND=RP), DIMENSION(n), INTENT(IN) :: z !! Nodes
      REAL(KIND=RP)              , INTENT(IN) :: x !! Eval Pt
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: j
      REAL(KIND=RP) :: eLP
      
!                                                                       
      IF(k == 1)     THEN 
         eLP = (x - z(2))/(z(k) - z(2)) 
         DO  j = 3,n 
            eLP = eLP*(x - z(j))/(z(k) - z(j))
         END DO
      ELSE 
         eLP = (x - z(1))/(z(k) - z(1)) 
         DO  j = 2,k-1 
            eLP = eLP*(x - z(j))/(z(k) - z(j))
         END DO
         DO j = k+1,n 
            eLP = eLP*(x - z(j))/(z(k) - z(j))
         END DO
      END IF
      EvaluateLagrangePoly= eLp
!                                                                       
      END FUNCTION EvaluateLagrangePoly
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute at x the derivative of the lagrange polynomial L_k of
!!    degree n-1 whose zeros are given by the z(i)
!---------------------------------------------------------------------
!
      FUNCTION EvaluateLagrangePolyDeriv(k,x,n,z)
!
      IMPLICIT NONE
      REAL(KIND=RP) :: EvaluateLagrangePolyDeriv
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER, INTENT(IN)                         :: k !! Which poly
      INTEGER, INTENT(IN)                         :: n !! Order+1
      REAL(KIND=RP), DIMENSION(n), INTENT(IN) :: z !! Nodes
      REAL(KIND=RP)              , INTENT(IN) :: x !! Eval Pt
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: l,m
      REAL(KIND=RP) :: hp,poly
!                                                                       
      hp = 0.0_RP
      DO l = 1,n 
         if(l == k)     CYCLE
         poly = 1.0_RP
         DO m = 1,n 
            if(m == l)     CYCLE
            if(m == k)     CYCLE 
            poly = poly*(x - z(m))/(z(k) - z(m))
         END DO
         hp = hp + poly/(z(k) - z(l)) 
      END DO
      EvaluateLagrangePolyDeriv = hp 
!                                                                       
      END FUNCTION EvaluateLagrangePolyDeriv

!
!  **********     
   END MODULE FacePatchClass
!  **********     