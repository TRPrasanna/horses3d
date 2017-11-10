!
!////////////////////////////////////////////////////////////////////////
!
!      MappedGeometry.f95
!      Created: 2008-06-19 15:58:02 -0400 
!      By: David Kopriva  
!
!      Modification history:
!        2008-06-19: Created by David Kopriva
!        XXXX-XX-XX: Gonzalo Rubio implemented cross-product metrics
!        2017-05-05: Andrés Rueda implemented polynomial anisotropy
!      Contains:
!         ALGORITHM 101: MappedGeometryClass
!         ALGORITHM 102: ConstructMappedGeometry
!
!////////////////////////////////////////////////////////////////////////
!
Module MappedGeometryClass 
   USE SMConstants
   USE TransfiniteMapClass
   USE NodalStorageClass
   IMPLICIT NONE
!
!     ---------
!     Constants
!     ---------
!
      integer, parameter :: EFRONT = 1, EBACK = 2, EBOTTOM = 3
      integer, parameter :: ERIGHT = 4, ETOP = 5, ELEFT = 6
!
!     -----
!     Class
!     -----
!
      TYPE MappedGeometry
            INTEGER                                         :: Nx, Ny, Nz                    ! Polynomial order
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: jGradXi, jGradEta, jGradZeta  ! 
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: x                             ! Position of points in absolute coordinates
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: xb                            ! 
            REAL(KIND=RP), DIMENSION(:,:,:)   , ALLOCATABLE :: jacobian, scal                ! 
            REAL(KIND=RP), DIMENSION(:,:,:,:) , ALLOCATABLE :: normal                        ! Normal vectors on the nodes that are on the boundaries of the element (allocated with maximum order)
            
            CONTAINS
            
            PROCEDURE :: construct => ConstructMappedGeometry
            PROCEDURE :: destruct  => DestructMappedGeometry
      END TYPE MappedGeometry
      
      type MappedGeometryFace
         real(kind=RP), dimension(:,:,:),   allocatable :: x
         real(kind=RP), dimension(:,:)  , allocatable :: scal   ! |ja^i|: Normalization term of the normal vectors on a face
         real(kind=RP), dimension(:,:,:), allocatable :: normal ! normal vector on a face
         
         contains
            procedure :: construct => ConstructMappedGeometryFace
            procedure :: destruct  => DestructMappedGeometryFace
      end type MappedGeometryFace
      
      LOGICAL       :: useCrossProductMetrics = .true. ! A switch for debugging purposes. Cross product metrics are fine (and more precise) for 2D geometries... But not for 3D.
                                                        ! Before changing, read: Kopriva, David A. "Metric identities and the discontinuous spectral element method on curvilinear meshes." Journal of Scientific Computing 26.3 (2006): 301-327.
!
!  ========
   CONTAINS 
!  ========
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ConstructMappedGeometry( self, spAxi, spAeta, spAzeta, mapper )
      IMPLICIT NONE
!
!      ---------
!      Arguments
!      ---------
!
      CLASS(MappedGeometry)  , intent(inout) :: self
      TYPE(TransfiniteHexMap), intent(in)    :: mapper
      TYPE(NodalStorage)     , intent(in)    :: spAxi
      TYPE(NodalStorage)     , intent(in)    :: spAeta
      TYPE(NodalStorage)     , intent(in)    :: spAzeta
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER       :: Nx, Ny, Nz, Nmax
      INTEGER       :: i, j, k
      REAL(KIND=RP) :: nrm
      REAL(KIND=RP) :: grad_x(3,3), jGrad(3)
!
!     -----------
!     Allocations
!     -----------
!
      Nx        = spAxi   % N
      Ny        = spAeta  % N
      Nz        = spAzeta % N
      Nmax      = MAX(Nx,Ny,Nz)
      self % Nx = Nx
      self % Ny = Ny
      self % Nz = Nz
      
      ALLOCATE( self % JGradXi  (3,0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % JGradEta (3,0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % JGradZeta(3,0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % jacobian   (0:Nx,0:Ny,0:Nz) )
      ALLOCATE( self % x        (3,0:Nx,0:Ny,0:Nz)    )
      
      ALLOCATE( self % xb    (3,0:Nmax,0:Nmax,6)     )
      ALLOCATE( self % normal(3,0:Nmax,0:Nmax,6) )
      ALLOCATE( self % scal    (0:Nmax,0:Nmax,6)     )
!
!     --------------------------
!     Compute interior locations
!     --------------------------
!
!     TODO: this could be inconsistent if curveOrder > N. Otherwise I think is ok.
      DO k = 0, Nz
         DO j= 0, Ny       
            DO i = 0,Nx 
               self % x(:,i,j,k) = mapper %  transfiniteMapAt([spAxi % x(i), spAeta % x(j), spAzeta % x(k)])
            END DO
         END DO
      END DO
!
!     ----------------------
!     Compute face locations
!     ----------------------
!
      ! y-z planes
      DO j = 0, Nz
         DO i = 0, Ny
            self % xb(:,i,j,ELEFT)   = mapper % transfiniteMapAt([-1.0_RP    , spAeta % x(i), spAzeta % x(j)])
            self % xb(:,i,j,ERIGHT)  = mapper % transfiniteMapAt([ 1.0_RP    , spAeta % x(i), spAzeta % x(j)])
         END DO
      END DO 
      
      ! x-y planes
      DO j = 0, Ny
         DO i = 0, Nx
            self % xb(:,i,j,EBOTTOM) = mapper % transfiniteMapAt([spAxi % x(i), spAeta % x(j),    -1.0_RP   ])
            self % xb(:,i,j,ETOP)    = mapper % transfiniteMapAt([spAxi % x(i), spAeta % x(j),     1.0_RP   ])
         END DO
      END DO 
      
      ! x-z planes
      DO j = 0, Nz
         DO i = 0, Nx
            self % xb(:,i,j,EBACK)   = mapper % transfiniteMapAt([spAxi % x(i),  1.0_RP     , spAzeta % x(j)  ])
            self % xb(:,i,j,EFRONT)  = mapper % transfiniteMapAt([spAxi % x(i), -1.0_RP     , spAzeta % x(j)  ])
         END DO
      END DO 
!
!     ------------
!     Metric terms
!     ------------
!
!
!     ------------------------------------------------------------
!     If the faces are straight, the CrossProductForm is OK
!     If there are curved faces, the Conservative form is required
!     see Kopriva 2006
!     ------------------------------------------------------------
!
      IF (useCrossProductMetrics .OR. isHex8(mapper)) THEN 
      
         CALL computeMetricTermsCrossProductForm(self, spAxi, spAeta, spAzeta, mapper)
!
!     ----------------
!     Boundary Normals - Must be evaluated at the boundaries!
!     ----------------
!
         ! y-z planes
         DO j = 0, Nz
            DO i = 0, Ny
!
!           ---------
!           Left face
!           ---------
!
               grad_x = mapper % metricDerivativesAt([-1.0_RP    , spAeta % x(i), spAzeta % x(j)])
               CALL vCross(grad_x(:,2), grad_x(:,3), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,ELEFT) = -jGrad/nrm
               self % scal(i,j,ELEFT)     = nrm
!
!           ----------
!           Right face
!           ----------
!
               grad_x = mapper % metricDerivativesAt([ 1.0_RP    , spAeta % x(i), spAzeta % x(j)])
               CALL vCross(grad_x(:,2), grad_x(:,3), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,ERIGHT) = jGrad/nrm
               self % scal(i,j,ERIGHT)     = nrm     
            END DO
         END DO
         
         ! x-y planes
         DO j = 0, Ny
            DO i = 0, Nx
!
!           -----------
!           bottom face
!           -----------
!
               grad_x = mapper % metricDerivativesAt([spAxi % x(i), spAeta % x(j),    -1.0_RP   ])
               CALL vCross(grad_x(:,1), grad_x(:,2), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,EBOTTOM) = -jGrad/nrm
               self % scal(i,j,EBOTTOM)     = nrm
!
!           --------
!           top face
!           --------
!
               grad_x = mapper % metricDerivativesAt([spAxi % x(i), spAeta % x(j),     1.0_RP   ])
               CALL vCross(grad_x(:,1), grad_x(:,2), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,ETOP) = jGrad/nrm
               self % scal(i,j,ETOP)     = nrm
            END DO
         END DO
         
         ! x-z planes
         DO j = 0, Nz
            DO i = 0, Nx
!
!           ----------
!           front face
!           ----------
!
               grad_x = mapper % metricDerivativesAt([spAxi % x(i), -1.0_RP     , spAzeta % x(j)  ])
               CALL vCross(grad_x(:,3), grad_x(:,1), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,EFRONT) = -jGrad/nrm
               self % scal(i,j,EFRONT)     = nrm
!
!           ---------
!           back face
!           ---------
!
               grad_x = mapper % metricDerivativesAt([spAxi % x(i),  1.0_RP     , spAzeta % x(j)  ])
               CALL vCross(grad_x(:,3), grad_x(:,1), jGrad)
               nrm = NORM2(jGrad)
               self % normal(:,i,j,EBACK) = jGrad/nrm
               self % scal(i,j,EBACK)     = nrm
              
            END DO
         END DO
         
      ELSE
         
         CALL computeMetricTermsConservativeForm(self, spAxi, spAeta, spAzeta, mapper)
      
      ENDIF
      
   END SUBROUTINE ConstructMappedGeometry
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructMappedGeometry(self)
         IMPLICIT NONE 
         CLASS(MappedGeometry) :: self
         DEALLOCATE( self % jGradXi, self % jGradEta, self % jGradZeta, self % jacobian )
         DEALLOCATE( self % x)
         DEALLOCATE( self % xb )
         DEALLOCATE( self % normal, self % scal )
      END SUBROUTINE DestructMappedGeometry
!
!//////////////////////////////////////////////////////////////////////// 
!
!  -----------------------------------------------------------------------------------
!  Computation of the metric terms on a face: TODO only the Left element (rotation 0)
!  -----------------------------------------------------------------------------------
   subroutine ConstructMappedGeometryFace(self, Nf, spA, geom, hexMap, side)
      use PhysicsStorage
      implicit none
      class(MappedGeometryFace), intent(inout)  :: self
      integer,                   intent(in)     :: Nf(2)
      type(NodalStorage),        intent(in)     :: spA(2)
      type(MappedGeometry),      intent(in)     :: geom
      type(TransfiniteHexMap),   intent(in)     :: hexMap
      integer,                   intent(in)     :: side
!
!     ---------------
!     Local variables
!     ---------------
!
      integer        :: i, j, l
      real(kind=RP)  :: x_xi(NDIM,0:Nf(1),0:Nf(2))
      real(kind=RP)  :: x_eta(NDIM,0:Nf(1),0:Nf(2))
      
      
      allocate( self % x(NDIM, 0:Nf(1), 0:Nf(2)))
      allocate( self % scal(0:Nf(1), 0:Nf(2)))
      allocate( self % normal(NDIM, 0:Nf(1), 0:Nf(2)))
!
!     First step: get the surface coordinates directly from the mapping
!     -----------------------------------------------------------------
      select case(side)
         case(ELEFT)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([-1.0_RP      , spA(1) % x(i), spA(2) % x(j) ])
            end do ; end do
         
         case(ERIGHT)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([ 1.0_RP      , spA(1) % x(i), spA(2) % x(j) ])
            end do ; end do
         
         case(EBOTTOM)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i), spA(2) % x(j),    -1.0_RP    ])
            end do ; end do
            
         case(ETOP)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i), spA(2) % x(j),     1.0_RP    ])
            end do ; end do
            
         case(EFRONT)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i), -1.0_RP      , spA(2) % x(j) ])
            end do ; end do
            
         case(EBACK)
            do j = 0, Nf(2) ; do i = 0, Nf(1)
               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i),  1.0_RP      , spA(2) % x(j) ])
            end do ; end do
            
      end select
!
!     Get the mappings interpolant derivatives
!     ----------------------------------------
      x_xi = 0.0_RP
      do j = 0, Nf(2) ; do i = 0, Nf(1) ; do l = 0, Nf(1)  
         x_xi(:,i,j) = x_xi(:,i,j) + spA(1) % D(i,l) * self % x(:,l,j)
      end do              ; end do              ; end do

      x_eta = 0.0_RP
      do j = 0, Nf(2) ; do i = 0, Nf(1) ; do l = 0, Nf(2)  
         x_eta(:,i,j) = x_eta(:,i,j) + spA(2) % D(j,l) * self % x(:,i,l)
      end do              ; end do              ; end do
!
!     Compute the metric terms
!     ------------------------
      if ( useCrossProductMetrics ) then

         do j = 0, Nf(2) ; do i = 0, Nf(1)
            call vCross(x_xi(:,i,j),x_eta(:,i,j), self % normal(:,i,j))
            self % scal(i,j) = norm2(self % normal(:,i,j))
            self % normal(:,i,j) = self % normal(:,i,j) / self % scal(i,j)
         end do              ; end do

      else
      
print*, "Face conservative metrics not implemented"
      end if

      if ( (side .eq. ELEFT) .or. (side .eq. EBACK) .or. (side .eq. EBOTTOM)) self % normal = -self % normal

   end subroutine ConstructMappedGeometryFace

!   subroutine ConstructMappedGeometryFace(self,N,eLgeom,eLSide, spA, hexMap)
!      use PolynomialInterpAndDerivsModule
!      implicit none
!      class(MappedGeometryFace), intent(inout) :: self      !<> 
!      integer                  , intent(in)    :: N(2)      !<  Polynomial orders on face (xi,eta)
!      type(MappedGeometry)     , intent(in)    :: eLgeom    !<  Mapped geometry of the element on the left
!      integer                  , intent(in)    :: eLSide    !<  Side number of this face relative to element on the left
!      type(NodalStorage)       , intent(in)    :: spA(2)    !<  Nodal storage on the face 1:xi, 2:eta
!      TYPE(TransfiniteHexMap)  , intent(in)    :: hexMap
!      !-------------------------------------------------------------------
!      ! Variables for the computation of the metric terms. TODO: They will have to be allocated with the order of the map!
!      real(kind=RP) :: xiDerMat (0:N(1),0:N(1))          ! Derivative matrix in direcion xi (of the face)
!      real(kind=RP) :: etaDerMat(0:N(2),0:N(2))          ! Derivative matrix in direcion eta (of the face)
!      real(kind=RP) :: grad_x   (0:N(1),0:N(2),3,2)      ! Gradient matrices
!      real(kind=RP) :: jGrad  (3,0:N(1),0:N(2))          ! Variable containing the metric term ja^i on the surface, where i is the direction normal to the surface
!      real(kind=RP) :: tArray   (0:N(1),0:N(2))          ! Temporal array
!      real(kind=RP) :: dArray   (0:N(1),0:N(2))          ! Temporal array
!      ! Additional local variables
!      real(kind=RP) :: nrm          ! Norm of the contravariant vector (times J)
!      integer       :: i, j   ! Coordinate counters
!      integer       :: DIMi   ! Dimension counter
!      integer, dimension(0:4) :: iCycle = (/3,1,2,3,1/)
!      !-------------------------------------------------------------------
!      
!!
!!     Allocate quantities
!!     -------------------
!      allocate ( self % x      (3,0:N(1),0:N(2)) )
!      allocate ( self % scal    (0:N(1),0:N(2)) )
!      allocate ( self % normal(3,0:N(1),0:N(2)) )
!      
!!
!!     Compute coordinates of nodes on face
!!     ------------------------------------
!      select case(eLSide)
!         case(ELEFT)
!            do j = 0, N(2) ; do i = 0, N(1)
!               self % x(:,i,j) = hexMap % transfiniteMapAt([-1.0_RP      , spA(1) % x(i), spA(2) % x(j) ])
!            end do ; end do
!         
!         case(ERIGHT)
!            do j = 0, N(2) ; do i = 0, N(1)
!               self % x(:,i,j) = hexMap % transfiniteMapAt([ 1.0_RP      , spA(1) % x(i), spA(2) % x(j) ])
!            end do ; end do
!         
!         case(EBOTTOM)
!            do j = 0, N(2) ; do i = 0, N(1)
!               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i), spA(2) % x(j),    -1.0_RP    ])
!            end do ; end do
!            
!         case(ETOP)
!            do j = 0, N(2) ; do i = 0, N(1)
!               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i), spA(2) % x(j),     1.0_RP    ])
!            end do ; end do
!            
!         case(EFRONT)
!            do j = 0, N(2) ; do i = 0, N(1)
!               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i), -1.0_RP      , spA(2) % x(j) ])
!            end do ; end do
!            
!         case(EBACK)
!            do j = 0, N(2) ; do i = 0, N(1)
!               self % x(:,i,j) = hexMap % transfiniteMapAt([spA(1) % x(i),  1.0_RP      , spA(2) % x(j) ])
!            end do ; end do
!            
!      end select
!
!      
!      CALL PolynomialDerivativeMatrix( N(1), spA(1) % x , xiDerMat  )
!      CALL PolynomialDerivativeMatrix( N(2), spA(2) % x , etaDerMat )
!      ! Note that the matrices don't have to be transposed for MMMultiply2DX 
!      
!!
!!     Calculate face Derivatives
!!     --------------------------
!!
!      do DIMi = 1,3
!         call MMMultiply2D1(xiDerMat ,N,self % x(DIMi,0:,0:),grad_x(0:,0:,DIMi,1))
!         call MMMultiply2D2(etaDerMat,N,self % x(DIMi,0:,0:),grad_x(0:,0:,DIMi,2))
!      end do
!      
!      IF (useCrossProductMetrics .OR. isHex8(hexMap)) THEN 
!         
!         do j = 0, N(2)
!            o i = 0, N(1)
!               CALL vCross(grad_x(i,j,:,1), grad_x(i,j,:,2), jGrad(:,i,j))
!               nrm = NORM2(jGrad(:,i,j))
!               self % normal(:,i,j) = jGrad(:,i,j)/nrm
!               self % scal(i,j)     = nrm
!            end do
!         end do
!         
!!
!!        Flip normal direction according to side position AND cross product anticommutativity
!!        ------------------------------------------------------------------------------------
!         
!         if (eLSide == ELEFT .or. eLSide == EBOTTOM .or. eLSide == EBACK) self % normal = -self % normal
!         
!      else
!!
!!        Compute the metric terms
!!        ------------------------
!!
!!        First term:
!         
!         DO DIMi = 1, 3
!            DO j = 0, N(2)
!               DO i = 0, N(1)
!                  tArray(i,j) = self % x(iCycle(DIMi-1),i,j) * grad_x(i,j,iCycle(DIMi+1),1)
!               END DO
!            END DO
!            CALL MMMultiply2D2( etaDerMat, N, tArray, dArray )
!            jGrad(DIMi,0:,0:) = dArray
!         END DO
!!
!!        Second term:
!         
!         DO DIMi = 1,3
!            DO j = 0, N(2)
!               DO i = 0, N(1)
!                  tArray(i,j) = self % x(iCycle(DIMi-1),i,j)*grad_x(i,j,iCycle(DIMi+1),2)
!               END DO
!            END DO
!            CALL MMMultiply2D1( xiDerMat , N, tArray, dArray )
!            jGrad(DIMi,0:,0:) = jGrad(DIMi,0:,0:) - dArray
!         END DO
!
!   !!//////////////////////////////////////////////////
!   !!    TODO: Interpolate back to Gauss points on face
!   !!//////////////////////////////////////////////////
!
!   !
!   !     Get the normal vector and the scaling term
!   !     ------------------------------------------
!         
!         do j = 0, N(2)
!            do i = 0, N(1)
!               nrm = NORM2( jGrad(:,i,j) )
!               self % normal(:,i,j) = jGrad(:,i,j)/nrm
!               self % scal    (i,j) = nrm 
!            end do
!         end do
!         
!!
!!        Flip normal direction according to side position AND order of operations (first and second term)
!!        ------------------------------------------------------------------------------------------------
!         
!         if (eLSide == ELEFT .or. eLSide == EBOTTOM .or. eLSide == EBACK) self % normal = -self % normal
!         
!         
!      end if
!   end subroutine ConstructMappedGeometryFace
!
!//////////////////////////////////////////////////////////////////////// 
!
      subroutine DestructMappedGeometryFace(self)
         implicit none
         !-------------------------------------------------------------------
         class(MappedGeometryFace), intent(inout) :: self
         !-------------------------------------------------------------------
         
         deallocate (self % x    )
         deallocate (self % scal  )
         deallocate (self % normal)
         
      end subroutine DestructMappedGeometryFace
!
!//////////////////////////////////////////////////////////////////////// 
!  
      SUBROUTINE computeMetricTermsConservativeForm(self, spAxi, spAeta, spAzeta, mapper)  
         use PolynomialInterpAndDerivsModule
!
!     -----------------------------------------------
!     Compute the metric terms in conservative form 
!     -----------------------------------------------
!
      IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(MappedGeometry)   , intent(inout) :: self
         TYPE(NodalStorage)     , intent(in)    :: spAxi
         TYPE(NodalStorage)     , intent(in)    :: spAeta
         TYPE(NodalStorage)     , intent(in)    :: spAzeta
         TYPE(TransfiniteHexMap) :: mapper      
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: xiCL(0:spAxi % N), etaCL(0:spAeta % N),zetaCL(0:spAzeta % N)
      REAL(KIND=RP) :: wXi (0:spAxi % N), wEta (0:spAeta % N),wZeta (0:spAzeta % N)
      REAL(KIND=RP) :: xi(3)
      
      REAL(KIND=RP) :: grad_x(3, 3, 0:spAxi % N, 0:spAeta % N, 0:spAzeta % N)
      REAL(KIND=RP) :: xGauss(3, 0:spAxi % N, 0:spAeta % N, 0:spAzeta % N)
      
      REAL(KIND=RP) :: xiDermat  (0:spAxi % N, 0:spAxi % N)
      REAL(KIND=RP) :: etaDerMat (0:spAeta % N, 0:spAeta % N)
      REAL(KIND=RP) :: zetaDerMat(0:spAzeta % N, 0:spAzeta % N)
      
      REAL(KIND=RP) :: IdentityMatrix(0:max(spAxi % N,spAeta % N,spAzeta % N), 0:max(spAxi % N,spAeta % N,spAzeta % N)) 
      
      REAL(KIND=RP) :: tArray(0:spAxi % N, 0:spAeta % N, 0:spAzeta % N)
      REAL(KIND=RP) :: dArray(0:spAxi % N, 0:spAeta % N, 0:spAzeta % N)
      REAL(KIND=RP) :: vArray(0:spAxi % N, 0:spAeta % N, 0:spAzeta % N)
      
      REAL(KIND=RP) :: jGradXi  (3, 0:spAxi % N, 0:spAeta % N, 0:spAzeta % N)
      REAL(KIND=RP) :: jGradEta (3, 0:spAxi % N, 0:spAeta % N, 0:spAzeta % N)
      REAL(KIND=RP) :: jGradZeta(3, 0:spAxi % N, 0:spAeta % N, 0:spAzeta % N)
      
!~      REAL(KIND=RP) :: jacXi  (0:spAxi % N, 0:spAeta % N, 0:spAzeta % N)
!~      REAL(KIND=RP) :: jacEta (0:spAxi % N, 0:spAeta % N, 0:spAzeta % N)
!~      REAL(KIND=RP) :: jacZeta(0:spAxi % N, 0:spAeta % N, 0:spAzeta % N)
      
      REAL(KIND=RP) :: jGrad (3)
      REAL(KIND=RP) :: nrm
      
      INTEGER :: i, j, k, m, l, n
      
      INTEGER :: polOrder(3)
      
      REAL(KIND=RP), ALLOCATABLE :: xiInterpMat  (:,:)
      REAL(KIND=RP), ALLOCATABLE :: etaInterpMat (:,:)
      REAL(KIND=RP), ALLOCATABLE :: zetaInterpMat(:,:)      
!
!     ---------------------------
!     A convenience mapping array
!     ---------------------------
!
      INTEGER, DIMENSION(0:4) :: iCycle = (/3,1,2,3,1/)
      
      polOrder(:) = [spAxi % N, spAeta % N, spAzeta % N]     
!
!     -------------------------------------------
!     Compute the mesh on the Chebyshev Lobatto 
!     grid and compute the gradients on that mesh
!     -------------------------------------------
!
      
      xiCL   = (/ ( -cos(n*PI/polOrder(1)),n=0, polOrder(1)) /)
      etaCL  = (/ ( -cos(n*PI/polOrder(2)),n=0, polOrder(2)) /)
      zetaCL = (/ ( -cos(n*PI/polOrder(3)),n=0, polOrder(3)) /)
      
      DO l = 0,polOrder(3)
         xi(3) = zetaCL(l)
         DO m = 0,polOrder(2)
            xi(2) = etaCL(m)
            DO n = 0,polOrder(1)
               xi(1) = xiCL(n)
               CALL GeneralHexGradAndMap( xi, xGauss(:,n,m,l), grad_x(:,:,n,m,l), mapper%corners, mapper % faces )
            END DO
         END DO
      END DO
      
      CALL PolynomialDerivativeMatrix( polOrder(1), xiCL  , xiDerMat )
      CALL PolynomialDerivativeMatrix( polOrder(2), etaCL , etaDerMat )
      CALL PolynomialDerivativeMatrix( polOrder(3), zetaCL, zetaDerMat )
      xiDerMat = TRANSPOSE(xiDerMat)     ! The matrices have to be transposed for MatrixMultiplyDeriv
      etaDerMat = TRANSPOSE(etaDerMat)
      zetaDerMat = TRANSPOSE(zetaDerMat)
!
!     -----------------------------------------------------
!     Now compute metric terms at each grid point
!     This computes the jGradXi terms in conservative form.
!     See the notes for the derivations.
!     -----------------------------------------------------
!
!     ----------
!     First term
!     ----------
!
      iLoop: DO i = 1,3
         jLoop : DO j = 1,3
            
            tArray = xGauss(iCycle(j-1),:,:,:)*grad_x(iCycle(j+1),iCycle(i+1),:,:,:)
            
            SELECT CASE (i)
               CASE (1)
                  DO l = 0, polOrder(2)
                     DO m = 0, polOrder(1)   
                        CALL MatrixMultiplyDeriv(tArray(m,l,:), dArray(m,l,:), zetaDerMat, polOrder(3), transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO 
                  jGradXi(j,:,:,:) = dArray
               CASE (2)
                  DO l = 0, polOrder(3)
                     DO m = 0, polOrder(2)   
                        CALL MatrixMultiplyDeriv(tArray(:,m,l), dArray(:,m,l), xiDerMat  , polOrder(1), transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO 
                  jGradEta(j,:,:,:) = dArray
               CASE (3)
                  DO l = 0, polOrder(3)
                     DO m = 0, polOrder(1)   
                        CALL MatrixMultiplyDeriv(tArray(m,:,l), dArray(m,:,l), etaDerMat , polOrder(2), transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO            
                  jGradZeta(j,:,:,:) = dArray
            END SELECT
         
         END DO jLoop
      END DO iLoop
!
!     -----------
!     Second term
!     -----------
!
      iLoop2: DO i = 1,3
         jLoop2 : DO j = 1,3
            
            tArray = xGauss(iCycle(j-1),:,:,:)*grad_x(iCycle(j+1),iCycle(i-1),:,:,:)
            
            SELECT CASE (i)
               CASE (1)
                  DO l = 0, polOrder(3)
                     DO m = 0, polOrder(1)   
                        CALL MatrixMultiplyDeriv(tArray(m,:,l), dArray(m,:,l), etaDerMat , polOrder(2), transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO 
                  jGradXi(j,:,:,:) = jGradXi (j,:,:,:) - dArray
               CASE (2)
                  DO l = 0, polOrder(2)
                     DO m = 0, polOrder(1)   
                        CALL MatrixMultiplyDeriv(tArray(m,l,:), dArray(m,l,:), zetaDerMat, polOrder(3), transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO 
                  jGradEta(j,:,:,:) = jGradEta(j,:,:,:) - dArray
               CASE (3)
                  DO l = 0, polOrder(3)
                     DO m = 0, polOrder(2)   
                        CALL MatrixMultiplyDeriv(tArray(:,m,l), dArray(:,m,l), xiDerMat  , polOrder(1), transp = MXV_DIRECT)
                     ENDDO 
                  ENDDO 
                  jGradZeta(j,:,:,:) = jGradZeta(j,:,:,:) - dArray
            END SELECT
         
         END DO jLoop2
      END DO iLoop2    
      
      
!
!     ------------------------------------
!     Interpolate back onto the Gauss grid
!     ------------------------------------
!
      ALLOCATE( xiInterpMat  (0:polOrder(1),0:polOrder(1)) )
      ALLOCATE( etaInterpMat (0:polOrder(2),0:polOrder(2)) )
      ALLOCATE( zetaInterpMat(0:polOrder(3),0:polOrder(3)) )
     
      CALL BarycentricWeights( polOrder(1), xiCL  , wXi  )
      CALL PolynomialInterpolationMatrix( polOrder(1), polOrder(1), xiCL  , wXi  , spAxi % x  , xiInterpmat)
      CALL BarycentricWeights( polOrder(2), etaCL , wEta )
      CALL PolynomialInterpolationMatrix( polOrder(2), polOrder(2), etaCL , wEta , spAeta % x , etaInterpmat)
      CALL BarycentricWeights( polOrder(3), zetaCL, wZeta)
      CALL PolynomialInterpolationMatrix( polOrder(3), polOrder(3), zetaCL, wZeta, spAzeta % x, zetaInterpmat)      
      
      DO k = 1,3      
         
         tArray(:,:,:) = jGradXi(k,:,:,:)
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )
         jGradXi(k,:,:,:) = vArray(:,:,:)
         
         tArray(:,:,:) = jGradEta(k,:,:,:)
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )
         jGradEta(k,:,:,:) = vArray(:,:,:)
         
         tArray(:,:,:) = jGradZeta(k,:,:,:)
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )
         jGradZeta(k,:,:,:) = vArray(:,:,:)
         
      END DO
      
      self % jGradXi   = jGradXi
      self % jGradEta  = jGradEta
      self % jGradZeta = jGradZeta              
      
!
!     ----------------------------------
!     Compute the jacobian at each point
!     ----------------------------------
!
      DO l = 0,polOrder(3)
         DO m = 0,polOrder(2)
            DO n = 0,polOrder(1)
               tArray(n,m,l) = jacobian3D( grad_x(:,1,n,m,l), grad_x(:,2,n,m,l), grad_x(:,3,n,m,l) )
            END DO
         END DO
      END DO
      CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, etaInterpMat, zetaInterpMat )
      self % jacobian = vArray
      
!
!     ----------------------------
!     Boundary Normals Computation
!          TODO: To be moved to the faces! (different routine)
!     ----------------------------
!
!
!     ----------------------------------------
!     Interpolate to the faces. This is done 
!     by performing a 3D interpolation with 
!     Lobatto nodes in the direction of 
!     interpolation and Gauss nodes in the
!     other two. 
!     ----------------------------------------
      
      IdentityMatrix = 0.0_RP
      do i = 0, max(spAxi % N,spAeta % N,spAzeta % N)
         IdentityMatrix(i,i) = 1.0_RP         
      enddo 
      
      CALL BarycentricWeights( polOrder(1), spAxi % x  , wXi  )
      CALL PolynomialInterpolationMatrix( polOrder(1), polOrder(1), spAxi % x  , wXi  , xiCL  , xiInterpmat  )
      CALL BarycentricWeights( polOrder(2), spAeta % x , wEta )
      CALL PolynomialInterpolationMatrix( polOrder(2), polOrder(2), spAeta % x , wEta , etaCL , etaInterpmat )
      CALL BarycentricWeights( polOrder(3), spAzeta % x, wZeta)
      CALL PolynomialInterpolationMatrix( polOrder(3), polOrder(3), spAzeta % x, wZeta, zetaCL, zetaInterpmat)   
      
      
      DO k = 1,3
         tArray(:,:,:) = jGradXi(k,:,:,:)      ! TODO, change by self % jGradXi
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat   , IdentityMatrix(0:polOrder(2),0:polOrder(2)), &
                                                                                 IdentityMatrix(0:polOrder(3),0:polOrder(3)) )
         jGradXi(k,:,:,:) = vArray(:,:,:)
         
         
         tArray(:,:,:) = jGradEta(k,:,:,:)
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, IdentityMatrix(0:polOrder(1),0:polOrder(1)), etaInterpMat  , &
                                                                 IdentityMatrix(0:polOrder(3),0:polOrder(3)) )
         jGradEta(k,:,:,:) = vArray(:,:,:)
         
         tArray(:,:,:) = jGradZeta(k,:,:,:)
         CALL Interp3DArray( polOrder, tArray, polOrder, vArray, IdentityMatrix(0:polOrder(1),0:polOrder(1)), &
                                                                 IdentityMatrix(0:polOrder(2),0:polOrder(2)), zetaInterpMat  )
         jGradZeta(k,:,:,:) = vArray(:,:,:)
         
      END DO
      
!
!     ---------------------------
!     Evaluation of the normals
!     To Do: divide by the sign 
!     of the Jacobian (see notes)
!        -> Sign of the Jacobian contribution outcommented. Is it really possible for the Jacobian to be negative??
!     ---------------------------
      
!~      tArray(:,:,:) = self % jacobian(:,:,:)
!~      CALL Interp3DArray( polOrder, tArray, polOrder, vArray, xiInterpmat, IdentityMatrix, IdentityMatrix )
!~      jacXi(:,:,:) = vArray(:,:,:)
      
!~      tArray(:,:,:) = self % jacobian(:,:,:)
!~      CALL Interp3DArray( polOrder, tArray, polOrder, vArray, IdentityMatrix, etaInterpMat, IdentityMatrix )
!~      jacEta(:,:,:) = vArray(:,:,:)
      
!~      tArray(:,:,:) = self % jacobian(:,:,:)
!~      CALL Interp3DArray( polOrder, tArray, polOrder, vArray, IdentityMatrix, IdentityMatrix, zetaInterpMat )
!~      jacZeta(:,:,:) = vArray(:,:,:)

      ! y-z planes
      do j = 0, polOrder(3)
         do i = 0, polOrder(2)
!
!           ---------
!           Left face
!           ---------
!
            jGrad(:) = jGradXi (:,0,i,j)
            nrm = NORM2(jGrad)
            self % normal(:,i,j,ELEFT) = -jGrad/nrm !* sign(1._RP, jacXi(0,i,j))
            self % scal(i,j,ELEFT)     = nrm 
!
!           ----------
!           Right face
!           ----------
!
            jGrad(:) = jGradXi(:,polOrder(1),i,j)
            nrm = NORM2(jGrad)
            self % normal(:,i,j,ERIGHT) = jGrad/nrm !* sign(1._RP, jacXi(polOrder(1),i,j))
            self % scal(i,j,ERIGHT)     = nrm 
         end do
      end do
      
      ! x-y planes
      do j = 0, polOrder(2)
         do i = 0, polOrder(1)
!
!
!           -----------
!           bottom face
!           -----------
!
            jGrad(:) = jGradZeta(:,i,j,0)
            nrm = NORM2(jGrad)
            self % normal(:,i,j,EBOTTOM) = -jGrad/nrm !* sign(1._RP, jacZeta(i,j,0))
            self % scal(i,j,EBOTTOM)     = nrm 
!
!           --------
!           top face
!           --------
!
            jGrad(:) = jGradZeta(:,i,j,polOrder(3))
            nrm = NORM2(jGrad)
            self % normal(:,i,j,ETOP) = jGrad/nrm !* sign(1._RP, jacZeta(i,j,polOrder(1)))
            self % scal(i,j,ETOP)     = nrm 
         end do
      end do
            
      ! x-z planes
      do j = 0, polOrder(3)
         do i = 0, polOrder(1)
!
!
!           ----------
!           front face
!           ----------
!
            jGrad(:) = jGradEta(:,i,0,j)
            nrm = NORM2(jGrad)
            self % normal(:,i,j,EFRONT) = -jGrad/nrm !* sign(1._RP, jacEta(i,0,j))
            self % scal(i,j,EFRONT)     = nrm 
!
!           ---------
!           back face
!           ---------
!
            jGrad(:) = jGradEta(:,i,polOrder(2),j)
            nrm = NORM2(jGrad)
            self % normal(:,i,j,EBACK) = jGrad/nrm !* sign(1._RP, jacEta(i,polOrder(1),j))
            self % scal(i,j,EBACK)     = nrm 
           
         end do
      end do
     
      DEALLOCATE(xiInterpMat)
      DEALLOCATE(etaInterpMat)
      DEALLOCATE(zetaInterpMat)
     
      END SUBROUTINE computeMetricTermsConservativeForm
!
!//////////////////////////////////////////////////////////////////////// 
! 
!~       SUBROUTINE computeMetricTermsConservativeDirectly(self, spAxi, spAeta, spAzeta, mapper)  
!~ !
!~ !        -------------------------------------------------------------
!~ !        Compute the metric terms by explicitly writing out the terms:
!~ !
!~ !        Ja1 = [(I(Y_eta*Z))_zeta - (I(Y_zeta*Z))_eta] x^ +
!~ !              [(I(Z_eta*X))_zeta - (I(Z_zeta*X))_eta] y^ +
!~ !              [(I(X_eta*Y))_zeta - (I(X_zeta*Y))_eta] z^ 
!~ !
!~ !        Ja2 = [(I(Y_zeta*Z))_xi - (I(Y_xi*Z))_zeta] x^ +
!~ !              [(I(Z_zeta*X))_xi - (I(Z_xi*X))_zeta] y^ +
!~ !              [(I(X_zeta*Y))_xi - (I(X_xi*Y))_zeta] z^ 
!~ !
!~ !        Ja3 = [(I(Y_xi*Z))_eta - (I(Y_eta*Z))_xi] x^ +
!~ !              [(I(Z_xi*X))_eta - (I(Z_eta*X))_xi] y^ +
!~ !              [(I(X_xi*Y))_eta - (I(X_eta*Y))_xi] z^ 
!~ !
!~ !        -------------------------------------------------------------
!~ !
!~          IMPLICIT NONE  
!~ !
!~ !        ---------
!~ !        Arguments
!~ !        ---------
!~ !
!~          TYPE(MappedGeometry)    :: self
!~          TYPE(NodalStorage)     , intent(in)    :: spAxi
!~          TYPE(NodalStorage)     , intent(in)    :: spAeta
!~          TYPE(NodalStorage)     , intent(in)    :: spAzeta
!~          TYPE(TransfiniteHexMap) :: mapper
!~ !
!~ !        ---------------
!~ !        Local Variables
!~ !        ---------------
!~ !
!~          INTEGER       :: i,j,k, N
!~          REAL(KIND=RP) :: grad_x(3,3)
!~          N = spA % N
!~ !
!~ !        --------------------
!~ !        Compute the Jacobian
!~ !        --------------------
!~ !
!~          DO k = 0, N
!~             DO j = 0,N
!~                DO i = 0,N
!~                   grad_x                 = mapper % metricDerivativesAt([spAxi % x(i), spAeta % x(j), spAzeta % x(k)])
!~                   self % jacobian(i,j,k) = jacobian3D(a1 = grad_x(:,1),a2 = grad_x(:,2),a3 = grad_x(:,3))
!~                END DO   
!~             END DO   
!~          END DO  

         
         
!~       END SUBROUTINE computeMetricTermsConservativeDirectly
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE computeMetricTermsCrossProductForm(self, spAxi, spAeta, spAzeta, mapper)       
!
!     -----------------------------------------------
!     Compute the metric terms in cross product form 
!     -----------------------------------------------
!
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(MappedGeometry)   , intent(inout) :: self
         TYPE(NodalStorage)     , intent(in)    :: spAxi
         TYPE(NodalStorage)     , intent(in)    :: spAeta
         TYPE(NodalStorage)     , intent(in)    :: spAzeta
         TYPE(TransfiniteHexMap), intent(in)    :: mapper
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER       :: i,j,k
         INTEGER       :: Nx, Ny, Nz
         REAL(KIND=RP) :: grad_x(3,3)         

         Nx = spAxi % N
         Ny = spAeta % N
         Nz = spAzeta % N
         
         DO k = 0, Nz
            DO j = 0,Ny
               DO i = 0,Nx
                  grad_x = mapper % metricDerivativesAt([spAxi % x(i), spAeta % x(j), spAzeta % x(k)])
                 
                  CALL vCross( grad_x(:,2), grad_x(:,3), self % jGradXi  (:,i,j,k))
                  CALL vCross( grad_x(:,3), grad_x(:,1), self % jGradEta (:,i,j,k))
                  CALL vCross( grad_x(:,1), grad_x(:,2), self % jGradZeta(:,i,j,k))

                  self % jacobian(i,j,k) = jacobian3D(a1 = grad_x(:,1),a2 = grad_x(:,2),a3 = grad_x(:,3))
               END DO   
            END DO   
         END DO  

      END SUBROUTINE computeMetricTermsCrossProductForm
!
!///////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!     Returns the jacobian of the transformation computed from
!!     the three co-variant coordinate vectors.
!-------------------------------------------------------------------------------
!
      FUNCTION jacobian3D(a1,a2,a3)
!
      USE SMConstants
      IMPLICIT NONE

      REAL(KIND=RP)               :: jacobian3D
      REAL(KIND=RP), DIMENSION(3) :: a1,a2,a3,v
!
      CALL vCross(a2,a3,v)
      jacobian3D = vDot(a1,v)

      END FUNCTION jacobian3D
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!    Returns in result the cross product u x v
!-------------------------------------------------------------------------------
!
      SUBROUTINE vCross(u,v,result)
!
      IMPLICIT NONE
      
      REAL(KIND=RP), DIMENSION(3) :: u,v,result

      result(1) = u(2)*v(3) - v(2)*u(3)
      result(2) = u(3)*v(1) - v(3)*u(1)
      result(3) = u(1)*v(2) - v(1)*u(2)

      END SUBROUTINE vCross
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!    Returns the dot product u.v
!-------------------------------------------------------------------------------
!
      FUNCTION vDot(u,v)
!
      IMPLICIT NONE
      
      REAL(KIND=RP)               :: vDot
      REAL(KIND=RP), DIMENSION(3) :: u,v

      vDot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

      END FUNCTION vDot
!
!///////////////////////////////////////////////////////////////////////////////
!
!-------------------------------------------------------------------------------
!!    Returns the 2-norm of u
!-------------------------------------------------------------------------------
!
      FUNCTION vNorm(u)
!
      IMPLICIT NONE
      
      REAL(KIND=RP)               :: vNorm
      REAL(KIND=RP), DIMENSION(3) :: u

      vNorm = SQRT(u(1)*u(1) + u(2)*u(2) + u(3)*u(3))

      END FUNCTION vNorm
!
!///////////////////////////////////////////////////////////////////////////////
!
      
END Module MappedGeometryClass
