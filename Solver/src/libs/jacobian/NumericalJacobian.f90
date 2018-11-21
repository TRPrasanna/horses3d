!
!//////////////////////////////////////////////////////
!
!   @File: NumericalJacobian.f90
!   @Author: Andrés Rueda (am.rueda@upm.es) 
!   @Created: Tue Mar 31 17:05:00 2017
!   @Last revision date: Wed Nov 21 19:34:08 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 1c6c630e4fbb918c0c9a98d0bfd4d0b73101e65d
!
!//////////////////////////////////////////////////////
!
!      Routines for computing the Jacobian matrix numerically using the colorings technique
!  ! TODO1: Implement as a class with a destructor to prevent memory leaking
!  ! TODO2: Add MPI communication
!////////////////////////////////////////////////////////////////////////
module NumericalJacobian
   use SMConstants
   use MatrixClass
   use ColorsClass            , only: Colors_t
   use HexMeshClass           , only: HexMesh, Neighbor_t, NUM_OF_NEIGHBORS
   use DGSEMClass
   use ElementClass
   use Jacobian               , only: JACEPS, local2ijk, Look_for_neighbour
   use PhysicsStorage
   use Utilities              , only: Qsort
   use StorageClass           , only: SolutionStorage_t
   use IntegerDataLinkedList  , only: IntegerDataLinkedList_t
   implicit none
   
   private
   public NumericalJacobian_Compute
   
!
!  Module variables
!  -> TODO: They will have to be moved to the class definition or to other types in the future
!  *******************
   type(Neighbor_t), allocatable :: nbr(:)  ! Neighbors information
   type(Colors_t)               :: ecolors
   
   integer        , allocatable :: ndofelm(:)               ! Number of degrees of freedom for each element
   integer        , allocatable :: firstIdx(:)              ! relative position in Jacobian for each element 
   integer        , allocatable :: used(:)                  ! array containing index of elements whose contributions to Jacobian has already been considered (TODO: replace with integer linked list)
   integer                      :: usedctr                  ! counter to fill positions of used
   integer        , allocatable :: irow_0(:), irow(:)       ! Variables to store the row indexes to fill the Jacobian
   integer                      :: num_of_neighbor_levels   ! Number of neighboring levels that affect one element's column of the Jacobian
   integer                      :: max_num_of_neighbors     ! Maximum number of neighboring elements that affect one element's column of the Jacobian
   
contains
   subroutine NumericalJacobian_Compute(sem, nEqn, nGradEqn, t, Matrix, ComputeTimeDerivative, PINFO, eps_in )
      use StopwatchClass
      !-------------------------------------------------------------------
      type(DGSem),                intent(inout), target  :: sem
      integer,                    intent(in)             :: nEqn, nGradEqn
      real(kind=RP),              intent(in)             :: t
      class(Matrix_t)          ,  intent(inout)          :: Matrix
      procedure(ComputeTimeDerivative_f)                 :: ComputeTimeDerivative      !   
      logical,                    OPTIONAL               :: PINFO                      !<? Print information?
      real(kind=RP),              optional               :: eps_in
      !-------------------------------------------------------------------
      integer                                            :: nelm
      integer                                            :: thiscolor, thiselmidx, thiselm         ! specific counters
      integer                                            :: thisdof                           ! specific counters
      integer                                            :: ielm, felm                      
      integer, save                                      :: nnz, totalnnz
      integer                           , save           :: maxndofel
      integer, allocatable, dimension(:), save           :: Nx, Ny, Nz                             ! Polynomial orders
      integer, allocatable, dimension(:), save           :: ndofcol                                ! Maximum number of degrees of freedom in each color        
      integer, allocatable                               :: cols(:)
      integer, allocatable                               :: rows(:)
      integer, allocatable                               :: diag(:)
      real(kind=RP), allocatable, save                   :: Q0(:)
      real(kind=RP), allocatable, save                   :: QDot0(:)
      
      integer :: i, j ! General counters
      integer, dimension(4)                              :: ijkl                                   ! Indexes to locate certain degree of freedom i,j,k...l:equation number
      real(kind=RP), save                                :: eps                                    ! Perturbation magnitude
      
      logical, save                                      :: isfirst = .TRUE.
#if (!defined(NAVIERSTOKES))
      logical                                            :: computeGradients = .true.
#endif
      !-------------------------------------------------------------------
      
!
!     --------------------------------------------------------------------
!     Initialize variables that will be used throughout all the simulation
!     --------------------------------------------------------------------
!
      IF (isfirst) call Stopwatch % CreateNewEvent("Numerical Jacobian construction")
      call Stopwatch % Start("Numerical Jacobian construction")
      
      if (isfirst) then   
         nelm = size(sem % mesh % elements)
!
!        Define the number of needed neighbors
!        -> TODO: Define according to physics and discretization
!        -------------------------------------------------------
#if defined(CAHNHILLIARD)
         num_of_neighbor_levels = 2
#elif defined(NAVIERSTOKES)
         if (flowIsNavierStokes) then
            num_of_neighbor_levels = 2
         else
            num_of_neighbor_levels = 1
         end if
#else
         num_of_neighbor_levels = 2
#endif
         
!
!        Initialize the colorings structure
!        ----------------------------------
         allocate(nbr(nelm))
         CALL Look_for_neighbour(nbr, sem % mesh)
         call ecolors % construct(nbr, num_of_neighbor_levels)
         
!
!        Allocate storage
!        ----------------
         allocate(ndofelm(nelm), firstIdx(nelm+1))
         allocate(Nx(nelm), Ny(nelm), Nz(nelm))
         firstIdx(1) = 1
         
         do i = 1, nelm
            Nx(i) = sem % mesh % elements(i) % Nxyz(1)
            Ny(i) = sem % mesh % elements(i) % Nxyz(2)
            Nz(i) = sem % mesh % elements(i) % Nxyz(3)

            ndofelm(i) = nEqn * (Nx(i)+1)*(Ny(i)+1)*(Nz(i)+1)
            firstIdx(i+1) = firstIdx(i) + ndofelm(i)
         end do         

         maxndofel = MAXVAL(ndofelm)                                             ! TODO: if there's p-adaptation, this value has to be recomputed
!
!        -------------------
!        Row position arrays
!        -------------------
!
         allocate(irow  (maxndofel))
         allocate(irow_0(maxndofel))
         
         irow_0(1:maxndofel) = (/ (i, i=0,maxndofel-1) /)
         
!
!        ---------------------------------------------------------------------------------
!        Get the maximum number of neighbors ["of neighbors" * (num_of_neighbor_levels-1)] 
!        that are needed for the Jacobian computation (mesh dependent)
!        ---------------------------------------------------------------------------------
!
         max_num_of_neighbors = 0 ! Initialize to minimum possible value
         do i=1, nelm
            max_num_of_neighbors = max (getNumOfNeighbors (i, num_of_neighbor_levels), max_num_of_neighbors)
         end do
         
!
!        ---------------------------------------------------------------
!        Allocate the used array that will contain the information about
!        which neighbor elements were already used in the numerical
!        computation of the Jacobian matrix entries
!        -> The neighbors (including itself) and a last entry that will be 0 always (boundary index)
!        ---------------------------------------------------------------
!
         allocate ( used(max_num_of_neighbors+1) )
!
!        -------------------------------------------------------------------------
!        Set max number of nonzero values expected in a row of the Jacobian matrix    TODO: if there's p-adaptation, this has to be recomputed
!              Assumes Legendre-Gauss quadrature and neglects zero values in each 
!                 block (taken into account later when assembling)
!              For Legendre-Gauss-Lobatto much less entries are expected (a node on the
!                 interface has more cols than an interior node)
!              IMPORTANT: These numbers assume conforming meshes!
!        -------------------------------------------------------------------------
!
         nnz = maxndofel * max_num_of_neighbors
!
!        --------------------------------------------------------------
!        Compute the maximum number of degrees of freedom in each color               TODO: if there's p-adaptation, this has to be recomputed
!        --------------------------------------------------------------
!
         allocate(ndofcol(ecolors % num_of_colors))
         ndofcol = 0
         DO thiscolor = 1 , ecolors % num_of_colors
            ielm = ecolors%bounds(thiscolor)             
            felm = ecolors%bounds(thiscolor+1)
            DO thiselmidx = ielm, felm-1              !perturbs a dof in all elements within current color
               thiselm = ecolors%elmnts(thiselmidx)
               ndofcol(thiscolor) = MAX(ndofcol(thiscolor),ndofelm(thiselm))
            END DO
         END DO
         
         allocate(QDot0(size(sem % mesh % storage % QDot)))
         allocate(Q0   (size(sem % mesh % storage % QDot)))
         
         ! All initializations done!
         isfirst = .FALSE.
      end if !(isfirst)
!
!     ---------------------------------------------
!     Set value of eps (currently using Mettot et al. approach with L2 norm because it seems to work)
!        See:
!           > Mettot, Clément, Florent Renac, and Denis Sipp. "Computation of eigenvalue sensitivity to base flow modifications in a discrete framework: Application to open-loop control." Journal of Computational Physics 269 (2014): 234-258.
!           > Knoll, Dana A., and David E. Keyes. "Jacobian-free Newton–Krylov methods: a survey of approaches and applications." Journal of Computational Physics 193.2 (2004): 357-397.
!     --------------------------------------------
!
      if (present(eps_in)) then
         eps = eps_in
      else
         call sem % mesh % storage % local2GlobalQ (sem % NDOF)
         associate (Q => sem % mesh % storage % Q)
         eps = SQRT(EPSILON(eps))*(NORM2(Q)+1._RP)
         end associate
      end if
!
!     ---------------------------
!     Preallocate Jacobian matrix
!     ---------------------------
!
      select type(Matrix_p => Matrix)
         type is(DenseBlockDiagMatrix_t)
            call Matrix_p % Preallocate(nnzs=ndofelm) ! Constructing with block size
            CALL Matrix % Reset
         type is(CSRMat_t)
!~             call GetRowsAndColsVector(sem, nEqn, Matrix_p % numRows, totalnnz, firstIdx, rows, cols, diag)
!~             call Matrix_p % PreAllocateWithStructure(totalnnz, rows, cols, diag) 
            call Matrix_p % Preallocate
         class default ! Construct with nonzeros in each row
            call Matrix % Preallocate(nnz)
            CALL Matrix % Reset
      end select
      
#if defined(CAHNHILLIARD)
      CALL ComputeTimeDerivative( sem % mesh, sem % particles, t, CTD_ONLY_CH_LIN )
#else
      CALL ComputeTimeDerivative( sem % mesh, sem % particles, t, CTD_IGNORE_MODE )
#endif
!
!     Save base state in Q0 and QDot0
!     -------------------------------

#if defined(CAHNHILLIARD)
      call sem % mesh % SetStorageToEqn(2)
#endif
      
      call sem % mesh % storage % local2GlobalQdot (sem % NDOF)
      call sem % mesh % storage % local2GlobalQ    (sem % NDOF)
      QDot0 = sem % mesh % storage % QDot
      Q0    = sem % mesh % storage % Q
!
!     ------------------------------------------
!     Compute numerical Jacobian using colorings
!     ------------------------------------------
!

!
!     Go through every color to obtain its elements' contribution to the Jacobian
!     ***************************************************************************
      do thiscolor = 1 , ecolors % num_of_colors
         ielm = ecolors%bounds(thiscolor)             ! Initial element of the color
         felm = ecolors%bounds(thiscolor+1)           ! Final element of the color
!         
!        Iterate through the DOFs in thiscolor
!           ( Computes one column for each dof within an elment )
!        ********************************************************
         do thisdof = 1, ndofcol(thiscolor)
            
!           Perturb the current degree of freedom in all elements within current color
!           --------------------------------------------------------------------------
            DO thiselmidx = ielm, felm-1              
               thiselm = ecolors%elmnts(thiselmidx)
               IF (ndofelm(thiselm)<thisdof) CYCLE    ! Do nothing if the DOF exceeds the NDOF of thiselm
               
               ijkl = local2ijk(thisdof,nEqn,Nx(thiselm),Ny(thiselm),Nz(thiselm))
               
               sem%mesh%elements(thiselm)% storage % Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) = &
                                                   sem%mesh%elements(thiselm)% storage % Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) + eps 
            ENDDO
!
!           Compute the time derivative
!           ---------------------------
#if defined(CAHNHILLIARD)
            CALL ComputeTimeDerivative( sem % mesh, sem % particles, t, CTD_ONLY_CH_LIN )
#else
            CALL ComputeTimeDerivative( sem % mesh, sem % particles, t, CTD_IGNORE_MODE )
#endif
            call sem % mesh % storage % local2GlobalQdot (sem %NDOF)
            sem % mesh % storage % QDot = (sem % mesh % storage % QDot - QDot0) / eps
            call sem % mesh % storage % global2LocalQdot
            
!
!           Add the contributions to the Jacobian
!           -------------------------------------
            do thiselmidx = ielm, felm-1
               thiselm = ecolors%elmnts(thiselmidx)
               IF (ndofelm(thiselm)<thisdof) CYCLE
               ! Redifine used array and counter
               used    = 0
               usedctr = 1
               
               call AssignColToJacobian(Matrix, sem % mesh % storage, thiselm, thiselm, thisdof, num_of_neighbor_levels)
               
            END DO           
!
!           Restore original values for Q (TODO: this can be improved)
!           ----------------------------------------------------------
            sem % mesh % storage % Q = Q0
            call sem % mesh % storage % global2LocalQ
         ENDDO
      ENDDO
      
      CALL Matrix % Assembly()                             ! Matrix A needs to be assembled before being used
      call Matrix % SpecifyBlockInfo(firstIdx,ndofelm)
      
      call Stopwatch % Pause("Numerical Jacobian construction")
      IF (PRESENT(PINFO)) THEN
         IF (PINFO) PRINT*, "Numerical Jacobian construction: ", Stopwatch % ElapsedTime("Numerical Jacobian construction"), "seconds"
      ENDIF
      call Stopwatch % Reset("Numerical Jacobian construction")
                
   end subroutine NumericalJacobian_Compute
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------------------------------------------
!  Assign to Jacobian the column that corresponds to the element "eID" and the degree of freedom "thisdof" 
!  taking into account the contribution of the neighbors [(depth-1) * "of neighbors"]
!  -------------------------------------------------------------------------------------------------------
   recursive subroutine AssignColToJacobian(Matrix, storage, eID, eIDn, thisdof, depth)
      implicit none
      !-arguments---------------------------------------
      class(Matrix_t)                , intent(inout) :: Matrix  !<> Jacobian Matrix
      type(SolutionStorage_t), target, intent(in)    :: storage !<  storage
      integer                        , intent(in)    :: eID     !<  Element ID 
      integer                        , intent(in)    :: eIDn    !<  ID of the element, whose neighbors' contributions are added
      integer                        , intent(in)    :: thisdof !<  Current degree of freedom
      integer                        , intent(in)    :: depth   !<  Amount of neighbors to visit
      !-local-variables---------------------------------
      integer :: elmnbr                    ! Neighbor element index
      integer :: i                         ! Counter
      integer :: ndof                      ! Number of degrees of freedom of element
      integer :: icol                      ! Current column
      real(kind=RP), pointer :: pbuffer(:) ! Buffer to point to an element's Qdot
      !-------------------------------------------------
      
      if ( (eID  == 0) .or. (eIDn == 0) ) return
!
!     Go through all the neighbors
!     ----------------------------
      do i = 1, size(nbr(eIDn) % elmnt)
         elmnbr = nbr(eIDn) % elmnt(i) 
      
         if (.NOT. any(used == elmnbr)) THEN
            ndof   = ndofelm(elmnbr)
            pbuffer(1:ndof) => storage % elements(elmnbr) % QDot  !maps Qdot array into a 1D pointer
            
            irow = irow_0 + firstIdx(elmnbr)                      !generates the row indices vector
            where ( abs(pbuffer(1:ndof)) .LT. JACEPS) irow = -1   !MatSetvalues will ignore entries with irow=-1
            icol = firstIdx(eID) + thisdof - 1  
            call Matrix % SetColumn(ndof, irow(1:ndof), icol, pbuffer(1:ndof) )
            
            used(usedctr) = elmnbr
            usedctr = usedctr + 1
         end if
         
         if (depth > 1) call AssignColToJacobian(Matrix, storage, eID, elmnbr, thisdof, depth-1)
         
      end do
      
   end subroutine AssignColToJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------------------------------------
!  Returns the number of neighbors [(depth-1) * "of neighbors"] for a specific element (counting itself)
!  -----------------------------------------------------------------------------------------------------
   function getNumOfNeighbors (eID, depth) result(num)
      implicit none
      !-arguments---------------------------------------
      integer                      , intent(in)              :: eID     !<  Element ID 
      integer                      , intent(in)              :: depth   !<  Amount of neighbors to visit
      integer           :: num
      !-local-variables---------------------------------
      type(IntegerDataLinkedList_t) :: neighborsList
      !-------------------------------------------------
      
      num = 0
      if (eID == 0) return 
      
!
!     Create list of already counted elements
!     ---------------------------------------
      neighborsList = IntegerDataLinkedList_t(.FALSE.)
      
!
!     Add neighbors to list
!     ---------------------
      call addNeighborsToList(eID,depth,neighborsList)
      num = neighborsList % no_of_entries
      
!
!     destruct list of already counted elements
!     -----------------------------------------
      call neighborsList % destruct
      
   end function getNumOfNeighbors
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   recursive subroutine addNeighborsToList(eID, depth, neighborsList)
      implicit none
      !-arguments---------------------------------------
      integer                      , intent(in)    :: eID     !<  Element ID 
      integer                      , intent(in)    :: depth   !<  Amount of neighbors to visit
      type(IntegerDataLinkedList_t), intent(inout) :: neighborsList
      !-local-variables---------------------------------
      integer :: elmnbr                    ! Neighbor element index
      integer :: i                         ! Counter
      !-------------------------------------------------
      
      do i = 1, NUM_OF_NEIGHBORS + 1
         elmnbr = nbr(eID) % elmnt(i)
         
         if (elmnbr == 0) cycle
            
         call neighborsList % add(elmnbr)
         if (depth > 1) call addNeighborsToList (elmnbr, depth - 1, neighborsList)
         
      end do
   end subroutine addNeighborsToList
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
!  ! To be deprecated(?)
   subroutine GetRowsAndColsVector(sem, nEqn, nRows, nnz, firstIdx, rows, cols, diag)
      implicit none
      class(DGSEM)         :: sem
      integer, intent(in)  :: nEqn, nRows
      integer, intent(in)  :: firstIdx(1:sem % mesh % no_of_elements)
      integer, intent(out)  :: nnz
      integer, allocatable, intent(out) :: rows(:)
      integer, allocatable, intent(out) :: cols(:)
      integer, allocatable, intent(out) :: diag(:)
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                 :: eID, i, j, k, counter, csr_pos, ieq, nID
      integer                 :: ii, jj, kk, iieq
      integer                 :: pos_i, pos_j
      integer                 :: lb, ub
      integer, dimension(26)  :: neighbours
  
!
!     *********************
!     First loop to get nnz: TODO could I already set rows here?
!     *********************
!
!!$omp parallel private(counter, neighbours, i, j)
!!$omp do reduction(+:nnz)
      nnz = 0
      do eID = 1, sem % mesh % no_of_elements
         associate(e => sem % mesh % elements(eID))
!
!     1/ For each element, get its neighbours
!        ------------------------------------
         counter = 0
         neighbours = -1

         do i = 1, size(nbr(eID) % elmnt)
            if ( nbr(eID) % elmnt(i) .le. 0 ) cycle
            do j = 1, size(nbr(nbr(eID) % elmnt(i)) % elmnt)
               if (nbr(nbr(eID) % elmnt(i)) % elmnt(j) .le. 0) cycle

               if ( .not. any(neighbours .eq. nbr(nbr(eID) % elmnt(i)) % elmnt(j)) ) then
                  counter = counter + 1 
                  neighbours(counter) = nbr(nbr(eID) % elmnt(i)) % elmnt(j)
               end if
            end do
         end do

         do i = 1, counter
            associate(eL => sem % mesh % elements(neighbours(i)))
            nnz = nnz + nEqn*(eL % Nxyz(1)+1)*(eL % Nxyz(2)+1)*(eL % Nxyz(3)+1)*nEqn*(e % Nxyz(1)+1)*(e % Nxyz(2)+1)*(e % Nxyz(3)+1)
            end associate
         end do
         end associate
      end do
!!$omp end do
!!$omp end parallel

      allocate(rows(1:nRows+1))
      allocate(cols(1:nnz))
      allocate(diag(1:nRows))
!
!     ****************************
!     We need to set rows and cols
!     ****************************
!
      csr_pos = 1
      do eID = 1, sem % mesh % no_of_elements
         associate(e => sem % mesh % elements(eID))
!
!     1/ For each element, get its neighbours
!        ------------------------------------
         counter = 0
         neighbours = -1

         do i = 1, size(nbr(eID) % elmnt)
            if ( nbr(eID) % elmnt(i) .le. 0 ) cycle
            do j = 1, size(nbr(nbr(eID) % elmnt(i)) % elmnt)
               if (nbr(nbr(eID) % elmnt(i)) % elmnt(j) .le. 0) cycle

               if ( .not. any(neighbours .eq. nbr(nbr(eID) % elmnt(i)) % elmnt(j)) ) then
                  counter = counter + 1 
                  neighbours(counter) = nbr(nbr(eID) % elmnt(i)) % elmnt(j)
               end if
            end do
         end do
         call Qsort(neighbours(1:counter))

         pos_i = firstIdx(eID)

         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) ; do ieq = 1, nEqn
            rows(pos_i) = csr_pos
            pos_i = pos_i + 1 
            do nID = 1, counter
               associate(neigh => sem % mesh % elements(neighbours(nID)))

               pos_j = firstIdx(neighbours(nID))

               do kk = 0, neigh % Nxyz(3) ; do jj = 0, neigh % Nxyz(2) ; do ii = 0, neigh % Nxyz(1) ; do iieq = 1, nEqn
                    
                  cols(csr_pos) = pos_j

                  pos_j = pos_j + 1
                  csr_pos = csr_pos + 1 
               end do                  ; end do                ; end do                 ; end do
               end associate
            end do   
         end do                ; end do                ; end do                ; end do

         end associate
      end do

      rows(nRows+1) = nnz+1
!
!     **************************
!     Get the diagonal positions
!     **************************
!
      do i = 1, nRows
         lb = rows(i)
         ub = rows(i+1)-1

         pos_i = minloc(abs(cols(lb:ub)-i),dim=1)
         diag(i) = lb + pos_i - 1
      end do

   end subroutine GetRowsAndColsVector
end module NumericalJacobian
