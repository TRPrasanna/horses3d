module Solution2PltModule
   use SMConstants
   use SolutionFile
   implicit none

   private
   public   Solution2Plt

#define PRECISION_FORMAT "(E13.5)"

   contains
      subroutine Solution2Plt(meshName, solutionName, fixedOrder, basis, Nout)
         use getTask
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: basis
         logical,          intent(in)     :: fixedOrder
         integer,          intent(in)     :: Nout(3)
   
         select case ( basis )

         case(EXPORT_GAUSS)

            if ( fixedOrder ) then
               call Solution2Plt_GaussPoints_FixedOrder(meshName, solutionName, Nout)
   
            else
               call Solution2Plt_GaussPoints(meshName, solutionName)

            end if

         case(EXPORT_HOMOGENEOUS)
            
            call Solution2Plt_Homogeneous(meshName, solutionName, Nout)

         end select

      end subroutine Solution2Plt
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!     Gauss Points procedures
!     -----------------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Solution2Plt_GaussPoints(meshName, solutionName)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         interface
            character(len=LINE_LENGTH) function getFileName(inputLine)
               use SMConstants
               implicit none
               character(len=*), intent(in)     :: inputLine
            end function getFileName
         end interface
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                    :: mesh
         character(len=LINE_LENGTH)      :: meshPltName
         character(len=LINE_LENGTH)      :: solutionFile
         character(len=1024)             :: title
         integer                         :: no_of_elements, eID
         integer                         :: fid
         integer                         :: Nmesh(4), Nsol(4)
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
         no_of_elements = mesh % no_of_elements
!
!        Transform zones to the output variables
!        ---------------------------------------
         do eID = 1, no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = e % Nsol
!
!           Construct spectral basis
!           ------------------------
            call addNewSpectralBasis(spA, e % Nmesh)
            call addNewSpectralBasis(spA, e % Nsol)
!
!           Project mesh and solution
!           -------------------------
            call ProjectStorageGaussPoints(e, spA(e % Nmesh(1), e % Nmesh(2), e % Nmesh(3)), spA(e % Nsol(1), e % Nsol(2), e % Nsol(3)))

            end associate
         end do
!
!        Write the solution file name
!        ----------------------------
         solutionFile = trim(getFileName(solutionName)) // ".tec"
!
!        Create the file
!        ---------------
         open(newunit = fid, file = trim(solutionFile), action = "write", status = "unknown")
!
!        Add the title
!        -------------
         write(title,'(A,A,A,A,A)') '"Generated from ',trim(meshName),' and ',trim(solutionName),'"'
         write(fid,'(A,A)') "TITLE = ", trim(title)
!
!        Add the variables
!        -----------------
         write(fid,'(A,A)') 'VARIABLES = "x","y","z"', trim(getOutputVariablesLabel())
!
!        Write each element zone
!        -----------------------
         do eID = 1, no_of_elements
            associate ( e => mesh % elements(eID) )
!
!           Write the tecplot file
!           ----------------------
            call WriteElementToTecplot(fid, e, mesh % refs) 
            end associate
         end do
!
!        Close the file
!        --------------
         close(fid)
      
      end subroutine Solution2Plt_GaussPoints

      subroutine ProjectStorageGaussPoints(e, spAM, spAS)
         use Storage
         use NodalStorageClass
         use ProlongMeshAndSolution
         implicit none
         type(Element_t)     :: e
         type(NodalStorage),  intent(in)  :: spAM
         type(NodalStorage),  intent(in)  :: spAS
         
         e % Nout = e % Nsol
         if ( all(e % Nmesh .eq. e % Nout) ) then
            e % xOut(1:,0:,0:,0:) => e % x

         else
            allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            call prolongMeshToGaussPoints(e, spAM, spAS)

         end if

         e % Qout(0:,0:,0:,1:) => e % Q

      end subroutine ProjectStorageGaussPoints
!
!//////////////////////////////////////////////////////////////////////////////////
!
!     Gauss points with fixed order procedures
!     ----------------------------------------
!
!//////////////////////////////////////////////////////////////////////////////////
!
      subroutine Solution2Plt_GaussPoints_FixedOrder(meshName, solutionName, Nout)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: Nout(3)
         interface
            character(len=LINE_LENGTH) function getFileName(inputLine)
               use SMConstants
               implicit none
               character(len=*), intent(in)     :: inputLine
            end function getFileName
         end interface
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                               :: mesh
         character(len=LINE_LENGTH)                 :: meshPltName
         character(len=LINE_LENGTH)                 :: solutionFile
         character(len=1024)                        :: title
         integer                                    :: no_of_elements, eID
         integer                                    :: fid
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
!
!        Allocate the output spectral basis
!        ----------------------------------
         call spA(Nout(1), Nout(2), Nout(3)) % Construct(Nout(1), Nout(2), Nout(3))
!
!        Write each element zone
!        -----------------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = Nout
!
!           Construct spectral basis
!           ------------------------
            call addNewSpectralBasis(spA, e % Nmesh)
            call addNewSpectralBasis(spA, e % Nsol)
!
!           Construct interpolation matrices
!           --------------------------------
            associate( spAout => spA(Nout(1), Nout(2), Nout(3)) )
            call addNewInterpolationMatrix(Tset, e % Nsol(1), spA(e % Nsol(1), e % Nsol(2), e % Nsol(3)), e % Nout(1), spAout % xi)
            call addNewInterpolationMatrix(Tset, e % Nsol(2), spA(e % Nsol(1), e % Nsol(2), e % Nsol(3)), e % Nout(2), spAout % eta)
            call addNewInterpolationMatrix(Tset, e % Nsol(3), spA(e % Nsol(1), e % Nsol(2), e % Nsol(3)), e % Nout(3), spAout % zeta)
            end associate
!
!           Perform interpolation
!           ---------------------
            call ProjectStorageGaussPoints_FixedOrder(e, spA(e % Nmesh(1), e % Nmesh(2), e % Nmesh(3)), &
                                                            spA(e % Nsol(1), e % Nsol(2), e % Nsol(3)), &
                                                            spA(e % Nout(1), e % Nout(2), e % Nout(3)), &
                                                                    Tset(e % Nout(1), e % Nsol(1)) % T, &
                                                                    Tset(e % Nout(2), e % Nsol(2)) % T, &
                                                                    Tset(e % Nout(3), e % Nsol(3)) % T )

            end associate
         end do
!
!        Write the solution file name
!        ----------------------------
         solutionFile = trim(getFileName(solutionName)) // ".tec"
!
!        Create the file
!        ---------------
         open(newunit = fid, file = trim(solutionFile), action = "write", status = "unknown")
!
!        Add the title
!        -------------
         write(title,'(A,A,A,A,A)') '"Generated from ',trim(meshName),' and ',trim(solutionName),'"'
         write(fid,'(A,A)') "TITLE = ", trim(title)
!
!        Add the variables
!        -----------------
         write(fid,'(A,A)') 'VARIABLES = "x","y","z"', trim(getOutputVariablesLabel())
!
!        Write elements
!        --------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )

            call WriteElementToTecplot(fid, e, mesh % refs)
            end associate
         end do

!
!        Close the file
!        --------------
         close(fid)

      end subroutine Solution2Plt_GaussPoints_FixedOrder

      subroutine ProjectStorageGaussPoints_FixedOrder(e, spAM, spAS, spAout, Tx, Ty, Tz)
         use Storage
         use NodalStorageClass
         use ProlongMeshAndSolution
         implicit none
         type(Element_t)     :: e
         type(NodalStorage),  intent(in)  :: spAM
         type(NodalStorage),  intent(in)  :: spAS
         type(NodalStorage),  intent(in)  :: spAout
         real(kind=RP),       intent(in)  :: Tx(0:e % Nout(1), 0:e % Nsol(1))
         real(kind=RP),       intent(in)  :: Ty(0:e % Nout(2), 0:e % Nsol(2))
         real(kind=RP),       intent(in)  :: Tz(0:e % Nout(3), 0:e % Nsol(3))
!
!        Project mesh
!        ------------         
         if ( all(e % Nmesh .eq. e % Nout) ) then
            e % xOut(1:,0:,0:,0:) => e % x

         else
            allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
            call prolongMeshToGaussPoints(e, spAM, spAout)

         end if
!
!        Project the solution
!        -----------------------------------------------------
         if ( all( e % Nsol .eq. e % Nout ) ) then
            e % Qout(0:,0:,0:,1:) => e % Q
   
         else
            allocate( e % Qout(0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3), 1:5) )
            call prolongSolutionToGaussPoints(e % Nsol, e % Q, e % Nout, e % Qout, Tx, Ty, Tz)

         end if


      end subroutine ProjectStorageGaussPoints_FixedOrder
!
!////////////////////////////////////////////////////////////////////////////
!
!     Homogeneous procedures
!     ----------------------
!
!////////////////////////////////////////////////////////////////////////////
!
      subroutine Solution2Plt_Homogeneous(meshName, solutionName, Nout)
         use Storage
         use NodalStorageClass
         use SharedSpectralBasis
         use OutputVariables
         implicit none  
         character(len=*), intent(in)     :: meshName
         character(len=*), intent(in)     :: solutionName
         integer,          intent(in)     :: Nout(3)
         interface
            character(len=LINE_LENGTH) function getFileName(inputLine)
               use SMConstants
               implicit none
               character(len=*), intent(in)     :: inputLine
            end function getFileName
         end interface
!
!        ---------------
!        Local variables
!        ---------------
!
         type(Mesh_t)                               :: mesh
         character(len=LINE_LENGTH)                 :: meshPltName
         character(len=LINE_LENGTH)                 :: solutionFile
         character(len=1024)                        :: title
         integer                                    :: no_of_elements, eID
         integer                                    :: fid
         real(kind=RP)                              :: xi(0:Nout(1)), eta(0:Nout(2)), zeta(0:Nout(3))
         integer                                    :: i
!
!        Read the mesh and solution data
!        -------------------------------
         call mesh % ReadMesh(meshName)
         call mesh % ReadSolution(SolutionName)
!
!        Set homogeneous nodes
!        ---------------------
         xi   = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(1),i=0,Nout(1)) /), (/ Nout(1)+1 /) )
         eta  = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(2),i=0,Nout(2)) /), (/ Nout(2)+1 /) )
         zeta = RESHAPE( (/ (-1.0_RP + 2.0_RP*i/Nout(3),i=0,Nout(3)) /), (/ Nout(3)+1 /) )
!
!        Write each element zone
!        -----------------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            e % Nout = Nout
!
!           Construct spectral basis for both mesh and solution
!           ---------------------------------------------------
            call addNewSpectralBasis(spA, e % Nmesh)
            call addNewSpectralBasis(spA, e % Nsol)
!
!           Construct interpolation matrices for the mesh
!           ---------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nmesh(1), spA(e % Nmesh(1), e % Nmesh(2), e % Nmesh(3)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nmesh(2), spA(e % Nmesh(1), e % Nmesh(2), e % Nmesh(3)), e % Nout(2), eta)
            call addNewInterpolationMatrix(Tset, e % Nmesh(3), spA(e % Nmesh(1), e % Nmesh(2), e % Nmesh(3)), e % Nout(3), zeta)

!
!           Construct interpolation matrices for the solution
!           -------------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nsol(1), spA(e % Nsol(1), e % Nsol(2), e % Nsol(3)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nsol(2), spA(e % Nsol(1), e % Nsol(2), e % Nsol(3)), e % Nout(2), eta)
            call addNewInterpolationMatrix(Tset, e % Nsol(3), spA(e % Nsol(1), e % Nsol(2), e % Nsol(3)), e % Nout(3), zeta)
!
!           Perform interpolation
!           ---------------------
            call ProjectStorageHomogeneousPoints(e, Tset(e % Nout(1), e % Nmesh(1)) % T, &
                                                    Tset(e % Nout(2), e % Nmesh(2)) % T, &
                                                    Tset(e % Nout(3), e % Nmesh(3)) % T, &
                                                     Tset(e % Nout(1), e % Nsol(1)) % T, &
                                                     Tset(e % Nout(2), e % Nsol(2)) % T, &
                                                     Tset(e % Nout(3), e % Nsol(3)) % T )


            end associate
         end do
!
!        Write the solution file name
!        ----------------------------
         solutionFile = trim(getFileName(solutionName)) // ".tec"
!
!        Create the file
!        ---------------
         open(newunit = fid, file = trim(solutionFile), action = "write", status = "unknown")
!
!        Add the title
!        -------------
         write(title,'(A,A,A,A,A)') '"Generated from ',trim(meshName),' and ',trim(solutionName),'"'
         write(fid,'(A,A)') "TITLE = ", trim(title)
!
!        Add the variables
!        -----------------
         write(fid,'(A,A)') 'VARIABLES = "x","y","z"', trim(getOutputVariablesLabel())
!
!        Write elements
!        --------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )

            call WriteElementToTecplot(fid, e, mesh % refs)
            end associate
         end do

!
!        Close the file
!        --------------
         close(fid)

      end subroutine Solution2Plt_Homogeneous

      subroutine ProjectStorageHomogeneousPoints(e, TxMesh, TyMesh, TzMesh, TxSol, TySol, TzSol)
         use Storage
         use NodalStorageClass
         implicit none
         type(Element_t)     :: e
         real(kind=RP),       intent(in)  :: TxMesh(0:e % Nout(1), 0:e % Nmesh(1))
         real(kind=RP),       intent(in)  :: TyMesh(0:e % Nout(2), 0:e % Nmesh(2))
         real(kind=RP),       intent(in)  :: TzMesh(0:e % Nout(3), 0:e % Nmesh(3))
         real(kind=RP),       intent(in)  :: TxSol(0:e % Nout(1), 0:e % Nsol(1))
         real(kind=RP),       intent(in)  :: TySol(0:e % Nout(2), 0:e % Nsol(2))
         real(kind=RP),       intent(in)  :: TzSol(0:e % Nout(3), 0:e % Nsol(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, iVar, l, m, n
!
!        Project mesh
!        ------------         
         allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % xOut = 0.0_RP

         do n = 0, e % Nmesh(3) ; do m = 0, e % Nmesh(2) ; do l = 0, e % Nmesh(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % xOut(:,i,j,k) = e % xOut(:,i,j,k) + e % x(:,l,m,n) * TxMesh(i,l) * TyMesh(j,m) * TzMesh(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do

!
!        Project the solution
!        --------------------
         allocate( e % Qout(0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3), 1:5) )
         e % Qout = 0.0_RP

         do iVar = 1, 5
            do n = 0, e % Nsol(3) ; do m = 0, e % Nsol(2) ; do l = 0, e % Nsol(1)
               do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
                  e % Qout(i,j,k,iVar) = e % Qout(i,j,k,iVar) + e % Q(l,m,n,iVar) * TxSol(i,l) * TySol(j,m) * TzSol(k,n)
               end do            ; end do            ; end do
            end do            ; end do            ; end do
         end do

      end subroutine ProjectStorageHomogeneousPoints
!
!/////////////////////////////////////////////////////////////////////////////
!
!     Write solution
!     --------------
!
!/////////////////////////////////////////////////////////////////////////////
!
      subroutine WriteElementToTecplot(fid,e,refs)
         use Storage
         use NodalStorageClass
         use prolongMeshAndSolution
         use OutputVariables
         use SolutionFile
         implicit none
         integer,            intent(in) :: fid
         type(Element_t),    intent(in) :: e 
         real(kind=RP),      intent(in) :: refs(NO_OF_SAVED_REFS)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: i,j,k,var
         real(kind=RP)              :: outputVars(0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3),1:no_of_outputVariables)
         real(kind=RP)              :: outputVarsReord(1:no_of_outputVariables,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3))
         character(len=LINE_LENGTH) :: formatout
!
!        Get output variables
!        --------------------
         call ComputeOutputVariables(e % Nout, e % Qout, outputVars, refs)

         do i = 1 , no_of_outputVariables
            outputVarsReord(i,:,:,:) = outputVars(:,:,:,i)
         end do
!
!        Write variables
!        ---------------        
         write(fid,'(A,I0,A,I0,A,I0,A)') "ZONE I=",e % Nout(1)+1,", J=",e % Nout(2)+1, &
                                            ", K=",e % Nout(3)+1,", F=POINT"

         formatout = getFormat()

         do k = 0, e % Nout(3)   ; do j = 0, e % Nout(2)    ; do i = 0, e % Nout(1)
            write(fid,trim(formatout)) e % xOut(:,i,j,k), outputVarsReord(:,i,j,k)
         end do               ; end do                ; end do

      end subroutine WriteElementToTecplot

      character(len=LINE_LENGTH) function getFormat()
         use OutputVariables
         implicit none

         getFormat = ""

         write(getFormat,'(A,I0,A,A)') "(",3+no_of_outputVariables,PRECISION_FORMAT,")"

      end function getFormat
      
end module Solution2PltModule
