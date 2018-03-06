module ViscousMethodClass
   use SMConstants
   use MeshTypes
   use Physics
   use PhysicsStorage
   use MPI_Face_Class
   use DGSEMClass, only: BCState_FCN
   implicit none
!
   private
   public   ViscousMethod_t, BaseClass_ComputeGradient

   type ViscousMethod_t
      contains
         procedure      :: Initialize                => BaseClass_Initialize
         procedure      :: ComputeGradient           => BaseClass_ComputeGradient
         procedure      :: ComputeInnerFluxes        => BaseClass_ComputeInnerFluxes
         procedure      :: ComputeInnerFluxJacobian  => BaseClass_ComputeInnerFluxJacobian
         procedure      :: ComputeInnerFluxesWithSGS => BaseClass_ComputeInnerFluxesWithSGS
         procedure      :: RiemannSolver             => BaseClass_RiemannSolver
         procedure      :: RiemannSolverWithSGS      => BaseClass_RiemannSolverWithSGS
         procedure      :: Describe                  => BaseClass_Describe
   end type ViscousMethod_t
!
!  ========
   contains
!  ========
!
      subroutine BaseClass_Initialize(self, controlVariables)
         use FTValueDictionaryClass
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         use PhysicsStorage
         implicit none
         class(ViscousMethod_t)                :: self
         class(FTValueDictionary),  intent(in) :: controlVariables

      end subroutine BaseClass_Initialize

      subroutine BaseClass_Describe(self)
         implicit none
         class(ViscousMethod_t), intent(in)  :: self

      end subroutine BaseClass_Describe

      subroutine BaseClass_ComputeGradient( self , mesh , time , externalStateProcedure)
!
!        *****************************************************
!           BaseClass computes Local Gradients by default
!              Do not change.. Used by ComputeTimeDerivativeIsolated
!        *****************************************************
!           
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(ViscousMethod_t), intent(in) :: self
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(BCState_FCN)           :: externalStateProcedure
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID

!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % ComputeLocalGradient
         end do
!$omp end do

      end subroutine BaseClass_ComputeGradient

      subroutine BaseClass_ComputeInnerFluxes( self , e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         implicit none
         class(ViscousMethod_t) ,  intent (in)   :: self
         type(Element)                           :: e
         real(kind=RP)           ,  intent (out) :: contravariantFlux(1:N_EQN, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         contravariantFlux = 0.0_RP

      end subroutine BaseClass_ComputeInnerFluxes
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------
!     Subroutine to get the Jacobian of the contravariant fluxes
!     -> dFdQ (:,:,i,j,k,dim)
!                 |     |
!              jac|coord|flux in cartesian direction dim 
!     ----------------------------------------------------------
      subroutine BaseClass_ComputeInnerFluxJacobian( self, e, dFdQ) 
         use ElementClass
         use Physics
         use PhysicsStorage
         implicit none
         !--------------------------------------------
         class(ViscousMethod_t), intent(in)  :: self
         type(Element)         , intent(in)  :: e
         real(kind=RP)         , intent(out) :: dFdQ( NCONS, NCONS, NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3), NDIM )
         !--------------------------------------------
         real(kind=RP), DIMENSION(NCONS,NCONS,NDIM)  :: df_dgradq, dg_dgradq, dh_dgradq
         integer                                     :: i,j,k
         !--------------------------------------------
#if defined(NAVIERSTOKES)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            
            call ViscousJacobian(e % storage % Q(:,i,j,k), df_dgradq, dg_dgradq, dh_dgradq)
            
            
            dFdQ(:,:,:,i,j,k,IX) = e % geom % jGradXi  (1,i,j,k) * df_dgradq + &
                                   e % geom % jGradXi  (2,i,j,k) * dg_dgradq + &
                                   e % geom % jGradXi  (3,i,j,k) * dh_dgradq 

            dFdQ(:,:,:,i,j,k,IY) = e % geom % jGradEta (1,i,j,k) * df_dgradq + &
                                   e % geom % jGradEta (2,i,j,k) * dg_dgradq + &
                                   e % geom % jGradEta (3,i,j,k) * dh_dgradq 

            dFdQ(:,:,:,i,j,k,IZ) = e % geom % jGradZeta(1,i,j,k) * df_dgradq + &
                                   e % geom % jGradZeta(2,i,j,k) * dg_dgradq + &
                                   e % geom % jGradZeta(3,i,j,k) * dh_dgradq 
         end do                ; end do                ; end do
#endif
      end subroutine BaseClass_ComputeInnerFluxJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_ComputeInnerFluxesWithSGS( self , e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         implicit none
         class(ViscousMethod_t) ,  intent (in)   :: self
         type(Element)                           :: e
         real(kind=RP)           ,  intent (out) :: contravariantFlux(1:N_EQN, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         contravariantFlux = 0.0_RP

      end subroutine BaseClass_ComputeInnerFluxesWithSGS

      subroutine BaseClass_RiemannSolver ( self, f, QLeft, QRight, U_xLeft, U_yLeft, U_zLeft, U_xRight, U_yRight, U_zRight, &
                                           nHat, dWall, flux )
         use SMConstants
         use PhysicsStorage
         use FaceClass
         implicit none
         class(ViscousMethod_t)               :: self
         class(Face),   intent(in)            :: f
         real(kind=RP), dimension(N_EQN)      :: QLeft
         real(kind=RP), dimension(N_EQN)      :: QRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_xLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_yLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_zLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_xRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_yRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_zRight
         real(kind=RP), dimension(NDIM)       :: nHat
         real(kind=RP)                        :: dWall
         real(kind=RP), dimension(N_EQN)      :: flux
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         flux = 0.0_RP

      end subroutine BaseClass_RiemannSolver

      subroutine BaseClass_RiemannSolverWithSGS ( self, f, QLeft, QRight, U_xLeft, U_yLeft, U_zLeft, U_xRight, U_yRight, U_zRight, &
                                                  nHat, dWall, flux )
         use SMConstants
         use PhysicsStorage
         use FaceClass
         implicit none
         class(ViscousMethod_t)               :: self
         class(Face),   intent(in)            :: f
         real(kind=RP), dimension(N_EQN)      :: QLeft
         real(kind=RP), dimension(N_EQN)      :: QRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_xLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_yLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_zLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_xRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_yRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_zRight
         real(kind=RP), dimension(NDIM)       :: nHat
         real(kind=RP)                        :: dWall
         real(kind=RP), dimension(N_EQN)      :: flux
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         flux = 0.0_RP

      end subroutine BaseClass_RiemannSolverWithSGS
end module ViscousMethodClass
