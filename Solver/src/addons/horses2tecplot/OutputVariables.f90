!
!//////////////////////////////////////////////////////
!
!   @File:    OutputVariables.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sat Oct 14 20:44:38 2017
!   @Last revision date: Fri Oct 27 18:22:27 2017
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 44bfd062343e60e3001d20111ce1a76e97e9df5b
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module OutputVariables
!
!//////////////////////////////////////////////////////////////////////////////////
!
!        This module selects user-defined output variables for the 
!     tecplot file. The user may add extra variables, this can be done by
!     following the steps:
!
!           * Increase NO_OF_VARIABLES in 1.
!           * Adding a new ID for the variable NEWVAR_V (following the given order).
!           * Adding a new key for the variable, NEWVARKey.
!           * Adding the new key to the keys array.
!           * Adding the "select case" procedure that computes the output
!                 variable.
!
!//////////////////////////////////////////////////////////////////////////////////
!
   use SMConstants
   use PhysicsStorage
   use Headers

   private
   public   no_of_outputVariables
   public   getOutputVariables, ComputeOutputVariables, getOutputVariablesLabel

   integer, parameter   :: STR_VAR_LEN = 16
   integer, parameter   :: NO_OF_VARIABLES = 32
   integer, parameter   :: NO_OF_INVISCID_VARIABLES = 16
!
!  ***************************
!  Variables without gradients
!  ***************************
!
   integer, parameter :: Q_V    = 1
   integer, parameter :: RHO_V  = 2
   integer, parameter :: U_V    = 3
   integer, parameter :: V_V    = 4
   integer, parameter :: W_V    = 5
   integer, parameter :: P_V    = 6
   integer, parameter :: T_V    = 7
   integer, parameter :: Mach_V = 8
   integer, parameter :: S_V    = 9
   integer, parameter :: Vabs_V = 10
   integer, parameter :: Vvec_V = 11
   integer, parameter :: Ht_V   = 12
   integer, parameter :: RHOU_V = 13
   integer, parameter :: RHOV_V = 14
   integer, parameter :: RHOW_V = 15
   integer, parameter :: RHOE_V = 16
!
!  ************************
!  Variables with gradients
!  ************************
!
   integer, parameter :: GRADV_V = 17
   integer, parameter :: UX_V = 18
   integer, parameter :: VX_V = 19
   integer, parameter :: WX_V = 20
   integer, parameter :: UY_V = 21
   integer, parameter :: VY_V = 22
   integer, parameter :: WY_V = 23
   integer, parameter :: UZ_V = 24
   integer, parameter :: VZ_V = 25
   integer, parameter :: WZ_V = 26
   integer, parameter :: OMEGA_V = 27
   integer, parameter :: OMEGAX_V = 28
   integer, parameter :: OMEGAY_V = 29
   integer, parameter :: OMEGAZ_V = 30
   integer, parameter :: OMEGAABS_V = 31
   integer, parameter :: QCRIT_V = 32

   character(len = STR_VAR_LEN), parameter  :: QKey    = "Q"
   character(len = STR_VAR_LEN), parameter  :: RHOKey  = "rho"
   character(len = STR_VAR_LEN), parameter  :: UKey    = "u"
   character(len = STR_VAR_LEN), parameter  :: VKey    = "v"
   character(len = STR_VAR_LEN), parameter  :: WKey    = "w"
   character(len = STR_VAR_LEN), parameter  :: PKey    = "p"
   character(len = STR_VAR_LEN), parameter  :: TKey    = "T"
   character(len = STR_VAR_LEN), parameter  :: MachKey = "Mach"
   character(len = STR_VAR_LEN), parameter  :: SKey    = "s"
   character(len = STR_VAR_LEN), parameter  :: VabsKey = "Vabs"
   character(len = STR_VAR_LEN), parameter  :: VvecKey = "V"
   character(len = STR_VAR_LEN), parameter  :: HtKey   = "Ht"
   character(len = STR_VAR_LEN), parameter  :: RHOUKey = "rhou"
   character(len = STR_VAR_LEN), parameter  :: RHOVKey = "rhov"
   character(len = STR_VAR_LEN), parameter  :: RHOWKey = "rhow"
   character(len = STR_VAR_LEN), parameter  :: RHOEKey = "rhoe"
   character(len = STR_VAR_LEN), parameter  :: gradVKey = "gradV"
   character(len = STR_VAR_LEN), parameter  :: uxKey = "u_x"
   character(len = STR_VAR_LEN), parameter  :: vxKey = "v_x"
   character(len = STR_VAR_LEN), parameter  :: wxKey = "w_x"
   character(len = STR_VAR_LEN), parameter  :: uyKey = "u_y"
   character(len = STR_VAR_LEN), parameter  :: vyKey = "v_y"
   character(len = STR_VAR_LEN), parameter  :: wyKey = "w_y"
   character(len = STR_VAR_LEN), parameter  :: uzKey = "u_z"
   character(len = STR_VAR_LEN), parameter  :: vzKey = "v_z"
   character(len = STR_VAR_LEN), parameter  :: wzKey = "w_z"
   character(len = STR_VAR_LEN), parameter  :: omegaKey = "omega"
   character(len = STR_VAR_LEN), parameter  :: omegaxKey = "omega_x"
   character(len = STR_VAR_LEN), parameter  :: omegayKey = "omega_y"
   character(len = STR_VAR_LEN), parameter  :: omegazKey = "omega_z"
   character(len = STR_VAR_LEN), parameter  :: omegaAbsKey = "omega_abs"
   character(len = STR_VAR_LEN), parameter  :: QCriterionKey = "Qcrit"
   
   

   character(len=STR_VAR_LEN), dimension(32), parameter  :: variableNames = (/ QKey, RHOKey, UKey, VKey, WKey, &
                                                                            PKey, TKey, MachKey, SKey, VabsKey, &
                                                                            VvecKey, HtKey, RHOUKey, RHOVKey, RHOWKey, &
                                                                            RHOEKey, gradVKey, uxKey, vxKey, wxKey, &
                                                                            uyKey, vyKey, wyKey, uzKey, vzKey, wzKey, &
                                                                            omegaKey, omegaxKey, omegayKey, omegazKey, &
                                                                            omegaAbsKey, QCriterionKey /)
                                                               

   integer                :: no_of_outputVariables
   integer, allocatable   :: outputVariableNames(:)
   logical                :: outScale

   contains
!
!/////////////////////////////////////////////////////////////////////////////
!
      subroutine getOutputVariables(flag_name)
         implicit none
         character(len=*), intent(in)  :: flag_name
!
!        ---------------
!        Local variables
!        ---------------
!
         logical                       :: flagPresent
         integer                       :: pos, pos2, i, preliminarNoOfVariables
         integer, allocatable          :: preliminarVariables(:)
         character(len=LINE_LENGTH)    :: flag   
         character(len=STR_VAR_LEN)   :: inputVar

         flagPresent = .false.
         outScale    = .true.

         do i = 1, command_argument_count()
            call get_command_argument(i, flag)
            pos = index(trim(flag),trim(flag_name))

            if ( pos .ne. 0 ) then
               flagPresent = .true.
               exit
            end if 
!
!           Also, check if the dimensionless version is requested
!           -----------------------------------------------------
            pos = index(trim(flag),"--dimensionless")
            if ( pos .ne. 0 ) then
               outScale = .false.
            end if
         end do
!
!        ***********************************************************
!              Read the variables. They are first loaded into
!           a "preliminar" variables, prior to introduce them in 
!           the real outputVariables. This is because some of
!           the output variables lead to multiple variables (e.g.
!           Q or V).
!        ***********************************************************
!
         if ( .not. flagPresent ) then
!
!           Default: export Q
!           -------
            preliminarNoOfVariables = 1
            allocate( preliminarVariables(preliminarNoOfVariables) )
            preliminarVariables(1) = Q_V
   
         else
!
!           Get the variables from the command argument
!           -------------------------------------------
            pos = index(trim(flag),"=")

            if ( pos .eq. 0 ) then
               print*, 'Missing "=" operator in --output-variables flag'
               errorMessage(STD_OUT)
               stop
            end if
!
!           Prepare to read the variable names
!           ----------------------------------
            preliminarNoOfVariables = getNoOfCommas(trim(flag)) + 1
            allocate( preliminarVariables(preliminarNoOfVariables) )
            
            if ( preliminarNoOfVariables .eq. 1 ) then
               read(flag(pos+1:len_trim(flag)),*) inputVar
               preliminarVariables(1) = outputVariableForName(trim(inputVar))
            else
               do i = 1, preliminarNoOfVariables-1
                  pos2 = index(trim(flag(pos+1:)),",") + pos
                  read(flag(pos+1:pos2),*) inputVar
                  preliminarVariables(i) = outputVariableForName(trim(inputVar))
                  pos = pos2
               end do
            
               pos = index(trim(flag),",",BACK=.true.)
               preliminarVariables(preliminarNoOfVariables) = outputVariableForName(flag(pos+1:))
               
            end if
         end if
!
!        *******************************************************
!        Convert the preliminary variables into output variables
!        *******************************************************
!
         no_of_outputVariables = 0    
         
         do i = 1, preliminarNoOfVariables
            no_of_outputVariables = no_of_outputVariables + outputVariablesForVariable(preliminarVariables(i))
         end do

         allocate( outputVariableNames(no_of_outputVariables) )

         pos = 1
         do i = 1, preliminarNoOfVariables
            pos2 = pos + outputVariablesForVariable(preliminarVariables(i)) - 1
            call outputVariablesForPreliminarVariable(preliminarVariables(i), outputVariableNames(pos:pos2) )
      
            pos = pos + outputVariablesForVariable(preliminarVariables(i))
         end do
!
!        *****************************
!        Describe the output variables
!        *****************************
!
         write(STD_OUT,'(/)')
         call Section_Header("Output variables")
         write(STD_OUT,'(/)')
         call Subsection_Header("Selected output variables")
         do i = 1, no_of_outputVariables
            write(STD_OUT,'(30X,A,A)') "* ",trim(variableNames(outputVariableNames(i)))
         end do

         if ( outScale ) then
            write(STD_OUT,'(30X,A,A)') "-> Variables are exported with dimensions."
         else
            write(STD_OUT,'(30X,A,A)') "-> Dimensionless mode."
         end if
         

      end subroutine getOutputVariables

      subroutine ComputeOutputVariables(N, e, output, refs, hasGradients)
         use SolutionFile
         use Storage
         implicit none
         integer, intent(in)          :: N(3)
         class(Element_t), intent(in) :: e
         real(kind=RP), intent(out)   :: output(1:no_of_outputVariables,0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(in)    :: refs(NO_OF_SAVED_REFS)
         logical,       intent(in)    :: hasGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: var, i, j, k
         real(kind=RP) :: Sym, Asym
         
         do var = 1, no_of_outputVariables
            if ( hasGradients .or. (outputVariableNames(var) .le. NO_OF_INVISCID_VARIABLES ) ) then
               associate ( Q   => e % Qout, &
                           U_x => e % U_xout, & 
                           U_y => e % U_yout, & 
                           U_z => e % U_zout )

               select case (outputVariableNames(var))
   
               case(RHO_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * output(var,:,:,:)
   
               case(U_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHOU,i,j,k) / Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)
   
               case(V_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHOV,i,j,k) / Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)
   
               case(W_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHOW,i,j,k) / Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)
   
               case(P_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = (refs(GAMMA_REF) - 1.0_RP)*(Q(IRHOE,i,j,k) - 0.5_RP*&
                                       ( POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))) /Q(IRHO,i,j,k)) 
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * POW2(refs(V_REF)) * output(var,:,:,:)
   
               case(T_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = (refs(GAMMA_REF) - 1.0_RP) * refs(GAMMA_REF) * POW2(refs(MACH_REF)) * (Q(IRHOE,i,j,k) - 0.5_RP*&
                                       ( POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))) /POW2(Q(IRHO,i,j,k))) 
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(T_REF) * output(var,:,:,:)
   
               case(MACH_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))/POW2(Q(IRHO,i,j,k))     ! Vabs**2
                     output(var,i,j,k) = sqrt( output(var,i,j,k) / ( refs(GAMMA_REF)*(refs(GAMMA_REF)-1.0_RP)*(Q(IRHOE,i,j,k)/Q(IRHO,i,j,k)-0.5_RP * output(var,i,j,k)) ) )
                  end do         ; end do         ; end do
   
               case(S_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = refs(GAMMA_REF) * (refs(GAMMA_REF) - 1.0_RP)*(Q(IRHOE,i,j,k) - 0.5_RP * &
                                       ( POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))) /Q(IRHO,i,j,k)) / (Q(IRHO,i,j,k)**refs(GAMMA_REF))
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF)*POW2(refs(V_REF))/refs(RHO_REF)**refs(GAMMA_REF) * output(var,:,:,:)
      
               case(Vabs_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = sqrt(POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k)))/Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(V_REF) * output(var,:,:,:)
   
               case(Ht_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = refs(GAMMA_REF)*Q(IRHOE,i,j,k) - 0.5_RP*(refs(GAMMA_REF)-1.0_RP)*&
                                       ( POW2(Q(IRHOU,i,j,k)) + POW2(Q(IRHOV,i,j,k)) + POW2(Q(IRHOW,i,j,k))) /Q(IRHO,i,j,k) 
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * POW2(refs(V_REF)) * output(var,:,:,:)
   
               case(RHOU_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHOU,i,j,k) 
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * refs(V_REF) * output(var,:,:,:)
   
               case(RHOV_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHOV,i,j,k) 
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * refs(V_REF) * output(var,:,:,:)
         
               case(RHOW_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHOW,i,j,k) 
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * refs(V_REF) * output(var,:,:,:)
   
               case(RHOE_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = Q(IRHOE,i,j,k) 
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = refs(RHO_REF) * POW2(refs(V_REF)) * output(var,:,:,:)
!
!
!              ******************
!              Gradient variables   
!              ******************
!
               case(UX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_x(1,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)
               
               case(VX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_x(2,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)
               
               case(WX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_x(3,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)
               
               case(UY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_y(1,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)
               
               case(VY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_y(2,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)
               
               case(WY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_y(3,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)
               
               case(UZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_z(1,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)
               
               case(VZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_z(2,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)
               
               case(WZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = U_z(3,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)

               case(OMEGAX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = ( U_y(3,i,j,k) - U_z(2,i,j,k) )
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)

               case(OMEGAY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = ( U_z(1,i,j,k) - U_x(3,i,j,k) )
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)

               case(OMEGAZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = ( U_x(2,i,j,k) - U_y(1,i,j,k) )
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)

               case(OMEGAABS_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output(var,i,j,k) = sqrt(  POW2( U_y(3,i,j,k) - U_z(2,i,j,k) ) &
                                              + POW2( U_z(1,i,j,k) - U_x(3,i,j,k) ) &
                                              + POW2( U_x(2,i,j,k) - U_y(1,i,j,k) ) )
                  end do         ; end do         ; end do
                  if ( outScale ) output(var,:,:,:) = POW2(refs(V_REF)) * output(var,:,:,:)

               case(QCRIT_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     Sym =   POW2( U_x(1,i,j,k) ) + POW2( U_y(2,i,j,k) ) + POW2( U_z(3,i,j,k) )  &
                           + 2.0_RP *( POW2( 0.5_RP * (U_x(2,i,j,k) + U_y(1,i,j,k)) ) +          &
                                       POW2( 0.5_RP * (U_x(3,i,j,k) + U_z(1,i,j,k)) ) +          &
                                       POW2( 0.5_RP * (U_y(3,i,j,k) + U_z(2,i,j,k)) ) )

                     Asym =   2.0_RP *( POW2( 0.5_RP * (U_x(2,i,j,k) - U_y(1,i,j,k)) ) +        &
                                        POW2( 0.5_RP * (U_x(3,i,j,k) - U_z(1,i,j,k)) ) +        &
                                        POW2( 0.5_RP * (U_y(3,i,j,k) - U_z(2,i,j,k)) ) )

                     output(var,i,j,k) = 0.5_RP*( Asym - Sym )
                  end do            ; end do            ; end do
               
               end select
               end associate
   
            else
               output(var,:,:,:) = 0.0_RP

            end if
         end do

   
      end subroutine ComputeOutputVariables

      character(len=1024) function getOutputVariablesLabel()
         implicit none
!
!        ---------------
!        Local variables         
!        ---------------
!
         integer  :: iVar

         getOutputVariablesLabel = ""
         do iVar = 1, no_of_outputVariables
            write(getOutputVariablesLabel,'(A,A,A,A)') trim(getOutputVariablesLabel), ',"',trim(variableNames(outputVariableNames(iVar))),'"'
         end do

      end function getOutputVariablesLabel

      integer function outputVariablesForVariable(iVar)
!
!        ************************************************
!           This subroutine specifies if a variable
!           implies more than one variable, e.g.,
!           Q = [rho,rhou,rhov,rhow,rhoe], or
!           V = [u,v,w]
!        ************************************************
!
         implicit none
         integer,    intent(in)     :: iVar
   
         select case(iVar)
   
         case(Q_V)
            outputVariablesForVariable = 5

         case(Vvec_V)
            outputVariablesForVariable = 3

         case(gradV_V)
            outputVariablesForVariable = 9

         case(omega_V)
            outputVariablesForVariable = 3

         case default
            outputVariablesForVariable = 1
      
         end select

      end function outputVariablesForVariable

      subroutine OutputVariablesForPreliminarVariable(iVar, output)
         implicit none
         integer, intent(in)     :: iVar
         integer, intent(out)    :: output(:)

         select case(iVar)

         case(Q_V)
            output = (/RHO_V, RHOU_V, RHOV_V, RHOW_V, RHOE_V/)

         case(Vvec_V)
            output = (/U_V, V_V, W_V/)

         case(gradV_V)
            output = (/UX_V, VX_V, WX_V, UY_V, VY_V, WY_V, UZ_V, VZ_V, WZ_V/)

         case(omega_V)
            output = (/OMEGAX_V, OMEGAY_V, OMEGAZ_V/)

         case default
            output = iVar

         end select

      end subroutine OutputVariablesForPreliminarVariable

      integer function outputVariableForName(variableName)
         implicit none
         character(len=*),    intent(in)     :: variableName

         do outputVariableForName = 1, NO_OF_VARIABLES
            if ( trim(variableNames(outputVariableForName)) .eq. trim(variableName) ) return
         end do
!
!        Return "-1" if not found
!        ------------------------
         outputVariableForName = -1
      
      end function outputVariableForName

      integer function getNoOfCommas(input_line)
         implicit none
         character(len=*), intent(in)     :: input_line
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i

         getNoOfCommas = 0
         do i = 1, len_trim(input_line)
            if ( input_line(i:i) .eq. "," ) getNoOfCommas = getNoOfCommas + 1 
         end do

      end function getNoOfCommas

end module OutputVariables