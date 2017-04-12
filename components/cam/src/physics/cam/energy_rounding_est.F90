module rounding_tests

CONTAINS
subroutine energy_rounding_est( state, ptend_all, hdtime )

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use constituents,  only: pcnst, cnst_get_ind
  use ppgrid,        only: pcols, pver
  use physconst,     only: latvap
  use physics_types, only: physics_state, physics_ptend, physics_ptend_init

   ! --------------- !
   ! Input Auguments !
   ! --------------- !

   type(physics_state), intent(in)    :: state                    ! Physics state variables                 [vary]
   real(r8),            intent(in)    :: hdtime                   ! Host model timestep                     [s]

   ! ---------------------- !
   ! Output Auguments !
   ! ---------------------- !
   type(physics_ptend), intent(out)   :: ptend_all                 ! package tendencies

   !----------
  !type(physics_state) :: state1                ! Local copy of state variable
  !type(physics_ptend) :: ptend_loc             ! Local tendency from processes, added up to return as ptend_all

   integer :: ncol
   integer :: ixcldliq
   logical :: lq(pcnst)

   real(r8) :: zfrac_evap
   real(r8) :: zfac
   real(r8) :: zdqvdt(pcols,pver)

   !--------
   ncol = state%ncol

   ! Initialize ptend_all.
   ! Among all the tracers, only vapor and liquid mass concentrations will be changes.
   ! Dry static energy (temperature) will also be changed. Winds will not be touched

   lq(:) = .FALSE.
   lq(1) = .TRUE.

   call cnst_get_ind('CLDLIQ',ixcldliq)
   lq(ixcldliq) = .TRUE.

   call physics_ptend_init(ptend_all, state%psetcols, 'energy_conserv_rounding_estimate', ls=.true., lq=lq )

   ! We will evaporate a certain amount of the cloud liquid

   zfrac_evap = 0.02_r8   ! 5%

   zdqvdt(:ncol,:pver) = zfrac_evap * state%q(:ncol,:pver,ixcldliq) /hdtime

   ptend_all%q(:ncol,:pver,1)        =  zdqvdt 
   ptend_all%q(:ncol,:pver,ixcldliq) = -zdqvdt 

   ! But for evaporative cooling, we will deliberately introduce an error by using a zfac that is not unity. 

   zfac = 1.1_r8

   ptend_all%s(:ncol,:pver) = zfac * latvap*(-zdqvdt(:ncol,:pver))

end subroutine energy_rounding_est

end module rounding_tests
