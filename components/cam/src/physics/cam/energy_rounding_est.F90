module rounding_tests

implicit none

CONTAINS

!---------------------------------------------------------------------------------
! This subroutine is meant to evaluate the rounding error related to the total
! energy fixer used at the CAM-CLUBB interface.
!
! An artificial "physical process" is implemented which evaporates a certain
! fraction of the cloud condensate. A deliberate error is implemented in the
! calculation of evaporative cooling, which introduced energy conservation
! error. A total energy fixer is then applied to restore the conservation.
! The intention is then to monitor the energy conservation error printed out by
! the check_energy_chng subroutine (called one level up by tphysbc) to see whether
! any conservation error will be detected and if so, how large.
!---------------------------------------------------------------------------------
subroutine energy_rounding_est( state, ptend_all, hdtime, evap_frac, stend_relerr )

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use constituents,  only: pcnst, cnst_get_ind
  use ppgrid,        only: pcols, pver
  use physconst,     only: latvap,gravit
  use physics_types, only: physics_state, physics_ptend, physics_ptend_init, &
                           physics_state_copy, physics_ptend_copy, physics_update

   ! --------------- !
   ! Input Auguments !
   ! --------------- !

   type(physics_state), intent(in)    :: state                    ! Physics state variables
   real(r8),            intent(in)    :: hdtime                   ! Host model timestep, unit: s
   real(r8),            intent(in)    :: evap_frac
   real(r8),            intent(in)    :: stend_relerr

   ! ---------------------- !
   ! Output Auguments !
   ! ---------------------- !
   type(physics_ptend), intent(out)   :: ptend_all                 ! package tendencies

   !----------
   type(physics_state) :: state1                ! Local copy of state variable
   type(physics_ptend) :: ptend_loc             ! Local tendency from processes, added up to return as ptend_all

   integer :: ncol, ixcldliq, ktop, kbot, i, k
   logical :: lq(pcnst)

   real(r8) :: zfac
   real(r8) :: zdqvdt(pcols,pver)

   real(r8) :: energy_tot_before(pcols)
   real(r8) :: energy_tot_after(pcols)

   real(r8) :: esmall = 0._r8, zcorrection

   !--------
   ncol = state%ncol

   ! Initialize ptend_all.
   ! Among all the tracers, only vapor and liquid mass concentrations will be changes.
   ! Dry static energy (temperature) will also be changed. Winds will not be touched

   !! This whole routine is about fake parameterization for testing
   !! The first part is with an introduction of an error in energy conservation
   lq(:) = .FALSE.
   lq(1) = .TRUE.

   call cnst_get_ind('CLDLIQ',ixcldliq)
   lq(ixcldliq) = .TRUE.
!!$
   call physics_ptend_init(ptend_all, state%psetcols, 'energy_conserv_rounding_estimate', ls=.true., lq=lq )

   ! We will evaporate a certain amount of the cloud liquid

   zdqvdt(:ncol,:pver) = evap_frac * state%q(:ncol,:pver,ixcldliq) /hdtime

   !! assigned values to the tendencies which are strictly the part of interface, not really calculated by clubb
   !! if there is real clubb then these tendencies are translated from clubb
   !! all ptends are the tendencies
   ptend_all%q(:ncol,:pver,1)        =  zdqvdt 
   ptend_all%q(:ncol,:pver,ixcldliq) = -zdqvdt 

   ! But for evaporative cooling, we will deliberately introduce an error by using a zfac that is not unity. 

   zfac = 1._r8 + stend_relerr

   !! assigned values to the tendencies which are strictly the part of interface, not really calculated by clubb
   !! if there is real clubb then these tendencies are translated from clubb
   ptend_all%s(:ncol,:pver) = zfac * latvap*(-zdqvdt(:ncol,:pver))

   !! ***** this marks the end of the first part that introduced error in energy conservation

   !! The part below, the remaining part fixes the energy error 

   !---------------------------------
   ! ENERGY FIXER
   !---------------------------------
   ! Calculate total energy before state update

   call column_total_energy( state, energy_tot_before )


   ! Calculate total energy after state update

   call physics_state_copy( state, state1 )
   call physics_ptend_copy( ptend_all, ptend_loc)
   call physics_update( state1, ptend_loc, hdtime )

   call column_total_energy( state1, energy_tot_after )

   ! Adjust static energy to restore energy conservation.

   do i = 1,ncol

      ! First determine in which layers the process was active.

      ktop = 0
      kbot = 0

      do k = 1,pver
         if ( abs(ptend_all%q(i,k,ixcldliq)).gt.esmall ) then
             ktop = k
             exit
         end if
      end do

      do k = pver,1,-1
         if ( abs(ptend_all%q(i,k,ixcldliq)).gt.esmall ) then
             kbot = k
             exit
         end if
      end do

      ! Now evenly distribute the energy error to layers ktop to kbot

     !print*, 'column ',i, 'ktop, kbot = ',ktop, kbot
      if ( ktop > 0 .and. kbot > 0 .and. kbot >= ktop ) then

         !!This is the energy correction part
         zcorrection =  ( energy_tot_before(i) - energy_tot_after(i) ) &
                       *gravit/sum(state%pdel(i,ktop:kbot))/hdtime
         ptend_all%s(i,ktop:kbot) = ptend_all%s(i,ktop:kbot) + zcorrection
        !print*, 'zcorrection = ', zcorrection 

      end if

   end do

end subroutine energy_rounding_est


!-----------------------------------------------------------------------------
! Subroutine for calculating the column integrated total frozen static energy
!-----------------------------------------------------------------------------
subroutine column_total_energy( state, energy_tot )

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use ppgrid,        only: pcols, pver
  use physconst,     only: latvap, latice, gravit
  use constituents,  only: cnst_get_ind

  type(physics_state), intent(in) :: state
  real(r8),intent(out) :: energy_tot(pcols)

  real(r8) :: ke(pcols)
  real(r8) :: se(pcols)
  real(r8) :: wv(pcols)
  real(r8) :: wl(pcols)
  real(r8) :: zrhodz

  integer :: i, k, ncol
  integer :: ixq, ixcldliq

  ncol = state%ncol

  call cnst_get_ind('Q',ixq)
  call cnst_get_ind('CLDLIQ',ixcldliq)

  ke  = 0._r8
  se  = 0._r8
  wv  = 0._r8
  wl  = 0._r8

  do k=1,pver
    do i=1,ncol
      zrhodz = state%pdel(i,k)/gravit
      se(i) = se(i) + state%s(i,k) *zrhodz
      ke(i) = ke(i) + 0.5_r8*(state%u(i,k)**2+state%v(i,k)**2) *zrhodz
      wv(i) = wv(i) + state%q(i,k,ixq) *zrhodz
      wl(i) = wl(i) + state%q(i,k,ixcldliq) *zrhodz
    enddo
  enddo

  energy_tot(:ncol) = se(:ncol) + ke(:ncol) + (latvap+latice)*wv(:ncol) + latice*wl(:ncol)

end subroutine column_total_energy

!---------------------------------
end module rounding_tests
