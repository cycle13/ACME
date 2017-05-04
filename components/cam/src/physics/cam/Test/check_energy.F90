module check_energy

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to check 
!   1. vertically integrated total energy and water conservation for each
!      column within the physical parameterizations
!
!   2. global mean total energy conservation between the physics output state
!      and the input state on the next time step.
!
!   3. add a globally uniform heating term to account for any change of total energy in 2.
!
! Author: Byron Boville  Oct 31, 2002
!         
! Modifications:
!   2003-03  Boville  Add global energy check and fixer.        
!   2016-08  Kai Zhang (kai.zhang@pnnl.gov) & Phil Rasch 
!               1. Modifications to check_energy_chng (water conservation check part) for
!                  MG2 (now includes prognostic rain and snow)    
!               2. Better printout information for energy and water conservation check 
!               3. Additional water conservation check utilities 
!---------------------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver, begchunk, endchunk
  use spmd_utils,      only: masterproc
  
  use physconst,       only: gravit, latvap, latice, cpair, cpairv
  use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_ptend_init
  use constituents,    only: cnst_get_ind, pcnst, cnst_name, cnst_get_type_byind
  use cam_logfile,     only: iulog
  use cam_abortutils,  only: endrun 

  implicit none
  private

! Public types:
  public check_tracers_data

! Public methods
  public :: check_energy_timestep_init  ! timestep initialization of energy integrals and cumulative boundary fluxes
  public :: check_energy_chng      ! check changes in integrals against cumulative boundary fluxes



! Private module data

  logical  :: print_energy_errors = .true.

  real(r8) :: teout_glob           ! global mean energy of output state
  real(r8) :: teinp_glob           ! global mean energy of input state
  real(r8) :: tedif_glob           ! global mean energy difference
  real(r8) :: psurf_glob           ! global mean surface pressure
  real(r8) :: ptopb_glob           ! global mean top boundary pressure
  real(r8) :: heat_glob            ! global mean heating rate

! Physics buffer indices
  
  integer  :: teout_idx  = 0       ! teout index in physics buffer 
  integer  :: dtcore_idx = 0       ! dtcore index in physics buffer 

  type check_tracers_data
     real(r8) :: tracer(pcols,pcnst)       ! initial vertically integrated total (kinetic + static) energy
     real(r8) :: tracer_tnd(pcols,pcnst)   ! cumulative boundary flux of total energy
     integer :: count(pcnst)               ! count of values with significant imbalances
  end type check_tracers_data


!===============================================================================
contains
!===============================================================================


!================================================================================================


  subroutine check_energy_timestep_init(state, tend, col_type)
!-----------------------------------------------------------------------
! Compute initial values of energy and water integrals, 
! zero cumulative tendencies
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state),   intent(inout)    :: state
    type(physics_tend ),   intent(inout)    :: tend
    integer, optional                       :: col_type  ! Flag inidicating whether using grid or subcolumns
!---------------------------Local storage-------------------------------

    real(r8) :: ke(state%ncol)                     ! vertical integral of kinetic energy
    real(r8) :: se(state%ncol)                     ! vertical integral of static energy
    real(r8) :: wv(state%ncol)                     ! vertical integral of water (vapor)
    real(r8) :: wl(state%ncol)                     ! vertical integral of water (liquid)
    real(r8) :: wi(state%ncol)                     ! vertical integral of water (ice)

!!$    real(r8),allocatable :: cpairv_loc(:,:,:)

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    real(r8) :: wr(state%ncol)                     ! vertical integral of rain
    real(r8) :: ws(state%ncol)                     ! vertical integral of snow
    integer :: ixrain
    integer :: ixsnow
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain, abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow, abort=.false.)

!!$    ! cpairv_loc needs to be allocated to a size which matches state and ptend
!!$    ! If psetcols == pcols, cpairv is the correct size and just copy into cpairv_loc
!!$    ! If psetcols > pcols and all cpairv match cpair, then assign the constant cpair
!!$
!!$    if (state%psetcols == pcols) then
!!$       allocate (cpairv_loc(state%psetcols,pver,begchunk:endchunk))
!!$       cpairv_loc(:,:,:) = cpairv(:,:,:)
!!$    else if (state%psetcols > pcols .and. all(cpairv(:,:,:) == cpair)) then
!!$       allocate(cpairv_loc(state%psetcols,pver,begchunk:endchunk))
!!$       cpairv_loc(:,:,:) = cpair
!!$    else
!!$       call endrun('check_energy_timestep_init: cpairv is not allowed to vary when subcolumns are turned on')
!!$    end if

! Compute vertical integrals of dry static energy and water (vapor, liquid, ice)
    ke = 0._r8
    se = 0._r8
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    wr = 0._r8
    ws = 0._r8
    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5_r8*(state%u(i,k)**2 + state%v(i,k)**2)*state%pdel(i,k)/gravit
          se(i) = se(i) + state%s(i,k         )*state%pdel(i,k)/gravit
!!! cam6  se(i) = se(i) +         state%t(i,k)*cpairv_loc(i,k,lchnk)*state%pdel(i,k)/gravit
          wv(i) = wv(i) + state%q(i,k,1       )*state%pdel(i,k)/gravit
       end do
    end do
!!! cam6    do i = 1, ncol
!!! cam6       se(i) = se(i) + state%phis(i)*state%ps(i)/gravit
!!! cam6    end do

    ! Don't require cloud liq/ice to be present.  Allows for adiabatic/ideal phys.
    if (ixcldliq > 1  .and.  ixcldice > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit
             wi(i) = wi(i) + state%q(i,k,ixcldice)*state%pdel(i,k)/gravit
          end do
       end do
    end if

    if (ixrain   > 1  .and.  ixsnow   > 1 ) then
       do k = 1, pver
          do i = 1, ncol
             wr(i) = wr(i) + state%q(i,k,ixrain)*state%pdel(i,k)/gravit
             ws(i) = ws(i) + state%q(i,k,ixsnow)*state%pdel(i,k)/gravit
          end do
       end do
    end if


! Compute vertical integrals of frozen static energy and total water.
    do i = 1, ncol
!!     state%te_ini(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       state%te_ini(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*( wl(i) + wr(i) ) 
       state%tw_ini(i) = wv(i) + wl(i) + wi(i) + wr(i) + ws(i) 

       state%te_cur(i) = state%te_ini(i)
       state%tw_cur(i) = state%tw_ini(i)
    end do

! zero cummulative boundary fluxes 
    tend%te_tnd(:ncol) = 0._r8
    tend%tw_tnd(:ncol) = 0._r8

    state%count = 0

! initialize physics buffer
!!$    if (is_first_step()) then
!!$       call pbuf_set_field(pbuf, teout_idx, state%te_ini, col_type=col_type)
!!$    end if


!!$    deallocate(cpairv_loc)

  end subroutine check_energy_timestep_init

!===============================================================================

  subroutine check_energy_chng(state, tend, name, nstep, ztodt,        &
       flx_vap, flx_cnd, flx_ice, flx_sen)

!-----------------------------------------------------------------------
! Check that the energy and water change matches the boundary fluxes
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

!!$    use cam_history,       only: outfld

    type(physics_state)    , intent(inout) :: state
    type(physics_tend )    , intent(inout) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: flx_vap(:)          ! (pcols) - boundary flux of vapor         (kg/m2/s)
    real(r8), intent(in   ) :: flx_cnd(:)          ! (pcols) -boundary flux of liquid+ice    (m/s) (precip?)
    real(r8), intent(in   ) :: flx_ice(:)          ! (pcols) -boundary flux of ice           (m/s) (snow?)
    real(r8), intent(in   ) :: flx_sen(:)          ! (pcols) -boundary flux of sensible heat (w/m2)

!******************** BAB ******************************************************
!******* Note that the precip and ice fluxes are in precip units (m/s). ********
!******* I would prefer to have kg/m2/s.                                ********
!******* I would also prefer liquid (not total) and ice fluxes          ********
!*******************************************************************************

!---------------------------Local storage-------------------------------

    real(r8) :: te_xpd(state%ncol)                 ! expected value (f0 + dt*boundary_flux)
    real(r8) :: te_dif(state%ncol)                 ! energy of input state - original energy
    real(r8) :: te_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: te_rer(state%ncol)                 ! relative error in energy column
    real(r8) :: te_err(state%ncol)                 ! absolute error in energy column

    real(r8) :: tw_xpd(state%ncol)                 ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tw_dif(state%ncol)                 ! tw_inp - original water
    real(r8) :: tw_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: tw_rer(state%ncol)                 ! relative error in water column
    real(r8) :: tw_err(state%ncol)                 ! absolute error in water column

    real(r8) :: ke(state%ncol)                     ! vertical integral of kinetic energy
    real(r8) :: se(state%ncol)                     ! vertical integral of static energy
    real(r8) :: wv(state%ncol)                     ! vertical integral of water (vapor)
    real(r8) :: wl(state%ncol)                     ! vertical integral of water (liquid)
    real(r8) :: wi(state%ncol)                     ! vertical integral of water (ice)

    real(r8) :: te(state%ncol)                     ! vertical integral of total energy
    real(r8) :: tw(state%ncol)                     ! vertical integral of total water

    real(r8),allocatable :: cpairv_loc(:,:,:)

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    real(r8) :: wr(state%ncol)                     ! vertical integral of rain
    real(r8) :: ws(state%ncol)                     ! vertical integral of snow
    integer :: ixrain
    integer :: ixsnow
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain, abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow, abort=.false.)

!!$    ! cpairv_loc needs to be allocated to a size which matches state and ptend
!!$    ! If psetcols == pcols, cpairv is the correct size and just copy into cpairv_loc
!!$    ! If psetcols > pcols and all cpairv match cpair, then assign the constant cpair
!!$
!!$    if (state%psetcols == pcols) then
!!$       allocate (cpairv_loc(state%psetcols,pver,begchunk:endchunk))
!!$       cpairv_loc(:,:,:) = cpairv(:,:,:)
!!$    else if (state%psetcols > pcols .and. all(cpairv(:,:,:) == cpair)) then
!!$       allocate(cpairv_loc(state%psetcols,pver,begchunk:endchunk))
!!$       cpairv_loc(:,:,:) = cpair
!!$    else
!!$       call endrun('check_energy_chng: cpairv is not allowed to vary when subcolumns are turned on')
!!$    end if

    ! Compute vertical integrals of dry static energy and water (vapor, liquid, ice)
    ke = 0._r8
    se = 0._r8
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    wr = 0._r8
    ws = 0._r8
    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5_r8*(state%u(i,k)**2 + state%v(i,k)**2)*state%pdel(i,k)/gravit
          se(i) = se(i) + state%s(i,k         )*state%pdel(i,k)/gravit
!!!cam6   se(i) = se(i) + state%t(i,k)*cpairv_loc(i,k,lchnk)*state%pdel(i,k)/gravit
          wv(i) = wv(i) + state%q(i,k,1       )*state%pdel(i,k)/gravit
       end do
    end do
!!!cam6    do i = 1, ncol
!!!cam6       se(i) = se(i) + state%phis(i)*state%ps(i)/gravit
!!!cam6    end do

    ! Don't require cloud liq/ice to be present.  Allows for adiabatic/ideal phys.
    if (ixcldliq > 1  .and.  ixcldice > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit
             wi(i) = wi(i) + state%q(i,k,ixcldice)*state%pdel(i,k)/gravit
          end do
       end do
    end if

    if (ixrain   > 1  .and.  ixsnow   > 1 ) then
       do k = 1, pver
          do i = 1, ncol
             wr(i) = wr(i) + state%q(i,k,ixrain)*state%pdel(i,k)/gravit
             ws(i) = ws(i) + state%q(i,k,ixsnow)*state%pdel(i,k)/gravit
          end do
       end do
    end if

    ! Compute vertical integrals of frozen static energy and total water.
    do i = 1, ncol
!!     te(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       te(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*( wl(i) + wr(i) )
       tw(i) = wv(i) + wl(i) + wi(i) + wr(i) + ws(i)
    end do

    ! compute expected values and tendencies
    do i = 1, ncol
       ! change in static energy and total water
       te_dif(i) = te(i) - state%te_cur(i)
       tw_dif(i) = tw(i) - state%tw_cur(i)

       ! expected tendencies from boundary fluxes for last process
       te_tnd(i) = flx_vap(i)*(latvap+latice) - (flx_cnd(i) - flx_ice(i))*1000._r8*latice + flx_sen(i)
       tw_tnd(i) = flx_vap(i) - flx_cnd(i) *1000._r8

       ! cummulative tendencies from boundary fluxes
       tend%te_tnd(i) = tend%te_tnd(i) + te_tnd(i)
       tend%tw_tnd(i) = tend%tw_tnd(i) + tw_tnd(i)

       ! expected new values from previous state plus boundary fluxes
       te_xpd(i) = state%te_cur(i) + te_tnd(i)*ztodt
       tw_xpd(i) = state%tw_cur(i) + tw_tnd(i)*ztodt

       ! absolute error, expected value - input state / previous state 
       te_err(i) = te_xpd(i) - te(i)

       ! relative error, expected value - input state / previous state 
       te_rer(i) = (te_xpd(i) - te(i)) / state%te_cur(i)
    end do

    ! absolute error for total water (allow for dry atmosphere)
    tw_err = 0._r8
    tw_err(:ncol) = tw_xpd(:ncol) - tw(:ncol)

    ! relative error for total water (allow for dry atmosphere)
    tw_rer = 0._r8
    where (state%tw_cur(:ncol) > 0._r8) 
       tw_rer(:ncol) = (tw_xpd(:ncol) - tw(:ncol)) / state%tw_cur(:ncol)
    end where

    ! error checking
    ! the relative error threshold for the water budget has been reduced to 1.e-10
    ! to avoid messages generated by QNEG3 calls
    ! PJR- change to identify if error in energy or water 
    if (print_energy_errors .and. ( name.eq."clubb_tend" .or. name.eq."energy_rounding_est" ) ) then
       if (any(abs(te_rer(1:ncol)) > TOL )) then
          write(iulog,"(a,a,i8,i5)") "significant en_conservation errors from process, nstep, chunk ", &
                name, nstep, lchnk
          write(iulog,"(a5,2a10,6a15)") ' i', 'lat', 'lon', 'en', 'en_from_flux', 'diff', 'exptd diff', 'rerr', 'cum_diff'
	  do i = 1,ncol
             if (abs(te_rer(i)) > 1.E-14_r8 ) then 
                state%count = state%count + 1
                write(iulog,"(i5,2f10.2,6e15.7)") i, state%lat(i), state%lon(i), te(i),te_xpd(i),te_dif(i),  &
                      te_tnd(i)*ztodt,te_rer(i), tend%te_tnd(i)*ztodt
             endif
          end do
       end if

     ! if (any(abs(tw_rer(1:ncol)) > 1.E-10_r8)) then
     !    write(iulog,"(a,a,i8,i5)") "significant w_conservation errors from process, nstep, chunk ", &
     !         name, nstep, lchnk
     !    write(iulog,"(a5,2a10,6a15)") ' i', 'lat', 'lon', 'tw', 'tw_from_flux', 'diff', 'exptd diff', 'rerr', 'cum_diff'
     !    do i = 1, ncol
     !       if ( abs(tw_rer(i)) > 1.E-10_r8) then
     !          state%count = state%count + 1
     !          write(iulog,"(i5,2f10.2,6e15.7)") i, state%lat(i), state%lon(i),tw(i),tw_xpd(i),tw_dif(i), &
     !               tw_tnd(i)*ztodt, tw_rer(i), tend%tw_tnd(i)*ztodt
     !       end if
     !    end do
     ! end if
    end if

    ! copy new value to state
    do i = 1, ncol
       state%te_cur(i) = te(i)
       state%tw_cur(i) = tw(i)
    end do

!!$    deallocate(cpairv_loc)

  end subroutine check_energy_chng

end module check_energy

