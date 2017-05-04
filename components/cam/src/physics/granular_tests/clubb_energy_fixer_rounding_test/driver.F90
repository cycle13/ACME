  Program TestEnergyFixer

    use shr_kind_mod,  only: r8=>shr_kind_r8
    use physpkg,       only: phys_init,phys_final
    use physics_types, only: physics_state, physics_ptend,physics_tend,physics_ptend_scale, physics_update
  
    use rounding_tests, only: energy_rounding_est
    use check_energy,   only: check_energy_chng
    use constituents, only: cnst_add
    use physconst, only: mwdry, cpair, mwh2o, cpwv
    use ppgrid, only : begchunk, endchunk

    integer :: idummy

    type(physics_state), pointer       :: phys_state(:)
    type(physics_tend ), pointer       :: phys_tend(:)   
    type(physics_ptend )               :: phys_ptend 

    real(r8) :: hdtime, evap_frac, stend_relerr
    integer :: cld_macmic_num_steps=6
    integer :: ncol=PCOLS
    real(r8),dimension(PCOLS) :: zero
    integer :: nstep=STEP
    integer :: ichunk
    
    hdtime=1800.0_r8
    evap_frac = 0.2_r8
    stend_relerr = 0.1_r8
    zero=0.0
    begchunk = BCHNK
    endchunk = ECHNK

    ! Initialize tracer indices

    call cnst_add('Q',      mwh2o, cpwv,  1.E-12_r8, idummy, longname='Specific humidity')
    call cnst_add('CLDLIQ', mwdry, cpair, 0._r8,     idummy, longname='Grid box averaged cloud liquid amount')
    call cnst_add('CLDICE', mwdry, cpair, 0._r8,     idummy, longname='Grid box averaged cloud ice amount')

    ! Allocate memory for state and tend vectors

    call phys_init(phys_state, phys_tend,nstep)

    ! The actual test

    do ichunk=begchunk,endchunk
       call energy_rounding_est(phys_state(ichunk),phys_ptend,hdtime,evap_frac, stend_relerr )
       call physics_ptend_scale(phys_ptend, 1._r8/cld_macmic_num_steps, ncol)
       call physics_update(phys_state(ichunk), phys_ptend, hdtime, phys_tend(ichunk))
       call check_energy_chng(phys_state(ichunk), phys_tend(ichunk), "energy_rounding_est", 1, hdtime, &
            zero,zero,zero,zero) 
    end do

  end Program TestEnergyFixer
