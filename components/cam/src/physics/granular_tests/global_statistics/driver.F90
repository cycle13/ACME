  Program TestEnergyFixer

    use shr_kind_mod,  only: r8=>shr_kind_r8
    use physpkg,       only: phys_init,phys_final
    use physics_types, only: physics_state, physics_ptend,physics_tend
  
    use constituents, only: cnst_add
    use physconst, only: mwdry, cpair, mwh2o, cpwv
    use ppgrid, only : begchunk, endchunk

    use global_statistics, only: tp_statistics, get_chunk_stat, get_domain_stat, &
                                 LESS_THAN, GREATER_THAN

    implicit none

    integer :: idummy, mm, ilev, nchnk

    type(physics_state), pointer       :: phys_state(:)
    type(physics_tend ), pointer       :: phys_tend(:)   

   !type(tp_statistics)          :: global_stat
   !type(tp_statistics)          :: domain_stat
    type(tp_statistics), pointer :: chunk_stat(:)
    type(tp_statistics)          :: domain_stat

    integer :: ncol = PCOLS-1
    integer :: nstep = STEP

    integer :: nfld = PCNST
    character(len=20) :: fldname(PCNST)     ! field names

    integer :: ichunk

    logical :: l_print_always = .true.
   
    begchunk = BCHNK
    endchunk = ECHNK

    ! Initialize tracer indices

    call cnst_add('Q',      mwh2o, cpwv,  1.E-12_r8, idummy, longname='Specific humidity')
    call cnst_add('CLDLIQ', mwdry, cpair, 0._r8,     idummy, longname='Grid box averaged cloud liquid amount')
    call cnst_add('CLDICE', mwdry, cpair, 0._r8,     idummy, longname='Grid box averaged cloud ice amount')

    ! Register fields for calculation global statistics. 
    ! Later, this should be included in a subroutine, and the method 
    ! will likely change (-> like pbuf).

    fldname(1) = 'Q'
    fldname(2) = 'CLDLIQ'
    fldname(3) = 'CLDICE'

    allocate(chunk_stat(begchunk:endchunk))

    mm = 2

    domain_stat%procedure_name = "test" 
    domain_stat%field_name = fldname(mm) 

    chunk_stat(:)%procedure_name = domain_stat%procedure_name 
    chunk_stat(:)%field_name     = domain_stat%field_name 

    domain_stat%threshold = 1.E-12_r8
    domain_stat%stat_type = GREATER_THAN

    chunk_stat(:)%threshold = domain_stat%threshold 
    chunk_stat(:)%stat_type = domain_stat%stat_type

    ! Allocate memory for state and tend vectors

    call phys_init(phys_state, phys_tend,nstep)

    ! The actual test

    do ichunk=begchunk,endchunk

       write(*,*) 'chunk ',ichunk

       chunk_stat(ichunk)%extreme_chnk = ichunk

       ilev = 64
       call get_chunk_stat( ncol, fldname(mm), phys_state(ichunk)%q(:ncol,ilev,mm),         &! intent(in)
                            phys_state(ichunk)%lat, phys_state(ichunk)%lon, l_print_always, &! intent(in)
                            chunk_stat(ichunk) )                                             ! intent(inout)
    end do

    nchnk = endchunk-begchunk+1
    write(*,*) 'entire domain '
    call get_domain_stat( l_print_always, nchnk, chunk_stat, domain_stat ) !intent in, in, in, inout

  end Program TestEnergyFixer
