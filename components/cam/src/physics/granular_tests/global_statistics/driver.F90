  Program test_global_statistics

    use shr_kind_mod, only: r8=>shr_kind_r8
    use ppgrid,       only: begchunk, endchunk, pver
    use physics_types,only: physics_state, physics_tend
    use physpkg,      only: phys_init,phys_final
    use constituents, only: cnst_add, cnst_name
    use physconst,    only: mwdry, cpair, mwh2o, cpwv
    use cam_abortutils,only: endrun

    use global_statistics, only: tp_statistics, SMALLER_THAN, GREATER_THAN, &
                                 add_stat_field, get_stat_field_idx, &
                                 get_chunk_stat, get_domain_stat

    implicit none

    integer :: idummy, mm, icol
    integer :: itr, istat, icnst

    type(physics_state), pointer :: phys_state(:) => null()     ! shape: (begchunk:endchunk)
    type(physics_tend ), pointer :: phys_tend(:)  => null()     ! shape: (begchunk:endchunk)

    type(tp_statistics), pointer :: chunk_stat(:,:) => null()   ! shape: (begchunk:endchunk,nfld)
    type(tp_statistics), pointer :: domain_stat(:)  => null()   ! shape: (nfld)

    integer :: ncol = PCOLS-1  ! deliberately making ncol and pcols different to make sure 
                               ! the code works in such situations
    integer :: nstep = STEP
    integer :: nchnk, ichnk
    logical :: l_print_always = .true.

    integer :: n_tot_cnt_in_chunk (BCHNK:ECHNK,PCNST)
    integer :: n_tot_cnt_in_domain(PCNST)

    !----------------------
   
    begchunk = BCHNK
    endchunk = ECHNK

    write(*,*)
    write(*,*) 'Active domain size: ',pver,' levels, ',ncol,' columns, ',(endchunk-begchunk+1),' chunks.'
    write(*,*)

    ! Initialize tracer indices

    call cnst_add('Q',      mwh2o, cpwv,  1.E-12_r8, idummy, longname='Specific humidity')
    call cnst_add('CLDLIQ', mwdry, cpair, 0._r8,     idummy, longname='Grid box averaged cloud liquid amount')
    call cnst_add('CLDICE', mwdry, cpair, 0._r8,     idummy, longname='Grid box averaged cloud ice amount')

    ! Register fields for calculating global statistics.
    ! This has to be done before 'call phys_init' which allocates memory for 
    ! chunk_stat and domain_stat.

    call add_stat_field('Q','test_part_1',GREATER_THAN,1.E-4_r8)
    call add_stat_field('Q','test_part_2',SMALLER_THAN,1.E-4_r8)

    call add_stat_field('CLDLIQ','test_part_1',GREATER_THAN,1.E-9_r8)
    call add_stat_field('CLDLIQ','test_part_2',SMALLER_THAN,1.E-9_r8)

    call add_stat_field('CLDICE','test_part_1',SMALLER_THAN,5._r8)
    call add_stat_field('CLDICE','test_part_2',GREATER_THAN,5._r8)

    ! Allocate memory for state, tend, and stat vectors; read in initial conditions.

    call phys_init(phys_state, phys_tend, chunk_stat, domain_stat, nstep)

    ! Re-assign values to CLDICE to facilitate verification
    do ichnk=begchunk,endchunk
    do icol = 1,ncol
       phys_state(ichnk)%q(icol,:,3) = icol*1._r8
    end do
    end do

    ! Test functionalities for getting global statistics 

    ! The following chunk loop mimics the corresponding loop in phys_run1 (in which tphysbc is called)
    do ichnk=begchunk,endchunk

       write(*,*) '-------------'
       write(*,*) 'chunk ',ichnk

       do icnst = 1,PCNST

         n_tot_cnt_in_chunk(ichnk,icnst) = 0

         itr = icnst
         call get_stat_field_idx(cnst_name(icnst),'test_part_1',istat)
         call get_chunk_stat( ncol, pver, phys_state(ichnk)%q(:ncol,:,itr),                 &! intent(in)
                              phys_state(ichnk)%lat, phys_state(ichnk)%lon, l_print_always, &! intent(in)
                              chunk_stat(ichnk,istat) )                                      ! intent(inout)

         n_tot_cnt_in_chunk(ichnk,icnst) = n_tot_cnt_in_chunk(ichnk,icnst) + chunk_stat(ichnk,istat)%count

         call get_stat_field_idx(cnst_name(icnst),'test_part_2',istat)
         call get_chunk_stat( ncol, pver, phys_state(ichnk)%q(:ncol,:,itr),                 &! intent(in)
                              phys_state(ichnk)%lat, phys_state(ichnk)%lon, l_print_always, &! intent(in)
                              chunk_stat(ichnk,istat) )                                      ! intent(inout)

         n_tot_cnt_in_chunk(ichnk,icnst) = n_tot_cnt_in_chunk(ichnk,icnst) + chunk_stat(ichnk,istat)%count

       end do
    end do

    if (any(n_tot_cnt_in_chunk/=ncol*pver)) then
       write(*,*) n_tot_cnt_in_chunk
       call endrun('Test error in chunk_stat.')
    end if

    ! After all chunks have been processed, get the domain statistics

    nchnk = endchunk-begchunk+1
    write(*,*) '-----------------------------'
    write(*,*) 'entire domain on this process'

    do icnst = 1,PCNST

       n_tot_cnt_in_domain(icnst) = 0

       itr = icnst
       call get_stat_field_idx(cnst_name(icnst),'test_part_1',istat)
       call get_domain_stat( l_print_always, nchnk, chunk_stat(:,istat), domain_stat(istat) ) !intent: 3xin, inout

       n_tot_cnt_in_domain(icnst) = n_tot_cnt_in_domain(icnst) + domain_stat(istat)%count

       call get_stat_field_idx(cnst_name(icnst),'test_part_2',istat)
       call get_domain_stat( l_print_always, nchnk, chunk_stat(:,istat), domain_stat(istat) ) !intent: 3xin, inout

       n_tot_cnt_in_domain(icnst) = n_tot_cnt_in_domain(icnst) + domain_stat(istat)%count
    end do

    !----------------------------------
    ! Check if the results are correct
    !----------------------------------
    ! For each constituent, n_tot_cnt_in_domain should match the total number of cells in the domain

    if (any(n_tot_cnt_in_domain/=ncol*nchnk*pver)) then
       call endrun('Test error in domain_stat.')
    end if

    ! For each consitituent and stat_type, the total count for the domain shound match the sum of 
    ! the counts from all chunks.

    if (any(domain_stat(:)%count/=sum(chunk_stat(:,:)%count,dim=1))) then
       write(*,*) domain_stat(:)%count
       write(*,*) sum( chunk_stat(:,:)%count, dim=1 )
       call endrun('Test error. chunk_stat and domain_stat do not match.')
    end if

    print*, '============================'
    print*, ' Test finished correctly.'
    print*, '============================'

   !Below are intended errors. The code should stop in ENDRUN.
   !call get_stat_field_idx('Q','test_part_3',istat)
   !call get_stat_field_idx('U','test_part_1',istat)

  end Program test_global_statistics
