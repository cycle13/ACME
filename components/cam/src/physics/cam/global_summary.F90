module global_summary
!------------------------------------------------------------------------------------
! Description:

! This module provides utilities for getting statistical summaries for global fields.
! For example, the module can be used to identify negative values in a tracer 
! mixing ratio, or energy conservation errors exceeding a pre-defined threshold.
!
! The basic idea is to build a list of fields for which such summaries will be provided 
! during  model integration. For each field, first identify the violation within each 
! chunk of the physics grid; get a total count, and note down the extreme value and 
! its location (chunk/column index, lat, lon, and vertical level index if applicable).
! Then the total count and the extreme among all chunks on a single  MPI process ("domain")
! are obtained. lastly, the domain summaries are collected  by the master process to 
! provide a global summary.
! 
! Each "field" on the list is identified by a field name and a procedure name
! (see components of the derived type tp_stat_smry). For example, total energy error
! from deep convection, total energy error from cloud microphysics, and negative 
! cloud droplet mass concentration from deep convection would be considered 3 distinct
! fields.
!
! History:
!
! First version by Hui Wan (PNNL, 2017-05)
!------------------------------------------------------------------------------------

  use shr_kind_mod,   only: r8=>SHR_KIND_R8
  use shr_kind_mod,   only: shortchar=>SHR_KIND_CS, longchar=>SHR_KIND_CL
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog

  implicit none
  private

  public tp_stat_smry
  public global_smry_init
  public add_smry_field
  public get_smry_field_idx
  public get_chunk_smry
  public get_global_smry

  character(len=shortchar),private,parameter :: THIS_MODULE = 'global_summary'

  !-------------------------------------------
  ! Types of summary supported by this module

  integer,public,parameter :: SMALLER_THAN     = -1
  integer,public,parameter :: GREATER_EQ       =  1
  integer,public,parameter :: ABS_SMALLER_THAN = -2
  integer,public,parameter :: ABS_GREATER_EQ   =  2

  !----------------
  ! Data structure 

  type tp_stat_smry

    character(len=shortchar) :: procedure_name   ! a procedure could be any piece of code
    character(len=shortchar) :: field_name       ! the physical quantity to be evaluated

    integer  :: smry_type  ! one of the summary types defined above
    real(r8) :: threshold  ! threshold specified by developer/user
    integer  :: count = 0  ! total number of cells with values exceeding threshold

    ! extreme value and its location

    real(r8) :: extreme_val
    real(r8) :: extreme_lat  = -999._r8
    real(r8) :: extreme_lon  = -999._r8
    integer  :: extreme_chnk = -999
    integer  :: extreme_col  = -999
    integer  :: extreme_lev  = -999

  end type tp_stat_smry
  !-------------------------------

  integer,parameter       :: max_number_of_smry_fields = 1000
  integer,public          :: current_number_of_smry_fields = 0
  logical                 :: l_smry_arrays_allocated = .false.
  character(len=longchar) :: msg 

#ifdef UNIT_TEST
  logical :: l_print_always = .true.    ! always print message in log file 
                                        ! (even when there are no
                                        ! values exeeding threshold)
#else
  logical :: l_print_always = .false. 
#endif

  !-------------------------------------------------------------------
  ! The variable that contain a list of fields (on all processes)
  ! and the global summary (on the master proc only)

  type(tp_stat_smry) :: global_smry_1d(max_number_of_smry_fields)

  !-------------------------------
  interface get_chunk_smry
    module procedure get_chunk_smry_1_lev_real   ! for fields that do not have a vertical distribution
    module procedure get_chunk_smry_m_lev_real   ! for fields with multiple vertical levels
  end interface get_chunk_smry



contains

  !--------------------------------------------------------------------------------------------
  ! Description: 
  !  The subroutine registers a new field for getting global summary. It is expected to be
  !  called during the initialization of various parameterizations.
  !--------------------------------------------------------------------------------------------
  subroutine add_smry_field( fldname, procname, stattype, threshold, fldidx )

    character(len=*), intent(in)   :: fldname
    character(len=*), intent(in)   :: procname
    real(r8)                       :: threshold
    integer                        :: stattype
    integer, intent(out), optional :: fldidx

    integer :: ii

    ! During the model initialization, we first register all the fields that we would
    ! like to get global summary for, then allocate memory for the chunk/domain summary
    ! arrays. No new fields can be added after that allocation.

    if (l_smry_arrays_allocated) then
        msg = trim(THIS_MODULE)//': subroutine add_smry_field should not be called'//&
              ' after global_smry_init has been called.'
        call endrun(trim(msg))
    end if

    ! Check if the same field from the same procedure has already been registered.

    do ii = 1,current_number_of_smry_fields
       if (trim(global_smry_1d(ii)%field_name) == trim(fldname) .and. &
           trim(global_smry_1d(ii)%procedure_name) == trim(procname)  ) then
           msg = trim(THIS_MODULE)//': field '//trim(fldname)//' from procedure '// &
                 trim(procname)//' has already been added to list of statistical summary.'
           call endrun(trim(msg))
       end if
    end do

    ! Now add a new field

    current_number_of_smry_fields = current_number_of_smry_fields + 1

    if (current_number_of_smry_fields.gt.max_number_of_smry_fields) then
       msg = trim(THIS_MODULE)//': max. No. of fields exceeded when attempting to add '//trim(fldname)
       call endrun(trim(msg))
    end if

    ii = current_number_of_smry_fields
    global_smry_1d(ii)%field_name     = trim(fldname)
    global_smry_1d(ii)%procedure_name = trim(procname)
    global_smry_1d(ii)%threshold      = threshold
    global_smry_1d(ii)%smry_type      = stattype

    if (present(fldidx)) fldidx = ii 
 
  end subroutine add_smry_field

  !--------------------------------------------------------------------------------
  ! Description: 
  !   Allocate memory for the chunk/domain summary variables, and copy meta data 
  !   from global_smry_1d. This subroutine needs to be called during model
  !   after all "add_smry_field" calls. 
  !--------------------------------------------------------------------------------
  subroutine global_smry_init( chunk_smry_2d, domain_smry_1d, begchunk, endchunk )

    type(tp_stat_smry), pointer ::  chunk_smry_2d(:,:)
    type(tp_stat_smry), pointer :: domain_smry_1d(:)
    integer, intent(in) :: begchunk, endchunk

    integer :: ierr, ichnk

    ! Sanity check

    if (l_smry_arrays_allocated) then
        msg = trim(THIS_MODULE)//': attempting to call global_smry_init multiple times.'
        call endrun(trim(msg))
    end if

    !---------------------------------------
    ! Initialize array for domain summaries 
    !---------------------------------------
    ! Allocate memory

    allocate(domain_smry_1d(current_number_of_smry_fields), stat=ierr)
    if( ierr /= 0 ) then
       write(msg,*) trim(THIS_MODULE)//': domain_smry allocation error = ',ierr
       call endrun(trim(msg))
    end if

    ! Copy metadata

    domain_smry_1d(1:current_number_of_smry_fields)%field_name = &
    global_smry_1d(1:current_number_of_smry_fields)%field_name

    domain_smry_1d(1:current_number_of_smry_fields)%procedure_name = &
    global_smry_1d(1:current_number_of_smry_fields)%procedure_name

    domain_smry_1d(1:current_number_of_smry_fields)%smry_type = &
    global_smry_1d(1:current_number_of_smry_fields)%smry_type

    domain_smry_1d(1:current_number_of_smry_fields)%threshold = &
    global_smry_1d(1:current_number_of_smry_fields)%threshold

    !--------------------------------------
    ! Initialize array for chunk summaries 
    !--------------------------------------
    ! Allocate memory

    allocate(chunk_smry_2d(begchunk:endchunk,current_number_of_smry_fields), stat=ierr)
    if( ierr /= 0 ) then
       write(msg,*) trim(THIS_MODULE)//': chunk_smry allocation error = ',ierr
       call endrun(trim(msg))
    end if

    ! Copy metadata

    do ichnk = begchunk,endchunk 

       chunk_smry_2d(ichnk,1:current_number_of_smry_fields)%field_name = &
      global_smry_1d(      1:current_number_of_smry_fields)%field_name

       chunk_smry_2d(ichnk,1:current_number_of_smry_fields)%procedure_name = &
      global_smry_1d(      1:current_number_of_smry_fields)%procedure_name

       chunk_smry_2d(ichnk,1:current_number_of_smry_fields)%smry_type = &
      global_smry_1d(      1:current_number_of_smry_fields)%smry_type

       chunk_smry_2d(ichnk,1:current_number_of_smry_fields)%threshold = &
      global_smry_1d(      1:current_number_of_smry_fields)%threshold

       chunk_smry_2d(ichnk,1:current_number_of_smry_fields)%extreme_chnk = ichnk
    end do

    ! Set flag for sanity check later

    l_smry_arrays_allocated = .true.

  end subroutine global_smry_init

  !---------------------------------------------------------------------
  ! Description:
  !   Find a field on the global summary list and return the index.
  !---------------------------------------------------------------------
  subroutine get_smry_field_idx( fldname, procname, fldidx )

    character(len=*), intent(in)   :: fldname
    character(len=*), intent(in)   :: procname
    integer, intent(out)           :: fldidx

    integer :: ii
    logical :: found_field

    found_field = .false.
    do ii = 1,current_number_of_smry_fields 
       if (trim(global_smry_1d(ii)%field_name) == trim(fldname) .and. &
           trim(global_smry_1d(ii)%procedure_name) == trim(procname)  ) then
          found_field = .true.
          fldidx = ii
          exit
       end if
    end do

    if (.not.found_field) then
       write(msg,*) trim(THIS_MODULE)//': get_smry_field_idx, did not find '// &
                     trim(fldname)//' from '//trim(procname)//' on the list.'
       call endrun(trim(msg))
    end if
 
  end subroutine get_smry_field_idx

  !---------------------------------------------------------------------------------------
  ! Description:
  !   Identify values in a field that exceed a threshold. The total number of such values 
  !   and the location of the extreme are noted down for later use.
  !   This subroutine is meant to operate on all columns (or a subset of columns) 
  !   in a single chunk of CAM's physics grid. The particular incarnation of the 
  !   subroutine deals with fields that have multiple vertical levels. 
  !---------------------------------------------------------------------------------------
  subroutine get_chunk_smry_m_lev_real( ncol, nlev, array, lat, lon, chunk_smry )

    integer,           intent(in)    :: ncol              ! number of columns packed in array
    integer,           intent(in)    :: nlev              ! number of vertical levels
    real(r8),          intent(in)    :: array(ncol,nlev)  ! array of values to be checked
    real(r8),          intent(in)    :: lat(ncol)
    real(r8),          intent(in)    :: lon(ncol)
    type(tp_stat_smry),intent(inout) :: chunk_smry

    ! Local variables

    integer  :: iflag(ncol,nlev) 
    integer  :: idx(2)
    character(len=shortchar) :: smry_type_char
  
    !-------------------------------------------------------------------------
    ! Calculate the total number of grid cells with value exceeding threshold
    ! and identify location of the extremem value.
  
    iflag(:,:) = 0

    SELECT CASE (chunk_smry%smry_type)
    CASE (GREATER_EQ)
      smry_type_char = '>='
      where( array .ge. chunk_smry%threshold ) iflag = 1
      idx = maxloc( array )

    CASE (SMALLER_THAN)
      smry_type_char = '<'
      where( array .lt. chunk_smry%threshold ) iflag = 1
      idx = minloc( array )

    CASE (ABS_GREATER_EQ)
      smry_type_char = 'ABS >='
      WHERE( abs(array) .ge. chunk_smry%threshold ) iflag = 1
      idx = maxloc( abs(array) )

    CASE (ABS_SMALLER_THAN)
      smry_type_char = 'ABS < '
      WHERE( abs(array) .lt. chunk_smry%threshold ) iflag = 1
      idx = minloc( abs(array) )

    END SELECT

    ! Total number of values exceeding threshold

    chunk_smry%count = sum( iflag )

    ! The extreme value

    chunk_smry%extreme_val  = array(idx(1),idx(2))
    chunk_smry%extreme_col  =       idx(1)
    chunk_smry%extreme_lev  =       idx(2)
    chunk_smry%extreme_lat  =   lat(idx(1))
    chunk_smry%extreme_lon  =   lon(idx(1))
  
    ! Send message to log file
  
    if (l_print_always) then
       write(iulog,"(a,i8,a,e15.7,a,e15.7)") &
       '**** '//trim(chunk_smry%field_name)//' from '//trim(chunk_smry%procedure_name)//':', &
       chunk_smry%count, ' values '//trim(smry_type_char), chunk_smry%threshold, &
       ', extreme is ', chunk_smry%extreme_val
    end if
  
  end subroutine get_chunk_smry_m_lev_real

  !---------------------------------------------------------------------------------------
  ! Description:
  !   Identify values in a field that exceed a threshold. The total number of such values 
  !   and the location of the extreme are noted down for later use.
  !   This subroutine is meant to operate on all columns (or a subset of columns) 
  !   in a single chunk of CAM's physics grid. The particular incarnation of the 
  !   subroutine deals with fields that do not have a vertical distribution (e.g. surface
  !   fluxes and vertical integrals). 
  !---------------------------------------------------------------------------------------
  subroutine get_chunk_smry_1_lev_real( ncol, array, lat, lon, chunk_smry )

    integer,           intent(in)    :: ncol              ! number of columns packed in array
    real(r8),          intent(in)    :: array(ncol)       ! array of values to be checked
    real(r8),          intent(in)    :: lat(ncol)
    real(r8),          intent(in)    :: lon(ncol)
    type(tp_stat_smry),intent(inout) :: chunk_smry

    ! Local variables

    integer  :: iflag(ncol) 
    integer  :: idx(1)
    character(len=shortchar) :: smry_type_char
  
    !-----------------------------------------------------------------------
    ! Calculate the total number of columns with value exceeding threshold
    ! and identify location of the extremem value.
  
    iflag(:) = 0

    SELECT CASE (chunk_smry%smry_type)
    CASE (GREATER_EQ)
      smry_type_char = '>='
      where( array .ge. chunk_smry%threshold ) iflag = 1
      idx = maxloc( array )

    CASE (SMALLER_THAN)
      smry_type_char = '<'
      where( array .lt. chunk_smry%threshold ) iflag = 1
      idx = minloc( array )

    CASE (ABS_GREATER_EQ)
      smry_type_char = 'ABS >='
      WHERE( abs(array) .ge. chunk_smry%threshold ) iflag = 1
      idx = maxloc( abs(array) )

    CASE (ABS_SMALLER_THAN)
      smry_type_char = 'ABS <'
      WHERE( abs(array) .lt. chunk_smry%threshold ) iflag = 1
      idx = minloc( abs(array) )

    END SELECT

    ! Total number of values exceeding threshold

    chunk_smry%count = sum( iflag )

    ! The extreme value

    chunk_smry%extreme_val  = array(idx(1))
    chunk_smry%extreme_col  =       idx(1)
    chunk_smry%extreme_lat  =   lat(idx(1))
    chunk_smry%extreme_lon  =   lon(idx(1))
  
    ! Send message to log file
  
    if (l_print_always) then
       write(iulog,"(a,i8,a,e15.7,a,e15.7)") &
       '**** '//trim(chunk_smry%field_name)//' from '//trim(chunk_smry%procedure_name)//':', &
       chunk_smry%count, ' values '//trim(smry_type_char), chunk_smry%threshold, &
       ', extreme is ', chunk_smry%extreme_val
    end if
  
  end subroutine get_chunk_smry_1_lev_real

  !---------------------------------------------------------------------------------------
  ! Description:
  !   Assuming the values exceeding a threshold have been identified for each chunk 
  !   handled by the current MPI process, this subroutine gets a total count of 
  !   violations in the current MPI process ("domain"), and identify the extreme value 
  !   among all chunks. Each call of this subroutine handles one field.
  !---------------------------------------------------------------------------------------
  subroutine get_domain_smry( chunk_smry_of_all_chunks, domain_smry )

    type(tp_stat_smry), intent(inout) :: chunk_smry_of_all_chunks(:)  ! shape: (nchnk)
    type(tp_stat_smry), intent(inout) :: domain_smry

    ! Local variables

    integer  :: idx(1)
    integer  :: ichnk
    character(len=shortchar) :: smry_type_char
  
    !------------------------------------------------
    ! Get a total count of values exceeding threshold

    domain_smry%count = sum( chunk_smry_of_all_chunks(:)%count )

    !-------------------
    ! Locate the extreme 
  
    SELECT CASE (domain_smry%smry_type)
    CASE (GREATER_EQ)
      idx = maxloc( chunk_smry_of_all_chunks(:)%extreme_val )
      smry_type_char = '>='

    CASE (ABS_GREATER_EQ)
      idx = maxloc( abs(chunk_smry_of_all_chunks(:)%extreme_val) )
      smry_type_char = 'ABS >='

    CASE (SMALLER_THAN)
      idx = minloc( chunk_smry_of_all_chunks(:)%extreme_val )
      smry_type_char = '<'

    CASE (ABS_SMALLER_THAN)
      idx = minloc( abs(chunk_smry_of_all_chunks(:)%extreme_val) )
      smry_type_char = 'ABS <'
    END SELECT

    ichnk = idx(1)
    domain_smry%extreme_val  = chunk_smry_of_all_chunks(ichnk)%extreme_val
    domain_smry%extreme_lat  = chunk_smry_of_all_chunks(ichnk)%extreme_lat
    domain_smry%extreme_lon  = chunk_smry_of_all_chunks(ichnk)%extreme_lon 
    domain_smry%extreme_lev  = chunk_smry_of_all_chunks(ichnk)%extreme_lev 
    domain_smry%extreme_col  = chunk_smry_of_all_chunks(ichnk)%extreme_col
    domain_smry%extreme_chnk = chunk_smry_of_all_chunks(ichnk)%extreme_chnk 

    ! Send message to log file
  
    if (l_print_always) then
       write(iulog,"(a,i8,a,e15.7,a,e15.7,a,3(a,i4),2(a,f8.2))") &
         trim(domain_smry%field_name)//' from '//trim(domain_smry%procedure_name)//':', &
         domain_smry%count, ' values '//trim(smry_type_char), domain_smry%threshold, &
         ', extreme is ',domain_smry%extreme_val,' at ',&
         '  chnk ',domain_smry%extreme_chnk, &
         ', col. ',domain_smry%extreme_col, &
         ', lev. ',domain_smry%extreme_lev, &
         ', lat = ',domain_smry%extreme_lat, &
         ', lon = ',domain_smry%extreme_lon 
    end if
  
  end subroutine get_domain_smry

  !------------------------------------------------------------------------------------------
  ! Description:
  !   Assuming the values exceeding a threshold have been identified for each chunk, 
  !   this subroutine first calls get_domain_smry to get a summary for each MPI process,
  !   then collect information from all MPI processes and get the global summary.
  !   This subroutine is meant to handle all fields on the summary list.
  !------------------------------------------------------------------------------------------
  subroutine get_global_smry( chunk_smry_2d, domain_smry_1d, nstep)

#ifdef SPMD
    use mpishorthand, only: mpir8, mpiint, mpicom
    use spmd_utils,   only: masterproc, npes
    use cam_logfile,  only: iulog
#endif

    type(tp_stat_smry)  ::  chunk_smry_2d(:,:)  ! shape: (nchunk,nfld)
    type(tp_stat_smry)  :: domain_smry_1d(:)    ! shape: (nfld)
    integer,intent(in)   :: nstep               ! model time step

    integer :: ii

#ifdef SPMD
    integer           :: sndrcvcnt   ! number of element for single send/receive
    integer,parameter :: nreal = 3   ! number of real variables to pack for mpigather
    integer,parameter :: nintg = 4   ! number of integer variables to pack for mpigather

    real(r8) :: real_array          (nreal,current_number_of_smry_fields)
    real(r8) :: real_array_gathered (nreal,current_number_of_smry_fields,0:npes-1)

    integer  :: intg_array          (nintg,current_number_of_smry_fields)
    integer  :: intg_array_gathered (nintg,current_number_of_smry_fields,0:npes-1)

    character(len=shortchar) :: smry_type_char
    integer  :: idx(1), ipe
#endif

    !--------------------------------------------
    ! Get domain summaries for each MPI process
    !--------------------------------------------
    do ii = 1,current_number_of_smry_fields
       call get_domain_smry( chunk_smry_2d(:,ii), domain_smry_1d(ii) ) !intent: 3xin, inout
    end do

    !--------------------------------------------
    ! Get global summaries
    !--------------------------------------------
#ifdef SPMD
    ! Pack arrays for MPI communication

    real_array(1,:) = domain_smry_1d(:)%extreme_val
    real_array(2,:) = domain_smry_1d(:)%extreme_lat
    real_array(3,:) = domain_smry_1d(:)%extreme_lon

    intg_array(1,:) = domain_smry_1d(:)%extreme_lev
    intg_array(2,:) = domain_smry_1d(:)%extreme_col
    intg_array(3,:) = domain_smry_1d(:)%extreme_chnk
    intg_array(4,:) = domain_smry_1d(:)%count

    ! Master process gathers info from all processes

    sndrcvcnt = nreal*current_number_of_smry_fields
    call mpigather( real_array, sndrcvcnt, mpir8,  real_array_gathered, sndrcvcnt, mpir8,  0, mpicom )

    sndrcvcnt = nintg*current_number_of_smry_fields
    call mpigather( intg_array, sndrcvcnt, mpiint, intg_array_gathered, sndrcvcnt, mpiint, 0, mpicom )

    ! Add the counts and locate the extreme values

    if (masterproc) then

      write(iulog,*) '**** Global summary at step ',nstep,' ****'
      do ii = 1,current_number_of_smry_fields

       SELECT CASE (global_smry_1d(ii)%smry_type)
       CASE( GREATER_EQ )
         smry_type_char = '>='
         idx = maxloc( real_array_gathered(1,ii,:) )

       CASE( ABS_GREATER_EQ)
         smry_type_char = 'ABS >='
         idx = maxloc( abs(real_array_gathered(1,ii,:)) )

       CASE( SMALLER_THAN)
         smry_type_char = '<'
         idx = minloc( real_array_gathered(1,ii,:) )

       CASE( ABS_SMALLER_THAN)
         smry_type_char = 'ABS <'
         idx = minloc( abs(real_array_gathered(1,ii,:)) )
       END SELECT

       ipe = idx(1)

       global_smry_1d(ii)%extreme_val = real_array_gathered(1,ii,ipe)
       global_smry_1d(ii)%extreme_lat = real_array_gathered(2,ii,ipe)
       global_smry_1d(ii)%extreme_lon = real_array_gathered(3,ii,ipe)

       global_smry_1d(ii)%extreme_lev  = intg_array_gathered(1,ii,ipe)
       global_smry_1d(ii)%extreme_col  = intg_array_gathered(2,ii,ipe)
       global_smry_1d(ii)%extreme_chnk = intg_array_gathered(3,ii,ipe)
       global_smry_1d(ii)%count        = sum(intg_array_gathered(4,ii,:))

       ! Send message to log file

       write(iulog,'(a,i8,a,e15.7,a,e15.7,a,3(a,i4),2(a,f8.2))')    &
             trim(global_smry_1d(ii)%field_name)//' from '//trim(global_smry_1d(ii)%procedure_name)//':', &
             global_smry_1d(ii)%count, ' values '//trim(smry_type_char), global_smry_1d(ii)%threshold, &
             ', extreme is ', global_smry_1d(ii)%extreme_val, ' at ',&
             '  chnk ',global_smry_1d(ii)%extreme_chnk, &
             ', col. ',global_smry_1d(ii)%extreme_col, &
             ', lev. ',global_smry_1d(ii)%extreme_lev, &
             ', lat = ',global_smry_1d(ii)%extreme_lat, &
             ', lon = ',global_smry_1d(ii)%extreme_lon

      end do
      write(iulog,*) '**** End of global summary at step ',nstep,' ****'

    end if
#endif

  end subroutine get_global_smry

end module global_summary
