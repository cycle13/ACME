module global_statistics

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none
  private

  public SMALLER_THAN, GREATER_THAN, ABS_SMALLER_THAN, ABS_GREATER_THAN
  public tp_statistics
  public global_stat_init
  public add_stat_field
  public get_stat_field_idx
  public get_chunk_stat
  public get_domain_stat

  character(len=128), parameter :: THIS_MODULE = 'global_statistics'

  integer, parameter :: SMALLER_THAN     = -1
  integer, parameter :: GREATER_THAN     =  1
  integer, parameter :: ABS_SMALLER_THAN = -2
  integer, parameter :: ABS_GREATER_THAN =  2

!-------------------------------
  type tp_statistics

    character(len=128) :: procedure_name
    character(len=128) :: field_name

    integer  :: stat_type
    real(r8) :: threshold
    real(r8) :: stat_value

    integer  :: count = 0

    real(r8) :: extreme_lat  = -999._r8
    real(r8) :: extreme_lon  = -999._r8
    integer  :: extreme_chnk = -999
    integer  :: extreme_col  = -999

  end type tp_statistics
!-------------------------------

  integer,parameter  :: max_number_of_stat_fields = 1000
  type(tp_statistics) :: global_stat(max_number_of_stat_fields)
  integer :: current_number_of_stat_fields = 0

  character(len=256) :: msg 

contains

  subroutine add_stat_field( fldname, procname, stattype, threshold, fldidx )

    character(len=*), intent(in)   :: fldname
    character(len=*), intent(in)   :: procname
    real(r8)                       :: threshold
    integer                        :: stattype
    integer, intent(out), optional :: fldidx

    integer :: ii

    ! First check if the same field from the same procedure has already been registered.
    do ii = 1,current_number_of_stat_fields
       if (trim(global_stat(ii)%field_name) == trim(fldname) .and. &
           trim(global_stat(ii)%procedure_name) == trim(procname)  ) then
           msg = trim(THIS_MODULE)//': field '//trim(fldname)//' from procedure '// &
                 trim(procname)//' has already been added to list of statistics.'
           call endrun(trim(msg))
       end if
    end do

    ! Now add a new field

    current_number_of_stat_fields = current_number_of_stat_fields + 1

    if (current_number_of_stat_fields.gt.max_number_of_stat_fields) then
       msg = trim(THIS_MODULE)//': max. No. of fields exceeded when attempting to add '//trim(fldname)
       call endrun(trim(msg))
    end if

    ii = current_number_of_stat_fields
    global_stat(ii)%field_name     = trim(fldname)
    global_stat(ii)%procedure_name = trim(procname)
    global_stat(ii)%threshold      = threshold
    global_stat(ii)%stat_type      = stattype

    if (present(fldidx)) fldidx = ii 
 
  end subroutine add_stat_field

  !---------------------------------------------------------------------
  subroutine global_stat_init( chunk_stat, domain_stat, begchunk, endchunk )

    type(tp_statistics), pointer ::  chunk_stat(:,:)
    type(tp_statistics), pointer :: domain_stat(:)
    integer, intent(in) :: begchunk, endchunk

    integer :: ierr, ichnk

    ! Initialize domain_stat (1D array):
    ! Allocate memory

    allocate(domain_stat(current_number_of_stat_fields), stat=ierr)
    if( ierr /= 0 ) then
       write(msg,*) 'global_statistics: global_stat_init allocation error = ',ierr
       call endrun(trim(msg))
    end if

    ! Copy metadata

    domain_stat(1:current_number_of_stat_fields)%field_name = &
    global_stat(1:current_number_of_stat_fields)%field_name

    domain_stat(1:current_number_of_stat_fields)%procedure_name = &
    global_stat(1:current_number_of_stat_fields)%procedure_name

    domain_stat(1:current_number_of_stat_fields)%stat_type = &
    global_stat(1:current_number_of_stat_fields)%stat_type

    domain_stat(1:current_number_of_stat_fields)%threshold = &
    global_stat(1:current_number_of_stat_fields)%threshold

    ! Initialize chunk_stat (2D array):
    ! Allocate memory

    allocate(chunk_stat(begchunk:endchunk,current_number_of_stat_fields), stat=ierr)
    if( ierr /= 0 ) then
       write(msg,*) 'global_statistics: global_stat_init allocation error = ',ierr
       call endrun(trim(msg))
    end if

    ! Copy metadata

    do ichnk = begchunk,endchunk 

       chunk_stat(ichnk,1:current_number_of_stat_fields)%field_name = &
      global_stat(      1:current_number_of_stat_fields)%field_name

       chunk_stat(ichnk,1:current_number_of_stat_fields)%procedure_name = &
      global_stat(      1:current_number_of_stat_fields)%procedure_name

       chunk_stat(ichnk,1:current_number_of_stat_fields)%stat_type = &
      global_stat(      1:current_number_of_stat_fields)%stat_type

       chunk_stat(ichnk,1:current_number_of_stat_fields)%threshold = &
      global_stat(      1:current_number_of_stat_fields)%threshold
    
    end do

  end subroutine global_stat_init

  !---------------------------------------------------------------------
  subroutine get_stat_field_idx( fldname, procname, fldidx )

    character(len=*), intent(in)   :: fldname
    character(len=*), intent(in)   :: procname
    integer, intent(out)           :: fldidx

    integer :: ii
    logical :: found_field

    found_field = .false.
    do ii = 1,current_number_of_stat_fields 
       if (trim(global_stat(ii)%field_name) == trim(fldname) .and. &
           trim(global_stat(ii)%procedure_name) == trim(procname)  ) then
          found_field = .true.
          fldidx = ii
          exit
       end if
    end do

    if (.not.found_field) then
       write(msg,*) 'global_statistics:get_stat_field_idx, did not find '// &
                     trim(fldname)//' from '//trim(procname)//' on the list.'
       call endrun(trim(msg))
    end if
 
  end subroutine get_stat_field_idx

  !---------------------------------------------------------------------
  subroutine get_chunk_stat( ichnk, ncol, array, lat, lon, l_print_always,   &!intent(in)
                             chunk_stat ) ! intent(inout)

  !---------------------------------------------------------------------------------------
  ! Description:
  !
  ! Written by:
  !   Hui Wan (PNNL, 2017-05)
  !---------------------------------------------------------------------------------------
    use shr_kind_mod,  only: r8=>shr_kind_r8
    use cam_logfile,   only: iulog
  
    implicit none
  
    integer,          intent(in) :: ichnk 
    integer,          intent(in) :: ncol              ! number of columns packed in array
    real(r8),         intent(in) :: array(ncol)       ! array of values to be checked
                                                      ! occurrence will be reported
    real(r8),         intent(in) :: lat(ncol)
    real(r8),         intent(in) :: lon(ncol)
    logical,          intent(in) :: l_print_always    ! always print message in log file
                                                      ! (even when there are no
                                                      ! values exeeding tolerance)
    type(tp_statistics), intent(inout) :: chunk_stat

    ! Local variables

    integer  :: iflag(ncol) 
    integer  :: icol(1)
    logical  :: l_print               ! print message in log file
  
    !--------------------------------------------------------------------------------------
    ! Calculate the total number of columns with value exceeding threshold
    ! then get the index of the extremem value.
  
    iflag(:) = 0

    SELECT CASE (chunk_stat%stat_type)
    CASE (GREATER_THAN)
      where( array(:) .ge. chunk_stat%threshold ) iflag = 1
      icol = maxloc( array )

    CASE (SMALLER_THAN)
      where( array(:) .lt. chunk_stat%threshold ) iflag = 1
      icol = minloc( array )

    CASE (ABS_GREATER_THAN)
      WHERE( abs(array(:)) .ge. chunk_stat%threshold ) iflag = 1
      icol = maxloc( abs(array) )

    CASE (ABS_SMALLER_THAN)
      WHERE( abs(array(:)) .lt. chunk_stat%threshold ) iflag = 1
      icol = minloc( abs(array) )

    END SELECT

    ! Total number of values exceeding tolerance

    chunk_stat%count = sum( iflag )

    ! The extreme value

    chunk_stat%stat_value   = array(icol(1))
    chunk_stat%extreme_col  =       icol(1)
    chunk_stat%extreme_lat  =   lat(icol(1))
    chunk_stat%extreme_lon  =   lon(icol(1))
    chunk_stat%extreme_chnk =  ichnk
  
    ! Send message to log file
  
      l_print = l_print_always                                     &! always print
           .or. ( .not.l_print_always .and. (chunk_stat%count>0) ) ! found large values
  
      if (l_print) then
         write(iulog,"(a,i8,a,e15.7,a,e15.7)") &
               "*** Procedure "//trim(chunk_stat%procedure_name)//', field '//trim(chunk_stat%field_name)//": ", &
               chunk_stat%count, ' values exceeding ', &
               chunk_stat%threshold, &
               ', extreme value is ', chunk_stat%stat_value
      end if
  
  end subroutine get_chunk_stat


  subroutine get_domain_stat( l_print_always, nchnk, chunk_stat, domain_stat ) ! in, inout, in

  !---------------------------------------------------------------------------------------
  ! Description:
  !
  ! Written by:
  !   Hui Wan (PNNL, 2017-05)
  !---------------------------------------------------------------------------------------
    use shr_kind_mod,  only: r8=>shr_kind_r8
    use cam_logfile,   only: iulog
  
    implicit none
  
    logical, intent(in) :: l_print_always    ! always print message in log file
                                             ! (even when there are no
                                             ! values exeeding tolerance)
    integer, intent(in) :: nchnk

    type(tp_statistics), intent(inout) ::  chunk_stat(nchnk)
    type(tp_statistics), intent(inout) :: domain_stat

    ! Local variables

    integer  :: iflag(nchnk) 
    integer  :: icol(1)
    logical  :: l_print               ! print message in log file
  
    !--------------------------------------------------------------------------------------
    ! Calculate the total number of columns with value exceeding threshold
    ! then get the index of the extremem value.
  
    iflag(:) = 0

    SELECT CASE (domain_stat%stat_type)
    CASE (GREATER_THAN, ABS_GREATER_THAN)
      icol = maxloc( chunk_stat(:)%stat_value )
    CASE (SMALLER_THAN, ABS_SMALLER_THAN)
      icol = minloc( chunk_stat(:)%stat_value )
    END SELECT

    ! Total number of values exceeding tolerance

    domain_stat%count = sum( chunk_stat(:)%count )

    ! The extreme value

    domain_stat%stat_value  = chunk_stat(icol(1))%stat_value
    domain_stat%extreme_col = chunk_stat(icol(1))%extreme_col 
    domain_stat%extreme_lat = chunk_stat(icol(1))%extreme_lat
    domain_stat%extreme_lon = chunk_stat(icol(1))%extreme_lon 
  
    ! Send message to log file
  
      l_print = l_print_always                                     &! always print
           .or. ( .not.l_print_always .and. (domain_stat%count>0) ) ! found large values
  
      if (l_print) then
         write(iulog,"(a,i8,a,e15.7,a,e15.7)") &
               "*** Procedure "//trim(domain_stat%procedure_name)//', field '//trim(domain_stat%field_name)//": ", &
               domain_stat%count, ' values exceeding ', &
               domain_stat%threshold, &
               ', extreme value is ', domain_stat%stat_value
      end if
  
  end subroutine get_domain_stat

end module global_statistics
