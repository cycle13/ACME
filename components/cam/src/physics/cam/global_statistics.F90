module global_statistics

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none
  private

  public tp_statistics
  public global_stat_init
  public add_stat_field
  public get_stat_field_idx
  public get_chunk_stat
  public get_domain_stat

  !------------------
  ! Constants

  integer,public,parameter :: SMALLER_THAN     = -1
  integer,public,parameter :: GREATER_THAN     =  1
  integer,public,parameter :: ABS_SMALLER_THAN = -2
  integer,public,parameter :: ABS_GREATER_THAN =  2

  character(len=128),private,parameter :: THIS_MODULE = 'global_statistics'

  !-------------------------------
  ! Derived type

  type tp_statistics

    character(len=128) :: procedure_name
    character(len=128) :: field_name

    integer  :: stat_type
    real(r8) :: threshold
    integer  :: count = 0

    real(r8) :: extreme_val
    real(r8) :: extreme_lat  = -999._r8
    real(r8) :: extreme_lon  = -999._r8
    integer  :: extreme_chnk = -999
    integer  :: extreme_col  = -999
    integer  :: extreme_lev  = -999

  end type tp_statistics
  !-------------------------------

  interface get_chunk_stat
    module procedure get_chunk_stat_1d_real
    module procedure get_chunk_stat_2d_real
  end interface get_chunk_stat

  integer,parameter  :: max_number_of_stat_fields = 1000
  type(tp_statistics) :: global_stat(max_number_of_stat_fields)
  integer,public :: current_number_of_stat_fields = 0

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
       write(msg,*) 'global_statistics: domain_stat allocation error = ',ierr
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
       write(msg,*) 'global_statistics: chunk_stat allocation error = ',ierr
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

       chunk_stat(ichnk,1:current_number_of_stat_fields)%extreme_chnk = ichnk
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
  !---------------------------------------------------------------------------------------
  ! Description:
  !   Subroutine for finding values in array that exceed threshold.
  !   The total number of such values and the location of the extreme values are saved 
  !   for later use.
  !
  ! Written by:
  !   Hui Wan (PNNL, 2017-05)
  !---------------------------------------------------------------------------------------
  subroutine get_chunk_stat_2d_real( ncol, pver, array, lat, lon, l_print_always,   &!intent(in)
                                     chunk_stat ) ! intent(inout)


    use shr_kind_mod,  only: r8=>shr_kind_r8
    use cam_logfile,   only: iulog
  
    implicit none
  
    integer,          intent(in) :: ncol              ! number of columns packed in array
    integer,          intent(in) :: pver              ! number of vertical levels
    real(r8),         intent(in) :: array(ncol,pver)  ! array of values to be checked
                                                      ! occurrence will be reported
    real(r8),         intent(in) :: lat(ncol)
    real(r8),         intent(in) :: lon(ncol)
    logical,          intent(in) :: l_print_always    ! always print message in log file
                                                      ! (even when there are no
                                                      ! values exeeding threshold)
    type(tp_statistics), intent(inout) :: chunk_stat

    ! Local variables

    integer  :: iflag(ncol,pver) 
    integer  :: idx(2)
    logical  :: l_print               ! print message in log file
    character(len=16) :: stat_type_char
  
    !--------------------------------------------------------------------------------------
    ! Calculate the total number of columns with value exceeding threshold
    ! then get the index of the extremem value.
  
    iflag(:,:) = 0

    SELECT CASE (chunk_stat%stat_type)
    CASE (GREATER_THAN)
      stat_type_char = '>='
      where( array .ge. chunk_stat%threshold ) iflag = 1
      idx = maxloc( array )

    CASE (SMALLER_THAN)
      stat_type_char = '<'
      where( array .lt. chunk_stat%threshold ) iflag = 1
      idx = minloc( array )

    CASE (ABS_GREATER_THAN)
      stat_type_char = 'ABS >='
      WHERE( abs(array) .ge. chunk_stat%threshold ) iflag = 1
      idx = maxloc( abs(array) )

    CASE (ABS_SMALLER_THAN)
      stat_type_char = 'ABS < '
      WHERE( abs(array) .lt. chunk_stat%threshold ) iflag = 1
      idx = minloc( abs(array) )

    END SELECT

    ! Total number of values exceeding threshold

    chunk_stat%count = sum( iflag )

    ! The extreme value

    chunk_stat%extreme_val  = array(idx(1),idx(2))
    chunk_stat%extreme_col  =       idx(1)
    chunk_stat%extreme_lev  =       idx(2)
    chunk_stat%extreme_lat  =   lat(idx(1))
    chunk_stat%extreme_lon  =   lon(idx(1))
  
    ! Send message to log file
  
      l_print = l_print_always                                     &! always print
           .or. ( .not.l_print_always .and. (chunk_stat%count>0) ) ! found large values
  
      if (l_print) then
         write(iulog,"(a,i8,a,e15.7,a,e15.7)") &
               "*** Procedure "//trim(chunk_stat%procedure_name)// &
               ', field '//trim(chunk_stat%field_name)//": ", &
               chunk_stat%count, ' values '//trim(stat_type_char), &
               chunk_stat%threshold, &
               ', extreme value is ', chunk_stat%extreme_val
      end if
  
  end subroutine get_chunk_stat_2d_real
  !-------------

  subroutine get_chunk_stat_1d_real( ncol, array, lat, lon, l_print_always,   &!intent(in)
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
  
    integer,          intent(in) :: ncol              ! number of columns packed in array
    real(r8),         intent(in) :: array(ncol)       ! array of values to be checked
                                                      ! occurrence will be reported
    real(r8),         intent(in) :: lat(ncol)
    real(r8),         intent(in) :: lon(ncol)
    logical,          intent(in) :: l_print_always    ! always print message in log file
                                                      ! (even when there are no
                                                      ! values exeeding threshold)
    type(tp_statistics), intent(inout) :: chunk_stat

    ! Local variables

    integer  :: iflag(ncol) 
    integer  :: icol(1)
    logical  :: l_print               ! print message in log file
    character(len=16) :: stat_type_char
  
    !--------------------------------------------------------------------------------------
    ! Calculate the total number of columns with value exceeding threshold
    ! then get the index of the extremem value.
  
    iflag(:) = 0

    SELECT CASE (chunk_stat%stat_type)
    CASE (GREATER_THAN)
      stat_type_char = '>='
      where( array .ge. chunk_stat%threshold ) iflag = 1
      icol = maxloc( array )

    CASE (SMALLER_THAN)
      stat_type_char = '<'
      where( array .lt. chunk_stat%threshold ) iflag = 1
      icol = minloc( array )

    CASE (ABS_GREATER_THAN)
      stat_type_char = 'ABS >='
      WHERE( abs(array) .ge. chunk_stat%threshold ) iflag = 1
      icol = maxloc( abs(array) )

    CASE (ABS_SMALLER_THAN)
      stat_type_char = 'ABS <'
      WHERE( abs(array) .lt. chunk_stat%threshold ) iflag = 1
      icol = minloc( abs(array) )

    END SELECT

    ! Total number of values exceeding threshold

    chunk_stat%count = sum( iflag )

    ! The extreme value

    chunk_stat%extreme_val  = array(icol(1))
    chunk_stat%extreme_col  =       icol(1)
    chunk_stat%extreme_lat  =   lat(icol(1))
    chunk_stat%extreme_lon  =   lon(icol(1))
  
    ! Send message to log file
  
      l_print = l_print_always                                     &! always print
           .or. ( .not.l_print_always .and. (chunk_stat%count>0) ) ! found large values
  
      if (l_print) then
         write(iulog,"(a,i8,a,e15.7,a,e15.7)") &
               "*** Procedure "//trim(chunk_stat%procedure_name)// &
               ', field '//trim(chunk_stat%field_name)//": ", &
               chunk_stat%count, ' values '//trim(stat_type_char), &
               chunk_stat%threshold, &
               ', extreme value is ', chunk_stat%extreme_val
      end if
  
  end subroutine get_chunk_stat_1d_real


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
                                             ! values exeeding threshold)
    integer, intent(in) :: nchnk

    type(tp_statistics), intent(inout) ::  chunk_stat(nchnk)
    type(tp_statistics), intent(inout) :: domain_stat

    ! Local variables

    integer  :: iflag(nchnk) 
    integer  :: idx(1)
    integer  :: ichnk
    logical  :: l_print                   ! print message in log file
    character(len=16) :: stat_type_char
  
    !--------------------------------------------------------------------------------------
    ! Calculate the total number of columns with value exceeding threshold
    ! then get the index of the extremem value.
  
    iflag(:) = 0

    SELECT CASE (domain_stat%stat_type)
    CASE (GREATER_THAN)
      idx = maxloc( chunk_stat(:)%extreme_val )
      stat_type_char = '>='

    CASE (ABS_GREATER_THAN)
      idx = maxloc( chunk_stat(:)%extreme_val )
      stat_type_char = 'ABS >='

    CASE (SMALLER_THAN)
      idx = minloc( chunk_stat(:)%extreme_val )
      stat_type_char = '<'

    CASE (ABS_SMALLER_THAN)
      idx = minloc( chunk_stat(:)%extreme_val )
      stat_type_char = 'ABS <'
    END SELECT

    ! Total number of values exceeding threshold

    domain_stat%count = sum( chunk_stat(:)%count )

    ! The extreme value

    ichnk = idx(1)
    domain_stat%extreme_val  = chunk_stat(ichnk)%extreme_val
    domain_stat%extreme_lat  = chunk_stat(ichnk)%extreme_lat
    domain_stat%extreme_lon  = chunk_stat(ichnk)%extreme_lon 
    domain_stat%extreme_lev  = chunk_stat(ichnk)%extreme_lev 
    domain_stat%extreme_col  = chunk_stat(ichnk)%extreme_col
    domain_stat%extreme_chnk = chunk_stat(ichnk)%extreme_chnk 

    ! Send message to log file
  
      l_print = l_print_always                                     &! always print
           .or. ( .not.l_print_always .and. (domain_stat%count>0) ) ! found large values
  
      if (l_print) then
         write(iulog,"(a,i8,a,e15.7,a,e15.7,a,3(a,i4),2(a,f8.2))") &
               "*** Procedure "//trim(domain_stat%procedure_name)// &
               ', field '//trim(domain_stat%field_name)//": ", &
               domain_stat%count, ' values '//trim(stat_type_char), &
               domain_stat%threshold, &
               ', extreme value is ', domain_stat%extreme_val, &
               ' at ',&
               '  chnk ',domain_stat%extreme_chnk, &
               ', col. ',domain_stat%extreme_col, &
               ', lev. ',domain_stat%extreme_lev, &
               ', lat = ',domain_stat%extreme_lat, &
               ', lon = ',domain_stat%extreme_lon 
      end if
  
  end subroutine get_domain_stat

  subroutine get_global_stat( chunk_stat_2d, domain_stat, begchunk, endchunk )

#ifdef SPMD
    use mpishortland
#endif

    type(tp_statistics), pointer ::  chunk_stat_2d(:,:)
    type(tp_statistics), pointer :: domain_stat(:)
    integer,intent(in) :: begchunk, endchunk

    integer :: ii, nchnk
    logical :: l_print_always = .true.

#ifdef SPMD
    integer  :: scount
    real(r8) :: real_array         (3,current_number_of_stat_fields)
    real(r8) :: real_array_gathered(3,current_number_of_stat_fields,0:npes-1)
    integer  ::  int_array         (4,current_number_of_stat_fields)
    integer  ::  int_array_gathered(4,current_number_of_stat_fields,0:npes-1)
#endif

    nchnk = endchunk - begchunk + 1

    do ii = 1,current_number_of_stat_fields
       call get_domain_stat( l_print_always, nchnk, chunk_stat_2d(:,ii), domain_stat(ii) ) !intent: 3xin, inout
    end do

#ifdef SPMD
    real_array(1,:) = domain_stat(:)%extreme_val
    real_array(2,:) = domain_stat(:)%extreme_lat
    real_array(3,:) = domain_stat(:)%extreme_lon

     int_array(1,:) = domain_stat(:)%extreme_lev
     int_array(2,:) = domain_stat(:)%extreme_col
     int_array(3,:) = domain_stat(:)%extreme_chnk
     int_array(4,:) = domain_stat(:)%count

    scount = 3*current_number_of_stat_fields
    call MPI_GATHER(real_array, scount, mpir8,  real_array_gathered, scount, mpir8,  root, comm, ierr)

    scount = 4*current_number_of_stat_fields
    call MPI_GATHER( int_array, scount, mpiint,  int_array_gathered, scount, mpiint, root, comm, ierr)

    if (masterproc) then
    do ii = 1,current_number_of_stat_fields

       SELECT CASE (global_stat(ii)%stat_type)
       CASE( GREATER_THAN, ABS_GREATHER_THAN)
         idx = maxloc( real_array_gathered(1,ii,:) )
       CASE( SMALLER_THAN, ABS_SMALLER_THAN)
         idx = minloc( real_array_gathered(1,ii,:) )
       END SELECT

       global_stat(ii)%extreme_val = real_array_gathered(1,ii,idx(1))
       global_stat(ii)%extreme_lat = real_array_gathered(2,ii,idx(1))
       global_stat(ii)%extreme_lon = real_array_gathered(3,ii,idx(1))

       global_stat(ii)%extreme_lev  = int_array_gathered(1,ii,idx(1))
       global_stat(ii)%extreme_col  = int_array_gathered(2,ii,idx(1))
       global_stat(ii)%extreme_chnk = int_array_gathered(3,ii,idx(1))

       global_stat(ii)%count = sum(int_array_gathered(4,ii,:))

       write(iulog,*) 
    end do
    end if
#endif

  end subroutine get_global_stat

end module global_statistics
