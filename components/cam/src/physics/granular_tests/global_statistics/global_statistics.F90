module global_statistics

  use shr_kind_mod, only: r8=>shr_kind_r8

  implicit none
 !private
  public

! Public types:

! public tp_statistics

! Public procedures

! public stat_init
! public report_large_values
!

  integer, parameter :: LESS_THAN        = -1
  integer, parameter :: GREATER_THAN     =  1
  integer, parameter :: ABS_LESS_THAN    = -2
  integer, parameter :: ABS_GREATER_THAN =  2

!-------------------------------
  type tp_statistics

    character(len=128) :: procedure_name
    character(len=128) :: field_name

    integer :: stat_type

    real(r8) :: stat_value
    real(r8) :: threshold

    integer  :: count

    real(r8) :: extreme_lat
    real(r8) :: extreme_lon
    integer  :: extreme_chnk
    integer  :: extreme_col

  end type tp_statistics
!-------------------------------

contains

  subroutine get_chunk_stat( ncol, fldname, array, lat, lon, l_print_always,   &!intent(in)
                             domain_stat ) ! intent(inout)

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
    character(len=*), intent(in) :: fldname           ! field names
    real(r8),         intent(in) :: array(ncol)       ! array of values to be checked
                                                      ! occurrence will be reported
    real(r8),         intent(in) :: lat(ncol)
    real(r8),         intent(in) :: lon(ncol)
    logical,          intent(in) :: l_print_always    ! always print message in log file
                                                      ! (even when there are no
                                                      ! values exeeding tolerance)
    type(tp_statistics), intent(inout) :: domain_stat

    ! Local variables

    integer  :: iflag(ncol) 
    integer  :: icol(1)
    logical  :: l_print               ! print message in log file
  
    !--------------------------------------------------------------------------------------
    ! Calculate the total number of columns with value exceeding threshold
    ! then get the index of the extremem value.
  
    iflag(:) = 0

    SELECT CASE (domain_stat%stat_type)
    CASE (GREATER_THAN)
      where( array(:) .gt. domain_stat%threshold ) iflag = 1
      icol = maxloc( array )

    CASE (LESS_THAN)
      where( array(:) .lt. domain_stat%threshold ) iflag = 1
      icol = minloc( array )

    CASE (ABS_GREATER_THAN)
      WHERE( abs(array(:)) .gt. domain_stat%threshold ) iflag = 1
      icol = maxloc( abs(array) )

    CASE (ABS_LESS_THAN)
      WHERE( abs(array(:)) .lt. domain_stat%threshold ) iflag = 1
      icol = minloc( abs(array) )

    END SELECT

    ! Total number of values exceeding tolerance

    domain_stat%count = sum( iflag )

    ! The extreme value

    domain_stat%stat_value  = array(icol(1))
    domain_stat%extreme_col =       icol(1)
    domain_stat%extreme_lat =   lat(icol(1))
    domain_stat%extreme_lon =   lon(icol(1))
  
    ! Send message to log file
  
      l_print = l_print_always                                     &! always print
           .or. ( .not.l_print_always .and. (domain_stat%count>0) ) ! found large values
  
      if (l_print) then
         write(iulog,"(a,i8,a,e15.7,a,e15.7)") &
               "*** Procedure "//trim(domain_stat%procedure_name)//', field '//trim(domain_stat%field_name)//": ", &
               domain_stat%count, ' values exceeding ', &
               domain_stat%threshold, &
               ', extreme value is ', domain_stat%stat_value
         write(iulog,*)
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
  
    logical,          intent(in) :: l_print_always    ! always print message in log file
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
    CASE (LESS_THAN, ABS_LESS_THAN)
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
         write(iulog,*)
      end if
  
  end subroutine get_domain_stat

end module global_statistics
