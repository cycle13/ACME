module glc_elevclass_mod

  !---------------------------------------------------------------------
  !
  ! Purpose:
  !
  ! This module contains data and routines for operating on GLC elevation classes.

#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod
  use seq_comm_mct, only : logunit
  use shr_log_mod, only : errMsg => shr_log_errMsg

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: glc_elevclass_init            ! initialize GLC elevation class data
  public :: glc_elevclass_clean           ! deallocate memory allocated here
  public :: glc_get_num_elevation_classes ! get the number of elevation classes
  public :: glc_get_elevation_class       ! get the elevation class index for a given elevation
  public :: glc_get_elevclass_bounds      ! get the boundaries of all elevation classes
  public :: glc_mean_elevation_virtual    ! get the mean elevation of a virtual elevation class
  public :: glc_elevclass_as_string       ! returns a string corresponding to a given elevation class
  public :: glc_all_elevclass_strings     ! returns an array of strings for all elevation classes
  public :: glc_errcode_to_string         ! convert an error code into a string describing the error

  interface glc_elevclass_init
     module procedure glc_elevclass_init_default
     module procedure glc_elevclass_init_override
  end interface glc_elevclass_init


  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  ! Possible error code values
  integer, parameter, public :: GLC_ELEVCLASS_ERR_NONE = 0      ! err_code indicating no error
  integer, parameter, public :: GLC_ELEVCLASS_ERR_UNDEFINED = 1 ! err_code indicating elevation classes have not been defined
  integer, parameter, public :: GLC_ELEVCLASS_ERR_TOO_LOW = 2   ! err_code indicating topo below lowest elevation class
  integer, parameter, public :: GLC_ELEVCLASS_ERR_TOO_HIGH = 3  ! err_code indicating topo above highest elevation class

  ! String length for glc elevation classes represented as strings
  integer, parameter, public :: GLC_ELEVCLASS_STRLEN = 2

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! number of elevation classes
  integer :: glc_nec

  ! upper elevation limit of each class (m)
  ! indexing starts at 0, with topomax(0) giving the lower elevation limit of EC 1
  real(r8), allocatable :: topomax(:)


contains

  !-----------------------------------------------------------------------
  subroutine glc_elevclass_init_default(my_glc_nec)
    !
    ! !DESCRIPTION:
    ! Initialize GLC elevation class data to default boundaries, based on given glc_nec
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: my_glc_nec  ! number of GLC elevation classes
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_elevclass_init'
    !-----------------------------------------------------------------------

    glc_nec = my_glc_nec
    allocate(topomax(0:glc_nec))

    select case (glc_nec)
    case(0)
       ! do nothing
    case(1)
       topomax = [0._r8, 10000._r8]
    case(3)
       topomax = [0._r8,  1000._r8,  2000._r8, 10000._r8]
    case(5)
       topomax = [0._r8,   500._r8,  1000._r8,  1500._r8, 2000._r8, 10000._r8]
    case(10)
       topomax = [0._r8,   200._r8,   400._r8,   700._r8,  1000._r8,  1300._r8,  &
                          1600._r8,  2000._r8,  2500._r8,  3000._r8, 10000._r8]
    case(36)
       topomax = [  0._r8,   200._r8,   400._r8,   600._r8,   800._r8,  &
                 1000._r8,  1200._r8,  1400._r8,  1600._r8,  1800._r8,  &
                 2000._r8,  2200._r8,  2400._r8,  2600._r8,  2800._r8,  &
                 3000._r8,  3200._r8,  3400._r8,  3600._r8,  3800._r8,  &
                 4000._r8,  4200._r8,  4400._r8,  4600._r8,  4800._r8,  &
                 5000._r8,  5200._r8,  5400._r8,  5600._r8,  5800._r8,  &
                 6000._r8,  6200._r8,  6400._r8,  6600._r8,  6800._r8,  &
                 7000._r8, 10000._r8]
    case default
       write(logunit,*) subname,' ERROR: unknown glc_nec: ', glc_nec
       call shr_sys_abort(subname//' ERROR: unknown glc_nec')
    end select

  end subroutine glc_elevclass_init_default

  !-----------------------------------------------------------------------
  subroutine glc_elevclass_init_override(my_glc_nec, my_topomax)
    !
    ! !DESCRIPTION:
    ! Initialize GLC elevation class data to the given elevation class boundaries.
    !
    ! The input, my_topomax, should have (my_glc_nec + 1) elements.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: my_glc_nec      ! number of GLC elevation classes
    real(r8), intent(in) :: my_topomax(0:) ! elevation class boundaries (m)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_elevclass_init_override'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(my_topomax) == (/my_glc_nec/)), __FILE__, __LINE__)

    glc_nec = my_glc_nec
    allocate(topomax(0:glc_nec))
    topomax = my_topomax

  end subroutine glc_elevclass_init_override

  !-----------------------------------------------------------------------
  subroutine glc_elevclass_clean()
    !
    ! !DESCRIPTION:
    ! Deallocate memory allocated in this module

    character(len=*), parameter :: subname = 'glc_elevclass_clean'
    !-----------------------------------------------------------------------

    if (allocated(topomax)) then
       deallocate(topomax)
    end if
    glc_nec = 0

  end subroutine glc_elevclass_clean

  !-----------------------------------------------------------------------
  function glc_get_num_elevation_classes() result(num_elevation_classes)
    !
    ! !DESCRIPTION:
    ! Get the number of GLC elevation classes
    !
    ! !ARGUMENTS:
    integer :: num_elevation_classes  ! function result
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_get_num_elevation_classes'
    !-----------------------------------------------------------------------

    num_elevation_classes = glc_nec

  end function glc_get_num_elevation_classes

  !-----------------------------------------------------------------------
  subroutine glc_get_elevation_class(topo, elevation_class, err_code)
    !
    ! !DESCRIPTION:
    ! Get the elevation class index associated with a given topographic height.
    !
    ! The returned elevation_class will be between 1 and num_elevation_classes, if this
    ! topographic height is contained in an elevation class. In this case, err_code will
    ! be GLC_ELEVCLASS_ERR_NONE (no error).
    !
    ! If there are no elevation classes defined, the returned value will be 0, and
    ! err_code will be GLC_ELEVCLASS_ERR_UNDEFINED
    !
    ! If this topographic height is below the lowest elevation class, the returned value
    ! will be 1, and err_code will be GLC_ELEVCLASS_ERR_TOO_LOW.
    !
    ! If this topographic height is above the highest elevation class, the returned value
    ! will be (num_elevation_classes), and err_code will be GLC_ELEVCLASS_ERR_TOO_HIGH.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: topo ! topographic height (m)
    integer, intent(out) :: elevation_class  ! elevation class index
    integer, intent(out) :: err_code ! error code (see above for possible codes)
    !
    ! !LOCAL VARIABLES:
    integer :: ec  ! temporary elevation class

    character(len=*), parameter :: subname = 'glc_get_elevation_class'
    !-----------------------------------------------------------------------

    if (glc_nec < 1) then
       elevation_class = 0
       err_code = GLC_ELEVCLASS_ERR_UNDEFINED
    else if (topo < topomax(0)) then
       elevation_class = 1
       err_code = GLC_ELEVCLASS_ERR_TOO_LOW
    else if (topo >= topomax(glc_nec)) then
       elevation_class = glc_nec
       err_code = GLC_ELEVCLASS_ERR_TOO_HIGH
    else
       err_code = GLC_ELEVCLASS_ERR_NONE
       elevation_class = 0
       do ec = 1, glc_nec
          if (topo >= topomax(ec - 1) .and. topo < topomax(ec)) then
             elevation_class = ec
             exit
          end if
       end do

       SHR_ASSERT(elevation_class > 0, subname//' elevation class was not assigned')
    end if

  end subroutine glc_get_elevation_class

  !-----------------------------------------------------------------------
  function glc_get_elevclass_bounds() result(elevclass_bounds)
    !
    ! !DESCRIPTION:
    ! Get the boundaries of all elevation classes.
    !
    ! This returns an array of size glc_nec+1, since it contains both the lower and upper
    ! bounds of each elevation class.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: elevclass_bounds(0:glc_nec)  ! function result
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_get_elevclass_bounds'
    !-----------------------------------------------------------------------

    elevclass_bounds(:) = topomax(:)

  end function glc_get_elevclass_bounds


  !-----------------------------------------------------------------------
  function glc_mean_elevation_virtual(elevation_class) result(mean_elevation)
    !
    ! !DESCRIPTION:
    ! Returns the mean elevation of a virtual elevation class
    !
    ! !ARGUMENTS:
    real(r8) :: mean_elevation  ! function result
    integer, intent(in) :: elevation_class
    !
    ! !LOCAL VARIABLES:
    integer :: resulting_elevation_class
    integer :: err_code

    character(len=*), parameter :: subname = 'glc_mean_elevation_virtual'
    !-----------------------------------------------------------------------

    if (elevation_class == 0) then
       ! Bare land "elevation class"
       mean_elevation = 0._r8
    else
       if (elevation_class < glc_nec) then
          ! Normal case
          mean_elevation = (topomax(elevation_class - 1) + topomax(elevation_class)) / 2._r8
       else if (elevation_class == glc_nec) then
          ! In the top elevation class; in this case, assignment of a "mean" elevation is
          ! somewhat arbitrary (because we expect the upper bound of the top elevation
          ! class to be very high).

          if (glc_nec > 1) then
             mean_elevation = 2._r8 * topomax(elevation_class - 1) - topomax(elevation_class - 2)
          else
             ! entirely arbitrary
             mean_elevation = 1000._r8
          end if
       else
          write(logunit,*) subname,' ERROR: elevation class out of bounds: ', elevation_class
          call shr_sys_abort(subname // ' ERROR: elevation class out of bounds')
       end if
    end if

    ! Ensure that the resulting elevation is within the given elevation class
    if (elevation_class > 0) then
       call glc_get_elevation_class(mean_elevation, resulting_elevation_class, err_code)
       if (err_code /= GLC_ELEVCLASS_ERR_NONE) then
          write(logunit,*) subname, ' ERROR: generated elevation that results in an error'
          write(logunit,*) 'when trying to determine the resulting elevation class'
          write(logunit,*) glc_errcode_to_string(err_code)
          write(logunit,*) 'elevation_class, mean_elevation = ', elevation_class, mean_elevation
          call shr_sys_abort(subname // ' ERROR: generated elevation that results in an error')
       else if (resulting_elevation_class /= elevation_class) then
          write(logunit,*) subname, ' ERROR: generated elevation outside the given elevation class'
          write(logunit,*) 'elevation_class, mean_elevation, resulting_elevation_class = ', &
               elevation_class, mean_elevation, resulting_elevation_class
          call shr_sys_abort(subname // ' ERROR: generated elevation outside the given elevation class')
       end if
    end if

  end function glc_mean_elevation_virtual


  !-----------------------------------------------------------------------
  function glc_elevclass_as_string(elevation_class) result(ec_string)
    !
    ! !DESCRIPTION:
    ! Returns a string corresponding to a given elevation class.
    !
    ! This string can be used as a suffix for fields in MCT attribute vectors.
    !
    ! ! NOTE(wjs, 2015-01-19) This function doesn't fully belong in this module, since it
    ! doesn't refer to the data stored in this module. However, I can't think of a more
    ! appropriate place for it.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=GLC_ELEVCLASS_STRLEN) :: ec_string  ! function result
    integer, intent(in) :: elevation_class
    !
    ! !LOCAL VARIABLES:
    character(len=16) :: format_string

    character(len=*), parameter :: subname = 'glc_elevclass_as_string'
    !-----------------------------------------------------------------------

    ! e.g., for GLC_ELEVCLASS_STRLEN = 2, format_string will be '(i2.2)'
    write(format_string,'(a,i0,a,i0,a)') '(i', GLC_ELEVCLASS_STRLEN, '.', GLC_ELEVCLASS_STRLEN, ')'

    write(ec_string,trim(format_string)) elevation_class
  end function glc_elevclass_as_string

  !-----------------------------------------------------------------------
  function glc_all_elevclass_strings(include_zero) result(ec_strings)
    !
    ! !DESCRIPTION:
    ! Returns an array of strings corresponding to all elevation classes from 1 to glc_nec
    !
    ! If include_zero is present and true, then includes elevation class 0 - so goes from
    ! 0 to glc_nec
    !
    ! These strings can be used as suffixes for fields in MCT attribute vectors.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=GLC_ELEVCLASS_STRLEN), allocatable :: ec_strings(:)  ! function result
    logical, intent(in), optional :: include_zero   ! if present and true, include elevation class 0 (default is false)
    !
    ! !LOCAL VARIABLES:
    logical :: l_include_zero  ! local version of optional include_zero argument
    integer :: lower_bound
    integer :: i

    character(len=*), parameter :: subname = 'glc_all_elevclass_strings'
    !-----------------------------------------------------------------------

    if (present(include_zero)) then
       l_include_zero = include_zero
    else
       l_include_zero = .false.
    end if

    if (l_include_zero) then
       lower_bound = 0
    else
       lower_bound = 1
    end if

    allocate(ec_strings(lower_bound:glc_nec))
    do i = lower_bound, glc_nec
       ec_strings(i) = glc_elevclass_as_string(i)
    end do

  end function glc_all_elevclass_strings


  !-----------------------------------------------------------------------
  function glc_errcode_to_string(err_code) result(err_string)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=256) :: err_string  ! function result
    integer, intent(in) :: err_code   ! error code (one of the GLC_ELEVCLASS_ERR* values)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'glc_errcode_to_string'
    !-----------------------------------------------------------------------

    select case (err_code)
    case (GLC_ELEVCLASS_ERR_NONE)
       err_string = '(no error)'
    case (GLC_ELEVCLASS_ERR_UNDEFINED)
       err_string = 'Elevation classes have not yet been defined'
    case (GLC_ELEVCLASS_ERR_TOO_LOW)
       err_string = 'Topographic height below the lower bound of the lowest elevation class'
    case (GLC_ELEVCLASS_ERR_TOO_HIGH)
       err_string = 'Topographic height above the upper bound of the highest elevation class'
    case default
       err_string = 'UNKNOWN ERROR'
    end select

  end function glc_errcode_to_string


end module glc_elevclass_mod
