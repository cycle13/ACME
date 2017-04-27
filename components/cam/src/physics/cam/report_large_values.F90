subroutine report_large_values( ncol, nfld, fldname, array, tol, string )

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use cam_logfile,   only: iulog

  implicit none

  integer,          intent(in) :: ncol              ! number of columns packed in array
  integer,          intent(in) :: nfld              ! number of fields  packed in array
  character(len=*), intent(in) :: fldname(nfld)     ! field names
  real(r8),         intent(in) :: array(ncol,nfld)  ! array of values to be checked
  real(r8),         intent(in) :: tol               ! tolerance beyond which the number of occurance 
                                                ! will be reported
  character(len=*), intent(in) :: string

  ! Local variables

  integer  :: ifld                  ! loop index
  integer  :: iflag(ncol)           ! for flagging large values. 1 = above tol, 0 = below tol
  integer  :: n_large_local(nfld)   ! number of values above tol on this process
  real(r8) :: z_largest_local(nfld) ! largest value 

  ! calculate the total number of columns at which the value exceeds tolerance

  do ifld = 1,nfld
    iflag(:) = 0
    where( abs(array(:,ifld)) > tol ) iflag = 1
    n_large_local(ifld) = sum( iflag )
  end do

  ! get the largest values

  z_largest_local(:) = maxval( abs(array), dim=1 )


  ! send message to log file

  do ifld = 1,nfld
    write(iulog,"(a,i8,a,e15.7,a,e15.7)") string//', subr. report_large_values, field '//fldname(ifld), &
                                         n_large_local(ifld), ' values exceeding ', tol, &
                                         ', largest was ', z_largest_local(ifld)
  end do

end subroutine report_large_values
