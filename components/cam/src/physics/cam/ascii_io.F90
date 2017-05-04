module ascii_io

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use constituents,  only: pcnst, cnst_name 
  use ppgrid,        only: pcols, pver
  use physics_types, only: physics_state
  use cam_abortutils,only: endrun
  use cam_logfile,   only : iulog

  implicit none

CONTAINS

!---------------------------------------------------------------------------------
! Write out components of model state vector to ASCII file
!---------------------------------------------------------------------------------
subroutine write_state_ascii( nstep, state )

   integer, intent(in) :: nstep
   type(physics_state), intent(in) :: state       ! Physics state variables

   integer :: ncol, lchnk, ii, kk , m
   integer :: funit, ierr
   character(len=128) :: filename

   !--------
   ! chunk index will be part of the output file name

   lchnk = state%lchnk

   ! output will include all active columns in this chunk

   ncol  = state%ncol

   ! open a file

   funit = 1234
   write(filename,'(2(a,i5.5),a)') 'nstep_',nstep,'_chunk_',lchnk,'_state.asc'

   open(unit=funit,file=trim(filename),access='sequential',action='write',form='formatted',iostat=ierr)
   if (ierr/=0) then
      write(iulog,*) 'Failed to open file '//trim(filename)//' for ASCII output of model state.' 
      call endrun
   end if

   ! write out components of state

   do ii = 1,ncol

      write(funit,*) 'lat ',state%lat(ii),' lon ',state%lon(ii)

      write(funit,*) 'U'
      do kk = 1,pver
         write(funit,*) kk, state%u(ii,kk)
      end do

      write(funit,*) 'V'
      do kk = 1,pver
         write(funit,*) kk, state%v(ii,kk)
      end do

      write(funit,*) 'S'
      do kk = 1,pver
         write(funit,*) kk, state%s(ii,kk)
      end do

      write(funit,*) 'PDEL'
      do kk = 1,pver
         write(funit,*) kk, state%pdel(ii,kk)
      end do

      ! tracers
     !do m=1,pcnst
      do m=1,3

         write(funit,*) trim(cnst_name(m))
         do kk = 1,pver
            write(funit,*) kk, state%q(ii,kk,m)
         end do

         write(iulog,*) 'nstep ',nstep,', lat ',state%lat(ii),', lon ',state%lon(ii),': ', &
                        'max ',trim(cnst_name(m)),' = ',maxval(state%q(ii,:,m))

      end do

   end do  ! column loop

   close(unit=funit)

end subroutine write_state_ascii

end module ascii_io
