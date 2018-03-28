! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains a collection of utilities related to timings in ELSI.
!!
module ELSI_TIMINGS

   use ELSI_CONSTANTS, only: STR_LEN,UNSET,UNSET_STR,N_TIMINGS
   use ELSI_DATATYPE,  only: elsi_handle,elsi_timings_handle
   use ELSI_IO,        only: elsi_say
   use ELSI_PRECISION, only: i4,r8

   implicit none

   private

   public :: elsi_get_time
   public :: elsi_init_timings
   public :: elsi_add_timing
   public :: elsi_print_timings
   public :: elsi_finalize_timings

contains

!>
!! This routine gets the current wallclock time.
!! (Taken from FHI-aims with permission of copyright holders)
!!
subroutine elsi_get_time(wtime)

   implicit none

   real(kind=r8), intent(out) :: wtime

   character(len=8)  :: cdate
   character(len=10) :: ctime
   integer(kind=i4)  :: val
   integer(kind=i4)  :: int_year
   real(kind=r8)     :: year
   real(kind=r8)     :: month
   real(kind=r8)     :: day
   real(kind=r8)     :: hour
   real(kind=r8)     :: minute
   real(kind=r8)     :: second
   real(kind=r8)     :: millisecond

   character(len=40), parameter :: caller = "elsi_get_time"

   call date_and_time(cdate,ctime)

   read(cdate(1:4),'(I4)') val

   int_year = val
   year     = real(val,kind=r8)-2009.0_r8 ! 2009 is an arbitrary zero

   day = year*365+floor(year/4.0_r8)

   read(cdate(5:6),'(I2)') val

   val = val-1

   do while(val > 0)
      if(val == 1) then
         day = day+31
      elseif(val == 2) then
         if(mod(int_year,4) == 0) then
            day = day+29
         else
            day = day+28
         endif
      elseif(val == 3) then
         day = day+31
      elseif(val == 4) then
         day = day+30
      elseif(val == 5) then
         day = day+31
      elseif(val == 6) then
         day = day+30
      elseif(val == 7) then
         day = day+31
      elseif(val == 8) then
         day = day+31
      elseif(val == 9) then
         day = day+30
      elseif(val == 10) then
         day = day+31
      elseif(val == 11) then
         day = day+30
      endif

      val = val-1
   enddo

   read(cdate(7:8),'(I2)') val
   day = day+real(val,kind=r8)-1

   read(ctime(1:2),'(I2)') val
   hour = real(val,kind=r8)

   read(ctime(3:4),'(I2)') val
   minute = real(val,kind=r8)

   read(ctime(5:6),'(I2)') val
   second = real(val,kind=r8)

   read(ctime(8:10),'(I3)') val
   millisecond = real(val,kind=r8)

   wtime = day*24.0_r8*3600.0_r8+hour*3600.0_r8+minute*60.0_r8+second+&
              millisecond*0.001_r8

end subroutine

!>
!! This routine initializes a timing handle.
!!
subroutine elsi_init_timings(t_h,set_label)

   implicit none

   type(elsi_timings_handle),  intent(inout)        :: t_h
   character(len=*),           intent(in), optional :: set_label

   character(len=40), parameter :: caller = "elsi_init_timings"

   t_h%size_timings = N_TIMINGS
   t_h%n_timings    = 0
   t_h%user_tag     = UNSET_STR

   if(present(set_label)) then
      t_h%set_label = set_label
   else
      t_h%set_label = UNSET_STR
   endif

   allocate(t_h%times(t_h%size_timings))
   allocate(t_h%elsi_tags(t_h%size_timings))
   allocate(t_h%user_tags(t_h%size_timings))

end subroutine

!>
!! This routine resizes arrays in the timing handle by a factor of two.
!!
subroutine elsi_resize_timing_arrays(t_h)

   implicit none

   type(elsi_timings_handle), intent(inout) :: t_h

   integer(kind=i4) :: i_timing

   real(kind=r8),          allocatable :: tmp_times(:)
   character(len=STR_LEN), allocatable :: tmp_elsi_tags(:)
   character(len=STR_LEN), allocatable :: tmp_user_tags(:)

   character(len=40), parameter :: caller = "elsi_resize_timing_arrays"

   allocate(tmp_times(t_h%size_timings))
   allocate(tmp_elsi_tags(t_h%size_timings))
   allocate(tmp_user_tags(t_h%size_timings))

   ! Save old timings
   do i_timing = 1,t_h%size_timings
      tmp_times(i_timing)     = t_h%times(i_timing)
      tmp_elsi_tags(i_timing) = t_h%elsi_tags(i_timing)
      tmp_user_tags(i_timing) = t_h%user_tags(i_timing)
   enddo

   ! Resize arrays by factor of two
   deallocate(t_h%times)
   deallocate(t_h%elsi_tags)
   deallocate(t_h%user_tags)

   t_h%size_timings = 2*t_h%size_timings

   allocate(t_h%times(t_h%size_timings))
   allocate(t_h%elsi_tags(t_h%size_timings))
   allocate(t_h%user_tags(t_h%size_timings))

   ! Fill in arrays with old timings
   do i_timing = 1,t_h%size_timings/2
      t_h%times(i_timing)     = tmp_times(i_timing)
      t_h%elsi_tags(i_timing) = tmp_elsi_tags(i_timing)
      t_h%user_tags(i_timing) = tmp_user_tags(i_timing)
   enddo

   deallocate(tmp_times)
   deallocate(tmp_elsi_tags)
   deallocate(tmp_user_tags)

end subroutine

!>
!! This routine adds provided timing to current handle.
!!
subroutine elsi_add_timing(t_h,time,elsi_tag,user_tag_in,iter_in)

   implicit none

   type(elsi_timings_handle), intent(inout)        :: t_h
   real(kind=r8),             intent(in)           :: time
   character(len=*),          intent(in)           :: elsi_tag
   character(len=*),          intent(in), optional :: user_tag_in
   integer(kind=i4),          intent(in), optional :: iter_in

   character(len=STR_LEN) :: user_tag
   integer(kind=i4)               :: iter

   character(len=40), parameter :: caller = "elsi_add_timings"

   if(present(iter_in)) then
      iter = iter_in
   else
      t_h%n_timings = t_h%n_timings+1
      iter          = t_h%n_timings
   endif

   ! If more timings than maximum array size, resize arrays
   do while(iter >= t_h%size_timings)
     call elsi_resize_timing_arrays(t_h)
   enddo

   if(present(user_tag_in)) then
      user_tag = user_tag_in
   else
      user_tag = t_h%user_tag
   endif

   t_h%times(iter)     = time
   t_h%elsi_tags(iter) = elsi_tag
   t_h%user_tags(iter) = user_tag

end subroutine

!>
!! This routine prints out timings collected so far.
!!
subroutine elsi_print_timings(e_h,t_h)

   implicit none

   type(elsi_handle),         intent(in) :: e_h
   type(elsi_timings_handle), intent(in) :: t_h

   integer(kind=i4)   :: iter
   real(kind=r8)      :: tmp
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_timings"

   call elsi_say(e_h,t_h%set_label)
   call elsi_say(e_h,"   #  wall_clock   [s]     elsi_tag             user_tag")

   do iter = 1,min(3,t_h%n_timings)
      write(info_str,"(I4,2X,F12.3,9X,A,1X,A)") iter,t_h%times(iter),&
         t_h%elsi_tags(iter),t_h%user_tags(iter)
      call elsi_say(e_h,info_str)
   enddo

   if(t_h%n_timings > 0) then
      tmp = maxval(t_h%times(1:t_h%n_timings))
      write(info_str,"(6X,F12.3,9X,A,18X,A)") tmp,"MAX",UNSET_STR
      call elsi_say(e_h,info_str)

      tmp = minval(t_h%times(1:t_h%n_timings))
      write(info_str,"(6X,F12.3,9X,A,18X,A)") tmp,"MIN",UNSET_STR
      call elsi_say(e_h,info_str)

      tmp = sum(t_h%times(1:t_h%n_timings))/t_h%n_timings
      write(info_str,"(6X,F12.3,9X,A,14X,A)") tmp,"AVERAGE",UNSET_STR
      call elsi_say(e_h,info_str)
   endif

end subroutine

!>
!! This routine deallocates the timings handle.
!!
subroutine elsi_finalize_timings(t_h)

   implicit none

   type(elsi_timings_handle), intent(inout) :: t_h

   character(len=40), parameter :: caller = "elsi_finalize_timings"

   if(allocated(t_h%times)) then
      deallocate(t_h%times)
   endif
   if(allocated(t_h%elsi_tags)) then
      deallocate(t_h%elsi_tags)
   endif
   if(allocated(t_h%user_tags)) then
      deallocate(t_h%user_tags)
   endif

   t_h%n_timings    = 0
   t_h%size_timings = UNSET
   t_h%user_tag     = UNSET_STR
   t_h%set_label    = UNSET_STR

end subroutine

end module ELSI_TIMINGS
