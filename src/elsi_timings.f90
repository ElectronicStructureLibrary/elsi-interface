! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains a collection of utilities related to timings in ELSI.
!!
module ELSI_TIMINGS

   use ELSI_CONSTANTS, only: STR_LEN,UNSET,UNSET_STR,N_TIMINGS,TIME_LEN,HUMAN,&
                             NO_COMMA,COMMA_BEFORE,COMMA_AFTER,OUTPUT_EV,&
                             REAL_DATA
   use ELSI_DATATYPE,  only: elsi_handle,elsi_io_handle,elsi_timings_handle
   use ELSI_IO,        only: elsi_say,elsi_say_setting,elsi_print_versioning,&
                             elsi_print_solver_settings,&
                             elsi_print_matrix_settings,&
                             elsi_print_handle_summary,elsi_append_string,&
                             elsi_truncate_string,elsi_start_json_record,&
                             elsi_finish_json_record
   use ELSI_PRECISION, only: i4,r8
   use ELSI_UTILS,     only: elsi_get_solver_tag,elsi_get_datetime_rfc3339

   implicit none

   private

   public :: elsi_get_time
   public :: elsi_init_timings
   public :: elsi_print_timings
   public :: elsi_finalize_timings
   public :: elsi_process_timing

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

!>
!! This routine acts on the timing information, both to add it to the solver
!! timing summary and print it to the timing file.
!!
subroutine elsi_process_timing(e_h,output_type,data_type,solver_used,dt0,t0)

   implicit none

   type(elsi_handle),       intent(inout) :: e_h
   integer(kind=i4),        intent(in)    :: output_type
   integer(kind=i4),        intent(in)    :: data_type
   integer(kind=i4),        intent(in)    :: solver_used
   character(len=TIME_LEN), intent(in)    :: dt0
   real(kind=r8),           intent(in)    :: t0

   character(len=STR_LEN) :: solver_tag
   character(len=STR_LEN) :: elsi_tag
   character(len=STR_LEN) :: user_tag
   real(kind=r8)          :: t1
   real(kind=r8)          :: t_total
   integer(kind=i4)       :: tmp_int
   integer(kind=i4)       :: comma_json_save
   integer(kind=i4)       :: iter
   type(elsi_io_handle)   :: io_h

   character(len=40), parameter :: caller = "elsi_process_timing"

   io_h = e_h%timings_file

   call elsi_get_time(t1)

   tmp_int    = e_h%solver
   e_h%solver = solver_used
   t_total    = t1-t0

   ! Output information about this solver invocation
   call elsi_get_solver_tag(e_h,data_type,solver_tag)
   call elsi_add_timing(e_h%timings,t_total,solver_tag)

   if(e_h%output_timings) then
      iter = e_h%timings%n_timings

      ! Avoid comma at the end of the last entry
      comma_json_save = io_h%comma_json

      if(e_h%timings%n_timings == 1) then
         io_h%comma_json = NO_COMMA
      else
         io_h%comma_json = COMMA_BEFORE
      endif

      elsi_tag = adjustr(trim(solver_tag))
      user_tag = adjustr(trim(e_h%timings%user_tags(iter)))

      call elsi_print_timing(e_h,io_h,output_type,data_type,dt0,t_total,iter,&
              elsi_tag,user_tag)

      io_h%comma_json = comma_json_save
   endif

   e_h%solver = tmp_int

end subroutine

!>
!! This routine prints timing and relevant information.
!!
subroutine elsi_print_timing(e_h,io_h,output_type,data_type,dt0,t_total,iter,&
              elsi_tag,user_tag)

   implicit none

   type(elsi_handle),       intent(inout) :: e_h
   type(elsi_io_handle),    intent(inout) :: io_h
   integer(kind=i4),        intent(in)    :: output_type
   integer(kind=i4),        intent(in)    :: data_type
   character(len=TIME_LEN), intent(in)    :: dt0
   real(kind=r8),           intent(in)    :: t_total
   integer(kind=i4),        intent(in)    :: iter
   character(len=*),        intent(in)    :: elsi_tag
   character(len=*),        intent(in)    :: user_tag

   integer(kind=i4)        :: comma_json_save
   character(len=TIME_LEN) :: dt_record
   character(len=200)      :: info_str

   character(len=40), parameter :: caller = "elsi_print_timing"

   call elsi_get_datetime_rfc3339(dt_record)

   ! Print out patterned header, versioning information, and timing details
   if(io_h%file_format == HUMAN) then
      call elsi_say(e_h,"-------------------------------------------------------------------------",io_h)
      write(info_str,"(A,I10)") "Start of ELSI Solver Iteration ",iter
      call elsi_say(e_h,info_str,io_h)
      call elsi_say(e_h,"",io_h)
      call elsi_print_versioning(e_h,io_h)
      call elsi_say(e_h,"",io_h)
      call elsi_say(e_h,"Timing Details",io_h)
      call elsi_append_string(io_h%prefix,"  ")
      if(output_type == OUTPUT_EV) then
         call elsi_say_setting(e_h,"Output Type","EIGENSOLUTION",io_h)
      else
         call elsi_say_setting(e_h,"Output Type","DENSITY MATRIX",io_h)
      endif
      if(data_type == REAL_DATA) then
         call elsi_say_setting(e_h,"Data Type","REAL",io_h)
      else
         call elsi_say_setting(e_h,"Data Type","COMPLEX",io_h)
      endif
      call elsi_say_setting(e_h,"ELSI Tag",elsi_tag,io_h)
      call elsi_say_setting(e_h,"User Tag",user_tag,io_h)
      call elsi_say_setting(e_h,"Timing (s)",t_total,io_h)
      call elsi_truncate_string(io_h%prefix,2)
   else
      call elsi_start_json_record(e_h,io_h%comma_json==COMMA_BEFORE,io_h)
      io_h%comma_json = COMMA_AFTER ! Add commas behind all records before final
      call elsi_print_versioning(e_h,io_h)
      call elsi_say_setting(e_h,"iteration",iter,io_h)
      if(output_type == OUTPUT_EV) then
         call elsi_say_setting(e_h,"output_type","EIGENSOLUTION",io_h)
      else
         call elsi_say_setting(e_h,"output_type","DENSITY MATRIX",io_h)
      endif
      if(data_type == REAL_DATA) then
         call elsi_say_setting(e_h,"data_type","REAL",io_h)
      else
         call elsi_say_setting(e_h,"data_type","COMPLEX",io_h)
      endif
      call elsi_say_setting(e_h,"elsi_tag",elsi_tag,io_h)
      call elsi_say_setting(e_h,"user_tag",user_tag,io_h)
      call elsi_say_setting(e_h,"start_datetime",dt0,io_h)
      call elsi_say_setting(e_h,"record_datetime",dt_record,io_h)
      call elsi_say_setting(e_h,"total_time",t_total,io_h)
   endif

   ! Print out handle summary
   if(io_h%file_format == HUMAN) then
      call elsi_say(e_h,"",io_h)
   endif
   call elsi_print_handle_summary(e_h,io_h)

   ! Print out matrix storage format settings
   if(io_h%file_format == HUMAN) then
      call elsi_say(e_h,"",io_h)
   endif
   call elsi_print_matrix_settings(e_h,io_h)

   ! Print out solver settings
   if(io_h%file_format == HUMAN) then
      call elsi_say(e_h,"",io_h)
   endif
   io_h%comma_json = NO_COMMA ! Final record in this scope
   call elsi_print_solver_settings(e_h,io_h)

   ! Print out patterned footer
   io_h%comma_json = comma_json_save
   if(io_h%file_format == HUMAN) then
      call elsi_say(e_h,"",io_h)
      write(info_str,"(A,I10)") "End of ELSI Solver Iteration   ",iter
      call elsi_say(e_h,info_str,io_h)
      call elsi_say(e_h,"-------------------------------------------------------------------------",io_h)
      call elsi_say(e_h,"",io_h)
   else
      call elsi_finish_json_record(e_h,io_h%comma_json==COMMA_AFTER,io_h)
   endif

end subroutine

end module ELSI_TIMINGS
