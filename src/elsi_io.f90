! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module performs IO to stdout and files.
!!
module ELSI_IO

   use ELSI_CONSTANTS, only: UNSET,UNSET_STR,HUMAN,JSON,MULTI_PROC,SINGLE_PROC,&
                             ELPA_SOLVER,SIPS_SOLVER,OMM_SOLVER,PEXSI_SOLVER,&
                             DMP_SOLVER,BLACS_DENSE,PEXSI_CSC,SIESTA_CSC,&
                             COMMA_BEFORE,COMMA_AFTER,NO_COMMA,LOG_NAME,&
                             TIME_LEN,STR_LEN,UUID_LEN,OUTPUT_EV,REAL_DATA
   use ELSI_DATATYPE,  only: elsi_handle,elsi_io_handle
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi,mpi_character
   use ELSI_PRECISION, only: r8,i4,i8

   implicit none

   private

   public :: elsi_say
   public :: elsi_say_setting
   public :: elsi_init_io
   public :: elsi_reset_io_handle
   public :: elsi_io_add_entry
   public :: elsi_print_handle_summary
   public :: elsi_print_versioning
   public :: elsi_print_solver_settings
   public :: elsi_print_matrix_settings
   public :: elsi_append_string
   public :: elsi_truncate_string
   public :: elsi_open_json_file
   public :: elsi_close_json_file
   public :: elsi_start_json_record
   public :: elsi_finish_json_record
   public :: elsi_get_time
   public :: elsi_get_datetime_rfc3339

   interface elsi_say_setting
      module procedure elsi_say_setting_i4,&
                       elsi_say_setting_r8,&
                       elsi_say_setting_log,&
                       elsi_say_setting_str
   end interface

contains

!>
!! This routine prints a message.
!!
subroutine elsi_say(e_h,info_str,io_h)

   implicit none

   type(elsi_handle),    intent(in) :: e_h
   character(len=*),     intent(in) :: info_str
   type(elsi_io_handle), intent(in) :: io_h

   character(len=40), parameter :: caller = "elsi_say"

   if(io_h%print_info .and. e_h%myid_all == 0) then
      if(allocated(io_h%prefix)) then
         write(io_h%print_unit,"(2A)") io_h%prefix,trim(info_str)
      else
         write(io_h%print_unit,"(A)") trim(info_str)
      endif
   endif

end subroutine

!>
!! This routine initializes a handle for reading and writing to files.
!!
subroutine elsi_init_io(io_h,print_unit,file_name,file_format,print_info,&
              prefix,comma_json)

   implicit none

   type(elsi_io_handle), intent(out)          :: io_h
   integer(kind=i4),     intent(in)           :: print_unit
   character(len=*),     intent(in), optional :: file_name
   integer(kind=i4),     intent(in), optional :: file_format
   logical,              intent(in), optional :: print_info
   character(len=*),     intent(in), optional :: prefix
   integer(kind=i4),     intent(in), optional :: comma_json

   character(len=40), parameter :: caller = "elsi_init_io"

   ! For safety
   call elsi_reset_io_handle(io_h)

   io_h%handle_init = .true.
   io_h%print_unit  = print_unit

   if(present(file_name)) then
      io_h%file_name = file_name
   else
      io_h%file_name = LOG_NAME
   endif

   if(present(file_format)) then
      io_h%file_format = file_format
   else
      io_h%file_format = HUMAN
   endif

   if(present(print_info)) then
      io_h%print_info = print_info
   else
      io_h%print_info = .true.
   endif

   if(present(prefix)) then
      io_h%prefix = prefix
   else
      if(allocated(io_h%prefix)) then
         deallocate(io_h%prefix)
      endif
   endif

   if(present(comma_json)) then
      io_h%comma_json = comma_json
   else
      io_h%comma_json = UNSET
   endif

end subroutine

!>
!! This routine resets a handle.
!!
subroutine elsi_reset_io_handle(io_h)

   implicit none

   type(elsi_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_reset_io_handle"

   io_h%handle_init = .false.
   io_h%file_name   = UNSET_STR
   io_h%file_format = UNSET
   io_h%print_info  = .false.
   io_h%print_unit  = UNSET
   io_h%comma_json  = UNSET
   io_h%n_records   = 0
   io_h%user_tag    = UNSET_STR
   io_h%elsi_tag    = UNSET_STR

   if(allocated(io_h%prefix)) then
      deallocate(io_h%prefix)
   endif

end subroutine

!>
!! This routine adds an entry to the log file.
!!
subroutine elsi_io_add_entry(e_h,output_type,data_type,solver_used,dt0,t0)

   implicit none

   type(elsi_handle),       intent(inout) :: e_h
   integer(kind=i4),        intent(in)    :: output_type
   integer(kind=i4),        intent(in)    :: data_type
   integer(kind=i4),        intent(in)    :: solver_used
   character(len=TIME_LEN), intent(in)    :: dt0
   real(kind=r8),           intent(in)    :: t0

   real(kind=r8)           :: t1
   real(kind=r8)           :: t_total
   integer(kind=i4)        :: tmp_int
   integer(kind=i4)        :: comma_save
   character(len=STR_LEN)  :: solver_tag
   character(len=STR_LEN)  :: elsi_tag
   character(len=STR_LEN)  :: user_tag
   character(len=TIME_LEN) :: dt_record

   character(len=40), parameter :: caller = "elsi_io_add_entry"

   if(e_h%log_file%print_info) then
      call elsi_get_time(t1)
      call elsi_get_solver_tag(e_h,data_type,solver_tag)

      tmp_int                = e_h%solver
      e_h%solver             = solver_used
      t_total                = t1-t0
      e_h%log_file%n_records = e_h%log_file%n_records+1
      comma_save             = e_h%log_file%comma_json

      if(e_h%log_file%n_records == 1) then
         e_h%log_file%comma_json = NO_COMMA

         if(e_h%myid_all == 0) then
            call elsi_open_json_file(e_h,e_h%log_file%print_unit,&
                    e_h%log_file%file_name,.true.,e_h%log_file)
         else
            e_h%log_file%print_unit = UNSET
            e_h%log_file%file_name  = UNSET_STR
         endif

         if(.not. e_h%uuid_exists) then
            call elsi_gen_uuid(e_h)
            e_h%uuid_exists = .true.
         endif
      else
         e_h%log_file%comma_json = COMMA_BEFORE
      endif

      elsi_tag = adjustr(trim(solver_tag))
      user_tag = adjustr(trim(e_h%log_file%user_tag))

      call elsi_get_datetime_rfc3339(dt_record)

      call elsi_start_json_record(e_h,e_h%log_file%comma_json==COMMA_BEFORE,&
              e_h%log_file)
      e_h%log_file%comma_json = COMMA_AFTER
      call elsi_print_versioning(e_h,e_h%log_file)
      call elsi_say_setting(e_h,"iteration",e_h%log_file%n_records,e_h%log_file)
      if(output_type == OUTPUT_EV) then
         call elsi_say_setting(e_h,"output_type","EIGENSOLUTION",e_h%log_file)
      else
         call elsi_say_setting(e_h,"output_type","DENSITY MATRIX",e_h%log_file)
      endif
      if(data_type == REAL_DATA) then
         call elsi_say_setting(e_h,"data_type","REAL",e_h%log_file)
      else
         call elsi_say_setting(e_h,"data_type","COMPLEX",e_h%log_file)
      endif
      call elsi_say_setting(e_h,"elsi_tag",elsi_tag,e_h%log_file)
      call elsi_say_setting(e_h,"user_tag",user_tag,e_h%log_file)
      call elsi_say_setting(e_h,"start_datetime",dt0,e_h%log_file)
      call elsi_say_setting(e_h,"record_datetime",dt_record,e_h%log_file)
      call elsi_say_setting(e_h,"total_time",t_total,e_h%log_file)

      call elsi_print_handle_summary(e_h,e_h%log_file)
      call elsi_print_matrix_settings(e_h,e_h%log_file)

      e_h%log_file%comma_json = NO_COMMA

      call elsi_print_solver_settings(e_h,e_h%log_file)

      e_h%log_file%comma_json = comma_save

      call elsi_finish_json_record(e_h,e_h%log_file%comma_json==COMMA_AFTER,&
              e_h%log_file)

      e_h%log_file%comma_json = comma_save
      e_h%solver = tmp_int
   endif

end subroutine

!>
!! This routine prints the state of the handle.
!!
subroutine elsi_print_handle_summary(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   real(kind=r8)      :: sparsity
   integer(kind=i4)   :: comma_save
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_handle_summary"

   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Physical Properties"
      call elsi_say(e_h,info_str,io_h)

      call elsi_append_string(io_h%prefix,"  ")
      call elsi_say_setting(e_h,"Number of electrons",e_h%n_electrons,io_h)
      if(e_h%parallel_mode == MULTI_PROC) then
         call elsi_say_setting(e_h,"Number of spins",e_h%n_spins,io_h)
         call elsi_say_setting(e_h,"Number of k-points",e_h%n_kpts,io_h)
      endif
      if(e_h%solver == ELPA_SOLVER .or. e_h%solver == SIPS_SOLVER) then
         call elsi_say_setting(e_h,"Number of states",e_h%n_states,io_h)
      endif
      call elsi_truncate_string(io_h%prefix,2)

      write(info_str,"(A)") ""
      call elsi_say(e_h,info_str,io_h)
      write(info_str,"(A)") "Matrix Properties"
      call elsi_say(e_h,info_str,io_h)

      call elsi_append_string(io_h%prefix,"  ")
      if(e_h%matrix_format == BLACS_DENSE) then
         call elsi_say_setting(e_h,"Matrix format","BLACS_DENSE",io_h)
      elseif(e_h%matrix_format == PEXSI_CSC) then
         call elsi_say_setting(e_h,"Matrix format","PEXSI_CSC",io_h)
      elseif(e_h%matrix_format == SIESTA_CSC) then
         call elsi_say_setting(e_h,"Matrix format","SIESTA_CSC",io_h)
      endif
      call elsi_say_setting(e_h,"Number of basis functions",e_h%n_basis,io_h)
      if(e_h%parallel_mode == MULTI_PROC) then
         sparsity = 1.0_r8-(1.0_r8*e_h%nnz_g/e_h%n_basis/e_h%n_basis)
         call elsi_say_setting(e_h,"Matrix sparsity",sparsity,io_h)
      endif
      call elsi_truncate_string(io_h%prefix,2)

      write(info_str,"(A)") ""
      call elsi_say(e_h,info_str,io_h)
      write(info_str,"(A)") "Computational Details"
      call elsi_say(e_h,info_str,io_h)

      call elsi_append_string(io_h%prefix,"  ")
      if(e_h%parallel_mode == MULTI_PROC) then
         call elsi_say_setting(e_h,"Parallel mode","MULTI_PROC",io_h)
      elseif(e_h%parallel_mode == SINGLE_PROC) then
         call elsi_say_setting(e_h,"Parallel mode","SINGLE_PROC",io_h)
      endif
      call elsi_say_setting(e_h,"Number of MPI tasks",e_h%n_procs,io_h)
      if(e_h%solver == ELPA_SOLVER) then
         call elsi_say_setting(e_h,"Solver requested","ELPA",io_h)
      elseif(e_h%solver == OMM_SOLVER) then
         call elsi_say_setting(e_h,"Solver requested","libOMM",io_h)
      elseif(e_h%solver == PEXSI_SOLVER) then
         call elsi_say_setting(e_h,"Solver requested","PEXSI",io_h)
      elseif(e_h%solver == SIPS_SOLVER) then
         call elsi_say_setting(e_h,"Solver requested","SIPS",io_h)
      elseif(e_h%solver == DMP_SOLVER) then
         call elsi_say_setting(e_h,"Solver requested","DMP",io_h)
      else
         call elsi_stop(e_h,"Unsupported solver.",caller)
      endif
      call elsi_say_setting(e_h,"Number of ELSI calls",e_h%n_elsi_calls,io_h)
      call elsi_truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      comma_save      = io_h%comma_json
      io_h%comma_json = COMMA_AFTER

      call elsi_say_setting(e_h,"n_electrons",e_h%n_electrons,io_h)
      if(e_h%parallel_mode == MULTI_PROC) then
         call elsi_say_setting(e_h,"n_spin",e_h%n_spins,io_h)
         call elsi_say_setting(e_h,"n_kpts",e_h%n_kpts,io_h)
      endif
      if(e_h%solver == ELPA_SOLVER .or. e_h%solver == SIPS_SOLVER) then
         call elsi_say_setting(e_h,"n_states",e_h%n_states,io_h)
      endif

      if(e_h%matrix_format == BLACS_DENSE) then
         call elsi_say_setting(e_h,"matrix_format","BLACS_DENSE",io_h)
      elseif(e_h%matrix_format == PEXSI_CSC) then
         call elsi_say_setting(e_h,"matrix_format","PEXSI_CSC",io_h)
      elseif(e_h%matrix_format == SIESTA_CSC) then
         call elsi_say_setting(e_h,"matrix_format","SIESTA_CSC",io_h)
      endif
      call elsi_say_setting(e_h,"n_basis",e_h%n_basis,io_h)
      if(e_h%parallel_mode == MULTI_PROC) then
         sparsity = 1.0_r8-(1.0_r8*e_h%nnz_g/e_h%n_basis/e_h%n_basis)
         call elsi_say_setting(e_h,"sparsity",sparsity,io_h)
      endif

      if(e_h%parallel_mode == MULTI_PROC) then
         call elsi_say_setting(e_h,"parallel_mode","MULTI_PROC",io_h)
      elseif(e_h%parallel_mode == SINGLE_PROC) then
         call elsi_say_setting(e_h,"parallel_mode","SINGLE_PROC",io_h)
      endif
      io_h%comma_json = comma_save
      call elsi_say_setting(e_h,"n_procs",e_h%n_procs,io_h)
      if(e_h%solver == ELPA_SOLVER) then
         call elsi_say_setting(e_h,"solver","ELPA",io_h)
      elseif(e_h%solver == OMM_SOLVER) then
         call elsi_say_setting(e_h,"solver","libOMM",io_h)
      elseif(e_h%solver == PEXSI_SOLVER) then
         call elsi_say_setting(e_h,"solver","PEXSI",io_h)
      elseif(e_h%solver == SIPS_SOLVER) then
         call elsi_say_setting(e_h,"solver","SIPS",io_h)
      elseif(e_h%solver == DMP_SOLVER) then
         call elsi_say_setting(e_h,"solver","DMP",io_h)
      else
         call elsi_stop(e_h,"Unsupported solver.",caller)
      endif
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

end subroutine

!>
!! This routine prints versioning information.
!!
subroutine elsi_print_versioning(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   integer(kind=i4)   :: comma_save
   logical            :: COMMIT_MODIFIED
   character(len=10)  :: DATE_STAMP
   character(len=40)  :: COMMIT
   character(len=8)   :: COMMIT_ABBREV
   character(len=40)  :: COMMIT_MSG_ABBREV
   character(len=40)  :: HOSTNAME
   character(len=20)  :: DATETIME
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_versioning"

   call elsi_version_info(DATE_STAMP,COMMIT,COMMIT_ABBREV,COMMIT_MODIFIED,&
           COMMIT_MSG_ABBREV,HOSTNAME,DATETIME)

   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "ELSI Versioning Information"
      call elsi_say(e_h,info_str,io_h)
      call elsi_append_string(io_h%prefix,"  ")

      call elsi_say_setting(e_h,"Date stamp",trim(DATE_STAMP),io_h)
      call elsi_say_setting(e_h,"Git commit (abbrev.)",trim(COMMIT_ABBREV),io_h)
      call elsi_say_setting(e_h,"Git commit modified?",COMMIT_MODIFIED,io_h)

      call elsi_truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      comma_save      = io_h%comma_json
      io_h%comma_json = COMMA_AFTER

      call elsi_say_setting(e_h,"data_source","ELSI",io_h)
      call elsi_say_setting(e_h,"date_stamp",DATE_STAMP,io_h)
      call elsi_say_setting(e_h,"git_commit",COMMIT,io_h)
      call elsi_say_setting(e_h,"git_commit_modified",COMMIT_MODIFIED,io_h)
      call elsi_say_setting(e_h,"git_message_abbrev",COMMIT_MSG_ABBREV,io_h)
      call elsi_say_setting(e_h,"source_created_on_hostname",HOSTNAME,io_h)
      call elsi_say_setting(e_h,"source_created_at_datetime",DATETIME,io_h)
      call elsi_say_setting(e_h,"calling_code",e_h%caller,io_h)

      io_h%comma_json = comma_save
      call elsi_say_setting(e_h,"uuid",e_h%uuid,io_h)
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

end subroutine

!>
!! This routine prints out settings for the current (user-indicated) solver.
!!
subroutine elsi_print_solver_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_print_solver_settings"

   select case(e_h%solver)
   case(DMP_SOLVER)
      call elsi_print_dmp_settings(e_h,io_h)
   case(ELPA_SOLVER)
      call elsi_print_elpa_settings(e_h,io_h)
   case(OMM_SOLVER)
      call elsi_print_omm_settings(e_h,io_h)
   case(PEXSI_SOLVER)
      call elsi_print_pexsi_settings(e_h,io_h)
   case(SIPS_SOLVER)
      call elsi_print_sips_settings(e_h,io_h)
   case default
      call elsi_stop(e_h,"Unsupported solver.",caller)
   end select

end subroutine

!>
!! This routine prints out settings for DMP.
!!
subroutine elsi_print_dmp_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   integer(kind=i4)   :: comma_save
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_dmp_settings"

   comma_save      = io_h%comma_json
   io_h%comma_json = COMMA_AFTER

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Solver Settings (DMP)"
   else
      write(info_str,"(A)") '"solver_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call elsi_append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"dmp_n_states",e_h%dmp_n_states,io_h)
   call elsi_say_setting(e_h,"dmp_method",e_h%dmp_method,io_h)
   call elsi_say_setting(e_h,"dmp_max_power",e_h%dmp_max_power,io_h)
   call elsi_say_setting(e_h,"dmp_max_iter",e_h%dmp_max_iter,io_h)
   io_h%comma_json = NO_COMMA
   call elsi_say_setting(e_h,"dmp_tol",e_h%dmp_tol,io_h)
   call elsi_truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_save

   if(io_h%file_format == JSON) then
      if(io_h%comma_json == COMMA_AFTER) then
         write(info_str,"(A)") '},'
      else
         write(info_str,"(A)") '}'
      endif

      call elsi_say(e_h,info_str,io_h)
   endif

end subroutine

!>
!! This routine prints out settings for ELPA.
!!
subroutine elsi_print_elpa_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   integer(kind=i4)   :: comma_save
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_elpa_settings"

   comma_save      = io_h%comma_json
   io_h%comma_json = COMMA_AFTER

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Solver Settings (ELPA)"
   else
      write(info_str,"(A)") '"solver_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call elsi_append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"elpa_solver",e_h%elpa_solver,io_h)
   io_h%comma_json = NO_COMMA
   call elsi_say_setting(e_h,"elpa_n_states",e_h%n_states,io_h)
   call elsi_truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_save

   if(io_h%file_format == JSON) then
      if(io_h%comma_json == COMMA_AFTER) then
         write(info_str,"(A)") '},'
      else
         write(info_str,"(A)") '}'
      endif

      call elsi_say(e_h,info_str,io_h)
   endif

end subroutine

!>
!! This routine prints out settings for libOMM.
!!
subroutine elsi_print_omm_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   integer(kind=i4)   :: comma_save
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_omm_settings"

   comma_save      = io_h%comma_json
   io_h%comma_json = COMMA_AFTER

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Solver Settings (libOMM)"
   else
      write(info_str,"(A)") '"solver_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

  ! Settings
   call elsi_append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"omm_n_states",e_h%omm_n_states,io_h)
   call elsi_say_setting(e_h,"omm_n_elpa",e_h%omm_n_elpa,io_h)
   call elsi_say_setting(e_h,"omm_flavor",e_h%omm_flavor,io_h)
   io_h%comma_json = NO_COMMA
   call elsi_say_setting(e_h,"omm_tol",e_h%omm_tol,io_h)
   call elsi_truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_save

   if(io_h%file_format == JSON) then
      if(io_h%comma_json == COMMA_AFTER) then
         write(info_str,"(A)") '},'
      else
         write(info_str,"(A)") '}'
      endif

      call elsi_say(e_h,info_str,io_h)
   endif

end subroutine

!>
!! This routine prints out settings for PEXSI.
!!
subroutine elsi_print_pexsi_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   integer(kind=i4)   :: comma_save
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_pexsi_settings"

   comma_save      = io_h%comma_json
   io_h%comma_json = COMMA_AFTER

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Solver Settings (PEXSI)"
   else
      write(info_str,"(A)") '"solver_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call elsi_append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"pexsi_np_per_pole",e_h%pexsi_np_per_pole,io_h)
   call elsi_say_setting(e_h,"pexsi_np_per_point",e_h%pexsi_np_per_point,io_h)
   call elsi_say_setting(e_h,"pexsi_n_prow_pexsi",e_h%pexsi_n_prow,io_h)
   call elsi_say_setting(e_h,"pexsi_n_pcol_pexsi",e_h%pexsi_n_pcol,io_h)
   call elsi_say_setting(e_h,"pexsi_delta_e",e_h%pexsi_options%deltaE,io_h)
   call elsi_say_setting(e_h,"pexsi_gap",e_h%pexsi_options%gap,io_h)
   call elsi_say_setting(e_h,"pexsi_n_pole",e_h%pexsi_options%numPole,io_h)
   call elsi_say_setting(e_h,"pexsi_n_point",e_h%pexsi_options%nPoints,io_h)
   io_h%comma_json = NO_COMMA
   call elsi_say_setting(e_h,"pexsi_np_symbfact",e_h%pexsi_options%npSymbFact,&
           io_h)
   call elsi_truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_save

   if(io_h%file_format == JSON) then
      if(io_h%comma_json == COMMA_AFTER) then
         write(info_str,"(A)") '},'
      else
         write(info_str,"(A)") '}'
      endif

      call elsi_say(e_h,info_str,io_h)
   endif

end subroutine

!>
!! This routine prints out settings for SIPS.
!!
subroutine elsi_print_sips_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   integer(kind=i4)   :: comma_save
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_sips_settings"

   comma_save      = io_h%comma_json
   io_h%comma_json = COMMA_AFTER

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Solver Settings (SIPS)"
   else
      write(info_str,"(A)") '"solver_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call elsi_append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"sips_n_states",e_h%n_states,io_h)
   call elsi_say_setting(e_h,"sips_n_elpa",e_h%sips_n_elpa,io_h)
   call elsi_say_setting(e_h,"sips_n_slices",e_h%sips_n_slices,io_h)
   call elsi_say_setting(e_h,"sips_slice_type",e_h%sips_slice_type,io_h)
   call elsi_say_setting(e_h,"sips_np_per_slice",e_h%sips_np_per_slice,io_h)
   call elsi_say_setting(e_h,"sips_buffer",e_h%sips_buffer,io_h)
   io_h%comma_json = NO_COMMA
   call elsi_say_setting(e_h,"sips_first_ev",e_h%sips_first_ev,io_h)
   call elsi_truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_save

   if(io_h%file_format == JSON) then
      if(io_h%comma_json == COMMA_AFTER) then
         write(info_str,"(A)") '},'
      else
         write(info_str,"(A)") '}'
      endif

      call elsi_say(e_h,info_str,io_h)
   endif

end subroutine

!>
!! This routine prints out settings for the matirx format.
!!
subroutine elsi_print_matrix_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_print_matrix_settings"

   if(e_h%matrix_format == BLACS_DENSE) then
      call elsi_print_den_settings(e_h,io_h)
   else
      call elsi_print_csc_settings(e_h,io_h)
   endif

end subroutine

!>
!! This routine prints out settings for the dense matrix format.
!!
subroutine elsi_print_den_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   integer(kind=i4)   :: comma_save
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_den_settings"

   comma_save      = io_h%comma_json
   io_h%comma_json = COMMA_AFTER

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Dense Matrix Format Settings"
   else
      write(info_str,"(A)") '"matrix_format_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call elsi_append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"blk_row",e_h%blk_row,io_h)
   call elsi_say_setting(e_h,"blk_col",e_h%blk_col,io_h)
   call elsi_say_setting(e_h,"n_prow",e_h%n_prow,io_h)
   call elsi_say_setting(e_h,"n_pcol",e_h%n_pcol,io_h)
   io_h%comma_json = NO_COMMA
   call elsi_say_setting(e_h,"blacs_ready",e_h%blacs_ready,io_h)
   call elsi_truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_save

   if(io_h%file_format == JSON) then
      if(io_h%comma_json == COMMA_AFTER) then
         write(info_str,"(2A)") '},'
      else
         write(info_str,"(2A)") '}'
      endif

      call elsi_say(e_h,info_str,io_h)
   endif

end subroutine

!>
!! This routine prints out settings for the sparse matrix format.
!!
subroutine elsi_print_csc_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   integer(kind=i4)   :: comma_save
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_csc_settings"

   comma_save      = io_h%comma_json
   io_h%comma_json = COMMA_AFTER

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Sparse Matrix Format Settings"
   else
      write(info_str,"(A)") '"matrix_format_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call elsi_append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"zero_def",e_h%zero_def,io_h)
   call elsi_say_setting(e_h,"blk_sp2",e_h%blk_sp2,io_h)
   call elsi_say_setting(e_h,"pexsi_csc_ready",e_h%pexsi_csc_ready,io_h)
   io_h%comma_json = NO_COMMA
   call elsi_say_setting(e_h,"siesta_csc_ready",e_h%siesta_csc_ready,io_h)
   call elsi_truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_save

   if(io_h%file_format == JSON) then
      if(io_h%comma_json == COMMA_AFTER) then
         write(info_str,"(2A)") '},'
      else
         write(info_str,"(2A)") '}'
      endif

      call elsi_say(e_h,info_str,io_h)
   endif

end subroutine

!>
!! This routine prints out an integer-type setting.
!!
subroutine elsi_say_setting_i4(e_h,label,setting,io_h)

   implicit none

   type(elsi_handle),    intent(in) :: e_h
   character(len=*),     intent(in) :: label
   integer(kind=i4),     intent(in) :: setting
   type(elsi_io_handle), intent(in) :: io_h

   character(len=28) :: label_ljust
   character(len=20) :: int_string

   character(len=40), parameter :: caller = "elsi_say_setting_i4"

   write(int_string,'(I20)') setting

   label_ljust = label

   if(io_h%print_info .and. e_h%myid_all == 0) then
      if(io_h%file_format == HUMAN) then
         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,I25)") io_h%prefix,label_ljust,&
               " : ",setting
         else
            write(io_h%print_unit,"(A28,A3,I25)") label_ljust," : ",setting
         endif
      elseif(io_h%file_format == JSON) then
         if(io_h%comma_json == COMMA_AFTER) then
            if(allocated(io_h%prefix)) then
               write(io_h%print_unit,"(A)") io_h%prefix // '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(int_string)) // ","
            else
               write(io_h%print_unit,"(A)") '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(int_string)) // ","
            endif
         else
            if(allocated(io_h%prefix)) then
               write(io_h%print_unit,"(A)") io_h%prefix // '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(int_string))
            else
               write(io_h%print_unit,"(A)") '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(int_string))
            endif
         endif
      else
         call elsi_stop(e_h,"Unsupported output format.",caller)
      endif
   endif

end subroutine

!>
!! This routine prints out a real-type setting.
!!
subroutine elsi_say_setting_r8(e_h,label,setting,io_h)

   implicit none

   type(elsi_handle),    intent(in) :: e_h
   character(len=*),     intent(in) :: label
   real(kind=r8),        intent(in) :: setting
   type(elsi_io_handle), intent(in) :: io_h

   character(len=28) :: label_ljust
   character(len=20) :: real_string

   character(len=40), parameter :: caller = "elsi_say_setting_r8"

   write(real_string,'(E20.8)') setting

   label_ljust = label

   if(io_h%print_info .and. e_h%myid_all == 0) then
      if(io_h%file_format == HUMAN) then
         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,E25.8)") io_h%prefix,label_ljust,&
               " : ",setting
         else
            write(io_h%print_unit,"(A28,A3,E25.8)") label_ljust," : ",setting
         endif
      elseif(io_h%file_format == JSON) then
         if(io_h%comma_json == COMMA_AFTER) then
            if(allocated(io_h%prefix)) then
               write(io_h%print_unit,"(A)") io_h%prefix // '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(real_string)) // ","
            else
               write(io_h%print_unit,"(A)") '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(real_string)) // ","
            endif
         else
            if(allocated(io_h%prefix)) then
               write(io_h%print_unit,"(A)") io_h%prefix // '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(real_string))
            else
               write(io_h%print_unit,"(A)") '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(real_string))
            endif
         endif
      else
         call elsi_stop(e_h,"Unsupported output format.",caller)
      endif
   endif

end subroutine

!>
!! This routine prints out a logical-type setting.
!!
subroutine elsi_say_setting_log(e_h,label,setting,io_h)

   implicit none

   type(elsi_handle),    intent(in) :: e_h
   character(len=*),     intent(in) :: label
   logical,              intent(in) :: setting
   type(elsi_io_handle), intent(in) :: io_h

   character(len=28) :: label_ljust
   character(len=20) :: log_string

   character(len=40), parameter :: caller = "elsi_say_setting_log"

   if(io_h%file_format == HUMAN) then
      if(setting) then
         log_string = "                TRUE"
      else
         log_string = "               FALSE"
      endif
   elseif(io_h%file_format == JSON) then
      if(setting) then
         log_string = "                true"
      else
         log_string = "               false"
      endif
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

   label_ljust = label

   if(io_h%print_info .and. e_h%myid_all == 0) then
      if(io_h%file_format == HUMAN) then
         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,A25)") io_h%prefix,label_ljust,&
               " : ",log_string
         else
            write(io_h%print_unit,"(A28,A3,A25)") label_ljust," : ",log_string
         endif
      elseif(io_h%file_format == JSON) then
         if(io_h%comma_json == COMMA_AFTER) then
            if(allocated(io_h%prefix)) then
               write(io_h%print_unit,"(A)") io_h%prefix // '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(log_string)) // ","
            else
               write(io_h%print_unit,"(A)") '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(log_string)) // ","
            endif
         else
            if(allocated(io_h%prefix)) then
               write(io_h%print_unit,"(A)") io_h%prefix // '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(log_string))
            else
               write(io_h%print_unit,"(A)") '"' // &
                  trim(adjustl(label_ljust)) // '": ' // &
                  trim(adjustl(log_string))
            endif
         endif
      else
         call elsi_stop(e_h,"Unsupported output format.",caller)
      endif
   endif

end subroutine

!>
!! This routine prints out a string-type setting.
!!
subroutine elsi_say_setting_str(e_h,label,setting,io_h)

   implicit none

   type(elsi_handle),    intent(in) :: e_h
   character(len=*),     intent(in) :: label
   character(len=*),     intent(in) :: setting
   type(elsi_io_handle), intent(in) :: io_h

   character(len=28) :: label_ljust

   character(len=40), parameter :: caller = "elsi_say_setting_str"

   label_ljust = label

   if(io_h%print_info .and. e_h%myid_all == 0) then
      if(io_h%file_format == HUMAN) then
         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,A25)") io_h%prefix,label_ljust,&
               " : ",setting
         else
            write(io_h%print_unit,"(A28,A3,A25)") label_ljust," : ",setting
         endif
      elseif(io_h%file_format == JSON) then
         if(io_h%comma_json == COMMA_AFTER) then
            if(allocated(io_h%prefix)) then
               write(io_h%print_unit,"(A)") io_h%prefix // '"' // &
                  trim(adjustl(label_ljust)) // '": "' // &
                  trim(adjustl(setting)) // '",'
            else
               write(io_h%print_unit,"(A)") '"' // &
                  trim(adjustl(label_ljust)) // '": "' // &
                  trim(adjustl(setting)) // '",'
            endif
         else
            if(allocated(io_h%prefix)) then
               write(io_h%print_unit,"(A)") io_h%prefix // '"' // &
                  trim(adjustl(label_ljust)) // '": "' // &
                  trim(adjustl(setting)) // '"'
            else
               write(io_h%print_unit,"(A)") '"' // &
                  trim(adjustl(label_ljust)) // '": "' // &
                  trim(adjustl(setting)) // '"'
            endif
         endif
      else
         call elsi_stop(e_h,"Unsupported output format.",caller)
      endif
   endif

end subroutine

!>
!! This routine generates a new (dynamic) string with another string appended to
!! the end.
!!
subroutine elsi_append_string(l_string,r_string)

   implicit none

   character(len=:), intent(inout), allocatable :: l_string
   character(len=*), intent(in)                 :: r_string

   character(len=:), allocatable :: tmp_string

   character(len=40), parameter :: caller = "elsi_append_string"

   if(allocated(l_string)) then
      tmp_string = l_string // r_string
   else
      tmp_string = r_string
   endif

   if(allocated(l_string)) then
      deallocate(l_string)
   endif

   l_string = tmp_string

   deallocate(tmp_string)

end subroutine

!>
!! This routine generates a new string with the indicated number of characters
!! removed from the end.
!!
subroutine elsi_truncate_string(l_string,n_chars_to_remove)

   implicit none

   character(len=:), intent(inout), allocatable :: l_string
   integer(kind=i4), intent(in)                 :: n_chars_to_remove

   integer(kind=i4) :: size_new_string

   character(len=:), allocatable :: tmp_string

   character(len=40), parameter :: caller = "elsi_truncate_string"

   if(allocated(l_string)) then
      size_new_string = len(l_string)-n_chars_to_remove
   else
      return
   endif

   if(size_new_string < 1) then
      deallocate(l_string)
      return
   endif

   tmp_string = l_string(1:size_new_string)

   deallocate(l_string)

   l_string = tmp_string

   deallocate(tmp_string)

end subroutine

!>
!! This routine generates a file IO handle for the JSON file and opens the file
!!
subroutine elsi_open_json_file(e_h,print_unit,file_name,opening_bracket,io_h)

   implicit none

   type(elsi_handle),    intent(in)  :: e_h
   integer(kind=i4),     intent(in)  :: print_unit
   character(len=*),     intent(in)  :: file_name
   logical,              intent(in)  :: opening_bracket
   type(elsi_io_handle), intent(out) :: io_h

   character(len=40), parameter :: caller = "elsi_open_json_file"

   call elsi_init_io(io_h,print_unit,file_name,JSON,.true.,"",COMMA_AFTER)

   open(unit=io_h%print_unit,file=io_h%file_name)

   if(opening_bracket) then
      call elsi_say(e_h,"[",io_h)
      call elsi_append_string(io_h%prefix,"  ")
   endif

end subroutine

!>
!! This routine closes the JSON file and tears down the file IO handle.
!!
subroutine elsi_close_json_file(e_h,closing_bracket,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   logical,              intent(in)    :: closing_bracket
   type(elsi_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_close_json_file"

   if(io_h%file_format /= JSON) then
      call elsi_stop(e_h,"This routine requires a file handle in JSON format.",&
              caller)
   endif

   if(closing_bracket) then
      call elsi_truncate_string(io_h%prefix,2)
      call elsi_say(e_h,"]",io_h)
   endif

   close(io_h%print_unit)

   call elsi_reset_io_handle(io_h)

end subroutine

!>
!! This routine starts a new record in the JSON file.
!!
subroutine elsi_start_json_record(e_h,comma_before,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   logical,              intent(in)    :: comma_before
   type(elsi_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_start_json_record"

   if(io_h%file_format /= JSON) then
      call elsi_stop(e_h,"This routine requires a file handle in JSON format.",&
              caller)
   endif

   if(comma_before) then
      call elsi_say(e_h,',{',io_h)
   else
      call elsi_say(e_h,'{',io_h)
   endif

   call elsi_append_string(io_h%prefix,"  ")

end subroutine

!>
!! This routine finishes the current record in the JSON file.
!!
subroutine elsi_finish_json_record(e_h,comma_after,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   logical,              intent(in)    :: comma_after
   type(elsi_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_finish_json_record"

   if(io_h%file_format /= JSON) then
      call elsi_stop(e_h,"This routine requires a file handle in JSON format.",&
              caller)
   endif

   call elsi_truncate_string(io_h%prefix,2)

   if(comma_after) then
      call elsi_say(e_h,'},',io_h)
   else
      call elsi_say(e_h,'}',io_h)
   endif

end subroutine

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
!! This routine returns the current date and time, formatted as an RFC3339
!! string with time zone offset.
!!
subroutine elsi_get_datetime_rfc3339(dt)

   implicit none

   character(len=TIME_LEN), intent(out) :: dt

   integer(kind=i4) :: datetime(8)
   integer(kind=i4) :: tmp_int
   character(len=4) :: year
   character(len=2) :: month
   character(len=2) :: day
   character(len=2) :: hour
   character(len=2) :: minute
   character(len=2) :: second
   character(len=3) :: millisecond
   character(len=1) :: timezone_sign
   character(len=2) :: timezone_hour
   character(len=2) :: timezone_min

   character(len=40), parameter :: caller = "elsi_get_datetime_rfc3339"

   call date_and_time(values=datetime)

   ! Get year
   if(datetime(1) < 10) then
      write(year,'(A3,I1)') "000",datetime(1)
   elseif(datetime(1) < 100) then
      write(year,'(A2,I2)') "00",datetime(1)
   elseif(datetime(1) < 1000) then
      write(year,'(A1,I3)') "0",datetime(1)
   else
      write(year,'(I4)') datetime(1)
   endif

   ! Get month
   if(datetime(2) < 10) then
      write(month,'(A1,I1)') "0",datetime(2)
   else
      write(month,'(I2)') datetime(2)
   endif

   ! Get day
   if(datetime(3) < 10) then
      write(day,'(A1,I1)') "0",datetime(3)
   else
      write(day,'(I2)') datetime(3)
   endif

   ! Get hour
   if(datetime(5) < 10) then
      write(hour,'(A1,I1)') "0",datetime(5)
   else
      write(hour,'(I2)') datetime(5)
   endif

   ! Get minute
   if(datetime(6) < 10) then
      write(minute,'(A1,I1)') "0",datetime(6)
   else
      write(minute,'(I2)') datetime(6)
   endif

   ! Get second
   if(datetime(7) < 10) then
      write(second,'(A1,I1)') "0",datetime(7)
   else
      write(second,'(I2)') datetime(7)
   endif

   ! Get millisecond
   if(datetime(8) < 10) then
      write(millisecond,'(A2,I1)') "00",datetime(8)
   elseif(datetime(8) < 100) then
      write(millisecond,'(A1,I2)') "0",datetime(8)
   else
      write(millisecond,'(I3)') datetime(8)
   endif

   ! Get time zone sign (ahead or behind UTC)
   if(datetime(4) < 0) then
      timezone_sign = "-"
      datetime(4)   = -1*datetime(4)
   else
      timezone_sign = "+"
   endif

   ! Get timezone minutes
   tmp_int = mod(datetime(4),60)
   if(tmp_int < 10) then
      write(timezone_min,'(A1,I1)') "0",tmp_int
   else
      write(timezone_min,'(I2)') tmp_int
   endif

   ! Get timezone hours
   tmp_int = datetime(4)/60
   if(tmp_int < 10) then
      write(timezone_hour,'(A1,I1)') "0",tmp_int
   else
      write(timezone_hour,'(I2)') tmp_int
   endif

   write(dt,'(A4,A1,A2,A1,A2,A1,A2,A1,A2,A1,A2,A1,A3,A1,A2,A1,A2)')&
      year,"-",month,"-",day,"T",hour,":",minute,":",second,".",millisecond,&
      timezone_sign,timezone_hour,":",timezone_min

end subroutine

!>
!! This routine generates a string identifying the current solver.
!!
subroutine elsi_get_solver_tag(e_h,data_type,tag)

   implicit none

   type(elsi_handle),      intent(in)  :: e_h
   integer(kind=i4),       intent(in)  :: data_type
   character(len=STR_LEN), intent(out) :: tag

   character(len=40), parameter :: caller = "elsi_get_solver_tag"

   if(data_type == REAL_DATA) then
      select case(e_h%solver)
      case(ELPA_SOLVER)
         if(e_h%parallel_mode == SINGLE_PROC) then
            tag = "LAPACK_REAL"
         else ! MULTI_PROC
            if(e_h%elpa_solver == 1) then
               tag = "ELPA_1STAGE_REAL"
            else
               tag = "ELPA_2STAGE_REAL"
            endif
         endif
      case(OMM_SOLVER)
         tag = "LIBOMM_REAL"
      case(PEXSI_SOLVER)
         tag = "PEXSI_REAL"
      case(SIPS_SOLVER)
         tag = "SIPS_REAL"
      case(DMP_SOLVER)
         tag = "DMP_REAL"
      end select
   else
      select case(e_h%solver)
      case(ELPA_SOLVER)
         if(e_h%parallel_mode == SINGLE_PROC) then
            tag = "LAPACK_CMPLX"
         else ! MULTI_PROC
            if(e_h%elpa_solver == 1) then
               tag = "ELPA_1STAGE_CMPLX"
            else
               tag = "ELPA_2STAGE_CMPLX"
            endif
         endif
      case(OMM_SOLVER)
         tag = "LIBOMM_CMPLX"
      case(PEXSI_SOLVER)
         tag = "PEXSI_CMPLX"
      case(SIPS_SOLVER)
         tag = "SIPS_CMPLX"
      case(DMP_SOLVER)
         tag = "DMP_CMPLX"
      end select
   endif

end subroutine

!>
!! Generate a UUID (unique identifier) in RFC 4122 format.
!! (Taken from FHI-aims with permission of copyright holders)
!!
subroutine elsi_gen_uuid(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   integer(kind=i4) :: ii3
   integer(kind=i4) :: ii4
   integer(kind=i4) :: i_entry
   integer(kind=i4) :: ierr
   real(kind=r8)    :: rr(8)
   character(len=3) :: ss3
   character(len=4) :: ss4

   character(len=40), parameter :: caller = "elsi_gen_uuid"

   call elsi_init_random_seed()

   ii3 = 4095
   ii4 = 65535

   call random_number(rr)

   do i_entry=1,8
      write(ss3,"(Z3.3)") transfer(int(rr(i_entry)*ii3),16)
      write(ss4,"(Z4.4)") transfer(int(rr(i_entry)*ii4),16)

      if(i_entry == 1) then
         write(e_h%uuid,'(A)') ss4
      elseif(i_entry == 2) then
         write(e_h%uuid,'(2A)') trim(e_h%uuid),ss4
      elseif(i_entry == 3) then
         write(e_h%uuid,'(3A)') trim(e_h%uuid),'-',ss4
      elseif(i_entry == 4) then
         write(e_h%uuid,'(3A)') trim(e_h%uuid),'-4',ss3
      elseif(i_entry == 5) then
         write(e_h%uuid,'(4A)') trim(e_h%uuid),'-A',ss3,'-'
      else
         write(e_h%uuid,'(2A)') trim(e_h%uuid),ss4
      endif
   enddo

   if(e_h%parallel_mode == MULTI_PROC) then
      call MPI_Bcast(e_h%uuid,UUID_LEN,mpi_character,0,e_h%mpi_comm_all,ierr)
      call elsi_check_mpi(e_h,"MPI_Bcast",ierr,caller)
   endif

end subroutine

!>
!! Linear congruential generator.
!! (Taken from FHI-aims with permission of copyright holders)
!!
integer(kind=i4) function lcg(s)

   implicit none

   integer(kind=i8) :: s

   if(s == 0) then
      s = 104729
   else
      s = mod(s,int(huge(0_2)*2,kind=i8))
   endif

   s   = mod(s*int(huge(0_2),kind=i8),int(huge(0_2)*2,kind=i8))
   lcg = int(mod(s,int(huge(0),kind=i8)),kind(0))

end function

!>
!! Set the seed in the built-in random_seed subroutine using the system clock
!! modified by lcg().
!! (Taken from FHI-aims with permission of copyright holders)
!!
subroutine elsi_init_random_seed()

   implicit none

   integer(kind=i4) :: i
   integer(kind=i4) :: n
   integer(kind=i4) :: dt(8)
   integer(kind=i8) :: t

   integer(kind=i4), allocatable :: seed(:)

   character(len=40), parameter :: caller = "elsi_init_random_seed"

   call random_seed(size=n)
   call date_and_time(values=dt)

   t = (dt(1)-1970)*365*24*60*60*1000+dt(2)*31*24*60*60*1000+&
          dt(3)*24*60*60*1000+dt(5)*60*60*1000+dt(6)*60*1000+dt(7)*1000+dt(8)

   allocate(seed(n))

   ! Writting the array with seeds
   do i = 1,n
      seed(i) = lcg(t)
   enddo

   call random_seed(put=seed)

   deallocate(seed)

end subroutine

end module ELSI_IO
