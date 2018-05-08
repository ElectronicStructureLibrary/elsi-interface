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
                             BLACS_DENSE,PEXSI_CSC,SIESTA_CSC,STR_LEN,&
                             ELPA_SOLVER,PEXSI_SOLVER,SIPS_SOLVER,OMM_SOLVER,&
                             DMP_SOLVER,FILENAME_LEN,DATETIME_LEN
   use ELSI_DATATYPE,  only: elsi_handle,elsi_io_handle
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi
   use ELSI_PRECISION, only: r8,i4,i8
   use ELSI_JSON,      only: elsi_append_string,elsi_truncate_string,&
                             elsi_write_json_name_value,elsi_reset_json_handle,&
                             elsi_start_json_name_object,elsi_start_json_array,&
                             elsi_finish_json_object,elsi_open_json_file,&
                             elsi_start_json_object,elsi_get_datetime_rfc3339

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
   public :: elsi_get_time

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
subroutine elsi_say(io_h,info_str)

   implicit none

   type(elsi_io_handle), intent(in) :: io_h
   character(len=*),     intent(in) :: info_str

   character(len=40), parameter :: caller = "elsi_say"

   if(io_h%print_info) then
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
subroutine elsi_init_io(io_h,print_unit,file_name,file_format,print_info,prefix)

   implicit none

   type(elsi_io_handle), intent(out) :: io_h
   integer(kind=i4),     intent(in)  :: print_unit
   character(len=*),     intent(in)  :: file_name
   integer(kind=i4),     intent(in)  :: file_format
   logical,              intent(in)  :: print_info
   character(len=*),     intent(in)  :: prefix

   character(len=40), parameter :: caller = "elsi_init_io"

   ! For safety
   call elsi_reset_io_handle(io_h)
   call elsi_reset_json_handle(io_h%json_h)

   io_h%handle_init = .true.
   io_h%print_unit  = print_unit
   io_h%file_name   = file_name
   io_h%file_format = file_format
   io_h%print_info  = print_info
   io_h%prefix      = prefix

end subroutine

!>
!! This routine resets a handle.
!!
subroutine elsi_reset_io_handle(io_h)

   implicit none

   type(elsi_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_reset_io_handle"

   call elsi_reset_json_handle(io_h%json_h)

   io_h%handle_init = .false.
   io_h%file_name   = UNSET_STR
   io_h%file_format = UNSET
   io_h%print_info  = .false.
   io_h%print_unit  = UNSET
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
subroutine elsi_io_add_entry(e_h,dt0,t0,caller)

   implicit none

   type(elsi_handle),           intent(inout) :: e_h
   character(len=DATETIME_LEN), intent(in)    :: dt0
   real(kind=r8),               intent(in)    :: t0
   character(len=*),            intent(in)    :: caller

   real(kind=r8)               :: t1
   real(kind=r8)               :: t_total
   integer(kind=i4)            :: log_unit
   character(len=FILENAME_LEN) :: log_name
   character(len=STR_LEN)      :: solver_tag
   character(len=STR_LEN)      :: elsi_tag
   character(len=STR_LEN)      :: user_tag
   character(len=DATETIME_LEN) :: dt_record

   if(e_h%log_file%print_info .and. e_h%myid_all == 0) then
      if(.not. e_h%log_file%handle_init) then
         log_unit = e_h%log_file%print_unit
         log_name = e_h%log_file%file_name

         call elsi_open_json_file(e_h%log_file%json_h,log_unit,log_name)
         call elsi_start_json_array(e_h%log_file%json_h)

         if(.not. e_h%uuid_exists) then
            call elsi_gen_uuid(e_h)
            e_h%uuid_exists = .true.
         endif

         e_h%log_file%handle_init = .true.
      endif

      call elsi_get_time(t1)
      call elsi_get_solver_tag(e_h,solver_tag)
      call elsi_get_datetime_rfc3339(dt_record)

      t_total                = t1-t0
      e_h%log_file%n_records = e_h%log_file%n_records+1
      elsi_tag               = adjustr(trim(solver_tag))
      user_tag               = adjustr(trim(e_h%log_file%user_tag))

      call elsi_start_json_object(e_h%log_file%json_h)

      call elsi_print_versioning(e_h,e_h%log_file)
      call elsi_say_setting(e_h%log_file,"iteration",e_h%log_file%n_records)
      if(caller(6:6) == "e") then
         call elsi_say_setting(e_h%log_file,"output_type","EIGENSOLUTION")
      else
         call elsi_say_setting(e_h%log_file,"output_type","DENSITY MATRIX")
      endif
      if(caller(9:9) == "r") then
         call elsi_say_setting(e_h%log_file,"data_type","REAL")
      else
         call elsi_say_setting(e_h%log_file,"data_type","COMPLEX")
      endif
      call elsi_say_setting(e_h%log_file,"elsi_tag",trim(adjustl(elsi_tag)))
      call elsi_say_setting(e_h%log_file,"user_tag",trim(adjustl(user_tag)))
      call elsi_say_setting(e_h%log_file,"start_datetime",dt0)
      call elsi_say_setting(e_h%log_file,"record_datetime",dt_record)
      call elsi_say_setting(e_h%log_file,"total_time",t_total)

      call elsi_print_handle_summary(e_h,e_h%log_file)
      call elsi_print_matrix_settings(e_h,e_h%log_file)

      call elsi_print_solver_settings(e_h,e_h%log_file)

      call elsi_finish_json_object(e_h%log_file%json_h)
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
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_handle_summary"

   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Physical Properties"
      call elsi_say(io_h,info_str)

      call elsi_append_string(io_h%prefix,"  ")
      call elsi_say_setting(io_h,"Number of electrons",e_h%n_electrons)
      if(e_h%parallel_mode == MULTI_PROC) then
         call elsi_say_setting(io_h,"Number of spins",e_h%n_spins)
         call elsi_say_setting(io_h,"Number of k-points",e_h%n_kpts)
      endif
      if(e_h%solver == ELPA_SOLVER .or. e_h%solver == SIPS_SOLVER) then
         call elsi_say_setting(io_h,"Number of states",e_h%n_states)
      endif
      call elsi_truncate_string(io_h%prefix,2)

      write(info_str,"(A)") ""
      call elsi_say(io_h,info_str)
      write(info_str,"(A)") "Matrix Properties"
      call elsi_say(io_h,info_str)

      call elsi_append_string(io_h%prefix,"  ")
      if(e_h%matrix_format == BLACS_DENSE) then
         call elsi_say_setting(io_h,"Matrix format","BLACS_DENSE")
      elseif(e_h%matrix_format == PEXSI_CSC) then
         call elsi_say_setting(io_h,"Matrix format","PEXSI_CSC")
      elseif(e_h%matrix_format == SIESTA_CSC) then
         call elsi_say_setting(io_h,"Matrix format","SIESTA_CSC")
      endif
      call elsi_say_setting(io_h,"Number of basis functions",e_h%n_basis)
      if(e_h%parallel_mode == MULTI_PROC) then
         sparsity = 1.0_r8-(1.0_r8*e_h%nnz_g/e_h%n_basis/e_h%n_basis)
         call elsi_say_setting(io_h,"Matrix sparsity",sparsity)
      endif
      call elsi_truncate_string(io_h%prefix,2)

      write(info_str,"(A)") ""
      call elsi_say(io_h,info_str)
      write(info_str,"(A)") "Computational Details"
      call elsi_say(io_h,info_str)

      call elsi_append_string(io_h%prefix,"  ")
      if(e_h%parallel_mode == MULTI_PROC) then
         call elsi_say_setting(io_h,"Parallel mode","MULTI_PROC")
      elseif(e_h%parallel_mode == SINGLE_PROC) then
         call elsi_say_setting(io_h,"Parallel mode","SINGLE_PROC")
      endif
      call elsi_say_setting(io_h,"Number of MPI tasks",e_h%n_procs)
      if(e_h%solver == ELPA_SOLVER) then
         call elsi_say_setting(io_h,"Solver requested","ELPA")
      elseif(e_h%solver == OMM_SOLVER) then
         call elsi_say_setting(io_h,"Solver requested","libOMM")
      elseif(e_h%solver == PEXSI_SOLVER) then
         call elsi_say_setting(io_h,"Solver requested","PEXSI")
      elseif(e_h%solver == SIPS_SOLVER) then
         call elsi_say_setting(io_h,"Solver requested","SLEPc_SIPs")
      elseif(e_h%solver == DMP_SOLVER) then
         call elsi_say_setting(io_h,"Solver requested","DMP")
      else
         call elsi_stop(e_h,"Unsupported solver.",caller)
      endif
      call elsi_say_setting(io_h,"Number of ELSI calls",e_h%n_elsi_calls)
      call elsi_truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      call elsi_say_setting(io_h,"n_electrons",e_h%n_electrons)
      if(e_h%parallel_mode == MULTI_PROC) then
         call elsi_say_setting(io_h,"n_spin",e_h%n_spins)
         call elsi_say_setting(io_h,"n_kpts",e_h%n_kpts)
      endif
      if(e_h%solver == ELPA_SOLVER .or. e_h%solver == SIPS_SOLVER) then
         call elsi_say_setting(io_h,"n_states",e_h%n_states)
      endif

      if(e_h%matrix_format == BLACS_DENSE) then
         call elsi_say_setting(io_h,"matrix_format","BLACS_DENSE")
      elseif(e_h%matrix_format == PEXSI_CSC) then
         call elsi_say_setting(io_h,"matrix_format","PEXSI_CSC")
      elseif(e_h%matrix_format == SIESTA_CSC) then
         call elsi_say_setting(io_h,"matrix_format","SIESTA_CSC")
      endif
      call elsi_say_setting(io_h,"n_basis",e_h%n_basis)
      if(e_h%parallel_mode == MULTI_PROC) then
         sparsity = 1.0_r8-(1.0_r8*e_h%nnz_g/e_h%n_basis/e_h%n_basis)
         call elsi_say_setting(io_h,"sparsity",sparsity)
      endif

      if(e_h%parallel_mode == MULTI_PROC) then
         call elsi_say_setting(io_h,"parallel_mode","MULTI_PROC")
      elseif(e_h%parallel_mode == SINGLE_PROC) then
         call elsi_say_setting(io_h,"parallel_mode","SINGLE_PROC")
      endif
      call elsi_say_setting(io_h,"n_procs",e_h%n_procs)
      if(e_h%solver == ELPA_SOLVER) then
         call elsi_say_setting(io_h,"solver","ELPA")
      elseif(e_h%solver == OMM_SOLVER) then
         call elsi_say_setting(io_h,"solver","libOMM")
      elseif(e_h%solver == PEXSI_SOLVER) then
         call elsi_say_setting(io_h,"solver","PEXSI")
      elseif(e_h%solver == SIPS_SOLVER) then
         call elsi_say_setting(io_h,"solver","SLEPc_SIPs")
      elseif(e_h%solver == DMP_SOLVER) then
         call elsi_say_setting(io_h,"solver","DMP")
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
      call elsi_say(io_h,info_str)
      call elsi_append_string(io_h%prefix,"  ")

      call elsi_say_setting(io_h,"Date stamp",trim(adjustl(DATE_STAMP)))
      call elsi_say_setting(io_h,"Git commit (abbrev.)",&
              trim(adjustl(COMMIT_ABBREV)))
      call elsi_say_setting(io_h,"Git commit modified?",COMMIT_MODIFIED)

      call elsi_truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      call elsi_write_json_name_value(io_h%json_h,"data_source","ELSI")
      call elsi_write_json_name_value(io_h%json_h,"date_stamp",&
              trim(adjustl(DATE_STAMP)))
      call elsi_write_json_name_value(io_h%json_h,"git_commit",&
              trim(adjustl(COMMIT)))
      call elsi_write_json_name_value(io_h%json_h,"git_commit_modified",&
              COMMIT_MODIFIED)
      call elsi_write_json_name_value(io_h%json_h,"git_message_abbrev",&
              trim(adjustl(COMMIT_MSG_ABBREV)))
      call elsi_write_json_name_value(io_h%json_h,"source_created_on_hostname",&
              trim(adjustl(HOSTNAME)))
      call elsi_write_json_name_value(io_h%json_h,"source_created_at_datetime",&
              trim(adjustl(DATETIME)))
      call elsi_write_json_name_value(io_h%json_h,"calling_code",&
              trim(adjustl(e_h%caller)))
      call elsi_write_json_name_value(io_h%json_h,"uuid",trim(adjustl(e_h%uuid)))
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

   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_dmp_settings"

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Solver Settings (DMP)"
      call elsi_say(io_h,info_str)
      call elsi_append_string(io_h%prefix,"  ")
   elseif(io_h%file_format == JSON) then
      call elsi_start_json_name_object(io_h%json_h, "solver_settings")
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

   ! Settings
   call elsi_say_setting(io_h,"dmp_n_states",e_h%dmp_n_states)
   call elsi_say_setting(io_h,"dmp_method",e_h%dmp_method)
   call elsi_say_setting(io_h,"dmp_max_power",e_h%dmp_max_power)
   call elsi_say_setting(io_h,"dmp_max_iter",e_h%dmp_max_iter)
   call elsi_say_setting(io_h,"dmp_tol",e_h%dmp_tol)

   ! Footer
   if(io_h%file_format == HUMAN) then
      call elsi_truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      call elsi_finish_json_object(io_h%json_h)
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

end subroutine

!>
!! This routine prints out settings for ELPA.
!!
subroutine elsi_print_elpa_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_elpa_settings"

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Solver Settings (ELPA)"
      call elsi_say(io_h,info_str)
      call elsi_append_string(io_h%prefix,"  ")
   elseif(io_h%file_format == JSON) then
      call elsi_start_json_name_object(io_h%json_h, "solver_settings")
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

   ! Settings
   call elsi_say_setting(io_h,"elpa_solver",e_h%elpa_solver)
   call elsi_say_setting(io_h,"elpa_n_states",e_h%n_states)
   call elsi_say_setting(io_h,"elpa_gpu",e_h%elpa_gpu)
   call elsi_say_setting(io_h,"elpa_gpu_kernels",e_h%elpa_gpu_kernels)

   ! Footer
   if(io_h%file_format == HUMAN) then
      call elsi_truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      call elsi_finish_json_object(io_h%json_h)
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

end subroutine

!>
!! This routine prints out settings for libOMM.
!!
subroutine elsi_print_omm_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_omm_settings"

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Solver Settings (libOMM)"
      call elsi_say(io_h,info_str)
      call elsi_append_string(io_h%prefix,"  ")
   elseif(io_h%file_format == JSON) then
      call elsi_start_json_name_object(io_h%json_h, "solver_settings")
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

  ! Settings
   call elsi_say_setting(io_h,"omm_n_states",e_h%omm_n_states)
   call elsi_say_setting(io_h,"omm_n_elpa",e_h%omm_n_elpa)
   call elsi_say_setting(io_h,"omm_flavor",e_h%omm_flavor)
   call elsi_say_setting(io_h,"omm_tol",e_h%omm_tol)

   ! Footer
   if(io_h%file_format == HUMAN) then
      call elsi_truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      call elsi_finish_json_object(io_h%json_h)
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

end subroutine

!>
!! This routine prints out settings for PEXSI.
!!
subroutine elsi_print_pexsi_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_pexsi_settings"

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Solver Settings (PEXSI)"
      call elsi_say(io_h,info_str)
      call elsi_append_string(io_h%prefix,"  ")
   elseif(io_h%file_format == JSON) then
      call elsi_start_json_name_object(io_h%json_h, "solver_settings")
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

   ! Settings
   call elsi_say_setting(io_h,"pexsi_n_pole",e_h%pexsi_options%numPole)
   call elsi_say_setting(io_h,"pexsi_n_point",e_h%pexsi_options%nPoints)
   call elsi_say_setting(io_h,"pexsi_np_per_pole",e_h%pexsi_np_per_pole)
   call elsi_say_setting(io_h,"pexsi_np_per_point",e_h%pexsi_np_per_point)
   call elsi_say_setting(io_h,"pexsi_n_prow_pexsi",e_h%pexsi_n_prow)
   call elsi_say_setting(io_h,"pexsi_n_pcol_pexsi",e_h%pexsi_n_pcol)
   call elsi_say_setting(io_h,"pexsi_temperature",e_h%pexsi_options%temperature)
   call elsi_say_setting(io_h,"pexsi_delta_e",e_h%pexsi_options%deltaE)
   call elsi_say_setting(io_h,"pexsi_gap",e_h%pexsi_options%gap)
   call elsi_say_setting(io_h,"pexsi_mu_min",e_h%pexsi_options%muMin0)
   call elsi_say_setting(io_h,"pexsi_mu_max",e_h%pexsi_options%muMax0)
   call elsi_say_setting(io_h,"pexsi_do_inertia",&
           e_h%pexsi_options%isInertiaCount)
   call elsi_say_setting(io_h,"pexsi_intertia_tol",&
           e_h%pexsi_options%muInertiaTolerance)
   call elsi_say_setting(io_h,"pexsi_tol",&
           e_h%pexsi_options%numElectronPEXSITolerance)
   call elsi_say_setting(io_h,"pexsi_np_symbfact",e_h%pexsi_options%npSymbFact)
   call elsi_say_setting(io_h,"pexsi_reordering",e_h%pexsi_options%ordering)

   ! Footer
   if(io_h%file_format == HUMAN) then
      call elsi_truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      call elsi_finish_json_object(io_h%json_h)
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

end subroutine

!>
!! This routine prints out settings for SIPS.
!!
subroutine elsi_print_sips_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_sips_settings"

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Solver Settings (SLEPc-SIPs)"
      call elsi_say(io_h,info_str)
      call elsi_append_string(io_h%prefix,"  ")
   elseif(io_h%file_format == JSON) then
      call elsi_start_json_name_object(io_h%json_h, "solver_settings")
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

   ! Settings
   call elsi_say_setting(io_h,"sips_n_states",e_h%n_states)
   call elsi_say_setting(io_h,"sips_n_elpa",e_h%sips_n_elpa)
   call elsi_say_setting(io_h,"sips_slice_type",e_h%sips_slice_type)
   call elsi_say_setting(io_h,"sips_n_slices",e_h%sips_n_slices)
   call elsi_say_setting(io_h,"sips_np_per_slice",e_h%sips_np_per_slice)
   call elsi_say_setting(io_h,"sips_buffer",e_h%sips_buffer)

   ! Footer
   if(io_h%file_format == HUMAN) then
      call elsi_truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      call elsi_finish_json_object(io_h%json_h)
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
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

   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_den_settings"

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Dense Matrix Format Settings"
      call elsi_say(io_h,info_str)
      call elsi_append_string(io_h%prefix,"  ")
   elseif(io_h%file_format == JSON) then
      call elsi_start_json_name_object(io_h%json_h, "matrix_format_settings")
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

   ! Settings
   call elsi_say_setting(io_h,"blk_row",e_h%blk_row)
   call elsi_say_setting(io_h,"blk_col",e_h%blk_col)
   call elsi_say_setting(io_h,"n_prow",e_h%n_prow)
   call elsi_say_setting(io_h,"n_pcol",e_h%n_pcol)
   call elsi_say_setting(io_h,"blacs_ready",e_h%blacs_ready)

   ! Footer
   if(io_h%file_format == HUMAN) then
      call elsi_truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      call elsi_finish_json_object(io_h%json_h)
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

end subroutine

!>
!! This routine prints out settings for the sparse matrix format.
!!
subroutine elsi_print_csc_settings(e_h,io_h)

   implicit none

   type(elsi_handle),    intent(in)    :: e_h
   type(elsi_io_handle), intent(inout) :: io_h

   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_csc_settings"

   ! Header
   if(io_h%file_format == HUMAN) then
      write(info_str,"(A)") "Sparse Matrix Format Settings"
      call elsi_say(io_h,info_str)
      call elsi_append_string(io_h%prefix,"  ")
   elseif(io_h%file_format == JSON) then
      call elsi_start_json_name_object(io_h%json_h, "matrix_format_settings")
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

   ! Settings
   call elsi_say_setting(io_h,"zero_def",e_h%zero_def)
   call elsi_say_setting(io_h,"blk_sp2",e_h%blk_sp2)
   call elsi_say_setting(io_h,"pexsi_csc_ready",e_h%pexsi_csc_ready)
   call elsi_say_setting(io_h,"siesta_csc_ready",e_h%siesta_csc_ready)

   ! Footer
   if(io_h%file_format == HUMAN) then
      call elsi_truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      call elsi_finish_json_object(io_h%json_h)
   else
      call elsi_stop(e_h,"Unsupported output format.",caller)
   endif

end subroutine

!>
!! This routine prints out an integer-type setting.
!!
subroutine elsi_say_setting_i4(io_h,label,setting)

   implicit none

   type(elsi_io_handle), intent(inout) :: io_h
   character(len=*),     intent(in)    :: label
   integer(kind=i4),     intent(in)    :: setting

   character(len=28) :: label_ljust
   character(len=20) :: int_string

   character(len=40), parameter :: caller = "elsi_say_setting_i4"

   if(io_h%print_info) then
      if(io_h%file_format == HUMAN) then
         write(int_string,"(I20)") setting
         label_ljust = label

         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,I25)") io_h%prefix,label_ljust,&
               " : ",setting
         else
            write(io_h%print_unit,"(A28,A3,I25)") label_ljust," : ",setting
         endif
      elseif(io_h%file_format == JSON) then
         call elsi_write_json_name_value(io_h%json_h, label, setting)
      endif
   endif

end subroutine

!>
!! This routine prints out a real-type setting.
!!
subroutine elsi_say_setting_r8(io_h,label,setting)

   implicit none

   type(elsi_io_handle), intent(inout) :: io_h
   character(len=*),     intent(in)    :: label
   real(kind=r8),        intent(in)    :: setting

   character(len=28) :: label_ljust
   character(len=20) :: real_string

   character(len=40), parameter :: caller = "elsi_say_setting_r8"

   if(io_h%print_info) then
      if(io_h%file_format == HUMAN) then
         write(real_string,"(E20.8)") setting
         label_ljust = label

         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,E25.8)") io_h%prefix,label_ljust,&
               " : ",setting
         else
            write(io_h%print_unit,"(A28,A3,E25.8)") label_ljust," : ",setting
         endif
      elseif(io_h%file_format == JSON) then
         call elsi_write_json_name_value(io_h%json_h, label, setting)
      endif
   endif

end subroutine

!>
!! This routine prints out a logical-type setting.
!!
subroutine elsi_say_setting_log(io_h,label,setting)

   implicit none

   type(elsi_io_handle), intent(inout) :: io_h
   character(len=*),     intent(in)    :: label
   logical,              intent(in)    :: setting

   character(len=28) :: label_ljust
   character(len=20) :: log_string

   character(len=40), parameter :: caller = "elsi_say_setting_log"

   if(io_h%print_info) then
      if(io_h%file_format == HUMAN) then
         if(setting) then
            log_string = "                TRUE"
         else
            log_string = "               FALSE"
         endif
         label_ljust = label

         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,A25)") io_h%prefix,label_ljust,&
               " : ",log_string
         else
            write(io_h%print_unit,"(A28,A3,A25)") label_ljust," : ",log_string
         endif
      elseif(io_h%file_format == JSON) then
         call elsi_write_json_name_value(io_h%json_h, label, setting)
      endif
   endif

end subroutine

!>
!! This routine prints out a string-type setting.
!!
subroutine elsi_say_setting_str(io_h,label,setting)

   implicit none

   type(elsi_io_handle), intent(inout) :: io_h
   character(len=*),     intent(in)    :: label
   character(len=*),     intent(in)    :: setting

   character(len=28) :: label_ljust

   character(len=40), parameter :: caller = "elsi_say_setting_str"

   if(io_h%print_info) then
      if(io_h%file_format == HUMAN) then
         label_ljust = label

         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,A25)") io_h%prefix,label_ljust,&
               " : ",setting
         else
            write(io_h%print_unit,"(A28,A3,A25)") label_ljust," : ",setting
         endif
      elseif(io_h%file_format == JSON) then
         call elsi_write_json_name_value(io_h%json_h, label, setting)
      endif
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

   read(cdate(1:4),"(I4)") val

   int_year = val
   year     = real(val,kind=r8)-2009.0_r8 ! 2009 is an arbitrary zero
   day      = year*365+floor(year/4.0_r8)

   read(cdate(5:6),"(I2)") val

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

   read(cdate(7:8),"(I2)") val
   day = day+real(val,kind=r8)-1

   read(ctime(1:2),"(I2)") val
   hour = real(val,kind=r8)

   read(ctime(3:4),"(I2)") val
   minute = real(val,kind=r8)

   read(ctime(5:6),"(I2)") val
   second = real(val,kind=r8)

   read(ctime(8:10),"(I3)") val
   millisecond = real(val,kind=r8)

   wtime = day*24.0_r8*3600.0_r8+hour*3600.0_r8+minute*60.0_r8+second+&
              millisecond*0.001_r8

end subroutine

!>
!! This routine generates a string identifying the current solver.
!!
subroutine elsi_get_solver_tag(e_h,tag)

   implicit none

   type(elsi_handle),      intent(in)  :: e_h
   character(len=STR_LEN), intent(out) :: tag

   character(len=40), parameter :: caller = "elsi_get_solver_tag"

   select case(e_h%solver)
   case(ELPA_SOLVER)
      if(e_h%parallel_mode == SINGLE_PROC) then
         tag = "LAPACK"
      else
         tag = "ELPA"
      endif
   case(OMM_SOLVER)
      if(e_h%n_elsi_calls <= e_h%omm_n_elpa) then
         tag = "ELPA"
      else
         tag = "LIBOMM"
      endif
   case(PEXSI_SOLVER)
      tag = "PEXSI"
   case(SIPS_SOLVER)
      if(e_h%n_elsi_calls <= e_h%sips_n_elpa) then
         tag = "ELPA"
      else
         tag = "SIPS"
      endif
   case(DMP_SOLVER)
      tag = "DMP"
   end select

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
         write(e_h%uuid,"(A)") ss4
      elseif(i_entry == 2) then
         write(e_h%uuid,"(2A)") trim(e_h%uuid),ss4
      elseif(i_entry == 3) then
         write(e_h%uuid,"(3A)") trim(e_h%uuid),"-",ss4
      elseif(i_entry == 4) then
         write(e_h%uuid,"(3A)") trim(e_h%uuid),"-4",ss3
      elseif(i_entry == 5) then
         write(e_h%uuid,"(4A)") trim(e_h%uuid),"-A",ss3,"-"
      else
         write(e_h%uuid,"(2A)") trim(e_h%uuid),ss4
      endif
   enddo

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