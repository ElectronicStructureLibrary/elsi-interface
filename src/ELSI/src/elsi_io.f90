! Copyright (c) 2015-2018, the ELSI team. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!  * Neither the name of the "ELectronic Structure Infrastructure" project nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!>
!! This module performs IO to stdout and files.
!!
module ELSI_IO

   use ELSI_CONSTANTS, only: UNSET,UNSET_STRING,HUMAN_READ,JSON,MULTI_PROC,&
                             SINGLE_PROC,ELPA_SOLVER,SIPS_SOLVER,OMM_SOLVER,&
                             PEXSI_SOLVER,CHESS_SOLVER,DMP_SOLVER,BLACS_DENSE,&
                             PEXSI_CSC,COMMA_AFTER,NO_COMMA
   use ELSI_DATATYPE,  only: elsi_handle,elsi_file_io_handle
   use ELSI_MPI,       only: elsi_stop
   use ELSI_PRECISION, only: r8,i4
   use ELSI_VERSION,   only: RELEASE_DATE, GIT_COMMIT, GIT_COMMIT_ABBREV, &
                             GIT_COMMIT_WAS_MODIFIED, GIT_COMMIT_MSG_ABBREV, &
                             SOURCE_HOSTNAME, SOURCE_LOCAL_DATE, &
                             SOURCE_LOCAL_TIME, SOURCE_DATETIME

   implicit none

   private

   public :: elsi_say
   public :: elsi_say_setting
   public :: elsi_init_file_io
   public :: elsi_reset_file_io_handle
   public :: elsi_finalize_file_io
   public :: elsi_print_handle_summary
   public :: elsi_print_versioning
   public :: elsi_print_settings
   public :: elsi_print_solver_settings
   public :: elsi_print_chess_settings
   public :: elsi_print_dmp_settings
   public :: elsi_print_elpa_settings
   public :: elsi_print_omm_settings
   public :: elsi_print_pexsi_settings
   public :: elsi_print_sips_settings
   public :: elsi_print_matrix_format_settings
   public :: elsi_print_blacs_dense_settings
   public :: elsi_print_pexsi_csc_settings
   public :: append_string
   public :: truncate_string
   public :: elsi_open_json_file
   public :: elsi_close_json_file
   public :: elsi_start_json_record
   public :: elsi_finish_json_record

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
subroutine elsi_say(e_h,info_str,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   character(len=*),          intent(in)           :: info_str
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   type(elsi_file_io_handle) :: io_h

   character(len=40), parameter :: caller = "elsi_init_io"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   if(io_h%print_info .and. e_h%myid_all == 0) then
      if(allocated(io_h%prefix)) then
         write(io_h%print_unit,"(A,A)") io_h%prefix,trim(info_str)
      else
         write(io_h%print_unit,"(A)") trim(info_str)
      endif
   endif

end subroutine

!>
!! This routine initializes a handle for reading and writing to files.
!!
subroutine elsi_init_file_io(io_h,print_unit,file_name,file_format,print_info,&
              prefix,comma_json)

   implicit none

   type(elsi_file_io_handle), intent(out)          :: io_h
   integer(kind=i4),          intent(in)           :: print_unit
   character(len=*),          intent(in), optional :: file_name
   integer(kind=i4),          intent(in), optional :: file_format
   logical,                   intent(in), optional :: print_info
   character(len=*),          intent(in), optional :: prefix
   integer(kind=i4),          intent(in), optional :: comma_json

   character(len=40), parameter :: caller = "elsi_init_io"

   ! For safety
   call elsi_reset_file_io_handle(io_h)

   io_h%handle_init = .true.
   io_h%print_unit  = print_unit

   if(present(file_name)) then
      io_h%file_name = file_name
   else
      io_h%file_name = UNSET_STRING
   endif

   if(present(file_format)) then
      io_h%file_format = file_format
   else
      io_h%file_format = HUMAN_READ
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
!! This routine finalizes a file io handle. Note this subroutine does NOT close
!! file units.
!!
subroutine elsi_finalize_file_io(e_h,io_h)

   implicit none

   type(elsi_handle),         intent(in)    :: e_h
   type(elsi_file_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_finalize_file_io"

   call elsi_check_file_io_handle(e_h,io_h,caller)
   call elsi_reset_file_io_handle(io_h)

end subroutine

!>
!! This routine checks whether a handle has been properly initialized for
!! reading and writing to files.
!!
subroutine elsi_check_file_io_handle(e_h,io_h,caller)

   implicit none

   type(elsi_handle),         intent(in) :: e_h
   type(elsi_file_io_handle), intent(in) :: io_h
   character(len=*),          intent(in) :: caller

   if(.not. io_h%handle_init) then
      call elsi_stop(" Invalid handle! Not initialized.",e_h,caller)
   endif

end subroutine

!>
!! This routine resets a handle.
!!
subroutine elsi_reset_file_io_handle(io_h)

   implicit none

   type(elsi_file_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_reset_file_io_handle"

   io_h%handle_init = .false.
   io_h%print_unit  = UNSET
   io_h%file_name   = UNSET_STRING
   io_h%file_format = UNSET
   io_h%print_info  = .false.
   io_h%comma_json  = UNSET

   if(allocated(io_h%prefix)) then
      deallocate(io_h%prefix)
   endif

end subroutine

!>
!! This routine prints the state of the handle
!!
subroutine elsi_print_handle_summary(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   real(kind=r8)             :: sparsity
   integer(kind=i4)          :: comma_json_save
   type(elsi_file_io_handle) :: io_h
   character(len=200)        :: info_str

   character(len=40), parameter :: caller = "elsi_print_handle_summary"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   if(io_h%file_format == HUMAN_READ) then
      write(info_str,"(A)") "Physical Properties"
      call elsi_say(e_h,info_str,io_h)

      call append_string(io_h%prefix,"  ")
      call elsi_say_setting(e_h,"Number of electrons",e_h%n_electrons,io_h)
      if(e_h%parallel_mode == MULTI_PROC) then
         call elsi_say_setting(e_h,"Number of spins",e_h%n_spins,io_h)
         call elsi_say_setting(e_h,"Number of k-points",e_h%n_kpts,io_h)
      endif
      if(e_h%solver == ELPA_SOLVER .or. e_h%solver == SIPS_SOLVER) then
         call elsi_say_setting(e_h,"Number of states",e_h%n_states,io_h)
      endif
      call truncate_string(io_h%prefix,2)

      write(info_str,"(A)") ""
      call elsi_say(e_h,info_str,io_h)
      write(info_str,"(A)") "Matrix Properties"
      call elsi_say(e_h,info_str,io_h)

      call append_string(io_h%prefix,"  ")
      if(e_h%matrix_format == BLACS_DENSE) then
         call elsi_say_setting(e_h,"Matrix format","BLACS_DENSE",io_h)
      elseif(e_h%matrix_format == PEXSI_CSC) then
         call elsi_say_setting(e_h,"Matrix format","PEXSI_CSC",io_h)
      endif
      call elsi_say_setting(e_h,"Number of basis functions",e_h%n_basis,io_h)
      if(e_h%parallel_mode == MULTI_PROC) then
         sparsity = 1.0_r8-(1.0_r8*e_h%nnz_g/e_h%n_basis/e_h%n_basis)
         call elsi_say_setting(e_h,"Matrix sparsity",sparsity,io_h)
      endif
      call truncate_string(io_h%prefix,2)

      write(info_str,"(A)") ""
      call elsi_say(e_h,info_str,io_h)
      write(info_str,"(A)") "Computational Details"
      call elsi_say(e_h,info_str,io_h)

      call append_string(io_h%prefix,"  ")
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
      elseif(e_h%solver == CHESS_SOLVER) then
         call elsi_say_setting(e_h,"Solver requested","CheSS",io_h)
      elseif(e_h%solver == SIPS_SOLVER) then
         call elsi_say_setting(e_h,"Solver requested","SIPs",io_h)
      elseif(e_h%solver == DMP_SOLVER) then
         call elsi_say_setting(e_h,"Solver requested","DMP",io_h)
      else
         call elsi_stop(" Unsupported solver.",e_h,caller)
      endif
      call truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      comma_json_save = io_h%comma_json
      io_h%comma_json = COMMA_AFTER ! Add commas behind all records before final

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
      io_h%comma_json = comma_json_save ! Final record, restore comma_json
      call elsi_say_setting(e_h,"n_procs",e_h%n_procs,io_h)
      if(e_h%solver == ELPA_SOLVER) then
         call elsi_say_setting(e_h,"solver","ELPA",io_h)
      elseif(e_h%solver == OMM_SOLVER) then
         call elsi_say_setting(e_h,"solver","libOMM",io_h)
      elseif(e_h%solver == PEXSI_SOLVER) then
         call elsi_say_setting(e_h,"solver","PEXSI",io_h)
      elseif(e_h%solver == CHESS_SOLVER) then
         call elsi_say_setting(e_h,"solver","CheSS",io_h)
      elseif(e_h%solver == SIPS_SOLVER) then
         call elsi_say_setting(e_h,"solver","SIPs",io_h)
      elseif(e_h%solver == DMP_SOLVER) then
         call elsi_say_setting(e_h,"solver","DMP",io_h)
      else
         call elsi_stop(" Unsupported solver.",e_h,caller)
      endif
   else
      call elsi_stop(" Unsupported output format.",e_h,caller)
   endif

end subroutine

!>
!! This routine prints versioning information
!!
subroutine elsi_print_versioning(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   integer(kind=i4)          :: comma_json_save
   type(elsi_file_io_handle) :: io_h
   character(len=200)        :: info_str

   character(len=40), parameter :: caller = "elsi_print_versioning"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   if(io_h%file_format == HUMAN_READ) then
      write(info_str,"(A)") "ELSI Versioning Information:"
      call elsi_say(e_h,info_str,io_h)
      call append_string(io_h%prefix,"  ")

      call elsi_say_setting(e_h,"ELSI release date",trim(RELEASE_DATE),io_h)
      call elsi_say_setting(e_h,"ELSI git commit (abbrev.)",trim(GIT_COMMIT_ABBREV),io_h)
      call elsi_say_setting(e_h,"Was git commit modified?",GIT_COMMIT_WAS_MODIFIED,io_h)
      call elsi_say_setting(e_h,"git commit message (abbrev.)",trim(GIT_COMMIT_MSG_ABBREV),io_h)
      call elsi_say_setting(e_h,"Source created on hostname",trim(SOURCE_HOSTNAME),io_h)
      call elsi_say_setting(e_h,"Source created at local date",trim(SOURCE_LOCAL_DATE),io_h)
      call elsi_say_setting(e_h,"Source created at local time",trim(SOURCE_LOCAL_TIME),io_h)
      call elsi_say_setting(e_h,"Name of code calling ELSI",trim(e_h%calling_code),io_h)
      call elsi_say_setting(e_h,"Version of code calling ELSI",trim(e_h%calling_code_ver),io_h)
      call elsi_say_setting(e_h,"UUID for this run",trim(e_h%uuid),io_h)

      call truncate_string(io_h%prefix,2)
   elseif(io_h%file_format == JSON) then
      comma_json_save = io_h%comma_json
      io_h%comma_json = COMMA_AFTER ! Add commas behind all records before final

      call elsi_say_setting(e_h,"data_source","ELSI",io_h)
      call elsi_say_setting(e_h,"release_date",RELEASE_DATE,io_h)
      call elsi_say_setting(e_h,"git_commit",GIT_COMMIT,io_h)
      call elsi_say_setting(e_h,"git_commit_modified",GIT_COMMIT_WAS_MODIFIED,io_h)
      call elsi_say_setting(e_h,"git_message_abbrev",GIT_COMMIT_MSG_ABBREV,io_h)
      call elsi_say_setting(e_h,"source_created_on_hostname",SOURCE_HOSTNAME,io_h)
      call elsi_say_setting(e_h,"source_created_at_datetime",SOURCE_DATETIME,io_h)
      call elsi_say_setting(e_h,"calling_code",e_h%calling_code,io_h)
      call elsi_say_setting(e_h,"calling_code_version",e_h%calling_code_ver,io_h)

      io_h%comma_json = comma_json_save ! Final record, restore comma_json
      call elsi_say_setting(e_h,"uuid",e_h%uuid,io_h)
   else
      call elsi_stop(" Unsupported output format.",e_h,caller)
   endif

end subroutine

!>
!! This routine prints ELSI settings.
!!
subroutine elsi_print_settings(e_h)

   implicit none

   type(elsi_handle), intent(in) :: e_h !< Handle

   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_print_settings"

   select case(e_h%solver)
   case(CHESS_SOLVER)
      call elsi_say(e_h,"  CheSS settings:")

      write(info_str,"('  | Error function decay length ',E10.2)")&
         e_h%chess_erf_decay
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Lower bound of decay length ',E10.2)")&
         e_h%chess_erf_min
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Upper bound of decay length ',E10.2)")&
         e_h%chess_erf_max
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Lower bound of H eigenvalue ',E10.2)")&
         e_h%chess_ev_ham_min
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Upper bound of H eigenvalue ',E10.2)")&
         e_h%chess_ev_ham_max
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Lower bound of S eigenvalue ',E10.2)")&
         e_h%chess_ev_ovlp_min
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Upper bound of S eigenvalue ',E10.2)") &
         e_h%chess_ev_ovlp_max
      call elsi_say(e_h,info_str)
   case(ELPA_SOLVER)
      if(e_h%parallel_mode == MULTI_PROC) then
         call elsi_say(e_h,"  ELPA settings:")

         write(info_str,"('  | ELPA solver ',I10)") e_h%elpa_solver
         call elsi_say(e_h,info_str)
      endif
   case(OMM_SOLVER)
      call elsi_say(e_h,"  libOMM settings:")

      write(info_str,"('  | Number of ELPA steps       ',I10)") e_h%omm_n_elpa
      call elsi_say(e_h,info_str)

      write(info_str,"('  | OMM minimization flavor    ',I10)") e_h%omm_flavor
      call elsi_say(e_h,info_str)

      write(info_str,"('  | OMM minimization tolerance ',E10.2)") e_h%omm_tol
      call elsi_say(e_h,info_str)
   case(PEXSI_SOLVER)
      call elsi_say(e_h,"  PEXSI settings:")

      write(info_str,"('  | Electron temperature       ',E10.2)")&
         e_h%pexsi_options%temperature
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Spectral gap               ',F10.3)")&
         e_h%pexsi_options%gap
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Spectral width             ',F10.3)")&
         e_h%pexsi_options%deltaE
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Number of poles            ',I10)")&
         e_h%pexsi_options%numPole
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Number of mu points        ',I10)")&
         e_h%pexsi_options%nPoints
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Lower bound of mu          ',E10.2)")&
         e_h%pexsi_options%muMin0
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Upper bound of mu          ',E10.2)")&
         e_h%pexsi_options%muMax0
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Inertia counting tolerance ',E10.2)")&
         e_h%pexsi_options%muInertiaTolerance
      call elsi_say(e_h,info_str)

      write(info_str,"('  | MPI tasks for symbolic     ',I10)")&
         e_h%pexsi_options%npSymbFact
      call elsi_say(e_h,info_str)
   case(SIPS_SOLVER)
      call elsi_say(e_h,"  SIPs settings:")

      write(info_str,"('  | Number of ELPA steps ',I10)") e_h%sips_n_elpa
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Slice buffer         ',E10.2)") e_h%sips_buffer
      call elsi_say(e_h,info_str)

      write(info_str,"('  | First eigensolution  ',I10)") e_h%sips_first_ev
      call elsi_say(e_h,info_str)
   case(DMP_SOLVER)
      call elsi_say(e_h,"  DMP settings:")

      write(info_str,"('  | Purification method              ',I10)")&
         e_h%dmp_method
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Max number of purification steps ',I10)")&
         e_h%dmp_max_iter
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Convergence tolerance            ',E10.2)")&
         e_h%dmp_tol
      call elsi_say(e_h,info_str)
   end select

end subroutine

!>
!! This routine prints out settings for the current (user-indicated) solver.
!!
subroutine elsi_print_solver_settings(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   type(elsi_file_io_handle) :: io_h

   character(len=40), parameter :: caller = "elsi_print_solver_settings"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   select case(e_h%solver)
   case(CHESS_SOLVER)
      call elsi_print_chess_settings(e_h,io_h)
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
      call elsi_stop(" Unsupported solver.",e_h,caller)
   end select

end subroutine

!>
!! This routine prints out settings for CheSS.
!!
subroutine elsi_print_chess_settings(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   integer(kind=i4)          :: comma_json_save
   type(elsi_file_io_handle) :: io_h
   character(len=200)        :: info_str

   character(len=40), parameter :: caller = "elsi_print_chess_settings"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   comma_json_save = io_h%comma_json
   io_h%comma_json = COMMA_AFTER ! Add commas behind all records before final

   ! Header
   if(io_h%file_format == HUMAN_READ) then
      write(info_str,"(A)") "Solver Settings (CheSS)"
   else
      write(info_str,"(A)") '"solver_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"chess_erf_decay",e_h%chess_erf_decay,io_h)
   call elsi_say_setting(e_h,"chess_erf_min",e_h%chess_erf_min,io_h)
   call elsi_say_setting(e_h,"chess_erf_max",e_h%chess_erf_max,io_h)
   call elsi_say_setting(e_h,"chess_ev_ham_min",e_h%chess_ev_ham_min,io_h)
   call elsi_say_setting(e_h,"chess_ev_ham_max",e_h%chess_ev_ham_max,io_h)
   call elsi_say_setting(e_h,"chess_ev_ovlp_min",e_h%chess_ev_ovlp_min,io_h)
   call elsi_say_setting(e_h,"chess_ev_ovlp_max",e_h%chess_ev_ovlp_max,io_h)
   io_h%comma_json = NO_COMMA ! Final record in this scope
   call elsi_say_setting(e_h,"chess_beta",e_h%chess_beta,io_h)
   call truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_json_save ! Final record, restore comma_json

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
!! This routine prints out settings for DMP.
!!
subroutine elsi_print_dmp_settings(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   integer(kind=i4)          :: comma_json_save
   type(elsi_file_io_handle) :: io_h
   character(len=200)        :: info_str

   character(len=40), parameter :: caller = "elsi_print_dmp_settings"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   comma_json_save = io_h%comma_json
   io_h%comma_json = COMMA_AFTER ! Add commas behind all records before final

   ! Header
   if(io_h%file_format == HUMAN_READ) then
      write(info_str,"(A)") "Solver Settings (DMP)"
   else
      write(info_str,"(A)") '"solver_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"dmp_n_states",e_h%dmp_n_states,io_h)
   call elsi_say_setting(e_h,"dmp_method",e_h%dmp_method,io_h)
   call elsi_say_setting(e_h,"dmp_max_power",e_h%dmp_max_power,io_h)
   call elsi_say_setting(e_h,"dmp_max_iter",e_h%dmp_max_iter,io_h)
   io_h%comma_json = NO_COMMA ! Final record in this scope
   call elsi_say_setting(e_h,"dmp_tol",e_h%dmp_tol,io_h)
   call truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_json_save ! Final record, restore comma_json

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
subroutine elsi_print_elpa_settings(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   integer(kind=i4)          :: comma_json_save
   type(elsi_file_io_handle) :: io_h
   character(len=200)        :: info_str

   character(len=40), parameter :: caller = "elsi_print_elpa_settings"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   comma_json_save = io_h%comma_json
   io_h%comma_json = COMMA_AFTER ! Add commas behind all records before final

   ! Header
   if(io_h%file_format == HUMAN_READ) then
      write(info_str,"(A)") "Solver Settings (ELPA)"
   else
      write(info_str,"(A)") '"solver_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"elpa_solver",e_h%elpa_solver,io_h)
   call elsi_say_setting(e_h,"elpa_n_states",e_h%n_states,io_h)
   io_h%comma_json = NO_COMMA ! Final record in this scope
   call elsi_say_setting(e_h,"elpa_n_single",e_h%elpa_n_single,io_h)
   call truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_json_save ! Final record, restore comma_json

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
subroutine elsi_print_omm_settings(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   integer(kind=i4)          :: comma_json_save
   type(elsi_file_io_handle) :: io_h
   character(len=200)        :: info_str

   character(len=40), parameter :: caller = "elsi_print_omm_settings"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   comma_json_save = io_h%comma_json
   io_h%comma_json = COMMA_AFTER ! Add commas behind all records before final

   ! Header
   if(io_h%file_format == HUMAN_READ) then
      write(info_str,"(A)") "Solver Settings (libOMM)"
   else
      write(info_str,"(A)") '"solver_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

  ! Settings
   call append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"omm_n_states",e_h%omm_n_states,io_h)
   call elsi_say_setting(e_h,"omm_n_elpa",e_h%omm_n_elpa,io_h)
   call elsi_say_setting(e_h,"omm_flavor",e_h%omm_flavor,io_h)
   io_h%comma_json = NO_COMMA ! Final record in this scope
   call elsi_say_setting(e_h,"omm_tol",e_h%omm_tol,io_h)
   call truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_json_save ! Final record, restore comma_json

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
subroutine elsi_print_pexsi_settings(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   integer(kind=i4)          :: comma_json_save
   type(elsi_file_io_handle) :: io_h
   character(len=200)        :: info_str

   character(len=40), parameter :: caller = "elsi_print_pexsi_settings"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   comma_json_save = io_h%comma_json
   io_h%comma_json = COMMA_AFTER ! Add commas behind all records before final

   ! Header
   if(io_h%file_format == HUMAN_READ) then
      write(info_str,"(A)") "Solver Settings (PEXSI)"
   else
      write(info_str,"(A)") '"solver_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"pexsi_np_per_pole",e_h%pexsi_np_per_pole,io_h)
   call elsi_say_setting(e_h,"pexsi_np_per_point",e_h%pexsi_np_per_point,io_h)
   call elsi_say_setting(e_h,"pexsi_n_prow_pexsi",e_h%pexsi_n_prow,io_h)
   call elsi_say_setting(e_h,"pexsi_n_pcol_pexsi",e_h%pexsi_n_pcol,io_h)
   call elsi_say_setting(e_h,"pexsi_delta_e",e_h%pexsi_options%deltaE,io_h)
   call elsi_say_setting(e_h,"pexsi_gap",e_h%pexsi_options%gap,io_h)
   call elsi_say_setting(e_h,"pexsi_n_pole",e_h%pexsi_options%numPole,io_h)
   call elsi_say_setting(e_h,"pexsi_n_point",e_h%pexsi_options%nPoints,io_h)
   io_h%comma_json = NO_COMMA ! Final record in this scope
   call elsi_say_setting(e_h,"pexsi_np_symbfact",e_h%pexsi_options%npSymbFact,&
           io_h)
   call truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_json_save ! Final record, restore comma_json

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
!! This routine prints out settings for SIPs.
!!
subroutine elsi_print_sips_settings(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   integer(kind=i4)          :: comma_json_save
   type(elsi_file_io_handle) :: io_h
   character(len=200)        :: info_str

   character(len=40), parameter :: caller = "elsi_print_sips_settings"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   comma_json_save = io_h%comma_json
   io_h%comma_json = COMMA_AFTER ! Add commas behind all records before final

   ! Header
   if(io_h%file_format == HUMAN_READ) then
      write(info_str,"(A)") "Solver Settings (SIPs)"
   else
      write(info_str,"(A)") '"solver_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"sips_n_states",e_h%n_states,io_h)
   call elsi_say_setting(e_h,"sips_n_elpa",e_h%sips_n_elpa,io_h)
   call elsi_say_setting(e_h,"sips_n_slices",e_h%sips_n_slices,io_h)
   call elsi_say_setting(e_h,"sips_np_per_slice",e_h%sips_np_per_slice,io_h)
   call elsi_say_setting(e_h,"sips_buffer",e_h%sips_buffer,io_h)
   io_h%comma_json = NO_COMMA ! Final record in this scope
   call elsi_say_setting(e_h,"sips_first_ev",e_h%sips_first_ev,io_h)
   call truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_json_save ! Final record, restore comma_json

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
!! This routine prints out settings for the current (user-indicated) solver.
!!
subroutine elsi_print_matrix_format_settings(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   type(elsi_file_io_handle) :: io_h

   character(len=40), parameter :: caller = "elsi_print_matrix_format_settings"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   select case(e_h%matrix_format)
   case(BLACS_DENSE)
      call elsi_print_blacs_dense_settings(e_h,io_h)
   case(PEXSI_CSC)
      call elsi_print_pexsi_csc_settings(e_h,io_h)
   case default
      call elsi_stop(" Unsupported matrix storage format.",e_h,caller)
   end select

end subroutine

!>
!! This routine prints out settings for the BLACS_DENSE matrix storage format.
!!
subroutine elsi_print_blacs_dense_settings(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   integer(kind=i4)          :: comma_json_save
   type(elsi_file_io_handle) :: io_h
   character(len=200)        :: info_str

   character(len=40), parameter :: caller = "elsi_print_blacs_dense_settings"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   comma_json_save = io_h%comma_json
   io_h%comma_json = COMMA_AFTER ! Add commas behind all records before final

   ! Header
   if(io_h%file_format == HUMAN_READ) then
      write(info_str,"(A)") "Matrix Storage Format Settings (BLACS_DENSE)"
   else
      write(info_str,"(A)") '"matrix_format_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"blk_row",e_h%blk_row,io_h)
   call elsi_say_setting(e_h,"blk_col",e_h%blk_col,io_h)
   call elsi_say_setting(e_h,"n_prow",e_h%n_prow,io_h)
   call elsi_say_setting(e_h,"n_pcol",e_h%n_pcol,io_h)
   io_h%comma_json = NO_COMMA ! Final record in this scope
   call elsi_say_setting(e_h,"blacs_ready",e_h%blacs_ready,io_h)
   call truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_json_save ! Final record, restore comma_json

   if(io_h%file_format == JSON) then
      if(io_h%comma_json == COMMA_AFTER) then
         write(info_str,"(A,A)") '},'
      else
         write(info_str,"(A,A)") '}'
      endif

      call elsi_say(e_h,info_str,io_h)
   endif

end subroutine

!>
!! This routine prints out settings for the PEXSI_CSC matrix storage format.
!!
subroutine elsi_print_pexsi_csc_settings(e_h,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   integer(kind=i4)          :: comma_json_save
   type(elsi_file_io_handle) :: io_h
   character(len=200)        :: info_str

   character(len=40), parameter :: caller = "elsi_print_pexsi_csc_settings"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   comma_json_save = io_h%comma_json
   io_h%comma_json = COMMA_AFTER ! Add commas behind all records before final

   ! Header
   if(io_h%file_format == HUMAN_READ) then
      write(info_str,"(A)") "Matrix Storage Format Settings (PESXI_CSC)"
   else
      write(info_str,"(A)") '"matrix_format_settings": {'
   endif

   call elsi_say(e_h,info_str,io_h)

   ! Settings
   call append_string(io_h%prefix,"  ")
   call elsi_say_setting(e_h,"nnz_g",e_h%nnz_g,io_h)
   call elsi_say_setting(e_h,"zero_def",e_h%zero_def,io_h)
   io_h%comma_json = NO_COMMA ! Final record in this scope
   call elsi_say_setting(e_h,"sparsity_ready",e_h%sparsity_ready,io_h)
   call truncate_string(io_h%prefix,2)

   ! Footer (only for JSON)
   io_h%comma_json = comma_json_save ! Final record, restore comma_json

   if(io_h%file_format == JSON) then
      if(io_h%comma_json == COMMA_AFTER) then
         write(info_str,"(A,A)") '},'
      else
         write(info_str,"(A,A)") '}'
      endif

      call elsi_say(e_h,info_str,io_h)
   endif

end subroutine

!>
!! This module procedure prints out ELSI settings in a systematic fashion.
!!
subroutine elsi_say_setting_i4(e_h,label,setting,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   character(len=*),          intent(in)           :: label
   integer(kind=i4),          intent(in)           :: setting
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   character(len=28) :: label_ljust
   character(len=20) :: int_string

   type(elsi_file_io_handle) :: io_h

   character(len=40), parameter :: caller = "elsi_say_setting_i4"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   write(int_string,'(I20)') setting

   label_ljust = label ! Store the label string in fixed-length character array

   if(io_h%print_info .and. e_h%myid_all == 0) then
      if(io_h%file_format == HUMAN_READ) then
         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,I40)") io_h%prefix,label_ljust,&
               " : ",setting
         else
            write(io_h%print_unit,"(A28,A3,I40)") label_ljust," : ",setting
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
         call elsi_stop(" Unsupported output format.",e_h,caller)
      endif
   endif

end subroutine

subroutine elsi_say_setting_r8(e_h,label,setting,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   character(len=*),          intent(in)           :: label
   real(kind=r8),             intent(in)           :: setting
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   character(len=28) :: label_ljust
   character(len=20) :: real_string

   type(elsi_file_io_handle) :: io_h

   character(len=40), parameter :: caller = "elsi_say_setting_r8"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   write(real_string,'(E20.8)') setting

   label_ljust = label ! Store the label string in fixed-length character array

   if(io_h%print_info .and. e_h%myid_all == 0) then
      if(io_h%file_format == HUMAN_READ) then
         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,E40.8)") io_h%prefix,label_ljust,&
               " : ",setting
         else
            write(io_h%print_unit,"(A28,A3,E40.8)") label_ljust," : ",setting
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
         call elsi_stop(" Unsupported output format.",e_h,caller)
      endif
   endif

end subroutine

subroutine elsi_say_setting_log(e_h,label,setting,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   character(len=*),          intent(in)           :: label
   logical,                   intent(in)           :: setting
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   character(len=28) :: label_ljust
   character(len=20) :: log_string

   type(elsi_file_io_handle) :: io_h

   character(len=40), parameter :: caller = "elsi_say_setting_log"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   if(io_h%file_format == HUMAN_READ) then
      if(setting) then
         log_string = "                TRUE"
      else
         log_string = "               FALSE" ! TIME PARADOX
      endif
   elseif(io_h%file_format == JSON) then
      ! By convention, JSON strings are lower-case
      if(setting) then
         log_string = "                true"
      else
         log_string = "               false" ! TIME PARADOX
      endif
   else
      call elsi_stop(" Unsupported output format.",e_h,caller)
   endif

   label_ljust = label ! Store the label string in fixed-length character array

   if(io_h%print_info .and. e_h%myid_all == 0) then
      if(io_h%file_format == HUMAN_READ) then
         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,A40)") io_h%prefix,label_ljust,&
               " : ",log_string
         else
            write(io_h%print_unit,"(A28,A3,A40)") label_ljust," : ",log_string
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
         call elsi_stop(" Unsupported output format.",e_h,caller)
      endif
   endif

end subroutine

subroutine elsi_say_setting_str(e_h,label,setting,io_h_in)

   implicit none

   type(elsi_handle),         intent(in)           :: e_h
   character(len=*),          intent(in)           :: label
   character(len=*),          intent(in)           :: setting
   type(elsi_file_io_handle), intent(in), optional :: io_h_in

   character(len=28) :: label_ljust

   type(elsi_file_io_handle) :: io_h

   character(len=40), parameter :: caller = "elsi_say_setting_str"

   if(present(io_h_in)) then
      io_h = io_h_in
   else
      io_h = e_h%stdio
   endif

   label_ljust = label ! Store the label string in fixed-length character array

   if(io_h%print_info .and. e_h%myid_all == 0) then
      if(io_h%file_format == HUMAN_READ) then
         if(allocated(io_h%prefix)) then
            write(io_h%print_unit,"(A,A28,A3,A40)") io_h%prefix,label_ljust,&
               " : ",setting
         else
            write(io_h%print_unit,"(A28,A3,A40)") label_ljust," : ",setting
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
         call elsi_stop(" Unsupported output format.",e_h,caller)
      endif
   endif

end subroutine

!>
!! This routine generates a new (dynamic) string with another string appended to
!! the end. Whitespace is preserved deliberately.
!!
subroutine append_string(l_string,r_string)

   implicit none

   character(len=:), intent(inout), allocatable :: l_string
   character(len=*), intent(in)                 :: r_string

   character(len=:), allocatable :: temp_string

   ! Create a temporary character array holding the new character array
   if(allocated(l_string)) then
      temp_string = l_string // r_string
   else
      temp_string = r_string
   endif

   ! Now deallocate the old character array and replace with the new one
   if(allocated(l_string)) then
      deallocate(l_string)
   endif

   l_string = temp_string

   deallocate(temp_string)

end subroutine

!>
!! This routine generates a new string with the indicated number of characters
!! removed from the end.
!!
subroutine truncate_string(l_string,n_chars_to_remove)

   implicit none

   character(len=:), intent(inout), allocatable :: l_string
   integer(kind=i4), intent(in)                 :: n_chars_to_remove

   integer(kind=i4) :: size_new_string

   character(len=:), allocatable :: temp_string

   ! Find size of new character array
   if(allocated(l_string)) then
      size_new_string = len(l_string)-n_chars_to_remove
   else
      return
   endif

   if(size_new_string < 1) then
      deallocate(l_string)
      return
   endif

   ! Create a temporary character array holding the new character array
   temp_string = l_string(1:size_new_string)

   ! Now deallocate the old character array and replace with the new one
   deallocate(l_string)

   l_string = temp_string

   deallocate(temp_string)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines for streamlining the patterned output of JSON files with ELSI.   !
! These subroutines create the relevant file IO handle and write the various   !
! opening/closing brackets for the JSON arrays and records.  The actual        !
! name/value pairs contained in the records are written by elsi_say_setting,   !
! which is not included in this set of subroutines as it is also used for      !
! human-readable output.                                                       !
! These subroutines require the ELSI handle only to specify the MPI layout.    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! This routine generates a file IO handle for the JSON file and opens the file
!!
subroutine elsi_open_json_file(e_h,print_unit,file_name,opening_bracket,io_h)

   implicit none

   type(elsi_handle),         intent(in)  :: e_h
   integer(kind=i4),          intent(in)  :: print_unit
   character(len=*),          intent(in)  :: file_name
   logical,                   intent(in)  :: opening_bracket
   type(elsi_file_io_handle), intent(out) :: io_h

   character(len=40), parameter :: caller = "elsi_open_json_file"

   call elsi_init_file_io(io_h,print_unit,file_name,JSON,.true.,"",COMMA_AFTER)

   open(unit=io_h%print_unit,file=io_h%file_name)

   if(opening_bracket) then
      call elsi_say(e_h,"[",io_h)
      call append_string(io_h%prefix,"  ")
   endif

end subroutine

!>
!! This routine closes the JSON file and tears down the file IO handle
!!
subroutine elsi_close_json_file(e_h,closing_bracket,io_h)

   implicit none

   type(elsi_handle),         intent(in)    :: e_h
   logical,                   intent(in)    :: closing_bracket
   type(elsi_file_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_close_json_file"

   if(io_h%file_format /= JSON) then
      call elsi_stop(" This routine requires a file handle using JSON format.",&
              e_h,caller)
   endif

   ! Closing bracket to signify end of JSON array
   if(closing_bracket) then
      call truncate_string(io_h%prefix,2)
      call elsi_say(e_h,"]",io_h)
   endif

   close(io_h%print_unit)

   call elsi_finalize_file_io(e_h,io_h)

end subroutine

!>
!! This routine starts a new record in the JSON file
!!
subroutine elsi_start_json_record(e_h,comma_before,io_h)

   implicit none

   type(elsi_handle),         intent(in)    :: e_h
   logical,                   intent(in)    :: comma_before
   type(elsi_file_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_start_json_record"

   if(io_h%file_format /= JSON) then
      call elsi_stop(" This routine requires a file handle using JSON format.",&
              e_h,caller)
   endif

   if(comma_before) then
      call elsi_say(e_h,',{',io_h)
   else
      call elsi_say(e_h,'{',io_h)
   endif

   call append_string(io_h%prefix,"  ")

end subroutine

!>
!! This routine finishes the current record in the JSON file
!!
subroutine elsi_finish_json_record(e_h,comma_after,io_h)

   implicit none

   type(elsi_handle),         intent(in)    :: e_h
   logical,                   intent(in)    :: comma_after
   type(elsi_file_io_handle), intent(inout) :: io_h

   character(len=40), parameter :: caller = "elsi_finish_json_record"

   if(io_h%file_format /= JSON) then
      call elsi_stop(" This routine requires a file handle using JSON format.",&
              e_h,caller)
   endif

   call truncate_string(io_h%prefix,2)

   if(comma_after) then
      call elsi_say(e_h,'},',io_h)
   else
      call elsi_say(e_h,'}',io_h)
   endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                              End of JSON Code                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ELSI_IO
