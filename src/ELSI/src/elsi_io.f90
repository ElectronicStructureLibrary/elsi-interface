! Copyright (c) 2015-2017, the ELSI team. All rights reserved.
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
!! This module performs IO to stdout and files.  For matrix IO, see the ELSI_MATIO
!! module.
!!
module ELSI_IO

   use ELSI_CONSTANTS
   use ELSI_DATATYPE
   use ELSI_MPI,       only: elsi_stop
   use ELSI_PRECISION, only: r8,i4

   implicit none

   private

   ! WPH:  I'm using "IO" to describe these subroutines for now, but at the time
   !       of this writing (25 December 2017), we only support output to non-matrix
   !       files.  I expect this will change eventually.

   ! Core IO Subroutines
   public :: elsi_say
   public :: elsi_say_setting
   ! IO for ELSI handle settings
   public :: elsi_print_handle_summary
   ! IO for solver settings 
   public :: elsi_print_settings
   public :: elsi_print_solver_settings
   public :: elsi_print_chess_settings
   public :: elsi_print_dmp_settings
   public :: elsi_print_elpa_settings
   public :: elsi_print_omm_settings
   public :: elsi_print_pexsi_settings
   public :: elsi_print_sips_settings
   ! IO for matrix storage format settings
   public :: elsi_print_matrix_format_settings
   public :: elsi_print_blacs_dense_settings
   public :: elsi_print_pexsi_csc_settings

   interface elsi_say_setting
      module procedure elsi_say_setting_i4,&
                       elsi_say_setting_r8,&
                       elsi_say_setting_log,&
                       elsi_say_setting_str
   end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                           CORE IO SUBROUTINES                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! This routine prints a message.
!!
subroutine elsi_say(e_h,info_str,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: info_str !< Message to print
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to print to

   integer(kind=i4) :: my_unit

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   if(e_h%print_info) then
      if(e_h%myid_all == 0) then
         write(my_unit,"(A)") trim(info_str)
      endif
   endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         IO FOR HANDLE SETTINGS                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! This routine prints the state of the handle
!!
subroutine elsi_print_handle_summary(e_h,prefix,use_unit)

   implicit none

   type(elsi_handle),          intent(in) :: e_h      !< Handle
   character(len=*),           intent(in) :: prefix   !< Prefix for every line
   integer(kind=i4), optional, intent(in) :: use_unit

   real(kind=r8)    :: sparsity
   character*200    :: info_str
   integer(kind=i4) :: my_unit

   character*40, parameter :: caller = "elsi_print_handle_summary"

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   write(info_str,"(A,A)") prefix,          "Physical Properties"
   call elsi_say(e_h,info_str,my_unit)

   call elsi_say_setting(e_h,prefix,        "  Number of electrons",e_h%n_electrons,my_unit)
   if(e_h%parallel_mode == MULTI_PROC) then
      call elsi_say_setting(e_h,prefix,     "  Number of spins",e_h%n_spins,my_unit)
      call elsi_say_setting(e_h,prefix,     "  Number of k-points",e_h%n_kpts,my_unit)
   endif
   if(e_h%solver == ELPA_SOLVER .or. e_h%solver == SIPS_SOLVER) then
      call elsi_say_setting(e_h,prefix,     "  Number of states",e_h%n_states,my_unit)
   endif

   write(info_str,"(A,A)") prefix,          ""
   call elsi_say(e_h,info_str,my_unit)
   write(info_str,"(A,A)") prefix,          "Matrix Properties"
   call elsi_say(e_h,info_str,my_unit)

   if(e_h%matrix_format == BLACS_DENSE) then
      call elsi_say_setting(e_h,prefix,     "  Matrix format","BLACS_DENSE",my_unit)
   elseif(e_h%matrix_format == PEXSI_CSC) then
      call elsi_say_setting(e_h,prefix,     "  Matrix format","PEXSI_CSC",my_unit)
   endif
   call elsi_say_setting(e_h,prefix,        "  Number of basis functions",e_h%n_basis,my_unit)
   if(e_h%parallel_mode == MULTI_PROC) then
      sparsity = 1.0_r8-(1.0_r8*e_h%nnz_g/e_h%n_basis/e_h%n_basis)
      call elsi_say_setting(e_h,prefix,     "  Matrix sparsity",sparsity,my_unit)
   endif

   write(info_str,"(A,A)") prefix,          ""
   call elsi_say(e_h,info_str,my_unit)
   write(info_str,"(A,A)") prefix,          "Computational Details"
   call elsi_say(e_h,info_str,my_unit)

   if(e_h%parallel_mode == MULTI_PROC) then
      call elsi_say_setting(e_h,prefix,     "  Parallel mode","MULTI_PROC",my_unit)
   elseif(e_h%parallel_mode == SINGLE_PROC) then
      call elsi_say_setting(e_h,prefix,     "  Parallel mode","SINGLE_PROC",my_unit)
   endif
   call elsi_say_setting(e_h,prefix,        "  Number of MPI tasks",e_h%n_procs,my_unit)
   if(e_h%solver == ELPA_SOLVER) then
      call elsi_say_setting(e_h,prefix,     "  Solver requested","ELPA",my_unit)
   elseif(e_h%solver == OMM_SOLVER) then
      call elsi_say_setting(e_h,prefix,     "  Solver requested","libOMM",my_unit)
   elseif(e_h%solver == PEXSI_SOLVER) then
      call elsi_say_setting(e_h,prefix,     "  Solver requested","PEXSI",my_unit)
   elseif(e_h%solver == CHESS_SOLVER) then
      call elsi_say_setting(e_h,prefix,     "  Solver requested","CheSS",my_unit)
   elseif(e_h%solver == SIPS_SOLVER) then
      call elsi_say_setting(e_h,prefix,     "  Solver requested","SIPs",my_unit)
   elseif(e_h%solver == DMP_SOLVER) then
      call elsi_say_setting(e_h,prefix,     "  Solver requested","DMP",my_unit)
   endif

end subroutine

!>
!! This routine prints ELSI settings.
!! NOTE:  This subroutine is functionally identical to elsi_print_solver_settings
!!        and should be deprecated in favor of it.
!!
subroutine elsi_print_settings(e_h)

   implicit none

   type(elsi_handle), intent(in) :: e_h !< Handle

   character*200 :: info_str

   character*40, parameter :: caller = "elsi_print_settings"

   select case(e_h%solver)
   case(CHESS_SOLVER)
      call elsi_say(e_h,"  CheSS settings:")

      write(info_str,"('  | Error function decay length ',E10.2)") e_h%erf_decay
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Lower bound of decay length ',E10.2)")&
         e_h%erf_decay_min
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Upper bound of decay length ',E10.2)")&
         e_h%erf_decay_max
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Lower bound of H eigenvalue ',E10.2)")&
         e_h%ev_ham_min
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Upper bound of H eigenvalue ',E10.2)")&
         e_h%ev_ham_max
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Lower bound of S eigenvalue ',E10.2)")&
         e_h%ev_ovlp_min
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Upper bound of S eigenvalue ',E10.2)") &
         e_h%ev_ovlp_max
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

      write(info_str,"('  | OMM minimization tolerance ',E10.2)") e_h%min_tol
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

      write(info_str,"('  | Slicing method            ',I10)")&
         e_h%slicing_method
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Lower bound of eigenvalue ',E10.2)") e_h%ev_min
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Upper bound of eigenvalue ',E10.2)") e_h%ev_max
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Inertia counting          ',I10)")&
         e_h%inertia_option
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Left bound option         ',I10)") e_h%unbound
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Slice buffer              ',E10.2)")&
         e_h%slice_buffer
      call elsi_say(e_h,info_str)
   case(DMP_SOLVER)
      call elsi_say(e_h,"  DMP settings:")

      write(info_str,"('  | Purification method              ',I10)")&
         e_h%dmp_method
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Max number of purification steps ',I10)")&
         e_h%max_dmp_iter
      call elsi_say(e_h,info_str)

      write(info_str,"('  | Convergence tolerance            ',E10.2)")&
         e_h%dmp_tol
      call elsi_say(e_h,info_str)
   end select

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         IO FOR SOLVER SETTINGS                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! This routine prints out settings for the current (user-indicated) solver.
!!
subroutine elsi_print_solver_settings(e_h,prefix,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix for every line
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to write to

   integer(kind=i4) :: my_unit

   character*40, parameter :: caller = "elsi_print_solver_settings"

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   select case(e_h%solver)
   case(CHESS_SOLVER)
      call elsi_print_chess_settings(e_h,prefix,my_unit)
   case(DMP_SOLVER)
      call elsi_print_dmp_settings(e_h,prefix,my_unit)
   case(ELPA_SOLVER)
      call elsi_print_elpa_settings(e_h,prefix,my_unit)
   case(OMM_SOLVER)
      call elsi_print_omm_settings(e_h,prefix,my_unit)
   case(PEXSI_SOLVER)
      call elsi_print_pexsi_settings(e_h,prefix,my_unit)
   case(SIPS_SOLVER)
      call elsi_print_sips_settings(e_h,prefix,my_unit)
   case default
      call elsi_stop(" Unsupported solver.",e_h,caller)
   end select

end subroutine


!>
!! This routine prints out settings for CheSS.
!!
subroutine elsi_print_chess_settings(e_h,prefix,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix for every line
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to write to

   integer(kind=i4) :: my_unit
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_print_chess_settings"

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   write(info_str,"(A,A)")   prefix,"Solver Settings (CheSS)"
   call elsi_say(e_h,info_str,my_unit)
   call elsi_say_setting(e_h,prefix,"  erf_decay",e_h%erf_decay,my_unit)
   call elsi_say_setting(e_h,prefix,"  erf_decay_min",e_h%erf_decay_min,my_unit)
   call elsi_say_setting(e_h,prefix,"  erf_decay_max",e_h%erf_decay_max,my_unit)
   call elsi_say_setting(e_h,prefix,"  ev_ham_min",e_h%ev_ham_min,my_unit)
   call elsi_say_setting(e_h,prefix,"  ev_ham_max",e_h%ev_ham_max,my_unit)
   call elsi_say_setting(e_h,prefix,"  ev_ovlp_min",e_h%ev_ovlp_min,my_unit)
   call elsi_say_setting(e_h,prefix,"  ev_ovlp_max",e_h%ev_ovlp_max,my_unit)
   call elsi_say_setting(e_h,prefix,"  beta",e_h%beta,my_unit)
   call elsi_say_setting(e_h,prefix,"  chess_started",e_h%chess_started,my_unit)

end subroutine

!>
!! This routine prints out settings for DMP.
!!
subroutine elsi_print_dmp_settings(e_h,prefix,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix for every line
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to write to

   integer(kind=i4) :: my_unit
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_print_dmp_settings"

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   write(info_str,"(A,A)")   prefix,"Solver Settings (DMP)"
   call elsi_say(e_h,info_str,my_unit)
   call elsi_say_setting(e_h,prefix,"  n_states_dmp",e_h%n_states_dmp,my_unit)
   call elsi_say_setting(e_h,prefix,"  dmp_method",e_h%dmp_method,my_unit)
   call elsi_say_setting(e_h,prefix,"  max_power_iter",e_h%max_power_iter,my_unit)
   call elsi_say_setting(e_h,prefix,"  max_dmp_iter",e_h%max_dmp_iter,my_unit)
   call elsi_say_setting(e_h,prefix,"  dmp_tol",e_h%dmp_tol,my_unit)
   call elsi_say_setting(e_h,prefix,"  ne_dmp",e_h%ne_dmp,my_unit)

end subroutine

!>
!! This routine prints out settings for ELPA.
!!
subroutine elsi_print_elpa_settings(e_h,prefix,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix for every line
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to write to

   integer(kind=i4) :: my_unit
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_print_elpa_settings"

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   write(info_str,"(A,A)")   prefix,"Solver Settings (ELPA)"
   call elsi_say(e_h,info_str,my_unit)
   call elsi_say_setting(e_h,prefix,"  elpa_solver",e_h%elpa_solver,my_unit)
   call elsi_say_setting(e_h,prefix,"  n_states",e_h%n_states,my_unit)
   call elsi_say_setting(e_h,prefix,"  n_single_steps",e_h%n_single_steps,my_unit)
   call elsi_say_setting(e_h,prefix,"  elpa_output",e_h%elpa_output,my_unit)
   call elsi_say_setting(e_h,prefix,"  elpa_started",e_h%elpa_started,my_unit)

end subroutine

!>
!! This routine prints out settings for libOMM.
!!
subroutine elsi_print_omm_settings(e_h,prefix,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix for every line
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to write to

   integer(kind=i4) :: my_unit
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_print_omm_settings"

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   write(info_str,"(A,A)")   prefix,"Solver Settings (libOMM)"
   call elsi_say(e_h,info_str,my_unit)
   call elsi_say_setting(e_h,prefix,"  n_states_omm",e_h%n_states_omm,my_unit)
   call elsi_say_setting(e_h,prefix,"  omm_n_elpa",e_h%omm_n_elpa,my_unit)
   call elsi_say_setting(e_h,prefix,"  new_ovlp",e_h%new_ovlp,my_unit)
   call elsi_say_setting(e_h,prefix,"  coeff_ready",e_h%coeff_ready,my_unit)
   call elsi_say_setting(e_h,prefix,"  omm_flavor",e_h%omm_flavor,my_unit)
   call elsi_say_setting(e_h,prefix,"  scale_kinetic",e_h%scale_kinetic,my_unit)
   call elsi_say_setting(e_h,prefix,"  calc_ed",e_h%calc_ed,my_unit)
   call elsi_say_setting(e_h,prefix,"  eta",e_h%eta,my_unit)
   call elsi_say_setting(e_h,prefix,"  min_tol",e_h%min_tol,my_unit)
   call elsi_say_setting(e_h,prefix,"  omm_output",e_h%omm_output,my_unit)
   call elsi_say_setting(e_h,prefix,"  do_dealloc",e_h%do_dealloc,my_unit)
   call elsi_say_setting(e_h,prefix,"  use_psp",e_h%use_psp,my_unit)

end subroutine

!>
!! This routine prints out settings for PEXSI.
!!
subroutine elsi_print_pexsi_settings(e_h,prefix,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix for every line
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to write to

   integer(kind=i4) :: my_unit
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_print_pexsi_settings"

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   write(info_str,"(A,A)")   prefix,"Solver Settings (PEXSI)"
   call elsi_say(e_h,info_str,my_unit)
   call elsi_say_setting(e_h,prefix,"  np_per_pole",e_h%np_per_pole,my_unit)
   call elsi_say_setting(e_h,prefix,"  np_per_point",e_h%np_per_point,my_unit)
   call elsi_say_setting(e_h,prefix,"  n_prow_pexsi",e_h%n_prow_pexsi,my_unit)
   call elsi_say_setting(e_h,prefix,"  n_pcol_pexsi",e_h%n_pcol_pexsi,my_unit)
   call elsi_say_setting(e_h,prefix,"  ne_pexsi",e_h%ne_pexsi,my_unit)
   call elsi_say_setting(e_h,prefix,"  pexsi_started",e_h%pexsi_started,my_unit)

   ! The following are specific to PEXSI parallelization and vary from task to task.
   ! I'm leaving them in, but commented out, should we ever need debug-level output.
   ! This would required patterned IO be written.
   !call elsi_say_setting(e_h,prefix,"my_prow_pexsi",e_h%my_prow_pexsi,my_unit)
   !call elsi_say_setting(e_h,prefix,"my_pcol_pexsi",e_h%my_pcol_pexsi,my_unit)
   !call elsi_say_setting(e_h,prefix,"my_point",e_h%my_point,my_unit)
   !call elsi_say_setting(e_h,prefix,"myid_point",e_h%myid_point,my_unit)
   !call elsi_say_setting(e_h,prefix,"comm_among_pole",e_h%comm_among_pole,my_unit)
   !call elsi_say_setting(e_h,prefix,"comm_in_pole",e_h%comm_in_pole,my_unit)
   !call elsi_say_setting(e_h,prefix,"comm_among_point",e_h%comm_among_point,my_unit)
   !call elsi_say_setting(e_h,prefix,"comm_in_point",e_h%comm_in_point,my_unit)
end subroutine

!>
!! This routine prints out settings for SIPs.
!!
subroutine elsi_print_sips_settings(e_h,prefix,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix for every line
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to write to

   integer(kind=i4) :: my_unit
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_print_sips_settings"

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   write(info_str,"(A,A)")   prefix,"Solver Settings (SIPs)"
   call elsi_say(e_h,info_str,my_unit)
   call elsi_say_setting(e_h,prefix,"  n_states",e_h%n_states,my_unit)
   call elsi_say_setting(e_h,prefix,"  sips_n_elpa",e_h%sips_n_elpa,my_unit)
   call elsi_say_setting(e_h,prefix,"  np_per_slice",e_h%np_per_slice,my_unit)
   call elsi_say_setting(e_h,prefix,"  n_inertia_steps",e_h%n_inertia_steps,my_unit)
   call elsi_say_setting(e_h,prefix,"  slicing_method",e_h%slicing_method,my_unit)
   call elsi_say_setting(e_h,prefix,"  inertia_option",e_h%inertia_option,my_unit)
   call elsi_say_setting(e_h,prefix,"  unbound",e_h%unbound,my_unit)
   call elsi_say_setting(e_h,prefix,"  n_slices",e_h%n_slices,my_unit)
   call elsi_say_setting(e_h,prefix,"  interval(1)",e_h%interval(1),my_unit)
   call elsi_say_setting(e_h,prefix,"  interval(2)",e_h%interval(2),my_unit)
   call elsi_say_setting(e_h,prefix,"  slice_buffer",e_h%slice_buffer,my_unit)
   call elsi_say_setting(e_h,prefix,"  ev_min",e_h%ev_min,my_unit)
   call elsi_say_setting(e_h,prefix,"  ev_max",e_h%ev_max,my_unit)
   call elsi_say_setting(e_h,prefix,"  sips_started",e_h%sips_started,my_unit)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  IO FOR MATRIX STORGE FORMAT SETTINGS                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! This routine prints out settings for the current (user-indicated) solver.
!!
subroutine elsi_print_matrix_format_settings(e_h,prefix,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix for every line
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to write to

   integer(kind=i4) :: my_unit

   character*40, parameter :: caller = "elsi_print_matrix_format_settings"

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   select case(e_h%matrix_format)
   case(BLACS_DENSE)
      call elsi_print_blacs_dense_settings(e_h,prefix,my_unit)
   case(PEXSI_CSC)
      call elsi_print_pexsi_csc_settings(e_h,prefix,my_unit)
   case default
      call elsi_stop(" Unsupported matrix storage format.",e_h,caller)
   end select

end subroutine

!>
!! This routine prints out settings for the BLACS_DENSE matrix storage format.
!!
subroutine elsi_print_blacs_dense_settings(e_h,prefix,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix for every line
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to write to

   integer(kind=i4) :: my_unit
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_print_blacs_dense_settings"

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   write(info_str,"(A,A)")   prefix,"Matrix Storage Format Settings (BLACS_DENSE)"
   call elsi_say(e_h,info_str,my_unit)
   call elsi_say_setting(e_h,prefix,"  blk_row",e_h%blk_row,my_unit)
   call elsi_say_setting(e_h,prefix,"  blk_col",e_h%blk_col,my_unit)
   call elsi_say_setting(e_h,prefix,"  n_prow",e_h%n_prow,my_unit)
   call elsi_say_setting(e_h,prefix,"  n_pcol",e_h%n_pcol,my_unit)
   call elsi_say_setting(e_h,prefix,"  blacs_ready",e_h%blacs_ready,my_unit)

end subroutine

!>
!! This routine prints out settings for the PEXSI_CSC matrix storage format.
!!
subroutine elsi_print_pexsi_csc_settings(e_h,prefix,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix for every line
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to write to

   integer(kind=i4) :: my_unit
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_print_pexsi_csc_settings"

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   write(info_str,"(A,A)")   prefix,"Matrix Storage Format Settings (PESXI_CSC)"
   call elsi_say(e_h,info_str,my_unit)
   call elsi_say_setting(e_h,prefix,"  nnz_g",e_h%nnz_g,my_unit)
   call elsi_say_setting(e_h,prefix,"  zero_def",e_h%zero_def,my_unit)
   call elsi_say_setting(e_h,prefix,"  sparsity_ready",e_h%sparsity_ready,my_unit)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                            ELSI_SAY_SETTINGS                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! This module procedure is used to print out ELSI settings in a systematic fashion.
!! TODO:  Generate formatting strings on-the-fly so that we can rid of hard-coded
!!        constants outside of ELSI_CONSTANTS
!!

subroutine elsi_say_setting_i4(e_h,prefix,label,setting,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix before label
   character(len=*),            intent(in) :: label    !< Label for setting to print
   integer(kind=i4),            intent(in) :: setting  !< Value for setting to print
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to print to

   integer(kind=i4)  :: my_unit
   character(len=27) :: label_ljust

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   label_ljust = label ! Store the label string in a fixed-length character array 
                       ! so that it is right-justified when output.

   if(e_h%print_info) then
      if(e_h%myid_all == 0) then
         write(my_unit,"(A,A27,A3,I20)") prefix, label_ljust, " : ", setting
      endif
   endif

end subroutine

subroutine elsi_say_setting_r8(e_h,prefix,label,setting,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix before label
   character(len=*),            intent(in) :: label    !< Label for setting to print
   real(kind=r8),               intent(in) :: setting  !< Value for setting to print
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to print to

   integer(kind=i4)  :: my_unit
   character(len=27) :: label_ljust

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   label_ljust = label ! Store the label string in a fixed-length character array 
                       ! so that it is right-justified when output.

   if(e_h%print_info) then
      if(e_h%myid_all == 0) then
         write(my_unit,"(A,A27,A3,E20.8)") prefix, label_ljust, " : ", setting
      endif
   endif

end subroutine

subroutine elsi_say_setting_log(e_h,prefix,label,setting,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix before label
   character(len=*),            intent(in) :: label    !< Label for setting to print
   logical,                     intent(in) :: setting  !< Value for setting to print
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to print to

   integer(kind=i4)  :: my_unit
   character(len=27) :: label_ljust
   character(len=20) :: truth

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   label_ljust = label ! Store the label string in a fixed-length character array 
                       ! so that it is right-justified when output.

   if(setting) then
      truth = "                TRUE"
   else
      truth = "               FALSE" ! TIME PARADOX
   endif

   if(e_h%print_info) then
      if(e_h%myid_all == 0) then
         write(my_unit,"(A,A27,A3,A20)") prefix, label_ljust, " : ", truth
      endif
   endif

end subroutine

subroutine elsi_say_setting_str(e_h,prefix,label,setting,use_unit)

   implicit none

   type(elsi_handle),           intent(in) :: e_h      !< Handle
   character(len=*),            intent(in) :: prefix   !< Prefix before label
   character(len=*),            intent(in) :: label    !< Label for setting to print
   character(len=*),            intent(in) :: setting  !< Value for setting to print
   integer(kind=i4),  optional, intent(in) :: use_unit !< Unit to print to

   integer(kind=i4)  :: my_unit
   character(len=27) :: label_ljust

   if(present(use_unit)) then
      my_unit = use_unit
   else
      my_unit = e_h%print_unit
   endif

   label_ljust = label ! Store the label string in a fixed-length character array 
                       ! so that it is left-justified when output.

   if(e_h%print_info) then
      if(e_h%myid_all == 0) then
         write(my_unit,"(A,A27,A3,A20)") prefix, label_ljust, " : ", setting
      endif
   endif

end subroutine

end module ELSI_IO
