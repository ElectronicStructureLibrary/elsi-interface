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
!! This module contains a collection of basic utility routines.
!!
module ELSI_UTILS

   use ELSI_CONSTANTS
   use ELSI_DATATYPE
   use ELSI_PRECISION, only: i4,r8

   implicit none

   include "mpif.h"

   public :: elsi_say
   public :: elsi_stop
   public :: elsi_check
   public :: elsi_check_handle
   public :: elsi_get_global_row
   public :: elsi_get_global_col
   public :: elsi_get_local_nnz
   public :: elsi_trace_mat
   public :: elsi_trace_mat_mat
   public :: elsi_init_timer
   public :: elsi_get_time

   interface elsi_get_local_nnz
      module procedure elsi_get_local_nnz_real,&
                       elsi_get_local_nnz_complex
   end interface

   interface elsi_trace_mat
      module procedure elsi_trace_mat_real,&
                       elsi_trace_mat_complex
   end interface

   interface elsi_trace_mat_mat
      module procedure elsi_trace_mat_mat_real,&
                       elsi_trace_mat_mat_complex
   end interface

contains

!>
!! This routine prints a message.
!!
subroutine elsi_say(e_h,info_str)

   implicit none

   type(elsi_handle), intent(in) :: e_h      !< Handle
   character(len=*),  intent(in) :: info_str !< Message to print

   if(e_h%print_info) then
      if(e_h%myid_all == 0) then
         write(e_h%print_unit,"(A)") trim(info_str)
      endif
   endif

end subroutine

!>
!! Clean shutdown in case of errors.
!!
subroutine elsi_stop(info,e_h,caller)

   implicit none

   character(len=*),  intent(in) :: info   !< Error message
   type(elsi_handle), intent(in) :: e_h    !< Handle
   character(len=*),  intent(in) :: caller !< Caller

   character*800    :: info_str
   integer(kind=i4) :: mpierr

   if(e_h%global_mpi_ready) then
      write(info_str,"(A,I7,5A)") "**Error! MPI task ",e_h%myid_all," in ",&
         trim(caller),": ",trim(info)," Exiting..."
      write(e_h%print_unit,"(A)") trim(info_str)

      if(e_h%n_procs_all > 1) then
         call MPI_Abort(e_h%mpi_comm_all,0,mpierr)
      endif
   elseif(e_h%mpi_ready) then
      write(info_str,"(A,I7,5A)") "**Error! MPI task ",e_h%myid," in ",&
         trim(caller),": ",trim(info)," Exiting..."
      write(e_h%print_unit,"(A)") trim(info_str)

      if(e_h%n_procs > 1) then
         call MPI_Abort(e_h%mpi_comm,0,mpierr)
      endif
   else
      write(info_str,"(5A)") "**Error! ",trim(caller),": ",trim(info),&
         " Exiting..."
      write(e_h%print_unit,"(A)") trim(info_str)
   endif

   stop

end subroutine

!>
!! This routine resets an ELSI handle.
!!
subroutine elsi_reset_handle(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*40, parameter :: caller = "elsi_reset_handle"

   e_h%handle_ready     = .false.
   e_h%solver           = UNSET
   e_h%data_type        = UNSET
   e_h%matrix_format    = UNSET
   e_h%uplo             = FULL_MAT
   e_h%parallel_mode    = UNSET
   e_h%print_info       = .false.
   e_h%print_mem        = .false.
   e_h%print_unit       = 6
   e_h%n_elsi_calls     = 0
   e_h%myid             = UNSET
   e_h%myid_all         = UNSET
   e_h%n_procs          = UNSET
   e_h%n_procs_all      = UNSET
   e_h%mpi_comm         = UNSET
   e_h%mpi_comm_all     = UNSET
   e_h%mpi_comm_row     = UNSET
   e_h%mpi_comm_col     = UNSET
   e_h%mpi_ready        = .false.
   e_h%global_mpi_ready = .false.
   e_h%blacs_ctxt       = UNSET
   e_h%sc_desc          = UNSET
   e_h%blk_row          = UNSET
   e_h%blk_col          = UNSET
   e_h%n_prow           = UNSET
   e_h%n_pcol           = UNSET
   e_h%my_prow          = UNSET
   e_h%my_pcol          = UNSET
   e_h%n_lrow           = UNSET
   e_h%n_lcol           = UNSET
   e_h%blacs_ready      = .false.
   e_h%nnz_g            = UNSET
   e_h%nnz_l            = UNSET
   e_h%nnz_l_sp         = UNSET
   e_h%n_lcol_sp        = UNSET
   e_h%zero_def         = 1.0e-15_r8
   e_h%sparsity_ready   = .false.
   e_h%ovlp_is_unit     = .false.
   e_h%ovlp_is_sing     = .false.
   e_h%no_sing_check    = .false.
   e_h%sing_tol         = 1.0e-5_r8
   e_h%stop_sing        = .false.
   e_h%n_nonsing        = UNSET
   e_h%n_electrons      = 0.0_r8
   e_h%mu               = 0.0_r8
   e_h%n_basis          = UNSET
   e_h%n_spins          = 1
   e_h%n_kpts           = 1
   e_h%n_states         = UNSET
   e_h%n_states_solve   = UNSET
   e_h%i_spin           = 1
   e_h%i_kpt            = 1
   e_h%i_weight         = 1.0_r8
   e_h%energy_hdm       = 0.0_r8
   e_h%energy_sedm      = 0.0_r8
   e_h%broaden_scheme   = 0
   e_h%broaden_width    = 1.0e-2_r8
   e_h%occ_tolerance    = 1.0e-13_r8
   e_h%max_mu_steps     = 100
   e_h%spin_degen       = 0.0_r8
   e_h%spin_is_set      = .false.
   e_h%mu_ready         = .false.
   e_h%edm_ready_real   = .false.
   e_h%edm_ready_cmplx  = .false.
   e_h%elpa_solver      = UNSET
   e_h%n_single_steps   = UNSET
   e_h%elpa_output      = .false.
   e_h%elpa_started     = .false.
   e_h%n_states_omm     = UNSET
   e_h%omm_n_elpa       = UNSET
   e_h%new_overlap      = .true.
   e_h%coeff_ready      = .false.
   e_h%omm_flavor       = UNSET
   e_h%scale_kinetic    = 0.0_r8
   e_h%calc_ed          = .false.
   e_h%eta              = 0.0_r8
   e_h%min_tol          = 1.0e-12_r8
   e_h%omm_output       = .false.
   e_h%do_dealloc       = .false.
   e_h%use_psp          = .false.
   e_h%np_per_pole      = UNSET
   e_h%np_per_point     = UNSET
   e_h%my_prow_pexsi    = UNSET
   e_h%my_pcol_pexsi    = UNSET
   e_h%n_prow_pexsi     = UNSET
   e_h%n_pcol_pexsi     = UNSET
   e_h%my_point         = UNSET
   e_h%myid_point       = UNSET
   e_h%comm_among_pole  = UNSET
   e_h%comm_in_pole     = UNSET
   e_h%comm_among_point = UNSET
   e_h%comm_in_point    = UNSET
   e_h%ne_pexsi         = 0.0_r8
   e_h%pexsi_started    = .false.
   e_h%erf_decay        = 0.0_r8
   e_h%erf_decay_min    = 0.0_r8
   e_h%erf_decay_max    = 0.0_r8
   e_h%ev_ham_min       = 0.0_r8
   e_h%ev_ham_max       = 0.0_r8
   e_h%ev_ovlp_min      = 0.0_r8
   e_h%ev_ovlp_max      = 0.0_r8
   e_h%beta             = 0.0_r8
   e_h%chess_started    = .false.
   e_h%sips_n_elpa      = UNSET
   e_h%np_per_slice     = UNSET
   e_h%n_inertia_steps  = UNSET
   e_h%slicing_method   = UNSET
   e_h%inertia_option   = UNSET
   e_h%unbound          = UNSET
   e_h%n_slices         = UNSET
   e_h%interval         = 0.0_r8
   e_h%slice_buffer     = 0.0_r8
   e_h%ev_min           = 0.0_r8
   e_h%ev_max           = 0.0_r8
   e_h%sips_started     = .false.
   e_h%clock_rate       = UNSET

end subroutine

!>
!! This routine guarantees that there are no mutually conflicting parameters.
!!
subroutine elsi_check(e_h,caller)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   character(len=*),  intent(in)    :: caller !< Caller

   ! General check of solver, parallel mode, data type, matrix format
   if(e_h%solver < 0 .or. e_h%solver >= N_SOLVERS) then
      call elsi_stop(" Unsupported solver.",e_h,caller)
   endif

   if(e_h%parallel_mode < 0 .or. e_h%parallel_mode >= N_PARALLEL_MODES) then
      call elsi_stop(" Unsupported parallel mode.",e_h,caller)
   endif

   if(e_h%data_type < 0 .or. e_h%data_type >= 2) then
      call elsi_stop(" Unsupported matirx data type.",e_h,caller)
   endif

   if(e_h%matrix_format < 0 .or. e_h%matrix_format >= N_MATRIX_FORMATS) then
      call elsi_stop(" Unsupported matirx format.",e_h,caller)
   endif

   if(e_h%uplo /= FULL_MAT) then
      if(e_h%matrix_format /= BLACS_DENSE .or. &
         e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(" Triangular matrix input only supported with BLACS_"//&
                 "DENSE matrix format and MULTI_PROC parallel mode.",e_h,caller)
      endif

      if(e_h%uplo /= UT_MAT .and. e_h%uplo /= LT_MAT) then
         call elsi_stop(" Invalid choice of uplo.",e_h,caller)
      endif
   endif

   ! Spin and k-point
   if(e_h%n_spins*e_h%n_kpts > 1) then
      if(.not. e_h%global_mpi_ready) then
         call elsi_stop(" Spin/k-point calculations require a global MPI"//&
                 " communicator.",e_h,caller)
      endif
   endif

   if(.not. e_h%global_mpi_ready) then
      e_h%mpi_comm_all = e_h%mpi_comm
      e_h%n_procs_all  = e_h%n_procs
      e_h%myid_all     = e_h%myid
   endif

   if(e_h%parallel_mode == MULTI_PROC) then
      if(.not. e_h%mpi_ready) then
         call elsi_stop(" MULTI_PROC parallel mode requires MPI.",e_h,caller)
      endif
   endif

   if(e_h%matrix_format == BLACS_DENSE) then
      if(.not. e_h%blacs_ready .and. e_h%parallel_mode /= SINGLE_PROC) then
         call elsi_stop(" BLACS_DENSE matrix format requires BLACS.",e_h,caller)
      endif
   else
      if(.not. e_h%sparsity_ready) then
         call elsi_stop(" PEXSI_CSC matrix format requires a sparsity"//&
                 " pattern.",e_h,caller)
      endif
   endif

   ! Specific check for each solver
   select case(e_h%solver)
   case(AUTO)
      call elsi_stop(" Solver auto-selection not yet available.",e_h,caller)
   case(ELPA_SOLVER)
      ! Nothing
   case(OMM_SOLVER)
      if(e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(" libOMM solver requires MULTI_PROC parallel mode.",&
                 e_h,caller)
      endif
   case(PEXSI_SOLVER)
      if(e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(" PEXSI solver requires MULTI_PROC parallel mode.",e_h,&
                 caller)
      endif

      if(e_h%n_basis < e_h%np_per_pole) then
         call elsi_stop(" For this number of MPI tasks, the matrix size is"//&
                 " too small to use PEXSI.",e_h,caller)
      endif

      if(e_h%np_per_pole == UNSET) then
         if(mod(e_h%n_procs,e_h%pexsi_options%numPole*&
            e_h%pexsi_options%nPoints) /= 0) then
            call elsi_stop(" To use PEXSI, the total number of MPI tasks"//&
                    " must be a multiple of the number of MPI tasks per pole"//&
                    " times the number of mu points.",e_h,caller)
         endif
      else
         if(mod(e_h%n_procs,e_h%np_per_pole*&
            e_h%pexsi_options%nPoints) /= 0) then
            call elsi_stop(" To use PEXSI, the total number of MPI tasks"//&
                    " must be a multiple of the number of MPI tasks per pole"//&
                    " times the number of mu points.",e_h,caller)
         endif

         if(e_h%np_per_pole*e_h%pexsi_options%numPole*&
            e_h%pexsi_options%nPoints < e_h%n_procs) then
            call elsi_stop(" Specified number of MPI tasks per pole is too"//&
                    " small for the total number of MPI tasks.",e_h,caller)
         endif
      endif
   case(CHESS_SOLVER)
      call elsi_say(e_h,"  ATTENTION! CheSS is EXPERIMENTAL.")

      if(e_h%n_basis < e_h%n_procs) then
         call elsi_stop(" For this number of MPI tasks, the matrix size is"//&
                 " too small to use CheSS.",e_h,caller)
      endif

      if(e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(" CheSS solver requires MULTI_PROC parallel mode.",e_h,&
                 caller)
      endif

      if(e_h%ovlp_is_unit) then
         call elsi_stop(" CheSS solver with an identity overlap matrix not"//&
                 " yet available.",e_h,caller)
      endif
   case(SIPS_SOLVER)
      call elsi_say(e_h,"  ATTENTION! SIPs is EXPERIMENTAL.")

      if(e_h%n_basis < e_h%n_procs) then
         call elsi_stop(" For this number of MPI tasks, the matrix size is"//&
                 " too small to use SIPs.",e_h,caller)
      endif

      if(e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(" SIPs solver requires MULTI_PROC parallel mode.",e_h,&
                 caller)
      endif
   case(DMP_SOLVER)
      call elsi_say(e_h,"  ATTENTION! DMP is EXPERIMENTAL.")

      if(e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(" DMP solver requires MULTI_PROC parallel mode.",e_h,&
                 caller)
      endif
   case default
      call elsi_stop(" Unsupported solver.",e_h,caller)
   end select

end subroutine

!>
!! This routine checks whether a handle has been properly initialized. This is
!! called at the very beginning of all public-facing subroutines except the one
!! that initializes a handle, elsi_init.
!!
subroutine elsi_check_handle(e_h,caller)

   implicit none

   type(elsi_handle), intent(in) :: e_h    !< Handle
   character(len=*),  intent(in) :: caller !< Caller

   if(.not. e_h%handle_ready) then
      call elsi_stop(" Invalid handle! Not initialized.",e_h,caller)
   endif

end subroutine

!>
!! This routine computes the global row index based on the local row index.
!!
subroutine elsi_get_global_row(e_h,g_id,l_id)

   implicit none

   type(elsi_handle), intent(in)  :: e_h  !< Handle
   integer(kind=i4),  intent(in)  :: l_id !< Local index
   integer(kind=i4),  intent(out) :: g_id !< Global index

   integer(kind=i4) :: block
   integer(kind=i4) :: idx

   character*40, parameter :: caller = "elsi_get_global_row"

   block = (l_id-1)/e_h%blk_row
   idx   = l_id-block*e_h%blk_row
   g_id  = e_h%my_prow*e_h%blk_row+block*e_h%blk_row*e_h%n_prow+idx

end subroutine

!>
!! This routine computes the global column index based on the local column index.
!!
subroutine elsi_get_global_col(e_h,g_id,l_id)

   implicit none

   type(elsi_handle), intent(in)  :: e_h  !< Handle
   integer(kind=i4),  intent(in)  :: l_id !< Local index
   integer(kind=i4),  intent(out) :: g_id !< Global index

   integer(kind=i4) :: block
   integer(kind=i4) :: idx

   character*40, parameter :: caller = "elsi_get_global_col"

   block = (l_id-1)/e_h%blk_col
   idx   = l_id-block*e_h%blk_col
   g_id  = e_h%my_pcol*e_h%blk_col+block*e_h%blk_col*e_h%n_pcol+idx

end subroutine

!>
!! This routine counts the local number of non_zero elements.
!!
subroutine elsi_get_local_nnz_real(e_h,mat,n_row,n_col,nnz)

   implicit none

   type(elsi_handle), intent(in)  :: e_h              !< Handle
   real(kind=r8),     intent(in)  :: mat(n_row,n_col) !< Local matrix
   integer(kind=i4),  intent(in)  :: n_row            !< Local rows
   integer(kind=i4),  intent(in)  :: n_col            !< Local cols
   integer(kind=i4),  intent(out) :: nnz              !< Number of non-zero

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col

   character*40, parameter :: caller = "elsi_get_local_nnz_real"

   nnz = 0

   do i_col = 1,n_col
      do i_row = 1,n_row
         if(abs(mat(i_row,i_col)) > e_h%zero_def) then
            nnz = nnz+1
         endif
      enddo
   enddo

end subroutine

!>
!! This routine counts the local number of non_zero elements.
!!
subroutine elsi_get_local_nnz_complex(e_h,mat,n_row,n_col,nnz)

   implicit none

   type(elsi_handle), intent(in)  :: e_h              !< Handle
   complex(kind=r8),  intent(in)  :: mat(n_row,n_col) !< Local matrix
   integer(kind=i4),  intent(in)  :: n_row            !< Local rows
   integer(kind=i4),  intent(in)  :: n_col            !< Local cols
   integer(kind=i4),  intent(out) :: nnz              !< Number of non-zero

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col

   character*40, parameter :: caller = "elsi_get_local_nnz_complex"

   nnz = 0

   do i_col = 1,n_col
      do i_row = 1,n_row
         if(abs(mat(i_row,i_col)) > e_h%zero_def) then
            nnz = nnz+1
         endif
      enddo
   enddo

end subroutine

!>
!! This routine computes the trace of a matrix. The size of the matrix is
!! restricted to be identical to Hamiltonian.
!!
subroutine elsi_trace_mat_real(e_h,mat,trace)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                        !< Handle
   real(kind=r8),     intent(in)    :: mat(e_h%n_lrow,e_h%n_lcol) !< Matrix
   real(kind=r8),     intent(out)   :: trace                      !< Trace(mat)

   integer(kind=i4) :: i
   integer(kind=i4) :: mpierr
   real(kind=r8)    :: l_trace ! Local result

   character*40, parameter :: caller = "elsi_trace_mat_real"

   do i = 1,e_h%n_basis
      if(e_h%loc_row(i) > 0 .and. e_h%loc_col(i) > 0) then
         l_trace = l_trace + mat(e_h%loc_row(i),e_h%loc_col(i))
      endif
   enddo

   call MPI_Allreduce(l_trace,trace,1,mpi_real8,mpi_sum,e_h%mpi_comm,mpierr)

end subroutine

!>
!! This routine computes the trace of a matrix. The size of the matrix is
!! restricted to be identical to Hamiltonian.
!!
subroutine elsi_trace_mat_complex(e_h,mat,trace)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                        !< Handle
   complex(kind=r8),  intent(in)    :: mat(e_h%n_lrow,e_h%n_lcol) !< Matrix
   complex(kind=r8),  intent(out)   :: trace                      !< Trace(mat)

   integer(kind=i4) :: i
   integer(kind=i4) :: mpierr
   complex(kind=r8) :: l_trace ! Local result

   character*40, parameter :: caller = "elsi_trace_mat_complex"

   do i = 1,e_h%n_basis
      if(e_h%loc_row(i) > 0 .and. e_h%loc_col(i) > 0) then
         l_trace = l_trace + mat(e_h%loc_row(i),e_h%loc_col(i))
      endif
   enddo

   call MPI_Allreduce(l_trace,trace,1,mpi_complex16,mpi_sum,e_h%mpi_comm,mpierr)

end subroutine

!>
!! This routine computes the trace of the product of two matrices. The size of
!! the two matrices is restricted to be identical to Hamiltonian.
!!
subroutine elsi_trace_mat_mat_real(e_h,mat1,mat2,trace)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                         !< Handle
   real(kind=r8),     intent(in)    :: mat1(e_h%n_lrow,e_h%n_lcol) !< Matrix
   real(kind=r8),     intent(in)    :: mat2(e_h%n_lrow,e_h%n_lcol) !< Matrix
   real(kind=r8),     intent(out)   :: trace                       !< Trace(mat1*mat2)

   real(kind=r8)    :: l_trace ! Local result
   integer(kind=i4) :: mpierr

   real(kind=r8), external :: ddot

   character*40, parameter :: caller = "elsi_trace_mat_mat_real"

   l_trace = ddot(e_h%n_lrow*e_h%n_lcol,mat1,1,mat2,1)

   call MPI_Allreduce(l_trace,trace,1,mpi_real8,mpi_sum,e_h%mpi_comm,mpierr)

end subroutine

!>
!! This routine computes the trace of the product of two matrices. The size of
!! the two matrices is restricted to be identical to Hamiltonian.
!!
subroutine elsi_trace_mat_mat_complex(e_h,mat1,mat2,trace)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                         !< Handle
   complex(kind=r8),  intent(in)    :: mat1(e_h%n_lrow,e_h%n_lcol) !< Matrix
   complex(kind=r8),  intent(in)    :: mat2(e_h%n_lrow,e_h%n_lcol) !< Matrix
   complex(kind=r8),  intent(out)   :: trace                       !< Trace(mat1*mat2)

   complex(kind=r8) :: l_trace ! Local result
   integer(kind=i4) :: mpierr

   complex(kind=r8), external :: zdotu

   character*40, parameter :: caller = "elsi_trace_mat_mat_complex"

   l_trace = zdotu(e_h%n_lrow*e_h%n_lcol,mat1,1,mat2,1)

   call MPI_Allreduce(l_trace,trace,1,mpi_complex16,mpi_sum,e_h%mpi_comm,mpierr)

end subroutine

!>
!! This routine initializes the timer.
!!
subroutine elsi_init_timer(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: initial_time
   integer(kind=i4) :: clock_max

   character*40, parameter :: caller = "elsi_init_timer"

   call system_clock(initial_time,e_h%clock_rate,clock_max)

end subroutine

!>
!! This routine gets the current wallclock time.
!!
subroutine elsi_get_time(e_h,wtime)

   implicit none

   type(elsi_handle), intent(in)  :: e_h   !< Handle
   real(kind=r8),     intent(out) :: wtime !< Time

   integer(kind=i4) :: tics

   character*40, parameter :: caller = "elsi_get_time"

   call system_clock(tics)

   wtime = 1.0_r8*tics/e_h%clock_rate

end subroutine

end module ELSI_UTILS
