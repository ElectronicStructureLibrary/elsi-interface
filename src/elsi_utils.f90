! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains a collection of basic utility routines.
!!
module ELSI_UTILS

   use ELSI_CONSTANTS, only: UNSET,FULL_MAT,N_SOLVERS,N_PARALLEL_MODES,&
                             N_MATRIX_FORMATS,MULTI_PROC,SINGLE_PROC,&
                             BLACS_DENSE,PEXSI_CSC,SIESTA_CSC,AUTO,ELPA_SOLVER,&
                             OMM_SOLVER,PEXSI_SOLVER,CHESS_SOLVER,SIPS_SOLVER,&
                             DMP_SOLVER,UNSET_STR,JSON,REAL_DATA,&
                             CMPLX_DATA,STR_LEN,TIME_LEN,&
                             UUID_LEN,COMMA_AFTER,HUMAN
   use ELSI_DATATYPE,  only: elsi_handle
   use ELSI_IO,        only: elsi_say,elsi_init_file_io,elsi_open_json_file
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,&
                             mpi_complex16,mpi_character
   use ELSI_PRECISION, only: i4,i8,r8

   implicit none

   private

   public :: elsi_check
   public :: elsi_check_handle
   public :: elsi_reset_handle
   public :: elsi_get_global_row
   public :: elsi_get_global_col
   public :: elsi_get_global_col_sp2
   public :: elsi_get_local_nnz_real
   public :: elsi_get_local_nnz_cmplx
   public :: elsi_get_solver_tag
   public :: elsi_ready_handle
   public :: elsi_trace_mat_real
   public :: elsi_trace_mat_cmplx
   public :: elsi_trace_mat_mat_real
   public :: elsi_trace_mat_mat_cmplx
   public :: elsi_get_datetime_rfc3339

contains

!>
!! This routine resets an ELSI handle.
!!
subroutine elsi_reset_handle(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   character(len=40), parameter :: caller = "elsi_reset_handle"

   e_h%handle_init            = .false.
   e_h%handle_ready           = .false.
   e_h%solver                 = UNSET
   e_h%matrix_format          = UNSET
   e_h%uplo                   = FULL_MAT
   e_h%parallel_mode          = UNSET
   e_h%print_mem              = .false.
   e_h%n_elsi_calls           = 0
   e_h%myid                   = UNSET
   e_h%myid_all               = UNSET
   e_h%n_procs                = UNSET
   e_h%n_procs_all            = UNSET
   e_h%mpi_comm               = UNSET
   e_h%mpi_comm_all           = UNSET
   e_h%mpi_comm_row           = UNSET
   e_h%mpi_comm_col           = UNSET
   e_h%mpi_ready              = .false.
   e_h%global_mpi_ready       = .false.
   e_h%blacs_ctxt             = UNSET
   e_h%sc_desc                = UNSET
   e_h%blk_row                = UNSET
   e_h%blk_col                = UNSET
   e_h%n_prow                 = UNSET
   e_h%n_pcol                 = UNSET
   e_h%my_prow                = UNSET
   e_h%my_pcol                = UNSET
   e_h%n_lrow                 = UNSET
   e_h%n_lcol                 = UNSET
   e_h%nnz_l                  = UNSET
   e_h%blacs_ready            = .false.
   e_h%nnz_g                  = UNSET
   e_h%nnz_l_sp               = UNSET
   e_h%n_lcol_sp              = UNSET
   e_h%zero_def               = 1.0e-15_r8
   e_h%nnz_l_sp1              = UNSET
   e_h%n_lcol_sp1             = UNSET
   e_h%pexsi_csc_ready        = .false.
   e_h%nnz_l_sp2              = UNSET
   e_h%n_lcol_sp2             = UNSET
   e_h%blk_sp2                = UNSET
   e_h%siesta_csc_ready       = .false.
   e_h%ovlp_is_unit           = .false.
   e_h%ovlp_is_sing           = .false.
   e_h%check_sing             = .true.
   e_h%sing_tol               = 1.0e-5_r8
   e_h%stop_sing              = .false.
   e_h%n_nonsing              = UNSET
   e_h%n_electrons            = 0.0_r8
   e_h%n_basis                = UNSET
   e_h%n_spins                = 1
   e_h%n_kpts                 = 1
   e_h%n_states               = UNSET
   e_h%n_states_solve         = UNSET
   e_h%i_spin                 = 1
   e_h%i_kpt                  = 1
   e_h%i_weight               = 1.0_r8
   e_h%spin_degen             = 0.0_r8
   e_h%energy_hdm             = 0.0_r8
   e_h%energy_sedm            = 0.0_r8
   e_h%mu                     = 0.0_r8
   e_h%ts                     = 0.0_r8
   e_h%broaden_scheme         = 0
   e_h%broaden_width          = 1.0e-2_r8
   e_h%occ_tolerance          = 1.0e-13_r8
   e_h%max_mu_steps           = 100
   e_h%mp_order               = 1
   e_h%spin_is_set            = .false.
   e_h%mu_ready               = .false.
   e_h%ts_ready               = .false.
   e_h%edm_ready_real         = .false.
   e_h%edm_ready_cmplx        = .false.
   e_h%elpa_solver            = UNSET
   e_h%elpa_n_single          = UNSET
   e_h%elpa_output            = .false.
   e_h%elpa_started           = .false.
   e_h%omm_n_states           = UNSET
   e_h%omm_n_elpa             = UNSET
   e_h%omm_flavor             = UNSET
   e_h%omm_ev_shift           = 0.0_r8
   e_h%omm_tol                = 1.0e-12_r8
   e_h%omm_output             = .false.
   e_h%pexsi_np_per_pole      = UNSET
   e_h%pexsi_np_per_point     = UNSET
   e_h%pexsi_my_prow          = UNSET
   e_h%pexsi_my_pcol          = UNSET
   e_h%pexsi_n_prow           = UNSET
   e_h%pexsi_n_pcol           = UNSET
   e_h%pexsi_my_point         = UNSET
   e_h%pexsi_myid_point       = UNSET
   e_h%pexsi_comm_among_pole  = UNSET
   e_h%pexsi_comm_in_pole     = UNSET
   e_h%pexsi_comm_among_point = UNSET
   e_h%pexsi_comm_in_point    = UNSET
   e_h%pexsi_ne               = 0.0_r8
   e_h%pexsi_started          = .false.
   e_h%sips_n_elpa            = UNSET
   e_h%sips_np_per_slice      = UNSET
   e_h%sips_n_slices          = UNSET
   e_h%sips_slice_type        = UNSET
   e_h%sips_first_ev          = UNSET
   e_h%sips_buffer            = 0.0_r8
   e_h%sips_interval          = 0.0_r8
   e_h%sips_inertia_tol       = 0.0_r8
   e_h%sips_do_inertia        = .false.
   e_h%sips_started           = .false.
   e_h%dmp_n_states           = UNSET
   e_h%dmp_method             = UNSET
   e_h%dmp_max_power          = UNSET
   e_h%dmp_max_iter           = UNSET
   e_h%dmp_ev_ham_max         = 0.0_r8
   e_h%dmp_ev_ham_min         = 0.0_r8
   e_h%dmp_tol                = 1e-8_r8
   e_h%dmp_ne                 = 0.0_r8
   e_h%output_timings_file    = .false.
   e_h%solver_timings_unit    = UNSET
   e_h%solver_timings_name    = UNSET_STR
   e_h%calling_code           = UNSET_STR
   e_h%calling_code_ver       = UNSET_STR
   e_h%uuid                   = UNSET_STR
   e_h%uuid_exists            = .false.

end subroutine

!>
!! This routine guarantees that there are no mutually conflicting parameters.
!!
subroutine elsi_check(e_h,caller)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   character(len=*),  intent(in)    :: caller

   ! General check of solver, parallel mode, matrix format
   if(e_h%solver < 0 .or. e_h%solver >= N_SOLVERS) then
      call elsi_stop(" Unsupported solver.",e_h,caller)
   endif

   if(e_h%parallel_mode < 0 .or. e_h%parallel_mode >= N_PARALLEL_MODES) then
      call elsi_stop(" Unsupported parallel mode.",e_h,caller)
   endif

   if(e_h%matrix_format < 0 .or. e_h%matrix_format >= N_MATRIX_FORMATS) then
      call elsi_stop(" Unsupported matirx format.",e_h,caller)
   endif

   if(e_h%uplo /= FULL_MAT) then
      call elsi_stop(" Triangular matrix input not yet supported.",e_h,caller)
   endif

   ! Spin and k-point
   if(e_h%n_spins*e_h%n_kpts > 1) then
      if(.not. e_h%global_mpi_ready) then
         call elsi_stop(" Spin/k-point calculations require a global MPI"//&
                 " communicator.",e_h,caller)
      endif
   endif

   if(.not. e_h%spin_is_set) then
      if(e_h%n_spins == 2) then
         e_h%spin_degen = 1.0_r8
      else
         e_h%spin_degen = 2.0_r8
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
         call elsi_stop(" BLACS not properly set up.",e_h,caller)
      endif
   elseif(e_h%matrix_format == SIESTA_CSC) then
      if(e_h%blk_sp2 == UNSET .or. .not. e_h%siesta_csc_ready) then
         call elsi_stop(" SIESTA_CSC format not properly set up.",e_h,caller)
      endif
   elseif(e_h%matrix_format == PEXSI_CSC) then
      if(.not. e_h%pexsi_csc_ready) then
         call elsi_stop(" PEXSI_CSC format not properly set up.",e_h,caller)
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

      if(e_h%n_basis < e_h%pexsi_np_per_pole) then
         call elsi_stop(" For this number of MPI tasks, the matrix size is"//&
                 " too small to use PEXSI.",e_h,caller)
      endif

      if(e_h%pexsi_np_per_pole == UNSET) then
         if(mod(e_h%n_procs,e_h%pexsi_options%numPole*&
            e_h%pexsi_options%nPoints) /= 0) then
            call elsi_stop(" To use PEXSI, the total number of MPI tasks"//&
                    " must be a multiple of the number of MPI tasks per pole"//&
                    " times the number of mu points.",e_h,caller)
         endif
      else
         if(mod(e_h%n_procs,e_h%pexsi_np_per_pole*&
            e_h%pexsi_options%nPoints) /= 0) then
            call elsi_stop(" To use PEXSI, the total number of MPI tasks"//&
                    " must be a multiple of the number of MPI tasks per pole"//&
                    " times the number of mu points.",e_h,caller)
         endif

         if(e_h%pexsi_np_per_pole*e_h%pexsi_options%numPole*&
            e_h%pexsi_options%nPoints < e_h%n_procs) then
            call elsi_stop(" Specified number of MPI tasks per pole is too"//&
                    " small for the total number of MPI tasks.",e_h,caller)
         endif
      endif
   case(CHESS_SOLVER)
      call elsi_stop(" CheSS solver not yet available.",e_h,caller)
   case(SIPS_SOLVER)
      if(e_h%n_basis < e_h%n_procs) then
         call elsi_stop(" For this number of MPI tasks, the matrix size is"//&
                 " too small to use SIPS.",e_h,caller)
      endif

      if(e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(" SIPS solver requires MULTI_PROC parallel mode.",e_h,&
                 caller)
      endif

      if(e_h%n_spins > 1) then
         call elsi_stop(" Spin-polarized case not yet supported with SIPS.",&
                 e_h,caller)
      endif

      if(e_h%n_kpts > 1) then
         call elsi_stop(" k-points not yet supported with SIPS.",e_h,caller)
      endif
   case(DMP_SOLVER)
      if(e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(" DMP solver requires MULTI_PROC parallel mode.",e_h,&
                 caller)
      endif

      if(e_h%n_spins > 1) then
         call elsi_stop(" Spin-polarized case not yet supported with DMP.",e_h,&
                 caller)
      endif

      if(e_h%n_kpts > 1) then
         call elsi_stop(" k-points not yet supported with DMP.",e_h,caller)
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

   type(elsi_handle), intent(in) :: e_h
   character(len=*),  intent(in) :: caller

   if(.not. e_h%handle_init) then
      call elsi_stop(" Invalid handle! Not initialized.",e_h,caller)
   endif

end subroutine

!>
!! This subroutine signifies that a handle is believed to be ready to be used.
!! There are certain tasks (such as IO) that are only possible once the user has
!! specified sufficient information about the problem, but they may call the
!! relevant mutators in any order, making it difficult to pin down exactly when
!! certain initialization-like tasks can be done. Thus, this subroutine is
!! called at the beginning of all public-facing subroutines in which ELSI is
!! executing a task that would require it to have been sufficiently initialized
!! by the user: currently, only solvers. This does not apply to mutators,
!! initalization, finalization, or matrix IO.
!!
subroutine elsi_ready_handle(e_h,caller)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   character(len=*),  intent(in)    :: caller

   character(len=UUID_LEN) :: uuid_temp

   call elsi_check_handle(e_h,caller)

   if(.not. e_h%handle_ready) then
      call elsi_check(e_h,caller)

      e_h%handle_ready   = .true.

      ! Perform initialization-like tasks which require MPI

      ! First, solver timings
      if(e_h%output_timings_file) then
         if(e_h%myid_all == 0) then
            ! Open the JSON file to output timings to (the .false. statement
            ! selects human-readable output, kept for debugging purposes)
            if(.true.) then
               call elsi_open_json_file(e_h,e_h%solver_timings_unit,&
                       e_h%solver_timings_name,.true.,e_h%timings_file)
            else
               call elsi_init_file_io(e_h%timings_file,e_h%solver_timings_unit,&
                       e_h%solver_timings_name,HUMAN,.true.,"",COMMA_AFTER)
               open(unit=e_h%timings_file%print_unit,&
                  file=e_h%timings_file%file_name)
            endif
         else
            ! De-initialize e_h%myid_all /= 0
            e_h%timings_file%print_unit = UNSET
            e_h%timings_file%file_name  = UNSET_STR
         endif
      endif

      ! Generate UUID
      if (.not. e_h%uuid_exists) then
         call elsi_gen_uuid(e_h)
         e_h%uuid_exists = .true.
      endif
   endif

end subroutine

!>
!! This routine computes the global row index based on the local row index.
!!
subroutine elsi_get_global_row(e_h,g_id,l_id)

   implicit none

   type(elsi_handle), intent(in)  :: e_h
   integer(kind=i4),  intent(in)  :: l_id
   integer(kind=i4),  intent(out) :: g_id

   integer(kind=i4) :: blk
   integer(kind=i4) :: idx

   character(len=40), parameter :: caller = "elsi_get_global_row"

   blk  = (l_id-1)/e_h%blk_row
   idx  = l_id-blk*e_h%blk_row
   g_id = e_h%my_prow*e_h%blk_row+blk*e_h%blk_row*e_h%n_prow+idx

end subroutine

!>
!! This routine computes the global column index based on the local column index.
!!
subroutine elsi_get_global_col(e_h,g_id,l_id)

   implicit none

   type(elsi_handle), intent(in)  :: e_h
   integer(kind=i4),  intent(in)  :: l_id
   integer(kind=i4),  intent(out) :: g_id

   integer(kind=i4) :: blk
   integer(kind=i4) :: idx

   character(len=40), parameter :: caller = "elsi_get_global_col"

   blk  = (l_id-1)/e_h%blk_col
   idx  = l_id-blk*e_h%blk_col
   g_id = e_h%my_pcol*e_h%blk_col+blk*e_h%blk_col*e_h%n_pcol+idx

end subroutine

!>
!! This routine computes the global column index based on the local column index.
!!
subroutine elsi_get_global_col_sp2(e_h,g_id,l_id)

   implicit none

   type(elsi_handle), intent(in)  :: e_h
   integer(kind=i4),  intent(in)  :: l_id
   integer(kind=i4),  intent(out) :: g_id

   integer(kind=i4) :: blk
   integer(kind=i4) :: idx

   character(len=40), parameter :: caller = "elsi_get_global_col_sp2"

   blk  = (l_id-1)/e_h%blk_sp2
   idx  = l_id-blk*e_h%blk_sp2
   g_id = e_h%myid*e_h%blk_sp2+blk*e_h%blk_sp2*e_h%n_procs+idx

end subroutine

!>
!! This routine counts the local number of non_zero elements.
!!
subroutine elsi_get_local_nnz_real(e_h,mat,n_row,n_col,nnz)

   implicit none

   type(elsi_handle), intent(in)  :: e_h
   real(kind=r8),     intent(in)  :: mat(n_row,n_col)
   integer(kind=i4),  intent(in)  :: n_row
   integer(kind=i4),  intent(in)  :: n_col
   integer(kind=i4),  intent(out) :: nnz

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col

   character(len=40), parameter :: caller = "elsi_get_local_nnz_real"

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
subroutine elsi_get_local_nnz_cmplx(e_h,mat,n_row,n_col,nnz)

   implicit none

   type(elsi_handle), intent(in)  :: e_h
   complex(kind=r8),  intent(in)  :: mat(n_row,n_col)
   integer(kind=i4),  intent(in)  :: n_row
   integer(kind=i4),  intent(in)  :: n_col
   integer(kind=i4),  intent(out) :: nnz

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col

   character(len=40), parameter :: caller = "elsi_get_local_nnz_cmplx"

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

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(in)    :: mat(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(out)   :: trace

   integer(kind=i4) :: i
   integer(kind=i4) :: ierr
   real(kind=r8)    :: l_trace ! Local result

   character(len=40), parameter :: caller = "elsi_trace_mat_real"

   l_trace = 0.0_r8

   do i = 1,e_h%n_basis
      if(e_h%loc_row(i) > 0 .and. e_h%loc_col(i) > 0) then
         l_trace = l_trace + mat(e_h%loc_row(i),e_h%loc_col(i))
      endif
   enddo

   call MPI_Allreduce(l_trace,trace,1,mpi_real8,mpi_sum,e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

end subroutine

!>
!! This routine computes the trace of a matrix. The size of the matrix is
!! restricted to be identical to Hamiltonian.
!!
subroutine elsi_trace_mat_cmplx(e_h,mat,trace)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(in)    :: mat(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(out)   :: trace

   integer(kind=i4) :: i
   integer(kind=i4) :: ierr
   complex(kind=r8) :: l_trace ! Local result

   character(len=40), parameter :: caller = "elsi_trace_mat_cmplx"

   l_trace = 0.0_r8

   do i = 1,e_h%n_basis
      if(e_h%loc_row(i) > 0 .and. e_h%loc_col(i) > 0) then
         l_trace = l_trace + mat(e_h%loc_row(i),e_h%loc_col(i))
      endif
   enddo

   call MPI_Allreduce(l_trace,trace,1,mpi_complex16,mpi_sum,e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

end subroutine

!>
!! This routine computes the trace of the product of two matrices. The size of
!! the two matrices is restricted to be identical to Hamiltonian.
!!
subroutine elsi_trace_mat_mat_real(e_h,mat1,mat2,trace)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(in)    :: mat1(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(in)    :: mat2(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(out)   :: trace

   real(kind=r8)    :: l_trace ! Local result
   integer(kind=i4) :: ierr

   real(kind=r8), external :: ddot

   character(len=40), parameter :: caller = "elsi_trace_mat_mat_real"

   l_trace = ddot(e_h%n_lrow*e_h%n_lcol,mat1,1,mat2,1)

   call MPI_Allreduce(l_trace,trace,1,mpi_real8,mpi_sum,e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

end subroutine

!>
!! This routine computes the trace of the product of two matrices. The size of
!! the two matrices is restricted to be identical to Hamiltonian.
!!
subroutine elsi_trace_mat_mat_cmplx(e_h,mat1,mat2,trace)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(in)    :: mat1(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(in)    :: mat2(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(out)   :: trace

   complex(kind=r8) :: l_trace ! Local result
   integer(kind=i4) :: ierr

   complex(kind=r8), external :: zdotu

   character(len=40), parameter :: caller = "elsi_trace_mat_mat_cmplx"

   l_trace = zdotu(e_h%n_lrow*e_h%n_lcol,mat1,1,mat2,1)

   call MPI_Allreduce(l_trace,trace,1,mpi_complex16,mpi_sum,e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

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
            elseif(e_h%elpa_solver == 2) then
               tag = "ELPA_2STAGE_REAL"
            else
               call elsi_stop(" Unsupported ELPA flavor.",e_h,caller)
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
      case default
         call elsi_stop(" Unsupported solver.",e_h,caller)
      end select
   elseif(data_type == CMPLX_DATA) then
      select case(e_h%solver)
      case(ELPA_SOLVER)
         if(e_h%parallel_mode == SINGLE_PROC) then
            tag = "LAPACK_CMPLX"
         else ! MULTI_PROC
            if(e_h%elpa_solver == 1) then
               tag = "ELPA_1STAGE_CMPLX"
            elseif(e_h%elpa_solver == 2) then
               tag = "ELPA_2STAGE_CMPLX"
            else
               call elsi_stop(" Unsupported ELPA flavor.",e_h,caller)
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
      case default
         call elsi_stop(" Unsupported solver.",e_h,caller)
      end select
   else
      call elsi_stop(" Unsupported data type.",e_h,caller)
   endif

end subroutine

!>
!! This routine returns the current date and time, formatted as an RFC3339
!! string with time zone offset.
!!
subroutine elsi_get_datetime_rfc3339(datetime_rfc3339)

   implicit none

   character(len=TIME_LEN), intent(out) :: datetime_rfc3339

   integer(kind=i4) :: datetime(8)
   integer(kind=i4) :: temp_int
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
   temp_int = mod(datetime(4),60)
   if(temp_int < 10) then
      write(timezone_min,'(A1,I1)') "0",temp_int
   else
      write(timezone_min,'(I2)') temp_int
   endif

   ! Get timezone hours
   temp_int = datetime(4)/60
   if(temp_int < 10) then
      write(timezone_hour,'(A1,I1)') "0",temp_int
   else
      write(timezone_hour,'(I2)') temp_int
   endif

   write(datetime_rfc3339,&
      '(A4,A1,A2,A1,A2,A1,A2,A1,A2,A1,A2,A1,A3,A1,A2,A1,A2)')&
      year,"-",month,"-",day,"T",hour,":",minute,":",second,".",millisecond,&
      timezone_sign,timezone_hour,":",timezone_min

end subroutine

!>
!! Generate a UUID (unique identifier) in RFC 4122 format. The random seed must
!! have been set.
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

end module ELSI_UTILS
