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
                             DMP_SOLVER,UNSET_STR
   use ELSI_DATATYPE,  only: elsi_handle
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,&
                             mpi_complex16,mpi_comm_self
   use ELSI_PRECISION, only: i4,r8

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
   public :: elsi_trace_mat_real
   public :: elsi_trace_mat_cmplx
   public :: elsi_trace_mat_mat_real
   public :: elsi_trace_mat_mat_cmplx

contains

!>
!! This routine resets an ELSI handle.
!!
subroutine elsi_reset_handle(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   character(len=40), parameter :: caller = "elsi_reset_handle"

   e_h%handle_init            = .false.
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
   e_h%edm_ready_real         = .false.
   e_h%edm_ready_cmplx        = .false.
   e_h%elpa_solver            = UNSET
   e_h%elpa_n_single          = UNSET
   e_h%elpa_gpu               = .false.
   e_h%elpa_output            = .false.
   e_h%elpa_started           = .false.
   e_h%omm_n_states           = UNSET
   e_h%omm_n_elpa             = UNSET
   e_h%omm_flavor             = UNSET
   e_h%omm_ev_shift           = 0.0_r8
   e_h%omm_tol                = 1.0e-12_r8
   e_h%omm_output             = .false.
   e_h%omm_started            = .false.
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
   e_h%dmp_started            = .false.
   e_h%caller                 = UNSET_STR
   e_h%uuid                   = UNSET_STR
   e_h%uuid_exists            = .false.

end subroutine

!>
!! This routine guarantees that there are no unsupported or mutually conflicting
!! parameters before running actual calculations.
!!
subroutine elsi_check(e_h,caller)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   character(len=*),  intent(in)    :: caller

   ! General check of solver, parallel mode, matrix format
   if(e_h%solver < 0 .or. e_h%solver >= N_SOLVERS) then
      call elsi_stop(e_h,"Unsupported solver.",caller)
   endif

   if(e_h%parallel_mode < 0 .or. e_h%parallel_mode >= N_PARALLEL_MODES) then
      call elsi_stop(e_h,"Unsupported parallel mode.",caller)
   endif

   if(e_h%matrix_format < 0 .or. e_h%matrix_format >= N_MATRIX_FORMATS) then
      call elsi_stop(e_h,"Unsupported matirx format.",caller)
   endif

   if(e_h%uplo /= FULL_MAT) then
      call elsi_stop(e_h,"Triangular matrix input not yet supported.",caller)
   endif

   ! Spin and k-point
   if(e_h%n_spins*e_h%n_kpts > 1) then
      if(.not. e_h%global_mpi_ready) then
         call elsi_stop(e_h,"Spin/k-point calculations require a global MPI"//&
                 " communicator.",caller)
      endif
   endif

   if(.not. e_h%spin_is_set) then
      if(e_h%n_spins == 2) then
         e_h%spin_degen = 1.0_r8
      else
         e_h%spin_degen = 2.0_r8
      endif
   endif

   if(.not. e_h%mpi_ready) then
      e_h%mpi_comm = mpi_comm_self
      e_h%n_procs  = 1
      e_h%myid     = 0
   endif

   if(.not. e_h%global_mpi_ready) then
      e_h%mpi_comm_all = e_h%mpi_comm
      e_h%n_procs_all  = e_h%n_procs
      e_h%myid_all     = e_h%myid
   endif

   if(e_h%myid_all /= 0) then
      e_h%stdio%print_info    = .false.
      e_h%log_file%print_info = .false.
   endif

   if(e_h%parallel_mode == MULTI_PROC) then
      if(.not. e_h%mpi_ready) then
         call elsi_stop(e_h,"MULTI_PROC parallel mode requires MPI.",caller)
      endif
   endif

   if(e_h%matrix_format == BLACS_DENSE) then
      if(.not. e_h%blacs_ready .and. e_h%parallel_mode /= SINGLE_PROC) then
         call elsi_stop(e_h,"BLACS not properly set up.",caller)
      endif
   elseif(e_h%matrix_format == SIESTA_CSC) then
      if(e_h%blk_sp2 == UNSET .or. .not. e_h%siesta_csc_ready) then
         call elsi_stop(e_h,"SIESTA_CSC not properly set up.",caller)
      endif
   elseif(e_h%matrix_format == PEXSI_CSC) then
      if(.not. e_h%pexsi_csc_ready) then
         call elsi_stop(e_h,"PEXSI_CSC not properly set up.",caller)
      endif
   endif

   ! Specific check for each solver
   select case(e_h%solver)
   case(AUTO)
      call elsi_stop(e_h,"Solver auto-selection not yet available.",caller)
   case(ELPA_SOLVER)
      ! Nothing
   case(OMM_SOLVER)
      if(e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(e_h,"libOMM solver requires MULTI_PROC parallel mode.",&
                 caller)
      endif
   case(PEXSI_SOLVER)
      if(e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(e_h,"PEXSI solver requires MULTI_PROC parallel mode.",&
                 caller)
      endif

      if(e_h%n_basis < e_h%pexsi_np_per_pole) then
         call elsi_stop(e_h,"For this number of MPI tasks, the matrix size"//&
                 " is too small to use PEXSI.",caller)
      endif

      if(e_h%pexsi_np_per_pole == UNSET) then
         if(mod(e_h%n_procs,e_h%pexsi_options%numPole*&
            e_h%pexsi_options%nPoints) /= 0) then
            call elsi_stop(e_h,"To use PEXSI, the total number of MPI tasks"//&
                    " must be a multiple of the number of MPI tasks per pole"//&
                    " times the number of mu points.",caller)
         endif
      else
         if(mod(e_h%n_procs,e_h%pexsi_np_per_pole*&
            e_h%pexsi_options%nPoints) /= 0) then
            call elsi_stop(e_h,"To use PEXSI, the total number of MPI tasks"//&
                    " must be a multiple of the number of MPI tasks per pole"//&
                    " times the number of mu points.",caller)
         endif

         if(e_h%pexsi_np_per_pole*e_h%pexsi_options%numPole*&
            e_h%pexsi_options%nPoints < e_h%n_procs) then
            call elsi_stop(e_h,"Specified number of MPI tasks per pole is"//&
                    " too small for the total number of MPI tasks.",caller)
         endif
      endif
   case(CHESS_SOLVER)
      call elsi_stop(e_h,"CheSS solver not yet supported.",caller)
   case(SIPS_SOLVER)
      if(e_h%n_basis < e_h%n_procs) then
         call elsi_stop(e_h,"For this number of MPI tasks, the matrix size"//&
                 " is too small to use SIPS.",caller)
      endif

      if(e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(e_h,"SIPS solver requires MULTI_PROC parallel mode.",&
                 caller)
      endif

      if(e_h%n_spins > 1) then
         call elsi_stop(e_h,"Spin-polarized case not yet supported with SIPS.",&
                 caller)
      endif

      if(e_h%n_kpts > 1) then
         call elsi_stop(e_h,"k-points not yet supported with SIPS.",caller)
      endif
   case(DMP_SOLVER)
      if(e_h%parallel_mode /= MULTI_PROC) then
         call elsi_stop(e_h,"DMP solver requires MULTI_PROC parallel mode.",&
                 caller)
      endif

      if(e_h%n_spins > 1) then
         call elsi_stop(e_h,"Spin-polarized case not yet supported with DMP.",&
                 caller)
      endif

      if(e_h%n_kpts > 1) then
         call elsi_stop(e_h,"k-points not yet supported with DMP.",caller)
      endif
   case default
      call elsi_stop(e_h,"Unsupported solver.",caller)
   end select

end subroutine

!>
!! This routine checks whether a handle has been properly initialized.
!!
subroutine elsi_check_handle(e_h,caller)

   implicit none

   type(elsi_handle), intent(in) :: e_h
   character(len=*),  intent(in) :: caller

   if(.not. e_h%handle_init) then
      call elsi_stop(e_h,"Invalid handle! Not initialized.",caller)
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

end module ELSI_UTILS
