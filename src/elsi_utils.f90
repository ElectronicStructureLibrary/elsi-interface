! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains a collection of basic utility routines.
!!
module ELSI_UTILS

   use ELSI_CONSTANTS, only: UNSET,UT_MAT,LT_MAT,N_SOLVERS,N_PARALLEL_MODES,&
                             N_MATRIX_FORMATS,MULTI_PROC,SINGLE_PROC,&
                             BLACS_DENSE,PEXSI_CSC,SIESTA_CSC,AUTO,ELPA_SOLVER,&
                             OMM_SOLVER,PEXSI_SOLVER,CHESS_SOLVER,SIPS_SOLVER,&
                             NTPOLY_SOLVER
   use ELSI_DATATYPE,  only: elsi_param_t,elsi_basic_t
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,&
                             mpi_complex16,mpi_comm_self
   use ELSI_PRECISION, only: i4,r8

   implicit none

   private

   public :: elsi_check
   public :: elsi_check_init
   public :: elsi_reset_param
   public :: elsi_reset_basic
   public :: elsi_get_gid
   public :: elsi_get_lid
   public :: elsi_get_nnz
   public :: elsi_trace_mat
   public :: elsi_trace_mat_mat
   public :: elsi_set_full_mat

   interface elsi_get_nnz
      module procedure elsi_get_nnz_real
      module procedure elsi_get_nnz_cmplx
   end interface

   interface elsi_trace_mat
      module procedure elsi_trace_mat_real
      module procedure elsi_trace_mat_cmplx
   end interface

   interface elsi_trace_mat_mat
      module procedure elsi_trace_mat_mat_real
      module procedure elsi_trace_mat_mat_cmplx
   end interface

   interface elsi_set_full_mat
      module procedure elsi_set_full_mat_real
      module procedure elsi_set_full_mat_cmplx
   end interface

contains

!>
!! This routine resets ELSI runtime parameters.
!!
subroutine elsi_reset_param(ph)

   implicit none

   type(elsi_param_t), intent(out) :: ph

   character(len=40), parameter :: caller = "elsi_reset_handle"

   ph%solver                 = UNSET
   ph%matrix_format          = UNSET
   ph%parallel_mode          = UNSET
   ph%n_calls                = 0
   ph%ovlp_is_unit           = .false.
   ph%ovlp_is_sing           = .false.
   ph%check_sing             = .true.
   ph%sing_tol               = 1.0e-5_r8
   ph%stop_sing              = .false.
   ph%n_good                 = UNSET
   ph%n_electrons            = 0.0_r8
   ph%n_basis                = UNSET
   ph%n_spins                = 1
   ph%n_kpts                 = 1
   ph%n_states               = UNSET
   ph%n_states_solve         = UNSET
   ph%i_spin                 = 1
   ph%i_kpt                  = 1
   ph%i_weight               = 1.0_r8
   ph%spin_degen             = 0.0_r8
   ph%spin_is_set            = .false.
   ph%ebs                    = 0.0_r8
   ph%edm_ready_real         = .false.
   ph%edm_ready_cmplx        = .false.
   ph%mu                     = 0.0_r8
   ph%ts                     = 0.0_r8
   ph%mu_scheme              = 0
   ph%mu_width               = 1.0e-2_r8
   ph%mu_tol                 = 1.0e-13_r8
   ph%mu_max_steps           = 100
   ph%mu_mp_order            = 1
   ph%elpa_solver            = 2
   ph%elpa_n_single          = 0
   ph%elpa_comm_row          = UNSET
   ph%elpa_comm_col          = UNSET
   ph%elpa_gpu               = .false.
   ph%elpa_gpu_kernels       = .false.
   ph%elpa_autotune          = .true.
   ph%elpa_output            = .false.
   ph%elpa_started           = .false.
   ph%omm_n_lrow             = UNSET
   ph%omm_n_states           = UNSET
   ph%omm_n_elpa             = 5
   ph%omm_flavor             = 0
   ph%omm_desc               = UNSET
   ph%omm_tol                = 1.0e-12_r8
   ph%omm_output             = .false.
   ph%omm_started            = .false.
   ph%pexsi_np_per_pole      = UNSET
   ph%pexsi_np_per_point     = UNSET
   ph%pexsi_my_prow          = UNSET
   ph%pexsi_my_pcol          = UNSET
   ph%pexsi_n_prow           = UNSET
   ph%pexsi_n_pcol           = UNSET
   ph%pexsi_my_point         = UNSET
   ph%pexsi_myid_point       = UNSET
   ph%pexsi_comm_intra_pole  = UNSET
   ph%pexsi_comm_inter_pole  = UNSET
   ph%pexsi_comm_inter_point = UNSET
   ph%pexsi_ne               = 0.0_r8
   ph%pexsi_started          = .false.
   ph%sips_n_elpa            = 0
   ph%sips_np_per_slice      = UNSET
   ph%sips_n_slices          = UNSET
   ph%sips_slice_type        = 2
   ph%sips_first_ev          = 1
   ph%sips_buffer            = 1.0e-2_r8
   ph%sips_interval(1)       = -2.0_r8
   ph%sips_interval(2)       = 2.0_r8
   ph%sips_inertia_tol       = 1.0e-3_r8
   ph%sips_do_inertia        = .true.
   ph%sips_started           = .false.
   ph%nt_n_group             = 1
   ph%nt_method              = 0
   ph%nt_max_iter            = 100
   ph%nt_tol                 = 1.0e-8_r8
   ph%nt_filter              = 1.0e-15_r8
   ph%nt_output              = .false.
   ph%nt_started             = .false.

end subroutine

!>
!! This routine resets ELSI basic information.
!!
subroutine elsi_reset_basic(bh)

   implicit none

   type(elsi_basic_t), intent(out) :: bh

   character(len=40), parameter :: caller = "elsi_reset_basic"

   bh%print_info       = 0
   bh%print_unit       = 6
   bh%print_json       = 0
   bh%json_init        = .false.
   bh%user_tag         = ""
   bh%uuid             = ""
   bh%uuid_ready       = .false.
   bh%myid             = UNSET
   bh%myid_all         = UNSET
   bh%n_procs          = UNSET
   bh%n_procs_all      = UNSET
   bh%comm             = UNSET
   bh%comm_all         = UNSET
   bh%mpi_ready        = .false.
   bh%mpi_all_ready    = .false.
   bh%blacs_ctxt       = UNSET
   bh%desc             = UNSET
   bh%blk              = UNSET
   bh%n_prow           = UNSET
   bh%n_pcol           = UNSET
   bh%my_prow          = UNSET
   bh%my_pcol          = UNSET
   bh%n_lrow           = UNSET
   bh%n_lcol           = UNSET
   bh%nnz_l            = UNSET
   bh%blacs_ready      = .false.
   bh%nnz_g            = UNSET
   bh%nnz_l_sp         = UNSET
   bh%n_lcol_sp        = UNSET
   bh%def0             = 1.0e-15_r8
   bh%nnz_l_sp1        = UNSET
   bh%n_lcol_sp1       = UNSET
   bh%pexsi_csc_ready  = .false.
   bh%nnz_l_sp2        = UNSET
   bh%n_lcol_sp2       = UNSET
   bh%blk_sp2          = UNSET
   bh%siesta_csc_ready = .false.

end subroutine

!>
!! This routine guarantees that there are no unsupported or mutually conflicting
!! parameters before running actual calculations.
!!
subroutine elsi_check(ph,bh,caller)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   character(len=*),   intent(in)    :: caller

   ! General check of solver, parallel mode, matrix format
   if(ph%solver < 0 .or. ph%solver >= N_SOLVERS) then
      call elsi_stop(bh,"Unsupported solver.",caller)
   endif

   if(ph%parallel_mode < 0 .or. ph%parallel_mode >= N_PARALLEL_MODES) then
      call elsi_stop(bh,"Unsupported parallel mode.",caller)
   endif

   if(ph%matrix_format < 0 .or. ph%matrix_format >= N_MATRIX_FORMATS) then
      call elsi_stop(bh,"Unsupported matirx format.",caller)
   endif

   ! Spin
   if(ph%n_spins > 1 .and. .not. bh%mpi_all_ready) then
      call elsi_stop(bh,"Calculations with two spin channels require a"//&
              " global MPI communicator.",caller)
   endif

   ! k-point
   if(ph%n_kpts > 1 .and. .not. bh%mpi_all_ready) then
      call elsi_stop(bh,"Calculations with multiple k-points require a"//&
              " global MPI communicator.",caller)
   endif

   if(.not. ph%spin_is_set) then
      if(ph%n_spins == 2) then
         ph%spin_degen = 1.0_r8
      else
         ph%spin_degen = 2.0_r8
      endif
   endif

   if(.not. bh%mpi_ready) then
      bh%comm    = mpi_comm_self
      bh%n_procs = 1
      bh%myid    = 0
   endif

   if(.not. bh%mpi_all_ready) then
      bh%comm_all    = bh%comm
      bh%n_procs_all = bh%n_procs
      bh%myid_all    = bh%myid
   endif

   if(bh%myid_all /= 0) then
      bh%print_info = 0
   endif

   if(ph%parallel_mode == MULTI_PROC) then
      if(.not. bh%mpi_ready) then
         call elsi_stop(bh,"MULTI_PROC parallel mode requires MPI.",caller)
      endif
   endif

   if(ph%matrix_format == BLACS_DENSE) then
      if(.not. bh%blacs_ready .and. ph%parallel_mode /= SINGLE_PROC) then
         call elsi_stop(bh,"BLACS matrix format not properly set up.",caller)
      endif
   elseif(ph%matrix_format == SIESTA_CSC) then
      if(.not. bh%siesta_csc_ready) then
         call elsi_stop(bh,"SIESTA_CSC matrix format not properly set up.",&
                 caller)
      endif

      if(bh%blk_sp2 == UNSET) then
         call elsi_stop(bh,"Block size should be set for SIESTA_CSC matrix"//&
                 " format.",caller)
      endif
   elseif(ph%matrix_format == PEXSI_CSC) then
      if(.not. bh%pexsi_csc_ready) then
         call elsi_stop(bh,"PEXSI_CSC matrix format not properly set up.",&
                 caller)
      endif

      if(ph%solver == PEXSI_SOLVER .and. ph%pexsi_np_per_pole == UNSET) then
         call elsi_stop(bh,"Number of MPI tasks per pole should be set for"//&
                 " PEXSI_CSC matrix format and PEXSI solver.",caller)
      endif
   endif

   ! Specific check for each solver
   select case(ph%solver)
   case(ELPA_SOLVER)
      ! Nothing
   case(OMM_SOLVER)
      if(ph%parallel_mode /= MULTI_PROC) then
         call elsi_stop(bh,"libOMM requires MULTI_PROC parallel mode.",caller)
      endif
   case(PEXSI_SOLVER)
      if(ph%parallel_mode /= MULTI_PROC) then
         call elsi_stop(bh,"PEXSI requires MULTI_PROC parallel mode.",caller)
      endif

      if(mod(bh%n_procs,ph%pexsi_options%nPoints) /= 0) then
         call elsi_stop(bh,"To use PEXSI, number of mu points must be a"//&
                 " divisor of number of MPI tasks.",caller)
      endif

      if(ph%pexsi_np_per_pole /= UNSET) then
         if(mod(bh%n_procs,ph%pexsi_np_per_pole*&
            ph%pexsi_options%nPoints) /= 0) then
            call elsi_stop(bh,"To use PEXSI, specified number of MPI tasks"//&
                    " per pole times number of mu points must be a divisor"//&
                    " of number of MPI tasks.",caller)
         endif

         if(ph%pexsi_np_per_pole*ph%pexsi_options%numPole*&
            ph%pexsi_options%nPoints < bh%n_procs) then
            call elsi_stop(bh,"Specified number of MPI tasks per pole too"//&
                    " small for this number of MPI tasks.",caller)
         endif
      endif
   case(CHESS_SOLVER)
      call elsi_stop(bh,"CheSS not yet supported.",caller)
   case(SIPS_SOLVER)
      if(ph%n_basis < bh%n_procs) then
         call elsi_stop(bh,"Matrix size too small to use SLEPc-SIPs with"//&
                 " this number of MPI tasks.",caller)
      endif

      if(ph%parallel_mode /= MULTI_PROC) then
         call elsi_stop(bh,"SLEPc-SIPs requires MULTI_PROC parallel mode.",&
                 caller)
      endif

      if(ph%n_spins > 1) then
         call elsi_stop(bh,"Calculations with two spin channels not yet"//&
                 " supported with SLEPc-SIPs.",caller)
      endif

      if(ph%n_kpts > 1) then
         call elsi_stop(bh,"Calculations with multiple k-points not yet"//&
                 " supported with SLEPc-SIPs.",caller)
      endif
   case(NTPOLY_SOLVER)
      if(ph%parallel_mode /= MULTI_PROC) then
         call elsi_stop(bh,"NTPoly requires MULTI_PROC parallel mode.",caller)
      endif

      if(ph%n_spins > 1) then
         call elsi_stop(bh,"Calculations with two spin channels not yet"//&
                 " supported with NTPoly.",caller)
      endif

      if(ph%n_kpts > 1) then
         call elsi_stop(bh,"Calculations with multiple k-points not yet"//&
                 " supported with NTPoly.",caller)
      endif
   end select

end subroutine

!>
!! This routine checks whether a handle has been properly initialized.
!!
subroutine elsi_check_init(bh,init,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   logical,            intent(in) :: init
   character(len=*),   intent(in) :: caller

   if(.not. init) then
      call elsi_stop(bh,"Invalid handle! Not initialized.",caller)
   endif

end subroutine

!>
!! This routine gets the global index from the local index of a 2D block-cyclic
!! distribution.
!!
subroutine elsi_get_gid(myid,n_procs,blk,lid,gid)

   implicit none

   integer(kind=i4), intent(in)  :: myid
   integer(kind=i4), intent(in)  :: n_procs
   integer(kind=i4), intent(in)  :: blk
   integer(kind=i4), intent(in)  :: lid
   integer(kind=i4), intent(out) :: gid

   character(len=40), parameter :: caller = "elsi_get_gid"

   gid = lid+myid*blk+((lid-1)/blk)*blk*(n_procs-1)

end subroutine

!>
!! This routine gets the local index from the global index of a 2D block-cyclic
!! distribution.
!!
subroutine elsi_get_lid(n_procs,blk,gid,lid)

   implicit none

   integer(kind=i4), intent(in)  :: n_procs
   integer(kind=i4), intent(in)  :: blk
   integer(kind=i4), intent(in)  :: gid
   integer(kind=i4), intent(out) :: lid

   character(len=40), parameter :: caller = "elsi_get_lid"

   lid = (gid-1)/(n_procs*blk)*blk+mod((gid-1),blk)+1

end subroutine

!>
!! This routine counts the number of non_zero elements in a matrix.
!!
subroutine elsi_get_nnz_real(def0,mat,n_row,n_col,nnz)

   implicit none

   real(kind=r8),    intent(in)  :: def0
   real(kind=r8),    intent(in)  :: mat(n_row,n_col)
   integer(kind=i4), intent(in)  :: n_row
   integer(kind=i4), intent(in)  :: n_col
   integer(kind=i4), intent(out) :: nnz

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col

   character(len=40), parameter :: caller = "elsi_get_nnz_real"

   nnz = 0

   do i_col = 1,n_col
      do i_row = 1,n_row
         if(abs(mat(i_row,i_col)) > def0) then
            nnz = nnz+1
         endif
      enddo
   enddo

end subroutine

!>
!! This routine counts the number of non_zero elements in a matrix.
!!
subroutine elsi_get_nnz_cmplx(def0,mat,n_row,n_col,nnz)

   implicit none

   real(kind=r8),    intent(in)  :: def0
   complex(kind=r8), intent(in)  :: mat(n_row,n_col)
   integer(kind=i4), intent(in)  :: n_row
   integer(kind=i4), intent(in)  :: n_col
   integer(kind=i4), intent(out) :: nnz

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col

   character(len=40), parameter :: caller = "elsi_get_nnz_cmplx"

   nnz = 0

   do i_col = 1,n_col
      do i_row = 1,n_row
         if(abs(mat(i_row,i_col)) > def0) then
            nnz = nnz+1
         endif
      enddo
   enddo

end subroutine

!>
!! This routine computes the trace of a matrix. The size of the matrix is
!! restricted to be identical to Hamiltonian.
!!
subroutine elsi_trace_mat_real(ph,bh,row_map,col_map,mat,trace)

   implicit none

   type(elsi_param_t), intent(in)  :: ph
   type(elsi_basic_t), intent(in)  :: bh
   integer(kind=i4),   intent(in)  :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)  :: col_map(ph%n_basis)
   real(kind=r8),      intent(in)  :: mat(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out) :: trace

   integer(kind=i4) :: i
   integer(kind=i4) :: ierr
   real(kind=r8)    :: l_trace ! Local result

   character(len=40), parameter :: caller = "elsi_trace_mat_real"

   l_trace = 0.0_r8

   do i = 1,ph%n_basis
      if(row_map(i) > 0 .and. col_map(i) > 0) then
         l_trace = l_trace + mat(row_map(i),col_map(i))
      endif
   enddo

   call MPI_Allreduce(l_trace,trace,1,mpi_real8,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

end subroutine

!>
!! This routine computes the trace of a matrix. The size of the matrix is
!! restricted to be identical to Hamiltonian.
!!
subroutine elsi_trace_mat_cmplx(ph,bh,row_map,col_map,mat,trace)

   implicit none

   type(elsi_param_t), intent(in)  :: ph
   type(elsi_basic_t), intent(in)  :: bh
   integer(kind=i4),   intent(in)  :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)  :: col_map(ph%n_basis)
   complex(kind=r8),   intent(in)  :: mat(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(out) :: trace

   integer(kind=i4) :: i
   integer(kind=i4) :: ierr
   complex(kind=r8) :: l_trace ! Local result

   character(len=40), parameter :: caller = "elsi_trace_mat_cmplx"

   l_trace = 0.0_r8

   do i = 1,ph%n_basis
      if(row_map(i) > 0 .and. col_map(i) > 0) then
         l_trace = l_trace + mat(row_map(i),col_map(i))
      endif
   enddo

   call MPI_Allreduce(l_trace,trace,1,mpi_complex16,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

end subroutine

!>
!! This routine computes the trace of the product of two matrices. The size of
!! the two matrices is restricted to be identical to Hamiltonian.
!!
subroutine elsi_trace_mat_mat_real(bh,mat1,mat2,trace)

   implicit none

   type(elsi_basic_t), intent(in)  :: bh
   real(kind=r8),      intent(in)  :: mat1(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)  :: mat2(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out) :: trace

   real(kind=r8)    :: l_trace ! Local result
   integer(kind=i4) :: ierr

   real(kind=r8), external :: ddot

   character(len=40), parameter :: caller = "elsi_trace_mat_mat_real"

   l_trace = ddot(bh%n_lrow*bh%n_lcol,mat1,1,mat2,1)

   call MPI_Allreduce(l_trace,trace,1,mpi_real8,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

end subroutine

!>
!! This routine computes the trace of the product of two matrices. The size of
!! the two matrices is restricted to be identical to Hamiltonian.
!!
subroutine elsi_trace_mat_mat_cmplx(bh,mat1,mat2,trace)

   implicit none

   type(elsi_basic_t), intent(in)  :: bh
   complex(kind=r8),   intent(in)  :: mat1(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(in)  :: mat2(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(out) :: trace

   complex(kind=r8) :: l_trace ! Local result
   integer(kind=i4) :: ierr

   complex(kind=r8), external :: zdotc

   character(len=40), parameter :: caller = "elsi_trace_mat_mat_cmplx"

   l_trace = zdotc(bh%n_lrow*bh%n_lcol,mat1,1,mat2,1)

   call MPI_Allreduce(l_trace,trace,1,mpi_complex16,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

end subroutine

!>
!! This routine symmetrizes an upper or lower triangular matrix. The size of
!! the matrix should be the same as the Hamiltonian matrix.
!!
subroutine elsi_set_full_mat_real(ph,bh,uplo,row_map,col_map,mat)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(in)    :: bh
   integer(kind=i4),   intent(in)    :: uplo
   integer(kind=i4),   intent(in)    :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)    :: col_map(ph%n_basis)
   real(kind=r8),      intent(inout) :: mat(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i
   integer(kind=i4) :: j

   real(kind=r8), allocatable :: tmp_real(:,:)

   character(len=40), parameter :: caller = "elsi_set_full_mat_real"

   call elsi_allocate(bh,tmp_real,bh%n_lrow,bh%n_lcol+2*bh%blk,"tmp_real",&
           caller)

   call pdtran(ph%n_basis,ph%n_basis,1.0_r8,mat,1,1,bh%desc,0.0_r8,tmp_real,1,&
           1,bh%desc)

   if(uplo == UT_MAT) then ! Upper triangular
      do j = 1,ph%n_basis-1
         if(col_map(j) > 0) then
            do i = j+1,ph%n_basis
               if(row_map(i) > 0) then
                  mat(row_map(i),col_map(j)) = tmp_real(row_map(i),col_map(j))
               endif
            enddo
         endif
      enddo
   elseif(uplo == LT_MAT) then ! Lower triangular
      do j = 2,ph%n_basis
         if(col_map(j) > 0) then
            do i = 1,j-1
               if(row_map(i) > 0) then
                  mat(row_map(i),col_map(j)) = tmp_real(row_map(i),col_map(j))
               endif
            enddo
         endif
      enddo
   endif

   call elsi_deallocate(bh,tmp_real,"tmp_real")

end subroutine

!>
!! This routine symmetrizes an upper or lower triangular matrix. The size of
!! the matrix should be the same as the Hamiltonian matrix.
!!
subroutine elsi_set_full_mat_cmplx(ph,bh,uplo,row_map,col_map,mat)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(in)    :: bh
   integer(kind=i4),   intent(in)    :: uplo
   integer(kind=i4),   intent(in)    :: row_map(ph%n_basis)
   integer(kind=i4),   intent(in)    :: col_map(ph%n_basis)
   complex(kind=r8),   intent(inout) :: mat(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i
   integer(kind=i4) :: j

   complex(kind=r8), allocatable :: tmp_cmplx(:,:)

   character(len=40), parameter :: caller = "elsi_set_full_mat_cmplx"

   call elsi_allocate(bh,tmp_cmplx,bh%n_lrow,bh%n_lcol+2*bh%blk,"tmp_cmplx",&
           caller)

   call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),mat,1,1,bh%desc,&
           (0.0_r8,0.0_r8),tmp_cmplx,1,1,bh%desc)

   if(uplo == UT_MAT) then ! Upper triangular
      do j = 1,ph%n_basis-1
         if(col_map(j) > 0) then
            do i = j+1,ph%n_basis
               if(row_map(i) > 0) then
                  mat(row_map(i),col_map(j)) = tmp_cmplx(row_map(i),col_map(j))
               endif
            enddo
         endif
      enddo
   elseif(uplo == LT_MAT) then ! Lower triangular
      do j = 2,ph%n_basis
         if(col_map(j) > 0) then
            do i = 1,j-1
               if(row_map(i) > 0) then
                  mat(row_map(i),col_map(j)) = tmp_cmplx(row_map(i),col_map(j))
               endif
            enddo
         endif
      enddo
   endif

   call elsi_deallocate(bh,tmp_cmplx,"tmp_cmplx")

   ! Make diagonal real
   do j = 1,ph%n_basis
      if(col_map(j) > 0 .and. row_map(j) > 0) then
         mat(row_map(j),col_map(j)) = real(mat(row_map(j),col_map(j)),kind=r8)
      endif
   enddo

end subroutine

end module ELSI_UTILS
