! Copyright (c) 2015-2020, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Contain a collection of utility routines.
!!
module ELSI_UTIL

   use ELSI_CONSTANT, only: UNSET,N_SOLVERS,N_PARALLEL_MODES,N_MATRIX_FORMATS,&
       MULTI_PROC,SINGLE_PROC,BLACS_DENSE,PEXSI_CSC,SIESTA_CSC,GENERIC_COO,&
       AUTO_SOLVER,ELPA_SOLVER,OMM_SOLVER,PEXSI_SOLVER,EIGENEXA_SOLVER,&
       SIPS_SOLVER,NTPOLY_SOLVER,MAGMA_SOLVER,BSEPACK_SOLVER,UT_MAT,GET_DM
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI
   use ELSI_OUTPUT, only: elsi_say,elsi_get_time
   use ELSI_PRECISION, only: i4,r8

   implicit none

   private

   public :: elsi_check
   public :: elsi_check_init
   public :: elsi_check_err
   public :: elsi_reset_param
   public :: elsi_reset_basic
   public :: elsi_get_gid
   public :: elsi_get_lid
   public :: elsi_reduce_energy
   public :: elsi_set_full_mat
   public :: elsi_build_dm_edm
   public :: elsi_gram_schmidt

   interface elsi_set_full_mat
      module procedure elsi_set_full_mat_real
      module procedure elsi_set_full_mat_cmplx
   end interface

   interface elsi_build_dm_edm
      module procedure elsi_build_dm_edm_real
      module procedure elsi_build_dm_edm_cmplx
   end interface

   interface elsi_gram_schmidt
      module procedure elsi_gram_schmidt_real
      module procedure elsi_gram_schmidt_cmplx
   end interface

contains

!>
!! Reset ELSI runtime parameters.
!!
subroutine elsi_reset_param(ph)

   implicit none

   type(elsi_param_t), intent(out) :: ph

   character(len=*), parameter :: caller = "elsi_reset_handle"

   ph%solver = UNSET
   ph%matrix_format = UNSET
   ph%parallel_mode = UNSET
   ph%n_calls = 0
   ph%n_calls_all = 0
   ph%save_ovlp = .false.
   ph%unit_ovlp = .false.
   ph%ill_ovlp = .false.
   ph%ill_check = .true.
   ph%ill_tol = 1.0e-5_r8
   ph%n_good = UNSET
   ph%ovlp_ev_min = 0.0_r8
   ph%ovlp_ev_max = 0.0_r8
   ph%n_electrons = 0.0_r8
   ph%n_basis = UNSET
   ph%n_spins = 1
   ph%n_kpts = 1
   ph%n_states = UNSET
   ph%n_states_solve = UNSET
   ph%i_spin = 1
   ph%i_kpt = 1
   ph%i_wt = 1.0_r8
   ph%spin_degen = 0.0_r8
   ph%spin_is_set = .false.
   ph%ebs = 0.0_r8
   ph%energy_gap = 0.0_r8
   ph%spectrum_width = 1.0e3_r8
   ph%dimensionality = 3
   ph%extrapolation = 0
   ph%edm_ready = .false.
   ph%eval_ready = .false.
   ph%evec_ready = .false.
   ph%occ_ready = .false.
   ph%mu = 0.0_r8
   ph%ts = 0.0_r8
   ph%mu_scheme = 0
   ph%mu_width = 1.0e-2_r8
   ph%mu_tol = 1.0e-13_r8
   ph%mu_max_steps = 100
   ph%mu_mp_order = 1
   ph%sparsity_mask = 2
   ph%first_blacs_to_ntpoly = .true.
   ph%first_blacs_to_pexsi = .true.
   ph%first_generic_to_blacs = .true.
   ph%first_generic_to_ntpoly = .true.
   ph%first_generic_to_pexsi = .true.
   ph%first_siesta_to_blacs = .true.
   ph%first_siesta_to_ntpoly = .true.
   ph%first_siesta_to_pexsi = .true.
   ph%first_sips_to_blacs = .true.
   ph%first_sips_to_ntpoly = .true.
   ph%elpa_comm_row = UNSET
   ph%elpa_comm_col = UNSET
   ph%elpa_solver = UNSET
   ph%elpa_n_single = 0
   ph%elpa_gpu = UNSET
   ph%elpa_autotune = UNSET
   ph%elpa_n_lrow = UNSET
   ph%elpa_n_lcol = UNSET
   ph%elpa_first = .true.
   ph%elpa_started = .false.
   ph%omm_n_lrow = UNSET
   ph%omm_n_elpa = 5
   ph%omm_flavor = 0
   ph%omm_desc(:) = UNSET
   ph%omm_tol = 1.0e-12_r8
   ph%omm_output = .false.
   ph%omm_first = .true.
   ph%omm_started = .false.
   ph%pexsi_np_per_pole = UNSET
   ph%pexsi_np_per_point = UNSET
   ph%pexsi_my_prow = UNSET
   ph%pexsi_my_pcol = UNSET
   ph%pexsi_n_prow = UNSET
   ph%pexsi_n_pcol = UNSET
   ph%pexsi_my_point = UNSET
   ph%pexsi_myid_point = UNSET
   ph%pexsi_comm_intra_pole = UNSET
   ph%pexsi_comm_inter_pole = UNSET
   ph%pexsi_comm_inter_point = UNSET
   ph%pexsi_ne = 0.0_r8
   ph%pexsi_mu_min = 0.0_r8
   ph%pexsi_mu_max = 0.0_r8
   ph%pexsi_first = .true.
   ph%pexsi_started = .false.
   ph%exa_n_lrow = UNSET
   ph%exa_n_lcol = UNSET
   ph%exa_my_prow = UNSET
   ph%exa_my_pcol = UNSET
   ph%exa_n_prow = UNSET
   ph%exa_n_pcol = UNSET
   ph%exa_method = 2
   ph%exa_blk_fwd = 48
   ph%exa_blk_bkwd = 128
   ph%exa_first = .true.
   ph%exa_started = .false.
   ph%sips_n_elpa = 1
   ph%sips_n_slices = UNSET
   ph%sips_buffer = 1.0e-1_r8
   ph%sips_interval(:) = 0.0_r8
   ph%sips_inertia_tol = 1.0e-3_r8
   ph%sips_do_inertia = .true.
   ph%sips_first = .true.
   ph%sips_started = .false.
   ph%nt_n_prow = UNSET
   ph%nt_n_pcol = UNSET
   ph%nt_method = 2
   ph%nt_isr = 5
   ph%nt_max_iter = 100
   ph%nt_tol = 1.0e-8_r8
   ph%nt_filter = 1.0e-15_r8
   ph%nt_output = .false.
   ph%nt_first = .true.
   ph%nt_started = .false.
   ph%magma_solver = 1
   ph%magma_n_gpus = 1
   ph%magma_started = .false.
   ph%bse_n_lrow = UNSET
   ph%bse_n_lcol = UNSET
   ph%bse_desc(:) = UNSET

end subroutine

!>
!! Reset ELSI basic information.
!!
subroutine elsi_reset_basic(bh)

   implicit none

   type(elsi_basic_t), intent(out) :: bh

   character(len=*), parameter :: caller = "elsi_reset_basic"

   bh%print_info = 0
   bh%print_unit = 6
   bh%print_json = 0
   bh%json_init = .false.
   bh%myid = UNSET
   bh%myid_all = UNSET
   bh%n_procs = UNSET
   bh%n_procs_all = UNSET
   bh%comm = UNSET
   bh%comm_all = UNSET
   bh%mpi_ready = .false.
   bh%mpi_all_ready = .false.
   bh%blacs_ctxt = UNSET
   bh%desc(:) = UNSET
   bh%blk = UNSET
   bh%n_prow = UNSET
   bh%n_pcol = UNSET
   bh%my_prow = UNSET
   bh%my_pcol = UNSET
   bh%n_lrow = UNSET
   bh%n_lcol = UNSET
   bh%nnz_l = UNSET
   bh%blacs_ready = .false.
   bh%nnz_g = UNSET
   bh%nnz_l_sp = UNSET
   bh%n_lcol_sp = UNSET
   bh%def0 = 1.0e-15_r8
   bh%nnz_l_sp1 = UNSET
   bh%n_lcol_sp1 = UNSET
   bh%pexsi_csc_ready = .false.
   bh%nnz_l_sp2 = UNSET
   bh%n_lcol_sp2 = UNSET
   bh%blk_sp2 = UNSET
   bh%siesta_csc_ready = .false.
   bh%nnz_l_sp3 = UNSET
   bh%generic_coo_ready = .false.

end subroutine

!>
!! Ensure there are no unsupported or mutually conflicting parameters before
!! running actual calculations.
!!
subroutine elsi_check(ph,bh,caller)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   character(len=*), intent(in) :: caller

   integer(kind=i4) :: solver_enabled
   character(len=200) :: msg

   ! General check of solver, parallel mode, matrix format
   if(ph%solver < 0 .or. ph%solver >= N_SOLVERS) then
      write(msg,"(A)") "Unsupported solver"
      call elsi_stop(bh,msg,caller)
   end if

   if(ph%parallel_mode < 0 .or. ph%parallel_mode >= N_PARALLEL_MODES) then
      write(msg,"(A)") "Unsupported parallel mode"
      call elsi_stop(bh,msg,caller)
   end if

   if(ph%matrix_format < 0 .or. ph%matrix_format >= N_MATRIX_FORMATS) then
      write(msg,"(A)") "Unsupported matrix format"
      call elsi_stop(bh,msg,caller)
   end if

   ! Spin
   if(ph%n_spins < 1 .or. ph%n_spins > 2) then
      write(msg,"(A)") "Number of spin channels should be 1 or 2"
      call elsi_stop(bh,msg,caller)
   end if

   if(ph%n_spins == 2 .and. .not. bh%mpi_all_ready) then
      write(msg,"(A)") "Two spin channels requested, but global MPI"//&
         " communicator not set up"
      call elsi_stop(bh,msg,caller)
   end if

   if(.not. ph%spin_is_set) then
      if(ph%n_spins == 2) then
         ph%spin_degen = 1.0_r8
      else
         ph%spin_degen = 2.0_r8
      end if
   end if

   ! k-point
   if(ph%n_kpts < 1) then
      write(msg,"(A)") "Number of k-points cannot be smaller than 1"
      call elsi_stop(bh,msg,caller)
   end if

   if(ph%n_kpts > 1 .and. .not. bh%mpi_all_ready) then
      write(msg,"(A)") "Multiple k-points requested, but global MPI"//&
         " communicator not set up"
      call elsi_stop(bh,msg,caller)
   end if

   if(.not. bh%mpi_ready) then
      bh%comm = MPI_COMM_SELF
      bh%n_procs = 1
      bh%myid = 0
   end if

   if(.not. bh%mpi_all_ready) then
      bh%comm_all = bh%comm
      bh%n_procs_all = bh%n_procs
      bh%myid_all = bh%myid
   end if

   if(bh%myid_all /= 0) then
      bh%print_info = 0
   end if

   if(ph%parallel_mode == MULTI_PROC) then
      if(.not. bh%mpi_ready) then
         write(msg,"(A)") "MULTI_PROC parallel mode requires MPI"
         call elsi_stop(bh,msg,caller)
      end if
   end if

   select case(ph%matrix_format)
   case(BLACS_DENSE)
      if(.not. bh%blacs_ready .and. ph%parallel_mode /= SINGLE_PROC) then
         write(msg,"(A)") "BLACS_DENSE matrix format requested but not"//&
            " properly set up"
         call elsi_stop(bh,msg,caller)
      end if
   case(SIESTA_CSC)
      if(.not. bh%siesta_csc_ready) then
         write(msg,"(A)") "SIESTA_CSC matrix format requested but not"//&
            " properly set up"
         call elsi_stop(bh,msg,caller)
      end if

      if(bh%blk_sp2 == UNSET) then
         write(msg,"(A)") "SIESTA_CSC matrix format requested but block size"//&
            " not set"
         call elsi_stop(bh,msg,caller)
      end if
   case(PEXSI_CSC)
      if(.not. bh%pexsi_csc_ready) then
         write(msg,"(A)") "PEXSI_CSC matrix format requested but not"//&
            " properly set up"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%solver == PEXSI_SOLVER .and. ph%pexsi_np_per_pole == UNSET) then
         write(msg,"(A)") "PEXSI_CSC matrix format requested but number of"//&
            " MPI tasks per pole not set"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%solver == AUTO_SOLVER) then
         write(msg,"(A)") "Solver automatic selection not implemented for"//&
            " PEXSI_CSC matrix format"
         call elsi_stop(bh,msg,caller)
      end if
   case(GENERIC_COO)
      if(.not. bh%generic_coo_ready) then
         write(msg,"(A)") "GENERIC_COO matrix format requested but not"//&
            " properly set up"
         call elsi_stop(bh,msg,caller)
      end if
   end select

   if(ph%unit_ovlp) then
      ph%save_ovlp = .false.
   end if

   ! Specific check for each solver
   select case(ph%solver)
   case(ELPA_SOLVER)
      if(ph%ill_check .and. ph%elpa_autotune > 0) then
         write(msg,"(A)") "ELPA autotuning not compatible with"//&
            " ill-conditioning check"
         call elsi_stop(bh,msg,caller)
      end if
   case(OMM_SOLVER)
      if(ph%parallel_mode /= MULTI_PROC) then
         write(msg,"(A)") "libOMM requires MULTI_PROC parallel mode"
         call elsi_stop(bh,msg,caller)
      end if
   case(PEXSI_SOLVER)
      call elsi_get_pexsi_enabled(solver_enabled)

      if(solver_enabled == 0) then
         write(msg,"(A)") "PEXSI is not enabled in this ELSI installation"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%parallel_mode /= MULTI_PROC) then
         write(msg,"(A)") "PEXSI requires MULTI_PROC parallel mode"
         call elsi_stop(bh,msg,caller)
      end if

      if(mod(bh%n_procs,ph%pexsi_options%nPoints) /= 0) then
         write(msg,"(A)") "To use PEXSI, number of mu points must be a"//&
            " divisor of total number of MPI tasks"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%pexsi_np_per_pole /= UNSET) then
         if(mod(bh%n_procs,ph%pexsi_np_per_pole*ph%pexsi_options%nPoints) /= 0)&
            then
            write(msg,"(A)") "To use PEXSI, number of MPI tasks per pole"//&
               " times number of mu points must be a divisor of total number"//&
               " of MPI tasks"
            call elsi_stop(bh,msg,caller)
         end if

         if(ph%pexsi_np_per_pole*ph%pexsi_options%numPole&
            *ph%pexsi_options%nPoints < bh%n_procs) then
            write(msg,"(A)") "Number of MPI tasks per pole too small"
            call elsi_stop(bh,msg,caller)
         end if
      end if

      if(ph%pexsi_options%method == 1) then
         if(ph%pexsi_options%numPole > 120&
            .or. mod(ph%pexsi_options%numPole,10) /= 0) then
            write(msg,"(A)") "Number of poles with this PEXSI method must be"//&
               " chosen from 10, 20, 30, ..., 120"
            call elsi_stop(bh,msg,caller)
         end if
      else
         if(ph%pexsi_options%numPole > 40&
            .or. mod(ph%pexsi_options%numPole,5) /= 0) then
            write(msg,"(A)") "Number of poles with this PEXSI method must be"//&
               " chosen from 5, 10, 15, ..., 40"
            call elsi_stop(bh,msg,caller)
         end if
      end if
   case(EIGENEXA_SOLVER)
      call elsi_get_eigenexa_enabled(solver_enabled)

      if(solver_enabled == 0) then
         write(msg,"(A)") "EigenEXA is not enabled in this ELSI installation"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%parallel_mode /= MULTI_PROC) then
         write(msg,"(A)") "EigenExa requires MULTI_PROC parallel mode"
         call elsi_stop(bh,msg,caller)
      end if
   case(SIPS_SOLVER)
      call elsi_get_sips_enabled(solver_enabled)

      if(solver_enabled == 0) then
         write(msg,"(A)") "SLEPc-SIPs is not enabled in this ELSI installation"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%n_basis < bh%n_procs) then
         write(msg,"(A)") "Number of MPI tasks too large"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%parallel_mode /= MULTI_PROC) then
         write(msg,"(A)") "SLEPc-SIPs requires MULTI_PROC parallel mode"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%n_spins > 1) then
         write(msg,"(A)") "Two spin channels not supported with SLEPc-SIPs"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%n_kpts > 1) then
         write(msg,"(A)") "Multiple k-points not supported with SLEPc-SIPs"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%sips_n_slices /= UNSET) then
         if(mod(bh%n_procs,ph%sips_n_slices) /= 0) then
            write(msg,"(A)") "To use SLEPc-SIPs, number of slices must be a"//&
               " divisor of total number of MPI tasks"
            call elsi_stop(bh,msg,caller)
         end if
      end if
   case(NTPOLY_SOLVER)
      if(ph%parallel_mode /= MULTI_PROC) then
         write(msg,"(A)") "NTPoly requires MULTI_PROC parallel mode"
         call elsi_stop(bh,msg,caller)
      end if
   case(MAGMA_SOLVER)
      call elsi_get_magma_enabled(solver_enabled)

      if(solver_enabled == 0) then
         write(msg,"(A)") "MAGMA is not enabled in this ELSI installation"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%parallel_mode /= SINGLE_PROC) then
         write(msg,"(A)") "MAGMA requires SINGLE_PROC parallel mode"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%matrix_format /= BLACS_DENSE) then
         write(msg,"(A)") "MAGMA requires BLACS_DENSE matrix format"
         call elsi_stop(bh,msg,caller)
      end if
   case(BSEPACK_SOLVER)
      call elsi_get_bsepack_enabled(solver_enabled)

      if(solver_enabled == 0) then
         write(msg,"(A)") "BSEPACK is not enabled in this ELSI installation"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%parallel_mode /= MULTI_PROC) then
         write(msg,"(A)") "BSEPACK requires MULTI_PROC parallel mode"
         call elsi_stop(bh,msg,caller)
      end if

      if(ph%matrix_format /= BLACS_DENSE) then
         write(msg,"(A)") "BSEPACK requires BLACS_DENSE matrix format"
         call elsi_stop(bh,msg,caller)
      end if
   end select

end subroutine

!>
!! Check if a handle has been properly initialized.
!!
subroutine elsi_check_init(bh,init,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   logical, intent(in) :: init
   character(len=*), intent(in) :: caller

   character(len=200) :: msg

   if(.not. init) then
      write(msg,"(A)") "Handle not initialized"
      call elsi_stop(bh,msg,caller)
   end if

end subroutine

!>
!! Check if error code is 0.
!!
subroutine elsi_check_err(bh,routine,ierr,caller)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   character(len=*), intent(in) :: routine
   integer(kind=i4), intent(in) :: ierr
   character(len=*), intent(in) :: caller

   character(len=200) :: msg

   if(ierr /= 0) then
      write(msg,"(2A)") trim(routine)," failed"

      call elsi_stop(bh,msg,caller)
   end if

end subroutine

!>
!! Get the global index from the local index of a block-cyclic distribution.
!!
subroutine elsi_get_gid(myid,n_procs,blk,lid,gid)

   implicit none

   integer(kind=i4), intent(in) :: myid
   integer(kind=i4), intent(in) :: n_procs
   integer(kind=i4), intent(in) :: blk
   integer(kind=i4), intent(in) :: lid
   integer(kind=i4), intent(out) :: gid

   character(len=*), parameter :: caller = "elsi_get_gid"

   gid = lid+myid*blk+((lid-1)/blk)*blk*(n_procs-1)

end subroutine

!>
!! Get the local index from the global index of a block-cyclic distribution.
!!
subroutine elsi_get_lid(n_procs,blk,gid,lid)

   implicit none

   integer(kind=i4), intent(in) :: n_procs
   integer(kind=i4), intent(in) :: blk
   integer(kind=i4), intent(in) :: gid
   integer(kind=i4), intent(out) :: lid

   character(len=*), parameter :: caller = "elsi_get_lid"

   lid = (gid-1)/(n_procs*blk)*blk+mod((gid-1),blk)+1

end subroutine

!>
!! Reduce energy over spin channels and k-points.
!!
subroutine elsi_reduce_energy(ph,bh,energy)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(inout) :: energy

   real(kind=r8) :: tmp
   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_reduce_energy"

   energy = energy*ph%i_wt

   if(ph%n_spins*ph%n_kpts > 1) then
      if(bh%myid /= 0) then
         energy = 0.0_r8
      end if

      call MPI_Allreduce(energy,tmp,1,MPI_REAL8,MPI_SUM,bh%comm_all,ierr)

      call elsi_check_err(bh,"MPI_Allreduce",ierr,caller)

      energy = tmp
   end if

end subroutine

!>
!! Symmetrize an upper or lower triangular matrix.
!!
subroutine elsi_set_full_mat_real(ph,bh,uplo,mat)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: uplo
   real(kind=r8), intent(inout) :: mat(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: g_row
   integer(kind=i4) :: g_col

   real(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_set_full_mat_real"

   if(uplo == UT_MAT) then
      do j = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,j,g_col)

         do i = 1,bh%n_lrow
            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i,g_row)

            if(g_row == g_col) then
               mat(i,j) = 0.5_r8*mat(i,j)
            else if(g_row > g_col) then
               mat(i,j) = 0.0_r8
            end if
         end do
      end do
   else
      do j = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,j,g_col)

         do i = 1,bh%n_lrow
            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i,g_row)

            if(g_row == g_col) then
               mat(i,j) = 0.5_r8*mat(i,j)
            else if(g_row < g_col) then
               mat(i,j) = 0.0_r8
            end if
         end do
      end do
   end if

   ! Allocate slightly more to work around crashes in some PBLAS
   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol+2*bh%blk,"tmp",caller)

   call pdtran(ph%n_basis,ph%n_basis,1.0_r8,mat,1,1,bh%desc,0.0_r8,tmp,1,1,&
        bh%desc)

   mat(:,:) = mat+tmp(:,1:bh%n_lcol)

   call elsi_deallocate(bh,tmp,"tmp")

end subroutine

!>
!! Symmetrize an upper or lower triangular matrix.
!!
subroutine elsi_set_full_mat_cmplx(ph,bh,uplo,mat)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: uplo
   complex(kind=r8), intent(inout) :: mat(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: g_row
   integer(kind=i4) :: g_col

   complex(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_set_full_mat_cmplx"

   if(uplo == UT_MAT) then
      do j = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,j,g_col)

         do i = 1,bh%n_lrow
            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i,g_row)

            if(g_row == g_col) then
               mat(i,j) = (0.5_r8,0.0_r8)*mat(i,j)
            else if(g_row > g_col) then
               mat(i,j) = (0.0_r8,0.0_r8)
            end if
         end do
      end do
   else
      do j = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,j,g_col)

         do i = 1,bh%n_lrow
            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i,g_row)

            if(g_row == g_col) then
               mat(i,j) = (0.5_r8,0.0_r8)*mat(i,j)
            else if(g_row < g_col) then
               mat(i,j) = (0.0_r8,0.0_r8)
            end if
         end do
      end do
   end if

   ! Allocate slightly more to work around crashes in some PBLAS
   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol+2*bh%blk,"tmp",caller)

   call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),mat,1,1,bh%desc,&
        (0.0_r8,0.0_r8),tmp,1,1,bh%desc)

   mat(:,:) = mat+tmp(:,1:bh%n_lcol)

   call elsi_deallocate(bh,tmp,"tmp")

   ! Make diagonal real
   do j = 1,bh%n_lcol
      call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,j,g_col)

      do i = 1,bh%n_lrow
         call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i,g_row)

         if(g_row == g_col) then
            mat(i,j) = real(mat(i,j),kind=r8)
         end if
      end do
   end do

end subroutine

!>
!! Construct density matrix or energy-weighted density matrix from eigenvectors.
!! Factor contains occupation numbers (GET_DM) or (minus) occupation numbers
!! multiplied with eigenvalues (GET_EDM).
!!
subroutine elsi_build_dm_edm_real(ph,bh,factor,evec,dm,which)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: factor(ph%n_states)
   real(kind=r8), intent(in) :: evec(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: dm(bh%n_lrow,bh%n_lcol)
   integer(kind=i4), intent(in) :: which

   real(kind=r8) :: alpha
   real(kind=r8) :: dummy(1,1)
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: i
   integer(kind=i4) :: gid
   integer(kind=i4) :: max_state
   integer(kind=i4) :: ierr
   logical :: use_elpa_mult
   character(len=200) :: msg

   real(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_build_dm_edm_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)

   dm(:,:) = 0.0_r8

   if(which == GET_DM) then
      alpha = 1.0_r8
   else
      alpha = -1.0_r8
   end if

   ! Methfessel-Paxton or energy density matrix may have negative factors
   if(any(factor < 0.0_r8)) then
      do i = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i,gid)

         if(gid <= ph%n_states_solve) then
            tmp(:,i) = evec(:,i)*factor(gid)
         end if
      end do

      call pdgemm("N","T",ph%n_basis,ph%n_basis,ph%n_states_solve,alpha,tmp,1,&
           1,bh%desc,evec,1,1,bh%desc,0.0_r8,dm,1,1,bh%desc)
   else if(any(factor > 0.0_r8)) then
      do i = 1,ph%n_states
         if(factor(i) > 0.0_r8) then
            max_state = i
         end if
      end do

      do i = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i,gid)

         if(gid <= max_state) then
            tmp(:,i) = evec(:,i)*sqrt(factor(gid))
         end if
      end do

      if(associated(ph%elpa_aux)) then
         call ph%elpa_aux%set("multiply_at_a",1,ierr)

         ! ELPA routine only faster on GPUs
         if(ierr /= 0 .or. ph%elpa_gpu == 0 .or. ph%solver /= ELPA_SOLVER&
            .or. bh%blk*(max(bh%n_prow,bh%n_pcol)-1) >= max_state) then
            use_elpa_mult = .false.
         else
            use_elpa_mult = .true.
         end if
      else
         use_elpa_mult = .false.
      end if

      if(use_elpa_mult) then
         call pdtran(ph%n_basis,ph%n_basis,1.0_r8,tmp,1,1,bh%desc,0.0_r8,dm,1,&
              1,bh%desc)

         call ph%elpa_aux%hermitian_multiply("N","U",max_state,dm,dummy,&
              bh%n_lrow,bh%n_lcol,tmp,bh%n_lrow,bh%n_lcol,ierr)

         call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

         dm(:,:) = alpha*tmp
      else
         call pdsyrk("U","N",ph%n_basis,max_state,alpha,tmp,1,1,bh%desc,0.0_r8,&
              dm,1,1,bh%desc)
      end if

      if(associated(ph%elpa_aux)) then
         call ph%elpa_aux%set("multiply_at_a",0,ierr)
      end if

      call elsi_set_full_mat(ph,bh,UT_MAT,dm)
   end if

   call elsi_deallocate(bh,tmp,"tmp")

   call elsi_get_time(t1)

   if(which == GET_DM) then
      write(msg,"(A)") "Finished density matrix calculation"
   else
      write(msg,"(A)") "Finished energy density matrix calculation"
   end if

   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Construct density matrix or energy-weighted density matrix from eigenvectors.
!! Factor contains occupation numbers (GET_DM) or (minus) occupation numbers
!! multiplied with eigenvalues (GET_EDM).
!!
subroutine elsi_build_dm_edm_cmplx(ph,bh,factor,evec,dm,which)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: factor(ph%n_states)
   complex(kind=r8), intent(in) :: evec(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: dm(bh%n_lrow,bh%n_lcol)
   integer(kind=i4), intent(in) :: which

   complex(kind=r8) :: alpha
   complex(kind=r8) :: dummy(1,1)
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: i
   integer(kind=i4) :: gid
   integer(kind=i4) :: max_state
   integer(kind=i4) :: ierr
   logical :: use_elpa_mult
   character(len=200) :: msg

   complex(kind=r8), allocatable :: tmp(:,:)

   character(len=*), parameter :: caller = "elsi_build_dm_edm_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp,bh%n_lrow,bh%n_lcol,"tmp",caller)

   dm(:,:) = (0.0_r8,0.0_r8)

   if(which == GET_DM) then
      alpha = (1.0_r8,0.0_r8)
   else
      alpha = (-1.0_r8,0.0_r8)
   end if

   if(any(factor /= 0.0_r8)) then
      do i = 1,ph%n_states
         if(factor(i) /= 0.0_r8) then
            max_state = i
         end if
      end do

      do i = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i,gid)

         if(gid <= max_state) then
            tmp(:,i) = evec(:,i)*sqrt(cmplx(factor(gid),kind=r8))
         end if
      end do

      if(associated(ph%elpa_aux)) then
         call ph%elpa_aux%set("multiply_at_a",1,ierr)

         ! ELPA routine only faster on GPUs
         if(ierr /= 0 .or. ph%elpa_gpu == 0 .or. ph%solver /= ELPA_SOLVER&
            .or. bh%blk*(max(bh%n_prow,bh%n_pcol)-1) >= max_state) then
            use_elpa_mult = .false.
         else
            use_elpa_mult = .true.
         end if
      else
         use_elpa_mult = .false.
      end if

      if(use_elpa_mult) then
         call pztranc(ph%n_basis,ph%n_basis,(1.0_r8,0.0_r8),tmp,1,1,bh%desc,&
              (0.0_r8,0.0_r8),dm,1,1,bh%desc)

         call ph%elpa_aux%hermitian_multiply("N","U",max_state,dm,dummy,&
              bh%n_lrow,bh%n_lcol,tmp,bh%n_lrow,bh%n_lcol,ierr)

         call elsi_check_err(bh,"ELPA matrix multiplication",ierr,caller)

         dm(:,:) = alpha*tmp
      else
         call pzherk("U","N",ph%n_basis,max_state,alpha,tmp,1,1,bh%desc,&
              (0.0_r8,0.0_r8),dm,1,1,bh%desc)
      end if

      if(associated(ph%elpa_aux)) then
         call ph%elpa_aux%set("multiply_at_a",0,ierr)
      end if

      call elsi_set_full_mat(ph,bh,UT_MAT,dm)
   end if

   call elsi_deallocate(bh,tmp,"tmp")

   call elsi_get_time(t1)

   if(which == GET_DM) then
      write(msg,"(A)") "Finished density matrix calculation"
   else
      write(msg,"(A)") "Finished energy density matrix calculation"
   end if

   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Orthonormalize eigenvectors with respect to an overlap matrix using a
!! modified Gram-Schmidt algorithm.
!!
subroutine elsi_gram_schmidt_real(ph,bh,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: norm
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: j
   integer(kind=i4) :: lid
   integer(kind=i4) :: i_done
   integer(kind=i4) :: n_block
   character(len=200) :: msg

   real(kind=r8), allocatable :: tmp1(:)
   real(kind=r8), allocatable :: tmp2(:,:)
   real(kind=r8), allocatable :: tmp3(:,:)

   character(len=*), parameter :: caller = "elsi_gram_schmidt_real"

   call elsi_get_time(t0)

   n_block = bh%n_pcol*bh%blk

   call elsi_allocate(bh,tmp1,bh%n_lrow,"tmp1",caller)
   call elsi_allocate(bh,tmp2,bh%n_lrow,bh%n_lcol,"tmp2",caller)
   call elsi_allocate(bh,tmp3,bh%n_lrow,bh%n_lcol,"tmp3",caller)

   call pdsymm("L","U",ph%n_basis,ph%n_states_solve,1.0_r8,ovlp,1,1,bh%desc,&
        evec,1,1,bh%desc,0.0_r8,tmp2,1,1,bh%desc)

   i_done = 0

   do j = 1,ph%n_states_solve
      if(j > i_done+1) then
         ! Dot product of evec(j) with evec(i_done+1..j-1)
         call pdgemv("T",ph%n_basis,j-1-i_done,1.0_r8,tmp2,1,i_done+1,bh%desc,&
              evec,1,j,bh%desc,1,0.0_r8,tmp1,1,1,bh%desc,1)

         ! Orthogonalize
         call pdgemv("N",ph%n_basis,j-1-i_done,-1.0_r8,evec,1,i_done+1,bh%desc,&
              tmp1,1,1,bh%desc,1,1.0_r8,evec,1,j,bh%desc,1)

         call pdgemv("N",ph%n_basis,j-1-i_done,-1.0_r8,tmp2,1,i_done+1,bh%desc,&
              tmp1,1,1,bh%desc,1,1.0_r8,tmp2,1,j,bh%desc,1)
      end if

      ! Normalize
      norm = 0.0_r8

      call pddot(ph%n_basis,norm,tmp2,1,j,bh%desc,1,evec,1,j,bh%desc,1)

      if(mod((j-1)/bh%blk,bh%n_pcol) == bh%my_pcol) then
         call elsi_get_lid(bh%n_pcol,bh%blk,j,lid)

         evec(:,lid) = evec(:,lid)/sqrt(norm)
         tmp2(:,lid) = tmp2(:,lid)/sqrt(norm)
      end if

      if(j-i_done == n_block .and. j < ph%n_states_solve) then
         ! Dot product of evec(i_done+1..j) with evec(j+1..n_states_solve)
         call pdgemm("T","N",n_block,ph%n_states_solve-j,ph%n_basis,1.0_r8,&
              tmp2,1,i_done+1,bh%desc,evec,1,j+1,bh%desc,0.0_r8,tmp3,1,j+1,&
              bh%desc)

         ! Orthogonalize
         call pdgemm("N","N",ph%n_basis,ph%n_states_solve-j,n_block,-1.0_r8,&
              evec,1,i_done+1,bh%desc,tmp3,1,j+1,bh%desc,1.0_r8,evec,1,j+1,&
              bh%desc)

         i_done = i_done+n_block
      end if
   end do

   call elsi_deallocate(bh,tmp1,"tmp1")
   call elsi_deallocate(bh,tmp2,"tmp2")
   call elsi_deallocate(bh,tmp3,"tmp3")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished Gram-Schmidt orthonormalization"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Orthonormalize eigenvectors with respect to an overlap matrix using a
!! modified Gram-Schmidt algorithm.
!!
subroutine elsi_gram_schmidt_cmplx(ph,bh,ovlp,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(in) :: ovlp(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: evec(bh%n_lrow,bh%n_lcol)

   complex(kind=r8) :: norm
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: j
   integer(kind=i4) :: lid
   integer(kind=i4) :: i_done
   integer(kind=i4) :: n_block
   character(len=200) :: msg

   complex(kind=r8), allocatable :: tmp1(:)
   complex(kind=r8), allocatable :: tmp2(:,:)
   complex(kind=r8), allocatable :: tmp3(:,:)

   character(len=*), parameter :: caller = "elsi_new_ev_cmplx"

   call elsi_get_time(t0)

   n_block = bh%n_pcol*bh%blk

   call elsi_allocate(bh,tmp1,bh%n_lrow,"tmp1",caller)
   call elsi_allocate(bh,tmp2,bh%n_lrow,bh%n_lcol,"tmp2",caller)
   call elsi_allocate(bh,tmp3,bh%n_lrow,bh%n_lcol,"tmp3",caller)

   call pzhemm("L","U",ph%n_basis,ph%n_states_solve,(1.0_r8,0.0_r8),ovlp,1,1,&
        bh%desc,evec,1,1,bh%desc,(0.0_r8,0.0_r8),tmp2,1,1,bh%desc)

   i_done = 0

   do j = 1,ph%n_states_solve
      if(j > i_done+1) then
         ! Dot product of evec(j) with evec(i_done+1..j-1)
         call pzgemv("C",ph%n_basis,j-1-i_done,(1.0_r8,0.0_r8),tmp2,1,i_done+1,&
              bh%desc,evec,1,j,bh%desc,1,(0.0_r8,0.0_r8),tmp1,1,1,bh%desc,1)

         ! Orthogonalize
         call pzgemv("N",ph%n_basis,j-1-i_done,(-1.0_r8,0.0_r8),evec,1,&
              i_done+1,bh%desc,tmp1,1,1,bh%desc,1,(1.0_r8,0.0_r8),evec,1,j,&
              bh%desc,1)

         call pzgemv("N",ph%n_basis,j-1-i_done,(-1.0_r8,0.0_r8),tmp2,1,&
              i_done+1,bh%desc,tmp1,1,1,bh%desc,1,(1.0_r8,0.0_r8),tmp2,1,j,&
              bh%desc,1)
      end if

      ! Normalize
      norm = (0.0_r8,0.0_r8)

      call pzdotc(ph%n_basis,norm,tmp2,1,j,bh%desc,1,evec,1,j,bh%desc,1)

      if(mod((j-1)/bh%blk,bh%n_pcol) == bh%my_pcol) then
         call elsi_get_lid(bh%n_pcol,bh%blk,j,lid)

         evec(:,lid) = evec(:,lid)/sqrt(real(norm,kind=r8))
         tmp2(:,lid) = tmp2(:,lid)/sqrt(real(norm,kind=r8))
      end if

      if(j-i_done == n_block .and. j < ph%n_states_solve) then
         ! Dot product of evec(i_done+1..j) with evec(j+1..n_states_solve)
         call pzgemm("C","N",n_block,ph%n_states_solve-j,ph%n_basis,&
              (1.0_r8,0.0_r8),tmp2,1,i_done+1,bh%desc,evec,1,j+1,bh%desc,&
              (0.0_r8,0.0_r8),tmp3,1,j+1,bh%desc)

         ! Orthogonalize
         call pzgemm("N","N",ph%n_basis,ph%n_states_solve-j,n_block,&
              (-1.0_r8,0.0_r8),evec,1,i_done+1,bh%desc,tmp3,1,j+1,bh%desc,&
              (1.0_r8,0.0_r8),evec,1,j+1,bh%desc)

         i_done = i_done+n_block
      end if
   end do

   call elsi_deallocate(bh,tmp1,"tmp1")
   call elsi_deallocate(bh,tmp2,"tmp2")
   call elsi_deallocate(bh,tmp3,"tmp3")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished Gram-Schmidt orthonormalization"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

end module ELSI_UTIL
