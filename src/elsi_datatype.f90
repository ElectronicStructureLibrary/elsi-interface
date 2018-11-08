! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module contains variables accessible in ELSI and related modules.
!!
module ELSI_DATATYPE

   use, intrinsic :: ISO_C_BINDING
   use ELSI_PRECISION, only: r8,i4
   use ELPA, only: elpa_t,elpa_autotune_t
   use F_PPEXSI_INTERFACE, only: f_ppexsi_options
   use FORTJSON, only: fjson_handle
   use MATRIXSWITCH, only: matrix
   use NTPOLY, only: Permutation_t,Matrix_ps,SolverParameters_t,ProcessGrid_t

   implicit none

   private

   type, public :: elsi_basic_t

      ! IO
      integer(kind=i4) :: print_info
      integer(kind=i4) :: print_unit
      integer(kind=i4) :: print_json
      logical :: json_init = .false.
      character(len=20) :: user_tag
      character(len=36) :: uuid
      logical :: uuid_ready = .false.

      ! MPI
      integer(kind=i4) :: myid
      integer(kind=i4) :: myid_all
      integer(kind=i4) :: n_procs
      integer(kind=i4) :: n_procs_all
      integer(kind=i4) :: comm
      integer(kind=i4) :: comm_all
      logical :: mpi_ready = .false.
      logical :: mpi_all_ready = .false.

      ! BLACS
      integer(kind=i4) :: blacs_ctxt
      integer(kind=i4) :: desc(9)
      integer(kind=i4) :: blk
      integer(kind=i4) :: n_prow
      integer(kind=i4) :: n_pcol
      integer(kind=i4) :: my_prow
      integer(kind=i4) :: my_pcol
      integer(kind=i4) :: n_lrow
      integer(kind=i4) :: n_lcol
      integer(kind=i4) :: nnz_l ! Local number of nonzeros
      logical :: blacs_ready = .false.

      ! Sparse matrix information (common)
      integer(kind=i4) :: nnz_g ! Global number of nonzeros
      integer(kind=i4) :: nnz_l_sp ! Local number of nonzeros
      integer(kind=i4) :: n_lcol_sp ! Local number of columns
      real(kind=r8) :: def0 ! Zero threshold

      ! Sparse matrix information (1D block)
      integer(kind=i4) :: nnz_l_sp1 ! Local number of nonzeros
      integer(kind=i4) :: n_lcol_sp1 ! Local number of columns
      logical :: pexsi_csc_ready = .false.

      ! Sparse matrix information (1D block-cyclic)
      integer(kind=i4) :: nnz_l_sp2 ! Local number of nonzeros
      integer(kind=i4) :: n_lcol_sp2 ! Local number of columns
      integer(kind=i4) :: blk_sp2
      logical :: siesta_csc_ready = .false.

   end type

   type, public :: elsi_param_t

      ! General info
      integer(kind=i4) :: solver
      integer(kind=i4) :: matrix_format
      integer(kind=i4) :: parallel_mode
      integer(kind=i4) :: n_calls ! Number of calls in this geometry step
      integer(kind=i4) :: n_calls_all ! Total number of calls

      ! Overlap
      logical :: unit_ovlp
      logical :: ill_ovlp
      logical :: ill_check
      real(kind=r8) :: ill_tol
      logical :: ill_abort ! Always stop if overlap is ill-conditioned?
      integer(kind=i4) :: n_good ! Number of non-ill-conditioned basis functions

      ! Physics
      real(kind=r8) :: n_electrons
      integer(kind=i4) :: n_basis
      integer(kind=i4) :: n_spins
      integer(kind=i4) :: n_kpts
      integer(kind=i4) :: n_states
      integer(kind=i4) :: n_states_solve
      integer(kind=i4) :: i_spin ! Local spin index
      integer(kind=i4) :: i_kpt ! Local k-point index
      real(kind=r8) :: i_wt ! Weight of local k-point
      real(kind=r8) :: spin_degen
      logical :: spin_is_set = .false.
      real(kind=r8) :: ebs ! Band structure energy
      logical :: edm_ready_real = .false.
      logical :: edm_ready_cmplx = .false.

      ! Chemical potential
      real(kind=r8) :: mu ! Fermi level
      real(kind=r8) :: ts ! Entropy
      integer(kind=i4) :: mu_scheme
      real(kind=r8) :: mu_width
      real(kind=r8) :: mu_tol
      integer(kind=i4) :: mu_max_steps
      integer(kind=i4) :: mu_mp_order

      ! Matrix redistribution
      logical :: first_blacs_to_ntpoly
      logical :: first_blacs_to_pexsi
      logical :: first_blacs_to_sips
      logical :: first_siesta_to_blacs
      logical :: first_siesta_to_pexsi
      logical :: first_sips_to_blacs
      logical :: first_sips_to_ntpoly

      ! ELPA
      integer(kind=i4) :: elpa_solver
      integer(kind=i4) :: elpa_n_single
      integer(kind=i4) :: elpa_comm_row
      integer(kind=i4) :: elpa_comm_col
      logical :: elpa_gpu
      logical :: elpa_gpu_kernels
      logical :: elpa_autotune
      logical :: elpa_output
      logical :: elpa_first
      logical :: elpa_started = .false.
      class(elpa_t), pointer :: elpa_aux
      class(elpa_t), pointer :: elpa_solve
      class(elpa_autotune_t), pointer :: elpa_tune

      ! libOMM
      integer(kind=i4) :: omm_n_lrow
      integer(kind=i4) :: omm_n_states
      integer(kind=i4) :: omm_n_elpa
      integer(kind=i4) :: omm_flavor
      integer(kind=i4) :: omm_desc(9)
      real(kind=r8) :: omm_tol
      logical :: omm_output
      logical :: omm_first
      logical :: omm_started = .false.

      ! PEXSI
      integer(kind=i4) :: pexsi_np_per_pole
      integer(kind=i4) :: pexsi_np_per_point
      integer(kind=i4) :: pexsi_my_prow
      integer(kind=i4) :: pexsi_my_pcol
      integer(kind=i4) :: pexsi_n_prow
      integer(kind=i4) :: pexsi_n_pcol
      integer(kind=i4) :: pexsi_my_point
      integer(kind=i4) :: pexsi_myid_point
      integer(kind=i4) :: pexsi_comm_intra_pole
      integer(kind=i4) :: pexsi_comm_inter_pole
      integer(kind=i4) :: pexsi_comm_inter_point
      real(kind=r8) :: pexsi_ne
      logical :: pexsi_first
      logical :: pexsi_started = .false.
      integer(kind=c_intptr_t) :: pexsi_plan
      type(f_ppexsi_options) :: pexsi_options

      ! SLEPc-SIPs
      integer(kind=i4) :: sips_n_elpa
      integer(kind=i4) :: sips_np_per_slice
      integer(kind=i4) :: sips_n_slices
      integer(kind=i4) :: sips_slice_type
      real(kind=r8) :: sips_buffer ! Buffer for adjusting interval
      real(kind=r8) :: sips_interval(2)
      real(kind=r8) :: sips_inertia_tol
      logical :: sips_do_inertia
      logical :: sips_first
      logical :: sips_started = .false.

      ! NTPoly
      integer(kind=i4) :: nt_n_group
      integer(kind=i4) :: nt_method
      integer(kind=i4) :: nt_isr ! Method to find S^(-1/2)
      integer(kind=i4) :: nt_max_iter
      real(kind=r8) :: nt_tol
      real(kind=r8) :: nt_filter
      logical :: nt_output
      logical :: nt_first
      logical :: nt_started = .false.
      type(Matrix_ps) :: nt_ham
      type(Matrix_ps) :: nt_ovlp
      type(Matrix_ps) :: nt_dm
      type(SolverParameters_t) :: nt_options
      type(Permutation_t) :: nt_perm
      type(ProcessGrid_t) :: nt_pgrid

   end type

   type, public :: elsi_handle

      type(elsi_basic_t) :: bh
      type(elsi_param_t) :: ph
      type(fjson_handle) :: jh

      ! Dense
      real(kind=r8), allocatable :: ham_real_den(:,:)
      complex(kind=r8), allocatable :: ham_cmplx_den(:,:)
      real(kind=r8), allocatable :: ovlp_real_den(:,:)
      complex(kind=r8), allocatable :: ovlp_cmplx_den(:,:)
      real(kind=r8), allocatable :: eval(:)
      real(kind=r8), allocatable :: evec_real(:,:)
      complex(kind=r8), allocatable :: evec_cmplx(:,:)
      real(kind=r8), allocatable :: dm_real_den(:,:)
      complex(kind=r8), allocatable :: dm_cmplx_den(:,:)

      ! Sparse
      real(kind=r8), allocatable :: ham_real_csc(:)
      complex(kind=r8), allocatable :: ham_cmplx_csc(:)
      real(kind=r8), allocatable :: ovlp_real_csc(:)
      complex(kind=r8), allocatable :: ovlp_cmplx_csc(:)
      real(kind=r8), allocatable :: dm_real_csc(:)
      complex(kind=r8), allocatable :: dm_cmplx_csc(:)
      integer(kind=i4), allocatable :: row_ind_sp1(:)
      integer(kind=i4), allocatable :: col_ptr_sp1(:)
      integer(kind=i4), allocatable :: row_ind_sp2(:)
      integer(kind=i4), allocatable :: col_ptr_sp2(:)

      ! Auxiliary
      real(kind=r8), allocatable :: ovlp_real_copy(:,:)
      complex(kind=r8), allocatable :: ovlp_cmplx_copy(:,:)
      real(kind=r8), allocatable :: occ(:,:,:)
      integer(kind=i4), allocatable :: row_map(:)
      integer(kind=i4), allocatable :: col_map(:)
      real(kind=r8), allocatable :: omm_c_real(:,:)
      complex(kind=r8), allocatable :: omm_c_cmplx(:,:)
      real(kind=r8), allocatable :: pexsi_ne_vec(:)

      logical :: handle_init = .false.

   end type

   type, public :: elsi_rw_handle

      type(elsi_basic_t) :: bh

      integer(kind=i4) :: rw_task
      integer(kind=i4) :: parallel_mode
      integer(kind=i4) :: matrix_format
      real(kind=r8) :: n_electrons
      integer(kind=i4) :: n_basis
      integer(kind=i4) :: header_user(8)
      logical :: handle_init = .false.

   end type

end module ELSI_DATATYPE
