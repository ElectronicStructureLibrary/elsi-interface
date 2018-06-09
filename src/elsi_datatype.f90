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
   use ELSI_CONSTANTS,     only: STR_LEN,UUID_LEN
   use ELSI_PRECISION,     only: r8,i4
   use ELPA,               only: elpa_t,elpa_autotune_t
   use F_PPEXSI_INTERFACE, only: f_ppexsi_options
   use FORTJSON,           only: fjson_handle
   use MATRIXSWITCH,       only: matrix

   implicit none

   private

   type, public :: elsi_basic_t

      ! IO
      integer(kind=i4)        :: print_info
      integer(kind=i4)        :: print_unit
      integer(kind=i4)        :: print_json
      logical                 :: json_init = .false.
      character(len=STR_LEN)  :: user_tag
      character(len=UUID_LEN) :: uuid
      logical                 :: uuid_ready = .false.

      ! MPI
      integer(kind=i4) :: myid
      integer(kind=i4) :: myid_all
      integer(kind=i4) :: n_procs
      integer(kind=i4) :: n_procs_all
      integer(kind=i4) :: comm
      integer(kind=i4) :: comm_all
      logical          :: mpi_ready = .false.
      logical          :: mpi_all_ready = .false.

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
      logical          :: blacs_ready = .false.

      ! Sparse matrix information (common)
      integer(kind=i4) :: nnz_g     ! Global number of nonzeros
      integer(kind=i4) :: nnz_l_sp  ! Local number of nonzeros
      integer(kind=i4) :: n_lcol_sp ! Local number of columns
      real(kind=r8)    :: def0

      ! Sparse matrix information (1D block)
      integer(kind=i4) :: nnz_l_sp1  ! Local number of nonzeros
      integer(kind=i4) :: n_lcol_sp1 ! Local number of columns
      logical          :: pexsi_csc_ready = .false.

      ! Sparse matrix information (1D block-cyclic)
      integer(kind=i4) :: nnz_l_sp2  ! Local number of nonzeros
      integer(kind=i4) :: n_lcol_sp2 ! Local number of columns
      integer(kind=i4) :: blk_sp2
      logical          :: siesta_csc_ready = .false.

   end type

   type, public :: elsi_param_t

      ! General info
      integer(kind=i4) :: solver
      integer(kind=i4) :: matrix_format
      integer(kind=i4) :: parallel_mode
      integer(kind=i4) :: n_calls

      ! Overlap
      logical          :: ovlp_is_unit
      logical          :: ovlp_is_sing ! Is overlap singular?
      logical          :: check_sing   ! Check overlap singularity?
      real(kind=r8)    :: sing_tol     ! Overlap singularity tolerance
      logical          :: stop_sing    ! Always stop if overlap is singular?
      integer(kind=i4) :: n_good       ! Number of nonsingular basis functions

      ! Physics
      real(kind=r8)    :: n_electrons
      integer(kind=i4) :: n_basis
      integer(kind=i4) :: n_spins
      integer(kind=i4) :: n_kpts
      integer(kind=i4) :: n_states
      integer(kind=i4) :: n_states_solve
      integer(kind=i4) :: i_spin
      integer(kind=i4) :: i_kpt
      real(kind=r8)    :: i_weight
      real(kind=r8)    :: spin_degen
      logical          :: spin_is_set = .false.
      real(kind=r8)    :: ebs ! Band structure energy
      logical          :: edm_ready_real = .false.
      logical          :: edm_ready_cmplx = .false.

      ! Chemical potential
      real(kind=r8)    :: mu
      real(kind=r8)    :: ts ! Entropy
      integer(kind=i4) :: mu_scheme
      real(kind=r8)    :: mu_width
      real(kind=r8)    :: mu_tol
      integer(kind=i4) :: mu_max_steps
      integer(kind=i4) :: mu_mp_order

      ! ELPA
      integer(kind=i4)                :: elpa_solver
      integer(kind=i4)                :: elpa_n_single
      integer(kind=i4)                :: elpa_comm_row
      integer(kind=i4)                :: elpa_comm_col
      logical                         :: elpa_gpu
      logical                         :: elpa_gpu_kernels
      logical                         :: elpa_output
      logical                         :: elpa_started = .false.
      class(elpa_t),          pointer :: elpa_main
      class(elpa_autotune_t), pointer :: elpa_tune

      ! libOMM
      integer(kind=i4) :: omm_n_lrow
      integer(kind=i4) :: omm_n_states ! Number of states used in libOMM
      integer(kind=i4) :: omm_n_elpa   ! Number of ELPA steps
      integer(kind=i4) :: omm_flavor   ! 0 (basic) or 2 (Cholesky)
      integer(kind=i4) :: omm_desc(9)
      real(kind=r8)    :: omm_tol
      logical          :: omm_output
      logical          :: omm_started = .false.

      ! PEXSI
      integer(kind=i4)         :: pexsi_np_per_pole
      integer(kind=i4)         :: pexsi_np_per_point
      integer(kind=i4)         :: pexsi_my_prow
      integer(kind=i4)         :: pexsi_my_pcol
      integer(kind=i4)         :: pexsi_n_prow
      integer(kind=i4)         :: pexsi_n_pcol
      integer(kind=i4)         :: pexsi_my_point
      integer(kind=i4)         :: pexsi_myid_point
      integer(kind=i4)         :: pexsi_comm_in_pole
      integer(kind=i4)         :: pexsi_comm_among_point
      real(kind=r8)            :: pexsi_ne
      logical                  :: pexsi_started = .false.
      integer(kind=c_intptr_t) :: pexsi_plan
      type(f_ppexsi_options)   :: pexsi_options

      ! SIPS
      integer(kind=i4) :: sips_n_elpa
      integer(kind=i4) :: sips_np_per_slice
      integer(kind=i4) :: sips_n_slices
      integer(kind=i4) :: sips_slice_type
      integer(kind=i4) :: sips_first_ev ! Index of first eigenvalue to compute
      real(kind=r8)    :: sips_buffer   ! Buffer for adjusting interval
      real(kind=r8)    :: sips_interval(2)
      real(kind=r8)    :: sips_inertia_tol
      logical          :: sips_do_inertia
      logical          :: sips_started = .false.

      ! DMP
      integer(kind=i4) :: dmp_n_states
      integer(kind=i4) :: dmp_method    ! 0 (TRS) or 1 (canonical)
      integer(kind=i4) :: dmp_max_power ! Maximum number of power iterations
      integer(kind=i4) :: dmp_max_iter  ! Maximum number of purification steps
      real(kind=r8)    :: dmp_ev_ham_max
      real(kind=r8)    :: dmp_ev_ham_min
      real(kind=r8)    :: dmp_tol
      logical          :: dmp_started = .false.

   end type

   type, public :: elsi_handle

      type(elsi_basic_t) :: bh
      type(elsi_param_t) :: ph
      type(fjson_handle) :: jh

      ! Dense
      real(kind=r8),    allocatable :: ham_real_den(:,:)
      complex(kind=r8), allocatable :: ham_cmplx_den(:,:)
      real(kind=r8),    allocatable :: ovlp_real_den(:,:)
      complex(kind=r8), allocatable :: ovlp_cmplx_den(:,:)
      real(kind=r8),    allocatable :: eval(:)
      real(kind=r8),    allocatable :: evec_real(:,:)
      complex(kind=r8), allocatable :: evec_cmplx(:,:)
      real(kind=r8),    allocatable :: dm_real_den(:,:)
      complex(kind=r8), allocatable :: dm_cmplx_den(:,:)

      ! Sparse
      real(kind=r8),    allocatable :: ham_real_csc(:)
      complex(kind=r8), allocatable :: ham_cmplx_csc(:)
      real(kind=r8),    allocatable :: ovlp_real_csc(:)
      complex(kind=r8), allocatable :: ovlp_cmplx_csc(:)
      real(kind=r8),    allocatable :: dm_real_csc(:)
      complex(kind=r8), allocatable :: dm_cmplx_csc(:)
      integer(kind=i4), allocatable :: row_ind_sp1(:)
      integer(kind=i4), allocatable :: col_ptr_sp1(:)
      integer(kind=i4), allocatable :: row_ind_sp2(:)
      integer(kind=i4), allocatable :: col_ptr_sp2(:)

      ! Auxiliary
      real(kind=r8),    allocatable :: ham_real_copy(:,:)
      real(kind=r8),    allocatable :: ovlp_real_copy(:,:)
      complex(kind=r8), allocatable :: ovlp_cmplx_copy(:,:)
      real(kind=r8),    allocatable :: ovlp_real_inv(:,:)
      real(kind=r8),    allocatable :: occ(:,:,:)
      integer(kind=i4), allocatable :: row_map(:)
      integer(kind=i4), allocatable :: col_map(:)
      real(kind=r8),    allocatable :: omm_c_real(:,:)
      complex(kind=r8), allocatable :: omm_c_cmplx(:,:)
      real(kind=r8),    allocatable :: pexsi_ne_vec(:)
      real(kind=r8),    allocatable :: dmp_vec1(:)
      real(kind=r8),    allocatable :: dmp_vec2(:)

      logical :: handle_init = .false.

   end type

   type, public :: elsi_rw_handle

      type(elsi_basic_t) :: bh

      integer(kind=i4) :: rw_task ! 0 (READ_FILE) or 1 (WRITE_FILE)
      integer(kind=i4) :: parallel_mode
      integer(kind=i4) :: matrix_format
      real(kind=r8)    :: n_electrons
      integer(kind=i4) :: n_basis
      integer(kind=i4) :: header_user(8)
      logical          :: handle_init = .false.

   end type

end module ELSI_DATATYPE
