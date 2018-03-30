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
   use ELSI_CONSTANTS,     only: FILENAME_LEN,STR_LEN,UUID_LEN
   use ELSI_PRECISION,     only: r8,i4
   use F_PPEXSI_INTERFACE, only: f_ppexsi_options
   use MATRIXSWITCH,       only: matrix

   implicit none

   private

   type, public :: elsi_io_handle

      logical                       :: handle_init ! Is this a valid handle?
      integer(kind=i4)              :: print_unit  ! Unit to print to
      character(len=FILENAME_LEN)   :: file_name   ! Name of file
      integer(kind=i4)              :: file_format ! Human-readable, JSON, etc.?
      logical                       :: print_info  ! Are we actually printing?
      character(len=:), allocatable :: prefix      ! Prefix for each line
      integer(kind=i4)              :: comma_json  ! Comma placement in JSON

   end type

   type, public :: elsi_timings_handle

      integer(kind=i4)       :: size_timings ! Array dimension
      integer(kind=i4)       :: n_timings
      character(len=STR_LEN) :: user_tag
      character(len=STR_LEN) :: set_label    ! Timing set identifier

      real(kind=r8),          allocatable :: times(:)
      character(len=STR_LEN), allocatable :: elsi_tags(:)
      character(len=STR_LEN), allocatable :: user_tags(:)

   end type

   type, public :: elsi_handle

      ! ELPA
      real(kind=r8),    allocatable :: ham_real_elpa(:,:)
      complex(kind=r8), allocatable :: ham_cmplx_elpa(:,:)
      real(kind=r8),    allocatable :: ovlp_real_elpa(:,:)
      complex(kind=r8), allocatable :: ovlp_cmplx_elpa(:,:)
      real(kind=r8),    allocatable :: eval_elpa(:)
      real(kind=r8),    allocatable :: evec_real_elpa(:,:)
      complex(kind=r8), allocatable :: evec_cmplx_elpa(:,:)
      real(kind=r8),    allocatable :: dm_real_elpa(:,:)
      complex(kind=r8), allocatable :: dm_cmplx_elpa(:,:)
      real(kind=r8),    allocatable :: occ_num(:,:,:)
      real(kind=r8),    allocatable :: eval_all(:,:,:) ! All eigenvalues
      real(kind=r8),    allocatable :: k_weight(:)

      ! libOMM
      type(Matrix)                  :: ham_omm
      type(Matrix)                  :: ovlp_omm
      type(Matrix)                  :: c_omm   ! Coefficient matrix
      type(Matrix)                  :: dm_omm
      type(Matrix)                  :: tdm_omm ! Kinetic energy matrix

      ! PESXI
      real(kind=r8),    allocatable :: ham_real_pexsi(:)
      complex(kind=r8), allocatable :: ham_cmplx_pexsi(:)
      real(kind=r8),    allocatable :: ovlp_real_pexsi(:)
      complex(kind=r8), allocatable :: ovlp_cmplx_pexsi(:)
      real(kind=r8),    allocatable :: dm_real_pexsi(:)
      complex(kind=r8), allocatable :: dm_cmplx_pexsi(:)
      integer(kind=i4), allocatable :: row_ind_pexsi(:)
      integer(kind=i4), allocatable :: col_ptr_pexsi(:)
      real(kind=r8),    allocatable :: ne_vec_pexsi(:)

      ! SIPS
      real(kind=r8),    allocatable :: evec_real_sips(:,:)
      complex(kind=r8), allocatable :: evec_cmplx_sips(:,:)

      ! DMP
      real(kind=r8),    allocatable :: ham_real_dmp(:,:)
      real(kind=r8),    allocatable :: ovlp_real_dmp(:,:)
      real(kind=r8),    allocatable :: dm_real_dmp(:,:)
      real(kind=r8),    allocatable :: ovlp_real_inv_dmp(:,:)
      real(kind=r8),    allocatable :: evec1_dmp(:)
      real(kind=r8),    allocatable :: evec2_dmp(:)

      ! Auxiliary
      real(kind=r8),    allocatable :: ham_real_copy(:,:)
      real(kind=r8),    allocatable :: ovlp_real_copy(:,:)
      complex(kind=r8), allocatable :: ovlp_cmplx_copy(:,:)
      integer(kind=i4), allocatable :: loc_row(:)
      integer(kind=i4), allocatable :: loc_col(:)
      integer(kind=i4), allocatable :: row_ind_sp2(:)
      integer(kind=i4), allocatable :: col_ptr_sp2(:)

      ! Is this a valid handle?
      logical          :: handle_init    = .false.
      ! Is this handle ready to be used?
      logical          :: handle_ready   = .false.

      ! Solver (AUTO=0,ELPA=1,OMM=2,PEXSI=3,CHESS=4,SIPS=5,DMP=6)
      integer(kind=i4) :: solver

      ! data_type has been removed, as it is not a property of the handle
      ! (a given instance of the handle can solve real or complex problems)

      ! Matrix format (BLACS_DENSE=0,PEXSI_CSC=1,SIESTA_CSC=2)
      integer(kind=i4) :: matrix_format

      ! Is input matrix triangular? (FULL_MAT=0,UT_MAT=1,LT_MAT=2)
      integer(kind=i4) :: uplo

      ! Parallel mode (SINGLE_PROC=0,MULTI_PROC=1)
      integer(kind=i4) :: parallel_mode

      ! Output
      logical          :: print_mem

      ! Number of ELSI being called
      integer(kind=i4) :: n_elsi_calls

      ! MPI
      integer(kind=i4) :: myid
      integer(kind=i4) :: myid_all
      integer(kind=i4) :: n_procs
      integer(kind=i4) :: n_procs_all
      integer(kind=i4) :: mpi_comm
      integer(kind=i4) :: mpi_comm_all
      integer(kind=i4) :: mpi_comm_row
      integer(kind=i4) :: mpi_comm_col
      logical          :: mpi_ready
      logical          :: global_mpi_ready

      ! BLACS
      integer(kind=i4) :: blacs_ctxt
      integer(kind=i4) :: sc_desc(9)
      integer(kind=i4) :: blk_row
      integer(kind=i4) :: blk_col
      integer(kind=i4) :: n_prow
      integer(kind=i4) :: n_pcol
      integer(kind=i4) :: my_prow
      integer(kind=i4) :: my_pcol
      integer(kind=i4) :: n_lrow
      integer(kind=i4) :: n_lcol
      integer(kind=i4) :: nnz_l ! Local number of nonzeros
      logical          :: blacs_ready

      ! Sparse matrix information (common)
      integer(kind=i4) :: nnz_g     ! Global number of nonzeros
      integer(kind=i4) :: nnz_l_sp  ! Local number of nonzeros
      integer(kind=i4) :: n_lcol_sp ! Local number of columns
      real(kind=r8)    :: zero_def

      ! Sparse matrix information (1D block)
      integer(kind=i4) :: nnz_l_sp1  ! Local number of nonzeros
      integer(kind=i4) :: n_lcol_sp1 ! Local number of columns
      logical          :: pexsi_csc_ready

      ! Sparse matrix information (1D block-cyclic)
      integer(kind=i4) :: nnz_l_sp2  ! Local number of nonzeros
      integer(kind=i4) :: n_lcol_sp2 ! Local number of columns
      integer(kind=i4) :: blk_sp2
      logical          :: siesta_csc_ready

      ! Overlap
      logical          :: ovlp_is_unit
      logical          :: ovlp_is_sing ! Is overlap singular?
      logical          :: check_sing   ! Check overlap singularity?
      real(kind=r8)    :: sing_tol     ! Overlap singularity tolerance
      logical          :: stop_sing    ! Always stop if overlap is singular?
      integer(kind=i4) :: n_nonsing    ! Number of nonsingular basis functions

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
      real(kind=r8)    :: energy_hdm
      real(kind=r8)    :: energy_sedm
      real(kind=r8)    :: mu
      real(kind=r8)    :: ts ! Entropy

      ! Chemical potential
      integer(kind=i4) :: broaden_scheme
      real(kind=r8)    :: broaden_width
      real(kind=r8)    :: occ_tolerance
      integer(kind=i4) :: max_mu_steps
      integer(kind=i4) :: mp_order
      logical          :: spin_is_set
      logical          :: mu_ready
      logical          :: ts_ready
      logical          :: edm_ready_real
      logical          :: edm_ready_cmplx

      ! ELPA
      integer(kind=i4) :: elpa_solver
      logical          :: elpa_output
      logical          :: elpa_started = .false.

      ! libOMM
      integer(kind=i4) :: omm_n_states ! Number of states used in libOMM
      integer(kind=i4) :: omm_n_elpa   ! Number of ELPA steps
      integer(kind=i4) :: omm_flavor   ! 0 = Basic
                                       ! 2 = Cholesky already performed
      real(kind=r8)    :: omm_ev_shift ! Eigenspectrum shift parameter
      real(kind=r8)    :: omm_tol      ! Tolerance for minimization
      logical          :: omm_output

      ! PEXSI
      integer(kind=i4) :: pexsi_np_per_pole
      integer(kind=i4) :: pexsi_np_per_point
      integer(kind=i4) :: pexsi_my_prow
      integer(kind=i4) :: pexsi_my_pcol
      integer(kind=i4) :: pexsi_n_prow
      integer(kind=i4) :: pexsi_n_pcol
      integer(kind=i4) :: pexsi_my_point
      integer(kind=i4) :: pexsi_myid_point
      integer(kind=i4) :: pexsi_comm_among_pole
      integer(kind=i4) :: pexsi_comm_in_pole
      integer(kind=i4) :: pexsi_comm_among_point
      integer(kind=i4) :: pexsi_comm_in_point
      real(kind=r8)    :: pexsi_ne ! Number of electrons computed by PEXSI
      logical          :: pexsi_started = .false.
      integer(kind=c_intptr_t) :: pexsi_plan
      type(f_ppexsi_options)   :: pexsi_options

      ! SIPS
      integer(kind=i4) :: sips_n_elpa
      integer(kind=i4) :: sips_np_per_slice
      integer(kind=i4) :: sips_n_slices
      integer(kind=i4) :: sips_slice_type
      integer(kind=i4) :: sips_first_ev
      real(kind=r8)    :: sips_buffer     ! Adjust interval
      real(kind=r8)    :: sips_interval(2)
      real(kind=r8)    :: sips_inertia_tol
      logical          :: sips_do_inertia
      logical          :: sips_started = .false.

      ! DMP
      integer(kind=i4) :: dmp_n_states  ! Number of states used in DMP
      integer(kind=i4) :: dmp_method    ! 0 = Trace correcting
                                        ! 1 = Canonical
      integer(kind=i4) :: dmp_max_power ! Maximum number of power iterations
      integer(kind=i4) :: dmp_max_iter  ! Maximum number of purification steps
      real(kind=r8)    :: dmp_ev_ham_max
      real(kind=r8)    :: dmp_ev_ham_min
      real(kind=r8)    :: dmp_tol       ! Tolerance for purification
      real(kind=r8)    :: dmp_ne        ! Number of electrons computed by DMP

      ! ELSI IO files
      type(elsi_io_handle) :: stdio
      type(elsi_io_handle) :: timings_file

      ! Timer and timings
      type(elsi_timings_handle)   :: timings
      logical                     :: output_timings
      integer(kind=i4)            :: timings_unit
      character(len=FILENAME_LEN) :: timings_name

      ! Versioning
      character(len=STR_LEN)  :: caller
      character(len=UUID_LEN) :: uuid
      logical                 :: uuid_exists

   end type

   type, public :: elsi_rw_handle

      ! Is this a valid handle?
      logical          :: handle_init    = .false.

      ! Reading and writing task (READ_FILE=0,WRITE_FILE=1)
      integer(kind=i4) :: rw_task

      ! Parallel mode (SINGLE_PROC=0,MULTI_PROC=1)
      integer(kind=i4) :: parallel_mode

      ! Matrix format (BLACS_DENSE=0,PEXSI_CSC=1,SIESTA_CSC=2)
      integer(kind=i4) :: matrix_format

      ! Output
      logical          :: print_info
      logical          :: print_mem
      integer(kind=i4) :: print_unit

      ! MPI
      integer(kind=i4) :: myid
      integer(kind=i4) :: n_procs
      integer(kind=i4) :: mpi_comm
      logical          :: mpi_ready

      ! BLACS
      integer(kind=i4) :: blacs_ctxt
      integer(kind=i4) :: blk
      integer(kind=i4) :: n_lrow
      integer(kind=i4) :: n_lcol
      logical          :: blacs_ready

      ! Sparse matrix information
      integer(kind=i4) :: nnz_g     ! Global number of nonzeros
      integer(kind=i4) :: nnz_l_sp  ! Local number of nonzeros
      integer(kind=i4) :: n_lcol_sp ! Local number of columns
      real(kind=r8)    :: zero_def

      ! Physics
      real(kind=r8)    :: n_electrons
      integer(kind=i4) :: n_basis

      ! User header
      integer(kind=i4) :: header_user(8)

   end type

end module ELSI_DATATYPE
