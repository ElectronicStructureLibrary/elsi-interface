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
!! This module contains variables accessible in ELSI and related modules.
!!
module ELSI_DATATYPE

   use, intrinsic :: ISO_C_BINDING
   use ELSI_CONSTANTS, only: FULL_MAT,UNSET
   use ELSI_PRECISION, only: r8,i4
   use FOE_BASE, only: foe_data
   use F_PPEXSI_INTERFACE, only: f_ppexsi_options
   use MATRIXSWITCH, only: matrix
   use SPARSEMATRIX_BASE, only: matrices,sparse_matrix

   implicit none

   logical :: print_info = .false.
   logical :: print_mem  = .false.

   type :: elsi_handle

      ! Pointers used when input format compatible with chosen solver
      real(kind=r8),    pointer :: ham_real(:,:)       ! Real Hamiltonian
      complex(kind=r8), pointer :: ham_complex(:,:)    ! Complex Hamiltonian
      real(kind=r8),    pointer :: ovlp_real(:,:)      ! Real overlap
      complex(kind=r8), pointer :: ovlp_complex(:,:)   ! Complex overlap
      real(kind=r8),    pointer :: eval(:)             ! Eigenvalues
      real(kind=r8),    pointer :: evec_real(:,:)      ! Real eigenvectors
      complex(kind=r8), pointer :: evec_complex(:,:)   ! Complex eigenvectors
      real(kind=r8),    pointer :: dm_real(:,:)        ! Real density matrix
      complex(kind=r8), pointer :: dm_complex(:,:)     ! Complex density matrix
      real(kind=r8),    pointer :: ham_real_ccs(:)     ! Real Hamiltonian
      complex(kind=r8), pointer :: ham_complex_ccs(:)  ! Complex Hamiltonian
      real(kind=r8),    pointer :: ovlp_real_ccs(:)    ! Real overlap
      complex(kind=r8), pointer :: ovlp_complex_ccs(:) ! Complex overlap
      real(kind=r8),    pointer :: dm_real_ccs(:)      ! Real density matrix
      complex(kind=r8), pointer :: dm_complex_ccs(:)   ! Complex density matrix
      integer(kind=i4), pointer :: row_ind_ccs(:)      ! Row index
      integer(kind=i4), pointer :: col_ptr_ccs(:)      ! Column pointer

      ! Allocatables
      ! ELPA
      real(kind=r8),    allocatable :: ham_real_elpa(:,:)     ! Real Hamiltonian
      complex(kind=r8), allocatable :: ham_complex_elpa(:,:)  ! Complex Hamiltonian
      real(kind=r8),    allocatable :: ovlp_real_elpa(:,:)    ! Real overlap
      complex(kind=r8), allocatable :: ovlp_complex_elpa(:,:) ! Complex overlap
      real(kind=r8),    allocatable :: eval_elpa(:)           ! Eigenvalues
      real(kind=r8),    allocatable :: evec_real_elpa(:,:)    ! Real eigenvectors
      complex(kind=r8), allocatable :: evec_complex_elpa(:,:) ! Complex eigenvectors
      real(kind=r8),    allocatable :: dm_real_elpa(:,:)      ! Real density matrix
      complex(kind=r8), allocatable :: dm_complex_elpa(:,:)   ! Complex density matrix
      real(kind=r8),    allocatable :: occ_num(:,:,:)         ! Occupation numbers
      real(kind=r8),    allocatable :: eval_all(:,:,:)        ! Eigenvalues
      real(kind=r8),    allocatable :: k_weight(:)            ! K-point weights
      integer(kind=i4), allocatable :: local_row(:)           ! Local-global index mapping
      integer(kind=i4), allocatable :: local_col(:)           ! Local-global index mapping

      ! libOMM
      type(Matrix)                  :: ham_omm               ! Hamiltonian
      type(Matrix)                  :: ovlp_omm              ! Overlap
      type(Matrix)                  :: coeff_omm             ! Coefficient matrix
      type(Matrix)                  :: dm_omm                ! Density matrix
      type(Matrix)                  :: tdm_omm               ! Kinetic energy density matrix
      real(kind=r8),    allocatable :: ovlp_real_omm(:,:)    ! Copy of real overlap
      complex(kind=r8), allocatable :: ovlp_complex_omm(:,:) ! Copy of complex overlap

      ! PESXI
      real(kind=r8),    allocatable :: ham_real_pexsi(:)     ! Sparse real Hamiltonian
      complex(kind=r8), allocatable :: ham_complex_pexsi(:)  ! Sparse complex Hamiltonian
      real(kind=r8),    allocatable :: ovlp_real_pexsi(:)    ! Sparse real overlap
      complex(kind=r8), allocatable :: ovlp_complex_pexsi(:) ! Sparse complex overlap
      real(kind=r8),    allocatable :: dm_real_pexsi(:)      ! Sparse real density matrix
      complex(kind=r8), allocatable :: dm_complex_pexsi(:)   ! Sparse complex density matrix
      integer(kind=i4), allocatable :: row_ind_pexsi(:)      ! Row index
      integer(kind=i4), allocatable :: col_ptr_pexsi(:)      ! Column pointer
      real(kind=r8),    allocatable :: ne_vec(:)             ! Number of electrons at all points

      ! CheSS
      real(kind=r8),    allocatable :: ham_real_chess(:)  ! Sparse real Hamiltonian
      real(kind=r8),    allocatable :: ovlp_real_chess(:) ! Sparse real overlap
      integer(kind=i4), allocatable :: row_ind_chess(:)   ! Row index
      integer(kind=i4), allocatable :: col_ptr_chess(:)   ! Column pointer
      integer(kind=i4), allocatable :: row_ind_buffer(:)  ! Row index of sparsity buffer
      integer(kind=i4), allocatable :: col_ptr_buffer(:)  ! Column pointer of sparsity buffer
      type(matrices)                :: ham_chess          ! Hamiltonian
      type(matrices)                :: ovlp_chess         ! Overlap
      type(matrices)                :: dm_chess           ! Density matrix
      type(matrices)                :: edm_chess          ! Energy density matrix
      type(matrices)                :: ovlp_inv_sqrt(1)   ! (1/overlap)^(-1/2)
      type(sparse_matrix)           :: sparse_mat(2)      ! Sparsity pattern etc

      ! SIPs
      real(kind=r8),    allocatable :: ham_real_sips(:)  ! Sparse real Hamiltonian
      real(kind=r8),    allocatable :: ovlp_real_sips(:) ! Sparse real overlap
      integer(kind=i4), allocatable :: row_ind_sips(:)   ! Row index
      integer(kind=i4), allocatable :: col_ptr_sips(:)   ! Column pointer
      real(kind=r8),    allocatable :: slices(:)         ! Slices

      ! Is this a valid handle?
      logical :: handle_ready = .false.

      ! Solver (AUTO=0,ELPA=1,LIBOMM=2,PEXSI=3,CHESS=4,SIPS=5)
      integer(kind=i4) :: solver = UNSET

      ! Real or complex data (REAL_VALUES=0,COMPLEX_VALUES=1)
      integer(kind=i4) :: matrix_data_type = UNSET

      ! Matrix storage format (BLACS_DENSE=0,PEXSI_CSC=1)
      integer(kind=i4) :: matrix_format = UNSET

      ! Is input matrix triangular? (FULL_MAT=0,UT_MAT=1,LT_MAT=2)
      integer(kind=i4) :: uplo = FULL_MAT

      ! Parallel mode (SINGLE_PROC=0,MULTI_PROC=1)
      integer(kind=i4) :: parallel_mode = UNSET

      ! Number of ELSI being called
      integer(kind=i4) :: n_elsi_calls = 0

      ! Block size in case of block-cyclic distribution
      integer(kind=i4) :: n_b_rows = UNSET
      integer(kind=i4) :: n_b_cols = UNSET

      ! Processor grid
      integer(kind=i4) :: n_p_rows = UNSET
      integer(kind=i4) :: n_p_cols = UNSET

      ! Local matrix size
      integer(kind=i4) :: n_l_rows = UNSET
      integer(kind=i4) :: n_l_cols = UNSET

      ! MPI
      integer(kind=i4) :: myid = UNSET
      integer(kind=i4) :: myid_all = UNSET
      integer(kind=i4) :: n_procs = UNSET
      integer(kind=i4) :: n_procs_all = UNSET
      integer(kind=i4) :: mpi_comm = UNSET
      integer(kind=i4) :: mpi_comm_all = UNSET
      integer(kind=i4) :: mpi_comm_row = UNSET
      integer(kind=i4) :: mpi_comm_col = UNSET
      integer(kind=i4) :: my_p_row = UNSET
      integer(kind=i4) :: my_p_col = UNSET
      logical :: mpi_ready = .false.
      logical :: global_mpi_ready = .false.

      ! BLACS
      integer(kind=i4) :: blacs_ctxt = UNSET
      integer(kind=i4) :: sc_desc(9) = UNSET
      logical :: blacs_ready = .false.

      ! Sparse matrix information
      integer(kind=i4) :: nnz_g = UNSET            ! Global number of nonzeros
      integer(kind=i4) :: nnz_l = UNSET            ! Local number of nonzeros
      integer(kind=i4) :: n_l_cols_sp = UNSET      ! Local number of columns
      integer(kind=i4) :: nnz_l_sp = UNSET         ! Local number of nonzeros
      real(kind=r8) :: zero_threshold = 1.0e-15_r8 ! Threshold to define numerical zero
      logical :: sparsity_ready = .false.          ! Is sparsity pattern set by user?

      ! Overlap
      logical :: ovlp_is_unit = .false.     ! Is overlap matrix unit?
      logical :: ovlp_is_sing = .false.     ! Is overlap matrix singular?
      logical :: no_sing_check = .false.    ! Disable checking for singular overlap?
      real(kind=r8) :: sing_tol = 1.0e-5_r8 ! Eigenfunctions of overlap with eigenvalues
                                            ! smaller than this value will be removed
      logical :: stop_sing = .false.        ! Always stop if overlap is singular?
      integer(kind=i4) :: n_nonsing = UNSET ! Number of nonsingular basis functions

      ! Physics
      real(kind=r8) :: n_electrons = 0.0_r8         ! Number of electrons
      real(kind=r8) :: mu = 0.0_r8                  ! Chemical potential
      integer(kind=i4) :: n_basis = UNSET           ! Number of basis functions
      integer(kind=i4) :: n_spins = 1               ! Number of spin channels
      integer(kind=i4) :: n_kpts = 1                ! Number of k-points
      integer(kind=i4) :: n_states = UNSET          ! Number of states
      integer(kind=i4) :: n_occupied_states = UNSET ! Number of occupied states
      integer(kind=i4) :: i_spin = 1
      integer(kind=i4) :: i_kpt = 1
      real(kind=r8) :: i_weight = 1.0_r8
      real(kind=r8) :: energy_hdm = 0.0_r8
      real(kind=r8) :: energy_sedm = 0.0_r8
      real(kind=r8) :: free_energy = 0.0_r8

      ! Chemical potential
      integer(kind=i4) :: broadening_scheme = 0
      real(kind=r8) :: broadening_width = 1.0e-2_r8
      real(kind=r8) :: occ_tolerance = 1.0e-13_r8   ! Accuracy of occupation numbers
      integer(kind=i4) :: max_mu_steps = 100        ! Maximum number of steps to find mu
      real(kind=r8) :: spin_degen = 0.0_r8          ! Spin degeneracy
      logical :: spin_is_set = .false.
      logical :: mu_ready = .false.
      logical :: edm_ready_real = .false.
      logical :: edm_ready_complex = .false.

      ! ELPA
      integer(kind=i4) :: elpa_solver = UNSET ! ELPA 1-stage or 2-stage solver
      logical :: elpa_output = .false.        ! Output level

      ! libOMM
      integer(kind=i4) :: n_states_omm = UNSET ! Number of states used in libOMM
      integer(kind=i4) :: n_elpa_steps = UNSET ! Use ELPA eigenvectors as initial guess
      logical :: new_overlap = .true.          ! Is a new overlap matrix provided?
      logical :: coeff_ready = .false.         ! Is coefficient matrix initialized?
      integer(kind=i4) :: omm_flavor = UNSET   ! 0 = Basic
                                               ! 1 = Cholesky factorisation (not supported)
                                               ! 2 = Cholesky already performed
                                               ! 3 = Preconditioning (not supported)
      real(kind=r8) :: scale_kinetic = 0.0_r8  ! Factor to scale kinetic energy matrix
      logical :: calc_ed = .false.             ! Calculate energy weighted density matrix?
      real(kind=r8) :: eta = 0.0_r8            ! Eigenspectrum shift parameter
      real(kind=r8) :: min_tol = 0.0_r8        ! Tolerance for minimization
      logical :: omm_output = .false.
      logical :: do_dealloc = .false.
      logical :: use_psp = .false.             ! Use PSP

      ! PEXSI
      integer(kind=i4) :: n_p_per_pole = UNSET
      integer(kind=i4) :: n_p_per_point = UNSET
      integer(kind=i4) :: my_p_row_pexsi = UNSET
      integer(kind=i4) :: my_p_col_pexsi = UNSET
      integer(kind=i4) :: n_p_rows_pexsi = UNSET
      integer(kind=i4) :: n_p_cols_pexsi = UNSET
      integer(kind=i4) :: my_point = UNSET
      integer(kind=i4) :: myid_point = UNSET
      integer(kind=i4) :: comm_among_pole = UNSET
      integer(kind=i4) :: comm_in_pole = UNSET
      integer(kind=i4) :: comm_among_point = UNSET
      integer(kind=i4) :: comm_in_point = UNSET
      real(kind=r8)    :: ne_pexsi ! Number of electrons computed by PEXSI
      logical :: pexsi_started = .false.
      integer(kind=c_intptr_t) :: pexsi_plan
      type(f_ppexsi_options)   :: pexsi_options

      ! CheSS
      type(foe_data) :: foe_obj
      type(foe_data) :: ice_obj
      real(kind=r8) :: erf_decay = 0.0_r8     ! Initial guess of error function decay length
      real(kind=r8) :: erf_decay_min = 0.0_r8 ! Lower bound of decay length
      real(kind=r8) :: erf_decay_max = 0.0_r8 ! Upper bound of decay length
      real(kind=r8) :: ev_ham_min = 0.0_r8
      real(kind=r8) :: ev_ham_max = 0.0_r8
      real(kind=r8) :: ev_ovlp_min = 0.0_r8
      real(kind=r8) :: ev_ovlp_max = 0.0_r8
      real(kind=r8) :: beta = 0.0_r8          ! A patameter used to estimate eigenspectrum
      logical :: chess_started = .false.

      ! SIPs
      integer(kind=i4) :: n_p_per_slice = UNSET
      integer(kind=i4) :: n_inertia_steps = UNSET ! Number of inertia counting steps
      integer(kind=i4) :: slicing_method = UNSET  ! 0 = Equally spaced subintervals
                                                  ! 1 = K-means after equally spaced subintervals
                                                  ! 2 = Equally populated subintervals
                                                  ! 3 = K-means after equally populated subintervals
      integer(kind=i4) :: inertia_option = UNSET  ! Extra inertia computations before solve?
                                                  ! 0 = No
                                                  ! 1 = Yes
      integer(kind=i4) :: unbound = UNSET         ! How to bound left side
                                                  ! 0 = Bound interval
                                                  ! 1 = -infinity
      integer(kind=i4) :: n_slices = UNSET        ! Number of slices
      real(kind=r8) :: interval(2) = 0.0_r8       ! Global interval to search eigenvalues
      real(kind=r8) :: slice_buffer = 0.0_r8      ! Small buffer to expand interval
      logical :: sips_started = .false.

      integer(kind=i4) :: clock_rate = UNSET

   end type

end module ELSI_DATATYPE