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
module ELSI_DIMENSIONS

   use iso_c_binding
   use ELSI_CONSTANTS, only: FULL_MAT,UNSET
   use ELSI_PRECISION, only: r8,i4
   use f_ppexsi_interface, only: f_ppexsi_options
   use MatrixSwitch, only: matrix

   implicit none

   logical :: print_info = .false.

   type :: elsi_handle

      !> Pointers used when input format compatible with chosen solver
      real(kind=r8),    pointer :: ham_real(:,:)       !< Real Hamiltonian
      complex(kind=r8), pointer :: ham_complex(:,:)    !< Complex Hamiltonian
      real(kind=r8),    pointer :: ovlp_real(:,:)      !< Real overlap
      complex(kind=r8), pointer :: ovlp_complex(:,:)   !< Complex overlap
      real(kind=r8),    pointer :: eval(:)             !< Eigenvalues
      real(kind=r8),    pointer :: evec_real(:,:)      !< Real eigenvectors
      complex(kind=r8), pointer :: evec_complex(:,:)   !< Complex eigenvectors
      real(kind=r8),    pointer :: den_mat(:,:)        !< Density matrix
      real(kind=r8),    pointer :: ham_real_ccs(:)     !< Real Hamiltonian
      complex(kind=r8), pointer :: ham_complex_ccs(:)  !< Complex Hamiltonian
      real(kind=r8),    pointer :: ovlp_real_ccs(:)    !< Real overlap
      complex(kind=r8), pointer :: ovlp_complex_ccs(:) !< Complex overlap
      real(kind=r8),    pointer :: den_mat_ccs(:)      !< Density matrix
      integer(kind=i4), pointer :: row_ind_ccs(:)      !< Row index
      integer(kind=i4), pointer :: col_ptr_ccs(:)      !< Column pointer

      !> Allocatables used when input format incompatible with chosen solver
      !> ELPA
      real(kind=r8),    allocatable :: ham_real_elpa(:,:)     !< Real Hamiltonian
      complex(kind=r8), allocatable :: ham_complex_elpa(:,:)  !< Complex Hamiltonian
      real(kind=r8),    allocatable :: ovlp_real_elpa(:,:)    !< Real overlap
      complex(kind=r8), allocatable :: ovlp_complex_elpa(:,:) !< Complex overlap
      real(kind=r8),    allocatable :: eval_elpa(:)           !< Eigenvalues
      real(kind=r8),    allocatable :: evec_real_elpa(:,:)    !< Real eigenvectors
      complex(kind=r8), allocatable :: evec_complex_elpa(:,:) !< Complex eigenvectors
      real(kind=r8),    allocatable :: den_mat_elpa(:,:)      !< Density matrix
      real(kind=r8),    allocatable :: occ_elpa(:)            !< Occupation numbers
      integer(kind=i4), allocatable :: local_row(:)
      integer(kind=i4), allocatable :: local_col(:)

      !> libOMM
      type(Matrix)               :: ham_omm            !< Hamiltonian
      type(Matrix)               :: ovlp_omm           !< Overlap
      type(Matrix)               :: coeff_omm          !< Coefficient matrix
      type(Matrix)               :: den_mat_omm        !< Density matrix
      type(Matrix)               :: t_den_mat_omm      !< Kinetic energy density matrix
      real(kind=r8), allocatable :: ovlp_real_omm(:,:) !< Copy of overlap matrix

      !> PESXI
      real(kind=r8),    allocatable :: ham_real_pexsi(:)     !< Sparse real Hamiltonian
      complex(kind=r8), allocatable :: ham_complex_pexsi(:)  !< Sparse complex Hamiltonian
      real(kind=r8),    allocatable :: ovlp_real_pexsi(:)    !< Sparse real overlap
      complex(kind=r8), allocatable :: ovlp_complex_pexsi(:) !< Sparse complex overlap
      real(kind=r8),    allocatable :: den_mat_pexsi(:)      !< Sparse density matrix
      real(kind=r8),    allocatable :: e_den_mat_pexsi(:)    !< Sparse energy density matrix
      real(kind=r8),    allocatable :: f_den_mat_pexsi(:)    !< Sparse free energy density matrix
      integer(kind=i4), allocatable :: row_ind_pexsi(:)      !< Row index
      integer(kind=i4), allocatable :: col_ptr_pexsi(:)      !< Column pointer

      !> SIPs
      real(kind=r8),    allocatable :: ham_real_sips(:)  !< Sparse real Hamiltonian
      real(kind=r8),    allocatable :: ovlp_real_sips(:) !< Sparse real overlap
      integer(kind=i4), allocatable :: row_ind_sips(:)   !< Row index
      integer(kind=i4), allocatable :: col_ptr_sips(:)   !< Column pointer
      real(kind=r8),    allocatable :: slices(:)         !< Slices
      real(kind=r8),    allocatable :: shifts(:)         !< Shifts
      integer(kind=i4), allocatable :: inertias(:)       !< Inertia count at each shift

      !> Is this a valid handle? (!!!)
      logical :: handle_initialized = .false.

      !> Solver (AUTO=0,ELPA=1,LIBOMM=2,PEXSI=3,CHESS=4,SIPS=5)
      integer(kind=i4) :: solver = UNSET

      !> Real or complex data (REAL_VALUES=0,COMPLEX_VALUES=1)
      integer(kind=i4) :: matrix_data_type = UNSET

      !> Matrix storage format (BLACS_DENSE=0,PEXSI_CSC=1)
      integer(kind=i4) :: matrix_storage_format = UNSET

      !> Is input matrix triangular? (FULL_MAT=0,UT_MAT=1,LT_MAT=2)
      integer(kind=i4) :: uplo = FULL_MAT

      !> Parallel mode (SINGLE_PROC=0,MULTI_PROC=1)
      integer(kind=i4) :: parallel_mode = UNSET

      !> Number of ELSI being called
      integer(kind=i4) :: n_elsi_calls = 0

      !> Global matrix size
      integer(kind=i4) :: n_g_size = UNSET
  
      !> Block size in case of block-cyclic distribution
      integer(kind=i4) :: n_b_rows = UNSET
      integer(kind=i4) :: n_b_cols = UNSET

      !> Processor grid
      integer(kind=i4) :: n_p_rows = UNSET
      integer(kind=i4) :: n_p_cols = UNSET
  
      !> Local matrix size
      integer(kind=i4) :: n_l_rows = UNSET
      integer(kind=i4) :: n_l_cols = UNSET

      !> MPI
      integer(kind=i4) :: myid = UNSET
      integer(kind=i4) :: n_procs = UNSET
      integer(kind=i4) :: mpi_comm = UNSET
      logical          :: mpi_is_setup = .false.

      !> BLACS
      integer(kind=i4) :: blacs_ctxt = UNSET
      integer(kind=i4) :: sc_desc(9) = UNSET
      integer(kind=i4) :: mpi_comm_row = UNSET
      integer(kind=i4) :: mpi_comm_col = UNSET
      integer(kind=i4) :: my_p_row = UNSET
      integer(kind=i4) :: my_p_col = UNSET
      logical          :: blacs_is_setup = .false.

      !> Sparse matrix information
      integer(kind=i4) :: nnz_g = UNSET               !< Global number of nonzeros
      integer(kind=i4) :: nnz_l = UNSET               !< Local number of nonzeros
      real(kind=r8)    :: zero_threshold = 1.0e-15_r8 !< Threshold to define numerical zero

      !> Overlap
      logical :: overlap_is_unit = .false.               !< Is overlap matrix unit?
      logical :: overlap_is_singular = .false.           !< Is overlap matrix singular?
      logical :: no_singularity_check = .false.          !< Disable checking for singular overlap?
      real(kind=r8) :: singularity_tolerance = 1.0e-5_r8 !< Eigenfunctions of overlap matrix with
                                                         !! eigenvalues smaller than this value
                                                         !! will be removed to avoid singularity
      logical :: stop_singularity = .false.              !< Always stop if overlap is singular?
      integer(kind=i4) :: n_nonsingular = UNSET          !< Number of nonsingular basis functions

      !> Physics
      real(kind=r8) :: n_electrons = 0.0_r8         !< Number of electrons in system
      integer(kind=i4) :: n_states = UNSET          !< Number of total states
      integer(kind=i4) :: n_occupied_states = UNSET !< Number of occupied states

      !> Chemical potential
      integer(kind=i4) :: broadening_scheme = 0     !< Broadening scheme for occupation numbers
      real(kind=r8) :: broadening_width = 1.0e-2_r8 !< Broadening width for occupation numbers
      real(kind=r8) :: occ_tolerance = 1.0e-13_r8   !< Maximum allowed difference between actual number
                                                    !! of electrons and the number computed by ELSI
      integer(kind=i4) :: max_mu_steps = 100        !< Maximum number of steps to find the chemical potential
      real(kind=r8) :: spin_degen = 0.0_r8          !< Spin degeneracy

      !> ELPA
      logical :: elpa_one_always = .false. !< Always use 1-stage solver
      logical :: elpa_two_always = .false. !< Always use 2-stage solver

      !> libOMM
      integer(kind=i4) :: n_elpa_steps = UNSET   !< Use ELPA eigenvectors as initial guess
      logical :: new_overlap = .true.            !< Is a new overlap matrix provided?
      logical :: coeff_initialized = .false.     !< Is coefficient matrix initialized?
      real(kind=r8) :: total_energy = 0.0_r8     !< Energy of the system
      integer(kind=i4) :: omm_flavor = UNSET     !< How to perform the calculation
                                                 !! 0 = Basic
                                                 !! 1 = Cholesky factorisation (not supported)
                                                 !! 2 = Cholesky already performed
                                                 !! 3 = Preconditioning (not supported)
      real(kind=r8) :: scale_kinetic = 0.0_r8    !< Scaling of the kinetic energy matrix
      logical :: calc_ed = .false.               !< Calculate energy weighted density matrix?
      real(kind=r8) :: eta = 0.0_r8              !< Eigenspectrum shift parameter
      real(kind=r8) :: min_tol = 1.0e-12_r8      !< Tolerance for minimization
      integer(kind=i4) :: nk_times_nspin = UNSET !< n_k_points * n_spin
      integer(kind=i4) :: i_k_spin = UNSET       !< Combined k_point spin index
      logical :: omm_verbose = .false.           !< Output level
      logical :: do_dealloc = .false.            !< Deallocate internal storage?
      logical :: use_psp = .false.               !< Use pspBLAS sparse linear algebra?

      !> PEXSI
      integer(kind=i4) :: my_p_row_pexsi = UNSET
      integer(kind=i4) :: my_p_col_pexsi = UNSET
      integer(kind=i4) :: n_b_rows_pexsi = UNSET
      integer(kind=i4) :: n_b_cols_pexsi = UNSET
      integer(kind=i4) :: n_p_rows_pexsi = UNSET
      integer(kind=i4) :: n_p_cols_pexsi = UNSET
      integer(kind=i4) :: n_l_rows_pexsi = UNSET
      integer(kind=i4) :: n_l_cols_pexsi = UNSET
      integer(kind=i4) :: n_p_per_pole_pexsi = UNSET !< Number of processors per pole
      integer(kind=i4) :: nnz_l_pexsi = UNSET        !< Local number of nonzeros in PEXSI distribution
      logical :: sparsity_pattern_ready = .false.    !< Is sparsity pattern set by user?
      logical :: n_p_per_pole_ready = .false.        !< Is number of processors per pole set by user?
      logical :: small_pexsi_tol = .false.           !< Is user-defined tolerance smaller than default?
      logical :: pexsi_started = .false.             !< Is PEXSI started?

      real(kind=r8)            :: final_pexsi_tol = 1.0e-2_r8
      integer(kind=c_intptr_t) :: pexsi_plan
      type(f_ppexsi_options)   :: pexsi_options
      integer(kind=i4)         :: pexsi_info = UNSET
      integer(kind=i4)         :: pexsi_output_file_index = UNSET
      real(kind=r8)            :: mu_pexsi = 0.0_r8          !< Chemical potential computed by PEXSI
      real(kind=r8)            :: n_electrons_pexsi = 0.0_r8 !< Number of electrons computed by PEXSI
      real(kind=r8)            :: mu_min_inertia = 0.0_r8
      real(kind=r8)            :: mu_max_inertia = 0.0_r8
      integer(kind=i4)         :: n_total_inertia_iter = UNSET
      integer(kind=i4)         :: n_total_pexsi_iter = UNSET
      real(kind=r8)            :: e_tot_h = 0.0_r8
      real(kind=r8)            :: e_tot_s = 0.0_r8
      real(kind=r8)            :: f_tot =0.0_r8

      !> SIPs
      integer(kind=i4) :: n_b_rows_sips = UNSET
      integer(kind=i4) :: n_b_cols_sips = UNSET
      integer(kind=i4) :: n_l_rows_sips = UNSET
      integer(kind=i4) :: n_l_cols_sips = UNSET
      integer(kind=i4) :: nnz_l_sips = UNSET         !< Local number of nonzeros in SIPs distribution
      integer(kind=i4) :: n_p_per_slice_sips = UNSET !< Number of processors per slice
      integer(kind=i4) :: n_inertia_steps = UNSET    !< Number of inertia counting steps
      integer(kind=i4) :: n_solve_steps = UNSET      !< Number of solution steps
      integer(kind=i4) :: slicing_method = UNSET     !< Type of slices
                                                     !! 0 = Equally spaced subintervals
                                                     !! 1 = K-meaans after equally spaced subintervals
                                                     !! 2 = Equally populated subintervals
                                                     !! 3 = K-means after equally populated subintervals
      integer(kind=i4) :: inertia_option = UNSET     !< Extra inertia computations before solve?
                                                     !! 0 = No
                                                     !! 1 = Yes
      integer(kind=i4) :: unbound = UNSET            !< How to bound the left side of the interval
                                                     !! 0 = Bound interval
                                                     !! 1 = -infinity
      integer(kind=i4) :: n_slices = UNSET           !< Number of slices
      real(kind=r8) :: interval(2) = 0.0_r8          !< Global interval to search eigenvalues
      real(kind=r8) :: slice_buffer = 0.0_r8         !< Small buffer to expand the eigenvalue interval
      logical :: sips_started = .false.              !< Is SIPs started?

      !> Timers
      real(kind=r8) :: t_generalized_evp
      real(kind=r8) :: t_generalized_evp_start
      real(kind=r8) :: t_redistribution
      real(kind=r8) :: t_redistribution_start
      real(kind=r8) :: t_transform_evp
      real(kind=r8) :: t_transform_evp_start
      real(kind=r8) :: t_back_transform_ev
      real(kind=r8) :: t_back_transform_ev_start
      real(kind=r8) :: t_singularity_check
      real(kind=r8) :: t_singularity_check_start
      real(kind=r8) :: t_standard_evp
      real(kind=r8) :: t_standard_evp_start
      real(kind=r8) :: t_density_matrix
      real(kind=r8) :: t_density_matrix_start
      real(kind=r8) :: t_cholesky
      real(kind=r8) :: t_cholesky_start

   end type

end module ELSI_DIMENSIONS
