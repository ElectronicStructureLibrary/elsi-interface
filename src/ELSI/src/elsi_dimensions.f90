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
   use ELSI_PRECISION, only : dp
   use f_ppexsi_interface
   use MatrixSwitch

   implicit none

! ========= ARRAYS =========

   !> Pointers used when input format compatible with chosen solver
   real(kind=dp),    pointer :: ham_real(:,:)       !< Real Hamiltonian
   complex(kind=dp), pointer :: ham_complex(:,:)    !< Complex Hamiltonian
   real(kind=dp),    pointer :: ovlp_real(:,:)      !< Real overlap
   complex(kind=dp), pointer :: ovlp_complex(:,:)   !< Complex overlap
   real(kind=dp),    pointer :: eval(:)             !< Eigenvalues
   real(kind=dp),    pointer :: evec_real(:,:)      !< Real eigenvectors
   complex(kind=dp), pointer :: evec_complex(:,:)   !< Complex eigenvectors
   real(kind=dp),    pointer :: den_mat(:,:)        !< Density matrix
   real(kind=dp),    pointer :: ham_real_ccs(:)     !< Real Hamiltonian
   complex(kind=dp), pointer :: ham_complex_ccs(:)  !< Complex Hamiltonian
   real(kind=dp),    pointer :: ovlp_real_ccs(:)    !< Real overlap
   complex(kind=dp), pointer :: ovlp_complex_ccs(:) !< Complex overlap
   real(kind=dp),    pointer :: den_mat_ccs(:)      !< Density matrix
   integer,          pointer :: row_ind_ccs(:)      !< Row index
   integer,          pointer :: col_ptr_ccs(:)      !< Column pointer

   !> Allocatables used when input format incompatible with chosen solver
   !> ELPA
   real(kind=dp),    allocatable :: ham_real_elpa(:,:)     !< Real Hamiltonian
   complex(kind=dp), allocatable :: ham_complex_elpa(:,:)  !< Complex Hamiltonian
   real(kind=dp),    allocatable :: ovlp_real_elpa(:,:)    !< Real overlap
   complex(kind=dp), allocatable :: ovlp_complex_elpa(:,:) !< Complex overlap
   real(kind=dp),    allocatable :: eval_elpa(:)           !< Eigenvalues
   real(kind=dp),    allocatable :: evec_real_elpa(:,:)    !< Real eigenvectors
   complex(kind=dp), allocatable :: evec_complex_elpa(:,:) !< Complex eigenvectors
   real(kind=dp),    allocatable :: den_mat_elpa(:,:)      !< Density matrix
   real(kind=dp),    allocatable :: occ_elpa(:)            !< Occupation numbers

   !> libOMM
   type(Matrix)               :: ham_omm            !< Hamiltonian
   type(Matrix)               :: ovlp_omm           !< Overlap
   type(Matrix)               :: coeff_omm          !< Coefficient matrix
   type(Matrix)               :: den_mat_omm        !< Density matrix
   type(Matrix)               :: t_den_mat_omm      !< Kinetic energy density matrix
   real(kind=dp), allocatable :: ovlp_real_omm(:,:) !< Copy of overlap matrix

   !> PESXI
   real(kind=dp),    allocatable :: ham_real_pexsi(:)     !< Sparse real Hamiltonian
   complex(kind=dp), allocatable :: ham_complex_pexsi(:)  !< Sparse complex Hamiltonian
   real(kind=dp),    allocatable :: ovlp_real_pexsi(:)    !< Sparse real overlap
   complex(kind=dp), allocatable :: ovlp_complex_pexsi(:) !< Sparse complex overlap
   real(kind=dp),    allocatable :: den_mat_pexsi(:)      !< Sparse density matrix
   real(kind=dp),    allocatable :: e_den_mat_pexsi(:)    !< Sparse energy density matrix
   real(kind=dp),    allocatable :: f_den_mat_pexsi(:)    !< Sparse free energy density matrix
   integer,          allocatable :: row_ind_pexsi(:)      !< Row index
   integer,          allocatable :: col_ptr_pexsi(:)      !< Column pointer

   !> SIPs
   real(kind=dp), allocatable :: ham_real_sips(:)  !< Sparse real Hamiltonian
   real(kind=dp), allocatable :: ovlp_real_sips(:) !< Sparse real overlap
   integer,       allocatable :: row_ind_sips(:)   !< Row index
   integer,       allocatable :: col_ptr_sips(:)   !< Column pointer
   real(kind=dp), allocatable :: slices(:)         !< Slices
   real(kind=dp), allocatable :: shifts(:)         !< Shifts
   integer,       allocatable :: inertias(:)       !< Inertia count at each shift

   !> BLACS
   integer, allocatable :: local_row(:)
   integer, allocatable :: local_col(:)

! ========= PARAMETERS =========

   !> Solver (AUTO=0,ELPA=1,LIBOMM=2,PEXSI=3,CHESS=4,SIPS=5)
   integer :: method = -1

   !> Real or complex data (REAL_VALUES=0,COMPLEX_VALUES=1)
   integer :: mode = -1

   !> Matrix storage format (BLACS_DENSE=0,PEXSI_CSC=1)
   integer :: storage = -1

   !> Parallel mode (SINGLE_PROC=0,MULTI_PROC=1)
   integer :: parallelism = -1

   !> Print detailed ELSI info?
   logical :: print_info = .false.

   !> Number of ELSI being called
   integer :: n_elsi_calls

   !> Global matrix size
   integer :: n_g_size
  
   !> Block size in case of block-cyclic distribution
   integer :: n_b_rows
   integer :: n_b_cols

   !> Processor grid
   integer :: n_p_rows
   integer :: n_p_cols
  
   !> Local matrix size
   integer :: n_l_rows
   integer :: n_l_cols

   !> MPI
   integer :: myid
   integer :: n_procs
   integer :: mpi_comm_global
   integer :: mpierr
   logical :: mpi_is_setup = .false.

   !> BLACS
   integer :: blacs_ctxt
   integer :: sc_desc(9)
   integer :: mpi_comm_row
   integer :: mpi_comm_col
   integer :: my_p_row
   integer :: my_p_col
   integer :: blacs_info
   logical :: blacs_is_setup = .false.

   !> Sparse matrix information
   integer :: nnz_g                              !< Global number of nonzeros
   integer :: nnz_l                              !< Local number of nonzeros
   real(kind=dp)  :: zero_threshold = 1.0e-15_dp !< Threshold to define numerical zero

   !> Overlap
   logical :: overlap_is_unit = .false.                !< Is overlap matrix unit?
   logical :: overlap_is_singular = .false.            !< Is overlap matrix singular?
   logical :: no_singularity_check = .false.           !< Disable checking for singular overlap?
   real(kind=dp)  :: singularity_tolerance = 1.0e-5_dp !< Eigenfunctions of overlap matrix with
                                                       !! eigenvalues smaller than this value
                                                       !! will be removed to avoid singularity
   logical :: stop_singularity = .false.               !< Always stop if overlap is singular?
   integer :: n_nonsingular                            !< Number of nonsingular basis functions

   !> Physics
   real(kind=dp)  :: n_electrons              !< Number of electrons in system
   integer :: n_states                        !< Number of total states
   integer :: n_occupied_states               !< Number of occupied states
   real(kind=dp)  :: hartree = 27.21138602_dp !< Hartree to eV (source: Codata 2015)

   !> Chemical potential
   integer :: broaden_method = 0                !< Broadening scheme for occupation numbers
   real(kind=dp)  :: broaden_width = 1.0e-2_dp  !< Broadening width for occupation numbers
   real(kind=dp)  :: occ_tolerance = 1.0e-13_dp !< Maximum allowed difference between actual number
                                                !! of electrons and the number computed by ELSI
   integer :: max_mu_steps = 100                !< Maximum number of steps to find the chemical potential
   real(kind=dp)  :: spin_degen = 0.0_dp        !< Spin degeneracy

   !> ELPA
   logical :: elpa_one_always = .false. !< Always use 1-stage solver
   logical :: elpa_two_always = .false. !< Always use 2-stage solver

   !> libOMM
   integer :: n_elpa_steps         !< Use ELPA eigenvectors as initial guess
   logical :: new_overlap          !< Is a new overlap matrix provided?
   logical :: coeff_initialized    !< Is coefficient matrix initialized?
   real(kind=dp)  :: total_energy  !< Energy of the system
   integer :: omm_flavor = -1      !< How to perform the calculation
                                   !! 0 = Basic
                                   !! 1 = Cholesky factorisation
                                   !! 2 = Cholesky already performed
                                   !! 3 = Preconditioning
   real(kind=dp)  :: scale_kinetic !< Scaling of the kinetic energy matrix
   logical :: calc_ed = .false.    !< Calculate energy weighted density matrix?
   real(kind=dp)  :: eta           !< Eigenspectrum shift parameter
   real(kind=dp)  :: min_tol       !< Tolerance for minimization
   integer :: nk_times_nspin = -1  !< n_k_points * n_spin
   integer :: i_k_spin = -1        !< Combined k_point spin index
   logical :: omm_verbose          !< Output level
   logical :: do_dealloc           !< Deallocate internal storage?
   logical :: use_psp = .false.    !< Use pspBLAS sparse linear algebra?

   !> PEXSI
   integer :: my_p_row_pexsi
   integer :: my_p_col_pexsi
   integer :: n_b_rows_pexsi
   integer :: n_b_cols_pexsi
   integer :: n_p_rows_pexsi
   integer :: n_p_cols_pexsi
   integer :: n_l_rows_pexsi
   integer :: n_l_cols_pexsi
   integer :: n_p_per_pole_pexsi               !< Number of processors per pole
   integer :: nnz_l_pexsi                      !< Local number of nonzeros in PEXSI distribution
   logical :: sparsity_pattern_ready = .false. !< Is sparsity pattern set by user?
   logical :: n_p_per_pole_ready = .false.     !< Is number of processors per pole set by user?
   logical :: small_pexsi_tol = .false.        !< Is user-defined tolerance smaller than default?

   real(c_double)         :: final_pexsi_tol = 1.0e-2_dp !< Default final PEXSI tolerance
   integer(c_intptr_t)    :: pexsi_plan
   type(f_ppexsi_options) :: pexsi_options
   integer(c_int)         :: pexsi_info
   integer(c_int)         :: pexsi_output_file_index
   real(c_double)         :: mu_pexsi = 0.0_dp !< Chemical potential computed by PEXSI
   real(c_double)         :: n_electrons_pexsi !< Number of electrons computed by PEXSI
   real(c_double)         :: mu_min_inertia 
   real(c_double)         :: mu_max_inertia 
   integer(c_int)         :: n_total_inertia_iter
   integer(c_int)         :: n_total_pexsi_iter
   real(c_double)         :: e_tot_h
   real(c_double)         :: e_tot_s
   real(c_double)         :: f_tot

   !> SIPs
   integer :: n_b_rows_sips
   integer :: n_b_cols_sips
   integer :: n_l_rows_sips
   integer :: n_l_cols_sips
   integer :: nnz_l_sips         !< Local number of nonzeros in SIPs distribution
   integer :: n_p_per_slice_sips !< Number of processors per slice
   integer :: n_inertia_steps    !< Number of inertia counting steps
   integer :: n_solve_steps      !< Number of solution steps
   integer :: slicing_method     !< Type of slices
                                 !! 0 = Equally spaced subintervals
                                 !! 1 = K-meaans after equally spaced subintervals
                                 !! 2 = Equally populated subintervals
                                 !! 3 = K-means after equally populated subintervals
   integer :: inertia_option     !< Extra inertia computations before solve?
                                 !! 0 = No
                                 !! 1 = Yes
   integer :: unbound            !< How to bound the left side of the interval
                                 !! 0 = Bound interval
                                 !! 1 = -infinity
   integer :: n_slices           !< Number of slices
   real(kind=dp) :: interval(2)  !< Global interval to search eigenvalues
   real(kind=dp) :: slice_buffer !< Small buffer to expand the eigenvalue interval

! ========= ALIAS =========

   !> Method names
   enum, bind( C )
      enumerator :: AUTO, ELPA, LIBOMM, PEXSI, CHESS, SIPS
   end enum

   !> Real or complex data
   enum, bind( C )
      enumerator :: REAL_VALUES, COMPLEX_VALUES
   end enum

   !> Storage formats
   enum, bind( C )
      enumerator :: BLACS_DENSE, PEXSI_CSC
   end enum

   !> Parallel modes
   enum, bind( C )
      enumerator :: SINGLE_PROC, MULTI_PROC
   end enum

   !> Broadening type (used if ELPA is chosen to compute density matrix)
   enum, bind( C )
      enumerator :: GAUSSIAN, FERMI, METHFESSEL_PAXTON_0, METHFESSEL_PAXTON_1
   end enum

end module ELSI_DIMENSIONS
