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
   use f_ppexsi_interface
   use MatrixSwitch

   implicit none

! ========= MATRIX =========

   !> Pointers used when input format compatible with chosen solver
   real*8,     pointer :: ham_real(:,:)       !< Real Hamiltonian
   complex*16, pointer :: ham_complex(:,:)    !< Complex Hamiltonian
   real*8,     pointer :: ovlp_real(:,:)      !< Real overlap
   complex*16, pointer :: ovlp_complex(:,:)   !< Complex overlap
   real*8,     pointer :: eval(:)             !< Eigenvalues
   real*8,     pointer :: evec_real(:,:)      !< Real eigenvectors
   complex*16, pointer :: evec_complex(:,:)   !< Complex eigenvectors
   real*8,     pointer :: den_mat(:,:)        !< Density matrix
   real*8,     pointer :: ham_real_ccs(:)     !< Real Hamiltonian
   complex*16, pointer :: ham_complex_ccs(:)  !< Complex Hamiltonian
   real*8,     pointer :: ovlp_real_ccs(:)    !< Real overlap
   complex*16, pointer :: ovlp_complex_ccs(:) !< Complex overlap
   real*8,     pointer :: den_mat_ccs(:)      !< Density matrix
   integer,    pointer :: row_ind_ccs(:)      !< Row index
   integer,    pointer :: col_ptr_ccs(:)      !< Column pointer

   !> Allocatables used when input format incompatible with chosen solver
   real*8,     allocatable :: ham_real_elpa(:,:)     !< Real Hamiltonian
   complex*16, allocatable :: ham_complex_elpa(:,:)  !< Complex Hamiltonian
   real*8,     allocatable :: ovlp_real_elpa(:,:)    !< Real overlap
   complex*16, allocatable :: ovlp_complex_elpa(:,:) !< Complex overlap
   real*8,     allocatable :: eval_elpa(:)           !< Eigenvalues
   real*8,     allocatable :: evec_real_elpa(:,:)    !< Real eigenvectors
   complex*16, allocatable :: evec_complex_elpa(:,:) !< Complex eigenvectors
   real*8,     allocatable :: den_mat_elpa(:,:)      !< Density matrix
   real*8,     allocatable :: occ_elpa(:)            !< Occupation numbers

   !> libOMM
   type(Matrix) :: ham_omm       !< Hamiltonian
   type(Matrix) :: ovlp_omm      !< Overlap
   type(Matrix) :: coeff_omm     !< Coefficient matrix
   type(Matrix) :: den_mat_omm   !< Density matrix
   type(Matrix) :: t_den_mat_omm !< Kinetic energy density matrix

   !> PESXI
   real*8,     allocatable :: ham_real_pexsi(:)     !< Sparse real Hamiltonian
   complex*16, allocatable :: ham_complex_pexsi(:)  !< Sparse complex Hamiltonian
   real*8,     allocatable :: ovlp_real_pexsi(:)    !< Sparse real overlap
   complex*16, allocatable :: ovlp_complex_pexsi(:) !< Sparse complex overlap
   real*8,     allocatable :: den_mat_pexsi(:)      !< Sparse density matrix
   real*8,     allocatable :: e_den_mat_pexsi(:)    !< Sparse energy density matrix
   real*8,     allocatable :: f_den_mat_pexsi(:)    !< Sparse free energy density matrix
   integer,    allocatable :: row_ind_pexsi(:)      !< Row index
   integer,    allocatable :: col_ptr_pexsi(:)      !< Column pointer

   !> BLACS
   integer, allocatable :: local_row(:)
   integer, allocatable :: local_col(:)

! ========= PARAMETER =========

   !> Solver (AUTO=0,ELPA=1,LIBOMM=2,PEXSI=3,CHESS=4)
   integer :: method = -1

   !> Real or complex data (REAL_VALUES=0,COMPLEX_VALUES=1)
   integer :: mode = -1

   !> Matrix storage format (DENSE=0,CCS=1,CSC=2,CRS=3,CSR=4)
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

   !> PEXSI may use a different processor grid setup
   integer :: my_p_row_pexsi
   integer :: my_p_col_pexsi
   integer :: n_b_rows_pexsi
   integer :: n_b_cols_pexsi
   integer :: n_p_rows_pexsi
   integer :: n_p_cols_pexsi
   integer :: n_l_rows_pexsi
   integer :: n_l_cols_pexsi
   integer :: n_p_per_pole_pexsi

   !> Sparse matrix information
   integer :: nnz_g                  !< Global number of nonzeros
   integer :: nnz_l                  !< Local number of nonzeros
   integer :: nnz_l_pexsi            !< Local number of nonzeros in PEXSI distribution
   real*8  :: zero_threshold = 1d-13 !< Threshold to define numerical zero

   !> Overlap
   logical :: overlap_is_unit = .false.      !< Is overlap matrix unit?
   logical :: overlap_is_singular = .false.  !< Is overlap matrix singular?
   logical :: no_singularity_check = .false. !< Disable checking for singular overlap?
   real*8  :: singularity_tolerance = 1d-5   !< Eigenfunctions of overlap matrix with
                                             !! eigenvalues smaller than this value
                                             !! will be removed to avoid singularity
   logical :: stop_singularity = .false.     !< Always stop if overlap is singular?
   integer :: n_nonsingular                  !< Number of nonsingular basis functions

   !> Physics
   real*8  :: n_electrons         !< Number of electrons in system
   integer :: n_states            !< Number of total states
   integer :: n_occupied_states   !< Number of occupied states
   real*8  :: hartree = 27.211386 !< Hartree to eV (source: Codata 2015)

   !> ELPA
   logical :: elpa_one_always = .false. !< Always use 1-stage solver
   logical :: elpa_two_always = .false. !< Always use 2-stage solver
   integer :: broaden_method = 0        !< Broadening scheme for occupation numbers
   real*8  :: broaden_width = 1d-2      !< Broadening width for occupation numbers
   real*8  :: occ_tolerance = 1d-10     !< Maximum allowed difference between actual number
                                        !! of electrons and the number computed by ELSI
   integer :: max_mu_steps = 100        !< Maximum number of steps to find the chemical potential

   !> libOMM
   integer :: n_elpa_steps         !< Use ELPA eigenvectors as initial guess
   logical :: new_overlap          !< Is a new overlap matrix provided?
   logical :: coeff_initialized    !< Is coefficient matrix initialized?
   real*8  :: total_energy         !< Energy of the system
   integer :: omm_flavour = -1     !< How to perform the calculation
                                   !! 0 = Basic
                                   !! 1 = Cholesky factorisation
                                   !! 2 = Cholesky already performed
                                   !! 3 = Preconditioning
   real*8  :: scale_kinetic        !< Scaling of the kinetic energy matrix
   logical :: calc_ed = .false.    !< Calculate energy weighted density matrix?
   real*8  :: eta                  !< Eigenspectrum shift parameter
   real*8  :: min_tol              !< Tolerance for minimization
   integer :: nk_times_nspin = -1  !< n_k_points * n_spin
   integer :: i_k_spin = -1        !< Combined k_point spin index
   logical :: omm_verbose          !< Output level
   logical :: do_dealloc           !< Deallocate internal storage?
   logical :: use_psp = .false.    !< Use pspBLAS sparse linear algebra?

   !> PEXSI
   logical                :: pexsi_n_poles_set = .false.      !< Is number of PEXSI poles set by user?
   logical                :: sparsity_pattern_ready = .false. !< Is sparsity pattern set by user?
   logical                :: pole_parallelism = .true.        !< Use PEXSI pole parallelism?
   logical                :: small_pexsi_tol = .false.        !< Is user-defined tolerance smaller than default?
   real(c_double)         :: final_pexsi_tol = 1d-2           !< Default final PEXSI tolerance
   integer(c_intptr_t)    :: pexsi_plan
   type(f_ppexsi_options) :: pexsi_options
   integer(c_int)         :: pexsi_info
   integer(c_int)         :: pexsi_output_file_index
   real(c_double)         :: mu_pexsi = 0d0                   !< Chemical potential computed by PEXSI
   real(c_double)         :: n_electrons_pexsi                !< Number of electrons computed by PEXSI
   real(c_double)         :: mu_min_inertia 
   real(c_double)         :: mu_max_inertia 
   integer(c_int)         :: n_total_inertia_iter
   integer(c_int)         :: n_total_pexsi_iter
   real(c_double)         :: e_tot_h
   real(c_double)         :: e_tot_s
   real(c_double)         :: f_tot

! ========= ALIAS =========

   !> Method names
   enum, bind( C )
      ENUMERATOR :: AUTO, ELPA, LIBOMM, PEXSI, CHESS
   end enum

   !> Real or complex data
   enum, bind( C )
      ENUMERATOR :: REAL_VALUES, COMPLEX_VALUES
   end enum

   !> Storage formats
   enum, bind( C )
      ENUMERATOR :: DENSE, CCS, CSC, CRS, CSR
   end enum

   !> Parallel modes
   enum, bind( C )
      ENUMERATOR :: SINGLE_PROC, MULTI_PROC
   end enum

   !> Broadening type (used if ELPA is chosen to compute density matrix)
   enum, bind( C )
      ENUMERATOR :: GAUSSIAN, FERMI, METHFESSEL_PAXTON
   end enum

end module ELSI_DIMENSIONS