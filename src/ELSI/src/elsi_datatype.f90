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
   use ELSI_PRECISION,     only: r8,i4
   use FOE_BASE,           only: foe_data
   use F_PPEXSI_INTERFACE, only: f_ppexsi_options
   use MATRIXSWITCH,       only: matrix
   use SPARSEMATRIX_BASE,  only: matrices,sparse_matrix

   implicit none

   private

   type, public :: elsi_handle

      ! Pointers used when input format compatible with chosen solver
      real(kind=r8),    pointer :: ham_real(:,:)
      complex(kind=r8), pointer :: ham_cmplx(:,:)
      real(kind=r8),    pointer :: ovlp_real(:,:)
      complex(kind=r8), pointer :: ovlp_cmplx(:,:)
      real(kind=r8),    pointer :: eval(:)
      real(kind=r8),    pointer :: evec_real(:,:)
      complex(kind=r8), pointer :: evec_cmplx(:,:)
      real(kind=r8),    pointer :: dm_real(:,:)
      complex(kind=r8), pointer :: dm_cmplx(:,:)
      real(kind=r8),    pointer :: ham_real_ccs(:)
      complex(kind=r8), pointer :: ham_cmplx_ccs(:)
      real(kind=r8),    pointer :: ovlp_real_ccs(:)
      complex(kind=r8), pointer :: ovlp_cmplx_ccs(:)
      real(kind=r8),    pointer :: dm_real_ccs(:)
      complex(kind=r8), pointer :: dm_cmplx_ccs(:)
      integer(kind=i4), pointer :: row_ind_ccs(:)
      integer(kind=i4), pointer :: col_ptr_ccs(:)

      ! Allocatables
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
      real(kind=r8),    allocatable :: k_weight(:)     ! K-point weights

      ! libOMM
      type(Matrix)                  :: ham_omm
      type(Matrix)                  :: ovlp_omm
      type(Matrix)                  :: coeff   ! Coefficient matrix
      type(Matrix)                  :: dm_omm
      type(Matrix)                  :: tdm_omm ! Kinetic energy density matrix
      real(kind=r8),    allocatable :: ovlp_real_copy(:,:)
      complex(kind=r8), allocatable :: ovlp_cmplx_copy(:,:)

      ! PESXI
      real(kind=r8),    allocatable :: ham_real_pexsi(:)
      complex(kind=r8), allocatable :: ham_cmplx_pexsi(:)
      real(kind=r8),    allocatable :: ovlp_real_pexsi(:)
      complex(kind=r8), allocatable :: ovlp_cmplx_pexsi(:)
      real(kind=r8),    allocatable :: dm_real_pexsi(:)
      complex(kind=r8), allocatable :: dm_cmplx_pexsi(:)
      integer(kind=i4), allocatable :: row_ind_pexsi(:)
      integer(kind=i4), allocatable :: col_ptr_pexsi(:)
      real(kind=r8),    allocatable :: ne_vec(:) ! Electron count at mu points

      ! CheSS
      real(kind=r8),    allocatable :: ham_real_chess(:)
      real(kind=r8),    allocatable :: ovlp_real_chess(:)
      integer(kind=i4), allocatable :: row_ind_chess(:)
      integer(kind=i4), allocatable :: col_ptr_chess(:)
      integer(kind=i4), allocatable :: row_ind_buf(:)
      integer(kind=i4), allocatable :: col_ptr_buf(:)
      type(matrices)                :: ham_chess
      type(matrices)                :: ovlp_chess
      type(matrices)                :: dm_chess
      type(matrices)                :: edm_chess
      type(matrices)                :: ovlp_inv_sqrt(1) ! ovlp^(-1/2)
      type(sparse_matrix)           :: sparse_mat(2)

      ! SIPs
      real(kind=r8),    allocatable :: ham_real_sips(:)
      complex(kind=r8), allocatable :: ham_cmplx_sips(:)
      real(kind=r8),    allocatable :: ovlp_real_sips(:)
      complex(kind=r8), allocatable :: ovlp_cmplx_sips(:)
      real(kind=r8),    allocatable :: dm_real_sips(:)
      complex(kind=r8), allocatable :: dm_cmplx_sips(:)
      integer(kind=i4), allocatable :: row_ind_sips(:)
      integer(kind=i4), allocatable :: col_ptr_sips(:)
      real(kind=r8),    allocatable :: slices(:)

      ! Local-global index mapping
      integer(kind=i4), allocatable :: loc_row(:)
      integer(kind=i4), allocatable :: loc_col(:)

      ! Is this a valid handle?
      logical          :: handle_ready = .false.

      ! Solver (AUTO=0,ELPA=1,OMM=2,PEXSI=3,CHESS=4,SIPS=5)
      integer(kind=i4) :: solver

      ! Real or complex data (REAL_VALUES=0,COMPLEX_VALUES=1)
      integer(kind=i4) :: data_type

      ! Matrix format (BLACS_DENSE=0,PEXSI_CSC=1)
      integer(kind=i4) :: matrix_format

      ! Is input matrix triangular? (FULL_MAT=0,UT_MAT=1,LT_MAT=2)
      integer(kind=i4) :: uplo

      ! Parallel mode (SINGLE_PROC=0,MULTI_PROC=1)
      integer(kind=i4) :: parallel_mode

      ! Output
      logical          :: print_info
      logical          :: print_mem
      integer(kind=i4) :: print_unit

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
      logical          :: blacs_ready

      ! Sparse matrix information
      integer(kind=i4) :: nnz_g     ! Global number of nonzeros
      integer(kind=i4) :: nnz_l     ! Local number of nonzeros
      integer(kind=i4) :: nnz_l_sp  ! Local number of nonzeros
      integer(kind=i4) :: n_lcol_sp ! Local number of columns
      real(kind=r8)    :: zero_def
      logical          :: sparsity_ready

      ! Overlap
      logical          :: ovlp_is_unit
      logical          :: ovlp_is_sing  ! Is overlap singular?
      logical          :: no_sing_check ! Disable singularity check?
      real(kind=r8)    :: sing_tol      ! Overlap singularity tolerance
      logical          :: stop_sing     ! Always stop if overlap is singular?
      integer(kind=i4) :: n_nonsing     ! Number of nonsingular basis functions

      ! Physics
      real(kind=r8)    :: n_electrons
      real(kind=r8)    :: mu
      integer(kind=i4) :: n_basis
      integer(kind=i4) :: n_spins
      integer(kind=i4) :: n_kpts
      integer(kind=i4) :: n_states
      integer(kind=i4) :: n_states_solve
      integer(kind=i4) :: i_spin
      integer(kind=i4) :: i_kpt
      real(kind=r8)    :: i_weight
      real(kind=r8)    :: energy_hdm
      real(kind=r8)    :: energy_sedm

      ! Chemical potential
      integer(kind=i4) :: broaden_scheme
      real(kind=r8)    :: broaden_width
      real(kind=r8)    :: occ_tolerance
      integer(kind=i4) :: max_mu_steps
      real(kind=r8)    :: spin_degen
      logical          :: spin_is_set
      logical          :: mu_ready
      logical          :: edm_ready_real
      logical          :: edm_ready_cmplx

      ! ELPA
      integer(kind=i4) :: elpa_solver
      integer(kind=i4) :: n_single_steps
      logical          :: elpa_output
      logical          :: elpa_started = .false.

      ! libOMM
      integer(kind=i4) :: n_states_omm  ! Number of states used in libOMM
      integer(kind=i4) :: omm_n_elpa    ! Number of ELPA steps
      logical          :: new_overlap
      logical          :: coeff_ready   ! Is coefficient initialized?
      integer(kind=i4) :: omm_flavor    ! 0 = Basic
                                        ! 2 = Cholesky already performed
      real(kind=r8)    :: scale_kinetic ! Factor to scale kinetic energy matrix
      logical          :: calc_ed       ! Calculate energy density matrix?
      real(kind=r8)    :: eta           ! Eigenspectrum shift parameter
      real(kind=r8)    :: min_tol       ! Tolerance for minimization
      logical          :: omm_output
      logical          :: do_dealloc
      logical          :: use_psp

      ! PEXSI
      integer(kind=i4) :: np_per_pole
      integer(kind=i4) :: np_per_point
      integer(kind=i4) :: my_prow_pexsi
      integer(kind=i4) :: my_pcol_pexsi
      integer(kind=i4) :: n_prow_pexsi
      integer(kind=i4) :: n_pcol_pexsi
      integer(kind=i4) :: my_point
      integer(kind=i4) :: myid_point
      integer(kind=i4) :: comm_among_pole
      integer(kind=i4) :: comm_in_pole
      integer(kind=i4) :: comm_among_point
      integer(kind=i4) :: comm_in_point
      real(kind=r8)    :: ne_pexsi ! Number of electrons computed by PEXSI
      logical          :: pexsi_started = .false.
      integer(kind=c_intptr_t) :: pexsi_plan
      type(f_ppexsi_options)   :: pexsi_options

      ! CheSS
      type(foe_data)   :: foe_obj
      type(foe_data)   :: ice_obj
      real(kind=r8)    :: erf_decay     ! Error function decay length
      real(kind=r8)    :: erf_decay_min ! Lower bound of decay length
      real(kind=r8)    :: erf_decay_max ! Upper bound of decay length
      real(kind=r8)    :: ev_ham_min
      real(kind=r8)    :: ev_ham_max
      real(kind=r8)    :: ev_ovlp_min
      real(kind=r8)    :: ev_ovlp_max
      real(kind=r8)    :: beta          ! A patameter for eigenspectrum estimate
      logical          :: chess_started = .false.

      ! SIPs
      integer(kind=i4) :: sips_n_elpa    ! Number of ELPA steps
      integer(kind=i4) :: np_per_slice
      integer(kind=i4) :: n_inertia_steps
      integer(kind=i4) :: slicing_method ! 0 = Equally spaced
                                         ! 2 = Equally populated
                                         ! 3 = K-means + equally populated
      integer(kind=i4) :: inertia_option ! Do inertia counting?
                                         ! 0 = No
                                         ! 1 = Yes
      integer(kind=i4) :: unbound        ! How to bound left side?
                                         ! 0 = Bound interval
                                         ! 1 = -infinity
      integer(kind=i4) :: n_slices
      real(kind=r8)    :: interval(2)    ! Interval to search eigenvalues
      real(kind=r8)    :: slice_buffer   ! Small buffer to expand interval
      real(kind=r8)    :: ev_min         ! Lower bound of eigenvalue
      real(kind=r8)    :: ev_max         ! Upper bound of eigenvalue
      logical          :: sips_started = .false.

      ! Timer 
      integer(kind=i4) :: clock_rate

   end type

   type, public :: elsi_rw_handle

      ! Is this a valid handle?
      logical          :: handle_ready = .false.

      ! Reading and writing task (READ_FILE=0,WRITE_FILE=1)
      integer(kind=i4) :: rw_task

      ! Parallel mode (SINGLE_PROC=0,MULTI_PROC=1)
      integer(kind=i4) :: parallel_mode

      ! Matrix format (BLACS_DENSE=0,PEXSI_CSC=1)
      integer(kind=i4) :: matrix_format

      ! File format (DENSE_FILE=0,CSC_FILE=1)
      integer(kind=i4) :: file_format

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
      logical          :: sparsity_ready

      ! Physics
      real(kind=r8)    :: n_electrons
      integer(kind=i4) :: n_basis

      ! User header
      integer(kind=i4) :: header_user(8)

   end type

end module ELSI_DATATYPE
