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
!! This module provides routines for setting up and solving or circumventing
!! an Kohn-Sham eigenvalue problem using ELPA, libOMM, PEXSI, CheSS, SIPs.
!!
module ELSI

   use iso_c_binding
   use ELSI_CONSTANTS, only: ELPA,LIBOMM,PEXSI,CHESS,SIPS,REAL_VALUES,&
                             COMPLEX_VALUES,SINGLE_PROC,MULTI_PROC,UNSET
   use ELSI_DIMENSIONS, only: elsi_handle,print_info
   use ELSI_ELPA
   use ELSI_MATCONV, only: elsi_blacs_to_pexsi,elsi_blacs_to_sips,&
                           elsi_pexsi_to_blacs_dm,elsi_pexsi_to_blacs
   use ELSI_MU, only: elsi_compute_mu_and_occ
   use ELSI_OMM
   use ELSI_PEXSI
   use ELSI_PRECISION, only: r8,i4
   use ELSI_SIPS
   use ELSI_TIMERS
   use ELSI_UTILS
   use MatrixSwitch, only: m_allocate,ms_scalapack_setup

   implicit none
   private

   public :: elsi_handle             !< ELSI handle data type
   public :: elsi_init               !< Initialize
   public :: elsi_set_solver         !< Set solver from calling code
   public :: elsi_set_mpi            !< Set MPI from calling code
   public :: elsi_set_blacs          !< Set BLACS from calling code
   public :: elsi_set_csc            !< Set CSC sparsity pattern from calling code
   public :: elsi_customize          !< Override ELSI default
   public :: elsi_customize_elpa     !< Override ELPA default
   public :: elsi_customize_omm      !< Override libOMM default
   public :: elsi_customize_pexsi    !< Override PEXSI default
   public :: elsi_customize_mu       !< Override chemical potential determination
   public :: elsi_ev_real            !< Compute eigenvalues and eigenvectors
   public :: elsi_ev_complex         !< Compute eigenvalues and eigenvectors
   public :: elsi_dm_real            !< Compute density matrix
   public :: elsi_dm_complex         !< Compute density matrix
   public :: elsi_dm_real_sparse     !< Compute density matrix
   public :: elsi_compute_mu_and_occ !< Compute chemical potential and occupation numbers
   public :: elsi_collect_pexsi      !< Collect additional PEXSI results
   public :: elsi_finalize           !< Clean memory and print timings

   integer(kind=i4), external :: numroc

contains

!=====================
! ELSI tools:
!
!   elsi_init
!   elsi_set_mpi
!   elsi_set_blacs
!   elsi_set_csc
!   elsi_get_energy
!   elsi_finalize   
!=====================

!>
!! This routine initializes ELSI with the solver, parallel mode, matrix storage
!! format, number of basis functions (global size of the Hamiltonian matrix),
!! number of electrons, and number of states.
!!
subroutine elsi_init(elsi_h,solver,parallel_mode,matrix_storage_format,&
                     n_basis,n_electron,n_state)

   implicit none

   type(elsi_handle), intent(out) :: elsi_h                !< Handle of this ELSI instance
   integer(kind=i4),  intent(in)  :: solver                !< AUTO,ELPA,LIBOMM,PEXSI,CHESS,SIPS
   integer(kind=i4),  intent(in)  :: parallel_mode         !< SINGLE_PROC,MULTI_PROC
   integer(kind=i4),  intent(in)  :: matrix_storage_format !< BLACS_DENSE,PEXSI_CSC
   integer(kind=i4),  intent(in)  :: n_basis               !< Number of basis functions
   real(kind=r8),     intent(in)  :: n_electron            !< Number of electrons
   integer(kind=i4),  intent(in)  :: n_state               !< Number of states

   character*40, parameter :: caller = "elsi_init"

   elsi_h%handle_initialized    = .true.
   elsi_h%n_g_size              = n_basis
   elsi_h%n_nonsingular         = n_basis
   elsi_h%n_electrons           = n_electron
   elsi_h%solver                = solver
   elsi_h%matrix_storage_format = matrix_storage_format
   elsi_h%parallel_mode         = parallel_mode
   elsi_h%n_elsi_calls          = 0

   if(parallel_mode == SINGLE_PROC) then
      elsi_h%n_l_rows = n_basis
      elsi_h%n_l_cols = n_basis
      elsi_h%n_b_rows = n_basis
      elsi_h%n_b_cols = n_basis
   endif

   if(solver == LIBOMM) then
      ! Set number of occupied states for libOMM
      elsi_h%n_states = nint(elsi_h%n_electrons/2.0_r8)
      ! Set libOMM default settings
      call elsi_set_omm_default_options(elsi_h)
   else
      elsi_h%n_states = n_state
   endif

   if(solver == PEXSI) then
      ! Set PEXSI default settings
      call elsi_set_pexsi_default_options(elsi_h)
   endif

   if(solver == SIPS) then
      ! Set SIPs default settings
      call elsi_set_sips_default_options(elsi_h)
   endif

   call elsi_init_timers(elsi_h)

end subroutine

!>
!! This routine sets the solver.
!!
subroutine elsi_set_solver(elsi_h,solver)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance
   integer(kind=i4),  intent(in)    :: solver !< AUTO,ELPA,LIBOMM,PEXSI,CHESS,SIPS

   character*40, parameter :: caller = "elsi_set_solver"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%solver = solver

end subroutine

!>
!! Set MPI.
!!
subroutine elsi_set_mpi(elsi_h,mpi_comm)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h      !< Handle of this ELSI instance
   integer(kind=i4),  intent(in)    :: mpi_comm !< MPI communicator

   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_set_mpi"

   call elsi_check_handle(elsi_h,caller)

   if(elsi_h%parallel_mode == MULTI_PROC) then
      elsi_h%mpi_comm = mpi_comm

      call MPI_Comm_rank(mpi_comm,elsi_h%myid,mpierr)
      call MPI_Comm_size(mpi_comm,elsi_h%n_procs,mpierr)

      elsi_h%mpi_is_setup = .true.
   endif

end subroutine

!>
!! Set BLACS.
!!
subroutine elsi_set_blacs(elsi_h,blacs_ctxt,block_size)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h     !< Handle of this ELSI instance
   integer(kind=i4),  intent(in)    :: blacs_ctxt !< BLACS context
   integer(kind=i4),  intent(in)    :: block_size !< Block size

   integer(kind=i4) :: i,i_row,i_col
   integer(kind=i4) :: blacs_info

   character*40, parameter :: caller = "elsi_set_blacs"

   call elsi_check_handle(elsi_h,caller)

   if(elsi_h%parallel_mode == MULTI_PROC) then
      elsi_h%blacs_ctxt = blacs_ctxt
      elsi_h%n_b_rows = block_size
      elsi_h%n_b_cols = block_size

      ! Get processor grid information
      call blacs_gridinfo(elsi_h%blacs_ctxt,elsi_h%n_p_rows,elsi_h%n_p_cols,&
                          elsi_h%my_p_row,elsi_h%my_p_col)

      ! Get local size of matrix
      elsi_h%n_l_rows = numroc(elsi_h%n_g_size,elsi_h%n_b_rows,&
                               elsi_h%my_p_row,0,elsi_h%n_p_rows)
      elsi_h%n_l_cols = numroc(elsi_h%n_g_size,elsi_h%n_b_cols,&
                               elsi_h%my_p_col,0,elsi_h%n_p_cols)

      ! Get BLACS descriptor
      call descinit(elsi_h%sc_desc,elsi_h%n_g_size,elsi_h%n_g_size,&
                    elsi_h%n_b_rows,elsi_h%n_b_cols,0,0,elsi_h%blacs_ctxt,&
                    max(1,elsi_h%n_l_rows),blacs_info)

      ! Get ELPA communicators
      call elsi_get_elpa_comms(elsi_h)

      ! Compute global-local mapping
      call elsi_allocate(elsi_h,elsi_h%local_row,elsi_h%n_g_size,"local_row",caller)
      call elsi_allocate(elsi_h,elsi_h%local_col,elsi_h%n_g_size,"local_col",caller)

      i_row = 0
      i_col = 0

      do i = 1,elsi_h%n_g_size
         if(mod((i-1)/elsi_h%n_b_rows,elsi_h%n_p_rows) == elsi_h%my_p_row) then
            i_row = i_row+1
            elsi_h%local_row(i) = i_row
         endif
         if(mod((i-1)/elsi_h%n_b_cols,elsi_h%n_p_cols) == elsi_h%my_p_col) then
            i_col = i_col+1
            elsi_h%local_col(i) = i_col
         endif
      enddo

      ! Set up MatrixSwitch
      if(elsi_h%solver == LIBOMM) then
         call ms_scalapack_setup(elsi_h%mpi_comm,elsi_h%n_p_rows,'r',&
                                 elsi_h%n_b_rows,icontxt=elsi_h%blacs_ctxt)
      endif

      elsi_h%blacs_is_setup = .true.
   endif

end subroutine

!>
!! This routine sets the sparsity pattern.
!!
subroutine elsi_set_csc(elsi_h,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h            !< Handle of this ELSI instance
   integer(kind=i4),  intent(in)    :: nnz_g             !< Global number of nonzeros
   integer(kind=i4),  intent(in)    :: nnz_l             !< Local number of nonzeros
   integer(kind=i4),  intent(in)    :: n_l_cols          !< Local number of columns
   integer(kind=i4)                 :: row_ind(nnz_l)    !< Row index
   integer(kind=i4)                 :: col_ptr(n_l_cols) !< Column pointer

   character*40, parameter :: caller = "elsi_set_csc"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%nnz_g          = nnz_g
   elsi_h%nnz_l_pexsi    = nnz_l
   elsi_h%n_l_cols_pexsi = n_l_cols

   call elsi_set_row_ind(elsi_h,row_ind)
   call elsi_set_col_ptr(elsi_h,col_ptr)

   elsi_h%sparsity_pattern_ready = .true.

end subroutine

!>
!! This routine gets the energy.
!!
subroutine elsi_get_energy(elsi_h,energy)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance
   real(kind=r8),     intent(out)   :: energy !< Energy of the system

   ! Only spin-nonpolarized case is supported now.
   real(kind=r8), parameter :: n_spin = 2.0_r8
   integer(kind=i4) :: i_state
   character*200 :: info_str

   character*40, parameter :: caller = "elsi_get_energy"

   call elsi_check_handle(elsi_h,caller)

   select case (elsi_h%solver)
      case (ELPA)
         energy = 0.0_r8
         do i_state =1,elsi_h%n_states
            energy = energy+elsi_h%occ_elpa(i_state)*elsi_h%eval(i_state)
         enddo
      case (LIBOMM)
         energy = n_spin*elsi_h%total_energy
      case (PEXSI)
         energy = elsi_h%e_tot_H
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine finalizes ELSI.
!!
subroutine elsi_finalize(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle of this ELSI instance

   character*40, parameter :: caller = "elsi_finalize"

   call elsi_check_handle(elsi_h,caller)
   call elsi_final_print(elsi_h)
   call elsi_cleanup(elsi_h)

end subroutine

!=========================
! ELSI customize routines
!
!   elsi_customize
!   elsi_customize_omm
!   elsi_customize_pexsi
!   elsi_customize_elpa
!   elsi_customize_mu
!   elsi_collect_pexsi
!=========================

!>
!! This routine overrides ELSI default settings.
!!
subroutine elsi_customize(elsi_h,print_detail,overlap_is_unit,zero_threshold,&
                          no_singularity_check,singularity_tolerance,stop_singularity)

   implicit none

   type(elsi_handle), intent(inout)        :: elsi_h                 !< Handle of this ELSI instance
   logical,           intent(in), optional :: print_detail           !< Print detailed info?
   logical,           intent(in), optional :: overlap_is_unit        !< Is overlap matrix unit?
   real(kind=r8),     intent(in), optional :: zero_threshold         !< Threshold to define "zero"
   logical,           intent(in), optional :: no_singularity_check   !< Do not perform singularity check of overlap
   real(kind=r8),     intent(in), optional :: singularity_tolerance  !< Tolerance of overlap singularity
   logical,           intent(in), optional :: stop_singularity       !< Stop if overlap is singular

   character*40, parameter :: caller = "elsi_customize"

   call elsi_check_handle(elsi_h,caller)

   ! Print detailed ELSI information? [Default: .false.]
   if(present(print_detail)) &
      print_info = print_detail
   ! Is the overlap matrix unit? [Default: .false.]
   if(present(overlap_is_unit)) &
      elsi_h%overlap_is_unit = overlap_is_unit
   ! Threshold to define numerical zero [Default: 1e-13_r8]
   if(present(zero_threshold)) &
      elsi_h%zero_threshold = zero_threshold
   ! Disable checking for overlap singularity? [Default: .false.]
   if(present(no_singularity_check)) &
      elsi_h%no_singularity_check = no_singularity_check
   ! Eigenfunctions of overlap matrix with eigenvalues smaller than
   ! this value will be removed to avoid singularity [Default: 1e-5_r8]
   if(present(singularity_tolerance)) &
      elsi_h%singularity_tolerance = singularity_tolerance
   ! Always stop if overlap is singular? [Default: .false.]
   if(present(stop_singularity)) &
      elsi_h%stop_singularity = stop_singularity

end subroutine

!>
!! This routine overrides libOMM default settings.
!!
subroutine elsi_customize_omm(elsi_h,n_elpa_steps,omm_flavor,eigen_shift,&
                              omm_tolerance,use_pspblas,omm_output)

   implicit none

   type(elsi_handle), intent(inout)        :: elsi_h        !< Handle of this ELSI instance
   integer(kind=i4),  intent(in), optional :: n_elpa_steps  !< Number of ELPA steps before libOMM
   integer(kind=i4),  intent(in), optional :: omm_flavor    !< How to perform orbital minimization
   real(kind=r8),     intent(in), optional :: eigen_shift   !< Eigenspectrum shift parameter
   real(kind=r8),     intent(in), optional :: omm_tolerance !< Tolerance of minimization
   logical,           intent(in), optional :: use_pspblas   !< Use pspBLAS sparse linear algebra?
   logical,           intent(in), optional :: omm_output    !< Output details?

   character*40, parameter :: caller = "elsi_customize_omm"

   call elsi_check_handle(elsi_h,caller)

   ! Number of ELPA steps
   if(present(n_elpa_steps)) &
      elsi_h%n_elpa_steps = n_elpa_steps
   ! How to perform orbital minimization?
   if(present(omm_flavor)) &
      elsi_h%omm_flavor = omm_flavor
   ! Eigenspectrum shift parameter
   if(present(eigen_shift)) &
      elsi_h%eta = eigen_shift
   ! Tolerance for minimization
   if(present(omm_tolerance)) &
      elsi_h%min_tol = omm_tolerance
   ! Use pspBLAS sparse linear algebra?
   if(present(use_pspblas)) &
      elsi_h%use_psp = use_pspblas
   ! Output details?
   if(present(omm_output)) &
      elsi_h%omm_verbose = omm_output

   if(elsi_h%solver .ne. LIBOMM) then
      call elsi_statement_print("  The chosen solver is not libOMM."//&
                                " Ignore elsi_customize_omm call.",elsi_h)
   endif

end subroutine

!>
!! This routine overrides PEXSI default settings.
!!
subroutine elsi_customize_pexsi(elsi_h,temperature,gap,delta_e,n_poles,n_procs_per_pole,&
                                max_iteration,mu_min,mu_max,mu0,mu_inertia_tolerance,&
                                mu_inertia_expansion,mu_safeguard,n_electron_accuracy,&
                                matrix_type,is_symbolic_factorize,ordering,&
                                np_symbolic_factorize,verbosity)

   implicit none

   type(elsi_handle), intent(inout)        :: elsi_h                !< Handle of this ELSI instance
   real(kind=r8),     intent(in), optional :: temperature           !< Temperature, in the same unit as Hamiltonian
   real(kind=r8),     intent(in), optional :: gap                   !< Spectral gap, can be set to 0 in most cases
   real(kind=r8),     intent(in), optional :: delta_e               !< Upper bound for the spectral radius of S^(-1)H
   integer(kind=i4),  intent(in), optional :: n_poles               !< Number of poles
   integer(kind=i4),  intent(in), optional :: n_procs_per_pole      !< Number of processors for one pole
   integer(kind=i4),  intent(in), optional :: max_iteration         !< Maximum number of PEXSI iterations
   real(kind=r8),     intent(in), optional :: mu_min                !< Lower bound of chemical potential
   real(kind=r8),     intent(in), optional :: mu_max                !< Upper bound of chemical potential
   real(kind=r8),     intent(in), optional :: mu0                   !< Initial guess of chemical potential
   real(kind=r8),     intent(in), optional :: mu_inertia_tolerance  !< Tolerance of inertia counting
   real(kind=r8),     intent(in), optional :: mu_inertia_expansion  !< Expansion step size in inertia counting
   real(kind=r8),     intent(in), optional :: mu_safeguard          !< Safeguard to reinvoke inertia counting
   real(kind=r8),     intent(in), optional :: n_electron_accuracy   !< Accuracy of number of electrons
   integer(kind=i4),  intent(in), optional :: matrix_type           !< Type of input matrices
   integer(kind=i4),  intent(in), optional :: is_symbolic_factorize !< Perform symbolic factorization?
   integer(kind=i4),  intent(in), optional :: ordering              !< Ordering strategy for factorization and selected inversion
   integer(kind=i4),  intent(in), optional :: np_symbolic_factorize !< Number of processors for ParMETIS, only used if ordering=0
   integer(kind=i4),  intent(in), optional :: verbosity             !< Level of output info

   character*40, parameter :: caller = "elsi_customize_pexsi"

   call elsi_check_handle(elsi_h,caller)

   ! Temperature, in the same unit as H
   if(present(temperature)) &
      elsi_h%pexsi_options%temperature = temperature

   ! Spectral gap, can be set to 0 in most cases (default)
   if(present(gap)) &
      elsi_h%pexsi_options%gap = gap

   ! Upper bound for the spectral radius of S^(-1)H
   ! default: 10
   if(present(delta_e)) &
      elsi_h%pexsi_options%deltaE = delta_e

   ! Number of poles
   ! default: 40
   if(present(n_poles)) &
      elsi_h%pexsi_options%numPole = n_poles

   ! Number of processors for one pole
   ! default: decided from n_procs and n_poles
   if(present(n_procs_per_pole)) then
      if((mod(elsi_h%n_procs,n_procs_per_pole) == 0) .and. &
         (elsi_h%matrix_storage_format .ne. 0)) then
         elsi_h%n_p_per_pole_pexsi = n_procs_per_pole
         elsi_h%n_p_per_pole_ready = .true.
      else
         call elsi_stop("  The total number of MPI tasks must be a"//&
                        " multiple of the number of MPI tasks per"//&
                        " pole. Exiting...",elsi_h,caller)
      endif
   endif

   ! Maximum number of PEXSI iterations after each inertia
   ! counting procedure
   ! default: 3
   if(present(max_iteration)) &
      elsi_h%pexsi_options%maxPEXSIIter = max_iteration

   ! From the second step, initial guess of mu is from previous step
   if(elsi_h%n_elsi_calls == 0) then
      ! Initial guess of mu
      ! default: 0.0
      if(present(mu0)) &
         elsi_h%pexsi_options%mu0 = mu0
   endif

   ! Initial guess of lower bound for mu
   ! default: -10.0
   if(present(mu_min)) &
      elsi_h%pexsi_options%muMin0 = mu_min

   ! Initial guess of upper bound for mu
   ! default: 10.0
   if(present(mu_max)) &
      elsi_h%pexsi_options%muMax0 = mu_max

   ! Stopping criterion in terms of the chemical potential
   ! for the inertia counting procedure
   ! default: 0.05
   if(present(mu_inertia_tolerance)) &
      elsi_h%pexsi_options%muInertiaTolerance = mu_inertia_tolerance

   ! If the chemical potential is not in the initial interval,
   ! the interval is expanded by this value
   ! default: 0.3
   if(present(mu_inertia_expansion)) &
      elsi_h%pexsi_options%muInertiaExpansion = mu_inertia_expansion

   ! Safeguard criterion in terms of the chemical potential to
   ! reinvoke the inertia counting procedure
   ! default: 0.05
   if(present(mu_safeguard)) &
      elsi_h%pexsi_options%muPEXSISafeGuard = mu_safeguard

   ! Stopping criterion of the PEXSI iteration in terms of the
   ! number of electrons compared to the exact number
   ! default: 0.01
   if(present(n_electron_accuracy)) then
      elsi_h%pexsi_options%numElectronPEXSITolerance = n_electron_accuracy
      if(n_electron_accuracy < 1.0e-2_r8) then
         elsi_h%small_pexsi_tol = .true.
         elsi_h%final_pexsi_tol = n_electron_accuracy
      endif
   endif

   ! Type of input H and S matrices
   ! 0: real symmetric (default)
   ! 1: general complex
   if(present(matrix_type)) &
      elsi_h%pexsi_options%matrixType = matrix_type

   ! Whether to perform symbolic factorization
   ! default: 1
   if(present(is_symbolic_factorize)) &
      elsi_h%pexsi_options%isSymbolicFactorize = is_symbolic_factorize

   ! Ordering strategy for factorization and selected inversion
   ! 0: parallel ordering using ParMETIS
   ! 1: sequential ordering using METIS
   ! 2: multiple minimum degree ordering
   if(present(ordering)) &
      elsi_h%pexsi_options%ordering = ordering

   ! Number of processors for ParMETIS, only used if ordering=0
   if(present(np_symbolic_factorize)) &
      elsi_h%pexsi_options%npSymbFact = np_symbolic_factorize

   ! Level of output information
   ! 0: no output
   ! 1: basic output (default)
   ! 2: detailed output
   if(present(verbosity)) &
      elsi_h%pexsi_options%verbosity = verbosity

   if(elsi_h%solver .ne. PEXSI) then
      call elsi_statement_print("  The chosen solver is not PEXSI."//&
                                " Ignore elsi_customize_pexsi call.",elsi_h)
   endif

end subroutine

!>
!! This routine overrides ELPA default settings.
!!
subroutine elsi_customize_elpa(elsi_h,elpa_solver)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h      !< Handle of this ELSI instance
   integer(kind=i4),  intent(in)    :: elpa_solver !< Always use 1-stage or 2-stage solver

   character*40, parameter :: caller = "elsi_customize_elpa"

   call elsi_check_handle(elsi_h,caller)

   if(elpa_solver == 1) then
      elsi_h%elpa_one_always = .true.
      elsi_h%elpa_two_always = .false.
   elseif(elpa_solver == 2) then
      elsi_h%elpa_one_always = .false.
      elsi_h%elpa_two_always = .true.
   endif

   if(elsi_h%solver .ne. ELPA) then
      call elsi_statement_print("  The chosen solver is not ELPA."//&
                                " Ignore elsi_customize_elpa call.",elsi_h)
   endif

end subroutine

!>
!! This routine overrides ELSI default settings for the chemical potential
!! determination module.
!!
subroutine elsi_customize_mu(elsi_h,broadening_scheme,broadening_width,&
                             occ_accuracy,mu_max_steps,spin_degeneracy)

   implicit none

   type(elsi_handle), intent(inout)        :: elsi_h            !< Handle of this ELSI instance
   integer(kind=i4),  intent(in), optional :: broadening_scheme !< Broadening method in chemical potential determination
   real(kind=r8),     intent(in), optional :: broadening_width  !< Broadening width in chemical potential determination
   real(kind=r8),     intent(in), optional :: occ_accuracy      !< Accuracy in electron count (sum of occ)
   integer(kind=i4),  intent(in), optional :: mu_max_steps      !< Maximum number of steps to find the chemical potential
   real(kind=r8),     intent(in), optional :: spin_degeneracy   !< Spin degeneracy

   character*40, parameter :: caller = "elsi_customize_mu"

   call elsi_check_handle(elsi_h,caller)

   ! Broadening scheme to compute Fermi level [Default: GAUSSIAN]
   if(present(broadening_scheme)) &
      elsi_h%broadening_scheme = broadening_scheme
   ! Broadening width to compute Fermi level [Default: 1e-2_r8]
   if(present(broadening_width)) &
      elsi_h%broadening_width = broadening_width
   ! Accuracy for chemical potential determination [Default: 1e-10_r8]
   if(present(occ_accuracy)) &
      elsi_h%occ_tolerance = occ_accuracy
   ! Maximum steps to determine the chemical potential [Default: 100]
   if(present(mu_max_steps)) &
      elsi_h%max_mu_steps = mu_max_steps
   ! Spin degeneracy [Default: 2.0_r8/n_spin]
   if(present(spin_degeneracy)) &
      elsi_h%spin_degen = spin_degeneracy

end subroutine

!>
!! This routine collects results (other than density matrix)
!! after a PEXSI calculation.
!!
subroutine elsi_collect_pexsi(elsi_h,mu,edm,fdm)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h                  !< Handle of this ELSI instance
   real(kind=r8),     intent(out) :: mu                      !< Chemical potential
   real(kind=r8),     intent(out) :: edm(elsi_h%nnz_l_pexsi) !< Energy density matrix
   real(kind=r8),     intent(out) :: fdm(elsi_h%nnz_l_pexsi) !< Free energy density matrix

   character*40, parameter :: caller = "elsi_collect_pexsi"

   call elsi_check_handle(elsi_h,caller)

   mu = elsi_h%mu_pexsi
   edm(1:elsi_h%nnz_l_pexsi) = elsi_h%e_den_mat_pexsi(1:elsi_h%nnz_l_pexsi)
   fdm(1:elsi_h%nnz_l_pexsi) = elsi_h%f_den_mat_pexsi(1:elsi_h%nnz_l_pexsi)

end subroutine

!=======================
! ELSI solvers
!
!   elsi_ev_real
!   elsi_ev_complex
!   elsi_ev_real_sparse
!   elsi_dm_real
!   elsi_dm_complex
!   elsi_dm_real_sparse
!=======================

!>
!! This routine computes eigenvalues and eigenvectors.
!!
subroutine elsi_ev_real(elsi_h,H_in,S_in,e_val_out,e_vec_out)

   implicit none

   type(elsi_handle) :: elsi_h                                     !< Handle of this ELSI instance
   real(kind=r8)     :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols)      !< Hamiltonian
   real(kind=r8)     :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols)      !< Overlap
   real(kind=r8)     :: e_val_out(elsi_h%n_g_size)                 !< Eigenvalues
   real(kind=r8)     :: e_vec_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_real"

   call elsi_check_handle(elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! REAL case
   elsi_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   ! Here the only supported solver is ELPA
   select case (elsi_h%solver)
      case (ELPA)
         ! Set matrices
         call elsi_set_hamiltonian(elsi_h,H_in)
         if(.not. elsi_h%overlap_is_unit) then
            call elsi_set_overlap(elsi_h,S_in)
         endif
         call elsi_set_eigenvector(elsi_h,e_vec_out)
         call elsi_set_eigenvalue(elsi_h,e_val_out)

         ! Solve eigenvalue problem
         if(elsi_h%parallel_mode == SINGLE_PROC) then
            call elsi_solve_evp_elpa_sp(elsi_h)
         else ! Multi-proc
            call elsi_solve_evp_elpa(elsi_h)
         endif

      case (LIBOMM)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",elsi_h,caller)
      case (PEXSI)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",elsi_h,caller)
      case (CHESS)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_print_sips_options(elsi_h)

         ! Initialize SIPs
         call elsi_init_sips(elsi_h)

         ! Convert matrix format and distribution from BLACS to SIPs
         call elsi_blacs_to_sips(elsi_h,H_in,S_in)

         ! Set matrices
         call elsi_set_sparse_hamiltonian(elsi_h,elsi_h%ham_real_sips)
         if(.not. elsi_h%overlap_is_unit) then
            call elsi_set_sparse_overlap(elsi_h,elsi_h%ovlp_real_sips)
         endif
         call elsi_set_row_ind(elsi_h,elsi_h%row_ind_sips)
         call elsi_set_col_ptr(elsi_h,elsi_h%col_ptr_sips)
         call elsi_set_eigenvalue(elsi_h,e_val_out)
         call elsi_set_eigenvector(elsi_h,e_vec_out)

         call elsi_solve_evp_sips(elsi_h)

      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA solver to compute"//&
                        " eigenvalues and eigenvectors. Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET

end subroutine

!>
!! This routine computes eigenvalues and eigenvectors.
!!
subroutine elsi_ev_complex(elsi_h,H_in,S_in,e_val_out,e_vec_out)

   implicit none

   type(elsi_handle) :: elsi_h                                     !< Handle of this ELSI instance
   complex(kind=r8)  :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols)      !< Hamiltonian
   complex(kind=r8)  :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols)      !< Overlap
   real(kind=r8)     :: e_val_out(elsi_h%n_g_size)                 !< Eigenvalues
   complex(kind=r8)  :: e_vec_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_complex"

   call elsi_check_handle(elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! COMPLEX case
   elsi_h%matrix_data_type = COMPLEX_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   ! Here the only supported solver is ELPA
   select case (elsi_h%solver)
      case (ELPA)
         ! Set matrices
         call elsi_set_hamiltonian(elsi_h,H_in)
         if(.not. elsi_h%overlap_is_unit) then
            call elsi_set_overlap(elsi_h,S_in)
         endif
         call elsi_set_eigenvector(elsi_h,e_vec_out)
         call elsi_set_eigenvalue(elsi_h,e_val_out)

         ! Solve eigenvalue problem
         if(elsi_h%parallel_mode == SINGLE_PROC) then
            call elsi_solve_evp_elpa_sp(elsi_h)
         else ! Multi-proc
            call elsi_solve_evp_elpa(elsi_h)
         endif

      case (LIBOMM)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",elsi_h,caller)
      case (PEXSI)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",elsi_h,caller)
      case (CHESS)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA solver to compute"//&
                        " eigenvalues and eigenvectors. Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET

end subroutine

!>
!! This routine computes eigenvalues and eigenvectors.
!!
subroutine elsi_ev_real_sparse(elsi_h,H_in,S_in,e_val_out,e_vec_out)

   implicit none

   type(elsi_handle) :: elsi_h                                     !< Handle of this ELSI instance
   real(kind=r8)     :: H_in(elsi_h%nnz_l_pexsi)                   !< Hamiltonian
   real(kind=r8)     :: S_in(elsi_h%nnz_l_pexsi)                   !< Overlap
   real(kind=r8)     :: e_val_out(elsi_h%n_g_size)                 !< Eigenvalues
   real(kind=r8)     :: e_vec_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_real_sparse"

   call elsi_check_handle(elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! REAL case
   elsi_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   ! Here the only supported solver is ELPA
   select case (elsi_h%solver)
      case (ELPA)
         ! Convert matrix format and distribution from PEXSI to BLACS
         call elsi_pexsi_to_blacs(elsi_h,H_in,S_in)

         ! Set matrices
         call elsi_set_hamiltonian(elsi_h,elsi_h%ham_real_elpa)
         if(.not. elsi_h%overlap_is_unit) then
            call elsi_set_overlap(elsi_h,elsi_h%ovlp_real_elpa)
         endif
         call elsi_set_eigenvector(elsi_h,e_vec_out)
         call elsi_set_eigenvalue(elsi_h,e_val_out)

         ! Solve eigenvalue problem
         call elsi_solve_evp_elpa(elsi_h)

      case (LIBOMM)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",elsi_h,caller)
      case (PEXSI)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",elsi_h,caller)
      case (CHESS)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",elsi_h,caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA solver to compute"//&
                        " eigenvalues and eigenvectors. Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET

end subroutine

!>
!! This routine computes density matrix.
!!
subroutine elsi_dm_real(elsi_h,H_in,S_in,D_out,energy_out)

   implicit none

   type(elsi_handle) :: elsi_h                                 !< Handle of this ELSI instance
   real(kind=r8)     :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols)  !< Hamiltonian
   real(kind=r8)     :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols)  !< Overlap
   real(kind=r8)     :: D_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix
   real(kind=r8)     :: energy_out                             !< Energy

   character*40, parameter :: caller = "elsi_dm_real"

   call elsi_check_handle(elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! REAL case
   elsi_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   ! Solve eigenvalue problem
   select case (elsi_h%solver)
      case (ELPA)
         ! Set matrices
         if(.not. allocated(elsi_h%eval_elpa)) then
            call elsi_allocate(elsi_h,elsi_h%eval_elpa,elsi_h%n_g_size,&
                               "eval_elpa",caller)
         endif
         if(.not. allocated(elsi_h%evec_real_elpa)) then
            call elsi_allocate(elsi_h,elsi_h%evec_real_elpa,elsi_h%n_l_rows,&
                               elsi_h%n_l_cols,"evec_real_elpa",caller)
         endif

         call elsi_set_hamiltonian(elsi_h,H_in)
         if(.not. elsi_h%overlap_is_unit) then
            call elsi_set_overlap(elsi_h,S_in)
         endif
         call elsi_set_eigenvector(elsi_h,elsi_h%evec_real_elpa)
         call elsi_set_eigenvalue(elsi_h,elsi_h%eval_elpa)
         call elsi_set_density_matrix(elsi_h,D_out)

         call elsi_solve_evp_elpa(elsi_h)

         call elsi_compute_occ_elpa(elsi_h)
         call elsi_compute_dm_elpa(elsi_h)
         call elsi_get_energy(elsi_h,energy_out)

      case (LIBOMM)
         call elsi_print_omm_options(elsi_h)

         if(elsi_h%n_elsi_calls .le. elsi_h%n_elpa_steps) then
            if((elsi_h%n_elsi_calls == 1) .and. (elsi_h%omm_flavor == 0)) then
               ! Overlap will be destroyed by Cholesky
               call elsi_allocate(elsi_h,elsi_h%ovlp_real_omm,elsi_h%n_l_rows,&
                                  elsi_h%n_l_cols,"ovlp_real_omm",caller)
               elsi_h%ovlp_real_omm(1:elsi_h%n_l_rows,1:elsi_h%n_l_cols) = &
                  S_in(1:elsi_h%n_l_rows,1:elsi_h%n_l_cols)
            endif

            ! Compute libOMM initial guess by ELPA
            elsi_h%solver = ELPA

            ! Set matrices
            if(.not. allocated(elsi_h%eval_elpa)) then
               call elsi_allocate(elsi_h,elsi_h%eval_elpa,elsi_h%n_g_size,&
                                  "eval_elpa",caller)
            endif
            if(.not. allocated(elsi_h%evec_real_elpa)) then
               call elsi_allocate(elsi_h,elsi_h%evec_real_elpa,elsi_h%n_l_rows,&
                                  elsi_h%n_l_cols,"evec_real_elpa",caller)
            endif

            call elsi_set_hamiltonian(elsi_h,H_in)
            if(.not. elsi_h%overlap_is_unit) then
               call elsi_set_overlap(elsi_h,S_in)
            endif
            call elsi_set_eigenvector(elsi_h,elsi_h%evec_real_elpa)
            call elsi_set_eigenvalue(elsi_h,elsi_h%eval_elpa)
            call elsi_set_density_matrix(elsi_h,D_out)

            ! Solve by ELPA
            call elsi_solve_evp_elpa(elsi_h)
            call elsi_compute_occ_elpa(elsi_h)
            call elsi_compute_dm_elpa(elsi_h)
            call elsi_get_energy(elsi_h,energy_out)

            ! Switch back to libOMM here to guarantee elsi_customize_omm
            elsi_h%solver = LIBOMM

         else ! ELPA is done
            if(allocated(elsi_h%ovlp_real_omm)) then
               ! Retrieve overlap matrix that has been destroyed by Cholesky
               S_in(1:elsi_h%n_l_rows,1:elsi_h%n_l_cols) = &
                  elsi_h%ovlp_real_omm(1:elsi_h%n_l_rows,1:elsi_h%n_l_cols)
               deallocate(elsi_h%ovlp_real_omm)
            endif

            ! Set matrices
            call elsi_set_hamiltonian(elsi_h,H_in)
            if(.not. elsi_h%overlap_is_unit) then
               call elsi_set_overlap(elsi_h,S_in)
            endif
            call elsi_set_density_matrix(elsi_h,D_out)

            if(.not. elsi_h%coeff_omm%is_initialized) then
               call m_allocate(elsi_h%coeff_omm,elsi_h%n_states,&
                               elsi_h%n_g_size,"pddbc")
            endif

            ! Initialize coefficient matrix with ELPA eigenvectors if possible
            if((elsi_h%n_elpa_steps > 0) .and. &
               (elsi_h%n_elsi_calls == elsi_h%n_elpa_steps+1)) then
               ! D_out is used for temporary storage here
               D_out(1:elsi_h%n_l_rows,1:elsi_h%n_l_cols) = &
                     elsi_h%evec_real(1:elsi_h%n_l_rows,1:elsi_h%n_l_cols)
               ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
               call pdtran(elsi_h%n_g_size,elsi_h%n_g_size,1.0_r8,D_out,1,1,&
                           elsi_h%sc_desc,0.0_r8,elsi_h%evec_real,1,1,elsi_h%sc_desc)

               elsi_h%coeff_omm%dval(1:elsi_h%coeff_omm%iaux2(1),1:elsi_h%coeff_omm%iaux2(2)) = &
                  elsi_h%evec_real(1:elsi_h%coeff_omm%iaux2(1),1:elsi_h%coeff_omm%iaux2(2))

               ! ELPA matrices are no longer needed
               if(associated(elsi_h%ham_real))      nullify(elsi_h%ham_real)
               if(associated(elsi_h%ovlp_real))     nullify(elsi_h%ovlp_real)
               if(associated(elsi_h%evec_real))     nullify(elsi_h%evec_real)
               if(associated(elsi_h%den_mat))       nullify(elsi_h%den_mat)
               if(associated(elsi_h%eval))          nullify(elsi_h%eval)
               if(allocated(elsi_h%evec_real_elpa)) deallocate(elsi_h%evec_real_elpa)
               if(allocated(elsi_h%eval_elpa))      deallocate(elsi_h%eval_elpa)
               if(allocated(elsi_h%occ_elpa))       deallocate(elsi_h%occ_elpa)
            endif

            ! Continue the computation using libOMM
            call elsi_solve_evp_omm(elsi_h)

            elsi_h%den_mat_omm%dval = 2.0_r8*elsi_h%den_mat_omm%dval
            call elsi_get_energy(elsi_h,energy_out)
         endif

      case (PEXSI)
         call elsi_print_pexsi_options(elsi_h)

         ! PEXSI may use different process grid to achieve
         ! the efficient 2-level parallelization
         call elsi_init_pexsi(elsi_h)

         ! Convert matrix format and distribution from BLACS to PEXSI
         call elsi_blacs_to_pexsi(elsi_h,H_in,S_in)

         ! Set matrices
         call elsi_set_sparse_hamiltonian(elsi_h,elsi_h%ham_real_pexsi)
         if(.not. elsi_h%overlap_is_unit) then
            call elsi_set_sparse_overlap(elsi_h,elsi_h%ovlp_real_pexsi)
         endif
         call elsi_set_row_ind(elsi_h,elsi_h%row_ind_pexsi)
         call elsi_set_col_ptr(elsi_h,elsi_h%col_ptr_pexsi)
         if(.not. allocated(elsi_h%den_mat_pexsi)) then
            call elsi_allocate(elsi_h,elsi_h%den_mat_pexsi,elsi_h%nnz_l_pexsi,&
                               "den_mat_pexsi",caller)
         endif
         elsi_h%den_mat_pexsi = 0.0_r8
         call elsi_set_sparse_density_matrix(elsi_h,elsi_h%den_mat_pexsi)

         call elsi_solve_evp_pexsi(elsi_h)

         ! Convert 1D block CCS sparse density matrix to 2D
         ! block-cyclic dense format
         call elsi_pexsi_to_blacs_dm(elsi_h,D_out)
         call elsi_get_energy(elsi_h,energy_out)

      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case default
         call elsi_stop(" No supported solver has been chosen."//&
                        " Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET

end subroutine

!>
!! This routine computes density matrix.
!!
subroutine elsi_dm_complex(elsi_h,H_in,S_in,D_out,energy_out)

   implicit none

   type(elsi_handle) :: elsi_h                                 !< Handle of this ELSI instance
   complex(kind=r8)  :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols)  !< Hamiltonian
   complex(kind=r8)  :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols)  !< Overlap
   real(kind=r8)     :: D_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix
   real(kind=r8)     :: energy_out                             !< Energy

   character*40, parameter :: caller = "elsi_dm_complex"

   call elsi_check_handle(elsi_h,caller)

   call elsi_stop(" ELSI density matrix solver for complex case not yet available."//&
                  " Exiting...",elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! COMPLEX case
   elsi_h%matrix_data_type = COMPLEX_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   ! Solve eigenvalue problem
   select case (elsi_h%solver)
      case (ELPA)
         ! Set matrices
         if(.not. allocated(elsi_h%eval_elpa)) then
            call elsi_allocate(elsi_h,elsi_h%eval_elpa,elsi_h%n_g_size,"eval_elpa",caller)
         endif
         if(.not. allocated(elsi_h%evec_complex_elpa)) then
            call elsi_allocate(elsi_h,elsi_h%evec_complex_elpa,elsi_h%n_l_rows,&
                               elsi_h%n_l_cols,"evec_complex_elpa",caller)
         endif

         call elsi_set_hamiltonian(elsi_h,H_in)
         if(.not. elsi_h%overlap_is_unit) then
            call elsi_set_overlap(elsi_h,S_in)
         endif
         call elsi_set_eigenvector(elsi_h,elsi_h%evec_complex_elpa)
         call elsi_set_eigenvalue(elsi_h,elsi_h%eval_elpa)
         call elsi_set_density_matrix(elsi_h,D_out)

         call elsi_solve_evp_elpa(elsi_h)

         call elsi_compute_occ_elpa(elsi_h)
         call elsi_compute_dm_elpa(elsi_h)
         call elsi_get_energy(elsi_h,energy_out)

      case (LIBOMM)
         call elsi_print_omm_options(elsi_h)

         ! Set Hamiltonian and overlap matrices
         call elsi_set_hamiltonian(elsi_h,H_in)
         if(.not. elsi_h%overlap_is_unit) then
            call elsi_set_overlap(elsi_h,S_in)
         endif
         call elsi_set_density_matrix(elsi_h,D_out)

         if(.not. elsi_h%coeff_omm%is_initialized) then
            call m_allocate(elsi_h%coeff_omm,elsi_h%n_states,&
                            elsi_h%n_g_size,"pddbc")
         endif

         call elsi_solve_evp_omm(elsi_h)

         elsi_h%den_mat_omm%dval = 2.0_r8*elsi_h%den_mat_omm%dval
         call elsi_get_energy(elsi_h,energy_out)

      case (PEXSI)
         call elsi_stop(" PEXSI not yet implemented. Exiting...",elsi_h,caller)
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case default
         call elsi_stop(" No supported solver has been chosen."//&
                        " Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET

end subroutine

!>
!! This routine computes density matrix.
!!
subroutine elsi_dm_real_sparse(elsi_h,H_in,S_in,D_out,energy_out)

   implicit none

   type(elsi_handle) :: elsi_h                    !< Handle of this ELSI instance
   real(kind=r8)     :: H_in(elsi_h%nnz_l_pexsi)  !< Hamiltonian
   real(kind=r8)     :: S_in(elsi_h%nnz_l_pexsi)  !< Overlap
   real(kind=r8)     :: D_out(elsi_h%nnz_l_pexsi) !< Density matrix
   real(kind=r8)     :: energy_out                !< Energy

   character*40, parameter :: caller = "elsi_dm_real_sparse"

   call elsi_check_handle(elsi_h,caller)

   ! Update counter
   elsi_h%n_elsi_calls = elsi_h%n_elsi_calls+1

   ! REAL case
   elsi_h%matrix_data_type = REAL_VALUES

   ! Safety check
   call elsi_check(elsi_h,caller)

   ! Solve eigenvalue problem
   select case (elsi_h%solver)
      case (ELPA)
         call elsi_stop(" ELPA cannot elsi_h sparse matrices."//&
                        " Exiting...",elsi_h,caller)
      case (LIBOMM)
         call elsi_stop(" libOMM cannot elsi_h sparse matrices."//&
                        " Exiting...",elsi_h,caller)
      case (PEXSI)
         call elsi_print_pexsi_options(elsi_h)

         ! Set matrices
         call elsi_set_sparse_hamiltonian(elsi_h,H_in)
         if(.not. elsi_h%overlap_is_unit) then
            call elsi_set_sparse_overlap(elsi_h,S_in)
         endif
         call elsi_set_sparse_density_matrix(elsi_h,D_out)

         call elsi_init_pexsi(elsi_h)

         call elsi_solve_evp_pexsi(elsi_h)

         call elsi_get_energy(elsi_h,energy_out)

      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case default
         call elsi_stop(" No supported solver has been chosen."//&
                        " Exiting...",elsi_h,caller)
   end select

   elsi_h%matrix_data_type = UNSET

end subroutine

end module ELSI
