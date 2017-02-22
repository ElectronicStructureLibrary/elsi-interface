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
!! an eigenvalue problem using ELPA, libOMM, PEXSI, or CheSS.
!!

module ELSI

   use iso_c_binding
   use ELSI_DIMENSIONS
   use ELSI_TIMERS
   use ELSI_UTILS
   use ELSI_MATRIX_CONVERSION
   use ELSI_ELPA
   use ELSI_OMM
   use ELSI_PEXSI
   use MatrixSwitch

   implicit none
   private

   !> Public routines
   public :: elsi_init            !< Initialize
   public :: elsi_set_method      !< Select solver
   public :: elsi_set_mpi         !< Set MPI from calling code
   public :: elsi_set_blacs       !< Set BLACS from calling code
   public :: elsi_set_sparsity    !< Set sparsity pattern from calling code
   public :: elsi_customize       !< Override ELSI default
   public :: elsi_customize_elpa  !< Override ELPA default
   public :: elsi_customize_omm   !< Override libOMM default
   public :: elsi_customize_pexsi !< Override PEXSI default
   public :: elsi_ev_real         !< Compute eigenvalues and eigenvectors
   public :: elsi_ev_complex      !< Compute eigenvalues and eigenvectors
   public :: elsi_dm_real         !< Compute density matrix
   public :: elsi_dm_complex      !< Compute density matrix
   public :: elsi_dm_real_sparse  !< Compute density matrix
   public :: elsi_finalize        !< Clean memory and print timings

   integer, external :: numroc

contains

!=====================
! ELSI tools:
!
!   elsi_init
!   elsi_set_method
!   elsi_set_mode
!   elsi_set_storage
!   elsi_set_parallel
!   elsi_set_mpi
!   elsi_set_blacs
!   elsi_set_sparsity
!   elsi_get_energy
!   elsi_finalize   
!=====================

!>
!! This routine initializes ELSI with solver, parallelism, matrix format,
!! global matrix size, number of electrons, and number of states.
!!
subroutine elsi_init(solver,parallel_mode,matrix_format,matrix_size,&
                     n_electrons_in,n_states_in)

   implicit none

   integer, intent(in) :: solver         !< AUTO,ELPA,LIBOMM,PEXSI,CHESS
   integer, intent(in) :: parallel_mode  !< SINGLE_PROC,MULTI_PROC
   integer, intent(in) :: matrix_format  !< DENSE,CCS,CSC,CRS,CSR
   integer, intent(in) :: matrix_size    !< Global dimension of matrix
   real*8,  intent(in) :: n_electrons_in !< Number of electrons
   integer, intent(in) :: n_states_in    !< Number of states

   n_g_size = matrix_size
   n_nonsingular = matrix_size
   n_electrons = n_electrons_in

   call elsi_set_method(solver)
   call elsi_set_storage(matrix_format)
   call elsi_set_parallel(parallel_mode)

   if(solver == 2) then
      ! Set number of occupied states for libOMM
      n_states = NINT(n_electrons/2.0d0)
      ! Set libOMM default settings
      call elsi_set_omm_default_options()
   else
      n_states = n_states_in
   endif

   if(solver == 3) then
      ! Set PEXSI default settings
      call elsi_set_pexsi_default_options()
   endif

   n_elsi_calls = 0

   call elsi_init_timers()

end subroutine

!>
!! This routine sets the method.
!!
subroutine elsi_set_method(i_method)

   implicit none

   integer, intent(in) :: i_method !< AUTO,ELPA,LIBOMM,PEXSI,CHESS
   
   method = i_method

end subroutine

!>
!! This routine sets the mode (real or complex).
!!
subroutine elsi_set_mode(i_mode)

   implicit none

   integer, intent(in) :: i_mode !< REAL_VALUES,COMPLEX_VALUES

   mode = i_mode

end subroutine

!>
!! This routine sets the matrix storage format.
!!
subroutine elsi_set_storage(i_storage)

   implicit none

   integer, intent(in) :: i_storage !< DENSE,CCS,CSC,CRS,CSR

   storage = i_storage

end subroutine

!>
!! This routine sets the parallel mode.
!!
subroutine elsi_set_parallel(i_parallel)

   implicit none

   integer, intent(in) :: i_parallel !< SINGLE_PROC,MULTI_PROC

   parallelism = i_parallel

   if(i_parallel == 0) then ! SINGLE_PROC
      n_l_rows = n_g_size
      n_l_cols = n_g_size
      n_b_rows = n_g_size
      n_b_cols = n_g_size
   endif

end subroutine

!>
!! Set MPI.
!!
subroutine elsi_set_mpi(mpi_comm_global_in)

   implicit none
   include "mpif.h"

   integer, intent(in) :: mpi_comm_global_in !< global mpi communicator

   if(parallelism == 1) then ! MULTI_PROC
      mpi_comm_global = mpi_comm_global_in

      call MPI_Comm_rank(mpi_comm_global,myid,mpierr)
      call MPI_Comm_size(mpi_comm_global,n_procs,mpierr)

      mpi_is_setup = .true.
   endif

end subroutine

!>
!! Set BLACS.
!!
subroutine elsi_set_blacs(blacs_ctxt_in,n_b_rows_in,n_b_cols_in,&
                          n_p_rows_in,n_p_cols_in)

   implicit none

   integer, intent(in) :: blacs_ctxt_in !< BLACS context
   integer, intent(in) :: n_b_rows_in   !< Block size
   integer, intent(in) :: n_b_cols_in   !< Block size
   integer, intent(in) :: n_p_rows_in   !< Number of processes in row
   integer, intent(in) :: n_p_cols_in   !< Number of processes in column

   integer :: i,i_row,i_col
   integer :: blacs_info

   character*40, parameter :: caller = "elsi_set_blacs"

   if(parallelism == 1) then ! MULTI_PROC
      blacs_ctxt = blacs_ctxt_in
      n_b_rows = n_b_rows_in
      n_b_cols = n_b_cols_in
      n_p_rows = n_p_rows_in
      n_p_cols = n_p_cols_in
      call blacs_pcoord(blacs_ctxt,myid,my_p_row,my_p_col)
      n_l_rows = numroc(n_g_size,n_b_rows,my_p_row,0,n_p_rows)
      n_l_cols = numroc(n_g_size,n_b_cols,my_p_col,0,n_p_cols)

      call descinit(sc_desc,n_g_size,n_g_size,n_b_rows,n_b_cols,0,0,&
                    blacs_ctxt,MAX(1,n_l_rows),blacs_info)

      call elsi_get_elpa_comms(mpi_comm_row,mpi_comm_col)

      ! Compute global-local mapping
      call elsi_allocate(local_row,n_g_size,"local_row",caller)
      call elsi_allocate(local_col,n_g_size,"local_col",caller)

      i_row = 0
      i_col = 0

      do i = 1,n_g_size
         if(MOD((i-1)/n_b_rows,n_p_rows) == my_p_row) then
            i_row = i_row+1
            local_row(i) = i_row
         endif
         if(MOD((i-1)/n_b_cols,n_p_cols) == my_p_col) then
            i_col = i_col+1
            local_col(i) = i_col
         endif
      enddo

      if(method == LIBOMM) then
         call ms_scalapack_setup(mpi_comm_global,n_p_rows,'r',n_b_rows,&
                                 icontxt=blacs_ctxt)
      endif

      blacs_is_setup = .true.
   endif

end subroutine

!>
!! This routine sets the sparsity pattern.
!!
subroutine elsi_set_sparsity(nnz_g_in,nnz_l_in,nnz_l_cols_in,&
                             row_ind_in,col_ptr_in)

   implicit none

   integer, intent(in) :: nnz_g_in      !< Global number of nonzeros
   integer, intent(in) :: nnz_l_in      !< Local number of nonzeros
   integer, intent(in) :: nnz_l_cols_in !< Local number of columns
   integer, target     :: row_ind_in(*) !< Row index
   integer, target     :: col_ptr_in(*) !< Column pointer

   character*40, parameter :: caller = "elsi_set_sparsity"

   nnz_g        = nnz_g_in
   nnz_l        = nnz_l_in
   nnz_l_pexsi  = nnz_l_cols_in

   call elsi_set_row_ind(row_ind_in)
   call elsi_set_col_ptr(col_ptr_in)

   sparsity_pattern_ready = .true.

end subroutine

!>
!! This routine gets the energy.
!!
subroutine elsi_get_energy(energy_out)

   implicit none

   real*8, intent(out) :: energy_out !< Energy of the system

   ! Only spin-nonpolarized case is supported now.
   real*8, parameter :: n_spin = 2.0d0
   integer :: i_state
   character*200 :: info_str

   character*40, parameter :: caller = "elsi_get_energy"

   select case (method)
      case (ELPA)
         energy_out = 0.0d0
         do i_state =1,n_states
            energy_out = energy_out+occ_elpa(i_state)*eval(i_state)
         enddo
      case (LIBOMM)
         energy_out = n_spin*total_energy
      case (PEXSI)
         energy_out = e_tot_H
      case (CHESS)
         call elsi_stop(" CHESS: not yet implemented! Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

   write(info_str,"(A,F15.5,A)") "  | Energy = ",energy_out," Ha"
   call elsi_statement_print(info_str)
   write(info_str,"(A,F15.5,A)") "  |        = ",energy_out*hartree," eV"
   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine finalizes ELSI.
!!
subroutine elsi_finalize()

   implicit none

   call elsi_cleanup()
   call elsi_print_timers()

   n_elsi_calls = 0

end subroutine

!=========================
! ELSI customize routines
!
!   elsi_customize
!   elsi_customize_omm
!   elsi_customize_pexsi
!   elsi_customize_elpa
!=========================

!>
!! This routine overrides ELSI default settings.
!!
subroutine elsi_customize(print_detail,unit_overlap,hartree_to_ev,numerical_zero,&
                          no_check_singularity,singularity_threshold,&
                          force_stop_singularity,broadening_scheme,broadening_width,&
                          mu_accuracy,mu_max_steps)

   implicit none

   logical, intent(in), optional :: print_detail           !< Print detailed info?
   logical, intent(in), optional :: unit_overlap           !< Is overlap matrix unit?
   real*8,  intent(in), optional :: hartree_to_ev          !< Conversion constant between Ha and eV
   real*8,  intent(in), optional :: numerical_zero         !< Threshold to define "zero"
   logical, intent(in), optional :: no_check_singularity   !< Do not perform singularity check of overlap
   real*8,  intent(in), optional :: singularity_threshold  !< Tolerance of overlap singularity
   logical, intent(in), optional :: force_stop_singularity !< Stop if overlap is singular
   integer, intent(in), optional :: broadening_scheme      !< Broadening method in chemical potential determination
   real*8,  intent(in), optional :: broadening_width       !< Broadening width in chemical potential determination
   real*8,  intent(in), optional :: mu_accuracy            !< Tolerance in chemical potential determination
   integer, intent(in), optional :: mu_max_steps           !< Maximum number of steps to find the chemical potential

   ! Print detailed ELSI information? [Default: .false.]
   if(PRESENT(print_detail)) &
      print_info = print_detail
   ! Is the overlap matrix unit? [Default: .false.]
   if(PRESENT(unit_overlap)) &
      overlap_is_unit = unit_overlap
   ! User-defined value for Hartree [Default: 27.211386, Codata 2015]
   if(PRESENT(hartree_to_ev)) &
      hartree = hartree_to_ev
   ! Threshold to define numerical zero [Default: 1d-13]
   if(PRESENT(numerical_zero)) &
      zero_threshold = numerical_zero
   ! Disable checking for overlap singularity? [Default: .false.]
   if(PRESENT(no_check_singularity)) &
      no_singularity_check = no_check_singularity
   ! Eigenfunctions of overlap matrix with eigenvalues smaller than
   ! this value will be removed to avoid singularity [Default: 1d-5]
   if(PRESENT(singularity_threshold)) &
      singularity_tolerance = singularity_threshold
   ! Always stop if overlap is singular? [Default: .false.]
   if(PRESENT(force_stop_singularity)) &
      stop_singularity = force_stop_singularity
   ! Broadening scheme to compute Fermi level [Default: GAUSSIAN]
   if(PRESENT(broadening_scheme)) &
      broaden_method = broadening_scheme
   ! Broadening width to compute Fermi level [Default: 1d-2]
   if(PRESENT(broadening_width)) &
      broaden_width = broadening_width
   ! Accuracy for chemical potential determination [Default: 1d-10]
   if(PRESENT(mu_accuracy)) &
      occ_tolerance = mu_accuracy
   ! Maximum steps to determine the chemical potential [Default: 100]
   if(PRESENT(mu_max_steps)) &
      max_mu_steps = mu_max_steps

end subroutine

!>
!! This routine overrides libOMM default settings.
!!
subroutine elsi_customize_omm(n_elpa_steps_omm,eigenspectrum_shift,&
                              omm_tolerance,use_pspblas)

   implicit none

   integer, intent(in), optional :: n_elpa_steps_omm    !< Number of ELPA steps before libOMM
   real*8,  intent(in), optional :: eigenspectrum_shift !< Eigenspectrum shift parameter
   real*8,  intent(in), optional :: omm_tolerance       !< Tolerance of minimization
   logical, intent(in), optional :: use_pspblas         !< Use pspBLAS sparse linear algebra?

   ! Number of ELPA steps
   if(PRESENT(n_elpa_steps_omm)) &
      n_elpa_steps = n_elpa_steps_omm
   ! Eigenspectrum shift parameter
   if(PRESENT(eigenspectrum_shift)) &
      eta = eigenspectrum_shift
   ! Tolerance for minimization
   if(PRESENT(omm_tolerance)) &
      min_tol = omm_tolerance
   ! Use pspBLAS sparse linear algebra?
   if(PRESENT(use_pspblas)) &
      use_psp = use_pspblas

   if(method .ne. LIBOMM) then
      call elsi_statement_print("  The chosen method is not libOMM."//&
                                " Ignore elsi_customize_omm call.")
   endif

end subroutine

!>
!! This routine overrides PEXSI default settings.
!!
subroutine elsi_customize_pexsi(temperature,gap,delta_E,n_poles,max_iteration,&
                                mu_min,mu_max,mu0,mu_inertia_tolerance,&
                                mu_inertia_expansion,mu_safeguard,n_electron_accuracy,&
                                matrix_type,is_symbolic_factorize,ordering,&
                                np_symbolic_factorize,verbosity)

   implicit none

   real(c_double), intent(in), optional :: temperature           !< Temperature, in the same unit as Hamiltonian
   real(c_double), intent(in), optional :: gap                   !< Spectral gap, can be set to 0 in most cases
   real(c_double), intent(in), optional :: delta_E               !< Upper bound for the spectral radius of S^(-1)H
   integer(c_int), intent(in), optional :: n_poles               !< Number of poles
   integer(c_int), intent(in), optional :: max_iteration         !< Maximum number of PEXSI iterations
   real(c_double), intent(in), optional :: mu_min                !< Lower bound of chemical potential
   real(c_double), intent(in), optional :: mu_max                !< Upper bound of chemical potential
   real(c_double), intent(in), optional :: mu0                   !< Initial guess of chemical potential
   real(c_double), intent(in), optional :: mu_inertia_tolerance  !< Tolerance of inertia counting
   real(c_double), intent(in), optional :: mu_inertia_expansion  !< Expansion step size in inertia counting
   real(c_double), intent(in), optional :: mu_safeguard          !< Safeguard to reinvoke inertia counting
   real(c_double), intent(in), optional :: n_electron_accuracy   !< Accuracy of number of electrons
   integer(c_int), intent(in), optional :: matrix_type           !< Type of input matrices
   integer(c_int), intent(in), optional :: is_symbolic_factorize !< Perform symbolic factorization?
   integer(c_int), intent(in), optional :: ordering              !< Ordering strategy for factorization and selected inversion
   integer(c_int), intent(in), optional :: np_symbolic_factorize !< Number of processors for ParMETIS, only used if ordering=0
   integer(c_int), intent(in), optional :: verbosity             !< Level of output info

   ! Temperature, in the same unit as H
   ! default: 0.0019 = 300K
   if(PRESENT(temperature)) &
      pexsi_options%temperature = temperature

   ! Spectral gap, can be set to 0 in most cases (default)
   if(PRESENT(gap)) &
      pexsi_options%gap = gap

   ! Upper bound for the spectral radius of S^(-1)H
   ! default: 10
   if(PRESENT(delta_E)) &
      pexsi_options%deltaE = delta_E

   ! Number of poles
   ! default: 40
   if(PRESENT(n_poles)) then
      pexsi_options%numPole = n_poles
      pexsi_n_poles_set = .true.
   endif

   ! Maximum number of PEXSI iterations after each inertia
   ! counting procedure
   ! default: 3
   if(PRESENT(max_iteration)) &
      pexsi_options%maxPEXSIIter = max_iteration

   ! From the second step, initial guess of mu is from previous step
   if(n_elsi_calls == 0) then
      ! Initial guess of mu
      ! default: 0.0
      if(PRESENT(mu0)) &
         pexsi_options%mu0 = mu0
   endif

   ! Initial guess of lower bound for mu
   ! default: -10.0
   if(PRESENT(mu_min)) &
      pexsi_options%muMin0 = mu_min

   ! Initial guess of upper bound for mu
   ! default: 10.0
   if(PRESENT(mu_max)) &
      pexsi_options%muMax0 = mu_max

   ! Stopping criterion in terms of the chemical potential
   ! for the inertia counting procedure
   ! default: 0.05
   if(PRESENT(mu_inertia_tolerance)) &
      pexsi_options%muInertiaTolerance = mu_inertia_tolerance

   ! If the chemical potential is not in the initial interval,
   ! the interval is expanded by this value
   ! default: 0.3
   if(PRESENT(mu_inertia_expansion)) &
      pexsi_options%muInertiaExpansion = mu_inertia_expansion

   ! Safeguard criterion in terms of the chemical potential to
   ! reinvoke the inertia counting procedure
   ! default: 0.05
   if(PRESENT(mu_safeguard)) &
      pexsi_options%muPEXSISafeGuard = mu_safeguard

   ! Stopping criterion of the PEXSI iteration in terms of the
   ! number of electrons compared to the exact number
   ! default: 0.01
   if(PRESENT(n_electron_accuracy)) then
      pexsi_options%numElectronPEXSITolerance = n_electron_accuracy
      if(n_electron_accuracy < 1.0d-2) then
         small_pexsi_tol = .true.
         final_pexsi_tol = n_electron_accuracy
      endif
   endif

   ! Type of input H and S matrices
   ! 0: real symmetric (default)
   ! 1: general complex
   if(PRESENT(matrix_type)) &
      pexsi_options%matrixType = matrix_type

   ! Whether to perform symbolic factorization
   ! default: 1
   if(PRESENT(is_symbolic_factorize)) &
      pexsi_options%isSymbolicFactorize = is_symbolic_factorize

   ! Ordering strategy for factorization and selected inversion
   ! 0: parallel ordering using ParMETIS
   ! 1: sequential ordering using METIS
   ! 2: multiple minimum degree ordering
   if(PRESENT(ordering)) &
      pexsi_options%ordering = ordering

   ! Number of processors for ParMETIS, only used if ordering=0
   if(PRESENT(np_symbolic_factorize)) &
      pexsi_options%npSymbFact = np_symbolic_factorize

   ! Level of output information
   ! 0: no output
   ! 1: basic output (default)
   ! 2: detailed output
   if(PRESENT(verbosity)) &
      pexsi_options%verbosity = verbosity

   if(method .ne. PEXSI) then
      call elsi_statement_print("  The chosen method is not PEXSI."//&
                                " Ignore elsi_customize_pexsi call.")
   endif

end subroutine

!>
!! This routine overrides ELPA default settings.
!!
subroutine elsi_customize_elpa(elpa_solver)

   implicit none

   integer, intent(in) :: elpa_solver !< Always use 1-stage or 2-stage solver

   if(elpa_solver == 1) then
      elpa_one_always = .true.
      elpa_two_always = .false.
   elseif(elpa_solver == 2) then
      elpa_one_always = .false.
      elpa_two_always = .true.
   else
      elpa_one_always = .false.
      elpa_two_always = .false.
   endif

   if(method .ne. ELPA) then
      call elsi_statement_print("  The chosen method is not ELPA."//&
                                " Ignore elsi_customize_elpa call.")
   endif

end subroutine

!=======================
! ELSI solvers
!
!   elsi_ev_real
!   elsi_ev_complex
!   elsi_dm_real
!   elsi_dm_complex
!   elsi_dm_real_sparse
!=======================

!>
!! This routine computes eigenvalues and eigenvectors.
!!
subroutine elsi_ev_real(H_in,S_in,e_val_out,e_vec_out)

   implicit none

   real*8, target :: H_in(n_l_rows,*)      !< Hamiltonian
   real*8, target :: S_in(n_l_rows,*)      !< Overlap
   real*8, target :: e_val_out(n_g_size)   !< Eigenvalues
   real*8, target :: e_vec_out(n_l_rows,*) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_real"

   ! Safety check
   call elsi_check()

   ! Update counter
   n_elsi_calls = n_elsi_calls+1

   ! Here the only supported method is ELPA
   select case (method)
      case (ELPA)
         ! REAL case
         call elsi_set_mode(REAL_VALUES)

         ! Set matrices
         call elsi_set_hamiltonian(H_in)
         if(.not.overlap_is_unit) then
            call elsi_set_overlap(S_in)
         endif
         call elsi_set_eigenvector(e_vec_out)
         call elsi_set_eigenvalue(e_val_out)

         ! Solve eigenvalue problem
         if(parallelism == 0) then ! Single-proc
            call elsi_solve_evp_elpa_sp()
         else  ! Multi-proc
            call elsi_solve_evp_elpa()
         endif

      case (LIBOMM)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",caller)
      case (PEXSI)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",caller)
      case (CHESS)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA solver to compute"//&
                        " eigenvalues and eigenvectors. Exiting...",caller)
   end select

end subroutine

!>
!! This routine computes eigenvalues and eigenvectors.
!!
subroutine elsi_ev_complex(H_in,S_in,e_val_out,e_vec_out)

   implicit none

   complex*16, target :: H_in(n_l_rows,*)      !< Hamiltonian
   complex*16, target :: S_in(n_l_rows,*)      !< Overlap
   real*8,     target :: e_val_out(n_g_size)   !< Eigenvalues
   complex*16, target :: e_vec_out(n_l_rows,*) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_complex"

   ! Safety check
   call elsi_check()

   ! Update counter
   n_elsi_calls = n_elsi_calls+1

   ! Here the only supported method is ELPA
   select case (method)
      case (ELPA)
         ! COMPLEX case
         call elsi_set_mode(COMPLEX_VALUES)

         ! Set matrices
         call elsi_set_hamiltonian(H_in)
         if(.not.overlap_is_unit) then
            call elsi_set_overlap(S_in)
         endif
         call elsi_set_eigenvector(e_vec_out)
         call elsi_set_eigenvalue(e_val_out)

         ! Solve eigenvalue problem
         if(parallelism==0) then ! Single-proc
            call elsi_solve_evp_elpa_sp()
         else ! Multi-proc
            call elsi_solve_evp_elpa()
         endif

      case (LIBOMM)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",caller)
      case (PEXSI)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",caller)
      case (CHESS)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors."//&
                        " Choose ELPA if necessary. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA solver to compute"//&
                        " eigenvalues and eigenvectors. Exiting...",caller)
   end select

end subroutine

!>
!! This routine computes density matrix.
!!
subroutine elsi_dm_real(H_in,S_in,D_out,energy_out)

   implicit none

   real*8, target      :: H_in(n_l_rows,*)  !< Hamiltonian
   real*8, target      :: S_in(n_l_rows,*)  !< Overlap
   real*8, target      :: D_out(n_l_rows,*) !< Density matrix
   real*8, intent(out) :: energy_out        !< Energy

   character*40, parameter :: caller = "elsi_dm_real"

   ! Safety check
   call elsi_check()

   ! Update counter
   n_elsi_calls = n_elsi_calls+1

   ! REAL case
   call elsi_set_mode(REAL_VALUES)

   ! Solve eigenvalue problem
   select case (method)
      case (ELPA)
         ! Set matrices
         if(.not.allocated(eval_elpa)) then
            call elsi_allocate(eval_elpa,n_g_size,"eval_elpa",caller)
         endif
         if(.not.allocated(evec_real_elpa)) then
            call elsi_allocate(evec_real_elpa,n_l_rows,n_l_cols,"evec_real_elpa",caller)
         endif

         call elsi_set_hamiltonian(H_in)
         if(.not.overlap_is_unit) then
            call elsi_set_overlap(S_in)
         endif
         call elsi_set_eigenvector(evec_real_elpa)
         call elsi_set_eigenvalue(eval_elpa)
         call elsi_set_density_matrix(D_out)

         call elsi_solve_evp_elpa()

         call elsi_compute_occ_elpa()
         call elsi_compute_dm_elpa()
         call elsi_get_energy(energy_out)

      case (LIBOMM)
         call elsi_print_omm_options()

         if(n_elsi_calls .le. n_elpa_steps) then ! Compute libOMM initial guess by ELPA
            call elsi_set_method(ELPA)

            ! Set matrices
            if(.not.allocated(eval_elpa)) then
               call elsi_allocate(eval_elpa,n_g_size,"eval_elpa",caller)
            endif
            if(.not.allocated(evec_real_elpa)) then
               call elsi_allocate(evec_real_elpa,n_l_rows,n_l_cols,"evec_real_elpa",caller)
            endif

            call elsi_set_hamiltonian(H_in)
            if(.not.overlap_is_unit) then
               call elsi_set_overlap(S_in)
            endif
            call elsi_set_eigenvector(evec_real_elpa)
            call elsi_set_eigenvalue(eval_elpa)
            call elsi_set_density_matrix(D_out)

            ! Solve by ELPA
            call elsi_solve_evp_elpa()
            call elsi_compute_occ_elpa()
            call elsi_compute_dm_elpa()
            call elsi_get_energy(energy_out)

            ! Switch back to libOMM here to guarantee elsi_customize_omm
            call elsi_set_method(LIBOMM)

         else ! ELPA is done
            ! Set matrices
            call elsi_set_hamiltonian(H_in)
            if(.not.overlap_is_unit) then
               call elsi_set_overlap(S_in)
            endif
            call elsi_set_density_matrix(D_out)

            if(.not.coeff_omm%is_initialized) then
               call m_allocate(coeff_omm,n_states,n_g_size,"pddbc")
            endif

            ! Initialize coefficient matrix with ELPA eigenvectors if possible
            if(n_elpa_steps > 0 .and. n_elsi_calls == n_elpa_steps+1) then
               ! D_out is used for temporary storage here
               D_out(1:n_l_rows,1:n_l_cols) = evec_real(1:n_l_rows,1:n_l_cols)
               ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
               call pdtran(n_g_size,n_g_size,1.0d0,D_out,1,1,sc_desc,0.0d0,evec_real,1,1,sc_desc)

               coeff_omm%dval(1:coeff_omm%iaux2(1),1:coeff_omm%iaux2(2)) = &
                  evec_real(1:coeff_omm%iaux2(1),1:coeff_omm%iaux2(2))

               ! ELPA matrices are no longer needed
               if(ASSOCIATED(ham_real))      nullify(ham_real)
               if(ASSOCIATED(ovlp_real))     nullify(ovlp_real)
               if(ASSOCIATED(evec_real))     nullify(evec_real)
               if(ASSOCIATED(den_mat))       nullify(den_mat)
               if(ASSOCIATED(eval))          nullify(eval)
               if(ALLOCATED(evec_real_elpa)) deallocate(evec_real_elpa)
               if(ALLOCATED(eval_elpa))      deallocate(eval_elpa)
               if(ALLOCATED(occ_elpa))       deallocate(occ_elpa)
            endif

            ! Continue the computation using libOMM
            call elsi_solve_evp_omm()

            den_mat_omm%dval = 2.0d0*den_mat_omm%dval
            call elsi_get_energy(energy_out)
         endif

      case (PEXSI)
         call elsi_print_pexsi_options()

         ! PEXSI may use different process grid to achieve
         ! the efficient 2-level parallelization
         call elsi_init_pexsi()

         ! Convert matrix format and distribution from BLACS to PEXSI
         call elsi_blacs_to_pexsi(H_in,S_in)

         ! Set matices
         call elsi_set_sparse_hamiltonian(ham_real_pexsi)
         if(.not.overlap_is_unit) then
            call elsi_set_sparse_overlap(ovlp_real_pexsi)
         endif
         call elsi_set_row_ind(row_ind_pexsi)
         call elsi_set_col_ptr(col_ptr_pexsi)
         if(.not.ALLOCATED(den_mat_pexsi)) then
            call elsi_allocate(den_mat_pexsi,nnz_l_pexsi,"den_mat_pexsi",caller)
         endif
         den_mat_pexsi = 0.0d0
         call elsi_set_sparse_density_matrix(den_mat_pexsi)

         call elsi_solve_evp_pexsi()

         ! Convert 1D block CCS sparse density matrix to 2D
         ! block-cyclic dense format
         call elsi_pexsi_to_blacs_dm(D_out)
         call elsi_get_energy(energy_out)

         if(ASSOCIATED(ham_real_ccs))   nullify(ham_real_ccs)
         if(ASSOCIATED(ovlp_real_ccs))  nullify(ovlp_real_ccs)
         if(ASSOCIATED(den_mat_ccs))    nullify(den_mat_ccs)
         if(ASSOCIATED(row_ind_ccs))    nullify(row_ind_ccs)
         if(ASSOCIATED(col_ptr_ccs))    nullify(col_ptr_ccs)
         if(ALLOCATED(ham_real_pexsi))  deallocate(ham_real_pexsi)
         if(ALLOCATED(ovlp_real_pexsi)) deallocate(ovlp_real_pexsi)
         if(ALLOCATED(den_mat_pexsi))   deallocate(den_mat_pexsi)
         if(ALLOCATED(e_den_mat_pexsi)) deallocate(e_den_mat_pexsi)
         if(ALLOCATED(f_den_mat_pexsi)) deallocate(f_den_mat_pexsi)
         if(ALLOCATED(row_ind_pexsi))   deallocate(row_ind_pexsi)
         if(ALLOCATED(col_ptr_pexsi))   deallocate(col_ptr_pexsi)

         call f_ppexsi_plan_finalize(pexsi_plan,pexsi_info)

      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case default
         call elsi_stop(" No supported solver has been chosen."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine computes density matrix.
!!
subroutine elsi_dm_complex(H_in,S_in,D_out,energy_out)

   implicit none

   complex*16, target  :: H_in(n_l_rows,*)  !< Hamiltonian
   complex*16, target  :: S_in(n_l_rows,*)  !< Overlap
   real*8,     target  :: D_out(n_l_rows,*) !< Density matrix
   real*8, intent(out) :: energy_out        !< Energy

   character*40, parameter :: caller = "elsi_dm_complex"

   call elsi_stop(" ELSI density matrix solver for complex case not yet available."//&
                  " Exiting...",caller)

   ! Safety check
   call elsi_check()

   ! Update counter
   n_elsi_calls = n_elsi_calls+1

   ! COMPLEX case
   call elsi_set_mode(COMPLEX_VALUES)

   ! Solve eigenvalue problem
   select case (method)
      case (ELPA)
         ! Set matrices
         if(.not.allocated(eval_elpa)) then
            call elsi_allocate(eval_elpa,n_g_size,"eval_elpa",caller)
         endif
         if(.not.allocated(evec_complex_elpa)) then
            call elsi_allocate(evec_complex_elpa,n_l_rows,n_l_cols,"evec_complex_elpa",caller)
         endif

         call elsi_set_hamiltonian(H_in)
         if(.not.overlap_is_unit) then
            call elsi_set_overlap(S_in)
         endif
         call elsi_set_eigenvector(evec_complex_elpa)
         call elsi_set_eigenvalue(eval_elpa)
         call elsi_set_density_matrix(D_out)

         call elsi_solve_evp_elpa()

         call elsi_compute_occ_elpa()
         call elsi_compute_dm_elpa()
         call elsi_get_energy(energy_out)

      case (LIBOMM)
         call elsi_print_omm_options()

         ! Set Hamiltonian and overlap matrices
         call elsi_set_hamiltonian(H_in)
         if(.not.overlap_is_unit) then
            call elsi_set_overlap(S_in)
         endif
         call elsi_set_density_matrix(D_out)

         if(.not.coeff_omm%is_initialized) then
            call m_allocate(coeff_omm,n_states,n_g_size,"pddbc")
         endif

         call elsi_solve_evp_omm()

         den_mat_omm%dval = 2.0d0*den_mat_omm%dval
         call elsi_get_energy(energy_out)

      case (PEXSI)
         call elsi_stop(" PEXSI not yet implemented. Exiting...",caller)
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case default
         call elsi_stop(" No supported solver has been chosen."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine computes density matrix.
!!
subroutine elsi_dm_real_sparse(H_in,S_in,D_out,energy_out)

   implicit none

   real*8,  target      :: H_in(*)    !< Hamiltonian
   real*8,  target      :: S_in(*)    !< Overlap
   real*8,  target      :: D_out(*)   !< Density matrix
   real*8,  intent(out) :: energy_out !< Energy

   character*40, parameter :: caller = "elsi_dm_real_sparse"

   ! Safety check
   call elsi_check()

   ! Update counter
   n_elsi_calls = n_elsi_calls+1

   ! REAL case
   call elsi_set_mode(REAL_VALUES)

   ! Solve eigenvalue problem
   select case (method)
      case (ELPA)
         call elsi_stop(" ELPA cannot handle sparse matrices."//&
                        " Exiting...",caller)
      case (LIBOMM)
         call elsi_stop(" libOMM cannot handle sparse matrices."//&
                        " Exiting...",caller)
      case (PEXSI)
         call elsi_print_pexsi_options()

         ! Set matrices
         call elsi_set_sparse_hamiltonian(H_in)
         if(.not.overlap_is_unit) then
            call elsi_set_sparse_overlap(S_in)
         endif
         call elsi_set_sparse_density_matrix(D_out)

         call elsi_init_pexsi()

         call elsi_solve_evp_pexsi()

         call elsi_get_energy(energy_out)

         call f_ppexsi_plan_finalize(pexsi_plan,pexsi_info)

      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case default
         call elsi_stop(" No supported solver has been chosen."//&
                        " Exiting...",caller)
   end select

end subroutine

end module ELSI
