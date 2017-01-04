!Copyright (c) 2016, ELSI consortium
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in the
!   documentation and/or other materials provided with the distribution.
! * Neither the name of the ELSI project nor the
!   names of its contributors may be used to endorse or promote products
!   derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

   ! Public routines
   public :: elsi_init            !< Initialize
   public :: elsi_set_method      !< Select solver
   public :: elsi_set_mpi         !< Set MPI from calling code
   public :: elsi_set_blacs       !< Set BLACS from calling code
   public :: elsi_customize       !< Override ELSI default
   public :: elsi_customize_elpa  !< Override ELPA default
   public :: elsi_customize_omm   !< Override libOMM default
   public :: elsi_customize_pexsi !< Override PEXSI default
   public :: elsi_ev_real         !< Compute eigenvalues and eigenvectors
   public :: elsi_ev_complex      !< Compute eigenvalues and eigenvectors
   public :: elsi_dm_real         !< Compute density matrix
   public :: elsi_dm_complex      !< Compute density matrix
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
!   elsi_get_energy
!   elsi_get_dm
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
   integer, intent(in) :: parallel_mode  !< SERIAL,PARALLEL
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
      n_states = NINT(n_electrons/2d0)
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

      call elsi_get_elpa_comms(mpi_comm_global,my_p_row,my_p_col,&
                               mpi_comm_row,mpi_comm_col)

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
!! This routine gets the energy.
!!
subroutine elsi_get_energy(energy_out)

   implicit none

   real*8, intent(out) :: energy_out !< Energy of the system

   ! Only spin-nonpolarized case is supported now.
   real*8, parameter :: n_spin = 2d0
   integer :: i_state
   character*200 :: info_str

   character*40, parameter :: caller = "elsi_get_energy"

   select case (method)
      case (ELPA)
         energy_out = 0d0
         do i_state =1,n_states
            energy_out = energy_out+occ_elpa(i_state)*eigenvalues(i_state)
         enddo
      case (LIBOMM)
         energy_out = n_spin*total_energy
      case (PEXSI)
         energy_out = e_tot_H
      case (CHESS)
         call elsi_stop(" CHESS: not yet implemented! Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting...",caller)
   end select

   write(info_str,"(A,F15.5,A)") "  | Energy = ",energy_out," Ha"
   call elsi_statement_print(info_str)
   write(info_str,"(A,F15.5,A)") "  |        = ",energy_out*hartree," eV"
   call elsi_statement_print(info_str)

end subroutine

!>
!! This routine gets the density matrix.
!!
subroutine elsi_get_dm(D_out)

   implicit none

   real*8, intent(out) :: D_out(n_l_rows,n_l_cols) !< Density matrix

   character*40, parameter :: caller = "elsi_get_dm"

   select case (method)
      case (ELPA)
         call elsi_stop(" ELPA needs to compute density matrix from eigenvectors. "//&
                        " Exiting...",caller)
      case (LIBOMM)
         ! Note that here the convention of density matrix used in libOMM is
         ! changed to the one used in ELPA and PEXSI.
         D_out = 2d0*D_omm%dval
      case (PEXSI)
         call elsi_pexsi_to_blacs_dm(D_out)
      case (CHESS)
         call elsi_stop(" CHESS: not yet implemented! Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine finalizes ELSI.
!!
subroutine elsi_finalize()

   implicit none
   include "mpif.h"

   if(parallelism == 1) then
      call MPI_Barrier(mpi_comm_global,mpierr)
   endif

   call elsi_deallocate_matrices()
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
                          mu_accuracy,no_check_singularity,singularity_threshold,&
                          force_stop_singularity,broadening_scheme,broadening_width)

   implicit none

   logical, intent(in), optional :: print_detail
   logical, intent(in), optional :: unit_overlap
   real*8,  intent(in), optional :: hartree_to_ev
   real*8,  intent(in), optional :: numerical_zero
   real*8,  intent(in), optional :: mu_accuracy
   logical, intent(in), optional :: no_check_singularity
   real*8,  intent(in), optional :: singularity_threshold
   logical, intent(in), optional :: force_stop_singularity
   integer, intent(in), optional :: broadening_scheme
   real*8,  intent(in), optional :: broadening_width

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
   ! Accuracy for chemical potential determination [Default: 1d-10]
   if(PRESENT(mu_accuracy)) &
      occ_tolerance = mu_accuracy
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

end subroutine

!>
!! This routine overrides libOMM default settings.
!!
subroutine elsi_customize_omm(n_elpa_steps_omm,scale_kinetic_energy,calculate_ed,&
                              eigenspectrum_shift,omm_tolerance)

   implicit none

   integer, intent(in), optional :: n_elpa_steps_omm
   real*8,  intent(in), optional :: scale_kinetic_energy
   logical, intent(in), optional :: calculate_ed
   real*8,  intent(in), optional :: eigenspectrum_shift
   real*8,  intent(in), optional :: omm_tolerance

   ! Number of ELPA steps
   if(PRESENT(n_elpa_steps_omm)) &
      n_elpa_steps = n_elpa_steps_omm
   ! Scaling of kinetic energy matrix
   if(PRESENT(scale_kinetic_energy)) &
      scale_kinetic = scale_kinetic_energy
   ! Calculate energy weigthed density matrix?
   if(PRESENT(calculate_ed)) &
      calc_ed = calculate_ed
   ! Eigenspectrum shift parameter
   if(PRESENT(eigenspectrum_shift)) &
      eta = eigenspectrum_shift
   ! Tolerance for minimization
   if(PRESENT(omm_tolerance)) &
      min_tol = omm_tolerance

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

   real(c_double), intent(in), optional :: temperature
   real(c_double), intent(in), optional :: gap
   real(c_double), intent(in), optional :: delta_E
   integer(c_int), intent(in), optional :: n_poles
   integer(c_int), intent(in), optional :: max_iteration
   real(c_double), intent(in), optional :: mu_min
   real(c_double), intent(in), optional :: mu_max
   real(c_double), intent(in), optional :: mu0
   real(c_double), intent(in), optional :: mu_inertia_tolerance
   real(c_double), intent(in), optional :: mu_inertia_expansion
   real(c_double), intent(in), optional :: mu_safeguard
   real(c_double), intent(in), optional :: n_electron_accuracy
   integer(c_int), intent(in), optional :: matrix_type
   integer(c_int), intent(in), optional :: is_symbolic_factorize
   integer(c_int), intent(in), optional :: ordering
   integer(c_int), intent(in), optional :: np_symbolic_factorize
   integer(c_int), intent(in), optional :: verbosity

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
   if(PRESENT(n_poles)) &
      pexsi_options%numPole = n_poles

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
      if(n_electron_accuracy < 1d-2) then
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

!===================
! ELSI solvers
!
!   elsi_ev_real
!   elsi_ev_complex
!   elsi_dm_real
!   elsi_dm_complex
!===================

!>
!! This routine computes eigenvalues and eigenvectors.
!!
subroutine elsi_ev_real(H_in,S_in,e_val_out,e_vec_out)

   implicit none

   real*8, target, intent(in) :: H_in(n_l_rows,n_l_cols)      !< Hamiltonian
   real*8, target, intent(in) :: S_in(n_l_rows,n_l_cols)      !< Overlap
   real*8, intent(out)        :: e_val_out(n_states)          !< Eigenvalues
   real*8, intent(out)        :: e_vec_out(n_l_rows,n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_real"

   ! Update counter
   n_elsi_calls = n_elsi_calls+1

   ! Here the only supported method is ELPA
   select case (method)
      case (ELPA)
         ! REAL case
         call elsi_set_mode(REAL_VALUES)

         ! Allocate matrices
         call elsi_allocate_matrices()

         ! Set Hamiltonian and Overlap matrices
         call elsi_set_hamiltonian(H_in)
         if(.not.overlap_is_unit) then
            call elsi_set_overlap(S_in)
         endif

         ! Solve eigenvalue problem
         if(parallelism == 0) then ! Single-proc
            call elsi_solve_evp_elpa_sp()
         else  ! Multi-proc
            call elsi_solve_evp_elpa()
         endif

         call elsi_get_eigenvalues(e_val_out)
         call elsi_get_eigenvectors(e_vec_out)

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
         call elsi_stop(" No supported method has been chosen."//&
                        " Please choose ELPA to compute eigenpairs."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine computes eigenvalues and eigenvectors.
!!
subroutine elsi_ev_complex(H_in,S_in,e_val_out,e_vec_out)

   implicit none

   complex*16, target, intent(in) :: H_in(n_l_rows,n_l_cols)      !< Hamiltonian
   complex*16, target, intent(in) :: S_in(n_l_rows,n_l_cols)      !< Overlap
   real*8, intent(out)            :: e_val_out(n_states)          !< Eigenvalues
   complex*16, intent(out)        :: e_vec_out(n_l_rows,n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_ev_complex"

   ! Update counter
   n_elsi_calls = n_elsi_calls+1

   ! Here the only supported method is ELPA
   select case (method)
      case (ELPA)
         ! COMPLEX case
         call elsi_set_mode(COMPLEX_VALUES)

         ! Allocate matrices
         call elsi_allocate_matrices()

         ! Set Hamiltonian and Overlap matrices
         call elsi_set_hamiltonian(H_in)
         if(.not.overlap_is_unit) then
            call elsi_set_overlap(S_in)
         endif

         ! Solve eigenvalue problem
         if(parallelism==0) then ! Single-proc
            call elsi_solve_evp_elpa_sp()
         else ! Multi-proc
            call elsi_solve_evp_elpa()
         endif

         call elsi_get_eigenvalues(e_val_out)
         call elsi_get_eigenvectors(e_vec_out)

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
         call elsi_stop(" No supported method has been chosen."//&
                        " Please choose ELPA to compute eigenpairs."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine computes density matrix.
!!
subroutine elsi_dm_real(H_in,S_in,D_out,energy_out)

   implicit none

   real*8,  target, intent(in)   :: H_in(n_l_rows,n_l_cols)  !< Hamiltonian
   real*8,  target, intent(in)   :: S_in(n_l_rows,n_l_cols)  !< Overlap
   real*8,  intent(out)          :: D_out(n_l_rows,n_l_cols) !< Density matrix
   real*8,  intent(out)          :: energy_out               !< Energy

   character*40, parameter :: caller = "elsi_dm_real"

   ! Update counter
   n_elsi_calls = n_elsi_calls+1

   ! REAL case
   call elsi_set_mode(REAL_VALUES)

   ! Allocation of matrices
   call elsi_allocate_matrices()

   ! Solve eigenvalue problem
   select case (method)
      case (ELPA)
         ! Set Hamiltonian and overlap matrices
         call elsi_set_hamiltonian(H_in)
         if(.not.overlap_is_unit) then
            call elsi_set_overlap(S_in)
         endif

         call elsi_solve_evp_elpa()

         call elsi_compute_occ_elpa()
         call elsi_compute_dm_elpa(D_out)
         call elsi_get_energy(energy_out)

      case (LIBOMM)
         if(overlap_is_unit) then
            call elsi_stop(" Unit overlap in libOMM not yet implemented."//&
                           " Exiting...",caller)
         endif

         if(MOD(NINT(n_electrons),2) /= 0) then
            call elsi_stop(" The current implementation of libOMM does not"//&
                           " work with fractional occupation numbers. This"//&
                           " means number of electrons in non-spin-polarized"//&
                           " system cannot be odd. Exiting...",caller)
         endif

         call elsi_print_omm_options()

         if(n_elsi_calls .le. n_elpa_steps) then ! Compute libOMM initial guess by ELPA
            call elsi_set_method(ELPA)

            ! Allocate ELPA matrices
            call elsi_allocate_matrices()

            ! Set Hamiltonian and overlap matrices
            call elsi_set_hamiltonian(H_in)
            if(.not.overlap_is_unit) then
               call elsi_set_overlap(S_in)
            endif

            ! Solve by ELPA
            call elsi_solve_evp_elpa()
            call elsi_compute_occ_elpa()
            call elsi_compute_dm_elpa(D_out)
            call elsi_get_energy(energy_out)

            ! Switch back to libOMM here to guarantee elsi_customize_omm
            call elsi_set_method(LIBOMM)

         else ! ELPA is done
            ! Set Hamiltonian and overlap matrices
            call elsi_set_hamiltonian(H_in)
            if(.not.overlap_is_unit) then
               call elsi_set_overlap(S_in)
            endif

            ! Initialize coefficient matrix with ELPA eigenvectors if possible
            if(n_elpa_steps > 0 .and. n_elsi_calls == n_elpa_steps+1) then
               ! D_elpa is used for temporary storage here
               D_elpa = C_real
               ! libOMM coefficient matrix is the transpose of ELPA eigenvectors
               call pdtran(n_g_size,n_g_size,1d0,D_elpa,1,1,sc_desc,0d0,C_real,1,1,sc_desc)

               Coeff_omm%dval(1:Coeff_omm%iaux2(1),1:Coeff_omm%iaux2(2)) = &
                  C_real(1:Coeff_omm%iaux2(1),1:Coeff_omm%iaux2(2))

               ! ELPA matrices are no longer needed at this point
               if(ASSOCIATED(H_real))     nullify(H_real)
               if(ASSOCIATED(S_real))     nullify(S_real)
               if(ALLOCATED(C_real))      deallocate(C_real)
               if(ALLOCATED(eigenvalues)) deallocate(eigenvalues)
               if(ALLOCATED(D_elpa))      deallocate(D_elpa)
               if(ALLOCATED(occ_elpa))    deallocate(occ_elpa)
            endif

            ! Continue the computation using libOMM
            call elsi_solve_evp_omm()

            call elsi_get_dm(D_out)
            call elsi_get_energy(energy_out)
         endif

      case (PEXSI)
         if(n_g_size < n_procs) then
            call elsi_stop(" The (global) size of matrix is too small for"//&
                           " this number of processes. Exiting...",caller)
         endif

         call elsi_print_pexsi_options()

         ! PEXSI may use different process grid to achieve
         ! the efficient 2-level parallelization
         call elsi_init_pexsi()

         ! Convert matrix format and distribution from BLACS to PEXSI
         call elsi_blacs_to_pexsi(H_in,S_in)

         call elsi_solve_evp_pexsi()

         ! Convert 1D block CCS sparse density matrix to 2D
         ! block-cyclic dense format
         call elsi_get_dm(D_out)
         call elsi_get_energy(energy_out)

         if(ALLOCATED(H_real_pexsi))  deallocate(H_real_pexsi)
         if(ALLOCATED(S_real_pexsi))  deallocate(S_real_pexsi)
         if(ALLOCATED(D_pexsi))       deallocate(D_pexsi)
         if(ALLOCATED(ED_pexsi))      deallocate(ED_pexsi)
         if(ALLOCATED(FD_pexsi))      deallocate(FD_pexsi)
         if(ALLOCATED(row_ind_pexsi)) deallocate(row_ind_pexsi)
         if(ALLOCATED(col_ptr_pexsi)) deallocate(col_ptr_pexsi)

         call f_ppexsi_plan_finalize(pexsi_plan,pexsi_info)

      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case default
         call elsi_stop(" No supported method has been chosen."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine computes density matrix.
!!
subroutine elsi_dm_complex(H_in,S_in,D_out,energy_out)

   implicit none

   complex*16, target, intent(in) :: H_in(n_l_rows,n_l_cols)  !< Hamiltonian
   complex*16, target, intent(in) :: S_in(n_l_rows,n_l_cols)  !< Overlap
   real*8, intent(out)            :: D_out(n_l_rows,n_l_cols) !< Density matrix
   real*8, intent(out)            :: energy_out               !< Energy

   character*40, parameter :: caller = "elsi_dm_complex"

   call elsi_stop(" ELSI density matrix solver for complex case not yet available."//&
                  " Exiting...",caller)

   ! Update counter
   n_elsi_calls = n_elsi_calls+1

   ! COMPLEX case
   call elsi_set_mode(COMPLEX_VALUES)

   ! Allocation of matrices
   call elsi_allocate_matrices()

   ! Solve eigenvalue problem
   select case (method)
      case (ELPA)
         ! Set Hamiltonian and overlap matrices
         call elsi_set_hamiltonian(H_in)
         if(.not.overlap_is_unit) then
            call elsi_set_overlap(S_in)
         endif

         call elsi_solve_evp_elpa()

         call elsi_compute_occ_elpa()
         call elsi_compute_dm_elpa(D_out)
         call elsi_get_energy(energy_out)

      case (LIBOMM)
         if(overlap_is_unit) then
            call elsi_stop(" Unit overlap in libOMM not yet implemented."//&
                           " Exiting...",caller)
         endif

         if(MOD(NINT(n_electrons),2) /= 0) then
            call elsi_stop(" The current implementation of libOMM does not"//&
                           " work with fractional occupation numbers. This"//&
                           " means number of electrons in non-spin-polarized"//&
                           " system cannot be odd. Exiting...",caller)
         endif

         call elsi_print_omm_options()

         ! Set Hamiltonian and overlap matrices
         call elsi_set_hamiltonian(H_in)
         if(.not.overlap_is_unit) then
            call elsi_set_overlap(S_in)
         endif

         call elsi_solve_evp_omm()
         call elsi_get_dm(D_out)
         call elsi_get_energy(energy_out)

      case (PEXSI)
         call elsi_stop(" PEXSI not yet implemented. Exiting...",caller)
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case default
         call elsi_stop(" No supported method has been chosen."//&
                        " Exiting...",caller)
   end select

end subroutine

end module ELSI
