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
!! This module provides routines for modifying ELSI and solver parameters.
!!
module ELSI_MUTATOR

   use ELSI_CONSTANTS, only: ELPA,LIBOMM,PEXSI,CHESS,SIPS,REAL_VALUES,&
                             COMPLEX_VALUES,UNSET
   use ELSI_DATATYPE
   use ELSI_ELPA, only: elsi_compute_edm_elpa
   use ELSI_MATCONV, only: elsi_pexsi_to_blacs_dm,elsi_blacs_to_sips_dm
   use ELSI_OMM, only: elsi_compute_edm_omm
   use ELSI_PEXSI, only: elsi_compute_edm_pexsi
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS

   implicit none

   private

   ! Mutator
   public :: elsi_set_output
   public :: elsi_set_unit_ovlp
   public :: elsi_set_zero_def
   public :: elsi_set_sing_check
   public :: elsi_set_sing_tol
   public :: elsi_set_sing_stop
   public :: elsi_set_uplo
   public :: elsi_set_elpa_solver
   public :: elsi_set_omm_flavor
   public :: elsi_set_omm_n_elpa
   public :: elsi_set_omm_tol
   public :: elsi_set_omm_ev_shift
   public :: elsi_set_omm_psp
   public :: elsi_set_pexsi_n_mu
   public :: elsi_set_pexsi_n_pole
   public :: elsi_set_pexsi_np_per_pole
   public :: elsi_set_pexsi_np_symbo
   public :: elsi_set_pexsi_temp
   public :: elsi_set_pexsi_gap
   public :: elsi_set_pexsi_delta_e
   public :: elsi_set_pexsi_mu_min
   public :: elsi_set_pexsi_mu_max
   public :: elsi_set_pexsi_inertia_tol
   public :: elsi_set_chess_erf_decay
   public :: elsi_set_chess_erf_decay_min
   public :: elsi_set_chess_erf_decay_max
   public :: elsi_set_chess_ev_ham_min
   public :: elsi_set_chess_ev_ham_max
   public :: elsi_set_chess_ev_ovlp_min
   public :: elsi_set_chess_ev_ovlp_max
   public :: elsi_set_sips_slice_type
   public :: elsi_set_sips_n_slice
   public :: elsi_set_sips_left_bound
   public :: elsi_set_sips_slice_buf
   public :: elsi_set_mu_broaden_scheme
   public :: elsi_set_mu_broaden_width
   public :: elsi_set_mu_tol
   public :: elsi_set_mu_spin_degen
   public :: elsi_get_pexsi_mu_min
   public :: elsi_get_pexsi_mu_max
   public :: elsi_get_ovlp_sing
   public :: elsi_get_n_sing
   public :: elsi_get_mu
   public :: elsi_get_edm_real
   public :: elsi_get_edm_complex
   public :: elsi_get_edm_real_sparse
   public :: elsi_customize
   public :: elsi_customize_elpa
   public :: elsi_customize_omm
   public :: elsi_customize_pexsi
   public :: elsi_customize_sips
   public :: elsi_customize_mu

contains

!>
!! This routine overrides ELSI default settings.
!!
subroutine elsi_customize(elsi_h,print_detail,overlap_is_unit,zero_threshold,&
              no_singularity_check,singularity_tolerance,stop_singularity,uplo)

   implicit none

   type(elsi_handle), intent(inout)        :: elsi_h                !< Handle
   logical,           intent(in), optional :: print_detail          !< Print detailed info?
   logical,           intent(in), optional :: overlap_is_unit       !< Is overlap matrix unit?
   real(kind=r8),     intent(in), optional :: zero_threshold        !< Threshold to define "zero"
   logical,           intent(in), optional :: no_singularity_check  !< Do not perform singularity check
   real(kind=r8),     intent(in), optional :: singularity_tolerance !< Tolerance of overlap singularity
   logical,           intent(in), optional :: stop_singularity      !< Stop if overlap is singular
   integer,           intent(in), optional :: uplo                  !< Is input upper/lower triangular?

   character*40, parameter :: caller = "elsi_customize"

   call elsi_check_handle(elsi_h,caller)

   ! Is the overlap matrix unit? [Default: .false.]
   if(present(overlap_is_unit)) then
      elsi_h%ovlp_is_unit = overlap_is_unit
   endif

   ! Threshold to define numerical zero [Default: 1e-15_r8]
   if(present(zero_threshold)) then
      elsi_h%zero_threshold = zero_threshold
   endif

   ! Disable checking for overlap singularity? [Default: .false.]
   if(present(no_singularity_check)) then
      elsi_h%no_sing_check = no_singularity_check
   endif

   ! Eigenfunctions of overlap matrix with eigenvalues smaller than
   ! this value will be removed to avoid singularity [Default: 1e-5_r8]
   if(present(singularity_tolerance)) then
      elsi_h%sing_tol = singularity_tolerance
   endif

   ! Always stop if overlap is singular? [Default: .false.]
   if(present(stop_singularity)) then
      elsi_h%stop_sing = stop_singularity
   endif

   ! Is the input matrices upper or lower triangular?
   ! 0: FULL_MAT, 1: UT_MAT, 2: LT_MAT [Default: 0]
   if(present(uplo)) then
      elsi_h%uplo = uplo
   endif

end subroutine

!>
!! This routine overrides libOMM default settings.
!!
subroutine elsi_customize_omm(elsi_h,n_elpa_steps,omm_flavor,eigen_shift,&
              omm_tolerance,use_pspblas,omm_output)

   implicit none

   type(elsi_handle), intent(inout)        :: elsi_h        !< Handle
   integer(kind=i4),  intent(in), optional :: n_elpa_steps  !< Number of ELPA steps
   integer(kind=i4),  intent(in), optional :: omm_flavor    !< OMM method
   real(kind=r8),     intent(in), optional :: eigen_shift   !< Eigenspectrum shift
   real(kind=r8),     intent(in), optional :: omm_tolerance !< Tolerance of minimization
   logical,           intent(in), optional :: use_pspblas   !< Use PSP sparse linear algebra?
   logical,           intent(in), optional :: omm_output    !< Output details?

   character*40, parameter :: caller = "elsi_customize_omm"

   call elsi_check_handle(elsi_h,caller)

   ! Number of ELPA steps [Default: 6]
   if(present(n_elpa_steps)) then
      elsi_h%n_elpa_steps = n_elpa_steps
   endif

   ! How to perform orbital minimization? [Default: 0]
   if(present(omm_flavor)) then
      elsi_h%omm_flavor = omm_flavor
   endif

   ! Eigenspectrum shift parameter [Default: 0.0_r8]
   if(present(eigen_shift)) then
      elsi_h%eta = eigen_shift
   endif

   ! Tolerance for minimization [Default: 1.0e-10_r8]
   if(present(omm_tolerance)) then
      elsi_h%min_tol = omm_tolerance
   endif

   ! Use pspBLAS sparse linear algebra? [Default: .false.]
   if(present(use_pspblas)) then
      elsi_h%use_psp = .false.
   endif

end subroutine

!>
!! This routine overrides PEXSI default settings.
!!
subroutine elsi_customize_pexsi(elsi_h,temperature,gap,delta_e,n_poles,&
              n_procs_per_pole,max_iteration,mu_min,mu_max,mu0,&
              mu_inertia_tolerance,mu_inertia_expansion,mu_safeguard,&
              n_electron_accuracy,matrix_type,is_symbolic_factorize,&
              ordering,np_symbolic_factorize,verbosity)

   implicit none

   type(elsi_handle), intent(inout)        :: elsi_h                !< Handle
   real(kind=r8),     intent(in), optional :: temperature           !< Temperature
   real(kind=r8),     intent(in), optional :: gap                   !< Spectral gap
   real(kind=r8),     intent(in), optional :: delta_e               !< Upper bound of spectral radius of S^(-1)H
   integer(kind=i4),  intent(in), optional :: n_poles               !< Number of poles
   integer(kind=i4),  intent(in), optional :: n_procs_per_pole      !< Number of processes for one pole
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
   integer(kind=i4),  intent(in), optional :: ordering              !< Ordering strategy
   integer(kind=i4),  intent(in), optional :: np_symbolic_factorize !< Number of processes for symbolic factorization
   integer(kind=i4),  intent(in), optional :: verbosity             !< Level of output info

   character*40, parameter :: caller = "elsi_customize_pexsi"

   call elsi_check_handle(elsi_h,caller)

   ! Temperature, in the same unit as H
   if(present(temperature)) then
      elsi_h%pexsi_options%temperature = temperature
   endif

   ! Spectral gap, can be set to 0 in most cases [Default: 0.0_r8]
   if(present(gap)) then
      elsi_h%pexsi_options%gap = gap
   endif

   ! Upper bound for the spectral radius of S^(-1)H [Default: 10.0_r8]
   if(present(delta_e)) then
      elsi_h%pexsi_options%deltaE = delta_e
   endif

   ! Number of poles [Default: 20]
   if(present(n_poles)) then
      elsi_h%pexsi_options%numPole = n_poles
   endif

   ! Number of processes for one pole
   if(present(n_procs_per_pole)) then
      if(mod(elsi_h%n_procs,n_procs_per_pole) == 0) then
         elsi_h%n_p_per_pole = n_procs_per_pole
      else
         call elsi_stop("  The total number of MPI tasks must be a"//&
                        " multiple of the number of MPI tasks per"//&
                        " pole. Exiting...",elsi_h,caller)
      endif
   endif

   ! Initial guess of lower bound for mu [Default: -10.0_r8]
   if(present(mu_min)) then
      elsi_h%pexsi_options%muMin0 = mu_min
   endif

   ! Initial guess of upper bound for mu [Default: 10.0_r8]
   if(present(mu_max)) then
      elsi_h%pexsi_options%muMax0 = mu_max
   endif

   ! Stopping criterion in terms of the chemical potential
   ! for the inertia counting procedure [Default: 0.05_r8]
   if(present(mu_inertia_tolerance)) then
      elsi_h%pexsi_options%muInertiaTolerance = mu_inertia_tolerance
   endif

   ! Whether to perform symbolic factorization [Default: 1]
   if(present(is_symbolic_factorize)) then
      elsi_h%pexsi_options%isSymbolicFactorize = is_symbolic_factorize
   endif

   ! Ordering strategy for factorization and selected inversion [Default: 0]
   ! 0: parallel ordering using ParMETIS
   ! 1: sequential ordering using METIS
   ! 2: multiple minimum degree ordering
   if(present(ordering)) then
      elsi_h%pexsi_options%ordering = ordering
   endif

   ! Number of processes for ParMETIS, only used if ordering=0 [Default: 1]
   if(present(np_symbolic_factorize)) then
      elsi_h%pexsi_options%npSymbFact = np_symbolic_factorize
   endif

end subroutine

!>
!! This routine overrides ELPA default settings.
!!
subroutine elsi_customize_elpa(elsi_h,elpa_solver,elpa_output)

   implicit none

   type(elsi_handle), intent(inout)        :: elsi_h      !< Handle
   integer(kind=i4),  intent(in), optional :: elpa_solver !< 1-stage or 2-stage solver?
   logical,           intent(in), optional :: elpa_output !< Output details?

   character*40, parameter :: caller = "elsi_customize_elpa"

   call elsi_check_handle(elsi_h,caller)

   ! 1-stage or 2-stage solver? [Default: 2]
   if(present(elpa_solver)) then
      elsi_h%elpa_solver = elpa_solver
   endif

end subroutine

!>
!! This routine overrides SIPs default settings.
!!
subroutine elsi_customize_sips(elsi_h,slicing_method,n_slices,inertia_option,&
              unbound,slice_buffer)

   implicit none

   type(elsi_handle), intent(inout)        :: elsi_h         !< Handle
   integer(kind=i4),  intent(in), optional :: slicing_method !< Method of slicing
   integer(kind=i4),  intent(in), optional :: n_slices       !< Number of slices
   integer(kind=i4),  intent(in), optional :: inertia_option !< Inertia counting before solve?
   integer(kind=i4),  intent(in), optional :: unbound        !< Bound the left side of the interval?
   real(kind=r8),     intent(in), optional :: slice_buffer   !< Small buffer to expand the interval

   character*40, parameter :: caller = "elsi_customize_sips"

   ! Method of slicing [Default: 3]
   ! 0: Equally spaced subintervals
   ! 1: K-meaans after equally spaced subintervals
   ! 2: Equally populated subintervals
   ! 3: K-means after equally populated subintervals
   if(present(slicing_method)) then
      elsi_h%slicing_method = slicing_method
   endif

   ! Number of slices
   if(present(n_slices)) then
      if(mod(elsi_h%n_procs,n_slices) == 0) then
         elsi_h%n_slices = n_slices
         elsi_h%n_p_per_slice = elsi_h%n_procs/elsi_h%n_slices
      else
         call elsi_stop("  The total number of MPI tasks must be"//&
                        " a multiple of the number of slices."//&
                        " Exiting...",elsi_h,caller)
      endif
   endif

   ! Perform inertia computations before solve? [Default: 1]
   if(present(inertia_option)) then
      elsi_h%inertia_option = inertia_option
   endif

   ! Bound the left side of the interval [Default: 0]
   ! 0: Bound interval
   ! 1: -Infinity
   if(present(unbound)) then
      elsi_h%unbound = unbound
   endif

   ! Small buffer to expand the interval [Default: 0.1_r8]
   if(present(slice_buffer)) then
      elsi_h%slice_buffer = slice_buffer
   endif

end subroutine

!>
!! This routine overrides ELSI default settings for the chemical potential
!! determination module.
!!
subroutine elsi_customize_mu(elsi_h,broadening_scheme,broadening_width,&
              occ_accuracy,mu_max_steps,spin_degeneracy)

   implicit none

   type(elsi_handle), intent(inout)        :: elsi_h            !< Handle
   integer(kind=i4),  intent(in), optional :: broadening_scheme !< Broadening method
   real(kind=r8),     intent(in), optional :: broadening_width  !< Broadening width
   real(kind=r8),     intent(in), optional :: occ_accuracy      !< Accuracy in electron count
   integer(kind=i4),  intent(in), optional :: mu_max_steps      !< Maximum number of bisection steps
   real(kind=r8),     intent(in), optional :: spin_degeneracy   !< Spin degeneracy

   character*40, parameter :: caller = "elsi_customize_mu"

   call elsi_check_handle(elsi_h,caller)

   ! Broadening scheme to compute Fermi level [Default: GAUSSIAN]
   if(present(broadening_scheme)) then
      elsi_h%broaden_scheme = broadening_scheme
   endif

   ! Broadening width to compute Fermi level [Default: 1e-2_r8]
   if(present(broadening_width)) then
      elsi_h%broaden_width = broadening_width
   endif

   ! Accuracy for chemical potential determination [Default: 1e-10_r8]
   if(present(occ_accuracy)) then
      elsi_h%occ_tolerance = occ_accuracy
   endif

   ! Maximum steps to determine the chemical potential [Default: 100]
   if(present(mu_max_steps)) then
      elsi_h%max_mu_steps = mu_max_steps
   endif

   ! Spin degeneracy [Default: 2.0_r8/n_spin]
   if(present(spin_degeneracy)) then
      elsi_h%spin_degen = spin_degeneracy
      elsi_h%spin_is_set = .true.
   endif

end subroutine

!>
!! This routine sets ELSI output level.
!!
subroutine elsi_set_output(elsi_h,out_level)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h    !< Handle
   integer(kind=i4),  intent(in)    :: out_level !< Output level of ELSI

   character*40, parameter :: caller = "elsi_set_output"

   call elsi_check_handle(elsi_h,caller)

   if(out_level <= 0) then
      print_info = .false.
      print_mem  = .false.
      elsi_h%omm_output = .false.
      elsi_h%pexsi_options%verbosity = 0
      elsi_h%elpa_output = .false.
   elseif(out_level == 1) then
      print_info = .true.
      print_mem  = .false.
      elsi_h%omm_output = .false.
      elsi_h%pexsi_options%verbosity = 0
      elsi_h%elpa_output = .false.
   elseif(out_level == 2) then
      print_info = .true.
      print_mem  = .false.
      elsi_h%omm_output = .true.
      elsi_h%pexsi_options%verbosity = 2
      elsi_h%elpa_output = .true.
   else
      print_info = .true.
      print_mem  = .true.
      elsi_h%omm_output = .true.
      elsi_h%pexsi_options%verbosity = 2
      elsi_h%elpa_output = .true.
   endif

end subroutine

!>
!! This routine sets the overlap matrix to be identity.
!!
subroutine elsi_set_unit_ovlp(elsi_h,unit_ovlp)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h    !< Handle
   integer(kind=i4),  intent(in)    :: unit_ovlp !< Overlap is an identity matrix?

   character*40, parameter :: caller = "elsi_set_unit_ovlp"

   call elsi_check_handle(elsi_h,caller)

   if(unit_ovlp == 0) then
      elsi_h%ovlp_is_unit = .false.
   else
      elsi_h%ovlp_is_unit = .true.
   endif

end subroutine

!>
!! This routine sets the threshold to define "zero".
!!
subroutine elsi_set_zero_def(elsi_h,zero_def)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h   !< Handle
   real(kind=r8),     intent(in)    :: zero_def !< Numbers smaller than this will be discarded

   character*40, parameter :: caller = "elsi_set_zero_def"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%zero_threshold = zero_def

end subroutine

!>
!! This routine switches on/off the singularity check of the overlap matrix.
!!
subroutine elsi_set_sing_check(elsi_h,sing_check)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h     !< Handle
   integer(kind=i4),  intent(in)    :: sing_check !< Perform singularity check?

   character*40, parameter :: caller = "elsi_set_sing_check"

   call elsi_check_handle(elsi_h,caller)

   if(sing_check == 0) then
      elsi_h%no_sing_check = .true.
   else
      elsi_h%no_sing_check = .false.
   endif

end subroutine

!>
!! This routine sets the tolerance of the singularity check.
!!
subroutine elsi_set_sing_tol(elsi_h,sing_tol)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h   !< Handle
   real(kind=r8),     intent(in)    :: sing_tol !< Singularity tolerance

   character*40, parameter :: caller = "elsi_set_sing_tol"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%sing_tol = sing_tol

end subroutine

!>
!! This routine sets whether to stop in case of singular overlap matrix.
!!
subroutine elsi_set_sing_stop(elsi_h,sing_stop)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h    !< Handle
   integer(kind=i4),  intent(in)    :: sing_stop !< Stop if overlap is singular

   character*40, parameter :: caller = "elsi_set_sing_stop"

   call elsi_check_handle(elsi_h,caller)

   if(sing_stop == 0) then
      elsi_h%stop_sing = .false.
   else
      elsi_h%stop_sing = .true.
   endif

end subroutine

!>
!! This routine sets the input matrices to be full, upper triangular, or
!! lower triangular.
!!
subroutine elsi_set_uplo(elsi_h,uplo)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   integer(kind=i4),  intent(in)    :: uplo   !< Input matrix triangular?

   character*40, parameter :: caller = "elsi_set_uplo"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%uplo = uplo

end subroutine

!>
!! This routine sets the ELPA solver.
!!
subroutine elsi_set_elpa_solver(elsi_h,elpa_solver)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h      !< Handle
   integer(kind=i4),  intent(in)    :: elpa_solver !< Which ELPA solver?

   character*40, parameter :: caller = "elsi_set_elpa_solver"

   call elsi_check_handle(elsi_h,caller)

   if((elpa_solver < 1) .or. (elpa_solver > 2)) then
      call elsi_stop("  Invalid choice of elpa_solver. Please"//&
              " choose 1 (1-stage solver) or 2 (2-stage solver)."//&
              " Exiting...",elsi_h,caller)
   endif

   elsi_h%elpa_solver = elpa_solver

end subroutine

!>
!! This routine sets the flavor of libOMM.
!!
subroutine elsi_set_omm_flavor(elsi_h,omm_flavor)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h     !< Handle
   integer(kind=i4),  intent(in)    :: omm_flavor !< Which libOMM flavor?

   character*40, parameter :: caller = "elsi_set_omm_flavor"

   call elsi_check_handle(elsi_h,caller)

   if((omm_flavor /= 0) .and. (omm_flavor /= 2)) then
      call elsi_stop("  Invalid choice of omm_flavor. Please"//&
              " choose 0 (basic flavor) or 2 (Cholesky flavor)."//&
              " Exiting...",elsi_h,caller)
   endif

   elsi_h%omm_flavor = omm_flavor

end subroutine

!>
!! This routine sets the number of ELPA steps when using libOMM.
!!
subroutine elsi_set_omm_n_elpa(elsi_h,n_elpa)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   integer(kind=i4),  intent(in)    :: n_elpa !< Number of ELPA steps before libOMM

   character*40, parameter :: caller = "elsi_set_n_elpa"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%n_elpa_steps = n_elpa

end subroutine

!>
!! This routine sets the tolerance of OMM minimization.
!!
subroutine elsi_set_omm_tol(elsi_h,min_tol)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h  !< Handle
   real(kind=r8),     intent(in)    :: min_tol !< Tolerance of OMM minimization

   character*40, parameter :: caller = "elsi_set_omm_tol"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%min_tol = min_tol

end subroutine

!>
!! This routine sets the shift of the eigenspectrum.
!!
subroutine elsi_set_omm_ev_shift(elsi_h,ev_shift)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h   !< Handle
   real(kind=r8),     intent(in)    :: ev_shift !< Shift of the eigenspectrum

   character*40, parameter :: caller = "elsi_set_omm_ev_shift"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%eta = ev_shift

end subroutine

!>
!! This routine switches on/off the matrix multiplications using PSP.
!!
subroutine elsi_set_omm_psp(elsi_h,use_psp)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h  !< Handle
   integer(kind=i4),  intent(in)    :: use_psp !< Use pspBLAS?

   character*40, parameter :: caller = "elsi_set_omm_psp"

   call elsi_check_handle(elsi_h,caller)

   if(use_psp == 0) then
      elsi_h%use_psp = .false.
   else
      elsi_h%use_psp = .false.
   endif

end subroutine

!>
!! This routine sets the number of mu points when using PEXSI driver 2.
!!
subroutine elsi_set_pexsi_n_mu(elsi_h,n_mu)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   integer(kind=i4),  intent(in)    :: n_mu   !< Number of mu points

   character*40, parameter :: caller = "elsi_set_pexsi_n_mu"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%pexsi_options%nPoints = n_mu

end subroutine

!>
!! This routine sets the number of poles in the pole expansion.
!!
subroutine elsi_set_pexsi_n_pole(elsi_h,n_pole)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   integer(kind=i4),  intent(in)    :: n_pole !< Number of poles

   character*40, parameter :: caller = "elsi_set_pexsi_n_pole"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%pexsi_options%numPole = n_pole

end subroutine

!>
!! This routine sets the number of MPI tasks assigned for one pole.
!!
subroutine elsi_set_pexsi_np_per_pole(elsi_h,np_per_pole)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h      !< Handle
   integer(kind=i4),  intent(in)    :: np_per_pole !< Number of processes per pole

   character*40, parameter :: caller = "elsi_set_pexsi_np_per_pole"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%n_p_per_pole = np_per_pole

end subroutine

!>
!! This routine sets the number of MPI tasks for the symbolic factorization.
!!
subroutine elsi_set_pexsi_np_symbo(elsi_h,np_symbo)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h   !< Handle
   integer(kind=i4),  intent(in)    :: np_symbo !< Number of processes for symbolic factorization

   character*40, parameter :: caller = "elsi_set_pexsi_np_symbo"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%pexsi_options%npSymbFact = np_symbo

end subroutine

!>
!! This routine sets the temperature parameter in PEXSI.
!!
subroutine elsi_set_pexsi_temp(elsi_h,temp)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(in)    :: temp   !< Temperature

   character*40, parameter :: caller = "elsi_set_pexsi_temp"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%pexsi_options%temperature = temp

end subroutine

!>
!! This routine sets the spectral gap in PEXSI.
!!
subroutine elsi_set_pexsi_gap(elsi_h,gap)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(in)    :: gap    !< Gap

   character*40, parameter :: caller = "elsi_set_pexsi_gap"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%pexsi_options%gap = gap

end subroutine

!>
!! This routine sets the spectrum width in PEXSI.
!!
subroutine elsi_set_pexsi_delta_e(elsi_h,delta_e)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h  !< Handle
   real(kind=r8),     intent(in)    :: delta_e !< Spectrum width

   character*40, parameter :: caller = "elsi_set_pexsi_delta_e"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%pexsi_options%deltaE = delta_e

end subroutine

!>
!! This routine sets the lower bound of the chemical potential in PEXSI.
!!
subroutine elsi_set_pexsi_mu_min(elsi_h,mu_min)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(in)    :: mu_min !< Lower bound of mu

   character*40, parameter :: caller = "elsi_set_pexsi_mu_min"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%pexsi_options%muMin0 = mu_min

end subroutine

!>
!! This routine sets the upper bound of the chemical potential in PEXSI.
!!
subroutine elsi_set_pexsi_mu_max(elsi_h,mu_max)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(in)    :: mu_max !< Upper bound of mu

   character*40, parameter :: caller = "elsi_set_pexsi_mu_max"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%pexsi_options%muMax0 = mu_max

end subroutine

!>
!! This routine sets the tolerance of the estimation of the chemical
!! potential in the inertia counting procedure.
!!
subroutine elsi_set_pexsi_inertia_tol(elsi_h,inertia_tol)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h      !< Handle
   real(kind=r8),     intent(in)    :: inertia_tol !< Tolerance of chemical potential

   character*40, parameter :: caller = "elsi_set_pexsi_inertia_tol"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%pexsi_options%muInertiaTolerance = inertia_tol

end subroutine

!>
!! This routine sets the initial guess of the error function decay
!! length in CheSS.
!!
subroutine elsi_set_chess_erf_decay(elsi_h,decay)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(in)    :: decay  !< Decay length

   character*40, parameter :: caller = "elsi_set_chess_erf_decay"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%erf_decay = decay

end subroutine

!>
!! This routine sets the lower bound of the decay length in CheSS.
!!
subroutine elsi_set_chess_erf_decay_min(elsi_h,decay_min)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h    !< Handle
   real(kind=r8),     intent(in)    :: decay_min !< Minimum decay length

   character*40, parameter :: caller = "elsi_set_chess_erf_decay_min"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%erf_decay_min = decay_min

end subroutine

!>
!! This routine sets the upper bound of the decay length in CheSS.
!!
subroutine elsi_set_chess_erf_decay_max(elsi_h,decay_max)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h    !< Handle
   real(kind=r8),     intent(in)    :: decay_max !< Maximum decay length

   character*40, parameter :: caller = "elsi_set_chess_erf_decay_max"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%erf_decay_max = decay_max

end subroutine

!>
!! This routine sets the lower bound of the eigenvalues of the
!! Hamiltonian matrix.
!!
subroutine elsi_set_chess_ev_ham_min(elsi_h,ev_min)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(in)    :: ev_min !< Minimum eigenvalue

   character*40, parameter :: caller = "elsi_set_chess_ev_ham_min"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%ev_ham_min = ev_min

end subroutine

!>
!! This routine sets the upper bound of the eigenvalues of the
!! Hamiltonian matrix.
!!
subroutine elsi_set_chess_ev_ham_max(elsi_h,ev_max)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(in)    :: ev_max !< Maximum eigenvalue

   character*40, parameter :: caller = "elsi_set_chess_ev_ham_max"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%ev_ham_max = ev_max

end subroutine

!>
!! This routine sets the lower bound of the eigenvalues of the
!! overlap matrix.
!!
subroutine elsi_set_chess_ev_ovlp_min(elsi_h,ev_min)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(in)    :: ev_min !< Minimum eigenvalue

   character*40, parameter :: caller = "elsi_set_chess_ev_ovlp_min"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%ev_ovlp_min = ev_min

end subroutine

!>
!! This routine sets the upper bound of the eigenvalues of the
!! overlap matrix.
!!
subroutine elsi_set_chess_ev_ovlp_max(elsi_h,ev_max)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(in)    :: ev_max !< Maximum eigenvalue

   character*40, parameter :: caller = "elsi_set_chess_ev_ovlp_max"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%ev_ovlp_max = ev_max

end subroutine

!>
!! This routine sets the slicing method when using SIPs.
!!
subroutine elsi_set_sips_slice_type(elsi_h,slice_type)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h     !< Handle
   integer(kind=i4),  intent(in)    :: slice_type !< Method of slicing

   character*40, parameter :: caller = "elsi_set_sips_slice_type"

   call elsi_check_handle(elsi_h,caller)

   if((slice_type < 0) .or. (slice_type > 5)) then
      call elsi_stop("  Invalid choice of slice_type. Please"//&
              " choose 0, 1, 2, 3, or 4. Exiting...",elsi_h,caller)
   endif

   elsi_h%slicing_method = slice_type

end subroutine

!>
!! This routine sets the number of slices in SIPs.
!!
subroutine elsi_set_sips_n_slice(elsi_h,n_slice)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h  !< Handle
   integer(kind=i4),  intent(in)    :: n_slice !< Number of slices

   character*40, parameter :: caller = "elsi_set_sips_n_slice"

   call elsi_check_handle(elsi_h,caller)

   if(mod(elsi_h%n_procs,n_slice) == 0) then
      elsi_h%n_slices = n_slice
      elsi_h%n_p_per_slice = elsi_h%n_procs/n_slice
   else
      call elsi_stop("  The total number of MPI tasks must be"//&
              " a multiple of the number of slices."//&
              " Exiting...",elsi_h,caller)
   endif

end subroutine

!>
!! This routine sets the method to bound the left side of the eigenvalue
!! interval in SIPs.
!!
subroutine elsi_set_sips_left_bound(elsi_h,left_bound)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h     !< Handle
   integer(kind=i4),  intent(in)    :: left_bound !< How to bound the left side?

   character*40, parameter :: caller = "elsi_set_sips_left_bound"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%unbound = left_bound

end subroutine

!>
!! This routine sets a small buffer to expand the eigenvalue interval
!! in SIPs.
!!
subroutine elsi_set_sips_slice_buf(elsi_h,slice_buffer)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h       !< Handle
   real(kind=r8),     intent(in)    :: slice_buffer !< Buffer to expand the interval

   character*40, parameter :: caller = "elsi_set_sips_slice_buf"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%slice_buffer = slice_buffer

end subroutine

!>
!! This routine sets the broadening scheme to determine the chemical
!! potential and the occupation numbers.
!!
subroutine elsi_set_mu_broaden_scheme(elsi_h,broaden_scheme)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h         !< Handle
   integer(kind=i4),  intent(in)    :: broaden_scheme !< Broadening method

   character*40, parameter :: caller = "elsi_set_mu_broaden_method"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%broaden_scheme = broaden_scheme

end subroutine

!>
!! This routine sets the broadening width to determine the chemical
!! potential and the occupation numbers.
!!
subroutine elsi_set_mu_broaden_width(elsi_h,broaden_width)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h        !< Handle
   real(kind=r8),     intent(in)    :: broaden_width !< Broadening width

   character*40, parameter :: caller = "elsi_set_mu_broaden_width"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%broaden_width = broaden_width

end subroutine

!>
!! This routine sets the desired accuracy of the determination of the
!! chemical potential and the occupation numbers.
!!
subroutine elsi_set_mu_tol(elsi_h,mu_tol)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(in)    :: mu_tol !< Accuracy of mu

   character*40, parameter :: caller = "elsi_set_mu_tol"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%occ_tolerance = mu_tol

end subroutine

!>
!! This routine sets the spin degeneracy in the determination of the
!! chemical potential and the occupation numbers.
!!
subroutine elsi_set_mu_spin_degen(elsi_h,spin_degen)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h     !< Handle
   real(kind=r8),     intent(in)    :: spin_degen !< Spin degeneracy

   character*40, parameter :: caller = "elsi_set_mu_spin_degen"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%spin_degen = spin_degen
   elsi_h%spin_is_set = .true.

end subroutine

!>
!! This routine gets the lower bound of the chemical potential
!! returned by the inertia counting in PEXSI.
!!
subroutine elsi_get_pexsi_mu_min(elsi_h,mu_min)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(out)   :: mu_min !< Lower bound of mu

   character*40, parameter :: caller = "elsi_get_pexsi_mu_min"

   call elsi_check_handle(elsi_h,caller)

   mu_min = elsi_h%pexsi_options%muMin0

end subroutine

!>
!! This routine gets the upper bound of the chemical potential
!! returned by the inertia counting in PEXSI.
!!
subroutine elsi_get_pexsi_mu_max(elsi_h,mu_max)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(out)   :: mu_max !< Upper bound of mu

   character*40, parameter :: caller = "elsi_get_pexsi_mu_max"

   call elsi_check_handle(elsi_h,caller)

   mu_max = elsi_h%pexsi_options%muMax0

end subroutine

!>
!! This routine gets the result of the singularity check of the
!! overlap matrix.
!!
subroutine elsi_get_ovlp_sing(elsi_h,ovlp_sing)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h    !< Handle
   integer(kind=i4),  intent(out)   :: ovlp_sing !< Is overlap singular?

   character*40, parameter :: caller = "elsi_get_ovlp_sing"

   call elsi_check_handle(elsi_h,caller)

   if(elsi_h%ovlp_is_sing) then
      ovlp_sing = 1
   else
      ovlp_sing = 0
   endif

end subroutine

!>
!! This routine gets the number of basis functions that are removed
!! due to overlap singularity.
!!
subroutine elsi_get_n_sing(elsi_h,n_sing)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   integer(kind=i4),  intent(out)   :: n_sing !< Number of singular basis

   character*40, parameter :: caller = "elsi_get_n_sing"

   call elsi_check_handle(elsi_h,caller)

   n_sing = elsi_h%n_basis-elsi_h%n_nonsing

end subroutine

!>
!! This routine gets the chemical potential.
!!
subroutine elsi_get_mu(elsi_h,mu)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   real(kind=r8),     intent(out)   :: mu     !< Chemical potential

   character*40, parameter :: caller = "elsi_get_mu"

   call elsi_check_handle(elsi_h,caller)

   mu = elsi_h%mu

   if(.not. elsi_h%mu_ready) then
      call elsi_statement_print("  ATTENTION! The return value of mu may"//&
              " be 0, since it has not been computed.",elsi_h)
   endif

   elsi_h%mu_ready = .false.

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_real(elsi_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                 !< Handle
   real(kind=r8),     intent(out)   :: d_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Energy density matrix

   character*40, parameter :: caller = "elsi_get_edm_real"

   call elsi_check_handle(elsi_h,caller)

   if(elsi_h%edm_ready_real) then
      elsi_h%matrix_data_type = REAL_VALUES

      select case(elsi_h%solver)
      case(ELPA)
         call elsi_set_dm(elsi_h,d_out)

         call elsi_compute_edm_elpa(elsi_h)
      case(LIBOMM)
         call elsi_set_dm(elsi_h,d_out)

         call elsi_compute_edm_omm(elsi_h)

         elsi_h%dm_omm%dval = 2.0_r8*elsi_h%dm_omm%dval
      case(PEXSI)
         call elsi_compute_edm_pexsi(elsi_h)

         call elsi_pexsi_to_blacs_dm(elsi_h,d_out)
      case(CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case(SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case default
         call elsi_stop(" No supported solver has been chosen."//&
                 " Exiting...",elsi_h,caller)
      end select

      elsi_h%edm_ready_real = .false.
      elsi_h%matrix_data_type = UNSET
   else
      call elsi_stop(" Energy weighted density matrix has not been."//&
              " computed. Exiting...",elsi_h,caller)
   endif

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_real_sparse(elsi_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                 !< Handle
   real(kind=r8),     intent(out)   :: d_out(elsi_h%nnz_l_sp) !< Energy density matrix

   character*40, parameter :: caller = "elsi_get_edm_real_sparse"

   call elsi_check_handle(elsi_h,caller)

   if(elsi_h%edm_ready_real) then
      elsi_h%matrix_data_type = REAL_VALUES

      select case(elsi_h%solver)
      case(ELPA)
         call elsi_compute_edm_elpa(elsi_h)

         call elsi_blacs_to_sips_dm(elsi_h,d_out)
      case(LIBOMM)
         call elsi_compute_edm_omm(elsi_h)

         elsi_h%dm_omm%dval = 2.0_r8*elsi_h%dm_omm%dval

         call elsi_blacs_to_sips_dm(elsi_h,d_out)
      case(PEXSI)
         call elsi_set_sparse_dm(elsi_h,d_out)

         call elsi_compute_edm_pexsi(elsi_h)
      case(CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case(SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case default
         call elsi_stop(" No supported solver has been chosen."//&
                 " Exiting...",elsi_h,caller)
      end select

      elsi_h%edm_ready_real = .false.
      elsi_h%matrix_data_type = UNSET
   else
      call elsi_stop(" Energy weighted density matrix has not been."//&
              " computed. Exiting...",elsi_h,caller)
   endif

end subroutine

!>
!! This routine gets the energy-weighted density matrix.
!!
subroutine elsi_get_edm_complex(elsi_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                 !< Handle
   complex(kind=r8),  intent(out)   :: d_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Energy density matrix

   character*40, parameter :: caller = "elsi_get_edm_complex"

   call elsi_check_handle(elsi_h,caller)

   if(elsi_h%edm_ready_complex) then
      elsi_h%matrix_data_type = COMPLEX_VALUES

      select case(elsi_h%solver)
      case(ELPA)
         call elsi_set_dm(elsi_h,d_out)

         call elsi_compute_edm_elpa(elsi_h)
      case(LIBOMM)
         call elsi_set_dm(elsi_h,d_out)

         call elsi_compute_edm_omm(elsi_h)

         elsi_h%dm_omm%zval = (2.0_r8,0.0_r8)*elsi_h%dm_omm%zval
      case(PEXSI)
         call elsi_compute_edm_pexsi(elsi_h)

         call elsi_pexsi_to_blacs_dm(elsi_h,d_out)
      case(CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case(SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case default
         call elsi_stop(" No supported solver has been chosen."//&
                 " Exiting...",elsi_h,caller)
      end select

      elsi_h%edm_ready_complex = .false.
      elsi_h%matrix_data_type = UNSET
   else
      call elsi_stop(" Energy weighted density matrix has not been."//&
              " computed. Exiting...",elsi_h,caller)
   endif

end subroutine

end module ELSI_MUTATOR
