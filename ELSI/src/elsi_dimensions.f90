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
!! This module contains variables accessible in ELSI and related modules.
!!

module ELSI_DIMENSIONS

  use iso_c_binding
  use f_ppexsi_interface

  implicit none
  
  !> Solver (AUTO=0,ELPA=1,LIBOMM=2,PEXSI=3,CHESS=4)
  integer :: method = -1

  !> Real or complex data (REAL_VALUES=0,COMPLEX_VALUES=1)
  integer :: mode = -1

  !> Matrix storage format (BLOCK_CYCLIC=0)
  integer :: storage = -1

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
  ! PEXSI may use another BLACS setup
  integer :: my_p_row_pexsi
  integer :: my_p_col_pexsi
  integer :: n_b_rows_pexsi
  integer :: n_b_cols_pexsi
  integer :: n_p_rows_pexsi
  integer :: n_p_cols_pexsi
  integer :: n_l_rows_pexsi
  integer :: n_l_cols_pexsi

  !> Sparse matrix information
  integer :: nnz_g !< Global number of nonzeros
  integer :: nnz_l !< Local number of nonzeros
  integer :: nnz_l_pexsi !< Local number of nonzeros in PEXSI distribution
  real*8, parameter :: threshold = 1.0d-15 !< Threshold to define numerical zero

  !> Overlap
  logical :: overlap_is_unity = .false. !< Is overlap matrix unity?

  !> Physics
  real*8  :: n_electrons !< Number of electrons in system
  integer :: n_states !< Number of total states
  integer :: n_occupied_states !< Number of occupied states
  real*8, parameter :: hartree = 27.211386 !< Hartree to eV (source: Codata 2015)

  !> ELPA
  logical :: elpa_one_always = .false. !< Always use 1-stage solver
  logical :: elpa_two_always = .false. !< Always use 2-stage solver
  integer :: broadening_type = -1 !< Broadening scheme used to find occupation numbers
  real*8  :: broadening_width !< Broadening width used to find occupation numbers

  !> OMM
  logical :: omm_customized = .false. !< Has elsi_customize_omm been called?
  logical :: new_overlap !< Is a new overlap matrix provided?
  logical :: C_matrix_initialized !< Is coefficient matrix initialized?
  real*8  :: total_energy !< Energy of the system
  !> How to perform the calculation
  !! 0 = Basic
  !! 1 = Cholesky factorisation of S requested
  !! 2 = Cholesky already performed, U is provided in S
  !! 3 = Use preconditioning based on the energy density
  integer :: omm_flavour = -1
  real*8  :: scale_kinetic       !< Scaling of the kinetic energy matrix
  logical :: calc_ed = .False.   !< Calculate energy weighted density matrix?
  real*8  :: eta                 !< Eigenspectrum shift parameter
  real*8  :: min_tol             !< Tolerance for minimization
  integer :: nk_times_nspin = -1 !< n_k_points * n_spin
  integer :: i_k_spin = -1       !< Combined k_point spin index
  logical :: omm_verbose         !< Output level
  logical :: do_dealloc          !< Deallocate internal storage?

  !> PEXSI
  logical                :: pexsi_customized = .false.
  integer(c_intptr_t)    :: pexsi_plan
  type(f_ppexsi_options) :: pexsi_options
  integer(c_int)         :: pexsi_info
  integer(c_int)         :: pexsi_output_file_index
  real(c_double)         :: mu_pexsi
  real(c_double)         :: n_electrons_pexsi
  real(c_double)         :: mu_min_inertia
  real(c_double)         :: mu_max_inertia
  integer(c_int)         :: n_total_inertia_iter
  integer(c_int)         :: n_total_pexsi_iter
  real(c_double)         :: e_tot_h
  real(c_double)         :: e_tot_s
  real(c_double)         :: f_tot

  !> Method names
  enum, bind( C )
     enumerator :: AUTO, ELPA, LIBOMM, PEXSI, CHESS
  end enum

  !> Mode
  enum, bind( C )
     enumerator :: REAL_VALUES, COMPLEX_VALUES
  end enum

  !> Storage formats
  enum, bind( C )
     enumerator :: BLOCK_CYCLIC
  end enum

  !> Broadening type (used if ELPA is chosen to compute density matrix)
  enum, bind( C )
     enumerator :: GAUSSIAN, FERMI
  end enum

contains

!>
!! Set PEXSI variables to ELSI default.
!!
subroutine elsi_set_pexsi_default_options()

   implicit none
   include "mpif.h"

   ! Use the PEXSI Default options
   call f_ppexsi_set_default_options(pexsi_options)

end subroutine

!>
!! Print PEXSI settings.
!!
subroutine elsi_print_pexsi_options()

   implicit none
   include "mpif.h"

   character(LEN=4096) :: string_message

   if(myid == 0) then
      write(*,"(A)") "  PEXSI settings used in ELSI (in the same unit of Hamiltonian):"

      write(string_message, "(1X,' | Temperature ',F10.4)") &
            pexsi_options%temperature
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Spectral Gap ',F10.4)") &
            pexsi_options%gap
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Delta E ',F10.4)") &
            pexsi_options%deltaE
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Number of Poles ',I5)") &
            pexsi_options%numPole
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Use Inertia Count ',I2)") &
            pexsi_options%isInertiaCount
      write(*,'(A)') trim(string_message)
      
      write(string_message, "(1X,' | Max Pexsi Iterations ',I5)") &
            pexsi_options%maxPEXSIIter
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Minimal Mu ',F10.4)") &
            pexsi_options%muMin0
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Maximal Mu ',F10.4)") &
            pexsi_options%muMax0
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Start Mu ',F10.4)") &
            pexsi_options%mu0
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Tolerance Mu ',E10.1)") &
            pexsi_options%muInertiaTolerance
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Expansion Mu ',F10.4)") &
            pexsi_options%muInertiaExpansion
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Safeguard Mu ',F10.4)") &
            pexsi_options%muPexsiSafeGuard
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Tolerance Electrons ',E10.1)") &
            pexsi_options%numElectronPEXSITolerance
      write(*,'(A)') trim(string_message)
   endif

end subroutine

!>
!! Set OMM variables to ELSI default.
!!
subroutine elsi_set_omm_default_options()

   implicit none
   include "mpif.h"

   !< How do we perform the calculation
   !! 0 = Basic
   !! 1 = Cholesky factorisation of S requested
   !! 2 = Cholesky already performed, U is provided in S
   !! 3 = Use preconditioning based on the energy density
   omm_flavour = 2
   !< Scaling of the kinetic energy matrix
   scale_kinetic = 5d0
   !< Calculate the energy weighted density matrix
   calc_ed = .false.
   !< Eigenspectrum shift parameter
   eta = 0d0
   !< Tolerance for minimization
   min_tol = 1.0d-12
   !< n_k_points * n_spin
   nk_times_nspin = 1
   !< Combined k_point spin index
   i_k_spin = 1
   omm_verbose = .False.
   do_dealloc = .False.

end subroutine

!>
!! Print OMM settings.
!!
subroutine elsi_print_omm_options()

   implicit none
   include "mpif.h"

   character(LEN=4096) :: string_message

   if(myid == 0) then
      write(*,"(A)") "  OMM settings used in ELSI: "

      write(string_message, "(1X,' | Eta (H) ',F10.4)") eta
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Kinetic scaling (H) ',F10.4)") scale_kinetic
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Tolerance (H) ',E10.1)") min_tol
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | OMM Flavour ',I1)") omm_flavour
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Energy weighted densigy matrix? ',L1)") calc_ed
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Verbose? ',L1)") omm_verbose
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | Dealloc?  ',L1)") do_dealloc
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | k points x spin ',I1)") nk_times_nspin
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,' | index ',I1)") i_k_spin
      write(*,'(A)') trim(string_message)
   endif

end subroutine

end module ELSI_DIMENSIONS
