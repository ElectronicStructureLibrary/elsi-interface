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
   use ELPA1
   use ELPA2
   use MatrixSwitch
   use f_ppexsi_interface

   implicit none
   private

   ! ELSI
   integer, save :: n_elsi_calls = 0

   ! Pointers
   !< Real Hamiltonian
   real*8, pointer     :: H_real(:,:)
   !< Complex Hamiltonian
   complex*16, pointer :: H_complex(:,:)
   !< Real overlap
   real*8, pointer     :: S_real(:,:)
   !< Complex overlap
   complex*16, pointer :: S_complex(:,:)

   ! ELPA
   !< Eigenvalues
   real*8, allocatable     :: eigenvalues(:)
   !< Real eigenvectors
   real*8, allocatable     :: C_real(:,:)
   !< Complex eigenvectors
   complex*16, allocatable :: C_complex(:,:)
   !< Density matrix
   real*8, allocatable     :: D_elpa(:,:)
   !< Occupation numbers used by ELPA to construct density matrix
   real*8, allocatable     :: occ_elpa(:)

   ! OMM
   !< Hamiltonian
   type(Matrix) :: H_omm
   !< Overlap
   type(Matrix) :: S_omm
   !< Coefficient matrix
   type(Matrix) :: Coeff_omm
   !< Density matrix
   type(Matrix) :: D_omm
   !< Kinetic energy density matrix
   type(Matrix) :: T_omm

   ! PESXI
   !< Sparse real Hamiltonian
   real*8, allocatable :: H_real_pexsi(:)
   !< Sparse complex Hamiltonian
   complex*16, allocatable :: H_complex_pexsi(:)
   !< Sparse real overlap
   real*8, allocatable :: S_real_pexsi(:)
   !< Sparse complex overlap
   complex*16, allocatable :: S_complex_pexsi(:)
   !< Sparse density matrix
   real*8, allocatable :: D_pexsi(:)
   !< Sparse energy density matrix
   real*8, allocatable :: ED_pexsi(:)
   !< Sparse free energy density matrix
   real*8, allocatable :: FD_pexsi(:)
   !< Row index in CCS format
   integer, allocatable :: row_ind_pexsi(:)
   !< Column pointer in CCS format
   integer, allocatable :: col_ptr_pexsi(:)

   !< From ELSI_DIMENSIONS
   public :: AUTO, ELPA, LIBOMM, PEXSI, CHESS   !< Solvers
   public :: REAL_VALUES, COMPLEX_VALUES        !< Real or complex
   public :: DENSE, CCS, CSC, CRS, CSR          !< Storage formats
   public :: GAUSSIAN, FERMI, METHFESSEL_PAXTON !< Broadening schemes

   ! Public routines
   public :: elsi_init            !< Initialize
   public :: elsi_set_method      !< Select solver
   public :: elsi_set_storage     !< Select matrix storage format
   public :: elsi_set_mpi         !< Set MPI from calling code
   public :: elsi_set_blacs       !< Set BLACS from calling code
   public :: elsi_customize       !< Override ELSI default
   public :: elsi_customize_elpa  !< Override ELPA default
   public :: elsi_customize_omm   !< Override OMM default
   public :: elsi_customize_pexsi !< Override PEXSI default
   public :: elsi_ev_real         !< Compute eigenvalues and eigenvectors
   public :: elsi_ev_complex      !< Compute eigenvalues and eigenvectors
   public :: elsi_dm_real         !< Compute density matrix
   public :: elsi_dm_complex      !< Compute density matrix
   public :: elsi_finalize        !< Clean memory and print timings

   integer, external :: numroc

   interface elsi_set_hamiltonian
      module procedure elsi_set_real_hamiltonian,&
                       elsi_set_complex_hamiltonian
   end interface

   interface elsi_set_overlap
      module procedure elsi_set_real_overlap,&
                       elsi_set_complex_overlap
   end interface

   interface elsi_get_eigenvectors
      module procedure elsi_get_real_eigenvectors,&
                       elsi_get_complex_eigenvectors
   end interface

contains

!===============
! ELSI routines
!===============

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
        ! Set number of occupied states for OMM
        n_states = nint(n_electrons/2d0)
     else
        n_states = n_states_in
     endif

     call elsi_init_timers()

  end subroutine ! elsi_init

!>
!! This routine sets the method.
!!
  subroutine elsi_set_method(i_method)

     implicit none

     integer, intent(in) :: i_method !< AUTO,ELPA,OMM,PEXSI,CHESS,...
   
     method = i_method

  end subroutine ! elsi_set_method

!>
!! This routine sets the mode (real or complex).
!!
  subroutine elsi_set_mode(i_mode)

     implicit none

     integer, intent(in) :: i_mode !< REAL_VALUES,COMPLEX_VALUES

     mode = i_mode

  end subroutine ! elsi_set_mode

!>
!! This routine sets the matrix storage format.
!!
  subroutine elsi_set_storage(i_storage)

     implicit none

     integer, intent(in) :: i_storage !< DENSE,CCS,CSC,CRS,CSR,...

     storage = i_storage

  end subroutine ! elsi_set_storage

!>
!! This routine sets the parallel mode.
!!
  subroutine elsi_set_parallel(i_parallel)

     implicit none

     integer, intent(in) :: i_parallel !< SERIAL,PARALLEL

     parallelism = i_parallel

     if(i_parallel == 0) then ! Serial version
        n_l_rows = n_g_size
        n_l_cols = n_g_size
        n_b_rows = n_g_size
        n_b_cols = n_g_size
     endif

  end subroutine ! elsi_set_parallel

!>
!! This routine overrides ELSI default settings.
!!
  subroutine elsi_customize(overlap_is_unit_in,hartree_in,zero_threshold_in,&
                            mu_accuracy_in,singularity_tolerance_in)

     implicit none

     logical, intent(in), optional :: overlap_is_unit_in
     real*8,  intent(in), optional :: hartree_in
     real*8,  intent(in), optional :: zero_threshold_in
     real*8,  intent(in), optional :: mu_accuracy_in
     real*8,  intent(in), optional :: singularity_tolerance_in

     ! Is the overlap matrix unit? [Default: .false.]
     if(present(overlap_is_unit_in)) &
        overlap_is_unit = overlap_is_unit_in
     ! User-defined value for Hartree [Default: 27.211386, Codata 2015]
     if(present(hartree_in)) &
        hartree = hartree_in
     ! Threshold to define numerical zero [Default: 1d-13]
     if(present(zero_threshold_in)) &
        zero_threshold = zero_threshold_in
     ! Accuracy for chemical potential determination [Default: 1d-10]
     if(present(mu_accuracy_in)) &
        occ_tolerance = mu_accuracy_in
     ! Eigenfunctions of overlap matrix with eigenvalues smaller than
     ! this value will be removed to avoid singularity [Default: 1d-5]
     if(present(singularity_tolerance_in)) &
        singularity_tolerance = singularity_tolerance_in

  end subroutine ! elsi_customize

!>
!! Set MPI.
!!
  subroutine elsi_set_mpi(mpi_comm_global_in,n_procs_in,myid_in)

     implicit none
     include "mpif.h"

     integer, intent(in) :: myid_in            !< local process id
     integer, intent(in) :: n_procs_in         !< number of mpi processes
     integer, intent(in) :: mpi_comm_global_in !< global mpi communicator

     mpi_is_setup    = .true.
     mpi_comm_global = mpi_comm_global_in
     n_procs         = n_procs_in
     myid            = myid_in

  end subroutine

!>
!! Set BLACS.
!!
  subroutine elsi_set_blacs(blacs_ctxt_in,n_b_rows_in,n_b_cols_in,n_p_rows_in,&
                            n_p_cols_in,mpi_comm_row_in,mpi_comm_col_in)

     implicit none

     integer, intent(in) :: blacs_ctxt_in             !< BLACS context
     integer, intent(in) :: n_b_rows_in               !< Block size
     integer, intent(in) :: n_b_cols_in               !< Block size
     integer, intent(in) :: n_p_rows_in               !< Number of processes in row
     integer, intent(in) :: n_p_cols_in               !< Number of processes in column
     integer, intent(in), optional :: mpi_comm_row_in !< row communicatior for ELPA
     integer, intent(in), optional :: mpi_comm_col_in !< column communicatior for ELPA

     integer :: blacs_info
     character*40, parameter :: caller = "elsi_set_blacs"

     blacs_is_setup = .true.

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

     mpi_comm_row = mpi_comm_row_in
     mpi_comm_col = mpi_comm_col_in

     if(method == LIBOMM) then
        call ms_scalapack_setup(myid,n_procs,n_p_rows,'r',n_b_rows,icontxt=blacs_ctxt)
     endif
  
  end subroutine

!>
!! This routine prepares the matrices.
!!
  subroutine elsi_allocate_matrices()

     implicit none

     character*40, parameter :: caller = "elsi_allocate_matrices"

     select case (method)
        case (ELPA)
           if(.not.allocated(eigenvalues)) then
              call elsi_allocate(eigenvalues, n_g_size, "eigenvalues", caller)
           endif
           eigenvalues = 0d0
           select case (mode)
              case (COMPLEX_VALUES)
                 if(.not.allocated(C_complex)) then
                    call elsi_allocate(C_complex, n_l_rows, n_l_cols, "C_complex", caller)
                 endif
                 C_complex = CMPLX(0d0, 0d0)
              case (REAL_VALUES)
                 if(.not.allocated(C_real)) then
                    call elsi_allocate(C_real, n_l_rows, n_l_cols, "C_real", caller)
                 endif
                 C_real = 0d0
              case DEFAULT
                 call elsi_stop(" No mode has been chosen. "//&
                                " Please choose REAL_VALUES or COMPLEX_VALUES. ",&
                                caller)
           end select

        case (LIBOMM)
           select case (mode)
              case (COMPLEX_VALUES)
                 if(.not.D_omm%is_initialized) then
                    call m_allocate(D_omm, n_g_size, n_g_size, "pddbc")
                 endif
                 if(.not.Coeff_omm%is_initialized) then
                    call m_allocate(Coeff_omm, n_states, n_g_size, "pddbc")
                 endif
              case (REAL_VALUES)
                 if(.not.D_omm%is_initialized) then
                    call m_allocate(D_omm, n_g_size, n_g_size, "pddbc")
                 endif
                 if(.not.Coeff_omm%is_initialized) then
                    call m_allocate(Coeff_omm, n_states, n_g_size, "pddbc")
                 endif
              case DEFAULT
                 call elsi_stop(" No mode has been chosen. "//&
                                " Please choose REAL_VALUES or COMPLEX_VALUES. ",&
                                caller)
           end select

        case (PEXSI)
           ! Nothing to be done here
        case (CHESS)
           call elsi_stop(" CHESS not yet implemented. Exiting...", caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...", caller)
     end select
  end subroutine ! elsi_allocate_matrices

!>
!! This routine sets the real hamiltonian matrix.
!!
  subroutine elsi_set_real_hamiltonian(H_in)

     implicit none

     real*8, target, intent(in) :: H_in(n_l_rows, n_l_cols) !< Hamiltonian

     character*40, parameter :: caller = "elsi_set_real_hamiltonian"

     select case (method)
        case (ELPA)
           H_real => H_in
        case (LIBOMM)
           call m_register_pdbc(H_omm, H_in, sc_desc)
        case (PEXSI)
           call elsi_stop(" PEXSI not yet implemented. Exiting...", caller)
        case (CHESS)
           call elsi_stop(" CHESS not yet implemented. Exiting...", caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...", caller)
     end select

  end subroutine ! elsi_set_real_hamiltonian

!>
!! This routine sets the complex hamiltonian matrix.
!!
  subroutine elsi_set_complex_hamiltonian(H_in)

     implicit none

     complex*16, target, intent(in) :: H_in(n_l_rows, n_l_cols) !< Hamiltonian

     character*40, parameter :: caller = "elsi_set_complex_hamiltonian"

     select case (method)
        case (ELPA)
           H_complex => H_in
        case (LIBOMM)
           call m_register_pdbc(H_omm, H_in, sc_desc)
        case (PEXSI)
           call elsi_stop(" PEXSI not yet implemented. Exiting...", caller)
        case (CHESS)
           call elsi_stop(" CHESS not yet implemented. Exiting...", caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...", caller)
     end select

  end subroutine ! elsi_set_complex_hamiltonian

!>
!! This routine sets the real overlap matrix.
!!
  subroutine elsi_set_real_overlap(S_in)

     implicit none

     real*8, target, intent(in) :: S_in(n_l_rows, n_l_cols) !< Overlap

     character*40, parameter :: caller = "elsi_set_real_overlap"

     select case (method)
        case (ELPA)
           S_real => S_in
        case (LIBOMM)
           call m_register_pdbc(S_omm, S_in, sc_desc)
        case (PEXSI)
           call elsi_stop(" PEXSI not yet implemented. Exiting...", caller)
        case (CHESS)
           call elsi_stop(" CHESS not yet implemented. Exiting...", caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...", caller)
     end select

  end subroutine ! elsi_set_real_overlap

!>
!! This routine sets the complex overlap matrix.
!!
  subroutine elsi_set_complex_overlap(S_in)

     implicit none

     complex*16, target, intent(in) :: S_in(n_l_rows, n_l_cols) !< Overlap

     character*40, parameter :: caller = "elsi_set_complex_overlap"

     select case (method)
        case (ELPA)
           S_complex => S_in
        case (LIBOMM)
           call m_register_pdbc(S_omm, S_in, sc_desc)
        case (PEXSI)
           call elsi_stop(" PEXSI not yet implemented. Exiting...", caller)
        case (CHESS)
           call elsi_stop(" CHESS not yet implemented. Exiting...", caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...", caller)
     end select

  end subroutine ! elsi_set_complex_overlap

!>
!! This routine gets the energy.
!!
  subroutine elsi_get_energy(energy_out)

     implicit none

     real*8, intent(out) :: energy_out !< Energy of the system

     ! Only spin-nonpolarized case is supported now.
     real*8, parameter :: n_spin = 2d0

     integer :: i_state
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
           call elsi_stop(" CHESS: not yet implemented! Exiting...", caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...", caller)
     end select

     if(myid == 0) then
        write(*,"(A,F15.5,A)") "  | Energy = ", energy_out, " Ha"
        write(*,"(A,F15.5,A)") "  |        = ", energy_out*hartree, " eV"
     endif

  end subroutine ! elsi_get_energy

!>
!! This routine gets the density matrix.
!!
  subroutine elsi_get_dm(D_out)

     implicit none

     real*8, intent(out) :: D_out(n_l_rows, n_l_cols) !< Density matrix

     character*40, parameter :: caller = "elsi_get_dm"

     select case (method)
        case (ELPA)
           call elsi_stop(" ELPA needs to compute density matrix from eigenvectors. "//&
                          " Exiting...", caller)
        case (LIBOMM)
           ! Note that here the convention of density matrix used in OMM is changed
           ! to the one used in ELPA and PEXSI.
           D_out = 2d0*D_omm%dval
        case (PEXSI)
           call elsi_1dbccs_to_2dbcd_dm_pexsi(D_out)
        case (CHESS)
           call elsi_stop(" CHESS: not yet implemented! Exiting...", caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...", caller)
     end select

  end subroutine ! elsi_get_dm

!>
!! This routine deallocates the matrices.
!!
  subroutine elsi_deallocate_matrices()

     implicit none

     ! Nullify pointers
     if(associated(H_real))       nullify(H_real)
     if(associated(H_complex))    nullify(H_complex)
     if(associated(S_real))       nullify(S_real)
     if(associated(S_complex))    nullify(S_complex)
 
     ! Free Memory
     ! ELPA
     if(allocated(C_real))        deallocate(C_real)
     if(allocated(C_complex))     deallocate(C_complex)
     if(allocated(eigenvalues))   deallocate(eigenvalues)
     if(allocated(D_elpa))        deallocate(D_elpa)
     if(allocated(occ_elpa))      deallocate(occ_elpa)
     ! PEXSI
     if(allocated(H_real_pexsi))  deallocate(H_real_pexsi)
     if(allocated(S_real_pexsi))  deallocate(S_real_pexsi)
     if(allocated(D_pexsi))       deallocate(D_pexsi)
     if(allocated(ED_pexsi))      deallocate(ED_pexsi)
     if(allocated(FD_pexsi))      deallocate(FD_pexsi)
     if(allocated(row_ind_pexsi)) deallocate(row_ind_pexsi)
     if(allocated(col_ptr_pexsi)) deallocate(col_ptr_pexsi)
     ! OMM
     if(H_omm%is_initialized)     call m_deallocate(H_omm)
     if(S_omm%is_initialized)     call m_deallocate(S_omm)
     if(D_omm%is_initialized)     call m_deallocate(D_omm)
     if(Coeff_omm%is_initialized) call m_deallocate(Coeff_omm)
     if(T_omm%is_initialized)     call m_deallocate(T_omm)

  end subroutine ! elsi_deallocate_matrices

!>
!! This routine finalizes ELSI.
!!
  subroutine elsi_finalize()

     implicit none
     include "mpif.h"

     call MPI_Barrier(mpi_comm_global, mpierr)

     call elsi_deallocate_matrices()

     n_elsi_calls = 0

     call elsi_print_timers()

  end subroutine ! elsi_finalize

!========================
! ELSI routines for ELPA
!========================

!>
!! This routine gets the eigenvalues.
!!
  subroutine elsi_get_eigenvalues(e_val_out)

     implicit none

     real*8, intent(out) :: e_val_out(n_states) !< Eigenvalues

     character*40, parameter :: caller = "elsi_get_eigenvalues"

     select case (method)
        case (ELPA)
           e_val_out(1:n_states) = eigenvalues(1:n_states)
        case (LIBOMM)
           call elsi_stop(" OMM does not compute eigenvalues! Exiting...", caller)
        case (PEXSI)
           call elsi_stop(" PEXSI does not compute eigenvalues! Exiting...", caller)
        case (CHESS)
           call elsi_stop(" CHESS does not compute eigenvalues! Exiting...", caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...", caller)
     end select

  end subroutine ! elsi_get_eigenvalues

!>
!! This routine gets the eigenvectors.
!!
  subroutine elsi_get_real_eigenvectors(e_vec_out)

     implicit none

     real*8, intent(out) :: e_vec_out(n_l_rows, n_l_cols) !< Eigenvectors

     character*40, parameter :: caller = "elsi_get_real_eigenvectors"

     select case (method)
        case (ELPA)
           e_vec_out = C_real
        case (LIBOMM)
           call elsi_stop(" OMM does not compute eigenvectors! Exiting...", caller)
        case (PEXSI)
           call elsi_stop(" PEXSI does not compute eigenvectors! Exiting...", caller)
        case (CHESS)
           call elsi_stop(" CHESS does not compute eigenvectors! Exiting...", caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...", caller)
     end select

  end subroutine ! elsi_get_real_eigenvectors

!>
!! This routine gets the eigenvectors.
!!
  subroutine elsi_get_complex_eigenvectors(e_vec_out)

     implicit none

     complex*16, intent(out) :: e_vec_out(n_l_rows, n_l_cols) !< Eigenvectors

     character*40, parameter :: caller = "elsi_get_complex_eigenvectors"

     select case (method)
        case (ELPA)
           e_vec_out = C_complex
        case (LIBOMM)
           call elsi_stop(" OMM does not compute eigenvectors! Exiting...", caller)
        case (PEXSI)
           call elsi_stop(" PEXSI does not compute eigenvectors! Exiting...", caller)
        case (CHESS)
           call elsi_stop(" CHESS does not compute eigenvectors! Exiting...", caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...", caller)
     end select

  end subroutine ! elsi_get_complex_eigenvectors

!>
!! This routine computes the chemical potential and occupation numbers
!! from eigenvalues and eigenvectors.
!!
  subroutine elsi_compute_occ_elpa()

     implicit none

     real*8 :: mu !< Chemical potential
     real*8 :: e_low !< Lowest eigenvalue
     real*8 :: e_high !< Highest eigenvalue
     real*8 :: mu_lower !< Lower bound of chemical potential
     real*8 :: mu_upper !< Upper bound of chemical potential
     real*8 :: diff_ne_lower !< Difference in number of electrons on lower bound
     real*8 :: diff_ne_upper !< Difference in number of electrons on upper bound

     integer :: i_state !< State index
     integer :: n_steps !< Number of steps to find chemical potential interval
     integer, parameter :: max_steps = 100 !< Maximum number of steps

     character*40, parameter :: caller = "elsi_compute_occ_elpa"

     ! Determine the smallest and largest eivenvalues
     e_low = eigenvalues(1)
     e_high = eigenvalues(n_states)

     do i_state = 1,n_states
        if(eigenvalues(i_state) < e_low) e_low = eigenvalues(i_state)
        if(eigenvalues(i_state) > e_high) e_high = eigenvalues(i_state)
     enddo

     ! Determine the upper and lower bounds for chemical potential
     mu_lower = e_low

     if(e_low == e_high) then
        mu_upper = 0.0
     else
        mu_upper = e_high
     endif

     if(.not.allocated(occ_elpa)) then
         call elsi_allocate(occ_elpa, n_states, "occ_elpa", caller)
     endif
     occ_elpa = 0d0

     ! Compute the difference of number of electrons
     call elsi_get_ne(mu_lower, diff_ne_lower)
     call elsi_get_ne(mu_upper, diff_ne_upper)

     ! If diff_ne_lower*diff_ne_upper > 0, it means that the solution is
     ! not in this interval.
     ! Enlarge the interval towards both sides, then recheck the condition.
     n_steps = 0
     do while(diff_ne_lower*diff_ne_upper > 0)
        n_steps = n_steps+1
        if(n_steps > max_steps) then
           call elsi_stop(" Chemical potential not found in 100 iterations! "//&
                          " Exiting...", caller)
        endif

        mu_lower = mu_lower-0.5d0*abs(e_high-e_low)
        mu_upper = mu_upper+0.5d0*abs(e_high-e_low)

        call elsi_get_ne(mu_lower, diff_ne_lower)
        call elsi_get_ne(mu_upper, diff_ne_upper)
     enddo

     ! At this point we should have the correct interval for chemical potential.
     ! Use simple bisection algorithm to find the solution.
     call elsi_get_mu(mu_lower, mu_upper, mu)

  end subroutine ! elsi_compute_occ_elpa

!>
!! This routine computes the number of electrons using a given chemical potential,
!! and returns the difference in number of electrons. The occupation numbers will
!! be updated as well.
!!
  subroutine elsi_get_ne(mu_in,diff_ne_out)

     implicit none

     real*8,  intent(in)  :: mu_in       !< Input chemical potential
     real*8,  intent(out) :: diff_ne_out !< Difference in number of electrons

     real*8 :: invert_width !< 1/broadening_width
     real*8 :: max_exp !< Maximum possible exponent
     real*8 :: this_exp !< Exponent in this step
     real*8, parameter :: n_spin = 2d0 !< Non spin-polarized case supported only
     integer :: i_state !< State index

     character*40, parameter :: caller = "elsi_get_ne"

     invert_width = 1d0/broadening_width
     diff_ne_out = -n_electrons

     select case (broadening_type)
        case(GAUSSIAN)
           do i_state = 1,n_states
              occ_elpa(i_state) = &
                 n_spin*0.5d0*(1-erf((eigenvalues(i_state)-mu_in)*invert_width))
              diff_ne_out = diff_ne_out+occ_elpa(i_state)
           enddo

        case(FERMI)
           max_exp = maxexponent(mu_in)*log(2d0)
           do i_state = 1,n_states
              this_exp = (eigenvalues(i_state)-mu_in)*invert_width
              if(this_exp < max_exp) then
                 occ_elpa(i_state) = n_spin/(1+exp(this_exp))
                 diff_ne_out = diff_ne_out+occ_elpa(i_state)
              else ! Exponent in this step is larger than the largest possible exponent
                 occ_elpa(i_state) = 0d0
              endif
           enddo

        case DEFAULT
           call elsi_stop(" No supperted broadening type has been chosen. "//&
                          " Exiting...", caller)
     end select

  end subroutine ! elsi_get_ne

!>
!! This routine computes the chemical potential using bisection algorithm.
!!
  subroutine elsi_get_mu(mu_lower_in,mu_upper_in,mu_out)

     implicit none

     real*8, intent(in)  :: mu_lower_in !< Lower bound of chemical potential
     real*8, intent(in)  :: mu_upper_in !< Upper bound of chemical potential
     real*8, intent(out) :: mu_out      !< Solution of chemical potential

     real*8  :: mu_left !< Left bound of chemical potential interval
     real*8  :: mu_right !< Right bound of chemical potential interval
     real*8  :: mu_mid !< Middle point of chemical potential interval
     real*8  :: diff_left !< Difference in number of electrons on left bound
     real*8  :: diff_right !< Difference in number of electrons on right bound
     real*8  :: diff_mid !< Difference in number of electrons on middle point
     logical :: found_mu !< Is chemical potential found?
     integer :: n_steps !< Number of steps to find chemical potential
     integer, parameter :: max_steps = 100 !< Maximum steps to find chemical potential

     character*40, parameter :: caller = "elsi_get_mu"

     n_steps = 0
     found_mu = .false.

     mu_left = mu_lower_in
     mu_right = mu_upper_in

     do while(.not.found_mu)
        call elsi_get_ne(mu_left, diff_left)
        call elsi_get_ne(mu_right, diff_right)

        if(abs(diff_left) < occ_tolerance) then
           mu_out = mu_left
           found_mu = .true.
        elseif(abs(diff_right) < occ_tolerance) then
           mu_out = mu_right
           found_mu = .true.
        else
           mu_mid = 0.5d0*(mu_left+mu_right)

           n_steps = n_steps+1
           if(n_steps > max_steps) then
              call elsi_stop(" Chemical potential not found in 100 iterations! "//&
                             " Exiting...", caller)
           endif

           call elsi_get_ne(mu_mid, diff_mid)

           if(abs(diff_mid) < occ_tolerance) then
              mu_out = mu_mid
              found_mu = .true.
           elseif(diff_mid < 0) then
              mu_left = mu_mid
           elseif(diff_mid > 0) then
              mu_right = mu_mid
           endif
        endif
     enddo

     if(myid == 0) then
        write(*,"(A,F15.5,A)") "  | Chemical potential = ", mu_out, " Ha"
        write(*,"(A,F15.5,A)") "  |                    = ", mu_out*hartree, " eV"
     endif

  end subroutine ! elsi_get_mu

!>
!! This routine constructs the density matrix using eigenvectors from ELPA.
!!
  subroutine elsi_compute_dm_elpa(D_out)

     implicit none

     real*8, intent(out) :: D_out(n_l_rows, n_l_cols) !< Density matrix

     real*8, allocatable     :: tmp_real(:,:)    !< Real eigenvectors, temporary
     complex*16, allocatable :: tmp_complex(:,:) !< Complex eigenvectors, temporary

     real*8 :: D_out_tmp(n_l_rows, n_l_cols) !< Density matrix from imaginary 
                                             !< part of complex eigenvectors

     real*8, allocatable :: factor(:) !< Factor to construct density matrix

     integer, allocatable :: local_col(:)
     integer :: i_col, i
     integer :: l_row, l_col !< Local index
     integer :: g_row, g_col !< Global index

     character*40, parameter :: caller = "elsi_compute_dm"

     select case (method)
        case (ELPA)
           if(.not.allocated(D_elpa)) then
              call elsi_allocate(D_elpa,n_l_rows,n_l_cols,"D_elpa",caller)
           endif
           D_elpa = 0d0

           ! Map global columns to local
           call elsi_allocate(local_col,n_g_size,"local_col",caller)

           i_col = 0 ! local column counter

           do i = 1,n_g_size
              if(MOD((i-1)/n_b_cols,n_p_cols) == my_p_col) then
                 i_col = i_col+1
                 local_col(i) = i_col
              endif
           enddo

           select case (mode)
              case (REAL_VALUES)
                 ! Get eigenvectors into tmp_real
                 call elsi_allocate(tmp_real,n_l_rows,n_l_cols,"tmp_real",caller)
                 call elsi_get_eigenvectors(tmp_real)

                 ! Compute the factors used to construct density matrix
                 call elsi_allocate(factor,n_states,"factor",caller)
                 factor = 0d0

                 do i = 1,n_states
                    if(occ_elpa(i) > 0d0) then
                       factor(i) = SQRT(occ_elpa(i))
                    endif
                 enddo

                 do i = 1,n_states
                    if(factor(i) > 0d0) then
                       if(local_col(i) > 0) then
                          tmp_real(:,local_col(i)) = tmp_real(:,local_col(i))*factor(i)
                       endif
                    elseif(local_col(i) .ne. 0) then
                       tmp_real(:,local_col(i)) = 0d0
                    endif
                 enddo

                 ! Compute density matrix
                 D_out = 0d0

                 ! D_out = tmp_real*tmp_real^T
                 call pdsyrk('U','N',n_g_size,n_states,1d0,tmp_real,1,1,sc_desc,&
                             0d0,D_out,1,1,sc_desc)

              case (COMPLEX_VALUES)
                 ! Get eigenvectors into tmp_complex
                 call elsi_allocate(tmp_complex,n_l_rows,n_l_cols,"tmp_complex",caller)
                 call elsi_get_eigenvectors(tmp_complex)

                 ! Compute the factors used to construct density matrix
                 call elsi_allocate(factor,n_states,"factor",caller)
                 factor = 0d0

                 do i = 1,n_states
                    if(occ_elpa(i) > 0d0) then
                       factor(i) = SQRT(occ_elpa(i))
                    endif
                 enddo

                 do i = 1,n_states
                    if(factor(i) > 0d0) then
                       if(local_col(i) > 0) then
                          tmp_complex(:,local_col(i)) = tmp_complex(:,local_col(i))*factor(i)
                       endif
                    elseif(local_col(i) .ne. 0) then
                       tmp_complex(:,local_col(i)) = 0d0
                    endif
                 enddo

                 ! Compute density matrix
                 D_out = 0d0
                 D_out_tmp = 0d0

                 ! D_out = tmp_complex*tmp_complex' A' here means transpose of A
                 !call pzherk('U','N',n_g_size,n_states,(1d0,0d0),tmp_complex,1,1,sc_desc,&
                 !            (0d0,0d0),D_out_tmp,1,1,sc_desc)

                 call pdsyrk('U','N',n_g_size,n_states,1d0,real(tmp_complex),1,1,sc_desc,&
                             0d0,D_out,1,1,sc_desc)
                 call pdsyrk('U','N',n_g_size,n_states,1d0,aimag(tmp_complex),1,1,sc_desc,&
                             0d0,D_out_tmp,1,1,sc_desc)

                 D_out = D_out+D_out_tmp
           end select

           deallocate(local_col)
           deallocate(factor)
           if(allocated(tmp_real))    deallocate(tmp_real)
           if(allocated(tmp_complex)) deallocate(tmp_complex)

           ! Now D_out is an upper triangle matrix
           ! Set D_out to full
           call elsi_allocate(tmp_real,n_l_rows,n_l_cols,"tmp_real",caller)
           tmp_real = D_out

           ! D_out = D_out + tmp_real' = D_out + D_out'
           call pdtran(n_g_size,n_g_size,1d0,tmp_real,1,1,sc_desc,1d0,D_out,1,1,sc_desc)

           deallocate(tmp_real)

           do l_row = 1,n_l_rows
              call elsi_get_global_row(g_row,l_row)
              do l_col = l_row,n_l_cols
                 call elsi_get_global_col(g_col,l_col)
                 if(g_row == g_col) then
                    D_out(l_row,l_col) = 0.5d0*D_out(l_row,l_col)
                 endif
              enddo
           enddo

        case (LIBOMM)
           call elsi_stop(" LIBOMM does not compute density matrix from eigenvectors! "//&
                          " Exiting...",caller)
        case (PEXSI)
           call elsi_stop(" PEXSI does not compute density matrix from eigenvectors! "//&
                          " Exiting...",caller)
        case (CHESS)
           call elsi_stop(" CHESS does not compute density matrix from eigenvectors! "//&
                          " Exiting...",caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...",caller)
     end select

  end subroutine ! elsi_compute_dm_elpa

!> 
!! This routine transforms a generalized eigenvalue problem (Ac = Bcv)
!! to standard form (A'c' = c'v)
!!
!! Starting from Hv = eSv, we first perform a Cholesky decomposition of S
!! S = (U^T)U, resulting in Hv = e(U^T)Uv
!!
!! Using 1=U^-1U we define a new standard eigenvalue problem by
!! H(U^-1)(Uv) = e(U^T)(Uv) => ((U^-1)^T)H(U^-1)(Uv) = e(Uv)
!!
!! On exit, (U^-1) is stored in S, to be used for back-transformation
!!
  subroutine elsi_to_standard_evp()


     real*8 :: ev_sqrt
     real*8, allocatable :: ev_overlap(:)
     real*8, allocatable :: buffer_real(:,:)
     complex*16, allocatable :: buffer_complex(:,:)
     integer, allocatable :: local_col(:)
     integer :: i,i_col
     logical :: success

     character*40, parameter :: caller = "elsi_to_standard_evp"

     select case (method)
        case (ELPA)
           select case (mode)
              case (COMPLEX_VALUES)
                 call elsi_allocate(buffer_complex,n_l_rows,n_l_cols,"temp",caller)

                 if(n_elsi_calls == 1) then
                    ! Check if overlap matrix is singular
                    call elsi_statement_print("  Checking singularity for overlap matrix")

                    ! Use buffer_complex to store overlap matrix, otherwise it will
                    ! be destroyed by eigenvalue calculation
                    ! The nonsingular eigenvalues must be the first ones, so find
                    ! eigenvalues of negative overlap matrix
                    buffer_complex = -S_complex

                    call elsi_allocate(ev_overlap,n_g_size,"ev_overlap",caller)

                    ! Always use ELPA 2-stage solver here
                    success = solve_evp_complex_2stage_double(n_g_size,n_g_size,&
                                 buffer_complex,n_l_rows,ev_overlap,C_complex,&
                                 n_l_rows,n_b_rows,n_l_cols,mpi_comm_row,&
                                 mpi_comm_col,mpi_comm_global)

                    if(.not.success) then
                       call elsi_stop(" ELPA failed when solving eigenvalue problem. "//&
                                      " Exiting...", caller)
                    endif

                    ! Invert signs of eigenvalues
                    ev_overlap = -ev_overlap

                    ! Get the number of nonsingular eigenvalues
                    do i = 1,n_g_size
                       if(ev_overlap(i) < singularity_tolerance) exit
                    enddo
                    n_nonsingular = i-1

                    ! Stop if n_states is larger that n_nonsingular
                    if(n_nonsingular < n_states) then ! Too singular to continue
                       call elsi_stop(" Overlap matrix is singular. The number of"//&
                                      " basis functions after removing singularity"//&
                                      " is smaller than the number of states. Try to"//&
                                      " a) decrease the size of basis set, or b)"//&
                                      " decrease the number of states, or c) increase"//&
                                      " the tolerance of basis singularity."//&
                                      " Exiting...",caller)
                    elseif(n_nonsingular < n_g_size) then ! Singular
                       overlap_is_singular = .true.

                       call elsi_statement_print("  Overlap matrix is singular. This"//&
                                                 " may mean that a very large basis"//&
                                                 " set is in use. The calculation"//&
                                                 " will continue. However, please"//&
                                                 " note that running with a near-"//&
                                                 " singular basis set may lead to"//&
                                                 " completely wrong numerical results.")

                       if(myid == 0) then
                          write(*,"(A,I13)") "  | Number of basis functions reduced to: ",&
                                             n_nonsingular
                       endif

                       call elsi_statement_print("  Using scaled eigenvectors of"//&
                                                 " overlap matrix for transformation")

                       ! Overlap matrix is overwritten with scaled eigenvectors
                       ! Map global columns to local
                       call elsi_allocate(local_col,n_g_size,"local_col",caller)

                       i_col = 0 ! local column counter
                       do i = 1,n_g_size
                          if(MOD((i-1)/n_b_cols,n_p_cols) == my_p_col) then
                             i_col = i_col+1
                             local_col(i) = i_col
                          endif
                       enddo

                       ! Scale eigenvectors
                       do i = 1,n_nonsingular
                          ev_sqrt = SQRT(ev_overlap(i))
                          if(local_col(i) == 0) cycle
                          S_complex(:,local_col(i)) = C_complex(:,local_col(i))/ev_sqrt
                       enddo

                    else ! Nonsingular
                       overlap_is_singular = .false.

                       call elsi_statement_print("  Overlap matrix is nonsingular")
                       call elsi_statement_print("  Starting Cholesty decomposition")

                       ! Compute S = (U^T)U, U -> S
                       success = elpa_cholesky_complex_double(n_g_size,S_complex,n_l_rows,&
                                    n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
                       if(.not.success) then
                          call elsi_stop(" Cholesky decomposition failed.",caller)
                       endif

                       ! compute U^-1 -> S
                       success = elpa_invert_trm_complex_double(n_g_size,S_complex,n_l_rows,&
                                    n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
                       if(.not.success) then
                          call elsi_stop(" Matrix invertion failed.", caller)
                       endif
                    endif ! Singular overlap?
                 endif ! First call?

                 if(overlap_is_singular) then ! Use scaled eigenvectors
                    ! buffer_complex = H_complex * S_complex
                    call pzgemm('N','N',n_g_size,n_nonsingular,n_g_size,(1d0,0d0),&
                                H_complex,1,1,sc_desc,S_complex,1,1,sc_desc,(0d0,0d0),&
                                buffer_complex,1,1,sc_desc)

                    ! H_complex = (S_complex)^T * buffer_complex
                    call pzgemm('C','N',n_nonsingular,n_nonsingular,n_g_size,(1d0,0d0),&
                                S_complex,1,1,sc_desc,buffer_complex,1,1,sc_desc,&
                                (0d0,0d0),H_complex,1,1,sc_desc)

                 else ! Use cholesky
                    ! compute H(U^-1) -> buff
                    ! buffer_complex = H_complex * S_complex
                    call pzgemm('C','N',n_g_size,n_g_size,n_g_size,(1d0,0d0),&
                                S_complex,1,1,sc_desc,H_complex,1,1,sc_desc,&
                                (0d0,0d0),buffer_complex,1,1,sc_desc)

                    ! compute ((U^-1)^T)H by (H(U^-1))^T -> H
                    ! H_complex = (buffer_complex)^T
                    call pztranc(n_g_size,n_g_size,(1d0,0d0),buffer_complex,1,1,&
                                 sc_desc,(0d0,0d0),H_complex,1,1,sc_desc)

                    ! compute ((U^-1)^T)H(U^-1) -> H
                    buffer_complex = H_complex
                    ! H_complex = buffer_complex * S_complex
                    call pzgemm('C','N',n_g_size,n_g_size,n_g_size,(1d0,0d0),&
                                S_complex,1,1,sc_desc,buffer_complex,1,1,&
                                sc_desc,(0d0,0d0),H_complex,1,1,sc_desc)
                 endif

              case (REAL_VALUES)
                 call elsi_allocate(buffer_real,n_l_rows,n_l_cols,"temp",caller)

                 if(n_elsi_calls == 1) then
                    ! Check if overlap matrix is singular
                    call elsi_statement_print("  Checking singularity for overlap matrix")

                    ! Use buffer_real to store overlap matrix, otherwise it will be
                    ! destroyed by eigenvalue calculation
                    ! The nonsingular eigenvalues must be the first ones, so find
                    ! eigenvalues of negative overlap matrix
                    buffer_real = -S_real

                    call elsi_allocate(ev_overlap,n_g_size,"ev_overlap",caller)

                    ! Always use ELPA 2-stage solver here
                    success = solve_evp_real_2stage_double(n_g_size,n_g_size,buffer_real,&
                                 n_l_rows,ev_overlap,C_real,n_l_rows,n_b_rows,n_l_cols,&
                                 mpi_comm_row,mpi_comm_col,mpi_comm_global)

                    if(.not.success) then
                       call elsi_stop(" ELPA failed when solving eigenvalue problem. "//&
                                      " Exiting...", caller)
                    endif

                    ! Invert signs of eigenvalues
                    ev_overlap = -ev_overlap

                    ! Get the number of nonsingular eigenvalues
                    do i = 1,n_g_size
                       if(ev_overlap(i) < singularity_tolerance) exit
                    enddo
                    n_nonsingular = i-1

                    ! Stop if n_states is larger that n_nonsingular
                    if(n_nonsingular < n_states) then ! Too singular to continue
                       call elsi_stop(" Overlap matrix is singular. The number of"//&
                                      " basis functions after removing singularity"//&
                                      " is smaller than the number of states. Try to"//&
                                      " a) decrease the size of basis set, or b)"//&
                                      " decrease the number of states, or c) increase"//&
                                      " the tolerance of basis singularity."//&
                                      " Exiting...",caller)
                    elseif(n_nonsingular < n_g_size) then ! Singular
                       overlap_is_singular = .true.

                       call elsi_statement_print("  Overlap matrix is singular. This"//&
                                                 " may mean that a very large basis"//&
                                                 " set is in use. The calculation"//&
                                                 " will continue. However, please"//&
                                                 " note that running with a near-"//&
                                                 " singular basis set may lead to"//&
                                                 " completely wrong numerical results.")

                       if(myid == 0) then
                          write(*,"(A,I13)") "  | Number of basis functions reduced to: ",&
                                             n_nonsingular
                       endif

                       call elsi_statement_print("  Using scaled eigenvectors of"//&
                                                 " overlap matrix for transformation")

                       ! Overlap matrix is overwritten with scaled eigenvectors
                       ! Map global columns to local
                       call elsi_allocate(local_col,n_g_size,"local_col",caller)

                       i_col = 0 ! local column counter
                       do i = 1,n_g_size
                          if(MOD((i-1)/n_b_cols,n_p_cols) == my_p_col) then
                             i_col = i_col+1
                             local_col(i) = i_col
                          endif
                       enddo

                       ! Scale eigenvectors
                       do i = 1,n_nonsingular
                          ev_sqrt = SQRT(ev_overlap(i))
                          if(local_col(i) == 0) cycle
                          S_real(:,local_col(i)) = C_real(:,local_col(i))/ev_sqrt
                       enddo

                    else ! Nonsingular
                       overlap_is_singular = .false.

                       call elsi_statement_print("  Overlap matrix is nonsingular")
                       call elsi_statement_print("  Starting Cholesty decomposition")

                       ! Compute S = (U^T)U, U -> S
                       success = elpa_cholesky_real_double(n_g_size,S_real,n_l_rows,&
                                    n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
                       if(.not.success) then
                          call elsi_stop(" Cholesky decomposition failed.",caller)
                       endif

                       ! compute U^-1 -> S
                       success = elpa_invert_trm_real_double(n_g_size,S_real,n_l_rows,&
                                    n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
                       if(.not.success) then
                          call elsi_stop(" Matrix invertion failed.",caller)
                       endif
                    endif ! Singular overlap?
                 endif ! First call?

                 if(overlap_is_singular) then ! Use scaled eigenvectors
                    ! buffer_real = H_real * S_real
                    call pdgemm('N','N',n_g_size,n_nonsingular,n_g_size,1d0,H_real,1,&
                                1,sc_desc,S_real,1,1,sc_desc,0d0,buffer_real,1,1,sc_desc)

                    ! H_real = (S_real)^T * buffer_real
                    call pdgemm('T','N',n_nonsingular,n_nonsingular,n_g_size,1d0,S_real,&
                                1,1,sc_desc,buffer_real,1,1,sc_desc,0d0,H_real,1,1,sc_desc)

                 else ! Use Cholesky
                    ! compute H(U^-1) -> buff
                    ! buffer_real = H_real * S_real
                    call pdgemm('N','N',n_g_size,n_g_size,n_g_size,1d0,H_real,1,1,&
                                sc_desc,S_real,1,1,sc_desc,0d0,buffer_real,1,1,sc_desc)

                    ! compute ((U^-1)^T)H by (H(U^-1))^T -> H
                    ! H_real = (buffer_real)^T
                    call pdtran(n_g_size,n_g_size,1d0,buffer_real,1,1,sc_desc,0d0,&
                                H_real,1,1,sc_desc)

                    ! compute ((U^-1)^T)H(U^-1) -> H
                    buffer_real = H_real
                    ! H_real = buffer_real * S_real
                    call pdgemm('N','N',n_g_size,n_g_size,n_g_size,1d0,buffer_real,&
                                1,1,sc_desc,S_real,1,1,sc_desc,0d0,H_real,1,1,sc_desc)
                 endif

           end select

        case (LIBOMM)
           call elsi_stop(" OMM does not need to transform evp. Exiting...",caller)
        case (PEXSI)
           call elsi_stop(" PEXSI does not need to transform evp. Exiting...",caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...",caller)
     end select

     if(allocated(ev_overlap))     deallocate(ev_overlap)
     if(allocated(local_col))      deallocate(local_col)
     if(allocated(buffer_real))    deallocate(buffer_real)
     if(allocated(buffer_complex)) deallocate(buffer_complex)

  end subroutine ! elsi_to_standard_evp

!> 
!! This routine does the back-transformation of the eigenvectors in standard
!! form (A'c' = c'v) to the original generalized form (Ac = Bcv)
!!
!! v = (U^-1)v'
!!
  subroutine elsi_to_original_ev()

     real*8, allocatable :: buffer_real(:,:)
     complex*16, allocatable :: buffer_complex(:,:)

     character*40, parameter :: caller = "elsi_to_original_ev"

     select case (method)
        case (ELPA)
           select case (mode)
              case (COMPLEX_VALUES)
                 call elsi_allocate(buffer_complex,n_l_rows,n_l_cols,"temp",caller)
                 buffer_complex = C_complex

                 if(overlap_is_singular) then
                    ! Transform matrix is stored in S_complex after elsi_to_standard_evp
                    call pzgemm('N','N',n_g_size,n_states,n_nonsingular,(1d0,0d0),&
                                S_complex,1,1,sc_desc,buffer_complex,1,1,sc_desc,&
                                (0d0,0d0),C_complex,1,1,sc_desc)
                 else ! Nonsingular, use Cholesky
                    ! (U^-1) is stored in S_complex after elsi_to_standard_evp
                    ! C_complex = S_complex * C_complex = S_complex * buffer_complex
                    call pzgemm('N','N',n_g_size,n_states,n_g_size,(1d0,0d0),S_complex,&
                                1,1,sc_desc,buffer_complex,1,1,sc_desc,(0d0,0d0),&
                                C_complex,1,1,sc_desc)
                 endif

              case (REAL_VALUES)
                 call elsi_allocate(buffer_real,n_l_rows,n_l_cols,"temp",caller)
                 buffer_real = C_real

                 if(overlap_is_singular) then
                    ! Transform matrix is stored in S_real after elsi_to_standard_evp
                    call pdgemm('N','N',n_g_size,n_states,n_nonsingular,1d0,S_real,1,&
                                1,sc_desc,buffer_real,1,1,sc_desc,0d0,C_real,1,1,sc_desc)
                 else ! Nonsingular, use Cholesky
                    ! (U^-1) is stored in S_real after elsi_to_standard_evp
                    ! C_real = S_real * C_real = S_real * buffer_real
                    call pdgemm('N','N',n_g_size,n_states,n_g_size,1d0,S_real,1,1,&
                                sc_desc,buffer_real,1,1,sc_desc,0d0,C_real,1,1,sc_desc)
                 endif

           end select

        case (LIBOMM)
           call elsi_stop(" OMM does not have eigenvectors. Exiting...",caller)
        case (PEXSI)
           call elsi_stop(" PEXSI does not have eigenvectors. Exiting...",caller)
        case DEFAULT
           call elsi_stop(" No supported method has been chosen. "//&
                          " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                          " Exiting...",caller)
     end select

     if(allocated(buffer_real))    deallocate(buffer_real)
     if(allocated(buffer_complex)) deallocate(buffer_complex)

  end subroutine ! elsi_to_original_ev

!>
!! This routine interfaces to ELPA.
!!
  subroutine elsi_solve_evp_elpa()

     implicit none
     include "mpif.h"

     logical :: success
     logical :: two_step_solver

     character*40, parameter :: caller = "elsi_solve_evp_elpa"

     call elsi_start_solve_evp_time()

     ! Choose 1-stage or 2-stage solver
     if(elpa_one_always) then
        two_step_solver = .false.
     elseif(elpa_two_always) then
        two_step_solver = .true.
     elseif(n_g_size < 256) then
        two_step_solver = .false.
     else
        two_step_solver = .true.
     endif

     ! Transform to standard form
     if(.not.overlap_is_unit) then
        call elsi_statement_print("  Tansforming to standard evp")
        call elsi_to_standard_evp()
     endif

     ! Solve evp, return eigenvalues and eigenvectors
     if(two_step_solver) then ! 2-stage solver
        call elsi_statement_print("  Starting ELPA 2-stage solver")
        select case (mode)
           case (COMPLEX_VALUES)
              success = solve_evp_complex_2stage_double(n_nonsingular,n_states,&
                           H_complex,n_l_rows,eigenvalues,C_complex,n_l_rows,n_b_rows,&
                           n_l_cols,mpi_comm_row,mpi_comm_col,mpi_comm_global)
           case (REAL_VALUES)
              success = solve_evp_real_2stage_double(n_nonsingular,n_states,H_real,&
                           n_l_rows,eigenvalues,C_real,n_l_rows,n_b_rows,&
                           n_l_cols,mpi_comm_row,mpi_comm_col,mpi_comm_global)
        end select
     else ! 1-stage solver
        call elsi_statement_print("  Starting ELPA 1-stage solver")
        select case (mode)
           case (COMPLEX_VALUES)
              success = solve_evp_complex_1stage_double(n_nonsingular,n_states,&
                           H_complex,n_l_rows,eigenvalues,C_complex,n_l_rows,&
                           n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col)
           case (REAL_VALUES)
              success = solve_evp_real_1stage_double(n_nonsingular,n_states,H_real,&
                           n_l_rows,eigenvalues,C_real,n_l_rows,n_b_rows,n_l_cols,&
                           mpi_comm_row,mpi_comm_col)
        end select
     endif

     if(.not.success) then
        call elsi_stop(" ELPA failed when solving eigenvalue problem. "//&
                       " Exiting...",caller)
     endif

     ! Back-transform eigenvectors
     if(.not.overlap_is_unit) then
        call elsi_statement_print("  Transforming to original eigenvectors")
        call elsi_to_original_ev()
     endif

     call MPI_Barrier(mpi_comm_global,mpierr)
     call elsi_stop_solve_evp_time()

  end subroutine ! elsi_solve_evp_elpa

!>
!! This routine overrides ELPA default settings.
!!
  subroutine elsi_customize_elpa(elpa_solver)

     implicit none

     integer, intent(in) :: elpa_solver !< Always use 1-stage or 2-stage solver

     if(method == ELPA) then
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
     else
        call elsi_statement_print("  The chosen method is not ELPA."//&
                                  " Ignore elsi_customize_elpa call.")
     endif

  end subroutine ! elsi_customize_elpa

!=======================
! ELSI routines for OMM
!=======================

!>
!! This routine interfaces to libOMM.
!!
  subroutine elsi_solve_evp_omm()

     implicit none
     include "mpif.h"

     logical, save :: first_call = .true.
     logical :: success
     character*40, parameter :: caller = "elsi_solve_evp_omm"

     if(overlap_is_singular) then
        call elsi_stop(" libOMM cannot treat singular overlap matrix yet. "//&
                       " Exiting...", caller)
     endif

     call elsi_start_solve_evp_time()

     if(n_elsi_calls == 1) then
        C_matrix_initialized = .false.
        ! Cholesky
        select case (mode)
           case (COMPLEX_VALUES)
              ! Compute S = (U^T)U, U -> S
              success = elpa_cholesky_complex_double(n_g_size,S_omm%zval,n_l_rows,&
                           n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
           case (REAL_VALUES)
              ! Compute S = (U^T)U, U -> S
              success = elpa_cholesky_real_double(n_g_size,S_omm%dval,n_l_rows,&
                           n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
        end select
     else
        C_matrix_initialized = .true.
     endif

     if(first_call) then
        new_overlap = .true.
     else
        new_overlap = .false.
     endif

     if(n_elpa_steps > 0 .and. n_elsi_calls == n_elpa_steps+1) then
        ! Invert U^(-1), which is used in ELPA cholesky and stored in S, to get U for OMM
        select case (mode)
           case (COMPLEX_VALUES)
              success = elpa_invert_trm_complex_double(n_g_size,S_omm%zval,n_l_rows,&
                           n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
           case (REAL_VALUES)
              success = elpa_invert_trm_real_double(n_g_size,S_omm%dval,n_l_rows,&
                           n_b_rows,n_l_cols,mpi_comm_row,mpi_comm_col,.false.)
        end select
     endif

     ! Shift eigenvalue spectrum
!    call m_add(S_elsi,'N',H_elsi,-eta,1d0,"lap")

     select case (mode)
        case (COMPLEX_VALUES)
           call omm(n_g_size,n_states,H_omm,S_omm,new_overlap,total_energy,D_omm,&
                    calc_ED,eta,Coeff_omm,C_matrix_initialized,T_omm,scale_kinetic,&
                    omm_flavour,nk_times_nspin,i_k_spin,min_tol,omm_verbose,&
                    do_dealloc,"pzdbc","lap",myid)
        case (REAL_VALUES)
           call omm(n_g_size,n_states,H_omm,S_omm,new_overlap,total_energy,D_omm,&
                    calc_ED,eta,Coeff_omm,C_matrix_initialized,T_omm,scale_kinetic,&
                    omm_flavour,nk_times_nspin,i_k_spin,min_tol,omm_verbose,&
                    do_dealloc,"pddbc","lap",myid)
     end select

     first_call = .false.

     call MPI_Barrier(mpi_comm_global,mpierr)
     call elsi_stop_solve_evp_time()

  end subroutine ! elsi_solve_evp_omm

!>
!! This routine overrides libOMM default settings.
!!
  subroutine elsi_customize_omm(n_elpa_steps_in,scale_kinetic_in,calc_ed_in,&
                                eta_in,min_tol_in)

     implicit none

     integer, intent(in), optional :: n_elpa_steps_in
     real*8,  intent(in), optional :: scale_kinetic_in
     logical, intent(in), optional :: calc_ed_in
     real*8,  intent(in), optional :: eta_in
     real*8,  intent(in), optional :: min_tol_in

     if(method == LIBOMM) then
        ! Set default settings
        call elsi_set_omm_default_options()

        ! Number of ELPA steps
        if(present(n_elpa_steps_in)) n_elpa_steps = n_elpa_steps_in
        ! Scaling of kinetic energy matrix
        if(present(scale_kinetic_in)) scale_kinetic = scale_kinetic_in
        ! Calculate energy weigthed density matrix?
        if(present(calc_ed_in)) calc_ed = calc_ed_in
        ! Eigenspectrum shift parameter
        if(present(eta_in)) eta = eta_in
        ! Tolerance for minimization
        if(present(min_tol_in)) min_tol = min_tol_in

        omm_customized = .true.
        call elsi_print_omm_options()
     else
        call elsi_statement_print("  The chosen method is not OMM."//&
                                  " Ignore elsi_customize_omm call.")
     endif

  end subroutine ! elsi_customize_omm

!=========================
! ELSI routines for PEXSI
!=========================

!>
!! PEXSI processor grid setup.
!!
  subroutine elsi_init_pexsi()

     implicit none

     character*40, parameter :: caller = "elsi_init_pexsi"

     if(method == PEXSI) then
        if(mod(n_procs,pexsi_options%numPole) == 0) then
           n_p_per_pole_pexsi = n_procs/pexsi_options%numPole
           call elsi_statement_print("  PEXSI parallel over poles.")
           if(myid == 0) &
              write(*,"(A,I13)") "  | Number of MPI tasks per pole: ", &
                    n_p_per_pole_pexsi
        else
           n_p_per_pole_pexsi = n_procs
           call elsi_statement_print("  PEXSI not parallel over poles. High performance"//&
                                     " is expected with number of MPI tasks being a"//&
                                     " multiple of number of poles.")
        endif

        ! Set square-like process grid for selected inversion of each pole
        do n_p_rows_pexsi = NINT(SQRT(REAL(n_p_per_pole_pexsi))),2,-1
           if(mod(n_p_per_pole_pexsi,n_p_rows_pexsi) == 0) exit
        enddo

        n_p_cols_pexsi = n_p_per_pole_pexsi/n_p_rows_pexsi

        ! PEXSI process grid
        my_p_col_pexsi = mod(myid,n_p_per_pole_pexsi)
        my_p_row_pexsi = myid/n_p_per_pole_pexsi

        ! PEXSI uses a pure block distribution in the first process row
        n_b_rows_pexsi = n_g_size

        ! The last process holds all remaining columns
        n_b_cols_pexsi = FLOOR(1d0*n_g_size/n_p_per_pole_pexsi)
        if(my_p_col_pexsi == n_p_per_pole_pexsi-1) then
           n_b_cols_pexsi = n_g_size-(n_p_per_pole_pexsi-1)*n_b_cols_pexsi
        endif

        n_l_rows_pexsi = n_b_rows_pexsi
        n_l_cols_pexsi = n_b_cols_pexsi

        ! Only master process outputs
        if(myid == 0) then
           pexsi_output_file_index = 0
        else
           pexsi_output_file_index = -1
        endif

        pexsi_plan = f_ppexsi_plan_initialize(mpi_comm_global,n_p_rows_pexsi,&
                        n_p_cols_pexsi,pexsi_output_file_index,pexsi_info)

        if(pexsi_info /= 0) &
           call elsi_stop(" PEXSI plan initialization failed. Exiting...",caller)

     endif

  end subroutine ! elsi_init_pexsi

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by PEXSI.
!!
  subroutine elsi_2dbcd_to_1dbccs_hs_pexsi(H_in,S_in)

     implicit none
     include "mpif.h"

     real*8, intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian matrix to be converted
     real*8, intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap matrix to be converted

     integer :: i_row !< Row counter
     integer :: i_col !< Col counter
     integer :: i_val !< Value counter
     integer :: i_proc !< Process counter
     integer :: global_col_id !< Global column id
     integer :: global_row_id !< Global row id
     integer :: local_col_id !< Local column id in 1D block distribution
     integer :: local_row_id !< Local row id in 1D block distribution
     integer :: nnz_l_tmp
     integer :: mpi_comm_aux_pexsi
     integer :: mpierr

     integer, allocatable :: dest(:) !< Destination of each element
     real*8 :: matrix_aux(n_l_rows_pexsi, n_l_cols_pexsi)

     ! For the meaning of each array here, see documentation of MPI_Alltoallv
     real*8, allocatable  :: h_val_send_buffer(:) !< Send buffer for Hamiltonian
     real*8, allocatable  :: s_val_send_buffer(:) !< Send buffer for overlap
     integer, allocatable :: pos_send_buffer(:)   !< Send buffer for global 1D id
     integer :: send_count(n_procs) !< Number of elements to send to each processor
     integer :: send_displ(n_procs) !< Displacement from which to take the outgoing data
     integer :: send_displ_aux      !< Auxiliary variable used to set displacement

     real*8, allocatable  :: h_val_recv_buffer(:) !< Receive buffer for Hamiltonian
     real*8, allocatable  :: s_val_recv_buffer(:) !< Receive buffer for overlap
     integer, allocatable :: pos_recv_buffer(:)   !< Receive buffer for global 1D id
     integer :: recv_count(n_procs) !< Number of elements to receive from each processor
     integer :: recv_displ(n_procs) !< Displacement at which to place the incoming data
     integer :: recv_displ_aux      !< Auxiliary variable used to set displacement

     character*40, parameter :: caller = "elsi_2dbcd_to_1dbccs_hs_pexsi"

     call elsi_start_2dbc_to_1dccs_time()
     call elsi_statement_print("  Matrix conversion: 2D block-cyclic dense ==> 1D block CCS sparse")

     send_count = 0
     send_displ = 0
     recv_count = 0
     recv_displ = 0

     call elsi_get_local_nnz(H_in, n_l_rows, n_l_cols, nnz_l)

     call elsi_allocate(dest, nnz_l, "dest", caller)
     call elsi_allocate(pos_send_buffer, nnz_l, "pos_send_buffer", caller)
     call elsi_allocate(h_val_send_buffer, nnz_l, "h_val_send_buffer", caller)
     if(.not.overlap_is_unit) then
        call elsi_allocate(s_val_send_buffer, nnz_l, "s_val_send_buffer", caller)
     endif

     i_val = 0
     ! Compute destination and global 1D id
     do i_col = 1, n_l_cols
        do i_row = 1, n_l_rows
           if(abs(H_in(i_row, i_col)) > zero_threshold) then
              i_val = i_val + 1
              call elsi_get_global_col(global_col_id, i_col)
              call elsi_get_global_row(global_row_id, i_row)

              ! Compute destination
              dest(i_val) = FLOOR(1d0*(global_col_id-1)/FLOOR(1d0*n_g_size/n_p_per_pole_pexsi))
              ! The last process may take more
              if(dest(i_val) > (n_p_per_pole_pexsi-1)) dest(i_val) = n_p_per_pole_pexsi-1

              ! Compute the global id
              ! Pack global id and data into buffers
              pos_send_buffer(i_val) = (global_col_id-1)*n_g_size+global_row_id
              h_val_send_buffer(i_val) = H_in(i_row, i_col)
              if(.not.overlap_is_unit) then
                 s_val_send_buffer(i_val) = S_in(i_row, i_col)
              endif
          endif
       enddo
     enddo

     ! Set send_count
     do i_proc = 0, n_procs-1
        do i_val = 1, nnz_l
           if(dest(i_val) == i_proc) then
              send_count(i_proc+1) = send_count(i_proc+1)+1
           endif
        enddo
     enddo

     deallocate(dest)

     ! Set recv_count
     call MPI_Alltoall(send_count, 1, mpi_integer, recv_count, &
                       1, mpi_integer, mpi_comm_global, mpierr)

     ! Set local/global number of nonzero
     nnz_l_pexsi = sum(recv_count, 1)
     call MPI_Allreduce(nnz_l_pexsi, nnz_g, 1, mpi_integer, mpi_sum, &
                        mpi_comm_global, mpierr)

     ! At this point only processes in the first row in PEXSI process gird
     ! have correct nnz_l_pexsi
     if(n_p_per_pole_pexsi < n_procs) then
        call MPI_Comm_split(mpi_comm_global, my_p_col_pexsi, my_p_row_pexsi, &
                            mpi_comm_aux_pexsi, mpierr)

        call MPI_Allreduce(nnz_l_pexsi, nnz_l_tmp, 1, mpi_integer, mpi_sum, &
                           mpi_comm_aux_pexsi, mpierr)

        nnz_l_pexsi = nnz_l_tmp
     endif

     ! Set send and receive displacement
     send_displ_aux = 0
     recv_displ_aux = 0

     do i_proc = 0, n_procs-1
        send_displ(i_proc+1) = send_displ_aux
        send_displ_aux = send_displ_aux + send_count(i_proc+1)

        recv_displ(i_proc+1) = recv_displ_aux
        recv_displ_aux = recv_displ_aux + recv_count(i_proc+1)
     enddo

     call elsi_allocate(pos_recv_buffer, nnz_l_pexsi, "pos_recv_buffer", caller)
     call elsi_allocate(h_val_recv_buffer, nnz_l_pexsi, "h_val_recv_buffer", caller)
     if(.not.overlap_is_unit) then
        call elsi_allocate(s_val_recv_buffer, nnz_l_pexsi, "s_val_recv_buffer", caller)
     endif

     ! Send and receive the packed data
     call MPI_Alltoallv(h_val_send_buffer, send_count, send_displ, mpi_real8, &
                        h_val_recv_buffer, recv_count, recv_displ, mpi_real8, &
                        mpi_comm_global, mpierr)

     if(.not.overlap_is_unit) then
        call MPI_Alltoallv(s_val_send_buffer, send_count, send_displ, mpi_real8, &
                           s_val_recv_buffer, recv_count, recv_displ, mpi_real8, &
                           mpi_comm_global, mpierr)
     endif

     call MPI_Alltoallv(pos_send_buffer, send_count, send_displ, mpi_integer, &
                        pos_recv_buffer, recv_count, recv_displ, mpi_integer, &
                        mpi_comm_global, mpierr)

     deallocate(pos_send_buffer)
     deallocate(h_val_send_buffer)
     if(.not.overlap_is_unit) then
        deallocate(s_val_send_buffer)
     endif

     ! TODO: double check the new algorithm and get rid of matrix_aux
     matrix_aux = 0d0

     ! Unpack Hamiltonian on the first process row in PEXSI process grid
     if(my_p_row_pexsi == 0) then
        do i_val = 1, nnz_l_pexsi
           ! Compute global 2d id
           global_col_id = FLOOR(1d0*(pos_recv_buffer(i_val)-1)/n_g_size)+1
           global_row_id = MOD(pos_recv_buffer(i_val), n_g_size)
           if(global_row_id == 0) global_row_id = n_g_size

           ! Compute local 2d id
           local_col_id = global_col_id-myid*FLOOR(1d0*n_g_size/n_p_per_pole_pexsi)
           local_row_id = global_row_id

           ! Put value to correct position
           matrix_aux(local_row_id,local_col_id) = h_val_recv_buffer(i_val)
        enddo
     endif

     deallocate(h_val_recv_buffer)

     ! Allocate PEXSI matrices
     if(.not.allocated(H_real_pexsi)) &
        call elsi_allocate(H_real_pexsi, nnz_l_pexsi, "H_real_pexsi", caller)
     H_real_pexsi = 0d0

     if(.not.overlap_is_unit) then
        if(.not.allocated(S_real_pexsi)) &
           call elsi_allocate(S_real_pexsi, nnz_l_pexsi, "S_real_pexsi", caller)
        S_real_pexsi = 0d0
     endif

     if(.not.allocated(row_ind_pexsi)) &
        call elsi_allocate(row_ind_pexsi, nnz_l_pexsi, "row_ind_pexsi", caller)
     row_ind_pexsi = 0

     if(.not.allocated(col_ptr_pexsi)) &
        call elsi_allocate(col_ptr_pexsi, (n_l_cols_pexsi+1), "col_ptr_pexsi", caller)
     col_ptr_pexsi = 0

     ! Transform Hamiltonian: 1D block dense ==> 1D block sparse CCS
     if(my_p_row_pexsi == 0) then
        call elsi_dense_to_ccs(matrix_aux, n_l_rows_pexsi, n_l_cols_pexsi, &
                               nnz_l_pexsi, H_real_pexsi, row_ind_pexsi, col_ptr_pexsi)
     endif

     matrix_aux = 0d0

     if(.not.overlap_is_unit) then
        ! Unpack overlap on the first process row in PEXSI process grid
        if(my_p_row_pexsi == 0) then
           do i_val = 1, nnz_l_pexsi
              ! Compute global 2d id
              global_col_id = FLOOR(1d0*(pos_recv_buffer(i_val)-1)/n_g_size)+1
              global_row_id = MOD(pos_recv_buffer(i_val), n_g_size)
              if(global_row_id == 0) global_row_id = n_g_size

              ! Compute local 2d id
              local_col_id = global_col_id-myid*FLOOR(1d0*n_g_size/n_p_per_pole_pexsi)
              local_row_id = global_row_id

              ! Put value to correct position
              matrix_aux(local_row_id,local_col_id) = s_val_recv_buffer(i_val)
           enddo
        endif

        deallocate(s_val_recv_buffer)

        if(my_p_row_pexsi == 0) then
           ! Transform overlap: 1D block dense ==> 1D block sparse CCS
           call elsi_dense_to_ccs_by_pattern(matrix_aux, n_l_rows_pexsi, n_l_cols_pexsi, &
                                             nnz_l_pexsi, row_ind_pexsi, col_ptr_pexsi, S_real_pexsi)
        endif
     endif

     deallocate(pos_recv_buffer)

     call elsi_stop_2dbc_to_1dccs_time()

  end subroutine ! elsi_2dbcd_to_1dbccs_hs_pexsi

!>
!! This routine converts density matrix computed by PEXSI and stored
!! in 1D block distributed sparse CCS format to 2D block-cyclic
!! distributed dense format.
!!
  subroutine elsi_1dbccs_to_2dbcd_dm_pexsi(D_out)

     implicit none
     include "mpif.h"

     real*8, intent(out) :: D_out(n_l_rows,n_l_cols) !< Density matrix to be converted

     integer :: i_row         !< Row counter
     integer :: i_col         !< Col counter
     integer :: i_val         !< Value counter
     integer :: j_val         !< Value counter
     integer :: i_proc        !< Process counter
     integer :: global_col_id !< Global column id
     integer :: global_row_id !< Global row id
     integer :: local_col_id  !< Local column id in 1D block distribution
     integer :: local_row_id  !< Local row id in 1D block distribution
     integer :: proc_col_id   !< Column id in process grid
     integer :: proc_row_id   !< Row id in process grid
     integer :: mpierr

     integer, allocatable :: dest(:)      !< Destination of each element
     integer, allocatable :: global_id(:) !< Global 1d id

     ! For the meaning of each array here, see documentation of MPI_Alltoallv
     real*8, allocatable :: val_send_buffer(:)  !< Send buffer for value
     integer, allocatable :: pos_send_buffer(:) !< Send buffer for global 1D id
     integer :: send_count(n_procs) !< Number of elements to send to each processor
     integer :: send_displ(n_procs) !< Displacement from which to take the outgoing data
     integer :: send_displ_aux      !< Auxiliary variable used to set displacement

     real*8, allocatable  :: val_recv_buffer(:) !< Receive buffer for value
     integer, allocatable :: pos_recv_buffer(:) !< Receive buffer for global 1D id
     integer :: recv_count(n_procs) !< Number of elements to receive from each processor
     integer :: recv_displ(n_procs) !< Displacement at which to place the incoming data
     integer :: recv_displ_aux      !< Auxiliary variable used to set displacement

     character*40, parameter :: caller = "elsi_1dbccs_to_2dbcd_dm_pexsi"

     call elsi_start_1dccs_to_2dbc_time()
     call elsi_statement_print("  Matrix conversion: 1D block CCS sparse ==> 2D block-cyclic dense")

     send_count = 0
     send_displ = 0
     recv_count = 0
     recv_displ = 0

     if(my_p_row_pexsi == 0) then
        call elsi_allocate(global_id, nnz_l_pexsi, "global_id", caller)
        call elsi_allocate(dest, nnz_l_pexsi, "dest", caller)
        call elsi_allocate(val_send_buffer, nnz_l_pexsi, "val_send_buffer", caller)
        call elsi_allocate(pos_send_buffer, nnz_l_pexsi, "pos_send_buffer", caller)

        i_col = 0
        ! Compute destination and global 1D id
        do i_val = 1, nnz_l_pexsi
           if(i_val == col_ptr_pexsi(i_col+1) .and. i_col /= n_l_cols_pexsi) then
              i_col = i_col+1
           endif
           i_row = row_ind_pexsi(i_val)

           ! Compute global id
           global_row_id = i_row
           global_col_id = i_col+myid*FLOOR(1d0*n_g_size/n_p_per_pole_pexsi)
           global_id(i_val) = (global_col_id-1)*n_g_size+global_row_id

           ! Compute destination
           proc_row_id = MOD(FLOOR(1d0*(global_row_id-1)/n_b_rows), n_p_rows)
           proc_col_id = MOD(FLOOR(1d0*(global_col_id-1)/n_b_cols), n_p_cols)
           dest(i_val) = proc_col_id+proc_row_id*n_p_cols
        enddo

        j_val = 0

        ! Set send_count
        do i_proc = 0, n_procs-1
           do i_val = 1, nnz_l_pexsi
              if(dest(i_val) == i_proc) then
                 j_val = j_val+1
                 val_send_buffer(j_val) = D_pexsi(i_val)
                 pos_send_buffer(j_val) = global_id(i_val)
                 send_count(i_proc+1) = send_count(i_proc+1)+1
              endif
           enddo
        enddo

        deallocate(global_id)
        deallocate(dest)
     endif

     ! Set recv_count
     call MPI_Alltoall(send_count, 1, mpi_integer, recv_count, &
                       1, mpi_integer, mpi_comm_global, mpierr)

     nnz_l = sum(recv_count,1)

     ! Set send and receive displacement
     send_displ_aux = 0
     recv_displ_aux = 0

     do i_proc = 0, n_procs-1
        send_displ(i_proc+1) = send_displ_aux
        send_displ_aux = send_displ_aux+send_count(i_proc+1)

        recv_displ(i_proc+1) = recv_displ_aux
        recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
     enddo

     call elsi_allocate(val_recv_buffer, nnz_l, "val_recv_buffer", caller)
     call elsi_allocate(pos_recv_buffer, nnz_l, "pos_recv_buffer", caller)

     ! Send and receive the packed data
     call MPI_Alltoallv(val_send_buffer, send_count, send_displ, mpi_real8, &
                        val_recv_buffer, recv_count, recv_displ, mpi_real8, &
                        mpi_comm_global, mpierr)

     call MPI_Alltoallv(pos_send_buffer, send_count, send_displ, mpi_integer, &
                        pos_recv_buffer, recv_count, recv_displ, mpi_integer, &
                        mpi_comm_global, mpierr)

     if(my_p_row_pexsi == 0) then
        deallocate(val_send_buffer)
        deallocate(pos_send_buffer)
     endif

     D_out = 0d0

     ! Unpack density matrix
     do i_val = 1, nnz_l
        ! Compute global 2d id
        global_col_id = FLOOR(1d0*(pos_recv_buffer(i_val)-1)/n_g_size)+1
        global_row_id = MOD(pos_recv_buffer(i_val), n_g_size)
        if(global_row_id == 0) global_row_id = n_g_size

        ! Compute local 2d id
        local_row_id = FLOOR(1d0*(global_row_id-1)/(n_p_rows*n_b_rows))*n_b_rows&
                       +MOD((global_row_id-1), n_b_rows)+1
        local_col_id = FLOOR(1d0*(global_col_id-1)/(n_p_cols*n_b_cols))*n_b_cols&
                       +MOD((global_col_id-1), n_b_cols)+1

        ! Put value to correct position
        D_out(local_row_id, local_col_id) = val_recv_buffer(i_val)
     enddo

     deallocate(val_recv_buffer)
     deallocate(pos_recv_buffer)

     call elsi_stop_1dccs_to_2dbc_time()

  end subroutine ! elsi_1dbccs_to_2dbcd_dm_pexsi

!>
!! This routine interfaces to PEXSI.
!!
  subroutine elsi_solve_evp_pexsi()

     implicit none
     include "mpif.h"

     real*8, save :: this_pexsi_tol = 1d-2

     character*40, parameter :: caller = "elsi_solve_evp_pexsi"

     call elsi_start_solve_evp_time()

     if(small_pexsi_tol) then
        pexsi_options%numElectronPEXSITolerance = this_pexsi_tol
        if(myid == 0) &
           write(*,"(A,E10.1)") "  | Current tolerance of number of electrons: ",&
                 this_pexsi_tol
     endif

     if(.not.allocated(D_pexsi)) then
        call elsi_allocate(D_pexsi,nnz_l_pexsi,"D_pexsi",caller)
     endif
     D_pexsi = 0d0

     if(.not.allocated(ED_pexsi)) then
        call elsi_allocate(ED_pexsi,nnz_l_pexsi,"ED_pexsi",caller)
     endif
     ED_pexsi = 0d0

     if(.not.allocated(FD_pexsi)) then
        call elsi_allocate(FD_pexsi,nnz_l_pexsi,"FD_pexsi",caller)
     endif
     FD_pexsi = 0d0

     ! Load sparse matrices for PEXSI
     if(overlap_is_unit) then
        call f_ppexsi_load_real_hs_matrix(pexsi_plan,pexsi_options,n_g_size,nnz_g,&
                                          nnz_l_pexsi,n_l_cols_pexsi,col_ptr_pexsi,&
                                          row_ind_pexsi,H_real_pexsi,1,S_real_pexsi,&
                                          pexsi_info)
     else
        call f_ppexsi_load_real_hs_matrix(pexsi_plan,pexsi_options,n_g_size,nnz_g,&
                                          nnz_l_pexsi,n_l_cols_pexsi,col_ptr_pexsi,&
                                          row_ind_pexsi,H_real_pexsi,0,S_real_pexsi,&
                                          pexsi_info)
     endif

     if(pexsi_info /= 0) &
        call elsi_stop(" PEXSI not able to load H/S matrix. Exiting...",caller)

     ! Inertia counting is only performed in the first few steps
     if(n_elsi_calls > n_inertia_steps) then
        pexsi_options%isInertiaCount = 0
     else
        pexsi_options%isInertiaCount = 1
     endif

     ! Solve the eigenvalue problem
     call elsi_statement_print("  Launch PEXSI DFT driver")

     call f_ppexsi_dft_driver(pexsi_plan,pexsi_options,n_electrons,mu_pexsi,&
                              n_electrons_pexsi,mu_min_inertia,mu_max_inertia,&
                              n_total_inertia_iter,n_total_pexsi_iter,pexsi_info)
         
     if(pexsi_info /= 0) &
        call elsi_stop(" PEXSI DFT Driver not able to solve problem. Exiting...",caller)

     if(small_pexsi_tol) then
        if(abs(n_electrons-n_electrons_pexsi) < this_pexsi_tol) then
           if(1d-1*this_pexsi_tol > final_pexsi_tol) then
              this_pexsi_tol = 1d-1*this_pexsi_tol
           else
              this_pexsi_tol = final_pexsi_tol
           endif
        endif
     endif

     ! Get the results
     call f_ppexsi_retrieve_real_dft_matrix(pexsi_plan,D_pexsi,ED_pexsi,FD_pexsi,&
                                            e_tot_H,e_tot_S,f_tot,pexsi_info)

     if(pexsi_info /= 0) then
        call elsi_stop(" PEXSI not able to retrieve solution. Exiting...",caller)
     endif

     call MPI_Barrier(mpi_comm_global,mpierr)
     call elsi_stop_solve_evp_time()

  end subroutine ! elsi_solve_evp_pexsi

!>
!! This routine overrides PEXSI default settings.
!!
  subroutine elsi_customize_pexsi(temperature_in,gap_in,delta_E_in,n_poles_in,&
                                  n_inertia_steps_in,max_iteration_in,mu_min_in,&
                                  mu_max_in,mu0_in,mu_inertia_tolerance_in,&
                                  mu_inertia_expansion_in,mu_safeguard_in,&
                                  n_electron_tolerance_in,matrix_type_in,&
                                  is_symbolic_factorize_in,ordering_in,&
                                  row_ordering_in,np_symbolic_factorize_in,&
                                  symmetric_in,transpose_in,verbosity_in)

     implicit none

     real(c_double), intent(in), optional :: temperature_in
     real(c_double), intent(in), optional :: gap_in
     real(c_double), intent(in), optional :: delta_E_in
     integer(c_int), intent(in), optional :: n_poles_in
     integer(c_int), intent(in), optional :: n_inertia_steps_in
     integer(c_int), intent(in), optional :: max_iteration_in
     real(c_double), intent(in), optional :: mu_min_in
     real(c_double), intent(in), optional :: mu_max_in
     real(c_double), intent(in), optional :: mu0_in
     real(c_double), intent(in), optional :: mu_inertia_tolerance_in
     real(c_double), intent(in), optional :: mu_inertia_expansion_in
     real(c_double), intent(in), optional :: mu_safeguard_in
     real(c_double), intent(in), optional :: n_electron_tolerance_in
     integer(c_int), intent(in), optional :: matrix_type_in
     integer(c_int), intent(in), optional :: is_symbolic_factorize_in
     integer(c_int), intent(in), optional :: ordering_in
     integer(c_int), intent(in), optional :: row_ordering_in
     integer(c_int), intent(in), optional :: np_symbolic_factorize_in
     integer(c_int), intent(in), optional :: symmetric_in
     integer(c_int), intent(in), optional :: transpose_in
     integer(c_int), intent(in), optional :: verbosity_in

     if(method == PEXSI) then
        ! Set default settings
        call elsi_set_pexsi_default_options()

        ! Temperature, in the same unit as H
        ! default: 0.0019 = 300K
        if(present(temperature_in)) &
           pexsi_options%temperature = temperature_in

        ! Spectral gap, can be set to 0 in most cases (default)
        if(present(gap_in)) &
           pexsi_options%gap = gap_in

        ! Upper bound for the spectral radius of S^(-1)H
        ! default: 10
        if(present(delta_E_in)) &
           pexsi_options%deltaE = delta_E_in

        ! Number of poles
        ! default: 40
        if(present(n_poles_in)) &
           pexsi_options%numPole = n_poles_in

        ! Number of steps to perform inertia counting
        ! default: 3
        if(present(n_inertia_steps_in)) &
           n_inertia_steps = n_inertia_steps_in

        ! Maximum number of PEXSI iterations after each inertia
        ! counting procedure
        ! default: 3
        if(present(max_iteration_in)) &
           pexsi_options%maxPEXSIIter = max_iteration_in

        ! From the second step, initial guess of mu is from previous step
        if(n_elsi_calls == 0) then
           ! Initial guess of mu
           ! default: 0.0
           if(present(mu0_in)) &
              pexsi_options%mu0 = mu0_in
        endif

        ! Initial guess of lower bound for mu
        ! default: -10.0
        if(present(mu_min_in)) &
           pexsi_options%muMin0 = mu_min_in

        ! Initial guess of upper bound for mu
        ! default: 10.0
        if(present(mu_max_in)) &
           pexsi_options%muMax0 = mu_max_in

        ! Stopping criterion in terms of the chemical potential
        ! for the inertia counting procedure
        ! default: 0.05
        if(present(mu_inertia_tolerance_in)) &
           pexsi_options%muInertiaTolerance = mu_inertia_tolerance_in

        ! If the chemical potential is not in the initial interval,
        ! the interval is expanded by this value
        ! default: 0.3
        if(present(mu_inertia_expansion_in)) &
           pexsi_options%muInertiaExpansion = mu_inertia_expansion_in

        ! Safeguard criterion in terms of the chemical potential to
        ! reinvoke the inertia counting procedure
        ! default: 0.05
        if(present(mu_safeguard_in)) &
           pexsi_options%muPEXSISafeGuard = mu_safeguard_in

        ! Stopping criterion of the PEXSI iteration in terms of the
        ! number of electrons compared to the exact number
        ! default: 0.01
        if(present(n_electron_tolerance_in)) then
           pexsi_options%numElectronPEXSITolerance = n_electron_tolerance_in
           if(n_electron_tolerance_in < 1d-2) then
              small_pexsi_tol = .true.
              final_pexsi_tol = n_electron_tolerance_in
           endif
        endif

        ! Type of input H and S matrices
        ! 0: real symmetric (default)
        ! 1: general complex
        if(present(matrix_type_in)) &
           pexsi_options%matrixType = matrix_type_in

        ! Whether to perform symbolic factorization
        ! default: 1
        if(present(is_symbolic_factorize_in)) &
           pexsi_options%isSymbolicFactorize = is_symbolic_factorize_in

        ! Ordering strategy for factorization and selected inversion
        ! 0: parallel ordering using ParMETIS
        ! 1: sequential ordering using METIS
        ! 2: multiple minimum degree ordering
        if(present(ordering_in)) &
           pexsi_options%ordering = ordering_in

        ! Row permutation strategy for factorization and selected inversion
        ! 0: no row permutation
        ! 1: make diagonal entry larger than off diagonal
        if(present(row_ordering_in)) &
           pexsi_options%rowOrdering = row_ordering_in

        ! Number of processors for ParMETIS, only used if ordering=0
        if(present(np_symbolic_factorize_in)) &
           pexsi_options%npSymbFact = np_symbolic_factorize_in

        ! Matrix structure
        ! 0: unsymmetric
        ! 1: symmetric (default)
        if(present(symmetric_in)) &
           pexsi_options%symmetric = symmetric_in

        ! Transpose
        ! 0: factor non transposed matrix (default)
        ! 1: factor transposed matrix
        if(present(transpose_in)) &
           pexsi_options%transpose = transpose_in

        ! Level of output information
        ! 0: no output
        ! 1: basic output (default)
        ! 2: detailed output
        if(present(verbosity_in)) &
           pexsi_options%verbosity = verbosity_in

        pexsi_customized = .true.
        call elsi_print_pexsi_options()
     else
        call elsi_statement_print("  The chosen method is not PEXSI."//&
                                  " Ignore elsi_customize_pexsi call.")
     endif

  end subroutine ! elsi_customize_pexsi

!======================
! ELSI solver routines
!======================

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
           call elsi_solve_evp_elpa()
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

  end subroutine ! elsi_ev_real

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
           call elsi_solve_evp_elpa()
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

  end subroutine ! elsi_ev_complex

!>
!! This routine computes density matrix.
!!
  subroutine elsi_dm_real(H_in,S_in,D_out,energy_out,broadening_in,width_in)

     implicit none

     real*8,  target, intent(in)   :: H_in(n_l_rows,n_l_cols)  !< Hamiltonian
     real*8,  target, intent(in)   :: S_in(n_l_rows,n_l_cols)  !< Overlap
     real*8,  intent(out)          :: D_out(n_l_rows,n_l_cols) !< Density matrix
     real*8,  intent(out)          :: energy_out               !< Energy
     integer, intent(in), optional :: broadening_in            !< (For ELPA) Occupation broadening scheme
     real*8,  intent(in), optional :: width_in                 !< (For ELPA) Occupation broadening width

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
           if(present(broadening_in)) then
              if(present(width_in)) then
                 broadening_type = broadening_in
                 broadening_width = width_in
              else
                 call elsi_stop(" ELPA is chosen to compute density matrix."//&
                                " Please specify broadening width! Exiting...",caller)
              endif
           else
              call elsi_stop(" ELPA is chosen to compute density matrix."//&
                             " Please specify broadening type! Exiting...",caller)
           endif

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
              call elsi_stop(" Unit overlap in OMM not yet implemented."//&
                             " Exiting...",caller)
           endif

           if(mod(nint(n_electrons),2) /= 0) then
              call elsi_stop(" The current implementation of OMM does not"//&
                             " work with fractional occupation numbers. This"//&
                             " means number of electrons in non-spin-polarized"//&
                             " system cannot be odd. Exiting...",caller)
           endif

           if(.not.omm_customized) then
              call elsi_set_omm_default_options()
              call elsi_print_omm_options()
           endif

           if(n_elsi_calls .le. n_elpa_steps) then ! Compute OMM initial guess by ELPA
              call elsi_set_method(ELPA)
              ! Broadening type and width are not critical here
              broadening_type = 0 ! Gaussian
              broadening_width = 0.01

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

              ! Switch back to OMM here to guarantee elsi_customize_omm
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
                 ! OMM coefficient matrix is the transpose of ELPA eigenvectors
                 call pdtran(n_g_size,n_g_size,1d0,D_elpa,1,1,sc_desc,0d0,C_real,1,1,sc_desc)

                 Coeff_omm%dval(1:Coeff_omm%iaux2(1),1:Coeff_omm%iaux2(2)) = &
                    C_real(1:Coeff_omm%iaux2(1),1:Coeff_omm%iaux2(2))

                 ! ELPA matrices are no longer needed at this point
                 if(associated(H_real))     nullify(H_real)
                 if(associated(S_real))     nullify(S_real)
                 if(allocated(C_real))      deallocate(C_real)
                 if(allocated(eigenvalues)) deallocate(eigenvalues)
                 if(allocated(D_elpa))      deallocate(D_elpa)
                 if(allocated(occ_elpa))    deallocate(occ_elpa)
              endif

              ! Continue the computation using OMM
              call elsi_solve_evp_omm()

              call elsi_get_dm(D_out)
              call elsi_get_energy(energy_out)
           endif

        case (PEXSI)
           if(.not.pexsi_customized) then
              call elsi_set_pexsi_default_options()
              call elsi_print_pexsi_options()
           endif

           ! PEXSI may use different process grid to achieve
           ! the efficient 2-level parallelization
           call elsi_init_pexsi()

           if(n_g_size < n_p_per_pole_pexsi) then
              call elsi_stop(" The (global) size of matrix is too small for"//&
                             " this number of processes. Exiting...",caller)
           endif

           ! Convert 2D block-cyclic dense Hamiltonian and overlap
           ! matrices to 1D block CCS sparse format
           call elsi_2dbcd_to_1dbccs_hs_pexsi(H_in,S_in)

           call elsi_solve_evp_pexsi()

           ! Convert 1D block CCS sparse density matrix to 2D
           ! block-cyclic dense format
           call elsi_get_dm(D_out)
           call elsi_get_energy(energy_out)

           if(allocated(H_real_pexsi))  deallocate(H_real_pexsi)
           if(allocated(S_real_pexsi))  deallocate(S_real_pexsi)
           if(allocated(D_pexsi))       deallocate(D_pexsi)
           if(allocated(ED_pexsi))      deallocate(ED_pexsi)
           if(allocated(FD_pexsi))      deallocate(FD_pexsi)
           if(allocated(row_ind_pexsi)) deallocate(row_ind_pexsi)
           if(allocated(col_ptr_pexsi)) deallocate(col_ptr_pexsi)

           call f_ppexsi_plan_finalize(pexsi_plan, pexsi_info)

        case (CHESS)
           call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
        case default
           call elsi_stop(" No supported method has been chosen."//&
                          " Exiting...",caller)
     end select

  end subroutine ! elsi_dm_real

!>
!! This routine computes density matrix.
!!
  subroutine elsi_dm_complex(H_in,S_in,D_out,energy_out,broadening_in,width_in)

     implicit none

     complex*16, target, intent(in) :: H_in(n_l_rows,n_l_cols)  !< Hamiltonian
     complex*16, target, intent(in) :: S_in(n_l_rows,n_l_cols)  !< Overlap
     real*8, intent(out)            :: D_out(n_l_rows,n_l_cols) !< Density matrix
     real*8, intent(out)            :: energy_out               !< Energy
     integer, intent(in), optional  :: broadening_in            !< (For ELPA) Occupation broadening scheme
     real*8, intent(in), optional   :: width_in                 !< (For ELPA) Occupation broadening width

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
           if(present(broadening_in)) then
              if(present(width_in)) then
                 broadening_type = broadening_in
                 broadening_width = width_in
              else
                 call elsi_stop(" ELPA is chosen to compute density matrix."//&
                                " Please specify broadening width! Exiting...",caller)
              endif
           else
              call elsi_stop(" ELPA is chosen to compute density matrix."//&
                             " Please specify broadening type! Exiting...",caller)
           endif

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
              call elsi_stop(" Unit overlap in OMM not yet implemented."//&
                             " Exiting...",caller)
           endif

           if(mod(nint(n_electrons),2) /= 0) then
              call elsi_stop(" The current implementation of OMM does not"//&
                             " work with fractional occupation numbers. This"//&
                             " means number of electrons in non-spin-polarized"//&
                             " system cannot be odd. Exiting...",caller)
           endif

           if(.not.omm_customized) then
              call elsi_set_omm_default_options()
              call elsi_print_omm_options()
           endif

           ! Set Hamiltonian and overlap matrices
           call elsi_set_hamiltonian(H_in)
           if(.not.overlap_is_unit) then
              call elsi_set_overlap(S_in)
           endif

           call elsi_solve_evp_omm()
           call elsi_get_dm(D_out)
           call elsi_get_energy(energy_out)

        case (PEXSI)
           if(n_g_size < n_p_per_pole_pexsi) then
              call elsi_stop(" The (global) size of matrix is too small for"//&
                             " this number of processes. Exiting...",caller)
           endif

           if(.not.pexsi_customized) then
              call elsi_set_pexsi_default_options()
              call elsi_print_pexsi_options()
           endif

           call elsi_stop(" PEXSI not yet implemented. Exiting...",caller)

        case (CHESS)
           call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
        case default
           call elsi_stop(" No supported method has been chosen."//&
                          " Exiting...",caller)
     end select

  end subroutine ! elsi_dm_complex

end module ELSI
