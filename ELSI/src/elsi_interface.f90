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

!============
! ELSI tools
!============

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
      n_states = nint(n_electrons/2d0)
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
!! This routine overrides ELSI default settings.
!!
subroutine elsi_customize(unit_overlap,hartree_to_ev,numerical_zero,mu_accuracy,&
                          no_check_singularity,singularity_threshold,&
                          force_stop_singularity,broadening_scheme,broadening_width)

   implicit none

   logical, intent(in), optional :: unit_overlap
   real*8,  intent(in), optional :: hartree_to_ev
   real*8,  intent(in), optional :: numerical_zero
   real*8,  intent(in), optional :: mu_accuracy
   logical, intent(in), optional :: no_check_singularity
   real*8,  intent(in), optional :: singularity_threshold
   logical, intent(in), optional :: force_stop_singularity
   integer, intent(in), optional :: broadening_scheme
   real*8,  intent(in), optional :: broadening_width

   ! Is the overlap matrix unit? [Default: .false.]
   if(present(unit_overlap)) &
      overlap_is_unit = unit_overlap
   ! User-defined value for Hartree [Default: 27.211386, Codata 2015]
   if(present(hartree_to_ev)) &
      hartree = hartree_to_ev
   ! Threshold to define numerical zero [Default: 1d-13]
   if(present(numerical_zero)) &
      zero_threshold = numerical_zero
   ! Accuracy for chemical potential determination [Default: 1d-10]
   if(present(mu_accuracy)) &
      occ_tolerance = mu_accuracy
   ! Disable checking for overlap singularity? [Default: .false.]
   if(present(no_check_singularity)) &
      no_singularity_check = no_check_singularity
   ! Eigenfunctions of overlap matrix with eigenvalues smaller than
   ! this value will be removed to avoid singularity [Default: 1d-5]
   if(present(singularity_threshold)) &
      singularity_tolerance = singularity_threshold
   ! Always stop if overlap is singular? [Default: .false.]
   if(present(force_stop_singularity)) &
      stop_singularity = force_stop_singularity
   ! Broadening scheme to compute Fermi level [Default: GAUSSIAN]
   if(present(broadening_scheme)) &
      broaden_method = broadening_scheme
   ! Broadening width to compute Fermi level [Default: 1d-2]
   if(present(broadening_width)) &
      broaden_width = broadening_width

end subroutine

!>
!! Set MPI.
!!
subroutine elsi_set_mpi(mpi_comm_global_in,n_procs_in,myid_in)

   implicit none

   integer, intent(in) :: myid_in            !< local process id
   integer, intent(in) :: n_procs_in         !< number of mpi processes
   integer, intent(in) :: mpi_comm_global_in !< global mpi communicator

   if(parallelism == 1) then ! MULTI_PROC
      mpi_comm_global = mpi_comm_global_in
      n_procs = n_procs_in
      myid = myid_in
      mpi_is_setup = .true.
   endif

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

      mpi_comm_row = mpi_comm_row_in
      mpi_comm_col = mpi_comm_col_in

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

end subroutine

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
         ! Note that here the convention of density matrix used in libOMM is
         ! changed to the one used in ELPA and PEXSI.
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

!==============
! ELSI solvers
!==============

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

         if(mod(nint(n_electrons),2) /= 0) then
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
               if(associated(H_real))     nullify(H_real)
               if(associated(S_real))     nullify(S_real)
               if(allocated(C_real))      deallocate(C_real)
               if(allocated(eigenvalues)) deallocate(eigenvalues)
               if(allocated(D_elpa))      deallocate(D_elpa)
               if(allocated(occ_elpa))    deallocate(occ_elpa)
            endif

            ! Continue the computation using libOMM
            call elsi_solve_evp_omm()

            call elsi_get_dm(D_out)
            call elsi_get_energy(energy_out)
         endif

      case (PEXSI)
         call elsi_print_pexsi_options()

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

         if(mod(nint(n_electrons),2) /= 0) then
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
         if(n_g_size < n_p_per_pole_pexsi) then
            call elsi_stop(" The (global) size of matrix is too small for"//&
                           " this number of processes. Exiting...",caller)
         endif

         call elsi_print_pexsi_options()

         call elsi_stop(" PEXSI not yet implemented. Exiting...",caller)

      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case default
         call elsi_stop(" No supported method has been chosen."//&
                        " Exiting...",caller)
   end select

end subroutine

end module ELSI
