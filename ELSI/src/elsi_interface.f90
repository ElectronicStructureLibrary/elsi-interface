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

!> ELSI Interchange
!! This module provides routines for setting up and solving or circumventing
!! an eigenvalue problem using ELPA, libOMM, PEXSI, or CheSS.
!!

module ELSI

   use iso_c_binding
   use ELSI_DIMENSIONS
   use ELSI_TIMERS
   use ELSI_MPI_TOOLS
   use ELSI_HDF5_TOOLS
   use ELSI_MATRIX_CONVERSION
   use ELPA1
   use ELPA2
   use MatrixSwitch
   use f_ppexsi_interface

   implicit none
   private

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
   !< Sparse Hamiltonian
   real*8,  allocatable :: H_real_sparse(:)
   !< Sparse overlap
   real*8,  allocatable :: S_real_sparse(:)
   !< Sparse density matrix
   real*8,  allocatable :: D_real_sparse(:)
   !< Sparse energy density matrix
   real*8,  allocatable :: ED_real_sparse(:)
   !< Sparse free energy density matrix
   real*8,  allocatable :: FD_real_sparse(:)
   !< Sparse index array
   integer, allocatable :: sparse_index(:)
   !< Sparse pointer array
   integer, allocatable :: sparse_pointer(:)

   !< From ELSI_DIMENSIONS
   public :: AUTO, ELPA, LIBOMM, PEXSI, CHESS
   public :: REAL_VALUES, COMPLEX_VALUES
   public :: BLOCK_CYCLIC

   ! From ELSI_MPI_TOOLS
   public :: elsi_init_mpi
   public :: elsi_init_blacs
   public :: elsi_set_mpi
   public :: elsi_set_blacs
   public :: elsi_get_global_row
   public :: elsi_get_global_col

   ! Public routines
   ! In use
   public :: elsi_init            !< Initialize
   public :: elsi_customize_elpa  !< Override ELPA default
   public :: elsi_customize_omm   !< Override OMM default
   public :: elsi_customize_pexsi !< Override PEXSI default
   public :: elsi_ev_real         !< Compute eigenvalues and eigenvectors
   public :: elsi_ev_complex      !< Compute eigenvalues and eigenvectors
   public :: elsi_dm_real         !< Compute density matrix
   public :: elsi_dm_complex      !< Compute density matrix
   public :: elsi_ev_real_ms      !< MatrixSwitch version
   public :: elsi_finalize        !< Clean memory and print timings
   ! Legacy
   public :: elsi_init_problem
   public :: elsi_init_problem_from_file
   public :: elsi_set_method
   public :: elsi_set_mode
   public :: elsi_allocate_matrices
   public :: elsi_deallocate_matrices
   public :: elsi_set_hamiltonian
   public :: elsi_set_overlap
   public :: elsi_solve_evp
   public :: elsi_set_hamiltonian_element
   public :: elsi_symmetrize_hamiltonian
   public :: elsi_get_hamiltonian
   public :: elsi_set_overlap_element
   public :: elsi_symmetrize_overlap
   public :: elsi_get_overlap
   public :: elsi_write_evp
   public :: elsi_read_evp
   public :: scalapack_dense_to_pexsi_sparse

   ! In use
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

   ! Legacy
   interface elsi_set_hamiltonian_element
      module procedure elsi_set_real_hamiltonian_element,&
                       elsi_set_complex_hamiltonian_element
   end interface

   interface elsi_get_hamiltonian
      module procedure elsi_get_real_hamiltonian,&
                       elsi_get_complex_hamiltonian
   end interface

   interface elsi_set_overlap_element
      module procedure elsi_set_real_overlap_element,&
                       elsi_set_complex_overlap_element
   end interface

   interface elsi_get_overlap
      module procedure elsi_get_real_overlap,&
                       elsi_get_complex_overlap
   end interface

contains

!========
! In use
!========

!>
!!  This routine initializes ELSI with global matrix size, 
!!  number of states, and method.
!!
subroutine elsi_init(solver, matrix_format, matrix_size, number_of_states)

   implicit none

   integer, intent(in) :: solver           !< AUTO,ELPA,LIBOMM,PEXSI,CHESS,...
   integer, intent(in) :: matrix_format    !< BLOCK_CYCLIC,...
   integer, intent(in) :: matrix_size      !< Global dimension of matrix
   integer, intent(in) :: number_of_states !< Number of states

   if(solver == 0) stop

   call elsi_set_method(solver)
   call elsi_set_storage(matrix_format)

   n_g_rank = matrix_size
   n_states = number_of_states

   call elsi_init_timers()
   call elsi_start_total_time()

end subroutine

!>
!!  This routine sets the method.
!!
subroutine elsi_set_method(i_method)

   implicit none

   integer, intent(in) :: i_method !< AUTO,ELPA,OMM,PEXSI,CHESS,...
   
   method = i_method

end subroutine

!>
!!  This routine sets the mode (real or complex).
!!
subroutine elsi_set_mode(i_mode)

   implicit none

   integer, intent(in) :: i_mode !< REAL_VALUES,COMPLEX_VALUES

   mode = i_mode

end subroutine 

!>
!!  This routine sets the matrix storage format.
!!
subroutine elsi_set_storage(i_storage)

   implicit none

   integer, intent(in) :: i_storage !< BLOCK_CYCLIC,...

   storage = i_storage

end subroutine

!>
!!  This routine prepares the matrices.
!!
subroutine elsi_allocate_matrices()

   implicit none

   character*40, parameter :: caller = "elsi_allocate_matrices"

   select case (method)
      case (ELPA)
         call elsi_allocate(eigenvalues,n_g_rank,"eigenvalues",caller)

         select case (mode)
            case (COMPLEX_VALUES)
               call elsi_allocate(C_complex, n_l_rows, n_l_cols,&
                                  "C_complex", caller)
            case (REAL_VALUES)
               call elsi_allocate(C_real, n_l_rows, n_l_cols,&
                                  "C_real", caller)
            case DEFAULT
               call elsi_stop(" No mode has been chosen. "//&
                              " Please choose REAL_VALUES or COMPLEX_VALUES. ",&
                              caller)
         end select

      case (LIBOMM)
         select case (mode)
            case (COMPLEX_VALUES)
               call m_allocate(D_omm, n_g_rank, n_g_rank, "pddbc")
               call m_allocate(Coeff_omm, n_states, n_g_rank, "pddbc")
            case (REAL_VALUES)
               call m_allocate(D_omm, n_g_rank, n_g_rank, "pddbc")
               call m_allocate(Coeff_omm, n_states, n_g_rank, "pddbc")
            case DEFAULT
               call elsi_stop(" No mode has been chosen. "//&
                              " Please choose REAL_VALUES or COMPLEX_VALUES. ",&
                              caller)
         end select

      case (PEXSI)
         call elsi_stop(" PEXSI not yet implemented. Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting... ",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select
end subroutine

!>
!!  This routine sets the real hamiltonian matrix.
!!
subroutine elsi_set_real_hamiltonian(H_in)

   implicit none

   real*8, target, intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian

   character*40, parameter :: caller = "elsi_set_real_hamiltonian"

   select case (method)
      case (ELPA)
         H_real => H_in
      case (LIBOMM)
         call m_register_pdbc(H_omm,H_in,sc_desc)
      case (PEXSI)
         call elsi_stop(" PEXSI not yet implemented. Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting... ",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select

end subroutine

!>
!!  This routine sets the complex hamiltonian matrix.
!!
subroutine elsi_set_complex_hamiltonian(H_in)

   implicit none

   complex*16, target, intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian

   character*40, parameter :: caller = "elsi_set_complex_hamiltonian"

   select case (method)
      case (ELPA)
         H_complex => H_in
      case (LIBOMM)
         call m_register_pdbc(H_omm,H_in,sc_desc)
      case (PEXSI)
         call elsi_stop(" PEXSI not yet implemented. Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting... ",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select

end subroutine

!>
!!  This routine sets the real overlap matrix.
!!
subroutine elsi_set_real_overlap(S_in)

   implicit none

   real*8, target, intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap

   character*40, parameter :: caller = "elsi_set_real_overlap"

   select case (method)
      case (ELPA)
         S_real => S_in
      case (LIBOMM)
         call m_register_pdbc(S_omm,S_in,sc_desc)
      case (PEXSI)
         call elsi_stop(" PEXSI not yet implemented. Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting... ",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select

   overlap_is_unity = .false.

end subroutine

!>
!!  This routine sets the complex overlap matrix.
!!
subroutine elsi_set_complex_overlap(S_in)

   implicit none

   complex*16, target, intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap

   character*40, parameter :: caller = "elsi_set_complex_overlap"

   select case (method)
      case (ELPA)
         S_complex => S_in
      case (LIBOMM)
         call m_register_pdbc(S_omm,S_in,sc_desc)
      case (PEXSI)
         call elsi_stop(" PEXSI not yet implemented. Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting... ",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select

   overlap_is_unity = .false.

end subroutine

!>
!!  This routine gets the eigenvalues.
!!
subroutine elsi_get_eigenvalues(e_val_out)

   implicit none

   real*8, intent(out) :: e_val_out(n_states) !< Eigenvalues

   character*40, parameter :: caller = "elsi_get_eigenvalues"

   select case (method)
      case (ELPA)
         e_val_out(1:n_states) = eigenvalues(1:n_states)
      case (LIBOMM)
         call elsi_stop(" OMM does not compute eigenvalues! Exiting...",caller)
      case (PEXSI)
         call elsi_stop(" PEXSI does not compute eigenvalues! Exiting...",caller)
      case (CHESS)
         call elsi_stop(" CHESS does not compute eigenvalues! Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select

end subroutine

!>
!!  This routine gets the eigenvectors.
!!
subroutine elsi_get_real_eigenvectors(e_vec_out)

   implicit none

   real*8, intent(out) :: e_vec_out(n_l_rows,n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_get_real_eigenvectors"

   select case (method)
      case (ELPA)
         e_vec_out = C_real
      case (LIBOMM)
         call elsi_stop(" OMM does not compute eigenvectors! Exiting...",caller)
      case (PEXSI)
         call elsi_stop(" PEXSI does not compute eigenvectors! Exiting...",caller)
      case (CHESS)
         call elsi_stop(" CHESS does not compute eigenvectors! Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select

end subroutine

!>
!!  This routine gets the eigenvectors.
!!
subroutine elsi_get_complex_eigenvectors(e_vec_out)

   implicit none

   complex*16, intent(out) :: e_vec_out(n_l_rows,n_l_cols) !< Eigenvectors

   character*40, parameter :: caller = "elsi_get_complex_eigenvectors"

   select case (method)
      case (ELPA)
         e_vec_out = C_complex
      case (LIBOMM)
         call elsi_stop(" OMM does not compute eigenvectors! Exiting...",caller)
      case (PEXSI)
         call elsi_stop(" PEXSI does not compute eigenvectors! Exiting...",caller)
      case (CHESS)
         call elsi_stop(" CHESS does not compute eigenvectors! Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select

end subroutine

!>
!!  This routine gets the number of occupied states from occupation numbers.
!!
subroutine elsi_get_occupied_number(occ)

   implicit none

   real*8, intent(in) :: occ(n_states)

   integer :: occupied, i

   do i = 1,n_states
      if(occ(i) > 0d0) then
         occupied = i
      endif
   enddo

   n_states = occupied

end subroutine

!>
!!  This routine gets the energy.
!!
subroutine elsi_get_energy(energy_out)

   implicit none

   real*8, intent(out) :: energy_out

   character*40, parameter :: caller = "elsi_get_energy"

   select case (method)
      case (ELPA)
         energy_out = SUM(eigenvalues(1:n_states))
      case (LIBOMM)
         energy_out = total_energy
      case (PEXSI)
         call elsi_stop(" PEXSI: not yet implemented! Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" CHESS: not yet implemented! Exiting... ",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select

end subroutine

!>
!!  This routine gets the density matrix.
!!
subroutine elsi_get_dm(D_out)

   implicit none

   real*8, intent(out) :: D_out(n_l_rows,n_l_cols) !< Density matrix

   character*40, parameter :: caller = "elsi_get_dm"

   select case (method)
      case (ELPA)
         call elsi_stop(" ELPA needs to compute density matrix from eigenvectors. "//&
                        " Exiting... ",caller)
      case (LIBOMM)
         D_out = D_omm%dval
      case (PEXSI)
         call elsi_stop(" PEXSI: not yet implemented! Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" CHESS: not yet implemented! Exiting... ",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select

end subroutine

!>
!!  This routine constructs the density matrix using eigenvectors from ELPA.
!!
subroutine elsi_compute_dm_elpa(D_out,occ)

   implicit none

   real*8, intent(out) :: D_out(n_l_rows,n_l_cols) !< Density matrix
   real*8, intent(in) :: occ(n_states) !< Occupation number

   real*8, allocatable :: tmp_real(:,:) !< Real eigenvectors, temporary
   complex*16, allocatable :: tmp_complex(:,:) !< Complex eigenvectors, temporary
   real*8, allocatable :: factor(:) !< Factor to construct density matrix

   integer, allocatable :: local_col(:)
   integer :: i_col, i
   integer :: l_row, l_col ! Local index
   integer :: g_row, g_col !< Global index


   character*40, parameter :: caller = "elsi_compute_dm"

   select case (method)
      case (ELPA)
         ! Map global columns to local
         call elsi_allocate(local_col,n_g_rank,"local_col",caller)

         i_col = 0 ! local column counter

         do i = 1,n_g_rank
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
                  if(occ(i) > 0d0) then
                     factor(i) = sqrt(occ(i))
                  endif
               enddo

               do i = 1,n_states
                  if(factor(i) > 0d0) then
                     if(local_col(i) > 0) then
                        tmp_real(:,local_col(i)) = tmp_real(:,local_col(i)) * factor(i)
                     endif
                  elseif(local_col(i) .ne. 0) then
                     tmp_real(:,local_col(i)) = 0d0
                  endif
               enddo

               ! Compute density matrix
               D_out = 0d0

               call pdsyrk('U','N',n_g_rank,n_states,1.d0,&
                           tmp_real,1,1,sc_desc,0.d0,D_out,1,1,sc_desc)

            case (COMPLEX_VALUES)
               ! Get eigenvectors into tmp_complex
               call elsi_allocate(tmp_complex,n_l_rows,n_l_cols,"tmp_complex",caller)
               call elsi_get_eigenvectors(tmp_complex)

               ! Compute the factors used to construct density matrix
               call elsi_allocate(factor,n_states,"factor",caller)
               factor = 0d0

               do i = 1,n_states
                  if(occ(i) > 0d0) then
                     factor(i) = sqrt(occ(i))
                  endif
               enddo

               do i = 1,n_states
                  if(factor(i) > 0d0) then
                     if(local_col(i) > 0) then
                        tmp_complex(:,local_col(i)) = tmp_complex(:,local_col(i)) * factor(i)
                     endif
                  elseif(local_col(i) .ne. 0) then
                     tmp_complex(:,local_col(i)) = 0d0
                  endif
               enddo

               ! Compute density matrix
               D_out = 0d0

               call pzherk('U','N',n_g_rank,n_states,(1.d0,0.d0), &
                           tmp_complex,1,1,sc_desc,(0.d0,0.d0),D_out,1,1,sc_desc)

         end select

         deallocate(local_col)
         deallocate(factor)
         if(allocated(tmp_real))    deallocate(tmp_real)
         if(allocated(tmp_complex)) deallocate(tmp_complex)

         ! Now D_out is an upper triangle matrix
         ! Set D_out to full
         call elsi_allocate(tmp_real,n_l_rows,n_l_cols,"tmp_real",caller)
         tmp_real = D_out

         call pdtran(n_g_rank,n_g_rank,1.d0,tmp_real,1,1,sc_desc,1.d0,D_out,1,1,sc_desc)

         deallocate(tmp_real)

         do l_row = 1,n_l_rows
            call elsi_get_global_row(g_row,l_row)
            do l_col = l_row,n_l_cols
               call elsi_get_global_col(g_col,l_col)
               if(g_row == g_col) then
                  D_out(l_row,l_col) = 0.5d0 * D_out(l_row,l_col)
               endif
            enddo
         enddo

      case (LIBOMM)
         call elsi_stop(" LIBOMM does not compute density matrix from eigenvectors! "//&
                        " Exiting... ",caller)
      case (PEXSI)
         call elsi_stop(" PEXSI does not compute density matrix from eigenvectors! "//&
                        " Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" CHESS does not compute density matrix from eigenvectors! "//&
                        " Exiting... ",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select

end subroutine

!>
!!  This routine deallocates the matrices.
!!
subroutine elsi_deallocate_matrices()

   implicit none

   ! Nullify pointers
   if(associated(H_real))        nullify(H_real)
   if(associated(H_complex))     nullify(H_complex)
   if(associated(S_real))        nullify(S_real)
   if(associated(S_complex))     nullify(S_complex)
 
   ! Free Memory
   ! ELPA
   if(allocated(C_real))         deallocate(C_real)
   if(allocated(C_complex))      deallocate(C_complex)
   if(allocated(eigenvalues))    deallocate(eigenvalues)
   if(allocated(D_elpa))         deallocate(D_elpa)
   ! PEXSI
   if(allocated(H_real_sparse))  deallocate(H_real_sparse)
   if(allocated(S_real_sparse))  deallocate(S_real_sparse)
   if(allocated(D_real_sparse))  deallocate(D_real_sparse)
   if(allocated(ED_real_sparse)) deallocate(ED_real_sparse)
   if(allocated(FD_real_sparse)) deallocate(FD_real_sparse)
   if(allocated(sparse_index))   deallocate(sparse_index)
   if(allocated(sparse_pointer)) deallocate(sparse_pointer)
   ! OMM
   if(H_omm%is_initialized)      call m_deallocate(H_omm)
   if(S_omm%is_initialized)      call m_deallocate(S_omm)
   if(D_omm%is_initialized)      call m_deallocate(D_omm)
   ! Coefficient matrix will only be deallocated when finalizing ELSI
!   if(Coeff_omm%is_initialized)  call m_deallocate(Coeff_omm)
   if(T_omm%is_initialized)      call m_deallocate(T_omm)

end subroutine

!>
!!  This routine finalizes ELSI.
!!
subroutine elsi_finalize()

   implicit none
   include "mpif.h"

   call MPI_BARRIER(mpi_comm_global, mpierr)

   call elsi_deallocate_matrices()

   if(Coeff_omm%is_initialized) call m_deallocate(Coeff_omm)
   if(method == PEXSI) call f_ppexsi_plan_finalize(pexsi_plan, pexsi_info)
   
   call elsi_stop_total_time()
   call elsi_print_timers()

   if(.not.external_blacs) call elsi_finalize_blacs()
   if(.not.external_mpi)   call elsi_finalize_mpi()

end subroutine

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
subroutine elsi_to_standard_evp(cholesky)

   logical, intent(in) :: cholesky                 !< If .True. factorize Overlap

   logical :: success                              !< Success flag
   real*8,     allocatable :: buffer_real (:,:)    !< Real valued matrix buffer
   complex*16, allocatable :: buffer_complex (:,:) !< Complex valued matrix buffer

   character*40, parameter :: caller = "elsi_to_standard_evp"

   select case (method)
      case (ELPA)
         select case (mode)
            case (COMPLEX_VALUES)
               call elsi_allocate(buffer_complex, n_l_rows, n_l_cols, &
                                  "buffer_complex",caller)

               if(cholesky) then
                  call elsi_statement_print(" Starting Cholesty decomposition")
                  ! Compute S = (U^T)U, U -> S
                  call cholesky_complex(n_g_rank, S_complex, n_l_rows, &
                                        n_b_rows, n_l_cols, mpi_comm_row, &
                                        mpi_comm_col, .False., success)
                  if(.not.success) then
                     call elsi_stop("Cholesky decomposition failed.",caller)
                  endif

                  ! compute U^-1 -> S
                  call invert_trm_complex(n_g_rank, S_complex, n_l_rows, &
                                          n_b_rows, n_l_cols, mpi_comm_row, &
                                          mpi_comm_col, .False., success)
                  if(.not.success) then
                     call elsi_stop("Matrix invertion failed.",caller)
                  endif
               endif

               ! compute H(U^-1) -> buff
               call pzgemm('N','N', n_g_rank, n_g_rank, n_g_rank, 1.0d0, &
                           H_complex, 1, 1, sc_desc, S_complex, 1, 1, sc_desc, &
                           0.0d0, buffer_complex, 1, 1, sc_desc)

               ! compute ((U^-1)^T)H by (H(U^-1))^T -> H
               call pztranc(n_g_rank, n_g_rank, 1.d0, buffer_complex, 1, 1, &
                            sc_desc, 0.d0, H_complex, 1, 1, sc_desc)

               ! compute ((U^-1)^T)H(U^-1) -> H
               buffer_complex = H_complex
               call pzgemm('N','N', n_g_rank, n_g_rank, n_g_rank, 1.0d0, &
                           buffer_complex, 1, 1, sc_desc, S_complex, 1, 1, &
                           sc_desc, 0.0d0, H_complex, 1, 1, sc_desc)

            case (REAL_VALUES)
               call elsi_allocate(buffer_real, n_l_rows, n_l_cols, &
                                  "buffer_real",caller)

               if(cholesky) then
                  call elsi_statement_print(" Starting Cholesty decomposition")
                  ! Compute S = (U^T)U, U -> S
                  call cholesky_real(n_g_rank, S_real, n_l_rows, &
                                     n_b_rows, n_l_cols, mpi_comm_row, &
                                     mpi_comm_col, .False., success)
                  if(.not.success) then
                     call elsi_stop("Cholesky decomposition failed.",caller)
                  endif

                  ! compute U^-1 -> S
                  call invert_trm_real(n_g_rank, S_real, n_l_rows, &
                                       n_b_rows, n_l_cols, mpi_comm_row, &
                                       mpi_comm_col, .False., success)
                  if(.not.success) then
                     call elsi_stop("Matrix invertion failed.",caller)
                  endif
               endif

               ! compute H(U^-1) -> buff
               call pdgemm('N','N', n_g_rank, n_g_rank, n_g_rank, 1.0d0, &
                           H_real, 1, 1, sc_desc, S_real, 1, 1, sc_desc, &
                           0.0d0, buffer_real, 1, 1, sc_desc)

               ! compute ((U^-1)^T)H by (H(U^-1))^T -> H
               call pdtran(n_g_rank, n_g_rank, 1.d0, buffer_real, 1, 1, &
                           sc_desc, 0.d0, H_real, 1, 1, sc_desc)

               ! compute ((U^-1)^T)H(U^-1) -> H
               buffer_real = H_real
               call pdgemm('N','N', n_g_rank, n_g_rank, n_g_rank, 1.0d0, &
                           buffer_real, 1, 1, sc_desc, S_real, 1, 1, &
                           sc_desc, 0.0d0, H_real, 1, 1, sc_desc)
         end select

      case (LIBOMM)
         call elsi_stop(" OMM does not need to transform evp. Exiting...",caller)
      case (PEXSI)
         call elsi_stop(" PEXSI does not need to transform evp. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting... ",caller)
   end select

   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

end subroutine

!> 
!! This routine does the back-transformation of the eigenvectors in standard
!! form (A'c' = c'v) to the original generalized form (Ac = Bcv)
!!
!! v = (U^-1)v'
!!
subroutine elsi_to_original_ev()

   real*8, allocatable :: buffer_real(:,:)         !< Real valued matrix buffer
   complex*16, allocatable :: buffer_complex (:,:) !< Complex valued matrix buffer

   character*40, parameter :: caller = "elsi_to_original_ev"

   select case (method)
      case (ELPA)
         select case (mode)
            case (COMPLEX_VALUES)
               ! (U^-1) is stored in S_complex after
               ! elsi_to_standard_evp
               ! vectors_complex = S_complex * vectors_complex
               call elsi_allocate(buffer_complex,n_l_rows,n_l_cols,"temp",caller)
               buffer_complex = C_complex

               call pzgemm('N', 'N', n_g_rank, n_states, n_g_rank, &
                           1.0d0, S_complex, 1, 1, sc_desc, buffer_complex, 1, 1, &
                           sc_desc, 0.0d0, C_complex, 1, 1, sc_desc)

            case (REAL_VALUES)
               ! (U^-1) is stored in S_real after
               ! elsi_to_standard_evp
               ! C_real = S_real * C_real
               call elsi_allocate(buffer_real,n_l_rows,n_l_cols,"temp",caller)
               buffer_real = C_real

               ! method (a)
!               call pdtran(n_g_rank, n_g_rank, 1.d0, S_real, 1, 1, sc_desc, &
!                           0.d0, H_real, 1, 1, sc_desc)
!               call mult_at_b_real('L', 'N', n_g_rank, n_states, H_real, &
!                                   n_l_rows, buffer_real, n_l_rows, n_b_rows, &
!                                   mpi_comm_row, mpi_comm_col, C_real, n_l_rows)

               ! method (b)
               call pdgemm('N', 'N', n_g_rank, n_states, n_g_rank, &
                           1.0d0, S_real, 1, 1, sc_desc, buffer_real , 1, 1, &
                           sc_desc, 0.0d0, C_real, 1, 1, sc_desc)
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

end subroutine

!>
!!  This routine interfaces to ELPA.
!!
subroutine elsi_solve_evp_elpa(cholesky)

   implicit none
   include "mpif.h"

   logical, intent(in) :: cholesky !< Cholesky factorize overlap?

   logical :: success
   logical :: two_step_solver

   character*40, parameter :: caller = "elsi_solve_evp_elpa"

   call elsi_start_solve_evp_time()

   if(mode == REAL_VALUES .and. .not. associated(H_real)) then
      call elsi_stop(" Hamiltonian not found. Exiting...",caller)
   endif

   if(mode == COMPLEX_VALUES .and. .not. associated(H_complex)) then
      call elsi_stop(" Hamiltonian not found. Exiting...",caller)
   endif

   ! Choose 1-stage or 2-stage solver
   if(elpa_one_always) then
      two_step_solver = .false.
   elseif(elpa_two_always) then
      two_step_solver = .true.
   elseif(n_g_rank < 256) then
      two_step_solver = .false.
   else
      two_step_solver = .true.
   endif

   ! Transform to standard form
   if(.not. overlap_is_unity) then
      call elsi_statement_print(" Tansforming to standard evp")
      call elsi_to_standard_evp(cholesky)
   endif

   ! Solve evp, return eigenvalues and eigenvectors
   if(two_step_solver) then ! 2-stage solver
      call elsi_statement_print(" Starting ELPA 2-stage solver")
      select case (mode)
         case (COMPLEX_VALUES)
            success = solve_evp_complex_2stage(n_g_rank, n_states, &
                         H_complex, n_l_rows, eigenvalues, &
                         C_complex, n_l_rows, n_b_rows, n_l_cols, &
                         mpi_comm_row, mpi_comm_col, mpi_comm_global)
         case (REAL_VALUES)
            success = solve_evp_real_2stage(n_g_rank, n_states, H_real, &
                         n_l_rows, eigenvalues, C_real, n_l_rows, &
                         n_b_rows, n_l_cols, mpi_comm_row, mpi_comm_col, &
                         mpi_comm_global)
      end select
   else ! 1-stage solver
      call elsi_statement_print(" Starting ELPA 1-stage solver")
      select case (mode)
         case (COMPLEX_VALUES)
            success = solve_evp_complex(n_g_rank, n_states, H_complex, &
                         n_l_rows, eigenvalues, C_complex, n_l_rows, &
                         n_b_rows, n_l_cols, mpi_comm_row, mpi_comm_col)
         case (REAL_VALUES)
            success = solve_evp_real(n_g_rank, n_states, H_real, &
                         n_l_rows, eigenvalues, C_real, n_l_rows, &
                         n_b_rows, n_l_cols, mpi_comm_row, mpi_comm_col)
      end select
   endif

   if(.not. success) then
      call elsi_stop(" ELPA failed when solving eigenvalue problem. "//&
                     " Exiting...",caller)
   endif

   ! Back-transform eigenvectors
   if(.not. overlap_is_unity) then
      call elsi_statement_print(" Transforming to original eigenvectors")
      call elsi_to_original_ev()
   endif

   call MPI_BARRIER(mpi_comm_global,mpierr)
   call elsi_stop_solve_evp_time()

end subroutine

!>
!!  This routine overrides ELPA default settings.
!!
subroutine elsi_customize_elpa(elpa_solver)

   implicit none

   integer, intent(in) :: elpa_solver

   character*40, parameter :: caller = "elsi_customize_elpa"

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
      call elsi_stop(" The chosen method is not ELPA. Exiting... ", caller)
   endif

end subroutine

!>
!!  This routine interfaces to libOMM.
!!
subroutine elsi_solve_evp_omm(cholesky)

   implicit none
   include "mpif.h"

   logical, intent(in) :: cholesky !< Cholesky factorize overlap?

   logical :: success
   character*40, parameter :: caller = "elsi_solve_evp_omm"

   call elsi_start_solve_evp_time()

   if(.not.omm_customized) then
      call elsi_set_omm_default_options()
      call elsi_print_omm_options()
   endif

   if(cholesky) then
      new_overlap = .true.
      C_matrix_initialized = .false.

      ! Cholesky
      select case (mode)
         case (COMPLEX_VALUES)
            ! Compute S = (U^T)U, U -> S
            call cholesky_complex(n_g_rank, S_omm%zval, n_l_rows, n_b_rows, n_l_cols, &
                                  mpi_comm_row, mpi_comm_col, .false., success)
         case (REAL_VALUES)
            ! Compute S = (U^T)U, U -> S
            call cholesky_real(n_g_rank, S_omm%dval, n_l_rows, n_b_rows, n_l_cols, &
                               mpi_comm_row, mpi_comm_col, .false., success)
      end select
   else
      new_overlap = .false.
      C_matrix_initialized = .true.
   endif

   ! Shift eigenvalue spectrum
!   call m_add(S_elsi,'N',H_elsi,-eta,1d0,"lap")

   select case (mode)
      case (COMPLEX_VALUES)
         call omm(n_g_rank, n_states, H_omm, S_omm, new_overlap, total_energy, &
                  D_omm, calc_ED, eta, Coeff_omm, C_matrix_initialized, T_omm, &
                  scale_kinetic, omm_flavour, nk_times_nspin, i_k_spin, min_tol, &
                  omm_verbose, do_dealloc, "pzdbc", "lap", myid)
      case (REAL_VALUES)
         call omm(n_g_rank, n_states, H_omm, S_omm, new_overlap, total_energy, &
                  D_omm, calc_ED, eta, Coeff_omm, C_matrix_initialized, T_omm, &
                  scale_kinetic, omm_flavour, nk_times_nspin, i_k_spin, min_tol, &
                  omm_verbose, do_dealloc, "pddbc", "lap", myid)
   end select

   call MPI_BARRIER(mpi_comm_global,mpierr)
   call elsi_stop_solve_evp_time()

end subroutine

!>
!!  This routine overrides libOMM default settings.
!!
subroutine elsi_customize_omm(scale_kinetic_in,calc_ed_in,eta_in,&
                              min_tol_in,nk_times_nspin_in,i_k_spin_in,&
                              omm_verbose_in,do_dealloc_in)

   implicit none

   real*8,  intent(in), optional :: scale_kinetic_in
   logical, intent(in), optional :: calc_ed_in
   real*8,  intent(in), optional :: eta_in
   real*8,  intent(in), optional :: min_tol_in
   integer, intent(in), optional :: nk_times_nspin_in
   integer, intent(in), optional :: i_k_spin_in
   logical, intent(in), optional :: omm_verbose_in
   logical, intent(in), optional :: do_dealloc_in

   character*40, parameter :: caller = "elsi_customize_omm"

   if(method == LIBOMM) then
      ! Set default settings
      call elsi_set_omm_default_options()

      ! Scaling of kinetic energy matrix
      if(present(scale_kinetic_in)) scale_kinetic = scale_kinetic_in
      ! Calculate energy weigthed density matrix?
      if(present(calc_ed_in)) calc_ed = calc_ed_in
      ! Eigenspectrum shift parameter
      if(present(eta_in)) eta = eta_in
      ! Tolerance for minimization
      if(present(min_tol_in)) min_tol = min_tol_in
      ! n_k_points * n_spin
      if(present(nk_times_nspin_in)) nk_times_nspin = nk_times_nspin_in
      ! Index
      if(present(i_k_spin_in)) i_k_spin = i_k_spin_in
      ! OMM output?
      if(present(omm_verbose_in)) omm_verbose = omm_verbose_in
      ! Deallocate internal storage?
      if(present(do_dealloc_in)) do_dealloc = do_dealloc_in

      omm_customized = .true.
      call elsi_print_omm_options()
   else
      call elsi_stop(" The chosen method is not OMM. Exiting... ", caller)
   endif

end subroutine

!>
!!  This routine interfaces to PEXSI.
!!
subroutine elsi_solve_evp_pexsi()

   implicit none
   include "mpif.h"

   character*40, parameter :: caller = "elsi_solve_evp_pexsi"

   call elsi_start_solve_evp_time()

   if(.not.pexsi_customized) then
      call elsi_set_pexsi_default_options()
      call elsi_print_pexsi_options()
   endif

   call MPI_BARRIER(mpi_comm_global,mpierr)
   call elsi_stop_solve_evp_time()

end subroutine

!>
!!  This routine overrides PEXSI default settings.
!!
subroutine elsi_customize_pexsi(temperature_in, gap_in, delta_E_in, n_poles_in, &
                                is_inertia_count_in, max_iteration_in, mu_min0_in, &
                                mu_max0_in, mu0_in, mu_inertia_tolerance_in, &
                                mu_inertia_expansion_in, mu_safeguard_in, &
                                n_electron_tolerance_in, matrix_type_in, &
                                is_symbolic_factorize_in, ordering_in, &
                                row_ordering_in, np_symbolic_factorize_in, &
                                symmetric_in, transpose_in, verbosity_in)

   implicit none

   real(c_double), intent(in), optional :: temperature_in
   real(c_double), intent(in), optional :: gap_in
   real(c_double), intent(in), optional :: delta_E_in
   integer(c_int), intent(in), optional :: n_poles_in
   integer(c_int), intent(in), optional :: is_inertia_count_in
   integer(c_int), intent(in), optional :: max_iteration_in
   real(c_double), intent(in), optional :: mu_min0_in
   real(c_double), intent(in), optional :: mu_max0_in
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

   character*40, parameter :: caller = "elsi_customize_pexsi"

   if(method == PEXSI) then
      ! Set default settings
      call elsi_set_pexsi_default_options()

      !TODO: comments here
      if(present(temperature_in)) &
         pexsi_options%temperature = temperature_in
      if(present(gap_in)) &
         pexsi_options%gap = gap_in
      if(present(delta_E_in)) &
         pexsi_options%deltaE = delta_E_in
      if(present(n_poles_in)) &
         pexsi_options%numPole = n_poles_in
      if(present(is_inertia_count_in)) &
         pexsi_options%isInertiaCount = is_inertia_count_in
      if(present(max_iteration_in)) &
         pexsi_options%maxPEXSIIter = max_iteration_in
      if(present(mu_min0_in)) &
         pexsi_options%muMin0 = mu_min0_in
      if(present(mu_max0_in)) &
         pexsi_options%muMax0 = mu_max0_in
      if(present(mu0_in)) &
         pexsi_options%mu0 = mu0_in
      if(present(mu_inertia_tolerance_in)) &
         pexsi_options%muInertiaTolerance = mu_inertia_tolerance_in
      if(present(mu_inertia_expansion_in)) &
         pexsi_options%muInertiaExpansion = mu_inertia_expansion_in
      if(present(mu_safeguard_in)) &
         pexsi_options%muPEXSISafeGuard = mu_safeguard_in
      if(present(n_electron_tolerance_in)) &
         pexsi_options%numElectronPEXSITolerance = n_electron_tolerance_in
      if(present(matrix_type_in)) &
         pexsi_options%matrixType = matrix_type_in
      if(present(is_symbolic_factorize_in)) &
         pexsi_options%isSymbolicFactorize = is_symbolic_factorize_in
      if(present(ordering_in)) &
         pexsi_options%ordering = ordering_in
      if(present(row_ordering_in)) &
         pexsi_options%rowOrdering = row_ordering_in
      if(present(np_symbolic_factorize_in)) &
         pexsi_options%npSymbFact = np_symbolic_factorize_in
      if(present(symmetric_in)) &
         pexsi_options%symmetric = symmetric_in
      if(present(transpose_in)) &
         pexsi_options%transpose = transpose_in
      if(present(verbosity_in)) &
         pexsi_options%verbosity = verbosity_in

      pexsi_customized = .true.
      call elsi_print_pexsi_options()
   else
      call elsi_stop(" The chosen method is not PEXSI. Exiting... ", caller)
   endif

end subroutine

!>
!!  This routine computes eigenvalues and eigenvectors.
!!
subroutine elsi_ev_real(H_in, S_in, e_val_out, e_vec_out, need_cholesky)

   implicit none

   real*8, target, intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian
   real*8, target, intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap
   real*8, intent(out) :: e_val_out(n_states)            !< Eigenvalues
   real*8, intent(out) :: e_vec_out(n_l_rows,n_l_cols)   !< Eigenvectors
   logical, intent(inout) :: need_cholesky               !< Cholesky factorize overlap?

   character*40, parameter :: caller = "elsi_ev_real"

   ! Here the only supported method is ELPA
   select case (method)
      case (ELPA)
         ! REAL case
         call elsi_set_mode(REAL_VALUES)

         ! Allocate matrices
         call elsi_allocate_matrices()

         ! Set Hamiltonian and Overlap matrices
         call elsi_set_hamiltonian(H_in)
         call elsi_set_overlap(S_in)

         ! Solve eigenvalue problem
         call elsi_solve_evp_elpa(need_cholesky)
         call elsi_get_eigenvalues(e_val_out)
         call elsi_get_eigenvectors(e_vec_out)

         ! Deallocate
         call elsi_deallocate_matrices()

      case (LIBOMM)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors. "//&
                        " Choose ELPA if necessary. Exiting... ",caller)
      case (PEXSI)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors. "//&
                        " Choose ELPA if necessary. Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors. "//&
                        " Choose ELPA if necessary. Exiting... ",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA to compute eigenpairs. "//&
                        " Exiting... ",caller)
   end select

   need_cholesky = .false.

end subroutine

subroutine elsi_ev_real_ms(H_ms, S_ms, e_val_out, C_ms, need_cholesky)
   implicit none

   type(matrix) :: H_ms                       !< Hamiltonian
   type(matrix) :: S_ms                       !< Overlap
   real*8, intent(out) :: e_val_out(n_states) !< Eigenvalues
   type(matrix) :: C_ms                       !< Eigenvectors
   logical, intent(inout) :: need_cholesky    !< Cholesky factorize overlap?

   character*40, parameter :: caller = "elsi_ev_real"

   ! Here the only supported method is ELPA
   select case (method)
      case (ELPA)
         ! REAL case
         call elsi_set_mode(REAL_VALUES)

         ! Allocate matrices
         call elsi_allocate_matrices()

         ! Set Hamiltonian and Overlap matrices
         call elsi_set_hamiltonian(H_ms%dval)
         call elsi_set_overlap(S_ms%dval)

         ! Solve eigenvalue problem
         call elsi_solve_evp_elpa(need_cholesky)
         call elsi_get_eigenvalues(e_val_out)
         call elsi_get_eigenvectors(C_ms%dval)

         ! Deallocate
         call elsi_deallocate_matrices()

      case (LIBOMM)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors. "//&
                        " Choose ELPA if necessary. Exiting... ",caller)
      case (PEXSI)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors. "//&
                        " Choose ELPA if necessary. Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors. "//&
                        " Choose ELPA if necessary. Exiting... ",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA to compute eigenpairs. "//&
                        " Exiting... ",caller)
   end select

   need_cholesky = .false.

end subroutine

!>
!!  This routine computes eigenvalues and eigenvectors.
!!
subroutine elsi_ev_complex(H_in, S_in, e_val_out, e_vec_out, need_cholesky)

   implicit none

   complex*16, target, intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian
   complex*16, target, intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap
   real*8, intent(out) :: e_val_out(n_states)                !< Eigenvalues
   complex*16, intent(out) :: e_vec_out(n_l_rows,n_l_cols)   !< Eigenvectors
   logical, intent(inout) :: need_cholesky                   !< Cholesky factorize Overlap?

   character*40, parameter :: caller = "elsi_ev_complex"

   ! Here the only supported method is ELPA
   select case (method)
      case (ELPA)
         ! COMPLEX case
         call elsi_set_mode(COMPLEX_VALUES)

         ! Allocate matrices
         call elsi_allocate_matrices()

         ! Set Hamiltonian and Overlap matrices
         call elsi_set_hamiltonian(H_in)
         call elsi_set_overlap(S_in)

         ! Solve eigenvalue problem
         call elsi_solve_evp(need_cholesky)
         call elsi_get_eigenvalues(e_val_out)
         call elsi_get_eigenvectors(e_vec_out)

         ! Deallocate
         call elsi_deallocate_matrices()

      case (LIBOMM)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors. "//&
                        " Choose ELPA if necessary. Exiting... ",caller)
      case (PEXSI)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors. "//&
                        " Choose ELPA if necessary. Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" Only ELPA computes eigenvalues and eigenvectors. "//&
                        " Choose ELPA if necessary. Exiting... ",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA to compute eigenpairs. "//&
                        " Exiting... ",caller)
   end select

   need_cholesky = .false.

end subroutine

!>
!!  This routine computes density matrix.
!!
subroutine elsi_dm_real(H_in, S_in, D_out, energy_out, need_cholesky, occupation)

   implicit none

   real*8, target, intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian
   real*8, target, intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap
   real*8, intent(out) :: D_out(n_l_rows,n_l_cols)       !< Density matrix
   real*8, intent(out) :: energy_out                     !< Energy
   logical, intent(inout) :: need_cholesky               !< Cholesky factorize overlap?
   real*8,  intent(in), optional :: occupation(n_states) !< Occupation number

   character*40, parameter :: caller = "elsi_dm_real"

   ! REAL case
   call elsi_set_mode(REAL_VALUES)

   ! Allocation of matrices
   call elsi_allocate_matrices()

   ! Set Hamiltonian and overlap matrices
   call elsi_set_hamiltonian(H_in)
   call elsi_set_overlap(S_in)

   ! Solve eigenvalue problem
   select case (method)
      case (ELPA)
         call elsi_get_occupied_number(occupation)
         call elsi_solve_evp_elpa(need_cholesky)
         call elsi_allocate(D_elpa, n_l_rows, n_l_cols, "D_elpa", caller)
         call elsi_compute_dm_elpa(D_out, occupation)
         call elsi_get_energy(energy_out)
      case (LIBOMM)
         call elsi_solve_evp_omm(need_cholesky)
         call elsi_get_dm(D_out)
         call elsi_get_energy(energy_out)
      case (PEXSI)
         call elsi_stop(" PEXSI not yet implemented. Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting... ",caller)
      case default
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA to compute eigenpairs. "//&
                        " Exiting... ",caller)
   end select

   ! Deallocate
   call elsi_deallocate_matrices()

   need_cholesky = .false.

end subroutine

!>
!!  This routine computes density matrix.
!!
subroutine elsi_dm_complex(H_in, S_in, D_out, energy_out, need_cholesky, occupation)

   implicit none

   complex*16, target, intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian
   complex*16, target, intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap
   real*8, intent(out) :: D_out(n_l_rows,n_l_cols)           !< Density matrix
   real*8, intent(out) :: energy_out                         !< Energy
   logical, intent(inout) :: need_cholesky                   !< Cholesky factorize overlap?
   real*8,  intent(in), optional :: occupation(n_states)     !< Occupation number

   character*40, parameter :: caller = "elsi_dm_complex"

   ! COMPLEX case
   call elsi_set_mode(COMPLEX_VALUES)

   ! Allocation of matrices
   call elsi_allocate_matrices()

   ! Set Hamiltonian and Overlap matrices
   call elsi_set_hamiltonian(H_in)
   call elsi_set_overlap(S_in)

   ! Solve eigenvalue problem
   select case (method)
      case (ELPA)
         call elsi_get_occupied_number(occupation)
         call elsi_solve_evp_elpa(need_cholesky)
         call elsi_allocate(D_elpa, n_l_rows, n_l_cols, "D_elpa", caller)
         call elsi_compute_dm_elpa(D_out, occupation)
         call elsi_get_energy(energy_out)
      case (LIBOMM)
         call elsi_solve_evp_omm(need_cholesky)
         call elsi_get_dm(D_out)
         call elsi_get_energy(energy_out)
      case (PEXSI)
         call elsi_stop(" PEXSI not yet implemented. Exiting... ",caller)
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting... ",caller)
      case default
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA to compute eigenpairs. "//&
                        " Exiting... ",caller)
   end select

   ! Deallocate
   call elsi_deallocate_matrices()

   need_cholesky = .false.

end subroutine

!==============================================
! Legecy: the following part is no longer used.
!==============================================

!>
!!  This routine sets the matrix dimensions.
!!
subroutine elsi_init_problem(matrixsize, block_rows, block_cols)

   implicit none

   integer, intent(in) :: matrixsize  !< Global dimension of matrix
   integer, intent(in) :: block_rows  !< Block rows of matrix
   integer, intent(in) :: block_cols  !< Block cols of matrix

   n_g_rank = matrixsize
   n_b_rows = block_rows
   n_b_cols = block_cols

end subroutine

!>
!!  This routine sets the element (i_row,i_col) in real hamiltonian matrix.
!!
subroutine elsi_set_real_hamiltonian_element(element,i_row,i_col)

   implicit none

   integer, intent(in) :: i_row   !< Row position
   integer, intent(in) :: i_col   !< Col position
   real*8, intent(in)  :: element !< Value 

   character*40, parameter :: caller = "elsi_set_real_hamiltonian_element"
   
   select case (method)
      case (ELPA,LIBOMM)
         if(mode == REAL_VALUES) then
            if(associated(H_real)) then
               H_real(i_row,i_col) = element
            else 
               call elsi_stop("Hamiltonian not created/linked.",caller)
            endif
         else  
            call elsi_stop("Wrong mode: "//&
                 "Complex valued hamiltonian to be written in real storage",&
                 caller)
         endif
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",caller)
      case DEFAULT
         call elsi_stop("No method has been chosen. "//&
                        "Please choose method ELPA, LIBOMM, or PEXSI",&
                        caller)
   end select

end subroutine

!>
!!  This routine sets the element (i_row, i_col) in the complex
!!  hamiltonian matrix.
!!
subroutine elsi_set_complex_hamiltonian_element(element,i_row,i_col)

   implicit none

   integer, intent(in)     :: i_row   !< Row position
   integer, intent(in)     :: i_col   !< Col position
   complex*16, intent(in)  :: element !< Value

   character*40, parameter :: caller = "elsi_set_complex_hamiltonian_element"

   select case (method)
      case (ELPA,LIBOMM)
         if(mode == COMPLEX_VALUES) then
            H_complex(i_row,i_col) = element      
         else  
            write(*,'(2a)') "Wrong mode:",&
                  "Complex valued hamiltonian to be written in real storage"
            stop
         endif
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",caller)
      case DEFAULT
         call elsi_stop("No method has been chosen. "//&
                        "Please choose method ELPA, LIBOMM, or PEXSI",&
                        caller)
   end select

end subroutine

!>
!!  This routine sets the element (i_row, i_col) in the real
!!  overlap matrix.
!!
subroutine elsi_set_real_overlap_element(element,i_row,i_col)

   implicit none

   integer, intent(in) :: i_row   !< Row position
   integer, intent(in) :: i_col   !< Col position
   real*8, intent(in)  :: element !< Value

   character*40, parameter :: caller = "elsi_set_real_overlap_element"

   select case (method)
      case (ELPA,LIBOMM)
         if(mode == REAL_VALUES) then
            S_real(i_row,i_col) = element      
         else  
            write(*,'(2a)') "Wrong mode:",&
                  "Complex valued overlap to be written in real storage"
            stop
         endif
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",caller)
      case DEFAULT
         call elsi_stop("No method has been chosen. "//&
                        "Please choose method ELPA, LIBOMM, or PEXSI",&
                        caller)
   end select

   overlap_is_unity = .False.

end subroutine

!>
!!  This routine symmetrizes an upper or lower triangle hamiltonian.
!!
subroutine elsi_symmetrize_hamiltonian()

  implicit none

  real*8,     allocatable :: buffer_real   (:,:)
  complex*16, allocatable :: buffer_complex(:,:)
  complex*16, parameter   :: CONE = (1d0,0d0)

  integer :: l_row, l_col  !< Local matrix indices
  integer :: g_row, g_col  !< Global matrix indices

  character*40, parameter :: caller = "elsi_symmetrize_hamiltonian"

  select case (method)
     case (ELPA,LIBOMM)
        if(mode == REAL_VALUES) then
           call elsi_allocate(buffer_real, n_l_rows, n_l_cols,&
                              "buffer_real", caller)
           buffer_real(:,:) = H_real(:,:)
           call pdtran(n_g_rank, n_g_rank, 1.d0, buffer_real, 1, 1,&
                       sc_desc, 1.d0, H_real, 1, 1, sc_desc)
           deallocate(buffer_real)
            
           do l_row = 1, n_l_rows
              call elsi_get_global_row(g_row, l_row)
              do l_col = l_row, n_l_cols
                 call elsi_get_global_col(g_col, l_col)
                 if(g_row == g_col) then
                    H_real(l_row,l_col) = 0.5d0 * H_real(l_row,l_col)
                 endif
              enddo
           enddo

        else 
           call elsi_allocate(buffer_complex, n_l_rows, n_l_cols,&
                              "buffer_complex", caller)
           buffer_complex(:,:) = H_complex(:,:)
           call pztranc(n_g_rank, n_g_rank, CONE, buffer_complex, 1, 1,&
                        sc_desc, CONE, H_complex, 1, 1, sc_desc)
           deallocate(buffer_complex)

           do l_row = 1, n_l_rows
              call elsi_get_global_row(g_row, l_row)
              do l_col = l_row, n_l_cols
                 call elsi_get_global_col(g_col, l_col)
                 if(g_row == g_col) then
                    H_complex(l_row,l_col) = 0.5d0 * H_complex(l_row,l_col)
                 endif
              enddo
           enddo

        endif
     case (PEXSI)
        call elsi_stop("PEXSI not implemented yet!",caller)
     case DEFAULT
        call elsi_stop("No method has been chosen. "//&
                       "Please choose method ELPA, LIBOMM, or PEXSI",&
                       caller)
   end select

   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

end subroutine 

!>
!!  This routine symmetrizes an upper or lower triangle overlap.
!!
subroutine elsi_symmetrize_overlap()

  implicit none

  real*8,     allocatable :: buffer_real   (:,:)
  complex*16, allocatable :: buffer_complex(:,:)
  complex*16, parameter   :: CONE = (1d0,0d0)

  integer :: l_row, l_col  !< Local matrix indices
  integer :: g_row, g_col  !< Global matrix indices

  character*40, parameter :: caller = "elsi_symmetrize_overlap"
  
  select case (method)
     case (ELPA,LIBOMM)
        if(mode == REAL_VALUES) then
           call elsi_allocate(buffer_real, n_l_rows, n_l_cols,&
                              "buffer_real", caller)
           buffer_real(:,:) = S_real(:,:)
           call pdtran(n_g_rank, n_g_rank, 1.d0, buffer_real, 1, 1,&
                       sc_desc, 1.d0, S_real, 1, 1, sc_desc)
           deallocate(buffer_real)

           do l_row = 1, n_l_rows
              call elsi_get_global_row(g_row, l_row)
              do l_col = l_row, n_l_cols
                 call elsi_get_global_col(g_col, l_col)
                 if(g_row == g_col) then
                    S_real(l_row,l_col) = 0.5d0 * S_real(l_row,l_col)
                 endif
              enddo
           enddo

        else
           call elsi_allocate(buffer_complex, n_l_rows, n_l_cols,&
                              "buffer_complex", caller)
           buffer_complex(:,:) = S_complex(:,:)
           call pztranc(n_g_rank, n_g_rank, CONE, buffer_complex, 1, 1,&
                        sc_desc, CONE, S_complex, 1, 1, sc_desc)
           deallocate(buffer_complex)

           do l_row = 1, n_l_rows
              call elsi_get_global_row(g_row, l_row)
              do l_col = l_row, n_l_cols
                 call elsi_get_global_col(g_col, l_col)
                 if(g_row == g_col) then
                    S_complex(l_row,l_col) = 0.5d0 * S_complex(l_row,l_col)
                 endif
              enddo
           enddo
        endif
     case (PEXSI)
        call elsi_stop("PEXSI not implemented yet!",caller)
     case DEFAULT
        call elsi_stop("No method has been chosen. "//&
                       "Please choose method ELPA, LIBOMM, or PEXSI",&
                       caller)
   end select
   
   if(allocated(buffer_real))    deallocate(buffer_real)
   if(allocated(buffer_complex)) deallocate(buffer_complex)

end subroutine 

!>
!!  This routine sets one element (i_row, i_col) of the complex overlap matrix.
!!
subroutine elsi_set_complex_overlap_element(element,i_row,i_col)

   implicit none

   integer, intent(in)     :: i_row   !< Row position
   integer, intent(in)     :: i_col   !< Col position
   complex*16, intent(in)  :: element !< Value

   character*40, parameter :: caller = "elsi_set_complex_overlap_element"

   select case (method)
      case (ELPA,LIBOMM)
         if(mode == COMPLEX_VALUES) then
            S_complex(i_row,i_col) = element      
         else  
            write(*,'(2a)') "Wrong mode:",&
                  "Real valued overlap to be written in complex storage"
            stop
         endif
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",caller)
      case DEFAULT
         call elsi_stop("No method has been chosen. "//&
                        "Please choose method ELPA, LIBOMM, or PEXSI",&
                        caller)
   end select

   overlap_is_unity = .False.

end subroutine

!>
!!  This routine gets the real hamiltonian matrix.
!!
subroutine elsi_get_real_hamiltonian(h,n_rows,n_cols)

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   real*8, intent(out) :: h(n_rows,n_cols) !< Hamiltonian

   character*40, parameter :: caller = "elsi_get_real_hamiltonian"

   select case (method)
      case (ELPA,LIBOMM)
         if(mode == REAL_VALUES) then
            h(:,:) = H_real(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:",&
                  "Real valued hamiltonian to be written in complex storage"
            stop
         endif
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",caller)
      case DEFAULT
         call elsi_stop("No method has been chosen. "//&
                        "Please choose method ELPA, LIBOMM, or PEXSI",&
                        caller)
   end select

end subroutine

!>
!!  This routine gets the complex hamiltonian matrix.
!!
subroutine elsi_get_complex_hamiltonian(h,n_rows,n_cols)

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   complex*16, intent(out) :: h(n_rows,n_cols) !< Hamiltonian

   character*40, parameter :: caller = "elsi_get_coomplex_hamiltonian"

   select case (method)
      case (ELPA,LIBOMM)
         if(mode == COMPLEX_VALUES) then
            h(:,:) = H_complex(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:",&
                  "Complex valued hamiltonian to be written in real storage"
            stop
         endif
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",caller)
      case DEFAULT
         call elsi_stop("No method has been chosen. "//&
                        "Please choose method ELPA, LIBOMM, or PEXSI",&
                        caller)
   end select

end subroutine

!>
!!  This routine gets the real overlap matrix.
!!
subroutine elsi_get_real_overlap(s,n_rows,n_cols)

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   real*8, intent(out) :: s(n_rows,n_cols) !< Overlap

   character*40, parameter :: caller = "elsi_get_real_overlap"

   select case (method)
      case (ELPA,LIBOMM)
         if(mode == REAL_VALUES) then
            s(:,:) = S_real(:,:)    
         else  
            write(*,'(2a)') "Wrong mode:",&
                  "Real valued overlap to be written in complex storage"
            stop
         endif
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",caller)
      case DEFAULT
         call elsi_stop("No method has been chosen. "//&
                        "Please choose method ELPA, LIBOMM, or PEXSI",&
                        caller)
   end select

end subroutine

!>
!!  This routine gets the complex overlap matrix.
!!
subroutine elsi_get_complex_overlap(s,n_rows,n_cols)

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   complex*16, intent(out) :: s(n_rows,n_cols) !< Overlap

   character*40, parameter :: caller = "elsi_get_complex_overlap"

   select case (method)
      case (ELPA,LIBOMM)
         if(mode == COMPLEX_VALUES) then
            s(:,:) = S_complex(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:",&
                  "Complex valued overlap to be written in real storage"
            stop
         endif
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",caller)
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ",&
               "Please choose method ELPA, LIBOMM, or PEXSI"
         stop
   end select

end subroutine

!>
!!  This routine writes the complete eigenvalue problem into a file.
!!
subroutine elsi_write_evp(file_name)

   implicit none
   include "mpif.h"

   character(len=*), intent(in) :: file_name !< File name where to write

   integer :: file_id   !< HDF5 File identifier
   integer :: group_id  !< HDF5 Group identifier
   real*8, allocatable :: buffer(:,:) !< Read buffer for PEXSI

   character*40, parameter :: caller = "elsi_write_evp"
  
   call elsi_start_write_evp_time()

   if(method == PEXSI) then
      call elsi_allocate(buffer, n_l_rows, n_l_cols, "buffer", caller) 
   endif

   call hdf5_create_file(file_name, mpi_comm_global, mpi_info_null, file_id)

   ! The Hamiltonian
   call hdf5_create_group(file_id, "hamiltonian", group_id)
   
   ! Matrix dimension
   call hdf5_write_attribute(group_id, "n_matrix_rows", n_g_rank)
   call hdf5_write_attribute(group_id, "n_matrix_cols", n_g_rank)

   ! Hamiltonian Write
   call hdf5_get_scalapack_pattern()
   
   select case (method)
      case (ELPA,LIBOMM)
         call hdf5_write_matrix_parallel(group_id, "matrix", H_real)
      case (PEXSI)
         call elsi_ccs_to_dense(buffer, n_l_rows, n_l_cols, H_real_sparse,&
                                n_l_nonzero, sparse_index, sparse_pointer)
         call hdf5_write_matrix_parallel(group_id, "matrix", buffer)
      case DEFAULT
         call elsi_stop("No supported method has been chosen. "&
                        // "Please choose method ELPA, LIBOMM, or PEXSI",&
                        caller)
   end select

   call hdf5_close_group(group_id)

   ! The Overlap Matrix
   call hdf5_create_group(file_id, "overlap", group_id)
   
   ! Matrix dimension
   call hdf5_write_attribute(group_id, "n_matrix_rows", n_g_rank)
   call hdf5_write_attribute(group_id, "n_matrix_cols", n_g_rank)

   ! Overlap Write
   select case (method)
      case (ELPA,LIBOMM) 
         call hdf5_write_matrix_parallel(group_id, "matrix", S_real)
      case (PEXSI)
         call elsi_ccs_to_dense(buffer, n_l_rows, n_l_cols, S_real_sparse,&
                                n_l_nonzero, sparse_index, sparse_pointer)
         call hdf5_write_matrix_parallel(group_id, "matrix", buffer)
      case DEFAULT
         call elsi_stop("No supported method has been chosen. "&
                        // "Please choose method ELPA, LIBOMM, or PEXSI",&
                        caller)
   end select

   call hdf5_close_group(group_id)

   call hdf5_close_file(file_id)

   call elsi_stop_write_evp_time()

   if(allocated(buffer)) deallocate(buffer)

end subroutine

!>
!!  This routine reads the eigenvalue problem from a file.
!!
subroutine elsi_read_evp(file_name)

   implicit none
   include "mpif.h"
   
   character(len=*), intent(in) :: file_name !< File to open

   integer :: file_id  !< HDF5 File identifier
   integer :: group_id !< HDF5 Group identifier
   real*8, allocatable :: buffer(:,:) !< Read buffer for PEXSI
  
   character*40, parameter :: caller = "elsi_read_evp"

   ! For Pexsi we need to create a buffer 
   ! we convert it directly to the CCS format

   call elsi_start_read_evp_time()

   if(method == PEXSI) then
      call elsi_allocate(buffer, n_l_rows, n_l_cols, "buffer", caller) 
   endif

   call hdf5_open_file(file_name, mpi_comm_global, mpi_info_null, file_id)

   ! The Hamiltonian
   call hdf5_open_group(file_id, "hamiltonian", group_id)
   
   ! Hamiltonian Read
   call hdf5_get_scalapack_pattern()
   select case (method)
      case (ELPA,LIBOMM)
         call hdf5_read_matrix_parallel(group_id, "matrix", H_real)
      case (PEXSI)
         call hdf5_read_matrix_parallel(group_id, "matrix", buffer)
         call elsi_compute_N_nonzero(buffer,n_l_rows, n_l_cols)
         call elsi_allocate(H_real_sparse, n_l_nonzero, "H_real_sparse", caller)
         call elsi_allocate(sparse_index, n_l_nonzero, "sparse_index", caller)
         call elsi_allocate(sparse_pointer, n_l_cols+1, "sparse_pointer", caller)
         call elsi_dense_to_ccs(buffer, n_l_rows, n_l_cols, H_real_sparse,&
                                n_l_nonzero, sparse_index, sparse_pointer)
      case DEFAULT
         call elsi_stop("No supported method has been chosen. "&
                        // "Please choose method ELPA, LIBOMM, or PEXSI",&
                        caller)
   end select

   call hdf5_close_group(group_id)

   ! The Overlap Matrix
   call hdf5_open_group(file_id, "overlap", group_id)
   
   ! Overlap Read
   select case (method)
      case (ELPA,LIBOMM)
         call hdf5_read_matrix_parallel(group_id, "matrix", S_real)
      case (PEXSI)
         call hdf5_read_matrix_parallel(group_id, "matrix", buffer)
         call elsi_allocate(S_real_sparse, n_l_nonzero, "S_real_sparse", caller)
         call elsi_dense_to_ccs_by_pattern(buffer, n_l_rows, n_l_cols, S_real_sparse,&
                                           n_l_nonzero, sparse_index, sparse_pointer)
      case DEFAULT
         call elsi_stop("No supported method has been chosen. "&
                        // "Please choose method ELPA, LIBOMM, or PEXSI",&
                        caller)
         stop
   end select

   ! TODO Check if overlap is unity
   overlap_is_unity = .False.
   
   call hdf5_close_group(group_id)

   call hdf5_close_file(file_id)

   if(allocated(buffer)) deallocate(buffer)

   call elsi_stop_read_evp_time()

end subroutine

!>
!!  This routine sets the method of choice for solving the eigenvalue problem.
!!
subroutine elsi_init_problem_from_file(file_name, block_rows, block_cols)

   implicit none
   include "mpif.h"

   character(len=*), intent(in) :: file_name !< File to open
   integer, intent(in) :: block_rows  !< Block rows of matrix
   integer, intent(in) :: block_cols  !< Block cols of matrix

   integer :: file_id  !< HDF5 File identifier 
   integer :: group_id !< HDF5 Group identifier

   character*40, parameter :: caller = "elsi_init_problem_from_file"

   call elsi_start_read_evp_time()

   n_b_rows = block_rows
   n_b_cols = block_cols

   if(.not. mpi_is_setup) call elsi_stop("MPI needs to be setup first!",caller)

   call hdf5_open_file(file_name, mpi_comm_global, mpi_info_null, file_id)

   ! The Hamiltonian
   call hdf5_open_group(file_id, "hamiltonian", group_id)

   ! Matrix dimension
   call hdf5_read_attribute(group_id, "n_matrix_rows", n_g_rank)
   call hdf5_read_attribute(group_id, "n_matrix_cols", n_g_rank)

   call hdf5_close_group(group_id)
   call hdf5_close_file(file_id)
   call elsi_stop_read_evp_time()

end subroutine

!>
!!  This routine interfaces to the eigenvalue solvers.
!!
subroutine elsi_solve_evp(cholesky)

   implicit none
   include "mpif.h"

   logical, intent(in) :: cholesky !< If .True. factorize Overlap

   logical :: success
   logical :: two_step_solver
   real*8  :: val
   integer :: i,j
   integer :: i_task

   character*100, parameter :: caller = "elsi_solve_evp"

   call elsi_start_solve_evp_time()

   if(method == ELPA .or. method == LIBOMM) then
      if(mode == REAL_VALUES .and. .not. associated(H_real)) then
         call elsi_stop("Hamiltonian not created/linked.",caller)
      endif
   
      if(mode == COMPLEX_VALUES .and. .not. associated(H_complex)) then
         call elsi_stop("Hamiltonian not created/linked.",caller)
      endif
   endif

   if(method == PEXSI) then
      if(mode == REAL_VALUES .and. .not. allocated(H_real_sparse)) then
         call elsi_stop("Hamiltonian not created/linked.",caller)
      endif
   endif

   select case (method)
      case (ELPA)
         ! Choose 1-stage or 2-stage solver
         two_step_solver = .True.

         if(n_g_rank < 256) then
            two_step_solver = .False.
         endif

         ! Transform to standard form
         if(.not. overlap_is_unity) then
            call elsi_statement_print(" Tansforming to standard evp")
            call elsi_to_standard_evp(cholesky)
         endif

         ! Solve evp, return eigenvalues and eigenvectors
         if(two_step_solver) then ! 2-stage solver
            call elsi_statement_print(" Starting ELPA 2-stage solver")
            select case (mode)
               case (COMPLEX_VALUES)
                  success = solve_evp_complex_2stage( &
                            n_g_rank, n_states, H_complex, &
                            n_l_rows, eigenvalues, C_complex, &
                            n_l_rows, n_b_rows, n_l_cols, &
                            mpi_comm_row, mpi_comm_col, mpi_comm_global)
               case (REAL_VALUES)
                  success = solve_evp_real_2stage( &
                            n_g_rank, n_states, H_real, &
                            n_l_rows, eigenvalues, C_real, &
                            n_l_rows, n_b_rows, n_l_cols, &
                            mpi_comm_row, mpi_comm_col, mpi_comm_global)
            end select
         else ! 1-stage solver
            call elsi_statement_print(" Starting ELPA 1-stage solver")
            select case (mode)
               case (COMPLEX_VALUES)
                  success = solve_evp_complex( &
                            n_g_rank, n_states, H_complex, &
                            n_l_rows, eigenvalues, C_complex, &
                            n_l_rows, n_b_rows, n_l_cols, &
                            mpi_comm_row, mpi_comm_col)
               case (REAL_VALUES)
                  success = solve_evp_real( &
                            n_g_rank, n_states, H_real, &
                            n_l_rows, eigenvalues, C_real, &
                            n_l_rows, n_b_rows, n_l_cols, &
                            mpi_comm_row, mpi_comm_col)
            end select
         endif

         if(.not.success) then
            call elsi_stop("ELPA failed in solving eigenvalue problem.",caller)
         endif

         ! Back-transform eigenvectors
         if(.not. overlap_is_unity) then
            call elsi_statement_print(" Transforming to original eigenvectors")
            call elsi_to_original_ev()
         endif

      case (LIBOMM)

         call elsi_set_omm_default_options()

         if(cholesky) then
            new_overlap = .True.
            C_matrix_initialized = .False.
         else
            new_overlap = .False.
            C_matrix_initialized = .True.
         endif

         ! Shift eigenvalue spectrum
!         call m_add(OMM_S_matrix,'N',OMM_H_matrix,-eta,1d0,"lap")

         call omm(n_g_rank, n_states, H_omm, S_omm, new_overlap, total_energy, &
                  D_omm, calc_ED, eta, Coeff_omm, C_matrix_initialized, T_omm, &
                  scale_kinetic, omm_flavour, nk_times_nspin, i_k_spin, &
                  min_tol, omm_verbose, do_dealloc, "pddbc", "lap", myid)

      case (PEXSI)

         if(.not. allocated(D_real_sparse))  allocate(D_real_sparse(n_l_nonzero))
         if(.not. allocated(ED_real_sparse)) allocate(ED_real_sparse(n_l_nonzero))
         if(.not. allocated(FD_real_sparse)) allocate(FD_real_sparse(n_l_nonzero))

         ! Set the default options
         ! TODO: User interface is missing
         call elsi_set_pexsi_default_options()

         ! Load sparse matrices for PEXSI
         if(overlap_is_unity) then
            call f_ppexsi_load_real_symmetric_hs_matrix(pexsi_plan, pexsi_options, &
                 n_g_rank, n_g_nonzero, n_l_nonzero, n_l_cols, sparse_pointer, &
                 sparse_index, H_real_sparse, 1, S_real_sparse, pexsi_info)
         else
            call f_ppexsi_load_real_symmetric_hs_matrix(pexsi_plan, pexsi_options, &
                 n_g_rank, n_g_nonzero, n_l_nonzero, n_l_cols, sparse_pointer, &
                 sparse_index, H_real_sparse, 0, S_real_sparse, pexsi_info)
         endif

         if(pexsi_info /= 0) then
            call elsi_stop("PEXSI not able to load H/S matrix.",caller)
         endif

         ! Solve the eigenvalue problem
         call f_ppexsi_dft_driver(pexsi_plan, pexsi_options, n_electrons, &
              mu_Pexsi, n_electrons_pexsi, mu_min_inertia, mu_max_inertia, &
              n_total_inertia_iter, n_total_pexsi_iter, pexsi_info)
         
         if(pexsi_info /= 0) then
            call elsi_stop("PEXSI DFT Driver not able to solve problem.",caller)
         endif

         ! Get the results
         call f_ppexsi_retrieve_real_symmetric_dft_matrix(pexsi_plan, &
              D_real_sparse, ED_real_sparse, FD_real_sparse, e_tot_H, &
              e_tot_S, f_tot, pexsi_info)

         if(pexsi_info /= 0) then
            call elsi_stop("PEXSI not able to retrieve solution.",caller)
         endif

         if(myid == 0) then
            write(*,*) "Total energy (H*DM)         = ", e_tot_H
            write(*,*) "Total energy (S*EDM)        = ", e_tot_S
            write(*,*) "Total free energy           = ", f_tot
         endif

      case DEFAULT
         call elsi_stop("No method has been chosen. "//&
                        "Please choose method ELPA, OMM, or PEXSI",caller)
   end select

   call MPI_BARRIER(mpi_comm_global,mpierr)
   call elsi_stop_solve_evp_time()

end subroutine

!>
!!  This routine transforms the eigenvalue problem from scalapack dense to pexsi sparse
!!
subroutine scalapack_dense_to_pexsi_sparse(H_external, S_external, &
           n_l_rows_external, n_l_cols_external, mpi_comm_external, &
           blacs_ctxt_external, sc_desc_external)

   implicit none
   include "mpif.h"

   ! Arguments
   integer :: n_l_rows_external
   integer :: n_l_cols_external
   real*8  :: H_external(n_l_rows_external, n_l_cols_external) 
   real*8  :: S_external(n_l_rows_external, n_l_cols_external) 
   integer :: mpi_comm_external
   integer :: blacs_ctxt_external
   integer :: sc_desc_external(9)

   ! External functions
   integer, external :: numroc

   ! Local variables
   integer :: n_p_rows_external
   integer :: n_p_cols_external
   integer :: my_p_row_external
   integer :: my_p_col_external
   integer :: mpi_comm_row_external
   integer :: mpi_comm_col_external
   
   real*8,  allocatable :: buffer(:,:)
   integer :: n_buffer_nonzero
   integer :: n_buffer_bcast_nonzero
   integer :: id_sent
   real*8,  allocatable :: buffer_sparse(:)
   integer, allocatable :: buffer_sparse_index(:)
   integer, allocatable :: buffer_sparse_pointer(:)
   real*8,  allocatable :: buffer_bcast_sparse(:)
   integer, allocatable :: buffer_bcast_sparse_index(:)
   integer, allocatable :: buffer_bcast_sparse_pointer(:)
   integer :: sc_desc_buffer(9)
   integer :: n_l_buffer_rows
   integer :: n_l_buffer_cols
   integer :: n_b_buffer_rows
   integer :: n_b_buffer_cols
   
   integer :: n_l_cols_g
   integer :: n_cols_add
   integer :: blacs_col_offset
   integer :: pexsi_col_offset
   integer :: my_col_offset
   integer :: g_column
   integer :: offset_B
   integer :: offset
   integer :: n_elements
   integer :: i_proc, i_col, id, l_col
   
   integer, allocatable :: n_l_nonzero_column(:)

   character*100, parameter :: caller = "scalapack_dense_to_pexsi_sparse"

   call elsi_start_dist_pexsi_time()

   call elsi_statement_print("Redistribution of H and S for PEXSI Started")

   ! Caution: only n_g_nonzero is meaningful, 
   ! n_l_nonzero refers to wrong distribution
   call elsi_compute_N_nonzero(H_external,n_l_rows_external, n_l_cols_external)
   n_l_nonzero = -1

   call BLACS_Gridinfo(blacs_ctxt_external, n_p_rows_external, &
                       n_p_cols_external, my_p_row_external, &
                       my_p_col_external)

   ! Scalapack has a bug:
   if((n_p_rows_external == 1 .or. n_p_cols_external == 1) &
      .and. n_procs /= 1) then
      if(myid == 0) &
      print *, "We encountered an scalapack bug when working "//&
      "with prime process numbers and using pdtran to transform a matrix "  //&
      "to a different blacs transcriptor. We stop here and wait for a "     //&
      "scalapack patch. While this setup is an inefficient corner case, "   //&
      "restart your calculation choosing a process number which in a "      //&
      "square number in best case for optimal performance and refrain from "//&
      "using prime numbers."
      call MPI_ABORT(mpi_comm_external)
      stop
   endif

   call mpi_comm_split(mpi_comm_external, my_p_col_external, &
                       my_p_row_external, mpi_comm_row_external, mpierr)
   call mpi_comm_split(mpi_comm_external, my_p_row_external, &
                       my_p_col_external, mpi_comm_col_external, mpierr)
   
   ! Determine dimensions of buffer for transformation to PEXSI layout
   n_b_buffer_rows = n_g_rank
   n_b_buffer_cols = CEILING(1d0*n_g_rank / n_p_cols_external)
   n_l_buffer_rows = n_g_rank
   n_l_buffer_cols = numroc(n_g_rank, n_b_buffer_cols, my_p_col_external, &
                            0, n_p_cols_external)

   if(myid == 0) then
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| Initial Parallel Distribution:                        ')")
      write(*,"('|-------------------------------------------------------')")
      write(*,"('| Process grid                   : ',I5,' x ',I5)")&
         n_p_rows_external, n_p_cols_external
      write(*,"('| Buffer Dimension               : ',I5,' x ',I5)")&
         n_l_buffer_rows, n_b_buffer_cols
   endif

   ! Setup new descriptor for buffer based on the new dimensions
   call descinit(sc_desc_buffer, n_g_rank, n_g_rank, n_b_buffer_rows, n_b_buffer_cols, &
                 0, 0, blacs_ctxt_external, n_l_buffer_rows, blacs_info)

   ! PEXSI setup
   n_l_rows = n_g_rank
   n_l_cols_g = FLOOR(1d0*n_g_rank / n_procs)
   n_cols_add = n_g_rank - n_l_cols_g * n_procs
   
   if(myid == n_procs - 1) then
      n_l_cols = n_l_cols_g + n_cols_add
   else
      n_l_cols = n_l_cols_g
   endif

   ! The Hamiltonian
   call elsi_allocate(buffer, n_l_buffer_rows, n_l_buffer_cols, &
                      "buffer", caller)
   call pdtran(n_g_rank, n_g_rank, 1d0, H_external, 1, 1, &
               sc_desc_external, 0d0, buffer, 1, 1, sc_desc_buffer)

   ! Calculate number of nonzero elements on this process
   ! Here we set no the meaningful n_l_nonzero
   call elsi_allocate(n_l_nonzero_column, n_g_rank, &
                      "n_l_nonzero_column", caller)
   blacs_col_offset = 1 + my_p_col_external * n_b_buffer_cols
   call elsi_get_n_nonzero_column(buffer, n_l_buffer_rows, n_l_buffer_cols, &
                                  blacs_col_offset, n_l_nonzero_column)
   pexsi_col_offset = 1 + myid * n_l_cols_g
   n_l_nonzero = &
   sum(n_l_nonzero_column(pexsi_col_offset:pexsi_col_offset + n_l_cols - 1))

   !call elsi_vector_print(n_l_nonzero_column, n_g_rank, "n_l_nonzero_column")
   !call elsi_value_print(n_l_nonzero, "n_l_nonzero")

   call elsi_allocate(H_real_sparse, n_l_nonzero, "H_real_sparse", caller)
   call elsi_allocate(sparse_index, n_l_nonzero, "sparse_index", caller)
   call elsi_allocate(sparse_pointer, n_l_cols + 1, "sparse_pointer", caller)
   sparse_index = -1
   sparse_pointer = -1

   call elsi_get_local_N_nonzero(buffer, n_l_buffer_rows, &
                                 n_b_buffer_cols, n_buffer_nonzero)

   call elsi_allocate(buffer_sparse, n_buffer_nonzero, "buffer_sparse", caller)
   call elsi_allocate(buffer_sparse_index, n_buffer_nonzero, &
                      "buffer_sparse_index", caller)
   call elsi_allocate(buffer_sparse_pointer, n_b_buffer_cols + 1, &
                      "buffer_sparse_pointer", caller)
   call elsi_dense_to_ccs(buffer, n_l_buffer_rows, n_l_buffer_cols, &
                          buffer_sparse, n_buffer_nonzero, buffer_sparse_index, &
                          buffer_sparse_pointer)
   deallocate(buffer)

   !call elsi_vector_print(buffer_sparse, n_buffer_nonzero, "Hamiltonian sparse")
   !call elsi_vector_print(buffer_sparse_index, n_buffer_nonzero, "Hamiltonian sparse index")
   !call elsi_vector_print(buffer_sparse_pointer, n_b_buffer_cols+1, "Hamiltonian sparse pointer")

   do i_proc = 0, n_p_cols_external - 1
   
      n_buffer_bcast_nonzero = n_buffer_nonzero

      id_sent = i_proc

      !call elsi_value_print(id_sent,"id sent")

      call MPI_Bcast(n_buffer_bcast_nonzero, 1, MPI_INT, id_sent, &
                     mpi_comm_external, mpierr)
      
      call elsi_allocate(buffer_bcast_sparse, n_buffer_bcast_nonzero, &
                         "buffer_bcast_spase", caller)
      call elsi_allocate(buffer_bcast_sparse_index, n_buffer_bcast_nonzero, &
                         "buffer_bcast_sparse_index", caller)
      call elsi_allocate(buffer_bcast_sparse_pointer,n_b_buffer_cols + 1, &
                         "buffer_bcast_sparse_pointer", caller)

      if(myid == id_sent) then
         buffer_bcast_sparse         = buffer_sparse
         buffer_bcast_sparse_index   = buffer_sparse_index
         buffer_bcast_sparse_pointer = buffer_sparse_pointer
      else
         buffer_bcast_sparse         = 0
         buffer_bcast_sparse_index   = 0
         buffer_bcast_sparse_pointer = 0
      endif

      ! However, no inter BLACS is possible so we have to do some by hand
      ! communication and broadcast the results over all the processes
      call MPI_Bcast(buffer_bcast_sparse, n_buffer_bcast_nonzero, &
                     MPI_DOUBLE, id_sent, mpi_comm_external, mpierr)
      call MPI_Bcast(buffer_bcast_sparse_index, n_buffer_bcast_nonzero, &
                     MPI_INT, id_sent, mpi_comm_external, mpierr)  
      call MPI_Bcast(buffer_bcast_sparse_pointer, n_b_buffer_cols + 1, &
                     MPI_INT, id_sent, mpi_comm_external, mpierr)  
 
      !call elsi_vector_print(buffer_bcast_sparse_pointer, n_b_buffer_cols+1,&
      !      "Hamiltonian sparse pointer bcast")

      blacs_col_offset = 1 + i_proc * n_b_buffer_cols
      pexsi_col_offset = 1 + myid * n_l_cols_g
      my_col_offset    = 1 + pexsi_col_offset - blacs_col_offset

      ! Fill elements from buffer
      do i_col = my_col_offset, my_col_offset + n_l_cols - 1
        
         if(i_col > 0 .and. i_col <= n_b_buffer_cols) then
         
            ! Get the global column
            !print *, "process ", myid, " pexsi_col_offset ", pexsi_col_offset
            g_column = i_col + blacs_col_offset - 1
            !print *, "process ", myid, " g_column ", g_column
         
            ! Get the offset in the buffer for this column
            offset_B = buffer_bcast_sparse_pointer(i_col)
            !print *, "process ", myid, " offset_B ", offset_B

            ! How many elements are in this column
            n_elements = buffer_bcast_sparse_pointer(i_col+1) &
                         - buffer_bcast_sparse_pointer(i_col)
            !print *, "process ", myid, " n_elements ", n_elements

            ! Local offset
            l_col = i_col - my_col_offset + 1
            !print *, "process ", myid, " l_col ", l_col
         
            offset = sum(n_l_nonzero_column(pexsi_col_offset:g_column-1)) + 1
            !print *, "process ", myid, " offset ", offset
         
            ! Populate sparse matrix
            sparse_pointer(l_col) = offset
            !print *, "process ", myid, " filling ", offset, " to ", &
            !   offset+n_elements-1, " from ", offset_B, " to ", &
            !   offset_B+n_elements-1 

            H_real_sparse(offset:offset+n_elements-1) = &
            buffer_bcast_sparse(offset_B:offset_B+n_elements-1)
            sparse_index(offset:offset+n_elements-1) = &
            buffer_bcast_sparse_index(offset_B:offset_B+n_elements-1)

        endif

     enddo
      
     deallocate(buffer_bcast_sparse)
     deallocate(buffer_bcast_sparse_index)
     deallocate(buffer_bcast_sparse_pointer)

  enddo

     ! Last element in sparse pointer is n_l_nonzero + 1
     sparse_pointer(n_l_cols + 1) = n_l_nonzero + 1

     !call elsi_stop("Stop here","dense_to_pexsi")

     ! The Overlap Matrix
     call elsi_allocate(buffer, n_l_buffer_rows, n_b_buffer_cols, &
                        "buffer", caller)
     buffer = 0d0 
     call pdtran(n_g_rank, n_g_rank, &
                 1d0, S_external, 1, 1, sc_desc_external, &
                 0d0, buffer, 1, 1, sc_desc_buffer)

     call elsi_allocate(S_real_sparse, n_l_nonzero, "S_real_sparse", caller)
     call elsi_dense_to_ccs_by_pattern(buffer, n_l_buffer_rows, n_l_buffer_cols, &
                                       buffer_sparse, n_buffer_nonzero, buffer_sparse_index, &
                                       buffer_sparse_pointer)
     deallocate(buffer)
   
     do i_proc = 0, n_p_cols_external - 1
   
        n_buffer_bcast_nonzero = n_buffer_nonzero

        id_sent = i_proc

        call MPI_Bcast(n_buffer_bcast_nonzero, 1, MPI_INT, id_sent, &
                       mpi_comm_external, mpierr)

        call elsi_allocate(buffer_bcast_sparse, n_buffer_bcast_nonzero, &
                           "buffer_bcast_sparse", caller)
        call elsi_allocate(buffer_bcast_sparse_index, n_buffer_bcast_nonzero, &
                           "buffer_bcast_sparse_index", caller)
        call elsi_allocate(buffer_bcast_sparse_pointer, n_b_buffer_cols + 1, &
                           "buffer_bcast_sparse_pointer", caller)

        if(myid == id_sent) then
           buffer_bcast_sparse         = buffer_sparse
           buffer_bcast_sparse_index   = buffer_sparse_index
           buffer_bcast_sparse_pointer = buffer_sparse_pointer
        else
           buffer_bcast_sparse         = 0
           buffer_bcast_sparse_index   = 0
           buffer_bcast_sparse_pointer = 0
        endif

        ! However, no inter BLACS is possible so we have to do some by hand
        ! communication and broadcast the results over all the processes
        call MPI_Bcast(buffer_bcast_sparse, n_buffer_bcast_nonzero, &
                       MPI_DOUBLE, id_sent, mpi_comm_external, mpierr)
        call MPI_Bcast(buffer_bcast_sparse_index, n_buffer_bcast_nonzero, &
                       MPI_INT, id_sent, mpi_comm_external, mpierr)  
        call MPI_Bcast(buffer_bcast_sparse_pointer, n_b_buffer_cols + 1, &
                       MPI_INT, id_sent, mpi_comm_external, mpierr)  
  
        blacs_col_offset = 1 + i_proc * n_b_buffer_cols
        pexsi_col_offset = 1 + myid * n_l_cols_g
        my_col_offset    = 1 + pexsi_col_offset - blacs_col_offset

        ! Fill elements from buffer
        do i_col = my_col_offset, my_col_offset + n_l_cols - 1
        
           if(i_col > 0 .and. i_col <= n_b_buffer_cols) then

              ! Get the global column
              g_column = i_col + blacs_col_offset - 1
         
              ! Get the offset in the buffer for this column
              offset_B = buffer_bcast_sparse_pointer(i_col)

              ! How many elements are in this column
              n_elements = buffer_bcast_sparse_pointer(i_col+1) &
                           - buffer_bcast_sparse_pointer(i_col)

              ! Local offset
              l_col = i_col - my_col_offset + 1

              ! Populate sparse matrix
              offset = sparse_pointer(l_col)
              S_real_sparse(offset:offset+n_elements-1) = &
              buffer_bcast_sparse(offset_B:offset_B+n_elements-1)

           endif

        enddo
       
        deallocate(buffer_bcast_sparse)
        deallocate(buffer_bcast_sparse_index)
        deallocate(buffer_bcast_sparse_pointer)

   enddo

   deallocate(buffer_sparse)
   deallocate(buffer_sparse_index)
   deallocate(buffer_sparse_pointer)
   deallocate(n_l_nonzero_column)

   ! TODO Check if overlap is unity
   overlap_is_unity = .False.
   
   call elsi_stop_dist_pexsi_time()
   
   call elsi_statement_print("Redistribution of H and S for PEXSI Done")

end subroutine

end module ELSI
