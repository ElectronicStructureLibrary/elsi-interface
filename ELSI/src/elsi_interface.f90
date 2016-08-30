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
   public :: AUTO, ELPA, LIBOMM, PEXSI, CHESS
   public :: REAL_VALUES, COMPLEX_VALUES
   public :: BLOCK_CYCLIC

   ! From ELSI_MPI_TOOLS
   public :: elsi_set_mpi
   public :: elsi_set_blacs

   ! Public routines
   public :: elsi_init            !< Initialize
   public :: elsi_customize_elpa  !< Override ELPA default
   public :: elsi_customize_omm   !< Override OMM default
   public :: elsi_customize_pexsi !< Override PEXSI default
   public :: elsi_ev_real         !< Compute eigenvalues and eigenvectors
   public :: elsi_ev_complex      !< Compute eigenvalues and eigenvectors
   public :: elsi_dm_real         !< Compute density matrix
   public :: elsi_dm_complex      !< Compute density matrix
   public :: elsi_finalize        !< Clean memory and print timings

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
         if(.not.allocated(eigenvalues)) then
            call elsi_allocate(eigenvalues,n_g_rank,"eigenvalues",caller)
         else
            eigenvalues = 0d0
         endif
         select case (mode)
            case (COMPLEX_VALUES)
               if(.not.allocated(C_complex)) then
                  call elsi_allocate(C_complex,n_l_rows,n_l_cols,"C_complex",caller)
               else
                  C_complex = CMPLX(0d0,0d0)
               endif
            case (REAL_VALUES)
               if(.not.allocated(C_real)) then
                  call elsi_allocate(C_real,n_l_rows,n_l_cols,"C_real",caller)
               else
                  C_real = 0d0
               endif
            case DEFAULT
               call elsi_stop(" No mode has been chosen. "//&
                              " Please choose REAL_VALUES or COMPLEX_VALUES. ",&
                              caller)
         end select

      case (LIBOMM)
         select case (mode)
            case (COMPLEX_VALUES)
               if(.not.D_omm%is_initialized) &
                  call m_allocate(D_omm,n_g_rank,n_g_rank,"pddbc")
               if(.not.Coeff_omm%is_initialized) &
                  call m_allocate(Coeff_omm,n_states,n_g_rank,"pddbc")
            case (REAL_VALUES)
               if(.not.D_omm%is_initialized) &
                  call m_allocate(D_omm,n_g_rank,n_g_rank,"pddbc")
               if(.not.Coeff_omm%is_initialized) &
                  call m_allocate(Coeff_omm,n_states,n_g_rank,"pddbc")
            case DEFAULT
               call elsi_stop(" No mode has been chosen. "//&
                              " Please choose REAL_VALUES or COMPLEX_VALUES. ",&
                              caller)
         end select

      case (PEXSI)
         ! Nothing to be done here
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
         energy_out = e_tot_H
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
!!  This routine deallocates the matrices.
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

end subroutine

!>
!!  This routine finalizes ELSI.
!!
subroutine elsi_finalize()

   implicit none
   include "mpif.h"

   call MPI_BARRIER(mpi_comm_global, mpierr)

   call elsi_deallocate_matrices()

   if(method == PEXSI) call f_ppexsi_plan_finalize(pexsi_plan, pexsi_info)
   
   call elsi_stop_total_time()
   call elsi_print_timers()

end subroutine

!========================
! ELSI routines for ELPA
!========================

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
                  success = elpa_cholesky_complex_double(n_g_rank, S_complex, n_l_rows, &
                                        n_b_rows, n_l_cols, mpi_comm_row, &
                                        mpi_comm_col, .False.)
                  if(.not.success) then
                     call elsi_stop("Cholesky decomposition failed.",caller)
                  endif

                  ! compute U^-1 -> S
                  success =  elpa_invert_trm_complex_double(n_g_rank, S_complex, n_l_rows, &
                                          n_b_rows, n_l_cols, mpi_comm_row, &
                                          mpi_comm_col, .False.)
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
                  success = elpa_cholesky_real_double(n_g_rank, S_real, n_l_rows, &
                                     n_b_rows, n_l_cols, mpi_comm_row, &
                                     mpi_comm_col, .False.)
                  if(.not.success) then
                     call elsi_stop("Cholesky decomposition failed.",caller)
                  endif

                  ! compute U^-1 -> S
                  success = elpa_invert_trm_real_double(n_g_rank, S_real, n_l_rows, &
                                       n_b_rows, n_l_cols, mpi_comm_row, &
                                       mpi_comm_col, .False.)
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
            success = solve_evp_complex_2stage_double(n_g_rank, n_states, &
                         H_complex, n_l_rows, eigenvalues, &
                         C_complex, n_l_rows, n_b_rows, n_l_cols, &
                         mpi_comm_row, mpi_comm_col, mpi_comm_global)
         case (REAL_VALUES)
            success = solve_evp_real_2stage_double(n_g_rank, n_states, H_real, &
                         n_l_rows, eigenvalues, C_real, n_l_rows, &
                         n_b_rows, n_l_cols, mpi_comm_row, mpi_comm_col, &
                         mpi_comm_global)
      end select
   else ! 1-stage solver
      call elsi_statement_print(" Starting ELPA 1-stage solver")
      select case (mode)
         case (COMPLEX_VALUES)
            success = solve_evp_complex_1stage_double(n_g_rank, n_states, H_complex, &
                         n_l_rows, eigenvalues, C_complex, n_l_rows, &
                         n_b_rows, n_l_cols, mpi_comm_row, mpi_comm_col)
         case (REAL_VALUES)
            success = solve_evp_real_1stage_double(n_g_rank, n_states, H_real, &
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

!=======================
! ELSI routines for OMM
!=======================

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
            success = elpa_cholesky_complex_double(n_g_rank, S_omm%zval, n_l_rows, n_b_rows, n_l_cols, &
                                  mpi_comm_row, mpi_comm_col, .false.)
         case (REAL_VALUES)
            ! Compute S = (U^T)U, U -> S
            success = elpa_cholesky_real_double(n_g_rank, S_omm%dval, n_l_rows, n_b_rows, n_l_cols, &
                               mpi_comm_row, mpi_comm_col, .false.)
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

!=========================
! ELSI routines for PEXSI
!=========================

!>
!!  This routine converts Halmitonian and overlap matrix stored in
!!  2D block-cyclic distributed dense format to 1D block distributed
!!  sparse CCS format, which can be used as input by PEXSI.
!!
subroutine elsi_2dbc_to_1dcss_pexsi(H_in, S_in, cholesky)

   implicit none
   include 'mpif.h'

   real*8, intent(in) :: H_in(n_l_rows,n_l_cols)
   real*8, intent(in) :: S_in(n_l_rows,n_l_cols)
   logical, intent(in) :: cholesky

   real*8 :: matrix_aux(n_l_rows_pexsi,n_l_cols_pexsi)

   character*40, parameter :: caller = "elsi_2dbc_to_1dcss_pexsi"

   call elsi_start_2dbc_to_1dccs_time()
   call elsi_statement_print(" Matrix conversion: 2D block-cyclic dense ==> 1D CCS sparse")

   if(cholesky) then ! Number of non-zeros, row index, column pointer unknown
      ! Transform Hamiltonian: 2D block-cyclic dense ==> 1D block dense
      call elsi_2dbc_to_1db(H_in,matrix_aux)
      ! Compute n_l_nonzero and n_g_nonzero
      call elsi_get_global_n_nonzero(matrix_aux,n_l_rows_pexsi,n_l_cols_pexsi)

      ! Allocate PEXSI matrices
      call elsi_allocate(H_real_pexsi,n_l_nonzero,"H_real_pexsi",caller)
      call elsi_allocate(S_real_pexsi,n_l_nonzero,"S_real_pexsi",caller)
      call elsi_allocate(row_ind_pexsi,n_l_nonzero,"row_ind_pexsi",caller)
      call elsi_allocate(col_ptr_pexsi,(n_l_cols_pexsi+1),"col_ptr_pexsi",caller)

      ! Transform Hamiltonian: 1D block dense ==> 1D block sparse CCS
      call elsi_dense_to_ccs(matrix_aux,n_l_rows_pexsi,n_l_cols_pexsi,&
                             n_l_nonzero,H_real_pexsi,row_ind_pexsi,col_ptr_pexsi)

      ! Transform overlap: 2D block-cyclic dense ==> 1D block dense ==> 1D block sparse CCS
      call elsi_2dbc_to_1db(S_in,matrix_aux)
      call elsi_dense_to_ccs_by_pattern(matrix_aux,n_l_rows_pexsi,n_l_cols_pexsi,&
                                        n_l_nonzero,row_ind_pexsi,col_ptr_pexsi,S_real_pexsi)

   else ! Number of non-zeros, row index, column pointer known
      ! Transform Hamiltonian: 2D block-cyclic dense ==> 1D block dense ==> 1D block sparse CCS
      call elsi_2dbc_to_1db(H_in,matrix_aux)
      call elsi_dense_to_ccs_by_pattern(matrix_aux,n_l_rows_pexsi,n_l_cols_pexsi,&
                                        n_l_nonzero,row_ind_pexsi,col_ptr_pexsi,H_real_pexsi)

      ! Transform overlap: 2D block-cyclic dense ==> 1D block dense ==> 1D block sparse CCS
      call elsi_2dbc_to_1db(S_in,matrix_aux)
      call elsi_dense_to_ccs_by_pattern(matrix_aux,n_l_rows_pexsi,n_l_cols_pexsi,&
                                        n_l_nonzero,row_ind_pexsi,col_ptr_pexsi,S_real_pexsi)
   endif

   call elsi_stop_2dbc_to_1dccs_time()

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

   if(.not.allocated(D_pexsi)) then
      call elsi_allocate(D_pexsi,n_l_nonzero,"D_pexsi",caller)
   endif
   if(.not.allocated(ED_pexsi)) then 
      call elsi_allocate(ED_pexsi,n_l_nonzero,"ED_pexsi",caller)
   endif
   if(.not.allocated(FD_pexsi)) then 
      call elsi_allocate(FD_pexsi,n_l_nonzero,"FD_pexsi",caller)
   endif

   ! Load sparse matrices for PEXSI
   if(overlap_is_unity) then
      call f_ppexsi_load_real_symmetric_hs_matrix(pexsi_plan,pexsi_options,&
              n_g_rank,n_g_nonzero,n_l_nonzero,n_l_cols_pexsi,col_ptr_pexsi,&
              row_ind_pexsi,H_real_pexsi,1,S_real_pexsi,pexsi_info)
   else
      call f_ppexsi_load_real_symmetric_hs_matrix(pexsi_plan,pexsi_options,&
              n_g_rank,n_g_nonzero,n_l_nonzero,n_l_cols_pexsi,col_ptr_pexsi,&
              row_ind_pexsi,H_real_pexsi,0,S_real_pexsi,pexsi_info)
   endif

   if(pexsi_info /= 0) then
      call elsi_stop(" PEXSI not able to load H/S matrix. Exiting... ",caller)
   endif

   ! Solve the eigenvalue problem
   call elsi_statement_print(" Launch PEXSI DFT driver ")

   call f_ppexsi_dft_driver(pexsi_plan,pexsi_options,n_electrons,mu_pexsi,&
                            n_electrons_pexsi,mu_min_inertia,mu_max_inertia,&
                            n_total_inertia_iter,n_total_pexsi_iter,pexsi_info)
         
   if(pexsi_info /= 0) then
      call elsi_stop(" PEXSI DFT Driver not able to solve problem. Exiting... ",caller)
   endif

   ! Get the results
   call f_ppexsi_retrieve_real_symmetric_dft_matrix(pexsi_plan,D_pexsi,ED_pexsi,FD_pexsi,&
                                                    e_tot_H,e_tot_S,f_tot,pexsi_info)

   if(pexsi_info /= 0) then
      call elsi_stop(" PEXSI not able to retrieve solution. Exiting... ",caller)
   endif

   if(myid == 0) then
      write(*,*) "Total energy (H*DM)  = ", e_tot_H
      write(*,*) "Total energy (S*EDM) = ", e_tot_S
      write(*,*) "Total free energy    = ", f_tot
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

!======================
! ELSI solver routines
!======================

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
         call elsi_solve_evp_elpa(need_cholesky)
         call elsi_get_eigenvalues(e_val_out)
         call elsi_get_eigenvectors(e_vec_out)

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

   ! Solve eigenvalue problem
   select case (method)
      case (ELPA)
         ! Set Hamiltonian and overlap matrices
         call elsi_set_hamiltonian(H_in)
         call elsi_set_overlap(S_in)

         call elsi_get_occupied_number(occupation)
         call elsi_solve_evp_elpa(need_cholesky)
         call elsi_allocate(D_elpa,n_l_rows,n_l_cols,"D_elpa",caller)
         call elsi_compute_dm_elpa(D_out,occupation)
         call elsi_get_energy(energy_out)
      case (LIBOMM)
         ! Set Hamiltonian and overlap matrices
         call elsi_set_hamiltonian(H_in)
         call elsi_set_overlap(S_in)

         call elsi_solve_evp_omm(need_cholesky)
         call elsi_get_dm(D_out)
         call elsi_get_energy(energy_out)
      case (PEXSI)
         ! PEXSI may use different process grid to achieve
         ! the efficient 2-level parallelization
         call elsi_init_pexsi()

         ! Convert 2D block-cyclic dense Hamiltonian and overlap
         ! matrices to 1D block CCS sparse format
         call elsi_2dbc_to_1dcss_pexsi(H_in,S_in,need_cholesky)

         call elsi_solve_evp_pexsi()

         ! Convert 1D block CCS sparse density matrix to 2D
         ! block-cyclic dense format
!         call elsi_get_dm(D_out)
         call elsi_get_energy(energy_out)
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting... ",caller)
      case default
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA to compute eigenpairs. "//&
                        " Exiting... ",caller)
   end select

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
         ! Set Hamiltonian and overlap matrices
         call elsi_set_hamiltonian(H_in)
         call elsi_set_overlap(S_in)

         call elsi_get_occupied_number(occupation)
         call elsi_solve_evp_elpa(need_cholesky)
         call elsi_allocate(D_elpa,n_l_rows,n_l_cols,"D_elpa",caller)
         call elsi_compute_dm_elpa(D_out,occupation)
         call elsi_get_energy(energy_out)
      case (LIBOMM)
         ! Set Hamiltonian and overlap matrices
         call elsi_set_hamiltonian(H_in)
         call elsi_set_overlap(S_in)

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

   need_cholesky = .false.

end subroutine

end module ELSI
