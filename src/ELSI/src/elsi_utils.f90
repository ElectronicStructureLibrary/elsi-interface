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
!! This module contains utilities for ELSI, including array allocations and
!! printings for debugging.
!!

module ELSI_UTILS

   use iso_c_binding
   use ELSI_DIMENSIONS
   use MatrixSwitch

   implicit none
   private

   public :: elsi_set_hamiltonian
   public :: elsi_set_overlap
   public :: elsi_set_eigenvalue
   public :: elsi_set_eigenvector
   public :: elsi_set_density_matrix
   public :: elsi_set_row_ind
   public :: elsi_set_col_ptr
   public :: elsi_set_sparse_hamiltonian
   public :: elsi_set_sparse_overlap
   public :: elsi_set_sparse_density_matrix
   public :: elsi_statement_print
   public :: elsi_allocate
   public :: elsi_cleanup
   public :: elsi_stop
   public :: elsi_check

   interface elsi_set_hamiltonian
      module procedure elsi_set_real_hamiltonian,&
                       elsi_set_complex_hamiltonian
   end interface

   interface elsi_set_sparse_hamiltonian
      module procedure elsi_set_sparse_real_hamiltonian,&
                       elsi_set_sparse_complex_hamiltonian
   end interface

   interface elsi_set_overlap
      module procedure elsi_set_real_overlap,&
                       elsi_set_complex_overlap
   end interface

   interface elsi_set_sparse_overlap
      module procedure elsi_set_sparse_real_overlap,&
                       elsi_set_sparse_complex_overlap
   end interface

   interface elsi_set_eigenvector
      module procedure elsi_set_real_eigenvector,&
                       elsi_set_complex_eigenvector
   end interface

   interface elsi_allocate
      module procedure elsi_allocate_int_vector,&
                       elsi_allocate_real_vector,&
                       elsi_allocate_complex_vector,&
                       elsi_allocate_int_matrix,&
                       elsi_allocate_real_matrix,&
                       elsi_allocate_complex_matrix
   end interface

contains

!>
!! This routine prints a statement.
!!
subroutine elsi_statement_print(message)

   implicit none

   character(len=*), intent(in) :: message !< Message to print

   if(print_info) then
      if(myid == 0) then
         write(*,'(A)') TRIM(message)
      endif
   endif

end subroutine

!>
!! This routine allocates a real vector.
!!
subroutine elsi_allocate_real_vector(vector,n_elements,vectorname,caller)

   implicit none

   real*8, allocatable, intent(inout) :: vector(:) !< Data
   integer, intent(in) :: n_elements               !< Size
   character(len=*), intent(in) :: vectorname      !< Name
   character(len=*), intent(in) :: caller          !< Caller

   integer :: error
   character*200 :: message

   allocate(vector(n_elements),stat=error)

   if(error > 0) then 
      write(message,"(A,A)") "Error in allocating ",TRIM(vectorname)
      call elsi_stop(message,caller)
   endif

   vector = 0d0 

end subroutine

!>
!! This routine allocates an integer vector.
!!
subroutine elsi_allocate_int_vector(vector,n_elements,vectorname,caller)

   implicit none

   integer, allocatable, intent(inout) :: vector(:) !< Data
   integer, intent(in) :: n_elements                !< Size
   character(len=*), intent(in) :: vectorname       !< Name
   character(len=*), intent(in) :: caller           !< Caller

   integer :: error
   character*200 :: message

   allocate(vector(n_elements),stat=error)

   if(error > 0) then 
      write(message,"(A,A)") "Error in allocating ",TRIM(vectorname)
      call elsi_stop(message,caller)
   endif

   vector = 0

end subroutine

!>
!! This routine allocates a complex vector.
!!
subroutine elsi_allocate_complex_vector(vector,n_elements,vectorname,caller)

   implicit none

   complex*16, allocatable, intent(inout) :: vector(:) !< Data
   integer, intent(in) :: n_elements                   !< Size
   character(len=*), intent(in) :: vectorname          !< Name
   character(len=*), intent(in) :: caller              !< Caller

   integer :: error
   character*200 :: message

   allocate(vector(n_elements),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",TRIM(vectorname)
      call elsi_stop(message,caller)
   endif

   vector = CMPLX(0d0,0d0)

end subroutine

!>
!! This routine allocates a real matrix.
!!
subroutine elsi_allocate_real_matrix(matrix,n_rows,n_cols,matrixname,caller)

   implicit none

   real*8, allocatable, intent(inout) :: matrix(:,:) !< Data
   integer, intent(in) :: n_rows                     !< Size
   integer, intent(in) :: n_cols                     !< Size
   character(len=*), intent(in) :: matrixname        !< Name
   character(len=*), intent(in) :: caller            !< Caller

   integer :: error
   character*200 :: message

   allocate(matrix(n_rows,n_cols),stat=error)

   if(error > 0) then 
      write(message,"(A,A)") "Error in allocating ",TRIM(matrixname)
      call elsi_stop(message,caller)
   endif

   matrix = 0d0 

end subroutine

!>
!! This routine allocates an integer matrix.
!!
subroutine elsi_allocate_int_matrix(matrix,n_rows,n_cols,matrixname,caller)

   implicit none

   integer, allocatable, intent(inout) :: matrix(:,:) !< Data
   integer, intent(in) :: n_rows                      !< Size
   integer, intent(in) :: n_cols                      !< Size
   character(len=*), intent(in) :: matrixname         !< Name
   character(len=*), intent(in) :: caller             !< Caller

   integer :: error
   character*200 :: message

   allocate(matrix(n_rows,n_cols),stat=error)

   if(error > 0) then 
      write(message,"(A,A)") "Error in allocating ",TRIM(matrixname)
      call elsi_stop(message,caller)
   endif

   matrix = 0 

end subroutine

!>
!! This routine allocates a complex matrix.
!!
subroutine elsi_allocate_complex_matrix(matrix,n_rows,n_cols,matrixname,caller)

   implicit none

   complex*16, allocatable, intent(inout) :: matrix(:,:) !< Data
   integer, intent(in) :: n_rows                         !< Size
   integer, intent(in) :: n_cols                         !< Size
   character(len=*), intent(in) :: matrixname            !< Name
   character(len=*), intent(in) :: caller                !< Caller

   integer :: error
   character*200 :: message

   allocate(matrix(n_rows,n_cols),stat=error)

   if(error > 0) then 
      write(message,"(A,A)") "Error in allocating ",TRIM(matrixname)
      call elsi_stop(message,caller)
   endif

   matrix = CMPLX(0d0,0d0) 

end subroutine

!>
!! Clean shutdown in case of error.
!!
subroutine elsi_stop(message,caller)

   implicit none
   include "mpif.h"

   character(len=*), intent(in) :: message !< Error message
   character(len=*), intent(in) :: caller  !< Caller

   character(len=4096) :: string_message
   integer :: i_task

   do i_task = 0, n_procs - 1
      if(myid == i_task) then
         write(string_message,"(1X,'*** Proc',I5,' in ',A,': ',A)")&
               myid,TRIM(caller),TRIM(message)
         write(*,'(A)') TRIM(string_message)
      endif
      call MPI_Barrier(mpi_comm_global,mpierr)
   enddo

   if(n_procs > 1) then
      call MPI_Abort(mpi_comm_global,0,mpierr)
   endif

   stop

end subroutine

!>
!! This routine sets the real hamiltonian matrix.
!!
subroutine elsi_set_real_hamiltonian(H_in)

   implicit none

   real*8, target :: H_in(n_l_rows,n_l_cols) !< Real Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_real_hamiltonian"

   select case (method)
      case (ELPA)
         ham_real => H_in
      case (LIBOMM)
         call m_register_pdbc(ham_omm,H_in,sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the complex hamiltonian matrix.
!!
subroutine elsi_set_complex_hamiltonian(H_in)

   implicit none

   complex*16, target :: H_in(n_l_rows,n_l_cols) !< Complex Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_complex_hamiltonian"

   select case (method)
      case (ELPA)
         ham_complex => H_in
      case (LIBOMM)
         call m_register_pdbc(ham_omm,H_in,sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the sparse real Hamiltonian matrix.
!!
subroutine elsi_set_sparse_real_hamiltonian(H_in)

   implicit none

   real*8, target :: H_in(nnz_l) !< Sparse real Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_sparse_real_hamiltonian"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ham_real_ccs => H_in
      case (CHESS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the sparse complex Hamiltonian matrix.
!!
subroutine elsi_set_sparse_complex_hamiltonian(H_in)

   implicit none

   complex*16, target :: H_in(nnz_l) !< Sparse complex Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_sparse_complex_hamiltonian"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ham_complex_ccs => H_in         
      case (CHESS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the real overlap matrix.
!!
subroutine elsi_set_real_overlap(S_in)

   implicit none

   real*8, target :: S_in(n_l_rows,n_l_cols) !< Real overlap matrix

   character*40, parameter :: caller = "elsi_set_real_overlap"

   select case (method)
      case (ELPA)
         ovlp_real => S_in
      case (LIBOMM)
         call m_register_pdbc(ovlp_omm,S_in,sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the complex overlap matrix.
!!
subroutine elsi_set_complex_overlap(S_in)

   implicit none

   complex*16, target :: S_in(n_l_rows,n_l_cols) !< Complex overlap matrix

   character*40, parameter :: caller = "elsi_set_complex_overlap"

   select case (method)
      case (ELPA)
         ovlp_complex => S_in
      case (LIBOMM)
         call m_register_pdbc(ovlp_omm,S_in,sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the sparse real overlap matrix.
!!
subroutine elsi_set_sparse_real_overlap(S_in)

   implicit none

   real*8, target :: S_in(nnz_l) !< Sparse real overlap matrix

   character*40, parameter :: caller = "elsi_set_sparse_real_overlap"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ovlp_real_ccs => S_in
      case (CHESS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the sparse complex overlap matrix.
!!
subroutine elsi_set_sparse_complex_overlap(S_in)

   implicit none

   complex*16, target :: S_in(nnz_l) !< Sparse complex overlap matrix

   character*40, parameter :: caller = "elsi_set_sparse_complex_overlap"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ovlp_complex_ccs => S_in
      case (CHESS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the eigenvalues.
!!
subroutine elsi_set_eigenvalue(e_val_in)

   implicit none

   real*8, target :: e_val_in(n_g_size)  !< Eigenvalues

   character*40, parameter :: caller = "elsi_set_eigenvalue"

   select case (method)
      case (ELPA)
         eval => e_val_in
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the real eigenvectors.
!!
subroutine elsi_set_real_eigenvector(e_vec_in)

   implicit none
   
   real*8, target :: e_vec_in(n_l_rows,n_l_cols) !< Real eigenvectors

   character*40, parameter :: caller = "elsi_set_real_eigenvector"

   select case (method)
      case (ELPA)
         evec_real => e_vec_in
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the complex eigenvectors.
!!
subroutine elsi_set_complex_eigenvector(e_vec_in)

   implicit none

   complex*16, target :: e_vec_in(n_l_rows,n_l_cols) !< Complex eigenvectors

   character*40, parameter :: caller = "elsi_set_complex_eigenvector"

   select case (method)
      case (ELPA)
         evec_complex => e_vec_in
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the density matrix.
!!
subroutine elsi_set_density_matrix(D_in)

   implicit none

   real*8, target :: D_in(n_l_rows,n_l_cols) !< Density matrix

   character*40, parameter :: caller = "elsi_set_density_matrix"

   select case (method)
      case (ELPA)
         den_mat => D_in
      case (LIBOMM)
         call m_register_pdbc(den_mat_omm,D_in,sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the sparse density matrix.
!!
subroutine elsi_set_sparse_density_matrix(D_in)

   implicit none

   real*8, target :: D_in(nnz_l) !< Sparse density matrix

   character*40, parameter :: caller = "elsi_set_sparse_density_matrix"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         den_mat_ccs => D_in
      case (CHESS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the row index.
!!
subroutine elsi_set_row_ind(row_ind_in)

   implicit none

   integer, target :: row_ind_in(nnz_l) !< Row index

   character*40, parameter :: caller = "elsi_set_row_ind"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         row_ind_ccs => row_ind_in
      case (CHESS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the column pointer.
!!
subroutine elsi_set_col_ptr(col_ptr_in)

   implicit none

   integer, target :: col_ptr_in(nnz_l) !< Column pointer

   character*40, parameter :: caller = "elsi_set_cplumn_pointer"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         col_ptr_ccs => col_ptr_in
      case (CHESS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine frees memory.
!!
subroutine elsi_cleanup()

   implicit none

   ! Nullify pointers
   if(ASSOCIATED(ham_real))         nullify(ham_real)
   if(ASSOCIATED(ham_complex))      nullify(ham_complex)
   if(ASSOCIATED(ovlp_real))        nullify(ovlp_real)
   if(ASSOCIATED(ovlp_complex))     nullify(ovlp_complex)
   if(ASSOCIATED(evec_real))        nullify(evec_real)
   if(ASSOCIATED(evec_complex))     nullify(evec_complex)
   if(ASSOCIATED(eval))             nullify(eval)
   if(ASSOCIATED(den_mat))          nullify(den_mat)
   if(ASSOCIATED(ham_real_ccs))     nullify(ham_real_ccs)
   if(ASSOCIATED(ham_complex_ccs))  nullify(ham_complex_ccs)
   if(ASSOCIATED(ovlp_real_ccs))    nullify(ovlp_real_ccs)
   if(ASSOCIATED(ovlp_complex_ccs)) nullify(ovlp_complex_ccs)
   if(ASSOCIATED(den_mat_ccs))      nullify(den_mat_ccs)
   if(ASSOCIATED(row_ind_ccs))      nullify(row_ind_ccs)
   if(ASSOCIATED(col_ptr_ccs))      nullify(col_ptr_ccs)

   ! ELPA
   if(ALLOCATED(ham_real_elpa))     deallocate(ham_real_elpa)
   if(ALLOCATED(ham_complex_elpa))  deallocate(ham_complex_elpa)
   if(ALLOCATED(ovlp_real_elpa))    deallocate(ovlp_real_elpa)
   if(ALLOCATED(ovlp_complex_elpa)) deallocate(ovlp_complex_elpa)
   if(ALLOCATED(evec_real_elpa))    deallocate(evec_real_elpa)
   if(ALLOCATED(evec_complex_elpa)) deallocate(evec_complex_elpa)
   if(ALLOCATED(eval_elpa))         deallocate(eval_elpa)
   if(ALLOCATED(den_mat_elpa))      deallocate(den_mat_elpa)
   if(ALLOCATED(occ_elpa))          deallocate(occ_elpa)

   ! PEXSI
   if(ALLOCATED(ham_real_pexsi))    deallocate(ham_real_pexsi)
   if(ALLOCATED(ovlp_real_pexsi))   deallocate(ovlp_real_pexsi)
   if(ALLOCATED(den_mat_pexsi))     deallocate(den_mat_pexsi)
   if(ALLOCATED(e_den_mat_pexsi))   deallocate(e_den_mat_pexsi)
   if(ALLOCATED(f_den_mat_pexsi))   deallocate(f_den_mat_pexsi)
   if(ALLOCATED(row_ind_pexsi))     deallocate(row_ind_pexsi)
   if(ALLOCATED(col_ptr_pexsi))     deallocate(col_ptr_pexsi)

   ! libOMM
   if(ham_omm%is_initialized)       call m_deallocate(ham_omm)
   if(ovlp_omm%is_initialized)      call m_deallocate(ovlp_omm)
   if(den_mat_omm%is_initialized)   call m_deallocate(den_mat_omm)
   if(coeff_omm%is_initialized)     call m_deallocate(coeff_omm)
   if(t_den_mat_omm%is_initialized) call m_deallocate(t_den_mat_omm)

end subroutine

!>
!! This routine guarantees that there are no mutually conflicting parameters.
!!
subroutine elsi_check()

   implicit none

   character*40, parameter :: caller = "elsi_check"

   ! The solver
   if(method == 0) then
      call elsi_stop(" AUTO not yet available."//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",caller)
   else if(method == 1) then
      if(storage /= 0) then
         call elsi_stop(" ELPA has been chosen as the solver."//&
                        " Please choose DENSE matrix format."//&
                        " Exiting...",caller)
      endif
   else if(method == 2) then
      if(parallelism == 0) then
         call elsi_stop(" libOMM has been chosen as the solver."//&
                        " Please choose MULTI_PROC parallelism."//&
                        " Exiting...",caller)
      endif
      if(storage /= 0) then
         call elsi_stop(" libOMM has been chosen as the solver."//&
                        " Please choose DENSE matrix format."//&
                        " Exiting...",caller)
      endif
      if(MOD(NINT(n_electrons),2) /= 0) then
         call elsi_stop(" libOMM has been chosen as the solver."//&
                        " The current implementation of libOMM does not"//&
                        " support fractional occupation numbers. This"//&
                        " means number of electrons in non-spin-polarized"//&
                        " system cannot be odd. Exiting...",caller)
      endif
   else if(method == 3) then
      if(parallelism == 0) then
         call elsi_stop(" PEXSI has been chosen as the solver."//&
                        " Please choose MULTI_PROC parallelism."//&
                        " Exiting...",caller)
      endif
      if(n_g_size < n_procs) then
         call elsi_stop(" PEXSI has been chosen as the solver."//&
                        " The (global) size of matrix is too small for"//&
                        " this number of processes. Exiting...",caller)
      endif
   else if(method == 4) then
      call elsi_stop(" CHESS not yet available. "//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",caller)
   else
      call elsi_stop(" No supported solver has been chosen."//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",caller)
   endif

   ! The parallelism
   if(parallelism < 0 .or. parallelism > 1) then
      call elsi_stop(" No supported parallelism has been chosen."//&
                     " Please choose SINGLE_PROC or MULTI_PROC parallelism."//&
                     " Exiting...",caller)
   endif

   ! The matrix format
   if(storage >0 .and. storage < 5) then ! Sparse foramt
      if(.not.sparsity_pattern_ready) then
         call elsi_stop(" A sparse format has been chosen."//&
                        " Please set the sparsity pattern before calling the solver."//&
                        " Exiting...",caller)
      endif
   else if(storage /= 0) then
      call elsi_stop(" No supported format has been chosen."//&
                     " Please choose DENSE, CCS, CSC, CRS, or CSR format."//&
                     " Exiting...",caller)
   endif

end subroutine

end module ELSI_UTILS
