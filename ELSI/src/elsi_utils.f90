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
   public :: elsi_statement_print
   public :: elsi_allocate
   public :: elsi_allocate_matrices
   public :: elsi_deallocate_matrices
   public :: elsi_stop

   interface elsi_set_hamiltonian
      module procedure elsi_set_real_hamiltonian,&
                       elsi_set_complex_hamiltonian
   end interface

   interface elsi_set_overlap
      module procedure elsi_set_real_overlap,&
                       elsi_set_complex_overlap
   end interface

   interface elsi_allocate
      module procedure elsi_allocate_int_vector, &
                       elsi_allocate_real_vector, &
                       elsi_allocate_complex_vector, &
                       elsi_allocate_int_matrix, &
                       elsi_allocate_real_matrix, &
                       elsi_allocate_complex_matrix
   end interface

contains

!>
!! This routine prints a statement.
!!
subroutine elsi_statement_print(message)

   implicit none

   character(len=*), intent(in) :: message

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

   real*8, allocatable, intent(inout) :: vector(:)
   integer, intent(in) :: n_elements
   character(len=*), intent(in) :: vectorname
   character(len=*), intent(in) :: caller

   integer :: error

   character*200 :: message

   real*8 :: memory

   memory = 8d0*n_elements/2d0**20

   allocate(vector(n_elements),stat=error)

   if(error > 0) then 
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ",&
            TRIM(vectorname),",",memory," MB needed."
      call elsi_stop(message,caller)
   endif

   vector = 0d0 

end subroutine

!>
!! This routine allocates an integer vector.
!!
subroutine elsi_allocate_int_vector(vector,n_elements,vectorname,caller)

   implicit none

   integer, allocatable, intent(inout) :: vector(:)
   integer, intent(in) :: n_elements
   character(len=*), intent(in) :: vectorname
   character(len=*), intent(in) :: caller

   integer :: error

   character*200 :: message

   real*8 :: memory

   memory = 4d0*n_elements/2d0**20

   allocate(vector(n_elements),stat=error)

   if(error > 0) then 
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ",&
            TRIM(vectorname),",",memory," MB needed."
      call elsi_stop(message,caller)
   endif

   vector = 0

end subroutine

!>
!! This routine allocates a complex vector.
!!
subroutine elsi_allocate_complex_vector(vector,n_elements,vectorname,caller)

   implicit none

   complex*16, allocatable, intent(inout) :: vector(:)
   integer, intent(in) :: n_elements
   character(len=*), intent(in) :: vectorname
   character(len=*), intent(in) :: caller

   integer :: error

   character*200 :: message

   real*8 :: memory

   memory = 16d0*n_elements/2d0**20

   allocate(vector(n_elements),stat=error)

   if(error > 0) then
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ",&
            TRIM(vectorname),",",memory," MB needed."
      call elsi_stop(message,caller)
   endif

   vector = CMPLX(0d0,0d0)

end subroutine

!>
!! This routine allocates a real matrix.
!!
subroutine elsi_allocate_real_matrix(matrix,n_rows,n_cols,matrixname,caller)

   implicit none

   real*8, allocatable, intent(inout) :: matrix(:,:)
   integer, intent(in) :: n_rows
   integer, intent(in) :: n_cols
   character(len=*), intent(in) :: matrixname
   character(len=*), intent(in) :: caller

   integer :: error

   character*200 :: message

   real*8 :: memory

   memory = 8d0*n_rows*n_cols/2d0**20

   allocate(matrix(n_rows,n_cols),stat=error)

   if(error > 0) then 
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ",&
            TRIM(matrixname),",",memory," MB needed."
      call elsi_stop(message,caller)
   endif

   matrix = 0d0 

end subroutine

!>
!! This routine allocates an integer matrix.
!!
subroutine elsi_allocate_int_matrix(matrix,n_rows,n_cols,matrixname,caller)

   implicit none

   integer, allocatable, intent(inout) :: matrix(:,:)
   integer, intent(in) :: n_rows
   integer, intent(in) :: n_cols
   character(len=*), intent(in) :: matrixname
   character(len=*), intent(in) :: caller

   integer :: error

   character*200 :: message

   real*8 :: memory

   memory = 4d0*n_rows*n_cols/2d0**20

   allocate(matrix(n_rows,n_cols),stat=error)

   if(error > 0) then 
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ",&
            TRIM(matrixname),",",memory," MB needed."
      call elsi_stop(message,caller)
   endif

   matrix = 0 

end subroutine

!>
!! This routine allocates a complex matrix.
!!
subroutine elsi_allocate_complex_matrix(matrix,n_rows,n_cols,matrixname,caller)

   implicit none

   complex*16, allocatable, intent(inout) :: matrix(:,:)
   integer, intent(in) :: n_rows
   integer, intent(in) :: n_cols
   character(len=*), intent(in) :: matrixname
   character(len=*), intent(in) :: caller

   integer :: error

   character*200 :: message

   real*8 :: memory

   memory = 16d0*n_rows*n_cols/2d0**20

   allocate(matrix(n_rows,n_cols),stat=error)

   if(error > 0) then 
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ",&
            TRIM(matrixname),",",memory," MB needed."
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

   character(len=*), intent(in) :: message
   character(len=*), intent(in) :: caller

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
!! This routine prepares the matrices.
!!
subroutine elsi_allocate_matrices()

   implicit none

   character*40, parameter :: caller = "elsi_allocate_matrices"

   select case (method)
      case (ELPA)
         if(.not.ALLOCATED(eigenvalues)) then
            call elsi_allocate(eigenvalues,n_g_size,"eigenvalues",caller)
         endif
         eigenvalues = 0d0
         select case (mode)
            case (COMPLEX_VALUES)
               if(.not.ALLOCATED(C_complex)) then
                  call elsi_allocate(C_complex,n_l_rows,n_l_cols,"C_complex",caller)
               endif
               C_complex = CMPLX(0d0,0d0)
            case (REAL_VALUES)
               if(.not.ALLOCATED(C_real)) then
                  call elsi_allocate(C_real,n_l_rows,n_l_cols,"C_real",caller)
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
                  call m_allocate(D_omm,n_g_size,n_g_size,"pddbc")
               endif
               if(.not.Coeff_omm%is_initialized) then
                  call m_allocate(Coeff_omm,n_states,n_g_size,"pddbc")
               endif
            case (REAL_VALUES)
               if(.not.D_omm%is_initialized) then
                  call m_allocate(D_omm,n_g_size,n_g_size,"pddbc")
               endif
               if(.not.Coeff_omm%is_initialized) then
                  call m_allocate(Coeff_omm,n_states,n_g_size,"pddbc")
               endif
            case DEFAULT
               call elsi_stop(" No mode has been chosen. "//&
                              " Please choose REAL_VALUES or COMPLEX_VALUES. ",&
                              caller)
         end select

      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting...",caller)
   end select
end subroutine

!>
!! This routine sets the real hamiltonian matrix.
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
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the complex hamiltonian matrix.
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
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the real overlap matrix.
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
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the complex overlap matrix.
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
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported method has been chosen. "//&
                        " Please choose ELPA, LIBOMM, PEXSI, or CHESS. "//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine deallocates the matrices.
!!
subroutine elsi_deallocate_matrices()

   implicit none

   ! Nullify pointers
   if(ASSOCIATED(H_real))       nullify(H_real)
   if(ASSOCIATED(H_complex))    nullify(H_complex)
   if(ASSOCIATED(S_real))       nullify(S_real)
   if(ASSOCIATED(S_complex))    nullify(S_complex)

   ! Free Memory
   ! ELPA
   if(ALLOCATED(C_real))        deallocate(C_real)
   if(ALLOCATED(C_complex))     deallocate(C_complex)
   if(ALLOCATED(eigenvalues))   deallocate(eigenvalues)
   if(ALLOCATED(D_elpa))        deallocate(D_elpa)
   if(ALLOCATED(occ_elpa))      deallocate(occ_elpa)
   ! PEXSI
   if(ALLOCATED(H_real_pexsi))  deallocate(H_real_pexsi)
   if(ALLOCATED(S_real_pexsi))  deallocate(S_real_pexsi)
   if(ALLOCATED(D_pexsi))       deallocate(D_pexsi)
   if(ALLOCATED(ED_pexsi))      deallocate(ED_pexsi)
   if(ALLOCATED(FD_pexsi))      deallocate(FD_pexsi)
   if(ALLOCATED(row_ind_pexsi)) deallocate(row_ind_pexsi)
   if(ALLOCATED(col_ptr_pexsi)) deallocate(col_ptr_pexsi)
   ! libOMM
   if(H_omm%is_initialized)     call m_deallocate(H_omm)
   if(S_omm%is_initialized)     call m_deallocate(S_omm)
   if(D_omm%is_initialized)     call m_deallocate(D_omm)
   if(Coeff_omm%is_initialized) call m_deallocate(Coeff_omm)
   if(T_omm%is_initialized)     call m_deallocate(T_omm)

end subroutine

end module ELSI_UTILS
