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

   public :: elsi_value_print
   public :: elsi_vector_print
   public :: elsi_matrix_print
   public :: elsi_statement_print
   public :: elsi_allocate

   interface elsi_vector_print
      module procedure elsi_int_vector_print, &
                       elsi_real_vector_print
   end interface

   interface elsi_value_print
      module procedure elsi_int_value_print, &
                       elsi_real_value_print
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

subroutine elsi_matrix_print(matrix, n_rows, n_cols, matrixname)

   implicit none

   integer, intent(in) :: n_rows
   integer, intent(in) :: n_cols
   real*8, intent(in) :: matrix(n_rows,n_cols) 
   character(len=*), intent(in) :: matrixname 

   integer :: id, i_row, i_col
   
   if(myid == 0) print *, trim(matrixname)

   do id = 0, n_procs - 1
      if(myid == id) then
         do i_row = 1, n_rows
            do i_col = 1, n_cols
               write(*,"(A,I6,A,I6,A,I6,A,F13.6)") " Process ",id,&
               " Row ",i_row," Col ",i_col," Value ",matrix(i_row,i_col)
            enddo
         enddo
      endif
      call MPI_Barrier(mpi_comm_global,mpierr)
   enddo 

end subroutine

subroutine elsi_int_vector_print(vector, n_dim, vectorname)

   implicit none

   integer, intent(in) :: n_dim
   integer, intent(in) :: vector(n_dim)
   character(len=*), intent(in) :: vectorname

   integer :: id, i_dim

   if(myid == 0) print *, trim(vectorname)

   do id = 0, n_procs - 1
      if(myid == id) then
         do i_dim = 1, n_dim
            write(*,"(A,I6,A,I10,A,I13)") " Process ",id,&
            " Entry ",i_dim," Value ",vector(i_dim)
         enddo
      endif
      call MPI_Barrier(mpi_comm_global,mpierr)
   enddo

end subroutine

subroutine elsi_int_value_print(val, valname)

   implicit none

   integer :: val
   character(len=*), intent(in) :: valname

   integer :: id

   if(myid == 0) print *, trim(valname)

   do id = 0, n_procs - 1
      if(myid == id) then
         write(*,"(A,I6,A,I13)") " Process ",id," Value ",val
      endif
      call MPI_Barrier(mpi_comm_global,mpierr)
   enddo

end subroutine

subroutine elsi_real_value_print(val, valname)

   implicit none

   real*8 :: val
   character(len=*), intent(in) :: valname

   integer :: id

   if(myid == 0) print *, trim(valname)

   do id = 0, n_procs - 1
      if(myid == id) then
         write(*,"(A,I4,A,F13.6)") " Process ",id," Value ",val
      endif
      call MPI_Barrier(mpi_comm_global,mpierr)
   enddo

end subroutine

subroutine elsi_real_vector_print(vector, n_dim, vectorname)

   implicit none

   integer, intent(in) :: n_dim
   real*8, intent(in) :: vector(n_dim)
   character(len=*), intent(in) :: vectorname

   integer :: id, i_dim

   if(myid == 0) print *, trim(vectorname)

   do id = 0, n_procs - 1
      if(myid == id) then
         do i_dim = 1, n_dim
            write(*,"(A,I6,A,I10,A,F13.6)") " Process ",id, &
            " Entry ",i_dim," Value ",vector(i_dim)
         enddo
      endif
      call MPI_Barrier(mpi_comm_global,mpierr)
   enddo

end subroutine

subroutine elsi_statement_print(message)

   implicit none

   character(len=*), intent(in) :: message

   if(myid == 0) print *, trim(message)
   call MPI_Barrier(mpi_comm_global,mpierr)

end subroutine

subroutine elsi_allocate_real_vector(vector, n_elements, vectorname, caller)

   implicit none

   real*8, allocatable, intent(inout) :: vector(:)
   integer, intent(in) :: n_elements
   character(len=*), intent(in) :: vectorname
   character(len=*), intent(in) :: caller

   integer :: error
   integer :: id

   character*200 :: message

   real*8 :: memory

   memory = 8d0 * n_elements / 2d0**20

   allocate(vector(n_elements), stat = error)

   if(error > 0) then 
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ", &
            trim(vectorname), ", ", memory, " MB needed."
      call elsi_stop(message, caller)
   endif

   vector = 0d0 

end subroutine

subroutine elsi_allocate_int_vector(vector, n_elements, vectorname, caller)

   implicit none

   integer, allocatable, intent(inout) :: vector(:)
   integer, intent(in) :: n_elements
   character(len=*), intent(in) :: vectorname
   character(len=*), intent(in) :: caller

   integer :: error
   integer :: id

   character*200 :: message

   real*8 :: memory

   memory = 4d0 * n_elements / 2d0**20

   allocate(vector(n_elements), stat = error)

   if(error > 0) then 
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ", &
            trim(vectorname), ", ", memory, " MB needed."
      call elsi_stop(message, caller)
   endif

   vector = 0

end subroutine

subroutine elsi_allocate_complex_vector(vector, n_elements, vectorname, caller)

   implicit none

   complex*16, allocatable, intent(inout) :: vector(:)
   integer, intent(in) :: n_elements
   character(len=*), intent(in) :: vectorname
   character(len=*), intent(in) :: caller

   integer :: error
   integer :: id

   character*200 :: message

   real*8 :: memory

   memory = 16d0 * n_elements / 2d0**20

   allocate(vector(n_elements), stat = error)

   if(error > 0) then
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ", &
            trim(vectorname), ", ", memory, " MB needed."
      call elsi_stop(message, caller)
   endif

   vector = CMPLX(0d0,0d0)

end subroutine

subroutine elsi_allocate_real_matrix(matrix, n_rows, n_cols, matrixname, caller)

   implicit none

   real*8, allocatable, intent(inout) :: matrix(:,:)
   integer, intent(in) :: n_rows
   integer, intent(in) :: n_cols
   character(len=*), intent(in) :: matrixname
   character(len=*), intent(in) :: caller

   integer :: error
   integer :: id

   character*200 :: message

   real*8 :: memory

   memory = 8d0 * n_rows * n_cols / 2d0**20

   allocate(matrix(n_rows,n_cols), stat = error)

   if(error > 0) then 
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ", &
            trim(matrixname), ", ", memory, " MB needed."
      call elsi_stop(message, caller)
   endif

   matrix = 0d0 

end subroutine

subroutine elsi_allocate_int_matrix(matrix, n_rows, n_cols, matrixname, caller)

   implicit none

   integer, allocatable, intent(inout) :: matrix(:,:)
   integer, intent(in) :: n_rows
   integer, intent(in) :: n_cols
   character(len=*), intent(in) :: matrixname
   character(len=*), intent(in) :: caller

   integer :: error
   integer :: id

   character*200 :: message

   real*8 :: memory

   memory = 4d0 * n_rows * n_cols / 2d0**20

   allocate(matrix(n_rows,n_cols), stat = error)

   if(error > 0) then 
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ", &
            trim(matrixname), ", ", memory, " MB needed."
      call elsi_stop(message, caller)
   endif

   matrix = 0 

end subroutine

subroutine elsi_allocate_complex_matrix(matrix, n_rows, n_cols, matrixname, caller)

   implicit none

   complex*16, allocatable, intent(inout) :: matrix(:,:)
   integer, intent(in) :: n_rows
   integer, intent(in) :: n_cols
   character(len=*), intent(in) :: matrixname
   character(len=*), intent(in) :: caller

   integer :: error
   integer :: id

   character*200 :: message

   real*8 :: memory

   memory = 16d0 * n_rows * n_cols / 2d0**20

   allocate(matrix(n_rows,n_cols), stat = error)

   if(error > 0) then 
      write(message,"(A,A,A,F10.3,A)") "Insufficient memory to allocate ", &
            trim(matrixname), ", ", memory, " MB needed."
      call elsi_stop(message, caller)
   endif

   matrix = CMPLX(0d0,0d0) 

end subroutine

end module ELSI_UTILS
