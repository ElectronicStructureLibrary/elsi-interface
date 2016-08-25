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
!! This module contains wrapper functions to the MPI library.
!!

module ELSI_MPI_TOOLS

   use iso_c_binding
   use ELSI_DIMENSIONS
   use MatrixSwitch

   implicit none
   private

   !> Calculate the local matrix dimensions based on the BLACS grid
   integer, external :: numroc

   public :: elsi_init_pexsi
   public :: elsi_set_mpi
   public :: elsi_set_blacs
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
      module procedure elsi_allocate_real_vector, &
                       elsi_allocate_real_matrix, &
                       elsi_allocate_int_vector, &
                       elsi_allocate_int_matrix, &
                       elsi_allocate_complex_vector, &
                       elsi_allocate_complex_matrix
   end interface

contains

!>
!! Set MPI.
!!
subroutine elsi_set_mpi(mpi_comm_global_in, n_procs_in, myid_in)

   implicit none
   include 'mpif.h'

   integer,intent(in) :: myid_in            !< local process id
   integer,intent(in) :: n_procs_in         !< number of mpi processes
   integer,intent(in) :: mpi_comm_global_in !< global mpi communicator

   ! MPI Initialization

   external_mpi = .True.
   mpi_is_setup = .True.
   mpi_comm_global = mpi_comm_global_in
   n_procs         = n_procs_in
   myid            = myid_in

end subroutine

!>
!! PEXSI processor grid setup.
!!
subroutine elsi_init_pexsi()

   implicit none
   include "mpif.h"

   if(method == PEXSI) then

      ! Find balancing between expansion parallel and matrix inversion parallel
      if(mod(n_procs,40) == 0 .and. n_procs >= 640) then
         n_p_rows_pexsi = 40
      else if (mod(n_procs,20) == 0 .and. n_procs >= 320) then
         n_p_rows_pexsi = 20
      else if (mod(n_procs,10) == 0 .and. n_procs >= 90) then
         n_p_rows_pexsi = 10
      else if (mod(n_procs,5) == 0 .and. n_procs >= 20) then
         n_p_rows_pexsi = 5
      else if (mod(n_procs,4) == 0 .and. n_procs >= 16) then
         n_p_rows_pexsi = 4
      else if (mod(n_procs,2) == 0 .and. n_procs >= 2) then
         n_p_rows_pexsi = 2
      else
         n_p_rows_pexsi = 1
      endif

      n_p_cols_pexsi = n_procs/n_p_rows_pexsi

      ! position in process grid must not be used
      my_p_col_pexsi = -1
      my_p_row_pexsi = -1

      ! PEXSI needs a pure block distribution
      n_b_rows_pexsi = n_g_rank

      ! The last process holds all remaining columns
      n_b_cols_pexsi = FLOOR(1d0*n_g_rank/n_procs)
      if(myid == n_procs-1) then
         n_b_cols_pexsi = n_g_rank - (n_procs-1) * n_b_cols_pexsi
      endif

      n_l_rows_pexsi = n_b_rows_pexsi
      n_l_cols_pexsi = n_b_cols_pexsi

      ! TODO: doublecheck this
      ! Each master process for each pole should write an output 
      if(mod(myid,n_p_cols_pexsi*n_p_rows_pexsi) == 0) then
         pexsi_output_file_index = myid/(n_p_cols_pexsi*n_p_rows_pexsi)
      else
         pexsi_output_file_index = -1
      endif

      pexsi_plan = f_ppexsi_plan_initialize(mpi_comm_global,n_p_rows_pexsi,&
                      n_p_cols_pexsi,pexsi_output_file_index,pexsi_info)

      if(pexsi_info /= 0) then
         call elsi_stop(" PEXSI plan initialization failed. Exiting... ", &
                        "elsi_init_pexsi")
      endif

   endif

end subroutine

!>
!! Set BLACS. 
!!
subroutine elsi_set_blacs(blacs_ctxt_in, n_b_rows_in, n_b_cols_in, n_p_rows_in, &
                          n_p_cols_in, mpi_comm_row_in, mpi_comm_col_in)

   implicit none
   include "mpif.h"

   integer, intent(in) :: blacs_ctxt_in             !< BLACS context
   integer, intent(in) :: n_b_rows_in               !< Block size
   integer, intent(in) :: n_b_cols_in               !< Block size
   integer, intent(in) :: n_p_rows_in               !< Number of processes in row
   integer, intent(in) :: n_p_cols_in               !< Number of processes in column
   integer, intent(in), optional :: mpi_comm_row_in !< row communicatior for ELPA
   integer, intent(in), optional :: mpi_comm_col_in !< column communicatior for ELPA

   integer :: blacs_info

   external_blacs = .true.
   blacs_is_setup = .true.

   blacs_ctxt = blacs_ctxt_in 
   n_b_rows = n_b_rows_in
   n_b_cols = n_b_cols_in
   n_p_rows = n_p_rows_in 
   n_p_cols = n_p_cols_in
   call blacs_pcoord(blacs_ctxt,myid,my_p_row,my_p_col)
   n_l_rows = numroc(n_g_rank,n_b_rows,my_p_row,0,n_p_rows)
   n_l_cols = numroc(n_g_rank,n_b_cols,my_p_col,0,n_p_cols)
   call descinit(sc_desc,n_g_rank,n_g_rank,n_b_rows,n_b_cols,0,0,&
                 blacs_ctxt,MAX(1,n_l_rows),blacs_info)

   mpi_comm_row = mpi_comm_row_in
   mpi_comm_col = mpi_comm_col_in

   if(method == LIBOMM) then
      call ms_scalapack_setup(myid,n_procs,n_p_rows,'r',n_b_rows,icontxt=blacs_ctxt)
   endif

end subroutine

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
      call MPI_Barrier(mpi_comm_global, mpierr)
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

end module ELSI_MPI_TOOLS
