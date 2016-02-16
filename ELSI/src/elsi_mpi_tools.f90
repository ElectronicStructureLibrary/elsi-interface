!Copyright (c) 2015, ELSI consortium
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
!! This module contains wrapper functions to the MPI library 
!!

module ELSI_MPI_TOOLS

  use iso_c_binding

  use ELSI_DIMENSIONS
  use matrixswitch

  implicit none
  private

  !> Calculates the local matrix dimensions based on the BLACS grid
  integer, external :: numroc

  public :: elsi_initialize_mpi
  public :: elsi_initialize_blacs 
  public :: elsi_set_mpi
  public :: elsi_set_blacs 
  public :: elsi_get_global_row
  public :: elsi_get_global_col
  public :: elsi_get_global_dimensions
  public :: elsi_get_local_dimensions
  public :: elsi_get_processor_grid
  public :: elsi_get_comm_grid
  public :: elsi_get_myid
  public :: elsi_finalize_mpi
  public :: elsi_finalize_blacs
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
!! Initialize MPI
!!
subroutine elsi_initialize_mpi()

   implicit none
   include 'mpif.h'

   ! MPI Initialization

   external_mpi = .False.
   mpi_is_setup = .True.

   call mpi_init(mpierr)
   if (mpierr /= 0) then
      write(*,'(a)') "MPI: Failed to initialize MPI."
      stop
   end if
   call mpi_comm_rank(mpi_comm_world, myid, mpierr)
   if (mpierr /= 0) then
      write(*,'(a)') "MPI: Failed to determine MPI rank."
      stop
   end if

   mpi_comm_global = mpi_comm_world

   call mpi_comm_size(mpi_comm_global, n_procs, mpierr)
   if (mpierr /= 0) then
      write(*,'(a)') "MPI: Failed to determine MPI size."
      stop
   end if

end subroutine

!> 
!! Set MPI from external
!!
subroutine elsi_set_mpi(global_comm, n_procs_in, myid_in)

   implicit none
   include 'mpif.h'

   integer,intent(in) :: myid_in            !< local process id
   integer,intent(in) :: n_procs_in         !< number of mpi processes
   integer,intent(in) :: global_comm        !< global mpi communicator


   ! MPI Initialization

   external_mpi = .True.
   mpi_is_setup = .True.
   mpi_comm_global = global_comm
   n_procs         = n_procs_in
   myid            = myid_in
end subroutine


!> 
!! Finalize MPI
!!
subroutine elsi_finalize_mpi()

   implicit none
   include 'mpif.h'

   ! MPI Finalization
   mpi_is_setup = .False.

   call mpi_finalize(mpierr)
   if (mpierr /= 0) then
      write(*,'(a)') "MPI: Failed to finalize MPI."
      stop
   end if

end subroutine


!> 
!! Initialize BLACS Grid
!!
subroutine elsi_initialize_blacs()

   implicit none

   include "mpif.h"

   ! Exception how to deviate from blocksize for OMM
   integer :: exception(2) = (/0,0/)


   if (method == PEXSI) then

     external_blacs = .False.
     blacs_is_setup = .False.

     ! Find balancing between expansion parallel and matrix inversion parallel
     if (mod(n_procs, 40) == 0 .and. n_procs >= 640) then
        n_p_rows = 10
     else if (mod(n_procs, 20) == 0 .and. n_procs >= 320) then
        n_p_rows = 20
     else if (mod(n_procs, 10) == 0 .and. n_procs >= 90) then
        n_p_rows = 10
     else if (mod(n_procs, 5) == 0 .and. n_procs >= 20) then
        n_p_rows = 5
     else if (mod(n_procs, 4) == 0 .and. n_procs >= 16) then
        n_p_rows = 4
     else if (mod(n_procs, 2) == 0 .and. n_procs >= 2) then
        n_p_rows = 2
     else 
        n_p_rows = 1
     end if

     n_p_cols = n_procs / n_p_rows

     ! position in process grid must not be used
     my_p_col = -1
     my_p_row = -1

     ! PEXSI needs a pure block distribution.
     ! This is why we reset the blocksize here 
     n_b_rows = n_g_rank

     ! The last process takes all remaining columns
     n_b_cols = FLOOR(1d0*n_g_rank / n_procs)
     if (myid == n_procs - 1) then
        n_b_cols = n_g_rank - (n_procs-1) * n_b_cols
     end if

     n_l_rows = n_b_rows
     n_l_cols = n_b_cols
 
     ! Each master process for each pole should write an output 
     if (mod(myid,n_p_cols * n_p_rows) == 0) then
        pexsi_output_file_index = myid / (n_p_cols * n_p_rows)
     else
        pexsi_output_file_index = -1
     end if

     pexsi_plan = f_ppexsi_plan_initialize( mpi_comm_global, n_p_rows, &
           n_p_cols, pexsi_output_file_index, pexsi_info)
     if (pexsi_info /= 0) then
        call elsi_stop("Pexsi Plan initialization faild.",&
              "elsi_initialize_blacs")
     end if


   else
     
     external_blacs = .False.
     blacs_is_setup = .True.
     
     ! Define blockcyclic setup
     do n_p_cols = NINT(SQRT(REAL(n_procs))),2,-1
       if(mod(n_procs,n_p_cols) == 0 ) exit
     enddo

     n_p_rows = n_procs / n_p_cols

     ! Set up BLACS and MPI communicators

     blacs_ctxt = mpi_comm_global
     call BLACS_Gridinit( blacs_ctxt, 'C', n_p_rows, n_p_cols )
     call BLACS_Gridinfo( blacs_ctxt, n_p_rows, n_p_cols, my_p_row, my_p_col )

     call mpi_comm_split(mpi_comm_global,my_p_col,my_p_row,mpi_comm_row,mpierr)

     if (mpierr /= 0) then
       call elsi_stop("Failed to get row communicators.",&
           "elsi_initialize_blacs")
     end if
 
     call mpi_comm_split(mpi_comm_global,my_p_row,my_p_col,mpi_comm_col,mpierr)
   
     if (mpierr /= 0) then
       call elsi_stop("Failed to get column communicators.",&
           "elsi_initialize_blacs")
     end if

     n_l_rows = numroc(n_g_rank, n_b_rows, my_p_row, 0, n_p_rows)
     n_l_cols = numroc(n_g_rank, n_b_cols, my_p_col, 0, n_p_cols)

     n_l_nonzero = n_l_rows * n_l_cols
     n_g_nonzero = n_g_rank * n_g_rank

     call descinit( sc_desc, n_g_rank, n_g_rank, n_b_rows, n_b_cols, 0, 0, &
                   blacs_ctxt, MAX(1,n_l_rows), blacs_info )

     if (method == OMM_DENSE) then
        call ms_scalapack_setup (n_procs, n_p_rows, 'c', n_b_rows, exception,&
             blacs_ctxt)
     end if

   end if

  ! For DEBUG
  !call elsi_variable_status()

end subroutine

!> 
!! Set BLACS Grid from external
!!
subroutine elsi_set_blacs( blacs_ctxt_in, n_p_rows_in, n_p_cols_in, &
      my_p_row_in, my_p_col_in, mpi_comm_row_in, mpi_comm_col_in,&
      n_l_rows_in, n_l_cols_in, sc_desc_in)

   implicit none

   include "mpif.h"
   ! Exception how to deviate from blocksize for OMM
   integer :: exception(2) = (/0,0/)

  integer,intent(in) :: blacs_ctxt_in     !< local blacs context
  integer,intent(in) :: sc_desc_in(9)     !< local blacs context
  integer,intent(in) :: mpi_comm_row_in   !< row    communicatior
  integer,intent(in) :: mpi_comm_col_in   !< column communicatior
  integer,intent(in) :: my_p_row_in       !< process row    position
  integer,intent(in) :: my_p_col_in       !< process column position
  integer,intent(in) :: n_p_rows_in       !< Number of processes in row
  integer,intent(in) :: n_p_cols_in       !< Number of processes in column
  integer,intent(in) :: n_l_rows_in       !< Number of local rows
  integer,intent(in) :: n_l_cols_in       !< Number of local columns
   
  external_blacs = .True.
  blacs_is_setup = .False.

  blacs_ctxt = blacs_ctxt_in 
  n_p_rows = n_p_rows_in 
  n_p_cols = n_p_cols_in
  my_p_row = my_p_row_in
  my_p_col = my_p_col_in
  mpi_comm_row = mpi_comm_row_in
  mpi_comm_col = mpi_comm_col_in
  n_l_rows = n_l_rows_in
  n_l_cols = n_l_cols_in
  sc_desc = sc_desc_in

  if (method == OMM_DENSE) then
     call ms_scalapack_setup (n_procs, n_p_rows, "C", n_b_rows, exception,&
           blacs_ctxt)
  end if

  ! For DEBUG
  !call elsi_variable_status()

end subroutine


!> 
!! Finalize BLACS grid
!!
subroutine elsi_finalize_blacs()

   implicit none

   if(blacs_is_setup) then
     blacs_is_setup = .False.
     call blacs_gridexit(blacs_ctxt)
   end if

end subroutine

!> 
!! Computes global row index based on local row index
!!
subroutine elsi_get_global_row (global_idx, local_idx)

   implicit none

   integer, intent(out) :: global_idx  !< Global index 
   integer, intent(in)  :: local_idx   !< Local index

   integer :: block !< local block 
   integer :: idx   !< local index in block


   block = FLOOR( 1d0 * (local_idx - 1) / n_b_rows) 
   idx = local_idx - block * n_b_rows

   global_idx = my_p_row * n_b_rows + block * n_b_rows * n_p_rows + idx
  
end subroutine


!> 
!! Computes global column index based on local column index
!!
subroutine elsi_get_global_col (global_idx, local_idx)

   implicit none

   integer, intent(out) :: global_idx   !< global index
   integer, intent(in)  :: local_idx    !< local index

   integer :: block !< local block
   integer :: idx   !< local index in block


   block = FLOOR( 1d0 * (local_idx - 1) / n_b_cols) 
   idx = local_idx - block * n_b_cols

   global_idx = my_p_col * n_b_cols + block * n_p_cols * n_b_cols + idx
  
end subroutine

!> 
!! Gets global matrix specifications such as rank and block dimensions
!!
subroutine elsi_get_global_dimensions (g_dim, g_block_rows, g_block_cols)

   implicit none

   integer, intent(out) :: g_dim   !< global rank
   integer, intent(out) :: g_block_rows !< global blocksize
   integer, intent(out) :: g_block_cols !< global blocksize

   g_dim   = n_g_rank 
   g_block_rows = n_b_rows
   g_block_cols = n_b_cols

end subroutine

!> 
!! Gets the local matrix dimensions
!!
subroutine elsi_get_local_dimensions (rows, cols)

   implicit none

   integer, intent(out) :: rows   !< output local rows
   integer, intent(out) :: cols   !< output local cols

   rows  = n_l_rows 
   cols  = n_l_cols

end subroutine

!> 
!! Get the processor grid as specified by BLACS
!!
subroutine elsi_get_processor_grid (rows, cols)

   implicit none

   integer, intent(out) :: rows   !< output processor rows
   integer, intent(out) :: cols   !< output processor cols

   rows  = n_p_rows 
   cols  = n_p_cols

end subroutine


!> 
!! Gets the row and column communicators for ELPA
!!
subroutine elsi_get_comm_grid (rows, cols)

   implicit none

   integer, intent(out) :: rows   !< output row communicator
   integer, intent(out) :: cols   !< output col communicator

   rows  = mpi_comm_row 
   cols  = mpi_comm_col

end subroutine

!> 
!! Gets the process id 
!!

subroutine elsi_get_myid (id)

   implicit none

   integer, intent(out) :: id   !< output process id

   id = myid 

end subroutine

!> 
!! Debug printout of a matrix
!!

subroutine elsi_matrix_print (matrix, n_rows, n_cols, matrixname)

   implicit none

   integer, intent(in) :: n_rows
   integer, intent(in) :: n_cols
   real*8, intent(in) :: matrix(n_rows,n_cols) 
   character(len=*), intent(in) :: matrixname 

   integer :: id, i_row, i_col
   
   if (myid == 0) print *, trim(matrixname)

   do id = 0, n_procs - 1
     if (myid == id) then
        do i_row = 1, n_rows
           do i_col = 1, n_cols
               write(*,"(A,I6,A,I6,A,I6,F13.6)") " Process ", id, &
               " Row ", i_row, " Col ", i_col, " value", matrix(i_row,i_col)
           end do
        end do
     end if
     call MPI_Barrier(mpi_comm_global, mpierr)
   end do 

end subroutine

!> 
!! Debug printout of a vector
!!

subroutine elsi_int_vector_print (vector, n_dim, vectorname)

   implicit none

   integer, intent(in) :: n_dim
   integer, intent(in) :: vector(n_dim) 
   character(len=*), intent(in) :: vectorname 

   integer :: id, i_dim

   if (myid == 0) print *, trim(vectorname)

   do id = 0, n_procs - 1
     if (myid == id) then
        do i_dim = 1, n_dim
            write(*,"(A,I6,A,I10,A,I13)") " Process ", id, " Entry ", i_dim,&
           " value", vector(i_dim)
        end do
     end if
     call MPI_Barrier(mpi_comm_global, mpierr)
   end do 

end subroutine

subroutine elsi_int_value_print (val, valname)

   implicit none

   integer :: val
   character(len=*), intent(in) :: valname 

   integer :: id

   if (myid == 0) print *, trim(valname)

   do id = 0, n_procs - 1
     if (myid == id) then
        write(*,"(A,I6,A,I13)") " Process ", id, " value", val
     end if
     call MPI_Barrier(mpi_comm_global, mpierr)
   end do 

end subroutine

subroutine elsi_real_value_print (val, valname)

   implicit none

   real :: val
   character(len=*), intent(in) :: valname 

   integer :: id

   if (myid == 0) print *, trim(valname)

   do id = 0, n_procs - 1
     if (myid == id) then
       write(*,"(A,I4,A,F13.6)") " Process ", id, " value", val
     end if
     call MPI_Barrier(mpi_comm_global, mpierr)
   end do 

end subroutine



subroutine elsi_real_vector_print (vector, n_dim, vectorname)

   implicit none

   integer, intent(in) :: n_dim
   real*8, intent(in) :: vector(n_dim) 
   character(len=*), intent(in) :: vectorname 

   integer :: id, i_dim

   if (myid == 0) print *, trim(vectorname)

   do id = 0, n_procs - 1
     if (myid == id) then
        do i_dim = 1, n_dim
           write(*,"(A,I6,A,I10,A,F13.6)") " Process ", id, " Entry ", i_dim,&
           " value", vector(i_dim)
        end do
     end if
     call MPI_Barrier(mpi_comm_global, mpierr)
   end do 

end subroutine


!> 
!! Debug printout of a vector
!!

subroutine elsi_statement_print (message)

   implicit none

   character(len=*), intent(in) :: message 

   if (myid == 0) print *, trim(message)
   call MPI_Barrier(mpi_comm_global, mpierr)

end subroutine

!>
!! Elsi allocation routines with error handling
!!

subroutine elsi_allocate_real_vector (vector, n_elements, vectorname, caller)

   implicit none

   real*8, allocatable, intent(inout) :: vector(:)
   integer, intent(in) :: n_elements
   character(len=*), intent(in) :: vectorname
   character(len=*), intent(in) :: caller

   integer :: error

   allocate(vector(n_elements), stat = error)

   if (error > 0) call elsi_stop ("Insufficient memory to allocate " &
         // trim(vectorname), caller)

   vector = 0d0 

end subroutine

subroutine elsi_allocate_int_vector (vector, n_elements, vectorname, caller)

   implicit none

   integer, allocatable, intent(inout) :: vector(:)
   integer, intent(in) :: n_elements
   character(len=*), intent(in) :: vectorname
   character(len=*), intent(in) :: caller

   integer :: error

   allocate(vector(n_elements), stat = error)

   if (error > 0) call elsi_stop ("Insufficient memory to allocate " &
         // trim(vectorname), caller)

   vector = 0

end subroutine

subroutine elsi_allocate_complex_vector (vector, n_elements, vectorname, caller)

   implicit none

   complex*16, allocatable, intent(inout) :: vector(:)
   integer, intent(in) :: n_elements
   character(len=*), intent(in) :: vectorname
   character(len=*), intent(in) :: caller

   integer :: error

   allocate(vector(n_elements), stat = error)

   if (error > 0) call elsi_stop ("Insufficient memory to allocate " &
         // trim(vectorname), caller)

   vector = CMPLX(0d0,0d0)

end subroutine

subroutine elsi_allocate_real_matrix (matrix, n_rows, n_cols, matrixname,&
      caller)

   implicit none

   real*8, allocatable, intent(inout) :: matrix(:,:)
   integer, intent(in) :: n_rows
   integer, intent(in) :: n_cols
   character(len=*), intent(in) :: matrixname
   character(len=*), intent(in) :: caller

   integer :: error

   allocate(matrix(n_rows,n_cols), stat = error)

   if (error > 0) call elsi_stop ("Insufficient memory to allocate " &
         // trim(matrixname), caller)

   matrix = 0d0 

end subroutine

subroutine elsi_allocate_int_matrix (matrix, n_rows, n_cols, matrixname,&
      caller)

   implicit none

   integer, allocatable, intent(inout) :: matrix(:,:)
   integer, intent(in) :: n_rows
   integer, intent(in) :: n_cols
   character(len=*), intent(in) :: matrixname
   character(len=*), intent(in) :: caller

   integer :: error

   allocate(matrix(n_rows,n_cols), stat = error)

   if (error > 0) call elsi_stop ("Insufficient memory to allocate " &
         // trim(matrixname), caller)

   matrix = 0 

end subroutine

subroutine elsi_allocate_complex_matrix (matrix, n_rows, n_cols, matrixname,&
      caller)

   implicit none

   complex*16, allocatable, intent(inout) :: matrix(:,:)
   integer, intent(in) :: n_rows
   integer, intent(in) :: n_cols
   character(len=*), intent(in) :: matrixname
   character(len=*), intent(in) :: caller

   integer :: error

   allocate(matrix(n_rows,n_cols), stat = error)

   if (error > 0) call elsi_stop ("Insufficient memory to allocate " &
         // trim(matrixname), caller)

   matrix = CMPLX(0d0,0d0) 

end subroutine

end module ELSI_MPI_TOOLS
