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

!> ELSI Matrix conversion
!! This module converts dense matrices to sparse matrices and vice versa
!! and solving or circumventing an eigenvalue problem using ELPA, OMM, or PEXSI
!!

module ELSI_MATRIX_CONVERSION

  use iso_c_binding
  use ELSI_DIMENSIONS

  implicit none
  private

  !< Internal Storage
  real*8 :: threshhold = 1.0d-12  !< threshhold to define numerical zero

  public :: elsi_compute_N_nonzero
  public :: elsi_get_local_N_nonzero
  public :: elsi_get_n_nonzero_column
  public :: elsi_dense_to_ccs
  public :: elsi_dense_to_ccs_by_pattern
  public :: elsi_ccs_to_dense

  contains

!>
!!  This routine computes the number of non_zero elements column-wise for 
!!  the full matrix and return it as a vector of dimension n_g_rank 
subroutine elsi_get_n_nonzero_column(matrix, n_rows, n_cols, l_offset,n_nonzero)

   implicit none
   include 'mpif.h'

   real*8,   intent(in)  :: matrix(n_rows,n_cols)  !< local matrix
   integer,  intent(in)  :: n_rows !< local rows
   integer,  intent(in)  :: n_cols !< local cols
   integer,  intent(in)  :: l_offset !< offset in global array
   integer,  intent(out) :: n_nonzero (n_g_rank)
   
   integer, allocatable :: buffer(:)  !< buffer for MPI communication
   integer             :: i_col       !< col counter
   integer             :: nnz_in_col  !< non-zero element counter in column
   integer             :: sent_n_elements  !< number of elements to sent in MPI call

   integer, parameter :: max_len = 1e6

   n_nonzero = 0

   do i_col = 1, n_cols
     call elsi_get_local_N_nonzero(matrix(:,i_col), n_rows, 1, nnz_in_col)
     n_nonzero(l_offset + i_col - 1) = nnz_in_col
   end do

   allocate(buffer (min(n_g_rank,max_len)))

   do i_col = 1,n_g_rank,max_len

     if (n_g_rank - i_col + 1 < max_len) then
        sent_n_elements = n_g_rank - i_col + 1
     else
        sent_n_elements = max_len
     end if

     call MPI_AllReduce(n_nonzero(i_col), buffer, sent_n_elements, MPI_INTEGER,&
           MPI_SUM, mpi_comm_global, mpierr)
     n_nonzero(i_col:i_col-1+sent_n_elements) = buffer(1:sent_n_elements)
   
   end do

   deallocate(buffer)


end subroutine


!>
!!  This routine computes the number of non_zero elements for a local matrix
!!
subroutine elsi_get_local_N_nonzero(matrix, n_rows, n_cols, n_nonzero)

   implicit none
   include 'mpif.h'

   real*8,  intent(in)  :: matrix(n_rows,n_cols)  !< local matrix
   integer,  intent(in) :: n_rows !< local rows
   integer,  intent(in) :: n_cols !< local cols
   integer, intent(out) :: n_nonzero
   
   integer             :: i_row       !< row counter
   integer             :: i_col       !< col counter
   integer             :: i_nonzero   !< non-zero element counter

   n_nonzero = 0
   do i_col = 1, n_cols
     do i_row = 1, n_rows
       if (abs(matrix(i_row,i_col)) > threshhold) then
         n_nonzero = n_nonzero + 1
       end if
     end do
   end do

end subroutine

!>
!!  This routine computes the number of non_zero elements
!!
subroutine elsi_compute_N_nonzero(matrix,n_rows,n_cols)

   implicit none
   include 'mpif.h'

   integer, intent(in) :: n_rows !< number of rows
   integer, intent(in) :: n_cols !< number of columns
   real*8,  intent(in) :: matrix(n_rows,n_cols) !< local matrix

   ! Set the number of non-zero elements in the local matrix
   call elsi_get_local_N_nonzero(matrix,n_rows,n_cols,n_l_nonzero)

   ! Set the number of non_zero elements in the global matrix
   call MPI_ALLREDUCE(n_l_nonzero, n_g_nonzero, 1, mpi_integer, mpi_sum,&
         mpi_comm_global, mpierr)

end subroutine

!>
!!  This routine transforms a dense matrix to the 
!!  Compressed Column Storage (CCS) Format
!!
subroutine elsi_dense_to_ccs(matrix, n_rows, n_cols, val, nnz, row_ind, col_ptr)

   implicit none

   integer, intent(in) :: n_rows !< number of rows
   integer, intent(in) :: n_cols !< number of columns
   integer, intent(in) :: nnz    !< number of non zero elements
   real*8,  intent(in) :: matrix(n_rows,n_cols) !< local matrix
   real*8,  intent(out):: val(nnz)          !< values
   integer, intent(out):: row_ind(nnz)      !< row index
   integer, intent(out):: col_ptr(n_cols + 1)     !< column pointer

   integer             :: i_row       !< row counter
   integer             :: i_col       !< col counter
   integer             :: i_val       !< value counter
   logical             :: first       !< First encounter in col?

   integer             :: indexshift

   i_val = 0
   col_ptr = 0
   do i_col = 1, n_cols
     first = .true.
     do i_row = 1, n_rows
       if (abs(matrix(i_row,i_col)) > threshhold) then
         i_val = i_val + 1
         if (first) then
           col_ptr(i_col) = i_val
           first = .false.
         end if
         val(i_val) = matrix(i_row,i_col)
         row_ind(i_val) = i_row
       end if
     end do
   end do

   col_ptr(n_cols + 1) = i_val + 1

   if (i_val /= nnz) then
      call elsi_stop("Number of non zero elements differ",&
         "elsi_dense_to_ccs")
   end if


end subroutine

!>
!!  This routine transforms a dense matrix to the 
!!  Compressed Column Storage (CCS) Format based 
!!  on a given col_ind and row_ptr
!!
subroutine elsi_dense_to_ccs_by_pattern(matrix, n_rows, n_cols, val, nnz, row_ind, col_ptr)

   implicit none

   integer, intent(in) :: n_rows !< number of rows
   integer, intent(in) :: n_cols !< number of columns
   integer, intent(in) :: nnz    !< number of non zero elements
   real*8,  intent(in) :: matrix(n_rows,n_cols) !< local matrix
   real*8,  intent(out):: val(nnz)          !< values
   integer, intent(in) :: row_ind(nnz)      !< row index
   integer, intent(in) :: col_ptr(n_cols + 1)     !< column pointer

   integer             :: i_row       !< row counter
   integer             :: i_col       !< col counter
   integer             :: i_val       !< value counter

   i_col = 0
   do i_val = 1, nnz
     if (i_val == col_ptr(i_col + 1) .and. i_col /= n_cols) then
       i_col = i_col + 1
     end if
     i_row = row_ind(i_val)
     val(i_val) = matrix(i_row,i_col)
   end do


end subroutine

!>
!!  This routine transforms a CCS matrix back to dense form  
!!  on a given col_ind and row_ptr
!!
subroutine elsi_ccs_to_dense(matrix, n_rows, n_cols, val, nnz, row_ind, col_ptr)

   implicit none

   integer, intent(in) :: n_rows !< number of rows
   integer, intent(in) :: n_cols !< number of columns
   integer, intent(in) :: nnz    !< number of non zero elements
   real*8,  intent(out) :: matrix(n_rows,n_cols) !< local matrix
   real*8,  intent(in)  :: val(nnz)          !< values
   integer, intent(in)  :: row_ind(nnz)      !< row index
   integer, intent(in)  :: col_ptr(n_cols + 1)     !< column pointer

   integer             :: i_row       !< row counter
   integer             :: i_col       !< col counter
   integer             :: i_val       !< value counter

   matrix = 0d0
   i_col = 0
   do i_val = 1, nnz
     if (i_val == col_ptr(i_col + 1) .and. i_col /= n_cols) then
       i_col = i_col + 1
     end if
     i_row = row_ind(i_val)
     matrix(i_row,i_col) = val(i_val)
   end do

end subroutine

end module ELSI_MATRIX_CONVERSION
