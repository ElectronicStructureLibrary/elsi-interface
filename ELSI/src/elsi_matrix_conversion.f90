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
  public :: elsi_dense_to_crs
  public :: elsi_dense_to_ccs
  public :: elsi_dense_to_ccs_by_pattern
  public :: elsi_ccs_to_dense

  contains

!>
!!  This routine computes the number of non_zero elements
!!
subroutine elsi_compute_N_nonzero(matrix)

   implicit none
   include 'mpif.h'

   real*8,  intent(in) :: matrix(:,:) !< local matrix
   
   integer             :: i_row       !< row counter
   integer             :: i_col       !< col counter
   integer             :: i_nonzero   !< non-zero element counter

   i_nonzero = 0
   do i_row = 1, n_l_rows
     do i_col = 1, n_l_cols
       if (abs(matrix(i_row,i_col)) > threshhold) then
         i_nonzero = i_nonzero + 1
       end if
     end do
   end do

   ! Set the number of non-zero elements in the local matrix
   n_l_nonzero = i_nonzero

   ! Set the number of non_zero elements in the global matrix
   call MPI_ALLREDUCE(i_nonzero, n_g_nonzero, 1, mpi_integer, mpi_sum,&
         mpi_comm_global, mpierr)

end subroutine

!>
!!  This routine transforms a dense matrix to the 
!!  Compressed Row Storage (CRS) Format
!!
subroutine elsi_dense_to_crs(matrix, val, col_ind, row_ptr)

   implicit none

   real*8,  intent(in) :: matrix(n_l_rows,n_l_cols) !< local matrix
   real*8,  intent(out):: val(n_l_nonzero)          !< values
   integer, intent(out):: col_ind(n_l_nonzero)      !< column index
   integer, intent(out):: row_ptr(n_l_rows + 1)     !< row pointer

   integer             :: i_row       !< row counter
   integer             :: i_col       !< col counter
   integer             :: i_val       !< value counter
   logical             :: first       !< First encounter in row?

   i_val = 0
   row_ptr = 0
   do i_row = 1, n_l_rows
     first = .true.
     do i_col = 1, n_l_cols
       if (abs(matrix(i_row,i_col)) > threshhold) then
         i_val = i_val + 1
         if (first) then
           row_ptr(i_row) = i_val
           first = .false.
         end if
         val(i_val) = matrix(i_row,i_col)
         col_ind(i_val) = i_col
       end if
     end do
   end do

   row_ptr(n_l_rows + 1) = i_val

   if (i_val /= n_l_nonzero) then
      call elsi_stop("Number of non zero elements differ",&
         "elsi_sparse_to_crs")
   end if

end subroutine

!>
!!  This routine transforms a dense matrix to the 
!!  Compressed Column Storage (CCS) Format
!!
subroutine elsi_dense_to_ccs(matrix, val, row_ind, col_ptr)

   implicit none

   real*8,  intent(in) :: matrix(n_l_rows,n_l_cols) !< local matrix
   real*8,  intent(out):: val(n_l_nonzero)          !< values
   integer, intent(out):: row_ind(n_l_nonzero)      !< row index
   integer, intent(out):: col_ptr(n_l_cols + 1)     !< column pointer

   integer             :: i_row       !< row counter
   integer             :: i_col       !< col counter
   integer             :: i_val       !< value counter
   logical             :: first       !< First encounter in col?

   integer             :: indexshift

   i_val = 0
   col_ptr = 0
   do i_col = 1, n_l_cols
     first = .true.
     do i_row = 1, n_l_rows
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

   col_ptr(n_l_cols + 1) = i_val + 1

   if (i_val /= n_l_nonzero) then
      call elsi_stop("Number of non zero elements differ",&
         "elsi_sparse_to_ccs")
   end if

end subroutine

!>
!!  This routine transforms a dense matrix to the 
!!  Compressed Column Storage (CCS) Format based 
!!  on a given col_ind and row_ptr
!!
subroutine elsi_dense_to_ccs_by_pattern(matrix, val, row_ind, col_ptr)

   implicit none

   real*8,  intent(in) :: matrix(n_l_rows,n_l_cols) !< local matrix
   real*8,  intent(out):: val(n_l_nonzero)          !< values
   integer, intent(in) :: row_ind(n_l_nonzero)      !< row index
   integer, intent(in) :: col_ptr(n_l_cols + 1)     !< column pointer

   integer             :: i_row       !< row counter
   integer             :: i_col       !< col counter
   integer             :: i_val       !< value counter

   integer             :: indexshift

   if (method == PEXSI) then
      indexshift = 1
   else 
      indexshift = 0
   end if

   i_col = 0
   do i_val = 1, n_l_nonzero
     if (i_val == col_ptr(i_col + 1) .and. i_col /= n_l_cols) then
       i_col = i_col + 1
     end if
     i_row = row_ind(i_val)
     val(i_val) = matrix(i_row,i_col)
   end do


end subroutine

!>
!!  This routine transforms a dense matrix to the 
!!  Compressed Column Storage (CCS) Format based 
!!  on a given col_ind and row_ptr
!!
subroutine elsi_ccs_to_dense(matrix, val, row_ind, col_ptr)

   implicit none

   real*8,  intent(out) :: matrix(n_l_rows,n_l_cols) !< local matrix
   real*8,  intent(in)  :: val(n_l_nonzero)          !< values
   integer, intent(in)  :: row_ind(n_l_nonzero)      !< row index
   integer, intent(in)  :: col_ptr(n_l_cols + 1)     !< column pointer

   integer             :: i_row       !< row counter
   integer             :: i_col       !< col counter
   integer             :: i_val       !< value counter

   i_col = 0
   do i_val = 1, n_l_nonzero
     if (i_val == col_ptr(i_col + 1) .and. i_col /= n_l_cols) then
       i_col = i_col + 1
     end if
     i_row = row_ind(i_val)
     matrix(i_row,i_col) = val(i_val)
   end do

end subroutine

end module ELSI_MATRIX_CONVERSION
