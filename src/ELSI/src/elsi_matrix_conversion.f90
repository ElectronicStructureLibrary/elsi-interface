!Copyright (c) 2015-2017, ELSI consortium
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
!! This module converts matrix storage format.
!!

module ELSI_MATRIX_CONVERSION

   use iso_c_binding
   use ELSI_DIMENSIONS
   use ELSI_UTILS

   implicit none
   private

   public :: elsi_get_global_col
   public :: elsi_get_global_row
   public :: elsi_get_global_nnz
   public :: elsi_get_local_nnz
   public :: elsi_dense_to_ccs
   public :: elsi_dense_to_ccs_by_pattern
   public :: elsi_ccs_to_dense

contains

!> 
!! This routine computes global row index based on local row index.
!!
subroutine elsi_get_global_row(global_idx,local_idx)

   implicit none

   integer, intent(in)  :: local_idx  !< Local index
   integer, intent(out) :: global_idx !< Global index

   integer :: block !< Local block
   integer :: idx !< Local index in block

   block = FLOOR(1d0*(local_idx-1)/n_b_rows)
   idx = local_idx-block*n_b_rows

   global_idx = my_p_row*n_b_rows+block*n_b_rows*n_p_rows+idx

end subroutine

!> 
!! This routine computes global column index based on local column index.
!!
subroutine elsi_get_global_col(global_idx,local_idx)

   implicit none

   integer, intent(in)  :: local_idx  !< Local index
   integer, intent(out) :: global_idx !< Global index

   integer :: block !< Local block
   integer :: idx !< Local index in block

   block = FLOOR(1d0*(local_idx-1)/n_b_cols)
   idx = local_idx-block*n_b_cols

   global_idx = my_p_col*n_b_cols+block*n_b_cols*n_p_cols+idx

end subroutine

!>
!! This routine computes the (local) number of non_zero elements.
!!
subroutine elsi_get_local_nnz(matrix,n_rows,n_cols,nnz)

   implicit none

   real*8,  intent(in)  :: matrix(n_rows,n_cols) !< Local matrix
   integer, intent(in)  :: n_rows                !< Local rows
   integer, intent(in)  :: n_cols                !< Local cols
   integer, intent(out) :: nnz                   !< Number of non-zero elements

   integer :: i_row !< Row counter
   integer :: i_col !< Column counter
   integer :: i_nz !< Non-zero element counter

   nnz = 0

   do i_col = 1,n_cols
      do i_row = 1,n_rows
         if(ABS(matrix(i_row,i_col)) > zero_threshold) then
            nnz = nnz+1
         endif
      enddo
   enddo

end subroutine

!>
!! This routine computes the (global) number of non_zero elements.
!!
subroutine elsi_get_global_nnz(matrix, n_rows, n_cols)

   implicit none
   include "mpif.h"

   real*8,  intent(in) :: matrix(n_rows,n_cols) !< Local matrix
   integer, intent(in) :: n_rows                !< Number of rows
   integer, intent(in) :: n_cols                !< Number of columns

   ! Set the number of non-zero elements in the local matrix
   call elsi_get_local_nnz(matrix,n_rows,n_cols,nnz_l)

   ! Set the number of non_zero elements in the global matrix
   call MPI_Allreduce(nnz_l,nnz_g,1,mpi_integer,mpi_sum,&
                      mpi_comm_global,mpierr)

end subroutine

!>
!! This routine transforms a (local) dense matrix to sparse CCS format.
!!
subroutine elsi_dense_to_ccs(matrix,n_rows,n_cols,nnz,val,row_ind,col_ptr)

   implicit none

   real*8,  intent(in)  :: matrix(n_rows,n_cols) !< Local matrix
   integer, intent(in)  :: n_rows                !< Number of rows
   integer, intent(in)  :: n_cols                !< Number of columns
   integer, intent(in)  :: nnz                   !< Number of non-zero elements
   real*8,  intent(out) :: val(nnz)              !< Values
   integer, intent(out) :: row_ind(nnz)          !< Row index
   integer, intent(out) :: col_ptr(n_cols+1)     !< Column pointer

   integer :: i_row !< Row counter
   integer :: i_col !< Column counter
   integer :: i_val !< Value counter
   logical :: new_col !< Enter a new column?

   character*40, parameter :: caller = "elsi_dense_to_ccs"

   i_val = 0
   col_ptr = 0
   do i_col = 1,n_cols
      new_col = .true.
      do i_row = 1,n_rows
         if(ABS(matrix(i_row,i_col)) > zero_threshold) then
            i_val = i_val+1
            if(new_col) then
               col_ptr(i_col) = i_val
               new_col = .false.
            endif
            val(i_val) = matrix(i_row,i_col)
            row_ind(i_val) = i_row
         endif
      enddo
   enddo

   col_ptr(n_cols+1) = i_val+1

   if(i_val /= nnz) then
      call elsi_stop(" Number of nonzero differs. Exiting...", caller)
   endif

end subroutine

!>
!! This routine transforms a (local) dense matrix to sparse CCS
!! format based on given col_ind and row_ptr.
!!
subroutine elsi_dense_to_ccs_by_pattern(matrix,n_rows,n_cols,&
                                        nnz,row_ind,col_ptr,val)

   implicit none

   real*8,  intent(in)  :: matrix(n_rows,n_cols) !< Local matrix
   integer, intent(in)  :: n_rows                !< Number of rows
   integer, intent(in)  :: n_cols                !< Number of columns
   integer, intent(in)  :: nnz                   !< Number of non-zero elements
   integer, intent(in)  :: row_ind(nnz)          !< Row index
   integer, intent(in)  :: col_ptr(n_cols+1)     !< Column pointer
   real*8,  intent(out) :: val(nnz)              !< Values

   integer :: i_row !< Row counter
   integer :: i_col !< Column counter
   integer :: i_val !< Value counter

   i_col = 0
   do i_val = 1, nnz
      if(i_val == col_ptr(i_col+1) .and. i_col /= n_cols) then
         i_col = i_col+1
      endif
      i_row = row_ind(i_val)
      val(i_val) = matrix(i_row,i_col)
   enddo

end subroutine

!>
!! This routine transforms a (local) sparse CCS matrix to dense
!! format based on given col_ind and row_ptr.
!!
subroutine elsi_ccs_to_dense(matrix,n_rows,n_cols,val,nnz,row_ind,col_ptr)

   implicit none

   integer, intent(in)  :: n_rows                !< Number of rows
   integer, intent(in)  :: n_cols                !< Number of columns
   integer, intent(in)  :: nnz                   !< Number of non-zero elements
   real*8,  intent(out) :: matrix(n_rows,n_cols) !< Local matrix
   real*8,  intent(in)  :: val(nnz)              !< Values
   integer, intent(in)  :: row_ind(nnz)          !< Row index
   integer, intent(in)  :: col_ptr(n_cols+1)     !< Column pointer

   integer :: i_row !< Row counter
   integer :: i_col !< Col counter
   integer :: i_val !< Value counter

   matrix = 0d0
   i_col = 0
   do i_val = 1,nnz
      if(i_val == col_ptr(i_col+1) .and. i_col /= n_cols) then
         i_col = i_col+1
      endif
      i_row = row_ind(i_val)
      matrix(i_row,i_col) = val(i_val)
   enddo

end subroutine

end module ELSI_MATRIX_CONVERSION
