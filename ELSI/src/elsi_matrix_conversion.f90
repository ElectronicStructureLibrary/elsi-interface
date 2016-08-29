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

!> ELSI Matrix conversion
!! This module converts matrix storage format.
!!

module ELSI_MATRIX_CONVERSION

   use iso_c_binding
   use ELSI_DIMENSIONS
   use ELSI_MPI_TOOLS

   implicit none
   private

   real*8, parameter :: threshold = 1.0d-12 !< Threshold to define numerical zero.

   public :: elsi_get_global_col
   public :: elsi_get_global_row
   public :: elsi_get_global_n_nonzero
   public :: elsi_get_local_n_nonzero
   public :: elsi_get_n_nonzero_column
   public :: elsi_dense_to_ccs
   public :: elsi_dense_to_ccs_by_pattern
   public :: elsi_ccs_to_dense
   public :: elsi_2dbc_to_1db

contains

!> 
!! Computes global row index based on local row index.
!!
subroutine elsi_get_global_row(global_idx, local_idx)

   implicit none

   integer, intent(in) :: local_idx !< Local index
   integer, intent(out) :: global_idx !< Global index

   integer :: block !< Local block
   integer :: idx !< Local index in block

   block = FLOOR(1d0*(local_idx-1)/n_b_rows)
   idx = local_idx-block*n_b_rows

   global_idx = my_p_row*n_b_rows+block*n_b_rows*n_p_rows+idx

end subroutine

!> 
!! Computes global column index based on local column index.
!!
subroutine elsi_get_global_col(global_idx, local_idx)

   implicit none

   integer, intent(in) :: local_idx !< Local index
   integer, intent(out) :: global_idx !< Global index

   integer :: block !< Local block
   integer :: idx !< Local index in block

   block = FLOOR(1d0*(local_idx-1)/n_b_cols)
   idx = local_idx-block*n_b_cols

   global_idx = my_p_col*n_b_cols+block*n_b_cols*n_p_cols+idx

end subroutine

!>
!!  This routine computes the number of non_zero elements column-wise for
!!  the full matrix and returns a vector of dimension n_g_rank.
!!
subroutine elsi_get_n_nonzero_column(matrix, n_rows, n_cols, l_offset, n_nonzero)

   implicit none
   include 'mpif.h'

   real*8, intent(in) :: matrix(n_rows,n_cols) !< Local matrix
   integer, intent(in) :: n_rows !< Local rows
   integer, intent(in) :: n_cols !< Local cols
   integer, intent(in) :: l_offset !< Offset in global array
   integer, intent(out) :: n_nonzero(n_g_rank)

   integer, allocatable :: buffer(:) !< MPI buffer
   integer :: i_col !< Column counter
   integer :: nnz_in_col !< Non-zero element in column counter
   integer :: send_n_elements !< Number of elements to send by MPI
   integer, parameter :: max_len = 1e6

   n_nonzero = 0

   do i_col = 1, n_cols
      call elsi_get_local_n_nonzero(matrix(:,i_col), n_rows, 1, nnz_in_col)
      n_nonzero(l_offset+i_col-1) = nnz_in_col
   enddo

   allocate(buffer(min(n_g_rank,max_len)))

   do i_col = 1,n_g_rank,max_len
      if(n_g_rank-i_col+1 < max_len) then
         send_n_elements = n_g_rank-i_col+1
      else
         send_n_elements = max_len
      endif

      call MPI_AllReduce(n_nonzero(i_col), buffer, send_n_elements, MPI_INTEGER, &
                         MPI_SUM, mpi_comm_global, mpierr)
      n_nonzero(i_col:i_col-1+send_n_elements) = buffer(1:send_n_elements)
   enddo

   deallocate(buffer)

end subroutine

!>
!!  This routine computes the (local) number of non_zero elements.
!!
subroutine elsi_get_local_n_nonzero(matrix, n_rows, n_cols, n_nonzero)

   implicit none
   include 'mpif.h'

   real*8, intent(in) :: matrix(n_rows,n_cols) !< Local matrix
   integer, intent(in) :: n_rows !< Local rows
   integer, intent(in) :: n_cols !< Local cols
   integer, intent(out) :: n_nonzero !< Number of non-zero elements

   integer :: i_row !< Row counter
   integer :: i_col !< Column counter
   integer :: i_nonzero !< Non-zero element counter

   n_nonzero = 0

   do i_col = 1, n_cols
      do i_row = 1, n_rows
         if(abs(matrix(i_row,i_col)) > threshold) then
            n_nonzero = n_nonzero+1
         endif
      enddo
   enddo

end subroutine

!>
!!  This routine computes the (global) number of non_zero elements.
!!
subroutine elsi_get_global_n_nonzero(matrix, n_rows, n_cols)

   implicit none
   include 'mpif.h'

   real*8, intent(in) :: matrix(n_rows,n_cols) !< Local matrix
   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of columns

   ! Set the number of non-zero elements in the local matrix
   call elsi_get_local_n_nonzero(matrix, n_rows, n_cols, n_l_nonzero)

   ! Set the number of non_zero elements in the global matrix
   call MPI_ALLREDUCE(n_l_nonzero, n_g_nonzero, 1, mpi_integer, mpi_sum, &
                      mpi_comm_global, mpierr)

end subroutine

!>
!!  This routine transforms a (local) dense matrix to sparse CCS format.
!!
subroutine elsi_dense_to_ccs(matrix, n_rows, n_cols, nnz, val, row_ind, col_ptr)

   implicit none

   real*8, intent(in) :: matrix(n_rows,n_cols) !< Local matrix
   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of columns
   integer, intent(in) :: nnz !< Number of non-zero elements
   real*8, intent(out) :: val(nnz) !< Values
   integer, intent(out) :: row_ind(nnz) !< Row index
   integer, intent(out) :: col_ptr(n_cols+1) !< Column pointer

   integer :: i_row !< Row counter
   integer :: i_col !< Column counter
   integer :: i_val !< Value counter
   logical :: new_col !< Enter a new column?

   character*40, parameter :: caller = "elsi_dense_to_ccs"

   i_val = 0
   col_ptr = 0
   do i_col = 1, n_cols
      new_col = .true.
      do i_row = 1, n_rows
         if(abs(matrix(i_row,i_col)) > threshold) then
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
      call elsi_stop(" Number of nonzero differs. Exiting... ", caller)
   endif

end subroutine

!>
!!  This routine transforms a (local) dense matrix to sparse CCS
!!  format based on given col_ind and row_ptr.
!!
subroutine elsi_dense_to_ccs_by_pattern(matrix, n_rows, n_cols, &
                                        nnz, row_ind, col_ptr, val)

   implicit none

   real*8, intent(in) :: matrix(n_rows,n_cols) !< Local matrix
   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of columns
   integer, intent(in) :: nnz !< Number of non-zero elements
   integer, intent(in) :: row_ind(nnz) !< Row index
   integer, intent(in) :: col_ptr(n_cols+1) !< Column pointer
   real*8, intent(out) :: val(nnz) !< Values

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
!!  This routine transforms a (local) sparse CCS matrix to dense
!!  format based on given col_ind and row_ptr.
!!
subroutine elsi_ccs_to_dense(matrix, n_rows, n_cols, val, nnz, row_ind, col_ptr)

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of columns
   integer, intent(in) :: nnz !< Number of non-zero elements
   real*8, intent(out) :: matrix(n_rows,n_cols) !< Local matrix
   real*8, intent(in) :: val(nnz) !< Values
   integer, intent(in) :: row_ind(nnz) !< Row index
   integer, intent(in) :: col_ptr(n_cols+1) !< Column pointer

   integer :: i_row !< Row counter
   integer :: i_col !< Col counter
   integer :: i_val !< Value counter

   matrix = 0d0
   i_col = 0
   do i_val = 1, nnz
      if(i_val == col_ptr(i_col+1) .and. i_col /= n_cols) then
         i_col = i_col+1
      endif
      i_row = row_ind(i_val)
      matrix(i_row,i_col) = val(i_val)
   enddo

end subroutine

!> 
!!  This routine converts the distribution of a dense matrix from
!!  2D block-cyclic to 1D block. The converted 1D dense matrix can
!!  be easily converted to 1D block distributed sparse CCS format
!!  used in PEXSI by "elsi_bc_to_ccs".
!!
subroutine elsi_2dbc_to_1db(matrix_in, matrix_out)

   implicit none
   include 'mpif.h'

   real*8, intent(in) :: matrix_in(n_l_rows,n_l_cols)
   real*8, intent(out) :: matrix_out(n_l_rows_pexsi, n_l_cols_pexsi)

   integer :: i_row !< Row counter
   integer :: i_col !< Col counter
   integer :: i_val !< Value counter
   integer :: i_proc !< Process counter
   integer :: global_col_id !< Global column id
   integer :: global_row_id !< Global row id
   integer :: local_col_id !< Local column id in 1D block distribution
   integer :: local_row_id !< Local row id in 1D block distribution

   integer :: dest(n_l_rows*n_l_cols) !< Destination of each element
   integer :: global_pos(n_l_rows*n_l_cols) !< Global 1d id of each element

   ! For the meaning of each array here, see documentation of MPI_Alltoallv
   real*8 :: val_send_buffer(n_l_rows*n_l_cols) !< Send buffer for value
   integer :: pos_send_buffer(n_l_rows*n_l_cols) !< Send buffer for global 1D id
   integer :: send_count(n_procs) !< Number of elements to send to each processor
   integer :: send_displ(n_procs) !< Displacement from which to take the outgoing data
   integer :: send_displ_aux !< Auxiliary variable used to set displacement

   real*8 :: val_recv_buffer(n_l_rows_pexsi*n_l_cols_pexsi) !< Receive buffer for value
   integer :: pos_recv_buffer(n_l_rows_pexsi*n_l_cols_pexsi) !< Receive buffer for global 1D id
   integer :: recv_count(n_procs) !< Number of elements to receive from each processor
   integer :: recv_displ(n_procs) !< Displacement at which to place the incoming data
   integer :: recv_displ_aux !< Auxiliary variable used to set displacement

   integer :: offset
   integer :: mpierr

   character*40, parameter :: caller = "elsi_2dbc_to_1db_pexsi"

   send_count = 0
   send_displ = 0
   recv_count = 0
   recv_displ = 0

   ! Compute destination and global 1D id
   do i_col = 1, n_l_cols
      call elsi_get_global_col(global_col_id,i_col)
      do i_row = 1, n_l_rows
         call elsi_get_global_row(global_row_id,i_row)
         i_val = i_row+n_l_rows*(i_col-1)

         ! Compute destination to send the element
         dest(i_val) = FLOOR(1d0*(global_col_id-1)/FLOOR(1d0*n_g_rank/n_procs))
         ! The last process may take more
         if(dest(i_val) > (n_procs-1)) dest(i_val) = n_procs-1

         ! Compute the global 1d id of the element
         global_pos(i_val) = (global_col_id-1)*n_g_rank+global_row_id
     enddo
   enddo

   offset = 0

   ! Pack all data into send buffers, in the order of receiver
   do i_proc = 0, n_procs-1
      do i_val = 1, n_l_cols*n_l_rows
         if(dest(i_val) == i_proc) then
            ! Pack value
            local_row_id = mod(i_val,n_l_rows)
            if(local_row_id == 0) local_row_id = n_l_rows
            local_col_id = FLOOR(1d0*(i_val-1)/n_l_rows)+1

            ! TODO: this can be simplified due to symmetry!
            val_send_buffer(send_count(i_proc+1)+1+offset) = matrix_in(local_row_id,local_col_id)

            ! Pack global 1d id
            pos_send_buffer(send_count(i_proc+1)+1+offset) = global_pos(i_val)

            ! Set send_count
            send_count(i_proc+1) = send_count(i_proc+1)+1
         endif
      enddo
      offset = offset+send_count(i_proc+1)
   enddo

   ! Set recv_count
   call MPI_Alltoall(send_count, 1, mpi_integer, recv_count, &
                     1, mpi_integer, mpi_comm_global, mpierr)

   ! Set send and receive displacement
   send_displ_aux = 0
   recv_displ_aux = 0

   do i_proc = 0, n_procs-1
      send_displ(i_proc+1) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc+1)

      recv_displ(i_proc+1) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
   enddo

   ! Send and receive the packed data
   call MPI_Alltoallv(val_send_buffer, send_count, send_displ, mpi_real8, &
                      val_recv_buffer, recv_count, recv_displ, mpi_real8, &
                      mpi_comm_global, mpierr)

   call MPI_Alltoallv(pos_send_buffer, send_count, send_displ, mpi_integer, &
                      pos_recv_buffer, recv_count, recv_displ, mpi_integer, &
                      mpi_comm_global, mpierr)

   ! Unpack data
   do i_val = 1,n_l_rows_pexsi*n_l_cols_pexsi
      ! Compute global 2d id
      global_col_id = FLOOR(1d0*(pos_recv_buffer(i_val)-1)/n_g_rank)+1
      global_row_id = MOD(pos_recv_buffer(i_val),n_g_rank)
      if(global_row_id == 0) global_row_id = n_g_rank

      ! Compute local 2d id
      local_col_id = global_col_id-myid*FLOOR(1d0*n_g_rank/n_procs)
      local_row_id = global_row_id

      ! Put value to correct position
      matrix_out(local_row_id,local_col_id) = val_recv_buffer(i_val)
   enddo

end subroutine

end module ELSI_MATRIX_CONVERSION
