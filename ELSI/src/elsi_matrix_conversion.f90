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
   public :: elsi_2dbc_to_1db

contains

!> 
!! This routine computes global row index based on local row index.
!!
subroutine elsi_get_global_row(global_idx, local_idx)

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
subroutine elsi_get_global_col(global_idx, local_idx)

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
subroutine elsi_get_local_nnz(matrix, n_rows, n_cols, nnz)

   implicit none

   real*8,  intent(in)  :: matrix(n_rows,n_cols) !< Local matrix
   integer, intent(in)  :: n_rows                !< Local rows
   integer, intent(in)  :: n_cols                !< Local cols
   integer, intent(out) :: nnz                   !< Number of non-zero elements

   integer :: i_row !< Row counter
   integer :: i_col !< Column counter
   integer :: i_nz !< Non-zero element counter

   nnz = 0

   do i_col = 1, n_cols
      do i_row = 1, n_rows
         if(abs(matrix(i_row,i_col)) > zero_threshold) then
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
   call elsi_get_local_nnz(matrix, n_rows, n_cols, nnz_l)

   ! Set the number of non_zero elements in the global matrix
   call MPI_Allreduce(nnz_l, nnz_g, 1, mpi_integer, mpi_sum, &
                      mpi_comm_global, mpierr)

end subroutine

!>
!! This routine transforms a (local) dense matrix to sparse CCS format.
!!
subroutine elsi_dense_to_ccs(matrix, n_rows, n_cols, nnz, val, row_ind, col_ptr)

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
   do i_col = 1, n_cols
      new_col = .true.
      do i_row = 1, n_rows
         if(abs(matrix(i_row,i_col)) > zero_threshold) then
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
!! This routine transforms a (local) dense matrix to sparse CCS
!! format based on given col_ind and row_ptr.
!!
subroutine elsi_dense_to_ccs_by_pattern(matrix, n_rows, n_cols, &
                                        nnz, row_ind, col_ptr, val)

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
subroutine elsi_ccs_to_dense(matrix, n_rows, n_cols, val, nnz, row_ind, col_ptr)

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
   do i_val = 1, nnz
      if(i_val == col_ptr(i_col+1) .and. i_col /= n_cols) then
         i_col = i_col+1
      endif
      i_row = row_ind(i_val)
      matrix(i_row,i_col) = val(i_val)
   enddo

end subroutine

!> 
!! This routine converts the distribution of a dense matrix from
!! 2D block-cyclic to 1D block. The converted 1D dense matrix can
!! be easily converted to 1D block distributed sparse CCS format
!! used in PEXSI.
!!
subroutine elsi_2dbc_to_1db(matrix_in, matrix_out)

   implicit none
   include "mpif.h"

   real*8, intent(in)  :: matrix_in(n_l_rows,n_l_cols)
   real*8, intent(out) :: matrix_out(n_l_rows_pexsi, n_l_cols_pexsi)

   integer :: nnz_local1 !< Local number of nonzeros before redistribution
   integer :: nnz_local2 !< Local number of nonzeros after redistribution

   integer :: i_row !< Row counter
   integer :: i_col !< Col counter
   integer :: i_val !< Value counter
   integer :: i_proc !< Process counter
   integer :: global_col_id !< Global column id
   integer :: global_row_id !< Global row id
   integer :: local_col_id !< Local column id in 1D block distribution
   integer :: local_row_id !< Local row id in 1D block distribution

   integer, allocatable :: dest(:) !< Destination of each element

   ! For the meaning of each array here, see documentation of MPI_Alltoallv
   real*8, allocatable :: val_send_buffer(:) !< Send buffer for value
   integer, allocatable :: pos_send_buffer(:) !< Send buffer for global 1D id
   integer :: send_count(n_procs) !< Number of elements to send to each processor
   integer :: send_displ(n_procs) !< Displacement from which to take the outgoing data
   integer :: send_displ_aux !< Auxiliary variable used to set displacement

   real*8, allocatable :: val_recv_buffer(:) !< Receive buffer for value
   integer, allocatable :: pos_recv_buffer(:) !< Receive buffer for global 1D id
   integer :: recv_count(n_procs) !< Number of elements to receive from each processor
   integer :: recv_displ(n_procs) !< Displacement at which to place the incoming data
   integer :: recv_displ_aux !< Auxiliary variable used to set displacement

   integer :: mpierr

   character*40, parameter :: caller = "elsi_2dbc_to_1db_pexsi"

   send_count = 0
   send_displ = 0
   recv_count = 0
   recv_displ = 0

   call elsi_get_local_nnz(matrix_in,n_l_rows,n_l_cols,nnz_local1)

   call elsi_allocate(dest,nnz_local1,"dest",caller)
   call elsi_allocate(val_send_buffer,nnz_local1,"val_send_buffer",caller)
   call elsi_allocate(pos_send_buffer,nnz_local1,"pos_send_buffer",caller)

   i_val = 0
   ! Compute destination and global 1D id
   do i_col = 1, n_l_cols
      do i_row = 1, n_l_rows
         if(abs(matrix_in(i_row,i_col)) > zero_threshold) then
            i_val = i_val+1
            call elsi_get_global_col(global_col_id,i_col)
            call elsi_get_global_row(global_row_id,i_row)

            ! Compute destination
            dest(i_val) = FLOOR(1d0*(global_col_id-1)/FLOOR(1d0*n_g_size/n_procs))
            ! The last process may take more
            if(dest(i_val) > (n_procs-1)) dest(i_val) = n_procs-1

            ! Compute the global id
            ! Pack global id and data into buffers
            pos_send_buffer(i_val) = (global_col_id-1)*n_g_size+global_row_id
            val_send_buffer(i_val) = matrix_in(i_row,i_col)
        endif
     enddo
   enddo

   ! Set send_count
   do i_proc = 0, n_procs-1
      do i_val = 1, nnz_local1
         if(dest(i_val) == i_proc) then
            send_count(i_proc+1) = send_count(i_proc+1)+1
         endif
      enddo
   enddo

   deallocate(dest)

   ! Set recv_count
   call MPI_Alltoall(send_count, 1, mpi_integer, recv_count, &
                     1, mpi_integer, mpi_comm_global, mpierr)

   nnz_local2 = sum(recv_count,1)

   ! Set send and receive displacement
   send_displ_aux = 0
   recv_displ_aux = 0

   do i_proc = 0, n_procs-1
      send_displ(i_proc+1) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc+1)

      recv_displ(i_proc+1) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
   enddo

   call elsi_allocate(val_recv_buffer,nnz_local2,"val_recv_buffer",caller)
   call elsi_allocate(pos_recv_buffer,nnz_local2,"pos_recv_buffer",caller)

   ! Send and receive the packed data
   call MPI_Alltoallv(val_send_buffer, send_count, send_displ, mpi_real8, &
                      val_recv_buffer, recv_count, recv_displ, mpi_real8, &
                      mpi_comm_global, mpierr)

   call MPI_Alltoallv(pos_send_buffer, send_count, send_displ, mpi_integer, &
                      pos_recv_buffer, recv_count, recv_displ, mpi_integer, &
                      mpi_comm_global, mpierr)

   deallocate(val_send_buffer)
   deallocate(pos_send_buffer)

   matrix_out = 0d0

   ! Unpack data
   do i_val = 1,nnz_local2
      ! Compute global 2d id
      global_col_id = FLOOR(1d0*(pos_recv_buffer(i_val)-1)/n_g_size)+1
      global_row_id = MOD(pos_recv_buffer(i_val),n_g_size)
      if(global_row_id == 0) global_row_id = n_g_size

      ! Compute local 2d id
      local_col_id = global_col_id-myid*FLOOR(1d0*n_g_size/n_procs)
      local_row_id = global_row_id

      ! Put value to correct position
      matrix_out(local_row_id,local_col_id) = val_recv_buffer(i_val)
   enddo

   deallocate(val_recv_buffer)
   deallocate(pos_recv_buffer)

end subroutine

end module ELSI_MATRIX_CONVERSION
