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
!! This module performs parallel matrix IO.
!!
module ELSI_IO

   use, intrinsic :: ISO_C_BINDING
   use ELSI_CONSTANTS, only: SINGLE_PROC,HEADER_SIZE,BLACS_DENSE,PEXSI_CSC
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_SETUP
   use ELSI_MATCONV, only: elsi_pexsi_to_blacs_dm,elsi_blacs_to_sips_hs
   use ELSI_PRECISION, only: r8,i4,i8
   use ELSI_UTILS

   implicit none

   private

   public :: elsi_read_mat_dim
   public :: elsi_read_mat_dim_sparse
   public :: elsi_read_mat_real
   public :: elsi_read_mat_real_sparse
   public :: elsi_read_mat_complex
   public :: elsi_read_mat_complex_sparse
   public :: elsi_write_mat_real
   public :: elsi_write_mat_real_sparse
   public :: elsi_write_mat_complex
   public :: elsi_write_mat_complex_sparse

contains

!>
!! This routine reads the dimensions of a 2D block-cyclic dense matrix
!! from file.
!!
subroutine elsi_read_mat_dim(f_name,mpi_comm,blacs_ctxt,block_size,&
              n_electron,n_basis,n_l_rows,n_l_cols)

   implicit none

   character(*),     intent(in)  :: f_name     !< File name
   integer(kind=i4), intent(in)  :: mpi_comm   !< MPI communicator
   integer(kind=i4), intent(in)  :: blacs_ctxt !< BLACS context
   integer(kind=i4), intent(in)  :: block_size !< Block size
   real(kind=r8),    intent(out) :: n_electron !< Number of electrons
   integer(kind=i4), intent(out) :: n_basis    !< Matrix size
   integer(kind=i4), intent(out) :: n_l_rows   !< Local number of rows
   integer(kind=i4), intent(out) :: n_l_cols   !< Local number of columns

   integer(kind=i4) :: myid
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_p_rows
   integer(kind=i4) :: n_p_cols
   integer(kind=i4) :: my_p_row
   integer(kind=i4) :: my_p_col
   integer(kind=i8) :: offset

   integer(kind=i4), external :: numroc

   character*40, parameter :: caller = "elsi_read_mat_dim"

   call MPI_Comm_rank(mpi_comm,myid,mpierr)

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Read header
   if(myid == 0) then
      offset = 0

      call MPI_File_read_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Broadcast header
   call MPI_Bcast(header,HEADER_SIZE,mpi_integer4,0,mpi_comm,mpierr)

   n_basis    = header(2)
   n_electron = real(header(3),kind=r8)

   ! Get processor grid information
   call blacs_gridinfo(blacs_ctxt,n_p_rows,n_p_cols,my_p_row,my_p_col)

   ! Get local size of matrix
   n_l_rows = numroc(n_basis,block_size,my_p_row,0,n_p_rows)
   n_l_cols = numroc(n_basis,block_size,my_p_col,0,n_p_cols)

end subroutine

!>
!! This routine reads the dimensions of a 1D block CSC matrix from file.
!!
subroutine elsi_read_mat_dim_sparse(f_name,mpi_comm,n_electron,n_basis,&
              nnz_g,nnz_l,n_l_cols)

   implicit none

   character(*),     intent(in)  :: f_name     !< File name
   integer(kind=i4), intent(in)  :: mpi_comm   !< MPI communicator
   real(kind=r8),    intent(out) :: n_electron !< Number of electrons
   integer(kind=i4), intent(out) :: n_basis    !< Matrix size
   integer(kind=i4), intent(out) :: nnz_g      !< Global number of nonzeros
   integer(kind=i4), intent(out) :: nnz_l      !< Local number of nonzeros
   integer(kind=i4), intent(out) :: n_l_cols   !< Local number of columns

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_l_cols0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset

   integer(kind=i4), allocatable :: col_ptr(:)

   type(elsi_handle) :: io_h

   character*40, parameter :: caller = "elsi_read_mat_dim_sparse"

   call elsi_init(io_h,ELPA,MULTI_PROC,BLACS_DENSE,n_basis,0.0_r8,0)
   call elsi_set_mpi(io_h,mpi_comm)
   call elsi_set_mpi_global(io_h,mpi_comm)

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Read header
   if(io_h%myid == 0) then
      offset = 0

      call MPI_File_read_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Broadcast header
   call MPI_Bcast(header,HEADER_SIZE,mpi_integer4,0,mpi_comm,mpierr)

   n_basis    = header(2)
   n_electron = real(header(3),kind=r8)
   nnz_g      = header(4)

   ! Compute n_l_cols
   n_l_cols  = n_basis/io_h%n_procs
   n_l_cols0 = n_l_cols
   if(io_h%myid == io_h%n_procs-1) then
      n_l_cols = n_basis-(io_h%n_procs-1)*n_l_cols0
   endif

   call elsi_allocate(io_h,col_ptr,n_l_cols+1,"col_ptr",caller)

   ! Read column pointer
   offset = HEADER_SIZE*4+io_h%myid*n_l_cols0*4

   call MPI_File_read_at_all(f_handle,offset,col_ptr,n_l_cols+1,mpi_integer4,&
           mpi_status_ignore,mpierr)

   if(io_h%myid == io_h%n_procs-1) then
      col_ptr(n_l_cols+1) = nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = col_ptr(1)-1
   col_ptr  = col_ptr-prev_nnz

   ! Compute nnz_l
   nnz_l = col_ptr(n_l_cols+1)-col_ptr(1)

   call elsi_deallocate(io_h,col_ptr,"col_ptr")
   call elsi_cleanup(io_h)

end subroutine

!>
!! This routine reads a 2D block-cyclic dense matrix from file.
!!
subroutine elsi_read_mat_real(f_name,mpi_comm,blacs_ctxt,block_size,&
              n_basis,n_l_rows,n_l_cols,mat)

   implicit none

   character(*),     intent(in)  :: f_name                 !< File name
   integer(kind=i4), intent(in)  :: mpi_comm               !< MPI communicator
   integer(kind=i4), intent(in)  :: blacs_ctxt             !< BLACS context
   integer(kind=i4), intent(in)  :: block_size             !< Block size
   integer(kind=i4), intent(in)  :: n_basis                !< Matrix size
   integer(kind=i4), intent(in)  :: n_l_rows               !< Local number of rows
   integer(kind=i4), intent(in)  :: n_l_cols               !< Local number of columns
   real(kind=r8),    intent(out) :: mat(n_l_rows,n_l_cols) !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_l_cols_sp
   integer(kind=i4) :: n_l_cols0
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)
   real(kind=r8),    allocatable :: nnz_val(:)

   type(elsi_handle) :: io_h

   character*40, parameter :: caller = "elsi_read_mat_real"

   call elsi_init(io_h,PEXSI,MULTI_PROC,BLACS_DENSE,n_basis,0.0_r8,0)
   call elsi_set_mpi(io_h,mpi_comm)
   call elsi_set_mpi_global(io_h,mpi_comm)
   call elsi_set_blacs(io_h,blacs_ctxt,block_size)
   call elsi_get_time(io_h,t0)

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Read header
   if(io_h%myid == 0) then
      offset = 0

      call MPI_File_read_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Broadcast header
   call MPI_Bcast(header,HEADER_SIZE,mpi_integer4,0,mpi_comm,mpierr)

   nnz_g = header(4)

   ! Compute n_l_cols_sp
   n_l_cols_sp  = n_basis/io_h%n_procs
   n_l_cols0    = n_l_cols_sp
   if(io_h%myid == io_h%n_procs-1) then
      n_l_cols_sp = n_basis-(io_h%n_procs-1)*n_l_cols0
   endif

   call elsi_allocate(io_h,col_ptr,n_l_cols_sp+1,"col_ptr",caller)

   ! Read column pointer
   offset = HEADER_SIZE*4+io_h%myid*n_l_cols0*4

   call MPI_File_read_at_all(f_handle,offset,col_ptr,n_l_cols_sp+1,&
           mpi_integer4,mpi_status_ignore,mpierr)

   if(io_h%myid == io_h%n_procs-1) then
      col_ptr(n_l_cols_sp+1) = nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = col_ptr(1)-1
   col_ptr  = col_ptr-prev_nnz

   ! Compute nnz_l
   nnz_l = col_ptr(n_l_cols_sp+1)-col_ptr(1)

   call elsi_allocate(io_h,row_ind,nnz_l,"row_ind",caller)

   ! Read row index
   offset = HEADER_SIZE*4+n_basis*4+prev_nnz*4

   call MPI_File_read_at_all(f_handle,offset,row_ind,nnz_l,mpi_integer4,&
           mpi_status_ignore,mpierr)

   ! Read non-zero value
   offset = HEADER_SIZE*4+n_basis*4+nnz_g*4+prev_nnz*8

   call elsi_allocate(io_h,nnz_val,nnz_l,"nnz_val",caller)

   call MPI_File_read_at_all(f_handle,offset,nnz_val,nnz_l,mpi_real8,&
           mpi_status_ignore,mpierr)

   ! Close file
   call MPI_File_close(f_handle,mpierr)

   ! Redistribute matrix
   call elsi_set_csc(io_h,nnz_g,nnz_l,n_l_cols_sp,row_ind,col_ptr)
   call elsi_set_sparse_dm(io_h,nnz_val)

   io_h%n_elsi_calls   = 1
   io_h%my_p_row_pexsi = 0
   io_h%n_p_per_pole   = io_h%n_procs

   call elsi_pexsi_to_blacs_dm(io_h,mat)

   call elsi_deallocate(io_h,col_ptr,"col_ptr")
   call elsi_deallocate(io_h,row_ind,"row_ind")
   call elsi_deallocate(io_h,nnz_val,"nnz_val")

   call elsi_get_time(io_h,t1)

   write(info_str,"('  Finished reading matrix')")
   call elsi_statement_print(info_str,io_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,io_h)

   call elsi_cleanup(io_h)

end subroutine

!>
!! This routine reads a 1D block CSC matrix from file.
!!
subroutine elsi_read_mat_real_sparse(f_name,mpi_comm,n_basis,nnz_g,nnz_l,&
              n_l_cols,row_ind,col_ptr,mat)

   implicit none

   character(*),     intent(in)  :: f_name              !< File name
   integer(kind=i4), intent(in)  :: mpi_comm            !< MPI communicator
   integer(kind=i4), intent(in)  :: n_basis             !< Matrix size
   integer(kind=i4), intent(in)  :: nnz_g               !< Global number of nonzeros
   integer(kind=i4), intent(in)  :: nnz_l               !< Local number of nonzeros
   integer(kind=i4), intent(in)  :: n_l_cols            !< Local number of columns
   integer(kind=i4), intent(out) :: row_ind(nnz_l)      !< Row index
   integer(kind=i4), intent(out) :: col_ptr(n_l_cols+1) !< Column pointer
   real(kind=r8),    intent(out) :: mat(nnz_l)          !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_l_cols0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   type(elsi_handle) :: io_h

   character*40, parameter :: caller = "elsi_read_mat_real_sparse"

   call elsi_init(io_h,ELPA,MULTI_PROC,BLACS_DENSE,n_basis,0.0_r8,0)
   call elsi_set_mpi(io_h,mpi_comm)
   call elsi_set_mpi_global(io_h,mpi_comm)
   call elsi_get_time(io_h,t0)

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Compute n_l_cols0
   n_l_cols0 = n_basis/io_h%n_procs

   ! Read column pointer
   offset = HEADER_SIZE*4+io_h%myid*n_l_cols0*4

   call MPI_File_read_at_all(f_handle,offset,col_ptr,n_l_cols+1,mpi_integer4,&
           mpi_status_ignore,mpierr)

   if(io_h%myid == io_h%n_procs-1) then
      col_ptr(n_l_cols+1) = nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = col_ptr(1)-1
   col_ptr  = col_ptr-prev_nnz

   ! Read row index
   offset = HEADER_SIZE*4+n_basis*4+prev_nnz*4

   call MPI_File_read_at_all(f_handle,offset,row_ind,nnz_l,mpi_integer4,&
           mpi_status_ignore,mpierr)

   ! Read non-zero value
   offset = HEADER_SIZE*4+n_basis*4+nnz_g*4+prev_nnz*8

   call MPI_File_read_at_all(f_handle,offset,mat,nnz_l,mpi_real8,&
           mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(io_h,t1)

   write(info_str,"('  Finished reading matrix')")
   call elsi_statement_print(info_str,io_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,io_h)

   call elsi_cleanup(io_h)

end subroutine

!>
!! This routine reads a 2D block-cyclic dense matrix from file.
!!
subroutine elsi_read_mat_complex(f_name,mpi_comm,blacs_ctxt,block_size,&
              n_basis,n_l_rows,n_l_cols,mat)

   implicit none

   character(*),     intent(in)  :: f_name                 !< File name
   integer(kind=i4), intent(in)  :: mpi_comm               !< MPI communicator
   integer(kind=i4), intent(in)  :: blacs_ctxt             !< BLACS context
   integer(kind=i4), intent(in)  :: block_size             !< Block size
   integer(kind=i4), intent(in)  :: n_basis                !< Matrix size
   integer(kind=i4), intent(in)  :: n_l_rows               !< Local number of rows
   integer(kind=i4), intent(in)  :: n_l_cols               !< Local number of columns
   complex(kind=r8), intent(out) :: mat(n_l_rows,n_l_cols) !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_l_cols_sp
   integer(kind=i4) :: n_l_cols0
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)
   complex(kind=r8), allocatable :: nnz_val(:)

   type(elsi_handle) :: io_h

   character*40, parameter :: caller = "elsi_read_mat_complex"

   call elsi_init(io_h,PEXSI,MULTI_PROC,BLACS_DENSE,n_basis,0.0_r8,0)
   call elsi_set_mpi(io_h,mpi_comm)
   call elsi_set_mpi_global(io_h,mpi_comm)
   call elsi_set_blacs(io_h,blacs_ctxt,block_size)
   call elsi_get_time(io_h,t0)

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Read header
   if(io_h%myid == 0) then
      offset = 0

      call MPI_File_read_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Broadcast header
   call MPI_Bcast(header,HEADER_SIZE,mpi_integer4,0,mpi_comm,mpierr)

   nnz_g = header(4)

   ! Compute n_l_cols_sp
   n_l_cols_sp  = n_basis/io_h%n_procs
   n_l_cols0    = n_l_cols_sp
   if(io_h%myid == io_h%n_procs-1) then
      n_l_cols_sp = n_basis-(io_h%n_procs-1)*n_l_cols0
   endif

   call elsi_allocate(io_h,col_ptr,n_l_cols_sp+1,"col_ptr",caller)

   ! Read column pointer
   offset = HEADER_SIZE*4+io_h%myid*n_l_cols0*4

   call MPI_File_read_at_all(f_handle,offset,col_ptr,n_l_cols_sp+1,&
           mpi_integer4,mpi_status_ignore,mpierr)

   if(io_h%myid == io_h%n_procs-1) then
      col_ptr(n_l_cols_sp+1) = nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = col_ptr(1)-1
   col_ptr  = col_ptr-prev_nnz

   ! Compute nnz_l
   nnz_l = col_ptr(n_l_cols_sp+1)-col_ptr(1)

   call elsi_allocate(io_h,row_ind,nnz_l,"row_ind",caller)

   ! Read row index
   offset = HEADER_SIZE*4+n_basis*4+prev_nnz*4

   call MPI_File_read_at_all(f_handle,offset,row_ind,nnz_l,mpi_integer4,&
           mpi_status_ignore,mpierr)

   ! Read non-zero value
   offset = HEADER_SIZE*4+n_basis*4+nnz_g*4+prev_nnz*16

   call elsi_allocate(io_h,nnz_val,nnz_l,"nnz_val",caller)

   call MPI_File_read_at_all(f_handle,offset,nnz_val,nnz_l,mpi_complex16,&
           mpi_status_ignore,mpierr)

   ! Close file
   call MPI_File_close(f_handle,mpierr)

   ! Redistribute matrix
   call elsi_set_csc(io_h,nnz_g,nnz_l,n_l_cols_sp,row_ind,col_ptr)
   call elsi_set_sparse_dm(io_h,nnz_val)

   io_h%n_elsi_calls   = 1
   io_h%my_p_row_pexsi = 0
   io_h%n_p_per_pole   = io_h%n_procs

   call elsi_pexsi_to_blacs_dm(io_h,mat)

   call elsi_deallocate(io_h,col_ptr,"col_ptr")
   call elsi_deallocate(io_h,row_ind,"row_ind")
   call elsi_deallocate(io_h,nnz_val,"nnz_val")

   call elsi_get_time(io_h,t1)

   write(info_str,"('  Finished reading matrix')")
   call elsi_statement_print(info_str,io_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,io_h)

   call elsi_cleanup(io_h)

end subroutine

!>
!! This routine reads a 1D block CSC matrix from file.
!!
subroutine elsi_read_mat_complex_sparse(f_name,mpi_comm,n_basis,nnz_g,&
              nnz_l,n_l_cols,row_ind,col_ptr,mat)

   implicit none

   character(*),     intent(in)  :: f_name              !< File name
   integer(kind=i4), intent(in)  :: mpi_comm            !< MPI communicator
   integer(kind=i4), intent(in)  :: n_basis             !< Matrix size
   integer(kind=i4), intent(in)  :: nnz_g               !< Global number of nonzeros
   integer(kind=i4), intent(in)  :: nnz_l               !< Local number of nonzeros
   integer(kind=i4), intent(in)  :: n_l_cols            !< Local number of columns
   integer(kind=i4), intent(out) :: row_ind(nnz_l)      !< Row index
   integer(kind=i4), intent(out) :: col_ptr(n_l_cols+1) !< Column pointer
   complex(kind=r8), intent(out) :: mat(nnz_l)          !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_l_cols0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   type(elsi_handle) :: io_h

   character*40, parameter :: caller = "elsi_read_mat_complex_sparse"

   call elsi_init(io_h,ELPA,MULTI_PROC,BLACS_DENSE,n_basis,0.0_r8,0)
   call elsi_set_mpi(io_h,mpi_comm)
   call elsi_set_mpi_global(io_h,mpi_comm)
   call elsi_get_time(io_h,t0)

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Compute n_l_cols0
   n_l_cols0 = n_basis/io_h%n_procs

   ! Read column pointer
   offset = HEADER_SIZE*4+io_h%myid*n_l_cols0*4

   call MPI_File_read_at_all(f_handle,offset,col_ptr,n_l_cols+1,mpi_integer4,&
           mpi_status_ignore,mpierr)

   if(io_h%myid == io_h%n_procs-1) then
      col_ptr(n_l_cols+1) = nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = col_ptr(1)-1
   col_ptr  = col_ptr-prev_nnz

   ! Read row index
   offset = HEADER_SIZE*4+n_basis*4+prev_nnz*4

   call MPI_File_read_at_all(f_handle,offset,row_ind,nnz_l,mpi_integer4,&
           mpi_status_ignore,mpierr)

   ! Read non-zero value
   offset = HEADER_SIZE*4+n_basis*4+nnz_g*4+prev_nnz*16

   call MPI_File_read_at_all(f_handle,offset,mat,nnz_l,mpi_complex16,&
           mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(io_h,t1)

   write(info_str,"('  Finished reading matrix')")
   call elsi_statement_print(info_str,io_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,io_h)

   call elsi_cleanup(io_h)

end subroutine

!>
!! This routine writes a 2D block-cyclic dense matrix to file.
!!
subroutine elsi_write_mat_real(f_name,mpi_comm,blacs_ctxt,block_size,&
              n_electron,n_basis,n_l_rows,n_l_cols,mat)

   implicit none

   character(*),     intent(in) :: f_name                 !< File name
   integer(kind=i4), intent(in) :: mpi_comm               !< MPI communicator
   integer(kind=i4), intent(in) :: blacs_ctxt             !< BLACS context
   integer(kind=i4), intent(in) :: block_size             !< Block size
   real(kind=r8),    intent(in) :: n_electron             !< Number of electrons
   integer(kind=i4), intent(in) :: n_basis                !< Matrix size
   integer(kind=i4), intent(in) :: n_l_rows               !< Local number of rows
   integer(kind=i4), intent(in) :: n_l_cols               !< Local number of columns
   real(kind=r8),    intent(in) :: mat(n_l_rows,n_l_cols) !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_l_cols0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   type(elsi_handle) :: io_h

   character*40, parameter :: caller = "elsi_write_mat_real"

   call elsi_init(io_h,ELPA,MULTI_PROC,BLACS_DENSE,n_basis,0.0_r8,0)
   call elsi_set_mpi(io_h,mpi_comm)
   call elsi_set_mpi_global(io_h,mpi_comm)
   call elsi_set_blacs(io_h,blacs_ctxt,block_size)
   call elsi_get_time(io_h,t0)

   io_h%n_elsi_calls = 1
   io_h%n_l_cols_sp = io_h%n_basis/io_h%n_procs
   if(io_h%myid == io_h%n_procs-1) then
      io_h%n_l_cols_sp = io_h%n_basis-(io_h%n_procs-1)*io_h%n_l_cols_sp
   endif

   call elsi_blacs_to_sips_hs(io_h,mat,mat)

   ! Open file
   f_mode = mpi_mode_wronly+mpi_mode_create

   call MPI_File_open(mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Write header
   header(1) = PEXSI_CSC
   header(2) = n_basis
   header(3) = int(n_electron,kind=i4)
   header(4) = io_h%nnz_g

   if(io_h%myid == 0) then
      offset = 0

      call MPI_File_write_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Compute shift of column pointers
   prev_nnz = 0

   call MPI_Exscan(io_h%nnz_l_sp,prev_nnz,1,mpi_integer4,mpi_sum,mpi_comm,&
           mpierr)

   ! Shift column pointer
   io_h%col_ptr_sips = io_h%col_ptr_sips+prev_nnz

   ! Write column pointer
   n_l_cols0 = n_basis/io_h%n_procs
   offset = HEADER_SIZE*4+io_h%myid*n_l_cols0*4

   call MPI_File_write_at_all(f_handle,offset,io_h%col_ptr_sips,&
           io_h%n_l_cols_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Write row index
   offset = HEADER_SIZE*4+n_basis*4+prev_nnz*4

   call MPI_File_write_at_all(f_handle,offset,io_h%row_ind_sips,&
           io_h%nnz_l_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Write non-zero value
   offset = HEADER_SIZE*4+n_basis*4+io_h%nnz_g*4+prev_nnz*8

   call MPI_File_write_at_all(f_handle,offset,io_h%ham_real_sips,&
           io_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(io_h,t1)

   write(info_str,"('  Finished writing matrix')")
   call elsi_statement_print(info_str,io_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,io_h)

   call elsi_cleanup(io_h)

end subroutine

!>
!! This routine writes a 2D block-cyclic dense matrix to file.
!!
subroutine elsi_write_mat_complex(f_name,mpi_comm,blacs_ctxt,block_size,&
              n_electron,n_basis,n_l_rows,n_l_cols,mat)

   implicit none

   character(*),     intent(in) :: f_name                 !< File name
   integer(kind=i4), intent(in) :: mpi_comm               !< MPI communicator
   integer(kind=i4), intent(in) :: blacs_ctxt             !< BLACS context
   integer(kind=i4), intent(in) :: block_size             !< Block size
   real(kind=r8),    intent(in) :: n_electron             !< Number of electrons
   integer(kind=i4), intent(in) :: n_basis                !< Matrix size
   integer(kind=i4), intent(in) :: n_l_rows               !< Local number of rows
   integer(kind=i4), intent(in) :: n_l_cols               !< Local number of columns
   complex(kind=r8), intent(in) :: mat(n_l_rows,n_l_cols) !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_l_cols0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   type(elsi_handle) :: io_h

   character*40, parameter :: caller = "elsi_write_mat_complex"

   call elsi_init(io_h,ELPA,MULTI_PROC,BLACS_DENSE,n_basis,0.0_r8,0)
   call elsi_set_mpi(io_h,mpi_comm)
   call elsi_set_mpi_global(io_h,mpi_comm)
   call elsi_set_blacs(io_h,blacs_ctxt,block_size)
   call elsi_get_time(io_h,t0)

   io_h%n_elsi_calls = 1
   io_h%n_l_cols_sp = io_h%n_basis/io_h%n_procs
   if(io_h%myid == io_h%n_procs-1) then
      io_h%n_l_cols_sp = io_h%n_basis-(io_h%n_procs-1)*io_h%n_l_cols_sp
   endif

   call elsi_blacs_to_sips_hs(io_h,mat,mat)

   ! Open file
   f_mode = mpi_mode_wronly+mpi_mode_create

   call MPI_File_open(mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Write header
   header(1) = PEXSI_CSC
   header(2) = n_basis
   header(3) = int(n_electron,kind=i4)
   header(4) = io_h%nnz_g

   if(io_h%myid == 0) then
      offset = 0

      call MPI_File_write_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Compute shift of column pointers
   prev_nnz = 0

   call MPI_Exscan(io_h%nnz_l_sp,prev_nnz,1,mpi_integer4,mpi_sum,mpi_comm,&
           mpierr)

   ! Shift column pointer
   io_h%col_ptr_sips = io_h%col_ptr_sips+prev_nnz

   ! Write column pointer
   n_l_cols0 = n_basis/io_h%n_procs
   offset = HEADER_SIZE*4+io_h%myid*n_l_cols0*4

   call MPI_File_write_at_all(f_handle,offset,io_h%col_ptr_sips,&
           io_h%n_l_cols_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Write row index
   offset = HEADER_SIZE*4+n_basis*4+prev_nnz*4

   call MPI_File_write_at_all(f_handle,offset,io_h%row_ind_sips,&
           io_h%nnz_l_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Write non-zero value
   offset = HEADER_SIZE*4+n_basis*4+io_h%nnz_g*4+prev_nnz*16

   call MPI_File_write_at_all(f_handle,offset,io_h%ham_complex_sips,&
           io_h%nnz_l_sp,mpi_complex16,mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(io_h,t1)

   write(info_str,"('  Finished writing matrix')")
   call elsi_statement_print(info_str,io_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,io_h)

   call elsi_cleanup(io_h)

end subroutine

!>
!! This routine writes a 1D block CSC matrix to file.
!!
subroutine elsi_write_mat_real_sparse(f_name,mpi_comm,n_electron,n_basis,&
              nnz_g,nnz_l,n_l_cols,row_ind,col_ptr,mat)

   implicit none

   character(*),     intent(in) :: f_name              !< File name
   integer(kind=i4), intent(in) :: mpi_comm            !< MPI communicator
   real(kind=r8),    intent(in) :: n_electron          !< Number of electrons
   integer(kind=i4), intent(in) :: n_basis             !< Matrix size
   integer(kind=i4), intent(in) :: nnz_g               !< Global number of nonzeros
   integer(kind=i4), intent(in) :: nnz_l               !< Local number of nonzeros
   integer(kind=i4), intent(in) :: n_l_cols            !< Local number of columns
   integer(kind=i4), intent(in) :: row_ind(nnz_l)      !< Row index
   integer(kind=i4), intent(in) :: col_ptr(n_l_cols+1) !< Column pointer
   real(kind=r8),    intent(in) :: mat(nnz_l)          !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: prev_nnz
   integer(kind=i4) :: n_l_cols0
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: col_ptr_shift(:)

   type(elsi_handle) :: io_h

   character*40, parameter :: caller = "elsi_write_mat_real_sparse"

   call elsi_init(io_h,ELPA,MULTI_PROC,BLACS_DENSE,n_basis,0.0_r8,0)
   call elsi_set_mpi(io_h,mpi_comm)
   call elsi_set_mpi_global(io_h,mpi_comm)
   call elsi_get_time(io_h,t0)

   ! Open file
   f_mode = mpi_mode_wronly+mpi_mode_create

   call MPI_File_open(mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Write header
   header(1) = PEXSI_CSC
   header(2) = n_basis
   header(3) = int(n_electron,kind=i4)
   header(4) = nnz_g

   if(io_h%myid == 0) then
      offset = 0

      call MPI_File_write_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Compute shift of column pointers
   prev_nnz = 0

   call MPI_Exscan(nnz_l,prev_nnz,1,mpi_integer4,mpi_sum,mpi_comm,mpierr)

   ! Shift column pointer
   call elsi_allocate(io_h,col_ptr_shift,n_l_cols+1,"col_ptr_shift",caller)

   col_ptr_shift = col_ptr+prev_nnz

   ! Write column pointer
   n_l_cols0 = n_basis/io_h%n_procs
   offset = HEADER_SIZE*4+io_h%myid*n_l_cols0*4

   call MPI_File_write_at_all(f_handle,offset,col_ptr_shift,n_l_cols,&
           mpi_integer4,mpi_status_ignore,mpierr)

   call elsi_deallocate(io_h,col_ptr_shift,"col_ptr_shift")

   ! Write row index
   offset = HEADER_SIZE*4+n_basis*4+prev_nnz*4

   call MPI_File_write_at_all(f_handle,offset,row_ind,nnz_l,mpi_integer4,&
           mpi_status_ignore,mpierr)

   ! Write non-zero value
   offset = HEADER_SIZE*4+n_basis*4+nnz_g*4+prev_nnz*8

   call MPI_File_write_at_all(f_handle,offset,mat,nnz_l,mpi_real8,&
           mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(io_h,t1)

   write(info_str,"('  Finished writing matrix')")
   call elsi_statement_print(info_str,io_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,io_h)

   call elsi_cleanup(io_h)

end subroutine

!>
!! This routine writes a 1D block CSC matrix to file.
!!
subroutine elsi_write_mat_complex_sparse(f_name,mpi_comm,n_electron,&
              n_basis,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr,mat)

   implicit none

   character(*),     intent(in) :: f_name              !< File name
   integer(kind=i4), intent(in) :: mpi_comm            !< MPI communicator
   real(kind=r8),    intent(in) :: n_electron          !< Number of electrons
   integer(kind=i4), intent(in) :: n_basis             !< Matrix size
   integer(kind=i4), intent(in) :: nnz_g               !< Global number of nonzeros
   integer(kind=i4), intent(in) :: nnz_l               !< Local number of nonzeros
   integer(kind=i4), intent(in) :: n_l_cols            !< Local number of columns
   integer(kind=i4), intent(in) :: row_ind(nnz_l)      !< Row index
   integer(kind=i4), intent(in) :: col_ptr(n_l_cols+1) !< Column pointer
   complex(kind=r8), intent(in) :: mat(nnz_l)          !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: prev_nnz
   integer(kind=i4) :: n_l_cols0
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: col_ptr_shift(:)

   type(elsi_handle) :: io_h

   character*40, parameter :: caller = "elsi_write_mat_complex_sparse"

   call elsi_init(io_h,ELPA,MULTI_PROC,BLACS_DENSE,n_basis,0.0_r8,0)
   call elsi_set_mpi(io_h,mpi_comm)
   call elsi_set_mpi_global(io_h,mpi_comm)
   call elsi_get_time(io_h,t0)

   ! Open file
   f_mode = mpi_mode_wronly+mpi_mode_create

   call MPI_File_open(mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Write header
   header(1) = PEXSI_CSC
   header(2) = n_basis
   header(3) = int(n_electron,kind=i4)
   header(4) = nnz_g

   if(io_h%myid == 0) then
      offset = 0

      call MPI_File_write_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Compute shift of column pointers
   prev_nnz = 0

   call MPI_Exscan(nnz_l,prev_nnz,1,mpi_integer4,mpi_sum,mpi_comm,mpierr)

   ! Shift column pointer
   call elsi_allocate(io_h,col_ptr_shift,n_l_cols+1,"col_ptr_shift",caller)

   col_ptr_shift = col_ptr+prev_nnz

   ! Write column pointer
   n_l_cols0 = n_basis/io_h%n_procs
   offset = HEADER_SIZE*4+io_h%myid*n_l_cols0*4

   call MPI_File_write_at_all(f_handle,offset,col_ptr_shift,n_l_cols,&
           mpi_integer4,mpi_status_ignore,mpierr)

   call elsi_deallocate(io_h,col_ptr_shift,"col_ptr_shift")

   ! Write row index
   offset = HEADER_SIZE*4+n_basis*4+prev_nnz*4

   call MPI_File_write_at_all(f_handle,offset,row_ind,nnz_l,mpi_integer4,&
           mpi_status_ignore,mpierr)

   ! Write non-zero value
   offset = HEADER_SIZE*4+n_basis*4+nnz_g*4+prev_nnz*16

   call MPI_File_write_at_all(f_handle,offset,mat,nnz_l,mpi_complex16,&
           mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(io_h,t1)

   write(info_str,"('  Finished writing matrix')")
   call elsi_statement_print(info_str,io_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,io_h)

   call elsi_cleanup(io_h)

end subroutine

end module ELSI_IO
