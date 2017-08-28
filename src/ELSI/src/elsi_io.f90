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
   use ELSI_CONSTANTS, only: SINGLE_PROC,HEADER_SIZE
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_MATCONV
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS

   implicit none

   private

!   public :: elsi_write_matrix_real
   public :: elsi_write_matrix_real_sparse
!   public :: elsi_read_matrix_real
   public :: elsi_read_matrix_real_sparse

contains

!>
!! This routine reads a 1D block distributed CSC matrix from file.
!!
subroutine elsi_read_matrix_real_sparse(elsi_h,filename,id)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h   !< Handle
   character(*),      intent(in)    :: filename !< File name
   integer(kind=i4),  intent(in)    :: id       !< H,S,D,other

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: filehandle
   integer(kind=i4) :: filemode
   integer(kind=i4) :: header(3)
   integer(kind=i4) :: n_l_cols0
   integer(kind=i4) :: prev_nnz

   integer(kind=mpi_offset_kind) :: offset

   character*40, parameter :: caller = "elsi_read_matrix_real_sparse"

   ! Open file
   filemode = mpi_mode_rdonly

   call MPI_File_open(elsi_h%mpi_comm,filename,filemode,mpi_info_null,&
           filehandle,mpierr)

   ! Read header
   if(elsi_h%myid == 0) then
      offset = 0

      call MPI_File_read_at(filehandle,offset,header,3,mpi_integer4,&
              mpi_status_ignore,mpierr)

!      call MPI_File_read_at(filehandle,offset,n_basis,1,mpi_integer4,&
!              mpi_status_ignore,mpierr)
!      header(1) = n_basis
!      offset = 4
!      call MPI_File_read_at(filehandle,offset,n_electron,1,mpi_integer4,&
!              mpi_status_ignore,mpierr)
!      header(2) = n_electron
!      offset = 8
!      call MPI_File_read_at(filehandle,offset,nnz_g,1,mpi_integer4,&
!              mpi_status_ignore,mpierr)
!      header(3) = nnz_g
   endif

   ! Broadcast header
   call MPI_Bcast(header,3,mpi_integer4,0,elsi_h%mpi_comm,mpierr)

   ! Set ELSI parameters
   elsi_h%n_basis     = header(1)
   elsi_h%n_nonsing   = header(1)
   elsi_h%n_electrons = real(header(2),kind=r8)
   elsi_h%nnz_g       = header(3)

   if(elsi_h%parallel_mode == SINGLE_PROC) then
      elsi_h%n_l_rows    = header(1)
      elsi_h%n_l_cols    = header(1)
      elsi_h%n_b_rows    = header(1)
      elsi_h%n_b_cols    = header(1)
   endif

   ! Compute n_l_col
   elsi_h%n_l_cols_sp = elsi_h%n_basis/elsi_h%n_procs
   n_l_cols0 = elsi_h%n_l_cols_sp
   if(elsi_h%myid == elsi_h%n_procs-1) then
      elsi_h%n_l_cols_sp = elsi_h%n_basis-(elsi_h%n_procs-1)*n_l_cols0
   endif

   call elsi_allocate(elsi_h,elsi_h%col_ptr_sips,elsi_h%n_l_cols_sp+1,&
           "col_ptr_sips",caller)

   ! Read column pointer
   offset = HEADER_SIZE+elsi_h%myid*n_l_cols0*4

   call MPI_File_read_at_all(filehandle,offset,elsi_h%col_ptr_sips,&
           elsi_h%n_l_cols_sp+1,mpi_integer4,mpi_status_ignore,mpierr)

   if(elsi_h%myid == elsi_h%n_procs-1) then
      elsi_h%col_ptr_sips(elsi_h%n_l_cols_sp+1) = elsi_h%nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = elsi_h%col_ptr_sips(1)-1
   elsi_h%col_ptr_sips = elsi_h%col_ptr_sips-prev_nnz

   ! Compute nnz_l
   elsi_h%nnz_l_sp = elsi_h%col_ptr_sips(elsi_h%n_l_cols_sp+1)-&
                        elsi_h%col_ptr_sips(1)

   call elsi_allocate(elsi_h,elsi_h%row_ind_sips,elsi_h%nnz_l_sp,&
           "row_ind_sips",caller)

   ! Read row index
   offset = HEADER_SIZE+elsi_h%n_basis*4+prev_nnz*4

   call MPI_File_read_at_all(filehandle,offset,elsi_h%row_ind_sips,&
           elsi_h%nnz_l_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Set ELSI sparsity pattern up
   call elsi_set_row_ind(elsi_h,elsi_h%row_ind_sips)
   call elsi_set_col_ptr(elsi_h,elsi_h%col_ptr_sips)

   elsi_h%sparsity_ready = .true.

   ! Read non-zero value
   offset = HEADER_SIZE+elsi_h%n_basis*4+elsi_h%nnz_g*4+prev_nnz*8

   select case(id)
   case(0)
      call elsi_allocate(elsi_h,elsi_h%ham_real_sips,elsi_h%nnz_l_sp,&
              "ham_real_sips",caller)

      call MPI_File_read_at_all(filehandle,offset,elsi_h%ham_real_sips,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)

      call elsi_set_sparse_ham(elsi_h,elsi_h%ham_real_sips)
   case(1)
      call elsi_allocate(elsi_h,elsi_h%ovlp_real_sips,elsi_h%nnz_l_sp,&
              "ovlp_real_sips",caller)

      call MPI_File_read_at_all(filehandle,offset,elsi_h%ovlp_real_sips,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)

      call elsi_set_sparse_ovlp(elsi_h,elsi_h%ovlp_real_sips)
   case(2)
      call elsi_allocate(elsi_h,elsi_h%dm_real_sips,elsi_h%nnz_l_sp,&
              "dm_real_sips",caller)

      call MPI_File_read_at_all(filehandle,offset,elsi_h%dm_real_sips,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)

      call elsi_set_sparse_dm(elsi_h,elsi_h%dm_real_sips)
   case default
      call elsi_stop(" Matrix not supported. Exiting...",elsi_h,caller)
   end select

   call MPI_File_close(filehandle,mpierr)

end subroutine

!>
!! This routine writes a 1D block distributed CSC matrix to file.
!!
subroutine elsi_write_matrix_real_sparse(elsi_h,filename,id)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h   !< Handle
   character(*),      intent(in)    :: filename !< File name
   integer(kind=i4),  intent(in)    :: id       !< H,S,D,other

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: filehandle
   integer(kind=i4) :: filemode
   integer(kind=i4) :: header(3)
   integer(kind=i4) :: prev_nnz
   integer(kind=i4) :: n_l_cols0

   integer(kind=mpi_offset_kind) :: offset

   character*40, parameter :: caller = "elsi_write_matrix_real_sparse"

   ! Open file
   filemode = mpi_mode_wronly+mpi_mode_create

   call MPI_File_open(elsi_h%mpi_comm,filename,filemode,mpi_info_null,&
           filehandle,mpierr)

   ! Write header
   header(1) = elsi_h%n_basis
   header(2) = int(elsi_h%n_electrons,kind=i4)
   header(3) = elsi_h%nnz_g

   if(elsi_h%myid == 0) then
      offset = 0

      call MPI_File_write_at(filehandle,offset,header,3,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Compute shift of column pointers
   prev_nnz = 0

   call MPI_Exscan(elsi_h%nnz_l_sp,prev_nnz,1,mpi_integer4,mpi_sum,&
           elsi_h%mpi_comm,mpierr)

   ! Shift column pointer
   elsi_h%col_ptr_ccs = elsi_h%col_ptr_ccs+prev_nnz

   ! Write column pointer
   n_l_cols0 = elsi_h%n_basis/elsi_h%n_procs
   offset = HEADER_SIZE+elsi_h%myid*n_l_cols0*4

   call MPI_File_write_at_all(filehandle,offset,elsi_h%col_ptr_ccs,&
           elsi_h%n_l_cols_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Unshift column pointer
   elsi_h%col_ptr_ccs = elsi_h%col_ptr_ccs-prev_nnz

   ! Write row index
   offset = HEADER_SIZE+elsi_h%n_basis*4+prev_nnz*4

   call MPI_File_write_at_all(filehandle,offset,elsi_h%row_ind_ccs,&
           elsi_h%nnz_l_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Write non-zero value
   offset = HEADER_SIZE+elsi_h%n_basis*4+elsi_h%nnz_g*4+prev_nnz*8

   select case(id)
   case(0)
      call MPI_File_write_at_all(filehandle,offset,elsi_h%ham_real_ccs,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)
   case(1)
      call MPI_File_write_at_all(filehandle,offset,elsi_h%ovlp_real_ccs,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)
   case(2)
      call MPI_File_write_at_all(filehandle,offset,elsi_h%dm_real_ccs,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)
   case default
      call elsi_stop(" Matrix not supported. Exiting...",elsi_h,caller)
   end select

   call MPI_File_close(filehandle,mpierr)

end subroutine

end module ELSI_IO
