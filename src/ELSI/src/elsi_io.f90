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
   use ELSI_CONSTANTS, only: SINGLE_PROC,HEADER_SIZE,BLACS_DENSE,PEXSI_CSC,&
                             MATRIX_H,MATRIX_S,MATRIX_D
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_MATCONV
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS

   implicit none

   private

   public :: elsi_write_matrix_real
   public :: elsi_write_matrix_real_sparse
!   public :: elsi_write_matrix_complex
!   public :: elsi_write_matrix_complex_sparse
   public :: elsi_read_matrix_real
   public :: elsi_read_matrix_real_sparse
!   public :: elsi_read_matrix_complex
!   public :: elsi_read_matrix_complex_sparse
   public :: elsi_get_matrix_dim
   public :: elsi_get_matrix_dim_sparse
   public :: elsi_get_csc
   public :: elsi_get_matrix_real
   public :: elsi_get_matrix_real_sparse
   public :: elsi_get_matrix_complex
   public :: elsi_get_matrix_complex_sparse

contains

!>
!! This routine reads a 2D block-cyclic dense matrix from file.
!!
subroutine elsi_read_matrix_real(elsi_h,id)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   integer(kind=i4),  intent(in)    :: id     !< H,S,D

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: filehandle
   integer(kind=i4) :: filemode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_l_cols
   integer(kind=i4) :: n_l_cols0
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: prev_nnz

   integer(kind=mpi_offset_kind) :: offset
   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)
   real(kind=r8),    allocatable :: nnz_val(:)

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: k_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: local_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: local_row_id ! Local row id in 1D block distribution
   integer(kind=i4) :: proc_col_id ! Column id in process grid
   integer(kind=i4) :: proc_row_id ! Row id in process grid
   integer(kind=i4), allocatable :: global_col_id(:) ! Global column id
   integer(kind=i4), allocatable :: global_row_id(:) ! Global row id
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: val_send_buffer(:)
   integer(kind=i4), allocatable :: row_send_buffer(:)
   integer(kind=i4), allocatable :: col_send_buffer(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: val_recv_buffer(:)
   integer(kind=i4), allocatable :: row_recv_buffer(:)
   integer(kind=i4), allocatable :: col_recv_buffer(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4) :: send_displ_aux
   integer(kind=i4) :: recv_displ_aux

   character    :: label
   character*40 :: filename

   character*40, parameter :: caller = "elsi_read_matrix_real"

   call elsi_check_handle(elsi_h,caller)

   if(.not. elsi_h%mpi_ready) then
      call elsi_stop(" Parallel matrix I/O requires MPI. Exiting...",&
              elsi_h,caller)
   endif

   ! File name
   select case(id)
   case(MATRIX_H)
      label = "H"
   case(MATRIX_S)
      label = "S"
   case(MATRIX_D)
      label = "D"
   case default
      call elsi_stop(" Matrix not supported. Exiting...",elsi_h,caller)
   end select

   write(filename,"(A,A,I1.1,A,I5.5,A)") label,"_spin_",elsi_h%i_spin,&
      "_kpt_",elsi_h%i_kpt,".csc"

   ! Open file
   filemode = mpi_mode_rdonly

   call MPI_File_open(elsi_h%mpi_comm,filename,filemode,mpi_info_null,&
           filehandle,mpierr)

   ! Read header
   if(elsi_h%myid == 0) then
      offset = 0

      call MPI_File_read_at(filehandle,offset,header,HEADER_SIZE,&
              mpi_integer4,mpi_status_ignore,mpierr)
   endif

   ! Broadcast header
   call MPI_Bcast(header,HEADER_SIZE,mpi_integer4,0,elsi_h%mpi_comm,mpierr)

   ! Set ELSI parameters
   if(header(1) /= BLACS_DENSE) then
      call elsi_stop(" Not reading a dense matrix. Exiting...",elsi_h,caller)
   endif

   elsi_h%matrix_format = BLACS_DENSE
   elsi_h%n_basis       = header(2)
   elsi_h%n_nonsing     = header(2)
   elsi_h%n_electrons   = real(header(3),kind=r8)
   elsi_h%nnz_g         = header(4)

   ! Compute n_l_col
   n_l_cols  = elsi_h%n_basis/elsi_h%n_procs
   n_l_cols0 = n_l_cols
   if(elsi_h%myid == elsi_h%n_procs-1) then
      elsi_h%n_l_cols = elsi_h%n_basis-(elsi_h%n_procs-1)*n_l_cols0
   endif

   call elsi_allocate(elsi_h,col_ptr,n_l_cols+1,"col_ptr",caller)

   ! Read column pointer
   offset = HEADER_SIZE*4+elsi_h%myid*n_l_cols0*4

   call MPI_File_read_at_all(filehandle,offset,col_ptr,n_l_cols+1,&
           mpi_integer4,mpi_status_ignore,mpierr)

   if(elsi_h%myid == elsi_h%n_procs-1) then
      col_ptr(n_l_cols+1) = elsi_h%nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = col_ptr(1)-1
   col_ptr  = col_ptr-prev_nnz

   ! Compute nnz_l
   nnz_l = col_ptr(n_l_cols+1)-col_ptr(1)

   call elsi_allocate(elsi_h,row_ind,nnz_l,"row_ind",caller)

   ! Read row index
   offset = HEADER_SIZE*4+elsi_h%n_basis*4+prev_nnz*4

   call MPI_File_read_at_all(filehandle,offset,row_ind,nnz_l,mpi_integer4,&
           mpi_status_ignore,mpierr)

   ! Read non-zero value
   offset = HEADER_SIZE*4+elsi_h%n_basis*4+elsi_h%nnz_g*4+prev_nnz*8

   call elsi_allocate(elsi_h,nnz_val,nnz_l,"nnz_val",caller)

   call MPI_File_read_at_all(filehandle,offset,nnz_val,nnz_l,mpi_real8,&
           mpi_status_ignore,mpierr)

   ! Close file
   call MPI_File_close(filehandle,mpierr)

   ! Redistribute matrix
   call elsi_allocate(elsi_h,val_send_buffer,nnz_l,"val_send_buffer",caller)
   call elsi_allocate(elsi_h,row_send_buffer,nnz_l,"row_send_buffer",caller)
   call elsi_allocate(elsi_h,col_send_buffer,nnz_l,"col_send_buffer",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)
   call elsi_allocate(elsi_h,global_row_id,nnz_l,"global_row_id",caller)
   call elsi_allocate(elsi_h,global_col_id,nnz_l,"global_col_id",caller)
   call elsi_allocate(elsi_h,dest,nnz_l,"dest",caller)

   ! Compute destination and global 1D id
   i_col = 0

   do i_val = 1,nnz_l
      if((i_val == col_ptr(i_col+1)) .and. (i_col /= n_l_cols)) then
         i_col = i_col+1
      endif
      i_row = row_ind(i_val)

      ! Compute global id
      global_row_id = i_row
      global_col_id = i_col+elsi_h%myid*n_l_cols0

      ! Compute destination
      proc_row_id = mod((global_row_id-1)/elsi_h%n_b_rows,elsi_h%n_p_rows)
      proc_col_id = mod((global_col_id-1)/elsi_h%n_b_cols,elsi_h%n_p_cols)
      dest(i_val) = proc_col_id+proc_row_id*elsi_h%n_p_cols
   enddo

   call elsi_deallocate(elsi_h,row_ind,"row_ind")
   call elsi_deallocate(elsi_h,col_ptr,"col_ptr")

   j_val = 0
   k_val = nnz_l+1

   ! Set send_count
   do i_proc = 1,elsi_h%n_procs
      do i_val = 1,nnz_l
         if(dest(i_val) == i_proc-1) then
            j_val = j_val+1
            val_send_buffer(j_val) = nnz_val(i_val)
            row_send_buffer(j_val) = global_row_id(i_val)
            col_send_buffer(j_val) = global_col_id(i_val)
            send_count(i_proc) = send_count(i_proc)+1
         endif
      enddo
   enddo

   call elsi_deallocate(elsi_h,nnz_val,"nnz_val")
   call elsi_deallocate(elsi_h,global_row_id,"global_row_id")
   call elsi_deallocate(elsi_h,global_col_id,"global_col_id")
   call elsi_deallocate(elsi_h,dest,"dest")

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           elsi_h%mpi_comm,mpierr)

   elsi_h%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   send_displ_aux = 0
   recv_displ_aux = 0

   do i_proc = 1,elsi_h%n_procs
      send_displ(i_proc) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc)

      recv_displ(i_proc) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc)
   enddo

   ! Send and receive the packed data
   ! Value
   call elsi_allocate(elsi_h,val_recv_buffer,elsi_h%nnz_l,"val_recv_buffer",&
           caller)

   call MPI_Alltoallv(val_send_buffer,send_count,send_displ,mpi_real8,&
           val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,&
           mpierr)

   call elsi_deallocate(elsi_h,val_send_buffer,"val_send_buffer")

   ! Row ID
   call elsi_allocate(elsi_h,row_recv_buffer,elsi_h%nnz_l,"row_recv_buffer",&
           caller)

   call MPI_Alltoallv(row_send_buffer,send_count,send_displ,mpi_integer4,&
           row_recv_buffer,recv_count,recv_displ,mpi_integer4,&
           elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buffer,"row_send_buffer")

   ! Column ID
   call elsi_allocate(elsi_h,col_recv_buffer,elsi_h%nnz_l,"col_recv_buffer",&
           caller)

   call MPI_Alltoallv(col_send_buffer,send_count,send_displ,mpi_integer4,&
           col_recv_buffer,recv_count,recv_displ,mpi_integer4,&
           elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buffer,"col_send_buffer")
   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   ! Unpack matrix
   select case(id)
   case(MATRIX_H)
      call elsi_allocate(elsi_h,elsi_h%ham_real_elpa,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,"ham_real_elpa",caller)

      do i_val = 1,elsi_h%nnz_l
         ! Compute local 2d id
         local_row_id = (row_recv_buffer(i_val)-1)/&
                           (elsi_h%n_p_rows*elsi_h%n_b_rows)*elsi_h%n_b_rows+&
                           mod((row_recv_buffer(i_val)-1),elsi_h%n_b_rows)+1
         local_col_id = (col_recv_buffer(i_val)-1)/&
                        (elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols+&
                        mod((col_recv_buffer(i_val)-1),elsi_h%n_b_cols)+1

         ! Put value to correct position
         elsi_h%ham_real_elpa(local_row_id,local_col_id) = val_recv_buffer(i_val)
      enddo

      call elsi_set_ham(elsi_h,elsi_h%ham_real_elpa)
   case(MATRIX_S)
      call elsi_allocate(elsi_h,elsi_h%ovlp_real_elpa,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,"ovlp_real_elpa",caller)

      do i_val = 1,elsi_h%nnz_l
         ! Compute local 2d id
         local_row_id = (row_recv_buffer(i_val)-1)/&
                           (elsi_h%n_p_rows*elsi_h%n_b_rows)*elsi_h%n_b_rows+&
                           mod((row_recv_buffer(i_val)-1),elsi_h%n_b_rows)+1
         local_col_id = (col_recv_buffer(i_val)-1)/&
                        (elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols+&
                        mod((col_recv_buffer(i_val)-1),elsi_h%n_b_cols)+1

         ! Put value to correct position
         elsi_h%ovlp_real_elpa(local_row_id,local_col_id) = val_recv_buffer(i_val)
      enddo

      call elsi_set_ovlp(elsi_h,elsi_h%ovlp_real_elpa)
   case(MATRIX_D)
      call elsi_allocate(elsi_h,elsi_h%dm_real_elpa,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,"dm_real_elpa",caller)

      do i_val = 1,elsi_h%nnz_l
         ! Compute local 2d id
         local_row_id = (row_recv_buffer(i_val)-1)/&
                           (elsi_h%n_p_rows*elsi_h%n_b_rows)*elsi_h%n_b_rows+&
                           mod((row_recv_buffer(i_val)-1),elsi_h%n_b_rows)+1
         local_col_id = (col_recv_buffer(i_val)-1)/&
                        (elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols+&
                        mod((col_recv_buffer(i_val)-1),elsi_h%n_b_cols)+1

         ! Put value to correct position
         elsi_h%dm_real_elpa(local_row_id,local_col_id) = val_recv_buffer(i_val)
      enddo

      call elsi_set_dm(elsi_h,elsi_h%dm_real_elpa)
   end select

   call elsi_deallocate(elsi_h,val_recv_buffer,"val_recv_buffer")
   call elsi_deallocate(elsi_h,row_recv_buffer,"row_recv_buffer")
   call elsi_deallocate(elsi_h,col_recv_buffer,"col_recv_buffer")

end subroutine

!>
!! This routine reads a 1D block CSC matrix from file.
!!
subroutine elsi_read_matrix_real_sparse(elsi_h,id)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   integer(kind=i4),  intent(in)    :: id     !< H,S,D

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: filehandle
   integer(kind=i4) :: filemode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_l_cols0
   integer(kind=i4) :: prev_nnz

   integer(kind=mpi_offset_kind) :: offset

   character    :: label
   character*40 :: filename

   character*40, parameter :: caller = "elsi_read_matrix_real_sparse"

   call elsi_check_handle(elsi_h,caller)

   if(.not. elsi_h%mpi_ready) then
      call elsi_stop(" Parallel matrix I/O requires MPI. Exiting...",&
              elsi_h,caller)
   endif

   ! File name
   select case(id)
   case(MATRIX_H)
      label = "H"
   case(MATRIX_S)
      label = "S"
   case(MATRIX_D)
      label = "D"
   case default
      call elsi_stop(" Matrix not supported. Exiting...",elsi_h,caller)
   end select

   write(filename,"(A,A,I1.1,A,I5.5,A)") label,"_spin_",elsi_h%i_spin,&
      "_kpt_",elsi_h%i_kpt,".csc"

   ! Open file
   filemode = mpi_mode_rdonly

   call MPI_File_open(elsi_h%mpi_comm,filename,filemode,mpi_info_null,&
           filehandle,mpierr)

   ! Read header
   if(elsi_h%myid == 0) then
      offset = 0

      call MPI_File_read_at(filehandle,offset,header,HEADER_SIZE,&
              mpi_integer4,mpi_status_ignore,mpierr)
   endif

   ! Broadcast header
   call MPI_Bcast(header,HEADER_SIZE,mpi_integer4,0,elsi_h%mpi_comm,mpierr)

   ! Set ELSI parameters
   if(header(1) /= PEXSI_CSC) then
      call elsi_stop(" Not reading a sparse matrix. Exiting...",elsi_h,caller)
   endif

   elsi_h%matrix_format = PEXSI_CSC
   elsi_h%n_basis       = header(2)
   elsi_h%n_nonsing     = header(2)
   elsi_h%n_electrons   = real(header(3),kind=r8)
   elsi_h%nnz_g         = header(4)

   ! Compute n_l_col
   elsi_h%n_l_cols_sp = elsi_h%n_basis/elsi_h%n_procs
   n_l_cols0 = elsi_h%n_l_cols_sp
   if(elsi_h%myid == elsi_h%n_procs-1) then
      elsi_h%n_l_cols_sp = elsi_h%n_basis-(elsi_h%n_procs-1)*n_l_cols0
   endif

   call elsi_allocate(elsi_h,elsi_h%col_ptr_sips,elsi_h%n_l_cols_sp+1,&
           "col_ptr_sips",caller)

   ! Read column pointer
   offset = HEADER_SIZE*4+elsi_h%myid*n_l_cols0*4

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
   offset = HEADER_SIZE*4+elsi_h%n_basis*4+prev_nnz*4

   call MPI_File_read_at_all(filehandle,offset,elsi_h%row_ind_sips,&
           elsi_h%nnz_l_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Set ELSI sparsity pattern up
   call elsi_set_row_ind(elsi_h,elsi_h%row_ind_sips)
   call elsi_set_col_ptr(elsi_h,elsi_h%col_ptr_sips)

   elsi_h%sparsity_ready = .true.

   ! Read non-zero value
   offset = HEADER_SIZE*4+elsi_h%n_basis*4+elsi_h%nnz_g*4+prev_nnz*8

   select case(id)
   case(MATRIX_H)
      call elsi_allocate(elsi_h,elsi_h%ham_real_sips,elsi_h%nnz_l_sp,&
              "ham_real_sips",caller)

      call MPI_File_read_at_all(filehandle,offset,elsi_h%ham_real_sips,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)

      call elsi_set_sparse_ham(elsi_h,elsi_h%ham_real_sips)
   case(MATRIX_S)
      call elsi_allocate(elsi_h,elsi_h%ovlp_real_sips,elsi_h%nnz_l_sp,&
              "ovlp_real_sips",caller)

      call MPI_File_read_at_all(filehandle,offset,elsi_h%ovlp_real_sips,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)

      call elsi_set_sparse_ovlp(elsi_h,elsi_h%ovlp_real_sips)
   case(MATRIX_D)
      call elsi_allocate(elsi_h,elsi_h%dm_real_sips,elsi_h%nnz_l_sp,&
              "dm_real_sips",caller)

      call MPI_File_read_at_all(filehandle,offset,elsi_h%dm_real_sips,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)

      call elsi_set_sparse_dm(elsi_h,elsi_h%dm_real_sips)
   end select

   call MPI_File_close(filehandle,mpierr)

end subroutine

!>
!! This routine writes a 2D block-cyclic dense matrix to file.
!!
subroutine elsi_write_matrix_real(elsi_h,id)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   integer(kind=i4),  intent(in)    :: id     !< H,S,D

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: filehandle
   integer(kind=i4) :: filemode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_l_cols
   integer(kind=i4) :: n_l_cols0
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: prev_nnz

   integer(kind=mpi_offset_kind) :: offset
   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)
   real(kind=r8),    allocatable :: nnz_val(:)

   integer(kind=mpi_offset_kind) :: offset

   character    :: label
   character*40 :: filename

   character*40, parameter :: caller = "elsi_write_matrix_real"

   call elsi_check_handle(elsi_h,caller)

   if(.not. elsi_h%mpi_ready) then
      call elsi_stop(" Parallel matrix I/O requires MPI. Exiting...",&
              elsi_h,caller)
   endif

   ! Redistribute matrix
   ! TODO

   ! File name
   select case(id)
   case(MATRIX_H)
      label = "H"
   case(MATRIX_S)
      label = "S"
   case(MATRIX_D)
      label = "D"
   case default
      call elsi_stop(" Matrix not supported. Exiting...",elsi_h,caller)
   end select

   write(filename,"(A,A,I1.1,A,I5.5,A)") label,"_spin_",elsi_h%i_spin,&
      "_kpt_",elsi_h%i_kpt,".csc"

   ! Open file
   filemode = mpi_mode_wronly+mpi_mode_create

   call MPI_File_open(elsi_h%mpi_comm,filename,filemode,mpi_info_null,&
           filehandle,mpierr)

   ! Write header
   header(1) = BLACS_DENSE
   header(2) = elsi_h%n_basis
   header(3) = int(elsi_h%n_electrons,kind=i4)
   header(4) = elsi_h%nnz_g

   if(elsi_h%myid == 0) then
      offset = 0

      call MPI_File_write_at(filehandle,offset,header,HEADER_SIZE,&
              mpi_integer4,mpi_status_ignore,mpierr)
   endif

   ! Compute shift of column pointers
   prev_nnz = 0

   call MPI_Exscan(nnz_l,prev_nnz,1,mpi_integer4,mpi_sum,elsi_h%mpi_comm,&
           mpierr)

   ! Shift column pointer
   col_ptr = col_ptr+prev_nnz

   ! Write column pointer
   n_l_cols0 = elsi_h%n_basis/elsi_h%n_procs
   offset = HEADER_SIZE*4+elsi_h%myid*n_l_cols0*4

   call MPI_File_write_at_all(filehandle,offset,col_ptr,n_l_cols,&
           mpi_integer4,mpi_status_ignore,mpierr)

   ! Unshift column pointer
   col_ptr = col_ptr-prev_nnz

   ! Write row index
   offset = HEADER_SIZE*4+elsi_h%n_basis*4+prev_nnz*4

   call MPI_File_write_at_all(filehandle,offset,row_ind,nnz_l,&
           mpi_integer4,mpi_status_ignore,mpierr)

   ! Write non-zero value
   offset = HEADER_SIZE*4+elsi_h%n_basis*4+elsi_h%nnz_g*4+prev_nnz*8

   call MPI_File_write_at_all(filehandle,offset,nnz_val,nnz_l,mpi_real8,&
           mpi_status_ignore,mpierr)

   call MPI_File_close(filehandle,mpierr)

end subroutine

!>
!! This routine writes a 1D block CSC matrix to file.
!!
subroutine elsi_write_matrix_real_sparse(elsi_h,id)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   integer(kind=i4),  intent(in)    :: id     !< H,S,D

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: filehandle
   integer(kind=i4) :: filemode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: prev_nnz
   integer(kind=i4) :: n_l_cols0

   integer(kind=mpi_offset_kind) :: offset

   character    :: label
   character*40 :: filename

   character*40, parameter :: caller = "elsi_write_matrix_real_sparse"

   call elsi_check_handle(elsi_h,caller)

   if(.not. elsi_h%mpi_ready) then
      call elsi_stop(" Parallel matrix I/O requires MPI. Exiting...",&
              elsi_h,caller)
   endif

   ! File name
   select case(id)
   case(MATRIX_H)
      label = "H"
   case(MATRIX_S)
      label = "S"
   case(MATRIX_D)
      label = "D"
   case default
      call elsi_stop(" Matrix not supported. Exiting...",elsi_h,caller)
   end select

   write(filename,"(A,A,I1.1,A,I5.5,A)") label,"_spin_",elsi_h%i_spin,&
      "_kpt_",elsi_h%i_kpt,".csc"

   ! Open file
   filemode = mpi_mode_wronly+mpi_mode_create

   call MPI_File_open(elsi_h%mpi_comm,filename,filemode,mpi_info_null,&
           filehandle,mpierr)

   ! Write header
   header(1) = PEXSI_CSC
   header(2) = elsi_h%n_basis
   header(3) = int(elsi_h%n_electrons,kind=i4)
   header(4) = elsi_h%nnz_g

   if(elsi_h%myid == 0) then
      offset = 0

      call MPI_File_write_at(filehandle,offset,header,HEADER_SIZE,&
              mpi_integer4,mpi_status_ignore,mpierr)
   endif

   ! Compute shift of column pointers
   prev_nnz = 0

   call MPI_Exscan(elsi_h%nnz_l_sp,prev_nnz,1,mpi_integer4,mpi_sum,&
           elsi_h%mpi_comm,mpierr)

   ! Shift column pointer
   elsi_h%col_ptr_ccs = elsi_h%col_ptr_ccs+prev_nnz

   ! Write column pointer
   n_l_cols0 = elsi_h%n_basis/elsi_h%n_procs
   offset = HEADER_SIZE*4+elsi_h%myid*n_l_cols0*4

   call MPI_File_write_at_all(filehandle,offset,elsi_h%col_ptr_ccs,&
           elsi_h%n_l_cols_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Unshift column pointer
   elsi_h%col_ptr_ccs = elsi_h%col_ptr_ccs-prev_nnz

   ! Write row index
   offset = HEADER_SIZE*4+elsi_h%n_basis*4+prev_nnz*4

   call MPI_File_write_at_all(filehandle,offset,elsi_h%row_ind_ccs,&
           elsi_h%nnz_l_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Write non-zero value
   offset = HEADER_SIZE*4+elsi_h%n_basis*4+elsi_h%nnz_g*4+prev_nnz*8

   select case(id)
   case(MATRIX_H)
      if(.not. associated(elsi_h%ham_real_ccs)) then
         call elsi_stop(" Sparse Hamiltonian not available. Exiting...",&
                 elsi_h,caller)
      endif

      call MPI_File_write_at_all(filehandle,offset,elsi_h%ham_real_ccs,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)
   case(MATRIX_S)
      if(.not. associated(elsi_h%ovlp_real_ccs)) then
         call elsi_stop(" Sparse overlap not available. Exiting...",&
                 elsi_h,caller)
      endif

      call MPI_File_write_at_all(filehandle,offset,elsi_h%ovlp_real_ccs,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)
   case(MATRIX_D)
      if(.not. associated(elsi_h%dm_real_ccs)) then
         call elsi_stop(" Sparse density matrix not available. Exiting...",&
                 elsi_h,caller)
      endif

      call MPI_File_write_at_all(filehandle,offset,elsi_h%dm_real_ccs,&
              elsi_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)
   end select

   call MPI_File_close(filehandle,mpierr)

end subroutine

!>
!! This routine gets the dimensions of 2D block-cyclic dense matrices in ELSI.
!!
subroutine elsi_get_matrix_dim(elsi_h,n_basis,n_row_l,n_col_l)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h  !< Handle
   integer(kind=i4),  intent(out) :: n_basis !< Number of basis functions
   integer(kind=i4),  intent(out) :: n_row_l !< Number of local rows
   integer(kind=i4),  intent(out) :: n_col_l !< Number of local columns

   character*40, parameter :: caller = "elsi_get_matrix_dim"

   call elsi_check_handle(elsi_h,caller)

   n_basis = elsi_h%n_basis
   n_row_l = elsi_h%n_l_rows
   n_col_l = elsi_h%n_l_cols

end subroutine

!>
!! This routine gets the dimensions of 1D block CSC matrices in ELSI.
!!
subroutine elsi_get_matrix_dim_sparse(elsi_h,n_basis,nnz_g,nnz_l,n_col_l)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h  !< Handle
   integer(kind=i4),  intent(out) :: n_basis !< Number of basis functions
   integer(kind=i4),  intent(out) :: nnz_g   !< Global number of nonzeros
   integer(kind=i4),  intent(out) :: nnz_l   !< Local number of nonzeros
   integer(kind=i4),  intent(out) :: n_col_l !< Number of local columns

   character*40, parameter :: caller = "elsi_get_matrix_dim_sparse"

   call elsi_check_handle(elsi_h,caller)

   n_basis = elsi_h%n_basis
   nnz_g   = elsi_h%nnz_g
   nnz_l   = elsi_h%nnz_l_sp
   n_col_l = elsi_h%n_l_cols_sp

end subroutine

!>
!! This routine gets the row index and column pointer arrays of 1D block
!! CSC matrices in ELSI.
!!
subroutine elsi_get_csc(elsi_h,row_ind,col_ptr)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h                        !< Handle
   integer(kind=i4),  intent(out) :: row_ind(elsi_h%nnz_l_sp)      !< Row index
   integer(kind=i4),  intent(out) :: col_ptr(elsi_h%n_l_cols_sp+1) !< Column pointer

   character*40, parameter :: caller = "elsi_get_csc"

   call elsi_check_handle(elsi_h,caller)

   if(.not. associated(elsi_h%row_ind_ccs)) then
      call elsi_stop(" Row index not available. Exiting...",&
              elsi_h,caller)
   endif

   if(.not. associated(elsi_h%col_ptr_ccs)) then
      call elsi_stop(" Column pointer not available. Exiting...",&
              elsi_h,caller)
   endif

   row_ind = elsi_h%row_ind_ccs
   col_ptr = elsi_h%col_ptr_ccs

end subroutine

!>
!! This routine gets a 2D block-cyclic dense matrix.
!!
subroutine elsi_get_matrix_real(elsi_h,id,matrix_out)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h                                      !< Handle
   integer(kind=i4),  intent(in)  :: id                                          !< H,S,D
   real(kind=r8),     intent(out) :: matrix_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Output matrix

   character*40, parameter :: caller = "elsi_get_matrix_real"

   call elsi_check_handle(elsi_h,caller)

   select case(id)
   case(MATRIX_H)
      if(.not. associated(elsi_h%ham_real)) then
         call elsi_stop(" Hamiltonian not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%ham_real
   case(MATRIX_S)
      if(.not. associated(elsi_h%ovlp_real)) then
         call elsi_stop(" Overlap not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%ovlp_real
   case(MATRIX_D)
      if(.not. associated(elsi_h%dm_real)) then
         call elsi_stop(" Density matrix not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%dm_real
   case default
      call elsi_stop(" Matrix not supported. Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine gets a 1D block CSC matrix.
!!
subroutine elsi_get_matrix_real_sparse(elsi_h,id,matrix_out)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h                      !< Handle
   integer(kind=i4),  intent(in)  :: id                          !< H,S,D
   real(kind=r8),     intent(out) :: matrix_out(elsi_h%nnz_l_sp) !< Output matrix

   character*40, parameter :: caller = "elsi_get_matrix_real_sparse"

   call elsi_check_handle(elsi_h,caller)

   select case(id)
   case(MATRIX_H)
      if(.not. associated(elsi_h%ham_real_ccs)) then
         call elsi_stop(" Sparse Hamiltonian not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%ham_real_ccs
   case(MATRIX_S)
      if(.not. associated(elsi_h%ovlp_real_ccs)) then
         call elsi_stop(" Sparse overlap not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%ovlp_real_ccs
   case(MATRIX_D)
      if(.not. associated(elsi_h%dm_real_ccs)) then
         call elsi_stop(" Sparse density matrix not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%dm_real_ccs
   case default
      call elsi_stop(" Matrix not supported. Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine gets a 2D block-cyclic dense matrix.
!!
subroutine elsi_get_matrix_complex(elsi_h,id,matrix_out)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h                                      !< Handle
   integer(kind=i4),  intent(in)  :: id                                          !< H,S,D
   complex(kind=r8),  intent(out) :: matrix_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Output matrix

   character*40, parameter :: caller = "elsi_get_matrix_complex"

   call elsi_check_handle(elsi_h,caller)

   select case(id)
   case(MATRIX_H)
      if(.not. associated(elsi_h%ham_real)) then
         call elsi_stop(" Hamiltonian not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%ham_complex
   case(MATRIX_S)
      if(.not. associated(elsi_h%ovlp_real)) then
         call elsi_stop(" Overlap not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%ovlp_complex
   case(MATRIX_D)
      if(.not. associated(elsi_h%dm_real)) then
         call elsi_stop(" Density matrix not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%dm_complex
   case default
      call elsi_stop(" Matrix not supported. Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine gets a 1D block CSC matrix.
!!
subroutine elsi_get_matrix_complex_sparse(elsi_h,id,matrix_out)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h                      !< Handle
   integer(kind=i4),  intent(in)  :: id                          !< H,S,D
   complex(kind=i4),  intent(out) :: matrix_out(elsi_h%nnz_l_sp) !< Output matrix

   character*40, parameter :: caller = "elsi_get_matrix_complex_sparse"

   call elsi_check_handle(elsi_h,caller)

   select case(id)
   case(MATRIX_H)
      if(.not. associated(elsi_h%ham_real_ccs)) then
         call elsi_stop(" Sparse Hamiltonian not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%ham_complex_ccs
   case(MATRIX_S)
      if(.not. associated(elsi_h%ovlp_real_ccs)) then
         call elsi_stop(" Sparse overlap not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%ovlp_complex_ccs
   case(MATRIX_D)
      if(.not. associated(elsi_h%dm_real_ccs)) then
         call elsi_stop(" Sparse density matrix not available. Exiting...",&
                 elsi_h,caller)
      endif

      matrix_out = elsi_h%dm_complex_ccs
   case default
      call elsi_stop(" Matrix not supported. Exiting...",elsi_h,caller)
   end select

end subroutine

end module ELSI_IO
