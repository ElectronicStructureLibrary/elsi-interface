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
!! This module provides matrix conversion and redistribution routines.
!!
module ELSI_MATCONV

   use ELSI_CONSTANTS, only: ELPA,LIBOMM
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4,i8
   use ELSI_UTILS

   implicit none

   private

   public :: elsi_blacs_to_chess_hs
   public :: elsi_blacs_to_sips_dm
   public :: elsi_blacs_to_pexsi_hs
   public :: elsi_blacs_to_sips_hs
   public :: elsi_chess_to_blacs_dm
   public :: elsi_pexsi_to_blacs_dm
   public :: elsi_sips_to_blacs_hs

   interface elsi_blacs_to_chess_hs
      module procedure elsi_blacs_to_chess_hs_real
   end interface

   interface elsi_blacs_to_sips_dm
      module procedure elsi_blacs_to_sips_dm_real
   end interface

   interface elsi_blacs_to_pexsi_hs
      module procedure elsi_blacs_to_pexsi_hs_real,&
                       elsi_blacs_to_pexsi_hs_complex
   end interface

   interface elsi_blacs_to_sips_hs
      module procedure elsi_blacs_to_sips_hs_real,&
                       elsi_blacs_to_sips_hs_complex
   end interface

   interface elsi_chess_to_blacs_dm
      module procedure elsi_chess_to_blacs_dm_real
   end interface

   interface elsi_pexsi_to_blacs_dm
      module procedure elsi_pexsi_to_blacs_dm_real,&
                       elsi_pexsi_to_blacs_dm_complex
   end interface

   interface elsi_sips_to_blacs_hs
      module procedure elsi_sips_to_blacs_hs_real
   end interface

   interface swap
      module procedure swap_i8,&
                       swap_i4,&
                       swap_r8,&
                       swap_c16
   end interface

   interface downheap
      module procedure downheap_real_v1,&
                       downheap_real_v2,&
                       downheap_complex_v1,&
                       downheap_complex_v2
   end interface

   interface heapsort
      module procedure heapsort_real_v1,&
                       heapsort_real_v2,&
                       heapsort_complex_v1,&
                       heapsort_complex_v2
   end interface

contains

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by PEXSI.
!!
subroutine elsi_blacs_to_pexsi_hs_real(elsi_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     target        :: h_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian
   real(kind=r8),     target        :: s_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap

   real(kind=r8), pointer :: ref(:,:) ! Sparsity pattern

   integer(kind=i4) :: n_para_task
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: global_col_id
   integer(kind=i4) :: global_row_id
   integer(kind=i4) :: d1,d2,d11,d12,d21,d22 ! Columns in the intermediate stage
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: this_n_cols
   integer(kind=i4) :: min_id
   integer(kind=i4) :: nnz_l_pexsi_aux
   integer(kind=i4) :: tmp_int
   integer(kind=i8) :: tmp_long
   real(kind=r8)    :: tmp_real
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send_buf(:)
   real(kind=r8),    allocatable :: s_val_send_buf(:)
   integer(kind=i4), allocatable :: row_send_buf(:)
   integer(kind=i4), allocatable :: col_send_buf(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: h_val_recv_buf(:)
   real(kind=r8),    allocatable :: s_val_recv_buf(:)
   integer(kind=i4), allocatable :: row_recv_buf(:)
   integer(kind=i4), allocatable :: col_recv_buf(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: locat(:) ! Location of each global column
   integer(kind=i8), allocatable :: global_id(:) ! Global 1D id

   character*40, parameter :: caller = "elsi_blacs_to_pexsi_hs_real"

   call elsi_get_time(elsi_h,t0)

   n_para_task = elsi_h%n_procs/elsi_h%n_p_per_pole

   if(.not. elsi_h%ovlp_is_unit) then
      ref => s_in
   else
      ref => h_in
   endif

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz(elsi_h,s_in,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,elsi_h%nnz_l)

      if(.not. elsi_h%ovlp_is_unit) then
         call elsi_allocate(elsi_h,s_val_send_buf,elsi_h%nnz_l,&
                 "s_val_send_buf",caller)
      endif
   endif

   call elsi_allocate(elsi_h,locat,elsi_h%n_basis,"locat",caller)
   call elsi_allocate(elsi_h,row_send_buf,elsi_h%nnz_l,&
           "row_send_buf",caller)
   call elsi_allocate(elsi_h,col_send_buf,elsi_h%nnz_l,&
           "col_send_buf",caller)
   call elsi_allocate(elsi_h,h_val_send_buf,elsi_h%nnz_l,&
           "h_val_send_buf",caller)

   ! Compute d1,d2,d11,d12,d21,d22 (need explanation)
   d1  = elsi_h%n_basis/elsi_h%n_p_per_pole
   d2  = elsi_h%n_basis-(elsi_h%n_p_per_pole-1)*d1
   d11 = d1/n_para_task
   d12 = d1-(n_para_task-1)*d11
   d21 = d2/n_para_task
   d22 = d2-(n_para_task-1)*d21

   i_val = 0

   if(d11 > 0) then
      do i_proc = 0,elsi_h%n_procs-n_para_task-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d12
         else
            this_n_cols = d11
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = 0,elsi_h%n_procs-n_para_task-1
         if(1+mod(i_proc,n_para_task) <= d1) then
            this_n_cols = 1
         else
            this_n_cols = 0
         endif

         if(this_n_cols /= 0) then
            locat(i_val+1:i_val+this_n_cols) = i_proc
            i_val = i_val+this_n_cols
         endif
      enddo
   endif

   if(d21 > 0) then
      do i_proc = elsi_h%n_procs-n_para_task,elsi_h%n_procs-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d22
         else
            this_n_cols = d21
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = elsi_h%n_procs-n_para_task,elsi_h%n_procs-1
         if(1+mod(i_proc,n_para_task) <= d2) then
            this_n_cols = 1
         else
            this_n_cols = 0
         endif

         if(this_n_cols /= 0) then
            locat(i_val+1:i_val+this_n_cols) = i_proc
            i_val = i_val+this_n_cols
         endif
      enddo
   endif

   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         call elsi_get_global_col(elsi_h,global_col_id,i_col)

         ! Compute destination
         dest = min(locat(global_col_id),elsi_h%n_procs-1)

         do i_row = 1,elsi_h%n_l_rows
            if(abs(ref(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1

               call elsi_get_global_row(elsi_h,global_row_id,i_row)

               ! Pack global id and data into bufs
               row_send_buf(i_val) = global_row_id
               col_send_buf(i_val) = global_col_id
               h_val_send_buf(i_val) = h_in(i_row,i_col)
               s_val_send_buf(i_val) = s_in(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   else
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         call elsi_get_global_col(elsi_h,global_col_id,i_col)

         ! Compute destination
         dest = min(locat(global_col_id),elsi_h%n_procs-1)

         do i_row = 1,elsi_h%n_l_rows
            if(abs(ref(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1

               call elsi_get_global_row(elsi_h,global_row_id,i_row)

               ! Pack global id and data into bufs
               row_send_buf(i_val) = global_row_id
               col_send_buf(i_val) = global_col_id
               h_val_send_buf(i_val) = h_in(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   endif

   nullify(ref)

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,&
           1,mpi_integer4,elsi_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   nnz_l_pexsi_aux = sum(recv_count,1)
   call MPI_Allreduce(nnz_l_pexsi_aux,elsi_h%nnz_g,1,mpi_integer4,&
           mpi_sum,elsi_h%mpi_comm,mpierr)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   do i_proc = 2,elsi_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(elsi_h,row_recv_buf,nnz_l_pexsi_aux,&
           "row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,&
           mpi_integer4,row_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buf,"row_send_buf")

   ! Column id
   call elsi_allocate(elsi_h,col_recv_buf,nnz_l_pexsi_aux,&
           "col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,&
           mpi_integer4,col_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buf,"col_send_buf")

   ! Hamiltonian value
   call elsi_allocate(elsi_h,h_val_recv_buf,nnz_l_pexsi_aux,&
           "h_val_recv_buf",caller)

   call MPI_Alltoallv(h_val_send_buf,send_count,send_displ,&
           mpi_real8,h_val_recv_buf,recv_count,recv_displ,&
           mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_send_buf,"h_val_send_buf")

   ! Overlap value
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      call elsi_allocate(elsi_h,s_val_recv_buf,nnz_l_pexsi_aux,&
              "s_val_recv_buf",caller)

      call MPI_Alltoallv(s_val_send_buf,send_count,send_displ,&
              mpi_real8,s_val_recv_buf,recv_count,recv_displ,&
              mpi_real8,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,s_val_send_buf,"s_val_send_buf")
   endif

   call elsi_allocate(elsi_h,global_id,nnz_l_pexsi_aux,"global_id",caller)

   ! Compute global 1D id
   do i_val = 1,nnz_l_pexsi_aux
      global_id(i_val) = int(col_recv_buf(i_val)-1,kind=i8)*&
                            int(elsi_h%n_basis,kind=i8)+&
                            int(row_recv_buf(i_val),kind=i8)
   enddo

   ! Sort
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      call heapsort(nnz_l_pexsi_aux,global_id,h_val_recv_buf,&
              s_val_recv_buf,row_recv_buf,col_recv_buf)
   else ! Row and column id not needed
      call heapsort(nnz_l_pexsi_aux,global_id,h_val_recv_buf)
   endif

   call elsi_deallocate(elsi_h,global_id,"global_id")

   ! Set send_count, all data sent to the first pole
   send_count = 0
   send_count(elsi_h%myid/n_para_task+1) = nnz_l_pexsi_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,&
           1,mpi_integer4,elsi_h%mpi_comm,mpierr)

   if(elsi_h%n_elsi_calls == 1) then
      ! Only the first pole knows nnz_l_sp
      elsi_h%nnz_l_sp = sum(recv_count,1)

      call MPI_Bcast(elsi_h%nnz_l_sp,1,mpi_integer4,0,&
              elsi_h%comm_among_pole,mpierr)
   endif

   ! Set send and receive displacement
   do i_proc = 2,elsi_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   if(elsi_h%n_elsi_calls == 1) then
      ! Row id
      call elsi_allocate(elsi_h,row_send_buf,elsi_h%nnz_l_sp,&
              "row_send_buf",caller)

      call MPI_Alltoallv(row_recv_buf,send_count,send_displ,&
              mpi_integer4,row_send_buf,recv_count,recv_displ,&
              mpi_integer4,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,row_recv_buf,"row_recv_buf")

      ! Column id
      call elsi_allocate(elsi_h,col_send_buf,elsi_h%nnz_l_sp,&
              "col_send_buf",caller)

      call MPI_Alltoallv(col_recv_buf,send_count,send_displ,&
              mpi_integer4,col_send_buf,recv_count,recv_displ,&
              mpi_integer4,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,col_recv_buf,"col_recv_buf")

      if(.not. elsi_h%ovlp_is_unit) then
         ! Overlap value
         call elsi_allocate(elsi_h,elsi_h%ovlp_real_pexsi,elsi_h%nnz_l_sp,&
                 "ovlp_real_pexsi",caller)

         call MPI_Alltoallv(s_val_recv_buf,send_count,send_displ,&
                 mpi_real8,elsi_h%ovlp_real_pexsi,recv_count,recv_displ,&
                 mpi_real8,elsi_h%mpi_comm,mpierr)

         call elsi_deallocate(elsi_h,s_val_recv_buf,"s_val_recv_buf")

         call elsi_allocate(elsi_h,elsi_h%ham_real_pexsi,elsi_h%nnz_l_sp,&
                 "ham_real_pexsi",caller)
      else
         call elsi_allocate(elsi_h,elsi_h%ovlp_real_pexsi,1,"dummy",caller)
      endif
   endif

   elsi_h%ham_real_pexsi = 0.0_r8

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv_buf,send_count,send_displ,&
           mpi_real8,elsi_h%ham_real_pexsi,recv_count,recv_displ,&
           mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_recv_buf,"h_val_recv_buf")
   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,elsi_h%row_ind_pexsi,elsi_h%nnz_l_sp,&
              "row_ind_pexsi",caller)

      call elsi_allocate(elsi_h,elsi_h%col_ptr_pexsi,(elsi_h%n_l_cols_sp+1),&
              "col_ptr_pexsi",caller)

      ! Only the first pole computes row index and column pointer
      if(elsi_h%my_p_row_pexsi == 0) then
         ! Compute row index and column pointer
         i_col = col_send_buf(1)-1
         do i_val = 1,elsi_h%nnz_l_sp
            elsi_h%row_ind_pexsi(i_val) = row_send_buf(i_val)

            if(col_send_buf(i_val) > i_col) then
               i_col = i_col+1
               elsi_h%col_ptr_pexsi(i_col-col_send_buf(1)+1) = i_val
            endif
         enddo

         elsi_h%col_ptr_pexsi(elsi_h%n_l_cols_sp+1) = elsi_h%nnz_l_sp+1
      endif

      call elsi_deallocate(elsi_h,row_send_buf,"row_send_buf")
      call elsi_deallocate(elsi_h,col_send_buf,"col_send_buf")
   endif

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by PEXSI.
!!
subroutine elsi_blacs_to_pexsi_hs_complex(elsi_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   complex(kind=r8),  target        :: h_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian
   complex(kind=r8),  target        :: s_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap

   complex(kind=r8), pointer :: ref(:,:) ! Sparsity pattern

   integer(kind=i4) :: n_para_task
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: global_col_id
   integer(kind=i4) :: global_row_id
   integer(kind=i4) :: d1,d2,d11,d12,d21,d22 ! Columns in the intermediate stage
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: this_n_cols
   integer(kind=i4) :: min_id
   integer(kind=i4) :: nnz_l_pexsi_aux
   integer(kind=i4) :: tmp_int
   integer(kind=i8) :: tmp_long
   complex(kind=r8) :: tmp_complex
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: h_val_send_buf(:)
   complex(kind=r8), allocatable :: s_val_send_buf(:)
   integer(kind=i4), allocatable :: row_send_buf(:)
   integer(kind=i4), allocatable :: col_send_buf(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   complex(kind=r8), allocatable :: h_val_recv_buf(:)
   complex(kind=r8), allocatable :: s_val_recv_buf(:)
   integer(kind=i4), allocatable :: row_recv_buf(:)
   integer(kind=i4), allocatable :: col_recv_buf(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: locat(:) ! Location of each global column
   integer(kind=i8), allocatable :: global_id(:) ! Global 1D id

   character*40, parameter :: caller = "elsi_blacs_to_pexsi_hs_complex"

   call elsi_get_time(elsi_h,t0)

   n_para_task = elsi_h%n_procs/elsi_h%n_p_per_pole

   if(.not. elsi_h%ovlp_is_unit) then
      ref => s_in
   else
      ref => h_in
   endif

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz(elsi_h,s_in,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,elsi_h%nnz_l)

      if(.not. elsi_h%ovlp_is_unit) then
         call elsi_allocate(elsi_h,s_val_send_buf,elsi_h%nnz_l,&
                 "s_val_send_buf",caller)
      endif
   endif

   call elsi_allocate(elsi_h,locat,elsi_h%n_basis,"locat",caller)
   call elsi_allocate(elsi_h,row_send_buf,elsi_h%nnz_l,&
           "row_send_buf",caller)
   call elsi_allocate(elsi_h,col_send_buf,elsi_h%nnz_l,&
           "col_send_buf",caller)
   call elsi_allocate(elsi_h,h_val_send_buf,elsi_h%nnz_l,&
           "h_val_send_buf",caller)

   ! Compute d1,d2,d11,d12,d21,d22 (need explanation)
   d1  = elsi_h%n_basis/elsi_h%n_p_per_pole
   d2  = elsi_h%n_basis-(elsi_h%n_p_per_pole-1)*d1
   d11 = d1/n_para_task
   d12 = d1-(n_para_task-1)*d11
   d21 = d2/n_para_task
   d22 = d2-(n_para_task-1)*d21

   i_val = 0

   if(d11 > 0) then
      do i_proc = 0,elsi_h%n_procs-n_para_task-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d12
         else
            this_n_cols = d11
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = 0,elsi_h%n_procs-n_para_task-1
         if(1+mod(i_proc,n_para_task) <= d1) then
            this_n_cols = 1
         else
            this_n_cols = 0
         endif

         if(this_n_cols /= 0) then
            locat(i_val+1:i_val+this_n_cols) = i_proc
            i_val = i_val+this_n_cols
         endif
      enddo
   endif

   if(d21 > 0) then
      do i_proc = elsi_h%n_procs-n_para_task,elsi_h%n_procs-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d22
         else
            this_n_cols = d21
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = elsi_h%n_procs-n_para_task,elsi_h%n_procs-1
         if(1+mod(i_proc,n_para_task) <= d2) then
            this_n_cols = 1
         else
            this_n_cols = 0
         endif

         if(this_n_cols /= 0) then
            locat(i_val+1:i_val+this_n_cols) = i_proc
            i_val = i_val+this_n_cols
         endif
      enddo
   endif

   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         call elsi_get_global_col(elsi_h,global_col_id,i_col)

         ! Compute destination
         dest = min(locat(global_col_id),elsi_h%n_procs-1)

         do i_row = 1,elsi_h%n_l_rows
            if(abs(ref(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1

               call elsi_get_global_row(elsi_h,global_row_id,i_row)

               ! Pack global id and data into bufs
               row_send_buf(i_val) = global_row_id
               col_send_buf(i_val) = global_col_id
               h_val_send_buf(i_val) = h_in(i_row,i_col)
               s_val_send_buf(i_val) = s_in(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   else
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         call elsi_get_global_col(elsi_h,global_col_id,i_col)

         ! Compute destination
         dest = min(locat(global_col_id),elsi_h%n_procs-1)

         do i_row = 1,elsi_h%n_l_rows
            if(abs(ref(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1

               call elsi_get_global_row(elsi_h,global_row_id,i_row)

               ! Pack global id and data into bufs
               row_send_buf(i_val) = global_row_id
               col_send_buf(i_val) = global_col_id
               h_val_send_buf(i_val) = h_in(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   endif

   nullify(ref)

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,&
           1,mpi_integer4,elsi_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   nnz_l_pexsi_aux = sum(recv_count,1)
   call MPI_Allreduce(nnz_l_pexsi_aux,elsi_h%nnz_g,1,mpi_integer4,&
           mpi_sum,elsi_h%mpi_comm,mpierr)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   do i_proc = 2,elsi_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(elsi_h,row_recv_buf,nnz_l_pexsi_aux,&
           "row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,&
           mpi_integer4,row_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buf,"row_send_buf")

   ! Column id
   call elsi_allocate(elsi_h,col_recv_buf,nnz_l_pexsi_aux,&
           "col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,&
           mpi_integer4,col_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buf,"col_send_buf")

   ! Hamiltonian value
   call elsi_allocate(elsi_h,h_val_recv_buf,nnz_l_pexsi_aux,&
           "h_val_recv_buf",caller)

   call MPI_Alltoallv(h_val_send_buf,send_count,send_displ,&
           mpi_complex16,h_val_recv_buf,recv_count,recv_displ,&
           mpi_complex16,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_send_buf,"h_val_send_buf")

   ! Overlap value
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      call elsi_allocate(elsi_h,s_val_recv_buf,nnz_l_pexsi_aux,&
              "s_val_recv_buf",caller)

      call MPI_Alltoallv(s_val_send_buf,send_count,send_displ,&
              mpi_complex16,s_val_recv_buf,recv_count,recv_displ,&
              mpi_complex16,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,s_val_send_buf,"s_val_send_buf")
   endif

   call elsi_allocate(elsi_h,global_id,nnz_l_pexsi_aux,"global_id",caller)

   ! Compute global 1D id
   do i_val = 1,nnz_l_pexsi_aux
      global_id(i_val) = int(col_recv_buf(i_val)-1,kind=i8)*&
                            int(elsi_h%n_basis,kind=i8)+&
                            int(row_recv_buf(i_val),kind=i8)
   enddo

   ! Sort
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      call heapsort(nnz_l_pexsi_aux,global_id,h_val_recv_buf,&
              s_val_recv_buf,row_recv_buf,col_recv_buf)
   else ! Row and column id not needed
      call heapsort(nnz_l_pexsi_aux,global_id,h_val_recv_buf)
   endif

   call elsi_deallocate(elsi_h,global_id,"global_id")

   ! Set send_count, all data sent to the first pole
   send_count = 0
   send_count(elsi_h%myid/n_para_task+1) = nnz_l_pexsi_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,&
           1,mpi_integer4,elsi_h%mpi_comm,mpierr)

   if(elsi_h%n_elsi_calls == 1) then
      ! Only the first pole knows nnz_l_sp
      elsi_h%nnz_l_sp = sum(recv_count,1)

      call MPI_Bcast(elsi_h%nnz_l_sp,1,mpi_integer4,0,&
              elsi_h%comm_among_pole,mpierr)
   endif

   ! Set send and receive displacement
   do i_proc = 2,elsi_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   if(elsi_h%n_elsi_calls == 1) then
      ! Row id
      call elsi_allocate(elsi_h,row_send_buf,elsi_h%nnz_l_sp,&
              "row_send_buf",caller)

      call MPI_Alltoallv(row_recv_buf,send_count,send_displ,&
              mpi_integer4,row_send_buf,recv_count,recv_displ,&
              mpi_integer4,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,row_recv_buf,"row_recv_buf")

      ! Column id
      call elsi_allocate(elsi_h,col_send_buf,elsi_h%nnz_l_sp,&
              "col_send_buf",caller)

      call MPI_Alltoallv(col_recv_buf,send_count,send_displ,&
              mpi_integer4,col_send_buf,recv_count,recv_displ,&
              mpi_integer4,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,col_recv_buf,"col_recv_buf")

      if(.not. elsi_h%ovlp_is_unit) then
         ! Overlap value
         call elsi_allocate(elsi_h,elsi_h%ovlp_complex_pexsi,elsi_h%nnz_l_sp,&
                 "ovlp_complex_pexsi",caller)

         call MPI_Alltoallv(s_val_recv_buf,send_count,send_displ,&
                 mpi_complex16,elsi_h%ovlp_complex_pexsi,recv_count,&
                 recv_displ,mpi_complex16,elsi_h%mpi_comm,mpierr)

         call elsi_deallocate(elsi_h,s_val_recv_buf,"s_val_recv_buf")

         call elsi_allocate(elsi_h,elsi_h%ham_complex_pexsi,elsi_h%nnz_l_sp,&
                 "ham_complex_pexsi",caller)
      else
         call elsi_allocate(elsi_h,elsi_h%ovlp_complex_pexsi,1,"dummy",caller)
      endif
   endif

   elsi_h%ham_complex_pexsi = (0.0_r8,0.0_r8)

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv_buf,send_count,send_displ,&
           mpi_complex16,elsi_h%ham_complex_pexsi,recv_count,recv_displ,&
           mpi_complex16,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_recv_buf,"h_val_recv_buf")
   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,elsi_h%row_ind_pexsi,elsi_h%nnz_l_sp,&
              "row_ind_pexsi",caller)

      call elsi_allocate(elsi_h,elsi_h%col_ptr_pexsi,(elsi_h%n_l_cols_sp+1),&
              "col_ptr_pexsi",caller)

      ! Only the first pole computes row index and column pointer
      if(elsi_h%my_p_row_pexsi == 0) then
         ! Compute row index and column pointer
         i_col = col_send_buf(1)-1
         do i_val = 1,elsi_h%nnz_l_sp
            elsi_h%row_ind_pexsi(i_val) = row_send_buf(i_val)

            if(col_send_buf(i_val) > i_col) then
               i_col = i_col+1
               elsi_h%col_ptr_pexsi(i_col-col_send_buf(1)+1) = i_val
            endif
         enddo

         elsi_h%col_ptr_pexsi(elsi_h%n_l_cols_sp+1) = elsi_h%nnz_l_sp+1
      endif

      call elsi_deallocate(elsi_h,row_send_buf,"row_send_buf")
      call elsi_deallocate(elsi_h,col_send_buf,"col_send_buf")
   endif

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine converts density matrix computed by PEXSI and stored
!! in 1D block distributed sparse CCS format to 2D block-cyclic
!! distributed dense format.
!!
subroutine elsi_pexsi_to_blacs_dm_real(elsi_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                 !< Handle
   real(kind=r8),     intent(out)   :: d_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: local_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: local_row_id ! Local row id in 1D block distribution
   integer(kind=i4) :: proc_col_id ! Column id in process grid
   integer(kind=i4) :: proc_row_id ! Row id in process grid
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: val_send_buf(:)
   integer(kind=i4), allocatable :: row_send_buf(:)
   integer(kind=i4), allocatable :: col_send_buf(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: val_recv_buf(:)
   integer(kind=i4), allocatable :: row_recv_buf(:)
   integer(kind=i4), allocatable :: col_recv_buf(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: global_col_id(:)
   integer(kind=i4), allocatable :: global_row_id(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_dm_real"

   call elsi_get_time(elsi_h,t0)

   call elsi_allocate(elsi_h,val_send_buf,elsi_h%nnz_l_sp,&
           "val_send_buf",caller)
   call elsi_allocate(elsi_h,row_send_buf,elsi_h%nnz_l_sp,&
           "row_send_buf",caller)
   call elsi_allocate(elsi_h,col_send_buf,elsi_h%nnz_l_sp,&
           "col_send_buf",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   if(elsi_h%my_p_row_pexsi == 0) then
      call elsi_allocate(elsi_h,global_row_id,elsi_h%nnz_l_sp,&
              "global_row_id",caller)
      call elsi_allocate(elsi_h,global_col_id,elsi_h%nnz_l_sp,&
              "global_col_id",caller)
      call elsi_allocate(elsi_h,dest,elsi_h%nnz_l_sp,"dest",caller)

      i_col = 0
      ! Compute destination and global id
      do i_val = 1,elsi_h%nnz_l_sp
         if(i_val == elsi_h%col_ptr_ccs(i_col+1) .and. &
            i_col /= elsi_h%n_l_cols_sp) then
            i_col = i_col+1
         endif
         i_row = elsi_h%row_ind_ccs(i_val)

         ! Compute global id
         global_row_id(i_val) = i_row
         global_col_id(i_val) = i_col+elsi_h%myid*(elsi_h%n_basis/elsi_h%n_p_per_pole)

         ! Compute destination
         proc_row_id = mod((global_row_id(i_val)-1)/elsi_h%n_b_rows,elsi_h%n_p_rows)
         proc_col_id = mod((global_col_id(i_val)-1)/elsi_h%n_b_cols,elsi_h%n_p_cols)
         dest(i_val) = proc_col_id+proc_row_id*elsi_h%n_p_cols
      enddo

      j_val = 0
      ! Set send_count
      do i_proc = 1,elsi_h%n_procs
         do i_val = 1,elsi_h%nnz_l_sp
            if(dest(i_val) == i_proc-1) then
               j_val = j_val+1
               val_send_buf(j_val) = elsi_h%dm_real_ccs(i_val)
               row_send_buf(j_val) = global_row_id(i_val)
               col_send_buf(j_val) = global_col_id(i_val)
               send_count(i_proc) = send_count(i_proc)+1
            endif
         enddo
      enddo

      call elsi_deallocate(elsi_h,global_row_id,"global_row_id")
      call elsi_deallocate(elsi_h,global_col_id,"global_col_id")
      call elsi_deallocate(elsi_h,dest,"dest")
   endif

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,&
           1,mpi_integer4,elsi_h%mpi_comm,mpierr)

   elsi_h%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   do i_proc = 2,elsi_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Value
   call elsi_allocate(elsi_h,val_recv_buf,elsi_h%nnz_l,&
           "val_recv_buf",caller)

   call MPI_Alltoallv(val_send_buf,send_count,send_displ,&
           mpi_real8,val_recv_buf,recv_count,recv_displ,&
           mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,val_send_buf,"val_send_buf")

   ! Row index
   call elsi_allocate(elsi_h,row_recv_buf,elsi_h%nnz_l,&
           "row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,&
           mpi_integer4,row_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buf,"row_send_buf")

   ! Column index
   call elsi_allocate(elsi_h,col_recv_buf,elsi_h%nnz_l,&
           "col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,&
           mpi_integer4,col_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buf,"col_send_buf")
   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   d_out = 0.0_r8

   ! Unpack density matrix
   do i_val = 1,elsi_h%nnz_l
      ! Compute local 2d id
      local_row_id = (row_recv_buf(i_val)-1)/(elsi_h%n_p_rows*elsi_h%n_b_rows)*&
                        elsi_h%n_b_rows+mod((row_recv_buf(i_val)-1),elsi_h%n_b_rows)+1
      local_col_id = (col_recv_buf(i_val)-1)/(elsi_h%n_p_cols*elsi_h%n_b_cols)*&
                        elsi_h%n_b_cols+mod((col_recv_buf(i_val)-1),elsi_h%n_b_cols)+1

      ! Put value to correct position
      d_out(local_row_id,local_col_id) = val_recv_buf(i_val)
   enddo

   call elsi_deallocate(elsi_h,val_recv_buf,"val_recv_buf")
   call elsi_deallocate(elsi_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(elsi_h,col_recv_buf,"col_recv_buf")

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine converts density matrix computed by PEXSI and stored
!! in 1D block distributed sparse CCS format to 2D block-cyclic
!! distributed dense format.
!!
subroutine elsi_pexsi_to_blacs_dm_complex(elsi_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                 !< Handle
   complex(kind=r8),  intent(out)   :: d_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: local_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: local_row_id ! Local row id in 1D block distribution
   integer(kind=i4) :: proc_col_id ! Column id in process grid
   integer(kind=i4) :: proc_row_id ! Row id in process grid
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: val_send_buf(:)
   integer(kind=i4), allocatable :: row_send_buf(:)
   integer(kind=i4), allocatable :: col_send_buf(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   complex(kind=r8), allocatable :: val_recv_buf(:)
   integer(kind=i4), allocatable :: row_recv_buf(:)
   integer(kind=i4), allocatable :: col_recv_buf(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: global_col_id(:)
   integer(kind=i4), allocatable :: global_row_id(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_dm_complex"

   call elsi_get_time(elsi_h,t0)

   call elsi_allocate(elsi_h,val_send_buf,elsi_h%nnz_l_sp,&
           "val_send_buf",caller)
   call elsi_allocate(elsi_h,row_send_buf,elsi_h%nnz_l_sp,&
           "row_send_buf",caller)
   call elsi_allocate(elsi_h,col_send_buf,elsi_h%nnz_l_sp,&
           "col_send_buf",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   if(elsi_h%my_p_row_pexsi == 0) then
      call elsi_allocate(elsi_h,global_row_id,elsi_h%nnz_l_sp,&
              "global_row_id",caller)
      call elsi_allocate(elsi_h,global_col_id,elsi_h%nnz_l_sp,&
              "global_col_id",caller)
      call elsi_allocate(elsi_h,dest,elsi_h%nnz_l_sp,"dest",caller)

      i_col = 0
      ! Compute destination and global id
      do i_val = 1,elsi_h%nnz_l_sp
         if(i_val == elsi_h%col_ptr_ccs(i_col+1) .and. &
            i_col /= elsi_h%n_l_cols_sp) then
            i_col = i_col+1
         endif
         i_row = elsi_h%row_ind_ccs(i_val)

         ! Compute global id
         global_row_id(i_val) = i_row
         global_col_id(i_val) = i_col+elsi_h%myid*(elsi_h%n_basis/elsi_h%n_p_per_pole)

         ! Compute destination
         proc_row_id = mod((global_row_id(i_val)-1)/elsi_h%n_b_rows,elsi_h%n_p_rows)
         proc_col_id = mod((global_col_id(i_val)-1)/elsi_h%n_b_cols,elsi_h%n_p_cols)
         dest(i_val) = proc_col_id+proc_row_id*elsi_h%n_p_cols
      enddo

      j_val = 0
      ! Set send_count
      do i_proc = 1,elsi_h%n_procs
         do i_val = 1,elsi_h%nnz_l_sp
            if(dest(i_val) == i_proc-1) then
               j_val = j_val+1
               val_send_buf(j_val) = elsi_h%dm_complex_ccs(i_val)
               row_send_buf(j_val) = global_row_id(i_val)
               col_send_buf(j_val) = global_col_id(i_val)
               send_count(i_proc) = send_count(i_proc)+1
            endif
         enddo
      enddo

      call elsi_deallocate(elsi_h,global_row_id,"global_row_id")
      call elsi_deallocate(elsi_h,global_col_id,"global_col_id")
      call elsi_deallocate(elsi_h,dest,"dest")
   endif

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,&
           1,mpi_integer4,elsi_h%mpi_comm,mpierr)

   elsi_h%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   do i_proc = 2,elsi_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Value
   call elsi_allocate(elsi_h,val_recv_buf,elsi_h%nnz_l,&
           "val_recv_buf",caller)

   call MPI_Alltoallv(val_send_buf,send_count,send_displ,&
           mpi_complex16,val_recv_buf,recv_count,recv_displ,&
           mpi_complex16,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,val_send_buf,"val_send_buf")

   ! Row index
   call elsi_allocate(elsi_h,row_recv_buf,elsi_h%nnz_l,&
           "row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,&
           mpi_integer4,row_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buf,"row_send_buf")

   ! Column index
   call elsi_allocate(elsi_h,col_recv_buf,elsi_h%nnz_l,&
           "col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,&
           mpi_integer4,col_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buf,"col_send_buf")
   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   d_out = 0.0_r8

   ! Unpack density matrix
   do i_val = 1,elsi_h%nnz_l
      ! Compute local 2d id
      local_row_id = (row_recv_buf(i_val)-1)/(elsi_h%n_p_rows*elsi_h%n_b_rows)*&
                        elsi_h%n_b_rows+mod((row_recv_buf(i_val)-1),elsi_h%n_b_rows)+1
      local_col_id = (col_recv_buf(i_val)-1)/(elsi_h%n_p_cols*elsi_h%n_b_cols)*&
                        elsi_h%n_b_cols+mod((col_recv_buf(i_val)-1),elsi_h%n_b_cols)+1

      ! Put value to correct position
      d_out(local_row_id,local_col_id) = val_recv_buf(i_val)
   enddo

   call elsi_deallocate(elsi_h,val_recv_buf,"val_recv_buf")
   call elsi_deallocate(elsi_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(elsi_h,col_recv_buf,"col_recv_buf")

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by SIPs.
!!
subroutine elsi_blacs_to_sips_hs_real(elsi_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     target        :: h_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian
   real(kind=r8),     target        :: s_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap

   real(kind=r8), pointer :: ref(:,:) ! Sparsity pattern

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: global_col_id
   integer(kind=i4) :: global_row_id
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: tmp_int
   integer(kind=i8) :: tmp_long
   integer(kind=i4) :: min_id
   real(kind=r8)    :: tmp_real
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send_buf(:)
   real(kind=r8),    allocatable :: s_val_send_buf(:)
   integer(kind=i4), allocatable :: row_send_buf(:)
   integer(kind=i4), allocatable :: col_send_buf(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv_buf(:)
   integer(kind=i4), allocatable :: col_recv_buf(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: global_id(:)

   character*40, parameter :: caller = "elsi_blacs_to_sips_hs_real"

   call elsi_get_time(elsi_h,t0)

   if(.not. elsi_h%ovlp_is_unit) then
      ref => s_in
   else
      ref => h_in
   endif

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz(elsi_h,h_in,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,elsi_h%nnz_l)

      if(.not. elsi_h%ovlp_is_unit) then
         call elsi_allocate(elsi_h,s_val_send_buf,elsi_h%nnz_l,&
                 "s_val_send_buf",caller)
      endif
   endif

   call elsi_allocate(elsi_h,row_send_buf,elsi_h%nnz_l,&
           "row_send_buf",caller)
   call elsi_allocate(elsi_h,col_send_buf,elsi_h%nnz_l,&
           "col_send_buf",caller)
   call elsi_allocate(elsi_h,h_val_send_buf,elsi_h%nnz_l,&
           "h_val_send_buf",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         call elsi_get_global_col(elsi_h,global_col_id,i_col)

         ! Compute destination
         dest = (global_col_id-1)/(elsi_h%n_basis/elsi_h%n_procs)
         ! The last process may take more
         dest = min(dest,elsi_h%n_procs-1)

         do i_row = 1,elsi_h%n_l_rows
            if(abs(ref(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1

               call elsi_get_global_row(elsi_h,global_row_id,i_row)

               ! Pack global id and data into bufs
               row_send_buf(i_val) = global_row_id
               col_send_buf(i_val) = global_col_id
               h_val_send_buf(i_val) = h_in(i_row,i_col)
               s_val_send_buf(i_val) = s_in(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
           endif
        enddo
      enddo
   else
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         call elsi_get_global_col(elsi_h,global_col_id,i_col)

         ! Compute destination
         dest = (global_col_id-1)/(elsi_h%n_basis/elsi_h%n_procs)
         ! The last process may take more
         dest = min(dest,elsi_h%n_procs-1)

         do i_row = 1,elsi_h%n_l_rows
            if(abs(ref(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1

               call elsi_get_global_row(elsi_h,global_row_id,i_row)

               ! Pack global id and data into bufs
               row_send_buf(i_val) = global_row_id
               col_send_buf(i_val) = global_col_id
               h_val_send_buf(i_val) = h_in(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
           endif
        enddo
      enddo
   endif

   nullify(ref)

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,&
           1,mpi_integer4,elsi_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   elsi_h%nnz_l_sp = sum(recv_count,1)
   call MPI_Allreduce(elsi_h%nnz_l_sp,elsi_h%nnz_g,1,mpi_integer4,&
            mpi_sum,elsi_h%mpi_comm,mpierr)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   do i_proc = 2,elsi_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(elsi_h,row_recv_buf,elsi_h%nnz_l_sp,&
           "row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,&
           mpi_integer4,row_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buf,"row_send_buf")

   ! Column id
   call elsi_allocate(elsi_h,col_recv_buf,elsi_h%nnz_l_sp,&
           "col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,&
           mpi_integer4,col_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buf,"col_send_buf")

   ! Hamiltonian value
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,elsi_h%ham_real_sips,elsi_h%nnz_l_sp,&
              "ham_real_sips",caller)
   endif

   call MPI_Alltoallv(h_val_send_buf,send_count,send_displ,&
           mpi_real8,elsi_h%ham_real_sips,recv_count,recv_displ,&
           mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_send_buf,"h_val_send_buf")

   ! Overlap value
   if(elsi_h%n_elsi_calls == 1) then
      if(.not. elsi_h%ovlp_is_unit) then
         call elsi_allocate(elsi_h,elsi_h%ovlp_real_sips,elsi_h%nnz_l_sp,&
                 "ovlp_real_sips",caller)

         call MPI_Alltoallv(s_val_send_buf,send_count,send_displ,&
                 mpi_real8,elsi_h%ovlp_real_sips,recv_count,recv_displ,&
                 mpi_real8,elsi_h%mpi_comm,mpierr)

         call elsi_deallocate(elsi_h,s_val_send_buf,"s_val_send_buf")
      else
         call elsi_allocate(elsi_h,elsi_h%ovlp_real_sips,1,"dummy",caller)
      endif
   endif

   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   call elsi_allocate(elsi_h,global_id,elsi_h%nnz_l_sp,"global_id",caller)

   ! Compute global 1D id
   do i_val = 1,elsi_h%nnz_l_sp
      global_id(i_val) = int(col_recv_buf(i_val)-1,kind=i8)*&
                            int(elsi_h%n_basis,kind=i8)+&
                            int(row_recv_buf(i_val),kind=i8)
   enddo

   ! Sort
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      call heapsort(elsi_h%nnz_l_sp,global_id,elsi_h%ham_real_sips,&
              elsi_h%ovlp_real_sips,row_recv_buf,col_recv_buf)
   else
      call heapsort(elsi_h%nnz_l_sp,global_id,elsi_h%ham_real_sips)
   endif

   call elsi_deallocate(elsi_h,global_id,"global_id")

   ! Compute row index and column pointer
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,elsi_h%row_ind_sips,elsi_h%nnz_l_sp,&
              "row_ind_sips",caller)

      call elsi_allocate(elsi_h,elsi_h%col_ptr_sips,(elsi_h%n_l_cols_sp+1),&
              "col_ptr_sips",caller)

      i_col = col_recv_buf(1)-1
      do i_val = 1,elsi_h%nnz_l_sp
         elsi_h%row_ind_sips(i_val) = row_recv_buf(i_val)

         if(col_recv_buf(i_val) > i_col) then
            i_col = i_col+1
            elsi_h%col_ptr_sips(i_col-col_recv_buf(1)+1) = i_val
         endif
      enddo

      elsi_h%col_ptr_sips(elsi_h%n_l_cols_sp+1) = elsi_h%nnz_l_sp+1
   endif

   call elsi_deallocate(elsi_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(elsi_h,col_recv_buf,"col_recv_buf")

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by SIPs.
!!
subroutine elsi_blacs_to_sips_hs_complex(elsi_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   complex(kind=r8),  target        :: h_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian
   complex(kind=r8),  target        :: s_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap

   complex(kind=r8), pointer :: ref(:,:) ! Sparsity pattern

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: global_col_id
   integer(kind=i4) :: global_row_id
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: tmp_int
   integer(kind=i8) :: tmp_long
   integer(kind=i4) :: min_id
   complex(kind=r8) :: tmp_complex
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: h_val_send_buf(:)
   complex(kind=r8), allocatable :: s_val_send_buf(:)
   integer(kind=i4), allocatable :: row_send_buf(:)
   integer(kind=i4), allocatable :: col_send_buf(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv_buf(:)
   integer(kind=i4), allocatable :: col_recv_buf(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: global_id(:)

   character*40, parameter :: caller = "elsi_blacs_to_sips_hs_complex"

   call elsi_get_time(elsi_h,t0)

   if(.not. elsi_h%ovlp_is_unit) then
      ref => s_in
   else
      ref => h_in
   endif

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz(elsi_h,h_in,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,elsi_h%nnz_l)

      if(.not. elsi_h%ovlp_is_unit) then
         call elsi_allocate(elsi_h,s_val_send_buf,elsi_h%nnz_l,&
                 "s_val_send_buf",caller)
      endif
   endif

   call elsi_allocate(elsi_h,row_send_buf,elsi_h%nnz_l,&
           "row_send_buf",caller)
   call elsi_allocate(elsi_h,col_send_buf,elsi_h%nnz_l,&
           "col_send_buf",caller)
   call elsi_allocate(elsi_h,h_val_send_buf,elsi_h%nnz_l,&
           "h_val_send_buf",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         call elsi_get_global_col(elsi_h,global_col_id,i_col)

         ! Compute destination
         dest = (global_col_id-1)/(elsi_h%n_basis/elsi_h%n_procs)
         ! The last process may take more
         dest = min(dest,elsi_h%n_procs-1)

         do i_row = 1,elsi_h%n_l_rows
            if(abs(ref(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1

               call elsi_get_global_row(elsi_h,global_row_id,i_row)

               ! Pack global id and data into bufs
               row_send_buf(i_val) = global_row_id
               col_send_buf(i_val) = global_col_id
               h_val_send_buf(i_val) = h_in(i_row,i_col)
               s_val_send_buf(i_val) = s_in(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
           endif
        enddo
      enddo
   else
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         call elsi_get_global_col(elsi_h,global_col_id,i_col)

         ! Compute destination
         dest = (global_col_id-1)/(elsi_h%n_basis/elsi_h%n_procs)
         ! The last process may take more
         dest = min(dest,elsi_h%n_procs-1)

         do i_row = 1,elsi_h%n_l_rows
            if(abs(ref(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1

               call elsi_get_global_row(elsi_h,global_row_id,i_row)

               ! Pack global id and data into bufs
               row_send_buf(i_val) = global_row_id
               col_send_buf(i_val) = global_col_id
               h_val_send_buf(i_val) = h_in(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
           endif
        enddo
      enddo
   endif

   nullify(ref)

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,&
           1,mpi_integer4,elsi_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   elsi_h%nnz_l_sp = sum(recv_count,1)
   call MPI_Allreduce(elsi_h%nnz_l_sp,elsi_h%nnz_g,1,mpi_integer4,&
            mpi_sum,elsi_h%mpi_comm,mpierr)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   do i_proc = 2,elsi_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(elsi_h,row_recv_buf,elsi_h%nnz_l_sp,&
           "row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,&
           mpi_integer4,row_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buf,"row_send_buf")

   ! Column id
   call elsi_allocate(elsi_h,col_recv_buf,elsi_h%nnz_l_sp,&
           "col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,&
           mpi_integer4,col_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buf,"col_send_buf")

   ! Hamiltonian value
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,elsi_h%ham_complex_sips,elsi_h%nnz_l_sp,&
              "ham_complex_sips",caller)
   endif

   call MPI_Alltoallv(h_val_send_buf,send_count,send_displ,&
           mpi_complex16,elsi_h%ham_complex_sips,recv_count,recv_displ,&
           mpi_complex16,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_send_buf,"h_val_send_buf")

   ! Overlap value
   if(elsi_h%n_elsi_calls == 1) then
      if(.not. elsi_h%ovlp_is_unit) then
         call elsi_allocate(elsi_h,elsi_h%ovlp_complex_sips,elsi_h%nnz_l_sp,&
                 "ovlp_complex_sips",caller)

         call MPI_Alltoallv(s_val_send_buf,send_count,send_displ,&
                 mpi_complex16,elsi_h%ovlp_complex_sips,recv_count,recv_displ,&
                 mpi_complex16,elsi_h%mpi_comm,mpierr)

         call elsi_deallocate(elsi_h,s_val_send_buf,"s_val_send_buf")
      else
         call elsi_allocate(elsi_h,elsi_h%ovlp_complex_sips,1,"dummy",caller)
      endif
   endif

   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   call elsi_allocate(elsi_h,global_id,elsi_h%nnz_l_sp,"global_id",caller)

   ! Compute global 1D id
   do i_val = 1,elsi_h%nnz_l_sp
      global_id(i_val) = int(col_recv_buf(i_val)-1,kind=i8)*&
                            int(elsi_h%n_basis,kind=i8)+&
                            int(row_recv_buf(i_val),kind=i8)
   enddo

   ! Sort
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      call heapsort(elsi_h%nnz_l_sp,global_id,elsi_h%ham_complex_sips,&
              elsi_h%ovlp_complex_sips,row_recv_buf,col_recv_buf)
   else
      call heapsort(elsi_h%nnz_l_sp,global_id,elsi_h%ham_complex_sips)
   endif

   call elsi_deallocate(elsi_h,global_id,"global_id")

   ! Compute row index and column pointer
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,elsi_h%row_ind_sips,elsi_h%nnz_l_sp,&
              "row_ind_sips",caller)

      call elsi_allocate(elsi_h,elsi_h%col_ptr_sips,(elsi_h%n_l_cols_sp+1),&
              "col_ptr_sips",caller)

      i_col = col_recv_buf(1)-1
      do i_val = 1,elsi_h%nnz_l_sp
         elsi_h%row_ind_sips(i_val) = row_recv_buf(i_val)

         if(col_recv_buf(i_val) > i_col) then
            i_col = i_col+1
            elsi_h%col_ptr_sips(i_col-col_recv_buf(1)+1) = i_val
         endif
      enddo

      elsi_h%col_ptr_sips(elsi_h%n_l_cols_sp+1) = elsi_h%nnz_l_sp+1
   endif

   call elsi_deallocate(elsi_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(elsi_h,col_recv_buf,"col_recv_buf")

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 1D block distributed sparse CCS format to 2D block-cyclic
!! distributed dense format, which can be used as input by ELPA.
!!
subroutine elsi_sips_to_blacs_hs_real(elsi_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                !< Handle
   real(kind=r8),     intent(in)    :: h_in(elsi_h%nnz_l_sp) !< Hamiltonian
   real(kind=r8),     intent(in)    :: s_in(elsi_h%nnz_l_sp) !< Overlap

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: local_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: local_row_id ! Local row id in 1D block distribution
   integer(kind=i4) :: proc_col_id ! Column id in process grid
   integer(kind=i4) :: proc_row_id ! Row id in process grid
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send_buf(:)
   real(kind=r8),    allocatable :: s_val_send_buf(:)
   integer(kind=i4), allocatable :: row_send_buf(:)
   integer(kind=i4), allocatable :: col_send_buf(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: h_val_recv_buf(:)
   real(kind=r8),    allocatable :: s_val_recv_buf(:)
   integer(kind=i4), allocatable :: row_recv_buf(:)
   integer(kind=i4), allocatable :: col_recv_buf(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: global_row_id(:)
   integer(kind=i4), allocatable :: global_col_id(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character*40, parameter :: caller = "elsi_sips_to_blacs_hs_real"

   call elsi_get_time(elsi_h,t0)

   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      call elsi_allocate(elsi_h,s_val_send_buf,elsi_h%nnz_l_sp,&
              "s_val_send_buf",caller)
   endif

   call elsi_allocate(elsi_h,h_val_send_buf,elsi_h%nnz_l_sp,&
           "h_val_send_buf",caller)
   call elsi_allocate(elsi_h,row_send_buf,elsi_h%nnz_l_sp,&
           "pos_send_buf",caller)
   call elsi_allocate(elsi_h,col_send_buf,elsi_h%nnz_l_sp,&
           "pos_send_buf",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,&
           "send_count",caller)

   call elsi_allocate(elsi_h,global_row_id,elsi_h%nnz_l_sp,&
           "global_row_id",caller)
   call elsi_allocate(elsi_h,global_col_id,elsi_h%nnz_l_sp,&
           "global_col_id",caller)
   call elsi_allocate(elsi_h,dest,elsi_h%nnz_l_sp,"dest",caller)

   i_col = 0
   ! Compute destination and global id
   do i_val = 1,elsi_h%nnz_l_sp
      if(i_val == elsi_h%col_ptr_ccs(i_col+1) .and. &
         i_col /= elsi_h%n_l_cols_sp) then
         i_col = i_col+1
      endif
      i_row = elsi_h%row_ind_ccs(i_val)

      ! Compute global id
      global_row_id(i_val) = i_row
      global_col_id(i_val) = i_col+elsi_h%myid*(elsi_h%n_basis/elsi_h%n_procs)

      ! Compute destination
      proc_row_id = mod((global_row_id(i_val)-1)/elsi_h%n_b_rows,elsi_h%n_p_rows)
      proc_col_id = mod((global_col_id(i_val)-1)/elsi_h%n_b_cols,elsi_h%n_p_cols)
      dest(i_val) = proc_col_id+proc_row_id*elsi_h%n_p_cols
   enddo

   j_val = 0
   ! Set send_count
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      do i_proc = 1,elsi_h%n_procs
         do i_val = 1,elsi_h%nnz_l_sp
            if(dest(i_val) == i_proc-1) then
               j_val = j_val+1
               h_val_send_buf(j_val) = h_in(i_val)
               s_val_send_buf(j_val) = s_in(i_val)
               row_send_buf(j_val) = global_row_id(i_val)
               col_send_buf(j_val) = global_col_id(i_val)
               send_count(i_proc) = send_count(i_proc)+1
            endif
         enddo
      enddo
   else
      do i_proc = 1,elsi_h%n_procs
         do i_val = 1,elsi_h%nnz_l_sp
            if(dest(i_val) == i_proc-1) then
               j_val = j_val+1
               h_val_send_buf(j_val) = h_in(i_val)
               row_send_buf(j_val) = global_row_id(i_val)
               col_send_buf(j_val) = global_col_id(i_val)
               send_count(i_proc) = send_count(i_proc)+1
            endif
         enddo
      enddo
   endif

   call elsi_deallocate(elsi_h,global_row_id,"global_row_id")
   call elsi_deallocate(elsi_h,global_col_id,"global_col_id")
   call elsi_deallocate(elsi_h,dest,"dest")

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,&
           1,mpi_integer4,elsi_h%mpi_comm,mpierr)

   elsi_h%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   do i_proc = 2,elsi_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row index
   call elsi_allocate(elsi_h,row_recv_buf,elsi_h%nnz_l,&
           "row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,&
           mpi_integer4,row_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buf,"row_send_buf")

   ! Column index
   call elsi_allocate(elsi_h,col_recv_buf,elsi_h%nnz_l,&
           "col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,&
           mpi_integer4,col_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buf,"col_send_buf")

   ! Hamiltonian Value
   call elsi_allocate(elsi_h,h_val_recv_buf,elsi_h%nnz_l,&
           "h_val_recv_buf",caller)

   call MPI_Alltoallv(h_val_send_buf,send_count,send_displ,&
           mpi_real8,h_val_recv_buf,recv_count,recv_displ,&
           mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_send_buf,"h_val_send_buf")

   ! Overlap value
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      call elsi_allocate(elsi_h,s_val_recv_buf,elsi_h%nnz_l,&
              "s_val_recv_buf",caller)

      call MPI_Alltoallv(s_val_send_buf,send_count,send_displ,&
              mpi_real8,s_val_recv_buf,recv_count,recv_displ,&
              mpi_real8,elsi_h%mpi_comm,mpierr)
   endif

   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   ! Unpack matrix
   if((elsi_h%n_elsi_calls == 1) .and. (.not. elsi_h%ovlp_is_unit)) then
      call elsi_allocate(elsi_h,elsi_h%ham_real_elpa,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,"ham_real_elpa",caller)

      call elsi_allocate(elsi_h,elsi_h%ovlp_real_elpa,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,"ovlp_real_elpa",caller)

      do i_val = 1,elsi_h%nnz_l
         ! Compute local 2d id
         local_row_id = (row_recv_buf(i_val)-1)/&
                           (elsi_h%n_p_rows*elsi_h%n_b_rows)*elsi_h%n_b_rows+&
                           mod((row_recv_buf(i_val)-1),elsi_h%n_b_rows)+1
         local_col_id = (col_recv_buf(i_val)-1)/&
                           (elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols+&
                           mod((col_recv_buf(i_val)-1),elsi_h%n_b_cols)+1

         ! Put value to correct position
         elsi_h%ham_real_elpa(local_row_id,local_col_id) = h_val_recv_buf(i_val)
         elsi_h%ovlp_real_elpa(local_row_id,local_col_id) = s_val_recv_buf(i_val)
      enddo

      call elsi_deallocate(elsi_h,s_val_recv_buf,"s_val_recv_buf")
   else
      if(elsi_h%n_elsi_calls == 1) then
         call elsi_allocate(elsi_h,elsi_h%ham_real_elpa,elsi_h%n_l_rows,&
                 elsi_h%n_l_cols,"ham_real_elpa",caller)

         call elsi_allocate(elsi_h,elsi_h%ovlp_real_elpa,1,1,"dummy",caller)
      endif

      elsi_h%ham_real_elpa = 0.0_r8

      do i_val = 1,elsi_h%nnz_l
         ! Compute local 2d id
         local_row_id = (row_recv_buf(i_val)-1)/&
                           (elsi_h%n_p_rows*elsi_h%n_b_rows)*elsi_h%n_b_rows+&
                           mod((row_recv_buf(i_val)-1),elsi_h%n_b_rows)+1
         local_col_id = (col_recv_buf(i_val)-1)/&
                           (elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols+&
                           mod((col_recv_buf(i_val)-1),elsi_h%n_b_cols)+1

         ! Put value to correct position
         elsi_h%ham_real_elpa(local_row_id,local_col_id) = h_val_recv_buf(i_val)
      enddo
   endif

   call elsi_deallocate(elsi_h,h_val_recv_buf,"h_val_recv_buf")
   call elsi_deallocate(elsi_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(elsi_h,col_recv_buf,"col_recv_buf")

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine converts density matrix in 2D block-cyclic distributed
!! dense format to 1D block distributed sparse CCS format.
!!
subroutine elsi_blacs_to_sips_dm_real(elsi_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                 !< Handle
   real(kind=r8),     intent(out)   :: d_out(elsi_h%nnz_l_sp) !< Density matrix

   real(kind=r8), pointer :: ref(:,:) ! Sparsity pattern

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: local_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: local_row_id ! Local row id in 1D block distribution
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: nnz_l_pexsi_aux
   integer(kind=i4) :: n_l_cols_pexsi_aux
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: val_send_buf(:)
   integer(kind=i4), allocatable :: row_send_buf(:)
   integer(kind=i4), allocatable :: col_send_buf(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: val_recv_buf(:)
   integer(kind=i4), allocatable :: row_recv_buf(:)
   integer(kind=i4), allocatable :: col_recv_buf(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)

   character*40, parameter :: caller = "elsi_blacs_to_sips_dm_real"

   call elsi_get_time(elsi_h,t0)

   if(elsi_h%solver == ELPA) then
      ref => elsi_h%dm_real
   elseif(elsi_h%solver == LIBOMM) then
      ref => elsi_h%dm_omm%dval
   endif

   call elsi_get_local_nnz(elsi_h,ref,elsi_h%n_l_rows,elsi_h%n_l_cols,&
           elsi_h%nnz_l)

   call elsi_allocate(elsi_h,val_send_buf,elsi_h%nnz_l,&
           "val_send_buf",caller)
   call elsi_allocate(elsi_h,row_send_buf,elsi_h%nnz_l,&
           "row_send_buf",caller)
   call elsi_allocate(elsi_h,col_send_buf,elsi_h%nnz_l,&
           "col_send_buf",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0
   do i_col = 1,elsi_h%n_l_cols
      do i_row = 1,elsi_h%n_l_rows
         if(abs(ref(i_row,i_col)) > elsi_h%zero_threshold) then
            i_val = i_val+1

            call elsi_get_global_row(elsi_h,row_send_buf(i_val),i_row)
            call elsi_get_global_col(elsi_h,col_send_buf(i_val),i_col)

            ! Compute destination
            dest = (col_send_buf(i_val)-1)/(elsi_h%n_basis/elsi_h%n_procs)
            ! The last process may take more
            dest = min(dest,elsi_h%n_procs-1)

            ! Pack data
            val_send_buf(i_val) = ref(i_row,i_col)

            ! Set send_count
            send_count(dest+1) = send_count(dest+1)+1
         endif
      enddo
   enddo

   nullify(ref)

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,&
           1,mpi_integer4,elsi_h%mpi_comm,mpierr)

   ! Set local number of nonzero
   nnz_l_pexsi_aux = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   do i_proc = 2,elsi_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(elsi_h,row_recv_buf,nnz_l_pexsi_aux,&
           "row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,&
           mpi_integer4,row_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buf,"row_send_buf")

   ! Column id
   call elsi_allocate(elsi_h,col_recv_buf,nnz_l_pexsi_aux,&
           "col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,&
           mpi_integer4,col_recv_buf,recv_count,recv_displ,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buf,"col_send_buf")

   ! Density matrix value
   call elsi_allocate(elsi_h,val_recv_buf,nnz_l_pexsi_aux,&
           "val_recv_buf",caller)

   call MPI_Alltoallv(val_send_buf,send_count,send_displ,&
           mpi_real8,val_recv_buf,recv_count,recv_displ,&
           mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,val_send_buf,"val_send_buf")

   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   d_out = 0.0_r8

   if(elsi_h%myid == elsi_h%n_procs-1) then
      n_l_cols_pexsi_aux = elsi_h%n_basis-elsi_h%n_l_cols_sp
      n_l_cols_pexsi_aux = n_l_cols_pexsi_aux/(elsi_h%n_procs-1)
   else
      n_l_cols_pexsi_aux = elsi_h%n_l_cols_sp
   endif

   ! Unpack matrix
   do i_val = 1,nnz_l_pexsi_aux
      local_col_id = col_recv_buf(i_val)-elsi_h%myid*n_l_cols_pexsi_aux
      local_row_id = row_recv_buf(i_val)

      do j_val = elsi_h%col_ptr_ccs(local_col_id),elsi_h%col_ptr_ccs(local_col_id+1)-1
         if(elsi_h%row_ind_ccs(j_val) == local_row_id) then
            d_out(j_val) = val_recv_buf(i_val)
         endif
      enddo
   enddo

   call elsi_deallocate(elsi_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(elsi_h,col_recv_buf,"col_recv_buf")
   call elsi_deallocate(elsi_h,val_recv_buf,"val_recv_buf")

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine converts matrix format and distribution from BLACS to CheSS.
!!
subroutine elsi_blacs_to_chess_hs_real(elsi_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     intent(in)    :: h_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian
   real(kind=r8),     intent(in)    :: s_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character*200 :: info_str

   character*40, parameter :: caller = "elsi_blacs_to_chess_hs_real"

   call elsi_get_time(elsi_h,t0)

   ! First convert to SIPs 1D block distribution
   elsi_h%n_l_cols_sp = elsi_h%n_basis/elsi_h%n_procs

   ! The last process holds all remaining columns
   if(elsi_h%myid == elsi_h%n_procs-1) then
      elsi_h%n_l_cols_sp = elsi_h%n_basis-&
                              (elsi_h%n_procs-1)*elsi_h%n_l_cols_sp
   endif

   call elsi_blacs_to_sips_hs_real(elsi_h,h_in,s_in)

   ! Then get the global matrices
   call elsi_sips_to_chess_hs(elsi_h)

   elsi_h%nnz_l_sp = elsi_h%nnz_g
   elsi_h%n_l_cols_sp = elsi_h%n_basis

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!>
!! This routine gets global matrices from 1D block distributed
!! sparse CCS format.
!!
subroutine elsi_sips_to_chess_hs(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)

   integer(kind=i4) :: prev_nnz
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_sips_to_chess_hs"

   ! Set recv_count and recv_displ
   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   call MPI_Allgather(elsi_h%nnz_l_sp,1,mpi_integer4,recv_count,1,&
           mpi_integer4,elsi_h%mpi_comm,mpierr)

   do i_proc = 2,elsi_h%n_procs
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Get the global matrices
   if(elsi_h%n_elsi_calls == 1) then
      if(.not. elsi_h%ovlp_is_unit) then
         ! Overlap value
         call elsi_allocate(elsi_h,elsi_h%ovlp_real_chess,elsi_h%nnz_g,&
                 "ovlp_real_chess",caller)

         call MPI_Allgatherv(elsi_h%ovlp_real_sips,elsi_h%nnz_l_sp,&
                 mpi_real8,elsi_h%ovlp_real_chess,recv_count,recv_displ,&
                 mpi_real8,elsi_h%mpi_comm,mpierr)

         call elsi_deallocate(elsi_h,elsi_h%ovlp_real_sips,"ovlp_real_sips")
      else
         call elsi_allocate(elsi_h,elsi_h%ovlp_real_chess,1,"dummy",caller)
      endif

      ! Row index
      call elsi_allocate(elsi_h,elsi_h%row_ind_chess,elsi_h%nnz_g,&
              "row_ind_chess",caller)

      call MPI_Allgatherv(elsi_h%row_ind_sips,elsi_h%nnz_l_sp,&
              mpi_integer4,elsi_h%row_ind_chess,recv_count,recv_displ,&
              mpi_integer4,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,elsi_h%row_ind_sips,"row_ind_sips")

      call elsi_allocate(elsi_h,elsi_h%ham_real_chess,elsi_h%nnz_g,&
              "ham_real_chess",caller)
   endif

   ! Hamiltonian value
   call MPI_Allgatherv(elsi_h%ham_real_sips,elsi_h%nnz_l_sp,&
           mpi_real8,elsi_h%ham_real_chess,recv_count,recv_displ,&
           mpi_real8,elsi_h%mpi_comm,mpierr)

   if(elsi_h%n_elsi_calls == 1) then
      ! Set recv_count and recv_displ
      recv_count(1:elsi_h%n_procs-1) = elsi_h%n_basis/elsi_h%n_procs
      recv_count(elsi_h%n_procs) = elsi_h%n_basis-&
                                      (elsi_h%n_procs-1)*recv_count(1)

      do i_proc = 2,elsi_h%n_procs
         recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
      enddo

      ! Shift column pointers
      prev_nnz = 0

      call MPI_Exscan(elsi_h%nnz_l_sp,prev_nnz,1,mpi_integer4,mpi_sum,&
              elsi_h%mpi_comm,mpierr)

      elsi_h%col_ptr_sips = elsi_h%col_ptr_sips+prev_nnz

      ! Column pointer
      call elsi_allocate(elsi_h,elsi_h%col_ptr_chess,elsi_h%n_basis+1,&
              "col_ptr_chess",caller)

      call MPI_Allgatherv(elsi_h%col_ptr_sips,elsi_h%n_l_cols_sp,&
              mpi_integer4,elsi_h%col_ptr_chess,recv_count,recv_displ,&
              mpi_integer4,elsi_h%mpi_comm,mpierr)

      elsi_h%col_ptr_chess(elsi_h%n_basis+1) = elsi_h%nnz_g+1

      call elsi_deallocate(elsi_h,elsi_h%col_ptr_sips,"col_ptr_sips")
   endif

   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

end subroutine

!>
!! This routine converts density matrix computed by CheSS and stored
!! in sparse CCS format to 2D block-cyclic distributed dense format.
!!
subroutine elsi_chess_to_blacs_dm_real(elsi_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                 !< Handle
   real(kind=r8),     intent(out)   :: d_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: g_row
   integer(kind=i4) :: g_col

   character*40, parameter :: caller = "elsi_chess_to_blacs_dm_real"

   d_out = 0.0_r8

   do i_col = 1,elsi_h%n_l_cols
      call elsi_get_global_col(elsi_h,g_col,i_col)

      do i_val = elsi_h%col_ptr_ccs(g_col),elsi_h%col_ptr_ccs(g_col+1)-1
         g_row = elsi_h%row_ind_ccs(i_val)

         if(elsi_h%local_row(g_row) == 0) cycle

         i_row = (g_row-1)/(elsi_h%n_p_rows*elsi_h%n_b_rows)*&
                    elsi_h%n_b_rows+mod((g_row-1),elsi_h%n_b_rows)+1

         d_out(i_row,i_col) = elsi_h%dm_chess%matrix_compr(i_val)
      enddo
   enddo

end subroutine

!>
!! This routine sorts an integer(kind=i8) array by heapsort,
!! moves two real(kind=r8) arrays and two integer(kind=i4)
!! arrays accordingly.
!!
subroutine heapsort_real_v1(length,a_i8,b_r8,c_r8,d_i4,e_i4)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length) !< i8 array to be sorted
   real(kind=i8)   , intent(inout) :: b_r8(length) !< r8 array to be moved
   real(kind=i8)   , intent(inout) :: c_r8(length) !< r8 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length) !< i4 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i8,b_r8,c_r8,d_i4,e_i4,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i8,1,i)
      call swap(length,b_r8,1,i)
      call swap(length,c_r8,1,i)
      call swap(length,d_i4,1,i)
      call swap(length,e_i4,1,i)

      i = i-1

      call downheap(length,a_i8,b_r8,c_r8,d_i4,e_i4,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_real_v1(length,a_i8,b_r8,c_r8,d_i4,e_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length) !< i8 array to be sorted
   real(kind=i8)   , intent(inout) :: b_r8(length) !< r8 array to be moved
   real(kind=i8)   , intent(inout) :: c_r8(length) !< r8 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length) !< i4 array to be moved
   integer(kind=i4), intent(in)    :: top          !< Top of heap
   integer(kind=i4), intent(in)    :: bottom       !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i8(w+1) > a_i8(w)) then
            w = w+1
         endif
      endif

      if(a_i8(v) >= a_i8(w)) then
         return
      else
         call swap(length,a_i8,v,w)
         call swap(length,b_r8,v,w)
         call swap(length,c_r8,v,w)
         call swap(length,d_i4,v,w)
         call swap(length,e_i4,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine sorts an integer(kind=i8) array by heapsort,
!! moves a real(kind=r8) array accordingly.
!!
subroutine heapsort_real_v2(length,a_i8,b_r8)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length) !< i8 array to be sorted
   real(kind=i8)   , intent(inout) :: b_r8(length) !< r8 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i8,b_r8,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i8,1,i)
      call swap(length,b_r8,1,i)

      i = i-1

      call downheap(length,a_i8,b_r8,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_real_v2(length,a_i8,b_r8,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length       !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length) !< i8 array to be sorted
   real(kind=i8)   , intent(inout) :: b_r8(length) !< r8 array to be moved
   integer(kind=i4), intent(in)    :: top          !< Top of heap
   integer(kind=i4), intent(in)    :: bottom       !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i8(w+1) > a_i8(w)) then
            w = w+1
         endif
      endif

      if(a_i8(v) >= a_i8(w)) then
         return
      else
         call swap(length,a_i8,v,w)
         call swap(length,b_r8,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine sorts an integer(kind=i8) array by heapsort,
!! moves two complex(kind=r8) arrays and two integer(kind=i4)
!! arrays accordingly.
!!
subroutine heapsort_complex_v1(length,a_i8,b_c16,c_c16,d_i4,e_i4)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length)  !< i8 array to be sorted
   complex(kind=i8), intent(inout) :: b_c16(length) !< c16 array to be moved
   complex(kind=i8), intent(inout) :: c_c16(length) !< c16 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length)  !< i4 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i8,b_c16,c_c16,d_i4,e_i4,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i8,1,i)
      call swap(length,b_c16,1,i)
      call swap(length,c_c16,1,i)
      call swap(length,d_i4,1,i)
      call swap(length,e_i4,1,i)

      i = i-1

      call downheap(length,a_i8,b_c16,c_c16,d_i4,e_i4,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_complex_v1(length,a_i8,b_c16,c_c16,d_i4,e_i4,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length)  !< i8 array to be sorted
   complex(kind=i8), intent(inout) :: b_c16(length) !< c16 array to be moved
   complex(kind=i8), intent(inout) :: c_c16(length) !< c16 array to be moved
   integer(kind=i4), intent(inout) :: d_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(inout) :: e_i4(length)  !< i4 array to be moved
   integer(kind=i4), intent(in)    :: top           !< Top of heap
   integer(kind=i4), intent(in)    :: bottom        !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i8(w+1) > a_i8(w)) then
            w = w+1
         endif
      endif

      if(a_i8(v) >= a_i8(w)) then
         return
      else
         call swap(length,a_i8,v,w)
         call swap(length,b_c16,v,w)
         call swap(length,c_c16,v,w)
         call swap(length,d_i4,v,w)
         call swap(length,e_i4,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine sorts an integer(kind=i8) array by heapsort,
!! moves a complex(kind=r8) array accordingly.
!!
subroutine heapsort_complex_v2(length,a_i8,b_c16)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length)  !< i8 array to be sorted
   complex(kind=i8), intent(inout) :: b_c16(length) !< c16 array to be moved

   integer(kind=i4) :: top
   integer(kind=i4) :: i

   ! Heapify
   top = length/2

   do i = top,1,-1
      call downheap(length,a_i8,b_c16,i,length)
   enddo

   i = length

   do while(i > 1)
      call swap(length,a_i8,1,i)
      call swap(length,b_c16,1,i)

      i = i-1

      call downheap(length,a_i8,b_c16,1,i)
   enddo

end subroutine

!>
!! This routine restores the max-heap structure.
!!
subroutine downheap_complex_v2(length,a_i8,b_c16,top,bottom)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i8), intent(inout) :: a_i8(length)  !< i8 array to be sorted
   complex(kind=i8), intent(inout) :: b_c16(length) !< c16 array to be moved
   integer(kind=i4), intent(in)    :: top           !< Top of heap
   integer(kind=i4), intent(in)    :: bottom        !< Bottom of heap

   integer(kind=i4) :: v
   integer(kind=i4) :: w

   v = top
   w = 2*v

   do while (w <= bottom)
      if(w+1 <= bottom) then
         if(a_i8(w+1) > a_i8(w)) then
            w = w+1
         endif
      endif

      if(a_i8(v) >= a_i8(w)) then
         return
      else
         call swap(length,a_i8,v,w)
         call swap(length,b_c16,v,w)

         v = w
         w = 2*v
      endif
   enddo

end subroutine

!>
!! This routine swaps two numbers in an integer(kind=i8) array.
!!
subroutine swap_i8(length,array,i,j)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i8), intent(inout) :: array(length) !< i8 array
   integer(kind=i4), intent(in)    :: i             !< Index to be swapped
   integer(kind=i4), intent(in)    :: j             !< Index to be swapped

   integer(kind=i8) :: tmp

   tmp      = array(i)
   array(i) = array(j)
   array(j) = tmp

end subroutine

!>
!! This routine swaps two numbers in an integer(kind=i4) array.
!!
subroutine swap_i4(length,array,i,j)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   integer(kind=i4), intent(inout) :: array(length) !< i4 array
   integer(kind=i4), intent(in)    :: i             !< Index to be swapped
   integer(kind=i4), intent(in)    :: j             !< Index to be swapped

   integer(kind=i4) :: tmp

   tmp      = array(i)
   array(i) = array(j)
   array(j) = tmp

end subroutine

!>
!! This routine swaps two numbers in a real(kind=r8) array.
!!
subroutine swap_r8(length,array,i,j)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   real(kind=r8),    intent(inout) :: array(length) !< r8 array
   integer(kind=i4), intent(in)    :: i             !< Index to be swapped
   integer(kind=i4), intent(in)    :: j             !< Index to be swapped

   real(kind=r8) :: tmp

   tmp      = array(i)
   array(i) = array(j)
   array(j) = tmp

end subroutine

!>
!! This routine swaps two numbers in a complex(kind=r8) array.
!!
subroutine swap_c16(length,array,i,j)

   implicit none

   integer(kind=i4), intent(in)    :: length        !< Length of array
   complex(kind=r8), intent(inout) :: array(length) !< c16 array
   integer(kind=i4), intent(in)    :: i             !< Index to be swapped
   integer(kind=i4), intent(in)    :: j             !< Index to be swapped

   complex(kind=r8) :: tmp

   tmp      = array(i)
   array(i) = array(j)
   array(j) = tmp

end subroutine

end module ELSI_MATCONV
