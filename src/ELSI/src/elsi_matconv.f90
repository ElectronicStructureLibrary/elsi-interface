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

   use ELSI_CONSTANTS, only: ELPA_SOLVER,OMM_SOLVER,DMP_SOLVER
   use ELSI_DATATYPE
   use ELSI_MALLOC
   use ELSI_PRECISION, only: r8,i4,i8
   use ELSI_SORT
   use ELSI_UTILS

   implicit none

   private

   public :: elsi_blacs_to_pexsi_hs
   public :: elsi_pexsi_to_blacs_dm
   public :: elsi_blacs_to_sips_hs
   public :: elsi_sips_to_blacs_hs
   public :: elsi_blacs_to_sips_dm
   public :: elsi_blacs_to_chess_hs
   public :: elsi_chess_to_blacs_dm

   interface elsi_blacs_to_chess_hs
      module procedure elsi_blacs_to_chess_hs_real
   end interface

   interface elsi_blacs_to_sips_dm
      module procedure elsi_blacs_to_sips_dm_real,&
                       elsi_blacs_to_sips_dm_complex
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
      module procedure elsi_sips_to_blacs_hs_real,&
                       elsi_sips_to_blacs_hs_complex
   end interface

contains

!>
!! This routine converts Halmitonian and overlap matrixs stored in 2D block-
!! cyclic dense format to 1D block sparse CCS format, which can be used as input
!! by PEXSI.
!!
subroutine elsi_blacs_to_pexsi_hs_real(e_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout)        :: e_h                         !< Handle
   real(kind=r8),     intent(in),   target :: h_in(e_h%n_lrow,e_h%n_lcol) !< Hamiltonian
   real(kind=r8),     intent(in),   target :: s_in(e_h%n_lrow,e_h%n_lcol) !< Overlap

   real(kind=r8), pointer :: ref(:,:) ! Sparsity pattern

   integer(kind=i4) :: n_para_task
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: g_col_id
   integer(kind=i4) :: g_row_id
   integer(kind=i4) :: d1
   integer(kind=i4) :: d2
   integer(kind=i4) :: d11
   integer(kind=i4) :: d12
   integer(kind=i4) :: d21
   integer(kind=i4) :: d22
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: this_n_cols
   integer(kind=i4) :: nnz_l_aux
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

   call elsi_get_time(e_h,t0)

   n_para_task = e_h%n_procs/e_h%np_per_pole

   if(.not. e_h%ovlp_is_unit) then
      ref => s_in
   else
      ref => h_in
   endif

   if(e_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz(e_h,ref,e_h%n_lrow,e_h%n_lcol,e_h%nnz_l)

      if(.not. e_h%ovlp_is_unit) then
         call elsi_allocate(e_h,s_val_send_buf,e_h%nnz_l,"s_val_send_buf",&
                 caller)
      endif
   endif

   call elsi_allocate(e_h,locat,e_h%n_basis,"locat",caller)
   call elsi_allocate(e_h,row_send_buf,e_h%nnz_l,"row_send_buf",caller)
   call elsi_allocate(e_h,col_send_buf,e_h%nnz_l,"col_send_buf",caller)
   call elsi_allocate(e_h,h_val_send_buf,e_h%nnz_l,"h_val_send_buf",caller)

   ! Compute d1,d2,d11,d12,d21,d22 (need explanation)
   d1  = e_h%n_basis/e_h%np_per_pole
   d2  = e_h%n_basis-(e_h%np_per_pole-1)*d1
   d11 = d1/n_para_task
   d12 = d1-(n_para_task-1)*d11
   d21 = d2/n_para_task
   d22 = d2-(n_para_task-1)*d21

   i_val = 0

   if(d11 > 0) then
      do i_proc = 0,e_h%n_procs-n_para_task-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d12
         else
            this_n_cols = d11
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = 0,e_h%n_procs-n_para_task-1
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
      do i_proc = e_h%n_procs-n_para_task,e_h%n_procs-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d22
         else
            this_n_cols = d21
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = e_h%n_procs-n_para_task,e_h%n_procs-1
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

   call elsi_allocate(e_h,send_count,e_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0
   do i_col = 1,e_h%n_lcol
      call elsi_get_global_col(e_h,g_col_id,i_col)

      ! Compute destination
      dest = min(locat(g_col_id),e_h%n_procs-1)

      do i_row = 1,e_h%n_lrow
         if(abs(ref(i_row,i_col)) > e_h%zero_def) then
            i_val = i_val+1

            call elsi_get_global_row(e_h,g_row_id,i_row)

            ! Pack global id and data into bufs
            row_send_buf(i_val)   = g_row_id
            col_send_buf(i_val)   = g_col_id
            h_val_send_buf(i_val) = h_in(i_row,i_col)
            if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
               s_val_send_buf(i_val) = s_in(i_row,i_col)
            endif

            ! Set send_count
            send_count(dest+1) = send_count(dest+1)+1
         endif
      enddo
   enddo

   nullify(ref)
   call elsi_deallocate(e_h,locat,"locat")
   call elsi_allocate(e_h,recv_count,e_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   nnz_l_aux = sum(recv_count,1)
   call MPI_Allreduce(nnz_l_aux,e_h%nnz_g,1,mpi_integer4,mpi_sum,e_h%mpi_comm,&
           mpierr)

   ! Set send and receive displacement
   call elsi_allocate(e_h,send_displ,e_h%n_procs,"send_displ",caller)
   call elsi_allocate(e_h,recv_displ,e_h%n_procs,"recv_displ",caller)

   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(e_h,row_recv_buf,nnz_l_aux,"row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,mpi_integer4,&
           row_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,row_send_buf,"row_send_buf")

   ! Column id
   call elsi_allocate(e_h,col_recv_buf,nnz_l_aux,"col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,mpi_integer4,&
           col_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,col_send_buf,"col_send_buf")

   ! Hamiltonian value
   call elsi_allocate(e_h,h_val_recv_buf,nnz_l_aux,"h_val_recv_buf",caller)

   call MPI_Alltoallv(h_val_send_buf,send_count,send_displ,mpi_real8,&
           h_val_recv_buf,recv_count,recv_displ,mpi_real8,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,h_val_send_buf,"h_val_send_buf")

   ! Overlap value
   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_allocate(e_h,s_val_recv_buf,nnz_l_aux,"s_val_recv_buf",caller)

      call MPI_Alltoallv(s_val_send_buf,send_count,send_displ,mpi_real8,&
              s_val_recv_buf,recv_count,recv_displ,mpi_real8,e_h%mpi_comm,&
              mpierr)

      call elsi_deallocate(e_h,s_val_send_buf,"s_val_send_buf")
   endif

   call elsi_allocate(e_h,global_id,nnz_l_aux,"global_id",caller)

   ! Compute global 1D id
   do i_val = 1,nnz_l_aux
      global_id(i_val) = int(col_recv_buf(i_val)-1,kind=i8)*&
                            int(e_h%n_basis,kind=i8)+&
                            int(row_recv_buf(i_val),kind=i8)
   enddo

   ! Sort
   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_heapsort(nnz_l_aux,global_id,h_val_recv_buf,s_val_recv_buf,&
              row_recv_buf,col_recv_buf)
   else
      call elsi_heapsort(nnz_l_aux,global_id,h_val_recv_buf,row_recv_buf,&
              col_recv_buf)
   endif

   call elsi_deallocate(e_h,global_id,"global_id")

   ! Set send_count, all data sent to the first pole
   send_count = 0
   send_count(e_h%myid/n_para_task+1) = nnz_l_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   if(e_h%n_elsi_calls == 1) then
      ! Only the first pole knows nnz_l_sp
      e_h%nnz_l_sp = sum(recv_count,1)

      call MPI_Bcast(e_h%nnz_l_sp,1,mpi_integer4,0,e_h%comm_among_pole,mpierr)
   endif

   ! Set send and receive displacement
   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   if(e_h%n_elsi_calls == 1) then
      ! Row id
      call elsi_allocate(e_h,row_send_buf,e_h%nnz_l_sp,"row_send_buf",caller)

      call MPI_Alltoallv(row_recv_buf,send_count,send_displ,mpi_integer4,&
              row_send_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,&
              mpierr)

      call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")

      ! Column id
      call elsi_allocate(e_h,col_send_buf,e_h%nnz_l_sp,"col_send_buf",caller)

      call MPI_Alltoallv(col_recv_buf,send_count,send_displ,mpi_integer4,&
              col_send_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,&
              mpierr)

      call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")

      if(.not. e_h%ovlp_is_unit) then
         ! Overlap value
         call elsi_allocate(e_h,e_h%ovlp_real_pexsi,e_h%nnz_l_sp,&
                 "ovlp_real_pexsi",caller)

         call MPI_Alltoallv(s_val_recv_buf,send_count,send_displ,mpi_real8,&
                 e_h%ovlp_real_pexsi,recv_count,recv_displ,mpi_real8,&
                 e_h%mpi_comm,mpierr)

         call elsi_deallocate(e_h,s_val_recv_buf,"s_val_recv_buf")

         call elsi_allocate(e_h,e_h%ham_real_pexsi,e_h%nnz_l_sp,&
                 "ham_real_pexsi",caller)
      else
         call elsi_allocate(e_h,e_h%ovlp_real_pexsi,1,"dummy",caller)
      endif
   else
      call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")
      call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")
   endif

   e_h%ham_real_pexsi = 0.0_r8

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv_buf,send_count,send_displ,mpi_real8,&
           e_h%ham_real_pexsi,recv_count,recv_displ,mpi_real8,e_h%mpi_comm,&
           mpierr)

   call elsi_deallocate(e_h,h_val_recv_buf,"h_val_recv_buf")
   call elsi_deallocate(e_h,send_count,"send_count")
   call elsi_deallocate(e_h,recv_count,"recv_count")
   call elsi_deallocate(e_h,send_displ,"send_displ")
   call elsi_deallocate(e_h,recv_displ,"recv_displ")

   if(e_h%n_elsi_calls == 1) then
      call elsi_allocate(e_h,e_h%row_ind_pexsi,e_h%nnz_l_sp,"row_ind_pexsi",&
              caller)

      call elsi_allocate(e_h,e_h%col_ptr_pexsi,(e_h%n_lcol_sp+1),&
              "col_ptr_pexsi",caller)

      ! Only the first pole computes row index and column pointer
      if(e_h%my_prow_pexsi == 0) then
         ! Compute row index and column pointer
         i_col = col_send_buf(1)-1
         do i_val = 1,e_h%nnz_l_sp
            e_h%row_ind_pexsi(i_val) = row_send_buf(i_val)

            if(col_send_buf(i_val) > i_col) then
               i_col = i_col+1
               e_h%col_ptr_pexsi(i_col-col_send_buf(1)+1) = i_val
            endif
         enddo

         e_h%col_ptr_pexsi(e_h%n_lcol_sp+1) = e_h%nnz_l_sp+1
      endif

      call elsi_deallocate(e_h,row_send_buf,"row_send_buf")
      call elsi_deallocate(e_h,col_send_buf,"col_send_buf")
   endif

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 2D block-
!! cyclic dense format to 1D block sparse CCS format, which can be used as input
!! by PEXSI.
!!
subroutine elsi_blacs_to_pexsi_hs_complex(e_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout)        :: e_h                         !< Handle
   complex(kind=r8),  intent(in),   target :: h_in(e_h%n_lrow,e_h%n_lcol) !< Hamiltonian
   complex(kind=r8),  intent(in),   target :: s_in(e_h%n_lrow,e_h%n_lcol) !< Overlap

   complex(kind=r8), pointer :: ref(:,:) ! Sparsity pattern

   integer(kind=i4) :: n_para_task
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: g_col_id
   integer(kind=i4) :: g_row_id
   integer(kind=i4) :: d1
   integer(kind=i4) :: d2
   integer(kind=i4) :: d11
   integer(kind=i4) :: d12
   integer(kind=i4) :: d21
   integer(kind=i4) :: d22
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: this_n_cols
   integer(kind=i4) :: nnz_l_aux
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

   call elsi_get_time(e_h,t0)

   n_para_task = e_h%n_procs/e_h%np_per_pole

   if(.not. e_h%ovlp_is_unit) then
      ref => s_in
   else
      ref => h_in
   endif

   if(e_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz(e_h,ref,e_h%n_lrow,e_h%n_lcol,e_h%nnz_l)

      if(.not. e_h%ovlp_is_unit) then
         call elsi_allocate(e_h,s_val_send_buf,e_h%nnz_l,"s_val_send_buf",&
                 caller)
      endif
   endif

   call elsi_allocate(e_h,locat,e_h%n_basis,"locat",caller)
   call elsi_allocate(e_h,row_send_buf,e_h%nnz_l,"row_send_buf",caller)
   call elsi_allocate(e_h,col_send_buf,e_h%nnz_l,"col_send_buf",caller)
   call elsi_allocate(e_h,h_val_send_buf,e_h%nnz_l,"h_val_send_buf",caller)

   ! Compute d1,d2,d11,d12,d21,d22 (need explanation)
   d1  = e_h%n_basis/e_h%np_per_pole
   d2  = e_h%n_basis-(e_h%np_per_pole-1)*d1
   d11 = d1/n_para_task
   d12 = d1-(n_para_task-1)*d11
   d21 = d2/n_para_task
   d22 = d2-(n_para_task-1)*d21

   i_val = 0

   if(d11 > 0) then
      do i_proc = 0,e_h%n_procs-n_para_task-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d12
         else
            this_n_cols = d11
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = 0,e_h%n_procs-n_para_task-1
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
      do i_proc = e_h%n_procs-n_para_task,e_h%n_procs-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d22
         else
            this_n_cols = d21
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = e_h%n_procs-n_para_task,e_h%n_procs-1
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

   call elsi_allocate(e_h,send_count,e_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0
   do i_col = 1,e_h%n_lcol
      call elsi_get_global_col(e_h,g_col_id,i_col)

      ! Compute destination
      dest = min(locat(g_col_id),e_h%n_procs-1)

      do i_row = 1,e_h%n_lrow
         if(abs(ref(i_row,i_col)) > e_h%zero_def) then
            i_val = i_val+1

            call elsi_get_global_row(e_h,g_row_id,i_row)

            ! Pack global id and data into bufs
            row_send_buf(i_val)   = g_row_id
            col_send_buf(i_val)   = g_col_id
            h_val_send_buf(i_val) = h_in(i_row,i_col)
            if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
               s_val_send_buf(i_val) = s_in(i_row,i_col)
            endif

            ! Set send_count
            send_count(dest+1) = send_count(dest+1)+1
         endif
      enddo
   enddo

   nullify(ref)
   call elsi_deallocate(e_h,locat,"locat")
   call elsi_allocate(e_h,recv_count,e_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   nnz_l_aux = sum(recv_count,1)
   call MPI_Allreduce(nnz_l_aux,e_h%nnz_g,1,mpi_integer4,mpi_sum,e_h%mpi_comm,&
           mpierr)

   ! Set send and receive displacement
   call elsi_allocate(e_h,send_displ,e_h%n_procs,"send_displ",caller)
   call elsi_allocate(e_h,recv_displ,e_h%n_procs,"recv_displ",caller)

   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(e_h,row_recv_buf,nnz_l_aux,"row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,mpi_integer4,&
           row_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,row_send_buf,"row_send_buf")

   ! Column id
   call elsi_allocate(e_h,col_recv_buf,nnz_l_aux,"col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,mpi_integer4,&
           col_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,col_send_buf,"col_send_buf")

   ! Hamiltonian value
   call elsi_allocate(e_h,h_val_recv_buf,nnz_l_aux,"h_val_recv_buf",caller)

   call MPI_Alltoallv(h_val_send_buf,send_count,send_displ,mpi_complex16,&
           h_val_recv_buf,recv_count,recv_displ,mpi_complex16,e_h%mpi_comm,&
           mpierr)

   call elsi_deallocate(e_h,h_val_send_buf,"h_val_send_buf")

   ! Overlap value
   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_allocate(e_h,s_val_recv_buf,nnz_l_aux,"s_val_recv_buf",caller)

      call MPI_Alltoallv(s_val_send_buf,send_count,send_displ,mpi_complex16,&
              s_val_recv_buf,recv_count,recv_displ,mpi_complex16,e_h%mpi_comm,&
              mpierr)

      call elsi_deallocate(e_h,s_val_send_buf,"s_val_send_buf")
   endif

   call elsi_allocate(e_h,global_id,nnz_l_aux,"global_id",caller)

   ! Compute global 1D id
   do i_val = 1,nnz_l_aux
      global_id(i_val) = int(col_recv_buf(i_val)-1,kind=i8)*&
                            int(e_h%n_basis,kind=i8)+&
                            int(row_recv_buf(i_val),kind=i8)
   enddo

   ! Sort
   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_heapsort(nnz_l_aux,global_id,h_val_recv_buf,s_val_recv_buf,&
              row_recv_buf,col_recv_buf)
   else
      call elsi_heapsort(nnz_l_aux,global_id,h_val_recv_buf,row_recv_buf,&
              col_recv_buf)
   endif

   call elsi_deallocate(e_h,global_id,"global_id")

   ! Set send_count, all data sent to the first pole
   send_count = 0
   send_count(e_h%myid/n_para_task+1) = nnz_l_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   if(e_h%n_elsi_calls == 1) then
      ! Only the first pole knows nnz_l_sp
      e_h%nnz_l_sp = sum(recv_count,1)

      call MPI_Bcast(e_h%nnz_l_sp,1,mpi_integer4,0,e_h%comm_among_pole,mpierr)
   endif

   ! Set send and receive displacement
   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   if(e_h%n_elsi_calls == 1) then
      ! Row id
      call elsi_allocate(e_h,row_send_buf,e_h%nnz_l_sp,"row_send_buf",caller)

      call MPI_Alltoallv(row_recv_buf,send_count,send_displ,mpi_integer4,&
              row_send_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,&
              mpierr)

      call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")

      ! Column id
      call elsi_allocate(e_h,col_send_buf,e_h%nnz_l_sp,"col_send_buf",caller)

      call MPI_Alltoallv(col_recv_buf,send_count,send_displ,mpi_integer4,&
              col_send_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,&
              mpierr)

      call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")

      if(.not. e_h%ovlp_is_unit) then
         ! Overlap value
         call elsi_allocate(e_h,e_h%ovlp_cmplx_pexsi,e_h%nnz_l_sp,&
                 "ovlp_cmplx_pexsi",caller)

         call MPI_Alltoallv(s_val_recv_buf,send_count,send_displ,mpi_complex16,&
                 e_h%ovlp_cmplx_pexsi,recv_count,recv_displ,mpi_complex16,&
                 e_h%mpi_comm,mpierr)

         call elsi_deallocate(e_h,s_val_recv_buf,"s_val_recv_buf")

         call elsi_allocate(e_h,e_h%ham_cmplx_pexsi,e_h%nnz_l_sp,&
                 "ham_cmplx_pexsi",caller)
      else
         call elsi_allocate(e_h,e_h%ovlp_cmplx_pexsi,1,"dummy",caller)
      endif
   else
      call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")
      call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")
   endif

   e_h%ham_cmplx_pexsi = (0.0_r8,0.0_r8)

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv_buf,send_count,send_displ,mpi_complex16,&
           e_h%ham_cmplx_pexsi,recv_count,recv_displ,mpi_complex16,&
           e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,h_val_recv_buf,"h_val_recv_buf")
   call elsi_deallocate(e_h,send_count,"send_count")
   call elsi_deallocate(e_h,recv_count,"recv_count")
   call elsi_deallocate(e_h,send_displ,"send_displ")
   call elsi_deallocate(e_h,recv_displ,"recv_displ")

   if(e_h%n_elsi_calls == 1) then
      call elsi_allocate(e_h,e_h%row_ind_pexsi,e_h%nnz_l_sp,"row_ind_pexsi",&
              caller)

      call elsi_allocate(e_h,e_h%col_ptr_pexsi,(e_h%n_lcol_sp+1),&
              "col_ptr_pexsi",caller)

      ! Only the first pole computes row index and column pointer
      if(e_h%my_prow_pexsi == 0) then
         ! Compute row index and column pointer
         i_col = col_send_buf(1)-1
         do i_val = 1,e_h%nnz_l_sp
            e_h%row_ind_pexsi(i_val) = row_send_buf(i_val)

            if(col_send_buf(i_val) > i_col) then
               i_col = i_col+1
               e_h%col_ptr_pexsi(i_col-col_send_buf(1)+1) = i_val
            endif
         enddo

         e_h%col_ptr_pexsi(e_h%n_lcol_sp+1) = e_h%nnz_l_sp+1
      endif

      call elsi_deallocate(e_h,row_send_buf,"row_send_buf")
      call elsi_deallocate(e_h,col_send_buf,"col_send_buf")
   endif

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine converts density matrix computed by PEXSI, stored in 1D block
!! sparse CCS format to 2D block-cyclic dense format.
!!
subroutine elsi_pexsi_to_blacs_dm_real(e_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                          !< Handle
   real(kind=r8),     intent(out)   :: d_out(e_h%n_lrow,e_h%n_lcol) !< Density matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: l_row_id ! Local row id in 1D block distribution
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
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_dm_real"

   call elsi_get_time(e_h,t0)

   call elsi_allocate(e_h,val_send_buf,e_h%nnz_l_sp,"val_send_buf",caller)
   call elsi_allocate(e_h,row_send_buf,e_h%nnz_l_sp,"row_send_buf",caller)
   call elsi_allocate(e_h,col_send_buf,e_h%nnz_l_sp,"col_send_buf",caller)
   call elsi_allocate(e_h,send_count,e_h%n_procs,"send_count",caller)

   if(e_h%my_prow_pexsi == 0) then
      call elsi_allocate(e_h,dest,e_h%nnz_l_sp,"dest",caller)

      i_col = 0
      ! Compute destination and global id
      do i_val = 1,e_h%nnz_l_sp
         if(i_val == e_h%col_ptr_ccs(i_col+1) .and. i_col /= e_h%n_lcol_sp) then
            i_col = i_col+1
         endif
         i_row = e_h%row_ind_ccs(i_val)

         ! Compute global id
         row_send_buf(i_val) = i_row
         col_send_buf(i_val) = i_col+e_h%myid*(e_h%n_basis/e_h%np_per_pole)
         val_send_buf(i_val) = e_h%dm_real_ccs(i_val)

         ! Compute destination
         proc_row_id = mod((row_send_buf(i_val)-1)/e_h%blk_row,e_h%n_prow)
         proc_col_id = mod((col_send_buf(i_val)-1)/e_h%blk_col,e_h%n_pcol)
         dest(i_val) = proc_col_id+proc_row_id*e_h%n_pcol

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      enddo

      call elsi_heapsort(e_h%nnz_l_sp,dest,val_send_buf,row_send_buf,&
              col_send_buf)

      call elsi_deallocate(e_h,dest,"dest")
   endif

   call elsi_allocate(e_h,recv_count,e_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   e_h%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(e_h,send_displ,e_h%n_procs,"send_displ",caller)
   call elsi_allocate(e_h,recv_displ,e_h%n_procs,"recv_displ",caller)

   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Value
   call elsi_allocate(e_h,val_recv_buf,e_h%nnz_l,"val_recv_buf",caller)

   call MPI_Alltoallv(val_send_buf,send_count,send_displ,mpi_real8,&
           val_recv_buf,recv_count,recv_displ,mpi_real8,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,val_send_buf,"val_send_buf")

   ! Row index
   call elsi_allocate(e_h,row_recv_buf,e_h%nnz_l,"row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,mpi_integer4,&
           row_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,row_send_buf,"row_send_buf")

   ! Column index
   call elsi_allocate(e_h,col_recv_buf,e_h%nnz_l,"col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,mpi_integer4,&
           col_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,col_send_buf,"col_send_buf")
   call elsi_deallocate(e_h,send_count,"send_count")
   call elsi_deallocate(e_h,recv_count,"recv_count")
   call elsi_deallocate(e_h,send_displ,"send_displ")
   call elsi_deallocate(e_h,recv_displ,"recv_displ")

   d_out = 0.0_r8

   ! Unpack density matrix
   do i_val = 1,e_h%nnz_l
      ! Compute local 2d id
      l_row_id = (row_recv_buf(i_val)-1)/(e_h%n_prow*e_h%blk_row)*e_h%blk_row+&
                    mod((row_recv_buf(i_val)-1),e_h%blk_row)+1
      l_col_id = (col_recv_buf(i_val)-1)/(e_h%n_pcol*e_h%blk_col)*e_h%blk_col+&
                    mod((col_recv_buf(i_val)-1),e_h%blk_col)+1

      ! Put value to correct position
      d_out(l_row_id,l_col_id) = val_recv_buf(i_val)
   enddo

   call elsi_deallocate(e_h,val_recv_buf,"val_recv_buf")
   call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine converts density matrix computed by PEXSI, stored in 1D block
!! sparse CCS format to 2D block-cyclic dense format.
!!
subroutine elsi_pexsi_to_blacs_dm_complex(e_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                          !< Handle
   complex(kind=r8),  intent(out)   :: d_out(e_h%n_lrow,e_h%n_lcol) !< Density matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: l_row_id ! Local row id in 1D block distribution
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
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_dm_complex"

   call elsi_get_time(e_h,t0)

   call elsi_allocate(e_h,val_send_buf,e_h%nnz_l_sp,"val_send_buf",caller)
   call elsi_allocate(e_h,row_send_buf,e_h%nnz_l_sp,"row_send_buf",caller)
   call elsi_allocate(e_h,col_send_buf,e_h%nnz_l_sp,"col_send_buf",caller)
   call elsi_allocate(e_h,send_count,e_h%n_procs,"send_count",caller)

   if(e_h%my_prow_pexsi == 0) then
      call elsi_allocate(e_h,dest,e_h%nnz_l_sp,"dest",caller)

      i_col = 0
      ! Compute destination and global id
      do i_val = 1,e_h%nnz_l_sp
         if(i_val == e_h%col_ptr_ccs(i_col+1) .and. i_col /= e_h%n_lcol_sp) then
            i_col = i_col+1
         endif
         i_row = e_h%row_ind_ccs(i_val)

         ! Compute global id
         row_send_buf(i_val) = i_row
         col_send_buf(i_val) = i_col+e_h%myid*(e_h%n_basis/e_h%np_per_pole)
         val_send_buf(i_val) = e_h%dm_cmplx_ccs(i_val)

         ! Compute destination
         proc_row_id = mod((row_send_buf(i_val)-1)/e_h%blk_row,e_h%n_prow)
         proc_col_id = mod((col_send_buf(i_val)-1)/e_h%blk_col,e_h%n_pcol)
         dest(i_val) = proc_col_id+proc_row_id*e_h%n_pcol

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      enddo

      call elsi_heapsort(e_h%nnz_l_sp,dest,val_send_buf,row_send_buf,&
              col_send_buf)

      call elsi_deallocate(e_h,dest,"dest")
   endif

   call elsi_allocate(e_h,recv_count,e_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   e_h%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(e_h,send_displ,e_h%n_procs,"send_displ",caller)
   call elsi_allocate(e_h,recv_displ,e_h%n_procs,"recv_displ",caller)

   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Value
   call elsi_allocate(e_h,val_recv_buf,e_h%nnz_l,"val_recv_buf",caller)

   call MPI_Alltoallv(val_send_buf,send_count,send_displ,mpi_complex16,&
           val_recv_buf,recv_count,recv_displ,mpi_complex16,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,val_send_buf,"val_send_buf")

   ! Row index
   call elsi_allocate(e_h,row_recv_buf,e_h%nnz_l,"row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,mpi_integer4,&
           row_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,row_send_buf,"row_send_buf")

   ! Column index
   call elsi_allocate(e_h,col_recv_buf,e_h%nnz_l,"col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,mpi_integer4,&
           col_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,col_send_buf,"col_send_buf")
   call elsi_deallocate(e_h,send_count,"send_count")
   call elsi_deallocate(e_h,recv_count,"recv_count")
   call elsi_deallocate(e_h,send_displ,"send_displ")
   call elsi_deallocate(e_h,recv_displ,"recv_displ")

   d_out = 0.0_r8

   ! Unpack density matrix
   do i_val = 1,e_h%nnz_l
      ! Compute local 2d id
      l_row_id = (row_recv_buf(i_val)-1)/(e_h%n_prow*e_h%blk_row)*e_h%blk_row+&
                    mod((row_recv_buf(i_val)-1),e_h%blk_row)+1
      l_col_id = (col_recv_buf(i_val)-1)/(e_h%n_pcol*e_h%blk_col)*e_h%blk_col+&
                    mod((col_recv_buf(i_val)-1),e_h%blk_col)+1

      ! Put value to correct position
      d_out(l_row_id,l_col_id) = val_recv_buf(i_val)
   enddo

   call elsi_deallocate(e_h,val_recv_buf,"val_recv_buf")
   call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 2D block-
!! cyclic dense format to 1D block sparse CCS format, which can be used as input
!! by SIPs.
!!
subroutine elsi_blacs_to_sips_hs_real(e_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout)        :: e_h                         !< Handle
   real(kind=r8),     intent(in),   target :: h_in(e_h%n_lrow,e_h%n_lcol) !< Hamiltonian
   real(kind=r8),     intent(in),   target :: s_in(e_h%n_lrow,e_h%n_lcol) !< Overlap

   real(kind=r8), pointer :: ref(:,:) ! Sparsity pattern

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: g_col_id
   integer(kind=i4) :: g_row_id
   integer(kind=i4) :: dest ! Destination of an element
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

   call elsi_get_time(e_h,t0)

   if(e_h%sips_n_elpa == UNSET) then
      e_h%sips_n_elpa = 0
   endif

   if(.not. e_h%ovlp_is_unit) then
      ref => s_in
   else
      ref => h_in
   endif

   if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa) then
      call elsi_get_local_nnz(e_h,ref,e_h%n_lrow,e_h%n_lcol,e_h%nnz_l)

      if(.not. e_h%ovlp_is_unit) then
         call elsi_allocate(e_h,s_val_send_buf,e_h%nnz_l,"s_val_send_buf",&
                 caller)
      endif
   endif

   call elsi_allocate(e_h,row_send_buf,e_h%nnz_l,"row_send_buf",caller)
   call elsi_allocate(e_h,col_send_buf,e_h%nnz_l,"col_send_buf",caller)
   call elsi_allocate(e_h,h_val_send_buf,e_h%nnz_l,"h_val_send_buf",caller)
   call elsi_allocate(e_h,send_count,e_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0
   do i_col = 1,e_h%n_lcol
      call elsi_get_global_col(e_h,g_col_id,i_col)

      ! Compute destination
      dest = (g_col_id-1)/(e_h%n_basis/e_h%n_procs)
      ! The last process may take more
      dest = min(dest,e_h%n_procs-1)

      do i_row = 1,e_h%n_lrow
         if(abs(ref(i_row,i_col)) > e_h%zero_def) then
            i_val = i_val+1

            call elsi_get_global_row(e_h,g_row_id,i_row)

            ! Pack global id and data into bufs
            row_send_buf(i_val)   = g_row_id
            col_send_buf(i_val)   = g_col_id
            h_val_send_buf(i_val) = h_in(i_row,i_col)
            if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa &
               .and. .not. e_h%ovlp_is_unit) then
               s_val_send_buf(i_val) = s_in(i_row,i_col)
            endif

            ! Set send_count
            send_count(dest+1) = send_count(dest+1)+1
        endif
     enddo
   enddo

   nullify(ref)

   call elsi_allocate(e_h,recv_count,e_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   e_h%nnz_l_sp = sum(recv_count,1)
   call MPI_Allreduce(e_h%nnz_l_sp,e_h%nnz_g,1,mpi_integer4,mpi_sum,&
           e_h%mpi_comm,mpierr)

   ! Set send and receive displacement
   call elsi_allocate(e_h,send_displ,e_h%n_procs,"send_displ",caller)
   call elsi_allocate(e_h,recv_displ,e_h%n_procs,"recv_displ",caller)

   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(e_h,row_recv_buf,e_h%nnz_l_sp,"row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,mpi_integer4,&
           row_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,row_send_buf,"row_send_buf")

   ! Column id
   call elsi_allocate(e_h,col_recv_buf,e_h%nnz_l_sp,"col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,mpi_integer4,&
           col_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,col_send_buf,"col_send_buf")

   ! Hamiltonian value
   if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa) then
      call elsi_allocate(e_h,e_h%ham_real_sips,e_h%nnz_l_sp,"ham_real_sips",&
              caller)
   endif

   call MPI_Alltoallv(h_val_send_buf,send_count,send_displ,mpi_real8,&
           e_h%ham_real_sips,recv_count,recv_displ,mpi_real8,e_h%mpi_comm,&
           mpierr)

   call elsi_deallocate(e_h,h_val_send_buf,"h_val_send_buf")

   ! Overlap value
   if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa) then
      if(.not. e_h%ovlp_is_unit) then
         call elsi_allocate(e_h,e_h%ovlp_real_sips,e_h%nnz_l_sp,&
                 "ovlp_real_sips",caller)

         call MPI_Alltoallv(s_val_send_buf,send_count,send_displ,mpi_real8,&
                 e_h%ovlp_real_sips,recv_count,recv_displ,mpi_real8,&
                 e_h%mpi_comm,mpierr)

         call elsi_deallocate(e_h,s_val_send_buf,"s_val_send_buf")
      else
         call elsi_allocate(e_h,e_h%ovlp_real_sips,1,"dummy",caller)
      endif
   endif

   call elsi_deallocate(e_h,send_count,"send_count")
   call elsi_deallocate(e_h,recv_count,"recv_count")
   call elsi_deallocate(e_h,send_displ,"send_displ")
   call elsi_deallocate(e_h,recv_displ,"recv_displ")

   call elsi_allocate(e_h,global_id,e_h%nnz_l_sp,"global_id",caller)

   ! Compute global 1D id
   do i_val = 1,e_h%nnz_l_sp
      global_id(i_val) = int(col_recv_buf(i_val)-1,kind=i8)*&
                            int(e_h%n_basis,kind=i8)+&
                            int(row_recv_buf(i_val),kind=i8)
   enddo

   ! Sort
   if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa .and. .not. e_h%ovlp_is_unit) then
      call elsi_heapsort(e_h%nnz_l_sp,global_id,e_h%ham_real_sips,&
              e_h%ovlp_real_sips,row_recv_buf,col_recv_buf)
   else
      call elsi_heapsort(e_h%nnz_l_sp,global_id,e_h%ham_real_sips,row_recv_buf,&
              col_recv_buf)
   endif

   call elsi_deallocate(e_h,global_id,"global_id")

   ! Compute row index and column pointer
   if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa) then
      call elsi_allocate(e_h,e_h%row_ind_sips,e_h%nnz_l_sp,"row_ind_sips",&
              caller)

      call elsi_allocate(e_h,e_h%col_ptr_sips,(e_h%n_lcol_sp+1),&
              "col_ptr_sips",caller)

      i_col = col_recv_buf(1)-1
      do i_val = 1,e_h%nnz_l_sp
         e_h%row_ind_sips(i_val) = row_recv_buf(i_val)

         if(col_recv_buf(i_val) > i_col) then
            i_col = i_col+1
            e_h%col_ptr_sips(i_col-col_recv_buf(1)+1) = i_val
         endif
      enddo

      e_h%col_ptr_sips(e_h%n_lcol_sp+1) = e_h%nnz_l_sp+1
   endif

   call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 2D block-
!! cyclic dense format to 1D block sparse CCS format, which can be used as input
!! by SIPs.
!!
subroutine elsi_blacs_to_sips_hs_complex(e_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout)      :: e_h                         !< Handle
   complex(kind=r8),  intent(in), target :: h_in(e_h%n_lrow,e_h%n_lcol) !< Hamiltonian
   complex(kind=r8),  intent(in), target :: s_in(e_h%n_lrow,e_h%n_lcol) !< Overlap

   complex(kind=r8), pointer :: ref(:,:) ! Sparsity pattern

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: g_col_id
   integer(kind=i4) :: g_row_id
   integer(kind=i4) :: dest ! Destination of an element
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

   call elsi_get_time(e_h,t0)

   if(e_h%sips_n_elpa == UNSET) then
      e_h%sips_n_elpa = 0
   endif

   if(.not. e_h%ovlp_is_unit) then
      ref => s_in
   else
      ref => h_in
   endif

   if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa) then
      call elsi_get_local_nnz(e_h,ref,e_h%n_lrow,e_h%n_lcol,e_h%nnz_l)

      if(.not. e_h%ovlp_is_unit) then
         call elsi_allocate(e_h,s_val_send_buf,e_h%nnz_l,"s_val_send_buf",&
                 caller)
      endif
   endif

   call elsi_allocate(e_h,row_send_buf,e_h%nnz_l,"row_send_buf",caller)
   call elsi_allocate(e_h,col_send_buf,e_h%nnz_l,"col_send_buf",caller)
   call elsi_allocate(e_h,h_val_send_buf,e_h%nnz_l,"h_val_send_buf",caller)
   call elsi_allocate(e_h,send_count,e_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0
   do i_col = 1,e_h%n_lcol
      call elsi_get_global_col(e_h,g_col_id,i_col)

      ! Compute destination
      dest = (g_col_id-1)/(e_h%n_basis/e_h%n_procs)
      ! The last process may take more
      dest = min(dest,e_h%n_procs-1)

      do i_row = 1,e_h%n_lrow
         if(abs(ref(i_row,i_col)) > e_h%zero_def) then
            i_val = i_val+1

            call elsi_get_global_row(e_h,g_row_id,i_row)

            ! Pack global id and data into bufs
            row_send_buf(i_val)   = g_row_id
            col_send_buf(i_val)   = g_col_id
            h_val_send_buf(i_val) = h_in(i_row,i_col)
            if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa &
               .and. .not. e_h%ovlp_is_unit) then
               s_val_send_buf(i_val) = s_in(i_row,i_col)
            endif

            ! Set send_count
            send_count(dest+1) = send_count(dest+1)+1
        endif
     enddo
   enddo

   nullify(ref)

   call elsi_allocate(e_h,recv_count,e_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   e_h%nnz_l_sp = sum(recv_count,1)
   call MPI_Allreduce(e_h%nnz_l_sp,e_h%nnz_g,1,mpi_integer4,mpi_sum,&
           e_h%mpi_comm,mpierr)

   ! Set send and receive displacement
   call elsi_allocate(e_h,send_displ,e_h%n_procs,"send_displ",caller)
   call elsi_allocate(e_h,recv_displ,e_h%n_procs,"recv_displ",caller)

   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(e_h,row_recv_buf,e_h%nnz_l_sp,"row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,mpi_integer4,&
           row_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,row_send_buf,"row_send_buf")

   ! Column id
   call elsi_allocate(e_h,col_recv_buf,e_h%nnz_l_sp,"col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,mpi_integer4,&
           col_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,col_send_buf,"col_send_buf")

   ! Hamiltonian value
   if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa) then
      call elsi_allocate(e_h,e_h%ham_cmplx_sips,e_h%nnz_l_sp,"ham_cmplx_sips",&
              caller)
   endif

   call MPI_Alltoallv(h_val_send_buf,send_count,send_displ,mpi_complex16,&
           e_h%ham_cmplx_sips,recv_count,recv_displ,mpi_complex16,e_h%mpi_comm,&
           mpierr)

   call elsi_deallocate(e_h,h_val_send_buf,"h_val_send_buf")

   ! Overlap value
   if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa) then
      if(.not. e_h%ovlp_is_unit) then
         call elsi_allocate(e_h,e_h%ovlp_cmplx_sips,e_h%nnz_l_sp,&
                 "ovlp_cmplx_sips",caller)

         call MPI_Alltoallv(s_val_send_buf,send_count,send_displ,mpi_complex16,&
                 e_h%ovlp_cmplx_sips,recv_count,recv_displ,mpi_complex16,&
                 e_h%mpi_comm,mpierr)

         call elsi_deallocate(e_h,s_val_send_buf,"s_val_send_buf")
      else
         call elsi_allocate(e_h,e_h%ovlp_cmplx_sips,1,"dummy",caller)
      endif
   endif

   call elsi_deallocate(e_h,send_count,"send_count")
   call elsi_deallocate(e_h,recv_count,"recv_count")
   call elsi_deallocate(e_h,send_displ,"send_displ")
   call elsi_deallocate(e_h,recv_displ,"recv_displ")

   call elsi_allocate(e_h,global_id,e_h%nnz_l_sp,"global_id",caller)

   ! Compute global 1D id
   do i_val = 1,e_h%nnz_l_sp
      global_id(i_val) = int(col_recv_buf(i_val)-1,kind=i8)*&
                            int(e_h%n_basis,kind=i8)+&
                            int(row_recv_buf(i_val),kind=i8)
   enddo

   ! Sort
   if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa .and. .not. e_h%ovlp_is_unit) then
      call elsi_heapsort(e_h%nnz_l_sp,global_id,e_h%ham_cmplx_sips,&
              e_h%ovlp_cmplx_sips,row_recv_buf,col_recv_buf)
   else
      call elsi_heapsort(e_h%nnz_l_sp,global_id,e_h%ham_cmplx_sips,&
              row_recv_buf,col_recv_buf)
   endif

   call elsi_deallocate(e_h,global_id,"global_id")

   ! Compute row index and column pointer
   if(e_h%n_elsi_calls == 1+e_h%sips_n_elpa) then
      call elsi_allocate(e_h,e_h%row_ind_sips,e_h%nnz_l_sp,"row_ind_sips",&
              caller)

      call elsi_allocate(e_h,e_h%col_ptr_sips,(e_h%n_lcol_sp+1),&
              "col_ptr_sips",caller)

      i_col = col_recv_buf(1)-1
      do i_val = 1,e_h%nnz_l_sp
         e_h%row_ind_sips(i_val) = row_recv_buf(i_val)

         if(col_recv_buf(i_val) > i_col) then
            i_col = i_col+1
            e_h%col_ptr_sips(i_col-col_recv_buf(1)+1) = i_val
         endif
      enddo

      e_h%col_ptr_sips(e_h%n_lcol_sp+1) = e_h%nnz_l_sp+1
   endif

   call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 1D block
!! sparse CCS format to 2D block-cyclic dense format, which can be used as input
!! by ELPA.
!!
subroutine elsi_sips_to_blacs_hs_real(e_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                !< Handle
   real(kind=r8),     intent(in)    :: h_in(e_h%nnz_l_sp) !< Hamiltonian
   real(kind=r8),     intent(in)    :: s_in(e_h%nnz_l_sp) !< Overlap

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: l_row_id ! Local row id in 1D block distribution
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
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character*40, parameter :: caller = "elsi_sips_to_blacs_hs_real"

   call elsi_get_time(e_h,t0)

   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_allocate(e_h,s_val_send_buf,e_h%nnz_l_sp,"s_val_send_buf",&
              caller)
   endif

   call elsi_allocate(e_h,h_val_send_buf,e_h%nnz_l_sp,"h_val_send_buf",caller)
   call elsi_allocate(e_h,row_send_buf,e_h%nnz_l_sp,"row_send_buf",caller)
   call elsi_allocate(e_h,col_send_buf,e_h%nnz_l_sp,"col_send_buf",caller)
   call elsi_allocate(e_h,send_count,e_h%n_procs,"send_count",caller)
   call elsi_allocate(e_h,dest,e_h%nnz_l_sp,"dest",caller)

   i_col = 0
   ! Compute destination and global id
   do i_val = 1,e_h%nnz_l_sp
      if(i_val == e_h%col_ptr_ccs(i_col+1) .and. i_col /= e_h%n_lcol_sp) then
         i_col = i_col+1
      endif
      i_row = e_h%row_ind_ccs(i_val)

      ! Compute global id
      row_send_buf(i_val)   = i_row
      col_send_buf(i_val)   = i_col+e_h%myid*(e_h%n_basis/e_h%n_procs)
      h_val_send_buf(i_val) = h_in(i_val)
      if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
         s_val_send_buf(i_val) = s_in(i_val)
      endif

      ! Compute destination
      proc_row_id = mod((row_send_buf(i_val)-1)/e_h%blk_row,e_h%n_prow)
      proc_col_id = mod((col_send_buf(i_val)-1)/e_h%blk_col,e_h%n_pcol)
      dest(i_val) = proc_col_id+proc_row_id*e_h%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   enddo

   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_heapsort(e_h%nnz_l_sp,dest,h_val_send_buf,s_val_send_buf,&
              row_send_buf,col_send_buf)
   else
      call elsi_heapsort(e_h%nnz_l_sp,dest,h_val_send_buf,row_send_buf,&
              col_send_buf)
   endif

   call elsi_deallocate(e_h,dest,"dest")

   call elsi_allocate(e_h,recv_count,e_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   e_h%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(e_h,send_displ,e_h%n_procs,"send_displ",caller)
   call elsi_allocate(e_h,recv_displ,e_h%n_procs,"recv_displ",caller)

   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row index
   call elsi_allocate(e_h,row_recv_buf,e_h%nnz_l,"row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,mpi_integer4,&
           row_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,row_send_buf,"row_send_buf")

   ! Column index
   call elsi_allocate(e_h,col_recv_buf,e_h%nnz_l,"col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,mpi_integer4,&
           col_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,col_send_buf,"col_send_buf")

   ! Hamiltonian Value
   call elsi_allocate(e_h,h_val_recv_buf,e_h%nnz_l,"h_val_recv_buf",caller)

   call MPI_Alltoallv(h_val_send_buf,send_count,send_displ,mpi_real8,&
           h_val_recv_buf,recv_count,recv_displ,mpi_real8,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,h_val_send_buf,"h_val_send_buf")

   ! Overlap value
   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_allocate(e_h,s_val_recv_buf,e_h%nnz_l,"s_val_recv_buf",caller)

      call MPI_Alltoallv(s_val_send_buf,send_count,send_displ,mpi_real8,&
              s_val_recv_buf,recv_count,recv_displ,mpi_real8,e_h%mpi_comm,&
              mpierr)
   endif

   call elsi_deallocate(e_h,send_count,"send_count")
   call elsi_deallocate(e_h,recv_count,"recv_count")
   call elsi_deallocate(e_h,send_displ,"send_displ")
   call elsi_deallocate(e_h,recv_displ,"recv_displ")

   ! Unpack matrix
   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_allocate(e_h,e_h%ham_real_elpa,e_h%n_lrow,e_h%n_lcol,&
              "ham_real_elpa",caller)

      call elsi_allocate(e_h,e_h%ovlp_real_elpa,e_h%n_lrow,e_h%n_lcol,&
              "ovlp_real_elpa",caller)

      do i_val = 1,e_h%nnz_l
         ! Compute local 2d id
         l_row_id = (row_recv_buf(i_val)-1)/&
                       (e_h%n_prow*e_h%blk_row)*e_h%blk_row+&
                       mod((row_recv_buf(i_val)-1),e_h%blk_row)+1
         l_col_id = (col_recv_buf(i_val)-1)/&
                       (e_h%n_pcol*e_h%blk_col)*e_h%blk_col+&
                       mod((col_recv_buf(i_val)-1),e_h%blk_col)+1

         ! Put value to correct position
         e_h%ham_real_elpa(l_row_id,l_col_id)  = h_val_recv_buf(i_val)
         e_h%ovlp_real_elpa(l_row_id,l_col_id) = s_val_recv_buf(i_val)
      enddo

      call elsi_deallocate(e_h,s_val_recv_buf,"s_val_recv_buf")
   else
      if(e_h%n_elsi_calls == 1) then
         call elsi_allocate(e_h,e_h%ham_real_elpa,e_h%n_lrow,e_h%n_lcol,&
                 "ham_real_elpa",caller)

         call elsi_allocate(e_h,e_h%ovlp_real_elpa,1,1,"dummy",caller)
      endif

      e_h%ham_real_elpa = 0.0_r8

      do i_val = 1,e_h%nnz_l
         ! Compute local 2d id
         l_row_id = (row_recv_buf(i_val)-1)/&
                       (e_h%n_prow*e_h%blk_row)*e_h%blk_row+&
                       mod((row_recv_buf(i_val)-1),e_h%blk_row)+1
         l_col_id = (col_recv_buf(i_val)-1)/&
                       (e_h%n_pcol*e_h%blk_col)*e_h%blk_col+&
                       mod((col_recv_buf(i_val)-1),e_h%blk_col)+1

         ! Put value to correct position
         e_h%ham_real_elpa(l_row_id,l_col_id) = h_val_recv_buf(i_val)
      enddo
   endif

   call elsi_deallocate(e_h,h_val_recv_buf,"h_val_recv_buf")
   call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 1D block
!! sparse CCS format to 2D block-cyclic dense format, which can be used as input
!! by ELPA.
!!
subroutine elsi_sips_to_blacs_hs_complex(e_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                !< Handle
   complex(kind=r8),  intent(in)    :: h_in(e_h%nnz_l_sp) !< Hamiltonian
   complex(kind=r8),  intent(in)    :: s_in(e_h%nnz_l_sp) !< Overlap

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: l_row_id ! Local row id in 1D block distribution
   integer(kind=i4) :: proc_col_id ! Column id in process grid
   integer(kind=i4) :: proc_row_id ! Row id in process grid
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
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character*40, parameter :: caller = "elsi_sips_to_blacs_hs_complex"

   call elsi_get_time(e_h,t0)

   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_allocate(e_h,s_val_send_buf,e_h%nnz_l_sp,"s_val_send_buf",&
              caller)
   endif

   call elsi_allocate(e_h,h_val_send_buf,e_h%nnz_l_sp,"h_val_send_buf",caller)
   call elsi_allocate(e_h,row_send_buf,e_h%nnz_l_sp,"row_send_buf",caller)
   call elsi_allocate(e_h,col_send_buf,e_h%nnz_l_sp,"col_send_buf",caller)
   call elsi_allocate(e_h,send_count,e_h%n_procs,"send_count",caller)
   call elsi_allocate(e_h,dest,e_h%nnz_l_sp,"dest",caller)

   i_col = 0
   ! Compute destination and global id
   do i_val = 1,e_h%nnz_l_sp
      if(i_val == e_h%col_ptr_ccs(i_col+1) .and. i_col /= e_h%n_lcol_sp) then
         i_col = i_col+1
      endif
      i_row = e_h%row_ind_ccs(i_val)

      ! Compute global id
      row_send_buf(i_val)   = i_row
      col_send_buf(i_val)   = i_col+e_h%myid*(e_h%n_basis/e_h%n_procs)
      h_val_send_buf(i_val) = h_in(i_val)
      if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
         s_val_send_buf(i_val) = s_in(i_val)
      endif

      ! Compute destination
      proc_row_id = mod((row_send_buf(i_val)-1)/e_h%blk_row,e_h%n_prow)
      proc_col_id = mod((col_send_buf(i_val)-1)/e_h%blk_col,e_h%n_pcol)
      dest(i_val) = proc_col_id+proc_row_id*e_h%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   enddo

   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_heapsort(e_h%nnz_l_sp,dest,h_val_send_buf,s_val_send_buf,&
              row_send_buf,col_send_buf)
   else
      call elsi_heapsort(e_h%nnz_l_sp,dest,h_val_send_buf,row_send_buf,&
              col_send_buf)
   endif

   call elsi_deallocate(e_h,dest,"dest")

   call elsi_allocate(e_h,recv_count,e_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   e_h%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(e_h,send_displ,e_h%n_procs,"send_displ",caller)
   call elsi_allocate(e_h,recv_displ,e_h%n_procs,"recv_displ",caller)

   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row index
   call elsi_allocate(e_h,row_recv_buf,e_h%nnz_l,"row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,mpi_integer4,&
           row_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,row_send_buf,"row_send_buf")

   ! Column index
   call elsi_allocate(e_h,col_recv_buf,e_h%nnz_l,"col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,mpi_integer4,&
           col_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,col_send_buf,"col_send_buf")

   ! Hamiltonian Value
   call elsi_allocate(e_h,h_val_recv_buf,e_h%nnz_l,"h_val_recv_buf",caller)

   call MPI_Alltoallv(h_val_send_buf,send_count,send_displ,mpi_complex16,&
           h_val_recv_buf,recv_count,recv_displ,mpi_complex16,e_h%mpi_comm,&
           mpierr)

   call elsi_deallocate(e_h,h_val_send_buf,"h_val_send_buf")

   ! Overlap value
   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_allocate(e_h,s_val_recv_buf,e_h%nnz_l,"s_val_recv_buf",caller)

      call MPI_Alltoallv(s_val_send_buf,send_count,send_displ,mpi_complex16,&
              s_val_recv_buf,recv_count,recv_displ,mpi_complex16,e_h%mpi_comm,&
              mpierr)
   endif

   call elsi_deallocate(e_h,send_count,"send_count")
   call elsi_deallocate(e_h,recv_count,"recv_count")
   call elsi_deallocate(e_h,send_displ,"send_displ")
   call elsi_deallocate(e_h,recv_displ,"recv_displ")

   ! Unpack matrix
   if(e_h%n_elsi_calls == 1 .and. .not. e_h%ovlp_is_unit) then
      call elsi_allocate(e_h,e_h%ham_cmplx_elpa,e_h%n_lrow,e_h%n_lcol,&
              "ham_cmplx_elpa",caller)

      call elsi_allocate(e_h,e_h%ovlp_cmplx_elpa,e_h%n_lrow,e_h%n_lcol,&
              "ovlp_cmplx_elpa",caller)

      do i_val = 1,e_h%nnz_l
         ! Compute local 2d id
         l_row_id = (row_recv_buf(i_val)-1)/&
                       (e_h%n_prow*e_h%blk_row)*e_h%blk_row+&
                       mod((row_recv_buf(i_val)-1),e_h%blk_row)+1
         l_col_id = (col_recv_buf(i_val)-1)/&
                       (e_h%n_pcol*e_h%blk_col)*e_h%blk_col+&
                       mod((col_recv_buf(i_val)-1),e_h%blk_col)+1

         ! Put value to correct position
         e_h%ham_cmplx_elpa(l_row_id,l_col_id)  = h_val_recv_buf(i_val)
         e_h%ovlp_cmplx_elpa(l_row_id,l_col_id) = s_val_recv_buf(i_val)
      enddo

      call elsi_deallocate(e_h,s_val_recv_buf,"s_val_recv_buf")
   else
      if(e_h%n_elsi_calls == 1) then
         call elsi_allocate(e_h,e_h%ham_cmplx_elpa,e_h%n_lrow,e_h%n_lcol,&
                 "ham_cmplx_elpa",caller)

         call elsi_allocate(e_h,e_h%ovlp_cmplx_elpa,1,1,"dummy",caller)
      endif

      e_h%ham_cmplx_elpa = (0.0_r8,0.0_r8)

      do i_val = 1,e_h%nnz_l
         ! Compute local 2d id
         l_row_id = (row_recv_buf(i_val)-1)/&
                       (e_h%n_prow*e_h%blk_row)*e_h%blk_row+&
                       mod((row_recv_buf(i_val)-1),e_h%blk_row)+1
         l_col_id = (col_recv_buf(i_val)-1)/&
                       (e_h%n_pcol*e_h%blk_col)*e_h%blk_col+&
                       mod((col_recv_buf(i_val)-1),e_h%blk_col)+1

         ! Put value to correct position
         e_h%ham_cmplx_elpa(l_row_id,l_col_id) = h_val_recv_buf(i_val)
      enddo
   endif

   call elsi_deallocate(e_h,h_val_recv_buf,"h_val_recv_buf")
   call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine converts density matrix in 2D block-cyclic dense format to 1D
!! block sparse CCS format.
!!
subroutine elsi_blacs_to_sips_dm_real(e_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                 !< Handle
   real(kind=r8),     intent(out)   :: d_out(e_h%nnz_l_sp) !< Density matrix

   real(kind=r8), pointer :: ref(:,:) ! Sparsity pattern

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: l_row_id ! Local row id in 1D block distribution
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: nnz_l_aux
   integer(kind=i4) :: n_lcol_aux
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

   call elsi_get_time(e_h,t0)

   if(e_h%solver == ELPA_SOLVER .or. e_h%solver == DMP_SOLVER) then
      ref => e_h%dm_real
   elseif(e_h%solver == OMM_SOLVER) then
      ref => e_h%dm_omm%dval
   endif

   call elsi_get_local_nnz(e_h,ref,e_h%n_lrow,e_h%n_lcol,e_h%nnz_l)

   call elsi_allocate(e_h,val_send_buf,e_h%nnz_l,"val_send_buf",caller)
   call elsi_allocate(e_h,row_send_buf,e_h%nnz_l,"row_send_buf",caller)
   call elsi_allocate(e_h,col_send_buf,e_h%nnz_l,"col_send_buf",caller)
   call elsi_allocate(e_h,send_count,e_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0
   do i_col = 1,e_h%n_lcol
      do i_row = 1,e_h%n_lrow
         if(abs(ref(i_row,i_col)) > e_h%zero_def) then
            i_val = i_val+1

            call elsi_get_global_row(e_h,row_send_buf(i_val),i_row)
            call elsi_get_global_col(e_h,col_send_buf(i_val),i_col)

            ! Compute destination
            dest = (col_send_buf(i_val)-1)/(e_h%n_basis/e_h%n_procs)
            ! The last process may take more
            dest = min(dest,e_h%n_procs-1)

            ! Pack data
            val_send_buf(i_val) = ref(i_row,i_col)

            ! Set send_count
            send_count(dest+1) = send_count(dest+1)+1
         endif
      enddo
   enddo

   nullify(ref)

   call elsi_allocate(e_h,recv_count,e_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(e_h,send_displ,e_h%n_procs,"send_displ",caller)
   call elsi_allocate(e_h,recv_displ,e_h%n_procs,"recv_displ",caller)

   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(e_h,row_recv_buf,nnz_l_aux,"row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,mpi_integer4,&
           row_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,row_send_buf,"row_send_buf")

   ! Column id
   call elsi_allocate(e_h,col_recv_buf,nnz_l_aux,"col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,mpi_integer4,&
           col_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,col_send_buf,"col_send_buf")

   ! Density matrix value
   call elsi_allocate(e_h,val_recv_buf,nnz_l_aux,"val_recv_buf",caller)

   call MPI_Alltoallv(val_send_buf,send_count,send_displ,mpi_real8,&
           val_recv_buf,recv_count,recv_displ,mpi_real8,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,val_send_buf,"val_send_buf")

   call elsi_deallocate(e_h,send_count,"send_count")
   call elsi_deallocate(e_h,recv_count,"recv_count")
   call elsi_deallocate(e_h,send_displ,"send_displ")
   call elsi_deallocate(e_h,recv_displ,"recv_displ")

   d_out = 0.0_r8

   if(e_h%myid == e_h%n_procs-1) then
      n_lcol_aux = e_h%n_basis-e_h%n_lcol_sp
      n_lcol_aux = n_lcol_aux/(e_h%n_procs-1)
   else
      n_lcol_aux = e_h%n_lcol_sp
   endif

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      l_col_id = col_recv_buf(i_val)-e_h%myid*n_lcol_aux
      l_row_id = row_recv_buf(i_val)

      do j_val = e_h%col_ptr_ccs(l_col_id),e_h%col_ptr_ccs(l_col_id+1)-1
         if(e_h%row_ind_ccs(j_val) == l_row_id) then
            d_out(j_val) = val_recv_buf(i_val)
         endif
      enddo
   enddo

   call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")
   call elsi_deallocate(e_h,val_recv_buf,"val_recv_buf")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine converts density matrix in 2D block-cyclic dense format to 1D
!! block sparse CCS format.
!!
subroutine elsi_blacs_to_sips_dm_complex(e_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                 !< Handle
   complex(kind=r8),  intent(out)   :: d_out(e_h%nnz_l_sp) !< Density matrix

   complex(kind=r8), pointer :: ref(:,:) ! Sparsity pattern

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: l_row_id ! Local row id in 1D block distribution
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: nnz_l_aux
   integer(kind=i4) :: n_lcol_aux
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

   character*40, parameter :: caller = "elsi_blacs_to_sips_dm_complex"

   call elsi_get_time(e_h,t0)

   if(e_h%solver == ELPA_SOLVER .or. e_h%solver == DMP_SOLVER) then
      ref => e_h%dm_cmplx
   elseif(e_h%solver == OMM_SOLVER) then
      ref => e_h%dm_omm%zval
   endif

   call elsi_get_local_nnz(e_h,ref,e_h%n_lrow,e_h%n_lcol,e_h%nnz_l)

   call elsi_allocate(e_h,val_send_buf,e_h%nnz_l,"val_send_buf",caller)
   call elsi_allocate(e_h,row_send_buf,e_h%nnz_l,"row_send_buf",caller)
   call elsi_allocate(e_h,col_send_buf,e_h%nnz_l,"col_send_buf",caller)
   call elsi_allocate(e_h,send_count,e_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0
   do i_col = 1,e_h%n_lcol
      do i_row = 1,e_h%n_lrow
         if(abs(ref(i_row,i_col)) > e_h%zero_def) then
            i_val = i_val+1

            call elsi_get_global_row(e_h,row_send_buf(i_val),i_row)
            call elsi_get_global_col(e_h,col_send_buf(i_val),i_col)

            ! Compute destination
            dest = (col_send_buf(i_val)-1)/(e_h%n_basis/e_h%n_procs)
            ! The last process may take more
            dest = min(dest,e_h%n_procs-1)

            ! Pack data
            val_send_buf(i_val) = ref(i_row,i_col)

            ! Set send_count
            send_count(dest+1) = send_count(dest+1)+1
         endif
      enddo
   enddo

   nullify(ref)

   call elsi_allocate(e_h,recv_count,e_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(e_h,send_displ,e_h%n_procs,"send_displ",caller)
   call elsi_allocate(e_h,recv_displ,e_h%n_procs,"recv_displ",caller)

   do i_proc = 2,e_h%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(e_h,row_recv_buf,nnz_l_aux,"row_recv_buf",caller)

   call MPI_Alltoallv(row_send_buf,send_count,send_displ,mpi_integer4,&
           row_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,row_send_buf,"row_send_buf")

   ! Column id
   call elsi_allocate(e_h,col_recv_buf,nnz_l_aux,"col_recv_buf",caller)

   call MPI_Alltoallv(col_send_buf,send_count,send_displ,mpi_integer4,&
           col_recv_buf,recv_count,recv_displ,mpi_integer4,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,col_send_buf,"col_send_buf")

   ! Density matrix value
   call elsi_allocate(e_h,val_recv_buf,nnz_l_aux,"val_recv_buf",caller)

   call MPI_Alltoallv(val_send_buf,send_count,send_displ,mpi_complex16,&
           val_recv_buf,recv_count,recv_displ,mpi_complex16,e_h%mpi_comm,mpierr)

   call elsi_deallocate(e_h,val_send_buf,"val_send_buf")

   call elsi_deallocate(e_h,send_count,"send_count")
   call elsi_deallocate(e_h,recv_count,"recv_count")
   call elsi_deallocate(e_h,send_displ,"send_displ")
   call elsi_deallocate(e_h,recv_displ,"recv_displ")

   d_out = 0.0_r8

   if(e_h%myid == e_h%n_procs-1) then
      n_lcol_aux = e_h%n_basis-e_h%n_lcol_sp
      n_lcol_aux = n_lcol_aux/(e_h%n_procs-1)
   else
      n_lcol_aux = e_h%n_lcol_sp
   endif

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      l_col_id = col_recv_buf(i_val)-e_h%myid*n_lcol_aux
      l_row_id = row_recv_buf(i_val)

      do j_val = e_h%col_ptr_ccs(l_col_id),e_h%col_ptr_ccs(l_col_id+1)-1
         if(e_h%row_ind_ccs(j_val) == l_row_id) then
            d_out(j_val) = val_recv_buf(i_val)
         endif
      enddo
   enddo

   call elsi_deallocate(e_h,row_recv_buf,"row_recv_buf")
   call elsi_deallocate(e_h,col_recv_buf,"col_recv_buf")
   call elsi_deallocate(e_h,val_recv_buf,"val_recv_buf")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine converts matrix format and distribution from BLACS to CheSS.
!!
subroutine elsi_blacs_to_chess_hs_real(e_h,h_in,s_in)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                         !< Handle
   real(kind=r8),     intent(in)    :: h_in(e_h%n_lrow,e_h%n_lcol) !< Hamiltonian
   real(kind=r8),     intent(in)    :: s_in(e_h%n_lrow,e_h%n_lcol) !< Overlap

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character*200 :: info_str

   character*40, parameter :: caller = "elsi_blacs_to_chess_hs_real"

   call elsi_get_time(e_h,t0)

   ! First convert to SIPs 1D block distribution
   e_h%n_lcol_sp = e_h%n_basis/e_h%n_procs

   ! The last process holds all remaining columns
   if(e_h%myid == e_h%n_procs-1) then
      e_h%n_lcol_sp = e_h%n_basis-(e_h%n_procs-1)*e_h%n_lcol_sp
   endif

   call elsi_blacs_to_sips_hs_real(e_h,h_in,s_in)

   ! Then get the global matrices
   call elsi_sips_to_chess_hs(e_h)

   e_h%nnz_l_sp  = e_h%nnz_g
   e_h%n_lcol_sp = e_h%n_basis

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine gets global matrices from 1D block sparse CCS format.
!!
subroutine elsi_sips_to_chess_hs(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)

   integer(kind=i4) :: prev_nnz
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_sips_to_chess_hs"

   ! Set recv_count and recv_displ
   call elsi_allocate(e_h,recv_count,e_h%n_procs,"recv_count",caller)
   call elsi_allocate(e_h,recv_displ,e_h%n_procs,"recv_displ",caller)

   call MPI_Allgather(e_h%nnz_l_sp,1,mpi_integer4,recv_count,1,mpi_integer4,&
           e_h%mpi_comm,mpierr)

   do i_proc = 2,e_h%n_procs
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Get the global matrices
   if(e_h%n_elsi_calls == 1) then
      if(.not. e_h%ovlp_is_unit) then
         ! Overlap value
         call elsi_allocate(e_h,e_h%ovlp_real_chess,e_h%nnz_g,&
                 "ovlp_real_chess",caller)

         call MPI_Allgatherv(e_h%ovlp_real_sips,e_h%nnz_l_sp,mpi_real8,&
                 e_h%ovlp_real_chess,recv_count,recv_displ,mpi_real8,&
                 e_h%mpi_comm,mpierr)

         call elsi_deallocate(e_h,e_h%ovlp_real_sips,"ovlp_real_sips")
      else
         call elsi_allocate(e_h,e_h%ovlp_real_chess,1,"dummy",caller)
      endif

      ! Row index
      call elsi_allocate(e_h,e_h%row_ind_chess,e_h%nnz_g,"row_ind_chess",caller)

      call MPI_Allgatherv(e_h%row_ind_sips,e_h%nnz_l_sp,mpi_integer4,&
              e_h%row_ind_chess,recv_count,recv_displ,mpi_integer4,&
              e_h%mpi_comm,mpierr)

      call elsi_deallocate(e_h,e_h%row_ind_sips,"row_ind_sips")

      call elsi_allocate(e_h,e_h%ham_real_chess,e_h%nnz_g,"ham_real_chess",&
              caller)
   endif

   ! Hamiltonian value
   call MPI_Allgatherv(e_h%ham_real_sips,e_h%nnz_l_sp,mpi_real8,&
           e_h%ham_real_chess,recv_count,recv_displ,mpi_real8,e_h%mpi_comm,&
           mpierr)

   if(e_h%n_elsi_calls == 1) then
      ! Set recv_count and recv_displ
      recv_count(1:e_h%n_procs-1) = e_h%n_basis/e_h%n_procs
      recv_count(e_h%n_procs)     = e_h%n_basis-(e_h%n_procs-1)*recv_count(1)

      do i_proc = 2,e_h%n_procs
         recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
      enddo

      ! Shift column pointers
      prev_nnz = 0

      call MPI_Exscan(e_h%nnz_l_sp,prev_nnz,1,mpi_integer4,mpi_sum,&
              e_h%mpi_comm,mpierr)

      e_h%col_ptr_sips = e_h%col_ptr_sips+prev_nnz

      ! Column pointer
      call elsi_allocate(e_h,e_h%col_ptr_chess,e_h%n_basis+1,"col_ptr_chess",&
              caller)

      call MPI_Allgatherv(e_h%col_ptr_sips,e_h%n_lcol_sp,mpi_integer4,&
              e_h%col_ptr_chess,recv_count,recv_displ,mpi_integer4,&
              e_h%mpi_comm,mpierr)

      e_h%col_ptr_chess(e_h%n_basis+1) = e_h%nnz_g+1

      call elsi_deallocate(e_h,e_h%col_ptr_sips,"col_ptr_sips")
   endif

   call elsi_deallocate(e_h,recv_count,"recv_count")
   call elsi_deallocate(e_h,recv_displ,"recv_displ")

end subroutine

!>
!! This routine converts density matrix computed by CheSS and stored in sparse
!! CCS format to 2D block-cyclic dense format.
!!
subroutine elsi_chess_to_blacs_dm_real(e_h,d_out)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                          !< Handle
   real(kind=r8),     intent(out)   :: d_out(e_h%n_lrow,e_h%n_lcol) !< Density matrix

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: g_row
   integer(kind=i4) :: g_col

   character*40, parameter :: caller = "elsi_chess_to_blacs_dm_real"

   d_out = 0.0_r8

   do i_col = 1,e_h%n_lcol
      call elsi_get_global_col(e_h,g_col,i_col)

      do i_val = e_h%col_ptr_ccs(g_col),e_h%col_ptr_ccs(g_col+1)-1
         g_row = e_h%row_ind_ccs(i_val)

         if(e_h%loc_row(g_row) == 0) cycle

         i_row = (g_row-1)/(e_h%n_prow*e_h%blk_row)*e_h%blk_row+&
                    mod((g_row-1),e_h%blk_row)+1

         d_out(i_row,i_col) = e_h%dm_chess%matrix_compr(i_val)
      enddo
   enddo

end subroutine

end module ELSI_MATCONV
