! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides matrix conversion and redistribution routines.
!!
module ELSI_REDIST

   use ELSI_CONSTANTS, only: UNSET
   use ELSI_DATATYPE,  only: elsi_param_t,elsi_basic_t
   use ELSI_IO,        only: elsi_say,elsi_get_time
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,       only: elsi_check_mpi,mpi_sum,mpi_real8,mpi_complex16,&
                             mpi_integer4,mpi_comm_self
   use ELSI_PRECISION, only: r8,i4,i8
   use ELSI_SORT,      only: elsi_heapsort
   use ELSI_UTILS,     only: elsi_get_nnz,elsi_get_gid,elsi_get_lid

   implicit none

   private

   public :: elsi_blacs_to_pexsi_hs_dim
   public :: elsi_blacs_to_pexsi_hs
   public :: elsi_blacs_to_siesta_dm
   public :: elsi_blacs_to_sips_dm
   public :: elsi_blacs_to_sips_hs_dim
   public :: elsi_blacs_to_sips_hs
   public :: elsi_pexsi_to_blacs_dm
   public :: elsi_pexsi_to_siesta_dm
   public :: elsi_siesta_to_blacs_hs
   public :: elsi_siesta_to_pexsi_hs_dim
   public :: elsi_siesta_to_pexsi_hs
   public :: elsi_siesta_to_sips_hs_dim
   public :: elsi_siesta_to_sips_hs
   public :: elsi_sips_to_blacs_dm
   public :: elsi_sips_to_blacs_ev
   public :: elsi_sips_to_blacs_hs
   public :: elsi_sips_to_siesta_dm

   interface elsi_blacs_to_pexsi_hs_dim
      module procedure elsi_blacs_to_pexsi_hs_dim_real
      module procedure elsi_blacs_to_pexsi_hs_dim_cmplx
   end interface

   interface elsi_blacs_to_pexsi_hs
      module procedure elsi_blacs_to_pexsi_hs_real
      module procedure elsi_blacs_to_pexsi_hs_cmplx
   end interface

   interface elsi_blacs_to_siesta_dm
      module procedure elsi_blacs_to_siesta_dm_real
      module procedure elsi_blacs_to_siesta_dm_cmplx
   end interface

   interface elsi_blacs_to_sips_dm
      module procedure elsi_blacs_to_sips_dm_real
      module procedure elsi_blacs_to_sips_dm_cmplx
   end interface

   interface elsi_blacs_to_sips_hs_dim
      module procedure elsi_blacs_to_sips_hs_dim_real
      module procedure elsi_blacs_to_sips_hs_dim_cmplx
   end interface

   interface elsi_blacs_to_sips_hs
      module procedure elsi_blacs_to_sips_hs_real
      module procedure elsi_blacs_to_sips_hs_cmplx
   end interface

   interface elsi_pexsi_to_blacs_dm
      module procedure elsi_pexsi_to_blacs_dm_real
      module procedure elsi_pexsi_to_blacs_dm_cmplx
   end interface

   interface elsi_pexsi_to_siesta_dm
      module procedure elsi_pexsi_to_siesta_dm_real
      module procedure elsi_pexsi_to_siesta_dm_cmplx
   end interface

   interface elsi_siesta_to_blacs_hs
      module procedure elsi_siesta_to_blacs_hs_real
      module procedure elsi_siesta_to_blacs_hs_cmplx
   end interface

   interface elsi_siesta_to_pexsi_hs
      module procedure elsi_siesta_to_pexsi_hs_real
      module procedure elsi_siesta_to_pexsi_hs_cmplx
   end interface

   interface elsi_siesta_to_sips_hs
      module procedure elsi_siesta_to_sips_hs_real
      module procedure elsi_siesta_to_sips_hs_cmplx
   end interface

   interface elsi_sips_to_blacs_dm
      module procedure elsi_sips_to_blacs_dm_real
      module procedure elsi_sips_to_blacs_dm_cmplx
   end interface

   interface elsi_sips_to_blacs_ev
      module procedure elsi_sips_to_blacs_ev_real
   end interface

   interface elsi_sips_to_blacs_hs
      module procedure elsi_sips_to_blacs_hs_real
      module procedure elsi_sips_to_blacs_hs_cmplx
   end interface

   interface elsi_sips_to_siesta_dm
      module procedure elsi_sips_to_siesta_dm_real
      module procedure elsi_sips_to_siesta_dm_cmplx
   end interface

contains

!>
!! This routine gets local number of nonzero elements in matrices in 1D block
!! CSC format converted from 2D block-cyclic dense format.
!!
subroutine elsi_blacs_to_pexsi_hs_dim_real(ph,bh,ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8),      intent(in)    :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)    :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: myid_in_pole
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=40), parameter :: caller = "elsi_blacs_to_pexsi_hs_dim_real"

   call elsi_allocate(bh,dest,ph%pexsi_np_per_pole,"dest",caller)
   call elsi_allocate(bh,nnz,ph%pexsi_np_per_pole,"nnz",caller)

   if(.not. ph%ovlp_is_unit) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/ph%pexsi_np_per_pole)
         i_proc = min(i_proc,ph%pexsi_np_per_pole-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            endif
         enddo
      enddo
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/ph%pexsi_np_per_pole)
         i_proc = min(i_proc,ph%pexsi_np_per_pole-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            endif
         enddo
      enddo
   endif

   call MPI_Allreduce(dest,nnz,bh%n_procs,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   call MPI_Comm_rank(ph%pexsi_comm_in_pole,myid_in_pole,ierr)

   call elsi_check_mpi(bh,"MPI_Comm_rank",ierr,caller)

   bh%nnz_l_sp  = nnz(myid_in_pole+1)
   bh%nnz_l_sp1 = bh%nnz_l_sp

   call MPI_Allreduce(bh%nnz_l_sp,bh%nnz_g,1,mpi_integer4,mpi_sum,&
           ph%pexsi_comm_among_pole,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! This routine gets local number of nonzero elements in matrices in 1D block
!! CSC format converted from 2D block-cyclic dense format.
!!
subroutine elsi_blacs_to_pexsi_hs_dim_cmplx(ph,bh,ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8),   intent(in)    :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(in)    :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: myid_in_pole
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=40), parameter :: caller = "elsi_blacs_to_pexsi_hs_dim_cmplx"

   call elsi_allocate(bh,dest,ph%pexsi_np_per_pole,"dest",caller)
   call elsi_allocate(bh,nnz,ph%pexsi_np_per_pole,"nnz",caller)

   if(.not. ph%ovlp_is_unit) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/ph%pexsi_np_per_pole)
         i_proc = min(i_proc,ph%pexsi_np_per_pole-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            endif
         enddo
      enddo
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/ph%pexsi_np_per_pole)
         i_proc = min(i_proc,ph%pexsi_np_per_pole-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            endif
         enddo
      enddo
   endif

   call MPI_Allreduce(dest,nnz,bh%n_procs,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   call MPI_Comm_rank(ph%pexsi_comm_in_pole,myid_in_pole,ierr)

   call elsi_check_mpi(bh,"MPI_Comm_rank",ierr,caller)

   bh%nnz_l_sp  = nnz(myid_in_pole+1)
   bh%nnz_l_sp1 = bh%nnz_l_sp

   call MPI_Allreduce(bh%nnz_l_sp,bh%nnz_g,1,mpi_integer4,mpi_sum,&
           ph%pexsi_comm_among_pole,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 2D block-
!! cyclic dense format to 1D block CCS format.
!!
subroutine elsi_blacs_to_pexsi_hs_real(ph,bh,ham_den,ovlp_den,ham_csc,ovlp_csc,&
              row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8),      intent(in)    :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)    :: ovlp_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: ham_csc(bh%nnz_l_sp)
   real(kind=r8),      intent(out)   :: ovlp_csc(bh%nnz_l_sp)
   integer(kind=i4),   intent(out)   :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(out)   :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4)   :: n_para_task
   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: g_col
   integer(kind=i4)   :: g_row
   integer(kind=i4)   :: d1
   integer(kind=i4)   :: d2
   integer(kind=i4)   :: d11
   integer(kind=i4)   :: d12
   integer(kind=i4)   :: d21
   integer(kind=i4)   :: d22
   integer(kind=i4)   :: dest ! Destination of an element
   integer(kind=i4)   :: this_n_cols
   integer(kind=i4)   :: nnz_l_aux
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send(:)
   real(kind=r8),    allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: h_val_recv(:)
   real(kind=r8),    allocatable :: s_val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: locat(:) ! Location of each global column
   integer(kind=i8), allocatable :: global_id(:) ! Global 1D id

   character(len=40), parameter :: caller = "elsi_blacs_to_pexsi_hs_real"

   call elsi_get_time(t0)

   n_para_task = bh%n_procs/ph%pexsi_np_per_pole

   if(ph%n_calls == 1) then
      if(.not. ph%ovlp_is_unit) then
         call elsi_get_nnz(bh%def0,ovlp_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)

         call elsi_allocate(bh,s_val_send,bh%nnz_l,"s_val_send",caller)
      else
         call elsi_get_nnz(bh%def0,ham_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)
      endif
   endif

   call elsi_allocate(bh,locat,ph%n_basis,"locat",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l,"h_val_send",caller)

   ! Compute d1,d2,d11,d12,d21,d22 (need explanation)
   d1  = ph%n_basis/ph%pexsi_np_per_pole
   d2  = ph%n_basis-(ph%pexsi_np_per_pole-1)*d1
   d11 = d1/n_para_task
   d12 = d1-(n_para_task-1)*d11
   d21 = d2/n_para_task
   d22 = d2-(n_para_task-1)*d21

   i_val = 0

   if(d11 > 0) then
      do i_proc = 0,bh%n_procs-n_para_task-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d12
         else
            this_n_cols = d11
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = 0,bh%n_procs-n_para_task-1
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
      do i_proc = bh%n_procs-n_para_task,bh%n_procs-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d22
         else
            this_n_cols = d21
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = bh%n_procs-n_para_task,bh%n_procs-1
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

   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0

   if(.not. ph%ovlp_is_unit) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = min(locat(g_col),bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val)   = g_row
               col_send(i_val)   = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)
               if(ph%n_calls == 1) then
                  s_val_send(i_val) = ovlp_den(i_row,i_col)
               endif

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = min(locat(g_col),bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val)   = g_row
               col_send(i_val)   = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   endif

   call elsi_deallocate(bh,locat,"locat")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(bh,row_recv,nnz_l_aux,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column id
   call elsi_allocate(bh,col_recv,nnz_l_aux,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian value
   call elsi_allocate(bh,h_val_recv,nnz_l_aux,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_real8,h_val_recv,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_recv,nnz_l_aux,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,s_val_recv,&
              recv_count,recv_displ,mpi_real8,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   endif

   call elsi_allocate(bh,global_id,nnz_l_aux,"global_id",caller)

   ! Compute global 1D id
   global_id = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+&
                  int(row_recv,kind=i8)

   ! Sort
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_heapsort(nnz_l_aux,global_id,h_val_recv,s_val_recv,row_recv,&
              col_recv)
   else
      call elsi_heapsort(nnz_l_aux,global_id,h_val_recv,row_recv,col_recv)
   endif

   call elsi_deallocate(bh,global_id,"global_id")

   ! Set send_count, all data sent to 1st pole
   send_count = 0
   send_count(bh%myid/n_para_task+1) = nnz_l_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send and receive displacement
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   if(ph%n_calls == 1) then
      ! Row id
      call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)

      call MPI_Alltoallv(row_recv,send_count,send_displ,mpi_integer4,row_send,&
              recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,row_recv,"row_recv")

      ! Column id
      call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)

      call MPI_Alltoallv(col_recv,send_count,send_displ,mpi_integer4,col_send,&
              recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,col_recv,"col_recv")

      if(.not. ph%ovlp_is_unit) then
         ! Overlap value
         call MPI_Alltoallv(s_val_recv,send_count,send_displ,mpi_real8,&
                 ovlp_csc,recv_count,recv_displ,mpi_real8,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

         call elsi_deallocate(bh,s_val_recv,"s_val_recv")
      endif
   else
      call elsi_deallocate(bh,row_recv,"row_recv")
      call elsi_deallocate(bh,col_recv,"col_recv")
   endif

   ham_csc = 0.0_r8

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv,send_count,send_displ,mpi_real8,ham_csc,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   if(ph%n_calls == 1) then
      ! Only 1st pole computes row index and column pointer
      if(ph%pexsi_my_prow == 0) then
         ! Compute row index and column pointer
         i_col = col_send(1)-1
         do i_val = 1,bh%nnz_l_sp
            row_ind(i_val) = row_send(i_val)

            if(col_send(i_val) > i_col) then
               i_col = i_col+1
               col_ptr(i_col-col_send(1)+1) = i_val
            endif
         enddo

         col_ptr(bh%n_lcol_sp+1) = bh%nnz_l_sp+1
      endif

      call elsi_deallocate(bh,row_send,"row_send")
      call elsi_deallocate(bh,col_send,"col_send")
   endif

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 2D block-
!! cyclic dense format to 1D block CCS format.
!!
subroutine elsi_blacs_to_pexsi_hs_cmplx(ph,bh,ham_den,ovlp_den,ham_csc,&
              ovlp_csc,row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8),   intent(in)    :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(in)    :: ovlp_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(out)   :: ham_csc(bh%nnz_l_sp)
   complex(kind=r8),   intent(out)   :: ovlp_csc(bh%nnz_l_sp)
   integer(kind=i4),   intent(out)   :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(out)   :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4)   :: n_para_task
   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: g_col
   integer(kind=i4)   :: g_row
   integer(kind=i4)   :: d1
   integer(kind=i4)   :: d2
   integer(kind=i4)   :: d11
   integer(kind=i4)   :: d12
   integer(kind=i4)   :: d21
   integer(kind=i4)   :: d22
   integer(kind=i4)   :: dest ! Destination of an element
   integer(kind=i4)   :: this_n_cols
   integer(kind=i4)   :: nnz_l_aux
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: h_val_send(:)
   complex(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   complex(kind=r8), allocatable :: h_val_recv(:)
   complex(kind=r8), allocatable :: s_val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: locat(:) ! Location of each global column
   integer(kind=i8), allocatable :: global_id(:) ! Global 1D id

   character(len=40), parameter :: caller = "elsi_blacs_to_pexsi_hs_cmplx"

   call elsi_get_time(t0)

   n_para_task = bh%n_procs/ph%pexsi_np_per_pole

   if(ph%n_calls == 1) then
      if(.not. ph%ovlp_is_unit) then
         call elsi_get_nnz(bh%def0,ovlp_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)

         call elsi_allocate(bh,s_val_send,bh%nnz_l,"s_val_send",caller)
      else
         call elsi_get_nnz(bh%def0,ham_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)
      endif
   endif

   call elsi_allocate(bh,locat,ph%n_basis,"locat",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l,"h_val_send",caller)

   ! Compute d1,d2,d11,d12,d21,d22 (need explanation)
   d1  = ph%n_basis/ph%pexsi_np_per_pole
   d2  = ph%n_basis-(ph%pexsi_np_per_pole-1)*d1
   d11 = d1/n_para_task
   d12 = d1-(n_para_task-1)*d11
   d21 = d2/n_para_task
   d22 = d2-(n_para_task-1)*d21

   i_val = 0

   if(d11 > 0) then
      do i_proc = 0,bh%n_procs-n_para_task-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d12
         else
            this_n_cols = d11
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = 0,bh%n_procs-n_para_task-1
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
      do i_proc = bh%n_procs-n_para_task,bh%n_procs-1
         if(mod((i_proc+1),n_para_task) == 0) then
            this_n_cols = d22
         else
            this_n_cols = d21
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = bh%n_procs-n_para_task,bh%n_procs-1
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

   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0
   if(.not. ph%ovlp_is_unit) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = min(locat(g_col),bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val)   = g_row
               col_send(i_val)   = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)
               if(ph%n_calls == 1) then
                  s_val_send(i_val) = ovlp_den(i_row,i_col)
               endif

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = min(locat(g_col),bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val)   = g_row
               col_send(i_val)   = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   endif

   call elsi_deallocate(bh,locat,"locat")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(bh,row_recv,nnz_l_aux,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column id
   call elsi_allocate(bh,col_recv,nnz_l_aux,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian value
   call elsi_allocate(bh,h_val_recv,nnz_l_aux,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_complex16,&
           h_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_recv,nnz_l_aux,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
              s_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   endif

   call elsi_allocate(bh,global_id,nnz_l_aux,"global_id",caller)

   ! Compute global 1D id
   global_id = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+&
                  int(row_recv,kind=i8)

   ! Sort
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_heapsort(nnz_l_aux,global_id,h_val_recv,s_val_recv,row_recv,&
              col_recv)
   else
      call elsi_heapsort(nnz_l_aux,global_id,h_val_recv,row_recv,col_recv)
   endif

   call elsi_deallocate(bh,global_id,"global_id")

   ! Set send_count, all data sent to 1st pole
   send_count = 0
   send_count(bh%myid/n_para_task+1) = nnz_l_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send and receive displacement
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   if(ph%n_calls == 1) then
      ! Row id
      call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)

      call MPI_Alltoallv(row_recv,send_count,send_displ,mpi_integer4,row_send,&
              recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,row_recv,"row_recv")

      ! Column id
      call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)

      call MPI_Alltoallv(col_recv,send_count,send_displ,mpi_integer4,col_send,&
              recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,col_recv,"col_recv")

      if(.not. ph%ovlp_is_unit) then
         ! Overlap value
         call MPI_Alltoallv(s_val_recv,send_count,send_displ,mpi_complex16,&
                 ovlp_csc,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

         call elsi_deallocate(bh,s_val_recv,"s_val_recv")
      endif
   else
      call elsi_deallocate(bh,row_recv,"row_recv")
      call elsi_deallocate(bh,col_recv,"col_recv")
   endif

   ham_csc = (0.0_r8,0.0_r8)

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv,send_count,send_displ,mpi_complex16,ham_csc,&
           recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   if(ph%n_calls == 1) then
      ! Only 1st pole computes row index and column pointer
      if(ph%pexsi_my_prow == 0) then
         ! Compute row index and column pointer
         i_col = col_send(1)-1
         do i_val = 1,bh%nnz_l_sp
            row_ind(i_val) = row_send(i_val)

            if(col_send(i_val) > i_col) then
               i_col = i_col+1
               col_ptr(i_col-col_send(1)+1) = i_val
            endif
         enddo

         col_ptr(bh%n_lcol_sp+1) = bh%nnz_l_sp+1
      endif

      call elsi_deallocate(bh,row_send,"row_send")
      call elsi_deallocate(bh,col_send,"col_send")
   endif

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts density matrix computed by PEXSI, stored in 1D block
!! CCS format to 2D block-cyclic dense format.
!!
subroutine elsi_pexsi_to_blacs_dm_real(ph,bh,row_ind,col_ptr,dm_csc,dm_den)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   real(kind=r8),      intent(in)    :: dm_csc(bh%nnz_l_sp)
   real(kind=r8),      intent(out)   :: dm_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   integer(kind=i4)   :: proc_col_id
   integer(kind=i4)   :: proc_row_id
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character(len=40), parameter :: caller = "elsi_pexsi_to_blacs_dm_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,val_send,bh%nnz_l_sp,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(ph%pexsi_my_prow == 0) then
      call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)

      i_col = 0
      ! Compute destination and global id
      do i_val = 1,bh%nnz_l_sp
         if(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp) then
            i_col = i_col+1
         endif
         i_row = row_ind(i_val)

         ! Compute global id
         row_send(i_val) = i_row
         col_send(i_val) = i_col+bh%myid*(ph%n_basis/ph%pexsi_np_per_pole)
         val_send(i_val) = dm_csc(i_val)

         ! Compute destination
         proc_row_id = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
         proc_col_id = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
         dest(i_val) = proc_col_id+proc_row_id*bh%n_pcol

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      enddo

      ! Sort
      call elsi_heapsort(bh%nnz_l_sp,dest,val_send,row_send,col_send)

      call elsi_deallocate(bh,dest,"dest")
   endif

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Value
   call elsi_allocate(bh,val_recv,bh%nnz_l,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_real8,val_recv,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")

   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   dm_den = 0.0_r8

   ! Unpack density matrix
   do i_val = 1,bh%nnz_l
      ! Compute local 2d id
      call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
      call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

      ! Put value to correct position
      dm_den(l_row,l_col) = val_recv(i_val)
   enddo

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts density matrix computed by PEXSI, stored in 1D block
!! CCS format to 2D block-cyclic dense format.
!!
subroutine elsi_pexsi_to_blacs_dm_cmplx(ph,bh,row_ind,col_ptr,dm_csc,dm_den)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   complex(kind=r8),   intent(in)    :: dm_csc(bh%nnz_l_sp)
   complex(kind=r8),   intent(out)   :: dm_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   integer(kind=i4)   :: proc_col_id
   integer(kind=i4)   :: proc_row_id
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   complex(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character(len=40), parameter :: caller = "elsi_pexsi_to_blacs_dm_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,val_send,bh%nnz_l_sp,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(ph%pexsi_my_prow == 0) then
      call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)

      i_col = 0
      ! Compute destination and global id
      do i_val = 1,bh%nnz_l_sp
         if(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp) then
            i_col = i_col+1
         endif
         i_row = row_ind(i_val)

         ! Compute global id
         row_send(i_val) = i_row
         col_send(i_val) = i_col+bh%myid*(ph%n_basis/ph%pexsi_np_per_pole)
         val_send(i_val) = dm_csc(i_val)

         ! Compute destination
         proc_row_id = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
         proc_col_id = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
         dest(i_val) = proc_col_id+proc_row_id*bh%n_pcol

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      enddo

      ! Sort
      call elsi_heapsort(bh%nnz_l_sp,dest,val_send,row_send,col_send)

      call elsi_deallocate(bh,dest,"dest")
   endif

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Value
   call elsi_allocate(bh,val_recv,bh%nnz_l,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_complex16,val_recv,&
           recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")

   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   dm_den = (0.0_r8,0.0_r8)

   ! Unpack density matrix
   do i_val = 1,bh%nnz_l
      ! Compute local 2d id
      call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
      call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

      ! Put value to correct position
      dm_den(l_row,l_col) = val_recv(i_val)
   enddo

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine gets local number of nonzero elements in matrices in 1D block
!! CSC format converted from 2D block-cyclic dense format.
!!
subroutine elsi_blacs_to_sips_hs_dim_real(ph,bh,ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8),      intent(in)    :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)    :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=40), parameter :: caller = "elsi_blacs_to_sips_hs_dim_real"

   call elsi_allocate(bh,dest,bh%n_procs,"dest",caller)
   call elsi_allocate(bh,nnz,bh%n_procs,"nnz",caller)

   if(.not. ph%ovlp_is_unit) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/bh%n_procs)
         i_proc = min(i_proc,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            endif
         enddo
      enddo
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/bh%n_procs)
         i_proc = min(i_proc,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            endif
         enddo
      enddo
   endif

   call MPI_Allreduce(dest,nnz,bh%n_procs,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   bh%nnz_l_sp  = nnz(bh%myid+1)
   bh%nnz_l_sp1 = bh%nnz_l_sp

   call MPI_Allreduce(bh%nnz_l_sp,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! This routine gets local number of nonzero elements in matrices in 1D block
!! CSC format converted from 2D block-cyclic dense format.
!!
subroutine elsi_blacs_to_sips_hs_dim_cmplx(ph,bh,ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8),   intent(in)    :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(in)    :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=40), parameter :: caller = "elsi_blacs_to_sips_hs_dim_cmplx"

   call elsi_allocate(bh,dest,bh%n_procs,"dest",caller)
   call elsi_allocate(bh,nnz,bh%n_procs,"nnz",caller)

   if(.not. ph%ovlp_is_unit) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/bh%n_procs)
         i_proc = min(i_proc,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            endif
         enddo
      enddo
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/bh%n_procs)
         i_proc = min(i_proc,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            endif
         enddo
      enddo
   endif

   call MPI_Allreduce(dest,nnz,bh%n_procs,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   bh%nnz_l_sp  = nnz(bh%myid+1)
   bh%nnz_l_sp1 = bh%nnz_l_sp

   call MPI_Allreduce(bh%nnz_l_sp,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 2D block-
!! cyclic dense format to 1D block CCS format.
!!
subroutine elsi_blacs_to_sips_hs_real(ph,bh,ham_den,ovlp_den,ham_csc,ovlp_csc,&
              row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8),      intent(in)    :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(in)    :: ovlp_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: ham_csc(bh%nnz_l_sp)
   real(kind=r8),      intent(out)   :: ovlp_csc(bh%nnz_l_sp)
   integer(kind=i4),   intent(out)   :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(out)   :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: g_col
   integer(kind=i4)   :: g_row
   integer(kind=i4)   :: dest ! Destination of an element
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send(:)
   real(kind=r8),    allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: global_id(:)

   character(len=40), parameter :: caller = "elsi_blacs_to_sips_hs_real"

   call elsi_get_time(t0)

   if(ph%n_calls == 1+ph%sips_n_elpa) then
      if(.not. ph%ovlp_is_unit) then
         call elsi_get_nnz(bh%def0,ovlp_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)

         call elsi_allocate(bh,s_val_send,bh%nnz_l,"s_val_send",caller)
      else
         call elsi_get_nnz(bh%def0,ham_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)
      endif
   endif

   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l,"h_val_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0

   if(.not. ph%ovlp_is_unit) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = (g_col-1)/(ph%n_basis/bh%n_procs)
         ! The last process may take more
         dest = min(dest,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val)   = g_row
               col_send(i_val)   = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)
               if(ph%n_calls == 1+ph%sips_n_elpa) then
                  s_val_send(i_val) = ovlp_den(i_row,i_col)
               endif

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = (g_col-1)/(ph%n_basis/bh%n_procs)
         ! The last process may take more
         dest = min(dest,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val)   = g_row
               col_send(i_val)   = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   endif

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(bh,row_recv,bh%nnz_l_sp,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column id
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_real8,ham_csc,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%n_calls == 1+ph%sips_n_elpa .and. .not. ph%ovlp_is_unit) then
      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,ovlp_csc,&
              recv_count,recv_displ,mpi_real8,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   endif

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   call elsi_allocate(bh,global_id,bh%nnz_l_sp,"global_id",caller)

   ! Compute global 1D id
   global_id = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+&
                  int(row_recv,kind=i8)

   ! Sort
   if(ph%n_calls == 1+ph%sips_n_elpa .and. .not. ph%ovlp_is_unit) then
      call elsi_heapsort(bh%nnz_l_sp,global_id,ham_csc,ovlp_csc,row_recv,&
              col_recv)
   else
      call elsi_heapsort(bh%nnz_l_sp,global_id,ham_csc,row_recv,col_recv)
   endif

   call elsi_deallocate(bh,global_id,"global_id")

   ! Compute row index and column pointer
   if(ph%n_calls == 1+ph%sips_n_elpa) then
      i_col = col_recv(1)-1
      do i_val = 1,bh%nnz_l_sp
         row_ind(i_val) = row_recv(i_val)

         if(col_recv(i_val) > i_col) then
            i_col = i_col+1
            col_ptr(i_col-col_recv(1)+1) = i_val
         endif
      enddo

      col_ptr(bh%n_lcol_sp+1) = bh%nnz_l_sp+1
   endif

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 2D block-
!! cyclic dense format to 1D block CCS format.
!!
subroutine elsi_blacs_to_sips_hs_cmplx(ph,bh,ham_den,ovlp_den,ham_csc,ovlp_csc,&
              row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8),   intent(in)    :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(in)    :: ovlp_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(out)   :: ham_csc(bh%nnz_l_sp)
   complex(kind=r8),   intent(out)   :: ovlp_csc(bh%nnz_l_sp)
   integer(kind=i4),   intent(out)   :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(out)   :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: g_col
   integer(kind=i4)   :: g_row
   integer(kind=i4)   :: dest ! Destination of an element
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: h_val_send(:)
   complex(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: global_id(:)

   character(len=40), parameter :: caller = "elsi_blacs_to_sips_hs_cmplx"

   call elsi_get_time(t0)

   if(ph%n_calls == 1+ph%sips_n_elpa) then
      if(.not. ph%ovlp_is_unit) then
         call elsi_get_nnz(bh%def0,ovlp_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)

         call elsi_allocate(bh,s_val_send,bh%nnz_l,"s_val_send",caller)
      else
         call elsi_get_nnz(bh%def0,ham_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)
      endif
   endif

   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l,"h_val_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0

   if(.not. ph%ovlp_is_unit) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = (g_col-1)/(ph%n_basis/bh%n_procs)
         ! The last process may take more
         dest = min(dest,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val)   = g_row
               col_send(i_val)   = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)
               if(ph%n_calls == 1+ph%sips_n_elpa) then
                  s_val_send(i_val) = ovlp_den(i_row,i_col)
               endif

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = (g_col-1)/(ph%n_basis/bh%n_procs)
         ! The last process may take more
         dest = min(dest,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val)   = g_row
               col_send(i_val)   = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   endif

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(bh,row_recv,bh%nnz_l_sp,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column id
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_complex16,ham_csc,&
           recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%n_calls == 1+ph%sips_n_elpa .and. .not. ph%ovlp_is_unit) then
      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
              ovlp_csc,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   endif

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   call elsi_allocate(bh,global_id,bh%nnz_l_sp,"global_id",caller)

   ! Compute global 1D id
   global_id = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+&
                  int(row_recv,kind=i8)

   ! Sort
   if(ph%n_calls == 1+ph%sips_n_elpa .and. .not. ph%ovlp_is_unit) then
      call elsi_heapsort(bh%nnz_l_sp,global_id,ham_csc,ovlp_csc,row_recv,&
              col_recv)
   else
      call elsi_heapsort(bh%nnz_l_sp,global_id,ham_csc,row_recv,col_recv)
   endif

   call elsi_deallocate(bh,global_id,"global_id")

   ! Compute row index and column pointer
   if(ph%n_calls == 1+ph%sips_n_elpa) then
      i_col = col_recv(1)-1
      do i_val = 1,bh%nnz_l_sp
         row_ind(i_val) = row_recv(i_val)

         if(col_recv(i_val) > i_col) then
            i_col = i_col+1
            col_ptr(i_col-col_recv(1)+1) = i_val
         endif
      enddo

      col_ptr(bh%n_lcol_sp+1) = bh%nnz_l_sp+1
   endif

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 1D block
!! CCS format to 2D block-cyclic dense format.
!!
subroutine elsi_sips_to_blacs_hs_real(ph,bh,row_ind,col_ptr,ham_csc,ovlp_csc,&
              ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   real(kind=r8),      intent(in)    :: ham_csc(bh%nnz_l_sp)
   real(kind=r8),      intent(in)    :: ovlp_csc(bh%nnz_l_sp)
   real(kind=r8),      intent(out)   :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   integer(kind=i4)   :: proc_col_id
   integer(kind=i4)   :: proc_row_id
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send(:)
   real(kind=r8),    allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: h_val_recv(:)
   real(kind=r8),    allocatable :: s_val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character(len=40), parameter :: caller = "elsi_sips_to_blacs_hs_real"

   call elsi_get_time(t0)

   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp,"s_val_send",caller)
   endif

   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)
   call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)

   i_col = 0
   ! Compute destination and global id
   do i_val = 1,bh%nnz_l_sp
      if(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp) then
         i_col = i_col+1
      endif
      i_row = row_ind(i_val)

      ! Compute global id
      row_send(i_val)   = i_row
      col_send(i_val)   = i_col+bh%myid*(ph%n_basis/bh%n_procs)
      h_val_send(i_val) = ham_csc(i_val)
      if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
         s_val_send(i_val) = ovlp_csc(i_val)
      endif

      ! Compute destination
      proc_row_id = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
      proc_col_id = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
      dest(i_val) = proc_col_id+proc_row_id*bh%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   enddo

   ! Sort
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_heapsort(bh%nnz_l_sp,dest,h_val_send,s_val_send,row_send,&
              col_send)
   else
      call elsi_heapsort(bh%nnz_l_sp,dest,h_val_send,row_send,col_send)
   endif

   call elsi_deallocate(bh,dest,"dest")

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian Value
   call elsi_allocate(bh,h_val_recv,bh%nnz_l,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_real8,h_val_recv,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_recv,bh%nnz_l,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,s_val_recv,&
              recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)
   endif

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   ! Unpack matrix
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      ham_den  = 0.0_r8
      ovlp_den = 0.0_r8

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col)  = h_val_recv(i_val)
         ovlp_den(l_row,l_col) = s_val_recv(i_val)
      enddo

      call elsi_deallocate(bh,s_val_recv,"s_val_recv")
   else
      ham_den = 0.0_r8

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
      enddo
   endif

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 1D block
!! CCS format to 2D block-cyclic dense format.
!!
subroutine elsi_sips_to_blacs_hs_cmplx(ph,bh,row_ind,col_ptr,ham_csc,ovlp_csc,&
              ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   complex(kind=r8),   intent(in)    :: ham_csc(bh%nnz_l_sp)
   complex(kind=r8),   intent(in)    :: ovlp_csc(bh%nnz_l_sp)
   complex(kind=r8),   intent(out)   :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(out)   :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   integer(kind=i4)   :: proc_col_id
   integer(kind=i4)   :: proc_row_id
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: h_val_send(:)
   complex(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   complex(kind=r8), allocatable :: h_val_recv(:)
   complex(kind=r8), allocatable :: s_val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character(len=40), parameter :: caller = "elsi_sips_to_blacs_hs_cmplx"

   call elsi_get_time(t0)

   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp,"s_val_send",caller)
   endif

   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)
   call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)

   i_col = 0
   ! Compute destination and global id
   do i_val = 1,bh%nnz_l_sp
      if(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp) then
         i_col = i_col+1
      endif
      i_row = row_ind(i_val)

      ! Compute global id
      row_send(i_val)   = i_row
      col_send(i_val)   = i_col+bh%myid*(ph%n_basis/bh%n_procs)
      h_val_send(i_val) = ham_csc(i_val)
      if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
         s_val_send(i_val) = ovlp_csc(i_val)
      endif

      ! Compute destination
      proc_row_id = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
      proc_col_id = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
      dest(i_val) = proc_col_id+proc_row_id*bh%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   enddo

   ! Sort
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_heapsort(bh%nnz_l_sp,dest,h_val_send,s_val_send,row_send,&
              col_send)
   else
      call elsi_heapsort(bh%nnz_l_sp,dest,h_val_send,row_send,col_send)
   endif

   call elsi_deallocate(bh,dest,"dest")

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian Value
   call elsi_allocate(bh,h_val_recv,bh%nnz_l,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_complex16,&
           h_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_recv,bh%nnz_l,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
              s_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)
   endif

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   ! Unpack matrix
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      ham_den  = (0.0_r8,0.0_r8)
      ovlp_den = (0.0_r8,0.0_r8)

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col)  = h_val_recv(i_val)
         ovlp_den(l_row,l_col) = s_val_recv(i_val)
      enddo

      call elsi_deallocate(bh,s_val_recv,"s_val_recv")
   else
      ham_den = (0.0_r8,0.0_r8)

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
      enddo
   endif

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts eigenvectors stored in 1D block dense format to 2D
!! block-cyclic dense format.
!!
subroutine elsi_sips_to_blacs_ev_real(ph,bh,evec_sips,evec)

   implicit none

   type(elsi_param_t), intent(in)  :: ph
   type(elsi_basic_t), intent(in)  :: bh
   real(kind=r8),      intent(in)  :: evec_sips(bh%n_lcol_sp1,ph%n_states)
   real(kind=r8),      intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   integer(kind=i4)   :: proc_col_id
   integer(kind=i4)   :: proc_row_id
   integer(kind=i4)   :: nnz_before
   integer(kind=i4)   :: nnz_after
   integer(kind=i4)   :: n_lrow_aux
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character(len=40), parameter :: caller = "elsi_sips_to_blacs_ev_real"

   call elsi_get_time(t0)

   n_lrow_aux = ph%n_basis/bh%n_procs
   nnz_before = 0

   do i_col = 1,ph%n_states
      do i_row = 1,bh%n_lcol_sp
         if(abs(evec_sips(i_row,i_col)) > bh%def0) then
            nnz_before = nnz_before+1
         endif
      enddo
   enddo

   call elsi_allocate(bh,val_send,nnz_before,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_before,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_before,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)
   call elsi_allocate(bh,dest,nnz_before,"dest",caller)

   i_val = 0
   ! Compute destination and global id
   do i_col = 1,ph%n_states
      do i_row = 1,bh%n_lcol_sp
         if(abs(evec_sips(i_row,i_col)) > bh%def0) then
            i_val = i_val+1

            ! Compute global id
            col_send(i_val) = i_col+ph%sips_first_ev-1
            row_send(i_val) = bh%myid*n_lrow_aux+i_row
            val_send(i_val) = evec_sips(i_row,i_col)

            ! Compute destination
            proc_row_id = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
            proc_col_id = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
            dest(i_val) = proc_col_id+proc_row_id*bh%n_pcol

            ! Set send_count
            send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
         endif
      enddo
   enddo

   ! Sort
   call elsi_heapsort(nnz_before,dest,val_send,row_send,col_send)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   nnz_after = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Value
   call elsi_allocate(bh,val_recv,nnz_after,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_real8,val_recv,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")

   ! Row index
   call elsi_allocate(bh,row_recv,nnz_after,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,nnz_after,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   evec = 0.0_r8

   ! Unpack eigenvectors
   do i_val = 1,nnz_after
      ! Compute local 2d id
      call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
      call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

      ! Put value to correct position
      evec(l_row,l_col) = val_recv(i_val)
   enddo

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts density matrix in 2D block-cyclic dense format to 1D
!! block CCS format.
!!
subroutine elsi_blacs_to_sips_dm_real(ph,bh,row_ind,col_ptr,dm_den,dm_csc)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   real(kind=r8),      intent(in)    :: dm_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: dm_csc(bh%nnz_l_sp)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: j_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   integer(kind=i4)   :: dest ! Destination of an element
   integer(kind=i4)   :: nnz_l_aux
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)

   character(len=40), parameter :: caller = "elsi_blacs_to_sips_dm_real"

   call elsi_get_time(t0)

   call elsi_get_nnz(bh%def0,dm_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)

   call elsi_allocate(bh,val_send,bh%nnz_l,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0

   do i_col = 1,bh%n_lcol
      do i_row = 1,bh%n_lrow
         if(abs(dm_den(i_row,i_col)) > bh%def0) then
            i_val = i_val+1

            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,row_send(i_val))
            call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,col_send(i_val))

            ! Compute destination
            dest = (col_send(i_val)-1)/(ph%n_basis/bh%n_procs)
            ! The last process may take more
            dest = min(dest,bh%n_procs-1)

            ! Pack data
            val_send(i_val) = dm_den(i_row,i_col)

            ! Set send_count
            send_count(dest+1) = send_count(dest+1)+1
         endif
      enddo
   enddo

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(bh,row_recv,nnz_l_aux,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column id
   call elsi_allocate(bh,col_recv,nnz_l_aux,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Density matrix value
   call elsi_allocate(bh,val_recv,nnz_l_aux,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_real8,val_recv,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   dm_csc = 0.0_r8

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      l_col = col_recv(i_val)-(ph%n_basis/bh%n_procs)*bh%myid
      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_csc(j_val) = val_recv(i_val)
         endif
      enddo
   enddo

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts density matrix in 2D block-cyclic dense format to 1D
!! block CCS format.
!!
subroutine elsi_blacs_to_sips_dm_cmplx(ph,bh,row_ind,col_ptr,dm_den,dm_csc)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   complex(kind=r8),   intent(in)    :: dm_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(out)   :: dm_csc(bh%nnz_l_sp)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: j_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   integer(kind=i4)   :: dest ! Destination of an element
   integer(kind=i4)   :: nnz_l_aux
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   complex(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)

   character(len=40), parameter :: caller = "elsi_blacs_to_sips_dm_cmplx"

   call elsi_get_time(t0)

   call elsi_get_nnz(bh%def0,dm_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)

   call elsi_allocate(bh,val_send,bh%nnz_l,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   ! Compute destination and global id
   i_val = 0

   do i_col = 1,bh%n_lcol
      do i_row = 1,bh%n_lrow
         if(abs(dm_den(i_row,i_col)) > bh%def0) then
            i_val = i_val+1

            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,row_send(i_val))
            call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,col_send(i_val))

            ! Compute destination
            dest = (col_send(i_val)-1)/(ph%n_basis/bh%n_procs)
            ! The last process may take more
            dest = min(dest,bh%n_procs-1)

            ! Pack data
            val_send(i_val) = dm_den(i_row,i_col)

            ! Set send_count
            send_count(dest+1) = send_count(dest+1)+1
         endif
      enddo
   enddo

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(bh,row_recv,nnz_l_aux,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column id
   call elsi_allocate(bh,col_recv,nnz_l_aux,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Density matrix value
   call elsi_allocate(bh,val_recv,nnz_l_aux,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_complex16,val_recv,&
           recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   dm_csc = (0.0_r8,0.0_r8)

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      l_col = col_recv(i_val)-(ph%n_basis/bh%n_procs)*bh%myid
      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_csc(j_val) = val_recv(i_val)
         endif
      enddo
   enddo

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 1D block-
!! cyclic CCS format to 2D block-cyclic dense format.
!!
subroutine elsi_siesta_to_blacs_hs_real(ph,bh,row_ind,col_ptr,ham_csc,ovlp_csc,&
              ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   real(kind=r8),      intent(in)    :: ham_csc(bh%nnz_l_sp)
   real(kind=r8),      intent(in)    :: ovlp_csc(bh%nnz_l_sp)
   real(kind=r8),      intent(out)   :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   integer(kind=i4)   :: proc_col_id
   integer(kind=i4)   :: proc_row_id
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send(:)
   real(kind=r8),    allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: h_val_recv(:)
   real(kind=r8),    allocatable :: s_val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character(len=40), parameter :: caller = "elsi_siesta_to_blacs_hs_real"

   call elsi_get_time(t0)

   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp,"s_val_send",caller)
   endif

   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)
   call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)

   i_col = 0
   ! Compute destination and global id
   do i_val = 1,bh%nnz_l_sp
      if(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp2) then
         i_col = i_col+1
      endif
      i_row = row_ind(i_val)

      ! Compute global id
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,col_send(i_val))
      row_send(i_val)   = i_row
      h_val_send(i_val) = ham_csc(i_val)
      if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
         s_val_send(i_val) = ovlp_csc(i_val)
      endif

      ! Compute destination
      proc_row_id = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
      proc_col_id = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
      dest(i_val) = proc_col_id+proc_row_id*bh%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   enddo

   ! Sort
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_heapsort(bh%nnz_l_sp,dest,h_val_send,s_val_send,row_send,&
              col_send)
   else
      call elsi_heapsort(bh%nnz_l_sp,dest,h_val_send,row_send,col_send)
   endif

   call elsi_deallocate(bh,dest,"dest")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian Value
   call elsi_allocate(bh,h_val_recv,bh%nnz_l,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_real8,h_val_recv,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_recv,bh%nnz_l,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,s_val_recv,&
              recv_count,recv_displ,mpi_real8,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)
   endif

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   ! Unpack matrix
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      ham_den  = 0.0_r8
      ovlp_den = 0.0_r8

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col)  = h_val_recv(i_val)
         ovlp_den(l_row,l_col) = s_val_recv(i_val)
      enddo

      call elsi_deallocate(bh,s_val_recv,"s_val_recv")
   else
      ham_den = 0.0_r8

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
      enddo
   endif

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 1D block-
!! cyclic CCS format to 2D block-cyclic dense format.
!!
subroutine elsi_siesta_to_blacs_hs_cmplx(ph,bh,row_ind,col_ptr,ham_csc,&
              ovlp_csc,ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   complex(kind=r8),   intent(in)    :: ham_csc(bh%nnz_l_sp)
   complex(kind=r8),   intent(in)    :: ovlp_csc(bh%nnz_l_sp)
   complex(kind=r8),   intent(out)   :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(out)   :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   integer(kind=i4)   :: proc_col_id
   integer(kind=i4)   :: proc_row_id
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: h_val_send(:)
   complex(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   complex(kind=r8), allocatable :: h_val_recv(:)
   complex(kind=r8), allocatable :: s_val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character(len=40), parameter :: caller = "elsi_siesta_to_blacs_hs_cmplx"

   call elsi_get_time(t0)

   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp,"s_val_send",caller)
   endif

   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)
   call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)

   i_col = 0
   ! Compute destination and global id
   do i_val = 1,bh%nnz_l_sp
      if(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp2) then
         i_col = i_col+1
      endif
      i_row = row_ind(i_val)

      ! Compute global id
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,col_send(i_val))
      row_send(i_val)   = i_row
      h_val_send(i_val) = ham_csc(i_val)
      if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
         s_val_send(i_val) = ovlp_csc(i_val)
      endif

      ! Compute destination
      proc_row_id = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
      proc_col_id = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
      dest(i_val) = proc_col_id+proc_row_id*bh%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   enddo

   ! Sort
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_heapsort(bh%nnz_l_sp,dest,h_val_send,s_val_send,row_send,&
              col_send)
   else
      call elsi_heapsort(bh%nnz_l_sp,dest,h_val_send,row_send,col_send)
   endif

   call elsi_deallocate(bh,dest,"dest")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian Value
   call elsi_allocate(bh,h_val_recv,bh%nnz_l,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_complex16,&
           h_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_recv,bh%nnz_l,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
              s_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)
   endif

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   ! Unpack matrix
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      ham_den  = (0.0_r8,0.0_r8)
      ovlp_den = (0.0_r8,0.0_r8)

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col)  = h_val_recv(i_val)
         ovlp_den(l_row,l_col) = s_val_recv(i_val)
      enddo

      call elsi_deallocate(bh,s_val_recv,"s_val_recv")
   else
      ham_den = (0.0_r8,0.0_r8)

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
      enddo
   endif

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts density matrix in 2D block-cyclic dense format to 1D
!! block-cyclic CCS format.
!!
subroutine elsi_blacs_to_siesta_dm_real(bh,row_ind,col_ptr,dm_den,dm_csc)

   implicit none

   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   real(kind=r8),      intent(in)    :: dm_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8),      intent(out)   :: dm_csc(bh%nnz_l_sp)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: j_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   integer(kind=i4)   :: nnz_l_aux
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character(len=40), parameter :: caller = "elsi_blacs_to_siesta_dm_real"

   call elsi_get_time(t0)

   call elsi_get_nnz(bh%def0,dm_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)

   call elsi_allocate(bh,val_send,bh%nnz_l,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)
   call elsi_allocate(bh,dest,bh%nnz_l,"dest",caller)

   ! Compute destination and global id
   i_val = 0

   do i_col = 1,bh%n_lcol
      do i_row = 1,bh%n_lrow
         if(abs(dm_den(i_row,i_col)) > bh%def0) then
            i_val = i_val+1

            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,row_send(i_val))
            call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,col_send(i_val))

            ! Compute destination
            dest(i_val) = mod((col_send(i_val)-1)/bh%blk_sp2,bh%n_procs)

            ! Pack data
            val_send(i_val) = dm_den(i_row,i_col)

            ! Set send_count
            send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
         endif
      enddo
   enddo

   ! Sort
   call elsi_heapsort(bh%nnz_l,dest,val_send,row_send,col_send)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(bh,row_recv,nnz_l_aux,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column id
   call elsi_allocate(bh,col_recv,nnz_l_aux,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Density matrix value
   call elsi_allocate(bh,val_recv,nnz_l_aux,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_real8,val_recv,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   dm_csc = 0.0_r8

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      call elsi_get_lid(bh%n_procs,bh%blk_sp2,col_recv(i_val),l_col)

      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_csc(j_val) = val_recv(i_val)
         endif
      enddo
   enddo

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts density matrix in 2D block-cyclic dense format to 1D
!! block-cyclic CCS format.
!!
subroutine elsi_blacs_to_siesta_dm_cmplx(bh,row_ind,col_ptr,dm_den,dm_csc)

   implicit none

   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   complex(kind=r8),   intent(in)    :: dm_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8),   intent(out)   :: dm_csc(bh%nnz_l_sp)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: j_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   integer(kind=i4)   :: nnz_l_aux
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   complex(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character(len=40), parameter :: caller = "elsi_blacs_to_siesta_dm_cmplx"

   call elsi_get_time(t0)

   call elsi_get_nnz(bh%def0,dm_den,bh%n_lrow,bh%n_lcol,bh%nnz_l)

   call elsi_allocate(bh,val_send,bh%nnz_l,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)
   call elsi_allocate(bh,dest,bh%nnz_l,"dest",caller)

   ! Compute destination and global id
   i_val = 0

   do i_col = 1,bh%n_lcol
      do i_row = 1,bh%n_lrow
         if(abs(dm_den(i_row,i_col)) > bh%def0) then
            i_val = i_val+1

            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,row_send(i_val))
            call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,col_send(i_val))

            ! Compute destination
            dest(i_val) = mod((col_send(i_val)-1)/bh%blk_sp2,bh%n_procs)

            ! Pack data
            val_send(i_val) = dm_den(i_row,i_col)

            ! Set send_count
            send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
         endif
      enddo
   enddo

   ! Sort
   call elsi_heapsort(bh%nnz_l,dest,val_send,row_send,col_send)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(bh,row_recv,nnz_l_aux,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column id
   call elsi_allocate(bh,col_recv,nnz_l_aux,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Density matrix value
   call elsi_allocate(bh,val_recv,nnz_l_aux,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_complex16,val_recv,&
           recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   dm_csc = (0.0_r8,0.0_r8)

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      call elsi_get_lid(bh%n_procs,bh%blk_sp2,col_recv(i_val),l_col)

      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_csc(j_val) = val_recv(i_val)
         endif
      enddo
   enddo

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine gets local number of nonzero elements in matrices in 1D block
!! CSC format converted from 1D block-cyclic dense format.
!!
subroutine elsi_siesta_to_pexsi_hs_dim(ph,bh,col_ptr2)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: col_ptr2(bh%n_lcol_sp2+1)

   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: myid_in_pole
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=40), parameter :: caller = "elsi_siesta_to_pexsi_hs_dim"

   call elsi_allocate(bh,dest,ph%pexsi_np_per_pole,"dest",caller)
   call elsi_allocate(bh,nnz,ph%pexsi_np_per_pole,"nnz",caller)

   do i_col = 1,bh%n_lcol_sp2
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,g_col)

      i_proc = (g_col-1)/(ph%n_basis/ph%pexsi_np_per_pole)
      i_proc = min(i_proc,ph%pexsi_np_per_pole-1)

      dest(i_proc+1) = dest(i_proc+1)+col_ptr2(i_col+1)-col_ptr2(i_col)
   enddo

   call MPI_Allreduce(dest,nnz,bh%n_procs,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   call MPI_Comm_rank(ph%pexsi_comm_in_pole,myid_in_pole,ierr)

   call elsi_check_mpi(bh,"MPI_Comm_rank",ierr,caller)

   bh%nnz_l_sp1 = nnz(myid_in_pole+1)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 1D block-
!! cyclic CCS format to 1D block CSC format.
!!
subroutine elsi_siesta_to_pexsi_hs_real(ph,bh,ham_csc2,ovlp_csc2,row_ind2,&
              col_ptr2,ham_csc1,ovlp_csc1,row_ind1,col_ptr1)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8),      intent(in)    :: ham_csc2(bh%nnz_l_sp2)
   real(kind=r8),      intent(in)    :: ovlp_csc2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: col_ptr2(bh%nnz_l_sp2)
   real(kind=r8),      intent(out)   :: ham_csc1(bh%nnz_l_sp1)
   real(kind=r8),      intent(out)   :: ovlp_csc1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(out)   :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(out)   :: col_ptr1(bh%nnz_l_sp1)

   integer(kind=i4)   :: n_para_task
   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: dest ! Destination of an element
   integer(kind=i4)   :: n_lcol_aux
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send(:)
   real(kind=r8),    allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: global_id(:) ! Global 1D id

   character(len=40), parameter :: caller = "elsi_siesta_to_pexsi_hs_real"

   call elsi_get_time(t0)

   n_para_task = bh%n_procs/ph%pexsi_np_per_pole
   n_lcol_aux  = ph%n_basis/ph%pexsi_np_per_pole

   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp2,"s_val_send",caller)
   endif

   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp2,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp2,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp2,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_col = 0
   ! Compute destination and global id
   do i_val = 1,bh%nnz_l_sp2
      if(i_val == col_ptr2(i_col+1) .and. i_col /= bh%n_lcol_sp2) then
         i_col = i_col+1
      endif
      i_row = row_ind2(i_val)

      ! Compute global id
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,col_send(i_val))
      row_send(i_val)   = i_row
      h_val_send(i_val) = ham_csc2(i_val)
      if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
         s_val_send(i_val) = ovlp_csc2(i_val)
      endif

      ! Compute destination
      dest = (col_send(i_val)-1)/n_lcol_aux
      ! The last process may take more
      dest = min(dest,ph%pexsi_np_per_pole-1)

      ! Set send_count
      send_count(dest+1) = send_count(dest+1)+1
   enddo

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l_sp1,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp1,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian Value
   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_real8,ham_csc1,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,ovlp_csc1,&
              recv_count,recv_displ,mpi_real8,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   endif

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   call elsi_allocate(bh,global_id,bh%nnz_l_sp1,"global_id",caller)

   ! Compute global 1D id
   global_id = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+&
                  int(row_recv,kind=i8)

   ! Sort
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_heapsort(bh%nnz_l_sp1,global_id,ham_csc1,ovlp_csc1,row_recv,&
              col_recv)
   else
      call elsi_heapsort(bh%nnz_l_sp1,global_id,ham_csc1,row_recv,col_recv)
   endif

   call elsi_deallocate(bh,global_id,"global_id")

   if(ph%n_calls == 1) then
      ! Only 1st pole computes row index and column pointer
      if(ph%pexsi_my_prow == 0) then
         ! Compute row index and column pointer
         i_col = col_recv(1)-1
         do i_val = 1,bh%nnz_l_sp1
            row_ind1(i_val) = row_recv(i_val)

            if(col_recv(i_val) > i_col) then
               i_col = i_col+1
               col_ptr1(i_col-col_recv(1)+1) = i_val
            endif
         enddo

         col_ptr1(bh%n_lcol_sp1+1) = bh%nnz_l_sp1+1
      endif
   endif

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 1D block-
!! cyclic CCS format to 1D block CSC format.
!!
subroutine elsi_siesta_to_pexsi_hs_cmplx(ph,bh,ham_csc2,ovlp_csc2,row_ind2,&
              col_ptr2,ham_csc1,ovlp_csc1,row_ind1,col_ptr1)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8),   intent(in)    :: ham_csc2(bh%nnz_l_sp2)
   complex(kind=r8),   intent(in)    :: ovlp_csc2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: col_ptr2(bh%nnz_l_sp2)
   complex(kind=r8),   intent(out)   :: ham_csc1(bh%nnz_l_sp1)
   complex(kind=r8),   intent(out)   :: ovlp_csc1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(out)   :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(out)   :: col_ptr1(bh%nnz_l_sp1)

   integer(kind=i4)   :: n_para_task
   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: dest ! Destination of an element
   integer(kind=i4)   :: n_lcol_aux
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: h_val_send(:)
   complex(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: global_id(:) ! Global 1D id

   character(len=40), parameter :: caller = "elsi_siesta_to_pexsi_hs_cmplx"

   call elsi_get_time(t0)

   n_para_task = bh%n_procs/ph%pexsi_np_per_pole
   n_lcol_aux  = ph%n_basis/ph%pexsi_np_per_pole

   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp2,"s_val_send",caller)
   endif

   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp2,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp2,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp2,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_col = 0
   ! Compute destination and global id
   do i_val = 1,bh%nnz_l_sp2
      if(i_val == col_ptr2(i_col+1) .and. i_col /= bh%n_lcol_sp2) then
         i_col = i_col+1
      endif
      i_row = row_ind2(i_val)

      ! Compute global id
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,col_send(i_val))
      row_send(i_val)   = i_row
      h_val_send(i_val) = ham_csc2(i_val)
      if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
         s_val_send(i_val) = ovlp_csc2(i_val)
      endif

      ! Compute destination
      dest = (col_send(i_val)-1)/n_lcol_aux
      ! The last process may take more
      dest = min(dest,ph%pexsi_np_per_pole-1)

      ! Set send_count
      send_count(dest+1) = send_count(dest+1)+1
   enddo

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l_sp1,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp1,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian Value
   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_complex16,ham_csc1,&
           recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
              ovlp_csc1,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   endif

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   call elsi_allocate(bh,global_id,bh%nnz_l_sp1,"global_id",caller)

   ! Compute global 1D id
   global_id = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+&
                  int(row_recv,kind=i8)

   ! Sort
   if(ph%n_calls == 1 .and. .not. ph%ovlp_is_unit) then
      call elsi_heapsort(bh%nnz_l_sp1,global_id,ham_csc1,ovlp_csc1,row_recv,&
              col_recv)
   else
      call elsi_heapsort(bh%nnz_l_sp1,global_id,ham_csc1,row_recv,col_recv)
   endif

   call elsi_deallocate(bh,global_id,"global_id")

   if(ph%n_calls == 1) then
      ! Only 1st pole computes row index and column pointer
      if(ph%pexsi_my_prow == 0) then
         ! Compute row index and column pointer
         i_col = col_recv(1)-1
         do i_val = 1,bh%nnz_l_sp1
            row_ind1(i_val) = row_recv(i_val)

            if(col_recv(i_val) > i_col) then
               i_col = i_col+1
               col_ptr1(i_col-col_recv(1)+1) = i_val
            endif
         enddo

         col_ptr1(bh%n_lcol_sp1+1) = bh%nnz_l_sp1+1
      endif
   endif

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts density matrix computed by PEXSI, stored in 1D block
!! CCS format to 1D block-cyclic CCS format.
!!
subroutine elsi_pexsi_to_siesta_dm_real(ph,bh,row_ind1,col_ptr1,dm_csc1,&
              row_ind2,col_ptr2,dm_csc2)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)    :: col_ptr1(bh%n_lcol_sp1+1)
   real(kind=r8),      intent(in)    :: dm_csc1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)    :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: col_ptr2(bh%n_lcol_sp2+1)
   real(kind=r8),      intent(out)   :: dm_csc2(bh%nnz_l_sp2)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: j_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8),    allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character(len=40), parameter :: caller = "elsi_pexsi_to_siesta_dm_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,val_send,bh%nnz_l_sp1,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp1,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp1,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(ph%pexsi_my_prow == 0) then
      call elsi_allocate(bh,dest,bh%nnz_l_sp1,"dest",caller)

      i_col = 0
      ! Compute destination and global id
      do i_val = 1,bh%nnz_l_sp1
         if(i_val == col_ptr1(i_col+1) .and. i_col /= bh%n_lcol_sp1) then
            i_col = i_col+1
         endif
         i_row = row_ind1(i_val)

         ! Compute global id
         row_send(i_val) = i_row
         col_send(i_val) = i_col+bh%myid*(ph%n_basis/ph%pexsi_np_per_pole)
         val_send(i_val) = dm_csc1(i_val)

         ! Compute destination
         dest(i_val) = mod((col_send(i_val)-1)/bh%blk_sp2,bh%n_procs)

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      enddo

      ! Sort
      call elsi_heapsort(bh%nnz_l_sp1,dest,val_send,row_send,col_send)

      call elsi_deallocate(bh,dest,"dest")
   endif

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l_sp2 = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Value
   call elsi_allocate(bh,val_recv,bh%nnz_l_sp2,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_real8,val_recv,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")

   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l_sp2,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp2,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   dm_csc2 = 0.0_r8

   ! Unpack matrix
   do i_val = 1,bh%nnz_l_sp2
      call elsi_get_lid(bh%n_procs,bh%blk_sp2,col_recv(i_val),l_col)

      l_row = row_recv(i_val)

      do j_val = col_ptr2(l_col),col_ptr2(l_col+1)-1
         if(row_ind2(j_val) == l_row) then
            dm_csc2(j_val) = val_recv(i_val)
         endif
      enddo
   enddo

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts density matrix computed by PEXSI, stored in 1D block
!! CCS format to 1D block-cyclic CCS format.
!!
subroutine elsi_pexsi_to_siesta_dm_cmplx(ph,bh,row_ind1,col_ptr1,dm_csc1,&
              row_ind2,col_ptr2,dm_csc2)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)    :: col_ptr1(bh%n_lcol_sp1+1)
   complex(kind=r8),   intent(in)    :: dm_csc1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)    :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: col_ptr2(bh%n_lcol_sp2+1)
   complex(kind=r8),   intent(out)   :: dm_csc2(bh%nnz_l_sp2)

   integer(kind=i4)   :: ierr
   integer(kind=i4)   :: i_row
   integer(kind=i4)   :: i_col
   integer(kind=i4)   :: i_val
   integer(kind=i4)   :: j_val
   integer(kind=i4)   :: i_proc
   integer(kind=i4)   :: l_col ! Local column id in 1D block distribution
   integer(kind=i4)   :: l_row ! Local row id in 1D block distribution
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   complex(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   character(len=40), parameter :: caller = "elsi_pexsi_to_siesta_dm_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,val_send,bh%nnz_l_sp1,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp1,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp1,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(ph%pexsi_my_prow == 0) then
      call elsi_allocate(bh,dest,bh%nnz_l_sp1,"dest",caller)

      i_col = 0
      ! Compute destination and global id
      do i_val = 1,bh%nnz_l_sp1
         if(i_val == col_ptr1(i_col+1) .and. i_col /= bh%n_lcol_sp1) then
            i_col = i_col+1
         endif
         i_row = row_ind1(i_val)

         ! Compute global id
         row_send(i_val) = i_row
         col_send(i_val) = i_col+bh%myid*(ph%n_basis/ph%pexsi_np_per_pole)
         val_send(i_val) = dm_csc1(i_val)

         ! Compute destination
         dest(i_val) = mod((col_send(i_val)-1)/bh%blk_sp2,bh%n_procs)

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      enddo

      ! Sort
      call elsi_heapsort(bh%nnz_l_sp1,dest,val_send,row_send,col_send)

      call elsi_deallocate(bh,dest,"dest")
   endif

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
           bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l_sp2 = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   enddo

   ! Send and receive the packed data
   ! Value
   call elsi_allocate(bh,val_recv,bh%nnz_l_sp2,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_complex16,val_recv,&
           recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")

   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l_sp2,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp2,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   dm_csc2 = (0.0_r8,0.0_r8)

   ! Unpack matrix
   do i_val = 1,bh%nnz_l_sp2
      call elsi_get_lid(bh%n_procs,bh%blk_sp2,col_recv(i_val),l_col)

      l_row = row_recv(i_val)

      do j_val = col_ptr2(l_col),col_ptr2(l_col+1)-1
         if(row_ind2(j_val) == l_row) then
            dm_csc2(j_val) = val_recv(i_val)
         endif
      enddo
   enddo

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished matrix redistribution"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine converts density matrix computed by SIPS, stored in 1D block
!! CCS format to 2D block-cyclic dense format.
!!
subroutine elsi_sips_to_blacs_dm_real(ph,bh,row_ind,col_ptr,dm_csc,dm_den)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   real(kind=r8),      intent(in)    :: dm_csc(bh%nnz_l_sp)
   real(kind=r8),      intent(out)   :: dm_den(bh%n_lrow,bh%n_lcol)

   character(len=40), parameter :: caller = "elsi_sips_to_blacs_dm_real"

   ph%pexsi_my_prow     = 0
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_pexsi_to_blacs_dm_real(ph,bh,row_ind,col_ptr,dm_csc,dm_den)

end subroutine

!>
!! This routine converts density matrix computed by SIPS, stored in 1D block
!! CCS format to 2D block-cyclic dense format.
!!
subroutine elsi_sips_to_blacs_dm_cmplx(ph,bh,row_ind,col_ptr,dm_csc,dm_den)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp+1)
   complex(kind=r8),   intent(in)    :: dm_csc(bh%nnz_l_sp)
   complex(kind=r8),   intent(out)   :: dm_den(bh%n_lrow,bh%n_lcol)

   character(len=40), parameter :: caller = "elsi_sips_to_blacs_dm_cmplx"

   ph%pexsi_my_prow     = 0
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_pexsi_to_blacs_dm_cmplx(ph,bh,row_ind,col_ptr,dm_csc,dm_den)

end subroutine

!>
!! This routine converts density matrix computed by SIPS, stored in 1D block
!! CCS format to 1D block-cyclic CCS format.
!!
subroutine elsi_sips_to_siesta_dm_real(ph,bh,row_ind1,col_ptr1,dm_csc1,&
              row_ind2,col_ptr2,dm_csc2)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)    :: col_ptr1(bh%n_lcol_sp1+1)
   real(kind=r8),      intent(in)    :: dm_csc1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)    :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: col_ptr2(bh%n_lcol_sp2+1)
   real(kind=r8),      intent(out)   :: dm_csc2(bh%nnz_l_sp2)

   character(len=40), parameter :: caller = "elsi_sips_to_siesta_dm_real"

   ph%pexsi_my_prow      = 0
   ph%pexsi_np_per_pole  = bh%n_procs

   call elsi_pexsi_to_siesta_dm_real(ph,bh,row_ind1,col_ptr1,dm_csc1,row_ind2,&
           col_ptr2,dm_csc2)

end subroutine

!>
!! This routine converts density matrix computed by SIPS, stored in 1D block
!! CCS format to 1D block-cyclic CCS format.
!!
subroutine elsi_sips_to_siesta_dm_cmplx(ph,bh,row_ind1,col_ptr1,dm_csc1,&
              row_ind2,col_ptr2,dm_csc2)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)    :: col_ptr1(bh%n_lcol_sp1+1)
   complex(kind=r8),   intent(in)    :: dm_csc1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)    :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: col_ptr2(bh%n_lcol_sp2+1)
   complex(kind=r8),   intent(out)   :: dm_csc2(bh%nnz_l_sp2)

   character(len=40), parameter :: caller = "elsi_sips_to_siesta_dm_cmplx"

   ph%pexsi_my_prow     = 0
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_pexsi_to_siesta_dm_cmplx(ph,bh,row_ind1,col_ptr1,dm_csc1,row_ind2,&
           col_ptr2,dm_csc2)

end subroutine

!>
!! This routine gets local number of nonzero elements in matrices in 1D block
!! CSC format converted from 1D block-cyclic dense format.
!!
subroutine elsi_siesta_to_sips_hs_dim(ph,bh,col_ptr2)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4),   intent(in)    :: col_ptr2(bh%n_lcol_sp2+1)

   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=40), parameter :: caller = "elsi_siesta_to_sips_hs_dim"

   call elsi_allocate(bh,dest,bh%n_procs,"dest",caller)
   call elsi_allocate(bh,nnz,bh%n_procs,"nnz",caller)

   do i_col = 1,bh%n_lcol_sp2
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,g_col)

      i_proc = (g_col-1)/(ph%n_basis/bh%n_procs)
      i_proc = min(i_proc,bh%n_procs)

      dest(i_proc+1) = dest(i_proc+1)+col_ptr2(i_col+1)-col_ptr2(i_col)
   enddo

   call MPI_Allreduce(dest,nnz,bh%n_procs,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   bh%nnz_l_sp1 = nnz(bh%myid+1)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 1D block-
!! cyclic CCS format to 1D block CSC format.
!!
subroutine elsi_siesta_to_sips_hs_real(ph,bh,ham_csc2,ovlp_csc2,row_ind2,&
              col_ptr2,ham_csc1,ovlp_csc1,row_ind1,col_ptr1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8),      intent(in)    :: ham_csc2(bh%nnz_l_sp2)
   real(kind=r8),      intent(in)    :: ovlp_csc2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: col_ptr2(bh%nnz_l_sp2)
   real(kind=r8),      intent(out)   :: ham_csc1(bh%nnz_l_sp1)
   real(kind=r8),      intent(out)   :: ovlp_csc1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(out)   :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(out)   :: col_ptr1(bh%nnz_l_sp1)

   integer(kind=i4) :: n_calls_save

   character(len=40), parameter :: caller = "elsi_siesta_to_sips_hs_real"

   ph%pexsi_my_prow         = 0
   ph%pexsi_np_per_pole     = bh%n_procs
   ph%pexsi_comm_among_pole = mpi_comm_self
   n_calls_save                = ph%n_calls

   if(ph%n_calls == 1+ph%sips_n_elpa) then
      ph%n_calls = 1
   endif

   call elsi_siesta_to_pexsi_hs_real(ph,bh,ham_csc2,ovlp_csc2,row_ind2,&
           col_ptr2,ham_csc1,ovlp_csc1,row_ind1,col_ptr1)

   ph%n_calls = n_calls_save

end subroutine

!>
!! This routine converts Halmitonian and overlap matrixs stored in 1D block-
!! cyclic CCS format to 1D block CSC format.
!!
subroutine elsi_siesta_to_sips_hs_cmplx(ph,bh,ham_csc2,ovlp_csc2,row_ind2,&
              col_ptr2,ham_csc1,ovlp_csc1,row_ind1,col_ptr1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8),   intent(in)    :: ham_csc2(bh%nnz_l_sp2)
   complex(kind=r8),   intent(in)    :: ovlp_csc2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4),   intent(in)    :: col_ptr2(bh%nnz_l_sp2)
   complex(kind=r8),   intent(out)   :: ham_csc1(bh%nnz_l_sp1)
   complex(kind=r8),   intent(out)   :: ovlp_csc1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(out)   :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4),   intent(out)   :: col_ptr1(bh%nnz_l_sp1)

   integer(kind=i4) :: n_calls_save

   character(len=40), parameter :: caller = "elsi_siesta_to_sips_hs_cmplx"

   ph%pexsi_my_prow         = 0
   ph%pexsi_np_per_pole     = bh%n_procs
   ph%pexsi_comm_among_pole = mpi_comm_self
   n_calls_save                = ph%n_calls

   if(ph%n_calls == 1+ph%sips_n_elpa) then
      ph%n_calls = 1
   endif

   call elsi_siesta_to_pexsi_hs_cmplx(ph,bh,ham_csc2,ovlp_csc2,row_ind2,&
           col_ptr2,ham_csc1,ovlp_csc1,row_ind1,col_ptr1)

   ph%n_calls = n_calls_save

end subroutine

end module ELSI_REDIST
