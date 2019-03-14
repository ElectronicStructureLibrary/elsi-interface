! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides matrix conversion and redistribution routines.
!!
module ELSI_REDIST

   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_IO, only: elsi_say,elsi_get_time
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI, only: elsi_check_mpi,mpi_sum,mpi_real8,mpi_complex16,&
       mpi_integer4
   use ELSI_NTPOLY, only: Triplet_r,Triplet_c,TripletList_r,TripletList_c,&
       Matrix_ps,ConstructEmptyMatrix,FillMatrixFromTripletList,&
       GetMatrixTripletList,ConstructTripletList,AppendToTripletList,&
       DestructTripletList
   use ELSI_PRECISION, only: r8,i4,i8
   use ELSI_SORT, only: elsi_heapsort,elsi_permute,elsi_unpermute
   use ELSI_UTILS, only: elsi_get_nnz,elsi_get_gid,elsi_get_lid

   implicit none

   private

   public :: elsi_blacs_to_generic_dm
   public :: elsi_blacs_to_ntpoly_hs
   public :: elsi_blacs_to_pexsi_hs_dim
   public :: elsi_blacs_to_pexsi_hs
   public :: elsi_blacs_to_siesta_dm
   public :: elsi_blacs_to_sips_dm
   public :: elsi_blacs_to_sips_hs_dim
   public :: elsi_blacs_to_sips_hs
   public :: elsi_generic_to_blacs_hs
   public :: elsi_generic_to_ntpoly_hs
   public :: elsi_generic_to_pexsi_hs_dim
   public :: elsi_generic_to_pexsi_hs
   public :: elsi_generic_to_sips_hs_dim
   public :: elsi_generic_to_sips_hs
   public :: elsi_ntpoly_to_blacs_dm
   public :: elsi_ntpoly_to_generic_dm
   public :: elsi_ntpoly_to_siesta_dm
   public :: elsi_ntpoly_to_sips_dm
   public :: elsi_pexsi_to_blacs_dm
   public :: elsi_pexsi_to_generic_dm
   public :: elsi_pexsi_to_siesta_dm
   public :: elsi_siesta_to_blacs_hs
   public :: elsi_siesta_to_ntpoly_hs
   public :: elsi_siesta_to_pexsi_hs_dim
   public :: elsi_siesta_to_pexsi_hs
   public :: elsi_siesta_to_sips_hs_dim
   public :: elsi_siesta_to_sips_hs
   public :: elsi_sips_to_blacs_dm
   public :: elsi_sips_to_blacs_ev
   public :: elsi_sips_to_blacs_hs
   public :: elsi_sips_to_generic_dm
   public :: elsi_sips_to_ntpoly_hs
   public :: elsi_sips_to_siesta_dm

   interface elsi_blacs_to_generic_dm
      module procedure elsi_blacs_to_generic_dm_real
      module procedure elsi_blacs_to_generic_dm_cmplx
   end interface

   interface elsi_blacs_to_ntpoly_hs
      module procedure elsi_blacs_to_ntpoly_hs_real
      module procedure elsi_blacs_to_ntpoly_hs_cmplx
   end interface

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

   interface elsi_generic_to_blacs_hs
      module procedure elsi_generic_to_blacs_hs_real
      module procedure elsi_generic_to_blacs_hs_cmplx
   end interface

   interface elsi_generic_to_ntpoly_hs
      module procedure elsi_generic_to_ntpoly_hs_real
      module procedure elsi_generic_to_ntpoly_hs_cmplx
   end interface

   interface elsi_generic_to_pexsi_hs
      module procedure elsi_generic_to_pexsi_hs_real
      module procedure elsi_generic_to_pexsi_hs_cmplx
   end interface

   interface elsi_generic_to_sips_hs
      module procedure elsi_generic_to_sips_hs_real
      module procedure elsi_generic_to_sips_hs_cmplx
   end interface

   interface elsi_ntpoly_to_blacs_dm
      module procedure elsi_ntpoly_to_blacs_dm_real
      module procedure elsi_ntpoly_to_blacs_dm_cmplx
   end interface

   interface elsi_ntpoly_to_generic_dm
      module procedure elsi_ntpoly_to_generic_dm_real
      module procedure elsi_ntpoly_to_generic_dm_cmplx
   end interface

   interface elsi_ntpoly_to_siesta_dm
      module procedure elsi_ntpoly_to_siesta_dm_real
      module procedure elsi_ntpoly_to_siesta_dm_cmplx
   end interface

   interface elsi_ntpoly_to_sips_dm
      module procedure elsi_ntpoly_to_sips_dm_real
      module procedure elsi_ntpoly_to_sips_dm_cmplx
   end interface

   interface elsi_pexsi_to_blacs_dm
      module procedure elsi_pexsi_to_blacs_dm_real
      module procedure elsi_pexsi_to_blacs_dm_cmplx
   end interface

   interface elsi_pexsi_to_generic_dm
      module procedure elsi_pexsi_to_generic_dm_real
      module procedure elsi_pexsi_to_generic_dm_cmplx
   end interface

   interface elsi_pexsi_to_siesta_dm
      module procedure elsi_pexsi_to_siesta_dm_real
      module procedure elsi_pexsi_to_siesta_dm_cmplx
   end interface

   interface elsi_siesta_to_blacs_hs
      module procedure elsi_siesta_to_blacs_hs_real
      module procedure elsi_siesta_to_blacs_hs_cmplx
   end interface

   interface elsi_siesta_to_ntpoly_hs
      module procedure elsi_siesta_to_ntpoly_hs_real
      module procedure elsi_siesta_to_ntpoly_hs_cmplx
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

   interface elsi_sips_to_generic_dm
      module procedure elsi_sips_to_generic_dm_real
      module procedure elsi_sips_to_generic_dm_cmplx
   end interface

   interface elsi_sips_to_ntpoly_hs
      module procedure elsi_sips_to_ntpoly_hs_real
      module procedure elsi_sips_to_ntpoly_hs_cmplx
   end interface

   interface elsi_sips_to_siesta_dm
      module procedure elsi_sips_to_siesta_dm_real
      module procedure elsi_sips_to_siesta_dm_cmplx
   end interface

contains

!>
!! Get local number of nonzero elements in matrices in 1D block CSC format
!! converted from 2D block-cyclic dense format.
!!
subroutine elsi_blacs_to_pexsi_hs_dim_real(ph,bh,ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=*), parameter :: caller = "elsi_blacs_to_pexsi_hs_dim_real"

   call elsi_allocate(bh,dest,ph%pexsi_np_per_pole,"dest",caller)
   call elsi_allocate(bh,nnz,ph%pexsi_np_per_pole,"nnz",caller)

   if(.not. ph%unit_ovlp) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/ph%pexsi_np_per_pole)
         i_proc = min(i_proc,ph%pexsi_np_per_pole-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            end if
         end do
      end do
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/ph%pexsi_np_per_pole)
         i_proc = min(i_proc,ph%pexsi_np_per_pole-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            end if
         end do
      end do
   end if

   call MPI_Allreduce(dest,nnz,ph%pexsi_np_per_pole,mpi_integer4,mpi_sum,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   bh%nnz_l_sp = nnz(ph%pexsi_my_pcol+1)
   bh%nnz_l_sp1 = bh%nnz_l_sp

   call MPI_Allreduce(bh%nnz_l_sp,bh%nnz_g,1,mpi_integer4,mpi_sum,&
        ph%pexsi_comm_intra_pole,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! Get local number of nonzero elements in matrices in 1D block CSC format
!! converted from 2D block-cyclic dense format.
!!
subroutine elsi_blacs_to_pexsi_hs_dim_cmplx(ph,bh,ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=*), parameter :: caller = "elsi_blacs_to_pexsi_hs_dim_cmplx"

   call elsi_allocate(bh,dest,ph%pexsi_np_per_pole,"dest",caller)
   call elsi_allocate(bh,nnz,ph%pexsi_np_per_pole,"nnz",caller)

   if(.not. ph%unit_ovlp) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/ph%pexsi_np_per_pole)
         i_proc = min(i_proc,ph%pexsi_np_per_pole-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            end if
         end do
      end do
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/ph%pexsi_np_per_pole)
         i_proc = min(i_proc,ph%pexsi_np_per_pole-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            end if
         end do
      end do
   end if

   call MPI_Allreduce(dest,nnz,ph%pexsi_np_per_pole,mpi_integer4,mpi_sum,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   bh%nnz_l_sp = nnz(ph%pexsi_my_pcol+1)
   bh%nnz_l_sp1 = bh%nnz_l_sp

   call MPI_Allreduce(bh%nnz_l_sp,bh%nnz_g,1,mpi_integer4,mpi_sum,&
        ph%pexsi_comm_intra_pole,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 2D block-cyclic dense
!! format to 1D block CSC format.
!!
subroutine elsi_blacs_to_pexsi_hs_real(ph,bh,ham_den,ovlp_den,ham_sp,ovlp_sp,&
   row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: ovlp_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: ham_sp(bh%nnz_l_sp)
   real(kind=r8), intent(inout) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(inout) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(inout) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: n_group
   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: p_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: g_col
   integer(kind=i4) :: g_row
   integer(kind=i4) :: d1
   integer(kind=i4) :: d2
   integer(kind=i4) :: d11
   integer(kind=i4) :: d21
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: h_val_send(:)
   real(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: h_val_recv(:)
   real(kind=r8), allocatable :: s_val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of columns
   integer(kind=i4), allocatable :: perm(:)
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id

   character(len=*), parameter :: caller = "elsi_blacs_to_pexsi_hs_real"

   call elsi_get_time(t0)

   ! Compute destination of columns
   n_group = bh%n_procs/ph%pexsi_np_per_pole
   d1 = ph%n_basis/ph%pexsi_np_per_pole
   d2 = ph%n_basis-(ph%pexsi_np_per_pole-1)*d1
   d11 = max(d1/n_group,1)
   d21 = max(d2/n_group,1)

   call elsi_allocate(bh,dest,bh%n_lcol,"dest",caller)

   do i_col = 1,bh%n_lcol
      call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

      p_col = (g_col-1)/d1
      p_col = min(p_col,ph%pexsi_np_per_pole-1)

      if(p_col < ph%pexsi_np_per_pole-1) then
         p_row = mod(g_col-1,d1)/d11
      else
         p_row = (g_col-(ph%pexsi_np_per_pole-1)*d1-1)/d21
      end if

      p_row = min(p_row,n_group-1)

      dest(i_col) = p_row+p_col*n_group
   end do

   if(ph%first_blacs_to_pexsi) then
      if(.not. ph%unit_ovlp) then
         call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,ovlp_den,bh%nnz_l)

         call elsi_allocate(bh,s_val_send,bh%nnz_l,"s_val_send",caller)
      else
         call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,ham_den,bh%nnz_l)
      end if
   end if

   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l,"h_val_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_val = 0

   if(.not. ph%unit_ovlp) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               row_send(i_val) = g_row
               col_send(i_val) = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               if(ph%first_blacs_to_pexsi) then
                  s_val_send(i_val) = ovlp_den(i_row,i_col)
               end if

               ! Set send_count
               send_count(dest(i_col)+1) = send_count(dest(i_col)+1)+1
            end if
         end do
      end do
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               row_send(i_val) = g_row
               col_send(i_val) = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               ! Set send_count
               send_count(dest(i_col)+1) = send_count(dest(i_col)+1)+1
            end if
         end do
      end do
   end if

   call elsi_deallocate(bh,dest,"dest")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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
   if(ph%first_blacs_to_pexsi .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_recv,nnz_l_aux,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,s_val_recv,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   end if

   call elsi_allocate(bh,gid,nnz_l_aux,"gid",caller)
   call elsi_allocate(bh,perm,nnz_l_aux,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_recv,kind=i8)

   ! Sort
   call elsi_heapsort(nnz_l_aux,gid,perm)
   call elsi_permute(nnz_l_aux,perm,h_val_recv)
   call elsi_permute(nnz_l_aux,perm,row_recv)
   call elsi_permute(nnz_l_aux,perm,col_recv)

   if(ph%first_blacs_to_pexsi .and. .not. ph%unit_ovlp) then
      call elsi_permute(nnz_l_aux,perm,s_val_recv)
   end if

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")

   ! Set send_count, all data sent to 1st pole
   send_count = 0
   send_count(bh%myid/n_group+1) = nnz_l_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
   if(ph%first_blacs_to_pexsi) then
      ! Row id
      call MPI_Alltoallv(row_recv,send_count,send_displ,mpi_integer4,row_ind,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,row_recv,"row_recv")

      ! Column id
      call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)

      call MPI_Alltoallv(col_recv,send_count,send_displ,mpi_integer4,col_send,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,col_recv,"col_recv")

      if(.not. ph%unit_ovlp) then
         ! Overlap value
         call MPI_Alltoallv(s_val_recv,send_count,send_displ,mpi_real8,ovlp_sp,&
              recv_count,recv_displ,mpi_real8,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

         call elsi_deallocate(bh,s_val_recv,"s_val_recv")
      end if
   else
      call elsi_deallocate(bh,row_recv,"row_recv")
      call elsi_deallocate(bh,col_recv,"col_recv")
   end if

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv,send_count,send_displ,mpi_real8,ham_sp,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   if(ph%first_blacs_to_pexsi) then
      ! Only 1st pole computes row index and column pointer
      if(ph%pexsi_my_prow == 0) then
         col_ptr = 0
         col_ptr(bh%n_lcol_sp+1) = bh%nnz_l_sp+1

         do i_val = 1,bh%nnz_l_sp
            i_col = col_send(i_val)-d1*ph%pexsi_my_pcol
            col_ptr(i_col) = col_ptr(i_col)+1
         end do

         do i_col = bh%n_lcol_sp,1,-1
            col_ptr(i_col) = col_ptr(i_col+1)-col_ptr(i_col)
         end do
      end if

      call elsi_deallocate(bh,col_send,"col_send")

      call MPI_Bcast(row_ind,bh%nnz_l_sp,mpi_integer4,0,&
           ph%pexsi_comm_inter_pole,ierr)

      call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

      call MPI_Bcast(col_ptr,bh%n_lcol_sp+1,mpi_integer4,0,&
           ph%pexsi_comm_inter_pole,ierr)

      call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_blacs_to_pexsi = .false.

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 2D block-cyclic dense
!! format to 1D block CSC format.
!!
subroutine elsi_blacs_to_pexsi_hs_cmplx(ph,bh,ham_den,ovlp_den,ham_sp,ovlp_sp,&
   row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: ovlp_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: ham_sp(bh%nnz_l_sp)
   complex(kind=r8), intent(inout) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(inout) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(inout) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: n_group
   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: p_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: g_col
   integer(kind=i4) :: g_row
   integer(kind=i4) :: d1
   integer(kind=i4) :: d2
   integer(kind=i4) :: d11
   integer(kind=i4) :: d21
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

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
   integer(kind=i4), allocatable :: dest(:) ! Destination of columns
   integer(kind=i4), allocatable :: perm(:)
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id

   character(len=*), parameter :: caller = "elsi_blacs_to_pexsi_hs_cmplx"

   call elsi_get_time(t0)

   ! Compute destination of columns
   n_group = bh%n_procs/ph%pexsi_np_per_pole
   d1 = ph%n_basis/ph%pexsi_np_per_pole
   d2 = ph%n_basis-(ph%pexsi_np_per_pole-1)*d1
   d11 = max(d1/n_group,1)
   d21 = max(d2/n_group,1)

   call elsi_allocate(bh,dest,bh%n_lcol,"dest",caller)

   do i_col = 1,bh%n_lcol
      call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

      p_col = (g_col-1)/d1
      p_col = min(p_col,ph%pexsi_np_per_pole-1)

      if(p_col < ph%pexsi_np_per_pole-1) then
         p_row = mod(g_col-1,d1)/d11
      else
         p_row = (g_col-(ph%pexsi_np_per_pole-1)*d1-1)/d21
      end if

      p_row = min(p_row,n_group-1)

      dest(i_col) = p_row+p_col*n_group
   end do

   if(ph%first_blacs_to_pexsi) then
      if(.not. ph%unit_ovlp) then
         call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,ovlp_den,bh%nnz_l)

         call elsi_allocate(bh,s_val_send,bh%nnz_l,"s_val_send",caller)
      else
         call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,ham_den,bh%nnz_l)
      end if
   end if

   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l,"h_val_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_val = 0

   if(.not. ph%unit_ovlp) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               row_send(i_val) = g_row
               col_send(i_val) = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               if(ph%first_blacs_to_pexsi) then
                  s_val_send(i_val) = ovlp_den(i_row,i_col)
               end if

               ! Set send_count
               send_count(dest(i_col)+1) = send_count(dest(i_col)+1)+1
            end if
         end do
      end do
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               row_send(i_val) = g_row
               col_send(i_val) = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               ! Set send_count
               send_count(dest(i_col)+1) = send_count(dest(i_col)+1)+1
            end if
         end do
      end do
   end if

   call elsi_deallocate(bh,dest,"dest")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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
   if(ph%first_blacs_to_pexsi .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_recv,nnz_l_aux,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
           s_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   end if

   call elsi_allocate(bh,gid,nnz_l_aux,"gid",caller)
   call elsi_allocate(bh,perm,nnz_l_aux,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_recv,kind=i8)

   ! Sort
   call elsi_heapsort(nnz_l_aux,gid,perm)
   call elsi_permute(nnz_l_aux,perm,h_val_recv)
   call elsi_permute(nnz_l_aux,perm,row_recv)
   call elsi_permute(nnz_l_aux,perm,col_recv)

   if(ph%first_blacs_to_pexsi .and. .not. ph%unit_ovlp) then
      call elsi_permute(nnz_l_aux,perm,s_val_recv)
   end if

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")

   ! Set send_count, all data sent to 1st pole
   send_count = 0
   send_count(bh%myid/n_group+1) = nnz_l_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
   if(ph%first_blacs_to_pexsi) then
      ! Row id
      call MPI_Alltoallv(row_recv,send_count,send_displ,mpi_integer4,row_ind,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,row_recv,"row_recv")

      ! Column id
      call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)

      call MPI_Alltoallv(col_recv,send_count,send_displ,mpi_integer4,col_send,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,col_recv,"col_recv")

      if(.not. ph%unit_ovlp) then
         ! Overlap value
         call MPI_Alltoallv(s_val_recv,send_count,send_displ,mpi_complex16,&
              ovlp_sp,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

         call elsi_deallocate(bh,s_val_recv,"s_val_recv")
      end if
   else
      call elsi_deallocate(bh,row_recv,"row_recv")
      call elsi_deallocate(bh,col_recv,"col_recv")
   end if

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv,send_count,send_displ,mpi_complex16,ham_sp,&
        recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   if(ph%first_blacs_to_pexsi) then
      ! Only 1st pole computes row index and column pointer
      if(ph%pexsi_my_prow == 0) then
         col_ptr = 0
         col_ptr(bh%n_lcol_sp+1) = bh%nnz_l_sp+1

         do i_val = 1,bh%nnz_l_sp
            i_col = col_send(i_val)-d1*ph%pexsi_my_pcol
            col_ptr(i_col) = col_ptr(i_col)+1
         end do

         do i_col = bh%n_lcol_sp,1,-1
            col_ptr(i_col) = col_ptr(i_col+1)-col_ptr(i_col)
         end do
      end if

      call elsi_deallocate(bh,col_send,"col_send")

      call MPI_Bcast(row_ind,bh%nnz_l_sp,mpi_integer4,0,&
           ph%pexsi_comm_inter_pole,ierr)

      call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

      call MPI_Bcast(col_ptr,bh%n_lcol_sp+1,mpi_integer4,0,&
           ph%pexsi_comm_inter_pole,ierr)

      call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_blacs_to_pexsi = .false.

end subroutine

!>
!! Convert density matrix computed by PEXSI, stored in 1D block CSC format to
!! 2D block-cyclic dense format.
!!
subroutine elsi_pexsi_to_blacs_dm_real(ph,bh,dm_sp,row_ind,col_ptr,dm_den)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   real(kind=r8), intent(out) :: dm_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: n_group
   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: p_row
   integer(kind=i4) :: row0
   integer(kind=i4) :: row1
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_pexsi_to_blacs_dm_real"

   call elsi_get_time(t0)

   n_group = bh%n_procs/ph%pexsi_np_per_pole
   row0 = ph%pexsi_my_prow*(ph%n_basis/n_group)+1
   row1 = (ph%pexsi_my_prow+1)*(ph%n_basis/n_group)

   if(ph%pexsi_my_prow == n_group-1) then
      row1 = ph%n_basis
   end if

   nnz_l_aux = 0

   do i_val = 1,bh%nnz_l_sp
      if(row_ind(i_val) >= row0 .and. row_ind(i_val) <= row1) then
         nnz_l_aux = nnz_l_aux+1
      end if
   end do

   call elsi_allocate(bh,dest,nnz_l_aux,"dest",caller)
   call elsi_allocate(bh,perm,nnz_l_aux,"perm",caller)
   call elsi_allocate(bh,val_send,nnz_l_aux,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_aux,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_aux,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_col = 0
   j_val = 0

   do i_val = 1,bh%nnz_l_sp
      do while(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp)
         i_col = i_col+1
      end do

      i_row = row_ind(i_val)

      if(row_ind(i_val) >= row0 .and. row_ind(i_val) <= row1) then
         j_val = j_val+1

         ! Compute global id
         row_send(j_val) = i_row
         col_send(j_val) = i_col+ph%pexsi_my_pcol&
            *(ph%n_basis/ph%pexsi_np_per_pole)
         val_send(j_val) = dm_sp(i_val)

         ! Compute destination
         p_row = mod((row_send(j_val)-1)/bh%blk,bh%n_prow)
         p_col = mod((col_send(j_val)-1)/bh%blk,bh%n_pcol)
         dest(j_val) = p_col+p_row*bh%n_pcol

         ! Set send_count
         send_count(dest(j_val)+1) = send_count(dest(j_val)+1)+1
      end if
   end do

   ! Sort
   call elsi_heapsort(nnz_l_aux,dest,perm)
   call elsi_permute(nnz_l_aux,perm,val_send)
   call elsi_permute(nnz_l_aux,perm,row_send)
   call elsi_permute(nnz_l_aux,perm,col_send)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Density matrix value
   call elsi_allocate(bh,val_recv,bh%nnz_l,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_real8,val_recv,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
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
   end do

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix computed by PEXSI, stored in 1D block CSC format to 2D
!! block-cyclic dense format.
!!
subroutine elsi_pexsi_to_blacs_dm_cmplx(ph,bh,dm_sp,row_ind,col_ptr,dm_den)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   complex(kind=r8), intent(out) :: dm_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: n_group
   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: p_row
   integer(kind=i4) :: row0
   integer(kind=i4) :: row1
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

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
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_pexsi_to_blacs_dm_cmplx"

   call elsi_get_time(t0)

   n_group = bh%n_procs/ph%pexsi_np_per_pole
   row0 = ph%pexsi_my_prow*(ph%n_basis/n_group)+1
   row1 = (ph%pexsi_my_prow+1)*(ph%n_basis/n_group)

   if(ph%pexsi_my_prow == n_group-1) then
      row1 = ph%n_basis
   end if

   nnz_l_aux = 0

   do i_val = 1,bh%nnz_l_sp
      if(row_ind(i_val) >= row0 .and. row_ind(i_val) <= row1) then
         nnz_l_aux = nnz_l_aux+1
      end if
   end do

   call elsi_allocate(bh,dest,nnz_l_aux,"dest",caller)
   call elsi_allocate(bh,perm,nnz_l_aux,"perm",caller)
   call elsi_allocate(bh,val_send,nnz_l_aux,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_aux,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_aux,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_col = 0
   j_val = 0

   do i_val = 1,bh%nnz_l_sp
      do while(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp)
         i_col = i_col+1
      end do

      i_row = row_ind(i_val)

      if(row_ind(i_val) >= row0 .and. row_ind(i_val) <= row1) then
         j_val = j_val+1

         ! Compute global id
         row_send(j_val) = i_row
         col_send(j_val) = i_col+ph%pexsi_my_pcol&
            *(ph%n_basis/ph%pexsi_np_per_pole)
         val_send(j_val) = dm_sp(i_val)

         ! Compute destination
         p_row = mod((row_send(j_val)-1)/bh%blk,bh%n_prow)
         p_col = mod((col_send(j_val)-1)/bh%blk,bh%n_pcol)
         dest(j_val) = p_col+p_row*bh%n_pcol

         ! Set send_count
         send_count(dest(j_val)+1) = send_count(dest(j_val)+1)+1
      end if
   end do

   ! Sort
   call elsi_heapsort(nnz_l_aux,dest,perm)
   call elsi_permute(nnz_l_aux,perm,val_send)
   call elsi_permute(nnz_l_aux,perm,row_send)
   call elsi_permute(nnz_l_aux,perm,col_send)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Density matrix value
   call elsi_allocate(bh,val_recv,bh%nnz_l,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_complex16,val_recv,&
        recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
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
   end do

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Get local number of nonzero elements in matrices in 1D block CSC format
!! converted from 2D block-cyclic dense format.
!!
subroutine elsi_blacs_to_sips_hs_dim_real(ph,bh,ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=*), parameter :: caller = "elsi_blacs_to_sips_hs_dim_real"

   call elsi_allocate(bh,dest,bh%n_procs,"dest",caller)
   call elsi_allocate(bh,nnz,bh%n_procs,"nnz",caller)

   if(.not. ph%unit_ovlp) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/bh%n_procs)
         i_proc = min(i_proc,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            end if
         end do
      end do
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/bh%n_procs)
         i_proc = min(i_proc,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            end if
         end do
      end do
   end if

   call MPI_Allreduce(dest,nnz,bh%n_procs,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   bh%nnz_l_sp = nnz(bh%myid+1)
   bh%nnz_l_sp1 = bh%nnz_l_sp

   call MPI_Allreduce(bh%nnz_l_sp,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! Get local number of nonzero elements in matrices in 1D block CSC format
!! converted from 2D block-cyclic dense format.
!!
subroutine elsi_blacs_to_sips_hs_dim_cmplx(ph,bh,ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=*), parameter :: caller = "elsi_blacs_to_sips_hs_dim_cmplx"

   call elsi_allocate(bh,dest,bh%n_procs,"dest",caller)
   call elsi_allocate(bh,nnz,bh%n_procs,"nnz",caller)

   if(.not. ph%unit_ovlp) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/bh%n_procs)
         i_proc = min(i_proc,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            end if
         end do
      end do
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         i_proc = (g_col-1)/(ph%n_basis/bh%n_procs)
         i_proc = min(i_proc,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               dest(i_proc+1) = dest(i_proc+1)+1
            end if
         end do
      end do
   end if

   call MPI_Allreduce(dest,nnz,bh%n_procs,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   bh%nnz_l_sp = nnz(bh%myid+1)
   bh%nnz_l_sp1 = bh%nnz_l_sp

   call MPI_Allreduce(bh%nnz_l_sp,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 2D block-cyclic dense
!! format to 1D block CSC format.
!!
subroutine elsi_blacs_to_sips_hs_real(ph,bh,ham_den,ovlp_den,ham_sp,ovlp_sp,&
   row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: ovlp_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: ham_sp(bh%nnz_l_sp)
   real(kind=r8), intent(inout) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(inout) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(inout) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: g_col
   integer(kind=i4) :: g_row
   integer(kind=i4) :: dest ! Destination of an element
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: h_val_send(:)
   real(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_blacs_to_sips_hs_real"

   call elsi_get_time(t0)

   if(ph%first_blacs_to_sips) then
      if(.not. ph%unit_ovlp) then
         call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,ovlp_den,bh%nnz_l)

         call elsi_allocate(bh,s_val_send,bh%nnz_l,"s_val_send",caller)
      else
         call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,ham_den,bh%nnz_l)
      end if
   end if

   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l,"h_val_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_val = 0

   if(.not. ph%unit_ovlp) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = (g_col-1)/(ph%n_basis/bh%n_procs)
         dest = min(dest,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val) = g_row
               col_send(i_val) = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               if(ph%first_blacs_to_sips) then
                  s_val_send(i_val) = ovlp_den(i_row,i_col)
               end if

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            end if
         end do
      end do
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = (g_col-1)/(ph%n_basis/bh%n_procs)
         dest = min(dest,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val) = g_row
               col_send(i_val) = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            end if
         end do
      end do
   end if

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
   ! Row id
   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_ind,&
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
   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_real8,ham_sp,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%first_blacs_to_sips .and. .not. ph%unit_ovlp) then
      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,ovlp_sp,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_ind,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,gid,perm)
   call elsi_permute(bh%nnz_l_sp,perm,ham_sp)
   call elsi_permute(bh%nnz_l_sp,perm,row_ind)
   call elsi_permute(bh%nnz_l_sp,perm,col_recv)

   if(ph%first_blacs_to_sips .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp,perm,ovlp_sp)
   end if

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")

   ! Compute row index and column pointer
   if(ph%first_blacs_to_sips) then
      col_ptr = 0
      col_ptr(bh%n_lcol_sp+1) = bh%nnz_l_sp+1

      do i_val = 1,bh%nnz_l_sp
         i_col = col_recv(i_val)-(ph%n_basis/bh%n_procs)*bh%myid
         col_ptr(i_col) = col_ptr(i_col)+1
      end do

      do i_col = bh%n_lcol_sp,1,-1
         col_ptr(i_col) = col_ptr(i_col+1)-col_ptr(i_col)
      end do
   end if

   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_blacs_to_sips = .false.

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 2D block-cyclic dense
!! format to 1D block CSC format.
!!
subroutine elsi_blacs_to_sips_hs_cmplx(ph,bh,ham_den,ovlp_den,ham_sp,ovlp_sp,&
   row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: ovlp_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: ham_sp(bh%nnz_l_sp)
   complex(kind=r8), intent(inout) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(inout) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(inout) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: g_col
   integer(kind=i4) :: g_row
   integer(kind=i4) :: dest ! Destination of an element
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: h_val_send(:)
   complex(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_blacs_to_sips_hs_cmplx"

   call elsi_get_time(t0)

   if(ph%first_blacs_to_sips) then
      if(.not. ph%unit_ovlp) then
         call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,ovlp_den,bh%nnz_l)

         call elsi_allocate(bh,s_val_send,bh%nnz_l,"s_val_send",caller)
      else
         call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,ham_den,bh%nnz_l)
      end if
   end if

   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l,"h_val_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_val = 0

   if(.not. ph%unit_ovlp) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = (g_col-1)/(ph%n_basis/bh%n_procs)
         dest = min(dest,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val) = g_row
               col_send(i_val) = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               if(ph%first_blacs_to_sips) then
                  s_val_send(i_val) = ovlp_den(i_row,i_col)
               end if

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            end if
         end do
      end do
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         ! Compute destination
         dest = (g_col-1)/(ph%n_basis/bh%n_procs)
         dest = min(dest,bh%n_procs-1)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               i_val = i_val+1

               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               ! Pack global id and data into bufs
               row_send(i_val) = g_row
               col_send(i_val) = g_col
               h_val_send(i_val) = ham_den(i_row,i_col)

               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            end if
         end do
      end do
   end if

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
   ! Row id
   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_ind,&
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
   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_complex16,ham_sp,&
        recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%first_blacs_to_sips .and. .not. ph%unit_ovlp) then
      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
           ovlp_sp,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_ind,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,gid,perm)
   call elsi_permute(bh%nnz_l_sp,perm,ham_sp)
   call elsi_permute(bh%nnz_l_sp,perm,row_ind)
   call elsi_permute(bh%nnz_l_sp,perm,col_recv)

   if(ph%first_blacs_to_sips .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp,perm,ovlp_sp)
   end if

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")

   ! Compute row index and column pointer
   if(ph%first_blacs_to_sips) then
      col_ptr = 0
      col_ptr(bh%n_lcol_sp+1) = bh%nnz_l_sp+1

      do i_val = 1,bh%nnz_l_sp
         i_col = col_recv(i_val)-(ph%n_basis/bh%n_procs)*bh%myid
         col_ptr(i_col) = col_ptr(i_col)+1
      end do

      do i_col = bh%n_lcol_sp,1,-1
         col_ptr(i_col) = col_ptr(i_col+1)-col_ptr(i_col)
      end do
   end if

   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_blacs_to_sips = .false.

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 1D block CSC format to 2D
!! block-cyclic dense format.
!!
subroutine elsi_sips_to_blacs_hs_real(ph,bh,ham_sp,ovlp_sp,row_ind,col_ptr,&
   ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   real(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   real(kind=r8), intent(out) :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: p_row
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: h_val_send(:)
   real(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: h_val_recv(:)
   real(kind=r8), allocatable :: s_val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_sips_to_blacs_hs_real"

   call elsi_get_time(t0)

   if(ph%first_sips_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp,"s_val_send",caller)
   end if

   call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_col = 0

   do i_val = 1,bh%nnz_l_sp
      do while(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp)
         i_col = i_col+1
      end do

      i_row = row_ind(i_val)

      ! Compute global id
      row_send(i_val) = i_row
      col_send(i_val) = i_col+bh%myid*(ph%n_basis/bh%n_procs)
      h_val_send(i_val) = ham_sp(i_val)

      if(ph%first_sips_to_blacs .and. .not. ph%unit_ovlp) then
         s_val_send(i_val) = ovlp_sp(i_val)
      end if

      ! Compute destination
      p_row = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
      p_col = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
      dest(i_val) = p_col+p_row*bh%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   end do

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,dest,perm)
   call elsi_permute(bh%nnz_l_sp,perm,h_val_send)
   call elsi_permute(bh%nnz_l_sp,perm,row_send)
   call elsi_permute(bh%nnz_l_sp,perm,col_send)

   if(ph%first_sips_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp,perm,s_val_send)
   end if

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send_recv and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Hamiltonian value
   call elsi_allocate(bh,h_val_recv,bh%nnz_l,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_real8,h_val_recv,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%first_sips_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_recv,bh%nnz_l,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,s_val_recv,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   ! Unpack matrix
   if(ph%first_sips_to_blacs .and. .not. ph%unit_ovlp) then
      ham_den = 0.0_r8
      ovlp_den = 0.0_r8

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
         ovlp_den(l_row,l_col) = s_val_recv(i_val)
      end do

      call elsi_deallocate(bh,s_val_recv,"s_val_recv")
   else
      ham_den = 0.0_r8

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
      end do
   end if

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_sips_to_blacs = .false.

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 1D block CSC format to 2D
!! block-cyclic dense format.
!!
subroutine elsi_sips_to_blacs_hs_cmplx(ph,bh,ham_sp,ovlp_sp,row_ind,col_ptr,&
   ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   complex(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   complex(kind=r8), intent(out) :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: p_row
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

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
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_sips_to_blacs_hs_cmplx"

   call elsi_get_time(t0)

   if(ph%first_sips_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp,"s_val_send",caller)
   end if

   call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_col = 0

   do i_val = 1,bh%nnz_l_sp
      do while(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp)
         i_col = i_col+1
      end do

      i_row = row_ind(i_val)

      ! Compute global id
      row_send(i_val) = i_row
      col_send(i_val) = i_col+bh%myid*(ph%n_basis/bh%n_procs)
      h_val_send(i_val) = ham_sp(i_val)

      if(ph%first_sips_to_blacs .and. .not. ph%unit_ovlp) then
         s_val_send(i_val) = ovlp_sp(i_val)
      end if

      ! Compute destination
      p_row = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
      p_col = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
      dest(i_val) = p_col+p_row*bh%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   end do

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,dest,perm)
   call elsi_permute(bh%nnz_l_sp,perm,h_val_send)
   call elsi_permute(bh%nnz_l_sp,perm,row_send)
   call elsi_permute(bh%nnz_l_sp,perm,col_send)

   if(ph%first_sips_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp,perm,s_val_send)
   end if

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Hamiltonian value
   call elsi_allocate(bh,h_val_recv,bh%nnz_l,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_complex16,&
        h_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%first_sips_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_recv,bh%nnz_l,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
           s_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   ! Unpack matrix
   if(ph%first_sips_to_blacs .and. .not. ph%unit_ovlp) then
      ham_den = (0.0_r8,0.0_r8)
      ovlp_den = (0.0_r8,0.0_r8)

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
         ovlp_den(l_row,l_col) = s_val_recv(i_val)
      end do

      call elsi_deallocate(bh,s_val_recv,"s_val_recv")
   else
      ham_den = (0.0_r8,0.0_r8)

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
      end do
   end if

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_sips_to_blacs = .false.

end subroutine

!>
!! Convert eigenvectors stored in 1D block dense format to 2D block-cyclic dense
!! format.
!!
subroutine elsi_sips_to_blacs_ev_real(ph,bh,evec_sips,evec)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: evec_sips(bh%n_lcol_sp1,ph%n_states)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: p_row
   integer(kind=i4) :: nnz_before
   integer(kind=i4) :: nnz_after
   integer(kind=i4) :: n_lrow_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_sips_to_blacs_ev_real"

   call elsi_get_time(t0)

   n_lrow_aux = ph%n_basis/bh%n_procs
   nnz_before = 0

   do i_col = 1,ph%n_states
      do i_row = 1,bh%n_lcol_sp
         if(abs(evec_sips(i_row,i_col)) > bh%def0) then
            nnz_before = nnz_before+1
         end if
      end do
   end do

   call elsi_allocate(bh,dest,nnz_before,"dest",caller)
   call elsi_allocate(bh,perm,nnz_before,"perm",caller)
   call elsi_allocate(bh,val_send,nnz_before,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_before,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_before,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_val = 0

   do i_col = 1,ph%n_states
      do i_row = 1,bh%n_lcol_sp
         if(abs(evec_sips(i_row,i_col)) > bh%def0) then
            i_val = i_val+1

            ! Compute global id
            col_send(i_val) = i_col
            row_send(i_val) = bh%myid*n_lrow_aux+i_row
            val_send(i_val) = evec_sips(i_row,i_col)

            ! Compute destination
            p_row = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
            p_col = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
            dest(i_val) = p_col+p_row*bh%n_pcol

            ! Set send_count
            send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
         end if
      end do
   end do

   ! Sort
   call elsi_heapsort(nnz_before,dest,perm)
   call elsi_permute(nnz_before,perm,val_send)
   call elsi_permute(nnz_before,perm,row_send)
   call elsi_permute(nnz_before,perm,col_send)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   nnz_after = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Eigenvector value
   call elsi_allocate(bh,val_recv,nnz_after,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_real8,val_recv,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
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
   end do

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix in 2D block-cyclic dense format to 1D block CSC
!! format.
!!
subroutine elsi_blacs_to_sips_dm_real(ph,bh,dm_den,dm_sp,row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: dm_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)

   character(len=*), parameter :: caller = "elsi_blacs_to_sips_dm_real"

   call elsi_get_time(t0)

   call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,dm_den,bh%nnz_l)

   call elsi_allocate(bh,val_send,bh%nnz_l,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_val = 0

   do i_col = 1,bh%n_lcol
      do i_row = 1,bh%n_lrow
         if(abs(dm_den(i_row,i_col)) > bh%def0) then
            i_val = i_val+1

            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,row_send(i_val))
            call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,col_send(i_val))

            ! Compute destination
            dest = (col_send(i_val)-1)/(ph%n_basis/bh%n_procs)
            dest = min(dest,bh%n_procs-1)

            ! Pack data
            val_send(i_val) = dm_den(i_row,i_col)

            ! Set send_count
            send_count(dest+1) = send_count(dest+1)+1
         end if
      end do
   end do

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   dm_sp = 0.0_r8

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      l_col = col_recv(i_val)-(ph%n_basis/bh%n_procs)*bh%myid
      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_sp(j_val) = val_recv(i_val)
         end if
      end do
   end do

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix in 2D block-cyclic dense format to 1D block CSC
!! format.
!!
subroutine elsi_blacs_to_sips_dm_cmplx(ph,bh,dm_den,dm_sp,row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: dm_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

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

   character(len=*), parameter :: caller = "elsi_blacs_to_sips_dm_cmplx"

   call elsi_get_time(t0)

   call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,dm_den,bh%nnz_l)

   call elsi_allocate(bh,val_send,bh%nnz_l,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_val = 0

   do i_col = 1,bh%n_lcol
      do i_row = 1,bh%n_lrow
         if(abs(dm_den(i_row,i_col)) > bh%def0) then
            i_val = i_val+1

            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,row_send(i_val))
            call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,col_send(i_val))

            ! Compute destination
            dest = (col_send(i_val)-1)/(ph%n_basis/bh%n_procs)
            dest = min(dest,bh%n_procs-1)

            ! Pack data
            val_send(i_val) = dm_den(i_row,i_col)

            ! Set send_count
            send_count(dest+1) = send_count(dest+1)+1
         end if
      end do
   end do

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   dm_sp = (0.0_r8,0.0_r8)

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      l_col = col_recv(i_val)-(ph%n_basis/bh%n_procs)*bh%myid
      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_sp(j_val) = val_recv(i_val)
         end if
      end do
   end do

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 1D block-cyclic CSC format
!! to 2D block-cyclic dense format.
!!
subroutine elsi_siesta_to_blacs_hs_real(ph,bh,ham_sp,ovlp_sp,row_ind,col_ptr,&
   ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   real(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   real(kind=r8), intent(out) :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: p_row
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: h_val_send(:)
   real(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: h_val_recv(:)
   real(kind=r8), allocatable :: s_val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_siesta_to_blacs_hs_real"

   call elsi_get_time(t0)

   if(ph%first_siesta_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp,"s_val_send",caller)
   end if

   call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_col = 0

   do i_val = 1,bh%nnz_l_sp
      do while(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp2)
         i_col = i_col+1
      end do

      i_row = row_ind(i_val)

      ! Compute global id
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,col_send(i_val))

      row_send(i_val) = i_row
      h_val_send(i_val) = ham_sp(i_val)

      if(ph%first_siesta_to_blacs .and. .not. ph%unit_ovlp) then
         s_val_send(i_val) = ovlp_sp(i_val)
      end if

      ! Compute destination
      p_row = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
      p_col = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
      dest(i_val) = p_col+p_row*bh%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   end do

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,dest,perm)
   call elsi_permute(bh%nnz_l_sp,perm,h_val_send)
   call elsi_permute(bh%nnz_l_sp,perm,row_send)
   call elsi_permute(bh%nnz_l_sp,perm,col_send)

   if(ph%first_siesta_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp,perm,s_val_send)
   end if

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Hamiltonian value
   call elsi_allocate(bh,h_val_recv,bh%nnz_l,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_real8,h_val_recv,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%first_siesta_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_recv,bh%nnz_l,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,s_val_recv,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   ! Unpack matrix
   if(ph%first_siesta_to_blacs .and. .not. ph%unit_ovlp) then
      ham_den = 0.0_r8
      ovlp_den = 0.0_r8

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
         ovlp_den(l_row,l_col) = s_val_recv(i_val)
      end do

      call elsi_deallocate(bh,s_val_recv,"s_val_recv")
   else
      ham_den = 0.0_r8

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
      end do
   end if

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_siesta_to_blacs = .false.

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 1D block-cyclic CSC format
!! to 2D block-cyclic dense format.
!!
subroutine elsi_siesta_to_blacs_hs_cmplx(ph,bh,ham_sp,ovlp_sp,row_ind,col_ptr,&
   ham_den,ovlp_den)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   complex(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   complex(kind=r8), intent(out) :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: p_row
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

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
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_siesta_to_blacs_hs_cmplx"

   call elsi_get_time(t0)

   if(ph%first_siesta_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp,"s_val_send",caller)
   end if

   call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_col = 0

   do i_val = 1,bh%nnz_l_sp
      do while(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp2)
         i_col = i_col+1
      end do

      i_row = row_ind(i_val)

      ! Compute global id
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,col_send(i_val))

      row_send(i_val) = i_row
      h_val_send(i_val) = ham_sp(i_val)

      if(ph%first_siesta_to_blacs .and. .not. ph%unit_ovlp) then
         s_val_send(i_val) = ovlp_sp(i_val)
      end if

      ! Compute destination
      p_row = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
      p_col = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
      dest(i_val) = p_col+p_row*bh%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   end do

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,dest,perm)
   call elsi_permute(bh%nnz_l_sp,perm,h_val_send)
   call elsi_permute(bh%nnz_l_sp,perm,row_send)
   call elsi_permute(bh%nnz_l_sp,perm,col_send)

   if(ph%first_siesta_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp,perm,s_val_send)
   end if

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Hamiltonian value
   call elsi_allocate(bh,h_val_recv,bh%nnz_l,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_complex16,&
        h_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%first_siesta_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_recv,bh%nnz_l,"s_val_recv",caller)

      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
           s_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   ! Unpack matrix
   if(ph%first_siesta_to_blacs .and. .not. ph%unit_ovlp) then
      ham_den = (0.0_r8,0.0_r8)
      ovlp_den = (0.0_r8,0.0_r8)

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
         ovlp_den(l_row,l_col) = s_val_recv(i_val)
      end do

      call elsi_deallocate(bh,s_val_recv,"s_val_recv")
   else
      ham_den = (0.0_r8,0.0_r8)

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
      end do
   end if

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_siesta_to_blacs = .false.

end subroutine

!>
!! Convert density matrix in 2D block-cyclic dense format to 1D block-cyclic CSC
!! format.
!!
subroutine elsi_blacs_to_siesta_dm_real(bh,dm_den,dm_sp,row_ind,col_ptr)

   implicit none

   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: dm_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_blacs_to_siesta_dm_real"

   call elsi_get_time(t0)

   call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,dm_den,bh%nnz_l)

   call elsi_allocate(bh,dest,bh%nnz_l,"dest",caller)
   call elsi_allocate(bh,perm,bh%nnz_l,"perm",caller)
   call elsi_allocate(bh,val_send,bh%nnz_l,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

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
         end if
      end do
   end do

   ! Sort
   call elsi_heapsort(bh%nnz_l,dest,perm)
   call elsi_permute(bh%nnz_l,perm,val_send)
   call elsi_permute(bh%nnz_l,perm,row_send)
   call elsi_permute(bh%nnz_l,perm,col_send)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   dm_sp = 0.0_r8

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      call elsi_get_lid(bh%n_procs,bh%blk_sp2,col_recv(i_val),l_col)

      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_sp(j_val) = val_recv(i_val)
         end if
      end do
   end do

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix in 2D block-cyclic dense format to 1D block-cyclic CSC
!! format.
!!
subroutine elsi_blacs_to_siesta_dm_cmplx(bh,dm_den,dm_sp,row_ind,col_ptr)

   implicit none

   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: dm_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

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
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_blacs_to_siesta_dm_cmplx"

   call elsi_get_time(t0)

   call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,dm_den,bh%nnz_l)

   call elsi_allocate(bh,dest,bh%nnz_l,"dest",caller)
   call elsi_allocate(bh,perm,bh%nnz_l,"perm",caller)
   call elsi_allocate(bh,val_send,bh%nnz_l,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

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
         end if
      end do
   end do

   ! Sort
   call elsi_heapsort(bh%nnz_l,dest,perm)
   call elsi_permute(bh%nnz_l,perm,val_send)
   call elsi_permute(bh%nnz_l,perm,row_send)
   call elsi_permute(bh%nnz_l,perm,col_send)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   dm_sp = (0.0_r8,0.0_r8)

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      call elsi_get_lid(bh%n_procs,bh%blk_sp2,col_recv(i_val),l_col)

      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_sp(j_val) = val_recv(i_val)
         end if
      end do
   end do

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Get local number of nonzero elements in matrices in 1D block CSC format
!! converted from 1D block-cyclic dense format.
!!
subroutine elsi_siesta_to_pexsi_hs_dim(ph,bh,col_ptr2)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4), intent(in) :: col_ptr2(bh%n_lcol_sp2+1)

   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=*), parameter :: caller = "elsi_siesta_to_pexsi_hs_dim"

   call elsi_allocate(bh,dest,ph%pexsi_np_per_pole,"dest",caller)
   call elsi_allocate(bh,nnz,ph%pexsi_np_per_pole,"nnz",caller)

   do i_col = 1,bh%n_lcol_sp2
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,g_col)

      i_proc = (g_col-1)/(ph%n_basis/ph%pexsi_np_per_pole)
      i_proc = min(i_proc,ph%pexsi_np_per_pole-1)

      dest(i_proc+1) = dest(i_proc+1)+col_ptr2(i_col+1)-col_ptr2(i_col)
   end do

   call MPI_Allreduce(dest,nnz,ph%pexsi_np_per_pole,mpi_integer4,mpi_sum,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   bh%nnz_l_sp1 = nnz(ph%pexsi_my_pcol+1)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 1D block-cyclic CSC format
!! to 1D block CSC format.
!!
subroutine elsi_siesta_to_pexsi_hs_real(ph,bh,ham_sp2,ovlp_sp2,row_ind2,&
   col_ptr2,ham_sp1,ovlp_sp1,row_ind1,col_ptr1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: ham_sp2(bh%nnz_l_sp2)
   real(kind=r8), intent(in) :: ovlp_sp2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: col_ptr2(bh%n_lcol_sp2+1)
   real(kind=r8), intent(out) :: ham_sp1(bh%nnz_l_sp1)
   real(kind=r8), intent(inout) :: ovlp_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: col_ptr1(bh%n_lcol_sp1+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: n_lcol_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: h_val_send(:)
   real(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_siesta_to_pexsi_hs_real"

   call elsi_get_time(t0)

   n_lcol_aux = ph%n_basis/ph%pexsi_np_per_pole

   if(ph%first_siesta_to_pexsi .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp2,"s_val_send",caller)
   end if

   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp2,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp2,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp2,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_col = 0

   do i_val = 1,bh%nnz_l_sp2
      do while(i_val == col_ptr2(i_col+1) .and. i_col /= bh%n_lcol_sp2)
         i_col = i_col+1
      end do

      i_row = row_ind2(i_val)

      ! Compute global id
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,col_send(i_val))

      row_send(i_val) = i_row
      h_val_send(i_val) = ham_sp2(i_val)

      if(ph%first_siesta_to_pexsi .and. .not. ph%unit_ovlp) then
         s_val_send(i_val) = ovlp_sp2(i_val)
      end if

      ! Compute destination
      dest = (col_send(i_val)-1)/n_lcol_aux
      dest = min(dest,ph%pexsi_np_per_pole-1)

      ! Set send_count
      send_count(dest+1) = send_count(dest+1)+1
   end do

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
   ! Row index
   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_ind1,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp1,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_real8,ham_sp1,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%first_siesta_to_pexsi .and. .not. ph%unit_ovlp) then
      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,ovlp_sp1,&
           recv_count,recv_displ,mpi_real8,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp1,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp1,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_ind1,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp1,gid,perm)
   call elsi_permute(bh%nnz_l_sp1,perm,ham_sp1)
   call elsi_permute(bh%nnz_l_sp1,perm,row_ind1)
   call elsi_permute(bh%nnz_l_sp1,perm,col_recv)

   if(ph%first_siesta_to_pexsi .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp1,perm,ovlp_sp1)
   end if

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")

   if(ph%first_siesta_to_pexsi) then
      ! Only 1st pole computes row index and column pointer
      if(ph%pexsi_my_prow == 0) then
         col_ptr1 = 0
         col_ptr1(bh%n_lcol_sp1+1) = bh%nnz_l_sp1+1

         do i_val = 1,bh%nnz_l_sp1
            i_col = col_recv(i_val)-n_lcol_aux*ph%pexsi_my_pcol
            col_ptr1(i_col) = col_ptr1(i_col)+1
         end do

         do i_col = bh%n_lcol_sp1,1,-1
            col_ptr1(i_col) = col_ptr1(i_col+1)-col_ptr1(i_col)
         end do
      end if
   end if

   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_siesta_to_pexsi = .false.

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 1D block-cyclic CSC format
!! to 1D block CSC format.
!!
subroutine elsi_siesta_to_pexsi_hs_cmplx(ph,bh,ham_sp2,ovlp_sp2,row_ind2,&
   col_ptr2,ham_sp1,ovlp_sp1,row_ind1,col_ptr1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: ham_sp2(bh%nnz_l_sp2)
   complex(kind=r8), intent(in) :: ovlp_sp2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: col_ptr2(bh%n_lcol_sp2+1)
   complex(kind=r8), intent(out) :: ham_sp1(bh%nnz_l_sp1)
   complex(kind=r8), intent(inout) :: ovlp_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: col_ptr1(bh%n_lcol_sp1+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: n_lcol_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: h_val_send(:)
   complex(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_siesta_to_pexsi_hs_cmplx"

   call elsi_get_time(t0)

   n_lcol_aux = ph%n_basis/ph%pexsi_np_per_pole

   if(ph%first_siesta_to_pexsi .and. .not. ph%unit_ovlp) then
      call elsi_allocate(bh,s_val_send,bh%nnz_l_sp2,"s_val_send",caller)
   end if

   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp2,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp2,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp2,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_col = 0

   do i_val = 1,bh%nnz_l_sp2
      do while(i_val == col_ptr2(i_col+1) .and. i_col /= bh%n_lcol_sp2)
         i_col = i_col+1
      end do

      i_row = row_ind2(i_val)

      ! Compute global id
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,col_send(i_val))

      row_send(i_val) = i_row
      h_val_send(i_val) = ham_sp2(i_val)

      if(ph%first_siesta_to_pexsi .and. .not. ph%unit_ovlp) then
         s_val_send(i_val) = ovlp_sp2(i_val)
      end if

      ! Compute destination
      dest = (col_send(i_val)-1)/n_lcol_aux
      dest = min(dest,ph%pexsi_np_per_pole-1)

      ! Set send_count
      send_count(dest+1) = send_count(dest+1)+1
   end do

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
   ! Row index
   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_ind1,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp1,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_complex16,ham_sp1,&
        recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   ! Overlap value
   if(ph%first_siesta_to_pexsi .and. .not. ph%unit_ovlp) then
      call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
           ovlp_sp1,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,s_val_send,"s_val_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp1,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp1,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_ind1,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp1,gid,perm)
   call elsi_permute(bh%nnz_l_sp1,perm,ham_sp1)
   call elsi_permute(bh%nnz_l_sp1,perm,row_ind1)
   call elsi_permute(bh%nnz_l_sp1,perm,col_recv)

   if(ph%first_siesta_to_pexsi .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp1,perm,ovlp_sp1)
   end if

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")

   if(ph%first_siesta_to_pexsi) then
      ! Only 1st pole computes row index and column pointer
      if(ph%pexsi_my_prow == 0) then
         col_ptr1 = 0
         col_ptr1(bh%n_lcol_sp1+1) = bh%nnz_l_sp1+1

         do i_val = 1,bh%nnz_l_sp1
            i_col = col_recv(i_val)-n_lcol_aux*ph%pexsi_my_pcol
            col_ptr1(i_col) = col_ptr1(i_col)+1
         end do

         do i_col = bh%n_lcol_sp1,1,-1
            col_ptr1(i_col) = col_ptr1(i_col+1)-col_ptr1(i_col)
         end do
      end if
   end if

   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_siesta_to_pexsi = .false.

end subroutine

!>
!! Convert density matrix computed by PEXSI, stored in 1D block CSC format to 1D
!! block-cyclic CSC format.
!!
subroutine elsi_pexsi_to_siesta_dm_real(ph,bh,dm_sp1,row_ind1,col_ptr1,dm_sp2,&
   row_ind2,col_ptr2)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: dm_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: col_ptr1(bh%n_lcol_sp1+1)
   real(kind=r8), intent(out) :: dm_sp2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: col_ptr2(bh%n_lcol_sp2+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_pexsi_to_siesta_dm_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,val_send,bh%nnz_l_sp1,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp1,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp1,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(ph%pexsi_my_prow == 0) then
      call elsi_allocate(bh,dest,bh%nnz_l_sp1,"dest",caller)
      call elsi_allocate(bh,perm,bh%nnz_l_sp1,"perm",caller)

      i_col = 0

      do i_val = 1,bh%nnz_l_sp1
         do while(i_val == col_ptr1(i_col+1) .and. i_col /= bh%n_lcol_sp1)
            i_col = i_col+1
         end do

         i_row = row_ind1(i_val)

         ! Compute global id
         row_send(i_val) = i_row
         col_send(i_val) = i_col+bh%myid*(ph%n_basis/ph%pexsi_np_per_pole)
         val_send(i_val) = dm_sp1(i_val)

         ! Compute destination
         dest(i_val) = mod((col_send(i_val)-1)/bh%blk_sp2,bh%n_procs)

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      end do

      ! Sort
      call elsi_heapsort(bh%nnz_l_sp1,dest,perm)
      call elsi_permute(bh%nnz_l_sp1,perm,val_send)
      call elsi_permute(bh%nnz_l_sp1,perm,row_send)
      call elsi_permute(bh%nnz_l_sp1,perm,col_send)

      call elsi_deallocate(bh,dest,"dest")
      call elsi_deallocate(bh,perm,"perm")
   end if

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Density matrix value
   call elsi_allocate(bh,val_recv,bh%nnz_l_sp2,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_real8,val_recv,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   dm_sp2 = 0.0_r8

   ! Unpack matrix
   do i_val = 1,bh%nnz_l_sp2
      call elsi_get_lid(bh%n_procs,bh%blk_sp2,col_recv(i_val),l_col)

      l_row = row_recv(i_val)

      do j_val = col_ptr2(l_col),col_ptr2(l_col+1)-1
         if(row_ind2(j_val) == l_row) then
            dm_sp2(j_val) = val_recv(i_val)
         end if
      end do
   end do

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix computed by PEXSI, stored in 1D block CSC format to 1D
!! block-cyclic CSC format.
!!
subroutine elsi_pexsi_to_siesta_dm_cmplx(ph,bh,dm_sp1,row_ind1,col_ptr1,dm_sp2,&
   row_ind2,col_ptr2)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: dm_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: col_ptr1(bh%n_lcol_sp1+1)
   complex(kind=r8), intent(out) :: dm_sp2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: col_ptr2(bh%n_lcol_sp2+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

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
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_pexsi_to_siesta_dm_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,val_send,bh%nnz_l_sp1,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp1,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp1,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(ph%pexsi_my_prow == 0) then
      call elsi_allocate(bh,dest,bh%nnz_l_sp1,"dest",caller)
      call elsi_allocate(bh,perm,bh%nnz_l_sp1,"perm",caller)

      i_col = 0

      do i_val = 1,bh%nnz_l_sp1
         do while(i_val == col_ptr1(i_col+1) .and. i_col /= bh%n_lcol_sp1)
            i_col = i_col+1
         end do

         i_row = row_ind1(i_val)

         ! Compute global id
         row_send(i_val) = i_row
         col_send(i_val) = i_col+bh%myid*(ph%n_basis/ph%pexsi_np_per_pole)
         val_send(i_val) = dm_sp1(i_val)

         ! Compute destination
         dest(i_val) = mod((col_send(i_val)-1)/bh%blk_sp2,bh%n_procs)

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      end do

      ! Sort
      call elsi_heapsort(bh%nnz_l_sp1,dest,perm)
      call elsi_permute(bh%nnz_l_sp1,perm,val_send)
      call elsi_permute(bh%nnz_l_sp1,perm,row_send)
      call elsi_permute(bh%nnz_l_sp1,perm,col_send)

      call elsi_deallocate(bh,dest,"dest")
      call elsi_deallocate(bh,perm,"perm")
   end if

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Density matrix value
   call elsi_allocate(bh,val_recv,bh%nnz_l_sp2,"val_recv",caller)

   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_complex16,val_recv,&
        recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   dm_sp2 = (0.0_r8,0.0_r8)

   ! Unpack matrix
   do i_val = 1,bh%nnz_l_sp2
      call elsi_get_lid(bh%n_procs,bh%blk_sp2,col_recv(i_val),l_col)

      l_row = row_recv(i_val)

      do j_val = col_ptr2(l_col),col_ptr2(l_col+1)-1
         if(row_ind2(j_val) == l_row) then
            dm_sp2(j_val) = val_recv(i_val)
         end if
      end do
   end do

   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")
   call elsi_deallocate(bh,val_recv,"val_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix computed by SLEPc-SIPs, stored in 1D block CSC format
!! to 2D block-cyclic dense format.
!!
subroutine elsi_sips_to_blacs_dm_real(ph,bh,dm_sp,row_ind,col_ptr,dm_den)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   real(kind=r8), intent(out) :: dm_den(bh%n_lrow,bh%n_lcol)

   character(len=*), parameter :: caller = "elsi_sips_to_blacs_dm_real"

   ph%pexsi_my_prow = 0
   ph%pexsi_my_pcol = bh%myid
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_pexsi_to_blacs_dm_real(ph,bh,dm_sp,row_ind,col_ptr,dm_den)

end subroutine

!>
!! Convert density matrix computed by SLEPc-SIPs, stored in 1D block CSC format
!! to 2D block-cyclic dense format.
!!
subroutine elsi_sips_to_blacs_dm_cmplx(ph,bh,dm_sp,row_ind,col_ptr,dm_den)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   complex(kind=r8), intent(out) :: dm_den(bh%n_lrow,bh%n_lcol)

   character(len=*), parameter :: caller = "elsi_sips_to_blacs_dm_cmplx"

   ph%pexsi_my_prow = 0
   ph%pexsi_my_pcol = bh%myid
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_pexsi_to_blacs_dm_cmplx(ph,bh,dm_sp,row_ind,col_ptr,dm_den)

end subroutine

!>
!! Convert density matrix computed by SLPEc-SIPs, stored in 1D block CSC format
!! to 1D block-cyclic CSC format.
!!
subroutine elsi_sips_to_siesta_dm_real(ph,bh,dm_sp1,row_ind1,col_ptr1,dm_sp2,&
   row_ind2,col_ptr2)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: dm_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: col_ptr1(bh%n_lcol_sp1+1)
   real(kind=r8), intent(out) :: dm_sp2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: col_ptr2(bh%n_lcol_sp2+1)

   character(len=*), parameter :: caller = "elsi_sips_to_siesta_dm_real"

   ph%pexsi_my_prow = 0
   ph%pexsi_my_pcol = bh%myid
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_pexsi_to_siesta_dm_real(ph,bh,dm_sp1,row_ind1,col_ptr1,dm_sp2,&
        row_ind2,col_ptr2)

end subroutine

!>
!! Convert density matrix computed by SLEPc-SIPs, stored in 1D block CSC format
!! to 1D block-cyclic CSC format.
!!
subroutine elsi_sips_to_siesta_dm_cmplx(ph,bh,dm_sp1,row_ind1,col_ptr1,dm_sp2,&
   row_ind2,col_ptr2)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: dm_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: col_ptr1(bh%n_lcol_sp1+1)
   complex(kind=r8), intent(out) :: dm_sp2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: col_ptr2(bh%n_lcol_sp2+1)

   character(len=*), parameter :: caller = "elsi_sips_to_siesta_dm_cmplx"

   ph%pexsi_my_prow = 0
   ph%pexsi_my_pcol = bh%myid
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_pexsi_to_siesta_dm_cmplx(ph,bh,dm_sp1,row_ind1,col_ptr1,dm_sp2,&
        row_ind2,col_ptr2)

end subroutine

!>
!! Get local number of nonzero elements in matrices in 1D block CSC format
!! converted from 1D block-cyclic dense format.
!!
subroutine elsi_siesta_to_sips_hs_dim(ph,bh,col_ptr2)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4), intent(in) :: col_ptr2(bh%n_lcol_sp2+1)

   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=*), parameter :: caller = "elsi_siesta_to_sips_hs_dim"

   call elsi_allocate(bh,dest,bh%n_procs,"dest",caller)
   call elsi_allocate(bh,nnz,bh%n_procs,"nnz",caller)

   do i_col = 1,bh%n_lcol_sp2
      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,g_col)

      i_proc = (g_col-1)/(ph%n_basis/bh%n_procs)
      i_proc = min(i_proc,bh%n_procs-1)

      dest(i_proc+1) = dest(i_proc+1)+col_ptr2(i_col+1)-col_ptr2(i_col)
   end do

   call MPI_Allreduce(dest,nnz,bh%n_procs,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   bh%nnz_l_sp1 = nnz(bh%myid+1)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 1D block-cyclic CSC format
!! to 1D block CSC format.
!!
subroutine elsi_siesta_to_sips_hs_real(ph,bh,ham_sp2,ovlp_sp2,row_ind2,&
   col_ptr2,ham_sp1,ovlp_sp1,row_ind1,col_ptr1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: ham_sp2(bh%nnz_l_sp2)
   real(kind=r8), intent(in) :: ovlp_sp2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: col_ptr2(bh%nnz_l_sp2)
   real(kind=r8), intent(out) :: ham_sp1(bh%nnz_l_sp1)
   real(kind=r8), intent(inout) :: ovlp_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: col_ptr1(bh%nnz_l_sp1)

   character(len=*), parameter :: caller = "elsi_siesta_to_sips_hs_real"

   ph%pexsi_my_prow = 0
   ph%pexsi_my_pcol = bh%myid
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_siesta_to_pexsi_hs_real(ph,bh,ham_sp2,ovlp_sp2,row_ind2,col_ptr2,&
        ham_sp1,ovlp_sp1,row_ind1,col_ptr1)

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in 1D block-cyclic CSC format
!! to 1D block CSC format.
!!
subroutine elsi_siesta_to_sips_hs_cmplx(ph,bh,ham_sp2,ovlp_sp2,row_ind2,&
   col_ptr2,ham_sp1,ovlp_sp1,row_ind1,col_ptr1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: ham_sp2(bh%nnz_l_sp2)
   complex(kind=r8), intent(in) :: ovlp_sp2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: row_ind2(bh%nnz_l_sp2)
   integer(kind=i4), intent(in) :: col_ptr2(bh%nnz_l_sp2)
   complex(kind=r8), intent(out) :: ham_sp1(bh%nnz_l_sp1)
   complex(kind=r8), intent(inout) :: ovlp_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: col_ptr1(bh%nnz_l_sp1)

   character(len=*), parameter :: caller = "elsi_siesta_to_sips_hs_cmplx"

   ph%pexsi_my_prow = 0
   ph%pexsi_my_pcol = bh%myid
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_siesta_to_pexsi_hs_cmplx(ph,bh,ham_sp2,ovlp_sp2,row_ind2,col_ptr2,&
        ham_sp1,ovlp_sp1,row_ind1,col_ptr1)

end subroutine

!>
!! Construct Halmitonian and overlep matrices in NTPoly format from matrices
!! stored in 2D block-cyclic dense format.
!!
subroutine elsi_blacs_to_ntpoly_hs_real(ph,bh,ham_den,ovlp_den,ham_nt,ovlp_nt)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(in) :: ovlp_den(bh%n_lrow,bh%n_lcol)
   type(Matrix_ps), intent(inout) :: ham_nt
   type(Matrix_ps), intent(inout) :: ovlp_nt

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_row
   integer(kind=i4) :: g_col
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   type(TripletList_r) :: ham_list
   type(TripletList_r) :: ovlp_list
   type(Triplet_r) :: coo

   character(len=*), parameter :: caller = "elsi_blacs_to_ntpoly_hs_real"

   call elsi_get_time(t0)

   if(ph%first_blacs_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call ConstructEmptyMatrix(ovlp_nt,ph%n_basis,ph%nt_pgrid,.false.)
         call ConstructTripletList(ovlp_list)
      end if

      call ConstructEmptyMatrix(ham_nt,ph%n_basis,ph%nt_pgrid,.false.)
   end if

   call ConstructTripletList(ham_list)

   if(.not. ph%unit_ovlp) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               coo%point_value = ham_den(i_row,i_col)
               coo%index_column = g_col
               coo%index_row = g_row

               call AppendToTripletList(ham_list,coo)

               if(ph%first_blacs_to_ntpoly) then
                  coo%point_value = ovlp_den(i_row,i_col)

                  call AppendToTripletList(ovlp_list,coo)
               end if
            end if
         end do
      end do
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               coo%point_value = ham_den(i_row,i_col)
               coo%index_column = g_col
               coo%index_row = g_row

               call AppendToTripletList(ham_list,coo)
            end if
         end do
      end do
   end if

   if(ph%first_blacs_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call MPI_Allreduce(ovlp_list%CurrentSize,bh%nnz_g,1,mpi_integer4,&
              mpi_sum,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

         call FillMatrixFromTripletList(ovlp_nt,ovlp_list)
         call DestructTripletList(ovlp_list)
      else
         call MPI_Allreduce(ham_list%CurrentSize,bh%nnz_g,1,mpi_integer4,&
              mpi_sum,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
      end if
   end if

   call FillMatrixFromTripletList(ham_nt,ham_list)
   call DestructTripletList(ham_list)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_blacs_to_ntpoly = .false.

end subroutine

!>
!! Construct Halmitonian and overlep matrices in NTPoly format from matrices
!! stored in 2D block-cyclic dense format.
!!
subroutine elsi_blacs_to_ntpoly_hs_cmplx(ph,bh,ham_den,ovlp_den,ham_nt,ovlp_nt)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(in) :: ovlp_den(bh%n_lrow,bh%n_lcol)
   type(Matrix_ps), intent(inout) :: ham_nt
   type(Matrix_ps), intent(inout) :: ovlp_nt

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_row
   integer(kind=i4) :: g_col
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   type(TripletList_c) :: ham_list
   type(TripletList_c) :: ovlp_list
   type(Triplet_c) :: coo

   character(len=*), parameter :: caller = "elsi_blacs_to_ntpoly_hs_cmplx"

   call elsi_get_time(t0)

   if(ph%first_blacs_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call ConstructEmptyMatrix(ovlp_nt,ph%n_basis,ph%nt_pgrid,.true.)
         call ConstructTripletList(ovlp_list)
      end if

      call ConstructEmptyMatrix(ham_nt,ph%n_basis,ph%nt_pgrid,.true.)
   end if

   call ConstructTripletList(ham_list)

   if(.not. ph%unit_ovlp) then
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         do i_row = 1,bh%n_lrow
            if(abs(ovlp_den(i_row,i_col)) > bh%def0) then
               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               coo%point_value = ham_den(i_row,i_col)
               coo%index_column = g_col
               coo%index_row = g_row

               call AppendToTripletList(ham_list,coo)

               if(ph%first_blacs_to_ntpoly) then
                  coo%point_value = ovlp_den(i_row,i_col)

                  call AppendToTripletList(ovlp_list,coo)
               end if
            end if
         end do
      end do
   else
      do i_col = 1,bh%n_lcol
         call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,g_col)

         do i_row = 1,bh%n_lrow
            if(abs(ham_den(i_row,i_col)) > bh%def0) then
               call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,g_row)

               coo%point_value = ham_den(i_row,i_col)
               coo%index_column = g_col
               coo%index_row = g_row

               call AppendToTripletList(ham_list,coo)
            end if
         end do
      end do
   end if

   if(ph%first_blacs_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call MPI_Allreduce(ovlp_list%CurrentSize,bh%nnz_g,1,mpi_integer4,&
              mpi_sum,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

         call FillMatrixFromTripletList(ovlp_nt,ovlp_list)
         call DestructTripletList(ovlp_list)
      else
         call MPI_Allreduce(ham_list%CurrentSize,bh%nnz_g,1,mpi_integer4,&
              mpi_sum,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
      end if
   end if

   call FillMatrixFromTripletList(ham_nt,ham_list)
   call DestructTripletList(ham_list)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_blacs_to_ntpoly = .false.

end subroutine

!>
!! Convert density matrix computed by NTPoly to 2D block-cyclic dense format.
!!
subroutine elsi_ntpoly_to_blacs_dm_real(ph,bh,dm_nt,dm_den)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(Matrix_ps), intent(inout) :: dm_nt
   real(kind=r8), intent(out) :: dm_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: p_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: nnz_l_nt
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   type(TripletList_r) :: dm_list

   character(len=*), parameter :: caller = "elsi_ntpoly_to_blacs_dm_real"

   call elsi_get_time(t0)

   call GetMatrixTripletList(dm_nt,dm_list)

   nnz_l_nt = dm_list%CurrentSize

   call elsi_allocate(bh,val_send,nnz_l_nt,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_nt,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_nt,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(bh%myid < bh%n_procs/ph%nt_n_layers) then
      call elsi_allocate(bh,dest,nnz_l_nt,"dest",caller)
      call elsi_allocate(bh,perm,nnz_l_nt,"perm",caller)

      do i_val = 1,nnz_l_nt
         ! Compute global id
         row_send(i_val) = dm_list%data(i_val)%index_row
         col_send(i_val) = dm_list%data(i_val)%index_column
         val_send(i_val) = dm_list%data(i_val)%point_value

         ! Compute destination
         p_row = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
         p_col = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
         dest(i_val) = p_col+p_row*bh%n_pcol

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      end do

      ! Sort
      call elsi_heapsort(nnz_l_nt,dest,perm)
      call elsi_permute(nnz_l_nt,perm,val_send)
      call elsi_permute(nnz_l_nt,perm,row_send)
      call elsi_permute(nnz_l_nt,perm,col_send)

      call elsi_deallocate(bh,dest,"dest")
      call elsi_deallocate(bh,perm,"perm")
   end if

   call DestructTripletList(dm_list)

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   dm_den = 0.0_r8

   ! Unpack density matrix
   do i_val = 1,nnz_l_aux
      ! Compute local 2d id
      call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
      call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

      ! Put value to correct position
      dm_den(l_row,l_col) = val_recv(i_val)
   end do

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix computed by NTPoly to 2D block-cyclic dense format.
!!
subroutine elsi_ntpoly_to_blacs_dm_cmplx(ph,bh,dm_nt,dm_den)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(Matrix_ps), intent(inout) :: dm_nt
   complex(kind=r8), intent(out) :: dm_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: p_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: nnz_l_nt
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

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
   integer(kind=i4), allocatable :: perm(:)

   type(TripletList_c) :: dm_list

   character(len=*), parameter :: caller = "elsi_ntpoly_to_blacs_dm_cmplx"

   call elsi_get_time(t0)

   call GetMatrixTripletList(dm_nt,dm_list)

   nnz_l_nt = dm_list%CurrentSize

   call elsi_allocate(bh,val_send,nnz_l_nt,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_nt,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_nt,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(bh%myid < bh%n_procs/ph%nt_n_layers) then
      call elsi_allocate(bh,dest,nnz_l_nt,"dest",caller)
      call elsi_allocate(bh,perm,nnz_l_nt,"perm",caller)

      do i_val = 1,nnz_l_nt
         ! Compute global id
         row_send(i_val) = dm_list%data(i_val)%index_row
         col_send(i_val) = dm_list%data(i_val)%index_column
         val_send(i_val) = dm_list%data(i_val)%point_value

         ! Compute destination
         p_row = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
         p_col = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
         dest(i_val) = p_col+p_row*bh%n_pcol

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      end do

      ! Sort
      call elsi_heapsort(nnz_l_nt,dest,perm)
      call elsi_permute(nnz_l_nt,perm,val_send)
      call elsi_permute(nnz_l_nt,perm,row_send)
      call elsi_permute(nnz_l_nt,perm,col_send)

      call elsi_deallocate(bh,dest,"dest")
      call elsi_deallocate(bh,perm,"perm")
   end if

   call DestructTripletList(dm_list)

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   dm_den = (0.0_r8,0.0_r8)

   ! Unpack density matrix
   do i_val = 1,nnz_l_aux
      ! Compute local 2d id
      call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
      call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

      ! Put value to correct position
      dm_den(l_row,l_col) = val_recv(i_val)
   end do

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Construct Halmitonian and overlep matrices in NTPoly format from matrices
!! stored in 1D block CSC format.
!!
subroutine elsi_sips_to_ntpoly_hs_real(ph,bh,ham_sp,ovlp_sp,row_ind,col_ptr,&
   ham_nt,ovlp_nt)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   real(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   type(Matrix_ps), intent(inout) :: ham_nt
   type(Matrix_ps), intent(inout) :: ovlp_nt

   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_col
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   type(TripletList_r) :: ham_list
   type(TripletList_r) :: ovlp_list
   type(Triplet_r) :: coo

   character(len=*), parameter :: caller = "elsi_sips_to_ntpoly_hs_real"

   call elsi_get_time(t0)

   if(ph%first_sips_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call ConstructEmptyMatrix(ovlp_nt,ph%n_basis,ph%nt_pgrid,.false.)
         call ConstructTripletList(ovlp_list)
      end if

      call ConstructEmptyMatrix(ham_nt,ph%n_basis,ph%nt_pgrid,.false.)
   end if

   call ConstructTripletList(ham_list)

   i_col = 0

   do i_val = 1,bh%nnz_l_sp
      do while(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp)
         i_col = i_col+1
      end do

      coo%point_value = ham_sp(i_val)
      coo%index_row = row_ind(i_val)
      coo%index_column = i_col+bh%myid*(ph%n_basis/bh%n_procs)

      call AppendToTripletList(ham_list,coo)

      if(ph%first_sips_to_ntpoly .and. .not. ph%unit_ovlp) then
         coo%point_value = ovlp_sp(i_val)

         call AppendToTripletList(ovlp_list,coo)
      end if
   end do

   if(ph%first_sips_to_ntpoly .and. .not. ph%unit_ovlp) then
      call FillMatrixFromTripletList(ovlp_nt,ovlp_list)
      call DestructTripletList(ovlp_list)
   end if

   call FillMatrixFromTripletList(ham_nt,ham_list)
   call DestructTripletList(ham_list)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_sips_to_ntpoly = .false.

end subroutine

!>
!! Construct Halmitonian and overlep matrices in NTPoly format from matrices
!! stored in 1D block CSC format.
!!
subroutine elsi_sips_to_ntpoly_hs_cmplx(ph,bh,ham_sp,ovlp_sp,row_ind,col_ptr,&
   ham_nt,ovlp_nt)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   complex(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   type(Matrix_ps), intent(inout) :: ham_nt
   type(Matrix_ps), intent(inout) :: ovlp_nt

   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_col
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   type(TripletList_c) :: ham_list
   type(TripletList_c) :: ovlp_list
   type(Triplet_c) :: coo

   character(len=*), parameter :: caller = "elsi_sips_to_ntpoly_hs_cmplx"

   call elsi_get_time(t0)

   if(ph%first_sips_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call ConstructEmptyMatrix(ovlp_nt,ph%n_basis,ph%nt_pgrid,.true.)
         call ConstructTripletList(ovlp_list)
      end if

      call ConstructEmptyMatrix(ham_nt,ph%n_basis,ph%nt_pgrid,.true.)
   end if

   call ConstructTripletList(ham_list)

   i_col = 0

   do i_val = 1,bh%nnz_l_sp
      do while(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp)
         i_col = i_col+1
      end do

      coo%point_value = ham_sp(i_val)
      coo%index_row = row_ind(i_val)
      coo%index_column = i_col+bh%myid*(ph%n_basis/bh%n_procs)

      call AppendToTripletList(ham_list,coo)

      if(ph%first_sips_to_ntpoly .and. .not. ph%unit_ovlp) then
         coo%point_value = ovlp_sp(i_val)

         call AppendToTripletList(ovlp_list,coo)
      end if
   end do

   if(ph%first_sips_to_ntpoly .and. .not. ph%unit_ovlp) then
      call FillMatrixFromTripletList(ovlp_nt,ovlp_list)
      call DestructTripletList(ovlp_list)
   end if

   call FillMatrixFromTripletList(ham_nt,ham_list)
   call DestructTripletList(ham_list)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_sips_to_ntpoly = .false.

end subroutine

!>
!! Convert density matrix computed by NTPoly to 1D block CSC format.
!!
subroutine elsi_ntpoly_to_sips_dm_real(ph,bh,dm_nt,dm_sp,row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(Matrix_ps), intent(inout) :: dm_nt
   real(kind=r8), intent(out) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: nnz_l_nt
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)

   type(TripletList_r) :: dm_list

   character(len=*), parameter :: caller = "elsi_ntpoly_to_sips_dm_real"

   call elsi_get_time(t0)

   call GetMatrixTripletList(dm_nt,dm_list)

   nnz_l_nt = dm_list%CurrentSize

   call elsi_allocate(bh,val_send,nnz_l_nt,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_nt,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_nt,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(bh%myid < bh%n_procs/ph%nt_n_layers) then
      do i_val = 1,nnz_l_nt
         ! Compute global id
         row_send(i_val) = dm_list%data(i_val)%index_row
         col_send(i_val) = dm_list%data(i_val)%index_column
         val_send(i_val) = dm_list%data(i_val)%point_value

         ! Compute destination
         dest = (col_send(i_val)-1)/(ph%n_basis/bh%n_procs)
         dest = min(dest,bh%n_procs-1)

         ! Set send_count
         send_count(dest+1) = send_count(dest+1)+1
      end do
   end if

   call DestructTripletList(dm_list)

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   dm_sp = 0.0_r8

   ! Unpack density matrix
   do i_val = 1,nnz_l_aux
      l_col = col_recv(i_val)-(ph%n_basis/bh%n_procs)*bh%myid
      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_sp(j_val) = val_recv(i_val)
         end if
      end do
   end do

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix computed by NTPoly to 1D block CSC format.
!!
subroutine elsi_ntpoly_to_sips_dm_cmplx(ph,bh,dm_nt,dm_sp,row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(Matrix_ps), intent(inout) :: dm_nt
   complex(kind=r8), intent(out) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: nnz_l_nt
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

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

   type(TripletList_r) :: dm_list

   character(len=*), parameter :: caller = "elsi_ntpoly_to_sips_dm_complex"

   call elsi_get_time(t0)

   call GetMatrixTripletList(dm_nt,dm_list)

   nnz_l_nt = dm_list%CurrentSize

   call elsi_allocate(bh,val_send,nnz_l_nt,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_nt,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_nt,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(bh%myid < bh%n_procs/ph%nt_n_layers) then
      do i_val = 1,nnz_l_nt
         ! Compute global id
         row_send(i_val) = dm_list%data(i_val)%index_row
         col_send(i_val) = dm_list%data(i_val)%index_column
         val_send(i_val) = dm_list%data(i_val)%point_value

         ! Compute destination
         dest = (col_send(i_val)-1)/(ph%n_basis/bh%n_procs)
         dest = min(dest,bh%n_procs-1)

         ! Set send_count
         send_count(dest+1) = send_count(dest+1)+1
      end do
   end if

   call DestructTripletList(dm_list)

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   dm_sp = (0.0_r8,0.0_r8)

   ! Unpack density matrix
   do i_val = 1,nnz_l_aux
      l_col = col_recv(i_val)-(ph%n_basis/bh%n_procs)*bh%myid
      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_sp(j_val) = val_recv(i_val)
         end if
      end do
   end do

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Construct Halmitonian and overlep matrices in NTPoly format from matrices
!! stored in 1D block-cyclic CSC format.
!!
subroutine elsi_siesta_to_ntpoly_hs_real(ph,bh,ham_sp,ovlp_sp,row_ind,col_ptr,&
   ham_nt,ovlp_nt)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   real(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   type(Matrix_ps), intent(inout) :: ham_nt
   type(Matrix_ps), intent(inout) :: ovlp_nt

   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   type(TripletList_r) :: ham_list
   type(TripletList_r) :: ovlp_list
   type(Triplet_r) :: coo

   character(len=*), parameter :: caller = "elsi_siesta_to_ntpoly_hs_real"

   call elsi_get_time(t0)

   if(ph%first_siesta_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call ConstructEmptyMatrix(ovlp_nt,ph%n_basis,ph%nt_pgrid,.false.)
         call ConstructTripletList(ovlp_list)
      end if

      call ConstructEmptyMatrix(ham_nt,ph%n_basis,ph%nt_pgrid,.false.)
   end if

   call ConstructTripletList(ham_list)

   i_col = 0

   do i_val = 1,bh%nnz_l_sp
      do while(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp)
         i_col = i_col+1
      end do

      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,g_col)

      coo%point_value = ham_sp(i_val)
      coo%index_row = row_ind(i_val)
      coo%index_column = g_col

      call AppendToTripletList(ham_list,coo)

      if(ph%first_siesta_to_ntpoly .and. .not. ph%unit_ovlp) then
         coo%point_value = ovlp_sp(i_val)

         call AppendToTripletList(ovlp_list,coo)
      end if
   end do

   if(ph%first_siesta_to_ntpoly .and. .not. ph%unit_ovlp) then
      call FillMatrixFromTripletList(ovlp_nt,ovlp_list)
      call DestructTripletList(ovlp_list)
   end if

   call FillMatrixFromTripletList(ham_nt,ham_list)
   call DestructTripletList(ham_list)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_siesta_to_ntpoly = .false.

end subroutine

!>
!! Construct Halmitonian and overlep matrices in NTPoly format from matrices
!! stored in 1D block-cyclic CSC format.
!!
subroutine elsi_siesta_to_ntpoly_hs_cmplx(ph,bh,ham_sp,ovlp_sp,row_ind,col_ptr,&
   ham_nt,ovlp_nt)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)
   complex(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   complex(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   type(Matrix_ps), intent(inout) :: ham_nt
   type(Matrix_ps), intent(inout) :: ovlp_nt

   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_col
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   type(TripletList_c) :: ham_list
   type(TripletList_c) :: ovlp_list
   type(Triplet_c) :: coo

   character(len=*), parameter :: caller = "elsi_siesta_to_ntpoly_hs_cmplx"

   call elsi_get_time(t0)

   if(ph%first_siesta_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call ConstructEmptyMatrix(ovlp_nt,ph%n_basis,ph%nt_pgrid,.true.)
         call ConstructTripletList(ovlp_list)
      end if

      call ConstructEmptyMatrix(ham_nt,ph%n_basis,ph%nt_pgrid,.true.)
   end if

   call ConstructTripletList(ham_list)

   i_col = 0

   do i_val = 1,bh%nnz_l_sp
      do while(i_val == col_ptr(i_col+1) .and. i_col /= bh%n_lcol_sp)
         i_col = i_col+1
      end do

      call elsi_get_gid(bh%myid,bh%n_procs,bh%blk_sp2,i_col,g_col)

      coo%point_value = ham_sp(i_val)
      coo%index_row = row_ind(i_val)
      coo%index_column = g_col

      call AppendToTripletList(ham_list,coo)

      if(ph%first_siesta_to_ntpoly .and. .not. ph%unit_ovlp) then
         coo%point_value = ovlp_sp(i_val)

         call AppendToTripletList(ovlp_list,coo)
      end if
   end do

   if(ph%first_siesta_to_ntpoly .and. .not. ph%unit_ovlp) then
      call FillMatrixFromTripletList(ovlp_nt,ovlp_list)
      call DestructTripletList(ovlp_list)
   end if

   call FillMatrixFromTripletList(ham_nt,ham_list)
   call DestructTripletList(ham_list)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_siesta_to_ntpoly = .false.

end subroutine

!>
!! Convert density matrix computed by NTPoly to 1D block-cyclic CSC format.
!!
subroutine elsi_ntpoly_to_siesta_dm_real(ph,bh,dm_nt,dm_sp,row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(Matrix_ps), intent(inout) :: dm_nt
   real(kind=r8), intent(out) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: nnz_l_nt
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: val_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   type(TripletList_r) :: dm_list

   character(len=*), parameter :: caller = "elsi_ntpoly_to_siesta_dm_real"

   call elsi_get_time(t0)

   call GetMatrixTripletList(dm_nt,dm_list)

   nnz_l_nt = dm_list%CurrentSize

   call elsi_allocate(bh,val_send,nnz_l_nt,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_nt,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_nt,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(bh%myid < bh%n_procs/ph%nt_n_layers) then
      call elsi_allocate(bh,dest,nnz_l_nt,"dest",caller)
      call elsi_allocate(bh,perm,nnz_l_nt,"perm",caller)

      do i_val = 1,nnz_l_nt
         ! Compute global id
         row_send(i_val) = dm_list%data(i_val)%index_row
         col_send(i_val) = dm_list%data(i_val)%index_column
         val_send(i_val) = dm_list%data(i_val)%point_value

         ! Compute destination
         dest(i_val) = mod((col_send(i_val)-1)/bh%blk_sp2,bh%n_procs)

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      end do

      ! Sort
      call elsi_heapsort(nnz_l_nt,dest,perm)
      call elsi_permute(nnz_l_nt,perm,val_send)
      call elsi_permute(nnz_l_nt,perm,row_send)
      call elsi_permute(nnz_l_nt,perm,col_send)

      call elsi_deallocate(bh,dest,"dest")
      call elsi_deallocate(bh,perm,"perm")
   end if

   call DestructTripletList(dm_list)

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   dm_sp = 0.0_r8

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      call elsi_get_lid(bh%n_procs,bh%blk_sp2,col_recv(i_val),l_col)

      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_sp(j_val) = val_recv(i_val)
         end if
      end do
   end do

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix computed by NTPoly to 1D block-cyclic CSC format.
!!
subroutine elsi_ntpoly_to_siesta_dm_cmplx(ph,bh,dm_nt,dm_sp,row_ind,col_ptr)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(Matrix_ps), intent(inout) :: dm_nt
   complex(kind=r8), intent(out):: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp+1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: nnz_l_nt
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

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
   integer(kind=i4), allocatable :: perm(:)

   type(TripletList_c) :: dm_list

   character(len=*), parameter :: caller = "elsi_ntpoly_to_siesta_dm_cmplx"

   call elsi_get_time(t0)

   call GetMatrixTripletList(dm_nt,dm_list)

   nnz_l_nt = dm_list%CurrentSize

   call elsi_allocate(bh,val_send,nnz_l_nt,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_nt,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_nt,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(bh%myid < bh%n_procs/ph%nt_n_layers) then
      call elsi_allocate(bh,dest,nnz_l_nt,"dest",caller)
      call elsi_allocate(bh,perm,nnz_l_nt,"perm",caller)

      do i_val = 1,nnz_l_nt
         ! Compute global id
         row_send(i_val) = dm_list%data(i_val)%index_row
         col_send(i_val) = dm_list%data(i_val)%index_column
         val_send(i_val) = dm_list%data(i_val)%point_value

         ! Compute destination
         dest(i_val) = mod((col_send(i_val)-1)/bh%blk_sp2,bh%n_procs)

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      end do

      ! Sort
      call elsi_heapsort(nnz_l_nt,dest,perm)
      call elsi_permute(nnz_l_nt,perm,val_send)
      call elsi_permute(nnz_l_nt,perm,row_send)
      call elsi_permute(nnz_l_nt,perm,col_send)

      call elsi_deallocate(bh,dest,"dest")
      call elsi_deallocate(bh,perm,"perm")
   end if

   call DestructTripletList(dm_list)

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set local number of nonzero
   nnz_l_aux = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   dm_sp = (0.0_r8,0.0_r8)

   ! Unpack matrix
   do i_val = 1,nnz_l_aux
      call elsi_get_lid(bh%n_procs,bh%blk_sp2,col_recv(i_val),l_col)

      l_row = row_recv(i_val)

      do j_val = col_ptr(l_col),col_ptr(l_col+1)-1
         if(row_ind(j_val) == l_row) then
            dm_sp(j_val) = val_recv(i_val)
         end if
      end do
   end do

   call elsi_deallocate(bh,val_recv,"val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in generic COO format to 2D
!! block-cyclic dense format.
!!
subroutine elsi_generic_to_blacs_hs_real(ph,bh,ham_sp,ovlp_sp,row_ind,col_ind,&
   ham_den,ovlp_den,map_den)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   real(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ind(bh%nnz_l_sp)
   real(kind=r8), intent(out) :: ham_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp_den(bh%n_lrow,bh%n_lcol)
   integer(kind=i4), intent(inout) :: map_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: p_row
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: h_val_send(:)
   real(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: map_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   real(kind=r8), allocatable :: h_val_recv(:)
   real(kind=r8), allocatable :: s_val_recv(:)
   integer(kind=i4), allocatable :: map_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_generic_to_blacs_hs_real"

   call elsi_get_time(t0)

   if(ph%first_generic_to_blacs) then
      if(.not. ph%unit_ovlp) then
         call elsi_allocate(bh,s_val_send,bh%nnz_l_sp,"s_val_send",caller)
      end if

      call elsi_allocate(bh,map_send,bh%nnz_l_sp,"map_send",caller)

      map_send = bh%myid
   end if

   call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   row_send = row_ind
   col_send = col_ind
   h_val_send = ham_sp

   if(ph%first_generic_to_blacs .and. .not. ph%unit_ovlp) then
      s_val_send = ovlp_sp
   end if

   do i_val = 1,bh%nnz_l_sp
      ! Compute destination
      p_row = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
      p_col = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
      dest(i_val) = p_col+p_row*bh%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   end do

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,dest,perm)
   call elsi_permute(bh%nnz_l_sp,perm,h_val_send)
   call elsi_permute(bh%nnz_l_sp,perm,row_send)
   call elsi_permute(bh%nnz_l_sp,perm,col_send)

   if(ph%first_generic_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp,perm,s_val_send)
   end if

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Hamiltonian value
   call elsi_allocate(bh,h_val_recv,bh%nnz_l,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_real8,h_val_recv,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   if(ph%first_generic_to_blacs) then
      ! Overlap value
      if(.not. ph%unit_ovlp) then
         call elsi_allocate(bh,s_val_recv,bh%nnz_l,"s_val_recv",caller)

         call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,&
              s_val_recv,recv_count,recv_displ,mpi_real8,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

         call elsi_deallocate(bh,s_val_send,"s_val_send")
      end if

      call elsi_allocate(bh,map_recv,bh%nnz_l,"map_recv",caller)

      call MPI_Alltoallv(map_send,send_count,send_displ,mpi_integer4,map_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,map_send,"map_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   ! Unpack matrix
   if(ph%first_generic_to_blacs) then
      ham_den = 0.0_r8
      map_den = -1

      if(.not. ph%unit_ovlp) then
         ovlp_den = 0.0_r8
      end if

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
         map_den(l_row,l_col) = map_recv(i_val)

         if(.not. ph%unit_ovlp) then
            ovlp_den(l_row,l_col) = s_val_recv(i_val)
         end if
      end do

      call elsi_deallocate(bh,map_recv,"map_recv")

      if(.not. ph%unit_ovlp) then
         call elsi_deallocate(bh,s_val_recv,"s_val_recv")
      end if
   else
      ham_den = 0.0_r8

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
      end do
   end if

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_generic_to_blacs = .false.

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in generic COO format to 2D
!! block-cyclic dense format.
!!
subroutine elsi_generic_to_blacs_hs_cmplx(ph,bh,ham_sp,ovlp_sp,row_ind,col_ind,&
   ham_den,ovlp_den,map_den)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   complex(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   complex(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ind(bh%nnz_l_sp)
   complex(kind=r8), intent(out) :: ham_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(inout) :: ovlp_den(bh%n_lrow,bh%n_lcol)
   integer(kind=i4), intent(inout) :: map_den(bh%n_lrow,bh%n_lcol)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: l_col
   integer(kind=i4) :: l_row
   integer(kind=i4) :: p_col
   integer(kind=i4) :: p_row
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: h_val_send(:)
   complex(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: map_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   complex(kind=r8), allocatable :: h_val_recv(:)
   complex(kind=r8), allocatable :: s_val_recv(:)
   integer(kind=i4), allocatable :: map_recv(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_generic_to_blacs_hs_cmplx"

   call elsi_get_time(t0)

   if(ph%first_generic_to_blacs) then
      if(.not. ph%unit_ovlp) then
         call elsi_allocate(bh,s_val_send,bh%nnz_l_sp,"s_val_send",caller)
      end if

      call elsi_allocate(bh,map_send,bh%nnz_l_sp,"map_send",caller)

      map_send = bh%myid
   end if

   call elsi_allocate(bh,dest,bh%nnz_l_sp,"dest",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   row_send = row_ind
   col_send = col_ind
   h_val_send = ham_sp

   if(ph%first_generic_to_blacs .and. .not. ph%unit_ovlp) then
      s_val_send = ovlp_sp
   end if

   do i_val = 1,bh%nnz_l_sp
      ! Compute destination
      p_row = mod((row_send(i_val)-1)/bh%blk,bh%n_prow)
      p_col = mod((col_send(i_val)-1)/bh%blk,bh%n_pcol)
      dest(i_val) = p_col+p_row*bh%n_pcol

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   end do

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,dest,perm)
   call elsi_permute(bh%nnz_l_sp,perm,h_val_send)
   call elsi_permute(bh%nnz_l_sp,perm,row_send)
   call elsi_permute(bh%nnz_l_sp,perm,col_send)

   if(ph%first_generic_to_blacs .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp,perm,s_val_send)
   end if

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   bh%nnz_l = sum(recv_count,1)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Hamiltonian value
   call elsi_allocate(bh,h_val_recv,bh%nnz_l,"h_val_recv",caller)

   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_complex16,&
        h_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   if(ph%first_generic_to_blacs) then
      ! Overlap value
      if(.not. ph%unit_ovlp) then
         call elsi_allocate(bh,s_val_recv,bh%nnz_l,"s_val_recv",caller)

         call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
              s_val_recv,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

         call elsi_deallocate(bh,s_val_send,"s_val_send")
      end if

      call elsi_allocate(bh,map_recv,bh%nnz_l,"map_recv",caller)

      call MPI_Alltoallv(map_send,send_count,send_displ,mpi_integer4,map_recv,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,map_send,"map_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")

   ! Unpack matrix
   if(ph%first_generic_to_blacs) then
      ham_den = (0.0_r8,0.0_r8)
      map_den = -1

      if(.not. ph%unit_ovlp) then
         ovlp_den = (0.0_r8,0.0_r8)
      end if

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
         map_den(l_row,l_col) = map_recv(i_val)

         if(.not. ph%unit_ovlp) then
            ovlp_den(l_row,l_col) = s_val_recv(i_val)
         end if
      end do

      call elsi_deallocate(bh,map_recv,"map_recv")

      if(.not. ph%unit_ovlp) then
         call elsi_deallocate(bh,s_val_recv,"s_val_recv")
      end if
   else
      ham_den = (0.0_r8,0.0_r8)

      do i_val = 1,bh%nnz_l
         ! Compute local 2d id
         call elsi_get_lid(bh%n_prow,bh%blk,row_recv(i_val),l_row)
         call elsi_get_lid(bh%n_pcol,bh%blk,col_recv(i_val),l_col)

         ! Put value to correct position
         ham_den(l_row,l_col) = h_val_recv(i_val)
      end do
   end if

   call elsi_deallocate(bh,h_val_recv,"h_val_recv")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_generic_to_blacs = .false.

end subroutine

!>
!! Construct Halmitonian and overlep matrices in NTPoly format from matrices
!! stored in generic COO format.
!!
subroutine elsi_generic_to_ntpoly_hs_real(ph,bh,ham_sp,ovlp_sp,row_ind,col_ind,&
   ham_nt,ovlp_nt,map_nt)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   real(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ind(bh%nnz_l_sp)
   type(Matrix_ps), intent(inout) :: ham_nt
   type(Matrix_ps), intent(inout) :: ovlp_nt
   type(Matrix_ps), intent(inout) :: map_nt

   integer(kind=i4) :: i_val
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   type(TripletList_r) :: ham_list
   type(TripletList_r) :: ovlp_list
   type(TripletList_r) :: map_list
   type(Triplet_r) :: coo

   character(len=*), parameter :: caller = "elsi_generic_to_ntpoly_hs_real"

   call elsi_get_time(t0)

   if(ph%first_generic_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call ConstructEmptyMatrix(ovlp_nt,ph%n_basis,ph%nt_pgrid,.false.)
         call ConstructTripletList(ovlp_list)
      end if

      call ConstructEmptyMatrix(ham_nt,ph%n_basis,ph%nt_pgrid,.false.)
      call ConstructEmptyMatrix(map_nt,ph%n_basis,ph%nt_pgrid,.false.)
      call ConstructTripletList(map_list)
   end if

   call ConstructTripletList(ham_list)

   do i_val = 1,bh%nnz_l_sp
      coo%point_value = ham_sp(i_val)
      coo%index_row = row_ind(i_val)
      coo%index_column = col_ind(i_val)

      call AppendToTripletList(ham_list,coo)

      if(ph%first_generic_to_ntpoly) then
         if(.not. ph%unit_ovlp) then
            coo%point_value = ovlp_sp(i_val)

            call AppendToTripletList(ovlp_list,coo)
         end if

         coo%point_value = real(bh%myid,kind=r8)

         call AppendToTripletList(map_list,coo)
      end if
   end do

   if(ph%first_generic_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call FillMatrixFromTripletList(ovlp_nt,ovlp_list)
         call DestructTripletList(ovlp_list)
      end if

      call FillMatrixFromTripletList(map_nt,map_list)
      call DestructTripletList(map_list)
   end if

   call FillMatrixFromTripletList(ham_nt,ham_list)
   call DestructTripletList(ham_list)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_generic_to_ntpoly = .false.

end subroutine

!>
!! Construct Halmitonian and overlep matrices in NTPoly format from matrices
!! stored in generic COO format.
!!
subroutine elsi_generic_to_ntpoly_hs_cmplx(ph,bh,ham_sp,ovlp_sp,row_ind,&
   col_ind,ham_nt,ovlp_nt,map_nt)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(in) :: ham_sp(bh%nnz_l_sp)
   complex(kind=r8), intent(in) :: ovlp_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: col_ind(bh%nnz_l_sp)
   type(Matrix_ps), intent(inout) :: ham_nt
   type(Matrix_ps), intent(inout) :: ovlp_nt
   type(Matrix_ps), intent(inout) :: map_nt

   integer(kind=i4) :: i_val
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   type(TripletList_r) :: ham_list
   type(TripletList_r) :: ovlp_list
   type(TripletList_r) :: map_list
   type(Triplet_r) :: coo

   character(len=*), parameter :: caller = "elsi_generic_to_ntpoly_hs_cmplx"

   call elsi_get_time(t0)

   if(ph%first_generic_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call ConstructEmptyMatrix(ovlp_nt,ph%n_basis,ph%nt_pgrid,.true.)
         call ConstructTripletList(ovlp_list)
      end if

      call ConstructEmptyMatrix(ham_nt,ph%n_basis,ph%nt_pgrid,.true.)
      call ConstructEmptyMatrix(map_nt,ph%n_basis,ph%nt_pgrid,.true.)
      call ConstructTripletList(map_list)
   end if

   call ConstructTripletList(ham_list)

   do i_val = 1,bh%nnz_l_sp
      coo%point_value = ham_sp(i_val)
      coo%index_row = row_ind(i_val)
      coo%index_column = col_ind(i_val)

      call AppendToTripletList(ham_list,coo)

      if(ph%first_generic_to_ntpoly) then
         if(.not. ph%unit_ovlp) then
            coo%point_value = ovlp_sp(i_val)

            call AppendToTripletList(ovlp_list,coo)
         end if

         coo%point_value = cmplx(bh%myid,kind=r8)

         call AppendToTripletList(map_list,coo)
      end if
   end do

   if(ph%first_generic_to_ntpoly) then
      if(.not. ph%unit_ovlp) then
         call FillMatrixFromTripletList(ovlp_nt,ovlp_list)
         call DestructTripletList(ovlp_list)
      end if

      call FillMatrixFromTripletList(map_nt,map_list)
      call DestructTripletList(map_list)
   end if

   call FillMatrixFromTripletList(ham_nt,ham_list)
   call DestructTripletList(ham_list)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_generic_to_ntpoly = .false.

end subroutine

!>
!! Get local number of nonzero elements in matrices in 1D block CSC format
!! converted from generic COO format.
!!
subroutine elsi_generic_to_pexsi_hs_dim(ph,bh,col_ind3)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4), intent(in) :: col_ind3(bh%nnz_l_sp3)

   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=*), parameter :: caller = "elsi_generic_to_pexsi_hs_dim"

   call elsi_allocate(bh,dest,ph%pexsi_np_per_pole,"dest",caller)
   call elsi_allocate(bh,nnz,ph%pexsi_np_per_pole,"nnz",caller)

   do i_val = 1,bh%nnz_l_sp3
      i_proc = (col_ind3(i_val)-1)/(ph%n_basis/ph%pexsi_np_per_pole)
      i_proc = min(i_proc,ph%pexsi_np_per_pole-1)

      dest(i_proc+1) = dest(i_proc+1)+1
   end do

   call MPI_Allreduce(dest,nnz,ph%pexsi_np_per_pole,mpi_integer4,mpi_sum,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   bh%nnz_l_sp1 = nnz(ph%pexsi_my_pcol+1)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in generic COO format to 1D
!! block CSC format.
!!
subroutine elsi_generic_to_pexsi_hs_real(ph,bh,ham_sp3,ovlp_sp3,row_ind3,&
   col_ind3,ham_sp1,ovlp_sp1,row_ind1,col_ptr1,map_sp1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: ham_sp3(bh%nnz_l_sp3)
   real(kind=r8), intent(in) :: ovlp_sp3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: row_ind3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: col_ind3(bh%nnz_l_sp3)
   real(kind=r8), intent(out) :: ham_sp1(bh%nnz_l_sp1)
   real(kind=r8), intent(inout) :: ovlp_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: col_ptr1(bh%n_lcol_sp1+1)
   integer(kind=i4), intent(inout) :: map_sp1(bh%nnz_l_sp1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: n_lcol_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: h_val_send(:)
   real(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: map_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_generic_to_pexsi_hs_real"

   call elsi_get_time(t0)

   n_lcol_aux = ph%n_basis/ph%pexsi_np_per_pole

   if(ph%first_generic_to_pexsi) then
      if(.not. ph%unit_ovlp) then
         call elsi_allocate(bh,s_val_send,bh%nnz_l_sp3,"s_val_send",caller)
      end if

      call elsi_allocate(bh,map_send,bh%nnz_l_sp3,"map_send",caller)

      map_send = bh%myid
   end if

   call elsi_allocate(bh,dest,bh%nnz_l_sp3,"dest",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp3,"perm",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp3,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp3,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp3,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   h_val_send = ham_sp3
   row_send = row_ind3
   col_send = col_ind3

   if(ph%first_generic_to_pexsi .and. .not. ph%unit_ovlp) then
      s_val_send = ovlp_sp3
   end if

   do i_val = 1,bh%nnz_l_sp3
      ! Compute destination
      dest(i_val) = (col_send(i_val)-1)/n_lcol_aux
      dest(i_val) = min(dest(i_val),ph%pexsi_np_per_pole-1)

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   end do

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp3,dest,perm)
   call elsi_permute(bh%nnz_l_sp3,perm,h_val_send)
   call elsi_permute(bh%nnz_l_sp3,perm,row_send)
   call elsi_permute(bh%nnz_l_sp3,perm,col_send)

   if(ph%first_generic_to_pexsi .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp3,perm,s_val_send)
   end if

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
   ! Row index
   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_ind1,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp1,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_real8,ham_sp1,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   if(ph%first_generic_to_pexsi) then
      if(.not. ph%unit_ovlp) then
         ! Overlap value
         call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_real8,&
              ovlp_sp1,recv_count,recv_displ,mpi_real8,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

         call elsi_deallocate(bh,s_val_send,"s_val_send")
      end if

      call MPI_Alltoallv(map_send,send_count,send_displ,mpi_integer4,map_sp1,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,map_send,"map_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp1,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp1,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_ind1,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp1,gid,perm)
   call elsi_permute(bh%nnz_l_sp1,perm,ham_sp1)
   call elsi_permute(bh%nnz_l_sp1,perm,row_ind1)
   call elsi_permute(bh%nnz_l_sp1,perm,col_recv)

   if(ph%first_generic_to_pexsi) then
      if(.not. ph%unit_ovlp) then
         call elsi_permute(bh%nnz_l_sp1,perm,ovlp_sp1)
      end if

      call elsi_permute(bh%nnz_l_sp1,perm,map_sp1)
   end if

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")

   if(ph%first_generic_to_pexsi) then
      ! Only 1st pole computes row index and column pointer
      if(ph%pexsi_my_prow == 0) then
         col_ptr1 = 0
         col_ptr1(bh%n_lcol_sp1+1) = bh%nnz_l_sp1+1

         do i_val = 1,bh%nnz_l_sp1
            i_col = col_recv(i_val)-n_lcol_aux*ph%pexsi_my_pcol
            col_ptr1(i_col) = col_ptr1(i_col)+1
         end do

         do i_col = bh%n_lcol_sp1,1,-1
            col_ptr1(i_col) = col_ptr1(i_col+1)-col_ptr1(i_col)
         end do
      end if
   end if

   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_generic_to_pexsi = .false.

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in generic COO format to 1D
!! block CSC format.
!!
subroutine elsi_generic_to_pexsi_hs_cmplx(ph,bh,ham_sp3,ovlp_sp3,row_ind3,&
   col_ind3,ham_sp1,ovlp_sp1,row_ind1,col_ptr1,map_sp1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(in) :: ham_sp3(bh%nnz_l_sp3)
   complex(kind=r8), intent(in) :: ovlp_sp3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: row_ind3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: col_ind3(bh%nnz_l_sp3)
   complex(kind=r8), intent(out) :: ham_sp1(bh%nnz_l_sp1)
   complex(kind=r8), intent(inout) :: ovlp_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: col_ptr1(bh%n_lcol_sp1+1)
   integer(kind=i4), intent(inout) :: map_sp1(bh%nnz_l_sp1)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: n_lcol_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: h_val_send(:)
   complex(kind=r8), allocatable :: s_val_send(:)
   integer(kind=i4), allocatable :: map_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_generic_to_pexsi_hs_cmplx"

   call elsi_get_time(t0)

   n_lcol_aux = ph%n_basis/ph%pexsi_np_per_pole

   if(ph%first_generic_to_pexsi) then
      if(.not. ph%unit_ovlp) then
         call elsi_allocate(bh,s_val_send,bh%nnz_l_sp3,"s_val_send",caller)
      end if

      call elsi_allocate(bh,map_send,bh%nnz_l_sp3,"map_send",caller)

      map_send = bh%myid
   end if

   call elsi_allocate(bh,dest,bh%nnz_l_sp3,"dest",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp3,"perm",caller)
   call elsi_allocate(bh,h_val_send,bh%nnz_l_sp3,"h_val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp3,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp3,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   h_val_send = ham_sp3
   row_send = row_ind3
   col_send = col_ind3

   if(ph%first_generic_to_pexsi .and. .not. ph%unit_ovlp) then
      s_val_send = ovlp_sp3
   end if

   do i_val = 1,bh%nnz_l_sp3
      ! Compute destination
      dest(i_val) = (col_send(i_val)-1)/n_lcol_aux
      dest(i_val) = min(dest(i_val),ph%pexsi_np_per_pole-1)

      ! Set send_count
      send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
   end do

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp3,dest,perm)
   call elsi_permute(bh%nnz_l_sp3,perm,h_val_send)
   call elsi_permute(bh%nnz_l_sp3,perm,row_send)
   call elsi_permute(bh%nnz_l_sp3,perm,col_send)

   if(ph%first_generic_to_pexsi .and. .not. ph%unit_ovlp) then
      call elsi_permute(bh%nnz_l_sp3,perm,s_val_send)
   end if

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
   ! Row index
   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_ind1,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp1,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_send,send_count,send_displ,mpi_complex16,ham_sp1,&
        recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,h_val_send,"h_val_send")

   if(ph%first_generic_to_pexsi) then
      if(.not. ph%unit_ovlp) then
         ! Overlap value
         call MPI_Alltoallv(s_val_send,send_count,send_displ,mpi_complex16,&
              ovlp_sp1,recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

         call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

         call elsi_deallocate(bh,s_val_send,"s_val_send")
      end if

      call MPI_Alltoallv(map_send,send_count,send_displ,mpi_integer4,map_sp1,&
           recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

      call elsi_deallocate(bh,map_send,"map_send")
   end if

   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp1,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp1,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_ind1,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp1,gid,perm)
   call elsi_permute(bh%nnz_l_sp1,perm,ham_sp1)
   call elsi_permute(bh%nnz_l_sp1,perm,row_ind1)
   call elsi_permute(bh%nnz_l_sp1,perm,col_recv)

   if(ph%first_generic_to_pexsi) then
      if(.not. ph%unit_ovlp) then
         call elsi_permute(bh%nnz_l_sp1,perm,ovlp_sp1)
      end if

      call elsi_permute(bh%nnz_l_sp1,perm,map_sp1)
   end if

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")

   if(ph%first_generic_to_pexsi) then
      ! Only 1st pole computes row index and column pointer
      if(ph%pexsi_my_prow == 0) then
         col_ptr1 = 0
         col_ptr1(bh%n_lcol_sp1+1) = bh%nnz_l_sp1+1

         do i_val = 1,bh%nnz_l_sp1
            i_col = col_recv(i_val)-n_lcol_aux*ph%pexsi_my_pcol
            col_ptr1(i_col) = col_ptr1(i_col)+1
         end do

         do i_col = bh%n_lcol_sp1,1,-1
            col_ptr1(i_col) = col_ptr1(i_col+1)-col_ptr1(i_col)
         end do
      end if
   end if

   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%first_generic_to_pexsi = .false.

end subroutine

!>
!! Get local number of nonzero elements in matrices in 1D block CSC format
!! converted from generic COO format.
!!
subroutine elsi_generic_to_sips_hs_dim(ph,bh,col_ind3)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   integer(kind=i4), intent(in) :: col_ind3(bh%nnz_l_sp3)

   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: ierr

   integer(kind=i4), allocatable :: dest(:)
   integer(kind=i4), allocatable :: nnz(:)

   character(len=*), parameter :: caller = "elsi_generic_to_sips_hs_dim"

   call elsi_allocate(bh,dest,bh%n_procs,"dest",caller)
   call elsi_allocate(bh,nnz,bh%n_procs,"nnz",caller)

   do i_val = 1,bh%nnz_l_sp3
      i_proc = (col_ind3(i_val)-1)/(ph%n_basis/bh%n_procs)
      i_proc = min(i_proc,bh%n_procs-1)

      dest(i_proc+1) = dest(i_proc+1)+1
   end do

   call MPI_Allreduce(dest,nnz,bh%n_procs,mpi_integer4,mpi_sum,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   bh%nnz_l_sp1 = nnz(bh%myid+1)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,nnz,"nnz")

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in generic COO format to 1D
!! block CSC format.
!!
subroutine elsi_generic_to_sips_hs_real(ph,bh,ham_sp3,ovlp_sp3,row_ind3,&
   col_ind3,ham_sp1,ovlp_sp1,row_ind1,col_ptr1,map_sp1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: ham_sp3(bh%nnz_l_sp3)
   real(kind=r8), intent(in) :: ovlp_sp3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: row_ind3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: col_ind3(bh%nnz_l_sp3)
   real(kind=r8), intent(out) :: ham_sp1(bh%nnz_l_sp1)
   real(kind=r8), intent(inout) :: ovlp_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: col_ptr1(bh%n_lcol_sp1+1)
   integer(kind=i4), intent(inout) :: map_sp1(bh%nnz_l_sp1)

   character(len=*), parameter :: caller = "elsi_generic_to_sips_hs_real"

   ph%pexsi_my_prow = 0
   ph%pexsi_my_pcol = bh%myid
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_generic_to_pexsi_hs_real(ph,bh,ham_sp3,ovlp_sp3,row_ind3,col_ind3,&
        ham_sp1,ovlp_sp1,row_ind1,col_ptr1,map_sp1)

end subroutine

!>
!! Convert Halmitonian and overlap matrices stored in generic COO format to 1D
!! block CSC format.
!!
subroutine elsi_generic_to_sips_hs_cmplx(ph,bh,ham_sp3,ovlp_sp3,row_ind3,&
   col_ind3,ham_sp1,ovlp_sp1,row_ind1,col_ptr1,map_sp1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(in) :: ham_sp3(bh%nnz_l_sp3)
   complex(kind=r8), intent(in) :: ovlp_sp3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: row_ind3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: col_ind3(bh%nnz_l_sp3)
   complex(kind=r8), intent(out) :: ham_sp1(bh%nnz_l_sp1)
   complex(kind=r8), intent(inout) :: ovlp_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: col_ptr1(bh%n_lcol_sp1+1)
   integer(kind=i4), intent(inout) :: map_sp1(bh%nnz_l_sp1)

   character(len=*), parameter :: caller = "elsi_generic_to_sips_hs_cmplx"

   ph%pexsi_my_prow = 0
   ph%pexsi_my_pcol = bh%myid
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_generic_to_pexsi_hs_cmplx(ph,bh,ham_sp3,ovlp_sp3,row_ind3,col_ind3,&
        ham_sp1,ovlp_sp1,row_ind1,col_ptr1,map_sp1)

end subroutine

!>
!! Convert density matrix in 2D block-cyclic dense format to generic COO format.
!!
subroutine elsi_blacs_to_generic_dm_real(ph,bh,dm_den,map_den,dm_sp,perm_sp)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: dm_den(bh%n_lrow,bh%n_lcol)
   integer(kind=i4), intent(in) :: map_den(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: perm_sp(bh%nnz_l_sp)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_row
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_blacs_to_generic_dm_real"

   call elsi_get_time(t0)

   nnz_l_aux = 0

   do i_col = 1,bh%n_lcol
      do i_row = 1,bh%n_lrow
         if(map_den(i_row,i_col) > -1) then
            nnz_l_aux = nnz_l_aux+1
         end if
      end do
   end do

   call elsi_allocate(bh,dest,nnz_l_aux,"dest",caller)
   call elsi_allocate(bh,perm,nnz_l_aux,"perm",caller)
   call elsi_allocate(bh,val_send,nnz_l_aux,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_aux,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_aux,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_val = 0

   do i_col = 1,bh%n_lcol
      do i_row = 1,bh%n_lrow
         if(map_den(i_row,i_col) > -1) then
            i_val = i_val+1

            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,row_send(i_val))
            call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,col_send(i_val))

            ! Compute destination
            dest(i_val) = map_den(i_row,i_col)

            ! Pack data
            val_send(i_val) = dm_den(i_row,i_col)

            ! Set send_count
            send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
         end if
      end do
   end do

   ! Sort
   call elsi_heapsort(nnz_l_aux,dest,perm)
   call elsi_permute(nnz_l_aux,perm,val_send)
   call elsi_permute(nnz_l_aux,perm,row_send)
   call elsi_permute(nnz_l_aux,perm,col_send)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Density matrix value
   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_real8,dm_sp,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_recv,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,gid,perm)
   call elsi_permute(bh%nnz_l_sp,perm,dm_sp)
   call elsi_unpermute(bh%nnz_l_sp,perm_sp,dm_sp)

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix in 2D block-cyclic dense format to generic COO format.
!!
subroutine elsi_blacs_to_generic_dm_cmplx(ph,bh,dm_den,map_den,dm_sp,perm_sp)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(in) :: dm_den(bh%n_lrow,bh%n_lcol)
   integer(kind=i4), intent(in) :: map_den(bh%n_lrow,bh%n_lcol)
   complex(kind=r8), intent(out) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: perm_sp(bh%nnz_l_sp)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_row
   integer(kind=i4) :: nnz_l_aux
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_blacs_to_generic_dm_cmplx"

   call elsi_get_time(t0)

   nnz_l_aux = 0

   do i_col = 1,bh%n_lcol
      do i_row = 1,bh%n_lrow
         if(map_den(i_row,i_col) > -1) then
            nnz_l_aux = nnz_l_aux+1
         end if
      end do
   end do

   call elsi_allocate(bh,dest,nnz_l_aux,"dest",caller)
   call elsi_allocate(bh,perm,nnz_l_aux,"perm",caller)
   call elsi_allocate(bh,val_send,nnz_l_aux,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_aux,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_aux,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   i_val = 0

   do i_col = 1,bh%n_lcol
      do i_row = 1,bh%n_lrow
         if(map_den(i_row,i_col) > -1) then
            i_val = i_val+1

            call elsi_get_gid(bh%my_prow,bh%n_prow,bh%blk,i_row,row_send(i_val))
            call elsi_get_gid(bh%my_pcol,bh%n_pcol,bh%blk,i_col,col_send(i_val))

            ! Compute destination
            dest(i_val) = map_den(i_row,i_col)

            ! Pack data
            val_send(i_val) = dm_den(i_row,i_col)

            ! Set send_count
            send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
         end if
      end do
   end do

   ! Sort
   call elsi_heapsort(nnz_l_aux,dest,perm)
   call elsi_permute(nnz_l_aux,perm,val_send)
   call elsi_permute(nnz_l_aux,perm,row_send)
   call elsi_permute(nnz_l_aux,perm,col_send)

   call elsi_deallocate(bh,dest,"dest")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Density matrix value
   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_complex16,dm_sp,&
        recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_recv,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,gid,perm)
   call elsi_permute(bh%nnz_l_sp,perm,dm_sp)
   call elsi_unpermute(bh%nnz_l_sp,perm_sp,dm_sp)

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix stored in 1D block CSC format to generic COO format.
!!
subroutine elsi_pexsi_to_generic_dm_real(ph,bh,dm_sp1,row_ind1,col_ptr1,&
   map_sp1,dm_sp3,perm_sp3)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: dm_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: col_ptr1(bh%n_lcol_sp1+1)
   integer(kind=i4), intent(in) :: map_sp1(bh%nnz_l_sp1)
   real(kind=r8), intent(out) :: dm_sp3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: perm_sp3(bh%nnz_l_sp3)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_proc
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_pexsi_to_generic_dm_real"

   call elsi_get_time(t0)

   call elsi_allocate(bh,val_send,bh%nnz_l_sp1,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp1,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp1,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(ph%pexsi_my_prow == 0) then
      call elsi_allocate(bh,dest,bh%nnz_l_sp1,"dest",caller)
      call elsi_allocate(bh,perm,bh%nnz_l_sp1,"perm",caller)

      i_col = 0

      do i_val = 1,bh%nnz_l_sp1
         do while(i_val == col_ptr1(i_col+1) .and. i_col /= bh%n_lcol_sp1)
            i_col = i_col+1
         end do

         i_row = row_ind1(i_val)

         ! Compute global id
         row_send(i_val) = i_row
         col_send(i_val) = i_col+bh%myid*(ph%n_basis/ph%pexsi_np_per_pole)
         val_send(i_val) = dm_sp1(i_val)

         ! Compute destination
         dest(i_val) = map_sp1(i_val)

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      end do

      ! Sort
      call elsi_heapsort(bh%nnz_l_sp1,dest,perm)
      call elsi_permute(bh%nnz_l_sp1,perm,val_send)
      call elsi_permute(bh%nnz_l_sp1,perm,row_send)
      call elsi_permute(bh%nnz_l_sp1,perm,col_send)

      call elsi_deallocate(bh,dest,"dest")
      call elsi_deallocate(bh,perm,"perm")
   end if

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l_sp3,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp3,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Density matrix value
   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_real8,dm_sp3,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp3,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp3,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_recv,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp3,gid,perm)
   call elsi_permute(bh%nnz_l_sp3,perm,dm_sp3)
   call elsi_unpermute(bh%nnz_l_sp3,perm_sp3,dm_sp3)

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix stored in 1D block CSC format to generic COO format.
!!
subroutine elsi_pexsi_to_generic_dm_cmplx(ph,bh,dm_sp1,row_ind1,col_ptr1,&
   map_sp1,dm_sp3,perm_sp3)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(in) :: dm_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: col_ptr1(bh%n_lcol_sp1+1)
   integer(kind=i4), intent(in) :: map_sp1(bh%nnz_l_sp1)
   complex(kind=r8), intent(out) :: dm_sp3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: perm_sp3(bh%nnz_l_sp3)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_proc
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: perm(:)

   character(len=*), parameter :: caller = "elsi_pexsi_to_generic_dm_cmplx"

   call elsi_get_time(t0)

   call elsi_allocate(bh,val_send,bh%nnz_l_sp1,"val_send",caller)
   call elsi_allocate(bh,row_send,bh%nnz_l_sp1,"row_send",caller)
   call elsi_allocate(bh,col_send,bh%nnz_l_sp1,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(ph%pexsi_my_prow == 0) then
      call elsi_allocate(bh,dest,bh%nnz_l_sp1,"dest",caller)
      call elsi_allocate(bh,perm,bh%nnz_l_sp1,"perm",caller)

      i_col = 0

      do i_val = 1,bh%nnz_l_sp1
         do while(i_val == col_ptr1(i_col+1) .and. i_col /= bh%n_lcol_sp1)
            i_col = i_col+1
         end do

         i_row = row_ind1(i_val)

         ! Compute global id
         row_send(i_val) = i_row
         col_send(i_val) = i_col+bh%myid*(ph%n_basis/ph%pexsi_np_per_pole)
         val_send(i_val) = dm_sp1(i_val)

         ! Compute destination
         dest(i_val) = map_sp1(i_val)

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      end do

      ! Sort
      call elsi_heapsort(bh%nnz_l_sp1,dest,perm)
      call elsi_permute(bh%nnz_l_sp1,perm,val_send)
      call elsi_permute(bh%nnz_l_sp1,perm,row_send)
      call elsi_permute(bh%nnz_l_sp1,perm,col_send)

      call elsi_deallocate(bh,dest,"dest")
      call elsi_deallocate(bh,perm,"perm")
   end if

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
   ! Row index
   call elsi_allocate(bh,row_recv,bh%nnz_l_sp3,"row_recv",caller)

   call MPI_Alltoallv(row_send,send_count,send_displ,mpi_integer4,row_recv,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,row_send,"row_send")

   ! Column index
   call elsi_allocate(bh,col_recv,bh%nnz_l_sp3,"col_recv",caller)

   call MPI_Alltoallv(col_send,send_count,send_displ,mpi_integer4,col_recv,&
        recv_count,recv_displ,mpi_integer4,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,col_send,"col_send")

   ! Density matrix value
   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_complex16,dm_sp3,&
        recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp3,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp3,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_recv,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp3,gid,perm)
   call elsi_permute(bh%nnz_l_sp3,perm,dm_sp3)
   call elsi_unpermute(bh%nnz_l_sp3,perm_sp3,dm_sp3)

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix stored in 1D block CSC format to generic COO format.
!!
subroutine elsi_sips_to_generic_dm_real(ph,bh,dm_sp1,row_ind1,col_ptr1,map_sp1,&
   dm_sp3,perm_sp3)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: dm_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: col_ptr1(bh%n_lcol_sp1+1)
   integer(kind=i4), intent(in) :: map_sp1(bh%nnz_l_sp1)
   real(kind=r8), intent(out) :: dm_sp3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: perm_sp3(bh%nnz_l_sp3)

   character(len=*), parameter :: caller = "elsi_sips_to_generic_dm_real"

   ph%pexsi_my_prow = 0
   ph%pexsi_my_pcol = bh%myid
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_pexsi_to_generic_dm_real(ph,bh,dm_sp1,row_ind1,col_ptr1,map_sp1,&
        dm_sp3,perm_sp3)

end subroutine

!>
!! Convert density matrix stored in 1D block CSC format to generic COO format.
!!
subroutine elsi_sips_to_generic_dm_cmplx(ph,bh,dm_sp1,row_ind1,col_ptr1,&
   map_sp1,dm_sp3,perm_sp3)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(in) :: dm_sp1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: row_ind1(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: col_ptr1(bh%n_lcol_sp1+1)
   integer(kind=i4), intent(in) :: map_sp1(bh%nnz_l_sp1)
   complex(kind=r8), intent(out) :: dm_sp3(bh%nnz_l_sp3)
   integer(kind=i4), intent(in) :: perm_sp3(bh%nnz_l_sp3)

   character(len=*), parameter :: caller = "elsi_sips_to_generic_dm_cmplx"

   ph%pexsi_my_prow = 0
   ph%pexsi_my_pcol = bh%myid
   ph%pexsi_np_per_pole = bh%n_procs

   call elsi_pexsi_to_generic_dm_cmplx(ph,bh,dm_sp1,row_ind1,col_ptr1,map_sp1,&
        dm_sp3,perm_sp3)

end subroutine

!>
!! Convert density matrix computed by NTPoly to generic COO format.
!!
subroutine elsi_ntpoly_to_generic_dm_real(ph,bh,dm_nt,map_nt,dm_sp,perm_sp)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(Matrix_ps), intent(inout) :: dm_nt
   type(Matrix_ps), intent(inout) :: map_nt
   real(kind=r8), intent(out) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: perm_sp(bh%nnz_l_sp)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: gid_map
   integer(kind=i4) :: gid_dm
   integer(kind=i4) :: nnz_l_nt
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   real(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: perm(:)

   type(TripletList_r) :: dm_list
   type(TripletList_r) :: map_list

   character(len=*), parameter :: caller = "elsi_ntpoly_to_generic_dm_real"

   call elsi_get_time(t0)

   call GetMatrixTripletList(dm_nt,dm_list)
   call GetMatrixTripletList(map_nt,map_list)

   nnz_l_nt = map_list%CurrentSize

   call elsi_allocate(bh,val_send,nnz_l_nt,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_nt,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_nt,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(bh%myid < bh%n_procs/ph%nt_n_layers) then
      call elsi_allocate(bh,dest,nnz_l_nt,"dest",caller)
      call elsi_allocate(bh,perm,nnz_l_nt,"perm",caller)

      ! Compute destination
      do i_val = 1,nnz_l_nt
         dest(i_val) = nint(map_list%data(i_val)%point_value,kind=i4)
         map_list%data(i_val)%point_value = 0.0_r8
      end do

      i_val = 1
      j_val = 1

      do while(i_val <= nnz_l_nt .and. j_val <= dm_list%CurrentSize)
         ! Compute global 1D id
         i_row = map_list%data(i_val)%index_row
         i_col = map_list%data(i_val)%index_column
         gid_map = int(i_col-1,kind=i8)*int(ph%n_basis,kind=i8)&
            +int(i_row,kind=i8)
         i_row = dm_list%data(j_val)%index_row
         i_col = dm_list%data(j_val)%index_column
         gid_dm = int(i_col-1,kind=i8)*int(ph%n_basis,kind=i8)&
            +int(i_row,kind=i8)

         if(gid_map == gid_dm) then
            map_list%data(i_val)%point_value = dm_list%data(j_val)%point_value
            i_val = i_val+1
            j_val = j_val+1
         else if(gid_map > gid_dm) then
            j_val = j_val+1
         else
            i_val = i_val+1
         end if
      end do

      do i_val = 1,nnz_l_nt
         ! Compute global id
         row_send(i_val) = map_list%data(i_val)%index_row
         col_send(i_val) = map_list%data(i_val)%index_column
         val_send(i_val) = map_list%data(i_val)%point_value

         ! Set send_count
        send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      end do

      ! Sort
      call elsi_heapsort(nnz_l_nt,dest,perm)
      call elsi_permute(nnz_l_nt,perm,val_send)
      call elsi_permute(nnz_l_nt,perm,row_send)
      call elsi_permute(nnz_l_nt,perm,col_send)

      call elsi_deallocate(bh,dest,"dest")
      call elsi_deallocate(bh,perm,"perm")
   end if

   call DestructTripletList(dm_list)
   call DestructTripletList(map_list)

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Density matrix value
   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_real8,dm_sp,&
        recv_count,recv_displ,mpi_real8,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_recv,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,gid,perm)
   call elsi_permute(bh%nnz_l_sp,perm,dm_sp)
   call elsi_unpermute(bh%nnz_l_sp,perm_sp,dm_sp)

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Convert density matrix computed by NTPoly to generic COO format.
!!
subroutine elsi_ntpoly_to_generic_dm_cmplx(ph,bh,dm_nt,map_nt,dm_sp,perm_sp)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(Matrix_ps), intent(inout) :: dm_nt
   type(Matrix_ps), intent(inout) :: map_nt
   complex(kind=r8), intent(out) :: dm_sp(bh%nnz_l_sp)
   integer(kind=i4), intent(in) :: perm_sp(bh%nnz_l_sp)

   integer(kind=i4) :: ierr
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: gid_map
   integer(kind=i4) :: gid_dm
   integer(kind=i4) :: nnz_l_nt
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   ! See documentation of MPI_Alltoallv
   complex(kind=r8), allocatable :: val_send(:)
   integer(kind=i4), allocatable :: row_send(:)
   integer(kind=i4), allocatable :: col_send(:)
   integer(kind=i4), allocatable :: send_count(:)
   integer(kind=i4), allocatable :: send_displ(:)
   integer(kind=i4), allocatable :: row_recv(:)
   integer(kind=i4), allocatable :: col_recv(:)
   integer(kind=i4), allocatable :: recv_count(:)
   integer(kind=i4), allocatable :: recv_displ(:)
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i8), allocatable :: gid(:) ! Global 1D id
   integer(kind=i4), allocatable :: perm(:)

   type(TripletList_r) :: dm_list
   type(TripletList_r) :: map_list

   character(len=*), parameter :: caller = "elsi_ntpoly_to_generic_dm_cmplx"

   call elsi_get_time(t0)

   call GetMatrixTripletList(dm_nt,dm_list)
   call GetMatrixTripletList(map_nt,map_list)

   nnz_l_nt = map_list%CurrentSize

   call elsi_allocate(bh,val_send,nnz_l_nt,"val_send",caller)
   call elsi_allocate(bh,row_send,nnz_l_nt,"row_send",caller)
   call elsi_allocate(bh,col_send,nnz_l_nt,"col_send",caller)
   call elsi_allocate(bh,send_count,bh%n_procs,"send_count",caller)

   if(bh%myid < bh%n_procs/ph%nt_n_layers) then
      call elsi_allocate(bh,dest,nnz_l_nt,"dest",caller)
      call elsi_allocate(bh,perm,nnz_l_nt,"perm",caller)

      ! Compute destination
      do i_val = 1,nnz_l_nt
         dest(i_val) = nint(map_list%data(i_val)%point_value,kind=i4)
         map_list%data(i_val)%point_value = (0.0_r8,0.0_r8)
      end do

      i_val = 1
      j_val = 1

      do while(i_val <= nnz_l_nt .and. j_val <= dm_list%CurrentSize)
         ! Compute global 1D id
         i_row = map_list%data(i_val)%index_row
         i_col = map_list%data(i_val)%index_column
         gid_map = int(i_col-1,kind=i8)*int(ph%n_basis,kind=i8)&
            +int(i_row,kind=i8)
         i_row = dm_list%data(j_val)%index_row
         i_col = dm_list%data(j_val)%index_column
         gid_dm = int(i_col-1,kind=i8)*int(ph%n_basis,kind=i8)&
            +int(i_row,kind=i8)

         if(gid_map == gid_dm) then
            map_list%data(i_val)%point_value = dm_list%data(j_val)%point_value
            i_val = i_val+1
            j_val = j_val+1
         else if(gid_map > gid_dm) then
            j_val = j_val+1
         else
            i_val = i_val+1
         end if
      end do

      do i_val = 1,nnz_l_nt
         ! Compute global id
         row_send(i_val) = map_list%data(i_val)%index_row
         col_send(i_val) = map_list%data(i_val)%index_column
         val_send(i_val) = map_list%data(i_val)%point_value

         ! Set send_count
         send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
      end do

      ! Sort
      call elsi_heapsort(nnz_l_nt,dest,perm)
      call elsi_permute(nnz_l_nt,perm,val_send)
      call elsi_permute(nnz_l_nt,perm,row_send)
      call elsi_permute(nnz_l_nt,perm,col_send)

      call elsi_deallocate(bh,dest,"dest")
      call elsi_deallocate(bh,perm,"perm")
   end if

   call DestructTripletList(dm_list)
   call DestructTripletList(map_list)

   call elsi_allocate(bh,recv_count,bh%n_procs,"recv_count",caller)
   call elsi_allocate(bh,send_displ,bh%n_procs,"send_displ",caller)
   call elsi_allocate(bh,recv_displ,bh%n_procs,"recv_displ",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer4,recv_count,1,mpi_integer4,&
        bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoall",ierr,caller)

   ! Set send_displ and recv_displ
   do i_proc = 2,bh%n_procs
      send_displ(i_proc) = sum(send_count(1:i_proc-1),1)
      recv_displ(i_proc) = sum(recv_count(1:i_proc-1),1)
   end do

   ! Redistribute packed data
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

   ! Density matrix value
   call MPI_Alltoallv(val_send,send_count,send_displ,mpi_complex16,dm_sp,&
        recv_count,recv_displ,mpi_complex16,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Alltoallv",ierr,caller)

   call elsi_deallocate(bh,val_send,"val_send")
   call elsi_deallocate(bh,send_count,"send_count")
   call elsi_deallocate(bh,recv_count,"recv_count")
   call elsi_deallocate(bh,send_displ,"send_displ")
   call elsi_deallocate(bh,recv_displ,"recv_displ")
   call elsi_allocate(bh,gid,bh%nnz_l_sp,"gid",caller)
   call elsi_allocate(bh,perm,bh%nnz_l_sp,"perm",caller)

   ! Compute global 1D id
   gid = int(col_recv-1,kind=i8)*int(ph%n_basis,kind=i8)+int(row_recv,kind=i8)

   ! Sort
   call elsi_heapsort(bh%nnz_l_sp,gid,perm)
   call elsi_permute(bh%nnz_l_sp,perm,dm_sp)
   call elsi_unpermute(bh%nnz_l_sp,perm_sp,dm_sp)

   call elsi_deallocate(bh,gid,"gid")
   call elsi_deallocate(bh,perm,"perm")
   call elsi_deallocate(bh,row_recv,"row_recv")
   call elsi_deallocate(bh,col_recv,"col_recv")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished matrix redistribution"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

end module ELSI_REDIST
