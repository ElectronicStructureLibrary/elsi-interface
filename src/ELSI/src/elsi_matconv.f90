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

   use ISO_C_BINDING
   use ELSI_CONSTANTS, only: ELPA,LIBOMM
   use ELSI_DIMENSIONS, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4,i8
   use ELSI_TIMERS
   use ELSI_UTILS

   implicit none
   private

   public :: elsi_blacs_to_pexsi_dm
   public :: elsi_blacs_to_pexsi_hs
   public :: elsi_blacs_to_sips_hs
   public :: elsi_pexsi_to_blacs_dm
   public :: elsi_pexsi_to_blacs_hs

contains

!>
!! This routine is a driver to convert matrix format and distribution
!! from BLACS to PEXSI.
!!
subroutine elsi_blacs_to_pexsi_hs(elsi_h,H_in,S_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     intent(inout) :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian matrix to be converted
   real(kind=r8),     intent(inout) :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap matrix to be converted

   character*40, parameter :: caller = "elsi_blacs_to_pexsi_hs"

   if(elsi_h%ovlp_is_unit) then
      !TODO
      call elsi_stop(" PEXSI with identity overlap matrix not yet available."//&
              " Exiting...",elsi_h,caller)
   else
      call elsi_set_full_mat(elsi_h,H_in)
      if(elsi_h%n_elsi_calls == 1) then
         call elsi_set_full_mat(elsi_h,S_in)
      endif

      if(elsi_h%n_basis < 46340) then
         call elsi_blacs_to_pexsi_hs_small(elsi_h,H_in,S_in)
      else ! use long integer
         call elsi_blacs_to_pexsi_hs_large(elsi_h,H_in,S_in)
      endif
   endif

end subroutine

!>
!! This routine is a driver to convert matrix format and distribution
!! from PEXSI to BLACS.
!!
subroutine elsi_pexsi_to_blacs_dm(elsi_h,D_out)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                 !< Handle
   real(kind=r8),     intent(out)   :: D_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix to be converted

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_dm"

   if(elsi_h%n_basis < 46340) then
      call elsi_pexsi_to_blacs_dm_small(elsi_h,D_out)
   else ! use long integer
      call elsi_pexsi_to_blacs_dm_large(elsi_h,D_out)
   endif

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by PEXSI.
!!
!! * PEXSI pole parallelism MUST be available. Data redistributed on
!!   the processes corresponding to the first pole.
!!
subroutine elsi_blacs_to_pexsi_hs_small(elsi_h,H_in,S_in)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     intent(in)    :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian matrix to be converted
   real(kind=r8),     intent(in)    :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap matrix to be converted

   integer(kind=i4) :: n_parallel_groups
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: global_col_id ! Global column id
   integer(kind=i4) :: global_row_id ! Global row id
   integer(kind=i4) :: d1,d2,d11,d12,d21,d22 ! Number of columns in the intermediate stage
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: this_n_cols
   integer(kind=i4) :: tmp_int
   integer(kind=i4) :: min_pos
   integer(kind=i4) :: min_id
   integer(kind=i4) :: nnz_l_pexsi_aux
   integer(kind=i4) :: mpi_comm_aux_pexsi
   integer(kind=i4), allocatable :: locat(:) ! Location of each global column
   real(kind=r8) :: tmp_real

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send_buffer(:) ! Send buffer for Hamiltonian
   real(kind=r8),    allocatable :: s_val_send_buffer(:) ! Send buffer for overlap
   integer(kind=i4), allocatable :: pos_send_buffer(:) ! Send buffer for global 1D id
   integer(kind=i4), allocatable :: send_count(:) ! Number of elements to send to each processor
   integer(kind=i4), allocatable :: send_displ(:) ! Displacement from which to take the outgoing data
   real(kind=r8),    allocatable :: h_val_recv_buffer(:) ! Receive buffer for Hamiltonian
   real(kind=r8),    allocatable :: s_val_recv_buffer(:) ! Receive buffer for overlap
   integer(kind=i4), allocatable :: pos_recv_buffer(:) ! Receive buffer for global 1D id
   integer(kind=i4), allocatable :: recv_count(:) ! Number of elements to receive from each processor
   integer(kind=i4), allocatable :: recv_displ(:) ! Displacement at which to place the incoming data
   integer(kind=i4) :: send_displ_aux ! Auxiliary variable used to set displacement
   integer(kind=i4) :: recv_displ_aux ! Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_blacs_to_pexsi_hs_small"

   call elsi_start_redistribution_time(elsi_h)

   n_parallel_groups = elsi_h%n_procs/elsi_h%n_p_per_pole

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz(elsi_h,S_in,elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%nnz_l)
      call elsi_allocate(elsi_h,s_val_send_buffer,elsi_h%nnz_l,"s_val_send_buffer",caller)
   endif

   call elsi_allocate(elsi_h,locat,elsi_h%n_basis,"locat",caller)
   call elsi_allocate(elsi_h,pos_send_buffer,elsi_h%nnz_l,"pos_send_buffer",caller)
   call elsi_allocate(elsi_h,h_val_send_buffer,elsi_h%nnz_l,"h_val_send_buffer",caller)

   ! Compute d1,d2,d11,d12,d21,d22 (need explanation)
   d1  = elsi_h%n_basis/elsi_h%n_p_per_pole
   d2  = elsi_h%n_basis-(elsi_h%n_p_per_pole-1)*d1
   d11 = d1/n_parallel_groups
   d12 = d1-(n_parallel_groups-1)*d11
   d21 = d2/n_parallel_groups
   d22 = d2-(n_parallel_groups-1)*d21

   i_val = 0

   if(d11 > 0) then
      do i_proc = 0,elsi_h%n_procs-n_parallel_groups-1
         if(mod((i_proc+1),n_parallel_groups) == 0) then
            this_n_cols = d12
         else
            this_n_cols = d11
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = 0,elsi_h%n_procs-n_parallel_groups-1
         if(1+mod(i_proc,n_parallel_groups) .le. d1) then
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
      do i_proc = elsi_h%n_procs-n_parallel_groups,elsi_h%n_procs-1
         if(mod((i_proc+1),n_parallel_groups) == 0) then
            this_n_cols = d22
         else
            this_n_cols = d21
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = elsi_h%n_procs-n_parallel_groups,elsi_h%n_procs-1
         if(1+mod(i_proc,n_parallel_groups) .le. d2) then
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

   ! Compute destination and global 1D id
   if(elsi_h%n_elsi_calls == 1) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(S_in(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest = min(locat(global_col_id),elsi_h%n_procs-1)
               ! Pack global id and data into buffers
               pos_send_buffer(i_val) = (global_col_id-1)*elsi_h%n_basis+global_row_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               s_val_send_buffer(i_val) = S_in(i_row,i_col)
               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   else
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(S_in(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest = min(locat(global_col_id),elsi_h%n_procs-1)
               ! Pack global id and data into buffers
               pos_send_buffer(i_val) = (global_col_id-1)*elsi_h%n_basis+global_row_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   endif

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   nnz_l_pexsi_aux = sum(recv_count,1)
   call MPI_Allreduce(nnz_l_pexsi_aux,elsi_h%nnz_g,1,mpi_integer,&
           mpi_sum,elsi_h%mpi_comm,mpierr)

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
   ! Position
   call elsi_allocate(elsi_h,pos_recv_buffer,nnz_l_pexsi_aux,"pos_recv_buffer",caller)

   call MPI_Alltoallv(pos_send_buffer,send_count,send_displ,mpi_integer,&
           pos_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,pos_send_buffer,"pos_send_buffer")

   ! Hamiltonian value
   call elsi_allocate(elsi_h,h_val_recv_buffer,nnz_l_pexsi_aux,"h_val_recv_buffer",caller)

   call MPI_Alltoallv(h_val_send_buffer,send_count,send_displ,mpi_real8,&
           h_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_send_buffer,"h_val_send_buffer")

   ! Overlap value
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,s_val_recv_buffer,nnz_l_pexsi_aux,"s_val_recv_buffer",caller)

      call MPI_Alltoallv(s_val_send_buffer,send_count,send_displ,mpi_real8,&
              s_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,s_val_send_buffer,"s_val_send_buffer")
   endif

   ! Sort
   if(elsi_h%n_elsi_calls == 1) then
      do i_val = 1,nnz_l_pexsi_aux
         min_id = minloc(pos_recv_buffer(i_val:nnz_l_pexsi_aux),1)+i_val-1

         tmp_int = pos_recv_buffer(i_val)
         pos_recv_buffer(i_val) = pos_recv_buffer(min_id)
         pos_recv_buffer(min_id) = tmp_int

         tmp_real = h_val_recv_buffer(i_val)
         h_val_recv_buffer(i_val) = h_val_recv_buffer(min_id)
         h_val_recv_buffer(min_id) = tmp_real

         tmp_real = s_val_recv_buffer(i_val)
         s_val_recv_buffer(i_val) = s_val_recv_buffer(min_id)
         s_val_recv_buffer(min_id) = tmp_real
      enddo
   else
      do i_val = 1,nnz_l_pexsi_aux
         min_id = minloc(pos_recv_buffer(i_val:nnz_l_pexsi_aux),1)+i_val-1

         tmp_int = pos_recv_buffer(i_val)
         pos_recv_buffer(i_val) = pos_recv_buffer(min_id)
         pos_recv_buffer(min_id) = tmp_int

         tmp_real = h_val_recv_buffer(i_val)
         h_val_recv_buffer(i_val) = h_val_recv_buffer(min_id)
         h_val_recv_buffer(min_id) = tmp_real
      enddo
   endif

   ! Set send_count, all data sent to the first pole
   send_count = 0
   send_count(elsi_h%myid/n_parallel_groups+1) = nnz_l_pexsi_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

   if(elsi_h%n_elsi_calls == 1) then
      elsi_h%nnz_l_pexsi = sum(recv_count,1)

      ! At this point only the first pole knows nnz_l_pexsi
      call MPI_Comm_split(elsi_h%mpi_comm,elsi_h%my_p_col_pexsi,elsi_h%my_p_row_pexsi,&
              mpi_comm_aux_pexsi,mpierr)

      call MPI_Bcast(elsi_h%nnz_l_pexsi,1,mpi_integer,0,mpi_comm_aux_pexsi,mpierr)
   endif

   ! Set send and receive displacement
   send_displ_aux = 0
   recv_displ_aux = 0

   do i_proc = 1,elsi_h%n_procs
      send_displ(i_proc) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc)
      recv_displ(i_proc) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc)
   enddo

   ! Allocate PEXSI matrices
   ! Only the Hamiltonian needs to be reset everytime
   if(.not. allocated(elsi_h%ham_real_pexsi)) &
      call elsi_allocate(elsi_h,elsi_h%ham_real_pexsi,elsi_h%nnz_l_pexsi,"ham_real_pexsi",caller)
   elsi_h%ham_real_pexsi = 0.0_r8

   if(.not. allocated(elsi_h%ovlp_real_pexsi)) &
      call elsi_allocate(elsi_h,elsi_h%ovlp_real_pexsi,elsi_h%nnz_l_pexsi,"ovlp_real_pexsi",caller)

   if(.not. allocated(elsi_h%row_ind_pexsi)) &
      call elsi_allocate(elsi_h,elsi_h%row_ind_pexsi,elsi_h%nnz_l_pexsi,"row_ind_pexsi",caller)

   if(.not. allocated(elsi_h%col_ptr_pexsi)) &
      call elsi_allocate(elsi_h,elsi_h%col_ptr_pexsi,(elsi_h%n_l_cols_pexsi+1),"col_ptr_pexsi",caller)

   ! Send and receive the packed data
   if(elsi_h%n_elsi_calls == 1) then
      ! Position
      call elsi_allocate(elsi_h,pos_send_buffer,elsi_h%nnz_l_pexsi,"pos_send_buffer",caller)

      call MPI_Alltoallv(pos_recv_buffer,send_count,send_displ,mpi_integer,&
              pos_send_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,pos_recv_buffer,"pos_recv_buffer")

      ! Overlap value
      call MPI_Alltoallv(s_val_recv_buffer,send_count,send_displ,mpi_real8,&
              elsi_h%ovlp_real_pexsi,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,s_val_recv_buffer,"pos_recv_buffer")
   endif

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv_buffer,send_count,send_displ,mpi_real8,&
           elsi_h%ham_real_pexsi,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_recv_buffer,"h_val_recv_buffer")
   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   ! Only the first pole computes row index and column pointer
   if(elsi_h%n_elsi_calls == 1) then
      if(elsi_h%my_p_row_pexsi == 0) then
         ! Compute row index and column pointer
         i_col = (pos_send_buffer(1)-1)/elsi_h%n_basis
         do i_val = 1,elsi_h%nnz_l_pexsi
            elsi_h%row_ind_pexsi(i_val) = mod(pos_send_buffer(i_val)-1,elsi_h%n_basis)+1
            if((pos_send_buffer(i_val)-1)/elsi_h%n_basis+1 > i_col) then
               i_col = i_col+1
               elsi_h%col_ptr_pexsi(i_col-(pos_send_buffer(1)-1)/elsi_h%n_basis) = i_val
            endif
         enddo

         elsi_h%col_ptr_pexsi(elsi_h%n_l_cols_pexsi+1) = elsi_h%nnz_l_pexsi+1
      endif

      call elsi_deallocate(elsi_h,pos_send_buffer,"pos_send_buffer")
   endif

   call elsi_stop_redistribution_time(elsi_h)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by PEXSI.
!!
!! * Integer(kind=8) is used to deal with large matrices.
!!
!! * PEXSI pole parallelism MUST be available. Data redistributed on
!!   the processes corresponding to the first pole.
!!
subroutine elsi_blacs_to_pexsi_hs_large(elsi_h,H_in,S_in)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     intent(in)    :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian matrix to be converted
   real(kind=r8),     intent(in)    :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap matrix to be converted

   integer(kind=i4) :: n_parallel_groups
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: global_col_id ! Global column id
   integer(kind=i4) :: global_row_id ! Global row id
   integer(kind=i4) :: d1,d2,d11,d12,d21,d22 ! Number of columns in the intermediate stage
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: this_n_cols
   integer(kind=i4) :: min_pos
   integer(kind=i4) :: min_id
   integer(kind=i4) :: nnz_l_pexsi_aux
   integer(kind=i4) :: mpi_comm_aux_pexsi
   integer(kind=i4) :: tmp_int
   integer(kind=i8) :: tmp_long
   integer(kind=i4), allocatable :: locat(:) ! Location of each global column
   integer(kind=i8), allocatable :: global_id(:) ! Global 1D id
   real(kind=r8) :: tmp_real

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send_buffer(:) ! Send buffer for Hamiltonian
   real(kind=r8),    allocatable :: s_val_send_buffer(:) ! Send buffer for overlap
   integer(kind=i4), allocatable :: row_send_buffer(:) ! Send buffer for global row id
   integer(kind=i4), allocatable :: col_send_buffer(:) ! Send buffer for global column id
   integer(kind=i4), allocatable :: send_count(:) ! Number of elements to send to each processor
   integer(kind=i4), allocatable :: send_displ(:) ! Displacement from which to take the outgoing data
   real(kind=r8),    allocatable :: h_val_recv_buffer(:) ! Receive buffer for Hamiltonian
   real(kind=r8),    allocatable :: s_val_recv_buffer(:) ! Receive buffer for overlap
   integer(kind=i4), allocatable :: row_recv_buffer(:) ! Receive buffer for global row id
   integer(kind=i4), allocatable :: col_recv_buffer(:) ! Receive buffer for global column id
   integer(kind=i4), allocatable :: recv_count(:) ! Number of elements to receive from each processor
   integer(kind=i4), allocatable :: recv_displ(:) ! Displacement at which to place the incoming data
   integer(kind=i4) :: recv_displ_aux ! Auxiliary variable used to set displacement
   integer(kind=i4) :: send_displ_aux ! Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_blacs_to_pexsi_hs_large"

   call elsi_start_redistribution_time(elsi_h)

   n_parallel_groups = elsi_h%n_procs/elsi_h%n_p_per_pole

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz(elsi_h,S_in,elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%nnz_l)
      call elsi_allocate(elsi_h,s_val_send_buffer,elsi_h%nnz_l,"s_val_send_buffer",caller)
   endif

   call elsi_allocate(elsi_h,locat,elsi_h%n_basis,"locat",caller)
   call elsi_allocate(elsi_h,row_send_buffer,elsi_h%nnz_l,"row_send_buffer",caller)
   call elsi_allocate(elsi_h,col_send_buffer,elsi_h%nnz_l,"col_send_buffer",caller)
   call elsi_allocate(elsi_h,h_val_send_buffer,elsi_h%nnz_l,"h_val_send_buffer",caller)

   ! Compute d1,d2,d11,d12,d21,d22 (need explanation)
   d1  = elsi_h%n_basis/elsi_h%n_p_per_pole
   d2  = elsi_h%n_basis-(elsi_h%n_p_per_pole-1)*d1
   d11 = d1/n_parallel_groups
   d12 = d1-(n_parallel_groups-1)*d11
   d21 = d2/n_parallel_groups
   d22 = d2-(n_parallel_groups-1)*d21

   i_val = 0

   if(d11 > 0) then
      do i_proc = 0,elsi_h%n_procs-n_parallel_groups-1
         if(mod((i_proc+1),n_parallel_groups) == 0) then
            this_n_cols = d12
         else
            this_n_cols = d11
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = 0,elsi_h%n_procs-n_parallel_groups-1
         if(1+mod(i_proc,n_parallel_groups) .le. d1) then
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
      do i_proc = elsi_h%n_procs-n_parallel_groups,elsi_h%n_procs-1
         if(mod((i_proc+1),n_parallel_groups) == 0) then
            this_n_cols = d22
         else
            this_n_cols = d21
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = elsi_h%n_procs-n_parallel_groups,elsi_h%n_procs-1
         if(1+mod(i_proc,n_parallel_groups) .le. d2) then
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
   if(elsi_h%n_elsi_calls == 1) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(S_in(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest = min(locat(global_col_id),elsi_h%n_procs-1)
               ! Pack global id and data into buffers
               row_send_buffer(i_val) = global_row_id
               col_send_buffer(i_val) = global_col_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               s_val_send_buffer(i_val) = S_in(i_row,i_col)
               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   else
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(S_in(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest = min(locat(global_col_id),elsi_h%n_procs-1)
               ! Pack global id and data into buffers
               row_send_buffer(i_val) = global_row_id
               col_send_buffer(i_val) = global_col_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
              ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   endif

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   nnz_l_pexsi_aux = sum(recv_count,1)
   call MPI_Allreduce(nnz_l_pexsi_aux,elsi_h%nnz_g,1,mpi_integer,&
           mpi_sum,elsi_h%mpi_comm,mpierr)

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
   ! Row id
   call elsi_allocate(elsi_h,row_recv_buffer,nnz_l_pexsi_aux,"row_recv_buffer",caller)

   call MPI_Alltoallv(row_send_buffer,send_count,send_displ,mpi_integer,&
           row_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buffer,"row_send_buffer")

   ! Column id
   call elsi_allocate(elsi_h,col_recv_buffer,nnz_l_pexsi_aux,"col_recv_buffer",caller)

   call MPI_Alltoallv(col_send_buffer,send_count,send_displ,mpi_integer,&
           col_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buffer,"col_send_buffer")

   ! Hamiltonian value
   call elsi_allocate(elsi_h,h_val_recv_buffer,nnz_l_pexsi_aux,"h_val_recv_buffer",caller)

   call MPI_Alltoallv(h_val_send_buffer,send_count,send_displ,mpi_real8,&
           h_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_send_buffer,"h_val_send_buffer")

   ! Overlap value
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,s_val_recv_buffer,nnz_l_pexsi_aux,"s_val_recv_buffer",caller)

      call MPI_Alltoallv(s_val_send_buffer,send_count,send_displ,mpi_real8,&
              s_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,s_val_send_buffer,"s_val_send_buffer")
   endif

   allocate(global_id(nnz_l_pexsi_aux))

   ! Compute global 1D id
   do i_val = 1,nnz_l_pexsi_aux
      global_id(i_val) = int(col_recv_buffer(i_val)-1,kind=i8)*int(elsi_h%n_basis,kind=i8)+&
                            int(row_recv_buffer(i_val),kind=i8)
   enddo

   ! Sort
   if(elsi_h%n_elsi_calls == 1) then
      do i_val = 1,nnz_l_pexsi_aux
         min_id = minloc(global_id(i_val:nnz_l_pexsi_aux),1)+i_val-1

         tmp_long = global_id(i_val)
         global_id(i_val) = global_id(min_id)
         global_id(min_id) = tmp_long

         tmp_int = row_recv_buffer(i_val)
         row_recv_buffer(i_val) = row_recv_buffer(min_id)
         row_recv_buffer(min_id) = tmp_int

         tmp_int = col_recv_buffer(i_val)
         col_recv_buffer(i_val) = col_recv_buffer(min_id)
         col_recv_buffer(min_id) = tmp_int

         tmp_real = h_val_recv_buffer(i_val)
         h_val_recv_buffer(i_val) = h_val_recv_buffer(min_id)
         h_val_recv_buffer(min_id) = tmp_real

         tmp_real = s_val_recv_buffer(i_val)
         s_val_recv_buffer(i_val) = s_val_recv_buffer(min_id)
         s_val_recv_buffer(min_id) = tmp_real
      enddo
   else ! Row and column id not needed
      do i_val = 1,nnz_l_pexsi_aux
         min_id = minloc(global_id(i_val:nnz_l_pexsi_aux),1)+i_val-1

         tmp_long = global_id(i_val)
         global_id(i_val) = global_id(min_id)
         global_id(min_id) = tmp_long

         tmp_real = h_val_recv_buffer(i_val)
         h_val_recv_buffer(i_val) = h_val_recv_buffer(min_id)
         h_val_recv_buffer(min_id) = tmp_real
      enddo
   endif

   deallocate(global_id)

   ! Set send_count, all data sent to the first pole
   send_count = 0
   send_count(elsi_h%myid/n_parallel_groups+1) = nnz_l_pexsi_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

   if(elsi_h%n_elsi_calls == 1) then
      elsi_h%nnz_l_pexsi = sum(recv_count,1)

      ! At this point only the first pole knows nnz_l_pexsi
      call MPI_Comm_split(elsi_h%mpi_comm,elsi_h%my_p_col_pexsi,elsi_h%my_p_row_pexsi,&
              mpi_comm_aux_pexsi,mpierr)

      call MPI_Bcast(elsi_h%nnz_l_pexsi,1,mpi_integer,0,mpi_comm_aux_pexsi,mpierr)
   endif

   ! Set send and receive displacement
   send_displ_aux = 0
   recv_displ_aux = 0

   do i_proc = 1,elsi_h%n_procs
      send_displ(i_proc) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc)
      recv_displ(i_proc) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc)
   enddo

   ! Allocate PEXSI matrices
   ! Only the Hamiltonian needs to be reset everytime
   if(.not. allocated(elsi_h%ham_real_pexsi)) &
      call elsi_allocate(elsi_h,elsi_h%ham_real_pexsi,elsi_h%nnz_l_pexsi,"ham_real_pexsi",caller)
   elsi_h%ham_real_pexsi = 0.0_r8

   if(.not. allocated(elsi_h%ovlp_real_pexsi)) &
      call elsi_allocate(elsi_h,elsi_h%ovlp_real_pexsi,elsi_h%nnz_l_pexsi,"ovlp_real_pexsi",caller)

   if(.not. allocated(elsi_h%row_ind_pexsi)) &
      call elsi_allocate(elsi_h,elsi_h%row_ind_pexsi,elsi_h%nnz_l_pexsi,"row_ind_pexsi",caller)

   if(.not. allocated(elsi_h%col_ptr_pexsi)) &
      call elsi_allocate(elsi_h,elsi_h%col_ptr_pexsi,(elsi_h%n_l_cols_pexsi+1),"col_ptr_pexsi",caller)

   ! Send and receive the packed data
   if(elsi_h%n_elsi_calls == 1) then
      ! Row id
      call elsi_allocate(elsi_h,row_send_buffer,elsi_h%nnz_l_pexsi,"row_send_buffer",caller)

      call MPI_Alltoallv(row_recv_buffer,send_count,send_displ,mpi_integer,&
              row_send_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,row_recv_buffer,"row_recv_buffer")

      ! Column id
      call elsi_allocate(elsi_h,col_send_buffer,elsi_h%nnz_l_pexsi,"col_send_buffer",caller)

      call MPI_Alltoallv(col_recv_buffer,send_count,send_displ,mpi_integer,&
              col_send_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,col_recv_buffer,"col_recv_buffer")

      ! Overlap value
      call MPI_Alltoallv(s_val_recv_buffer,send_count,send_displ,mpi_real8,&
              elsi_h%ovlp_real_pexsi,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,s_val_recv_buffer,"s_val_recv_buffer")
   endif

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv_buffer,send_count,send_displ,mpi_real8,&
           elsi_h%ham_real_pexsi,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_recv_buffer,"h_val_recv_buffer")
   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   ! Only the first pole computes row index and column pointer
   if(elsi_h%n_elsi_calls == 1) then
      if(elsi_h%my_p_row_pexsi == 0) then
         ! Compute row index and column pointer
         i_col = col_send_buffer(1)-1
         do i_val = 1,elsi_h%nnz_l_pexsi
            elsi_h%row_ind_pexsi(i_val) = row_send_buffer(i_val)
            if(col_send_buffer(i_val) > i_col) then
               i_col = i_col+1
               elsi_h%col_ptr_pexsi(i_col-col_send_buffer(1)+1) = i_val
            endif
         enddo

         elsi_h%col_ptr_pexsi(elsi_h%n_l_cols_pexsi+1) = elsi_h%nnz_l_pexsi+1
      endif

      call elsi_deallocate(elsi_h,row_send_buffer,"row_send_buffer")
      call elsi_deallocate(elsi_h,col_send_buffer,"col_send_buffer")
   endif

   call elsi_stop_redistribution_time(elsi_h)

end subroutine

!>
!! This routine converts density matrix computed by PEXSI and stored
!! in 1D block distributed sparse CCS format to 2D block-cyclic
!! distributed dense format.
!!
subroutine elsi_pexsi_to_blacs_dm_small(elsi_h,D_out)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h                                 !< Handle
   real(kind=r8),     intent(out)   :: D_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix to be converted

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: k_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: global_col_id ! Global column id
   integer(kind=i4) :: global_row_id ! Global row id
   integer(kind=i4) :: local_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: local_row_id ! Local row id in 1D block distribution
   integer(kind=i4) :: proc_col_id ! Column id in process grid
   integer(kind=i4) :: proc_row_id ! Row id in process grid
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: global_id(:) ! Global 1d id

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: val_send_buffer(:) ! Send buffer for value
   integer(kind=i4), allocatable :: pos_send_buffer(:) ! Send buffer for global 1D id
   integer(kind=i4), allocatable :: send_count(:) ! Number of elements to send to each processor
   integer(kind=i4), allocatable :: send_displ(:) ! Displacement from which to take the outgoing data
   real(kind=r8),    allocatable :: val_recv_buffer(:) ! Receive buffer for value
   integer(kind=i4), allocatable :: pos_recv_buffer(:) ! Receive buffer for global 1D id
   integer(kind=i4), allocatable :: recv_count(:) ! Number of elements to receive from each processor
   integer(kind=i4), allocatable :: recv_displ(:) ! Displacement at which to place the incoming data
   integer(kind=i4) :: send_displ_aux ! Auxiliary variable used to set displacement
   integer(kind=i4) :: recv_displ_aux ! Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_dm_small"

   call elsi_start_redistribution_time(elsi_h)

   call elsi_allocate(elsi_h,val_send_buffer,elsi_h%nnz_l_pexsi,"val_send_buffer",caller)
   call elsi_allocate(elsi_h,pos_send_buffer,elsi_h%nnz_l_pexsi,"pos_send_buffer",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   if(elsi_h%my_p_row_pexsi == 0) then
      call elsi_allocate(elsi_h,global_id,elsi_h%nnz_l_pexsi,"global_id",caller)
      call elsi_allocate(elsi_h,dest,elsi_h%nnz_l_pexsi,"dest",caller)

      i_col = 0
      ! Compute destination and global 1D id
      do i_val = 1,elsi_h%nnz_l_pexsi
         if(i_val == elsi_h%col_ptr_ccs(i_col+1) .and. i_col /= elsi_h%n_l_cols_pexsi) then
            i_col = i_col+1
         endif
         i_row = elsi_h%row_ind_ccs(i_val)

         ! Compute global id
         global_row_id = i_row
         global_col_id = i_col+elsi_h%myid*(elsi_h%n_basis/elsi_h%n_p_per_pole)
         global_id(i_val) = (global_col_id-1)*elsi_h%n_basis+global_row_id

         ! Compute destination
         proc_row_id = mod((global_row_id-1)/elsi_h%n_b_rows,elsi_h%n_p_rows)
         proc_col_id = mod((global_col_id-1)/elsi_h%n_b_cols,elsi_h%n_p_cols)
         dest(i_val) = proc_col_id+proc_row_id*elsi_h%n_p_cols
      enddo

      j_val = 0
      k_val = elsi_h%nnz_l_pexsi+1

      ! Set send_count
      do i_proc = 1,elsi_h%n_procs/2
         do i_val = 1,elsi_h%nnz_l_pexsi
            if(dest(i_val) == i_proc-1) then
               j_val = j_val+1
               val_send_buffer(j_val) = elsi_h%den_mat_ccs(i_val)
               pos_send_buffer(j_val) = global_id(i_val)
               send_count(i_proc) = send_count(i_proc)+1
            endif
            if(dest(i_val) == elsi_h%n_procs-i_proc) then
               k_val = k_val-1
               val_send_buffer(k_val) = elsi_h%den_mat_ccs(i_val)
               pos_send_buffer(k_val) = global_id(i_val)
               send_count(elsi_h%n_procs+1-i_proc) = &
                  send_count(elsi_h%n_procs+1-i_proc)+1
            endif
         enddo
      enddo

      call elsi_deallocate(elsi_h,global_id,"global_id")
      call elsi_deallocate(elsi_h,dest,"dest")
   endif

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

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
   call elsi_allocate(elsi_h,val_recv_buffer,elsi_h%nnz_l,"val_recv_buffer",caller)

   call MPI_Alltoallv(val_send_buffer,send_count,send_displ,mpi_real8,&
           val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,val_send_buffer,"val_send_buffer")

   ! Position
   call elsi_allocate(elsi_h,pos_recv_buffer,elsi_h%nnz_l,"pos_recv_buffer",caller)

   call MPI_Alltoallv(pos_send_buffer,send_count,send_displ,mpi_integer,&
           pos_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,pos_send_buffer,"pos_send_buffer")
   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   D_out = 0.0_r8

   ! Unpack density matrix
   do i_val = 1,elsi_h%nnz_l
      ! Compute global 2d id
      global_col_id = (pos_recv_buffer(i_val)-1)/elsi_h%n_basis+1
      global_row_id = mod(pos_recv_buffer(i_val)-1,elsi_h%n_basis)+1

      ! Compute local 2d id
      local_row_id = (global_row_id-1)/(elsi_h%n_p_rows*elsi_h%n_b_rows)*elsi_h%n_b_rows&
                        +mod((global_row_id-1),elsi_h%n_b_rows)+1
      local_col_id = (global_col_id-1)/(elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols&
                        +mod((global_col_id-1),elsi_h%n_b_cols)+1

      ! Put value to correct position
      D_out(local_row_id,local_col_id) = val_recv_buffer(i_val)
   enddo

   call elsi_deallocate(elsi_h,val_recv_buffer,"val_recv_buffer")
   call elsi_deallocate(elsi_h,pos_recv_buffer,"pos_recv_buffer")

   call elsi_stop_redistribution_time(elsi_h)

end subroutine

!>
!! This routine converts density matrix computed by PEXSI and stored
!! in 1D block distributed sparse CCS format to 2D block-cyclic
!! distributed dense format.
!!
!! * 2D index is used to deal with large matrices.
!!
subroutine elsi_pexsi_to_blacs_dm_large(elsi_h,D_out)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h                                 !< Handle
   real(kind=r8),     intent(out)   :: D_out(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Density matrix to be converted

   integer(kind=i4) :: mpierr
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
   real(kind=r8),    allocatable :: val_send_buffer(:)  ! Send buffer for value
   integer(kind=i4), allocatable :: row_send_buffer(:) ! Send buffer for global row id
   integer(kind=i4), allocatable :: col_send_buffer(:) ! Send buffer for global column id
   integer(kind=i4), allocatable :: send_count(:) ! Number of elements to send to each processor
   integer(kind=i4), allocatable :: send_displ(:) ! Displacement from which to take the outgoing data
   real(kind=r8),    allocatable :: val_recv_buffer(:) ! Receive buffer for value
   integer(kind=i4), allocatable :: row_recv_buffer(:) ! Receive buffer for global row id
   integer(kind=i4), allocatable :: col_recv_buffer(:) ! Receive buffer for global column id
   integer(kind=i4), allocatable :: recv_count(:) ! Number of elements to receive from each processor
   integer(kind=i4), allocatable :: recv_displ(:) ! Displacement at which to place the incoming data
   integer(kind=i4) :: send_displ_aux ! Auxiliary variable used to set displacement
   integer(kind=i4) :: recv_displ_aux ! Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_dm_large"

   call elsi_start_redistribution_time(elsi_h)

   call elsi_allocate(elsi_h,val_send_buffer,elsi_h%nnz_l_pexsi,"val_send_buffer",caller)
   call elsi_allocate(elsi_h,row_send_buffer,elsi_h%nnz_l_pexsi,"row_send_buffer",caller)
   call elsi_allocate(elsi_h,col_send_buffer,elsi_h%nnz_l_pexsi,"col_send_buffer",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   if(elsi_h%my_p_row_pexsi == 0) then
      call elsi_allocate(elsi_h,global_row_id,elsi_h%nnz_l_pexsi,"global_row_id",caller)
      call elsi_allocate(elsi_h,global_col_id,elsi_h%nnz_l_pexsi,"global_col_id",caller)
      call elsi_allocate(elsi_h,dest,elsi_h%nnz_l_pexsi,"dest",caller)

      i_col = 0
      ! Compute destination and global id
      do i_val = 1,elsi_h%nnz_l_pexsi
         if(i_val == elsi_h%col_ptr_ccs(i_col+1) .and. i_col /= elsi_h%n_l_cols_pexsi) then
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
      k_val = elsi_h%nnz_l_pexsi+1

      ! Set send_count
      do i_proc = 1,elsi_h%n_procs/2
         do i_val = 1,elsi_h%nnz_l_pexsi
            if(dest(i_val) == i_proc-1) then
               j_val = j_val+1
               val_send_buffer(j_val) = elsi_h%den_mat_ccs(i_val)
               row_send_buffer(j_val) = global_row_id(i_val)
               col_send_buffer(j_val) = global_col_id(i_val)
               send_count(i_proc) = send_count(i_proc)+1
            endif
            if(dest(i_val) == elsi_h%n_procs-i_proc) then
               k_val = k_val-1
               val_send_buffer(k_val) = elsi_h%den_mat_ccs(i_val)
               row_send_buffer(k_val) = global_row_id(i_val)
               col_send_buffer(k_val) = global_col_id(i_val)
               send_count(elsi_h%n_procs+1-i_proc) = &
                  send_count(elsi_h%n_procs+1-i_proc)+1
            endif
         enddo
      enddo

      call elsi_deallocate(elsi_h,global_row_id,"global_row_id")
      call elsi_deallocate(elsi_h,global_col_id,"global_col_id")
      call elsi_deallocate(elsi_h,dest,"dest")
   endif

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

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
   call elsi_allocate(elsi_h,val_recv_buffer,elsi_h%nnz_l,"val_recv_buffer",caller)

   call MPI_Alltoallv(val_send_buffer,send_count,send_displ,mpi_real8,&
           val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,val_send_buffer,"val_send_buffer")

   ! Row index
   call elsi_allocate(elsi_h,row_recv_buffer,elsi_h%nnz_l,"row_recv_buffer",caller)

   call MPI_Alltoallv(row_send_buffer,send_count,send_displ,mpi_integer,&
           row_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buffer,"row_send_buffer")

   ! Column index
   call elsi_allocate(elsi_h,col_recv_buffer,elsi_h%nnz_l,"col_recv_buffer",caller)

   call MPI_Alltoallv(col_send_buffer,send_count,send_displ,mpi_integer,&
           col_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buffer,"col_send_buffer")
   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   D_out = 0.0_r8

   ! Unpack density matrix
   do i_val = 1,elsi_h%nnz_l
      ! Compute local 2d id
      local_row_id = (row_recv_buffer(i_val)-1)/(elsi_h%n_p_rows*elsi_h%n_b_rows)*elsi_h%n_b_rows&
                        +mod((row_recv_buffer(i_val)-1),elsi_h%n_b_rows)+1
      local_col_id = (col_recv_buffer(i_val)-1)/(elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols&
                        +mod((col_recv_buffer(i_val)-1),elsi_h%n_b_cols)+1

      ! Put value to correct position
      D_out(local_row_id,local_col_id) = val_recv_buffer(i_val)
   enddo

   call elsi_deallocate(elsi_h,val_recv_buffer,"val_recv_buffer")
   call elsi_deallocate(elsi_h,row_recv_buffer,"row_recv_buffer")
   call elsi_deallocate(elsi_h,col_recv_buffer,"col_recv_buffer")

   call elsi_stop_redistribution_time(elsi_h)

end subroutine

!>
!! This routine is a driver to convert matrix format and distribution
!! from BLACS to SIPs.
!!
subroutine elsi_blacs_to_sips_hs(elsi_h,H_in,S_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     intent(inout) :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian matrix to be converted
   real(kind=r8),     intent(inout) :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap matrix to be converted

   character*40, parameter :: caller = "elsi_blacs_to_sips_hs"

   if(elsi_h%ovlp_is_unit) then
      call elsi_stop(" SIPs with identity overlap matrix not yet available."//&
              " Exiting...",elsi_h,caller)
   else
      call elsi_set_full_mat(elsi_h,H_in)
      if(elsi_h%n_elsi_calls == 1) then
         call elsi_set_full_mat(elsi_h,S_in)
      endif

      if(elsi_h%n_basis < 46340) then ! kind=4 integer works
         call elsi_blacs_to_sips_hs_small(elsi_h,H_in,S_in)
      else ! use kind=8 integer
         call elsi_blacs_to_sips_hs_large(elsi_h,H_in,S_in)
      endif
   endif

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by SIPs.
!!
subroutine elsi_blacs_to_sips_hs_small(elsi_h,H_in,S_in)

   implicit none

   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     intent(in)    :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian matrix to be converted
   real(kind=r8),     intent(in)    :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap matrix to be converted

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: global_col_id ! Global column id
   integer(kind=i4) :: global_row_id ! Global row id
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: tmp_int
   integer(kind=i4) :: min_pos
   integer(kind=i4) :: min_id
   real(kind=r8)    :: tmp_real

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send_buffer(:) ! Send buffer for Hamiltonian
   real(kind=r8),    allocatable :: s_val_send_buffer(:) ! Send buffer for overlap
   integer(kind=i4), allocatable :: pos_send_buffer(:) ! Send buffer for global 1D id
   integer(kind=i4), allocatable :: send_count(:) ! Number of elements to send to each processor
   integer(kind=i4), allocatable :: send_displ(:) ! Displacement from which to take the outgoing data
   real(kind=r8),    allocatable :: h_val_recv_buffer(:) ! Receive buffer for Hamiltonian
   real(kind=r8),    allocatable :: s_val_recv_buffer(:) ! Receive buffer for overlap
   integer(kind=i4), allocatable :: pos_recv_buffer(:) ! Receive buffer for global 1D id
   integer(kind=i4), allocatable :: recv_count(:) ! Number of elements to receive from each processor
   integer(kind=i4), allocatable :: recv_displ(:) ! Displacement at which to place the incoming data
   integer(kind=i4) :: send_displ_aux ! Auxiliary variable used to set displacement 
   integer(kind=i4) :: recv_displ_aux ! Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_blacs_to_sips_hs_small"

   call elsi_start_redistribution_time(elsi_h)

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz(elsi_h,H_in,elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%nnz_l)
      call elsi_allocate(elsi_h,s_val_send_buffer,elsi_h%nnz_l,"s_val_send_buffer",caller)
   endif

   call elsi_allocate(elsi_h,pos_send_buffer,elsi_h%nnz_l,"pos_send_buffer",caller)
   call elsi_allocate(elsi_h,h_val_send_buffer,elsi_h%nnz_l,"h_val_send_buffer",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   ! Compute destination and global 1D id
   if(elsi_h%n_elsi_calls == 1) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(S_in(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest = (global_col_id-1)/(elsi_h%n_basis/elsi_h%n_procs)
               ! The last process may take more
               dest = min(dest,elsi_h%n_procs-1)
               ! Compute the global id
               ! Pack global id and data into buffers
               pos_send_buffer(i_val) = (global_col_id-1)*elsi_h%n_basis+global_row_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               s_val_send_buffer(i_val) = S_in(i_row,i_col)
               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
           endif
        enddo
      enddo
   else
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(S_in(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest = (global_col_id-1)/(elsi_h%n_basis/elsi_h%n_procs)
               ! The last process may take more
               dest = min(dest,elsi_h%n_procs-1)
               ! Compute the global id
               ! Pack global id and data into buffers
               pos_send_buffer(i_val) = (global_col_id-1)*elsi_h%n_basis+global_row_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
           endif
        enddo
      enddo
   endif

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   elsi_h%nnz_l_sips = sum(recv_count,1)
   call MPI_Allreduce(elsi_h%nnz_l_sips,elsi_h%nnz_g,1,mpi_integer,&
            mpi_sum,elsi_h%mpi_comm,mpierr)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   send_displ_aux = 0
   recv_displ_aux = 0

   do i_proc = 0,elsi_h%n_procs-1
      send_displ(i_proc+1) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc+1)
      recv_displ(i_proc+1) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
   enddo

   ! Send and receive the packed data
   ! Position
   call elsi_allocate(elsi_h,pos_recv_buffer,elsi_h%nnz_l_sips,"pos_recv_buffer",caller)

   call MPI_Alltoallv(pos_send_buffer,send_count,send_displ,mpi_integer,&
           pos_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,pos_send_buffer,"pos_send_buffer")

   ! Hamiltonian value
   call elsi_allocate(elsi_h,h_val_recv_buffer,elsi_h%nnz_l_sips,"h_val_recv_buffer",caller)

   call MPI_Alltoallv(h_val_send_buffer,send_count,send_displ,mpi_real8,&
           h_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_send_buffer,"h_val_send_buffer")

   ! Overlap value
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,s_val_recv_buffer,elsi_h%nnz_l_sips,"s_val_recv_buffer",caller)

      call MPI_Alltoallv(s_val_send_buffer,send_count,send_displ,mpi_real8,&
              s_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,s_val_send_buffer,"s_val_send_buffer")
   endif

   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   ! Sort
   if(elsi_h%n_elsi_calls == 1) then
      do i_val = 1,elsi_h%nnz_l_sips
         min_id = minloc(pos_recv_buffer(i_val:elsi_h%nnz_l_sips),1)+i_val-1

         tmp_int = pos_recv_buffer(i_val)
         pos_recv_buffer(i_val) = pos_recv_buffer(min_id)
         pos_recv_buffer(min_id) = tmp_int

         tmp_real = h_val_recv_buffer(i_val)
         h_val_recv_buffer(i_val) = h_val_recv_buffer(min_id)
         h_val_recv_buffer(min_id) = tmp_real

         tmp_real = s_val_recv_buffer(i_val)
         s_val_recv_buffer(i_val) = s_val_recv_buffer(min_id)
         s_val_recv_buffer(min_id) = tmp_real
      enddo
   else
      do i_val = 1,elsi_h%nnz_l_sips
         min_id = minloc(pos_recv_buffer(i_val:elsi_h%nnz_l_sips),1)+i_val-1

         tmp_int = pos_recv_buffer(i_val)
         pos_recv_buffer(i_val) = pos_recv_buffer(min_id)
         pos_recv_buffer(min_id) = tmp_int

         tmp_real = h_val_recv_buffer(i_val)
         h_val_recv_buffer(i_val) = h_val_recv_buffer(min_id)
         h_val_recv_buffer(min_id) = tmp_real
      enddo
   endif

   ! Allocate SIPs matrices
   ! Only the Hamiltonian needs to be reset everytime
   if(.not. allocated(elsi_h%ham_real_sips)) &
      call elsi_allocate(elsi_h,elsi_h%ham_real_sips,elsi_h%nnz_l_sips,"ham_real_sips",caller)
   elsi_h%ham_real_sips = 0.0_r8

   if(.not. allocated(elsi_h%ovlp_real_sips)) &
      call elsi_allocate(elsi_h,elsi_h%ovlp_real_sips,elsi_h%nnz_l_sips,"ovlp_real_sips",caller)

   if(.not. allocated(elsi_h%row_ind_sips)) &
      call elsi_allocate(elsi_h,elsi_h%row_ind_sips,elsi_h%nnz_l_sips,"row_ind_sips",caller)

   if(.not. allocated(elsi_h%col_ptr_sips)) &
      call elsi_allocate(elsi_h,elsi_h%col_ptr_sips,(elsi_h%n_l_cols_sips+1),"col_ptr_sips",caller)

   elsi_h%ham_real_sips = h_val_recv_buffer
   call elsi_deallocate(elsi_h,h_val_recv_buffer,"h_val_recv_buffer")

   if(elsi_h%n_elsi_calls == 1) then
      elsi_h%ovlp_real_sips = s_val_recv_buffer
      call elsi_deallocate(elsi_h,s_val_recv_buffer,"s_val_recv_buffer")
   endif

   ! Compute row index and column pointer
   if(elsi_h%n_elsi_calls == 1) then
      i_col = (pos_recv_buffer(1)-1)/elsi_h%n_basis
      do i_val = 1,elsi_h%nnz_l_sips
         elsi_h%row_ind_sips(i_val) = mod(pos_recv_buffer(i_val)-1,elsi_h%n_basis)+1
         if((pos_recv_buffer(i_val)-1)/elsi_h%n_basis+1 > i_col) then
            i_col = i_col+1
            elsi_h%col_ptr_sips(i_col-(pos_recv_buffer(1)-1)/elsi_h%n_basis) = i_val
         endif
      enddo

      elsi_h%col_ptr_sips(elsi_h%n_l_cols_sips+1) = elsi_h%nnz_l_sips+1
   endif

   call elsi_deallocate(elsi_h,pos_recv_buffer,"pos_recv_buffer")

   ! Dummy
   elsi_h%nnz_l_pexsi = elsi_h%nnz_l_sips
   elsi_h%n_l_cols_pexsi = elsi_h%n_l_cols_sips

   call elsi_stop_redistribution_time(elsi_h)

end subroutine

!> 
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by SIPs.
!!
!! * Integer(kind=8) is used to deal with large matrices.
!!
subroutine elsi_blacs_to_sips_hs_large(elsi_h,H_in,S_in)

   implicit none

   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h                                !< Handle
   real(kind=r8),     intent(in)    :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian matrix to be converted
   real(kind=r8),     intent(in)    :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap matrix to be converted

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: global_col_id ! Global column id
   integer(kind=i4) :: global_row_id ! Global row id
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: tmp_int
   integer(kind=i8) :: tmp_long
   integer(kind=i4) :: min_pos
   integer(kind=i4) :: min_id
   real(kind=r8)    :: tmp_real

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send_buffer(:) ! Send buffer for Hamiltonian
   real(kind=r8),    allocatable :: s_val_send_buffer(:) ! Send buffer for overlap
   integer(kind=i4), allocatable :: row_send_buffer(:) ! Send buffer for global row id
   integer(kind=i4), allocatable :: col_send_buffer(:) ! Send buffer for global column id
   integer(kind=i4), allocatable :: send_count(:) ! Number of elements to send to each processor
   integer(kind=i4), allocatable :: send_displ(:) ! Displacement from which to take the outgoing data
   real(kind=r8),    allocatable :: h_val_recv_buffer(:) ! Receive buffer for Hamiltonian
   real(kind=r8),    allocatable :: s_val_recv_buffer(:) ! Receive buffer for overlap
   integer(kind=i4), allocatable :: row_recv_buffer(:) ! Receive buffer for global row id
   integer(kind=i4), allocatable :: col_recv_buffer(:) ! Receive buffer for global col id
   integer(kind=i4), allocatable :: recv_count(:) ! Number of elements to receive from each processor
   integer(kind=i4), allocatable :: recv_displ(:) ! Displacement at which to place the incoming data
   integer(kind=i4) :: send_displ_aux ! Auxiliary variable used to set displacement 
   integer(kind=i4) :: recv_displ_aux ! Auxiliary variable used to set displacement
   integer(kind=i8), allocatable :: global_id(:) ! Global 1D id

   character*40, parameter :: caller = "elsi_blacs_to_sips_hs_large"

   call elsi_start_redistribution_time(elsi_h)

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz(elsi_h,H_in,elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%nnz_l)
      call elsi_allocate(elsi_h,s_val_send_buffer,elsi_h%nnz_l,"s_val_send_buffer",caller)
   endif

   call elsi_allocate(elsi_h,row_send_buffer,elsi_h%nnz_l,"row_send_buffer",caller)
   call elsi_allocate(elsi_h,col_send_buffer,elsi_h%nnz_l,"col_send_buffer",caller)
   call elsi_allocate(elsi_h,h_val_send_buffer,elsi_h%nnz_l,"h_val_send_buffer",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   if(elsi_h%n_elsi_calls == 1) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(S_in(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest = (global_col_id-1)/(elsi_h%n_basis/elsi_h%n_procs)
               ! The last process may take more
               dest = min(dest,elsi_h%n_procs-1)
               ! Pack global id and data into buffers
               row_send_buffer(i_val) = global_row_id
               col_send_buffer(i_val) = global_col_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               s_val_send_buffer(i_val) = S_in(i_row,i_col)
               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
           endif
        enddo
      enddo
   else
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(S_in(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest = (global_col_id-1)/(elsi_h%n_basis/elsi_h%n_procs)
               ! The last process may take more
               dest = min(dest,elsi_h%n_procs-1)
               ! Pack global id and data into buffers
               row_send_buffer(i_val) = global_row_id
               col_send_buffer(i_val) = global_col_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
           endif
        enddo
      enddo
   endif

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

   ! Set local/global number of nonzero
   elsi_h%nnz_l_sips = sum(recv_count,1)
   call MPI_Allreduce(elsi_h%nnz_l_sips,elsi_h%nnz_g,1,mpi_integer,&
            mpi_sum,elsi_h%mpi_comm,mpierr)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   send_displ_aux = 0
   recv_displ_aux = 0

   do i_proc = 0,elsi_h%n_procs-1
      send_displ(i_proc+1) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc+1)
      recv_displ(i_proc+1) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(elsi_h,row_recv_buffer,elsi_h%nnz_l_sips,"row_recv_buffer",caller)

   call MPI_Alltoallv(row_send_buffer,send_count,send_displ,mpi_integer,&
           row_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buffer,"row_send_buffer")

   ! Column id
   call elsi_allocate(elsi_h,col_recv_buffer,elsi_h%nnz_l_sips,"col_recv_buffer",caller)

   call MPI_Alltoallv(col_send_buffer,send_count,send_displ,mpi_integer,&
           col_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buffer,"col_send_buffer")

   ! Hamiltonian value
   call elsi_allocate(elsi_h,h_val_recv_buffer,elsi_h%nnz_l_sips,"h_val_recv_buffer",caller)

   call MPI_Alltoallv(h_val_send_buffer,send_count,send_displ,mpi_real8,&
           h_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_send_buffer,"h_val_send_buffer")

   ! Overlap value
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,s_val_recv_buffer,elsi_h%nnz_l_sips,"s_val_recv_buffer",caller)

      call MPI_Alltoallv(s_val_send_buffer,send_count,send_displ,mpi_real8,&
              s_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

      call elsi_deallocate(elsi_h,s_val_send_buffer,"s_val_send_buffer")
   endif

   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   allocate(global_id(elsi_h%nnz_l_sips))

   ! Compute global 1D id
   do i_val = 1,elsi_h%nnz_l_sips
      global_id(i_val) = int(col_recv_buffer(i_val)-1,kind=i8)*int(elsi_h%n_basis,kind=i8)+&
                            int(row_recv_buffer(i_val),kind=i8)
   enddo

   ! Sort
   if(elsi_h%n_elsi_calls == 1) then
      do i_val = 1,elsi_h%nnz_l_sips
         min_id = minloc(global_id(i_val:elsi_h%nnz_l_sips),1)+i_val-1

         tmp_long = global_id(i_val)
         global_id(i_val) = global_id(min_id)
         global_id(min_id) = tmp_long

         tmp_int = row_recv_buffer(i_val)
         row_recv_buffer(i_val) = row_recv_buffer(min_id)
         row_recv_buffer(min_id) = tmp_int

         tmp_int = col_recv_buffer(i_val)
         col_recv_buffer(i_val) = col_recv_buffer(min_id)
         col_recv_buffer(min_id) = tmp_int

         tmp_real = h_val_recv_buffer(i_val)
         h_val_recv_buffer(i_val) = h_val_recv_buffer(min_id)
         h_val_recv_buffer(min_id) = tmp_real

         tmp_real = s_val_recv_buffer(i_val)
         s_val_recv_buffer(i_val) = s_val_recv_buffer(min_id)
         s_val_recv_buffer(min_id) = tmp_real
      enddo
   else
      do i_val = 1,elsi_h%nnz_l_sips
         min_id = minloc(global_id(i_val:elsi_h%nnz_l_sips),1)+i_val-1

         tmp_long = global_id(i_val)
         global_id(i_val) = global_id(min_id)
         global_id(min_id) = tmp_long

         tmp_real = h_val_recv_buffer(i_val)
         h_val_recv_buffer(i_val) = h_val_recv_buffer(min_id)
         h_val_recv_buffer(min_id) = tmp_real
      enddo
   endif

   deallocate(global_id)

   ! Allocate SIPs matrices
   ! Only the Hamiltonian needs to be reset everytime
   if(.not. allocated(elsi_h%ham_real_sips)) &
      call elsi_allocate(elsi_h,elsi_h%ham_real_sips,elsi_h%nnz_l_sips,"ham_real_sips",caller)
   elsi_h%ham_real_sips = 0.0_r8

   if(.not. allocated(elsi_h%ovlp_real_sips)) &
      call elsi_allocate(elsi_h,elsi_h%ovlp_real_sips,elsi_h%nnz_l_sips,"ovlp_real_sips",caller)

   if(.not. allocated(elsi_h%row_ind_sips)) &
      call elsi_allocate(elsi_h,elsi_h%row_ind_sips,elsi_h%nnz_l_sips,"row_ind_sips",caller)

   if(.not. allocated(elsi_h%col_ptr_sips)) &
      call elsi_allocate(elsi_h,elsi_h%col_ptr_sips,(elsi_h%n_l_cols_sips+1),"col_ptr_sips",caller)

   elsi_h%ham_real_sips = h_val_recv_buffer
   call elsi_deallocate(elsi_h,h_val_recv_buffer,"h_val_recv_buffer")

   if(elsi_h%n_elsi_calls == 1) then
      elsi_h%ovlp_real_sips = s_val_recv_buffer
      call elsi_deallocate(elsi_h,s_val_recv_buffer,"s_val_recv_buffer")
   endif

   ! Compute row index and column pointer
   if(elsi_h%n_elsi_calls == 1) then
      i_col = col_recv_buffer(1)-1
      do i_val = 1,elsi_h%nnz_l_sips
         elsi_h%row_ind_sips(i_val) = row_recv_buffer(i_val)
         if(col_recv_buffer(i_val) > i_col) then
            i_col = i_col+1
            elsi_h%col_ptr_sips(i_col-col_recv_buffer(1)+1) = i_val
         endif
      enddo

      elsi_h%col_ptr_sips(elsi_h%n_l_cols_sips+1) = elsi_h%nnz_l_sips+1
   endif

   call elsi_deallocate(elsi_h,row_recv_buffer,"row_recv_buffer")
   call elsi_deallocate(elsi_h,col_recv_buffer,"col_recv_buffer")

   ! Dummy
   elsi_h%nnz_l_pexsi = elsi_h%nnz_l_sips
   elsi_h%n_l_cols_pexsi = elsi_h%n_l_cols_sips

   call elsi_stop_redistribution_time(elsi_h)

end subroutine

!>
!! This routine is a driver to convert matrix format and distribution
!! from PEXSI to BLACS.
!!
subroutine elsi_pexsi_to_blacs_hs(elsi_h,H_in,S_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                   !< Handle
   real(kind=r8),     intent(in)    :: H_in(elsi_h%nnz_l_pexsi) !< Hamiltonian matrix to be converted
   real(kind=r8),     intent(in)    :: S_in(elsi_h%nnz_l_pexsi) !< Overlap matrix to be converted

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_hs"

   if(elsi_h%n_basis < 46340) then
      call elsi_pexsi_to_blacs_hs_small(elsi_h,H_in,S_in)
   else ! use long integer
      call elsi_pexsi_to_blacs_hs_large(elsi_h,H_in,S_in)
   endif

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 1D block distributed sparse CCS format to 2D block-cyclic
!! distributed dense format, which can be used as input by ELPA.
!!
subroutine elsi_pexsi_to_blacs_hs_small(elsi_h,H_in,S_in)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h                   !< Handle
   real(kind=r8),     intent(in)    :: H_in(elsi_h%nnz_l_pexsi) !< Hamiltonian matrix to be converted
   real(kind=r8),     intent(in)    :: S_in(elsi_h%nnz_l_pexsi) !< Overlap matrix to be converted

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: k_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: global_col_id ! Global column id
   integer(kind=i4) :: global_row_id ! Global row id
   integer(kind=i4) :: local_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: local_row_id ! Local row id in 1D block distribution
   integer(kind=i4) :: proc_col_id ! Column id in process grid
   integer(kind=i4) :: proc_row_id ! Row id in process grid
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element
   integer(kind=i4), allocatable :: global_id(:) ! Global 1d id

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send_buffer(:) ! Send buffer for H
   real(kind=r8),    allocatable :: s_val_send_buffer(:) ! Send buffer for S
   integer(kind=i4), allocatable :: pos_send_buffer(:) ! Send buffer for global 1D id
   integer(kind=i4), allocatable :: send_count(:) ! Number of elements to send to each processor
   integer(kind=i4), allocatable :: send_displ(:) ! Displacement from which to take the outgoing data
   real(kind=r8),    allocatable :: h_val_recv_buffer(:) ! Receive buffer for H
   real(kind=r8),    allocatable :: s_val_recv_buffer(:) ! Receive buffer for S
   integer(kind=i4), allocatable :: pos_recv_buffer(:) ! Receive buffer for global 1D id
   integer(kind=i4), allocatable :: recv_count(:) ! Number of elements to receive from each processor
   integer(kind=i4), allocatable :: recv_displ(:) ! Displacement at which to place the incoming data
   integer(kind=i4) :: send_displ_aux ! Auxiliary variable used to set displacement
   integer(kind=i4) :: recv_displ_aux ! Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_hs_small"

   call elsi_start_redistribution_time(elsi_h)

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,s_val_send_buffer,elsi_h%nnz_l_pexsi,"s_val_send_buffer",caller)
   endif

   call elsi_allocate(elsi_h,h_val_send_buffer,elsi_h%nnz_l_pexsi,"h_val_send_buffer",caller)
   call elsi_allocate(elsi_h,pos_send_buffer,elsi_h%nnz_l_pexsi,"pos_send_buffer",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   call elsi_allocate(elsi_h,global_id,elsi_h%nnz_l_pexsi,"global_id",caller)
   call elsi_allocate(elsi_h,dest,elsi_h%nnz_l_pexsi,"dest",caller)

   i_col = 0
   ! Compute destination and global 1D id
   do i_val = 1,elsi_h%nnz_l_pexsi
      if(i_val == elsi_h%col_ptr_ccs(i_col+1) .and. i_col /= elsi_h%n_l_cols_pexsi) then
         i_col = i_col+1
      endif
      i_row = elsi_h%row_ind_ccs(i_val)

      ! Compute global id
      global_row_id = i_row
      global_col_id = i_col+elsi_h%myid*(elsi_h%n_basis/elsi_h%n_procs)
      global_id(i_val) = (global_col_id-1)*elsi_h%n_basis+global_row_id

      ! Compute destination
      proc_row_id = mod((global_row_id-1)/elsi_h%n_b_rows,elsi_h%n_p_rows)
      proc_col_id = mod((global_col_id-1)/elsi_h%n_b_cols,elsi_h%n_p_cols)
      dest(i_val) = proc_col_id+proc_row_id*elsi_h%n_p_cols
   enddo

   j_val = 0
   k_val = elsi_h%nnz_l_pexsi+1

   ! Set send_count
   if(elsi_h%n_elsi_calls == 1) then
      do i_proc = 1,elsi_h%n_procs
         do i_val = 1,elsi_h%nnz_l_pexsi
            if(dest(i_val) == i_proc-1) then
               j_val = j_val+1
               h_val_send_buffer(j_val) = H_in(i_val)
               s_val_send_buffer(j_val) = S_in(i_val)
               pos_send_buffer(j_val) = global_id(i_val)
               send_count(i_proc) = send_count(i_proc)+1
            endif
         enddo
      enddo
   else
      do i_proc = 1,elsi_h%n_procs
         do i_val = 1,elsi_h%nnz_l_pexsi
            if(dest(i_val) == i_proc-1) then
               j_val = j_val+1
               h_val_send_buffer(j_val) = H_in(i_val)
               pos_send_buffer(j_val) = global_id(i_val)
               send_count(i_proc) = send_count(i_proc)+1
            endif
         enddo
      enddo
   endif

   call elsi_deallocate(elsi_h,global_id,"global_id")
   call elsi_deallocate(elsi_h,dest,"dest")

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

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
   ! Position
   call elsi_allocate(elsi_h,pos_recv_buffer,elsi_h%nnz_l,"pos_recv_buffer",caller)

   call MPI_Alltoallv(pos_send_buffer,send_count,send_displ,mpi_integer,&
           pos_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,pos_send_buffer,"pos_send_buffer")

   ! Hamiltonian Value
   call elsi_allocate(elsi_h,h_val_recv_buffer,elsi_h%nnz_l,"h_val_recv_buffer",caller)

   call MPI_Alltoallv(h_val_send_buffer,send_count,send_displ,mpi_real8,&
           h_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_send_buffer,"h_val_send_buffer")

   ! Overlap value
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,s_val_recv_buffer,elsi_h%nnz_l,"s_val_recv_buffer",caller)

      call MPI_Alltoallv(s_val_send_buffer,send_count,send_displ,mpi_real8,&
              s_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)
   endif

   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   ! Allocate ELPA matrices
   ! Only the Hamiltonian needs to be reset everytime
   if(.not. allocated(elsi_h%ham_real_elpa)) &
      call elsi_allocate(elsi_h,elsi_h%ham_real_elpa,elsi_h%n_l_rows,elsi_h%n_l_cols,&
              "ham_real_elpa",caller)
   elsi_h%ham_real_elpa = 0.0_r8

   if(.not. allocated(elsi_h%ovlp_real_elpa)) &
      call elsi_allocate(elsi_h,elsi_h%ovlp_real_elpa,elsi_h%n_l_rows,elsi_h%n_l_cols,&
              "ovlp_real_elpa",caller)

   ! Unpack matrix
   if(elsi_h%n_elsi_calls == 1) then
      do i_val = 1,elsi_h%nnz_l
         ! Compute global 2d id
         global_col_id = (pos_recv_buffer(i_val)-1)/elsi_h%n_basis+1
         global_row_id = mod(pos_recv_buffer(i_val)-1,elsi_h%n_basis)+1

         ! Compute local 2d id
         local_row_id = (global_row_id-1)/(elsi_h%n_p_rows*elsi_h%n_b_rows)*elsi_h%n_b_rows&
                           +mod((global_row_id-1),elsi_h%n_b_rows)+1
         local_col_id = (global_col_id-1)/(elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols&
                           +mod((global_col_id-1),elsi_h%n_b_cols)+1

         ! Put value to correct position
         elsi_h%ham_real_elpa(local_row_id,local_col_id) = h_val_recv_buffer(i_val)
         elsi_h%ovlp_real_elpa(local_row_id,local_col_id) = s_val_recv_buffer(i_val)
      enddo

      call elsi_deallocate(elsi_h,s_val_recv_buffer,"s_val_recv_buffer")
   else
      do i_val = 1,elsi_h%nnz_l
         ! Compute global 2d id
         global_col_id = (pos_recv_buffer(i_val)-1)/elsi_h%n_basis+1
         global_row_id = mod(pos_recv_buffer(i_val)-1,elsi_h%n_basis)+1

         ! Compute local 2d id
         local_row_id = (global_row_id-1)/(elsi_h%n_p_rows*elsi_h%n_b_rows)*elsi_h%n_b_rows&
                           +mod((global_row_id-1),elsi_h%n_b_rows)+1
         local_col_id = (global_col_id-1)/(elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols&
                           +mod((global_col_id-1),elsi_h%n_b_cols)+1

         ! Put value to correct position
         elsi_h%ham_real_elpa(local_row_id,local_col_id) = h_val_recv_buffer(i_val)
      enddo
   endif

   call elsi_deallocate(elsi_h,h_val_recv_buffer,"h_val_recv_buffer")
   call elsi_deallocate(elsi_h,pos_recv_buffer,"pos_recv_buffer")

   call elsi_stop_redistribution_time(elsi_h)

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 1D block distributed sparse CCS format to 2D block-cyclic
!! distributed dense format, which can be used as input by ELPA.
!!
!! * 2D index is used to deal with large matrices.
!!
subroutine elsi_pexsi_to_blacs_hs_large(elsi_h,H_in,S_in)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h                   !< Handle
   real(kind=r8),     intent(in)    :: H_in(elsi_h%nnz_l_pexsi) !< Hamiltonian matrix to be converted
   real(kind=r8),     intent(in)    :: S_in(elsi_h%nnz_l_pexsi) !< Overlap matrix to be converted

   integer(kind=i4) :: mpierr
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
   integer(kind=i4), allocatable :: global_row_id(:) ! Global row id
   integer(kind=i4), allocatable :: global_col_id(:) ! Global column id
   integer(kind=i4), allocatable :: dest(:) ! Destination of each element

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: h_val_send_buffer(:) ! Send buffer for H
   real(kind=r8),    allocatable :: s_val_send_buffer(:) ! Send buffer for S
   integer(kind=i4), allocatable :: row_send_buffer(:) ! Send buffer for global row id
   integer(kind=i4), allocatable :: col_send_buffer(:) ! Send buffer for global column id
   integer(kind=i4), allocatable :: send_count(:) ! Number of elements to send to each processor
   integer(kind=i4), allocatable :: send_displ(:) ! Displacement from which to take the outgoing data
   real(kind=r8),    allocatable :: h_val_recv_buffer(:) ! Receive buffer for H
   real(kind=r8),    allocatable :: s_val_recv_buffer(:) ! Receive buffer for S
   integer(kind=i4), allocatable :: row_recv_buffer(:) ! Receive buffer for global row id
   integer(kind=i4), allocatable :: col_recv_buffer(:) ! Receive buffer for global column id
   integer(kind=i4), allocatable :: recv_count(:) ! Number of elements to receive from each processor
   integer(kind=i4), allocatable :: recv_displ(:) ! Displacement at which to place the incoming data
   integer(kind=i4) :: send_displ_aux ! Auxiliary variable used to set displacement
   integer(kind=i4) :: recv_displ_aux ! Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_hs_large"

   call elsi_start_redistribution_time(elsi_h)

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,s_val_send_buffer,elsi_h%nnz_l_pexsi,"s_val_send_buffer",caller)
   endif

   call elsi_allocate(elsi_h,h_val_send_buffer,elsi_h%nnz_l_pexsi,"h_val_send_buffer",caller)
   call elsi_allocate(elsi_h,row_send_buffer,elsi_h%nnz_l_pexsi,"pos_send_buffer",caller)
   call elsi_allocate(elsi_h,col_send_buffer,elsi_h%nnz_l_pexsi,"pos_send_buffer",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   call elsi_allocate(elsi_h,global_row_id,elsi_h%nnz_l_pexsi,"global_row_id",caller)
   call elsi_allocate(elsi_h,global_col_id,elsi_h%nnz_l_pexsi,"global_col_id",caller)
   call elsi_allocate(elsi_h,dest,elsi_h%nnz_l_pexsi,"dest",caller)

   i_col = 0
   ! Compute destination and global id
   do i_val = 1,elsi_h%nnz_l_pexsi
      if(i_val == elsi_h%col_ptr_ccs(i_col+1) .and. i_col /= elsi_h%n_l_cols_pexsi) then
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
   k_val = elsi_h%nnz_l_pexsi+1

   ! Set send_count
   if(elsi_h%n_elsi_calls == 1) then
      do i_proc = 1,elsi_h%n_procs
         do i_val = 1,elsi_h%nnz_l_pexsi
            if(dest(i_val) == i_proc-1) then
               j_val = j_val+1
               h_val_send_buffer(j_val) = H_in(i_val)
               s_val_send_buffer(j_val) = S_in(i_val)
               row_send_buffer(j_val) = global_row_id(i_val)
               col_send_buffer(j_val) = global_col_id(i_val)
               send_count(i_proc) = send_count(i_proc)+1
            endif
         enddo
      enddo
   else
      do i_proc = 1,elsi_h%n_procs
         do i_val = 1,elsi_h%nnz_l_pexsi
            if(dest(i_val) == i_proc-1) then
               j_val = j_val+1
               h_val_send_buffer(j_val) = H_in(i_val)
               row_send_buffer(j_val) = global_row_id(i_val)
               col_send_buffer(j_val) = global_col_id(i_val)
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
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

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
   ! Row index
   call elsi_allocate(elsi_h,row_recv_buffer,elsi_h%nnz_l,"row_recv_buffer",caller)

   call MPI_Alltoallv(row_send_buffer,send_count,send_displ,mpi_integer,&
           row_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buffer,"row_send_buffer")

   ! Column index
   call elsi_allocate(elsi_h,col_recv_buffer,elsi_h%nnz_l,"col_recv_buffer",caller)

   call MPI_Alltoallv(col_send_buffer,send_count,send_displ,mpi_integer,&
           col_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buffer,"col_send_buffer")

   ! Hamiltonian Value
   call elsi_allocate(elsi_h,h_val_recv_buffer,elsi_h%nnz_l,"h_val_recv_buffer",caller)

   call MPI_Alltoallv(h_val_send_buffer,send_count,send_displ,mpi_real8,&
           h_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,h_val_send_buffer,"h_val_send_buffer")

   ! Overlap value
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,s_val_recv_buffer,elsi_h%nnz_l,"s_val_recv_buffer",caller)

      call MPI_Alltoallv(s_val_send_buffer,send_count,send_displ,mpi_real8,&
              s_val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)
   endif

   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   ! Allocate ELPA matrices
   ! Only the Hamiltonian needs to be reset everytime
   if(.not. allocated(elsi_h%ham_real_elpa)) &
      call elsi_allocate(elsi_h,elsi_h%ham_real_elpa,elsi_h%n_l_rows,elsi_h%n_l_cols,&
              "ham_real_elpa",caller)
   elsi_h%ham_real_elpa = 0.0_r8

   if(.not. allocated(elsi_h%ovlp_real_elpa)) &
      call elsi_allocate(elsi_h,elsi_h%ovlp_real_elpa,elsi_h%n_l_rows,elsi_h%n_l_cols,&
              "ovlp_real_elpa",caller)

   ! Unpack matrix
   if(elsi_h%n_elsi_calls == 1) then
      do i_val = 1,elsi_h%nnz_l
         ! Compute local 2d id
         local_row_id = (row_recv_buffer(i_val)-1)/(elsi_h%n_p_rows*elsi_h%n_b_rows)*elsi_h%n_b_rows&
                           +mod((row_recv_buffer(i_val)-1),elsi_h%n_b_rows)+1
         local_col_id = (col_recv_buffer(i_val)-1)/(elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols&
                           +mod((col_recv_buffer(i_val)-1),elsi_h%n_b_cols)+1

         ! Put value to correct position
         elsi_h%ham_real_elpa(local_row_id,local_col_id) = h_val_recv_buffer(i_val)
         elsi_h%ovlp_real_elpa(local_row_id,local_col_id) = s_val_recv_buffer(i_val)
      enddo

      call elsi_deallocate(elsi_h,s_val_recv_buffer,"s_val_recv_buffer")
   else
      do i_val = 1,elsi_h%nnz_l
         ! Compute local 2d id
         local_row_id = (row_recv_buffer(i_val)-1)/(elsi_h%n_p_rows*elsi_h%n_b_rows)*elsi_h%n_b_rows&
                           +mod((row_recv_buffer(i_val)-1),elsi_h%n_b_rows)+1
         local_col_id = (col_recv_buffer(i_val)-1)/(elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols&
                           +mod((col_recv_buffer(i_val)-1),elsi_h%n_b_cols)+1

         ! Put value to correct position
         elsi_h%ham_real_elpa(local_row_id,local_col_id) = h_val_recv_buffer(i_val)
      enddo
   endif

   call elsi_deallocate(elsi_h,h_val_recv_buffer,"h_val_recv_buffer")
   call elsi_deallocate(elsi_h,row_recv_buffer,"row_recv_buffer")
   call elsi_deallocate(elsi_h,col_recv_buffer,"col_recv_buffer")

   call elsi_stop_redistribution_time(elsi_h)

end subroutine

!>
!! This routine is a driver to convert matrix format and distribution
!! from BLACS to PEXSI.
!!
subroutine elsi_blacs_to_pexsi_dm(elsi_h,D_out)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h                    !< Handle
   real(kind=r8),     intent(out)   :: D_out(elsi_h%nnz_l_pexsi) !< Density matrix to be converted

   character*40, parameter :: caller = "elsi_blacs_to_pexsi_dm"

   if(elsi_h%n_basis < 46340) then
      call elsi_blacs_to_pexsi_dm_small(elsi_h,D_out)
   else ! use long integer
      call elsi_blacs_to_pexsi_dm_large(elsi_h,D_out)
   endif

end subroutine

!>
!! This routine converts density matrix in 2D block-cyclic distributed
!! dense format to 1D block distributed sparse CCS format.
!!
subroutine elsi_blacs_to_pexsi_dm_small(elsi_h,D_out)

   implicit none

   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h                    !< Handle
   real(kind=r8),     intent(out)   :: D_out(elsi_h%nnz_l_pexsi) !< Density matrix to be converted

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: i_val
   integer(kind=i4) :: j_val
   integer(kind=i4) :: i_proc
   integer(kind=i4) :: global_col_id ! Global column id
   integer(kind=i4) :: global_row_id ! Global row id
   integer(kind=i4) :: local_col_id ! Local column id in 1D block distribution
   integer(kind=i4) :: local_row_id ! Local row id in 1D block distribution
   integer(kind=i4) :: dest ! Destination of an element
   integer(kind=i4) :: nnz_l_pexsi_aux
   integer(kind=i4) :: n_l_cols_pexsi_aux

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: val_send_buffer(:) ! Send buffer for value
   integer(kind=i4), allocatable :: pos_send_buffer(:) ! Send buffer for global 1D id
   integer(kind=i4), allocatable :: send_count(:) ! Number of elements to send to each processor
   integer(kind=i4), allocatable :: send_displ(:) ! Displacement from which to take the outgoing data
   real(kind=r8),    allocatable :: val_recv_buffer(:) ! Receive buffer for value
   integer(kind=i4), allocatable :: pos_recv_buffer(:) ! Receive buffer for global 1D id
   integer(kind=i4), allocatable :: recv_count(:) ! Number of elements to receive from each processor
   integer(kind=i4), allocatable :: recv_displ(:) ! Displacement at which to place the incoming data
   integer(kind=i4) :: send_displ_aux ! Auxiliary variable used to set displacement 
   integer(kind=i4) :: recv_displ_aux ! Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_blacs_to_pexsi_dm_small"

   call elsi_start_redistribution_time(elsi_h)

   if(elsi_h%solver == ELPA) then
      call elsi_get_local_nnz(elsi_h,elsi_h%den_mat_elpa,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,elsi_h%nnz_l)
   elseif(elsi_h%solver == LIBOMM) then
      call elsi_get_local_nnz(elsi_h,elsi_h%den_mat_omm%dval,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,elsi_h%nnz_l)
   endif

   call elsi_allocate(elsi_h,val_send_buffer,elsi_h%nnz_l,"val_send_buffer",caller)
   call elsi_allocate(elsi_h,pos_send_buffer,elsi_h%nnz_l,"pos_send_buffer",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   ! Compute destination and global 1D id
   if(elsi_h%solver == ELPA) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(elsi_h%den_mat_elpa(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest = (global_col_id-1)/(elsi_h%n_basis/elsi_h%n_procs)
               ! The last process may take more
               dest = min(dest,elsi_h%n_procs-1)
               ! Compute the global id
               ! Pack global id and data into buffers
               pos_send_buffer(i_val) = (global_col_id-1)*elsi_h%n_basis+global_row_id
               val_send_buffer(i_val) = elsi_h%den_mat_elpa(i_row,i_col)
               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   elseif(elsi_h%solver == LIBOMM) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(elsi_h%den_mat_omm%dval(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest = (global_col_id-1)/(elsi_h%n_basis/elsi_h%n_procs)
               ! The last process may take more
               dest = min(dest,elsi_h%n_procs-1)
               ! Compute the global id
               ! Pack global id and data into buffers
               pos_send_buffer(i_val) = (global_col_id-1)*elsi_h%n_basis+global_row_id
               val_send_buffer(i_val) = elsi_h%den_mat_omm%dval(i_row,i_col)
               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   endif

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

   ! Set local number of nonzero
   nnz_l_pexsi_aux = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   send_displ_aux = 0
   recv_displ_aux = 0

   do i_proc = 0,elsi_h%n_procs-1
      send_displ(i_proc+1) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc+1)
      recv_displ(i_proc+1) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
   enddo

   ! Send and receive the packed data
   ! Position
   call elsi_allocate(elsi_h,pos_recv_buffer,nnz_l_pexsi_aux,"pos_recv_buffer",caller)

   call MPI_Alltoallv(pos_send_buffer,send_count,send_displ,mpi_integer,&
           pos_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,pos_send_buffer,"pos_send_buffer")

   ! Density matrix value
   call elsi_allocate(elsi_h,val_recv_buffer,nnz_l_pexsi_aux,"val_recv_buffer",caller)

   call MPI_Alltoallv(val_send_buffer,send_count,send_displ,mpi_real8,&
           val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,val_send_buffer,"val_send_buffer")

   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   D_out = 0.0_r8

   if(elsi_h%myid == elsi_h%n_procs-1) then
      n_l_cols_pexsi_aux = elsi_h%n_basis-elsi_h%n_l_cols_pexsi
      n_l_cols_pexsi_aux = n_l_cols_pexsi_aux/(elsi_h%n_procs-1)
   else
      n_l_cols_pexsi_aux = elsi_h%n_l_cols_pexsi
   endif

   ! Unpack matrix
   do i_val = 1,nnz_l_pexsi_aux
      global_col_id = (pos_recv_buffer(i_val)-1)/elsi_h%n_basis+1
      local_col_id = global_col_id-elsi_h%myid*n_l_cols_pexsi_aux
      local_row_id = mod(pos_recv_buffer(i_val)-1,elsi_h%n_basis)+1

      do j_val = elsi_h%col_ptr_ccs(local_col_id),elsi_h%col_ptr_ccs(local_col_id+1)-1
         if(elsi_h%row_ind_ccs(j_val) == local_row_id) then
            D_out(j_val) = val_recv_buffer(i_val)
!            exit
         endif
      enddo
   enddo

   call elsi_deallocate(elsi_h,pos_recv_buffer,"pos_recv_buffer")
   call elsi_deallocate(elsi_h,val_recv_buffer,"val_recv_buffer")

   call elsi_stop_redistribution_time(elsi_h)

end subroutine

!>
!! This routine converts density matrix in 2D block-cyclic distributed
!! dense format to 1D block distributed sparse CCS format.
!!
!! * 2D index is used to deal with large matrices.
!!
subroutine elsi_blacs_to_pexsi_dm_large(elsi_h,D_out)

   implicit none

   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h                    !< Handle
   real(kind=r8),     intent(out)   :: D_out(elsi_h%nnz_l_pexsi) !< Density matrix to be converted

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

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable :: val_send_buffer(:) ! Send buffer for value
   integer(kind=i4), allocatable :: row_send_buffer(:) ! Send buffer for global row id
   integer(kind=i4), allocatable :: col_send_buffer(:) ! Send buffer for global column id
   integer(kind=i4), allocatable :: send_count(:) ! Number of elements to send to each processor
   integer(kind=i4), allocatable :: send_displ(:) ! Displacement from which to take the outgoing data
   real(kind=r8),    allocatable :: val_recv_buffer(:) ! Receive buffer for value
   integer(kind=i4), allocatable :: row_recv_buffer(:) ! Receive buffer for global row id
   integer(kind=i4), allocatable :: col_recv_buffer(:) ! Receive buffer for global column id
   integer(kind=i4), allocatable :: recv_count(:) ! Number of elements to receive from each processor
   integer(kind=i4), allocatable :: recv_displ(:) ! Displacement at which to place the incoming data
   integer(kind=i4) :: send_displ_aux ! Auxiliary variable used to set displacement 
   integer(kind=i4) :: recv_displ_aux ! Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_blacs_to_pexsi_dm_large"

   call elsi_start_redistribution_time(elsi_h)

   if(elsi_h%solver == ELPA) then
      call elsi_get_local_nnz(elsi_h,elsi_h%den_mat_elpa,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,elsi_h%nnz_l)
   elseif(elsi_h%solver == LIBOMM) then
      call elsi_get_local_nnz(elsi_h,elsi_h%den_mat_omm%dval,elsi_h%n_l_rows,&
              elsi_h%n_l_cols,elsi_h%nnz_l)
   endif

   call elsi_allocate(elsi_h,val_send_buffer,elsi_h%nnz_l,"val_send_buffer",caller)
   call elsi_allocate(elsi_h,row_send_buffer,elsi_h%nnz_l,"row_send_buffer",caller)
   call elsi_allocate(elsi_h,col_send_buffer,elsi_h%nnz_l,"col_send_buffer",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   ! Compute destination and global id
   if(elsi_h%solver == ELPA) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(elsi_h%den_mat_elpa(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               ! Compute global id
               call elsi_get_global_col(elsi_h,col_send_buffer(i_val),i_col)
               call elsi_get_global_row(elsi_h,row_send_buffer(i_val),i_row)
               ! Compute destination
               dest = (col_send_buffer(i_val)-1)/(elsi_h%n_basis/elsi_h%n_procs)
               ! The last process may take more
               dest = min(dest,elsi_h%n_procs-1)
               ! Pack data
               val_send_buffer(i_val) = elsi_h%den_mat_elpa(i_row,i_col)
               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   elseif(elsi_h%solver == LIBOMM) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(elsi_h%den_mat_omm%dval(i_row,i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               ! Compute global id
               call elsi_get_global_col(elsi_h,col_send_buffer(i_val),i_col)
               call elsi_get_global_row(elsi_h,row_send_buffer(i_val),i_row)
               ! Compute destination
               dest = (col_send_buffer(i_val)-1)/(elsi_h%n_basis/elsi_h%n_procs)
               ! The last process may take more
               dest = min(dest,elsi_h%n_procs-1)
               ! Pack data
               val_send_buffer(i_val) = elsi_h%den_mat_omm%dval(i_row,i_col)
               ! Set send_count
               send_count(dest+1) = send_count(dest+1)+1
            endif
         enddo
      enddo
   endif

   call elsi_allocate(elsi_h,recv_count,elsi_h%n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
           1,mpi_integer,elsi_h%mpi_comm,mpierr)

   ! Set local number of nonzero
   nnz_l_pexsi_aux = sum(recv_count,1)

   ! Set send and receive displacement
   call elsi_allocate(elsi_h,send_displ,elsi_h%n_procs,"send_displ",caller)
   call elsi_allocate(elsi_h,recv_displ,elsi_h%n_procs,"recv_displ",caller)

   send_displ_aux = 0
   recv_displ_aux = 0

   do i_proc = 0,elsi_h%n_procs-1
      send_displ(i_proc+1) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc+1)
      recv_displ(i_proc+1) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(elsi_h,row_recv_buffer,nnz_l_pexsi_aux,"row_recv_buffer",caller)

   call MPI_Alltoallv(row_send_buffer,send_count,send_displ,mpi_integer,&
           row_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,row_send_buffer,"row_send_buffer")

   ! Column id
   call elsi_allocate(elsi_h,col_recv_buffer,nnz_l_pexsi_aux,"col_recv_buffer",caller)

   call MPI_Alltoallv(col_send_buffer,send_count,send_displ,mpi_integer,&
           col_recv_buffer,recv_count,recv_displ,mpi_integer,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,col_send_buffer,"col_send_buffer")

   ! Density matrix value
   call elsi_allocate(elsi_h,val_recv_buffer,nnz_l_pexsi_aux,"val_recv_buffer",caller)

   call MPI_Alltoallv(val_send_buffer,send_count,send_displ,mpi_real8,&
           val_recv_buffer,recv_count,recv_displ,mpi_real8,elsi_h%mpi_comm,mpierr)

   call elsi_deallocate(elsi_h,val_send_buffer,"val_send_buffer")

   call elsi_deallocate(elsi_h,send_count,"send_count")
   call elsi_deallocate(elsi_h,recv_count,"recv_count")
   call elsi_deallocate(elsi_h,send_displ,"send_displ")
   call elsi_deallocate(elsi_h,recv_displ,"recv_displ")

   D_out = 0.0_r8

   if(elsi_h%myid == elsi_h%n_procs-1) then
      n_l_cols_pexsi_aux = elsi_h%n_basis-elsi_h%n_l_cols_pexsi
      n_l_cols_pexsi_aux = n_l_cols_pexsi_aux/(elsi_h%n_procs-1)
   else
      n_l_cols_pexsi_aux = elsi_h%n_l_cols_pexsi
   endif

   ! Unpack matrix
   do i_val = 1,nnz_l_pexsi_aux
      local_col_id = col_recv_buffer(i_val)-elsi_h%myid*n_l_cols_pexsi_aux
      local_row_id = row_recv_buffer(i_val)

      do j_val = elsi_h%col_ptr_ccs(local_col_id),elsi_h%col_ptr_ccs(local_col_id+1)-1
         if(elsi_h%row_ind_ccs(j_val) == local_row_id) then
            D_out(j_val) = val_recv_buffer(i_val)
!            exit
         endif
      enddo
   enddo

   call elsi_deallocate(elsi_h,row_recv_buffer,"row_recv_buffer")
   call elsi_deallocate(elsi_h,col_recv_buffer,"col_recv_buffer")
   call elsi_deallocate(elsi_h,val_recv_buffer,"val_recv_buffer")

   call elsi_stop_redistribution_time(elsi_h)

end subroutine

end module ELSI_MATCONV
