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
!! This module provides interfaces to PEXSI.
!!

module ELSI_PEXSI

   use iso_c_binding
   use ELSI_DIMENSIONS
   use ELSI_TIMERS
   use ELSI_UTILS
   use f_ppexsi_interface

   implicit none
   private

   public :: elsi_init_pexsi
   public :: elsi_blacs_to_pexsi
   public :: elsi_pexsi_to_blacs_dm
   public :: elsi_solve_evp_pexsi
   public :: elsi_set_pexsi_default_options
   public :: elsi_print_pexsi_options

contains

!=========================
! ELSI routines for PEXSI
!=========================

!>
!! PEXSI processor grid setup.
!!
subroutine elsi_init_pexsi()

   implicit none

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_init_pexsi"

   if(n_elsi_calls == 1) then
      if(.not.n_p_per_pole_ready) then
         if(mod(n_procs,pexsi_options%numPole) == 0) then
            n_p_per_pole_pexsi = n_procs/pexsi_options%numPole

            call elsi_statement_print("  PEXSI parallel over poles.")
            write(info_str,"(A,I13)") "  | Number of MPI tasks per pole: ",&
               n_p_per_pole_pexsi
            call elsi_statement_print(info_str)
         else
            call elsi_stop("  PEXSI not parallel over poles. High"//&
                           " performance of PEXSI is expected if the"//&
                           " number of MPI tasks is a multiple of"//&
                           " the number of PEXSI poles. Please adjust"//&
                           " either the number of MPI tasks, or the"//&
                           " number of poles. Exiting...",caller)
         endif
      endif

      ! Set square-like process grid for selected inversion of each pole
      do n_p_rows_pexsi = nint(sqrt(real(n_p_per_pole_pexsi))),2,-1
         if(mod(n_p_per_pole_pexsi,n_p_rows_pexsi) == 0) exit
      enddo

      n_p_cols_pexsi = n_p_per_pole_pexsi/n_p_rows_pexsi

      ! PEXSI process grid
      my_p_col_pexsi = mod(myid,n_p_per_pole_pexsi)
      my_p_row_pexsi = myid/n_p_per_pole_pexsi

      ! PEXSI uses a pure block distribution in the first process row
      n_b_rows_pexsi = n_g_size

      ! The last process holds all remaining columns
      n_b_cols_pexsi = n_g_size/n_p_per_pole_pexsi
      if(my_p_col_pexsi == n_p_per_pole_pexsi-1) then
         n_b_cols_pexsi = n_g_size-(n_p_per_pole_pexsi-1)*n_b_cols_pexsi
      endif

      n_l_rows_pexsi = n_b_rows_pexsi
      n_l_cols_pexsi = n_b_cols_pexsi

      ! Only master process outputs
      if(myid == 0) then
         pexsi_output_file_index = 0
      else
         pexsi_output_file_index = -1
      endif

      pexsi_plan = f_ppexsi_plan_initialize(mpi_comm_global,n_p_rows_pexsi,&
                      n_p_cols_pexsi,pexsi_output_file_index,pexsi_info)

      if(pexsi_info /= 0) &
         call elsi_stop(" PEXSI plan initialization failed. Exiting...",caller)

   endif

end subroutine

!>
!! This routine is a driver to convert matrix format and distribution
!! from BLACS to PEXSI.
!!
subroutine elsi_blacs_to_pexsi(H_in,S_in)

   implicit none

   real*8, intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian matrix to be converted
   real*8, intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap matrix to be converted

   character*40, parameter :: caller = "elsi_blacs_to_pexsi"

   if(overlap_is_unit) then
      !TODO
      call elsi_stop(" PEXSI with identity overlap matrix not yet available."//&
                     " Exiting...",caller)
   else
      if(n_g_size < 46340) then ! kind=4 integer works
         call elsi_blacs_to_pexsi_hs_small(H_in,S_in)
      else ! use kind=8 integer
         call elsi_blacs_to_pexsi_hs_large(H_in,S_in)
      endif
   endif

end subroutine

!>
!! This routine is a driver to convert matrix format and distribution
!! from PEXSI to BLACS.
!!
subroutine elsi_pexsi_to_blacs_dm(D_out)

   implicit none

   real*8, intent(out) :: D_out(n_l_rows,n_l_cols) !< Density matrix to be converted

   if(n_g_size < 46340) then
      call elsi_pexsi_to_blacs_dm_small(D_out)
   else ! use kind=8 integer
      call elsi_pexsi_to_blacs_dm_large(D_out)
   endif

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by PEXSI.
!!
!! Usage:
!!
!! * PEXSI pole parallelism MUST be available. Data redistributed on
!!   the processors corresponding to the first pole.
!!
!! * Overlap is NOT an identity matrix. Both Hamiltonian and overlap
!!   are considered.
!!
!! * See also elsi_blacs_to_pexsi_hs_large
!!
subroutine elsi_blacs_to_pexsi_hs_small(H_in,S_in)

   implicit none
   include "mpif.h"

   real*8, intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian matrix to be converted
   real*8, intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap matrix to be converted
   integer :: i_row !< Row counter
   integer :: i_col !< Col counter
   integer :: i_val,j_val !< Value counter
   integer :: i_proc !< Process counter
   integer :: global_col_id !< Global column id
   integer :: global_row_id !< Global row id
   integer :: local_col_id !< Local column id in 1D block distribution
   integer :: local_row_id !< Local row id in 1D block distribution
   integer :: d1,d2,d11,d12,d21,d22 !< Number of columns in the intermediate stage
   integer :: this_n_cols
   integer :: tmp_int
   integer :: min_pos,min_id
   integer :: nnz_l_pexsi_aux,mpi_comm_aux_pexsi
   integer, allocatable :: dest(:) !< Destination of each element
   integer, allocatable :: locat(:) !< Location of each global column
   real*8 :: tmp_real

   ! See documentation of MPI_Alltoallv
   real*8, allocatable  :: h_val_send_buffer(:) !< Send buffer for Hamiltonian
   real*8, allocatable  :: s_val_send_buffer(:) !< Send buffer for overlap
   integer, allocatable :: pos_send_buffer(:)   !< Send buffer for global 1D id
   integer :: send_count(n_procs) !< Number of elements to send to each processor
   integer :: send_displ(n_procs) !< Displacement from which to take the outgoing data
   integer :: send_displ_aux      !< Auxiliary variable used to set displacement
   real*8, allocatable  :: h_val_recv_buffer(:) !< Receive buffer for Hamiltonian
   real*8, allocatable  :: s_val_recv_buffer(:) !< Receive buffer for overlap
   integer, allocatable :: pos_recv_buffer(:)   !< Receive buffer for global 1D id
   integer :: recv_count(n_procs) !< Number of elements to receive from each processor
   integer :: recv_displ(n_procs) !< Displacement at which to place the incoming data
   integer :: recv_displ_aux      !< Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_blacs_to_pexsi_hs_small"

   call elsi_start_redistribution_time()

   send_count = 0
   send_displ = 0
   recv_count = 0
   recv_displ = 0

   if(n_elsi_calls == 1) then
      call elsi_get_local_nnz(S_in,n_l_rows,n_l_cols,nnz_l)
      call elsi_allocate(s_val_send_buffer,nnz_l,"s_val_send_buffer",caller)
   endif

   call elsi_allocate(locat,n_g_size,"locat",caller)
   call elsi_allocate(dest,nnz_l,"dest",caller)
   call elsi_allocate(pos_send_buffer,nnz_l,"pos_send_buffer",caller)
   call elsi_allocate(h_val_send_buffer,nnz_l,"h_val_send_buffer",caller)

   ! Compute d1,d2,d11,d12,d21,d22 (need explanation)
   d1  = n_g_size/n_p_per_pole_pexsi
   d2  = n_g_size-(n_p_per_pole_pexsi-1)*d1
   d11 = d1/pexsi_options%numPole
   d12 = d1-(pexsi_options%numPole-1)*d11
   d21 = d2/pexsi_options%numPole
   d22 = d2-(pexsi_options%numPole-1)*d21

   i_val = 0

   if(d11 > 0) then
      do i_proc = 0,n_procs-pexsi_options%numPole-1
         if(mod((i_proc+1),pexsi_options%numPole) == 0) then
            this_n_cols = d12
         else
            this_n_cols = d11
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = 0,n_procs-pexsi_options%numPole-1
         if(1+mod(i_proc,pexsi_options%numPole) .le. d1) then
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
      do i_proc = n_procs-pexsi_options%numPole,n_procs-1
         if(mod((i_proc+1),pexsi_options%numPole) == 0) then
            this_n_cols = d22
         else
            this_n_cols = d21
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = n_procs-pexsi_options%numPole,n_procs-1
         if(1+mod(i_proc,pexsi_options%numPole) .le. d2) then
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

   ! Compute destination and global 1D id
   if(n_elsi_calls == 1) then
      i_val = 0
      do i_col = 1,n_l_cols
         do i_row = 1,n_l_rows
            if(abs(S_in(i_row,i_col)) > zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(global_col_id,i_col)
               call elsi_get_global_row(global_row_id,i_row)
               ! Compute destination
               dest(i_val) = locat(global_col_id)
               ! The last process may take more
               if(dest(i_val) > (n_procs-1)) dest(i_val) = n_procs-1
               ! Compute the global id
               ! Pack global id and data into buffers
               pos_send_buffer(i_val) = (global_col_id-1)*n_g_size+global_row_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               s_val_send_buffer(i_val) = S_in(i_row,i_col)
               ! Set send_count
               send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
            endif
         enddo
      enddo
   else
      i_val = 0
      do i_col = 1,n_l_cols
         do i_row = 1,n_l_rows
            if(abs(S_in(i_row,i_col)) > zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(global_col_id,i_col)
               call elsi_get_global_row(global_row_id,i_row)
               ! Compute destination
               dest(i_val) = locat(global_col_id)
               ! The last process may take more
               if(dest(i_val) > (n_procs-1)) dest(i_val) = n_procs-1
               ! Compute the global id
               ! Pack global id and data into buffers
               pos_send_buffer(i_val) = (global_col_id-1)*n_g_size+global_row_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               ! Set send_count
               send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
            endif
         enddo
      enddo
   endif

   deallocate(dest)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
                     1,mpi_integer,mpi_comm_global,mpierr)

   ! Set local/global number of nonzero
   nnz_l_pexsi_aux = sum(recv_count,1)
   call MPI_Allreduce(nnz_l_pexsi_aux,nnz_g,1,mpi_integer,mpi_sum,&
                      mpi_comm_global,mpierr)

   ! Set send and receive displacement
   send_displ_aux = 0
   recv_displ_aux = 0
   do i_proc = 0,n_procs-1
      send_displ(i_proc+1) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc+1)
      recv_displ(i_proc+1) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
   enddo

   ! Send and receive the packed data
   ! Position
   call elsi_allocate(pos_recv_buffer,nnz_l_pexsi_aux,"pos_recv_buffer",caller)

   call MPI_Alltoallv(pos_send_buffer,send_count,send_displ,mpi_integer,&
                      pos_recv_buffer,recv_count,recv_displ,mpi_integer,&
                      mpi_comm_global,mpierr)

   deallocate(pos_send_buffer)

   ! Hamiltonian value
   call elsi_allocate(h_val_recv_buffer,nnz_l_pexsi_aux,"h_val_recv_buffer",caller)

   call MPI_Alltoallv(h_val_send_buffer,send_count,send_displ,mpi_real8,&
                      h_val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                      mpi_comm_global,mpierr)

   deallocate(h_val_send_buffer)

   ! Overlap value
   if(n_elsi_calls == 1) then
      call elsi_allocate(s_val_recv_buffer,nnz_l_pexsi_aux,"s_val_recv_buffer",caller)

      call MPI_Alltoallv(s_val_send_buffer,send_count,send_displ,mpi_real8,&
                         s_val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                         mpi_comm_global,mpierr)

      deallocate(s_val_send_buffer)
   endif

   ! Unpack and reorder
   if(n_elsi_calls == 1) then
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
   send_count(myid/pexsi_options%numPole+1) = nnz_l_pexsi_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
                     1,mpi_integer,mpi_comm_global,mpierr)

   if(n_elsi_calls == 1) then
      nnz_l_pexsi = sum(recv_count,1)

      ! At this point only the first pole knows nnz_l_pexsi
      call MPI_Comm_split(mpi_comm_global,my_p_col_pexsi,my_p_row_pexsi,&
                          mpi_comm_aux_pexsi,mpierr)

      call MPI_Bcast(nnz_l_pexsi,1,mpi_integer,0,mpi_comm_aux_pexsi,mpierr)
   endif

   ! Set send and receive displacement
   send_displ_aux = 0
   recv_displ_aux = 0
   do i_proc = 0,n_procs-1
      send_displ(i_proc+1) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc+1)
      recv_displ(i_proc+1) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
   enddo

   ! Allocate PEXSI matrices
   ! Only the Hamiltonian needs to be reset everytime
   if(.not.allocated(ham_real_pexsi)) &
      call elsi_allocate(ham_real_pexsi,nnz_l_pexsi,"ham_real_pexsi",caller)
   ham_real_pexsi = 0.0d0

   if(.not.allocated(ovlp_real_pexsi)) &
      call elsi_allocate(ovlp_real_pexsi,nnz_l_pexsi,"ovlp_real_pexsi",caller)

   if(.not.allocated(row_ind_pexsi)) &
      call elsi_allocate(row_ind_pexsi,nnz_l_pexsi,"row_ind_pexsi",caller)

   if(.not.allocated(col_ptr_pexsi)) &
      call elsi_allocate(col_ptr_pexsi,(n_l_cols_pexsi+1),"col_ptr_pexsi",caller)

   ! Send and receive the packed data
   if(n_elsi_calls == 1) then
      ! Position
      call elsi_allocate(pos_send_buffer,nnz_l_pexsi,"pos_send_buffer",caller)

      call MPI_Alltoallv(pos_recv_buffer,send_count,send_displ,mpi_integer,&
                         pos_send_buffer,recv_count,recv_displ,mpi_integer,&
                         mpi_comm_global,mpierr)

      deallocate(pos_recv_buffer)

      ! Overlap value
      call MPI_Alltoallv(s_val_recv_buffer,send_count,send_displ,mpi_real8,&
                         ovlp_real_pexsi,recv_count,recv_displ,mpi_real8,&
                         mpi_comm_global,mpierr)

      deallocate(s_val_recv_buffer)
   endif

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv_buffer,send_count,send_displ,mpi_real8,&
                      ham_real_pexsi,recv_count,recv_displ,mpi_real8,&
                      mpi_comm_global,mpierr)

   deallocate(h_val_recv_buffer)

   ! Only the first pole computes row index and column pointer
   if(n_elsi_calls == 1) then
      if(my_p_row_pexsi == 0) then
         ! Compute row index and column pointer
         i_col = (pos_send_buffer(1)-1)/n_g_size
         do i_val = 1,nnz_l_pexsi
            row_ind_pexsi(i_val) = mod(pos_send_buffer(i_val)-1,n_g_size)+1
            if((pos_send_buffer(i_val)-1)/n_g_size+1 > i_col) then
               i_col = i_col+1
               col_ptr_pexsi(i_col-(pos_send_buffer(1)-1)/n_g_size) = i_val
            endif
         enddo

         col_ptr_pexsi(n_l_cols_pexsi+1) = nnz_l_pexsi+1
      endif

      deallocate(pos_send_buffer)
   endif

   call elsi_stop_redistribution_time()

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by PEXSI.
!!
!! Usage:
!!
!! * Integer(kind=8) is used to deal with large matrices.
!!
!! * PEXSI pole parallelism MUST be available. Data redistributed on
!!   the processors corresponding to the first pole.
!!
!! * Overlap is NOT an identity matrix. Both Hamiltonian and overlap
!!   are considered.
!!
!! * See also elsi_blacs_to_pexsi_hs_small.
!!
subroutine elsi_blacs_to_pexsi_hs_large(H_in,S_in)

   implicit none
   include "mpif.h"

   real*8, intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian matrix to be converted
   real*8, intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap matrix to be converted
   integer :: i_row !< Row counter
   integer :: i_col !< Col counter
   integer :: i_val,j_val !< Value counter
   integer :: i_proc !< Process counter
   integer :: global_col_id !< Global column id
   integer :: global_row_id !< Global row id
   integer :: local_col_id !< Local column id in 1D block distribution
   integer :: local_row_id !< Local row id in 1D block distribution
   integer :: d1,d2,d11,d12,d21,d22 !< Number of columns in the intermediate stage
   integer :: this_n_cols
   integer :: min_pos,min_id
   integer :: nnz_l_pexsi_aux,mpi_comm_aux_pexsi
   integer :: tmp_int
   integer(kind=8) :: tmp_long
   integer, allocatable :: dest(:) !< Destination of each element
   integer, allocatable :: locat(:) !< Location of each global column
   real*8 :: tmp_real

   ! See documentation of MPI_Alltoallv
   real*8, allocatable  :: h_val_send_buffer(:) !< Send buffer for Hamiltonian
   real*8, allocatable  :: s_val_send_buffer(:) !< Send buffer for overlap
   integer, allocatable :: row_send_buffer(:)   !< Send buffer for global row id
   integer, allocatable :: col_send_buffer(:)   !< Send buffer for global column id
   integer :: send_count(n_procs) !< Number of elements to send to each processor
   integer :: send_displ(n_procs) !< Displacement from which to take the outgoing data
   integer :: send_displ_aux      !< Auxiliary variable used to set displacement
   real*8, allocatable  :: h_val_recv_buffer(:) !< Receive buffer for Hamiltonian
   real*8, allocatable  :: s_val_recv_buffer(:) !< Receive buffer for overlap
   integer, allocatable :: row_recv_buffer(:)   !< Receive buffer for global row id
   integer, allocatable :: col_recv_buffer(:)   !< Receive buffer for global column id
   integer :: recv_count(n_procs) !< Number of elements to receive from each processor
   integer :: recv_displ(n_procs) !< Displacement at which to place the incoming data
   integer :: recv_displ_aux      !< Auxiliary variable used to set displacement
   integer(kind=8), allocatable :: global_id(:) !< Global 1D id

   character*40, parameter :: caller = "elsi_blacs_to_pexsi_hs_large"

   call elsi_start_redistribution_time()

   send_count = 0
   send_displ = 0
   recv_count = 0
   recv_displ = 0

   if(n_elsi_calls == 1) then
      call elsi_get_local_nnz(S_in,n_l_rows,n_l_cols,nnz_l)
      call elsi_allocate(s_val_send_buffer,nnz_l,"s_val_send_buffer",caller)
   endif

   call elsi_allocate(locat,n_g_size,"locat",caller)
   call elsi_allocate(dest,nnz_l,"dest",caller)
   call elsi_allocate(row_send_buffer,nnz_l,"row_send_buffer",caller)
   call elsi_allocate(col_send_buffer,nnz_l,"col_send_buffer",caller)
   call elsi_allocate(h_val_send_buffer,nnz_l,"h_val_send_buffer",caller)

   ! Compute d1,d2,d11,d12,d21,d22 (need explanation)
   d1  = n_g_size/n_p_per_pole_pexsi
   d2  = n_g_size-(n_p_per_pole_pexsi-1)*d1
   d11 = d1/pexsi_options%numPole
   d12 = d1-(pexsi_options%numPole-1)*d11
   d21 = d2/pexsi_options%numPole
   d22 = d2-(pexsi_options%numPole-1)*d21

   i_val = 0

   if(d11 > 0) then
      do i_proc = 0,n_procs-pexsi_options%numPole-1
         if(mod((i_proc+1),pexsi_options%numPole) == 0) then
            this_n_cols = d12
         else
            this_n_cols = d11
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = 0,n_procs-pexsi_options%numPole-1
         if(1+mod(i_proc,pexsi_options%numPole) .le. d1) then
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
      do i_proc = n_procs-pexsi_options%numPole,n_procs-1
         if(mod((i_proc+1),pexsi_options%numPole) == 0) then
            this_n_cols = d22
         else
            this_n_cols = d21
         endif

         locat(i_val+1:i_val+this_n_cols) = i_proc
         i_val = i_val+this_n_cols
      enddo
   else
      do i_proc = n_procs-pexsi_options%numPole,n_procs-1
         if(1+mod(i_proc,pexsi_options%numPole) .le. d2) then
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

   ! Compute destination and global id
   if(n_elsi_calls == 1) then
      i_val = 0
      do i_col = 1,n_l_cols
         do i_row = 1,n_l_rows
            if(abs(S_in(i_row,i_col)) > zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(global_col_id,i_col)
               call elsi_get_global_row(global_row_id,i_row)
               ! Compute destination
               dest(i_val) = locat(global_col_id)
               ! The last process may take more
               if(dest(i_val) > (n_procs-1)) dest(i_val) = n_procs-1
               ! Pack global id and data into buffers
               row_send_buffer(i_val) = global_row_id
               col_send_buffer(i_val) = global_col_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               s_val_send_buffer(i_val) = S_in(i_row,i_col)
               ! Set send_count
               send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
            endif
         enddo
      enddo
   else
      i_val = 0
      do i_col = 1,n_l_cols
         do i_row = 1,n_l_rows
            if(abs(S_in(i_row,i_col)) > zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(global_col_id,i_col)
               call elsi_get_global_row(global_row_id,i_row)
               ! Compute destination
               dest(i_val) = locat(global_col_id)
               ! The last process may take more
               if(dest(i_val) > (n_procs-1)) dest(i_val) = n_procs-1
               ! Pack global id and data into buffers
               row_send_buffer(i_val) = global_row_id
               col_send_buffer(i_val) = global_col_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
              ! Set send_count
               send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
            endif
         enddo
      enddo
   endif

   deallocate(dest)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
                     1,mpi_integer,mpi_comm_global,mpierr)

   ! Set local/global number of nonzero
   nnz_l_pexsi_aux = sum(recv_count,1)
   call MPI_Allreduce(nnz_l_pexsi_aux,nnz_g,1,mpi_integer,mpi_sum,&
                      mpi_comm_global,mpierr)

   ! Set send and receive displacement
   send_displ_aux = 0
   recv_displ_aux = 0
   do i_proc = 0,n_procs-1
      send_displ(i_proc+1) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc+1)
      recv_displ(i_proc+1) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
   enddo

   ! Send and receive the packed data
   ! Row id
   call elsi_allocate(row_recv_buffer,nnz_l_pexsi_aux,"row_recv_buffer",caller)

   call MPI_Alltoallv(row_send_buffer,send_count,send_displ,mpi_integer,&
                      row_recv_buffer,recv_count,recv_displ,mpi_integer,&
                      mpi_comm_global,mpierr)

   deallocate(row_send_buffer)

   ! Column id
   call elsi_allocate(col_recv_buffer,nnz_l_pexsi_aux,"col_recv_buffer",caller)

   call MPI_Alltoallv(col_send_buffer,send_count,send_displ,mpi_integer,&
                      col_recv_buffer,recv_count,recv_displ,mpi_integer,&
                      mpi_comm_global,mpierr)

   deallocate(col_send_buffer)

   ! Hamiltonian value
   call elsi_allocate(h_val_recv_buffer,nnz_l_pexsi_aux,"h_val_recv_buffer",caller)

   call MPI_Alltoallv(h_val_send_buffer,send_count,send_displ,mpi_real8,&
                      h_val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                      mpi_comm_global,mpierr)

   deallocate(h_val_send_buffer)

   ! Overlap value
   if(n_elsi_calls == 1) then
      call elsi_allocate(s_val_recv_buffer,nnz_l_pexsi_aux,"s_val_recv_buffer",caller)

      call MPI_Alltoallv(s_val_send_buffer,send_count,send_displ,mpi_real8,&
                         s_val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                         mpi_comm_global,mpierr)

      deallocate(s_val_send_buffer)
   endif

   allocate(global_id(nnz_l_pexsi_aux))

   ! Compute global 1D id
   do i_val = 1,nnz_l_pexsi_aux
      global_id(i_val) = int(col_recv_buffer(i_val)-1,kind=8)*int(n_g_size,kind=8)+&
                         int(row_recv_buffer(i_val),kind=8)
   enddo

   ! Reorder
   if(n_elsi_calls == 1) then
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
   send_count(myid/pexsi_options%numPole+1) = nnz_l_pexsi_aux

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
                     1,mpi_integer,mpi_comm_global,mpierr)

   if(n_elsi_calls == 1) then
      nnz_l_pexsi = sum(recv_count,1)

      ! At this point only the first pole knows nnz_l_pexsi
      call MPI_Comm_split(mpi_comm_global,my_p_col_pexsi,my_p_row_pexsi,&
                          mpi_comm_aux_pexsi,mpierr)

      call MPI_Bcast(nnz_l_pexsi,1,mpi_integer,0,mpi_comm_aux_pexsi,mpierr)
   endif

   ! Set send and receive displacement
   send_displ_aux = 0
   recv_displ_aux = 0
   do i_proc = 0,n_procs-1
      send_displ(i_proc+1) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc+1)
      recv_displ(i_proc+1) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
   enddo

   ! Allocate PEXSI matrices
   ! Only the Hamiltonian needs to be reset everytime
   if(.not.allocated(ham_real_pexsi)) &
      call elsi_allocate(ham_real_pexsi,nnz_l_pexsi,"ham_real_pexsi",caller)
   ham_real_pexsi = 0.0d0

   if(.not.allocated(ovlp_real_pexsi)) &
      call elsi_allocate(ovlp_real_pexsi,nnz_l_pexsi,"ovlp_real_pexsi",caller)

   if(.not.allocated(row_ind_pexsi)) &
      call elsi_allocate(row_ind_pexsi,nnz_l_pexsi,"row_ind_pexsi",caller)

   if(.not.allocated(col_ptr_pexsi)) &
      call elsi_allocate(col_ptr_pexsi,(n_l_cols_pexsi+1),"col_ptr_pexsi",caller)

   ! Send and receive the packed data
   if(n_elsi_calls == 1) then
      ! Row id
      call elsi_allocate(row_send_buffer,nnz_l_pexsi,"row_send_buffer",caller)

      call MPI_Alltoallv(row_recv_buffer,send_count,send_displ,mpi_integer,&
                         row_send_buffer,recv_count,recv_displ,mpi_integer,&
                         mpi_comm_global,mpierr)

      deallocate(row_recv_buffer)

      ! Column id
      call elsi_allocate(col_send_buffer,nnz_l_pexsi,"col_send_buffer",caller)

      call MPI_Alltoallv(col_recv_buffer,send_count,send_displ,mpi_integer,&
                         col_send_buffer,recv_count,recv_displ,mpi_integer,&
                         mpi_comm_global,mpierr)

      deallocate(col_recv_buffer)

      ! Overlap value
      call MPI_Alltoallv(s_val_recv_buffer,send_count,send_displ,mpi_real8,&
                         ovlp_real_pexsi,recv_count,recv_displ,mpi_real8,&
                         mpi_comm_global,mpierr)

      deallocate(s_val_recv_buffer)
   endif

   ! Hamiltonian value
   call MPI_Alltoallv(h_val_recv_buffer,send_count,send_displ,mpi_real8,&
                      ham_real_pexsi,recv_count,recv_displ,mpi_real8,&
                      mpi_comm_global,mpierr)

   deallocate(h_val_recv_buffer)

   ! Only the first pole computes row index and column pointer
   if(n_elsi_calls == 1) then
      if(my_p_row_pexsi == 0) then
         ! Compute row index and column pointer
         i_col = col_send_buffer(1)-1
         do i_val = 1,nnz_l_pexsi
            row_ind_pexsi(i_val) = row_send_buffer(i_val)
            if(col_send_buffer(i_val) > i_col) then
               i_col = i_col+1
               col_ptr_pexsi(i_col-col_send_buffer(1)+1) = i_val
            endif
         enddo

         col_ptr_pexsi(n_l_cols_pexsi+1) = nnz_l_pexsi+1
      endif

      deallocate(row_send_buffer)
      deallocate(col_send_buffer)
   endif

   call elsi_stop_redistribution_time()

end subroutine

!>
!! This routine converts density matrix computed by PEXSI and stored
!! in 1D block distributed sparse CCS format to 2D block-cyclic
!! distributed dense format.
!!
subroutine elsi_pexsi_to_blacs_dm_small(D_out)

   implicit none
   include "mpif.h"

   real*8, intent(out) :: D_out(n_l_rows,n_l_cols) !< Density matrix to be converted

   integer :: i_row         !< Row counter
   integer :: i_col         !< Col counter
   integer :: i_val         !< Value counter
   integer :: j_val         !< Value counter
   integer :: i_proc        !< Process counter
   integer :: global_col_id !< Global column id
   integer :: global_row_id !< Global row id
   integer :: local_col_id  !< Local column id in 1D block distribution
   integer :: local_row_id  !< Local row id in 1D block distribution
   integer :: proc_col_id   !< Column id in process grid
   integer :: proc_row_id   !< Row id in process grid

   integer, allocatable :: dest(:)      !< Destination of each element
   integer, allocatable :: global_id(:) !< Global 1d id

   ! See documentation of MPI_Alltoallv
   real*8, allocatable :: val_send_buffer(:)  !< Send buffer for value
   integer, allocatable :: pos_send_buffer(:) !< Send buffer for global 1D id
   integer :: send_count(n_procs) !< Number of elements to send to each processor
   integer :: send_displ(n_procs) !< Displacement from which to take the outgoing data
   integer :: send_displ_aux      !< Auxiliary variable used to set displacement

   real*8, allocatable  :: val_recv_buffer(:) !< Receive buffer for value
   integer, allocatable :: pos_recv_buffer(:) !< Receive buffer for global 1D id
   integer :: recv_count(n_procs) !< Number of elements to receive from each processor
   integer :: recv_displ(n_procs) !< Displacement at which to place the incoming data
   integer :: recv_displ_aux      !< Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_dm_small"

   call elsi_start_redistribution_time()

   send_count = 0
   send_displ = 0
   recv_count = 0
   recv_displ = 0

   call elsi_allocate(val_send_buffer,nnz_l_pexsi,"val_send_buffer",caller)
   call elsi_allocate(pos_send_buffer,nnz_l_pexsi,"pos_send_buffer",caller)

   if(my_p_row_pexsi == 0) then
      call elsi_allocate(global_id,nnz_l_pexsi,"global_id",caller)
      call elsi_allocate(dest,nnz_l_pexsi,"dest",caller)

      i_col = 0
      ! Compute destination and global 1D id
      do i_val = 1,nnz_l_pexsi
         if(i_val == col_ptr_ccs(i_col+1) .and. i_col /= n_l_cols_pexsi) then
            i_col = i_col+1
         endif
         i_row = row_ind_ccs(i_val)

         ! Compute global id
         global_row_id = i_row
         global_col_id = i_col+myid*(n_g_size/n_p_per_pole_pexsi)
         global_id(i_val) = (global_col_id-1)*n_g_size+global_row_id

         ! Compute destination
         proc_row_id = mod((global_row_id-1)/n_b_rows,n_p_rows)
         proc_col_id = mod((global_col_id-1)/n_b_cols,n_p_cols)
         dest(i_val) = proc_col_id+proc_row_id*n_p_cols
      enddo

      j_val = 0

      ! Set send_count
      do i_proc = 0,n_procs-1
         do i_val = 1,nnz_l_pexsi
            if(dest(i_val) == i_proc) then
               j_val = j_val+1
               val_send_buffer(j_val) = den_mat_ccs(i_val)
               pos_send_buffer(j_val) = global_id(i_val)
               send_count(i_proc+1) = send_count(i_proc+1)+1
            endif
         enddo
      enddo

      deallocate(global_id)
      deallocate(dest)
   endif

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
                     1,mpi_integer,mpi_comm_global,mpierr)

   nnz_l = sum(recv_count,1)

   ! Set send and receive displacement
   send_displ_aux = 0
   recv_displ_aux = 0

   do i_proc = 0,n_procs-1
      send_displ(i_proc+1) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc+1)

      recv_displ(i_proc+1) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc+1)
   enddo

   ! Send and receive the packed data
   ! Value
   call elsi_allocate(val_recv_buffer,nnz_l,"val_recv_buffer",caller)

   call MPI_Alltoallv(val_send_buffer,send_count,send_displ,mpi_real8,&
                      val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                      mpi_comm_global,mpierr)

   deallocate(val_send_buffer)

   ! Position
   call elsi_allocate(pos_recv_buffer,nnz_l,"pos_recv_buffer",caller)

   call MPI_Alltoallv(pos_send_buffer,send_count,send_displ,mpi_integer,&
                      pos_recv_buffer,recv_count,recv_displ,mpi_integer,&
                      mpi_comm_global,mpierr)

   deallocate(pos_send_buffer)

   D_out = 0.0d0

   ! Unpack density matrix
   do i_val = 1,nnz_l
      ! Compute global 2d id
      global_col_id = (pos_recv_buffer(i_val)-1)/n_g_size+1
      global_row_id = mod(pos_recv_buffer(i_val)-1,n_g_size)+1

      ! Compute local 2d id
      local_row_id = (global_row_id-1)/(n_p_rows*n_b_rows)*n_b_rows&
                     +mod((global_row_id-1),n_b_rows)+1
      local_col_id = (global_col_id-1)/(n_p_cols*n_b_cols)*n_b_cols&
                     +mod((global_col_id-1),n_b_cols)+1

      ! Put value to correct position
      D_out(local_row_id,local_col_id) = val_recv_buffer(i_val)
   enddo

   deallocate(val_recv_buffer)
   deallocate(pos_recv_buffer)

   call elsi_stop_redistribution_time()

end subroutine

!>
!! This routine converts density matrix computed by PEXSI and stored
!! in 1D block distributed sparse CCS format to 2D block-cyclic
!! distributed dense format.
!!
!! * 2D index is used to deal with large matrices.
!!
subroutine elsi_pexsi_to_blacs_dm_large(D_out)

   implicit none
   include "mpif.h"

   real*8, intent(out) :: D_out(n_l_rows,n_l_cols) !< Density matrix to be converted

   integer :: i_row         !< Row counter
   integer :: i_col         !< Col counter
   integer :: i_val         !< Value counter
   integer :: j_val         !< Value counter
   integer :: i_proc        !< Process counter
   integer :: local_col_id  !< Local column id in 1D block distribution
   integer :: local_row_id  !< Local row id in 1D block distribution
   integer :: proc_col_id   !< Column id in process grid
   integer :: proc_row_id   !< Row id in process grid

   integer, allocatable :: global_col_id(:) !< Global column id
   integer, allocatable :: global_row_id(:) !< Global row id
   integer, allocatable :: dest(:)          !< Destination of each element

   ! See documentation of MPI_Alltoallv
   real*8, allocatable :: val_send_buffer(:)  !< Send buffer for value
   integer, allocatable :: row_send_buffer(:) !< Send buffer for global row id
   integer, allocatable :: col_send_buffer(:) !< Send buffer for global column id
   integer :: send_count(n_procs) !< Number of elements to send to each processor
   integer :: send_displ(n_procs) !< Displacement from which to take the outgoing data
   integer :: send_displ_aux      !< Auxiliary variable used to set displacement

   real*8, allocatable  :: val_recv_buffer(:) !< Receive buffer for value
   integer, allocatable :: row_recv_buffer(:) !< Receive buffer for global row id
   integer, allocatable :: col_recv_buffer(:) !< Receive buffer for global column id
   integer :: recv_count(n_procs) !< Number of elements to receive from each processor
   integer :: recv_displ(n_procs) !< Displacement at which to place the incoming data
   integer :: recv_displ_aux      !< Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_pexsi_to_blacs_dm_large"

   call elsi_start_redistribution_time()

   send_count = 0
   send_displ = 0
   recv_count = 0
   recv_displ = 0

   call elsi_allocate(val_send_buffer,nnz_l_pexsi,"val_send_buffer",caller)
   call elsi_allocate(row_send_buffer,nnz_l_pexsi,"row_send_buffer",caller)
   call elsi_allocate(col_send_buffer,nnz_l_pexsi,"col_send_buffer",caller)

   if(my_p_row_pexsi == 0) then
      call elsi_allocate(global_row_id,nnz_l_pexsi,"global_row_id",caller)
      call elsi_allocate(global_col_id,nnz_l_pexsi,"global_col_id",caller)
      call elsi_allocate(dest,nnz_l_pexsi,"dest",caller)

      i_col = 0
      ! Compute destination and global id
      do i_val = 1,nnz_l_pexsi
         if(i_val == col_ptr_ccs(i_col+1) .and. i_col /= n_l_cols_pexsi) then
            i_col = i_col+1
         endif
         i_row = row_ind_ccs(i_val)

         ! Compute global id
         global_row_id(i_val) = i_row
         global_col_id(i_val) = i_col+myid*(n_g_size/n_p_per_pole_pexsi)

         ! Compute destination
         proc_row_id = mod((global_row_id(i_val)-1)/n_b_rows,n_p_rows)
         proc_col_id = mod((global_col_id(i_val)-1)/n_b_cols,n_p_cols)
         dest(i_val) = proc_col_id+proc_row_id*n_p_cols
      enddo

      j_val = 0

      ! Set send_count
      do i_proc = 0,n_procs-1
         do i_val = 1,nnz_l_pexsi
            if(dest(i_val) == i_proc) then
               j_val = j_val+1
               val_send_buffer(j_val) = den_mat_ccs(i_val)
               row_send_buffer(j_val) = global_row_id(i_val)
               col_send_buffer(j_val) = global_col_id(i_val)
               send_count(i_proc+1) = send_count(i_proc+1)+1
            endif
         enddo
      enddo

      deallocate(global_row_id)
      deallocate(global_col_id)
      deallocate(dest)
   endif

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
                     1,mpi_integer,mpi_comm_global,mpierr)

   nnz_l = sum(recv_count,1)

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
   ! Value
   call elsi_allocate(val_recv_buffer,nnz_l,"val_recv_buffer",caller)

   call MPI_Alltoallv(val_send_buffer,send_count,send_displ,mpi_real8,&
                      val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                      mpi_comm_global,mpierr)

   deallocate(val_send_buffer)

   ! Row index
   call elsi_allocate(row_recv_buffer,nnz_l,"row_recv_buffer",caller)

   call MPI_Alltoallv(row_send_buffer,send_count,send_displ,mpi_integer,&
                      row_recv_buffer,recv_count,recv_displ,mpi_integer,&
                      mpi_comm_global,mpierr)

   deallocate(row_send_buffer)

   ! Column index
   call elsi_allocate(col_recv_buffer,nnz_l,"col_recv_buffer",caller)

   call MPI_Alltoallv(col_send_buffer,send_count,send_displ,mpi_integer,&
                      col_recv_buffer,recv_count,recv_displ,mpi_integer,&
                      mpi_comm_global,mpierr)

   deallocate(col_send_buffer)

   D_out = 0.0d0

   ! Unpack density matrix
   do i_val = 1,nnz_l
      ! Compute local 2d id
      local_row_id = (row_recv_buffer(i_val)-1)/(n_p_rows*n_b_rows)*n_b_rows&
                     +mod((row_recv_buffer(i_val)-1),n_b_rows)+1
      local_col_id = (col_recv_buffer(i_val)-1)/(n_p_cols*n_b_cols)*n_b_cols&
                     +mod((col_recv_buffer(i_val)-1),n_b_cols)+1

      ! Put value to correct position
      D_out(local_row_id,local_col_id) = val_recv_buffer(i_val)
   enddo

   deallocate(val_recv_buffer)
   deallocate(row_recv_buffer)
   deallocate(col_recv_buffer)

   call elsi_stop_redistribution_time()

end subroutine

!>
!! This routine interfaces to PEXSI.
!!
subroutine elsi_solve_evp_pexsi()

   implicit none
   include "mpif.h"

   real*8, save :: this_pexsi_tol = 1.0d-2

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_solve_evp_pexsi"

   call elsi_start_density_matrix_time()

   if(small_pexsi_tol) then
      pexsi_options%numElectronPEXSITolerance = this_pexsi_tol

      write(info_str,"(A,E10.1)") "  | Current tolerance of number of electrons: ",&
         this_pexsi_tol
      call elsi_statement_print(info_str)
   endif

   if(n_elsi_calls == 1) then
      pexsi_options%isSymbolicFactorize = 1
   else
      pexsi_options%isSymbolicFactorize = 0
   endif

   if(.not.allocated(e_den_mat_pexsi)) then
      call elsi_allocate(e_den_mat_pexsi,nnz_l_pexsi,"e_den_mat_pexsi",caller)
   endif
   e_den_mat_pexsi = 0.0d0

   if(.not.allocated(f_den_mat_pexsi)) then
      call elsi_allocate(f_den_mat_pexsi,nnz_l_pexsi,"f_den_mat_pexsi",caller)
   endif
   f_den_mat_pexsi = 0.0d0

   ! Load sparse matrices for PEXSI
   if(overlap_is_unit) then
      call f_ppexsi_load_real_hs_matrix(pexsi_plan,pexsi_options,n_g_size,nnz_g,&
                                        nnz_l_pexsi,n_l_cols_pexsi,col_ptr_ccs,&
                                        row_ind_ccs,ham_real_ccs,1,ovlp_real_ccs,&
                                        pexsi_info)
   else
      call f_ppexsi_load_real_hs_matrix(pexsi_plan,pexsi_options,n_g_size,nnz_g,&
                                        nnz_l_pexsi,n_l_cols_pexsi,col_ptr_ccs,&
                                        row_ind_ccs,ham_real_ccs,0,ovlp_real_ccs,&
                                        pexsi_info)
   endif

   if(pexsi_info /= 0) &
      call elsi_stop(" PEXSI not able to load H/S matrix. Exiting...",caller)

   if(pexsi_options%isInertiaCount == 0) then
      call elsi_statement_print("  PEXSI inertia counting skipped")
   endif

   ! Solve the eigenvalue problem
   call elsi_statement_print("  Starting PEXSI density matrix solver")

   call f_ppexsi_dft_driver(pexsi_plan,pexsi_options,n_electrons,mu_pexsi,&
                            n_electrons_pexsi,mu_min_inertia,mu_max_inertia,&
                            n_total_inertia_iter,n_total_pexsi_iter,pexsi_info)
       
   if(pexsi_info /= 0) &
      call elsi_stop(" PEXSI DFT driver not able to solve problem. Exiting...",caller)

   ! Turn off inertia counting if chemical potential does not change a lot
   if(abs(mu_pexsi-pexsi_options%mu0) > 5.0d-3) then
      pexsi_options%isInertiaCount = 1
   else
      pexsi_options%isInertiaCount = 0
   endif

   ! Use chemical potential in this step as initial guess for next step
   pexsi_options%mu0 = mu_pexsi

   if(small_pexsi_tol) then
      if(abs(n_electrons-n_electrons_pexsi) < this_pexsi_tol) then
         if(1.0d-1*this_pexsi_tol > final_pexsi_tol) then
            this_pexsi_tol = 1.0d-1*this_pexsi_tol
         else
            this_pexsi_tol = final_pexsi_tol
         endif
      endif
   endif

   ! Get the results
   if((my_p_row_pexsi == 0) .or. (storage == BLACS_DENSE)) then
      call f_ppexsi_retrieve_real_dft_matrix(pexsi_plan,den_mat_ccs,&
              e_den_mat_pexsi,f_den_mat_pexsi,e_tot_H,e_tot_S,f_tot,pexsi_info)
   endif

   if(pexsi_info /= 0) then
      call elsi_stop(" PEXSI not able to retrieve solution. Exiting...",caller)
   endif

   call MPI_Barrier(mpi_comm_global,mpierr)
   call elsi_stop_density_matrix_time()

end subroutine

!>
!! Set PEXSI variables to ELSI default.
!!
subroutine elsi_set_pexsi_default_options()

   implicit none

   ! Use the PEXSI Default options
   call f_ppexsi_set_default_options(pexsi_options)

   ! Use 1 process in ParMETIS for symbolic factorization
   pexsi_options%npSymbFact = 1

end subroutine

!>
!! Print PEXSI settings.
!!
subroutine elsi_print_pexsi_options()

   implicit none

   character*200 :: info_str

   write(info_str,"(A)") "  PEXSI settings (in the same unit of Hamiltonian):"
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Temperature ',F10.4)") &
      pexsi_options%temperature
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Spectral gap ',F10.4)") &
      pexsi_options%gap
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Number of poles ',I5)") &
      pexsi_options%numPole
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Max PEXSI iterations ',I5)") &
      pexsi_options%maxPEXSIIter
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Lower bound of chemical potential ',F10.4)") &
      pexsi_options%muMin0
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Upper bound of chemical potential ',F10.4)") &
      pexsi_options%muMax0
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Initial guess of chemical potential ',F10.4)") &
      pexsi_options%mu0
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Tolerance of chemical potential ',E10.1)") &
      pexsi_options%muInertiaTolerance
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Safeguard of chemical potential ',F10.4)") &
      pexsi_options%muPexsiSafeGuard
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Tolerance of number of electrons ',E10.1)") &
      pexsi_options%numElectronPEXSITolerance
   call elsi_statement_print(info_str)

end subroutine

end module ELSI_PEXSI
