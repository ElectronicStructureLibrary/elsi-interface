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
!! This module provides interfaces to SIPs.
!!

module ELSI_SIPS

   use ELSI_PRECISION, only : dp
   use ELSI_DIMENSIONS
   use ELSI_TIMERS
   use ELSI_UTILS
   use m_qetsc

   implicit none
   private

   public :: elsi_init_sips
   public :: elsi_solve_evp_sips
   public :: elsi_set_sips_default_options
   public :: elsi_print_sips_options
   public :: elsi_blacs_to_sips

contains

!=========================
! ELSI routines for SIPs
!=========================

!>
!! This routine initializes SIPs solver.
!!
subroutine elsi_init_sips()

   implicit none

   integer :: i

   character*40, parameter :: caller = "elsi_init_sips"

   if(n_elsi_calls == 1) then
      call init_sips()
   endif

   ! Number of slices
   n_p_per_slice_sips = 1
   n_slices = 1
   ! n_p_per_slice_sips cannot be larger than 16
   do i = 1,5
      if((mod(n_procs,n_p_per_slice_sips) == 0) .and. &
         (n_procs/n_p_per_slice_sips .le. n_states)) then
         n_slices = n_procs/n_p_per_slice_sips
      endif

      n_p_per_slice_sips = n_p_per_slice_sips*2
   enddo

   n_p_per_slice_sips = n_procs/n_slices

   ! SIPs uses a pure block distribution
   n_b_rows_sips = n_g_size

   ! The last process holds all remaining columns
   n_b_cols_sips = n_g_size/n_procs
   if(myid == n_procs-1) then
      n_b_cols_sips = n_g_size-(n_procs-1)*n_b_cols_sips
   endif

   n_l_rows_sips = n_b_rows_sips
   n_l_cols_sips = n_b_cols_sips

   if(.not.allocated(slices)) then
      call elsi_allocate(slices,n_slices+1,"slices",caller)
   endif

   if(.not.allocated(inertias)) then
      call elsi_allocate(inertias,n_slices+1,"inertias",caller)
   endif

   if(.not.allocated(shifts)) then
      call elsi_allocate(shifts,n_slices+1,"shifts",caller)
   endif

end subroutine

!>
!! This routine is a driver to convert matrix format and distribution
!! from BLACS to SIPs.
!!
subroutine elsi_blacs_to_sips(H_in,S_in)

   implicit none

   real(kind=dp), intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian matrix to be converted
   real(kind=dp), intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap matrix to be converted

   character*40, parameter :: caller = "elsi_blacs_to_sips"

   if(overlap_is_unit) then
      call elsi_stop(" SIPs with identity overlap matrix not yet available."//&
                     " Exiting...",caller)
   else
      if(n_g_size < 46340) then ! kind=4 integer works
         call elsi_blacs_to_sips_hs_small(H_in,S_in)
      else ! use kind=8 integer
!         call elsi_blacs_to_sips_hs_large(H_in,S_in)
      endif
   endif

end subroutine

!>
!! This routine converts Halmitonian and overlap matrix stored in
!! 2D block-cyclic distributed dense format to 1D block distributed
!! sparse CCS format, which can be used as input by SIPs.
!!
subroutine elsi_blacs_to_sips_hs_small(H_in,S_in)

   implicit none

   include "mpif.h"

   real(kind=dp), intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian matrix to be converted
   real(kind=dp), intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap matrix to be converted

   integer :: i_row !< Row counter
   integer :: i_col !< Col counter
   integer :: i_val,j_val !< Value counter
   integer :: i_proc !< Process counter
   integer :: global_col_id !< Global column id
   integer :: global_row_id !< Global row id
   integer :: local_col_id !< Local column id in 1D block distribution
   integer :: local_row_id !< Local row id in 1D block distribution
   integer :: tmp_int
   integer :: min_pos,min_id
   real(kind=dp) :: tmp_real
   integer, allocatable :: dest(:) !< Destination of each element

   ! See documentation of MPI_Alltoallv
   real(kind=dp), allocatable  :: h_val_send_buffer(:) !< Send buffer for Hamiltonian
   real(kind=dp), allocatable  :: s_val_send_buffer(:) !< Send buffer for overlap
   integer, allocatable :: pos_send_buffer(:) !< Send buffer for global 1D id
   integer, allocatable :: send_count(:) !< Number of elements to send to each processor
   integer, allocatable :: send_displ(:) !< Displacement from which to take the outgoing data
   real(kind=dp), allocatable  :: h_val_recv_buffer(:) !< Receive buffer for Hamiltonian
   real(kind=dp), allocatable  :: s_val_recv_buffer(:) !< Receive buffer for overlap
   integer, allocatable :: pos_recv_buffer(:) !< Receive buffer for global 1D id
   integer, allocatable :: recv_count(:) !< Number of elements to receive from each processor
   integer, allocatable :: recv_displ(:) !< Displacement at which to place the incoming data
   integer :: send_displ_aux !< Auxiliary variable used to set displacement 
   integer :: recv_displ_aux !< Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_blacs_to_sips_hs_small"

   call elsi_start_redistribution_time()

   if(n_elsi_calls == 1) then
      call elsi_get_local_nnz(H_in,n_l_rows,n_l_cols,nnz_l)
      call elsi_allocate(s_val_send_buffer,nnz_l,"s_val_send_buffer",caller)
   endif

   call elsi_allocate(dest,nnz_l,"dest",caller)
   call elsi_allocate(pos_send_buffer,nnz_l,"pos_send_buffer",caller)
   call elsi_allocate(h_val_send_buffer,nnz_l,"h_val_send_buffer",caller)
   call elsi_allocate(send_count,n_procs,"send_count",caller)

   ! Compute destination and global 1D id
   if(n_elsi_calls == 1) then
      i_val = 0
      do i_col = 1,n_l_cols
         do i_row = 1,n_l_rows
            if(abs(S_in(i_row, i_col)) > zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(global_col_id,i_col)
               call elsi_get_global_row(global_row_id,i_row)
               ! Compute destination
               dest(i_val) = (global_col_id-1)/(n_g_size/n_procs)
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
            if(abs(S_in(i_row, i_col)) > zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(global_col_id,i_col)
               call elsi_get_global_row(global_row_id,i_row)
               ! Compute destination
               dest(i_val) = (global_col_id-1)/(n_g_size/n_procs)
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

   call elsi_allocate(recv_count,n_procs,"recv_count",caller)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
                     1,mpi_integer,mpi_comm_global,mpierr)

   ! Set local/global number of nonzero
   nnz_l_sips = sum(recv_count,1)
   call MPI_Allreduce(nnz_l_sips,nnz_g,1,mpi_integer,mpi_sum,&
                      mpi_comm_global,mpierr)

   ! Set send and receive displacement
   call elsi_allocate(send_displ,n_procs,"send_displ",caller)
   call elsi_allocate(recv_displ,n_procs,"recv_displ",caller)

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
   call elsi_allocate(pos_recv_buffer,nnz_l_sips,"pos_recv_buffer",caller)

   call MPI_Alltoallv(pos_send_buffer,send_count,send_displ,mpi_integer,&
                      pos_recv_buffer,recv_count,recv_displ,mpi_integer,&
                      mpi_comm_global,mpierr)

   deallocate(pos_send_buffer)

   ! Hamiltonian value
   call elsi_allocate(h_val_recv_buffer,nnz_l_sips,"h_val_recv_buffer",caller)

   call MPI_Alltoallv(h_val_send_buffer,send_count,send_displ,mpi_real8,&
                      h_val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                      mpi_comm_global,mpierr)

   deallocate(h_val_send_buffer)

   ! Overlap value
   if(n_elsi_calls == 1) then
      call elsi_allocate(s_val_recv_buffer,nnz_l_sips,"s_val_recv_buffer",caller)

      call MPI_Alltoallv(s_val_send_buffer,send_count,send_displ,mpi_real8,&
                         s_val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                         mpi_comm_global,mpierr)

      deallocate(s_val_send_buffer)
   endif

   deallocate(send_count)
   deallocate(recv_count)
   deallocate(send_displ)
   deallocate(recv_displ)

   ! Unpack and reorder
   if(n_elsi_calls == 1) then
      do i_val = 1,nnz_l_sips
         min_id = minloc(pos_recv_buffer(i_val:nnz_l_sips),1)+i_val-1

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
      do i_val = 1,nnz_l_sips
         min_id = minloc(pos_recv_buffer(i_val:nnz_l_sips),1)+i_val-1

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
   if(.not.allocated(ham_real_sips)) &
      call elsi_allocate(ham_real_sips,nnz_l_sips,"ham_real_sips",caller)
   ham_real_sips = 0.0_dp

   if(.not.allocated(ovlp_real_sips)) &
      call elsi_allocate(ovlp_real_sips,nnz_l_sips,"ovlp_real_sips",caller)

   if(.not.allocated(row_ind_sips)) &
      call elsi_allocate(row_ind_sips,nnz_l_sips,"row_ind_sips",caller)

   if(.not.allocated(col_ptr_sips)) &
      call elsi_allocate(col_ptr_sips,(n_l_cols_sips+1),"col_ptr_sips",caller)

   ham_real_sips = h_val_recv_buffer
   deallocate(h_val_recv_buffer)

   if(n_elsi_calls == 1) then
      ovlp_real_sips = s_val_recv_buffer
      deallocate(s_val_recv_buffer)
   endif

   ! Compute row index and column pointer
   if(n_elsi_calls == 1) then
      i_col = (pos_recv_buffer(1)-1)/n_g_size
      do i_val = 1,nnz_l_sips
         row_ind_sips(i_val) = mod(pos_recv_buffer(i_val)-1,n_g_size)+1
         if((pos_recv_buffer(i_val)-1)/n_g_size+1 > i_col) then
            i_col = i_col+1
            col_ptr_sips(i_col-(pos_recv_buffer(1)-1)/n_g_size) = i_val
         endif
      enddo

      col_ptr_sips(n_l_cols_sips+1) = nnz_l_sips+1
   endif

   deallocate(pos_recv_buffer)

   ! Dummy
   nnz_l_pexsi = nnz_l_sips
   n_l_cols_pexsi = n_l_cols_sips

   call elsi_stop_redistribution_time()

end subroutine

!>
!! This routine interfaces to SIPs via QETSC.
!!
subroutine elsi_solve_evp_sips()

   implicit none
   include "mpif.h"

   integer :: i_state
   integer :: i_row
   integer :: i_row2
   integer :: i_col
   integer :: g_row
   integer :: g_row2
   integer :: this_p_col
   real(kind=dp), allocatable :: tmp_real(:)

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_solve_evp_sips"

   call elsi_start_generalized_evp_time()

   write(info_str,"(1X,' | Number of slices ',I7)") n_slices
   call elsi_statement_print(info_str)

   ! Solve the eigenvalue problem
   call elsi_statement_print("  Starting SIPs eigensolver")

   ! Load H and S matrices
   call sips_load_ham(n_g_size,n_l_cols_sips,nnz_l_sips,row_ind_ccs,&
                      col_ptr_ccs,ham_real_ccs)

   call sips_load_ovlp(n_g_size,n_l_cols_sips,nnz_l_sips,row_ind_ccs,&
                       col_ptr_ccs,ovlp_real_ccs)

   ! Initialize an eigenvalue problem
   if(.not.overlap_is_unit) then
      call sips_set_evp(1)
   else
      call sips_set_evp(0)
   endif

   ! Estimate the lower and upper bounds of eigenvalues
   call sips_get_interval(interval)

   ! Compute slicing
   call sips_get_slicing(n_slices,slicing_method,unbound,interval,&
                         0.0_dp,0.0_dp,slices,n_states,eval)

   ! Solve eigenvalue problem
   call sips_solve_evp(n_states,n_slices,slices,n_solve_steps)

   ! Get eigenvalues
   call sips_get_eigenvalues(n_states,eval(1:n_states))

   ! Get and distribute eigenvectors
   call elsi_allocate(tmp_real,n_g_size,"tmp_real",caller)

   evec_real = 0.0_dp

   do i_state = 1,n_states
      call sips_get_eigenvectors(n_g_size,i_state,tmp_real)

      this_p_col = mod((i_state-1)/n_b_cols,n_procs)

      if(my_p_col == this_p_col) then
         i_col = (i_state-1)/(n_p_cols*n_b_cols)*n_b_cols&
                 +mod((i_state-1),n_b_cols)+1

         do i_row = 1,n_l_rows,n_b_rows
            i_row2 = i_row+n_b_rows-1

            call elsi_get_global_row(g_row,i_row)
            call elsi_get_global_row(g_row2,i_row2)

            if(g_row2 > n_g_size) then
               g_row2 = n_g_size
               i_row2 = i_row+g_row2-g_row
            endif

            evec_real(i_row:i_row2,i_col) = tmp_real(g_row:g_row2)
         enddo
      endif
   enddo

   deallocate(tmp_real)

   call MPI_Barrier(mpi_comm_global,mpierr)
   call elsi_stop_generalized_evp_time()

end subroutine

!>
!! Set SIPs variables to ELSI default.
!!
subroutine elsi_set_sips_default_options()

   implicit none

   character*40, parameter :: caller = "elsi_set_sips_default_options"

   !< Type of slices
   !! 0 = Equally spaced subintervals
   !! 1 = K-meaans after equally spaced subintervals
   !! 2 = Equally populated subintervals
   !! 3 = K-means after equally populated
   slicing_method = 0

   !< Extra inertia computations before solve?
   !! 0 = No
   !! 1 = Yes
   inertia_option = 0

   !< How to bound the left side of the interval
   !! 0 = Bounded
   !! 1 = -infinity
   unbound = 0

   !< Small buffer to expand the eigenvalue interval
   !! Smaller values improve performance if eigenvalue range known
   slice_buffer = 0.1_dp

end subroutine

!>
!! Print SIPs settings.
!!
subroutine elsi_print_sips_options()

   implicit none

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_print_sips_options"

   write(info_str,"(A)") "  SIPs settings:"
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Slicing method ',I2)") slicing_method
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Inertia option ',I2)") inertia_option
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Left bound ',I2)") unbound
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Slice buffer ',F10.4)") slice_buffer
   call elsi_statement_print(info_str)

end subroutine

end module ELSI_SIPS
