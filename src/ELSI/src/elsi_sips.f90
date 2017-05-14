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

   use ELSI_PRECISION, only: r8,i4
   use ELSI_CONSTANTS
   use ELSI_DIMENSIONS, only: elsi_handle
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
subroutine elsi_init_sips(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   integer(kind=i4) :: i

   character*40, parameter :: caller = "elsi_init_sips"

   if(elsi_h%n_elsi_calls == 1) then
      call init_sips()
   endif

   ! Number of slices
   elsi_h%n_p_per_slice_sips = 1
   elsi_h%n_slices = 1
   ! n_p_per_slice_sips cannot be larger than 16
   do i = 1,5
      if((mod(elsi_h%n_procs,elsi_h%n_p_per_slice_sips) == 0) .and. &
         (elsi_h%n_procs/elsi_h%n_p_per_slice_sips .le. elsi_h%n_states)) then
         elsi_h%n_slices = elsi_h%n_procs/elsi_h%n_p_per_slice_sips
      endif

      elsi_h%n_p_per_slice_sips = elsi_h%n_p_per_slice_sips*2
   enddo

   elsi_h%n_p_per_slice_sips = elsi_h%n_procs/elsi_h%n_slices

   ! SIPs uses a pure block distribution
   elsi_h%n_b_rows_sips = elsi_h%n_g_size

   ! The last process holds all remaining columns
   elsi_h%n_b_cols_sips = elsi_h%n_g_size/elsi_h%n_procs
   if(elsi_h%myid == elsi_h%n_procs-1) then
      elsi_h%n_b_cols_sips = elsi_h%n_g_size-(elsi_h%n_procs-1)*elsi_h%n_b_cols_sips
   endif

   elsi_h%n_l_rows_sips = elsi_h%n_b_rows_sips
   elsi_h%n_l_cols_sips = elsi_h%n_b_cols_sips

   if(.not. allocated(elsi_h%slices)) then
      call elsi_allocate(elsi_h,elsi_h%slices,elsi_h%n_slices+1,"slices",caller)
   endif

   if(.not. allocated(elsi_h%inertias)) then
      call elsi_allocate(elsi_h,elsi_h%inertias,elsi_h%n_slices+1,"inertias",caller)
   endif

   if(.not. allocated(elsi_h%shifts)) then
      call elsi_allocate(elsi_h,elsi_h%shifts,elsi_h%n_slices+1,"shifts",caller)
   endif

end subroutine

!>
!! This routine is a driver to convert matrix format and distribution
!! from BLACS to SIPs.
!!
subroutine elsi_blacs_to_sips(elsi_h,H_in,S_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   real(kind=r8),     intent(in)    :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian matrix to be converted
   real(kind=r8),     intent(in)    :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap matrix to be converted

   character*40, parameter :: caller = "elsi_blacs_to_sips"

   if(elsi_h%overlap_is_unit) then
      call elsi_stop(" SIPs with identity overlap matrix not yet available."//&
                     " Exiting...",elsi_h,caller)
   else
      if(elsi_h%n_g_size < 46340) then ! kind=4 integer works
         call elsi_blacs_to_sips_hs_small(elsi_h,H_in,S_in)
      else ! use kind=8 integer
!         call elsi_blacs_to_sips_hs_large(elsi_h,H_in,S_in)
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

   type(elsi_handle), intent(inout) :: elsi_h
   real(kind=r8),     intent(in)    :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Hamiltonian matrix to be converted
   real(kind=r8),     intent(in)    :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols) !< Overlap matrix to be converted

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: i_row !< Row counter
   integer(kind=i4) :: i_col !< Col counter
   integer(kind=i4) :: i_val,j_val !< Value counter
   integer(kind=i4) :: i_proc !< Process counter
   integer(kind=i4) :: global_col_id !< Global column id
   integer(kind=i4) :: global_row_id !< Global row id
   integer(kind=i4) :: local_col_id !< Local column id in 1D block distribution
   integer(kind=i4) :: local_row_id !< Local row id in 1D block distribution
   integer(kind=i4) :: tmp_int
   integer(kind=i4) :: min_pos,min_id
   real(kind=r8)    :: tmp_real
   integer(kind=i4), allocatable :: dest(:) !< Destination of each element

   ! See documentation of MPI_Alltoallv
   real(kind=r8),    allocatable  :: h_val_send_buffer(:) !< Send buffer for Hamiltonian
   real(kind=r8),    allocatable  :: s_val_send_buffer(:) !< Send buffer for overlap
   integer(kind=i4), allocatable :: pos_send_buffer(:) !< Send buffer for global 1D id
   integer(kind=i4), allocatable :: send_count(:) !< Number of elements to send to each processor
   integer(kind=i4), allocatable :: send_displ(:) !< Displacement from which to take the outgoing data
   real(kind=r8),    allocatable  :: h_val_recv_buffer(:) !< Receive buffer for Hamiltonian
   real(kind=r8),    allocatable  :: s_val_recv_buffer(:) !< Receive buffer for overlap
   integer(kind=i4), allocatable :: pos_recv_buffer(:) !< Receive buffer for global 1D id
   integer(kind=i4), allocatable :: recv_count(:) !< Number of elements to receive from each processor
   integer(kind=i4), allocatable :: recv_displ(:) !< Displacement at which to place the incoming data
   integer(kind=i4) :: send_displ_aux !< Auxiliary variable used to set displacement 
   integer(kind=i4) :: recv_displ_aux !< Auxiliary variable used to set displacement

   character*40, parameter :: caller = "elsi_blacs_to_sips_hs_small"

   call elsi_start_redistribution_time(elsi_h)

   if(elsi_h%n_elsi_calls == 1) then
      call elsi_get_local_nnz(elsi_h,H_in,elsi_h%n_l_rows,elsi_h%n_l_cols,elsi_h%nnz_l)
      call elsi_allocate(elsi_h,s_val_send_buffer,elsi_h%nnz_l,"s_val_send_buffer",caller)
   endif

   call elsi_allocate(elsi_h,dest,elsi_h%nnz_l,"dest",caller)
   call elsi_allocate(elsi_h,pos_send_buffer,elsi_h%nnz_l,"pos_send_buffer",caller)
   call elsi_allocate(elsi_h,h_val_send_buffer,elsi_h%nnz_l,"h_val_send_buffer",caller)
   call elsi_allocate(elsi_h,send_count,elsi_h%n_procs,"send_count",caller)

   ! Compute destination and global 1D id
   if(elsi_h%n_elsi_calls == 1) then
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(S_in(i_row, i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest(i_val) = (global_col_id-1)/(elsi_h%n_g_size/elsi_h%n_procs)
               ! The last process may take more
               if(dest(i_val) > (elsi_h%n_procs-1)) dest(i_val) = elsi_h%n_procs-1
               ! Compute the global id
               ! Pack global id and data into buffers
               pos_send_buffer(i_val) = (global_col_id-1)*elsi_h%n_g_size+global_row_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               s_val_send_buffer(i_val) = S_in(i_row,i_col)
               ! Set send_count
               send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
           endif
        enddo
      enddo
   else
      i_val = 0
      do i_col = 1,elsi_h%n_l_cols
         do i_row = 1,elsi_h%n_l_rows
            if(abs(S_in(i_row, i_col)) > elsi_h%zero_threshold) then
               i_val = i_val+1
               call elsi_get_global_col(elsi_h,global_col_id,i_col)
               call elsi_get_global_row(elsi_h,global_row_id,i_row)
               ! Compute destination
               dest(i_val) = (global_col_id-1)/(elsi_h%n_g_size/elsi_h%n_procs)
               ! The last process may take more
               if(dest(i_val) > (elsi_h%n_procs-1)) dest(i_val) = elsi_h%n_procs-1
               ! Compute the global id
               ! Pack global id and data into buffers
               pos_send_buffer(i_val) = (global_col_id-1)*elsi_h%n_g_size+global_row_id
               h_val_send_buffer(i_val) = H_in(i_row,i_col)
               ! Set send_count
               send_count(dest(i_val)+1) = send_count(dest(i_val)+1)+1
           endif
        enddo
      enddo
   endif

   deallocate(dest)

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
                      pos_recv_buffer,recv_count,recv_displ,mpi_integer,&
                      elsi_h%mpi_comm,mpierr)

   deallocate(pos_send_buffer)

   ! Hamiltonian value
   call elsi_allocate(elsi_h,h_val_recv_buffer,elsi_h%nnz_l_sips,"h_val_recv_buffer",caller)

   call MPI_Alltoallv(h_val_send_buffer,send_count,send_displ,mpi_real8,&
                      h_val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                      elsi_h%mpi_comm,mpierr)

   deallocate(h_val_send_buffer)

   ! Overlap value
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,s_val_recv_buffer,elsi_h%nnz_l_sips,"s_val_recv_buffer",caller)

      call MPI_Alltoallv(s_val_send_buffer,send_count,send_displ,mpi_real8,&
                         s_val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                         elsi_h%mpi_comm,mpierr)

      deallocate(s_val_send_buffer)
   endif

   deallocate(send_count)
   deallocate(recv_count)
   deallocate(send_displ)
   deallocate(recv_displ)

   ! Unpack and reorder
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
   deallocate(h_val_recv_buffer)

   if(elsi_h%n_elsi_calls == 1) then
      elsi_h%ovlp_real_sips = s_val_recv_buffer
      deallocate(s_val_recv_buffer)
   endif

   ! Compute row index and column pointer
   if(elsi_h%n_elsi_calls == 1) then
      i_col = (pos_recv_buffer(1)-1)/elsi_h%n_g_size
      do i_val = 1,elsi_h%nnz_l_sips
         elsi_h%row_ind_sips(i_val) = mod(pos_recv_buffer(i_val)-1,elsi_h%n_g_size)+1
         if((pos_recv_buffer(i_val)-1)/elsi_h%n_g_size+1 > i_col) then
            i_col = i_col+1
            elsi_h%col_ptr_sips(i_col-(pos_recv_buffer(1)-1)/elsi_h%n_g_size) = i_val
         endif
      enddo

      elsi_h%col_ptr_sips(elsi_h%n_l_cols_sips+1) = elsi_h%nnz_l_sips+1
   endif

   deallocate(pos_recv_buffer)

   ! Dummy
   elsi_h%nnz_l_pexsi = elsi_h%nnz_l_sips
   elsi_h%n_l_cols_pexsi = elsi_h%n_l_cols_sips

   call elsi_stop_redistribution_time(elsi_h)

end subroutine

!>
!! This routine interfaces to SIPs via QETSC.
!!
subroutine elsi_solve_evp_sips(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h

   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_row2
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_row
   integer(kind=i4) :: g_row2
   integer(kind=i4) :: this_p_col
   integer(kind=i4) :: mpierr
   real(kind=r8), allocatable :: tmp_real(:)

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_solve_evp_sips"

   call elsi_start_generalized_evp_time(elsi_h)

   write(info_str,"(1X,' | Number of slices ',I7)") elsi_h%n_slices
   call elsi_statement_print(info_str,elsi_h)

   ! Solve the eigenvalue problem
   call elsi_statement_print("  Starting SIPs eigensolver",elsi_h)

   ! Load H and S matrices
   call sips_load_ham(elsi_h%n_g_size,elsi_h%n_l_cols_sips,elsi_h%nnz_l_sips,&
                     elsi_h%row_ind_ccs,elsi_h%col_ptr_ccs,elsi_h%ham_real_ccs)

   call sips_load_ovlp(elsi_h%n_g_size,elsi_h%n_l_cols_sips,elsi_h%nnz_l_sips,&
                       elsi_h%row_ind_ccs,elsi_h%col_ptr_ccs,elsi_h%ovlp_real_ccs)

   ! Initialize an eigenvalue problem
   if(.not. elsi_h%overlap_is_unit) then
      call sips_set_evp(1)
   else
      call sips_set_evp(0)
   endif

   ! Estimate the lower and upper bounds of eigenvalues
   call sips_get_interval(elsi_h%interval)

   ! Compute slicing
   call sips_get_slicing(elsi_h%n_slices,elsi_h%slicing_method,elsi_h%unbound,&
                         elsi_h%interval,0.0_r8,0.0_r8,elsi_h%slices,&
                         elsi_h%n_states,elsi_h%eval)

   ! Solve eigenvalue problem
   call sips_solve_evp(elsi_h%n_states,elsi_h%n_slices,&
                       elsi_h%slices,elsi_h%n_solve_steps)

   ! Get eigenvalues
   call sips_get_eigenvalues(elsi_h%n_states,elsi_h%eval(1:elsi_h%n_states))

   ! Get and distribute eigenvectors
   call elsi_allocate(elsi_h,tmp_real,elsi_h%n_g_size,"tmp_real",caller)

   elsi_h%evec_real = 0.0_r8

   do i_state = 1,elsi_h%n_states
      call sips_get_eigenvectors(elsi_h%n_g_size,i_state,tmp_real)

      this_p_col = mod((i_state-1)/elsi_h%n_b_cols,elsi_h%n_procs)

      if(elsi_h%my_p_col == this_p_col) then
         i_col = (i_state-1)/(elsi_h%n_p_cols*elsi_h%n_b_cols)*elsi_h%n_b_cols&
                 +mod((i_state-1),elsi_h%n_b_cols)+1

         do i_row = 1,elsi_h%n_l_rows,elsi_h%n_b_rows
            i_row2 = i_row+elsi_h%n_b_rows-1

            call elsi_get_global_row(elsi_h,g_row,i_row)
            call elsi_get_global_row(elsi_h,g_row2,i_row2)

            if(g_row2 > elsi_h%n_g_size) then
               g_row2 = elsi_h%n_g_size
               i_row2 = i_row+g_row2-g_row
            endif

            elsi_h%evec_real(i_row:i_row2,i_col) = tmp_real(g_row:g_row2)
         enddo
      endif
   enddo

   deallocate(tmp_real)

   call MPI_Barrier(elsi_h%mpi_comm,mpierr)
   call elsi_stop_generalized_evp_time(elsi_h)

end subroutine

!>
!! Set SIPs variables to ELSI default.
!!
subroutine elsi_set_sips_default_options(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   character*40, parameter :: caller = "elsi_set_sips_default_options"

   !< Type of slices
   !! 0 = Equally spaced subintervals
   !! 1 = K-meaans after equally spaced subintervals
   !! 2 = Equally populated subintervals
   !! 3 = K-means after equally populated
   elsi_h%slicing_method = 0

   !< Extra inertia computations before solve?
   !! 0 = No
   !! 1 = Yes
   elsi_h%inertia_option = 0

   !< How to bound the left side of the interval
   !! 0 = Bounded
   !! 1 = -infinity
   elsi_h%unbound = 0

   !< Small buffer to expand the eigenvalue interval
   !! Smaller values improve performance if eigenvalue range known
   elsi_h%slice_buffer = 0.1_r8

end subroutine

!>
!! Print SIPs settings.
!!
subroutine elsi_print_sips_options(elsi_h)

   implicit none

   type(elsi_handle), intent(in) :: elsi_h

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_print_sips_options"

   write(info_str,"(A)") "  SIPs settings:"
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Slicing method ',I2)") elsi_h%slicing_method
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Inertia option ',I2)") elsi_h%inertia_option
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Left bound ',I2)") elsi_h%unbound
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Slice buffer ',F10.4)") elsi_h%slice_buffer
   call elsi_statement_print(info_str,elsi_h)

end subroutine

end module ELSI_SIPS
