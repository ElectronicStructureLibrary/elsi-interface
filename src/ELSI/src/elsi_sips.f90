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

   use ELSI_CONSTANTS
   use ELSI_DIMENSIONS, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4
   use ELSI_TIMERS
   use ELSI_UTILS
   use m_qetsc

   implicit none
   private

   public :: elsi_init_sips
   public :: elsi_solve_evp_sips
   public :: elsi_set_sips_default_options
   public :: elsi_print_sips_options

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
