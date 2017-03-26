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

   character*40, parameter :: caller = "elsi_init_sips"

   if(n_elsi_calls == 1) then
      call initialize_qetsc()
   endif

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

   real*8, intent(in) :: H_in(n_l_rows,n_l_cols) !< Hamiltonian matrix to be converted
   real*8, intent(in) :: S_in(n_l_rows,n_l_cols) !< Overlap matrix to be converted

   character*40, parameter :: caller = "elsi_blacs_to_sips"

   if(overlap_is_unit) then
      !TODO
      call elsi_stop(" SIPs with indentity overlap matrix not yet available."//&
                     " Exiting...",caller)
   else
      if(n_g_size < 46340) then ! kind=4 integer works
!         call elsi_blacs_to_sips_hs_small(H_in,S_in)
      else ! use kind=8 integer
!         call elsi_blacs_to_sips_hs_large(H_in,S_in)
      endif
   endif

end subroutine

!>
!! This routine interfaces to SIPs via QETSC.
!!
subroutine elsi_solve_evp_sips()

   implicit none
   include "mpif.h"

   integer :: istart
   integer :: iend

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_solve_evp_sips"

   call elsi_start_solve_evp_time()

   ! Solve the eigenvalue problem
   call elsi_statement_print("  Starting SIPs eigensolver")

   ! Load H and S matrices
   call load_elsi_ham(n_g_size,n_l_cols_sips,nnz_l_sips,row_ind_ccs,&
                      col_ptr_ccs,ham_real_ccs,istart,iend)

   call load_elsi_ovlp(istart,iend,n_l_cols_sips,nnz_l_sips,&
                       row_ind_ccs,col_ptr_ccs,ovlp_real_ccs)

   ! Initialize an eigenvalue problem
   if(.not.overlap_is_unit) then
      call set_eps(1)
   else
      call set_eps(0)
   endif

   ! Estimate the lower and upper bounds of eigenvalues
   interval = get_eps_interval()

   ! Compute slicing
   call compute_subintervals(n_slices,slicing_method,unbound,interval,&
                             0.0d0,0.0d0,slices)
   call set_eps_subintervals(n_slices,slices)

   ! Solve eigenvalue problem
   call solve_eps_check(n_states,n_slices,slices,n_solve_steps)

   ! Get results
   eval = get_eps_eigenvalues(n_states)

   call MPI_Barrier(mpi_comm_global,mpierr)
   call elsi_stop_solve_evp_time()

   ! TEST
   if(myid == 0) print *,eval
   call elsi_stop("End testing SIPs solver. Exiting...",caller)

end subroutine

!>
!! Set SIPs variables to ELSI default.
!!
subroutine elsi_set_sips_default_options()

   implicit none

   integer :: i

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

   !< Number of slices
   n_p_per_slice_sips = 1
   do i = 1,5
      if(mod(n_procs,n_p_per_slice_sips) == 0) then
         n_slices = n_procs/n_p_per_slice_sips
      endif
      n_p_per_slice_sips = n_p_per_slice_sips*2
   enddo

   !< Small buffer to expand the eigenvalue interval
   !! Smaller values improve performance if eigenvalue range known
   slice_buffer = 0.1d0

end subroutine

!>
!! Print SIPs settings.
!!
subroutine elsi_print_sips_options()

   implicit none

   character*200 :: info_str

   write(info_str,"(A)") "  SIPs settings:"
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Slicing method ',I2)") slicing_method
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Inertia option ',I2)") inertia_option
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Left bound ',I2)") unbound
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Number of slices ',I7)") n_slices
   call elsi_statement_print(info_str)

   write(info_str,"(1X,' | Slice buffer ',F10.4)") slice_buffer
   call elsi_statement_print(info_str)

end subroutine

end module ELSI_SIPS
