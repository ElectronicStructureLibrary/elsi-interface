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
   use m_qetsc_tools

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

!TODO: set_ham, set_ovlp, update_ham

!>
!! This routine initializes SIPs solver.
!!
subroutine elsi_init_sips()

   implicit none

   character*40, parameter :: caller = "elsi_init_sips"

   call initialize_qetsc()

   call elsi_allocate(slices,n_slices+1,"slices",caller)
   call elsi_allocate(inertias,n_slices+1,"inertias",caller)
   call elsi_allocate(shifts,n_slices+1,"shifts",caller)

   ! Initialize an eigenvalue problem
   call set_eps(ham_sips,ovlp_sips)

end subroutine

!>
!! This routine interfaces to SIPs via QETSC.
!!
subroutine elsi_solve_evp_sips()

   implicit none
   include "mpif.h"

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_solve_evp_sips"

   call elsi_start_solve_evp_time()

   ! Solve the eigenvalue problem
   call elsi_statement_print("  Starting SIPs eigensolver")

   ! Estimate the lower and upper bounds of eigenvalues
   interval = get_eps_interval()

   ! Expand by a smaller buffer
   ! Second SCF
!   interval(1) = eval(1) - slice_buffer
!   interval(2) = eval(n_states) + slice_buffer

   ! Compute slicing
   call compute_subintervals(n_slices,0,unbound,interval,&
                             0.0d0,0.0d0,slices)
   call set_eps_subintervals(n_slices,slices)

   if(inertia_option > 0 .and. n_slices > 1) then
      ! Computes inertias at shifts
      call run_eps_inertias_check(unbound,n_states,n_slices,slices,&
                                  shifts,inertias,n_inertia_counts)

      ! Computes estimated eigenvalues based on inertias
      call inertias_to_eigenvalues(n_slices+1,n_states,slice_buffer,&
                                   shifts,inertias,eval(1:n_states))

      ! Compute slicing
      call compute_subintervals(n_slices,subtype,unbound,interval,&
                                0.0d0,0.0d0,slices,eval(1:n_states))
      call set_eps_subintervals(n_slices,slices)
   endif

   ! Solve eigenvalue problem
   call solve_eps_check(n_states,n_slices,slices,n_solves)

   ! Get results
   eval = get_eps_eigenvalues(n_states)

!DEBUG
if(myid == 0) print *,eval
!DEBUG

   call MPI_Barrier(mpi_comm_global,mpierr)
   call elsi_stop_solve_evp_time()

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
      if(MOD(n_procs,n_p_per_slice_sips) == 0) then
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
