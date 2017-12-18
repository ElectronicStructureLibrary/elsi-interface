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

   use ELSI_CONSTANTS, only: UNSET
   use ELSI_DATATYPE
   use ELSI_MALLOC
   use ELSI_PRECISION, only: r8,i4
   use ELSI_TIMINGS,   only: elsi_get_time
   use ELSI_UTILS
   use M_QETSC

   implicit none

   private

   public :: elsi_set_sips_default
   public :: elsi_init_sips
   public :: elsi_solve_evp_sips
   public :: elsi_sips_to_blacs_ev

contains

!>
!! This routine initializes SIPs.
!!
subroutine elsi_init_sips(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*200 :: info_str

   character*40, parameter :: caller = "elsi_init_sips"

   if(e_h%n_elsi_calls == e_h%sips_n_elpa+1) then
      call initialize_qetsc()

      if(e_h%n_slices == UNSET) then
         ! TODO: Number of slices
         e_h%np_per_slice = 1
         e_h%n_slices     = e_h%n_procs
      endif

      ! 1D block distribution
      e_h%n_lcol_sp = e_h%n_basis/e_h%n_procs

      ! The last process holds all remaining columns
      if(e_h%myid == e_h%n_procs-1) then
         e_h%n_lcol_sp = e_h%n_basis-(e_h%n_procs-1)*e_h%n_lcol_sp
      endif

      call elsi_allocate(e_h,e_h%slices,e_h%n_slices+1,"slices",caller)

      e_h%sips_started = .true.
   endif

   write(info_str,"('  | Number of slices          ',I10)") e_h%n_slices
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine interfaces to SIPs via QETSC.
!!
subroutine elsi_solve_evp_sips(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: n_solve_steps
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   integer(kind=i4) :: mpierr
   character*200    :: info_str

   real(kind=r8),    allocatable :: shifts(:)
   integer(kind=i4), allocatable :: inertias(:)

   character*40, parameter :: caller = "elsi_solve_evp_sips"

   ! Solve the eigenvalue problem
   call elsi_say(e_h,"  Starting SIPs eigensolver")

   call elsi_get_time(e_h,t0)

   if(e_h%n_elsi_calls == e_h%sips_n_elpa+1) then
      ! Load H matrix
      call eps_load_ham(e_h%n_basis,e_h%n_lcol_sp,e_h%nnz_l_sp,&
              e_h%row_ind_ccs,e_h%col_ptr_ccs,e_h%ham_real_ccs)

      if(.not. e_h%ovlp_is_unit) then
         ! Load S matrix
         call eps_load_ovlp(e_h%n_basis,e_h%n_lcol_sp,e_h%nnz_l_sp,&
                 e_h%row_ind_ccs,e_h%col_ptr_ccs,e_h%ovlp_real_ccs)

         call set_eps(e_h%ev_min,e_h%ev_max,math,mats)
      else
         call set_eps(e_h%ev_min,e_h%ev_max,math)
      endif
   else ! n_elsi_calls > sips_n_elpa+1
      ! Update H matrix
      call eps_update_ham(e_h%n_basis,e_h%n_lcol_sp,e_h%nnz_l_sp,&
              e_h%row_ind_ccs,e_h%col_ptr_ccs,e_h%ham_real_ccs)

      call update_eps(e_h%n_slices)
   endif

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished loading matrices')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

   if(e_h%sips_n_elpa < 1 .and. e_h%n_elsi_calls == 1) then
      ! Estimate the lower and upper bounds of eigenvalues
      e_h%interval = get_eps_interval()

      ! Compute slicing
      call compute_subintervals(e_h%n_slices,0,e_h%unbound,e_h%interval,0.0_r8,&
              0.0_r8,e_h%slices)

      ! Run inertia counting
      if(e_h%inertia_option > 0 .and. e_h%n_slices > 1) then
         call elsi_get_time(e_h,t0)

         call elsi_allocate(e_h,inertias,e_h%n_slices+1,"inertias",caller)
         call elsi_allocate(e_h,shifts,e_h%n_slices+1,"shifts",caller)

         call run_eps_inertias_check(e_h%unbound,e_h%n_states,e_h%n_slices,&
                 e_h%slices,shifts,inertias,n_solve_steps)

         call inertias_to_eigenvalues(e_h%n_slices+1,e_h%n_states,&
                 e_h%slice_buffer,shifts,inertias,e_h%eval(1:e_h%n_states))

         call compute_subintervals(e_h%n_slices,e_h%slicing_method,e_h%unbound,&
                 e_h%interval,0.0_r8,0.0_r8,e_h%slices,e_h%eval(1:e_h%n_states))

         call elsi_deallocate(e_h,inertias,"inertias")
         call elsi_deallocate(e_h,shifts,"shifts")

         call elsi_get_time(e_h,t1)

         write(info_str,"('  Finished inertia counting')")
         call elsi_say(e_h,info_str)
         write(info_str,"('  | Time :',F10.3,' s')") t1-t0
         call elsi_say(e_h,info_str)
      endif
   else
      e_h%interval(1) = e_h%eval(1)-e_h%slice_buffer
      e_h%interval(2) = e_h%eval(e_h%n_states)+e_h%slice_buffer

      call compute_subintervals(e_h%n_slices,e_h%slicing_method,e_h%unbound,&
              e_h%interval,e_h%slice_buffer,1.0e-6_r8,e_h%slices,&
              e_h%eval(1:e_h%n_states))
   endif

   call set_eps_subintervals(e_h%n_slices,e_h%slices)

   call elsi_get_time(e_h,t0)

   ! Solve
   call solve_eps_check(e_h%n_states,e_h%n_slices,e_h%slices,n_solve_steps)

   ! Get eigenvalues
   e_h%eval(1:e_h%n_states) = get_eps_eigenvalues(e_h%n_states)

   call MPI_Barrier(e_h%mpi_comm,mpierr)

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished solving generalized eigenproblem')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine gets the eigenvectors computed by SIPs and distributes them in a
!! 2D block-cyclic fashion.
!!
subroutine elsi_sips_to_blacs_ev(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_row2
   integer(kind=i4) :: i_col
   integer(kind=i4) :: g_row
   integer(kind=i4) :: g_row2
   integer(kind=i4) :: this_pcol
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   real(kind=r8), allocatable :: tmp_real(:)

   character*40, parameter :: caller = "elsi_distribute_ev"

   call elsi_get_time(e_h,t0)

   call elsi_allocate(e_h,tmp_real,e_h%n_basis,"tmp_real",caller)

   e_h%evec_real = 0.0_r8

   do i_state = 1,e_h%n_states
      call get_eps_eigenvectors(e_h%n_basis,i_state,tmp_real)

      this_pcol = mod((i_state-1)/e_h%blk_col,e_h%n_pcol)

      if(e_h%my_pcol == this_pcol) then
         i_col = (i_state-1)/(e_h%n_pcol*e_h%blk_col)*e_h%blk_col+&
                    mod((i_state-1),e_h%blk_col)+1

         do i_row = 1,e_h%n_lrow,e_h%blk_row
            i_row2 = i_row+e_h%blk_row-1

            call elsi_get_global_row(e_h,g_row,i_row)
            call elsi_get_global_row(e_h,g_row2,i_row2)

            if(g_row2 > e_h%n_basis) then
               g_row2 = e_h%n_basis
               i_row2 = i_row+g_row2-g_row
            endif

            e_h%evec_real(i_row:i_row2,i_col) = tmp_real(g_row:g_row2)
         enddo
      endif
   enddo

   call elsi_deallocate(e_h,tmp_real,"tmp_real")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished matrix redistribution')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine sets default SIPs parameters.
!!
subroutine elsi_set_sips_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*40, parameter :: caller = "elsi_set_sips_default"

   ! How many steps of ELPA to run before SIPs
   e_h%sips_n_elpa = 0

   ! Type of slices
   ! 0 = Equally spaced subintervals
   ! 2 = Equally populated subintervals
   ! 3 = K-means after equally populated
   e_h%slicing_method = 2

   ! Extra inertia computations before solve?
   ! 0 = No
   ! 1 = Yes
   e_h%inertia_option = 0

   ! How to bound the left side of the interval
   ! 0 = Bounded
   ! 1 = -infinity
   e_h%unbound = 0

   ! Small buffer to expand the eigenvalue interval
   ! Smaller values improve performance if eigenvalue range known
   e_h%slice_buffer = 0.5_r8

   ! Lower bound of eigenvalue
   e_h%ev_min = -1.0e1_r8

   ! Upper bound of eigenvalue
   e_h%ev_max = 1.0e1_r8

end subroutine

end module ELSI_SIPS
