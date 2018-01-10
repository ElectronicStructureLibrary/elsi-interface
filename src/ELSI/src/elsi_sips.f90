! Copyright (c) 2015-2018, the ELSI team. All rights reserved.
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
   use ELSI_DATATYPE,  only: elsi_handle
   use ELSI_IO,        only: elsi_say
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,       only: elsi_check_mpi
   use ELSI_PRECISION, only: r8,i4
   use ELSI_TIMINGS,   only: elsi_get_time
   use ELSI_UTILS,     only: elsi_get_global_row
   use M_QETSC,        only: sips_initialize,sips_finalize,sips_load_ham_ovlp,&
                             sips_load_ham,sips_update_ham,sips_set_eps,&
                             sips_update_eps,sips_set_slices,sips_solve_eps,&
                             sips_get_eigenvalues,sips_get_eigenvectors

   implicit none

   private

   public :: elsi_set_sips_default
   public :: elsi_init_sips
   public :: elsi_solve_evp_sips_real
   public :: elsi_sips_to_blacs_ev_real

contains

!>
!! This routine initializes SIPs.
!! This does not change the state of the handle.
!!
subroutine elsi_init_sips(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   character*200 :: info_str

   character*40, parameter :: caller = "elsi_init_sips"

   if(e_h%n_elsi_calls == e_h%sips_n_elpa+1) then
      call sips_initialize()

      if(e_h%sips_n_slices == UNSET) then
         ! TODO: Number of slices
         e_h%sips_np_per_slice = 1
         e_h%sips_n_slices     = e_h%n_procs
      endif

      ! 1D block distribution
      e_h%n_lcol_sp = e_h%n_basis/e_h%n_procs

      ! The last process holds all remaining columns
      if(e_h%myid == e_h%n_procs-1) then
         e_h%n_lcol_sp = e_h%n_basis-(e_h%n_procs-1)*e_h%n_lcol_sp
      endif

      e_h%sips_started = .true.
   endif

   write(info_str,"('  | Number of slices     ',I10)") e_h%sips_n_slices
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine interfaces to SIPs via QETSC.
!!
subroutine elsi_solve_evp_sips_real(e_h,ham,ovlp,eval)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: ham(e_h%nnz_l_sp)
   real(kind=r8),     intent(inout) :: ovlp(e_h%nnz_l_sp)
   real(kind=r8),     intent(inout) :: eval(e_h%n_states)

   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   real(kind=r8)    :: lower
   real(kind=r8)    :: upper
   real(kind=r8)    :: max_diff
   integer(kind=i4) :: i
   integer(kind=i4) :: n_solved
   integer(kind=i4) :: n_grp
   integer(kind=i4) :: grp_ptr(100)
   integer(kind=i4) :: this
   integer(kind=i4) :: tmp
   integer(kind=i4) :: mpierr
   character*200    :: info_str

   real(kind=r8),    allocatable :: slices(:)
   real(kind=r8),    allocatable :: shifts(:)
   real(kind=r8),    allocatable :: eval_save(:)
   integer(kind=i4), allocatable :: inertias(:)

   character*40, parameter :: caller = "elsi_solve_evp_sips_real"

   ! Solve the eigenvalue problem
   call elsi_say(e_h,"  Starting SIPs eigensolver")

   call elsi_get_time(e_h,t0)

   if(e_h%n_elsi_calls == e_h%sips_n_elpa+1) then
      if(.not. e_h%ovlp_is_unit) then
         ! Load H and S
         call sips_load_ham_ovlp(e_h%n_basis,e_h%n_lcol_sp,e_h%nnz_l_sp,&
                 e_h%row_ind_sips,e_h%col_ptr_sips,ham,ovlp)

         call sips_set_eps(0)
      else
         ! Load H
         call sips_load_ham(e_h%n_basis,e_h%n_lcol_sp,e_h%nnz_l_sp,&
                 e_h%row_ind_sips,e_h%col_ptr_sips,ham)

         call sips_set_eps(1)
      endif
   else ! n_elsi_calls > sips_n_elpa+1
      ! Update H matrix
      call sips_update_ham(e_h%n_basis,e_h%n_lcol_sp,e_h%nnz_l_sp,&
              e_h%row_ind_sips,e_h%col_ptr_sips,ham)

      call sips_update_eps(e_h%sips_n_slices)
   endif

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished loading matrices')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

   call elsi_allocate(e_h,slices,e_h%sips_n_slices+1,"slices",caller)

   eval = eval+e_h%sips_ev_shift

   n_grp      = 1
   grp_ptr    = 0
   grp_ptr(1) = 1

   do i = 1,e_h%n_states-1
      if(eval(i+1)-eval(i) > 1.0_r8) then
         n_grp          = n_grp+1
         grp_ptr(n_grp) = i+1
      endif
   enddo

   if(e_h%sips_n_slices-n_grp+1 < n_grp) then
      lower = eval(1)-e_h%sips_buffer
      upper = eval(e_h%n_states)+e_h%sips_buffer

      call elsi_linspace(lower,upper,e_h%sips_n_slices+1,slices)
   else
      tmp = 0

      do i = 1,n_grp-1
         lower = eval(grp_ptr(i))-e_h%sips_buffer
         upper = eval(grp_ptr(i+1)-1)+e_h%sips_buffer
         this  = (grp_ptr(i+1)-grp_ptr(i))*(e_h%sips_n_slices-n_grp+1)&
                    /e_h%n_states
         this  = max(1,this)

         call elsi_linspace(lower,upper,this+1,slices(tmp+1:tmp+this+1))

         tmp = tmp+this+1
      enddo

      lower = eval(grp_ptr(n_grp))-e_h%sips_buffer
      upper = eval(e_h%n_states)+e_h%sips_buffer
      this  = e_h%sips_n_slices-tmp

      call elsi_linspace(lower,upper,this+1,slices(tmp+1:e_h%sips_n_slices+1))

      if(e_h%sips_do_inertia) then
         call elsi_allocate(e_h,shifts,e_h%sips_n_slices+1,"shifts",caller)
         call elsi_allocate(e_h,inertias,e_h%sips_n_slices+1,"inertias",caller)

!         call sips_get_inertia(e_h%sips_n_slices,slices,inertias)

         call elsi_deallocate(e_h,shifts,"shifts")
         call elsi_deallocate(e_h,inertias,"inertias")
      endif
   endif

   call sips_set_slices(e_h%sips_n_slices,slices)

   call elsi_get_time(e_h,t0)

   ! Solve
   call sips_solve_eps(e_h%n_states,n_solved)

   ! Get eigenvalues
   call elsi_allocate(e_h,eval_save,e_h%n_states,"eval_save",caller)
   eval_save = eval

   call sips_get_eigenvalues(e_h%n_states,eval(1:e_h%n_states))

   max_diff = maxval(abs(eval-eval_save))

   if(max_diff < e_h%sips_inertia_tol) then
      e_h%sips_do_inertia = .false.
   endif

! DEBUG
if(e_h%myid == 0) then
   print *,"VY: slices"
   do i = 1,e_h%sips_n_slices+1
      print *,slices(i)
   enddo
   print *

   print *,"VY: eigenvalues"
   do i = 1,n_grp-1
      print *,eval(grp_ptr(i))
      print *,eval(grp_ptr(i+1)-1)
   enddo
   print *,eval(e_h%n_states)
   print *
endif

   call MPI_Barrier(e_h%mpi_comm,mpierr)

   call elsi_check_mpi(e_h,"MPI_Barrier",mpierr,caller)

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
subroutine elsi_sips_to_blacs_ev_real(e_h,evec)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(out)   :: evec(e_h%n_lrow,e_h%n_lcol)

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

   character*40, parameter :: caller = "elsi_sips_to_blacs_ev_real"

   call elsi_get_time(e_h,t0)

   call elsi_allocate(e_h,tmp_real,e_h%n_basis,"tmp_real",caller)

   evec = 0.0_r8

   do i_state = 1,e_h%n_states
      call sips_get_eigenvectors(e_h%n_basis,i_state,tmp_real)

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

            evec(i_row:i_row2,i_col) = tmp_real(g_row:g_row2)
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
!! This routine gets a linearly spaced vector.
!!
subroutine elsi_linspace(lower,upper,length,vec)

   implicit none

   real(kind=r8),    intent(in)  :: lower
   real(kind=r8),    intent(in)  :: upper
   integer(kind=i4), intent(in)  :: length
   real(kind=r8),    intent(out) :: vec(length)

   character*40, parameter :: caller = "elsi_linspace"

   real(kind=r8)    :: step
   integer(kind=i4) :: i

   step = (upper-lower)/(length-1)

   do i = 1,length
      vec(i) = lower+(i-1)*step
   enddo

end subroutine

!>
!! This routine sets default SIPs parameters.
!!
subroutine elsi_set_sips_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   character*40, parameter :: caller = "elsi_set_sips_default"

   if(e_h%handle_ready) then
      e_h%handle_changed = .true.
   endif

   ! How many steps of ELPA to run before SIPs
   e_h%sips_n_elpa = 1

   ! Buffer to adjust global interval
   e_h%sips_buffer = 0.02_r8

   ! Shift of eigenspectrum between SCF steps
   e_h%sips_ev_shift = 0.0_r8

   ! Tolerance to stop inertia counting
   e_h%sips_inertia_tol = 1.0e-4_r8

   ! Do inertia counting
   e_h%sips_do_inertia = .true.

   ! Slice type (not used)
   e_h%sips_slice_type = 0

end subroutine

end module ELSI_SIPS
