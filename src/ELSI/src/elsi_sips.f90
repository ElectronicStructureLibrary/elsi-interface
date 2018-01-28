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
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi
   use ELSI_PRECISION, only: r8,i4
   use ELSI_TIMINGS,   only: elsi_get_time
   use M_QETSC,        only: sips_initialize,sips_load_ham_ovlp,sips_load_ham,&
                             sips_update_ham,sips_set_eps,sips_update_eps,&
                             sips_set_slices,sips_solve_eps,sips_get_inertias,&
                             sips_get_eigenvalues,sips_get_eigenvectors,&
                             sips_get_slices,sips_get_slices_from_inertias

   implicit none

   private

   public :: elsi_set_sips_default
   public :: elsi_init_sips
   public :: elsi_solve_evp_sips_real

contains

!>
!! This routine initializes SIPs.
!! This does not change the state of the handle.
!!
subroutine elsi_init_sips(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_init_sips"

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

      call elsi_allocate(e_h,e_h%evec_real_sips,e_h%n_lcol_sp,e_h%n_states,&
              "evec_real_sips",caller)

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

   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   real(kind=r8)      :: max_diff
   integer(kind=i4)   :: i
   integer(kind=i4)   :: n_solved
   integer(kind=i4)   :: ierr
   logical            :: inertia_ok
   character(len=200) :: info_str

   real(kind=r8),    allocatable :: eval_save(:)
   real(kind=r8),    allocatable :: slices(:)
   integer(kind=i4), allocatable :: inertias(:)

   character(len=40), parameter :: caller = "elsi_solve_evp_sips_real"

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

   if(e_h%n_elsi_calls == 1) then
      eval    = e_h%sips_interval(2)
      eval(1) = e_h%sips_interval(1)
   endif

   if(e_h%sips_do_inertia) then ! Inertia counting
      call elsi_get_time(e_h,t0)

      call elsi_allocate(e_h,inertias,e_h%sips_n_slices+1,"inertias",caller)

      inertia_ok         = .false.
      eval(1)            = eval(1)-e_h%sips_buffer
      eval(e_h%n_states) = eval(e_h%n_states)+e_h%sips_buffer

      do while(.not. inertia_ok)
         inertia_ok = .true.

         call sips_get_slices(0,e_h%n_states,e_h%sips_n_slices,0.0_r8,&
                 1.0e-5_r8,eval,slices)

         call sips_get_inertias(e_h%sips_n_slices,slices,inertias)

         if(inertias(1) > e_h%sips_first_ev-1) then
            eval(1)    = eval(1)-e_h%sips_buffer
            inertia_ok = .false.
         else
            do i = 1,e_h%sips_n_slices+1
               if(inertias(i+1) > e_h%sips_first_ev-1) then
                  eval(1) = slices(i)

                  exit
               endif
            enddo

            if(inertias(i) < e_h%sips_first_ev-1) then
               eval(1)    = eval(1)+e_h%sips_buffer
               inertia_ok = .false.
            endif
         endif

         if(inertias(e_h%sips_n_slices+1) < &
            e_h%n_states+e_h%sips_first_ev-1) then
            eval(e_h%n_states) = eval(e_h%n_states)+e_h%sips_buffer
            inertia_ok         = .false.
         else
            do i = e_h%sips_n_slices+1,1,-1
               if(inertias(i-1) < e_h%n_states+e_h%sips_first_ev-1) then
                  eval(e_h%n_states) = slices(i)

                  exit
               endif
            enddo
         endif
      enddo

      call sips_get_slices_from_inertias(e_h%n_states,e_h%sips_n_slices,&
              inertias,slices)

      ! DEBUG
      call sips_get_inertias(e_h%sips_n_slices,slices,inertias)

      if(e_h%myid == 0) then
         print *
         print *,"Final:"
         do i = 1,e_h%sips_n_slices+1
            print *,slices(i),":",inertias(i)
         enddo
         print *
      endif

      call elsi_deallocate(e_h,inertias,"inertias")

      call elsi_get_time(e_h,t1)

      write(info_str,"('  Finished inertia counting')")
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_say(e_h,info_str)
   else
      call sips_get_slices(e_h%sips_slice_type,e_h%n_states,e_h%sips_n_slices,&
              e_h%sips_inertia_tol*2,1.0e-5_r8,eval,slices)
   endif

   call sips_set_slices(e_h%sips_n_slices,slices)

   call elsi_get_time(e_h,t0)

   ! Solve
   call sips_solve_eps(e_h%n_states,n_solved)

   if(n_solved < e_h%n_states) then
      call elsi_stop(" SIPs solver failed.",e_h,caller)
   endif

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished solving generalized eigenproblem')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

   call elsi_get_time(e_h,t0)

   ! Get solutions
   call elsi_allocate(e_h,eval_save,e_h%n_states,"eval_save",caller)
   eval_save = eval

   call sips_get_eigenvalues(e_h%n_states,eval(1:e_h%n_states))
   call sips_get_eigenvectors(e_h%n_states,e_h%n_lcol_sp,e_h%evec_real_sips)

   max_diff = maxval(abs(eval_save-eval))

   if(max_diff > e_h%sips_inertia_tol) then
      e_h%sips_do_inertia = .true.
   else
      e_h%sips_do_inertia = .false.
   endif

   call elsi_deallocate(e_h,eval_save,"eval_save")
   call elsi_deallocate(e_h,slices,"slices")

   call MPI_Barrier(e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Barrier",ierr,caller)

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished retrieving eigensolutions')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine sets default SIPs parameters.
!!
subroutine elsi_set_sips_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   character(len=40), parameter :: caller = "elsi_set_sips_default"

   if(e_h%handle_ready) then
      e_h%handle_changed = .true.
   endif

   ! How many steps of ELPA to run before SIPs
   e_h%sips_n_elpa = 1

   ! Buffer to adjust interval
   e_h%sips_buffer = 0.02_r8

   ! Do inertia counting
   e_h%sips_do_inertia = .true.

   ! Criterion to stop inertia counting
   e_h%sips_inertia_tol = 0.001_r8

   ! Initial global interval
   e_h%sips_interval(1) = -2.0_r8
   e_h%sips_interval(2) = 2.0_r8

   ! Slice type
   e_h%sips_slice_type = 2

   ! Index of 1st eigensolution to be solved
   e_h%sips_first_ev = 1

end subroutine

end module ELSI_SIPS
