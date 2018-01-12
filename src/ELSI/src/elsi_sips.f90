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
   use ELSI_UTILS,     only: elsi_get_global_row
   use M_QETSC,        only: sips_initialize,sips_finalize,sips_load_ham_ovlp,&
                             sips_load_ham,sips_update_ham,sips_set_eps,&
                             sips_update_eps,sips_set_slices,sips_solve_eps,&
                             sips_get_eigenvalues,sips_get_eigenvectors,&
                             sips_get_inertias

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

   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   real(kind=r8)    :: lower
   real(kind=r8)    :: upper
   real(kind=r8)    :: prev_sum
   real(kind=r8)    :: new_sum
   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: this
   integer(kind=i4) :: tmp
   integer(kind=i4) :: n_iner_steps
   integer(kind=i4) :: n_solved
   integer(kind=i4) :: n_clst
   integer(kind=i4) :: clst(100,3) ! Eigenvalue clusters
   integer(kind=i4) :: ierr
   logical          :: slices_changed
   character*200    :: info_str

   real(kind=r8),    allocatable :: slices(:)
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

   n_clst    = 1
   clst      = 0
   clst(1,1) = 1

   do i = 1,e_h%n_states-1
      if(eval(i+1)-eval(i) > 1.0_r8) then
         n_clst         = n_clst+1
         clst(n_clst,1) = i+1
      endif
   enddo

   clst(n_clst+1,1) = e_h%n_states+1

   if(e_h%sips_n_slices-n_clst+1 < n_clst .or. n_clst == 1) then ! One slice
      lower = eval(1)-e_h%sips_buffer
      upper = eval(e_h%n_states)+e_h%sips_buffer

      call elsi_linspace(lower,upper,e_h%sips_n_slices+1,slices)

      if(e_h%sips_do_inertia) then ! Inertia counting
         call elsi_get_time(e_h,t0)

         call elsi_allocate(e_h,inertias,e_h%sips_n_slices+1,"inertias",caller)

         n_iner_steps   = 0
         slices_changed = .true.

         do while(n_iner_steps < 5 .and. slices_changed)
            n_iner_steps = n_iner_steps+1

            call sips_get_inertias(e_h%sips_n_slices,slices,inertias)

            slices_changed = .false.

            if(inertias(1) > 0) then
               slices(1)      = slices(1)-0.1_r8
               slices_changed = .true.
            endif
         enddo

         call elsi_deallocate(e_h,inertias,"inertias")

         lower = slices(1)
         upper = slices(e_h%sips_n_slices+1)

         call elsi_linspace(lower,upper,e_h%sips_n_slices+1,slices)

         call elsi_get_time(e_h,t1)

         write(info_str,"('  Finished inertia counting')")
         call elsi_say(e_h,info_str)
         write(info_str,"('  | Time :',F10.3,' s')") t1-t0
         call elsi_say(e_h,info_str)
      endif
   else
      tmp = 0

      do i = 1,n_clst-1
         this      = (clst(i+1,1)-clst(i,1))*(e_h%sips_n_slices-n_clst+1)/&
                        e_h%n_states
         this      = max(1,this)
         clst(i,2) = tmp+1
         clst(i,3) = tmp+this+1
         tmp       = tmp+this+1
      enddo

      clst(n_clst,2) = tmp+1
      clst(n_clst,3) = e_h%sips_n_slices+1

      do i = 1,n_clst
         lower = eval(clst(i,1))-e_h%sips_buffer
         upper = eval(clst(i+1,1)-1)+e_h%sips_buffer
         this  = clst(i,3)-clst(i,2)+1

         call elsi_linspace(lower,upper,this,slices(clst(i,2):clst(i,3)))
      enddo

      if(e_h%sips_do_inertia) then ! Inertia counting
         call elsi_get_time(e_h,t0)

         call elsi_allocate(e_h,inertias,e_h%sips_n_slices+1,"inertias",caller)

         n_iner_steps   = 0
         slices_changed = .true.

         do while(n_iner_steps < 5 .and. slices_changed)
            n_iner_steps = n_iner_steps+1

            call sips_get_inertias(e_h%sips_n_slices,slices,inertias)

            ! DEBUG
            if(e_h%myid == 0) then
               print *
               print *,"QETSc inertia counting iteration",n_iner_steps
               print *,"Shifts                  :    Inertias"
               do i = 1,e_h%sips_n_slices+1
                  print *,slices(i),":",inertias(i)
               enddo
               print *
            endif

            slices_changed = .false.

            do i = 1,n_clst
               do j = clst(i,2),clst(i,3)
                  if(inertias(j) > clst(i,1)-1) then
                     if(j == clst(i,2)) then
                        slices(clst(i,2)) = slices(clst(i,2))-0.1_r8
                        slices_changed    = .true.
                     elseif(j-1 /= clst(i,2)) then
                        slices(clst(i,2)) = slices(j-1)
                        slices_changed    = .true.
                     endif

                     exit
                  endif
               enddo

               do j = clst(i,3),clst(i,2),-1
                  if(inertias(j) < clst(i+1,1)-1) then
                     if(j == clst(i,3)) then
                        slices(clst(i,3)) = slices(clst(i,3))+0.1_r8
                        slices_changed    = .true.
                     elseif(j+1 /= clst(i,3)) then
                        slices(clst(i,3)) = slices(j+1)
                        slices_changed    = .true.
                     endif

                     exit
                  endif
               enddo

               lower = slices(clst(i,2))
               upper = slices(clst(i,3))
               this  = clst(i,3)-clst(i,2)+1

               call elsi_linspace(lower,upper,this,slices(clst(i,2):clst(i,3)))
            enddo
         enddo

         call elsi_deallocate(e_h,inertias,"inertias")

         call elsi_get_time(e_h,t1)

         write(info_str,"('  Finished inertia counting')")
         call elsi_say(e_h,info_str)
         write(info_str,"('  | Time :',F10.3,' s')") t1-t0
         call elsi_say(e_h,info_str)
      endif
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

   ! Get eigenvalues
   prev_sum = sum(eval(1:e_h%n_states))/e_h%n_states

   call sips_get_eigenvalues(e_h%n_states,eval(1:e_h%n_states))

   ! Adjust slice buffer
   new_sum         = sum(eval(1:e_h%n_states))/e_h%n_states
   e_h%sips_buffer = min(e_h%sips_buffer,abs(new_sum-prev_sum))

   ! Get eigenvectors
   call sips_get_eigenvectors(e_h%n_states,e_h%n_lcol_sp,e_h%evec_real_sips)

   call MPI_Barrier(e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Barrier",ierr,caller)

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished retrieving eigensolutions')")
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
   e_h%sips_inertia_tol = 1.0e-3_r8

   ! Do inertia counting
   e_h%sips_do_inertia = .true.

   ! Slice type (not used)
   e_h%sips_slice_type = 0

end subroutine

end module ELSI_SIPS
