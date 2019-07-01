! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Interface to SIPS.
!!
module ELSI_SIPS

   use ELSI_CONSTANT, only: UNSET
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI, only: elsi_stop
   use ELSI_OUTPUT, only: elsi_say,elsi_get_time
   use ELSI_PRECISION, only: r8,i4
   use M_SIPS, only: sips_initialize,sips_load_ham_ovlp,sips_load_ham,&
       sips_update_ham,sips_set_eps,sips_update_eps,sips_set_slices,&
       sips_get_slices,sips_get_inertias,sips_get_slices_from_inertias,&
       sips_solve_eps,sips_get_eigenvalues,sips_get_eigenvectors,sips_get_dm,&
       sips_get_edm,sips_finalize

   implicit none

   private

   public :: elsi_init_sips
   public :: elsi_cleanup_sips
   public :: elsi_solve_sips
   public :: elsi_build_dm_sips
   public :: elsi_build_edm_sips

   interface elsi_solve_sips
      module procedure elsi_solve_sips_real
   end interface

   interface elsi_build_dm_sips
      module procedure elsi_build_dm_sips_real
   end interface

   interface elsi_build_edm_sips
      module procedure elsi_build_edm_sips_real
   end interface

contains

!>
!! Initialize SLEPc-SIPs.
!!
subroutine elsi_init_sips(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh

   character(len=*), parameter :: caller = "elsi_init_sips"

   if(.not. ph%sips_started) then
      call sips_initialize()

      if(ph%sips_n_slices == UNSET) then
         ! TODO: Number of slices
         ph%sips_n_slices = bh%n_procs
      end if

      ! 1D block distribution
      bh%n_lcol_sp1 = ph%n_basis/bh%n_procs

      ! The last process holds all remaining columns
      if(bh%myid == bh%n_procs-1) then
         bh%n_lcol_sp1 = ph%n_basis-(bh%n_procs-1)*bh%n_lcol_sp1
      end if

      if(bh%n_lcol_sp == UNSET) then
         bh%n_lcol_sp = bh%n_lcol_sp1
      end if

      ph%sips_started = .true.
   end if

end subroutine

!>
!! Interface to SLEPc-SIPs.
!!
subroutine elsi_solve_sips_real(ph,bh,row_ind,col_ptr,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(inout) :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8), intent(in) :: ham(bh%nnz_l_sp1)
   real(kind=r8), intent(in) :: ovlp(bh%nnz_l_sp1)
   real(kind=r8), intent(inout) :: eval(ph%n_states)
   real(kind=r8), intent(out) :: evec(bh%n_lcol_sp1,ph%n_states)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   real(kind=r8) :: max_diff
   integer(kind=i4) :: i
   integer(kind=i4) :: n_solved
   logical :: inertia_ok
   character(len=200) :: msg

   real(kind=r8), allocatable :: eval_save(:)
   real(kind=r8), allocatable :: slices(:)
   integer(kind=i4), allocatable :: inertias(:)

   character(len=*), parameter :: caller = "elsi_solve_sips_real"

   write(msg,"(A)") "Starting SLEPc-SIPs eigensolver"
   call elsi_say(bh,msg)

   call elsi_get_time(t0)

   if(ph%sips_first) then
      if(.not. ph%unit_ovlp) then
         ! Load H and S
         call sips_load_ham_ovlp(ph%n_basis,bh%n_lcol_sp1,bh%nnz_l_sp1,row_ind,&
              col_ptr,ham,ovlp)

         call sips_set_eps(0)
      else
         ! Load H
         call sips_load_ham(ph%n_basis,bh%n_lcol_sp1,bh%nnz_l_sp1,row_ind,&
              col_ptr,ham)

         call sips_set_eps(1)
      end if
   else
      ! Update H matrix
      call sips_update_ham(ph%n_basis,bh%n_lcol_sp1,bh%nnz_l_sp1,row_ind,&
           col_ptr,ham)

      call sips_update_eps(ph%sips_n_slices)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished loading matrices"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   call elsi_allocate(bh,slices,ph%sips_n_slices+1,"slices",caller)

   if(ph%sips_first .and. ph%elpa_first) then
      eval = ph%sips_interval(2)
      eval(1) = ph%sips_interval(1)
   end if

   if(ph%sips_do_inertia) then
      call elsi_get_time(t0)

      call elsi_allocate(bh,inertias,ph%sips_n_slices+1,"inertias",caller)

      inertia_ok = .false.
      eval(1) = eval(1)-ph%sips_buffer
      eval(ph%n_states) = eval(ph%n_states)+ph%sips_buffer

      do while(.not. inertia_ok)
         inertia_ok = .true.

         call sips_get_slices(0,ph%n_states,ph%sips_n_slices,0.0_r8,1.0e-5_r8,&
              eval,slices)

         call sips_get_inertias(ph%sips_n_slices,slices,inertias)

         if(inertias(1) > 0) then
            eval(1) = eval(1)-ph%sips_buffer
            inertia_ok = .false.
         else
            do i = 1,ph%sips_n_slices
               if(inertias(i+1) > 0) then
                  eval(1) = slices(i)

                  if(i > 1) then
                     inertia_ok = .false.
                  end if

                  exit
               end if
            end do
         end if

         if(inertias(ph%sips_n_slices+1) < ph%n_states) then
            eval(ph%n_states) = eval(ph%n_states)+ph%sips_buffer
            inertia_ok = .false.
         else
            do i = ph%sips_n_slices+1,2,-1
               if(inertias(i-1) < ph%n_states) then
                  eval(ph%n_states) = slices(i)

                  if(i < ph%sips_n_slices+1) then
                     inertia_ok = .false.
                  end if

                  exit
               end if
            end do
         end if
      end do

      call sips_get_slices_from_inertias(ph%n_states,ph%sips_n_slices,inertias,&
           slices)

      call elsi_deallocate(bh,inertias,"inertias")

      call elsi_get_time(t1)

      write(msg,"(A)") "Finished inertia counting"
      call elsi_say(bh,msg)
      write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,msg)
   else
      call sips_get_slices(ph%sips_slice_type,ph%n_states,ph%sips_n_slices,&
           ph%sips_inertia_tol*2,1.0e-5_r8,eval,slices)
   end if

   call sips_set_slices(ph%sips_n_slices,slices)

   call elsi_get_time(t0)

   ! Solve
   call sips_solve_eps(n_solved)

   if(n_solved < ph%n_states) then
      write(msg,"(A)") "SLEPc-SIPs solver failed"
      call elsi_stop(bh,msg,caller)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished solving generalized eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   call elsi_get_time(t0)

   ! Get solutions
   call elsi_allocate(bh,eval_save,ph%n_states,"eval_save",caller)

   eval_save = eval

   call sips_get_eigenvalues(ph%n_states,eval(1:ph%n_states))
   call sips_get_eigenvectors(ph%n_states,bh%n_lcol_sp1,evec)

   max_diff = maxval(abs(eval_save-eval))

   if(max_diff > ph%sips_inertia_tol) then
      ph%sips_do_inertia = .true.
   else
      ph%sips_do_inertia = .false.
   end if

   call elsi_deallocate(bh,eval_save,"eval_save")
   call elsi_deallocate(bh,slices,"slices")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished retrieving eigensolutions"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%sips_first = .false.

end subroutine

!>
!! Construct the density matrix.
!!
subroutine elsi_build_dm_sips_real(ph,bh,row_ind,col_ptr,occ,dm)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(inout) :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8), intent(in) :: occ(ph%n_states)
   real(kind=r8), intent(out) :: dm(bh%nnz_l_sp1)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_build_dm_sips_real"

   call elsi_get_time(t0)

   call sips_get_dm(bh%n_lcol_sp1,bh%nnz_l_sp1,row_ind,col_ptr,ph%n_states,occ,&
        dm)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished density matrix calculation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Construct the energy-weighted density matrix.
!!
subroutine elsi_build_edm_sips_real(ph,bh,row_ind,col_ptr,occ,edm)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(inout) :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4), intent(inout) :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8), intent(in) :: occ(ph%n_states)
   real(kind=r8), intent(out) :: edm(bh%nnz_l_sp1)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_build_edm_sips_real"

   call elsi_get_time(t0)

   call sips_get_edm(bh%n_lcol_sp1,bh%nnz_l_sp1,row_ind,col_ptr,ph%n_states,&
        occ,edm)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished energy density matrix calculation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Clean up SLEPc-SIPs.
!!
subroutine elsi_cleanup_sips(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   character(len=*), parameter :: caller = "elsi_cleanup_sips"

   if(ph%sips_started) then
      call sips_finalize()
   end if

   ph%sips_first = .true.
   ph%sips_started = .false.

end subroutine

end module ELSI_SIPS
