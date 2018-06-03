! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to SIPS.
!!
module ELSI_SIPS

   use ELSI_CONSTANTS, only: UNSET
   use ELSI_DATATYPE,  only: elsi_param_t,elsi_basic_t
   use ELSI_IO,        only: elsi_say,elsi_get_time
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,       only: elsi_stop
   use ELSI_PRECISION, only: r8,i4
   use M_SIPS,         only: sips_initialize,sips_load_ham_ovlp,sips_load_ham,&
                             sips_update_ham,sips_set_eps,sips_update_eps,&
                             sips_set_slices,sips_solve_eps,sips_get_inertias,&
                             sips_get_eigenvalues,sips_get_eigenvectors,&
                             sips_get_slices,sips_get_slices_from_inertias,&
                             sips_get_dm,sips_get_edm,sips_finalize

   implicit none

   private

   public :: elsi_init_sips
   public :: elsi_cleanup_sips
   public :: elsi_solve_sips
   public :: elsi_compute_dm_sips
   public :: elsi_compute_edm_sips

   interface elsi_solve_sips
      module procedure elsi_solve_sips_real
   end interface

   interface elsi_compute_dm_sips
      module procedure elsi_compute_dm_sips_real
   end interface

   interface elsi_compute_edm_sips
      module procedure elsi_compute_edm_sips_real
   end interface

contains

!>
!! This routine initializes SIPS.
!!
subroutine elsi_init_sips(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh

   character(len=40), parameter :: caller = "elsi_init_sips"

   if(ph%n_calls == ph%sips_n_elpa+1) then
      call sips_initialize()

      if(ph%sips_n_slices == UNSET) then
         ! TODO: Number of slices
         ph%sips_np_per_slice = 1
         ph%sips_n_slices     = bh%n_procs
      endif

      ! 1D block distribution
      bh%n_lcol_sp1 = ph%n_basis/bh%n_procs

      ! The last process holds all remaining columns
      if(bh%myid == bh%n_procs-1) then
         bh%n_lcol_sp1 = ph%n_basis-(bh%n_procs-1)*bh%n_lcol_sp1
      endif

      if(bh%n_lcol_sp == UNSET) then
         bh%n_lcol_sp = bh%n_lcol_sp1
      endif

      ph%sips_started = .true.
   endif

end subroutine

!>
!! This routine interfaces to SIPS.
!!
subroutine elsi_solve_sips_real(ph,bh,row_ind,col_ptr,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8),      intent(in)    :: ham(bh%nnz_l_sp1)
   real(kind=r8),      intent(in)    :: ovlp(bh%nnz_l_sp1)
   real(kind=r8),      intent(inout) :: eval(ph%n_states)
   real(kind=r8),      intent(out)   :: evec(bh%n_lcol_sp1,ph%n_states)

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

   character(len=40), parameter :: caller = "elsi_solve_sips_real"

   ! Solve the eigenvalue problem
   write(info_str,"(2X,A)") "Starting SLEPc-SIPs eigensolver"
   call elsi_say(bh,info_str)

   call elsi_get_time(t0)

   if(ph%n_calls == ph%sips_n_elpa+1) then
      if(.not. ph%ovlp_is_unit) then
         ! Load H and S
         call sips_load_ham_ovlp(ph%n_basis,bh%n_lcol_sp1,bh%nnz_l_sp1,row_ind,&
                 col_ptr,ham,ovlp)

         call sips_set_eps(0)
      else
         ! Load H
         call sips_load_ham(ph%n_basis,bh%n_lcol_sp1,bh%nnz_l_sp1,row_ind,&
                 col_ptr,ham)

         call sips_set_eps(1)
      endif
   else ! n_calls > sips_n_elpa+1
      ! Update H matrix
      call sips_update_ham(ph%n_basis,bh%n_lcol_sp1,bh%nnz_l_sp1,row_ind,&
              col_ptr,ham)

      call sips_update_eps(ph%sips_n_slices)
   endif

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished loading matrices"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

   call elsi_allocate(bh,slices,ph%sips_n_slices+1,"slices",caller)

   if(ph%n_calls == 1) then
      eval    = ph%sips_interval(2)
      eval(1) = ph%sips_interval(1)
   endif

   if(ph%sips_do_inertia) then ! Inertia counting
      call elsi_get_time(t0)

      call elsi_allocate(bh,inertias,ph%sips_n_slices+1,"inertias",caller)

      inertia_ok        = .false.
      eval(1)           = eval(1)-ph%sips_buffer
      eval(ph%n_states) = eval(ph%n_states)+ph%sips_buffer

      do while(.not. inertia_ok)
         inertia_ok = .true.

         call sips_get_slices(0,ph%n_states,ph%sips_n_slices,0.0_r8,1.0e-5_r8,&
                 eval,slices)

         call sips_get_inertias(ph%sips_n_slices,slices,inertias)

         if(inertias(1) > ph%sips_first_ev-1) then
            eval(1)    = eval(1)-ph%sips_buffer
            inertia_ok = .false.
         else
            do i = 1,ph%sips_n_slices
               if(inertias(i+1) > ph%sips_first_ev-1) then
                  eval(1) = slices(i)

                  exit
               endif
            enddo

            if(inertias(i) < ph%sips_first_ev-1) then
               eval(1)    = eval(1)+ph%sips_buffer
               inertia_ok = .false.
            endif
         endif

         if(inertias(ph%sips_n_slices+1) < &
            ph%n_states+ph%sips_first_ev-1) then
            eval(ph%n_states) = eval(ph%n_states)+ph%sips_buffer
            inertia_ok        = .false.
         else
            do i = ph%sips_n_slices+1,2,-1
               if(inertias(i-1) < ph%n_states+ph%sips_first_ev-1) then
                  eval(ph%n_states) = slices(i)

                  exit
               endif
            enddo
         endif
      enddo

      call sips_get_slices_from_inertias(ph%n_states,ph%sips_n_slices,inertias,&
              slices)

      call elsi_deallocate(bh,inertias,"inertias")

      call elsi_get_time(t1)

      write(info_str,"(2X,A)") "Finished inertia counting"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,info_str)
   else
      call sips_get_slices(ph%sips_slice_type,ph%n_states,ph%sips_n_slices,&
              ph%sips_inertia_tol*2,1.0e-5_r8,eval,slices)
   endif

   call sips_set_slices(ph%sips_n_slices,slices)

   call elsi_get_time(t0)

   ! Solve
   call sips_solve_eps(ph%n_states,n_solved)

   if(n_solved < ph%n_states) then
      call elsi_stop(bh,"SLEPc-SIPs solver failed.",caller)
   endif

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished solving generalized eigenproblem"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

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
   endif

   call elsi_deallocate(bh,eval_save,"eval_save")
   call elsi_deallocate(bh,slices,"slices")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished retrieving eigensolutions"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine constructs the density matrix.
!!
subroutine elsi_compute_dm_sips_real(ph,bh,row_ind,col_ptr,occ,dm)

   implicit none

   type(elsi_param_t), intent(in)  :: ph
   type(elsi_basic_t), intent(in)  :: bh
   integer(kind=i4),   intent(in)  :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)  :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8),      intent(in)  :: occ(ph%n_states,ph%n_spins,ph%n_kpts)
   real(kind=r8),      intent(out) :: dm(bh%nnz_l_sp1)

   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_compute_dm_sips_real"

   call elsi_get_time(t0)

   call sips_get_dm(bh%n_lcol_sp1,bh%nnz_l_sp1,row_ind,col_ptr,ph%n_states,&
           occ(:,1,1),dm)

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished density matrix calculation"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine constructs the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_sips_real(ph,bh,row_ind,col_ptr,occ,edm)

   implicit none

   type(elsi_param_t), intent(in)  :: ph
   type(elsi_basic_t), intent(in)  :: bh
   integer(kind=i4),   intent(in)  :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)  :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8),      intent(in)  :: occ(ph%n_states,ph%n_spins,ph%n_kpts)
   real(kind=r8),      intent(out) :: edm(bh%nnz_l_sp1)

   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_compute_edm_sips_real"

   call elsi_get_time(t0)

   call sips_get_edm(bh%n_lcol_sp1,bh%nnz_l_sp1,row_ind,col_ptr,ph%n_states,&
           occ(:,1,1),edm)

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished energy density matrix calculation"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine cleans up SIPS.
!!
subroutine elsi_cleanup_sips(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   character(len=40), parameter :: caller = "elsi_cleanup_sips"

   if(ph%sips_started) then
      call sips_finalize()
   endif

   ph%sips_started = .false.

end subroutine

end module ELSI_SIPS
