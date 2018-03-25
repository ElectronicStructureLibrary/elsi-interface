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
   use ELSI_DATATYPE,  only: elsi_handle
   use ELSI_IO,        only: elsi_say
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,       only: elsi_stop,elsi_check_mpi
   use ELSI_PRECISION, only: r8,i4
   use ELSI_TIMINGS,   only: elsi_get_time
   use M_SIPS,         only: sips_initialize,sips_load_ham_ovlp,sips_load_ham,&
                             sips_update_ham,sips_set_eps,sips_update_eps,&
                             sips_set_slices,sips_solve_eps,sips_get_inertias,&
                             sips_get_eigenvalues,sips_get_eigenvectors,&
                             sips_get_slices,sips_get_slices_from_inertias,&
                             sips_get_dm,sips_get_edm

   implicit none

   private

   public :: elsi_set_sips_default
   public :: elsi_init_sips
   public :: elsi_solve_evp_sips_real
   public :: elsi_compute_dm_sips_real
   public :: elsi_compute_edm_sips_real

contains

!>
!! This routine initializes SIPS.
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
      e_h%n_lcol_sp1 = e_h%n_basis/e_h%n_procs

      ! The last process holds all remaining columns
      if(e_h%myid == e_h%n_procs-1) then
         e_h%n_lcol_sp1 = e_h%n_basis-(e_h%n_procs-1)*e_h%n_lcol_sp1
      endif

      call elsi_allocate(e_h,e_h%evec_real_sips,e_h%n_lcol_sp1,e_h%n_states,&
              "evec_real_sips",caller)

      if(e_h%n_lcol_sp == UNSET) then
         e_h%n_lcol_sp = e_h%n_lcol_sp1
      endif

      e_h%sips_started = .true.
   endif

   write(info_str,"('  | Number of slices     ',I10)") e_h%sips_n_slices
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine interfaces to SIPS.
!!
subroutine elsi_solve_evp_sips_real(e_h,ham,ovlp,eval)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: ham(e_h%nnz_l_sp1)
   real(kind=r8),     intent(inout) :: ovlp(e_h%nnz_l_sp1)
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
   call elsi_say(e_h,"  Starting SIPS eigensolver")

   call elsi_get_time(t0)

   if(e_h%n_elsi_calls == e_h%sips_n_elpa+1) then
      if(.not. e_h%ovlp_is_unit) then
         ! Load H and S
         call sips_load_ham_ovlp(e_h%n_basis,e_h%n_lcol_sp1,e_h%nnz_l_sp1,&
                 e_h%row_ind_pexsi,e_h%col_ptr_pexsi,ham,ovlp)

         call sips_set_eps(0)
      else
         ! Load H
         call sips_load_ham(e_h%n_basis,e_h%n_lcol_sp1,e_h%nnz_l_sp1,&
                 e_h%row_ind_pexsi,e_h%col_ptr_pexsi,ham)

         call sips_set_eps(1)
      endif
   else ! n_elsi_calls > sips_n_elpa+1
      ! Update H matrix
      call sips_update_ham(e_h%n_basis,e_h%n_lcol_sp1,e_h%nnz_l_sp1,&
              e_h%row_ind_pexsi,e_h%col_ptr_pexsi,ham)

      call sips_update_eps(e_h%sips_n_slices)
   endif

   call elsi_get_time(t1)

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
      call elsi_get_time(t0)

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
            do i = 1,e_h%sips_n_slices
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
            do i = e_h%sips_n_slices+1,2,-1
               if(inertias(i-1) < e_h%n_states+e_h%sips_first_ev-1) then
                  eval(e_h%n_states) = slices(i)

                  exit
               endif
            enddo
         endif
      enddo

      call sips_get_slices_from_inertias(e_h%n_states,e_h%sips_n_slices,&
              inertias,slices)

      call elsi_deallocate(e_h,inertias,"inertias")

      call elsi_get_time(t1)

      write(info_str,"('  Finished inertia counting')")
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_say(e_h,info_str)
   else
      call sips_get_slices(e_h%sips_slice_type,e_h%n_states,e_h%sips_n_slices,&
              e_h%sips_inertia_tol*2,1.0e-5_r8,eval,slices)
   endif

   call sips_set_slices(e_h%sips_n_slices,slices)

   call elsi_get_time(t0)

   ! Solve
   call sips_solve_eps(e_h%n_states,n_solved)

   if(n_solved < e_h%n_states) then
      call elsi_stop(" SIPS solver failed.",e_h,caller)
   endif

   call elsi_get_time(t1)

   write(info_str,"('  Finished solving generalized eigenproblem')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

   call elsi_get_time(t0)

   ! Get solutions
   call elsi_allocate(e_h,eval_save,e_h%n_states,"eval_save",caller)
   eval_save = eval

   call sips_get_eigenvalues(e_h%n_states,eval(1:e_h%n_states))
   call sips_get_eigenvectors(e_h%n_states,e_h%n_lcol_sp1,e_h%evec_real_sips)

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

   call elsi_get_time(t1)

   write(info_str,"('  Finished retrieving eigensolutions')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine constructs the density matrix.
!!
subroutine elsi_compute_dm_sips_real(e_h,dm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: dm(e_h%nnz_l_sp)

   character(len=40), parameter :: caller = "elsi_compute_dm_sips_real"

   call sips_get_dm(e_h%n_lcol_sp1,e_h%nnz_l_sp1,e_h%row_ind_pexsi,&
           e_h%col_ptr_pexsi,e_h%n_states,e_h%occ_num,dm)

end subroutine

!>
!! This routine constructs the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_sips_real(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: edm(e_h%nnz_l_sp)

   character(len=40), parameter :: caller = "elsi_compute_edm_sips_real"

   call sips_get_edm(e_h%n_lcol_sp1,e_h%nnz_l_sp1,e_h%row_ind_pexsi,&
           e_h%col_ptr_pexsi,e_h%n_states,e_h%occ_num,edm)

end subroutine

!>
!! This routine sets default SIPS parameters.
!!
subroutine elsi_set_sips_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   character(len=40), parameter :: caller = "elsi_set_sips_default"

   if(e_h%handle_ready) then
      e_h%handle_changed = .true.
   endif

   ! How many steps of ELPA to run before SIPS
   e_h%sips_n_elpa = 0

   ! Buffer to adjust interval
   e_h%sips_buffer = 0.01_r8

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
