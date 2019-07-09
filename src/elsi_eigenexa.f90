! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Interface to EigenExa.
!!
module ELSI_EIGENEXA

   use ELSI_CONSTANT, only: UNSET
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_ELPA, only: elsi_factor_ovlp_elpa,elsi_reduce_evp_elpa,&
       elsi_back_ev_elpa
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI, only: elsi_check_mpi,mpi_sum,mpi_integer4
   use ELSI_OUTPUT, only: elsi_say,elsi_get_time
   use ELSI_PRECISION, only: r8,i4
!   use ELSI_REDIST, only: elsi_blacs_to_eigenexa_h,elsi_eigenexa_to_blacs_ev
   use ELSI_UTIL, only: elsi_get_nnz
   use EIGEN_LIBS_MOD, only: eigen_init,eigen_get_procs,eigen_get_id,&
       eigen_get_matdims,eigen_s,eigen_sx,eigen_free

   implicit none

   private

   public :: elsi_init_eigenexa
   public :: elsi_cleanup_eigenexa
   public :: elsi_solve_eigenexa

   interface elsi_solve_eigenexa
      module procedure elsi_solve_eigenexa_real
   end interface

contains

!>
!! Initialize EigenExa.
!!
subroutine elsi_init_eigenexa(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh

   integer(kind=i4) :: n_procs
   integer(kind=i4) :: myid

   character(len=*), parameter :: caller = "elsi_init_eigenexa"

   if(.not. ph%exa_started) then
      call eigen_init(bh%comm)
      call eigen_get_procs(n_procs,ph%exa_n_prow,ph%exa_n_pcol)
      call eigen_get_id(myid,ph%exa_my_prow,ph%exa_my_pcol)
      call eigen_get_matdims(ph%n_basis,ph%exa_n_lrow,ph%exa_n_lcol)

      ph%exa_started = .true.
   end if

end subroutine

!>
!! Interface to EigenExa.
!!
subroutine elsi_solve_eigenexa_real(ph,bh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   real(kind=r8), intent(inout) :: ham(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(inout) :: ovlp(bh%n_lrow,bh%n_lcol)
   real(kind=r8), intent(out) :: eval(ph%n_basis)
   real(kind=r8), intent(out) :: evec(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: ierr
   character :: mode
   character(len=200) :: msg

   real(kind=r8), allocatable :: ham_exa(:,:)
   real(kind=r8), allocatable :: evec_exa(:,:)

   character(len=*), parameter :: caller = "elsi_solve_eigenexa_real"

   ! Compute sparsity
   if(bh%nnz_g == UNSET) then
      if(bh%nnz_l == UNSET) then
         call elsi_get_nnz(bh%def0,bh%n_lrow,bh%n_lcol,ham,bh%nnz_l)
      end if

      call MPI_Allreduce(bh%nnz_l,bh%nnz_g,1,mpi_integer4,mpi_sum,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   end if

   ! Transform to standard form
   if(.not. ph%unit_ovlp) then
      if(ph%exa_first) then
         call elsi_factor_ovlp_elpa(ph,bh,ovlp)
      end if

      call elsi_reduce_evp_elpa(ph,bh,ham,ovlp,evec)
   end if

   ! BLACS to EigenExa
   call elsi_allocate(bh,ham_exa,ph%exa_n_lrow,ph%exa_n_lcol,"ham_exa",caller)
   call elsi_allocate(bh,evec_exa,ph%exa_n_lrow,ph%exa_n_lcol,"evec_exa",caller)

!   call elsi_blacs_to_eigenexa_h(ph,bh,ham,ham_exa)

   ! Solve
   write(msg,"(A)") "Starting EigenExa eigensolver"
   call elsi_say(bh,msg)

   call elsi_get_time(t0)

   if(ph%n_states <= 0) then
      mode = "N"
   else
      mode = "A"
   end if

   if(ph%exa_method == 1) then
      call eigen_s(ph%n_basis,ph%n_basis,ham_exa,ph%exa_n_lrow,eval,evec_exa,&
           ph%exa_n_lrow,ph%exa_blk_fwd,ph%exa_blk_bkwd,mode)
   else
      call eigen_sx(ph%n_basis,ph%n_basis,ham_exa,ph%exa_n_lrow,eval,evec_exa,&
           ph%exa_n_lrow,ph%exa_blk_fwd,ph%exa_blk_bkwd,mode)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished solving standard eigenproblem"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ! EigenExa to BLACS
!   call elsi_eigenexa_to_blacs_ev(ph,bh,evec_exa,evec)

   call elsi_deallocate(bh,ham_exa,"ham_exa")
   call elsi_deallocate(bh,evec_exa,"evec_exa")

   ! Back-transform eigenvectors
   if(.not. ph%unit_ovlp) then
      call elsi_back_ev_elpa(ph,bh,ham,ovlp,evec)
   end if

   ph%exa_first = .false.

end subroutine

!>
!! Clean up EigenExa.
!!
subroutine elsi_cleanup_eigenexa(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   character(len=*), parameter :: caller = "elsi_cleanup_eigenexa"

   if(ph%exa_started) then
      call eigen_free()
   end if

   ph%exa_first = .true.
   ph%exa_started = .false.

end subroutine

end module ELSI_EIGENEXA
