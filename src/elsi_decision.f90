! Copyright (c) 2015-2020, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Decide which eigensolver or density matrix solver to use, when it is not
!! specified by the user.
!!
module ELSI_DECISION

   use ELSI_CONSTANT, only: AUTO_SOLVER,ELPA_SOLVER,PEXSI_SOLVER,NTPOLY_SOLVER,&
       UNSET
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_MPI
   use ELSI_OUTPUT, only: elsi_say
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTIL, only: elsi_check_err

   implicit none

   private

   public :: elsi_decide_ev
   public :: elsi_decide_dm

   interface elsi_decide_dm
      module procedure elsi_decide_dm_real
      module procedure elsi_decide_dm_cmplx
      module procedure elsi_decide_dm_sparse
   end interface

contains

!>
!! Decide which eigensolver to use.
!!
subroutine elsi_decide_ev(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_decide_ev"

   if(ph%solver == AUTO_SOLVER) then
      ph%solver = ELPA_SOLVER

      write(msg,"(A)") "ELPA selected"
      call elsi_say(bh,msg)
   end if

end subroutine

!>
!! Decide which density matrix solver to use.
!!
subroutine elsi_decide_dm_real(ph,bh,mat)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: mat(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: sparsity
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_decide_dm_real"

   if(ph%solver == AUTO_SOLVER) then
      if(ph%i_spin == 1 .and. ph%i_kpt == 1) then
         nnz_l = count(abs(mat) > bh%def0)

         call MPI_Allreduce(nnz_l,nnz_g,1,MPI_INTEGER4,MPI_SUM,bh%comm,ierr)

         call elsi_check_err(bh,"MPI_Allreduce",ierr,caller)

         sparsity = 1.0_r8-(1.0_r8*nnz_g/ph%n_basis/ph%n_basis)
      end if

      call MPI_Bcast(sparsity,1,MPI_REAL8,0,bh%comm_all,ierr)

      call elsi_check_err(bh,"MPI_Bcast",ierr,caller)

      call elsi_decide_dm_smart(ph,bh,sparsity)
   end if

end subroutine

!>
!! Decide which density matrix solver to use.
!!
subroutine elsi_decide_dm_cmplx(ph,bh,mat)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   complex(kind=r8), intent(in) :: mat(bh%n_lrow,bh%n_lcol)

   real(kind=r8) :: sparsity
   integer(kind=i4) :: nnz_g
   integer(kind=i4) :: nnz_l
   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_decide_dm_cmplx"

   if(ph%solver == AUTO_SOLVER) then
      if(ph%i_spin == 1 .and. ph%i_kpt == 1) then
         nnz_l = count(abs(mat) > bh%def0)

         call MPI_Allreduce(nnz_l,nnz_g,1,MPI_INTEGER4,MPI_SUM,bh%comm,ierr)

         call elsi_check_err(bh,"MPI_Allreduce",ierr,caller)

         sparsity = 1.0_r8-(1.0_r8*nnz_g/ph%n_basis/ph%n_basis)
      end if

      call MPI_Bcast(sparsity,1,MPI_REAL8,0,bh%comm_all,ierr)

      call elsi_check_err(bh,"MPI_Bcast",ierr,caller)

      call elsi_decide_dm_smart(ph,bh,sparsity)
   end if

end subroutine

!>
!! Decide which density matrix solver to use.
!!
subroutine elsi_decide_dm_sparse(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh

   real(kind=r8) :: sparsity
   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_decide_dm_sparse"

   if(ph%solver == AUTO_SOLVER) then
      sparsity = 1.0_r8-(1.0_r8*bh%nnz_g/ph%n_basis/ph%n_basis)

      call MPI_Bcast(sparsity,1,MPI_REAL8,0,bh%comm_all,ierr)

      call elsi_check_err(bh,"MPI_Bcast",ierr,caller)

      call elsi_decide_dm_smart(ph,bh,sparsity)
   end if

end subroutine

!>
!! If possible, quickly decide which density matrix solver to use.
!!
subroutine elsi_decide_dm_smart(ph,bh,sparsity)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: sparsity

   logical :: try_ntpoly
   logical :: try_pexsi
   integer(kind=i4) :: pexsi_enabled
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_decide_dm_smart"

   try_ntpoly = .true.
   try_pexsi = .true.

   ! Is PEXSI enabled?
   call elsi_get_pexsi_enabled(pexsi_enabled)

   if(pexsi_enabled == 0) then
      try_pexsi = .false.
   end if

   if(ph%n_basis < 20000 .or. sparsity < 0.95_r8) then
      try_pexsi = .false.
   end if

   if(mod(bh%n_procs,ph%pexsi_options%nPoints) /= 0) then
      try_pexsi = .false.
   end if

   if(ph%pexsi_np_per_pole /= UNSET) then
      if(mod(bh%n_procs,ph%pexsi_np_per_pole*ph%pexsi_options%nPoints) /= 0)&
         then
         try_pexsi = .false.
      end if

      if(ph%pexsi_np_per_pole*ph%pexsi_options%numPole*&
         ph%pexsi_options%nPoints < bh%n_procs) then
         try_pexsi = .false.
      end if
   end if

   if(ph%energy_gap < 0.5_r8) then
      try_ntpoly = .false.
   end if

   if(ph%n_basis < 50000 .or. sparsity < 0.98_r8) then
      try_ntpoly = .false.
   end if

   if(.not. try_pexsi .and. .not. try_ntpoly) then
      ph%solver = ELPA_SOLVER

      write(msg,"(A)") "ELPA selected"
      call elsi_say(bh,msg)
   end if

   if(try_pexsi .and. ph%dimensionality < 3) then
      ph%solver = PEXSI_SOLVER

      write(msg,"(A)") "PEXSI selected"
      call elsi_say(bh,msg)
   end if

   if(try_ntpoly .and. ph%n_basis > 100000 .and. sparsity > 0.99_r8) then
      ph%solver = NTPOLY_SOLVER

      write(msg,"(A)") "NTPoly selected"
      call elsi_say(bh,msg)
   end if

   if(ph%solver == AUTO_SOLVER) then
      ph%solver = ELPA_SOLVER
   end if

end subroutine

end module ELSI_DECISION
