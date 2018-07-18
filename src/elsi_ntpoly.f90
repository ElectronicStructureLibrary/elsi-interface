! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to NTPoly.
!!
module ELSI_NTPOLY

   use ELSI_CONSTANTS, only: UNSET,NTPOLY_TRS2,NTPOLY_TRS4,NTPOLY_HPCP
   use ELSI_DATATYPE,  only: elsi_param_t,elsi_basic_t
   use ELSI_IO,        only: elsi_say,elsi_get_time
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_PRECISION, only: r8,i4
   use NTPOLY,         only: TRS2,TRS4,HPCP,ScaleDistributedSparseMatrix,&
                             DistributedSparseMatrix_t,FillFromTripletList,&
                             ConstructEmptyDistributedSparseMatrix,&
                             FilterDistributedSparseMatrix,&
                             CopyDistributedSparseMatrix,GetTripletList,&
                             DestructDistributedSparseMatrix,InverseSquareRoot,&
                             IterativeSolverParameters_t,DestructPermutation,&
                             ConstructRandomPermutation,ConstructProcessGrid,&
                             Triplet_t,TripletList_t,ConstructTripletList,&
                             AppendToTripletList,DestructTripletList

   implicit none

   private

   public :: elsi_init_ntpoly
   public :: elsi_cleanup_ntpoly
   public :: elsi_solve_ntpoly
   public :: DistributedSparseMatrix_t
   public :: ConstructEmptyDistributedSparseMatrix
   public :: CopyDistributedSparseMatrix
   public :: FillFromTripletList
   public :: GetTripletList
   public :: FilterDistributedSparseMatrix
   public :: TripletList_t
   public :: ConstructTripletList
   public :: AppendToTripletList
   public :: DestructTripletList
   public :: Triplet_t

   interface elsi_solve_ntpoly
      module procedure elsi_solve_ntpoly_real
   end interface

contains

!>
!! This routine initializes NTPoly.
!!
subroutine elsi_init_ntpoly(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh

   integer(kind=i4) :: n_prow
   integer(kind=i4) :: n_pcol
   integer(kind=i4) :: np_per_group

   character(len=40), parameter :: caller = "elsi_init_ntpoly"

   if(.not. ph%nt_started) then
      np_per_group = bh%n_procs/ph%nt_n_group

      ! Set square-like process grid
      do n_prow = nint(sqrt(real(np_per_group,kind=r8))),2,-1
         if(mod(np_per_group,n_prow) == 0) then
            exit
         endif
      enddo

      n_pcol = np_per_group/n_prow

      call ConstructProcessGrid(bh%comm,n_prow,n_pcol,ph%nt_n_group)

      ph%nt_started = .true.
   endif

end subroutine

!>
!! This routine interfaces to NTPoly.
!!
subroutine elsi_solve_ntpoly_real(ph,bh,ham,ovlp,dm)

   implicit none

   type(elsi_param_t),              intent(inout) :: ph
   type(elsi_basic_t),              intent(in)    :: bh
   type(DistributedSparseMatrix_t), intent(in)    :: ham
   type(DistributedSparseMatrix_t), intent(inout) :: ovlp
   type(DistributedSparseMatrix_t), intent(inout) :: dm

   integer(kind=i4)   :: ne
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   type(DistributedSparseMatrix_t) :: ovlp_isq

   character(len=40), parameter :: caller = "elsi_solve_ntpoly_real"

   if(ph%n_calls == 1) then
      call elsi_get_time(t0)

      call ConstructRandomPermutation(ph%nt_perm,ovlp%logical_matrix_dimension)

      ph%nt_options = IterativeSolverParameters_t(ph%nt_tol,ph%nt_filter,&
                         ph%nt_max_iter,ph%nt_output,ph%nt_perm)

      call InverseSquareRoot(ovlp,ovlp_isq,ph%nt_options)
      call CopyDistributedSparseMatrix(ovlp_isq,ovlp)
      call DestructDistributedSparseMatrix(ovlp_isq)

      call elsi_get_time(t1)

      write(info_str,"(2X,A)") "Finished overlap matrix inversion"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,info_str)
   endif

   call elsi_get_time(t0)

   call ConstructEmptyDistributedSparseMatrix(dm,ph%n_basis)

   ne = int(ph%n_electrons,kind=i4)

   select case(ph%nt_method)
   case(NTPOLY_TRS2)
      call TRS2(ham,ovlp,ne,dm,ph%ebs,ph%mu,ph%nt_options)
   case(NTPOLY_TRS4)
      call TRS4(ham,ovlp,ne,dm,ph%ebs,ph%mu,ph%nt_options)
   case(NTPOLY_HPCP)
      call HPCP(ham,ovlp,ne,dm,ph%ebs,ph%mu,ph%nt_options)
   end select

   call ScaleDistributedSparseMatrix(dm,ph%spin_degen)

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished density matrix purification"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine cleans up NTPoly.
!!
subroutine elsi_cleanup_ntpoly(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   character(len=40), parameter :: caller = "elsi_cleanup_ntpoly"

   if(ph%nt_started) then
      call DestructDistributedSparseMatrix(ph%nt_ham)
      call DestructDistributedSparseMatrix(ph%nt_ovlp)
      call DestructDistributedSparseMatrix(ph%nt_dm)
      call DestructPermutation(ph%nt_perm)
   endif

   ph%nt_started = .false.

end subroutine

end module ELSI_NTPOLY
