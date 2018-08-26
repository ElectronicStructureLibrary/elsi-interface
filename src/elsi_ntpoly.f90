! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to NTPoly.
!!
module ELSI_NTPOLY

   use ELSI_CONSTANTS, only: UNSET,NTPOLY_PM,NTPOLY_TC2,NTPOLY_TRS4,NTPOLY_HPCP
   use ELSI_DATATYPE,  only: elsi_param_t,elsi_basic_t
   use ELSI_IO,        only: elsi_say,elsi_get_time
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,       only: elsi_check_mpi,mpi_logical
   use ELSI_PRECISION, only: r8,i4
   use NTPOLY,         only: PM,TRS2,TRS4,HPCP,EnergyDensityMatrix,&
                             DistributedSparseMatrix_t,FillFromTripletList,&
                             ConstructEmptyDistributedSparseMatrix,&
                             ScaleDistributedSparseMatrix,&
                             CopyDistributedSparseMatrix,&
                             FilterDistributedSparseMatrix,&
                             GetTripletList,DestructDistributedSparseMatrix,&
                             ConstructProcessGrid,DestructProcessGrid,&
                             ConstructRandomPermutation,DestructPermutation,&
                             InverseSquareRoot,IterativeSolverParameters_t,&
                             Triplet_t,TripletList_t,ConstructTripletList,&
                             AppendToTripletList,DestructTripletList

   implicit none

   private

   public :: elsi_init_ntpoly
   public :: elsi_cleanup_ntpoly
   public :: elsi_solve_ntpoly
   public :: elsi_compute_edm_ntpoly
   public :: DistributedSparseMatrix_t
   public :: ConstructEmptyDistributedSparseMatrix
   public :: FillFromTripletList
   public :: GetTripletList
   public :: TripletList_t
   public :: ConstructTripletList
   public :: AppendToTripletList
   public :: DestructTripletList
   public :: Triplet_t

   interface elsi_solve_ntpoly
      module procedure elsi_solve_ntpoly_real
   end interface

   interface elsi_compute_edm_ntpoly
      module procedure elsi_compute_edm_ntpoly_real
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
   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_init_ntpoly"

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

      call MPI_Bcast(ph%nt_output,1,mpi_logical,0,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

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

   type(DistributedSparseMatrix_t) :: ovlp_isr

   character(len=*), parameter :: caller = "elsi_solve_ntpoly_real"

   if(ph%n_calls == 1) then
      call elsi_get_time(t0)

      call ConstructRandomPermutation(ph%nt_perm,ovlp%logical_matrix_dimension)

      ph%nt_options = IterativeSolverParameters_t(ph%nt_tol,ph%nt_filter,&
                         ph%nt_max_iter,ph%nt_output,ph%nt_perm)

      ! Overlap seems more difficult to converge
      ph%nt_options%threshold = max(0.01_r8*ph%nt_filter,1.0e-15_r8)

      call InverseSquareRoot(ovlp,ovlp_isr,ph%nt_options,ph%nt_isr)
      call CopyDistributedSparseMatrix(ovlp_isr,ovlp)
      call DestructDistributedSparseMatrix(ovlp_isr)

      ph%nt_options%threshold = ph%nt_filter

      call elsi_get_time(t1)

      write(info_str,"(2X,A)") "Finished overlap matrix inverse square root"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,info_str)
   endif

   call elsi_get_time(t0)

   call ConstructEmptyDistributedSparseMatrix(dm,ph%n_basis)

   ne = int(ph%n_electrons,kind=i4)

   select case(ph%nt_method)
   case(NTPOLY_PM)
      call PM(ham,ovlp,ne,dm,ph%ebs,ph%mu,ph%nt_options)
   case(NTPOLY_TC2)
      call TRS2(ham,ovlp,ne,dm,ph%ebs,ph%mu,ph%nt_options)
   case(NTPOLY_TRS4)
      call TRS4(ham,ovlp,ne,dm,ph%ebs,ph%mu,ph%nt_options)
   case(NTPOLY_HPCP)
      call HPCP(ham,ovlp,ne,dm,ph%ebs,ph%mu,ph%nt_options)
   end select

   call ScaleDistributedSparseMatrix(dm,ph%spin_degen)

   if(ph%nt_filter < bh%def0) then
      call FilterDistributedSparseMatrix(dm,bh%def0)
   endif

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished density matrix purification"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine computes the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_ntpoly_real(ph,bh,ham,edm)

   implicit none

   type(elsi_param_t),              intent(in)    :: ph
   type(elsi_basic_t),              intent(in)    :: bh
   type(DistributedSparseMatrix_t), intent(inout) :: ham
   type(DistributedSparseMatrix_t), intent(inout) :: edm ! On entry: DM

   real(kind=r8)      :: factor
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   character(len=200) :: info_str

   type(DistributedSparseMatrix_t) :: tmp

   character(len=*), parameter :: caller = "elsi_compute_edm_ntpoly_real"

   call elsi_get_time(t0)

   factor = 1.0_r8/ph%spin_degen

   call ConstructEmptyDistributedSparseMatrix(tmp,ph%n_basis)
   call CopyDistributedSparseMatrix(edm,tmp)
   call EnergyDensityMatrix(ham,tmp,edm,ph%nt_options)
   call DestructDistributedSparseMatrix(tmp)
   call ScaleDistributedSparseMatrix(edm,factor)

   if(ph%nt_filter < bh%def0) then
      call FilterDistributedSparseMatrix(edm,bh%def0)
   endif

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished energy density matrix calculation"
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

   character(len=*), parameter :: caller = "elsi_cleanup_ntpoly"

   if(ph%nt_started) then
      call DestructDistributedSparseMatrix(ph%nt_ham)
      call DestructDistributedSparseMatrix(ph%nt_ovlp)
      call DestructDistributedSparseMatrix(ph%nt_dm)
      call DestructPermutation(ph%nt_perm)
      call DestructProcessGrid()
   endif

   ph%nt_started = .false.

end subroutine

end module ELSI_NTPOLY
