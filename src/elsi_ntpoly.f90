! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to NTPoly.
!!
module ELSI_NTPOLY

   use ELSI_CONSTANTS, only: NTPOLY_SOLVER,NTPOLY_PM,NTPOLY_TC2,NTPOLY_TRS4,&
       NTPOLY_HPCP
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_IO, only: elsi_say,elsi_get_time
   use ELSI_MPI, only: elsi_check_mpi,mpi_logical
   use ELSI_PRECISION, only: r8,i4
   use NTPOLY, only: PM,TRS2,TRS4,HPCP,EnergyDensityMatrix,LowdinExtrapolate,&
       Matrix_ps,ConstructEmptyMatrix,DestructMatrix,CopyMatrix,ScaleMatrix,&
       FillMatrixFromTripletList,GetMatrixTripletList,ProcessGrid_t,&
       ConstructNewProcessGrid,DestructProcessGrid,ConstructRandomPermutation,&
       DestructPermutation,InverseSquareRoot,SolverParameters_t,Triplet_r,&
       Triplet_c,TripletList_r,TripletList_c,ConstructTripletList,&
       AppendToTripletList,DestructTripletList,ActivateLogger,DeactivateLogger

   implicit none

   private

   public :: elsi_init_ntpoly
   public :: elsi_cleanup_ntpoly
   public :: elsi_solve_ntpoly
   public :: elsi_compute_edm_ntpoly
   public :: elsi_update_dm_ntpoly
   public :: Matrix_ps
   public :: ConstructEmptyMatrix
   public :: CopyMatrix
   public :: DestructMatrix
   public :: FillMatrixFromTripletList
   public :: GetMatrixTripletList
   public :: TripletList_r
   public :: TripletList_c
   public :: ConstructTripletList
   public :: AppendToTripletList
   public :: DestructTripletList
   public :: Triplet_r
   public :: Triplet_c
   public :: ProcessGrid_t

contains

!>
!! This routine initializes NTPoly.
!!
subroutine elsi_init_ntpoly(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh

   integer(kind=i4) :: n_prow
   integer(kind=i4) :: n_pcol
   integer(kind=i4) :: np_per_layer
   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_init_ntpoly"

   if(.not. ph%nt_started) then
      np_per_layer = bh%n_procs/ph%nt_n_layers

      ! Set 2D process grid
      do n_prow = nint(sqrt(real(np_per_layer,kind=r8)),kind=i4),1,-1
         if(mod(np_per_layer,n_prow) == 0) then
            n_pcol = np_per_layer/n_prow

            ! n_prow must be a multiple of n_pcol, or vice versa
            if(mod(max(n_prow,n_pcol),min(n_prow,n_pcol)) == 0) then
               exit
            end if
         end if
      end do

      call ConstructNewProcessGrid(ph%nt_pgrid,bh%comm,n_prow,n_pcol,&
           ph%nt_n_layers)

      ph%nt_n_prow = n_prow
      ph%nt_n_pcol = n_pcol

      call MPI_Bcast(ph%nt_output,1,mpi_logical,0,bh%comm,ierr)

      call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

      ph%nt_filter = max(ph%nt_filter,bh%def0)

      if(ph%nt_output .and. bh%myid_all == 0) then
         call ActivateLogger()
      end if

      ph%nt_started = .true.
   end if

end subroutine

!>
!! This routine interfaces to NTPoly.
!!
subroutine elsi_solve_ntpoly(ph,bh,ham,ovlp,dm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(Matrix_ps), intent(in) :: ham
   type(Matrix_ps), intent(inout) :: ovlp
   type(Matrix_ps), intent(inout) :: dm

   integer(kind=i4) :: ne
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   type(Matrix_ps) :: ovlp_isr

   character(len=*), parameter :: caller = "elsi_solve_ntpoly"

   if(ph%nt_first) then
      call elsi_get_time(t0)

      call DestructPermutation(ph%nt_perm)
      call ConstructRandomPermutation(ph%nt_perm,ovlp%logical_matrix_dimension,&
           ph%nt_pgrid)

      ph%nt_options = SolverParameters_t(ph%nt_tol,ph%nt_filter,ph%nt_max_iter,&
         ph%nt_output,ph%nt_perm)

      ! Overlap seems more difficult to converge
      ph%nt_options%threshold = max(0.01_r8*ph%nt_filter,1.0e-15_r8)

      call InverseSquareRoot(ovlp,ovlp_isr,ph%nt_options,ph%nt_isr)
      call CopyMatrix(ovlp_isr,ovlp)
      call DestructMatrix(ovlp_isr)

      ph%nt_options%threshold = ph%nt_filter

      call elsi_get_time(t1)

      write(msg,"(A)") "Finished overlap matrix inverse square root"
      call elsi_say(bh,msg)
      write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,msg)
   end if

   call elsi_get_time(t0)

   ne = nint(ph%n_electrons,kind=i4)

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

   call ScaleMatrix(dm,ph%spin_degen)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished density matrix purification"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   if(ph%decision_status == 1) then
      ph%decision_data(NTPOLY_SOLVER) = t1-t0
   end if

   ph%nt_first = .false.

end subroutine

!>
!! This routine computes the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_ntpoly(ph,bh,ham,edm)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(Matrix_ps), intent(inout) :: ham
   type(Matrix_ps), intent(inout) :: edm ! On entry: DM

   real(kind=r8) :: factor
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   type(Matrix_ps) :: tmp

   character(len=*), parameter :: caller = "elsi_compute_edm_ntpoly"

   call elsi_get_time(t0)

   factor = 1.0_r8/ph%spin_degen

   call CopyMatrix(edm,tmp)
   call EnergyDensityMatrix(ham,tmp,edm,ph%nt_filter)
   call DestructMatrix(tmp)
   call ScaleMatrix(edm,factor)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished energy density matrix calculation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! This routine extrapolates density matrix using Lowdin decomposition of the
!! old and new overlap matrices or 2nd order trace resetting purification.
!!
subroutine elsi_update_dm_ntpoly(ph,bh,ovlp0,ovlp1,dm0,dm1)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(Matrix_ps), intent(inout) :: ovlp0
   type(Matrix_ps), intent(inout) :: ovlp1
   type(Matrix_ps), intent(inout) :: dm0
   type(Matrix_ps), intent(inout) :: dm1

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_update_dm_ntpoly"

   call elsi_get_time(t0)

   if(ph%solver /= NTPOLY_SOLVER) then
      call DestructPermutation(ph%nt_perm)
      call ConstructRandomPermutation(ph%nt_perm,&
           ovlp0%logical_matrix_dimension,ph%nt_pgrid)

      ph%nt_options = SolverParameters_t(ph%nt_tol,ph%nt_filter,ph%nt_max_iter,&
         ph%nt_output,ph%nt_perm)
   end if

   call LowdinExtrapolate(dm0,ovlp0,ovlp1,dm1,ph%nt_options)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished density matrix extrapolation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! This routine cleans up NTPoly.
!!
subroutine elsi_cleanup_ntpoly(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   character(len=*), parameter :: caller = "elsi_cleanup_ntpoly"

   if(ph%nt_started) then
      call DeactivateLogger()
      call DestructMatrix(ph%nt_ham)
      call DestructMatrix(ph%nt_ovlp)
      call DestructMatrix(ph%nt_ovlp_copy)
      call DestructMatrix(ph%nt_dm)
      call DestructMatrix(ph%nt_map)
      call DestructPermutation(ph%nt_perm)
      call DestructProcessGrid(ph%nt_pgrid)
   end if

   ph%nt_first = .true.
   ph%nt_started = .false.

end subroutine

end module ELSI_NTPOLY
