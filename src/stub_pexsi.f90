! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide stub routines which are only compiled when the actual PEXSI is not
!! available.
!!
module F_PPEXSI_INTERFACE

   use, intrinsic :: ISO_C_BINDING
   use ELSI_PRECISION, only: r8,i4

   implicit none

   private

   public :: f_ppexsi_set_default_options
   public :: f_ppexsi_plan_initialize
   public :: f_ppexsi_load_real_hs_matrix
   public :: f_ppexsi_load_complex_hs_matrix
   public :: f_ppexsi_symbolic_factorize_real_symmetric_matrix
   public :: f_ppexsi_symbolic_factorize_complex_symmetric_matrix
   public :: f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix
   public :: f_ppexsi_inertia_count_real_matrix
   public :: f_ppexsi_inertia_count_complex_matrix
   public :: f_ppexsi_calculate_fermi_operator_real3
   public :: f_ppexsi_calculate_fermi_operator_complex
   public :: f_ppexsi_retrieve_real_dm
   public :: f_ppexsi_retrieve_complex_dm
   public :: f_ppexsi_retrieve_real_edm
   public :: f_ppexsi_retrieve_complex_edm
   public :: f_ppexsi_plan_finalize

   type, public :: f_ppexsi_options
      real(kind=r8) :: spin
      real(kind=r8) :: temperature
      real(kind=r8) :: gap
      real(kind=r8) :: deltaE
      integer(kind=i4) :: numPole
      integer(kind=i4) :: isInertiaCount
      integer(kind=i4) :: maxPEXSIIter
      real(kind=r8) :: muMin0
      real(kind=r8) :: muMax0
      real(kind=r8) :: mu0
      real(kind=r8) :: muInertiaTolerance
      real(kind=r8) :: muInertiaExpansion
      real(kind=r8) :: muPEXSISafeGuard
      real(kind=r8) :: numElectronPEXSITolerance
      integer(kind=i4) :: matrixType
      integer(kind=i4) :: isSymbolicFactorize
      integer(kind=i4) :: isConstructCommPattern
      integer(kind=i4) :: solver
      integer(kind=i4) :: symmetricStorage
      integer(kind=i4) :: ordering
      integer(kind=i4) :: rowOrdering
      integer(kind=i4) :: npSymbFact
      integer(kind=i4) :: symmetric
      integer(kind=i4) :: transpose
      integer(kind=i4) :: method
      integer(kind=i4) :: nPoints
      integer(kind=i4) :: verbosity
   end type

contains

subroutine f_ppexsi_set_default_options(options)

   implicit none

   type(f_ppexsi_options) :: options

   ! Do not stop here!
   options%numPole = 1
   options%nPoints = 1

end subroutine

function f_ppexsi_plan_initialize(fcomm,numProcRow,numProcCol,outputFileIndex,&
   info)

   implicit none

   integer(kind=i4) :: fcomm
   integer(kind=i4) :: numProcRow
   integer(kind=i4) :: numProcCol
   integer(kind=i4) :: outputFileIndex
   integer(kind=i4) :: info
   integer(kind=c_intptr_t) :: f_ppexsi_plan_initialize

   f_ppexsi_plan_initialize = int(0,kind=c_intptr_t)

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end function

subroutine f_ppexsi_load_real_hs_matrix(plan,options,nrows,nnz,nnzLocal,&
   numColLocal,colptrLocal,rowindLocal,HnzvalLocal,isSIdentity,SnzvalLocal,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   integer(kind=i4) :: nrows
   integer(kind=i4) :: nnz
   integer(kind=i4) :: nnzLocal
   integer(kind=i4) :: numColLocal
   integer(kind=i4) :: isSIdentity
   integer(kind=i4) :: colptrLocal(*)
   integer(kind=i4) :: rowindLocal(*)
   real(kind=r8) :: HnzvalLocal(*)
   real(kind=r8) :: SnzvalLocal(*)
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_load_complex_hs_matrix(plan,options,nrows,nnz,nnzLocal,&
   numColLocal,colptrLocal,rowindLocal,HnzvalLocal,isSIdentity,SnzvalLocal,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   integer(kind=i4) :: nrows
   integer(kind=i4) :: nnz
   integer(kind=i4) :: nnzLocal
   integer(kind=i4) :: numColLocal
   integer(kind=i4) :: isSIdentity
   integer(kind=i4) :: colptrLocal(*)
   integer(kind=i4) :: rowindLocal(*)
   complex(kind=r8) :: HnzvalLocal(*)
   complex(kind=r8) :: SnzvalLocal(*)
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_symbolic_factorize_real_symmetric_matrix(plan,options,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_symbolic_factorize_complex_symmetric_matrix(plan,options,&
   info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix(plan,options,&
   AnzvalLocal,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   complex(kind=r8) :: AnzvalLocal(*)
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_inertia_count_real_matrix(plan,options,numShift,shiftList,&
   inertiaList,info)

   implicit none

   integer(c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   integer(kind=i4) :: numShift
   real(kind=r8) :: shiftList(*)
   real(kind=r8) :: inertiaList(*)
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_inertia_count_complex_matrix(plan,options,numShift,&
   shiftList,inertiaList,info)

   implicit none

   integer(c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   integer(kind=i4) :: numShift
   real(kind=r8) :: shiftList(*)
   real(kind=r8) :: inertiaList(*)
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_calculate_fermi_operator_real3(plan,options,&
   numElectronExact,mu,numElectronPEXSI,info)

   implicit none

   integer(c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   real(kind=r8) :: numElectronExact
   real(kind=r8) :: mu
   real(kind=r8) :: numElectronPEXSI
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_calculate_fermi_operator_complex(plan,options,mu,&
   numElectronExact,numElectronPEXSI,numElectronDrvMuPEXSI,info)

   implicit none

   integer(c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   real(kind=r8) :: mu
   real(kind=r8) :: numElectronExact
   real(kind=r8) :: numElectronPEXSI
   real(kind=r8) :: numElectronDrvMuPEXSI
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_retrieve_real_dm(plan,DMnzvalLocal,totalEnergyH,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   real(kind=r8) :: DMnzvalLocal(*)
   real(kind=r8) :: totalEnergyH
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_retrieve_complex_dm(plan,DMnzvalLocal,totalEnergyH,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   complex(kind=r8) :: DMnzvalLocal(*)
   real(kind=r8) :: totalEnergyH
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_retrieve_real_edm(plan,options,EDMnzvalLocal,totalEnergyS,&
   info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   real(kind=r8) :: EDMnzvalLocal(*)
   real(kind=r8) :: totalEnergyS
   integer(kind=i4):: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_retrieve_complex_edm(plan,options,EDMnzvalLocal,&
   totalEnergyS,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   complex(kind=r8) :: EDMnzvalLocal(*)
   real(kind=r8) :: totalEnergyS
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

subroutine f_ppexsi_plan_finalize(plan,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   integer(kind=i4) :: info

   write(*,"(A)") "**Error! A PEXSI stub routine was called"
   stop

end subroutine

end module F_PPEXSI_INTERFACE
