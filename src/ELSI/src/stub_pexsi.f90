! Copyright (c) 2015-2017, the ELSI team. All rights reserved.
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
!! This module is only compiled when the actual PEXSI is not available,
!! to make the PEXSI part of ELSI compile.
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
   public :: f_ppexsi_calculate_edm_correction_real
   public :: f_ppexsi_calculate_edm_correction_complex
   public :: f_ppexsi_retrieve_real_dm
   public :: f_ppexsi_retrieve_complex_dm
   public :: f_ppexsi_retrieve_real_edm
   public :: f_ppexsi_retrieve_complex_edm
   public :: f_ppexsi_plan_finalize

   type, public, bind(C) :: f_ppexsi_options
      real(kind=r8)    :: temperature
      real(kind=r8)    :: gap
      real(kind=r8)    :: deltaE
      integer(kind=i4) :: numPole
      integer(kind=i4) :: isInertiaCount
      integer(kind=i4) :: maxPEXSIIter
      real(kind=r8)    :: muMin0
      real(kind=r8)    :: muMax0
      real(kind=r8)    :: mu0
      real(kind=r8)    :: muInertiaTolerance
      real(kind=r8)    :: muInertiaExpansion
      real(kind=r8)    :: muPEXSISafeGuard
      real(kind=r8)    :: numElectronPEXSITolerance
      integer(kind=i4) :: matrixType
      integer(kind=i4) :: isSymbolicFactorize
      integer(kind=i4) :: isConstructCommPattern
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

function f_ppexsi_plan_initialize(fcomm,numProcRow,numProcCol,outputFileIndex,info)

   implicit none

   integer(kind=i4)         :: fcomm
   integer(kind=i4)         :: numProcRow
   integer(kind=i4)         :: numProcCol
   integer(kind=i4)         :: outputFileIndex
   integer(kind=i4)         :: info
   integer(kind=c_intptr_t) :: f_ppexsi_plan_initialize

   f_ppexsi_plan_initialize = int(0,kind=c_intptr_t)

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end function

subroutine f_ppexsi_load_real_hs_matrix(plan,options,nrows,nnz,nnzLocal,&
              numColLocal,colptrLocal,rowindLocal,HnzvalLocal,isSIdentity,&
              SnzvalLocal,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options)   :: options
   integer(kind=i4)         :: nrows
   integer(kind=i4)         :: nnz
   integer(kind=i4)         :: nnzLocal
   integer(kind=i4)         :: numColLocal
   integer(kind=i4)         :: isSIdentity
   integer(kind=i4)         :: colptrLocal(*)
   integer(kind=i4)         :: rowindLocal(*)
   real(kind=r8)            :: HnzvalLocal(*)
   real(kind=r8)            :: SnzvalLocal(*)
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_load_complex_hs_matrix(plan,options,nrows,nnz,nnzLocal,&
              numColLocal,colptrLocal,rowindLocal,HnzvalLocal,isSIdentity,&
              SnzvalLocal,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options)   :: options
   integer(kind=i4)         :: nrows
   integer(kind=i4)         :: nnz
   integer(kind=i4)         :: nnzLocal
   integer(kind=i4)         :: numColLocal
   integer(kind=i4)         :: isSIdentity
   integer(kind=i4)         :: colptrLocal(*)
   integer(kind=i4)         :: rowindLocal(*)
   complex(kind=r8)         :: HnzvalLocal(*)
   complex(kind=r8)         :: SnzvalLocal(*)
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_symbolic_factorize_real_symmetric_matrix(plan,options,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options)   :: options
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_symbolic_factorize_complex_symmetric_matrix(plan,options,&
              info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options)   :: options
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix(plan,options,&
              AnzvalLocal,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options)   :: options
   complex(kind=r8)         :: AnzvalLocal(*)
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_inertia_count_real_matrix(plan,options,numShift,shiftList,&
              inertiaList,info)

   implicit none

   integer(c_intptr_t)    :: plan
   type(f_ppexsi_options) :: options
   integer(kind=i4)       :: numShift
   real(kind=r8)          :: shiftList(*)
   real(kind=r8)          :: inertiaList(*)
   integer(kind=i4)       :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_inertia_count_complex_matrix(plan,options,numShift,&
              shiftList,inertiaList,info)

   implicit none

   integer(c_intptr_t)    :: plan
   type(f_ppexsi_options) :: options
   integer(kind=i4)       :: numShift
   real(kind=r8)          :: shiftList(*)
   real(kind=r8)          :: inertiaList(*)
   integer(kind=i4)       :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_calculate_fermi_operator_real3(plan,options,mu,&
              numElectronExact,numElectronPEXSI,numElectronDrvMuPEXSI,info)

   implicit none

   integer(c_intptr_t)    :: plan
   type(f_ppexsi_options) :: options
   real(kind=r8)          :: mu
   real(kind=r8)          :: numElectronExact
   real(kind=r8)          :: numElectronPEXSI
   real(kind=r8)          :: numElectronDrvMuPEXSI
   integer(kind=i4)       :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_calculate_fermi_operator_complex(plan,options,mu,&
              numElectronExact,numElectronPEXSI,numElectronDrvMuPEXSI,info)

   implicit none

   integer(c_intptr_t)    :: plan
   type(f_ppexsi_options) :: options
   real(kind=r8)          :: mu
   real(kind=r8)          :: numElectronExact
   real(kind=r8)          :: numElectronPEXSI
   real(kind=r8)          :: numElectronDrvMuPEXSI
   integer(kind=i4)       :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_calculate_edm_correction_real(plan,options,info)

   implicit none

   integer(c_intptr_t)    :: plan
   type(f_ppexsi_options) :: options
   integer(kind=i4)       :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_calculate_edm_correction_complex(plan,options,info)

   implicit none

   integer(c_intptr_t)    :: plan
   type(f_ppexsi_options) :: options
   integer(kind=i4)       :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_retrieve_real_dm(plan,DMnzvalLocal,totalEnergyH,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   real(kind=r8)            :: DMnzvalLocal(*)
   real(kind=r8)            :: totalEnergyH
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_retrieve_complex_dm(plan,DMnzvalLocal,totalEnergyH,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   complex(kind=r8)         :: DMnzvalLocal(*)
   real(kind=r8)            :: totalEnergyH
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_retrieve_real_edm(plan,EDMnzvalLocal,totalEnergyS,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   real(kind=r8)            :: EDMnzvalLocal(*)
   real(kind=r8)            :: totalEnergyS
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_retrieve_complex_edm(plan,EDMnzvalLocal,totalEnergyS,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   complex(kind=r8)         :: EDMnzvalLocal(*)
   real(kind=r8)            :: totalEnergyS
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_plan_finalize(plan,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

end module F_PPEXSI_INTERFACE
