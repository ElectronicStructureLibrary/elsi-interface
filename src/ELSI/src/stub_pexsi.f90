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

module f_ppexsi_interface

   use iso_c_binding
   use ELSI_PRECISION, only: r8,i4

   implicit none

   public :: f_ppexsi_set_default_options
   public :: f_ppexsi_plan_initialize
   public :: f_ppexsi_load_real_hs_matrix
   public :: f_ppexsi_dft_driver
   public :: f_ppexsi_dft_driver3
   public :: f_ppexsi_retrieve_real_dft_matrix
   public :: f_ppexsi_retrieve_real_dft_matrix2
   public :: f_ppexsi_plan_finalize

   type, bind(C) :: f_ppexsi_options
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
      integer(kind=i4) :: verbosity
   end type f_ppexsi_options

contains

subroutine f_ppexsi_set_default_options(options)

   implicit none

   type(f_ppexsi_options) :: options

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

function f_ppexsi_plan_initialize(fcomm,numProcRow,numProcCol,outputFileIndex,info)

   implicit none

   integer                  :: fcomm
   integer(kind=i4)         :: numProcRow,numProcCol,outputFileIndex
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
   integer(kind=i4)         :: nrows,nnz,nnzLocal,numColLocal,isSIdentity
   integer(kind=i4)         :: colptrLocal(*),rowindLocal(*)
   real(kind=r8)            :: HnzvalLocal(*),SnzvalLocal(*)
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_dft_driver(plan,options,numElectronExact,muPEXSI,&
              numElectronPEXSI,muMinInertia,muMaxInertia,numTotalInertiaIter,&
              numTotalPEXSIIter,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options)   :: options
   real(kind=r8)            :: numElectronExact
   real(kind=r8)            :: muPEXSI,numElectronPEXSI
   real(kind=r8)            :: muMinInertia,muMaxInertia
   integer(kind=i4)         :: numTotalInertiaIter,numTotalPEXSIIter
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_dft_driver3(plan,options,numElectronExact,pexsi_driver,&
              n_mu_points,muPEXSI,numElectronPEXSI,numTotalInertiaIter,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   type(f_ppexsi_options)   :: options
   real(kind=r8)            :: numElectronExact
   integer(kind=i4)         :: pexsi_driver,n_mu_points
   real(kind=r8)            :: muPEXSI,numElectronPEXSI
   integer(kind=i4)         :: numTotalInertiaIter
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_retrieve_real_dft_matrix(plan,DMnzvalLocal,&
              EDMnzvalLocal,FDMnzvalLocal,totalEnergyH,totalEnergyS,&
              totalFreeEnergy,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   real(kind=r8)            :: DMnzvalLocal(*),EDMnzvalLocal(*),FDMnzvalLocal(*)
   real(kind=r8)            :: totalEnergyH,totalEnergyS,totalFreeEnergy
   integer(kind=i4)         :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_retrieve_real_dft_matrix2(plan,DMnzvalLocal,&
              EDMnzvalLocal,FDMnzvalLocal,totalEnergyH,totalEnergyS,&
              totalFreeEnergy,info)

   implicit none

   integer(kind=c_intptr_t) :: plan
   real(kind=r8)            :: DMnzvalLocal(*),EDMnzvalLocal(*),FDMnzvalLocal(*)
   real(kind=r8)            :: totalEnergyH,totalEnergyS,totalFreeEnergy
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

end module
