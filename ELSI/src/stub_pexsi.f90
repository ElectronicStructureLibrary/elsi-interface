!Copyright (c) 2016, ELSI consortium
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in the
!   documentation and/or other materials provided with the distribution.
! * Neither the name of the ELSI project nor the
!   names of its contributors may be used to endorse or promote products
!   derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! This module is only compiled when the actual PEXSI is not available, to make
! the PEXSI part of ELSI compile.

module f_ppexsi_interface

   use iso_c_binding

   implicit none

   public :: f_ppexsi_set_default_options
   public :: f_ppexsi_plan_initialize
   public :: f_ppexsi_load_real_hs_matrix
   public :: f_ppexsi_dft_driver
   public :: f_ppexsi_retrieve_real_dft_matrix
   public :: f_ppexsi_plan_finalize

   type, bind(C) :: f_ppexsi_options
      real(c_double) :: temperature
      real(c_double) :: gap
      real(c_double) :: deltaE
      integer(c_int) :: numPole
      integer(c_int) :: isInertiaCount
      integer(c_int) :: maxPEXSIIter
      real(c_double) :: muMin0
      real(c_double) :: muMax0
      real(c_double) :: mu0
      real(c_double) :: muInertiaTolerance
      real(c_double) :: muInertiaExpansion
      real(c_double) :: muPEXSISafeGuard
      real(c_double) :: numElectronPEXSITolerance
      integer(c_int) :: matrixType
      integer(c_int) :: isSymbolicFactorize
      integer(c_int) :: isConstructCommPattern
      integer(c_int) :: ordering
      integer(c_int) :: rowOrdering
      integer(c_int) :: npSymbFact
      integer(c_int) :: symmetric
      integer(c_int) :: transpose
      integer(c_int) :: verbosity
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

   integer :: fcomm
   integer(c_int) :: numProcRow, numProcCol,outputFileIndex
   integer(c_int) :: info
   integer(c_intptr_t) :: f_ppexsi_plan_initialize

end function

subroutine f_ppexsi_load_real_hs_matrix(plan,options,nrows,nnz,nnzLocal,&
              numColLocal,colptrLocal,rowindLocal,HnzvalLocal,isSIdentity,&
              SnzvalLocal,info)

   implicit none

   integer(c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   integer(c_int) :: nrows,nnz,nnzLocal,numColLocal,isSIdentity
   integer(c_int) :: colptrLocal(*),rowindLocal(*)
   real(c_double) :: HnzvalLocal(*),SnzvalLocal(*)
   integer(c_int) :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_dft_driver(plan,options,numElectronExact,muPEXSI,&
              numElectronPEXSI,muMinInertia,muMaxInertia,numTotalInertiaIter,&
              numTotalPEXSIIter,info)

   implicit none

   integer(c_intptr_t) :: plan
   type(f_ppexsi_options) :: options
   real(c_double) :: numElectronExact
   real(c_double) :: muPEXSI,numElectronPEXSI
   real(c_double) :: muMinInertia,muMaxInertia
   integer(c_int) :: numTotalInertiaIter,numTotalPEXSIIter
   integer(c_int) :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_ppexsi_retrieve_real_dft_matrix(plan,DMnzvalLocal,&
              EDMnzvalLocal,FDMnzvalLocal,totalEnergyH,totalEnergyS,&
              totalFreeEnergy,info)

   implicit none

   integer(c_intptr_t) :: plan
   real(c_double) :: DMnzvalLocal(*),EDMnzvalLocal(*),FDMnzvalLocal(*)
   real(c_double) :: totalEnergyH,totalEnergyS,totalFreeEnergy
   integer(c_int) :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine


subroutine f_ppexsi_plan_finalize(plan,info)
   implicit none

   integer(c_intptr_t) :: plan
   integer(c_int) :: info

   write(*,"(A)") " A PEXSI stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

end module
