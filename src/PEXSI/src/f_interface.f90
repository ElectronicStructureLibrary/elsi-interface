!	 Copyright (c) 2012 The Regents of the University of California,
!	 through Lawrence Berkeley National Laboratory.
!
!  Author: Lin Lin
!	 
!  This file is part of PEXSI. All rights reserved.
!
!	 Redistribution and use in source and binary forms, with or without
!	 modification, are permitted provided that the following conditions are met:
!
!	 (1) Redistributions of source code must retain the above copyright notice, this
!	 list of conditions and the following disclaimer.
!	 (2) Redistributions in binary form must reproduce the above copyright notice,
!	 this list of conditions and the following disclaimer in the documentation
!	 and/or other materials provided with the distribution.
!	 (3) Neither the name of the University of California, Lawrence Berkeley
!	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
!	 be used to endorse or promote products derived from this software without
!	 specific prior written permission.
!
!	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
!	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!	 You are under no obligation whatsoever to provide any bug fixes, patches, or
!	 upgrades to the features, functionality or performance of the source code
!	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
!	 available either publicly, or directly to Lawrence Berkeley National
!	 Laboratory, without imposing a separate written license agreement for such
!	 Enhancements, then you hereby grant the following license: a non-exclusive,
!	 royalty-free perpetual license to install, use, modify, prepare derivative
!	 works, incorporate into other computer software, distribute, and sublicense
!	 such enhancements or derivative works thereof, in binary and source code form.
!
!> @file f_interface.f90
!> @brief FORTRAN interface for %PEXSI library using ISO_C_BINDING.
!>
!> The ISO_C_BINDING feature is included in the FORTRAN 2003 standard, and is
!> implemented in most modern compilers.
!> 
!>
!> @note (From Alberto Garcia, 2013-04-10) Array arguments are *required* by the
!> standard to be explicit shape (e.g. a(3)), or assumed size (e.g. a(*)). This
!> avoids problems with the assumed shape specification (e.g. a(:)), which would
!> need to pass extra information in general. This permits the use of pointers
!> and array sections as actual arguments. The compiler would presumably make
!> on-the-fly copies in the case of non-contiguous data (a situation that should
!> be avoided on performance grounds).
!>
!> @see  c_pexsi_interface.h
!> @date 2014-04-01 Initial version
!> @date 2015-01-21 Interface with unsymmetric version of selected inversion.

! *********************************************************************
! Module for main PEXSI interface routines
! *********************************************************************

module f_ppexsi_interface
use, intrinsic :: iso_c_binding

! Struct for PPEXSIOptions
! NOTE: The order and the type of the parameters must be strictly the same as in PPEXSIOptions in 
! c_pexsi_interface.h
type, bind(C) :: f_ppexsi_options
  real(c_double)         :: temperature
  real(c_double)         :: gap
  real(c_double)         :: deltaE
  integer(c_int)         :: numPole
  integer(c_int)         :: isInertiaCount
  integer(c_int)         :: maxPEXSIIter
  real(c_double)         :: muMin0
  real(c_double)         :: muMax0
  real(c_double)         :: mu0
  real(c_double)         :: muInertiaTolerance
  real(c_double)         :: muInertiaExpansion
  real(c_double)         :: muPEXSISafeGuard
  real(c_double)         :: numElectronPEXSITolerance
  integer(c_int)         :: matrixType
  integer(c_int)         :: isSymbolicFactorize
  integer(c_int)         :: isConstructCommPattern
  integer(c_int)         :: solver
  integer(c_int)         :: ordering
  integer(c_int)         :: rowOrdering
  integer(c_int)         :: npSymbFact
  integer(c_int)         :: symmetric
  integer(c_int)         :: transpose
  integer(c_int)         :: method
  integer(c_int)         :: nPoints
  integer(c_int)         :: verbosity
end type f_ppexsi_options

interface  

  subroutine f_ppexsi_set_default_options(&
      options) &
      bind(C, Name="PPEXSISetDefaultOptions")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    type( f_ppexsi_options ), intent(out) :: options
  end subroutine 
 
  subroutine f_read_distsparsematrix_formatted_head (&
      filename,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      fcomm ) &
      bind(C, Name="f_read_distsparsematrix_formatted_head")
    use, intrinsic :: iso_c_binding
    implicit none
    character(c_char), intent(in) :: filename(*)
    integer(c_int), intent(out)        :: nrows, nnz, nnzLocal, numColLocal
    integer, intent(in)                :: fcomm
  end subroutine 

  subroutine f_read_distsparsematrix_formatted(&
      filename,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      colptrLocal,&
      rowindLocal,&
      nzvalLocal,&
      fcomm ) &
      bind(C, Name="f_read_distsparsematrix_formatted")
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char), intent(in) :: filename(*)
    integer(c_int), intent(in), value  :: nrows, nnz, nnzLocal, numColLocal
    integer(c_int), intent(out)        :: colptrLocal(*), rowindLocal(*)
    real(c_double), intent(out)        :: nzvalLocal(*)
    integer, intent(in)                :: fcomm
  end subroutine 

  subroutine f_read_distsparsematrix_head (&
      filename,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      fcomm ) &
      bind(C, Name="f_read_distsparsematrix_head")
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char), intent(in) :: filename(*)
    integer(c_int), intent(out)        :: nrows, nnz, nnzLocal, numColLocal
    integer, intent(in)                :: fcomm
  end subroutine 

  subroutine f_para_read_distsparsematrix(&
      filename,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      colptrLocal,&
      rowindLocal,&
      nzvalLocal,&
      fcomm ) &
      bind(C, Name="f_para_read_distsparsematrix")
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char), intent(in) :: filename(*)
    integer(c_int), intent(in), value  :: nrows, nnz, nnzLocal, numColLocal
    integer(c_int), intent(out)        :: colptrLocal(*), rowindLocal(*)
    real(c_double), intent(out)        :: nzvalLocal(*)
    integer, intent(in)                :: fcomm
  end subroutine 

  function f_ppexsi_plan_initialize(&
      fcomm,&
      numProcRow,&
      numProcCol,&
      outputFileIndex,&
      info) &
      bind(C, Name="f_ppexsi_plan_initialize")
    use, intrinsic :: iso_c_binding
    implicit none
    integer, intent(in)                :: fcomm
    integer(c_int), intent(in), value  :: numProcRow, numProcCol,outputFileIndex
    integer(c_int), intent(out)        :: info
    integer(c_intptr_t)                :: f_ppexsi_plan_initialize
  end function

  subroutine f_ppexsi_load_real_hs_matrix(&
      plan,&
      options,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      colptrLocal,&
      rowindLocal,&
      HnzvalLocal,&
      isSIdentity,&
      SnzvalLocal,&
      info) &
      bind(C, Name="PPEXSILoadRealHSMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    integer(c_int), value, intent(in)  :: &
      nrows, nnz, nnzLocal, numColLocal, isSIdentity
    integer(c_int), intent(in)  :: &
      colptrLocal(*), rowindLocal(*)
    real(c_double), intent(in)  :: &
      HnzvalLocal(*), SnzvalLocal(*)
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_load_complex_hs_matrix(&
      plan,&
      options,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      colptrLocal,&
      rowindLocal,&
      HnzvalLocal,&
      isSIdentity,&
      SnzvalLocal,&
      info) &
      bind(C, Name="PPEXSILoadComplexHSMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    integer(c_int), value, intent(in)  :: &
      nrows, nnz, nnzLocal, numColLocal, isSIdentity
    integer(c_int), intent(in)  :: &
      colptrLocal(*), rowindLocal(*)
    complex(c_double), intent(in)  :: &
      HnzvalLocal(*), SnzvalLocal(*)
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_symbolic_factorize_real_symmetric_matrix(&
      plan,&
      options,&
      info) &
      bind(C, Name="PPEXSISymbolicFactorizeRealSymmetricMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_symbolic_factorize_complex_symmetric_matrix(&
      plan,&
      options,&
      info) &
      bind(C, Name="PPEXSISymbolicFactorizeComplexSymmetricMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_symbolic_factorize_real_unsymmetric_matrix(&
      plan,&
      options,&
      AnzvalLocal,&
      info) &
      bind(C, Name="PPEXSISymbolicFactorizeRealUnsymmetricMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    real(c_double), intent(inout)          :: AnzvalLocal(*)
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix(&
      plan,&
      options,&
      AnzvalLocal,&
      info) &
      bind(C, Name="PPEXSISymbolicFactorizeComplexUnsymmetricMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    complex(c_double), intent(inout)       :: AnzvalLocal(*)
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_inertia_count_real_matrix(&
      plan,&
      options,&
      numShift,&
      shiftList,&
      inertiaList,&
      info) &
      bind(C, Name="PPEXSIInertiaCountRealMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    integer(c_int), value, intent(in)      :: numShift
    real(c_double), intent(in)             :: shiftList(*)
    real(c_double), intent(out)            :: inertiaList(*)
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_inertia_count_complex_matrix(&
      plan,&
      options,&
      numShift,&
      shiftList,&
      inertiaList,&
      info) &
      bind(C, Name="PPEXSIInertiaCountComplexMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    integer(c_int), value, intent(in)      :: numShift
    real(c_double), intent(in)             :: shiftList(*)
    real(c_double), intent(out)            :: inertiaList(*)
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_calculate_fermi_operator_real(&
      plan,&
      options,&
      mu,&
      numElectronExact,&
      numElectronPEXSI,&
      numElectronDrvMuPEXSI,&
      info) &
      bind(C, Name="PPEXSICalculateFermiOperatorReal")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    real(c_double), value, intent(in)      :: mu, numElectronExact
    real(c_double), intent(out)   :: numElectronPEXSI, numElectronDrvMuPEXSI
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_calculate_fermi_operator_real3(&
      plan,&
      options,&
      mu,&
      numElectronExact,&
      numElectronPEXSI,&
      numElectronDrvMuPEXSI,&
      info) &
      bind(C, Name="PPEXSICalculateFermiOperatorReal3")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value    :: plan
    type(f_ppexsi_options), value, intent(in) :: options
    real(c_double), value, intent(in)         :: mu, numElectronExact
    real(c_double), intent(out)               :: numElectronPEXSI, numElectronDrvMuPEXSI
    integer(c_int), intent(out)               :: info
  end subroutine

  subroutine f_ppexsi_calculate_fermi_operator_complex(&
      plan,&
      options,&
      mu,&
      numElectronExact,&
      numElectronPEXSI,&
      numElectronDrvMuPEXSI,&
      info) &
      bind(C, Name="PPEXSICalculateFermiOperatorComplex")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    real(c_double), value, intent(in)      :: mu, numElectronExact
    real(c_double), intent(out)   :: numElectronPEXSI, numElectronDrvMuPEXSI
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_selinv_real_symmetric_matrix(&
      plan,&
      options,&
      AnzvalLocal,&
      AinvnzvalLocal,&
      info) &
      bind(C, Name="PPEXSISelInvRealSymmetricMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    real(c_double), intent(in)  :: AnzvalLocal(*)
    real(c_double), intent(out) :: AinvnzvalLocal(*)
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_selinv_complex_symmetric_matrix(&
      plan,&
      options,&
      AnzvalLocal,&
      AinvnzvalLocal,&
      info) &
      bind(C, Name="PPEXSISelInvComplexSymmetricMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    real(c_double), intent(in)  :: AnzvalLocal(*)
    real(c_double), intent(out) :: AinvnzvalLocal(*)
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_selinv_real_unsymmetric_matrix(&
      plan,&
      options,&
      AnzvalLocal,&
      AinvnzvalLocal,&
      info) &
      bind(C, Name="PPEXSISelInvRealUnsymmetricMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    real(c_double), intent(in)  :: AnzvalLocal(*)
    real(c_double), intent(out) :: AinvnzvalLocal(*)
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_selinv_complex_unsymmetric_matrix(&
      plan,&
      options,&
      AnzvalLocal,&
      AinvnzvalLocal,&
      info) &
      bind(C, Name="PPEXSISelInvComplexUnsymmetricMatrix")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    real(c_double), intent(in)  :: AnzvalLocal(*)
    real(c_double), intent(out) :: AinvnzvalLocal(*)
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_dft_driver(&
      plan,&
      options,&
      numElectronExact,&
      muPEXSI,&
      numElectronPEXSI,&
      muMinInertia,&
      muMaxInertia,& 
      numTotalInertiaIter,&
      numTotalPEXSIIter,&
      info) &
      bind(C, Name="PPEXSIDFTDriver")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    real(c_double), value, intent(in)      :: numElectronExact
    real(c_double), intent(out)   :: muPEXSI, numElectronPEXSI
    real(c_double), intent(out)   :: muMinInertia, muMaxInertia
    integer(c_int), intent(out)   :: numTotalInertiaIter, numTotalPEXSIIter
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_dft_driver2_deprecate(&
      plan,&
      options,&
      numElectronExact,&
      muPEXSI,&
      numElectronPEXSI,&
      muMinInertia,&
      muMaxInertia,& 
      numTotalInertiaIter,&
      numTotalPEXSIIter,&
      info) &
      bind(C, Name="PPEXSIDFTDriver2_Deprecate")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    real(c_double), value, intent(in)      :: numElectronExact
    real(c_double), intent(out)   :: muPEXSI, numElectronPEXSI
    real(c_double), intent(out)   :: muMinInertia, muMaxInertia
    integer(c_int), intent(out)   :: numTotalInertiaIter, numTotalPEXSIIter
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_dft_driver2(&
      plan,&
      options,&
      numElectronExact,&
      method,&
      npoints,&
      muPEXSI,&
      numElectronPEXSI,&
      !muMinInertia,&
      !muMaxInertia,&
      numTotalInertiaIter,&
      info) &
      bind(C, Name="PPEXSIDFTDriver2")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t),    value, intent(in)    :: plan
    type(f_ppexsi_options),        intent(inout) :: options
    real(c_double),         value, intent(in)    :: numElectronExact
    integer(c_int),         value, intent(in)    :: method
    integer(c_int),         value, intent(in)    :: npoints
    real(c_double),                intent(out)   :: muPEXSI, numElectronPEXSI
    integer(c_int),                intent(out)   :: numTotalInertiaIter
    integer(c_int),                intent(out)   :: info
  end subroutine

  subroutine f_ppexsi_retrieve_real_dft_matrix(&
      plan,&
      DMnzvalLocal,&
      EDMnzvalLocal,&
      FDMnzvalLocal,&
      totalEnergyH,&
      totalEnergyS,&
      totalFreeEnergy,&
      info) &
      bind(C, Name="PPEXSIRetrieveRealDFTMatrix")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    real(c_double), intent(out)   :: DMnzvalLocal(*), EDMnzvalLocal(*), FDMnzvalLocal(*)
    real(c_double), intent(out)   :: totalEnergyH, totalEnergyS, totalFreeEnergy
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_retrieve_real_dft_matrix2(&
      plan,&
      DMnzvalLocal,&
      EDMnzvalLocal,&
      FDMnzvalLocal,&
      totalEnergyH,&
      totalEnergyS,&
      totalFreeEnergy,&
      info) &
      bind(C, Name="PPEXSIRetrieveRealDFTMatrix2")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), value, intent(in)  :: plan
    real(c_double),             intent(out) :: DMnzvalLocal(*), EDMnzvalLocal(*), FDMnzvalLocal(*)
    real(c_double),             intent(out) :: totalEnergyH, totalEnergyS, totalFreeEnergy
    integer(c_int),             intent(out) :: info
  end subroutine

  subroutine f_ppexsi_retrieve_complex_dft_matrix(&
      plan,&
      DMnzvalLocal,&
      EDMnzvalLocal,&
      FDMnzvalLocal,&
      totalEnergyH,&
      totalEnergyS,&
      totalFreeEnergy,&
      info) &
      bind(C, Name="PPEXSIRetrieveComplexDFTMatrix")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    real(c_double), intent(out)   :: DMnzvalLocal(*), EDMnzvalLocal(*), FDMnzvalLocal(*)
    real(c_double), intent(out)   :: totalEnergyH, totalEnergyS, totalFreeEnergy
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_plan_finalize(&
      plan,&
      info) &
      bind(C, Name="PPEXSIPlanFinalize")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_get_pole_dm(&
      zshift,&
      zweight,&
      Npole,&
      temp,&
      gap,&
      deltaE,&
      mu) &
      bind(C, Name="PPEXSIGetPoleDM")
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), intent(out)  :: zshift(*), zweight(*)
    integer(c_int), intent(in), value  :: Npole
    real(c_double), intent(in), value  :: temp, gap, deltaE, mu
  end subroutine 

  subroutine f_ppexsi_get_pole_edm(&
      zshift,&
      zweight,&
      Npole,&
      temp,&
      gap,&
      deltaE,&
      mu) &
      bind(C, Name="PPEXSIGetPoleEDM")
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), intent(out)  :: zshift(*), zweight(*)
    integer(c_int), intent(in), value  :: Npole
    real(c_double), intent(in), value  :: temp, gap, deltaE, mu
  end subroutine 

  subroutine f_ppexsi_get_pole_fdm(&
      zshift,&
      zweight,&
      Npole,&
      temp,&
      gap,&
      deltaE,&
      mu) &
      bind(C, Name="PPEXSIGetPoleFDM")
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), intent(out)  :: zshift(*), zweight(*)
    integer(c_int), intent(in), value  :: Npole
    real(c_double), intent(in), value  :: temp, gap, deltaE, mu
  end subroutine 

  subroutine f_ppexsi_calculate_edm_correction_complex(&
      plan,&
      options,&
      info) &
      bind(C, Name="PPEXSICalculateEDMCorrectionComplex")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    integer(c_int), intent(out)        :: info
  end subroutine     

  subroutine f_ppexsi_calculate_edm_correction_real(&
      plan,&
      options,&
      info) &
      bind(C, Name="PPEXSICalculateEDMCorrectionReal")
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), value, intent(in) :: options
    integer(c_int), intent(out)        :: info
  end subroutine     

  subroutine f_ppexsi_interpolate_dm_real(&
      plan,&
      options,&
      numElectronExact,&
      numElectronPEXSI,&
      NeVec,&
      muPEXSI, &
      info)&
      bind(C, Name="PPEXSIInterpolateDMReal") 
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), intent(inout) :: options
    real(c_double), intent(in)        :: numElectronExact 
    real(c_double), intent(in)        :: numElectronPEXSI
    real(c_double), intent(out)        :: NeVec(*)
    real(c_double), intent(out)        :: muPEXSI
    integer(c_int), intent(out)        :: info
  end subroutine     

  subroutine f_ppexsi_interpolate_dm_complex(&
      plan,&
      options,&
      numElectronExact,&
      numElectronPEXSI,&
      NeVec,&
      muPEXSI, &
      info)&
      bind(C, Name="PPEXSIInterpolateDMComplex") 
    use, intrinsic :: iso_c_binding
    import         :: f_ppexsi_options
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    type( f_ppexsi_options ), intent(inout) :: options
    real(c_double), intent(in)        :: numElectronExact 
    real(c_double), intent(in)        :: numElectronPEXSI
    real(c_double), intent(out)        :: NeVec(*)
    real(c_double), intent(out)        :: muPEXSI
    integer(c_int), intent(out)        :: info
  end subroutine     

  subroutine f_ppexsi_retrieve_real_dm_edm(&
      plan,&
      DMnzvalLocal,&
      EDMnzvalLocal,&
      nnzLocal,&
      info) &
      bind(C, Name="PPEXSIRetrieveRealDMEDM")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), value, intent(in)  :: plan
    real(c_double),             intent(out) :: DMnzvalLocal(*), EDMnzvalLocal(*)
    integer(c_int),             intent(out) :: nnzLocal
    integer(c_int),             intent(out) :: info
  end subroutine

! WPH: f_ppexsi_retrieve_complex_dft_matrix already binds to
!      PPEXSIRetrieveComplexDFTMatrix, and the Fortran 2003 standard
!      specifies that "the binding label of an entity must not be the
!      same as the binding label of another entity", thus this does
!      not compile with XLF.  Additionally, the interface does not match
!      PPEXSIRetrieveComplexDFTMatrix.
!      Since this subroutine is not being used with ELSI,commenting out.
!  subroutine f_ppexsi_retrieve_complex_dm_edm(&
!      plan,&
!      DMnzvalLocal,&
!      EDMnzvalLocal,&
!      nnzLocal,&
!      info) &
!      bind(C, Name="PPEXSIRetrieveComplexDFTMatrix")
!    use, intrinsic :: iso_c_binding
!    implicit none
!    integer(c_intptr_t), value, intent(in)  :: plan
!    real(c_double),             intent(out) :: DMnzvalLocal(*), EDMnzvalLocal(*)
!    integer(c_int),             intent(out) :: nnzLocal
!    integer(c_int),             intent(out) :: info
!  end subroutine 

   subroutine f_ppexsi_retrieve_real_dm(&
      plan,&
      DMnzvalLocal,&
      totalEnergyH,&
      info) &
      bind(C, Name="PPEXSIRetrieveRealDM")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    real(c_double), intent(out)   :: DMnzvalLocal(*)
    real(c_double), intent(out)   :: totalEnergyH
    integer(c_int), intent(out)   :: info
  end subroutine 

  subroutine f_ppexsi_retrieve_real_edm(&
      plan,&
      EDMnzvalLocal,&
      totalEnergyS,&
      info) &
      bind(C, Name="PPEXSIRetrieveRealEDM")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    real(c_double), intent(out)   :: EDMnzvalLocal(*)
    real(c_double), intent(out)   :: totalEnergyS
    integer(c_int), intent(out)   :: info
  end subroutine 

  subroutine f_ppexsi_retrieve_real_fdm(&
      plan,&
      FDMnzvalLocal,&
      totalFreeEnergy,&
      info) &
      bind(C, Name="PPEXSIRetrieveRealFDM")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    real(c_double), intent(out)   :: FDMnzvalLocal(*)
    real(c_double), intent(out)   :: totalFreeEnergy
    integer(c_int), intent(out)   :: info
  end subroutine 

  subroutine f_ppexsi_retrieve_complex_dm(&
      plan,&
      DMnzvalLocal,&
      totalEnergyH,&
      info) &
      bind(C, Name="PPEXSIRetrieveComplexDM")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    complex(c_double), intent(out)   :: DMnzvalLocal(*)
    real(c_double), intent(out)   :: totalEnergyH
    integer(c_int), intent(out)            :: info
  end subroutine 

  subroutine f_ppexsi_retrieve_complex_edm(&
      plan,&
      EDMnzvalLocal,&
      totalEnergyS,&
      info) &
      bind(C, Name="PPEXSIRetrieveComplexEDM")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    complex(c_double), intent(out)   :: EDMnzvalLocal(*)
    real(c_double), intent(out)   :: totalEnergyS
    integer(c_int), intent(out)   :: info
  end subroutine 

  subroutine f_ppexsi_retrieve_complex_fdm(&
      plan,&
      FDMnzvalLocal,&
      totalFreeEnergy,&
      info) &
      bind(C, Name="PPEXSIRetrieveComplexFDM")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_intptr_t), intent(in), value :: plan
    complex(c_double), intent(out)   :: FDMnzvalLocal(*)
    real(c_double), intent(out)   :: totalFreeEnergy
    integer(c_int), intent(out)   :: info
  end subroutine 

end interface

end module f_ppexsi_interface
