/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.

Author: Lin Lin

This file is part of PEXSI. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.
*/
/// @file interface.cpp
/// @brief Interface subroutines of PPEXSI that can be called by both C and FORTRAN.
///
/// This file will eventually merge with interface.cpp.
/// @date Original:      2013-02-03
/// @date Revision:      2014-01-01
/// @date Revision:      2014-03-07  Second generation interface.
#include "c_pexsi_interface.h"
#include "ppexsi.hpp"
#include "pexsi/blas.hpp"

// FIXME
// Error handling used in the C interface that is different from the
// throw/catch system.
#define iC(fun)  { int ierr=fun; if(ierr!=0) exit(1); }
#define iA(expr) { if((expr)==0) { std::cerr<<"wrong "<<__LINE__<<" in " <<__FILE__<<std::endl; std::cerr.flush(); exit(1); } }

using namespace PEXSI;

extern "C"
void ReadDistSparseMatrixFormattedHeadInterface (
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    MPI_Comm comm )
{
  int mpirank;  MPI_Comm_rank(comm, &mpirank);
  int mpisize;  MPI_Comm_size(comm, &mpisize);
  std::ifstream fin;
  if( mpirank == 0 ){
    fin.open(filename);
    if( !fin.good() ){
      ErrorHandling( "File cannot be openeded!" );
    }
    Int dummy;
    fin >> *size >> dummy;
    fin >> *nnz >> dummy;
  }

  MPI_Bcast( &(*size), 1, MPI_INT, 0, comm);
  MPI_Bcast( &(*nnz),  1, MPI_INT, 0, comm);

  IntNumVec  colptr(*size+1);
  if( mpirank == 0 ){
    Int* ptr = colptr.Data();
    for( Int i = 0; i < *size+1; i++ )
      fin >> *(ptr++);
  }

  MPI_Bcast(colptr.Data(), *size+1, MPI_INT, 0, comm);

  // Compute the number of columns on each processor
  IntNumVec numColLocalVec(mpisize);
  Int numColFirst;
  numColFirst = *size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = *size - numColFirst * (mpisize-1);
  // Modify the last entry

  *numColLocal = numColLocalVec[mpirank];

  *nnzLocal = colptr[mpirank * numColFirst + (*numColLocal)] -
    colptr[mpirank * numColFirst];

  // Close the file
  if( mpirank == 0 ){
    fin.close();
  }

  return;
}
// -----  end of function ReadDistSparseMatrixFormattedHeadInterface

extern "C"
void ReadDistSparseMatrixFormattedInterface(
    char*     filename,
    int       size,
    int       nnz,
    int       nnzLocal,
    int       numColLocal,
    int*      colptrLocal,
    int*      rowindLocal,
    double*   nzvalLocal,
    MPI_Comm  comm )
{
  DistSparseMatrix<Real> A;
  ReadDistSparseMatrixFormatted( filename, A, comm );
  iA( size == A.size );
  iA( nnz  == A.nnz  );
  iA( nnzLocal == A.nnzLocal );
  iA( numColLocal + 1 == A.colptrLocal.m() );

  blas::Copy( numColLocal+1, A.colptrLocal.Data(), 1,
      colptrLocal, 1 );

  blas::Copy( nnzLocal, A.rowindLocal.Data(), 1,
      rowindLocal, 1 );

  blas::Copy( nnzLocal, A.nzvalLocal.Data(), 1,
      nzvalLocal, 1 );

  return;
}
// -----  end of function ReadDistSparseMatrixFormattedInterface

extern "C"
void ReadDistSparseMatrixHeadInterface (
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    MPI_Comm comm )
{
  int mpirank;  MPI_Comm_rank(comm, &mpirank);
  int mpisize;  MPI_Comm_size(comm, &mpisize);
  std::ifstream fin;
  if( mpirank == 0 ){
    fin.open(filename);
    if( !fin.good() ){
      ErrorHandling( "File cannot be openeded!" );
    }
    fin.read((char*)size, sizeof(int));
    fin.read((char*)nnz,  sizeof(int));
  }

  MPI_Bcast( &(*size), 1, MPI_INT, 0, comm);
  MPI_Bcast( &(*nnz),  1, MPI_INT, 0, comm);

  IntNumVec  colptr(*size+1);
  if( mpirank == 0 ){
    Int tmp;
    fin.read((char*)&tmp, sizeof(int));

    if( tmp != (*size)+1 ){
      ErrorHandling( "colptr is not of the right size." );
    }

    Int* ptr = colptr.Data();
    fin.read((char*)ptr, (*size+1) * sizeof(int));
  }

  MPI_Bcast(colptr.Data(), *size+1, MPI_INT, 0, comm);

  // Compute the number of columns on each processor
  IntNumVec numColLocalVec(mpisize);
  Int numColFirst;
  numColFirst = *size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = *size - numColFirst * (mpisize-1);
  // Modify the last entry

  *numColLocal = numColLocalVec[mpirank];

  *nnzLocal = colptr[mpirank * numColFirst + (*numColLocal)] -
    colptr[mpirank * numColFirst];

  // Close the file
  if( mpirank == 0 ){
    fin.close();
  }

  return;
}
// -----  end of function ReadDistSparseMatrixHeadInterface

extern "C"
void ParaReadDistSparseMatrixInterface(
    char*     filename,
    int       size,
    int       nnz,
    int       nnzLocal,
    int       numColLocal,
    int*      colptrLocal,
    int*      rowindLocal,
    double*   nzvalLocal,
    MPI_Comm  comm )
{
  DistSparseMatrix<Real> A;
  ParaReadDistSparseMatrix( filename, A, comm );
  iA( size == A.size );
  iA( nnz  == A.nnz  );
  iA( nnzLocal == A.nnzLocal );
  iA( numColLocal + 1 == A.colptrLocal.m() );

  blas::Copy( numColLocal+1, A.colptrLocal.Data(), 1,
      colptrLocal, 1 );

  blas::Copy( nnzLocal, A.rowindLocal.Data(), 1,
      rowindLocal, 1 );

  blas::Copy( nnzLocal, A.nzvalLocal.Data(), 1,
      nzvalLocal, 1 );

  return;
}
// -----  end of function ReadDistSparseMatrixFormattedInterface

// *********************************************************************
// Second interface
// *********************************************************************

extern "C"
void PPEXSISetDefaultOptions(
    PPEXSIOptions*   options ){
  options->spin                  = 2.0;
  options->temperature           = 0.0019;   // 300K
  options->gap                   = 0.0;      // no gap
  options->deltaE                = 10.0;
  options->numPole               = 20;
  options->isInertiaCount        = 1;
  options->maxPEXSIIter          = 3;
  options->muMin0                = -10.0;
  options->muMax0                = +10.0;
  options->mu0                   =  0.0;  // For more control by the user
  options->muInertiaTolerance    = 0.05;
  options->muInertiaExpansion    = 0.3;
  options->muPEXSISafeGuard      = 0.05;
  options->numElectronPEXSITolerance = 1.0E-10;
  options->matrixType            = 0;
  options->isSymbolicFactorize   = 1;
  options->solver                = 0;
  options->symmetricStorage      = 0;
  options->ordering              = 0;
  options->rowOrdering           = 0;
  options->npSymbFact            = 1;
  options->symmetric             = 1;
  options->transpose             = 0;
  options->method                = 2;
  options->verbosity             = 1;
  options->nPoints               = 2;
}   // -----  end of function PPEXSISetDefaultOptions  -----

extern "C"
PPEXSIPlan PPEXSIPlanInitialize(
    MPI_Comm      comm,
    int           numProcRow,
    int           numProcCol,
    int           outputFileIndex,
    int*          info ){

  int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );

  *info = 0;
  PPEXSIData *ptrData;

  try{
    ptrData = new PPEXSIData( comm, numProcRow, numProcCol, outputFileIndex );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return reinterpret_cast<PPEXSIPlan>(ptrData);
}   // -----  end of function PPEXSIPlanInitialize  -----

extern "C"
void PPEXSILoadRealHSMatrix(
    PPEXSIPlan    plan,
    PPEXSIOptions options,
    int           nrows,
    int           nnz,
    int           nnzLocal,
    int           numColLocal,
    int*          colptrLocal,
    int*          rowindLocal,
    double*       HnzvalLocal,
    int           isSIdentity,
    double*       SnzvalLocal,
    int*          info ){

  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  *info = 0;

  try{
    reinterpret_cast<PPEXSIData*>(plan)->
      LoadRealMatrix(
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          HnzvalLocal,
          isSIdentity,
          SnzvalLocal,
          options.solver,
          options.verbosity );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return;
}   // -----  end of function PPEXSILoadRealHSMatrix  -----

extern "C"
void PPEXSILoadComplexHSMatrix(
    PPEXSIPlan    plan,
    PPEXSIOptions options,
    int           nrows,
    int           nnz,
    int           nnzLocal,
    int           numColLocal,
    int*          colptrLocal,
    int*          rowindLocal,
    double*       HnzvalLocal,
    int           isSIdentity,
    double*       SnzvalLocal,
    int*          info ){

  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  *info = 0;

  try{
    reinterpret_cast<PPEXSIData*>(plan)->
      LoadComplexMatrix(
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          reinterpret_cast<Complex*>(HnzvalLocal),
          isSIdentity,
          reinterpret_cast<Complex*>(SnzvalLocal),
          options.solver,
          options.verbosity );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return;
}   // -----  end of function PPEXSILoadComplexHSMatrix  -----

extern "C"
void PPEXSISymbolicFactorizeRealSymmetricMatrix(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int*              info ) {
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  *info = 0;

  try{
    std::string colPerm;
    switch (options.solver){
      case 0:
        {
          //Handle SuperLU ordering options
          switch (options.ordering){
            case 0:
              colPerm = "PARMETIS";
              break;
            case 1:
              colPerm = "METIS_AT_PLUS_A";
              break;
            case 2:
              colPerm = "MMD_AT_PLUS_A";
              break;
            case 3:
              colPerm = "NATURAL";
              break;
            default:
              ErrorHandling("Unsupported ordering strategy.");
              break;
          }
        }
        break;
#ifdef WITH_SYMPACK
      case 1:
        {
          //Handle symPACK ordering options
          switch (options.ordering){
            case 5:
              colPerm = "PTSCOTCH";
              break;
            case 6:
              colPerm = "SCOTCH";
              break;
            case 2:
              colPerm = "MMD";
              break;
            case 3:
              colPerm = "NATURAL";
              break;
            case 4:
              colPerm = "AMD";
              break;
            case 0:
              colPerm = "PARMETIS";
              break;
            case 1:
              colPerm = "METIS";
              break;
            default:
              ErrorHandling("Unsupported ordering strategy.");
              break;
          }
        }
        break;
#endif
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }

    reinterpret_cast<PPEXSIData*>(plan)->
      SymbolicFactorizeRealSymmetricMatrix(
          options.solver,
          options.symmetricStorage,
          colPerm,
          options.npSymbFact,
          options.verbosity );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSISymbolicFactorizeRealSymmetricMatrix  -----

extern "C"
void PPEXSISymbolicFactorizeRealUnsymmetricMatrix(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    int*              info ) {
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  *info = 0;

  try{
    std::string colPerm;
    switch (options.ordering){
      case 0:
        colPerm = "PARMETIS";
        break;
      case 1:
        colPerm = "METIS_AT_PLUS_A";
        break;
      case 2:
        colPerm = "MMD_AT_PLUS_A";
        break;
      case 3:
        colPerm = "NATURAL";
        break;
      default:
        ErrorHandling("Unsupported ordering strategy.");
    }

    std::string rowPerm;
    switch (options.rowOrdering){
      case 0:
        rowPerm = "NOROWPERM";
        break;
      case 1:
        rowPerm = "LargeDiag";
        break;
      default:
        ErrorHandling("Unsupported row ordering strategy.");
    }

    reinterpret_cast<PPEXSIData*>(plan)->
      SymbolicFactorizeRealUnsymmetricMatrix(
          options.solver,
          colPerm,
          rowPerm,
          options.npSymbFact,
          options.transpose,
          AnzvalLocal,
          options.verbosity );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSISymbolicFactorizeRealUnsymmetricMatrix  -----

extern "C"
void PPEXSISymbolicFactorizeComplexSymmetricMatrix(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int*              info ) {
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  *info = 0;

  try{
    std::string colPerm;
    switch (options.solver){
      case 0:
        {
          //Handle SuperLU ordering options
          switch (options.ordering){
            case 0:
              colPerm = "PARMETIS";
              break;
            case 1:
              colPerm = "METIS_AT_PLUS_A";
              break;
            case 2:
              colPerm = "MMD_AT_PLUS_A";
              break;
            case 3:
              colPerm = "NATURAL";
              break;
            default:
              ErrorHandling("Unsupported ordering strategy.");
              break;
          }
        }
        break;
#ifdef WITH_SYMPACK
      case 1:
        {
          //Handle symPACK ordering options
          switch (options.ordering){
            case 5:
              colPerm = "PTSCOTCH";
              break;
            case 6:
              colPerm = "SCOTCH";
              break;
            case 2:
              colPerm = "MMD";
              break;
            case 3:
              colPerm = "NATURAL";
              break;
            case 4:
              colPerm = "AMD";
              break;
            case 0:
              colPerm = "PARMETIS";
              break;
            case 1:
              colPerm = "METIS";
              break;
            default:
              ErrorHandling("Unsupported ordering strategy.");
              break;
          }
        }
        break;
#endif
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }

    reinterpret_cast<PPEXSIData*>(plan)->
      SymbolicFactorizeComplexSymmetricMatrix(
          options.solver,
          options.symmetricStorage,
          colPerm,
          options.npSymbFact,
          options.verbosity );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSISymbolicFactorizeRealSymmetricMatrix  -----

extern "C"
void PPEXSISymbolicFactorizeComplexUnsymmetricMatrix(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    int*              info ) {
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  *info = 0;

  try{
    std::string colPerm;
    switch (options.ordering){
      case 0:
        colPerm = "PARMETIS";
        break;
      case 1:
        colPerm = "METIS_AT_PLUS_A";
        break;
      case 2:
        colPerm = "MMD_AT_PLUS_A";
        break;
      case 3:
        colPerm = "NATURAL";
        break;
      default:
        ErrorHandling("Unsupported ordering strategy.");
    }

    std::string rowPerm;
    switch (options.rowOrdering){
      case 0:
        rowPerm = "NOROWPERM";
        break;
      case 1:
        rowPerm = "LargeDiag";
        break;
      default:
        ErrorHandling("Unsupported row ordering strategy.");
    }

    reinterpret_cast<PPEXSIData*>(plan)->
      SymbolicFactorizeComplexUnsymmetricMatrix(
          options.solver,
          colPerm,
          rowPerm,
          options.npSymbFact,
          options.transpose,
          AnzvalLocal,
          options.verbosity );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSISymbolicFactorizeRealUnsymmetricMatrix  -----

extern "C"
void PPEXSIInertiaCountRealMatrix(
    /* Input parameters */
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int               numShift,
    double*           shiftList,
    /* Output parameters */
    double*           inertiaList,
    int*              info ) {

  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  *info = 0;

  try{
    std::vector<Real>    shiftVec(numShift);
    std::vector<Real>    inertiaVec(numShift);
    for( Int i = 0; i < numShift; i++ ){
      shiftVec[i] = shiftList[i];
    }

    reinterpret_cast<PPEXSIData*>(plan)->
      CalculateNegativeInertiaReal(
          shiftVec,
          inertiaVec,
          options.solver,
          options.verbosity );

    for( Int i = 0; i < numShift; i++ ){
      inertiaList[i] = inertiaVec[i];
    }
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSIInertiaCountRealMatrix  -----

extern "C"
void PPEXSIInertiaCountComplexMatrix(
    /* Input parameters */
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int               numShift,
    double*           shiftList,
    /* Output parameters */
    double*           inertiaList,
    int*              info ) {

  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  *info = 0;

  try{
    std::vector<Real>    shiftVec(numShift);
    std::vector<Real>    inertiaVec(numShift);
    for( Int i = 0; i < numShift; i++ ){
      shiftVec[i] = shiftList[i];
    }

    reinterpret_cast<PPEXSIData*>(plan)->
      CalculateNegativeInertiaComplex(
          shiftVec,
          inertiaVec,
          options.solver,
          options.verbosity );

    for( Int i = 0; i < numShift; i++ ){
      inertiaList[i] = inertiaVec[i];
    }
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSIInertiaCountComplexMatrix  -----

extern "C"
void PPEXSICalculateFermiOperatorReal(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double            mu,
    double            numElectronExact,
    double*           numElectronPEXSI,
    double*           numElectronDrvMuPEXSI,
    int*              info )
{
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->CalculateFermiOperatorReal(
        options.numPole,
        options.temperature,
        options.gap,
        options.deltaE,
        mu,
        numElectronExact,
        options.numElectronPEXSITolerance,
        options.solver,
        options.verbosity,
        *numElectronPEXSI,
        *numElectronDrvMuPEXSI );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSICalculateFermiOperatorReal  -----

extern "C"
void PPEXSICalculateFermiOperatorReal3(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double            numElectronExact,
    double*           mu,
    double*           numElectronPEXSI,
    int*              info )
{
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->CalculateFermiOperatorReal3(
        options.numPole,
        options.temperature,
        options.gap,
        options.deltaE,
        numElectronExact,
        options.numElectronPEXSITolerance,
        options.solver,
        options.verbosity,
        *mu,
        *numElectronPEXSI,
        options.method,
        options.nPoints,
        options.spin);
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSICalculateFermiOperatorReal3  -----

extern "C"
void PPEXSICalculateFermiOperatorComplex(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double            mu,
    double            numElectronExact,
    double*           numElectronPEXSI,
    double*           numElectronDrvMuPEXSI,
    int*              info )
{
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->CalculateFermiOperatorComplex(
        options.numPole,
        options.temperature,
        options.gap,
        options.deltaE,
        mu,
        numElectronExact,
        options.numElectronPEXSITolerance,
        options.solver,
        options.verbosity,
        *numElectronPEXSI,
        *numElectronDrvMuPEXSI,
        options.method,
        options.nPoints,
        options.spin);
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSICalculateFermiOperatorComplex   -----

extern "C"
void PPEXSICalculateEDMCorrectionReal(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int*              info )
{
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->CalculateEDMCorrectionReal(
        options.numPole,
        options.solver,
        options.verbosity,
        options.nPoints,
        options.spin);
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSICalculateEDMCorrectionReal  -----

extern "C"
void PPEXSICalculateEDMCorrectionComplex(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int*              info )
{
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->CalculateEDMCorrectionComplex(
        options.numPole,
        options.solver,
        options.verbosity,
        options.nPoints,
        options.spin);
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSICalculateEDMCorrectionComplex  -----

extern "C"
void PPEXSIInterpolateDMReal(
    PPEXSIPlan        plan,
    PPEXSIOptions*    options,
    double            numElectronExact,
    double            numElectronPEXSI,
    double *          NeVec,
    double *          muPEXSI,
    int*              info )
{
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->InterpolateDMReal(
        numElectronExact,
        numElectronPEXSI,
        options->numElectronPEXSITolerance,
        options->nPoints,
        NeVec,
        options->muMin0,
        options->muMax0,
        *muPEXSI,
        options->method,
        options->verbosity);

  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PEXSIInterpolateDMReal  -----

extern "C"
void PPEXSIInterpolateDMComplex(
    PPEXSIPlan        plan,
    PPEXSIOptions*    options,
    double            numElectronExact,
    double            numElectronPEXSI,
    double *          NeVec,
    double *          muPEXSI,
    int*              info )
{
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->InterpolateDMComplex(
        numElectronExact,
        numElectronPEXSI,
        options->numElectronPEXSITolerance,
        options->nPoints,
        NeVec,
        options->muMin0,
        options->muMax0,
        *muPEXSI,
        options->method,
        options->verbosity);

  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PEXSIInterpolateDMComplex  -----

extern "C"
void PPEXSISelInvRealSymmetricMatrix (
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    double*           AinvnzvalLocal,
    int*              info )
{
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->SelInvRealSymmetricMatrix(
        options.solver,
        options.symmetricStorage,
        AnzvalLocal,
        options.verbosity,
        AinvnzvalLocal );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSISelInvRealSymmetricMatrix  -----

extern "C"
void PPEXSISelInvRealUnsymmetricMatrix (
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    double*           AinvnzvalLocal,
    int*              info )
{
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->SelInvRealUnsymmetricMatrix(
        options.solver,
        AnzvalLocal,
        options.verbosity,
        AinvnzvalLocal );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSISelInvRealUnsymmetricMatrix  -----

extern "C"
void PPEXSISelInvComplexSymmetricMatrix (
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    double*           AinvnzvalLocal,
    int*              info )
{
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->SelInvComplexSymmetricMatrix(
        options.solver,
        options.symmetricStorage,
        AnzvalLocal,
        options.verbosity,
        AinvnzvalLocal );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSISelInvComplexSymmetricMatrix  -----

extern "C"
void PPEXSISelInvComplexUnsymmetricMatrix (
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    double*           AinvnzvalLocal,
    int*              info )
{
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->SelInvComplexUnsymmetricMatrix(
        options.solver,
        AnzvalLocal,
        options.verbosity,
        AinvnzvalLocal );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }

  return ;
}		// -----  end of function PPEXSISelInvComplexUnsymmetricMatrix  -----

extern "C"
void PPEXSIDFTDriver(
    /* Input parameters */
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double            numElectronExact,
    /* Output parameters */
    double*           muPEXSI,
    double*           numElectronPEXSI,
    double*           muMinInertia,
    double*           muMaxInertia,
    int*              numTotalInertiaIter,
    int*              numTotalPEXSIIter,
    int*              info ){
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->DFTDriver(
        numElectronExact,
        options.temperature,
        options.gap,
        options.deltaE,
        options.numPole,
        options.isInertiaCount,
        options.maxPEXSIIter,
        options.muMin0,
        options.muMax0,
        options.mu0,
        options.muInertiaTolerance,
        options.muInertiaExpansion,
        options.muPEXSISafeGuard,
        options.numElectronPEXSITolerance,
        options.matrixType,
        options.isSymbolicFactorize,
        options.solver,
        options.symmetricStorage,
        options.ordering,
        options.npSymbFact,
        options.verbosity,
        *muPEXSI,
        *numElectronPEXSI,
        *muMinInertia,
        *muMaxInertia,
        *numTotalInertiaIter,
        *numTotalPEXSIIter );
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIDFTDriver  -----

extern "C"
void PPEXSIDFTDriver2(
    /* Input parameters */
    PPEXSIPlan        plan,
    PPEXSIOptions*    options,        // options is input and output
    double            numElectronExact,
    //int               method,
    //int               nPoints,
    /* Output parameters */
    double*           muPEXSI,
    double*           numElectronPEXSI,
    int*              numTotalInertiaIter,
    int*              info ){
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    reinterpret_cast<PPEXSIData*>(plan)->DFTDriver2(
        numElectronExact,
        options->temperature,
        options->gap,
        options->deltaE,
        options->numPole,
        options->muInertiaTolerance,
        options->numElectronPEXSITolerance,
        options->matrixType,
        options->isSymbolicFactorize,
        options->solver,
        options->symmetricStorage,
        options->ordering,
        options->npSymbFact,
        options->verbosity,
        *muPEXSI,
        *numElectronPEXSI,
        options->muMin0,
        options->muMax0,
        *numTotalInertiaIter,
        options->method,
        options->nPoints,
        options->spin);
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIDFTDriver2  -----

extern "C"
void PPEXSIRetrieveRealDFTMatrix(
    PPEXSIPlan        plan,
    double*      DMnzvalLocal,
    double*     EDMnzvalLocal,
    double*     FDMnzvalLocal,
    double*     totalEnergyH,
    double*     totalEnergyS,
    double*     totalFreeEnergy,
    int*              info ){
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();
  PPEXSIData* ptrData = reinterpret_cast<PPEXSIData*>(plan);

  try{
    Int nnzLocal = ptrData->RhoRealMat().nnzLocal;

    blas::Copy( nnzLocal, ptrData->RhoRealMat().nzvalLocal.Data(), 1,
        DMnzvalLocal, 1 );

    blas::Copy( nnzLocal, ptrData->EnergyDensityRealMat().nzvalLocal.Data(), 1,
        EDMnzvalLocal, 1 );

    blas::Copy( nnzLocal, ptrData->FreeEnergyDensityRealMat().nzvalLocal.Data(), 1,
        FDMnzvalLocal, 1 );

    *totalEnergyH = ptrData->TotalEnergyH();

    *totalEnergyS = ptrData->TotalEnergyS();

    *totalFreeEnergy = ptrData->TotalFreeEnergy();
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIRetrieveRealDFTMatrix  -----

extern "C"
void PPEXSIRetrieveComplexDFTMatrix(
    PPEXSIPlan        plan,
    double*      DMnzvalLocal,
    double*     EDMnzvalLocal,
    double*     FDMnzvalLocal,
    double*     totalEnergyH,
    double*     totalEnergyS,
    double*     totalFreeEnergy,
    int*              info ){
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();
  PPEXSIData* ptrData = reinterpret_cast<PPEXSIData*>(plan);

  try{
    Int nnzLocal = ptrData->RhoComplexMat().nnzLocal;

    blas::Copy( 2*nnzLocal,
        reinterpret_cast<double*>(ptrData->RhoComplexMat().nzvalLocal.Data()), 1,
        DMnzvalLocal, 1 );

    blas::Copy( 2*nnzLocal,
        reinterpret_cast<double*>(ptrData->EnergyDensityComplexMat().nzvalLocal.Data()), 1,
        EDMnzvalLocal, 1 );

    blas::Copy( 2*nnzLocal,
        reinterpret_cast<double*>(ptrData->FreeEnergyDensityComplexMat().nzvalLocal.Data()), 1,
        FDMnzvalLocal, 1 );

    *totalEnergyH = ptrData->TotalEnergyH();

    *totalEnergyS = ptrData->TotalEnergyS();

    *totalFreeEnergy = ptrData->TotalFreeEnergy();
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIRetrieveComplexDFTMatrix  -----

extern "C"
void PPEXSIPlanFinalize(
    PPEXSIPlan    plan,
    int*          info ){

  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();

  try{
    delete reinterpret_cast<PPEXSIData*>(plan);
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIPlanFinalize  -----

extern "C"
void PPEXSIGetPoleDM(
    double*       zshift,
    double*       zweight,
    int           Npole,
    double        temp,
    double        gap,
    double        deltaE,
    double        mu ){

  CpxNumVec zshiftVec(Npole), zweightVec(Npole);

  GetPoleDensity( zshiftVec.Data(), zweightVec.Data(),
      Npole, temp, gap, deltaE, mu );

  for( Int i = 0; i < Npole; i++ ){
    zshift[2*i]    = zshiftVec[i].real();
    zshift[2*i+1]  = zshiftVec[i].imag();
    zweight[2*i]   = zweightVec[i].real();
    zweight[2*i+1] = zweightVec[i].imag();
  }

  return;
}   // -----  end of function PPEXSIGetPoleDM  -----

extern "C"
void PPEXSIGetPoleEDM(
    double*       zshift,
    double*       zweight,
    int           Npole,
    double        temp,
    double        gap,
    double        deltaE,
    double        mu ){

  CpxNumVec zshiftVec(Npole), zweightVec(Npole);

  GetPoleForce( zshiftVec.Data(), zweightVec.Data(),
      Npole, temp, gap, deltaE, mu );

  for( Int i = 0; i < Npole; i++ ){
    zshift[2*i]    = zshiftVec[i].real();
    zshift[2*i+1]  = zshiftVec[i].imag();
    zweight[2*i]   = zweightVec[i].real();
    zweight[2*i+1] = zweightVec[i].imag();
  }

  return;
}   // -----  end of function PPEXSIGetPoleEDM  -----

extern "C"
void PPEXSIGetPoleFDM(
    double*       zshift,
    double*       zweight,
    int           Npole,
    double        temp,
    double        gap,
    double        deltaE,
    double        mu ){

  CpxNumVec zshiftVec(Npole), zweightVec(Npole);

  GetPoleHelmholtz( zshiftVec.Data(), zweightVec.Data(),
      Npole, temp, gap, deltaE, mu );

  for( Int i = 0; i < Npole; i++ ){
    zshift[2*i]    = zshiftVec[i].real();
    zshift[2*i+1]  = zshiftVec[i].imag();
    zweight[2*i]   = zweightVec[i].real();
    zweight[2*i+1] = zweightVec[i].imag();
  }

  return;
}   // -----  end of function PPEXSIGetPoleEDM  -----

// *********************************************************************
// FORTRAN interface
//
// NOTE:
//
// 1. Most of the interface routines do not explicitly depend on the MPI
// communicators, and therefore can be called directly through the
// ISO_C_BINDING feature as in f_interface.f90.
//
// 2. For routines use MPI communicators, the communicators from FORTRAN
// must be transferred to C communicators using f2c_comm, and the
// FORTRAN version is given below.
//
// 3. In ISO_C_BINDING, passing by value is allowed and it is not
// required that all arguments by given in the form of pointers.
//
// 4. The following routines need not be defined in the C header file.
// *********************************************************************

/// @brief Internal subroutine to convert FORTRAN communicator to C
extern "C"
MPI_Comm f2c_comm(MPI_Fint* Fcomm)
{
  return MPI_Comm_f2c((*Fcomm));
}  // -----  end of function f2c_comm -----

/// @brief FORTRAN interface for @ref ReadDistSparseMatrixFormattedHeadInterface.
extern "C"
void f_read_distsparsematrix_formatted_head (
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    MPI_Fint* Fcomm )
{
  ReadDistSparseMatrixFormattedHeadInterface(
      filename,
      size,
      nnz,
      nnzLocal,
      numColLocal,
      f2c_comm( Fcomm ) );

  return;
}  // -----  end of function f_read_distsparsematrix_formatted_head

/// @brief FORTRAN interface for @ref ReadDistSparseMatrixFormattedInterface.
extern "C"
void f_read_distsparsematrix_formatted (
    char*    filename,
    int      size,
    int      nnz,
    int      nnzLocal,
    int      numColLocal,
    int*     colptrLocal,
    int*     rowindLocal,
    double*  nzvalLocal,
    MPI_Fint* Fcomm )
{
  ReadDistSparseMatrixFormattedInterface(
      filename,
      size,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      nzvalLocal,
      f2c_comm( Fcomm ) );
  return;
} // -----  end of function f_read_distsparsematrix_formatted  -----

/// @brief FORTRAN interface for @ref ReadDistSparseMatrixHeadInterface.
extern "C"
void f_read_distsparsematrix_head (
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    MPI_Fint* Fcomm )
{
  ReadDistSparseMatrixHeadInterface(
      filename,
      size,
      nnz,
      nnzLocal,
      numColLocal,
      f2c_comm( Fcomm ) );

  return;
}  // -----  end of function f_read_distsparsematrix_head

/// @brief FORTRAN interface for @ref ParaReadDistSparseMatrixInterface.
extern "C"
void f_para_read_distsparsematrix (
    char*    filename,
    int      size,
    int      nnz,
    int      nnzLocal,
    int      numColLocal,
    int*     colptrLocal,
    int*     rowindLocal,
    double*  nzvalLocal,
    MPI_Fint* Fcomm )
{
  ParaReadDistSparseMatrixInterface(
      filename,
      size,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      nzvalLocal,
      f2c_comm( Fcomm ) );
  return;
} // -----  end of function f_para_read_distsparsematrix  -----

/// @brief FORTRAN interface for @ref PPEXSIPlanInitialize
extern "C"
PPEXSIPlan f_ppexsi_plan_initialize (
    MPI_Fint*     Fcomm,
    int           numProcRow,
    int           numProcCol,
    int           outputFileIndex,
    int*          info ){
  return PPEXSIPlanInitialize(
      f2c_comm( Fcomm ),
      numProcRow,
      numProcCol,
      outputFileIndex,
      info );
}		// -----  end of function f_ppexsi_plan_initialize  -----

extern "C"
void PPEXSIRetrieveComplexDMEDM(
    PPEXSIPlan  plan,
    double*     DMnzvalLocal,
    double*     EDMnzvalLocal,
    int*        nnz,
    int*        info ){
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();
  PPEXSIData* ptrData = reinterpret_cast<PPEXSIData*>(plan);

  try{
    Int nnzLocal = ptrData->RhoComplexMat().nnzLocal;

    blas::Copy( 2*nnzLocal,
        reinterpret_cast<double*>(ptrData->RhoComplexMat().nzvalLocal.Data()), 1,
        DMnzvalLocal, 1 );

    blas::Copy( 2*nnzLocal,
        reinterpret_cast<double*>(ptrData->EnergyDensityComplexMat().nzvalLocal.Data()), 1,
        EDMnzvalLocal, 1 );

    *nnz = nnzLocal;
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIRetrieveComplexDMEDM   -----

extern "C"
void PPEXSIRetrieveRealDMEDM(
    PPEXSIPlan  plan,
    double*     DMnzvalLocal,
    double*     EDMnzvalLocal,
    int*        nnz,
    int*        info ){
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();
  PPEXSIData* ptrData = reinterpret_cast<PPEXSIData*>(plan);

  try{

    Int nnzLocal = ptrData->RhoRealMat().nnzLocal;

    blas::Copy( nnzLocal, ptrData->RhoRealMat().nzvalLocal.Data(), 1,
        DMnzvalLocal, 1 );

    blas::Copy( nnzLocal, ptrData->EnergyDensityRealMat().nzvalLocal.Data(), 1,
        EDMnzvalLocal, 1 );

    *nnz = nnzLocal;
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIRetrieveRealDMEDM   -----

extern "C"
void PPEXSIRetrieveRealDM(
    PPEXSIPlan  plan,
    double*     DMnzvalLocal,
    double*     totalEnergyH,
    int*              info ){
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();
  PPEXSIData* ptrData = reinterpret_cast<PPEXSIData*>(plan);

  try{
    Int nnzLocal = ptrData->RhoRealMat().nnzLocal;

    blas::Copy( nnzLocal, ptrData->RhoRealMat().nzvalLocal.Data(), 1,
        DMnzvalLocal, 1 );
    *totalEnergyH = ptrData->TotalEnergyH();
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIRetrieveRealDM  -----

extern "C"
void PPEXSIRetrieveRealEDM(
    PPEXSIPlan  plan,
    PPEXSIOptions     options,
    double*     EDMnzvalLocal,
    double*     totalEnergyS,
    int*        info ){
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();
  PPEXSIData* ptrData = reinterpret_cast<PPEXSIData*>(plan);

  try{

    //FIXME method == 2
    if( options.method == 2){

      statusOFS << std::endl << " Warning, EDM correction before retrieve"
        << gridPole->mpirank << std::endl;
      reinterpret_cast<PPEXSIData*>(plan)->CalculateEDMCorrectionReal(
          options.numPole,
          options.solver,
          options.verbosity,
          options.nPoints,
          options.spin);
    }

    Int nnzLocal = ptrData->EnergyDensityRealMat().nnzLocal;

    blas::Copy( nnzLocal, ptrData->EnergyDensityRealMat().nzvalLocal.Data(), 1,
        EDMnzvalLocal, 1 );

    *totalEnergyS = ptrData->TotalEnergyS();
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIRetrieveRealEDM   -----

extern "C"
void PPEXSIRetrieveRealFDM(
    PPEXSIPlan  plan,
    PPEXSIOptions     options,
    double*     FDMnzvalLocal,
    double*     totalEnergyF,
    int*              info ){
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();
  PPEXSIData* ptrData = reinterpret_cast<PPEXSIData*>(plan);

  try{
    Int nnzLocal = ptrData->FreeEnergyDensityRealMat().nnzLocal;

    if( options.method == 2){
      blas::Copy( nnzLocal, ptrData->RhoRealMat().nzvalLocal.Data(), 1,
          FDMnzvalLocal, 1 );
    }
    else{
      blas::Copy( nnzLocal, ptrData->FreeEnergyDensityRealMat().nzvalLocal.Data(), 1,
          FDMnzvalLocal, 1 );
    }
    *totalEnergyF = ptrData->TotalFreeEnergy();
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIRetrieveRealFDM  -----

extern "C"
void PPEXSIRetrieveComplexDM(
    PPEXSIPlan        plan,
    double*      DMnzvalLocal,
    double*     totalEnergyH,
    int*              info ){
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();
  PPEXSIData* ptrData = reinterpret_cast<PPEXSIData*>(plan);

  try{
    Int nnzLocal = ptrData->RhoComplexMat().nnzLocal;

    blas::Copy( 2*nnzLocal,
        reinterpret_cast<double*>(ptrData->RhoComplexMat().nzvalLocal.Data()), 1,
        DMnzvalLocal, 1 );

    *totalEnergyH = ptrData->TotalEnergyH();
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIRetrieveComplexDM  -----

extern "C"
void PPEXSIRetrieveComplexEDM(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*     EDMnzvalLocal,
    double*     totalEnergyS,
    int*              info ){
  *info = 0;
  const GridType* gridPole =
    reinterpret_cast<PPEXSIData*>(plan)->GridPole();
  PPEXSIData* ptrData = reinterpret_cast<PPEXSIData*>(plan);

  try{
    //FIXME method == 2
    if( options.method == 2){

      statusOFS << std::endl << " Warning, EDM correction before retrieve"
        << gridPole->mpirank << std::endl;
      reinterpret_cast<PPEXSIData*>(plan)->CalculateEDMCorrectionComplex(
          options.numPole,
          options.solver,
          options.verbosity,
          options.nPoints,
          options.spin);
    }

    Int nnzLocal = ptrData->RhoComplexMat().nnzLocal;

    blas::Copy( 2*nnzLocal,
        reinterpret_cast<double*>(ptrData->EnergyDensityComplexMat().nzvalLocal.Data()), 1,
        EDMnzvalLocal, 1 );

    *totalEnergyS = ptrData->TotalEnergyS();
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank
      << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    *info = 1;
  }
  return;
}   // -----  end of function PPEXSIRetrieveComplexEDM -----

