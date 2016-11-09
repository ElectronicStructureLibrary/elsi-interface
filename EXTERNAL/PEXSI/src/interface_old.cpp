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
/// @file interface_old.cpp
/// @brief Interface subroutines of PPEXSI that can be called by both C and FORTRAN.
/// @date 2013-02-03
/// @date Last revision: 2014-01-01
#include "pexsi/c_pexsi_interface.h"
#include "ppexsi.hpp"
#include "pexsi/blas.hpp"

// FIXME
// Error handling used in the C interface that is different from the
// throw/catch system.
#define iC(fun)  { int ierr=fun; if(ierr!=0) exit(1); }
#define iA(expr) { if((expr)==0) { std::cerr<<"wrong "<<__LINE__<<" in " <<__FILE__<<std::endl; std::cerr.flush(); exit(1); } }

using namespace PEXSI;




// Derivative of the Fermi-Dirac distribution with respect to mu 
// (therefore there is no minus sign)
// fdDrvMu(x) = beta * exp(beta*x)/(1+exp(beta*z))^2
// Note: compared to fdDrvMu used in pole.cpp, the factor 2.0 is not
// present, because it is included in inertiaVec. 
inline Real fdDrvMu( Real x, double beta, double mu ){
  Real val;
  val = beta / ( 2.0 * ( 1.0 + std::cosh( beta * ( x - mu ) ) ) );
  return val;
}


// *********************************************************************
// C interface
// *********************************************************************

// Dummy interface for the purpose of testing the C interface.
// Not used in practice.
extern "C" 
void DummyInterface( MPI_Comm comm, int a )
{
  int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );
  if( mpirank == 0 ){
    std::cout << "Comm rank = " << mpirank << std::endl;
    std::cout << "Comm size = " << mpisize << std::endl;
    std::cout << "Dummy inteface is working and is outputing an integer " 
      << a << std::endl;
  }
  return;
}  // -----  end of function DummyInterface  ----- 

extern "C"
void ReadDistSparseMatrixFormattedHeadInterface (
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    MPI_Comm comm )
{
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
  std::ifstream fin;
  if( mpirank == 0 ){
    fin.open(filename);
    if( !fin.good() ){
      ErrorHandling( "File cannot be openeded!" );
    }
    Int dummy;
    fin >> *size >> dummy;
    fin >> *nnz;
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
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
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



extern "C"
void PPEXSIInertiaCountInterface(
    // Input parameters
    int           nrows,                        // Size of the matrix
    int           nnz,                          // Total number of nonzeros in H
    int           nnzLocal,                     // Number of nonzeros in H on this proc
    int           numColLocal,                  // Number of local columns for H
    int*          colptrLocal,                  // Colomn pointer in CSC format
    int*          rowindLocal,                  // Row index pointer in CSC format
    double*       HnzvalLocal,                  // Nonzero value of H in CSC format
    int           isSIdentity,                  // Whether S is an identity matrix. If so, the variable SnzvalLocal is omitted.
    double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
    double        temperature,                  // Temperature, in the same unit as H
    double        numElectronExact,             // Exact number of electrons
    double        muMin0,                       // Initial guess of lower bound for mu
    double        muMax0,                       // Initial guess of upper bound for mu
    int           numPole,                      // Number of shifts in computing the inertia, still called "Pole" for legacy reason
    int           maxIter,                      // Maximum number of iterations for computing the inertia
    double        numElectronTolerance,         // Stopping criterion of inertia count
    int           ordering,                     // SuperLUDIST ordering
    int           npPerPole,                    // Number of processors for each shift, still called "Pole" for legacy reason
    int           npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
    MPI_Comm	    comm,                         // Overall MPI communicator
    // Output parameters
    double*       muMinInertia,                 // Lower bound for mu after inertia count
    double*       muMaxInertia,                 // Upper bound for mu after inertia count
    double*       muLowerEdge,                  // Ne(muLowerEdge) = Ne - eps. For band gapped system
    double*       muUpperEdge,                  // Ne(muUpperEdge) = Ne + eps. For band gapped system
    int*          numIter,                      // Number of actual iterations for inertia count
    double*       muList,                       // The list of shifts
    double*       numElectronList,              // The number of electrons (finite temperature) corresponding to shifts
    int*          info                          // 0: successful exit.  1: unsuccessful
    )
{
  Int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );
  Real timeSta, timeEnd;
  // log files
  std::stringstream  ss;
  ss << "logPEXSI" << mpirank;
  // append to previous log files
  statusOFS.open( ss.str().c_str(), std::ios_base::app );

  try{

    if( mpisize % npPerPole != 0 ){
      std::ostringstream msg;
      msg 
        << "mpisize    = " << mpisize << std::endl
        << "npPerPole = " << npPerPole << std::endl
        << "mpisize is not divisible by npPerPole!" << std::endl;
      ErrorHandling( msg.str().c_str() );
    }

    if( numElectronTolerance < 1 ){
      std::ostringstream msg;
      msg 
        << "numElectronTolerance = " << numElectronTolerance 
        << ", which is less than 1. " <<
        "This is probably too tight for the purpose of inertia count." 
        << std::endl;
      ErrorHandling( msg.str().c_str() );
    }


    Int nprow = iround( std::sqrt( (double)npPerPole) );
    Int npcol = npPerPole / nprow;

    GridType gridPole( comm, mpisize / npPerPole, npPerPole );
    PPEXSIData pexsi( &gridPole, nprow, npcol );

    // Convert into H and S matrices
    DistSparseMatrix<Real> HMat, SMat;

    // The first row processors (size: npPerPole) read the matrix, and
    // then distribute the matrix to the rest of the processors.
    //
    // NOTE: The first row processor must have data for H/S.
    std::vector<char> sstr;
    Int sizeStm;
    if( MYROW( &gridPole ) == 0 ){
      std::stringstream sstm;

      HMat.size        = nrows;
      HMat.nnz         = nnz;
      HMat.nnzLocal    = nnzLocal;
      // The first row processor does not need extra copies of the index /
      // value of the matrix. 
      HMat.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
      HMat.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
      // H value
      HMat.nzvalLocal  = DblNumVec( nnzLocal,      false, HnzvalLocal );
      HMat.comm = gridPole.rowComm;

      // Serialization will copy the values regardless of the ownership
      serialize( HMat, sstm, NO_MASK );

      // S value
      if( isSIdentity ){
        SMat.size = 0;
        SMat.nnz  = 0;
        SMat.nnzLocal = 0;
        SMat.comm = HMat.comm; 
      }
      else{
        CopyPattern( HMat, SMat );
        SMat.comm = gridPole.rowComm;
        SMat.nzvalLocal  = DblNumVec( nnzLocal,      false, SnzvalLocal );
        serialize( SMat.nzvalLocal, sstm, NO_MASK );
      }


      sstr.resize( Size( sstm ) );
      sstm.read( &sstr[0], sstr.size() ); 	
      sizeStm = sstr.size();
    }

    MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole.colComm );

#if ( _DEBUGlevel_ >= 0 )
    statusOFS << "sizeStm = " << sizeStm << std::endl;
#endif

    if( MYROW( &gridPole ) != 0 ) sstr.resize( sizeStm );

    MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole.colComm );

    if( MYROW( &gridPole ) != 0 ){
      std::stringstream sstm;
      sstm.write( &sstr[0], sizeStm );
      deserialize( HMat, sstm, NO_MASK );
      // Communicator
      HMat.comm = gridPole.rowComm;
      if( isSIdentity ){
        SMat.size = 0;
        SMat.nnz  = 0;
        SMat.nnzLocal = 0;
        SMat.comm = HMat.comm;
      }
      else{
        CopyPattern( HMat, SMat );
        SMat.comm = gridPole.rowComm;
        deserialize( SMat.nzvalLocal, sstm, NO_MASK );
      }
    }
    sstr.clear();


#if ( _DEBUGlevel_ >= 0 )
    statusOFS << "H.size     = " << HMat.size     << std::endl;
    statusOFS << "H.nnzLocal = " << HMat.nnzLocal << std::endl;
    statusOFS << "H.nnz      = " << HMat.Nnz()    << std::endl;
    statusOFS << "S.size     = " << SMat.size     << std::endl;
    statusOFS << "S.nnzLocal = " << SMat.nnzLocal << std::endl;
    statusOFS << "S.nnz      = " << SMat.Nnz()    << std::endl;
#endif


    // Parameters
    std::string colPerm;
    switch (ordering){
      case 0:
        colPerm = "PARMETIS";
        break;
      case 1:
        colPerm = "METIS_AT_PLUS_A";
        break;
      case 2:
        colPerm = "MMD_AT_PLUS_A";
        break;
      default:
        ErrorHandling("Unsupported ordering strategy.");
    }



    // In inertia counts, "Pole" is the same as "Shifts" (since there is no actual pole here)
    // inertiaVec is the zero-temperature cumulative DOS.
    // inertiaFTVec is the finite temperature cumulative DOS, obtained from linear interpolation of
    // inertiaVec.
    Int  numShift = numPole;
    std::vector<Real>  shiftVec( numShift );
    std::vector<Real>  inertiaVec( numShift );
    std::vector<Real>  inertiaFTVec( numShift );

    Real muMin = muMin0;
    Real muMax = muMax0;
    Real muLow, muUpp;

    bool isConverged = false;
    Int  iter;
    for( iter = 1; iter <= maxIter; iter++ ){
      for( Int l = 0; l < numShift; l++ ){
        shiftVec[l] = muMin + l * (muMax - muMin) / (numShift-1);
      }


      GetTime( timeSta );
      pexsi.CalculateNegativeInertiaReal( 
          shiftVec,
          inertiaVec,
          HMat,
          SMat,
          colPerm,
          npSymbFact );

      GetTime( timeEnd );



      // Inertia is multiplied by 2.0 to reflect the doubly occupied
      // orbitals.
      for( Int l = 0; l < numShift; l++ ){
        inertiaVec[l] *= 2.0;
      }

      // Interpolation to compute the finite temperature cumulative DOS
      Int numInp = 1000;
      std::vector<Real>  shiftInpVec( numInp );
      std::vector<Real>  inertiaInpVec( numInp ); 
      std::vector<Real>  fdDrvVec( numInp );
      for( Int l = 0; l < numShift; l++ ){
        Real shiftInp0 = shiftVec[l] - 20 * temperature;
        Real shiftInp1 = shiftVec[l] + 20 * temperature;
        Real h = (shiftInp1 - shiftInp0) / (numInp-1);
        for( Int i = 0; i < numInp; i++ ){
          shiftInpVec[i] = shiftInp0 + h * i;
          // fdDrvMu(x) = beta * exp(beta*x)/(1+exp(beta*z))^2
          // Note: compared to fdDrvMu used in pole.cpp, the factor 2.0 is not
          // present, because it is included in inertiaVec. 
          fdDrvVec[i]    = 1.0 / ( 2.0 * temperature * (
                1.0 + std::cosh( ( shiftInpVec[i] - shiftVec[l] ) / temperature ) ) );
        }
        LinearInterpolation( shiftVec, inertiaVec, 
            shiftInpVec, inertiaInpVec );

        inertiaFTVec[l] = 0.0;
        for( Int i = 0; i < numInp; i++ ){
          inertiaFTVec[l] += fdDrvVec[i] * inertiaInpVec[i] * h;
        }

#if ( _DEBUGlevel_ >= 0 )
        statusOFS << std::setiosflags(std::ios::left) 
          << std::setw(LENGTH_VAR_NAME) << "Shift      = "
          << std::setw(LENGTH_VAR_DATA) << shiftVec[l]
          << std::setw(LENGTH_VAR_NAME) << "Inertia    = "
          << std::setw(LENGTH_VAR_DATA) << inertiaVec[l]
          << std::setw(LENGTH_VAR_NAME) << "InertiaFT  = "
          << std::setw(LENGTH_VAR_DATA) << inertiaFTVec[l]
          << std::endl;
#endif
      } // for (l)

#if ( _DEBUGlevel_ >= 0 )
      statusOFS << std::endl << "Time for iteration " 
        << iter << " of the inertia count is " 
        << timeEnd - timeSta << std::endl;
#endif

      // Convert the internal variables to output parameters
      // This early output is mainly for debugging purpose.
      {
        for( Int l = 0; l < numShift; l++ ){
          muList[l]          = shiftVec[l];
          numElectronList[l] = inertiaFTVec[l];
        }
      }


      if( inertiaFTVec[0] >= numElectronExact ||
          inertiaFTVec[numShift-1] <= numElectronExact ){
        std::ostringstream msg;
        msg 
          << "The solution is not in the prescribed (muMin, muMax) interval." << std::endl
          << "(muMin, muMax) ~ (" << shiftVec[0] << " , " << shiftVec[numShift-1] << " ) " << std::endl
          << "(Ne(muMin), Ne(muMax)) ~ (" << inertiaFTVec[0] << " , " << inertiaFTVec[numShift-1] 
          << " ) " << std::endl
          << "NeExact = " << numElectronExact << std::endl
          << "Try to increase numElectronTolerance or initial search interval for mu." << std::endl;
        ErrorHandling( msg.str().c_str() );
      }


      // Update muMin, muMax
      // NOTE: divide by 4 is just heuristics so that the interval does
      // not get to be too small for the PEXSI iteration

      // First find the smallest interval
      const Real EPS = 1e-1;
      Int idxMin = 0, idxMax = numShift-1;
      for( Int i = 0; i < numShift; i++ ){
        if( ( inertiaFTVec[i] < numElectronExact ) &&
            ( inertiaVec[i]   < numElectronExact - EPS ) ){
          idxMin = ( idxMin < i ) ? i : idxMin;
        }
        if( ( inertiaFTVec[i] > numElectronExact ) &&
            ( inertiaVec[i]   > numElectronExact + EPS ) ){
          idxMax = ( idxMax > i ) ? i : idxMax;
        }
      }

#if ( _DEBUGlevel_ >= 0 )
      statusOFS << "idxMin = " << idxMin << "inertiaVec = " << inertiaVec[idxMin] << std::endl;
      statusOFS << "idxMax = " << idxMax << "inertiaVec = " << inertiaVec[idxMax] << std::endl;
#endif



      if( inertiaFTVec[idxMax] - inertiaFTVec[idxMin] >= 
          numElectronTolerance ){
        // The smallest interval still contain too many eigenvalues
        muMin = shiftVec[idxMin];
        muMax = shiftVec[idxMax];
        Real NeMin, NeMax;
        NeMin = std::min( inertiaFTVec[idxMin], 
            numElectronExact - numElectronTolerance / 2 + EPS );
        NeMin = std::max( inertiaFTVec[0] + EPS, NeMin );
        NeMax = std::max( inertiaFTVec[idxMax],
            numElectronExact + numElectronTolerance / 2 - EPS );
        NeMax = std::min( inertiaFTVec[numShift-1] - EPS, NeMax );
        muMin = MonotoneRootFinding( shiftVec, inertiaFTVec, NeMin );
        muMax = MonotoneRootFinding( shiftVec, inertiaFTVec, NeMax );
      }
      else{
        // Root finding using linear interpolation
        Real NeMin, NeMax;
        NeMin = numElectronExact - numElectronTolerance / 2 + EPS;
        NeMin = std::max( inertiaFTVec[0] + EPS, NeMin );
        NeMax = numElectronExact + numElectronTolerance / 2 - EPS;
        NeMax = std::min( inertiaFTVec[numShift-1] - EPS, NeMax );
        muMin = MonotoneRootFinding( shiftVec, inertiaFTVec, NeMin );
        muMax = MonotoneRootFinding( shiftVec, inertiaFTVec, NeMax );
      }


      double mu0 = (muMin + muMax)/2.0;


      Real NeLower = numElectronExact - EPS;
      Real NeUpper = numElectronExact + EPS;



      muLow = MonotoneRootFinding( shiftVec, inertiaFTVec, NeLower );
      muUpp = MonotoneRootFinding( shiftVec, inertiaFTVec, NeUpper );

      // Convert the internal variables to output parameters
      // This early output is mainly for debugging purpose.
      {
        *muMinInertia = muMin;
        *muMaxInertia = muMax;
        *muLowerEdge  = muLow;
        *muUpperEdge  = muUpp;
        *numIter      = iter;
      }

      // Output some information after EACH STEP for dynamical adjustament.
      {
#if ( _DEBUGlevel_ >= 0 )
        statusOFS << std::endl << "After inertia count iteration " << iter
          << std::endl;
        Print( statusOFS, "muLowerEdge   = ", muLow );
        Print( statusOFS, "muUppperEdge  = ", muUpp );
        Print( statusOFS, "muMin         = ", muMin );
        Print( statusOFS, "muMax         = ", muMax );
        statusOFS << std::endl;
#endif
      }


      if( inertiaFTVec[numShift-1] - inertiaFTVec[0] < numElectronTolerance ){
        isConverged = true;
        break;
      }

    } // for (iter)


#if ( _DEBUGlevel_ >= 0 )
    if( isConverged ){
      statusOFS << std::endl << "Inertia count converged with " << iter
        << " iterations. N(muMax) - N(muMin) = " <<
        inertiaFTVec[numShift-1] - inertiaFTVec[0] << std::endl;
    }
    else {
      statusOFS << std::endl << "Inertia count did not converge with " << iter
        << " iterations. N(muMax) - N(muMin) = " <<
        inertiaFTVec[numShift-1] - inertiaFTVec[0] << std::endl;
    }
#endif

    // Convert the internal variables to output parameters
    *muMinInertia = muMin;
    *muMaxInertia = muMax;
    *muLowerEdge  = muLow;
    *muUpperEdge  = muUpp;
    *numIter      = iter;


    for( Int l = 0; l < numShift; l++ ){
      muList[l]          = shiftVec[l];
      numElectronList[l] = inertiaFTVec[l];
    }

    *info = 0;
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    //		std::cerr << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
    //			<< std::endl << e.what() << std::endl;
    *info = 1;
  }

  // Synchronize the info among all processors. 
  // If any processor gets error message, info = 1
  Int infoAll = 0;
  mpi::Allreduce( info, &infoAll, 1, MPI_MAX, comm  );
  *info = infoAll;

  statusOFS.close(); 
  return;
}  // -----  end of function PPEXSIInertiaCountInterface ----- 


extern "C"
void PPEXSIRawInertiaCountInterface(
    // Input parameters
    int           nrows,                        // Size of the matrix
    int           nnz,                          // Total number of nonzeros in H
    int           nnzLocal,                     // Number of nonzeros in H on this proc
    int           numColLocal,                  // Number of local columns for H
    int*          colptrLocal,                  // Colomn pointer in CSC format
    int*          rowindLocal,                  // Row index pointer in CSC format
    double*       HnzvalLocal,                  // Nonzero value of H in CSC format
    int           isSIdentity,                  // Whether S is an identity matrix. If so, the variable SnzvalLocal is omitted.
    double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
    double        muMin0,                       // Initial guess of lower bound for mu
    double        muMax0,                       // Initial guess of upper bound for mu
    int           numPole,                      // Number of shifts in computing the inertia, still called "Pole" for legacy reason
    int           ordering,                     // SuperLUDIST ordering
    int           npPerPole,                    // Number of processors for each shift, still called "Pole" for legacy reason
    int           npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
    MPI_Comm	    comm,                         // Overall MPI communicator
    // Output parameters
    double*       muList,                       // The list of shifts
    int   *       numElectronList,              // The number of electrons (finite temperature) corresponding to shifts
    int*          info                          // 0: successful exit.  1: unsuccessful
    )
{
  Int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );
  Real timeSta, timeEnd;
  // log files
  std::stringstream  ss;
  ss << "logPEXSI" << mpirank;
  // append to previous log files
  statusOFS.open( ss.str().c_str(), std::ios_base::app );

  try{

    if( mpisize % npPerPole != 0 ){
      std::ostringstream msg;
      msg 
        << "mpisize    = " << mpisize << std::endl
        << "npPerPole = " << npPerPole << std::endl
        << "mpisize is not divisible by npPerPole!" << std::endl;
      ErrorHandling( msg.str().c_str() );
    }

    Int nprow = iround( std::sqrt( (double)npPerPole) );
    Int npcol = npPerPole / nprow;

    GridType gridPole( comm, mpisize / npPerPole, npPerPole );
    PPEXSIData pexsi( &gridPole, nprow, npcol );

    // Convert into H and S matrices
    DistSparseMatrix<Real> HMat, SMat;

    // The first row processors (size: npPerPole) read the matrix, and
    // then distribute the matrix to the rest of the processors.
    //
    // NOTE: The first row processor must have data for H/S.
    std::vector<char> sstr;
    Int sizeStm;
    if( MYROW( &gridPole ) == 0 ){
      std::stringstream sstm;

      HMat.size        = nrows;
      HMat.nnz         = nnz;
      HMat.nnzLocal    = nnzLocal;
      // The first row processor does not need extra copies of the index /
      // value of the matrix. 
      HMat.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
      HMat.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
      // H value
      HMat.nzvalLocal  = DblNumVec( nnzLocal,      false, HnzvalLocal );
      HMat.comm = gridPole.rowComm;

      // Serialization will copy the values regardless of the ownership
      serialize( HMat, sstm, NO_MASK );

      // S value
      if( isSIdentity ){
        SMat.size = 0;
        SMat.nnz  = 0;
        SMat.nnzLocal = 0;
        SMat.comm = HMat.comm; 
      }
      else{
        CopyPattern( HMat, SMat );
        SMat.comm = gridPole.rowComm;
        SMat.nzvalLocal  = DblNumVec( nnzLocal,      false, SnzvalLocal );
        serialize( SMat.nzvalLocal, sstm, NO_MASK );
      }


      sstr.resize( Size( sstm ) );
      sstm.read( &sstr[0], sstr.size() ); 	
      sizeStm = sstr.size();
    }

    MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole.colComm );

#if ( _DEBUGlevel_ >= 0 )
    statusOFS << "sizeStm = " << sizeStm << std::endl;
#endif

    if( MYROW( &gridPole ) != 0 ) sstr.resize( sizeStm );

    MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole.colComm );

    if( MYROW( &gridPole ) != 0 ){
      std::stringstream sstm;
      sstm.write( &sstr[0], sizeStm );
      deserialize( HMat, sstm, NO_MASK );
      // Communicator
      HMat.comm = gridPole.rowComm;
      if( isSIdentity ){
        SMat.size = 0;
        SMat.nnz  = 0;
        SMat.nnzLocal = 0;
        SMat.comm = HMat.comm;
      }
      else{
        CopyPattern( HMat, SMat );
        SMat.comm = gridPole.rowComm;
        deserialize( SMat.nzvalLocal, sstm, NO_MASK );
      }
    }
    sstr.clear();


#if ( _DEBUGlevel_ >= 0 )
    statusOFS << "H.size     = " << HMat.size     << std::endl;
    statusOFS << "H.nnzLocal = " << HMat.nnzLocal << std::endl;
    statusOFS << "S.size     = " << SMat.size     << std::endl;
    statusOFS << "S.nnzLocal = " << SMat.nnzLocal << std::endl;
#endif


    // Parameters
    std::string colPerm;
    switch (ordering){
      case 0:
        colPerm = "PARMETIS";
        break;
      case 1:
        colPerm = "METIS_AT_PLUS_A";
        break;
      case 2:
        colPerm = "MMD_AT_PLUS_A";
        break;
      default:
        ErrorHandling("Unsupported ordering strategy.");
    }



    Int  numShift = numPole;
    std::vector<Real>  shiftVec( numShift );
    std::vector<Real>  inertiaVec( numShift );

    for( Int l = 0; l < numShift; l++ ){
      shiftVec[l] = muMin0 + l * (muMax0 - muMin0) / (numShift-1);
    }

    GetTime( timeSta );
    pexsi.CalculateNegativeInertiaReal( 
        shiftVec,
        inertiaVec,
        HMat,
        SMat,
        colPerm,
        npSymbFact );

    GetTime( timeEnd );

    // Inertia is multiplied by 2.0 to reflect the doubly occupied
    // orbitals.
    for( Int l = 0; l < numShift; l++ ){
      inertiaVec[l] *= 2.0;
    }

#if ( _DEBUGlevel_ >= 0 )
    statusOFS << std::endl << "Time for the inertia count is " <<
      timeEnd - timeSta << std::endl;
#endif

    // Convert the internal variables to output parameters
    for( Int l = 0; l < numShift; l++ ){
      muList[l]          = shiftVec[l];
      numElectronList[l] = inertiaVec[l];
    }

    *info = 0;
  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    //		std::cerr << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
    //			<< std::endl << e.what() << std::endl;
    *info = 1;
  }

  // Synchronize the info among all processors. 
  // If any processor gets error message, info = 1
  Int infoAll = 0;
  mpi::Allreduce( info, &infoAll, 1, MPI_MAX, comm  );
  *info = infoAll;

  statusOFS.close(); 
  return;
}  // -----  end of function PPEXSIRawInertiaCountInterface ----- 


extern "C" 
void PPEXSISolveInterface (
    // Input parameters
    int           nrows,                        // Size of the matrix
    int           nnz,                          // Total number of nonzeros in H
    int           nnzLocal,                     // Number of nonzeros in H on this proc
    int           numColLocal,                  // Number of local columns for H
    int*          colptrLocal,                  // Colomn pointer in CSC format
    int*          rowindLocal,                  // Row index pointer in CSC format
    double*       HnzvalLocal,                  // Nonzero value of H in CSC format
    int           isSIdentity,                  // Whether S is an identity matrix. If so, the variable SnzvalLocal is omitted.
    double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
    double        temperature,                  // Temperature, in the same unit as H
    double        numElectronExact,             // Exact number of electrons
    double        mu0,                          // Initial guess for mu
    double        muMin0,                       // Initial guess for lower bound of mu
    double        muMax0,                       // Initial guess for upper bound of mu
    double        gap,                          // Energy gap (lower bound)
    double        deltaE,                       // Spectral radius of S^{-1}H
    int           numPole,                      // Number of poles
    int           maxIter,                      // Maximum number of iterations for mu-iteration in PEXSI
    double        numElectronTolerance,         // Stopping criterion of PEXSI mu iteration.
    int           ordering,                     // SuperLUDIST ordering
    int           npPerPole,                    // Number of processors for each pole
    int           npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
    MPI_Comm	    comm,                         // Overall MPI communicator
    // Output parameters
    double*      DMnzvalLocal,                  // Nonzero value of density matrix in CSC format
    double*     EDMnzvalLocal,                  // Nonzero value of energy density matrix in CSC format
    double*     FDMnzvalLocal,                  // Nonzero value of free energy density matrix in CSC format
    double*       muPEXSI,                      // Final chemical potential
    double*       numElectronPEXSI,             // Computed number of electron at the final chemical potential
    double*       muMinPEXSI,                   // Final lower bound for mu.
    double*       muMaxPEXSI,                   // Final upper bound for mu
    int*          numIter,                      // Number of actual iterations for PEXSI
    double*       muList,                       // The history of mu
    double*       numElectronList,              // The history of number of electrons correspondig to mu
    double*       numElectronDrvList,           // The history of dN/dMu
    int*          info                          // 0: successful exit.  1: unsuccessful
    )
{
  Int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );
  Real timeSta, timeEnd;

  // log files
  std::stringstream  ss;
  ss << "logPEXSI" << mpirank;
  // append to previous log files
  statusOFS.open( ss.str().c_str(), std::ios_base::app );

  try{
    if( mpisize % npPerPole != 0 ){
      std::ostringstream msg;
      msg 
        << "mpisize    = " << mpisize << std::endl
        << "npPerPole  = " << npPerPole << std::endl
        << "mpisize is not divisible by npPerPole!" << std::endl;
      ErrorHandling( msg.str().c_str() );
    }


    Int nprow = iround( std::sqrt( (double)npPerPole) );
    Int npcol = npPerPole / nprow;

    GridType gridPole( comm, mpisize / npPerPole, npPerPole );
    PPEXSIData pexsi( &gridPole, nprow, npcol );

    // Convert into H and S matrices
    DistSparseMatrix<Real> HMat, SMat;

    // The first row processors (size: npPerPole) read the matrix, and
    // then distribute the matrix to the rest of the processors.
    //
    // NOTE: The first row processor must have data for H/S.
    std::vector<char> sstr;
    Int sizeStm;
    if( MYROW( &gridPole ) == 0 ){
      std::stringstream sstm;

      HMat.size        = nrows;
      HMat.nnz         = nnz;
      HMat.nnzLocal    = nnzLocal;
      // The first row processor does not need extra copies of the index /
      // value of the matrix. 
      HMat.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
      HMat.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
      // H value
      HMat.nzvalLocal  = DblNumVec( nnzLocal,      false, HnzvalLocal );
      HMat.comm = gridPole.rowComm;

      // Serialization will copy the values regardless of the ownership
      serialize( HMat, sstm, NO_MASK );

      // S value
      if( isSIdentity ){
        SMat.size = 0;
        SMat.nnz  = 0;
        SMat.nnzLocal = 0;
        SMat.comm = HMat.comm;
      }
      else{
        CopyPattern( HMat, SMat );
        SMat.comm = gridPole.rowComm;
        SMat.nzvalLocal  = DblNumVec( nnzLocal,      false, SnzvalLocal );
        serialize( SMat.nzvalLocal, sstm, NO_MASK );
      }


      sstr.resize( Size( sstm ) );
      sstm.read( &sstr[0], sstr.size() ); 	
      sizeStm = sstr.size();
    }

    MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole.colComm );

#if ( _DEBUGlevel_ >= 0 )
    statusOFS << "sizeStm = " << sizeStm << std::endl;
#endif

    if( MYROW( &gridPole ) != 0 ) sstr.resize( sizeStm );

    MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole.colComm );

    if( MYROW( &gridPole ) != 0 ){
      std::stringstream sstm;
      sstm.write( &sstr[0], sizeStm );
      deserialize( HMat, sstm, NO_MASK );
      // Communicator
      HMat.comm = gridPole.rowComm;
      if( isSIdentity ){
        SMat.size = 0;
        SMat.nnz  = 0;
        SMat.nnzLocal = 0;
        SMat.comm = HMat.comm;
      }
      else{
        CopyPattern( HMat, SMat );
        SMat.comm = gridPole.rowComm;
        deserialize( SMat.nzvalLocal, sstm, NO_MASK );
      }
    }
    sstr.clear();


#if ( _DEBUGlevel_ >= 0 )
    statusOFS << "H.size     = " << HMat.size     << std::endl;
    statusOFS << "H.nnzLocal = " << HMat.nnzLocal << std::endl;
    statusOFS << "S.size     = " << SMat.size     << std::endl;
    statusOFS << "S.nnzLocal = " << SMat.nnzLocal << std::endl;
#endif


    // Parameters
    bool isFreeEnergyDensityMatrix = true;
    bool isEnergyDensityMatrix     = true;
    bool isDerivativeTMatrix       = false;
    bool isConverged               = false; 	
    std::string colPerm;
    switch (ordering){
      case 0:
        colPerm = "PARMETIS";
        break;
      case 1:
        colPerm = "METIS_AT_PLUS_A";
        break;
      case 2:
        colPerm = "MMD_AT_PLUS_A";
        break;
      default:
        ErrorHandling("Unsupported ordering strategy.");
    }

    std::vector<Real>  muVec;
    std::vector<Real>  numElectronVec;
    std::vector<Real>  numElectronDrvVec;

    Real timeSolveSta, timeSolveEnd;

    Real muMin = muMin0; 
    Real muMax = muMax0;
    Real mu    = mu0;

    GetTime( timeSolveSta );
    pexsi.Solve( 
        numPole,
        temperature,
        numElectronExact,
        gap,
        deltaE,
        mu,
        muMin,
        muMax,
        HMat,
        SMat,
        maxIter,
        numElectronTolerance,
        colPerm,
        npSymbFact,
        isFreeEnergyDensityMatrix,
        isEnergyDensityMatrix,
        isDerivativeTMatrix,
        muVec,
        numElectronVec,
        numElectronDrvVec,
        isConverged	);
    GetTime( timeSolveEnd );

    Int muIter = muVec.size();

#if ( _DEBUGlevel_ >= 0 )
    PrintBlock( statusOFS, "Solve finished." );
    if( isConverged ){
      statusOFS << "PEXSI has converged with " << muIter << 
        " iterations" << std::endl;
    }
    else {
      statusOFS << "PEXSI did not converge with " << muIter << 
        " iterations" << std::endl;
    }
#endif


    // Convert the internal variables to output parameters

    // Update the density matrices for the first row processors (size: npPerPole) 
    if( MYROW( &gridPole ) == 0 ){
      blas::Copy( nnzLocal, pexsi.DensityMatrix().nzvalLocal.Data(), 
          1, DMnzvalLocal, 1 );

      blas::Copy( nnzLocal, pexsi.FreeEnergyDensityMatrix().nzvalLocal.Data(), 
          1, FDMnzvalLocal, 1 );

      blas::Copy( nnzLocal, pexsi.EnergyDensityMatrix().nzvalLocal.Data(), 
          1, EDMnzvalLocal, 1 );
    }

    *muPEXSI          = *(muVec.rbegin());
    *numElectronPEXSI = *(numElectronVec.rbegin());
    *muMinPEXSI       = muMin;
    *muMaxPEXSI       = muMax;
    *numIter          = muIter;
    for( Int i = 0; i < muIter; i++ ){
      muList[i]               = muVec[i];
      numElectronList[i]      = numElectronVec[i];
      numElectronDrvList[i]   = numElectronDrvVec[i];
    }


    MPI_Barrier( comm );


    // Compute the guess for T=0 chemical potential. Not used anymore
    //	*muZeroT = pexsi.EstimateZeroTemperatureChemicalPotential(
    //			temperature,
    //			*mu,
    //			SMat );
    //
    //	Print( statusOFS, "guess of mu(T=0) = ", 
    //			*muZeroT );

#if ( _DEBUGlevel_ >= 0 )
    Print( statusOFS, "Total time for PEXSI = ", 
        timeSolveEnd - timeSolveSta );
#endif

    *info = 0;
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    //		std::cerr << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
    //			<< std::endl << e.what() << std::endl;
    *info = 1;
  }

  // Synchronize the info among all processors. 
  // If any processor gets error message, info = 1
  Int infoAll = 0;
  mpi::Allreduce( info, &infoAll, 1, MPI_MAX, comm  );
  *info = infoAll;

  statusOFS.close();
  return;
}  // -----  end of function PPEXSISolveInterface ----- 

extern "C" 
void PPEXSISelInvInterface (
    // Input parameters
    int           nrows,                        // Size of the matrix
    int           nnz,                          // Total number of nonzeros in H
    int           nnzLocal,                     // Number of nonzeros in H on this proc
    int           numColLocal,                  // Number of local columns for H
    int*          colptrLocal,                  // Colomn pointer in CSC format
    int*          rowindLocal,                  // Row index pointer in CSC format
    double*       HnzvalLocal,                  // Nonzero value of H in CSC format
    int           isSIdentity,                  // Whether S is an identity matrix. If so, the variable SnzvalLocal is omitted.
    double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
    double*       zShift,                       // Shift. Real: zShift[0], Imag: zShift[1]    
    int           ordering,                     // SuperLUDIST ordering
    int           npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
    MPI_Comm	    comm,                         // Overall MPI communicator
    // Output parameters
    double*       AinvnzvalLocal,               // Nonzero value of Ainv = (H - z S)^{-1}
    int*          info                          // 0: successful exit.  1: unsuccessful
    )
{
  Int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );
  Real timeSta, timeEnd;


  // log files
  std::stringstream  ss;
  ss << "logPEXSI" << mpirank;
  // append to previous log files
  statusOFS.open( ss.str().c_str(), std::ios_base::app );


  try{
    Int nprow = iround( std::sqrt( (double)mpisize) );
    Int npcol = mpisize / nprow;
    if( mpisize != nprow * npcol || nprow != npcol ){
      ErrorHandling( "nprow == npcol is assumed in this test routine." );
    }

    SuperLUGrid<Scalar> g( comm, nprow, npcol );


    // Convert into H and S matrices
    DistSparseMatrix<Real> HMat, SMat;

    // The first row processors (size: npPerPole) read the matrix, and
    // then distribute the matrix to the rest of the processors.
    //
    // NOTE: The first row processor must have data for H/S.
    {
      HMat.size        = nrows;
      HMat.nnz         = nnz;
      HMat.nnzLocal    = nnzLocal;
      // The first row processor does not need extra copies of the index /
      // value of the matrix. 
      HMat.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
      HMat.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
      // H value
      HMat.nzvalLocal  = DblNumVec( nnzLocal,      false, HnzvalLocal );
      HMat.comm        = comm;

      // S value
      if( isSIdentity ){
        SMat.size = 0;
        SMat.nnz  = 0;
        SMat.nnzLocal = 0;
        SMat.comm = HMat.comm;
      }
      else{
        CopyPattern( HMat, SMat );
        SMat.comm = comm;
        SMat.nzvalLocal  = DblNumVec( nnzLocal,      false, SnzvalLocal );
      }
    }

    // Get the diagonal indices for H and save it n diagIdxLocal

    std::vector<Int>  diagIdxLocal;
    { 
      Int numColLocal      = HMat.colptrLocal.m() - 1;
      Int numColLocalFirst = HMat.size / mpisize;
      Int firstCol         = mpirank * numColLocalFirst;

      diagIdxLocal.clear();

      for( Int j = 0; j < numColLocal; j++ ){
        Int jcol = firstCol + j + 1;
        for( Int i = HMat.colptrLocal(j)-1; 
            i < HMat.colptrLocal(j+1)-1; i++ ){
          Int irow = HMat.rowindLocal(i);
          if( irow == jcol ){
            diagIdxLocal.push_back( i );
          }
        }
      } // for (j)
    }


    DistSparseMatrix<Complex>  AMat;
    AMat.size          = HMat.size;
    AMat.nnz           = HMat.nnz;
    AMat.nnzLocal      = HMat.nnzLocal;
    AMat.colptrLocal   = HMat.colptrLocal;
    AMat.rowindLocal   = HMat.rowindLocal;
    AMat.nzvalLocal.Resize( HMat.nnzLocal );
    AMat.comm          = comm;

    Complex *ptr0 = AMat.nzvalLocal.Data();
    Real *ptr1 = HMat.nzvalLocal.Data();
    Real *ptr2 = SMat.nzvalLocal.Data();
    Complex zshift = Complex(zShift[0], zShift[1]);

    if( SMat.size != 0 ){
      // S is not an identity matrix
      for( Int i = 0; i < HMat.nnzLocal; i++ ){
        AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift * SMat.nzvalLocal(i);
      }
    }
    else{
      // S is an identity matrix
      for( Int i = 0; i < HMat.nnzLocal; i++ ){
        AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
      }

      for( Int i = 0; i < diagIdxLocal.size(); i++ ){
        AMat.nzvalLocal( diagIdxLocal[i] ) -= zshift;
      }
    } // if (SMat.size != 0 )

    std::string colPerm;
    switch (ordering){
      case 0:
        colPerm = "PARMETIS";
        break;
      case 1:
        colPerm = "METIS_AT_PLUS_A";
        break;
      case 2:
        colPerm = "MMD_AT_PLUS_A";
        break;
      default:
        ErrorHandling("Unsupported ordering strategy.");
    }


    // *********************************************************************
    // Symbolic factorization 
    // *********************************************************************

    SuperLUOptions luOpt;
    luOpt.ColPerm = colPerm;
    luOpt.numProcSymbFact = npSymbFact;
    // TODO Introduce maxPipelineDepth as an adjustable parameter when needed.
    luOpt.maxPipelineDepth = -1;

    SuperLUMatrix<Scalar> luMat( g, luOpt );
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );

    luMat.SymbolicFactorize();
    luMat.DestroyAOnly();
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
    luMat.Distribute();
    luMat.NumericalFactorize();

    // *********************************************************************
    // Selected inversion
    // *********************************************************************

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "After numerical factorization." << std::endl;
#endif

    GridType g1( comm, nprow, npcol );
    SuperNodeType super;

    luMat.SymbolicToSuperNode( super );
    PMatrix<Scalar> PMloc( &g1, &super, &luOpt );
    luMat.LUstructToPMatrix( PMloc );


    // P2p communication version
    PMloc.ConstructCommunicationPattern();

    // Collective communication version
    //		PMloc.ConstructCommunicationPattern_Collectives();

    PMloc.PreSelInv();

    // Main subroutine for selected inversion
    // Use the broadcast pipelined version of SelInv.

    // P2p communication version
    PMloc.SelInv();

    // Collective communication version
    //		PMloc.SelInv_Collectives();

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "After selected inversion." << std::endl;
#endif

    DistSparseMatrix<Complex>  AinvMat;       // A^{-1} in DistSparseMatrix format

    PMloc.PMatrixToDistSparseMatrix2( AMat, AinvMat );
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "After conversion to DistSparseMatrix." << std::endl;
#endif

    // Convert the internal variables to output parameters

    blas::Copy( nnzLocal*2, reinterpret_cast<double*>(AinvMat.nzvalLocal.Data()), 
        1, AinvnzvalLocal, 1 );

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "After copying the matrix of Ainv." << std::endl;
#endif



    *info = 0;
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    //		std::cerr  << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
    //			<< std::endl << e.what() << std::endl;
    *info = 1;
  }

  // Synchronize the info among all processors. 
  // If any processor gets error message, info = 1
  Int infoAll = 0;
  mpi::Allreduce( info, &infoAll, 1, MPI_MAX, comm  );
  *info = infoAll;

  statusOFS.close();
  return;
}  // -----  end of function PPEXSISelInvInterface ----- 



extern "C" 
void PPEXSILocalDOSInterface (
    // Input parameters
    int           nrows,                        // Size of the matrix
    int           nnz,                          // Total number of nonzeros in H
    int           nnzLocal,                     // Number of nonzeros in H on this proc
    int           numColLocal,                  // Number of local columns for H
    int*          colptrLocal,                  // Colomn pointer in CSC format
    int*          rowindLocal,                  // Row index pointer in CSC format
    double*       HnzvalLocal,                  // Nonzero value of H in CSC format
    int           isSIdentity,                  // Whether S is an identity matrix. If so, the variable SnzvalLocal is omitted.
    double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
    double        Energy,                       // Real part of the shift
    double        eta,                          // Broadening parameter
    int           ordering,                     // SuperLUDIST ordering
    int           npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
    MPI_Comm	    comm,                         // Overall MPI communicator
    // Output parameters
    double*       localDOSnzvalLocal,           // Nonzero value of Im 1/pi (H - (E+ieta) S)^{-1}
    int*          info                          // 0: successful exit.  1: unsuccessful
    )
{
  Int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );
  Real timeSta, timeEnd;

  // log files
  std::stringstream  ss;
  ss << "logPEXSI" << mpirank;
  // append to previous log files
  statusOFS.open( ss.str().c_str(), std::ios_base::app );

  try{
    double zShift[2];
    zShift[0] = Energy;
    zShift[1] = eta;
    DblNumVec Ainvnzval( 2 * nnzLocal );
    SetValue( Ainvnzval, 0.0 );


    statusOFS.close();

    PPEXSISelInvInterface(
        nrows,
        nnz,
        nnzLocal,
        numColLocal,
        colptrLocal,
        rowindLocal,
        HnzvalLocal,
        isSIdentity,
        SnzvalLocal,
        zShift,
        ordering,
        npSymbFact,
        comm,	
        Ainvnzval.Data(),
        info );

    statusOFS.open( ss.str().c_str(), std::ios_base::app );


    if( *info ){
      ErrorHandling("Error in SelInv!.");
    }

    // Get the imaginary part
    blas::Copy( nnzLocal, Ainvnzval.Data()+1, 2, localDOSnzvalLocal, 1 );
    // Scale the imaginary part by 1/pi
    blas::Scal( nnzLocal, 1.0 / PI, localDOSnzvalLocal, 1 );
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    //		std::cerr  << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
    //			<< std::endl << e.what() << std::endl;
    *info = 1;
  }

  // Synchronize the info among all processors. 
  // If any processor gets error message, info = 1
  Int infoAll = 0;
  mpi::Allreduce( info, &infoAll, 1, MPI_MAX, comm  );
  *info = infoAll;

  statusOFS.close();
  return;
}  // -----  end of function PPEXSILocalDOSInterface  ----- 


extern "C"
void PSelInvComplexSymmetricInterface (
    int           nrows,                        
    int           nnz,                          
    int           nnzLocal,                     
    int           numColLocal,                  
    int*          colptrLocal,                  
    int*          rowindLocal,                  
    double*       AnzvalLocal,                  
    int           ordering,                
    int           npSymbFact,                   
    MPI_Comm	    comm,                         
    int           nprow,
    int           npcol,
    double*       AinvnzvalLocal,
    int*          info
    )
{
  Int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );
  Real timeSta, timeEnd;


  // log files
  std::stringstream  ss;
  ss << "logPEXSI" << mpirank;
  // append to previous log files
  statusOFS.open( ss.str().c_str(), std::ios_base::app );


  try{
    if( mpisize != nprow * npcol ){
      ErrorHandling( " mpisize != nprow * npcol ." );
    }

    SuperLUGrid<Complex> g( comm, nprow, npcol );

    // Convert into DistSparseMatrix
    DistSparseMatrix<Complex> AMat;

    {
      AMat.size        = nrows;
      AMat.nnz         = nnz;
      AMat.nnzLocal    = nnzLocal;
      // The first row processor does not need extra copies of the index /
      // value of the matrix. 
      AMat.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
      AMat.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
      // H value
      AMat.nzvalLocal  = CpxNumVec( nnzLocal,      false, 
          reinterpret_cast<Complex*>(AnzvalLocal) );
      AMat.comm        = comm;
    }

    std::string colPerm;
    switch (ordering){
      case 0:
        colPerm = "PARMETIS";
        break;
      case 1:
        colPerm = "METIS_AT_PLUS_A";
        break;
      case 2:
        colPerm = "MMD_AT_PLUS_A";
        break;
      default:
        ErrorHandling("Unsupported ordering strategy.");
    }


    // *********************************************************************
    // Symbolic factorization 
    // *********************************************************************

    SuperLUOptions luOpt;
    luOpt.ColPerm = colPerm;
    luOpt.numProcSymbFact = npSymbFact;
    // TODO Introduce maxPipelineDepth as an adjustable parameter when needed.
    luOpt.maxPipelineDepth = -1;

    SuperLUMatrix<Complex> luMat( g, luOpt );
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );

    luMat.SymbolicFactorize();
    luMat.DestroyAOnly();
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
    luMat.Distribute();
    luMat.NumericalFactorize();

    // *********************************************************************
    // Selected inversion
    // *********************************************************************

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "After numerical factorization." << std::endl;
#endif

    GridType g1( comm, nprow, npcol );
    SuperNodeType super;

    luMat.SymbolicToSuperNode( super );
    PMatrix<Complex> PMloc( &g1, &super, &luOpt );
    luMat.LUstructToPMatrix( PMloc );


    // P2p communication version
    PMloc.ConstructCommunicationPattern();

    // Collective communication version
    //		PMloc.ConstructCommunicationPattern_Collectives();

    PMloc.PreSelInv();

    // Main subroutine for selected inversion
    // Use the broadcast pipelined version of SelInv.

    // P2p communication version
    PMloc.SelInv();

    // Collective communication version
    //		PMloc.SelInv_Collectives();

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "After selected inversion." << std::endl;
#endif

    DistSparseMatrix<Complex>  AinvMat;       // A^{-1} in DistSparseMatrix format

    PMloc.PMatrixToDistSparseMatrix2( AMat, AinvMat );
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "After conversion to DistSparseMatrix." << std::endl;
#endif

    // Convert the internal variables to output parameters

    blas::Copy( nnzLocal*2, reinterpret_cast<double*>(AinvMat.nzvalLocal.Data()), 
        1, AinvnzvalLocal, 1 );

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "After copying the matrix of Ainv." << std::endl;
#endif



    *info = 0;
  }
  catch( std::exception& e ) {
    statusOFS << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
      << std::endl << e.what() << std::endl;
    //		std::cerr  << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
    //			<< std::endl << e.what() << std::endl;
    *info = 1;
  }

  // Synchronize the info among all processors. 
  // If any processor gets error message, info = 1
  Int infoAll = 0;
  mpi::Allreduce( info, &infoAll, 1, MPI_MAX, comm  );
  *info = infoAll;

  statusOFS.close();
  return;
}  // -----  end of function PSelInvComplexSymmetricInterface ----- 


// ********************************************************************* 
// FORTRAN interface
// 
// All FORTRAN interfaces call the C-interface subroutines above for
// actual computation.  
//
// NOTE: 
//
// 1. The FORTRAN communicators are converted to C communicators using 
// f2c_comm.
//
// 2. All string arguments should have the "correct input" from FORTRAN
// that is consistent with C-convention. This can be done by either
// using iso_c_binding, or using trim(string)//char(0) instead of string
// in FORTRAN.
//
// *********************************************************************

/// @brief Internal subroutine to convert FORTRAN communicator to C
extern "C" 
MPI_Comm f2c_comm(int *Fcomm)
{
  return MPI_Comm_f2c((MPI_Fint)(*Fcomm));
}  // -----  end of function f2c_comm ----- 



// Dummy interface for the purpose of testing the FORTRAN
// interface.

// Not used in practice.
extern "C" 
void FORTRAN(f_dummy_interface)( int* Fcomm, int* a ){
  DummyInterface( f2c_comm(Fcomm), *a );
  return;
}  // -----  end of function f_dummy_interface  ----- 

/// @brief FORTRAN interface for @ref ReadDistSparseMatrixFormattedHeadInterface.
extern "C"
void FORTRAN(f_read_distsparsematrix_formatted_head) (
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    int*     Fcomm )
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
void FORTRAN(f_read_distsparsematrix_formatted) (
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    int*     colptrLocal,
    int*     rowindLocal,
    double*  nzvalLocal,
    int*     Fcomm )
{
  ReadDistSparseMatrixFormattedInterface(
      filename,
      *size,
      *nnz,
      *nnzLocal,
      *numColLocal,
      colptrLocal,
      rowindLocal,
      nzvalLocal,
      f2c_comm( Fcomm ) );
  return;
} // -----  end of function f_read_distsparsematrix_formatted  ----- 

/// @brief FORTRAN interface for @ref ReadDistSparseMatrixHeadInterface.
extern "C"
void FORTRAN(f_read_distsparsematrix_head) (
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    int*     Fcomm )
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
void FORTRAN(f_para_read_distsparsematrix) (
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    int*     colptrLocal,
    int*     rowindLocal,
    double*  nzvalLocal,
    int*     Fcomm )
{
  ParaReadDistSparseMatrixInterface(
      filename,
      *size,
      *nnz,
      *nnzLocal,
      *numColLocal,
      colptrLocal,
      rowindLocal,
      nzvalLocal,
      f2c_comm( Fcomm ) );
  return;
} // -----  end of function f_para_read_distsparsematrix  ----- 



/// @brief FORTRAN interface for @ref PPEXSIInertiaCountInterface.
extern "C" 
void FORTRAN(f_ppexsi_inertiacount_interface)(
    // Input parameters
    int*          nrows,                        // Size of the matrix
    int*          nnz,                          // Total number of nonzeros in H
    int*          nnzLocal,                     // Number of nonzeros in H on this proc
    int*          numColLocal,                  // Number of local columns for H
    int*          colptrLocal,                  // Colomn pointer in CSC format
    int*          rowindLocal,                  // Row index pointer in CSC format
    double*       HnzvalLocal,                  // Nonzero value of H in CSC format
    int*          isSIdentity,                  // Whether S is an identity matrix. If so, the variable SnzvalLocal is omitted.
    double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
    double*       temperature,                  // Temperature, in the same unit as H
    double*       numElectronExact,             // Exact number of electrons
    double*       muMin0,                       // Initial guess of lower bound for mu
    double*       muMax0,                       // Initial guess of upper bound for mu
    int*          numPole,                      // Number of shifts in computing the inertia, still called "Pole" for legacy reason
    int*          maxIter,                      // Maximum number of iterations for computing the inertia
    double*       numElectronTolerance,         // Stopping criterion of inertia count
    int*          ordering,                     // SuperLUDIST ordering
    int*          npPerPole,                    // Number of processors for each shift, still called "Pole" for legacy reason
    int*          npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
    int*    	    Fcomm,                        // Overall MPI communicator
    // Output parameters
    double*       muMinInertia,                 // Lower bound for mu after inertia count
    double*       muMaxInertia,                 // Upper bound for mu after inertia count
    double*       muLowerEdge,                  // Ne(muLowerEdge) = Ne - eps. For band gapped system
    double*       muUpperEdge,                  // Ne(muUpperEdge) = Ne + eps. For band gapped system
    int*          numIter,                      // Number of actual iterations for inertia count
    double*       muList,                       // The list of shifts
    double*       numElectronList,              // The number of electrons (finite temperature) corresponding to shifts
    int*          info                          // 0: successful exit.  1: unsuccessful
    )
{
  PPEXSIInertiaCountInterface( 
      *nrows,
      *nnz,
      *nnzLocal,
      *numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      *isSIdentity,
      SnzvalLocal,
      *temperature,
      *numElectronExact,
      *muMin0,
      *muMax0,
      *numPole,
      *maxIter,
      *numElectronTolerance,
      *ordering,
      *npPerPole,
      *npSymbFact,
      f2c_comm(Fcomm),
      muMinInertia,
      muMaxInertia,
      muLowerEdge,
      muUpperEdge,
      numIter,
      muList,
      numElectronList,
      info );

  return;
} // -----  end of function f_ppexsi_inertiacount_interface  ----- 


/// @brief FORTRAN interface for @ref PPEXSIRawInertiaCountInterface.
extern "C" 
void FORTRAN(f_ppexsi_raw_inertiacount_interface)(
    // Input parameters
    int*          nrows,                        // Size of the matrix
    int*          nnz,                          // Total number of nonzeros in H
    int*          nnzLocal,                     // Number of nonzeros in H on this proc
    int*          numColLocal,                  // Number of local columns for H
    int*          colptrLocal,                  // Colomn pointer in CSC format
    int*          rowindLocal,                  // Row index pointer in CSC format
    double*       HnzvalLocal,                  // Nonzero value of H in CSC format
    int*          isSIdentity,                  // Whether S is an identity matrix. If so, the variable SnzvalLocal is omitted.
    double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
    double*       muMin0,                       // Initial guess of lower bound for mu
    double*       muMax0,                       // Initial guess of upper bound for mu
    int*          numPole,                      // Number of shifts in computing the inertia, still called "Pole" for legacy reason
    int*          ordering,                     // SuperLUDIST ordering
    int*          npPerPole,                    // Number of processors for each shift, still called "Pole" for legacy reason
    int*          npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
    int*    	    Fcomm,                        // Overall MPI communicator
    // Output parameters
    double*       muList,                       // The list of shifts
    int*          numElectronList,              // The number of electrons (finite temperature) corresponding to shifts
    int*          info                          // 0: successful exit.  1: unsuccessful
    )
{
  PPEXSIRawInertiaCountInterface( 
      *nrows,
      *nnz,
      *nnzLocal,
      *numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      *isSIdentity,
      SnzvalLocal,
      *muMin0,
      *muMax0,
      *numPole,
      *ordering,
      *npPerPole,
      *npSymbFact,
      f2c_comm(Fcomm),
      muList,
      numElectronList,
      info );

  return;
} // -----  end of function f_ppexsi_raw_inertiacount_interface  ----- 


/// @brief FORTRAN interface for @ref PPEXSISolveInterface.
extern "C" 
void FORTRAN(f_ppexsi_solve_interface)(
    // Input parameters
    int*          nrows,                        // Size of the matrix
    int*          nnz,                          // Total number of nonzeros in H
    int*          nnzLocal,                     // Number of nonzeros in H on this proc
    int*          numColLocal,                  // Number of local columns for H
    int*          colptrLocal,                  // Colomn pointer in CSC format
    int*          rowindLocal,                  // Row index pointer in CSC format
    double*       HnzvalLocal,                  // Nonzero value of H in CSC format
    int*          isSIdentity,                  // Whether S is an identity matrix. If so, the variable SnzvalLocal is omitted.
    double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
    double*       temperature,                  // Temperature, in the same unit as H
    double*       numElectronExact,             // Exact number of electrons
    double*       mu0,                          // Initial guess for mu
    double*       muMin0,                       // Initial guess for lower bound of mu
    double*       muMax0,                       // Initial guess for upper bound of mu
    double*       gap,                          // Energy gap (lower bound)
    double*       deltaE,                       // Spectral radius of S^{-1}H
    int*          numPole,                      // Number of poles
    int*          maxIter,                      // Maximum number of iterations for mu-iteration in PEXSI
    double*       numElectronTolerance,         // Stopping criterion of PEXSI mu iteration.
    int*          ordering,                     // SuperLUDIST ordering
    int*          npPerPole,                    // Number of processors for each pole
    int*          npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
    int*    	    Fcomm,                        // Overall MPI communicator
    // Output parameters
    double*      DMnzvalLocal,                  // Nonzero value of density matrix in CSC format
    double*     EDMnzvalLocal,                  // Nonzero value of energy density matrix in CSC format
    double*     FDMnzvalLocal,                  // Nonzero value of free energy density matrix in CSC format
    double*       muPEXSI,                      // Final chemical potential
    double*       numElectronPEXSI,             // Computed number of electron at the final chemical potential
    double*       muMinPEXSI,                   // Final lower bound for mu.
    double*       muMaxPEXSI,                   // Final upper bound for mu
    int*          numIter,                      // Number of actual iterations for PEXSI
    double*       muList,                       // The history of mu
    double*       numElectronList,              // The history of number of electrons correspondig to mu
    double*       numElectronDrvList,           // The history of dN/dMu
    int*          info                          // 0: successful exit.  1: unsuccessful
    )
{
  PPEXSISolveInterface( 
      *nrows,
      *nnz,
      *nnzLocal,
      *numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      *isSIdentity,
      SnzvalLocal,
      *temperature,
      *numElectronExact,
      *mu0,
      *muMin0,
      *muMax0,
      *gap,
      *deltaE,
      *numPole,
      *maxIter,
      *numElectronTolerance,
      *ordering,
      *npPerPole,
      *npSymbFact,
      f2c_comm(Fcomm),
      DMnzvalLocal,
      EDMnzvalLocal,
      FDMnzvalLocal,
      muPEXSI,
      numElectronPEXSI,
      muMinPEXSI,
      muMaxPEXSI,
      numIter,
      muList,
      numElectronList,
      numElectronDrvList,
      info);

  return;
} // -----  end of function f_ppexsi_solve_interface  ----- 


/// @brief FORTRAN interface for @ref PPEXSISelInvInterface.
extern "C" 
void FORTRAN(f_ppexsi_selinv_interface)(
    // Input parameters
    int*          nrows,                        // Size of the matrix
    int*          nnz,                          // Total number of nonzeros in H
    int*          nnzLocal,                     // Number of nonzeros in H on this proc
    int*          numColLocal,                  // Number of local columns for H
    int*          colptrLocal,                  // Colomn pointer in CSC format
    int*          rowindLocal,                  // Row index pointer in CSC format
    double*       HnzvalLocal,                  // Nonzero value of H in CSC format
    int*          isSIdentity,                  // Whether S is an identity matrix. If so, the variable SnzvalLocal is omitted.
    double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
    double*       zShift,                       // Shift. Real: zShift[0], Imag: zShift[1]    
    int*          ordering,                     // SuperLUDIST ordering
    int*          npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
    int*    	    Fcomm,                        // Overall MPI communicator
    // Output parameters
    double*       AinvnzvalLocal,               // Nonzero value of Ainv = (H - z S)^{-1}
    int*          info                          // 0: successful exit.  1: unsuccessful
    )
{
  PPEXSISelInvInterface(
      *nrows,
      *nnz,
      *nnzLocal,
      *numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      *isSIdentity,
      SnzvalLocal,
      zShift,
      *ordering,
      *npSymbFact,
      f2c_comm(Fcomm),
      AinvnzvalLocal,
      info );
} // -----  end of function f_ppexsi_selinv_interface  ----- 

/// @brief FORTRAN interface for @ref PPEXSILocalDOSInterface.
extern "C" 
void FORTRAN(f_ppexsi_localdos_interface)(
    // Input parameters
    int*          nrows,                        // Size of the matrix
    int*          nnz,                          // Total number of nonzeros in H
    int*          nnzLocal,                     // Number of nonzeros in H on this proc
    int*          numColLocal,                  // Number of local columns for H
    int*          colptrLocal,                  // Colomn pointer in CSC format
    int*          rowindLocal,                  // Row index pointer in CSC format
    double*       HnzvalLocal,                  // Nonzero value of H in CSC format
    int*          isSIdentity,                  // Whether S is an identity matrix. If so, the variable SnzvalLocal is omitted.
    double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
    double*       Energy,                       // Real part of the shift
    double*       eta,                          // Broadening parameter
    int*          ordering,                     // SuperLUDIST ordering
    int*          npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
    int*    	    Fcomm,                        // Overall MPI communicator
    // Output parameters
    double*       localDOSnzvalLocal,           // Nonzero value of Im 1/pi (H - (E+ieta) S)^{-1}
    int*          info                          // 0: successful exit.  1: unsuccessful
    )
{
  PPEXSILocalDOSInterface(
      *nrows,
      *nnz,
      *nnzLocal,
      *numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      *isSIdentity,
      SnzvalLocal,
      *Energy,
      *eta,
      *ordering,
      *npSymbFact,
      f2c_comm(Fcomm),
      localDOSnzvalLocal,
      info );
} // -----  end of function f_ppexsi_localdos_interface  ----- 


/// @brief FORTRAN interface for @ref PSelInvComplexSymmetricInterface.
extern "C" 
void FORTRAN(f_pselinv_complex_symmetric_interface)(
    // Input parameters
    int*          nrows,                        // Size of the matrix
    int*          nnz,                          // Total number of nonzeros in H
    int*          nnzLocal,                     // Number of nonzeros in H on this proc
    int*          numColLocal,                  // Number of local columns for H
    int*          colptrLocal,                  // Colomn pointer in CSC format
    int*          rowindLocal,                  // Row index pointer in CSC format
    double*       AnzvalLocal,                  // Nonzero value of A in CSC format.  A is complex.
    int*          ordering,                     // SuperLUDIST ordering
    int*          npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
    int*    	    Fcomm,                        // Overall MPI communicator
    int*          nprow,                        // Number of row processors
    int*          npcol,                        // Number of column processors
    // Output parameters
    double*       AinvnzvalLocal,               // Nonzero value of (A)^{-1}. A^{-1} is complex.
    int*          info                          // 0: successful exit.  1: unsuccessful
    )
{
  PSelInvComplexSymmetricInterface ( 
      *nrows,
      *nnz,
      *nnzLocal,
      *numColLocal,
      colptrLocal,
      rowindLocal,
      AnzvalLocal,
      *ordering,
      *npSymbFact,
      f2c_comm(Fcomm),
      *nprow,
      *npcol,
      AinvnzvalLocal,
      info
      );
} // -----  end of function f_pselinv_complex_symmetric_interface  ----- 



