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
/// @file ppexsi.cpp
/// @brief Implementation of the parallel %PEXSI.
/// @date Original:      2012-11-20  Initially started.
/// @date Revision:      2014-03-09  Second generation interface.
/// @date Revision:      2015-11-25  Update strategy with pole
/// expansion.

// FIXME
#define _DEBUGlevel_ 0

#include "ppexsi.hpp"
#include "pexsi/utility.hpp"

namespace PEXSI{

PPEXSIData::PPEXSIData	(
    MPI_Comm   comm,
    Int        numProcRow, 
    Int        numProcCol, 
    Int        outputFileIndex ){
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::PPEXSIData");
#endif

  Int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );

  Int npPerPole = numProcRow * numProcCol;
  if( mpisize % npPerPole != 0 ){
    std::ostringstream msg;
    msg 
      << "mpisize    = " << mpisize << std::endl
      << "npPerPole = " << npPerPole << std::endl
      << "mpisize is not divisible by npPerPole!" << std::endl;
#ifdef USE_ABORT
    abort();
#endif
    throw std::runtime_error( msg.str().c_str() );
  }

  gridPole_     = new GridType( comm, mpisize / npPerPole, npPerPole );
	gridSuperLUReal_     = new SuperLUGrid<Real>( 
      gridPole_->rowComm, numProcRow, numProcCol );
	gridSuperLUComplex_  = new SuperLUGrid<Complex>( 
      gridPole_->rowComm, numProcRow, numProcCol );

	gridSelInv_   = new GridType( gridPole_->rowComm, 
      numProcRow, numProcCol );

  // Start the log file. Append to previous log files
//#ifndef _RELEASE_
////  // All processors output
////  {
////    std::stringstream ss;
////    ss << "logPEXSI" << outputFileIndex;
////    statusOFS.open( ss.str().c_str(), std::ios_base::app );
////  }
//  // Only master processor output
//  {
//    if( mpirank == 0 ){
//      std::stringstream ss;
//      ss << "logPEXSI" << outputFileIndex;
//      statusOFS.open( ss.str().c_str(), std::ios_base::app );
//    }
//  }
//#else
//  // Only master processor output
//  {
//    if( mpirank == 0 ){
//      std::stringstream ss;
//      ss << "logPEXSI" << outputFileIndex;
//      statusOFS.open( ss.str().c_str(), std::ios_base::app );
//    }
//  }
//#endif
 
  // Gives an option of not outputing log files.
  if( outputFileIndex >= 0 ){
    std::stringstream ss;
    ss << "logPEXSI" << outputFileIndex;
    statusOFS.open( ss.str().c_str(), std::ios_base::app );
  }

  // Initialize the saved variables
  isMatrixLoaded_       = false;
  isRealSymmetricSymbolicFactorized_ = false;
  isComplexSymmetricSymbolicFactorized_ = false;

  // Initialize the empty matrices
  luRealMat_ = new SuperLUMatrix<Real>;
  luComplexMat_ = new SuperLUMatrix<Complex>;

  PMRealMat_ = new PMatrix<Real>;
  PMComplexMat_ = new PMatrix<Complex>;
  
  PMRealUnsymMat_ = new PMatrixUnsym<Real>;
  PMComplexUnsymMat_ = new PMatrixUnsym<Complex>;
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::PPEXSIData  ----- 


PPEXSIData::~PPEXSIData	(  )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::~PPEXSIData");
#endif
  if( luRealMat_ != NULL ){
    delete luRealMat_;
  }

  if( luComplexMat_ != NULL ){
    delete luComplexMat_;
  }

  if( PMRealMat_ != NULL ){
    delete PMRealMat_;
  }

  if( PMComplexMat_ != NULL ){
    delete PMComplexMat_;
  }

  if( PMRealUnsymMat_ != NULL ){
    delete PMRealUnsymMat_;
  }

  if( PMComplexUnsymMat_ != NULL ){
    delete PMComplexUnsymMat_;
  }


	if( gridSuperLUReal_ != NULL ){
		delete gridSuperLUReal_;
	}

	if( gridSuperLUComplex_ != NULL ){
		delete gridSuperLUComplex_;
	}

	if( gridSelInv_ != NULL ){
		delete gridSelInv_;
	}

  if( gridPole_    != NULL ){
    delete gridPole_;
  }


  // Close the log file
  statusOFS.close();

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::~PPEXSIData  ----- 


void
PPEXSIData::LoadRealSymmetricMatrix	(
    Int           nrows,                        
    Int           nnz,                          
    Int           nnzLocal,                     
    Int           numColLocal,                  
    Int*          colptrLocal,                  
    Int*          rowindLocal,                  
    Real*         HnzvalLocal,                  
    Int           isSIdentity,                  
    Real*         SnzvalLocal,
    Int           verbosity )
{
#ifndef _RELEASE_
  PushCallStack("PPEXSIData::LoadRealSymmetricMatrix");
#endif
  // Clear the previously saved information
  HRealMat_ = DistSparseMatrix<Real>();
  SRealMat_ = DistSparseMatrix<Real>();

  // Data communication
  std::vector<char> sstr;
  Int sizeStm;
  if( MYROW( gridPole_ ) == 0 ){
    std::stringstream sstm;

    HRealMat_.size        = nrows;
    HRealMat_.nnz         = nnz;
    HRealMat_.nnzLocal    = nnzLocal;
    // The first row processor does not need extra copies of the index /
    // value of the matrix. 
    HRealMat_.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
    HRealMat_.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
    // H value
    HRealMat_.nzvalLocal  = DblNumVec( nnzLocal,      false, HnzvalLocal );
    HRealMat_.comm        = gridPole_->rowComm;

    // Serialization will copy the values regardless of the ownership
    serialize( HRealMat_, sstm, NO_MASK );

    // S value
    if( isSIdentity ){
      SRealMat_.size = 0;
      SRealMat_.nnz  = 0;
      SRealMat_.nnzLocal = 0;
      SRealMat_.comm = HRealMat_.comm; 
    }
    else{
      CopyPattern( HRealMat_, SRealMat_ );
      SRealMat_.comm = HRealMat_.comm; 
      SRealMat_.nzvalLocal  = DblNumVec( nnzLocal,      false, SnzvalLocal );
      serialize( SRealMat_.nzvalLocal, sstm, NO_MASK );
    }

    sstr.resize( Size( sstm ) );
    sstm.read( &sstr[0], sstr.size() ); 	
    sizeStm = sstr.size();
  }

  MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole_->colComm );

  if( verbosity >= 2 ){
    statusOFS << "sizeStm = " << sizeStm << std::endl;
  }

  if( MYROW( gridPole_ ) != 0 ) sstr.resize( sizeStm );

  MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole_->colComm );

  if( MYROW( gridPole_ ) != 0 ){
    std::stringstream sstm;
    sstm.write( &sstr[0], sizeStm );
    deserialize( HRealMat_, sstm, NO_MASK );
    // Communicator
    HRealMat_.comm = gridPole_->rowComm;
    if( isSIdentity ){
      SRealMat_.size = 0;    // Means S is an identity matrix
      SRealMat_.nnz  = 0;
      SRealMat_.nnzLocal = 0;
      SRealMat_.comm = HRealMat_.comm;
    }
    else{
      CopyPattern( HRealMat_, SRealMat_ );
      SRealMat_.comm = HRealMat_.comm;
      deserialize( SRealMat_.nzvalLocal, sstm, NO_MASK );
    }
  }
  sstr.clear();


  if( verbosity >= 1 ){
    statusOFS << "H.size     = " << HRealMat_.size     << std::endl;
    statusOFS << "H.nnzLocal = " << HRealMat_.nnzLocal << std::endl;
    statusOFS << "S.size     = " << SRealMat_.size     << std::endl;
    statusOFS << "S.nnzLocal = " << SRealMat_.nnzLocal << std::endl;
    statusOFS << std::endl << std::endl;
  }


  // Record the index for the diagonal elements to handle the case if S
  // is identity.
	{
		Int numColLocal      = HRealMat_.colptrLocal.m() - 1;
		Int numColLocalFirst = HRealMat_.size / gridSelInv_->mpisize;
		Int firstCol         = gridSelInv_->mpirank * numColLocalFirst;
		
		diagIdxLocal_.clear();

		for( Int j = 0; j < numColLocal; j++ ){
			Int jcol = firstCol + j + 1;
			for( Int i = HRealMat_.colptrLocal(j)-1; 
				 	 i < HRealMat_.colptrLocal(j+1)-1; i++ ){
				Int irow = HRealMat_.rowindLocal(i);
				if( irow == jcol ){
					diagIdxLocal_.push_back( i );
				}
			}
		} // for (j)
	}

  isMatrixLoaded_ = true;

#ifndef _RELEASE_
  PopCallStack();
#endif

  return ;
}    	// -----  end of method PPEXSIData::LoadRealSymmetricMatrix  ----- 




void
PPEXSIData::LoadRealUnsymmetricMatrix	(
    Int           nrows,                        
    Int           nnz,                          
    Int           nnzLocal,                     
    Int           numColLocal,                  
    Int*          colptrLocal,                  
    Int*          rowindLocal,                  
    Real*         HnzvalLocal,                  
    Int           isSIdentity,                  
    Real*         SnzvalLocal,
    Int           verbosity )
{
#ifndef _RELEASE_
  PushCallStack("PPEXSIData::LoadRealUnsymmetricMatrix");
#endif
  // Clear the previously saved information
  HRealMat_ = DistSparseMatrix<Real>();
  SRealMat_ = DistSparseMatrix<Real>();

  // Data communication
  std::vector<char> sstr;
  Int sizeStm;
  if( MYROW( gridPole_ ) == 0 ){
    std::stringstream sstm;

    HRealMat_.size        = nrows;
    HRealMat_.nnz         = nnz;
    HRealMat_.nnzLocal    = nnzLocal;
    // The first row processor does not need extra copies of the index /
    // value of the matrix. 
    HRealMat_.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
    HRealMat_.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
    // H value
    HRealMat_.nzvalLocal  = DblNumVec( nnzLocal,      false, HnzvalLocal );
    HRealMat_.comm        = gridPole_->rowComm;

    // Serialization will copy the values regardless of the ownership
    serialize( HRealMat_, sstm, NO_MASK );

    // S value
    if( isSIdentity ){
      SRealMat_.size = 0;
      SRealMat_.nnz  = 0;
      SRealMat_.nnzLocal = 0;
      SRealMat_.comm = HRealMat_.comm; 
    }
    else{
      CopyPattern( HRealMat_, SRealMat_ );
      SRealMat_.comm = HRealMat_.comm; 
      SRealMat_.nzvalLocal  = DblNumVec( nnzLocal,      false, SnzvalLocal );
      serialize( SRealMat_.nzvalLocal, sstm, NO_MASK );
    }

    sstr.resize( Size( sstm ) );
    sstm.read( &sstr[0], sstr.size() ); 	
    sizeStm = sstr.size();
  }

  MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole_->colComm );

  if( verbosity >= 2 ){
    statusOFS << "sizeStm = " << sizeStm << std::endl;
  }

  if( MYROW( gridPole_ ) != 0 ) sstr.resize( sizeStm );

  MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole_->colComm );

  if( MYROW( gridPole_ ) != 0 ){
    std::stringstream sstm;
    sstm.write( &sstr[0], sizeStm );
    deserialize( HRealMat_, sstm, NO_MASK );
    // Communicator
    HRealMat_.comm = gridPole_->rowComm;
    if( isSIdentity ){
      SRealMat_.size = 0;    // Means S is an identity matrix
      SRealMat_.nnz  = 0;
      SRealMat_.nnzLocal = 0;
      SRealMat_.comm = HRealMat_.comm;
    }
    else{
      CopyPattern( HRealMat_, SRealMat_ );
      SRealMat_.comm = HRealMat_.comm;
      deserialize( SRealMat_.nzvalLocal, sstm, NO_MASK );
    }
  }
  sstr.clear();


  if( verbosity >= 1 ){
    statusOFS << "H.size     = " << HRealMat_.size     << std::endl;
    statusOFS << "H.nnzLocal = " << HRealMat_.nnzLocal << std::endl;
    statusOFS << "S.size     = " << SRealMat_.size     << std::endl;
    statusOFS << "S.nnzLocal = " << SRealMat_.nnzLocal << std::endl;
    statusOFS << std::endl << std::endl;
  }


  // Record the index for the diagonal elements to handle the case if S
  // is identity.
	{
		Int numColLocal      = HRealMat_.colptrLocal.m() - 1;
		Int numColLocalFirst = HRealMat_.size / gridSelInv_->mpisize;
		Int firstCol         = gridSelInv_->mpirank * numColLocalFirst;
		
		diagIdxLocal_.clear();

		for( Int j = 0; j < numColLocal; j++ ){
			Int jcol = firstCol + j + 1;
			for( Int i = HRealMat_.colptrLocal(j)-1; 
				 	 i < HRealMat_.colptrLocal(j+1)-1; i++ ){
				Int irow = HRealMat_.rowindLocal(i);
				if( irow == jcol ){
					diagIdxLocal_.push_back( i );
				}
			}
		} // for (j)
	}

  isMatrixLoaded_ = true;

#ifndef _RELEASE_
  PopCallStack();
#endif

  return ;
}    	// -----  end of method PPEXSIData::LoadRealUnsymmetricMatrix  ----- 







void
PPEXSIData::SymbolicFactorizeRealSymmetricMatrix	(
    std::string                    ColPerm,
    Int                            numProcSymbFact,
    Int                            verbosity )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::SymbolicFactorizeRealSymmetricMatrix");
#endif
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealSymmetricMatrix first." << std::endl;
#ifdef USE_ABORT
    abort();
#endif
throw std::runtime_error( msg.str().c_str() );
  }
  
  {
    if( verbosity >= 1 ){
      statusOFS << "Symbolic factorization for the real matrix."  << std::endl;
    }

    SuperLUMatrix<Real>&    luMat     = *luRealMat_;
    PMatrix<Real>&          PMloc     = *PMRealMat_;
    SuperNodeType&          super     = superReal_;

    // Clear the matrices first
    luMat = SuperLUMatrix<Real>();
    PMloc = PMatrix<Real>();


    luOpt_.ColPerm = ColPerm;
    luOpt_.numProcSymbFact = numProcSymbFact;
    selinvOpt_.maxPipelineDepth = -1;

    luMat.Setup( *gridSuperLUReal_, luOpt_ );  // SuperLU matrix.

    DistSparseMatrix<Real> AMat;
    CopyPattern( HRealMat_, AMat );

    SetValue( AMat.nzvalLocal, D_ZERO );          // Symbolic factorization does not need value

    Real timeSta, timeEnd;
    GetTime( timeSta );
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS << "Time for SuperMatrix conversion is " <<
        timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    GetTime( timeSta );
    luMat.SymbolicFactorize();
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for symbolic factorization is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    luMat.SymbolicToSuperNode( super );
    
    PMloc.Setup( gridSelInv_, &super , &selinvOpt_, &luOpt_ );
    GetTime( timeSta );
    luMat.LUstructToPMatrix( PMloc );
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for converting LUstruct to PMatrix is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    GetTime( timeSta );
    PMloc.ConstructCommunicationPattern();
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for constructing communication pattern is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }

    // Compute the number of nonzeros from PMatrix
    if( verbosity >= 1 ) {
      Int nnzLocal = PMloc.NnzLocal();
      statusOFS << "Number of local nonzeros (L+U) = " << nnzLocal << std::endl;
      LongInt nnz  = PMloc.Nnz();
      statusOFS << "Number of nonzeros (L+U)       = " << nnz << std::endl;
    }

    // Get ready for the factorization for another matrix using the same
    // sparsity pattern
    luMat.DestroyAOnly();

    if( verbosity >= 2 ){
      statusOFS << "perm: "    << std::endl << super.perm     << std::endl;
      statusOFS << "permInv: " << std::endl << super.permInv  << std::endl;
      statusOFS << "superIdx:" << std::endl << super.superIdx << std::endl;
      statusOFS << "superPtr:" << std::endl << super.superPtr << std::endl; 
    }
  }

  isRealSymmetricSymbolicFactorized_ = true;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::SymbolicFactorizeRealSymmetricMatrix  ----- 




void
PPEXSIData::SymbolicFactorizeRealUnsymmetricMatrix	(
    std::string                    ColPerm,
    std::string                    RowPerm,
    Int                            numProcSymbFact,
    Int                            Transpose,
    double*                        AnzvalLocal,                  
    Int                            verbosity )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::SymbolicFactorizeRealUnsymmetricMatrix");
#endif
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealUnsymmetricMatrix first." << std::endl;
#ifdef USE_ABORT
    abort();
#endif
throw std::runtime_error( msg.str().c_str() );
  }
  
  {
    if( verbosity >= 1 ){
      statusOFS << "Symbolic factorization for the real matrix."  << std::endl;
    }

    SuperLUMatrix<Real>&    luMat     = *luRealMat_;
    PMatrixUnsym<Real>&     PMloc     = *PMRealUnsymMat_;
    SuperNodeType&          super     = superReal_;

    // Clear the matrices first
    luMat = SuperLUMatrix<Real>();
    PMloc = PMatrixUnsym<Real>();


    luOpt_.ColPerm = ColPerm;
    luOpt_.RowPerm = RowPerm;
    luOpt_.numProcSymbFact = numProcSymbFact;
    selinvOpt_.maxPipelineDepth = -1;
    luOpt_.Symmetric = 0;
    luOpt_.Transpose = Transpose;

    luMat.Setup( *gridSuperLUReal_, luOpt_ );  // SuperLU matrix.

    DistSparseMatrix<Real> AMat;
    CopyPattern( HRealMat_, AMat );

    if(luOpt_.RowPerm == "LargeDiag"){
      if(AnzvalLocal != NULL){
        //NOT TRUE
        //SetValue( AMat.nzvalLocal, Z_ZERO );          // Symbolic factorization does not need value
        blas::Copy( AMat.nnzLocal, AnzvalLocal, 1, AMat.nzvalLocal.Data(), 1 );
      }
      else{
        std::ostringstream msg;
        msg  << std::endl
          << "LargeDiag requires the non zero values to be provided." << std::endl;
          throw std::runtime_error( msg.str().c_str() );

      }
    }



    Real timeSta, timeEnd;
    GetTime( timeSta );
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS << "Time for SuperMatrix conversion is " <<
        timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    GetTime( timeSta );
    luMat.SymbolicFactorize();
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for symbolic factorization is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    luMat.SymbolicToSuperNode( super );
    
    PMloc.Setup( gridSelInv_, &super , &selinvOpt_, &luOpt_ );
    GetTime( timeSta );
    luMat.LUstructToPMatrix( PMloc );
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for converting LUstruct to PMatrix is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    GetTime( timeSta );
    PMloc.ConstructCommunicationPattern();
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for constructing communication pattern is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }

    // Compute the number of nonzeros from PMatrix
    if( verbosity >= 1 ) {
      Int nnzLocal = PMloc.NnzLocal();
      statusOFS << "Number of local nonzeros (L+U) = " << nnzLocal << std::endl;
      LongInt nnz  = PMloc.Nnz();
      statusOFS << "Number of nonzeros (L+U)       = " << nnz << std::endl;
    }

    // Get ready for the factorization for another matrix using the same
    // sparsity pattern
    luMat.DestroyAOnly();

    if( verbosity >= 2 ){
      statusOFS << "perm: "    << std::endl << super.perm     << std::endl;
      statusOFS << "permInv: " << std::endl << super.permInv  << std::endl;
      statusOFS << "superIdx:" << std::endl << super.superIdx << std::endl;
      statusOFS << "superPtr:" << std::endl << super.superPtr << std::endl; 
    }
  }

  isRealUnsymmetricSymbolicFactorized_ = true;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::SymbolicFactorizeRealUnsymmetricMatrix  ----- 









void
PPEXSIData::SymbolicFactorizeComplexSymmetricMatrix	(
    std::string                    ColPerm,
    Int                            numProcSymbFact,
    Int                            verbosity )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::SymbolicFactorizeComplexSymmetricMatrix");
#endif
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealSymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }
  
  // Complex matrices
  {
    if( verbosity >= 1 ){
      statusOFS << "Symbolic factorization for the complex matrix."  << std::endl;
    }

    SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
    PMatrix<Complex>&          PMloc     = *PMComplexMat_;
    SuperNodeType&             super     = superComplex_;

    // Clear the matrices first
    luMat = SuperLUMatrix<Complex>();
    PMloc = PMatrix<Complex>();


    luOpt_.ColPerm = ColPerm;
    luOpt_.numProcSymbFact = numProcSymbFact;
    selinvOpt_.maxPipelineDepth = -1;

    luMat.Setup( *gridSuperLUComplex_, luOpt_ );  // SuperLU matrix.

    DistSparseMatrix<Complex> AMat;
    CopyPattern( HRealMat_, AMat );

    SetValue( AMat.nzvalLocal, Z_ZERO );          // Symbolic factorization does not need value

    Real timeSta, timeEnd;
    GetTime( timeSta );
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS << "Time for SuperMatrix conversion is " <<
        timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    GetTime( timeSta );
    luMat.SymbolicFactorize();
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for symbolic factorization is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    luMat.SymbolicToSuperNode( super );
    
    PMloc.Setup( gridSelInv_, &super , &selinvOpt_, &luOpt_ );
    GetTime( timeSta );
    luMat.LUstructToPMatrix( PMloc );
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for converting LUstruct to PMatrix is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    GetTime( timeSta );
    PMloc.ConstructCommunicationPattern();
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for constructing communication pattern is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }

    // Compute the number of nonzeros from PMatrix
    if( verbosity >= 1 ) {
      Int nnzLocal = PMloc.NnzLocal();
      statusOFS << "Number of local nonzeros (L+U) = " << nnzLocal << std::endl;
      LongInt nnz  = PMloc.Nnz();
      statusOFS << "Number of nonzeros (L+U)       = " << nnz << std::endl;
    }

    // Get ready for the factorization for another matrix using the same
    // sparsity pattern
    luMat.DestroyAOnly();

    if( verbosity >= 2 ){
      statusOFS << "perm: "    << std::endl << super.perm     << std::endl;
      statusOFS << "permInv: " << std::endl << super.permInv  << std::endl;
      statusOFS << "superIdx:" << std::endl << super.superIdx << std::endl;
      statusOFS << "superPtr:" << std::endl << super.superPtr << std::endl; 
    }
  }


  isComplexSymmetricSymbolicFactorized_ = true;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::SymbolicFactorizeComplexSymmetricMatrix  ----- 

void
PPEXSIData::SymbolicFactorizeComplexUnsymmetricMatrix	(
    std::string                    ColPerm,
    std::string                    RowPerm,
    Int                            numProcSymbFact,
    Int                            Transpose,
    double*                        AnzvalLocal,                  
    Int                            verbosity )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::SymbolicFactorizeComplexUnsymmetricMatrix");
#endif
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealUnsymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }
  
  // Complex matrices
  {
    if( verbosity >= 1 ){
      statusOFS << "Symbolic factorization for the complex matrix."  << std::endl;
    }

    SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
    PMatrixUnsym<Complex>&     PMloc     = *PMComplexUnsymMat_;
    SuperNodeType&             super     = superComplex_;

    // Clear the matrices first
    luMat = SuperLUMatrix<Complex>();
    PMloc = PMatrixUnsym<Complex>();


    luOpt_.ColPerm = ColPerm;
    luOpt_.RowPerm = RowPerm;
    luOpt_.numProcSymbFact = numProcSymbFact;
    selinvOpt_.maxPipelineDepth = -1;
    luOpt_.Symmetric = 0;
    luOpt_.Transpose = Transpose;

    luMat.Setup( *gridSuperLUComplex_, luOpt_ );  // SuperLU matrix.

    DistSparseMatrix<Complex> AMat;
    CopyPattern( HRealMat_, AMat );

    if(luOpt_.RowPerm == "LargeDiag"){
      if(AnzvalLocal != NULL){
        //NOT TRUE
        //SetValue( AMat.nzvalLocal, Z_ZERO );          // Symbolic factorization does not need value
        blas::Copy( 2*AMat.nnzLocal, AnzvalLocal, 1, 
            reinterpret_cast<double*>(AMat.nzvalLocal.Data()), 1 );
      }
      else{
        std::ostringstream msg;
        msg  << std::endl
          << "LargeDiag requires the non zero values to be provided." << std::endl;
          throw std::runtime_error( msg.str().c_str() );

      }
    }

    Real timeSta, timeEnd;
    GetTime( timeSta );
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS << "Time for SuperMatrix conversion is " <<
        timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    GetTime( timeSta );
    luMat.SymbolicFactorize();
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for symbolic factorization is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    luMat.SymbolicToSuperNode( super );
    
    PMloc.Setup( gridSelInv_, &super , &selinvOpt_, &luOpt_ );
    GetTime( timeSta );
    luMat.LUstructToPMatrix( PMloc );
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for converting LUstruct to PMatrix is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }
    GetTime( timeSta );
    PMloc.ConstructCommunicationPattern();
    GetTime( timeEnd );
    if( verbosity >= 1 ){
      statusOFS 
        << "Time for constructing communication pattern is " 
        << timeEnd - timeSta << " [s]" << std::endl << std::endl;
    }

    // Compute the number of nonzeros from PMatrix
    if( verbosity >= 1 ) {
      Int nnzLocal = PMloc.NnzLocal();
      statusOFS << "Number of local nonzeros (L+U) = " << nnzLocal << std::endl;
      LongInt nnz  = PMloc.Nnz();
      statusOFS << "Number of nonzeros (L+U)       = " << nnz << std::endl;
    }

    // Get ready for the factorization for another matrix using the same
    // sparsity pattern
    luMat.DestroyAOnly();

    if( verbosity >= 2 ){
      statusOFS << "perm: "    << std::endl << super.perm     << std::endl;
      statusOFS << "permInv: " << std::endl << super.permInv  << std::endl;
      statusOFS << "superIdx:" << std::endl << super.superIdx << std::endl;
      statusOFS << "superPtr:" << std::endl << super.superPtr << std::endl; 
    }
  }


  isComplexUnsymmetricSymbolicFactorized_ = true;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::SymbolicFactorizeComplexUnsymmetricMatrix  ----- 



void 
PPEXSIData::SelInvRealSymmetricMatrix(
    double*           AnzvalLocal,                  
    Int               verbosity,
    double*           AinvnzvalLocal )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::SelInvRealSymmetricMatrix");
#endif
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealSymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }

  if( isRealSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeRealSymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }
  
  // Only the processor group corresponding to the first pole participate
  if( MYROW( gridPole_ ) == 0 ){

    DistSparseMatrix<Real>& AMat      = shiftRealMat_;
    DistSparseMatrix<Real>& AinvMat   = shiftInvRealMat_;
    SuperLUMatrix<Real>&    luMat     = *luRealMat_;
    PMatrix<Real>&          PMloc     = *PMRealMat_;

    // Copy the pattern
    CopyPattern( HRealMat_, AMat );

    blas::Copy( AMat.nnzLocal, AnzvalLocal, 1, AMat.nzvalLocal.Data(), 1 );
    
    if( verbosity >= 2 ){
      statusOFS << "Before DistSparseMatrixToSuperMatrixNRloc." << std::endl;
    }
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
    if( verbosity >= 2 ){
      statusOFS << "After DistSparseMatrixToSuperMatrixNRloc." << std::endl;
    }

    Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

    GetTime( timeTotalFactorizationSta );

    // Data redistribution
    if( verbosity >= 2 ){
      statusOFS << "Before Distribute." << std::endl;
    }
    luMat.Distribute();
    if( verbosity >= 2 ){
      statusOFS << "After Distribute." << std::endl;
    }

    // Numerical factorization
    if( verbosity >= 2 ){
      statusOFS << "Before NumericalFactorize." << std::endl;
    }
    luMat.NumericalFactorize();
    if( verbosity >= 2 ){
      statusOFS << "After NumericalFactorize." << std::endl;
    }
    luMat.DestroyAOnly();

    GetTime( timeTotalFactorizationEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 
    }

    Real timeTotalSelInvSta, timeTotalSelInvEnd;
    GetTime( timeTotalSelInvSta );

    luMat.LUstructToPMatrix( PMloc );

    PMloc.PreSelInv();
    
    PMloc.SelInv();

    GetTime( timeTotalSelInvEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total selected inversion is " <<
        timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
    }


    Real timePostProcessingSta, timePostProcessingEnd;

    GetTime( timePostProcessingSta );

    PMloc.PMatrixToDistSparseMatrix( AMat, AinvMat );

    GetTime( timePostProcessingEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for postprocessing is " <<
        timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
    }

    // Return the data to AinvnzvalLocal
    blas::Copy( AMat.nnzLocal, AinvMat.nzvalLocal.Data(), 1, AinvnzvalLocal, 1 );
  }

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::SelInvRealSymmetricMatrix  ----- 

void 
PPEXSIData::SelInvRealUnsymmetricMatrix(
    double*           AnzvalLocal,                  
    Int               verbosity,
    double*           AinvnzvalLocal )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::SelInvRealUnsymmetricMatrix");
#endif
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealUnsymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }

  if( isRealUnsymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeRealUnsymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }
  
  // Only the processor group corresponding to the first pole participate
  if( MYROW( gridPole_ ) == 0 ){

    DistSparseMatrix<Real>& AMat      = shiftRealMat_;
    DistSparseMatrix<Real>& AinvMat   = shiftInvRealMat_;
    SuperLUMatrix<Real>&    luMat     = *luRealMat_;
    PMatrixUnsym<Real>&          PMloc     = *PMRealUnsymMat_;

    // Copy the pattern
    CopyPattern( HRealMat_, AMat );

    blas::Copy( AMat.nnzLocal, AnzvalLocal, 1, AMat.nzvalLocal.Data(), 1 );
    
    if( verbosity >= 2 ){
      statusOFS << "Before DistSparseMatrixToSuperMatrixNRloc." << std::endl;
    }
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
    if( verbosity >= 2 ){
      statusOFS << "After DistSparseMatrixToSuperMatrixNRloc." << std::endl;
    }

    Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

    GetTime( timeTotalFactorizationSta );

    // Data redistribution
    if( verbosity >= 2 ){
      statusOFS << "Before Distribute." << std::endl;
    }
    luMat.Distribute();
    if( verbosity >= 2 ){
      statusOFS << "After Distribute." << std::endl;
    }

    // Numerical factorization
    if( verbosity >= 2 ){
      statusOFS << "Before NumericalFactorize." << std::endl;
    }
    luMat.NumericalFactorize();
    if( verbosity >= 2 ){
      statusOFS << "After NumericalFactorize." << std::endl;
    }
    luMat.DestroyAOnly();

    GetTime( timeTotalFactorizationEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 
    }

    Real timeTotalSelInvSta, timeTotalSelInvEnd;
    GetTime( timeTotalSelInvSta );

    luMat.LUstructToPMatrix( PMloc );

    PMloc.PreSelInv();
    
    PMloc.SelInv();

    GetTime( timeTotalSelInvEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total selected inversion is " <<
        timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
    }


    Real timePostProcessingSta, timePostProcessingEnd;

    GetTime( timePostProcessingSta );

    PMloc.PMatrixToDistSparseMatrix( AMat, AinvMat );

    GetTime( timePostProcessingEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for postprocessing is " <<
        timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
    }

    // Return the data to AinvnzvalLocal
    blas::Copy( AMat.nnzLocal, AinvMat.nzvalLocal.Data(), 1, AinvnzvalLocal, 1 );
  }

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::SelInvRealUnsymmetricMatrix  ----- 


void 
PPEXSIData::SelInvComplexSymmetricMatrix(
    double*           AnzvalLocal,                  
    Int               verbosity,
    double*           AinvnzvalLocal )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::SelInvComplexSymmetricMatrix");
#endif
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealSymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }

  if( isComplexSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeComplexSymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }
  
  // Only the processor group corresponding to the first pole participate
  if( MYROW( gridPole_ ) == 0 ){

    DistSparseMatrix<Complex>& AMat      = shiftComplexMat_;
    DistSparseMatrix<Complex>& AinvMat   = shiftInvComplexMat_;
    SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
    PMatrix<Complex>&          PMloc     = *PMComplexMat_;

    // Copy the pattern
    CopyPattern( HRealMat_, AMat );

    blas::Copy( 2*AMat.nnzLocal, AnzvalLocal, 1, 
       reinterpret_cast<double*>(AMat.nzvalLocal.Data()), 1 );
    
    if( verbosity >= 2 ){
      statusOFS << "Before DistSparseMatrixToSuperMatrixNRloc." << std::endl;
    }
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
    if( verbosity >= 2 ){
      statusOFS << "After DistSparseMatrixToSuperMatrixNRloc." << std::endl;
    }

    Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

    GetTime( timeTotalFactorizationSta );

    // Data redistribution
    if( verbosity >= 2 ){
      statusOFS << "Before Distribute." << std::endl;
    }
    luMat.Distribute();
    if( verbosity >= 2 ){
      statusOFS << "After Distribute." << std::endl;
    }

    // Numerical factorization
    if( verbosity >= 2 ){
      statusOFS << "Before NumericalFactorize." << std::endl;
    }
    luMat.NumericalFactorize();
    if( verbosity >= 2 ){
      statusOFS << "After NumericalFactorize." << std::endl;
    }
    luMat.DestroyAOnly();

    GetTime( timeTotalFactorizationEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 
    }

    Real timeTotalSelInvSta, timeTotalSelInvEnd;
    GetTime( timeTotalSelInvSta );

    luMat.LUstructToPMatrix( PMloc );

    PMloc.PreSelInv();
    
    PMloc.SelInv();

    GetTime( timeTotalSelInvEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total selected inversion is " <<
        timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
    }


    Real timePostProcessingSta, timePostProcessingEnd;

    GetTime( timePostProcessingSta );

    PMloc.PMatrixToDistSparseMatrix( AMat, AinvMat );

    GetTime( timePostProcessingEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for postprocessing is " <<
        timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
    }

    // Return the data to AinvnzvalLocal
    blas::Copy( 2*AMat.nnzLocal, reinterpret_cast<double*>(AinvMat.nzvalLocal.Data()), 1, 
        AinvnzvalLocal, 1 );
  }

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::SelInvComplexSymmetricMatrix  ----- 

void 
PPEXSIData::SelInvComplexUnsymmetricMatrix(
    double*           AnzvalLocal,                  
    Int               verbosity,
    double*           AinvnzvalLocal )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::SelInvComplexUnsymmetricMatrix");
#endif
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealUnsymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }

  if( isComplexUnsymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeComplexUnsymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }
  
  // Only the processor group corresponding to the first pole participate
  if( MYROW( gridPole_ ) == 0 ){

    DistSparseMatrix<Complex>& AMat      = shiftComplexMat_;
    DistSparseMatrix<Complex>& AinvMat   = shiftInvComplexMat_;
    SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
    PMatrixUnsym<Complex>&          PMloc     = *PMComplexUnsymMat_;

    // Copy the pattern
    CopyPattern( HRealMat_, AMat );

    blas::Copy( 2*AMat.nnzLocal, AnzvalLocal, 1, 
       reinterpret_cast<double*>(AMat.nzvalLocal.Data()), 1 );
    
    if( verbosity >= 2 ){
      statusOFS << "Before DistSparseMatrixToSuperMatrixNRloc." << std::endl;
    }
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
    if( verbosity >= 2 ){
      statusOFS << "After DistSparseMatrixToSuperMatrixNRloc." << std::endl;
    }

    Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

    GetTime( timeTotalFactorizationSta );

    // Data redistribution
    if( verbosity >= 2 ){
      statusOFS << "Before Distribute." << std::endl;
    }
    luMat.Distribute();
    if( verbosity >= 2 ){
      statusOFS << "After Distribute." << std::endl;
    }

    // Numerical factorization
    if( verbosity >= 2 ){
      statusOFS << "Before NumericalFactorize." << std::endl;
    }
    luMat.NumericalFactorize();
    if( verbosity >= 2 ){
      statusOFS << "After NumericalFactorize." << std::endl;
    }
    luMat.DestroyAOnly();

    GetTime( timeTotalFactorizationEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 
    }

    Real timeTotalSelInvSta, timeTotalSelInvEnd;
    GetTime( timeTotalSelInvSta );

    luMat.LUstructToPMatrix( PMloc );

    PMloc.PreSelInv();
    
    PMloc.SelInv();

    GetTime( timeTotalSelInvEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total selected inversion is " <<
        timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
    }


    Real timePostProcessingSta, timePostProcessingEnd;

    GetTime( timePostProcessingSta );

    PMloc.PMatrixToDistSparseMatrix( AMat, AinvMat );

    GetTime( timePostProcessingEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for postprocessing is " <<
        timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
    }

    // Return the data to AinvnzvalLocal
    blas::Copy( 2*AMat.nnzLocal, reinterpret_cast<double*>(AinvMat.nzvalLocal.Data()), 1, 
        AinvnzvalLocal, 1 );
  }

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::SelInvComplexUnsymmetricMatrix  ----- 



void PPEXSIData::CalculateNegativeInertiaReal(
		const std::vector<Real>&       shiftVec, 
		std::vector<Real>&             inertiaVec,
    Int                            verbosity ){
#ifndef _RELEASE_
  PushCallStack("PPEXSIData::CalculateNegativeInertiaReal");
#endif

	// *********************************************************************
	// Initialize
	// *********************************************************************
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealSymmetricMatrix first." << std::endl;
#ifdef USE_ABORT
    abort();
#endif
throw std::runtime_error( msg.str().c_str() );
  }

  if( isRealSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeRealSymmetricMatrix first." << std::endl;
#ifdef USE_ABORT
    abort();
#endif
throw std::runtime_error( msg.str().c_str() );
  }


  // Rename for convenience.
  // The symbolic information should have already been there.
  DistSparseMatrix<Real>&  HMat     = HRealMat_;
  DistSparseMatrix<Real>&  SMat     = SRealMat_;
  DistSparseMatrix<Real>& AMat      = shiftRealMat_;  // A = H - \lambda  S
  SuperLUMatrix<Real>&    luMat     = *luRealMat_;
  PMatrix<Real>&          PMloc     = *PMRealMat_;

  CopyPattern( HMat, AMat );

	Real timeShiftSta, timeShiftEnd;
	
	Int numShift = shiftVec.size();
	std::vector<Real>  inertiaVecLocal(numShift);
	inertiaVec.resize(numShift);
	for(Int l = 0; l < numShift; l++){
		inertiaVecLocal[l] = 0.0;
		inertiaVec[l]      = 0.0;
	}

	for(Int l = 0; l < numShift; l++){
		if( MYROW( gridPole_ ) == PROW( l, gridPole_ ) ){

			GetTime( timeShiftSta );

      if( verbosity >= 1 ){
        statusOFS << "Shift " << l << " = " << shiftVec[l] 
          << " processing..." << std::endl;
      }

			if( SMat.size != 0 ){
				// S is not an identity matrix
				for( Int i = 0; i < HMat.nnzLocal; i++ ){
					AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - shiftVec[l] * SMat.nzvalLocal(i);
				}
			}
			else{
				// S is an identity matrix
				for( Int i = 0; i < HMat.nnzLocal; i++ ){
					AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
				}

				for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
					AMat.nzvalLocal( diagIdxLocal_[i] ) -= shiftVec[l];
				}
			} // if (SMat.size != 0 )


			// *********************************************************************
			// Factorization
			// *********************************************************************
			// Important: the distribution in pzsymbfact is going to mess up the
			// A matrix.  Recompute the matrix A here.
      if( verbosity >= 2 ){
        statusOFS << "Before DistSparseMatrixToSuperMatrixNRloc." << std::endl;
      }
			luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
      if( verbosity >= 2 ){
        statusOFS << "After DistSparseMatrixToSuperMatrixNRloc." << std::endl;
      }

			Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

			GetTime( timeTotalFactorizationSta );

			// Data redistribution
      if( verbosity >= 2 ){
        statusOFS << "Before Distribute." << std::endl;
      }
			luMat.Distribute();
      if( verbosity >= 2 ){
        statusOFS << "After Distribute." << std::endl;
      }

			// Numerical factorization
      if( verbosity >= 2 ){
        statusOFS << "Before NumericalFactorize." << std::endl;
      }
			luMat.NumericalFactorize();
      if( verbosity >= 2 ){
        statusOFS << "After NumericalFactorize." << std::endl;
      }
			luMat.DestroyAOnly();

			GetTime( timeTotalFactorizationEnd );

      if( verbosity >= 1 ){
        statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 
      }

			// *********************************************************************
			// Compute inertia
			// *********************************************************************
			Real timeInertiaSta, timeInertiaEnd;
			GetTime( timeInertiaSta );

			luMat.LUstructToPMatrix( PMloc );

			// Compute the negative inertia of the matrix.
			PMloc.GetNegativeInertia( inertiaVecLocal[l] );

			GetTime( timeInertiaEnd );

      if( verbosity >= 1 ){
        statusOFS << "Time for computing the inertia is " <<
          timeInertiaEnd  - timeInertiaSta << " [s]" << std::endl;
      }


		} // if I am in charge of this shift
	} // for(l)

	// Collect all the negative inertia together
	mpi::Allreduce( &inertiaVecLocal[0], &inertiaVec[0], numShift, 
			MPI_SUM, gridPole_->colComm );


#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::CalculateNegativeInertiaReal ----- 


// Main subroutine for the electronic structure calculation
void PPEXSIData::CalculateFermiOperatorReal(
    Int   numPole, 
    Real  temperature,
    Real  gap,
    Real  deltaE,
    Real  mu,
    Real  numElectronExact,
    Real  numElectronTolerance,
    Int   verbosity,
    Real& numElectron,
    Real& numElectronDrvMu ){

#ifndef _RELEASE_
  PushCallStack("PPEXSIData::CalculateFermiOperatorReal");
#endif
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealSymmetricMatrix first." << std::endl;
#ifdef USE_ABORT
    abort();
#endif
throw std::runtime_error( msg.str().c_str() );
  }

  if( isComplexSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeComplexSymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }

  // *********************************************************************
  // Check the input parameters
  // *********************************************************************
  if( numPole % 2 != 0 ){
#ifdef USE_ABORT
    abort();
#endif
throw std::logic_error( "Must be even number of poles!" );
  }

  // *********************************************************************
  // Initialize
  // *********************************************************************
  // Rename for convenience
  DistSparseMatrix<Real>&  HMat        = HRealMat_;
  DistSparseMatrix<Real>&  SMat        = SRealMat_;

  DistSparseMatrix<Real>& rhoMat       = rhoRealMat_;     
  DistSparseMatrix<Real>& rhoDrvMuMat  = rhoDrvMuRealMat_;
  DistSparseMatrix<Real>& rhoDrvTMat   = rhoDrvTRealMat_;
  DistSparseMatrix<Real>& hmzMat       = freeEnergyDensityRealMat_;
  DistSparseMatrix<Real>& frcMat       = energyDensityRealMat_;

  DistSparseMatrix<Complex>& AMat      = shiftComplexMat_;
  DistSparseMatrix<Complex>& AinvMat   = shiftInvComplexMat_;

  // The symbolic information should have already been there.
  SuperLUMatrix<Complex>& luMat        = *luComplexMat_;
  PMatrix<Complex>&       PMloc        = *PMComplexMat_;

  // 
  bool isFreeEnergyDensityMatrix = true;
  bool isEnergyDensityMatrix     = true;
  bool isDerivativeTMatrix       = false;

  // Copy the pattern
  CopyPattern( HMat, AMat );
  CopyPattern( HMat, rhoMat );
  CopyPattern( HMat, rhoDrvMuMat );
  if( isFreeEnergyDensityMatrix )
    CopyPattern( HMat, hmzMat );
  if( isEnergyDensityMatrix )
    CopyPattern( HMat, frcMat );
  if( isDerivativeTMatrix )
    CopyPattern( HMat, rhoDrvTMat );

  // Reinitialize the variables
  SetValue( rhoMat.nzvalLocal, 0.0 );
  SetValue( rhoDrvMuMat.nzvalLocal, 0.0 );
  if( isFreeEnergyDensityMatrix )
    SetValue( hmzMat.nzvalLocal, 0.0 );
  if( isEnergyDensityMatrix )
    SetValue( frcMat.nzvalLocal, 0.0 );
  if( isDerivativeTMatrix )
    SetValue( rhoDrvTMat.nzvalLocal, 0.0 );

  // Refine the pole expansion  
  // numPoleInput is the number of poles to be given to other parts of
  // the pole expansion, which is larger than or equal to numPole.
  Int numPoleInput;
  // poleIdx is a vector of size numPole.  Only poles with index in
  // poleIdx are used for actual computation. The rest of the poles are
  // discarded according to tolerance criterion
  //
  //   numElectronTolerance / numElectronExact / numPole
  //
  // FIXME The heuristics should be refined to give error estimate to
  // other quantities such as the energy.
  // FIXME The heuristics part should also be given in a separate
  // routine, and the input of this file does not need mu, gap etc.
  std::vector<Int>  poleIdx(numPole);
  {
    // Setup a grid from (mu - deltaE, mu + deltaE), and measure
    // the error bound in the L^infty sense on this grid.
    //
    // fdGrid:      Exact Fermi-Dirac function evaluated on xGrid
    // fdPoleGrid:  Fermi-Dirac function using pole expansion
    // evaluated on the grid.
    Int numX = 10000;
    std::vector<Real>    xGrid( numX );
    std::vector<Real>    fdGrid( numX );

    Real x0 = mu - deltaE;
    Real x1 = mu + deltaE;
    Real h  = (x1 - x0) / (numX - 1);
    Real ez;
    for( Int i = 0; i < numX; i++ ){
      xGrid[i]  = x0 + i * h;
      if( xGrid[i] - mu >= 0 ){
        ez = std::exp(- (xGrid[i] - mu) / temperature );
        fdGrid[i] = 2.0 * ez / (1.0 + ez);
      }
      else{
        ez = std::exp((xGrid[i] - mu) / temperature );
        fdGrid[i] = 2.0 / (1.0 + ez);
      }
    }


    numPoleInput = numPole;
    Real tol;
    Int  numPoleSignificant;

    Int poleIter = 0;
    do{
      // If the number of significant poles is less than numPole,
      // increase numPoleInput by 2 at a time and redo the
      // computation.
      if( poleIter > 0 )
        numPoleInput += 2;

      zshift_.resize( numPoleInput );
      zweightRho_.resize( numPoleInput );
      GetPoleDensity( &zshift_[0], &zweightRho_[0],
          numPoleInput, temperature, gap, deltaE, mu ); 

      std::vector<Complex>  zshiftTmp( numPoleInput );
      zweightForce_.resize( numPoleInput );
      GetPoleForce( &zshiftTmp[0], &zweightForce_[0],
          numPoleInput, temperature, gap, deltaE, mu ); 


      std::vector<Real>  maxMagPole(numPoleInput);
      for( Int l = 0; l < numPoleInput; l++ ){
        maxMagPole[l] = 0.0;
      }

      // Compute the approximation due to pole expansion, as well as
      // the maximum magnitude of each pole
      Complex cpxmag;
      Real    mag;
      numPoleSignificant = 0;
      tol = numElectronTolerance / numElectronExact / numPoleInput;
      for( Int l = 0; l < numPoleInput; l++ ){
        for( Int i = 0; i < numX; i++ ){
          cpxmag = zweightRho_[l] / ( xGrid[i] - zshift_[l] );
          mag    = cpxmag.imag();
          maxMagPole[l] = ( maxMagPole[l] >= mag ) ?  maxMagPole[l] : mag;
        }
        if( maxMagPole[l] > tol ){
          numPoleSignificant++;
        }	
      } // for (l)

      // Pick the most significant numPole poles and update poleIdx
      // Sort in DESCENDING order
      std::vector<Int>  sortIdx( numPoleInput );
      for( Int i = 0; i < sortIdx.size(); i++ ){
        sortIdx[i]      = i;
      }
      std::sort( sortIdx.begin(), sortIdx.end(), 
          IndexComp<std::vector<Real>& >( maxMagPole ) ) ;
      std::reverse( sortIdx.begin(), sortIdx.end() );

      for( Int l = 0; l < numPole; l++ ){
        poleIdx[l]      = sortIdx[l];
      }


      // Update poleIter
      poleIter++; 
    } while( numPoleSignificant < numPole );


    // Estimate the error of the number of electrons and the band energy
    // by assuming a flat density of states within a interval of size
    // deltaE, i.e. each unit interval contains 
    // HMat.size / deltaE 
    // number of electrons

    std::vector<Real>  fdPoleGrid( numX );
    std::vector<Real>  fdEPoleGrid( numX );
    std::vector<Real>  fdTimesEPoleGrid( numX );
    Real errAbs1, errorAbsMax1;
    Real errAbs2, errorAbsMax2;
    Real errorNumElectron, errorBandEnergy, errorBandEnergy2;
    Complex cpxmag1, cpxmag2;

    errorAbsMax1 = 0.0; 
    errorAbsMax2 = 0.0; 
    for( Int i = 0; i < numX; i++ ){
      fdPoleGrid[i] = 0.0;
      fdEPoleGrid[i] = 0.0;
      for( Int lidx = 0; lidx < numPoleInput; lidx++ ){
        Int l = lidx;
        cpxmag1 = zweightRho_[l] / ( xGrid[i] - zshift_[l] );
        cpxmag2 = zweightForce_[l] / ( xGrid[i] - zshift_[l] );
        fdPoleGrid[i] += cpxmag1.imag();
        fdTimesEPoleGrid[i] += cpxmag1.imag() * xGrid[i];
        fdEPoleGrid[i]      += cpxmag2.imag();
      }
      errAbs1 = std::abs( fdPoleGrid[i] - fdGrid[i] );
      errorAbsMax1 = ( errorAbsMax1 >= errAbs1 ) ? errorAbsMax1 : errAbs1;
      errAbs2 = std::abs( fdEPoleGrid[i] - fdTimesEPoleGrid[i] );
      errorAbsMax2 = ( errorAbsMax2 >= errAbs2 ) ? errorAbsMax2 : errAbs2;
    }

    errorNumElectron = errorAbsMax1 * HMat.size;
    errorBandEnergy  = 0.5 * deltaE * errorAbsMax1 * HMat.size;
    errorBandEnergy2 = errorAbsMax2 * HMat.size;

    if( verbosity >= 1 ){
      statusOFS 
        << std::endl 
        << "Estimated error of pole expansion by assuming a flat spectrum."
        << std::endl;

      // The estimation of energy using the difference of DM and EDM
      // seems to be more reliable
      Print( statusOFS, "Error of num electron             = ", errorNumElectron );
//      Print( statusOFS, "Error of band energy (DM only)    = ", errorBandEnergy );
      Print( statusOFS, "Error of band energy (DM and EDM) = ", errorBandEnergy2 );
      Print( statusOFS, "Required accuracy (num electron)  = ", numElectronTolerance );

      statusOFS << std::endl;

      if( errorNumElectron > numElectronTolerance ){
        statusOFS << "WARNING!!! " 
          << "Pole expansion may not be accurate enough to reach numElectronTolerance. " << std::endl
          << "Try to increase numPole or increase numElectronTolerance." << std::endl << std::endl;
      }
      statusOFS << "numPoleInput =" << numPoleInput << std::endl;
      statusOFS << "numPoleSignificant = " << numPoleSignificant << std::endl;
    }
  }

  // Initialize the number of electrons
  numElectron  = 0.0;

  //Initialize the pole expansion
  zweightRhoDrvMu_.resize( numPoleInput );

  GetPoleDensityDrvMu( &zshift_[0], &zweightRhoDrvMu_[0],
      numPoleInput, temperature, gap, deltaE, mu ); 

  if( isFreeEnergyDensityMatrix ){
    std::vector<Complex>  zshiftTmp( numPoleInput );
    zweightHelmholtz_.resize( numPoleInput );
    GetPoleHelmholtz( &zshiftTmp[0], &zweightHelmholtz_[0], 
        numPoleInput, temperature, gap, deltaE, mu ); 
  }

  if( isEnergyDensityMatrix ){
    std::vector<Complex>  zshiftTmp( numPoleInput );
    zweightForce_.resize( numPoleInput );
    GetPoleForce( &zshiftTmp[0], &zweightForce_[0],
        numPoleInput, temperature, gap, deltaE, mu ); 
  }

  if( isDerivativeTMatrix ){
    std::vector<Complex>  zshiftTmp( numPoleInput );
    zweightRhoDrvT_.resize( numPoleInput );
    GetPoleDensityDrvT( &zshiftTmp[0], &zweightRhoDrvT_[0],
        numPoleInput, temperature, gap, deltaE, mu ); 
  }

  if( verbosity >= 2 ){
    statusOFS << "zshift" << std::endl << zshift_ << std::endl;
    statusOFS << "zweightRho" << std::endl << zweightRho_ << std::endl;
  }

  // for each pole, perform LDLT factoriation and selected inversion
  Real timePoleSta, timePoleEnd;

  Int numPoleComputed = 0;
  for(Int lidx = 0; lidx < numPole; lidx++){
    if( MYROW( gridPole_ ) == PROW( lidx, gridPole_ ) ){

      Int l = poleIdx[lidx];

      GetTime( timePoleSta );

      if( verbosity >= 1 ){
        statusOFS << "Pole " << lidx << " processing..." << std::endl;
      }
      if( verbosity >= 2 ){
        statusOFS << "zshift           = " << zshift_[l] << std::endl;
        statusOFS	<< "zweightRho       = " << zweightRho_[l] << std::endl;
        statusOFS	<< "zweightRhoDrvMu  = " << zweightRhoDrvMu_[l] << std::endl;
        if( isFreeEnergyDensityMatrix )
          statusOFS << "zweightHelmholtz = " << zweightHelmholtz_[l] << std::endl;
        if( isEnergyDensityMatrix )
          statusOFS << "zweightForce     = " << zweightForce_[l] << std::endl;
        if( isDerivativeTMatrix )
          statusOFS << "zweightRhoDrvT   = " << zweightRhoDrvT_[l] << std::endl;
      }

      {
        AinvMat.Clear();


        numPoleComputed++;

        if( SMat.size != 0 ){
          // S is not an identity matrix
          for( Int i = 0; i < HMat.nnzLocal; i++ ){
            AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift_[l] * SMat.nzvalLocal(i);
          }
        }
        else{
          // S is an identity matrix
          for( Int i = 0; i < HMat.nnzLocal; i++ ){
            AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
          }

          for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
            AMat.nzvalLocal( diagIdxLocal_[i] ) -= zshift_[l];
          }
        } // if (SMat.size != 0 )


        // *********************************************************************
        // Factorization
        // *********************************************************************
        // Important: the distribution in pzsymbfact is going to mess up the
        // A matrix.  Recompute the matrix A here.
        if( verbosity >= 2 ){
          statusOFS << "Before DistSparseMatrixToSuperMatrixNRloc." << std::endl;
        }
        luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
        if( verbosity >= 2 ){
          statusOFS << "After DistSparseMatrixToSuperMatrixNRloc." << std::endl;
        }

        Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

        GetTime( timeTotalFactorizationSta );

        // Data redistribution
        if( verbosity >= 2 ){
          statusOFS << "Before Distribute." << std::endl;
        }
        luMat.Distribute();
        if( verbosity >= 2 ){
          statusOFS << "After Distribute." << std::endl;
        }

        // Numerical factorization
        if( verbosity >= 2 ){
          statusOFS << "Before NumericalFactorize." << std::endl;
        }
        luMat.NumericalFactorize();
        if( verbosity >= 2 ){
          statusOFS << "After NumericalFactorize." << std::endl;
        }
        luMat.DestroyAOnly();

        GetTime( timeTotalFactorizationEnd );

        if( verbosity >= 1 ){
          statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 
        }

        // *********************************************************************
        // Selected inversion
        // *********************************************************************
        Real timeTotalSelInvSta, timeTotalSelInvEnd;
        GetTime( timeTotalSelInvSta );

        luMat.LUstructToPMatrix( PMloc );


        // Collective communication version
        //          PMloc.ConstructCommunicationPattern_Collectives();

        PMloc.PreSelInv();

        // Main subroutine for selected inversion
        //
        // P2p communication version
        PMloc.SelInv();

        // Collective communication version
        //          PMloc.SelInv_Collectives();

        GetTime( timeTotalSelInvEnd );

        if( verbosity >= 1 ){
          statusOFS << "Time for total selected inversion is " <<
            timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
        }

        // *********************************************************************
        // Postprocessing
        // *********************************************************************

        Real timePostProcessingSta, timePostProcessingEnd;

        GetTime( timePostProcessingSta );

        PMloc.PMatrixToDistSparseMatrix( AMat, AinvMat );

        if( verbosity >= 2 ){
          statusOFS << "rhoMat.nnzLocal = " << rhoMat.nnzLocal << std::endl;
          statusOFS << "AinvMat.nnzLocal = " << AinvMat.nnzLocal << std::endl;
        }


        // Update the density matrix. The following lines are equivalent to
        //
        //				for( Int i = 0; i < rhoMat.nnzLocal; i++ ){
        //					rhoMat.nzvalLocal(i) += 
        //						zweightRho_[l].real() * AinvMat.nzvalLocal(i).imag() + 
        //						zweightRho_[l].imag() * AinvMat.nzvalLocal(i).real();
        //				}
        // 
        // But done more cache-efficiently with blas.
        Real* AinvMatRealPtr = (Real*)AinvMat.nzvalLocal.Data();
        Real* AinvMatImagPtr = AinvMatRealPtr + 1;
        blas::Axpy( rhoMat.nnzLocal, zweightRho_[l].real(), AinvMatImagPtr, 2, 
            rhoMat.nzvalLocal.Data(), 1 );
        blas::Axpy( rhoMat.nnzLocal, zweightRho_[l].imag(), AinvMatRealPtr, 2,
            rhoMat.nzvalLocal.Data(), 1 );

        // Derivative of the Fermi-Dirac with respect to mu
        blas::Axpy( rhoDrvMuMat.nnzLocal, zweightRhoDrvMu_[l].real(), AinvMatImagPtr, 2, 
            rhoDrvMuMat.nzvalLocal.Data(), 1 );
        blas::Axpy( rhoDrvMuMat.nnzLocal, zweightRhoDrvMu_[l].imag(), AinvMatRealPtr, 2,
            rhoDrvMuMat.nzvalLocal.Data(), 1 );

        // Free energy density matrix
        if( isFreeEnergyDensityMatrix ){
          blas::Axpy( hmzMat.nnzLocal, zweightHelmholtz_[l].real(), AinvMatImagPtr, 2,
              hmzMat.nzvalLocal.Data(), 1 );
          blas::Axpy( hmzMat.nnzLocal, zweightHelmholtz_[l].imag(), AinvMatRealPtr, 2,
              hmzMat.nzvalLocal.Data(), 1 );
        }

        // Energy density matrix
        if( isEnergyDensityMatrix ){
          blas::Axpy( frcMat.nnzLocal, zweightForce_[l].real(), AinvMatImagPtr, 2,
              frcMat.nzvalLocal.Data(), 1 );
          blas::Axpy( frcMat.nnzLocal, zweightForce_[l].imag(), AinvMatRealPtr, 2, 
              frcMat.nzvalLocal.Data(), 1 );
        }

        // Derivative of the Fermi-Dirac with respect to T
        if( isDerivativeTMatrix ){
          blas::Axpy( rhoDrvTMat.nnzLocal, zweightRhoDrvT_[l].real(), AinvMatImagPtr, 2, 
              rhoDrvTMat.nzvalLocal.Data(), 1 );
          blas::Axpy( rhoDrvTMat.nnzLocal, zweightRhoDrvT_[l].imag(), AinvMatRealPtr, 2,
              rhoDrvTMat.nzvalLocal.Data(), 1 );
        }

        GetTime( timePostProcessingEnd );

        if( verbosity >= 1 ){
          statusOFS << "Time for postprocessing is " <<
            timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
        }

      }
      GetTime( timePoleEnd );

      if( verbosity >= 1 ){
        statusOFS << "Time for pole " << lidx << " is " <<
          timePoleEnd - timePoleSta << " [s]" << std::endl << std::endl;
      }

    } // if I am in charge of this pole


  } // for(lidx)

  // Reduce the density matrix across the processor rows in gridPole_
  {
    DblNumVec nzvalRhoMatLocal = rhoMat.nzvalLocal;
    SetValue( rhoMat.nzvalLocal, 0.0 );

    mpi::Allreduce( nzvalRhoMatLocal.Data(), rhoMat.nzvalLocal.Data(),
        rhoMat.nnzLocal, MPI_SUM, gridPole_->colComm );
  }

  // Reduce the derivative of density matrix with respect to mu across
  // the processor rows in gridPole_ 
  {
    DblNumVec nzvalRhoDrvMuMatLocal = rhoDrvMuMat.nzvalLocal;
    SetValue( rhoDrvMuMat.nzvalLocal, 0.0 );

    mpi::Allreduce( nzvalRhoDrvMuMatLocal.Data(), rhoDrvMuMat.nzvalLocal.Data(),
        rhoDrvMuMat.nnzLocal, MPI_SUM, gridPole_->colComm );
  }

  // Reduce the free energy density matrix across the processor rows in gridPole_ 
  if( isFreeEnergyDensityMatrix ){
    DblNumVec nzvalHmzMatLocal = hmzMat.nzvalLocal;
    SetValue( hmzMat.nzvalLocal, 0.0 );

    mpi::Allreduce( nzvalHmzMatLocal.Data(), hmzMat.nzvalLocal.Data(),
        hmzMat.nnzLocal, MPI_SUM, gridPole_->colComm );
  }

  // Reduce the energy density matrix across the processor rows in gridPole_ 
  if( isEnergyDensityMatrix ){
    DblNumVec nzvalFrcMatLocal = frcMat.nzvalLocal;
    SetValue( frcMat.nzvalLocal, 0.0 );

    mpi::Allreduce( nzvalFrcMatLocal.Data(), frcMat.nzvalLocal.Data(),
        frcMat.nnzLocal, MPI_SUM, gridPole_->colComm );
  }

  // Reduce the derivative of density matrix with respect to T across
  // the processor rows in gridPole_ 
  if( isDerivativeTMatrix ){
    DblNumVec nzvalRhoDrvTMatLocal = rhoDrvTMat.nzvalLocal;
    SetValue( rhoDrvTMat.nzvalLocal, 0.0 );

    mpi::Allreduce( nzvalRhoDrvTMatLocal.Data(), rhoDrvTMat.nzvalLocal.Data(),
        rhoDrvTMat.nnzLocal, MPI_SUM, gridPole_->colComm );
  }

  // Compute the number of electrons
  // The number of electrons is computed by Tr[DM*S]
  {
    Real numElecLocal = 0.0;
    if( SMat.size != 0 ){
      // S is not an identity matrix
      numElecLocal = blas::Dot( SMat.nnzLocal, SMat.nzvalLocal.Data(),
          1, rhoMat.nzvalLocal.Data(), 1 );
    }
    else{
      // S is an identity matrix
      DblNumVec& nzval = rhoMat.nzvalLocal;
      for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
        numElecLocal += nzval(diagIdxLocal_[i]);
      }
    } // if ( SMat.size != 0 )
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "numElecLocal = " << numElecLocal << std::endl;
#endif

    mpi::Allreduce( &numElecLocal, &numElectron, 1, MPI_SUM, rhoMat.comm ); 
  }

  // Compute the derivative of the number of electrons with respect to mu
  // The number of electrons is computed by Tr[f'(H-\muS)*S]
  {
    Real numElecDrvLocal = 0.0;
    if( SMat.size != 0 ){
      // S is not an identity matrix
      numElecDrvLocal = blas::Dot( SMat.nnzLocal, SMat.nzvalLocal.Data(),
          1, rhoDrvMuMat.nzvalLocal.Data(), 1 );
    }
    else{
      // S is an identity matrix
      DblNumVec& nzval = rhoDrvMuMat.nzvalLocal;
      for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
        numElecDrvLocal += nzval(diagIdxLocal_[i]);
      }
    }

#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "numElecDrvLocal = " << numElecDrvLocal << std::endl;
#endif

    mpi::Allreduce( &numElecDrvLocal, &numElectronDrvMu, 1, MPI_SUM, rhoDrvMuMat.comm ); 
  }

#ifndef _RELEASE_
  PopCallStack();
#endif

  return ;
}    // -----  end of method PPEXSIData::CalculateFermiOperatorReal  ----- 


void
PPEXSIData::DFTDriver (
    Real       numElectronExact,
    Real       temperature,
    Real       gap,
    Real       deltaE,
    Int        numPole, 
    Int        isInertiaCount,
    Int        maxPEXSIIter,
    Real       muMin0,
    Real       muMax0,
    Real       mu0,
    Real       muInertiaTolerance,
    Real       muInertiaExpansion,
    Real       muPEXSISafeGuard,
    Real       numElectronPEXSITolerance,
    Int        matrixType,
    Int        isSymbolicFactorize,
    Int        ordering,
    Int        numProcSymbFact,
    Int        verbosity,
    Real&      muPEXSI,
    Real&      numElectronPEXSI,         
    Real&      muMinInertia,              
    Real&      muMaxInertia,             
    Int&       numTotalInertiaIter,   
    Int&       numTotalPEXSIIter )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::DFTDriver");
#endif
  numTotalInertiaIter = 0;
  numTotalPEXSIIter   = 0;

	Real timeSta, timeEnd;
  Real timeInertiaSta, timeInertiaEnd;
  Real timePEXSISta, timePEXSIEnd;
  Real timeTotalSta, timeTotalEnd;
  Real timeInertia = 0.0;
  Real timePEXSI   = 0.0;

  Real muLower;
  Real muUpper;

  // Initial setup
  muMinInertia = muMin0;
  muMaxInertia = muMax0;
  muPEXSI = mu0;        

  Real numElectronDrvMuPEXSI;

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
#ifdef USE_ABORT
      abort();
#endif
throw std::logic_error("Unsupported ordering strategy.");
  }

  if( matrixType != 0 ){
#ifdef USE_ABORT
    abort();
#endif
throw std::logic_error("Unsupported matrixType. The variable has to be 0.");
  }

  // Perform symbolic factorization first if required
  if( isSymbolicFactorize == true ){
    if( verbosity >= 1 ){
      statusOFS << "Perform symbolic factorization now." << std::endl;
    }
    if( matrixType == 0 ){
      SymbolicFactorizeRealSymmetricMatrix( 
          colPerm, 
          numProcSymbFact,
          verbosity );


      SymbolicFactorizeComplexSymmetricMatrix( 
          colPerm, 
          numProcSymbFact,
          verbosity );
    }
  }
  else{
    if( verbosity >= 1 ){
      statusOFS << "Skip symbolic factorization" << std::endl
        << "NOTE: This assumes that symbolic factorization has been done, "
        << "and the input H and S matrices have the same sparisty as "
        << "previously used." << std::endl;
    }
  }

  // Main loop: inertia count + PEXSI

  // Maximum number of total iterations
  const Int maxTotalInertiaIter = 10;
  // Whether the given interval contains the chemical potential
  bool  isBadBound = false;  
  // Whether the Newton update is in the trusted range
  bool  isDeltaMuSafe = true;
  // Whether the whole process has converged
  bool  isConverged = false;
   
  GetTime( timeTotalSta );
  while( numTotalInertiaIter < maxTotalInertiaIter ){

    // Perform inertia count
    if( isInertiaCount == 1 ){
      GetTime( timeInertiaSta );

      if( verbosity >= 1 ){
        PrintBlock( statusOFS, "Inertia counting phase" );
      }


      // Number of shifts is exactly determined by the number of
      // independent groups to minimize the cost
      // However, the minimum number of shifts is 10 to accelerate
      // convergence.
      Int numShift = std::max( gridPole_->numProcRow, 10 );
      std::vector<Real>  shiftVec( numShift );
      std::vector<Real>  inertiaVec( numShift );   // Zero temperature
      std::vector<Real>  inertiaFTVec( numShift ); // Finite temperature
      Int maxInertiaIter = std::max( 1, (Int)std::ceil( 
            std::log( (muMax0 - muMin0) / muInertiaTolerance ) /
            std::log( static_cast<Real>(numShift) ) ) ); 

      if( verbosity >= 1 ){
        statusOFS << std::endl
          << "From" << std::endl
          << "(muMin0, muMax0)   = " << "(" << muMinInertia << ", " << muMaxInertia
          << ")" << std::endl
          << "muInertiaTolerance = " << muInertiaTolerance << std::endl
          << "numShift           = " << numShift << std::endl
          << "we have " << std::endl
          << "maxInertiaIter     = " << maxInertiaIter << std::endl << std::endl;
      }

      for( Int iter = 1; iter <= maxInertiaIter; iter++ ){
        if( numTotalInertiaIter >= maxTotalInertiaIter ){
          std::ostringstream msg;
          msg  << std::endl
            << maxTotalInertiaIter 
            << " inertia counts have been proceeded." << std::endl
            << "Try to revise the initial interval for the chemical potential, "
            << "or increase muInertiaTolerance. " << std::endl;
#ifdef USE_ABORT
          abort();
#endif
          throw std::runtime_error( msg.str().c_str() );
        }

        numTotalInertiaIter++;

        for( Int l = 0; l < numShift; l++ ){
          shiftVec[l] = muMinInertia + l * (muMaxInertia - muMinInertia) / (numShift-1);
        }


        GetTime( timeSta );
        if( matrixType == 0 ){
          CalculateNegativeInertiaReal(
              shiftVec,
              inertiaVec,
              verbosity );
        }

        GetTime( timeEnd );

        // Inertia is multiplied by 2.0 to reflect the doubly occupied
        // orbitals.
        for( Int l = 0; l < numShift; l++ ){
          inertiaVec[l] *= 2.0;
        }

        // Using linear interpolation procedure to compute the finite
        // temperature cumulative DOS
        {
          Int numInp = 1000;   // Number of interpolation points
          Real shiftExpand = 20 * temperature;   // Expand the interval for interpolation
          std::vector<Real>  shiftInpVec( numInp );
          std::vector<Real>  inertiaInpVec( numInp ); 
          std::vector<Real>  fdDrvVec( numInp );
          for( Int l = 0; l < numShift; l++ ){
            Real shiftInp0 = shiftVec[l] - shiftExpand;
            Real shiftInp1 = shiftVec[l] + shiftExpand; 
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

            if( verbosity >= 1 ){
              statusOFS << std::setiosflags(std::ios::left) 
                << std::setw(LENGTH_VAR_NAME) << "Shift      = "
                << std::setw(LENGTH_VAR_DATA) << shiftVec[l]
                << std::setw(LENGTH_VAR_NAME) << "Inertia    = "
                << std::setw(LENGTH_VAR_DATA) << inertiaVec[l]
                << std::setw(LENGTH_VAR_NAME) << "InertiaFT  = "
                << std::setw(LENGTH_VAR_DATA) << inertiaFTVec[l]
                << std::endl;
            }
          } // for (l)
        }

        if( verbosity >= 1 ){
          statusOFS << std::endl << "Time for iteration " 
            << iter << " of the inertia count is " 
            << timeEnd - timeSta << std::endl;
        }

        // If the chemical potential does not fall into the
        // (muMinInertia, muMaxInertia) bracket, expand the interval.

        const Real EPS = 1e-1;

        isBadBound = false;

        if( inertiaFTVec[0] > numElectronExact ||
            inertiaVec[0] > numElectronExact - EPS ){
          isBadBound = true;
          muMaxInertia = muMaxInertia;
          muMinInertia = muMinInertia - muInertiaExpansion;
        }

        if( inertiaFTVec[numShift-1] < numElectronExact ||
            inertiaVec[numShift-1] < numElectronExact + EPS ){
          isBadBound = true;
          muMinInertia = muMinInertia;
          muMaxInertia = muMaxInertia + muInertiaExpansion;
        }

        if( isBadBound == true ){
          if( verbosity >= 1 ){
            statusOFS << std::endl << std::endl
              << "The solution is not in the provided interval." << std::endl
              << "(muMin, muMax) = ( " << shiftVec[0] << " , " << shiftVec[numShift-1] << " ) " << std::endl
              << "(Ne(muMin), Ne(muMax)) = ( " << inertiaFTVec[0] << " , " << inertiaFTVec[numShift-1] 
              << " ) " << std::endl
              << "NeExact = " << numElectronExact << std::endl
              << "Refine the interval to " << std::endl
              << "(muMin, muMax) = ( " << muMinInertia << " , " << muMaxInertia << " ) " << std::endl;
          }
          // Exit the loop
          break;
        }

        // Update muMin, muMax

        // First find the smallest interval
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

        if( verbosity >= 1 ){
          statusOFS << "idxMin = " << idxMin << ", inertiaVec = " << inertiaVec[idxMin] << std::endl;
          statusOFS << "idxMax = " << idxMax << ", inertiaVec = " << inertiaVec[idxMax] << std::endl;
        }


        // Heuristics to increase the span of mu a bit
//        muMinInertia = shiftVec[idxMin] - 5 * temperature;
//        muMaxInertia = shiftVec[idxMax] + 5 * temperature;
        muMinInertia = shiftVec[idxMin];
        muMaxInertia = shiftVec[idxMax];
        // Search instead for the band edges which is more stable
        muLower      = MonotoneRootFinding( shiftVec, inertiaFTVec, numElectronExact - 0.01);
        muUpper      = MonotoneRootFinding( shiftVec, inertiaFTVec, numElectronExact + 0.01 );
        muPEXSI =      (muLower + muUpper)/2.0;

        if( verbosity >= 1 ){
          statusOFS << "muLower = " << muLower << std::endl;
          statusOFS << "muUpper = " << muUpper << std::endl;
          statusOFS << "mu guessed by IC = " << muPEXSI << std::endl;
          statusOFS << "|ivec(idxMax)-Ne_exact| = " << 
                      std::abs( inertiaVec[idxMax] - numElectronExact ) << std::endl;
          statusOFS << "ivec(idxMax)-ivec(idxMin) = " << 
                       shiftVec[idxMax] - shiftVec[idxMin] << std::endl;
        }

        // Check convergence. Stop the inertia count after convergence.
        if( ( std::abs( inertiaVec[idxMax] - numElectronExact ) < EPS ) ||
            ( ( shiftVec[idxMax] - shiftVec[idxMin] ) < muInertiaTolerance ) ){
          isInertiaCount = 0;
        }
      } // for (iter)

      GetTime( timeInertiaEnd );
      timeInertia += timeInertiaEnd - timeInertiaSta;
    }

    // Immediately continue the inertia counting procedure
    if( isBadBound == true && isInertiaCount == 1 ){
      continue;
    }


    if( verbosity >= 1 ){
      PrintBlock( statusOFS, "PEXSI phase" );
    }

    // PEXSI iteration
    Real timeMuSta, timeMuEnd;
    isDeltaMuSafe = true;



    GetTime( timePEXSISta );
    for(Int iter = 1; iter <= maxPEXSIIter; iter++){
      GetTime( timeMuSta );

      numTotalPEXSIIter++;

      if( verbosity >= 1 ){
        statusOFS << "PEXSI Iteration " << iter 
          << ", mu = " << muPEXSI << std::endl;
      }

      if( matrixType == 0 ){
        CalculateFermiOperatorReal(
            numPole,
            temperature,
            gap,
            deltaE,
            muPEXSI,
            numElectronExact,
            numElectronPEXSITolerance,
            verbosity,
            numElectronPEXSI,
            numElectronDrvMuPEXSI );
      }

      GetTime( timeMuEnd );

      if( verbosity >= 1 ){
        statusOFS << std::endl;
        Print( statusOFS, "mu                          = ", muPEXSI );
        Print( statusOFS, "Computed number of electron = ", numElectronPEXSI );
        Print( statusOFS, "Exact number of electron    = ", numElectronExact );
        Print( statusOFS, "d Ne / d mu                 = ", numElectronDrvMuPEXSI );
        statusOFS << std::endl << "Time for mu iteration " << iter << " is " <<
          timeMuEnd - timeMuSta << " [s]" << std::endl << std::endl;
      }

      if( std::abs( numElectronPEXSI - numElectronExact ) < 
          numElectronPEXSITolerance ){
        isConverged = true;
        if( verbosity >= 0 ){
          statusOFS << "PEXSI has converged within " 
            << iter << " steps." << std::endl
            << "Exiting the PEXSI iteration now." << std::endl;
        }
        break;
      }
      else{

        // Update the chemical potential using Newton's method
        Real deltaMu = 
          -(numElectronPEXSI - numElectronExact) / numElectronDrvMuPEXSI;

        if( std::abs( deltaMu ) > muPEXSISafeGuard ){
          isDeltaMuSafe = false;
          isInertiaCount = 1;
          statusOFS << std::endl 
            << "deltaMu computed from Newton's method is " << deltaMu << std::endl
            << "and is larger than the trusted range " << muPEXSISafeGuard << std::endl
            << "Exiting PEXSI iteration now and reactivate the inertia counting scheme." 
            << std::endl;
          break;
        }
        else{
          muPEXSI += deltaMu;
        }
      }
    } // for (iter)
    GetTime( timePEXSIEnd );

    timePEXSI = timePEXSIEnd - timePEXSISta;

    // Exit no matter whether PEXSI has converged.
    if( isDeltaMuSafe == true ){
      break;
    }
  }

  GetTime( timeTotalEnd );
  
  // Compute the energy, and free energy
  if( matrixType == 0 ){
    // Energy computed from Tr[H*DM]
    {
      Real local = 0.0;
      local = blas::Dot( HRealMat_.nnzLocal, 
          HRealMat_.nzvalLocal.Data(),
          1, rhoRealMat_.nzvalLocal.Data(), 1 );
      mpi::Allreduce( &local, &totalEnergyH_, 1, MPI_SUM, 
          gridPole_->rowComm ); 
    }
  
    // Energy computed from Tr[S*EDM]
    {
      Real local = 0.0;
      if( SRealMat_.size != 0 ){
        local = blas::Dot( SRealMat_.nnzLocal, 
            SRealMat_.nzvalLocal.Data(),
            1, energyDensityRealMat_.nzvalLocal.Data(), 1 );
      }
      else{
        DblNumVec& nzval = energyDensityRealMat_.nzvalLocal;
        for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
          local += nzval(diagIdxLocal_[i]);
        }
      }

      mpi::Allreduce( &local, &totalEnergyS_, 1, MPI_SUM, 
          gridPole_->rowComm ); 
    }


    // Free energy 
    {
      Real local = 0.0;
      if( SRealMat_.size != 0 ){
        local = blas::Dot( SRealMat_.nnzLocal, 
            SRealMat_.nzvalLocal.Data(),
            1, freeEnergyDensityRealMat_.nzvalLocal.Data(), 1 );
      }
      else{
        DblNumVec& nzval = freeEnergyDensityRealMat_.nzvalLocal;
        for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
          local += nzval(diagIdxLocal_[i]);
        }
      }

      mpi::Allreduce( &local, &totalFreeEnergy_, 1, MPI_SUM, 
          gridPole_->rowComm ); 

      // Correction
      totalFreeEnergy_ += muPEXSI * numElectronPEXSI;
    }
  } // if( matrixType == 0 )

  if( verbosity == 1 ){
    if( isConverged == 1 ){
      statusOFS << std::endl << "PEXSI has converged!" << std::endl;
    }
    else{
      statusOFS << std::endl << "PEXSI did not converge!" << std::endl;
    }
    statusOFS << std::endl
      << "Total number of inertia counts       = " << numTotalInertiaIter << std::endl
      << "Total number of PEXSI iterations     = " << numTotalPEXSIIter << std::endl 
      << "Total time for inertia counts        = " << timeInertia << " [s]" << std::endl 
      << "Total time for PEXSI iterations      = " << timePEXSI   << " [s]" << std::endl
      << "Total time for the DFT driver        = " << timeTotalEnd - timeTotalSta   << " [s]" << std::endl
      << std::endl;
  }

  if( verbosity >= 1 ){
    statusOFS << "Final result " << std::endl;
    Print( statusOFS, "mu                          = ", muPEXSI );
    Print( statusOFS, "Computed number of electron = ", numElectronPEXSI );
    Print( statusOFS, "Exact number of electron    = ", numElectronExact );
    Print( statusOFS, "d Ne / d mu                 = ", numElectronDrvMuPEXSI );
		Print( statusOFS, "Total energy (H*DM)         = ", totalEnergyH_ );
		Print( statusOFS, "Total energy (S*EDM)        = ", totalEnergyS_ );
		Print( statusOFS, "Total free energy           = ", totalFreeEnergy_ );
    statusOFS << std::endl << std::endl;
  }


#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::DFTDriver  ----- 

void
PPEXSIData::DFTDriver2 (
    Real       numElectronExact,
    Real       temperature,
    Real       gap,
    Real       deltaE,
    Int        numPole, 
    Int        isInertiaCount,
    Real       muMin0,
    Real       muMax0,
    Real       mu0,
    Real       muInertiaTolerance,
    Real       muInertiaExpansion,
    Real       numElectronPEXSITolerance,
    Int        matrixType,
    Int        isSymbolicFactorize,
    Int        ordering,
    Int        numProcSymbFact,
    Int        verbosity,
    Real&      muPEXSI,
    Real&      numElectronPEXSI,         
    Real&      muMinInertia,              
    Real&      muMaxInertia,             
    Int&       numTotalInertiaIter ) 
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::DFTDriver2");
#endif
	Real timeSta, timeEnd;
  Real timeInertiaSta, timeInertiaEnd;
  Real timePEXSISta, timePEXSIEnd;
  Real timeTotalSta, timeTotalEnd;
  Real timeInertia = 0.0;
  Real timePEXSI   = 0.0;

  Real muLower;
  Real muUpper;

  // Initial setup
  muMinInertia = muMin0;
  muMaxInertia = muMax0;
  muPEXSI = mu0;        

  GetTime( timeTotalSta );

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
#ifdef USE_ABORT
      abort();
#endif
      throw std::logic_error("Unsupported ordering strategy.");
  }

  if( matrixType != 0 ){
#ifdef USE_ABORT
    abort();
#endif
    throw std::logic_error("Unsupported matrixType. The variable has to be 0.");
  }

  if( muInertiaTolerance < 4.0 * temperature ){
    statusOFS << "muInertiaTolerance cannot be smaller than 4*temperature = "  << 
      4.0 * temperature << std::endl;
    statusOFS << "Forcefully set muInertiaTolerance to 4*temperature." << std::endl;

    muInertiaTolerance = 4.0 * temperature;
  }


  // Perform symbolic factorization first if required
  if( isSymbolicFactorize == true ){
    if( verbosity >= 1 ){
      statusOFS << "Perform symbolic factorization now." << std::endl;
    }
    if( matrixType == 0 ){
      SymbolicFactorizeRealSymmetricMatrix( 
          colPerm, 
          numProcSymbFact,
          verbosity );


      SymbolicFactorizeComplexSymmetricMatrix( 
          colPerm, 
          numProcSymbFact,
          verbosity );
    }
  }
  else{
    if( verbosity >= 1 ){
      statusOFS << "Skip symbolic factorization" << std::endl
        << "NOTE: This assumes that symbolic factorization has been done, "
        << "and the input H and S matrices have the same sparisty as "
        << "previously used." << std::endl;
    }
  }
  

  // Inertia counting loop
  numTotalInertiaIter = 0;
  // Hard coded
  const Int maxTotalInertiaIter = 10;

  if( isInertiaCount == 1 ){
    if( verbosity >= 1 ){
      PrintBlock( statusOFS, "Inertia counting phase" );
    }
    bool  isBadBound = false;  

    while( numTotalInertiaIter < maxTotalInertiaIter ){
      GetTime( timeInertiaSta );

      // Number of shifts is exactly determined by the number of
      // independent groups to minimize the cost
      // However, the minimum number of shifts is 10 to accelerate
      // convergence.
      Int numShift = std::max( gridPole_->numProcRow, 10 );
      std::vector<Real>  shiftVec( numShift );
      std::vector<Real>  inertiaVec( numShift );   // Zero temperature
      std::vector<Real>  inertiaFTVec( numShift ); // Finite temperature
      Int maxInertiaIter = std::max( 1, (Int)std::ceil( 
            std::log( (muMaxInertia - muMinInertia) / muInertiaTolerance ) /
            std::log( static_cast<Real>(numShift) ) ) ); 

      if( verbosity >= 1 ){
        statusOFS << std::endl
          << "From" << std::endl
          << "(muMin0, muMax0)   = " << "(" << muMinInertia << ", " << muMaxInertia
          << ")" << std::endl
          << "muInertiaTolerance = " << muInertiaTolerance << std::endl
          << "numShift           = " << numShift << std::endl
          << "we have " << std::endl
          << "maxInertiaIter     = " << maxInertiaIter << std::endl << std::endl;
      }

      for( Int iter = 1; iter <= maxInertiaIter; iter++ ){
        if( numTotalInertiaIter >= maxTotalInertiaIter ){
          std::ostringstream msg;
          msg  << std::endl
            << maxTotalInertiaIter 
            << " inertia counts have been proceeded." << std::endl
            << "Try to revise the initial interval for the chemical potential, "
            << "or increase muInertiaTolerance. " << std::endl;
#ifdef USE_ABORT
          abort();
#endif
          throw std::runtime_error( msg.str().c_str() );
        }

        numTotalInertiaIter++;

        for( Int l = 0; l < numShift; l++ ){
          shiftVec[l] = muMinInertia + l * (muMaxInertia - muMinInertia) / (numShift-1);
        }


        GetTime( timeSta );
        if( matrixType == 0 ){
          CalculateNegativeInertiaReal(
              shiftVec,
              inertiaVec,
              verbosity );
        }

        GetTime( timeEnd );

        // Inertia is multiplied by 2.0 to reflect the doubly occupied
        // orbitals.  
        //
        // FIXME In the future numSpin should be introduced.
        for( Int l = 0; l < numShift; l++ ){
          inertiaVec[l] *= 2.0;
        }

        // Using linear interpolation procedure to compute the finite
        // temperature cumulative DOS
        {
          Int numInp = 1000;   // Number of interpolation points
          // Expand the interval for interpolation
          Real shiftExpand = 20 * temperature;   
          std::vector<Real>  shiftInpVec( numInp );
          std::vector<Real>  inertiaInpVec( numInp ); 
          std::vector<Real>  fdDrvVec( numInp );
          for( Int l = 0; l < numShift; l++ ){
            Real shiftInp0 = shiftVec[l] - shiftExpand;
            Real shiftInp1 = shiftVec[l] + shiftExpand; 
            Real h = (shiftInp1 - shiftInp0) / (numInp-1);
            for( Int i = 0; i < numInp; i++ ){
              shiftInpVec[i] = shiftInp0 + h * i;
              // fdDrvMu(x) = beta * exp(beta*x)/(1+exp(beta*z))^2
              // Note: compared to fdDrvMu used in pole.cpp, the
              // factor 2.0 is not present, because it is included in
              // inertiaVec. 
              fdDrvVec[i]    = 1.0 / ( 2.0 * temperature * 
                  ( 1.0 + std::cosh( ( shiftInpVec[i] - shiftVec[l] )
                                     / temperature ) ) ); }
              LinearInterpolation( shiftVec, inertiaVec, 
                  shiftInpVec, inertiaInpVec );

              inertiaFTVec[l] = 0.0;
              for( Int i = 0; i < numInp; i++ ){
                inertiaFTVec[l] += fdDrvVec[i] * inertiaInpVec[i] * h;
              }

              if( verbosity >= 1 ){
                statusOFS << std::setiosflags(std::ios::left) 
                  << std::setw(LENGTH_VAR_NAME) << "Shift      = "
                  << std::setw(LENGTH_VAR_DATA) << shiftVec[l]
                  << std::setw(LENGTH_VAR_NAME) << "Inertia    = "
                  << std::setw(LENGTH_VAR_DATA) << inertiaVec[l]
                  << std::setw(LENGTH_VAR_NAME) << "InertiaFT  = "
                  << std::setw(LENGTH_VAR_DATA) << inertiaFTVec[l]
                  << std::endl;
              }
          } // for (l)
        }

        if( verbosity >= 1 ){
          statusOFS << std::endl << "Time for iteration " 
            << iter << " of the inertia count is " 
            << timeEnd - timeSta << std::endl;
        }

        // If the chemical potential does not fall into the
        // (muMinInertia, muMaxInertia) bracket, expand the interval.

        const Real EPS = 1e-1;

        isBadBound = false;

        if( inertiaFTVec[0] > numElectronExact ||
            inertiaVec[0] > numElectronExact - EPS ){
          isBadBound = true;
          muMaxInertia = muMaxInertia;
          muMinInertia = muMinInertia - muInertiaExpansion;
        }

        if( inertiaFTVec[numShift-1] < numElectronExact ||
            inertiaVec[numShift-1] < numElectronExact + EPS ){
          isBadBound = true;
          muMinInertia = muMinInertia;
          muMaxInertia = muMaxInertia + muInertiaExpansion;
        }

        if( isBadBound == true ){
          if( verbosity >= 1 ){
            statusOFS << std::endl << std::endl
              << "The solution is not in the provided interval." << std::endl
              << "(muMin, muMax) = ( " << shiftVec[0] << " , " << shiftVec[numShift-1] << " ) " << std::endl
              << "(Ne(muMin), Ne(muMax)) = ( " << inertiaFTVec[0] << " , " << inertiaFTVec[numShift-1] 
              << " ) " << std::endl
              << "NeExact = " << numElectronExact << std::endl
              << "Refine the interval to " << std::endl
              << "(muMin, muMax) = ( " << muMinInertia << " , " << muMaxInertia << " ) " << std::endl;
          }
          // Exit the loop
          break;
        }

        // Update muMin, muMax

        // First find the smallest interval
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

        if( verbosity >= 1 ){
          statusOFS << "idxMin = " << idxMin << ", inertiaVec = " << inertiaVec[idxMin] << std::endl;
          statusOFS << "idxMax = " << idxMax << ", inertiaVec = " << inertiaVec[idxMax] << std::endl;
        }


        // Heuristics to increase the span of mu a bit
        // LL: 2/17/2016 NOTE This might end up with a wide interval that inertia counting never converges
//        muMinInertia = shiftVec[idxMin] - 3 * temperature;
//        muMaxInertia = shiftVec[idxMax] + 3 * temperature;
        // LL: 2/17/2016 This might end up with a interval that is too tight especialy at high temperature
//        muMinInertia = shiftVec[idxMin];
//        muMaxInertia = shiftVec[idxMax];
        // Search instead for the band edges which is more stable
        muLower      = MonotoneRootFinding( shiftVec, inertiaFTVec, numElectronExact - 0.01);
        muUpper      = MonotoneRootFinding( shiftVec, inertiaFTVec, numElectronExact + 0.01 );
        muPEXSI      = (muLower + muUpper)/2.0;

        // LL: 2/23/2016: New way for updating the interval to avoid too skewed interval
        muMinInertia = std::min(shiftVec[idxMin], muPEXSI - 2.0 * temperature);
        muMaxInertia = std::max(shiftVec[idxMax], muPEXSI + 2.0 * temperature);

        if( verbosity >= 1 ){
          statusOFS << "muLower = " << muLower << std::endl;
          statusOFS << "muUpper = " << muUpper << std::endl;
          statusOFS << "mu guessed by IC = " << muPEXSI << std::endl;
//          statusOFS << "|ivec(idxMax)-Ne_exact| = " << 
//            std::abs( inertiaVec[idxMax] - numElectronExact ) << std::endl;
//          statusOFS << "ivec(idxMax)-ivec(idxMin) = " << 
//            shiftVec[idxMax] - shiftVec[idxMin] << std::endl;
          statusOFS << "muMinInertia = " << muMinInertia << std::endl;
          statusOFS << "muMaxInertia = " << muMaxInertia << std::endl;
        }

        // Check convergence. Stop the inertia count after convergence.
        if( ( std::abs( inertiaVec[idxMax] - numElectronExact ) < EPS ) ||
            ( ( shiftVec[idxMax] - shiftVec[idxMin] ) < muInertiaTolerance ) ){
          break;
        }
      } // for (iter)

      GetTime( timeInertiaEnd );
      timeInertia += timeInertiaEnd - timeInertiaSta;
      // Immediately continue the inertia counting procedure

      if( isBadBound == false ){
        if( verbosity >= 1 ){
          statusOFS << std::endl << "Inertia counting has converged!" << std::endl;
        }
        break;
      }
    } // while (inertiaCount)
  } // if (isInertiaCount)

  // PEXSI phase.
  // No option for falling back to inertia counting
  bool  isConverged = false;
  {
    if( verbosity >= 1 ){
      PrintBlock( statusOFS, "PEXSI phase" );
    }

    GetTime( timePEXSISta );

    if( matrixType == 0 ){
      // The new routine not only computes the number of electrons at
      // mu, but also updates the value of mu using a chemical potential
      // update procedure by only perfoming the selected inversion once.
      //
      // Always perform one PEXSI iteration only and rely on the update
      // strategy to compute the chemical potential.
      //
      // 11/24/2015
      
      // muMinPEXSI / muMaxPEXSI are used to bound the searching
      // interval for the PEXSI update
      // 
      // The search range should be at most on the order of k_B T
      
      Real muWidth = std::min( muInertiaTolerance, 2.0*temperature );
      Real muMinPEXSI = muPEXSI - muWidth;
      Real muMaxPEXSI = muPEXSI + muWidth;

      CalculateFermiOperatorReal2(
          numPole,
          temperature,
          gap,
          deltaE,
          numElectronExact,
          numElectronPEXSITolerance,
          muMinPEXSI,
          muMaxPEXSI,
          verbosity,
          muPEXSI, 
          numElectronPEXSI, 
          isConverged );
    }
    
    GetTime( timePEXSIEnd );
    timePEXSI = timePEXSIEnd - timePEXSISta;
  }

  
  // Compute the energy, and free energy
  if( matrixType == 0 ){
    // Energy computed from Tr[H*DM]
    {
      Real local = 0.0;
      local = blas::Dot( HRealMat_.nnzLocal, 
          HRealMat_.nzvalLocal.Data(),
          1, rhoRealMat_.nzvalLocal.Data(), 1 );
      mpi::Allreduce( &local, &totalEnergyH_, 1, MPI_SUM, 
          gridPole_->rowComm ); 
    }
  
    // Energy computed from Tr[S*EDM]
    {
      Real local = 0.0;
      if( SRealMat_.size != 0 ){
        local = blas::Dot( SRealMat_.nnzLocal, 
            SRealMat_.nzvalLocal.Data(),
            1, energyDensityRealMat_.nzvalLocal.Data(), 1 );
      }
      else{
        DblNumVec& nzval = energyDensityRealMat_.nzvalLocal;
        for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
          local += nzval(diagIdxLocal_[i]);
        }
      }

      mpi::Allreduce( &local, &totalEnergyS_, 1, MPI_SUM, 
          gridPole_->rowComm ); 
    }


    // Free energy 
    {
      Real local = 0.0;
      if( SRealMat_.size != 0 ){
        local = blas::Dot( SRealMat_.nnzLocal, 
            SRealMat_.nzvalLocal.Data(),
            1, freeEnergyDensityRealMat_.nzvalLocal.Data(), 1 );
      }
      else{
        DblNumVec& nzval = freeEnergyDensityRealMat_.nzvalLocal;
        for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
          local += nzval(diagIdxLocal_[i]);
        }
      }

      mpi::Allreduce( &local, &totalFreeEnergy_, 1, MPI_SUM, 
          gridPole_->rowComm ); 

      // Correction
      totalFreeEnergy_ += muPEXSI * numElectronPEXSI;
    }
  } // if( matrixType == 0 )

  GetTime( timeTotalEnd );


  if( verbosity == 1 ){
    statusOFS << std::endl
      << "Total number of inertia counts       = " << numTotalInertiaIter << std::endl
      << "Total time for inertia count step    = " << timeInertia << " [s]" << std::endl 
      << "Total time for PEXSI step            = " << timePEXSI   << " [s]" << std::endl
      << "Total time for the DFT driver        = " << timeTotalEnd - timeTotalSta   << " [s]" << std::endl
      << std::endl;
  }

  if( verbosity >= 1 ){
    statusOFS << "Final result " << std::endl;
    Print( statusOFS, "mu                          = ", muPEXSI );
    Print( statusOFS, "Computed number of electron = ", numElectronPEXSI );
    Print( statusOFS, "Exact number of electron    = ", numElectronExact );
		Print( statusOFS, "Total energy (H*DM)         = ", totalEnergyH_ );
		Print( statusOFS, "Total energy (S*EDM)        = ", totalEnergyS_ );
		Print( statusOFS, "Total free energy           = ", totalFreeEnergy_ );
    statusOFS << std::endl << std::endl;
  }

  if( isConverged == false ){
    std::ostringstream msg;
    msg  << "PEXSI did not converge. " << std::endl;
    msg << "Aborting..." << std::endl;
#ifdef USE_ABORT
    abort();
#endif
    throw std::runtime_error( msg.str().c_str() );
  }

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::DFTDriver2  ----- 


// Main subroutine for the electronic structure calculation
// Use the same Green's functions at multiple mu so that PEXSI is only
// calculated once per SCF.
void PPEXSIData::CalculateFermiOperatorReal2(
    Int   numPole, 
    Real  temperature,
    Real  gap,
    Real  deltaE,
    Real  numElectronExact,
    Real  numElectronTolerance,
    Real  muMinPEXSI, 
    Real  muMaxPEXSI,
    Int   verbosity,
    Real& mu,
    Real& numElectron,
    bool& isPEXSIConverged){
#ifndef _RELEASE_
  PushCallStack("PPEXSIData::CalculateFermiOperatorReal2");
#endif
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealSymmetricMatrix first." << std::endl;
    #ifdef USE_ABORT
abort();
#endif
throw std::runtime_error( msg.str().c_str() );
  }

  if( isComplexSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeComplexSymmetricMatrix first." << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }

  // *********************************************************************
  // Check the input parameters
  // *********************************************************************
  if( numPole % 2 != 0 ){
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "Must be even number of poles!" );
  }

  // *********************************************************************
  // Initialize
  // *********************************************************************
  // Rename for convenience
  DistSparseMatrix<Real>&  HMat        = HRealMat_;
  DistSparseMatrix<Real>&  SMat        = SRealMat_;

  DistSparseMatrix<Real>& rhoMat       = rhoRealMat_;     
  DistSparseMatrix<Real>& rhoDrvMuMat  = rhoDrvMuRealMat_;
  DistSparseMatrix<Real>& rhoDrvTMat   = rhoDrvTRealMat_;
  DistSparseMatrix<Real>& hmzMat       = freeEnergyDensityRealMat_;
  DistSparseMatrix<Real>& frcMat       = energyDensityRealMat_;

  DistSparseMatrix<Complex>& AMat      = shiftComplexMat_;

  // The symbolic information should have already been there.
  SuperLUMatrix<Complex>& luMat        = *luComplexMat_;
  PMatrix<Complex>&       PMloc        = *PMComplexMat_;

  // 
  bool isFreeEnergyDensityMatrix = true;
  bool isEnergyDensityMatrix     = true;
  bool isDerivativeTMatrix       = false;


  // Copy the pattern
  CopyPattern( HMat, AMat );
  CopyPattern( HMat, rhoMat );
  CopyPattern( HMat, rhoDrvMuMat );
  if( isFreeEnergyDensityMatrix )
    CopyPattern( HMat, hmzMat );
  if( isEnergyDensityMatrix )
    CopyPattern( HMat, frcMat );
  if( isDerivativeTMatrix )
    CopyPattern( HMat, rhoDrvTMat );

  // Reinitialize the variables
  SetValue( rhoMat.nzvalLocal, 0.0 );
  SetValue( rhoDrvMuMat.nzvalLocal, 0.0 );
  if( isFreeEnergyDensityMatrix )
    SetValue( hmzMat.nzvalLocal, 0.0 );
  if( isEnergyDensityMatrix )
    SetValue( frcMat.nzvalLocal, 0.0 );
  if( isDerivativeTMatrix )
    SetValue( rhoDrvTMat.nzvalLocal, 0.0 );

  // Refine the pole expansion  
  // numPoleInput is the number of poles to be given to other parts of
  // the pole expansion, which is larger than or equal to numPole.
  Int numPoleInput;
  // poleIdx is a vector of size numPole.  Only poles with index in
  // poleIdx are used for actual computation. The rest of the poles are
  // discarded according to tolerance criterion
  //
  //   numElectronTolerance / numElectronExact / numPole
  //
  // FIXME The heuristics should be refined to give error estimate to
  // other quantities such as the energy.
  // FIXME The heuristics part should also be given in a separate
  // routine, and the input of this file does not need mu, gap etc.
  std::vector<Int>  poleIdx(numPole);
  {
    // Setup a grid from (mu - deltaE, mu + deltaE), and measure
    // the error bound in the L^infty sense on this grid.
    //
    // fdGrid:      Exact Fermi-Dirac function evaluated on xGrid
    // fdPoleGrid:  Fermi-Dirac function using pole expansion
    // evaluated on the grid.
    Int numX = 10000;
    std::vector<Real>    xGrid( numX );
    std::vector<Real>    fdGrid( numX );

    Real x0 = mu - deltaE;
    Real x1 = mu + deltaE;
    Real h  = (x1 - x0) / (numX - 1);
    Real ez;
    for( Int i = 0; i < numX; i++ ){
      xGrid[i]  = x0 + i * h;
      if( xGrid[i] - mu >= 0 ){
        ez = std::exp(- (xGrid[i] - mu) / temperature );
        fdGrid[i] = 2.0 * ez / (1.0 + ez);
      }
      else{
        ez = std::exp((xGrid[i] - mu) / temperature );
        fdGrid[i] = 2.0 / (1.0 + ez);
      }
    }


    numPoleInput = numPole;
    Real tol;
    Int  numPoleSignificant;

    Int poleIter = 0;
    do{
      // If the number of significant poles is less than numPole,
      // increase numPoleInput by 2 at a time and redo the
      // computation.
      if( poleIter > 0 )
        numPoleInput += 2;

      zshift_.resize( numPoleInput );
      zweightRho_.resize( numPoleInput );
      GetPoleDensity( &zshift_[0], &zweightRho_[0],
          numPoleInput, temperature, gap, deltaE, mu ); 

      std::vector<Complex>  zshiftTmp( numPoleInput );
      zweightForce_.resize( numPoleInput );
      GetPoleForce( &zshiftTmp[0], &zweightForce_[0],
          numPoleInput, temperature, gap, deltaE, mu ); 


      std::vector<Real>  maxMagPole(numPoleInput);
      for( Int l = 0; l < numPoleInput; l++ ){
        maxMagPole[l] = 0.0;
      }

      // Compute the approximation due to pole expansion, as well as
      // the maximum magnitude of each pole
      Complex cpxmag;
      Real    mag;
      numPoleSignificant = 0;
      tol = numElectronTolerance / numElectronExact / numPoleInput;
      for( Int l = 0; l < numPoleInput; l++ ){
        for( Int i = 0; i < numX; i++ ){
          cpxmag = zweightRho_[l] / ( xGrid[i] - zshift_[l] );
          mag    = cpxmag.imag();
          maxMagPole[l] = ( maxMagPole[l] >= mag ) ?  maxMagPole[l] : mag;
        }
        if( maxMagPole[l] > tol ){
          numPoleSignificant++;
        }	
      } // for (l)

      // Pick the most significant numPole poles and update poleIdx
      // Sort in DESCENDING order
      std::vector<Int>  sortIdx( numPoleInput );
      for( Int i = 0; i < sortIdx.size(); i++ ){
        sortIdx[i]      = i;
      }
      std::sort( sortIdx.begin(), sortIdx.end(), 
          IndexComp<std::vector<Real>& >( maxMagPole ) ) ;
      std::reverse( sortIdx.begin(), sortIdx.end() );

      for( Int l = 0; l < numPole; l++ ){
        poleIdx[l]      = sortIdx[l];
      }


      // Update poleIter
      poleIter++; 
    } while( numPoleSignificant < numPole );


    // Estimate the error of the number of electrons and the band energy
    // by assuming a flat density of states within a interval of size
    // deltaE, i.e. each unit interval contains 
    // HMat.size / deltaE 
    // number of electrons

    std::vector<Real>  fdPoleGrid( numX );
    std::vector<Real>  fdEPoleGrid( numX );
    std::vector<Real>  fdTimesEPoleGrid( numX );
    Real errAbs1, errorAbsMax1;
    Real errAbs2, errorAbsMax2;
    Real errorNumElectron, errorBandEnergy, errorBandEnergy2;
    Complex cpxmag1, cpxmag2;

    errorAbsMax1 = 0.0; 
    errorAbsMax2 = 0.0; 
    for( Int i = 0; i < numX; i++ ){
      fdPoleGrid[i] = 0.0;
      fdEPoleGrid[i] = 0.0;
      for( Int lidx = 0; lidx < numPoleInput; lidx++ ){
        Int l = lidx;
        cpxmag1 = zweightRho_[l] / ( xGrid[i] - zshift_[l] );
        cpxmag2 = zweightForce_[l] / ( xGrid[i] - zshift_[l] );
        fdPoleGrid[i] += cpxmag1.imag();
        fdTimesEPoleGrid[i] += cpxmag1.imag() * xGrid[i];
        fdEPoleGrid[i]      += cpxmag2.imag();
      }
      errAbs1 = std::abs( fdPoleGrid[i] - fdGrid[i] );
      errorAbsMax1 = ( errorAbsMax1 >= errAbs1 ) ? errorAbsMax1 : errAbs1;
      errAbs2 = std::abs( fdEPoleGrid[i] - fdTimesEPoleGrid[i] );
      errorAbsMax2 = ( errorAbsMax2 >= errAbs2 ) ? errorAbsMax2 : errAbs2;
    }

    errorNumElectron = errorAbsMax1 * HMat.size;
    errorBandEnergy  = 0.5 * deltaE * errorAbsMax1 * HMat.size;
    errorBandEnergy2 = errorAbsMax2 * HMat.size;

    if( verbosity >= 1 ){
      statusOFS 
        << std::endl 
        << "Estimated error of pole expansion by assuming a flat spectrum."
        << std::endl;

      // The estimation of energy using the difference of DM and EDM
      // seems to be more reliable
      Print( statusOFS, "Error of num electron             = ", errorNumElectron );
//      Print( statusOFS, "Error of band energy (DM only)    = ", errorBandEnergy );
      Print( statusOFS, "Error of band energy (DM and EDM) = ", errorBandEnergy2 );
      Print( statusOFS, "Required accuracy (num electron)  = ", numElectronTolerance );

      statusOFS << std::endl;

      if( errorNumElectron > numElectronTolerance ){
        statusOFS << "WARNING!!! " 
          << "Pole expansion may not be accurate enough to reach numElectronTolerance. " << std::endl
          << "Try to increase numPole or increase numElectronTolerance." << std::endl << std::endl;
      }
      statusOFS << "numPoleInput =" << numPoleInput << std::endl;
      statusOFS << "numPoleSignificant = " << numPoleSignificant << std::endl;
    }
  }

  // Initialize the number of electrons
  numElectron  = 0.0;

  //Initialize the pole expansion
  zweightRhoDrvMu_.resize( numPoleInput );

  GetPoleDensityDrvMu( &zshift_[0], &zweightRhoDrvMu_[0],
      numPoleInput, temperature, gap, deltaE, mu ); 

  if( verbosity >= 2 ){
    statusOFS << "zshift" << std::endl << zshift_ << std::endl;
    statusOFS << "zweightRho" << std::endl << zweightRho_ << std::endl;
  }

  // *********************************************************************
  // For each pole, perform LDLT factoriation and selected inversion
  // Store each Green's function for subsequent post-processing.
  //
  // This is the most time consuming part.
  // *********************************************************************

  Real timePoleSta, timePoleEnd;

  Int numPoleComputed = 0;

  // Store the Tr[Ainv*S] to compute the trace at multiple mu.
  CpxNumVec traceAinvSLocal(numPole);
  CpxNumVec traceAinvS(numPole);
  SetValue(traceAinvSLocal, Z_ZERO);
  SetValue(traceAinvS, Z_ZERO);
  // Need to store all Ainv matrices. If empty then Ainv corresponding
  // to this pole is not used.
  std::vector<DistSparseMatrix<Complex>> AinvMatVec(numPole);

  for(Int lidx = 0; lidx < numPole; lidx++){
    if( MYROW( gridPole_ ) == PROW( lidx, gridPole_ ) ){

      Int l = poleIdx[lidx];

      GetTime( timePoleSta );

      if( verbosity >= 1 ){
        statusOFS << "Pole " << lidx << " processing..." << std::endl;
      }
      if( verbosity >= 2 ){
        statusOFS << "zshift           = " << zshift_[l] << std::endl;
        statusOFS	<< "zweightRho       = " << zweightRho_[l] << std::endl;
        statusOFS	<< "zweightRhoDrvMu  = " << zweightRhoDrvMu_[l] << std::endl;
        if( isFreeEnergyDensityMatrix )
          statusOFS << "zweightHelmholtz = " << zweightHelmholtz_[l] << std::endl;
        if( isEnergyDensityMatrix )
          statusOFS << "zweightForce     = " << zweightForce_[l] << std::endl;
        if( isDerivativeTMatrix )
          statusOFS << "zweightRhoDrvT   = " << zweightRhoDrvT_[l] << std::endl;
      }

      {
        DistSparseMatrix<Complex>& AinvMat = AinvMatVec[lidx];

        AinvMat.Clear();

        numPoleComputed++;

        if( SMat.size != 0 ){
          // S is not an identity matrix
          for( Int i = 0; i < HMat.nnzLocal; i++ ){
            AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift_[l] * SMat.nzvalLocal(i);
          }
        }
        else{
          // S is an identity matrix
          for( Int i = 0; i < HMat.nnzLocal; i++ ){
            AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
          }

          for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
            AMat.nzvalLocal( diagIdxLocal_[i] ) -= zshift_[l];
          }
        } // if (SMat.size != 0 )


        // *********************************************************************
        // Factorization
        // *********************************************************************
        // Important: the distribution in pzsymbfact is going to mess up the
        // A matrix.  Recompute the matrix A here.
        if( verbosity >= 2 ){
          statusOFS << "Before DistSparseMatrixToSuperMatrixNRloc." << std::endl;
        }
        luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
        if( verbosity >= 2 ){
          statusOFS << "After DistSparseMatrixToSuperMatrixNRloc." << std::endl;
        }

        Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

        GetTime( timeTotalFactorizationSta );

        // Data redistribution
        if( verbosity >= 2 ){
          statusOFS << "Before Distribute." << std::endl;
        }
        luMat.Distribute();
        if( verbosity >= 2 ){
          statusOFS << "After Distribute." << std::endl;
        }

        // Numerical factorization
        if( verbosity >= 2 ){
          statusOFS << "Before NumericalFactorize." << std::endl;
        }
        luMat.NumericalFactorize();
        if( verbosity >= 2 ){
          statusOFS << "After NumericalFactorize." << std::endl;
        }
        luMat.DestroyAOnly();

        GetTime( timeTotalFactorizationEnd );

        if( verbosity >= 1 ){
          statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 
        }

        // *********************************************************************
        // Selected inversion
        // *********************************************************************
        Real timeTotalSelInvSta, timeTotalSelInvEnd;
        GetTime( timeTotalSelInvSta );

        luMat.LUstructToPMatrix( PMloc );


        // Collective communication version
        //          PMloc.ConstructCommunicationPattern_Collectives();

        PMloc.PreSelInv();

        // Main subroutine for selected inversion
        //
        // P2p communication version
        PMloc.SelInv();

        // Collective communication version
        //          PMloc.SelInv_Collectives();

        GetTime( timeTotalSelInvEnd );

        if( verbosity >= 1 ){
          statusOFS << "Time for total selected inversion is " <<
            timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
        }

        // *********************************************************************
        // Postprocessing
        // *********************************************************************

        Real timeConvertingSta, timeConvertingEnd;

        GetTime( timeConvertingSta );

        PMloc.PMatrixToDistSparseMatrix( AMat, AinvMat );

        // Compute Tr[Ainv*S]. 
        //
        // The number of electrons can be computed as
        //
        // Ne(mu) = sum_l imag( Tr[Ainv_{l}*S] * omega_l )
        {
          Complex trace = Z_ZERO;
          if( SMat.size != 0 ){
            // S is not an identity matrix
            DblNumVec& nzvalS = SMat.nzvalLocal;
            CpxNumVec& nzvalAinv = AinvMat.nzvalLocal;

            for( Int i = 0; i < SMat.nnzLocal; i++ ){
              trace += nzvalAinv(i) * nzvalS(i);
            }
          }
          else{
            // S is an identity matrix
            CpxNumVec& nzvalAinv = AinvMat.nzvalLocal;
            for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
              trace += nzvalAinv(diagIdxLocal_[i]);
            }
          } // if ( SMat.size != 0 )

          traceAinvSLocal[lidx] = trace;
        }

        GetTime( timeConvertingEnd );

        if( verbosity >= 1 ){
          statusOFS << "Time for converting the matrix from PMatrix format "
            << " to DistSparseMatrix fomat is " <<
            timeConvertingEnd - timeConvertingSta << " [s]" << std::endl;
        }

      }
      GetTime( timePoleEnd );

      if( verbosity >= 1 ){
        statusOFS << "Time for pole " << lidx << " is " <<
          timePoleEnd - timePoleSta << " [s]" << std::endl << std::endl;
      }
    } // if I am in charge of this pole
  } // for(lidx)

  // *********************************************************************
  // Reuse the Green's functions and update the chemical potential
  // 
  // Important to reuse poleIdx to get new weights of the pole expansion.
  // *********************************************************************
  
  // Reduce the information of trace of Ainv*S over ALL processors.
  mpi::Allreduce( traceAinvSLocal.Data(), traceAinvS.Data(), 
      numPole, MPI_SUM, gridPole_->comm ); 


  // A bisection strategy for finding mu.
  if( mu < muMinPEXSI || mu > muMaxPEXSI ){
    std::ostringstream msg;
    msg  << std::endl
      << "The chemical potential mu = " << mu  << std::endl
      << "This is outside the interval [muMinPEXSI,muMaxPEXSI] = " 
      << "[" << muMinPEXSI << ", " << muMaxPEXSI << "]" << std::endl;
    msg << "Aborting..." << std::endl;
//#ifdef USE_ABORT
//    abort();
//#endif
    throw std::runtime_error( msg.str().c_str() );
  }

  Real muInit = mu;
  Real numElectronInit;
  Real dmu = 0.0;
  Int maxBisectionIter = 50;
  isPEXSIConverged = false;

  std::vector<Complex> zshiftDummy( numPoleInput );
  // Overwrite zweightRho_ by the updated weights.
  for( Int il = 1; il < maxBisectionIter; il++ ){
    dmu = mu - muInit;
    GetPoleDensityUpdate( &zshiftDummy[0], &zweightRho_[0],
        numPoleInput, temperature, gap, deltaE, muInit, dmu ); 
    // Compute the number of electrons 
    numElectron = 0.0;
    for(Int lidx = 0; lidx < numPole; lidx++){
      Int l = poleIdx[lidx];
      numElectron += zweightRho_[l].real() * traceAinvS[lidx].imag() + 
        zweightRho_[l].imag() * traceAinvS[lidx].real();
    } // for(lidx)

    if( il == 1 )
      numElectronInit = numElectron;

    if( std::abs( numElectron - numElectronExact ) < 
        numElectronTolerance ){
      isPEXSIConverged = true;
      break;
    }
    
    // Bisection
    if( numElectron < numElectronExact ){
      muMinPEXSI = mu;
      mu = 0.5*(muMaxPEXSI+mu);
    }
    else{
      muMaxPEXSI = mu;
      mu = 0.5*(muMinPEXSI+mu);
    }
  } // for (bisection)

  if( isPEXSIConverged == true ){
    statusOFS << std::endl << "PEXSI has converged!" << std::endl;
  }
  else{
    statusOFS << std::endl 
      << "PEXSI did not converge." << std::endl
      << "The best guess is mu = " << mu << ", numElectron = " 
      << numElectron << std::endl
      << "You can consider the following options." << std::endl
      << "1. Make the pole expansion more accurate" << std::endl
      << "2. Make the inertia counting more accurate" << std::endl;
  }

  if( verbosity >= 1 ){
    statusOFS 
      << "mu (input)           = " << muInit << std::endl
      << "numElectron (input)  = " << numElectronInit << std::endl
      << "mu (output)          = " << mu << std::endl
      << "numElectron (output) = " << numElectron << std::endl;
  }

#if ( _DEBUGlevel_ >= 0 )
#endif


  // Compute the pole expansion for density etc using the update formula

  if( isFreeEnergyDensityMatrix ){
    zweightHelmholtz_.resize( numPoleInput );
    GetPoleHelmholtzUpdate( &zshiftDummy[0], &zweightHelmholtz_[0], 
        numPoleInput, temperature, gap, deltaE, muInit, dmu ); 
  }

  if( isEnergyDensityMatrix ){
    zweightForce_.resize( numPoleInput );
    GetPoleForceUpdate( &zshiftDummy[0], &zweightForce_[0],
        numPoleInput, temperature, gap, deltaE, muInit, dmu ); 
  }

  // *********************************************************************
  // Post-processing for obtaining the density matrix, and other derived
  // quantities at the corrected mu.
  // *********************************************************************

  // Update the density matrix. The following lines are equivalent to
  //
  //				for( Int i = 0; i < rhoMat.nnzLocal; i++ ){
  //					rhoMat.nzvalLocal(i) += 
  //						zweightRho_[l].real() * AinvMat.nzvalLocal(i).imag() + 
  //						zweightRho_[l].imag() * AinvMat.nzvalLocal(i).real();
  //				}
  // 
  // But done more cache-efficiently with blas.
  // 
  // Similarly for other quantities.
  
  for(Int lidx = 0; lidx < numPole; lidx++){
    if( MYROW( gridPole_ ) == PROW( lidx, gridPole_ ) ){
      Int l = poleIdx[lidx];
      Real* AinvMatRealPtr = (Real*)(AinvMatVec[lidx].nzvalLocal.Data());
      Real* AinvMatImagPtr = AinvMatRealPtr + 1;
      blas::Axpy( rhoMat.nnzLocal, zweightRho_[l].real(), AinvMatImagPtr, 2, 
          rhoMat.nzvalLocal.Data(), 1 );
      blas::Axpy( rhoMat.nnzLocal, zweightRho_[l].imag(), AinvMatRealPtr, 2,
          rhoMat.nzvalLocal.Data(), 1 );

      // Free energy density matrix
      if( isFreeEnergyDensityMatrix ){
        blas::Axpy( hmzMat.nnzLocal, zweightHelmholtz_[l].real(), AinvMatImagPtr, 2,
            hmzMat.nzvalLocal.Data(), 1 );
        blas::Axpy( hmzMat.nnzLocal, zweightHelmholtz_[l].imag(), AinvMatRealPtr, 2,
            hmzMat.nzvalLocal.Data(), 1 );
      }


      // Energy density matrix
      if( isEnergyDensityMatrix ){
        blas::Axpy( frcMat.nnzLocal, zweightForce_[l].real(), AinvMatImagPtr, 2,
            frcMat.nzvalLocal.Data(), 1 );
        blas::Axpy( frcMat.nnzLocal, zweightForce_[l].imag(), AinvMatRealPtr, 2, 
            frcMat.nzvalLocal.Data(), 1 );
      }


    } // if I am in charge of this pole
  } // for(lidx)

  // Reduce the density matrix across the processor rows in gridPole_
  {
    DblNumVec nzvalRhoMatLocal = rhoMat.nzvalLocal;
    SetValue( rhoMat.nzvalLocal, 0.0 );

    mpi::Allreduce( nzvalRhoMatLocal.Data(), rhoMat.nzvalLocal.Data(),
        rhoMat.nnzLocal, MPI_SUM, gridPole_->colComm );
  }

  // Reduce the free energy density matrix across the processor rows in gridPole_ 
  if( isFreeEnergyDensityMatrix ){
    DblNumVec nzvalHmzMatLocal = hmzMat.nzvalLocal;
    SetValue( hmzMat.nzvalLocal, 0.0 );

    mpi::Allreduce( nzvalHmzMatLocal.Data(), hmzMat.nzvalLocal.Data(),
        hmzMat.nnzLocal, MPI_SUM, gridPole_->colComm );
  }

  // Reduce the energy density matrix across the processor rows in gridPole_ 
  if( isEnergyDensityMatrix ){
    DblNumVec nzvalFrcMatLocal = frcMat.nzvalLocal;
    SetValue( frcMat.nzvalLocal, 0.0 );

    mpi::Allreduce( nzvalFrcMatLocal.Data(), frcMat.nzvalLocal.Data(),
        frcMat.nnzLocal, MPI_SUM, gridPole_->colComm );
  }


  // Free the space for the saved Green's functions
  AinvMatVec.clear();

#ifndef _RELEASE_
  PopCallStack();
#endif

  return ;
}    // -----  end of method PPEXSIData::CalculateFermiOperatorReal2  ----- 


} //  namespace PEXSI
