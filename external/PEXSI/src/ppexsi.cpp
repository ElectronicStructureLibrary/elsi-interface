/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.

Author: Lin Lin and Weile Jia

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
/// @date Revision:      2016-09-04  Update interface for unsymmetric
/// solvers.


#include "ppexsi.hpp"
#include "pexsi/utility.hpp"
#include "pexsi/getPole.hpp"

namespace PEXSI{

PPEXSIData::PPEXSIData    (
    MPI_Comm   comm,
    Int        numProcRow,
    Int        numProcCol,
    Int        outputFileIndex ){

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
    ErrorHandling( msg.str().c_str() );
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
  isRealUnsymmetricSymbolicFactorized_ = false;
  isComplexUnsymmetricSymbolicFactorized_ = false;

  // Initialize the empty matrices
  luRealMat_ = new SuperLUMatrix<Real>;
  luComplexMat_ = new SuperLUMatrix<Complex>;

#ifdef WITH_SYMPACK
  outputFileIndex_ = outputFileIndex;
  symPACKRealMat_ = nullptr;
  symPACKComplexMat_ = nullptr;
  symPACKOpt_.verbose=0;
  //    if( outputFileIndex >= 0 ){
  //      if(symPACK::logfileptr==NULL){
  //        //Initialize symPACK logfile
  //        std::stringstream suffix;
  //        suffix<<mpirank;
  //        symPACK::logfileptr = new symPACK::LogFile("status",suffix.str().c_str());
  //        symPACK::logfileptr->OFS()<<"********* LOGFILE OF P"<<mpirank<<" *********"<<std::endl;
  //        symPACK::logfileptr->OFS()<<"**********************************"<<std::endl;
  //      }
  //    }
  //
  //    symPACKRealMat_ = new symPACK::symPACKMatrix<Real>;
  //    symPACKComplexMat_ = new symPACK::symPACKMatrix<Complex>;
#endif

  PMRealMat_ = new PMatrix<Real>;
  PMComplexMat_ = new PMatrix<Complex>;

  PMRealUnsymMat_ = new PMatrixUnsym<Real>;
  PMComplexUnsymMat_ = new PMatrixUnsym<Complex>;

  return ;
}         // -----  end of method PPEXSIData::PPEXSIData  -----


PPEXSIData::~PPEXSIData    (  )
{
  if( luRealMat_ != NULL ){
    delete luRealMat_;
  }

  if( luComplexMat_ != NULL ){
    delete luComplexMat_;
  }


#ifdef WITH_SYMPACK
  if ( symPACKRealMat_ != NULL ){
    delete symPACKRealMat_;
  }

  if ( symPACKComplexMat_ != NULL ){
    delete symPACKComplexMat_;
  }

  delete symPACK::logfileptr;
#endif


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


  return ;
}         // -----  end of method PPEXSIData::~PPEXSIData  -----

void
PPEXSIData::LoadRealMatrix    (
    Int           nrows,
    Int           nnz,
    Int           nnzLocal,
    Int           numColLocal,
    Int*          colptrLocal,
    Int*          rowindLocal,
    Real*         HnzvalLocal,
    Int           isSIdentity,
    Real*         SnzvalLocal,
    Int           solver,
    Int           verbosity )
{
  // Clear the previously saved information
  HRealMat_ = DistSparseMatrix<Real>();
  SRealMat_ = DistSparseMatrix<Real>();
  DistSparseMatrix<Real>&  SMat    = SRealMat_;
  //#ifdef WITH_SYMPACK
  //      symmHRealMat_ = symPACK::DistSparseMatrix<Real>();
  //      symmSRealMat_ = symPACK::DistSparseMatrix<Real>();
  //#endif

  // Data communication
  //      switch (solver){
  //        case 0:
  //          {
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
      SMat.size = 0;
      SMat.nnz  = 0;
      SMat.nnzLocal = 0;
      SMat.comm = HRealMat_.comm;
    }
    else{
      CopyPattern( HRealMat_, SMat );
      SMat.comm = HRealMat_.comm;
      SMat.nzvalLocal  = DblNumVec( nnzLocal,      false, SnzvalLocal );
      serialize( SMat.nzvalLocal, sstm, NO_MASK );
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
      SMat.size = 0;    // Means S is an identity matrix
      SMat.nnz  = 0;
      SMat.nnzLocal = 0;
      SMat.comm = HRealMat_.comm;
    }
    else{
      CopyPattern( HRealMat_, SMat );
      SMat.comm = HRealMat_.comm;
      deserialize( SMat.nzvalLocal, sstm, NO_MASK );
    }
  }
  sstr.clear();


  if( verbosity >= 1 ){
    statusOFS << "H.size     = " << HRealMat_.size     << std::endl;
    statusOFS << "H.nnzLocal = " << HRealMat_.nnzLocal << std::endl;
    statusOFS << "S.size     = " << SMat.size     << std::endl;
    statusOFS << "S.nnzLocal = " << SMat.nnzLocal << std::endl;
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

  CopyPattern( HRealMat_, PatternMat_ );

#ifdef WITH_SYMPACK
  if(solver == 1){
    Convert(PatternMat_, symmPatternMat_);
    symmPatternMat_.ToLowerTriangular();
    //symmPatternMat_.GetLocalGraph().SetSorted(false);
    //symmPatternMat_.SortGraph();
  }
#endif
  //          }
  //          break;
  //#ifdef WITH_SYMPACK
  //        case 1:
  //          {
  //            std::vector<char> sstr;
  //            Int sizeStm;
  //            if( MYROW( gridPole_ ) == 0 ){
  //              std::stringstream sstm;
  //
  //              symmHRealMat_.comm        = gridPole_->rowComm;
  //              symmHRealMat_.size        = nrows;
  //              symmHRealMat_.nnz         = nnz;
  //              // The first row processor does not need extra copies of the index /
  //              // value of the matrix.
  //              symPACK::DistSparseMatrixGraph & Hgraph = symmHRealMat_.GetLocalGraph();
  //              Hgraph.size    = nrows;
  //              Hgraph.nnz    = nnz;
  //              Hgraph.SetComm(symmHRealMat_.comm);
  //              Hgraph.colptr.resize(numColLocal+1);
  //              std::copy(colptrLocal,colptrLocal+numColLocal+1,&Hgraph.colptr[0]);
  //              Hgraph.rowind.resize(nnzLocal+1);
  //              std::copy(rowindLocal,rowindLocal+nnzLocal,&Hgraph.rowind[0]);
  //              Hgraph.SetExpanded(true);
  //
  //              //          int mpisize,mpirank;
  //              //          MPI_Comm_size(symmHRealMat_.comm,&mpisize);
  //              //          MPI_Comm_rank(symmHRealMat_.comm,&mpirank);
  //              //
  //              //symPACK::gdb_lock();
  //              //          Int colPerProc = Hgraph.size / mpisize;
  //              //          Hgraph.vertexDist.resize(mpisize+1,colPerProc);
  //              //          Hgraph.vertexDist[0] = 1;
  //              //          std::partial_sum(Hgraph.vertexDist.begin(),Hgraph.vertexDist.end(),Hgraph.vertexDist.begin());
  //              //          Hgraph.vertexDist.back() = Hgraph.size+1;
  //              //symPACK::gdb_lock();
  //
  //              // H value
  //              symmHRealMat_.nzvalLocal.resize( nnzLocal);
  //              std::copy(HnzvalLocal,HnzvalLocal+nnzLocal,&symmHRealMat_.nzvalLocal[0]);
  //              //To Lower Triagular
  //              symmHRealMat_.ToLowerTriangular();
  //
  //              // Serialization will copy the values regardless of the ownership
  //              //TODO implement this
  //              PEXSI::serialize( symmHRealMat_, sstm, NO_MASK );
  //
  //              // S value
  //              if( isSIdentity ){
  //                symmSRealMat_.size = 0;
  //                symmSRealMat_.nnz  = 0;
  //                symmSRealMat_.GetLocalGraph().nnz = 0;
  //                symmSRealMat_.comm = symmHRealMat_.comm;
  //              }
  //              else{
  //                //CopyPattern( symmPatternMat_, symmSRealMat_ );
  //                PEXSI::CopyPattern( symmHRealMat_, symmSRealMat_ );
  //                symmSRealMat_.comm = symmHRealMat_.comm;
  //                symmSRealMat_.nzvalLocal.resize( nnzLocal);
  //                std::copy(SnzvalLocal,SnzvalLocal+nnzLocal,&symmSRealMat_.nzvalLocal[0]);
  //                //To Lower Triagular
  //                symmSRealMat_.ToLowerTriangular();
  //
  //                serialize( symmSRealMat_.nzvalLocal, sstm, NO_MASK );
  //              }
  //
  //              sstr.resize( Size( sstm ) );
  //              sstm.read( &sstr[0], sstr.size() );
  //              sizeStm = sstr.size();
  //            }
  //
  //            MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole_->colComm );
  //
  //            if( verbosity >= 2 ){
  //              statusOFS << "sizeStm = " << sizeStm << std::endl;
  //            }
  //
  //            if( MYROW( gridPole_ ) != 0 ) sstr.resize( sizeStm );
  //
  //            MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole_->colComm );
  //
  //            if( MYROW( gridPole_ ) != 0 ){
  //              std::stringstream sstm;
  //              sstm.write( &sstr[0], sizeStm );
  //              PEXSI::deserialize( symmHRealMat_, sstm, NO_MASK );
  //              // Communicator
  //              symmHRealMat_.comm = gridPole_->rowComm;
  //              if( isSIdentity ){
  //                symmSRealMat_.size = 0;
  //                symmSRealMat_.nnz  = 0;
  //                symmSRealMat_.GetLocalGraph().nnz = 0;
  //                symmSRealMat_.comm = symmHRealMat_.comm;
  //              }
  //              else{
  //                PEXSI::CopyPattern( symmHRealMat_, symmSRealMat_ );
  //                symmSRealMat_.comm = symmHRealMat_.comm;
  //                deserialize( symmSRealMat_.nzvalLocal, sstm, NO_MASK );
  //              }
  //            }
  //            sstr.clear();
  //
  //
  //            if( verbosity >= 1 ){
  //              statusOFS << "H.size     = " << symmHRealMat_.size     << std::endl;
  //              statusOFS << "H.nnzLocal = " << symmHRealMat_.GetLocalGraph().LocalEdgeCount() << std::endl;
  //              statusOFS << "S.size     = " << symmSRealMat_.size     << std::endl;
  //              statusOFS << "S.nnzLocal = " << symmSRealMat_.GetLocalGraph().LocalEdgeCount() << std::endl;
  //              statusOFS << std::endl << std::endl;
  //            }
  //
  //
  //            // Record the index for the diagonal elements to handle the case if S
  //            // is identity.
  //            {
  //              const symPACK::DistSparseMatrixGraph & Hgraph = symmHRealMat_.GetLocalGraph();
  //              Int numColLocal      = Hgraph.LocalVertexCount();
  //              Int firstCol         = Hgraph.LocalFirstVertex();
  //
  //              diagIdxLocal_.clear();
  //              diagIdxLocal_.reserve(symmHRealMat_.size);
  //              for( Int j = 0; j < numColLocal; j++ ){
  //                Int jcol = firstCol + j + 1;
  //                for( Int i = Hgraph.colptr[j]-1;
  //                    i < Hgraph.colptr[j+1]-1; i++ ){
  //                  Int irow = Hgraph.rowind[i];
  //                  if( irow == jcol ){
  //                    diagIdxLocal_.push_back( i );
  //                  }
  //                }
  //              } // for (j)
  //            }
  //
  //            isMatrixLoaded_ = true;
  //
  //            PEXSI::CopyPattern( symmHRealMat_, symmPatternMat_ );
  //            Convert(symmPatternMat_,PatternMat_);
  //
  //            Convert(symmSRealMat_,SRealMat_);
  //            Convert(symmHRealMat_,HRealMat_);
  //
  //
  //          }
  //          break;
  //#endif
  //        default:
  //          ErrorHandling("Unsupported solver.");
  //          break;
  //      }
  return ;
}        // -----  end of method PPEXSIData::LoadRealMatrix  -----

void
PPEXSIData::LoadComplexMatrix    (
    Int           nrows,
    Int           nnz,
    Int           nnzLocal,
    Int           numColLocal,
    Int*          colptrLocal,
    Int*          rowindLocal,
    Complex*      HnzvalLocal,
    Int           isSIdentity,
    Complex*      SnzvalLocal,
    Int           solver,
    Int           verbosity )
{
  // Clear the previously saved information
  HComplexMat_ = DistSparseMatrix<Complex>();
  SComplexMat_ = DistSparseMatrix<Complex>();
  DistSparseMatrix<Complex>&  SMat  = SComplexMat_;
  //#ifdef WITH_SYMPACK
  //      symmHComplexMat_ = symPACK::DistSparseMatrix<Complex>();
  //      symmSComplexMat_ = symPACK::DistSparseMatrix<Complex>();
  //#endif
  //
  //      // Data communication
  //      switch (solver) {
  //        case 0:
  //          {
  std::vector<char> sstr;
  Int sizeStm;
  if( MYROW( gridPole_ ) == 0 ){
    std::stringstream sstm;

    HComplexMat_.size        = nrows;
    HComplexMat_.nnz         = nnz;
    HComplexMat_.nnzLocal    = nnzLocal;
    // The first row processor does not need extra copies of the index /
    // value of the matrix.
    HComplexMat_.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
    HComplexMat_.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
    // H value
    HComplexMat_.nzvalLocal  = CpxNumVec( nnzLocal,      false, HnzvalLocal );
    HComplexMat_.comm        = gridPole_->rowComm;

    // Serialization will copy the values regardless of the ownership
    serialize( HComplexMat_, sstm, NO_MASK );

    // S value
    if( isSIdentity ){
      SMat.size = 0;
      SMat.nnz  = 0;
      SMat.nnzLocal = 0;
      SMat.comm = HComplexMat_.comm;
    }
    else{
      PEXSI::CopyPattern( HComplexMat_, SMat );
      SMat.comm = HComplexMat_.comm;
      SMat.nzvalLocal  = CpxNumVec( nnzLocal,      false, SnzvalLocal );
      serialize( SMat.nzvalLocal, sstm, NO_MASK );
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
    deserialize( HComplexMat_, sstm, NO_MASK );
    // Communicator
    HComplexMat_.comm = gridPole_->rowComm;
    if( isSIdentity ){
      SMat.size = 0;    // Means S is an identity matrix
      SMat.nnz  = 0;
      SMat.nnzLocal = 0;
      SMat.comm = HComplexMat_.comm;
    }
    else{
      CopyPattern( HComplexMat_, SMat );
      SMat.comm = HComplexMat_.comm;
      deserialize( SMat.nzvalLocal, sstm, NO_MASK );
    }
  }
  sstr.clear();


  if( verbosity >= 1 ){
    statusOFS << "H.size     = " << HComplexMat_.size     << std::endl;
    statusOFS << "H.nnzLocal = " << HComplexMat_.nnzLocal << std::endl;
    statusOFS << "S.size     = " << SMat.size     << std::endl;
    statusOFS << "S.nnzLocal = " << SMat.nnzLocal << std::endl;
    statusOFS << std::endl << std::endl;
  }


  // Record the index for the diagonal elements to handle the case if S
  // is identity.
  {
    Int numColLocal      = HComplexMat_.colptrLocal.m() - 1;
    Int numColLocalFirst = HComplexMat_.size / gridSelInv_->mpisize;
    Int firstCol         = gridSelInv_->mpirank * numColLocalFirst;

    diagIdxLocal_.clear();

    for( Int j = 0; j < numColLocal; j++ ){
      Int jcol = firstCol + j + 1;
      for( Int i = HComplexMat_.colptrLocal(j)-1;
          i < HComplexMat_.colptrLocal(j+1)-1; i++ ){
        Int irow = HComplexMat_.rowindLocal(i);
        if( irow == jcol ){
          diagIdxLocal_.push_back( i );
        }
      }
    } // for (j)
  }

  isMatrixLoaded_ = true;

  CopyPattern( HComplexMat_, PatternMat_ );

#ifdef WITH_SYMPACK
  if(solver == 1){
    Convert(PatternMat_, symmPatternMat_);
    symmPatternMat_.ToLowerTriangular();
    //symmPatternMat_.GetLocalGraph().SetSorted(false);
    //symmPatternMat_.SortGraph();
  }
#endif

  //          }
  //          break;
  //#ifdef WITH_SYMPACK
  //        case 1:
  //          {
  //            std::vector<char> sstr;
  //            Int sizeStm;
  //            if( MYROW( gridPole_ ) == 0 ){
  //              std::stringstream sstm;
  //
  //              symmHComplexMat_.comm        = gridPole_->rowComm;
  //              symmHComplexMat_.size        = nrows;
  //              symmHComplexMat_.nnz         = nnz;
  //              // The first row processor does not need extra copies of the index /
  //              // value of the matrix.
  //              symPACK::DistSparseMatrixGraph & Hgraph = symmHComplexMat_.GetLocalGraph();
  //              Hgraph.size    = nrows;
  //              Hgraph.nnz    = nnz;
  //              Hgraph.SetComm(symmHComplexMat_.comm);
  //              Hgraph.colptr.resize(numColLocal+1);
  //              std::copy(colptrLocal,colptrLocal+numColLocal+1,&Hgraph.colptr[0]);
  //              Hgraph.rowind.resize(nnzLocal+1);
  //              std::copy(rowindLocal,rowindLocal+nnzLocal,&Hgraph.rowind[0]);
  //              Hgraph.SetExpanded(true);
  //
  //              //          int mpisize,mpirank;
  //              //          MPI_Comm_size(symmHComplexMat_.comm,&mpisize);
  //              //          MPI_Comm_rank(symmHComplexMat_.comm,&mpirank);
  //              //
  //              //          Int colPerProc = Hgraph.size / mpisize;
  //              //          Hgraph.vertexDist.resize(mpisize+1,colPerProc);
  //              //          Hgraph.vertexDist[0] = 1;
  //              //          std::partial_sum(Hgraph.vertexDist.begin(),Hgraph.vertexDist.end(),Hgraph.vertexDist.begin());
  //              //          Hgraph.vertexDist.back() = Hgraph.size+1;
  //              //
  //              //
  //
  //              // H value
  //              symmHComplexMat_.nzvalLocal.resize( nnzLocal);
  //              std::copy(HnzvalLocal,HnzvalLocal+nnzLocal,&symmHComplexMat_.nzvalLocal[0]);
  //              //To Lower Triagular
  //              symmHComplexMat_.ToLowerTriangular();
  //
  //              // Serialization will copy the values regardless of the ownership
  //              //TODO implement this
  //              PEXSI::serialize( symmHComplexMat_, sstm, NO_MASK );
  //
  //              // S value
  //              if( isSIdentity ){
  //                symmSComplexMat_.size = 0;
  //                symmSComplexMat_.nnz  = 0;
  //                symmSComplexMat_.GetLocalGraph().nnz = 0;
  //                symmSComplexMat_.comm = symmHComplexMat_.comm;
  //              }
  //              else{
  //                //CopyPattern( symmPatternMat_, symmSComplexMat_ );
  //                PEXSI::CopyPattern( symmHComplexMat_, symmSComplexMat_ );
  //                symmSComplexMat_.comm = symmHComplexMat_.comm;
  //                symmSComplexMat_.nzvalLocal.resize( nnzLocal);
  //                std::copy(SnzvalLocal,SnzvalLocal+nnzLocal,&symmSComplexMat_.nzvalLocal[0]);
  //                //To Lower Triagular
  //                symmSComplexMat_.ToLowerTriangular();
  //
  //                serialize( symmSComplexMat_.nzvalLocal, sstm, NO_MASK );
  //              }
  //
  //              sstr.resize( Size( sstm ) );
  //              sstm.read( &sstr[0], sstr.size() );
  //              sizeStm = sstr.size();
  //            }
  //
  //            MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole_->colComm );
  //
  //            if( verbosity >= 2 ){
  //              statusOFS << "sizeStm = " << sizeStm << std::endl;
  //            }
  //
  //            if( MYROW( gridPole_ ) != 0 ) sstr.resize( sizeStm );
  //
  //            MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole_->colComm );
  //
  //            if( MYROW( gridPole_ ) != 0 ){
  //              std::stringstream sstm;
  //              sstm.write( &sstr[0], sizeStm );
  //              PEXSI::deserialize( symmHComplexMat_, sstm, NO_MASK );
  //              // Communicator
  //              symmHComplexMat_.comm = gridPole_->rowComm;
  //              if( isSIdentity ){
  //                symmSComplexMat_.size = 0;
  //                symmSComplexMat_.nnz  = 0;
  //                symmSComplexMat_.GetLocalGraph().nnz = 0;
  //                symmSComplexMat_.comm = symmHComplexMat_.comm;
  //              }
  //              else{
  //                PEXSI::CopyPattern( symmHComplexMat_, symmSComplexMat_ );
  //                symmSComplexMat_.comm = symmHComplexMat_.comm;
  //                deserialize( symmSComplexMat_.nzvalLocal, sstm, NO_MASK );
  //              }
  //            }
  //            sstr.clear();
  //
  //
  //            if( verbosity >= 1 ){
  //              statusOFS << "H.size     = " << symmHComplexMat_.size     << std::endl;
  //              statusOFS << "H.nnzLocal = " << symmHComplexMat_.GetLocalGraph().LocalEdgeCount() << std::endl;
  //              statusOFS << "S.size     = " << symmSComplexMat_.size     << std::endl;
  //              statusOFS << "S.nnzLocal = " << symmSComplexMat_.GetLocalGraph().LocalEdgeCount() << std::endl;
  //              statusOFS << std::endl << std::endl;
  //            }
  //
  //
  //            // Record the index for the diagonal elements to handle the case if S
  //            // is identity.
  //            {
  //              const symPACK::DistSparseMatrixGraph & Hgraph = symmHComplexMat_.GetLocalGraph();
  //              Int numColLocal      = Hgraph.LocalVertexCount();
  //              Int firstCol         = Hgraph.LocalFirstVertex();
  //
  //              diagIdxLocal_.clear();
  //              diagIdxLocal_.reserve(symmHComplexMat_.size);
  //              for( Int j = 0; j < numColLocal; j++ ){
  //                Int jcol = firstCol + j + 1;
  //                for( Int i = Hgraph.colptr[j]-1;
  //                    i < Hgraph.colptr[j+1]-1; i++ ){
  //                  Int irow = Hgraph.rowind[i];
  //                  if( irow == jcol ){
  //                    diagIdxLocal_.push_back( i );
  //                  }
  //                }
  //              } // for (j)
  //            }
  //
  //            isMatrixLoaded_ = true;
  //
  //            PEXSI::CopyPattern( symmHComplexMat_, symmPatternMat_ );
  //            Convert(symmPatternMat_,PatternMat_);
  //
  //            Convert(symmSComplexMat_,SComplexMat_);
  //            Convert(symmHComplexMat_,HComplexMat_);
  //          }
  //          break;
  //#endif
  //        default:
  //          ErrorHandling("Unsupported solver.");
  //          break;
  //      }

  return ;
}        // -----  end of method PPEXSIData::LoadComplexMatrix  -----

void
PPEXSIData::SymbolicFactorizeRealSymmetricMatrix    (
    Int                            solver,
    Int                            symmetricStorage,
    std::string                    ColPerm,
    Int                            numProcSymbFact,
    Int                            verbosity )
{
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  {
    if( verbosity >= 1 ){
      statusOFS << "Symbolic factorization for the real matrix."  << std::endl;
    }

    Real timeSta, timeEnd;

    PMatrix<Real>&          PMloc     = *PMRealMat_;
    SuperNodeType&          super     = superReal_;
    // Clear the PMatrix first
    PMloc = PMatrix<Real>();

    selinvOpt_.maxPipelineDepth = -1;
    selinvOpt_.symmetricStorage = symmetricStorage;
    factOpt_.ColPerm = ColPerm;

    switch (solver) {
      case 0:
        {
          // If we are using SuperLU
          SuperLUMatrix<Real>&    luMat     = *luRealMat_;
          // Clear the SuperLUMatrix first
          luMat = SuperLUMatrix<Real>();

          factOpt_.ColPerm = ColPerm;
          luOpt_.ColPerm = ColPerm;
          luOpt_.numProcSymbFact = numProcSymbFact;


          luMat.Setup( *gridSuperLUReal_, luOpt_ );  // SuperLU matrix.

          DistSparseMatrix<Real> AMat;
          CopyPattern( PatternMat_, AMat );
          SetValue( AMat.nzvalLocal, D_ZERO );          // Symbolic factorization does not need value


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
        }
        break;
#ifdef WITH_SYMPACK
      case 1:
        {
          if(symPACKRealMat_!=nullptr){
            delete symPACKRealMat_;
          }

          if( outputFileIndex_ >= 0 ){
            if(symPACK::logfileptr==NULL){
              auto mpirank = gridSelInv_->mpirank;
              //Initialize symPACK logfile
              std::stringstream suffix;
              suffix<<mpirank;
              symPACK::logfileptr = new symPACK::LogFile("status",suffix.str().c_str());
              symPACK::logfileptr->OFS()<<"********* LOGFILE OF P"<<mpirank<<" *********"<<std::endl;
              symPACK::logfileptr->OFS()<<"**********************************"<<std::endl;
            }
          }


          symPACKRealMat_ = new symPACK::symPACKMatrix<Real>();

          symPACK::symPACKOptions & optionsFact = symPACKOpt_;
          optionsFact.decomposition = symPACK::DecompositionType::LDL;
          optionsFact.orderingStr = ColPerm;
          optionsFact.MPIcomm = gridPole_->rowComm;
          optionsFact.NpOrdering = numProcSymbFact;

          symPACK::symPACKMatrix<Real>& symPACKMat = *symPACKRealMat_ ;
          symPACKMat.Init(optionsFact);

          //AMat.nzvalLocal.assign(AMat.nzvalLocal.size(), D_ZERO );          // Symbolic factorization does not need value
          symPACK::DistSparseMatrix<Real> ltAMat;
          Convert(PatternMat_, ltAMat);
          ltAMat.ToLowerTriangular();
          //              SetValue( AMat.nzvalLocal, D_ZERO );          // Symbolic factorization does not need value
          //ltAMat.nzvalLocal.assign(AMat.nzvalLocal.size(), D_ZERO );          // Symbolic factorization does not need value

          GetTime( timeSta );


          symPACKMat.SymbolicFactorization(ltAMat);
          GetTime( timeEnd );
          if( verbosity >= 1 ){
            statusOFS
              << "Time for symbolic factorization is "
              << timeEnd - timeSta << " [s]" << std::endl << std::endl;
          }

          symPACKMatrixToSuperNode( symPACKMat, super );
        }
        break;
#endif
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }



    PMloc.Setup( gridSelInv_, &super , &selinvOpt_, &factOpt_ );


    switch (solver) {
      case 0:
        {
          SuperLUMatrix<Real>&    luMat     = *luRealMat_;
          GetTime( timeSta );
          luMat.LUstructToPMatrix( PMloc );
          GetTime( timeEnd );
          if( verbosity >= 1 ){
            statusOFS
              << "Time for converting LUstruct to PMatrix is "
              << timeEnd - timeSta << " [s]" << std::endl << std::endl;
          }
        }
        break;
#ifdef WITH_SYMPACK
      case 1:
        {
          symPACK::symPACKMatrix<Real>& symPACKMat = *symPACKRealMat_ ;
          GetTime( timeSta );
          symPACKMatrixToPMatrix( symPACKMat, PMloc );
          GetTime( timeEnd );
          if( verbosity >= 1 ){
            statusOFS
              << "Time for converting symPACK Matrix to PMatrix is "
              << timeEnd - timeSta << " [s]" << std::endl << std::endl;
          }
        }
        break;
#endif
      default:
        ErrorHandling("Unsupported solver.");
        break;
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

    switch (solver) {
      case 0:
        {
          SuperLUMatrix<Real>&    luMat     = *luRealMat_;
          // Get ready for the factorization for another matrix using the same
          // sparsity pattern
          luMat.DestroyAOnly();
        }
        break;
#ifdef WITH_SYMPACK
      case 1:
        {
        }
        break;
#endif
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }


    if( verbosity >= 2 ){
      statusOFS << "perm: "    << std::endl << super.perm     << std::endl;
      statusOFS << "permInv: " << std::endl << super.permInv  << std::endl;
      statusOFS << "superIdx:" << std::endl << super.superIdx << std::endl;
      statusOFS << "superPtr:" << std::endl << super.superPtr << std::endl;
    }
  }

  isRealSymmetricSymbolicFactorized_ = true;


  return ;
}         // -----  end of method PPEXSIData::SymbolicFactorizeRealSymmetricMatrix  -----

void
PPEXSIData::SymbolicFactorizeRealUnsymmetricMatrix    (
    Int                            solver,
    std::string                    ColPerm,
    std::string                    RowPerm,
    Int                            numProcSymbFact,
    Int                            Transpose,
    double*                        AnzvalLocal,
    Int                            verbosity )
{
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  {
    if( verbosity >= 1 ){
      statusOFS << "Symbolic factorization for the real matrix."  << std::endl;
    }

    PMatrixUnsym<Real>&     PMloc     = *PMRealUnsymMat_;
    SuperNodeType&          super     = superReal_;
    // Clear the matrices first
    Real timeSta, timeEnd;
    PMloc = PMatrixUnsym<Real>();

    selinvOpt_.maxPipelineDepth = -1;
    factOpt_.ColPerm = ColPerm;
    factOpt_.RowPerm = RowPerm;
    factOpt_.Symmetric = 0;

    switch(solver){
      case 0:
        {
          SuperLUMatrix<Real>&    luMat     = *luRealMat_;
          // Clear the matrices first
          luMat = SuperLUMatrix<Real>();


          luOpt_.ColPerm = ColPerm;
          luOpt_.RowPerm = RowPerm;
          luOpt_.numProcSymbFact = numProcSymbFact;
          luOpt_.Symmetric = 0;
          luOpt_.Transpose = Transpose;

          luMat.Setup( *gridSuperLUReal_, luOpt_ );  // SuperLU matrix.

          DistSparseMatrix<Real> AMat;
          CopyPattern( PatternMat_, AMat );

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
              ErrorHandling( msg.str().c_str() );

            }
          }



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
        }
        break;
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }

    PMloc.Setup( gridSelInv_, &super , &selinvOpt_, &factOpt_ );

    switch (solver) {
      case 0:
        {
          SuperLUMatrix<Real>&    luMat     = *luRealMat_;
          GetTime( timeSta );
          luMat.LUstructToPMatrix( PMloc );
          GetTime( timeEnd );
          if( verbosity >= 1 ){
            statusOFS
              << "Time for converting LUstruct to PMatrix is "
              << timeEnd - timeSta << " [s]" << std::endl << std::endl;
          }
        }
        break;
      default:
        ErrorHandling("Unsupported solver.");
        break;
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

    switch (solver) {
      case 0:
        {
          SuperLUMatrix<Real>&    luMat     = *luRealMat_;
          // Get ready for the factorization for another matrix using the same
          // sparsity pattern
          luMat.DestroyAOnly();
        }
        break;
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }

    if( verbosity >= 2 ){
      statusOFS << "perm: "    << std::endl << super.perm     << std::endl;
      statusOFS << "permInv: " << std::endl << super.permInv  << std::endl;
      statusOFS << "superIdx:" << std::endl << super.superIdx << std::endl;
      statusOFS << "superPtr:" << std::endl << super.superPtr << std::endl;
    }
  }

  isRealUnsymmetricSymbolicFactorized_ = true;


  return ;
}         // -----  end of method PPEXSIData::SymbolicFactorizeRealUnsymmetricMatrix  -----

void
PPEXSIData::SymbolicFactorizeComplexSymmetricMatrix    (
    Int                            solver,
    Int                            symmetricStorage,
    std::string                    ColPerm,
    Int                            numProcSymbFact,
    Int                            verbosity )
{
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadComplexMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // Complex matrices
  {
    if( verbosity >= 1 ){
      statusOFS << "Symbolic factorization for the complex matrix."  << std::endl;
    }

    Real timeSta, timeEnd;
    PMatrix<Complex>&          PMloc     = *PMComplexMat_;
    SuperNodeType&             super     = superComplex_;

    // Clear the matrices first
    PMloc = PMatrix<Complex>();

    factOpt_.ColPerm = ColPerm;
    selinvOpt_.maxPipelineDepth = -1;
    selinvOpt_.symmetricStorage = symmetricStorage;

    switch (solver) {
      case 0:
        {
          SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
          // Clear the matrices first
          luMat = SuperLUMatrix<Complex>();

          luOpt_.ColPerm = ColPerm;
          luOpt_.numProcSymbFact = numProcSymbFact;

          luMat.Setup( *gridSuperLUComplex_, luOpt_ );  // SuperLU matrix.

          DistSparseMatrix<Complex> AMat;

          CopyPattern( PatternMat_, AMat );

          SetValue( AMat.nzvalLocal, Z_ZERO );          // Symbolic factorization does not need value

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

        }
        break;
#ifdef WITH_SYMPACK
      case 1:
        {
          if(symPACKComplexMat_!=nullptr){
            delete symPACKComplexMat_;
          }
          if( outputFileIndex_ >= 0 ){
            if(symPACK::logfileptr==NULL){
              //Initialize symPACK logfile
              auto mpirank = gridSelInv_->mpirank;
              std::stringstream suffix;
              suffix<<mpirank;
              symPACK::logfileptr = new symPACK::LogFile("status",suffix.str().c_str());
              symPACK::logfileptr->OFS()<<"********* LOGFILE OF P"<<mpirank<<" *********"<<std::endl;
              symPACK::logfileptr->OFS()<<"**********************************"<<std::endl;
            }
          }


          symPACKComplexMat_ = new symPACK::symPACKMatrix<Complex>();

          symPACK::symPACKOptions & optionsFact = symPACKOpt_;
          optionsFact.decomposition = symPACK::DecompositionType::LDL;
          optionsFact.orderingStr = ColPerm;
          optionsFact.MPIcomm = gridPole_->rowComm;
          optionsFact.NpOrdering = numProcSymbFact;

          symPACK::symPACKMatrix<Complex>& symPACKMat = *symPACKComplexMat_ ;
          symPACKMat.Init(optionsFact);

          DistSparseMatrix<Complex> AMat;
          CopyPattern(PatternMat_, AMat);
          symPACK::DistSparseMatrix<Complex> ltAMat;
          Convert(AMat, ltAMat);
          ltAMat.ToLowerTriangular();
          //SetValue( AMat.nzvalLocal, Z_ZERO );          // Symbolic factorization does not need value

          GetTime( timeSta );
          symPACKMat.SymbolicFactorization(ltAMat);
          GetTime( timeEnd );
          if( verbosity >= 1 ){
            statusOFS
              << "Time for symbolic factorization is "
              << timeEnd - timeSta << " [s]" << std::endl << std::endl;
          }

          symPACKMatrixToSuperNode( symPACKMat, super );
        }
        break;
#endif
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }


    PMloc.Setup( gridSelInv_, &super , &selinvOpt_, &factOpt_ );

    switch(solver){
      case 0:
        {
          SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
          GetTime( timeSta );
          luMat.LUstructToPMatrix( PMloc );
          GetTime( timeEnd );
          if( verbosity >= 1 ){
            statusOFS
              << "Time for converting LUstruct to PMatrix is "
              << timeEnd - timeSta << " [s]" << std::endl << std::endl;
          }
        }
        break;
#ifdef WITH_SYMPACK
      case 1:
        {
          symPACK::symPACKMatrix<Complex>& symPACKMat = *symPACKComplexMat_ ;
          GetTime( timeSta );
          symPACKMatrixToPMatrix( symPACKMat, PMloc );
          GetTime( timeEnd );
          if( verbosity >= 1 ){
            statusOFS
              << "Time for converting symPACK Matrix to PMatrix is "
              << timeEnd - timeSta << " [s]" << std::endl << std::endl;
          }
        }
        break;
#endif
      default:
        ErrorHandling("Unsupported solver.");
        break;
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


    switch (solver) {
      case 0:
        {
          SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
          // Get ready for the factorization for another matrix using the same
          // sparsity pattern
          luMat.DestroyAOnly();
        }
        break;
#ifdef WITH_SYMPACK
      case 1:
        {
        }
        break;
#endif
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }

    if( verbosity >= 2 ){
      statusOFS << "perm: "    << std::endl << super.perm     << std::endl;
      statusOFS << "permInv: " << std::endl << super.permInv  << std::endl;
      statusOFS << "superIdx:" << std::endl << super.superIdx << std::endl;
      statusOFS << "superPtr:" << std::endl << super.superPtr << std::endl;
    }
  }


  isComplexSymmetricSymbolicFactorized_ = true;


  return ;
}         // -----  end of method PPEXSIData::SymbolicFactorizeComplexSymmetricMatrix  -----

void
PPEXSIData::SymbolicFactorizeComplexUnsymmetricMatrix    (
    Int                            solver,
    std::string                    ColPerm,
    std::string                    RowPerm,
    Int                            numProcSymbFact,
    Int                            Transpose,
    double*                        AnzvalLocal,
    Int                            verbosity )
{
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadComplexMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // Complex matrices
  {
    if( verbosity >= 1 ){
      statusOFS << "Symbolic factorization for the complex matrix."  << std::endl;
    }

    Real timeSta, timeEnd;
    PMatrixUnsym<Complex>&     PMloc     = *PMComplexUnsymMat_;
    SuperNodeType&             super     = superComplex_;

    // Clear the matrices first
    PMloc = PMatrixUnsym<Complex>();

    selinvOpt_.maxPipelineDepth = -1;
    factOpt_.ColPerm = ColPerm;
    factOpt_.RowPerm = RowPerm;
    factOpt_.Symmetric = 0;

    switch(solver){
      case 0:
        {
          SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
          // Clear the matrices first
          luMat = SuperLUMatrix<Complex>();



          luOpt_.ColPerm = ColPerm;
          luOpt_.RowPerm = RowPerm;
          luOpt_.Symmetric = 0;
          luOpt_.numProcSymbFact = numProcSymbFact;
          luOpt_.Transpose = Transpose;

          luMat.Setup( *gridSuperLUComplex_, luOpt_ );  // SuperLU matrix.

          DistSparseMatrix<Complex> AMat;
          CopyPattern( PatternMat_, AMat );

          if(luOpt_.RowPerm == "LargeDiag"){
            if(AnzvalLocal != NULL){
              blas::Copy( 2*AMat.nnzLocal, AnzvalLocal, 1,
                  reinterpret_cast<double*>(AMat.nzvalLocal.Data()), 1 );
            }
            else{
              std::ostringstream msg;
              msg  << std::endl
                << "LargeDiag requires the non zero values to be provided." << std::endl;
              ErrorHandling( msg.str().c_str() );

            }
          }

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
        }
        break;
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }



    PMloc.Setup( gridSelInv_, &super , &selinvOpt_, &factOpt_ );



    switch(solver){
      case 0:
        {
          SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
          GetTime( timeSta );
          luMat.LUstructToPMatrix( PMloc );
          GetTime( timeEnd );
          if( verbosity >= 1 ){
            statusOFS
              << "Time for converting LUstruct to PMatrix is "
              << timeEnd - timeSta << " [s]" << std::endl << std::endl;
          }
        }
        break;
      default:
        ErrorHandling("Unsupported solver.");
        break;
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

    switch (solver) {
      case 0:
        {
          SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
          // Get ready for the factorization for another matrix using the same
          // sparsity pattern
          luMat.DestroyAOnly();
        }
        break;
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }

    if( verbosity >= 2 ){
      statusOFS << "perm: "    << std::endl << super.perm     << std::endl;
      statusOFS << "permInv: " << std::endl << super.permInv  << std::endl;
      statusOFS << "superIdx:" << std::endl << super.superIdx << std::endl;
      statusOFS << "superPtr:" << std::endl << super.superPtr << std::endl;
    }
  }


  isComplexUnsymmetricSymbolicFactorized_ = true;


  return ;
}         // -----  end of method PPEXSIData::SymbolicFactorizeComplexUnsymmetricMatrix  -----

void
PPEXSIData::SelInvRealSymmetricMatrix(
    Int               solver,
    Int               symmetricStorage,
    double*           AnzvalLocal,
    Int               verbosity,
    double*           AinvnzvalLocal )
{
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  if( isRealSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeRealSymmetricMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // Only the processor group corresponding to the first pole participate
  if( MYROW( gridPole_ ) == 0 ){

    Real timeTotalSelInvSta, timeTotalSelInvEnd;
    PMatrix<Real>&          PMloc     = *PMRealMat_;
    DistSparseMatrix<Real>& AMat      = shiftRealMat_;
    DistSparseMatrix<Real>& AinvMat   = shiftInvRealMat_;
    // Copy the pattern
    CopyPattern( PatternMat_, AMat );
    blas::Copy( AMat.nnzLocal, AnzvalLocal, 1, AMat.nzvalLocal.Data(), 1 );

    switch (solver) {
      case 0:
        {
          SuperLUMatrix<Real>&    luMat     = *luRealMat_;

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

          GetTime( timeTotalSelInvSta );

          luMat.LUstructToPMatrix( PMloc );
        }
        break;
#ifdef WITH_SYMPACK
      case 1:
        {
          symPACK::symPACKMatrix<Real>& symPACKMat = *symPACKRealMat_ ;
          symPACK::DistSparseMatrix<Real> ltAMat;
          if( verbosity >= 2 ){
            statusOFS << "Before ToLowerTriangular." << std::endl;
          }
          Convert(AMat,ltAMat);
          ltAMat.ToLowerTriangular();
          ltAMat.GetLocalGraph().SetSorted(false);
          ltAMat.SortGraph();
          if( verbosity >= 2 ){
            statusOFS << "After ToLowerTriangular." << std::endl;
          }


          Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

          GetTime( timeTotalFactorizationSta );

          // Data redistribution
          if( verbosity >= 2 ){
            statusOFS << "Before Distribute." << std::endl;
          }
          symPACKMat.DistributeMatrix(ltAMat);
          if( verbosity >= 2 ){
            statusOFS << "After Distribute." << std::endl;
          }

          // Numerical factorization
          if( verbosity >= 2 ){
            statusOFS << "Before NumericalFactorize." << std::endl;
          }
          symPACKMat.Factorize();
          // Numerical factorization
          if( verbosity >= 2 ){
            statusOFS << "After NumericalFactorize." << std::endl;
          }

          GetTime( timeTotalFactorizationEnd );

          if( verbosity >= 1 ){
            statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl;
          }

          //make a symmetric matrix out of that....
          GetTime( timeTotalSelInvSta );
          symPACKMatrixToPMatrix( symPACKMat, PMloc );
        }
        break;
#endif
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }

    PMloc.PreSelInv();

    PMloc.SelInv();

    GetTime( timeTotalSelInvEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total selected inversion is " <<
        timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
    }


    Real timePostProcessingSta, timePostProcessingEnd;

    GetTime( timePostProcessingSta );

    //TODO convert to symmAinvMat too
    PMloc.PMatrixToDistSparseMatrix( PatternMat_, AinvMat );

    GetTime( timePostProcessingEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for postprocessing is " <<
        timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
    }

    // Return the data to AinvnzvalLocal
    blas::Copy( AMat.nnzLocal, AinvMat.nzvalLocal.Data(), 1, AinvnzvalLocal, 1 );
  }


  return ;
}         // -----  end of method PPEXSIData::SelInvRealSymmetricMatrix  -----

void
PPEXSIData::SelInvRealUnsymmetricMatrix(
    Int               solver,
    double*           AnzvalLocal,
    Int               verbosity,
    double*           AinvnzvalLocal )
{
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  if( isRealUnsymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeRealUnsymmetricMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // Only the processor group corresponding to the first pole participate
  if( MYROW( gridPole_ ) == 0 ){

    PMatrixUnsym<Real>& PMloc         = *PMRealUnsymMat_;

    DistSparseMatrix<Real>& AMat      = shiftRealMat_;
    DistSparseMatrix<Real>& AinvMat   = shiftInvRealMat_;
    // Copy the pattern
    CopyPattern( PatternMat_, AMat );
    blas::Copy( AMat.nnzLocal, AnzvalLocal, 1, AMat.nzvalLocal.Data(), 1 );

    Real timeTotalSelInvSta, timeTotalSelInvEnd;
    switch (solver) {
      case 0:
        {
          SuperLUMatrix<Real>&    luMat     = *luRealMat_;


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

          GetTime( timeTotalSelInvSta );

          luMat.LUstructToPMatrix( PMloc );
        }
        break;
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }



    PMloc.PreSelInv();

    PMloc.SelInv();

    GetTime( timeTotalSelInvEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total selected inversion is " <<
        timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
    }


    Real timePostProcessingSta, timePostProcessingEnd;

    GetTime( timePostProcessingSta );

    PMloc.PMatrixToDistSparseMatrix( PatternMat_, AinvMat );

    GetTime( timePostProcessingEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for postprocessing is " <<
        timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
    }

    // Return the data to AinvnzvalLocal
    blas::Copy( AMat.nnzLocal, AinvMat.nzvalLocal.Data(), 1, AinvnzvalLocal, 1 );


  }


  return ;
}         // -----  end of method PPEXSIData::SelInvRealUnsymmetricMatrix  -----


void
PPEXSIData::SelInvComplexSymmetricMatrix(
    Int               solver,
    Int               symmetricStorage,
    double*           AnzvalLocal,
    Int               verbosity,
    double*           AinvnzvalLocal )
{
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadComplexMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  if( isComplexSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeComplexSymmetricMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // Only the processor group corresponding to the first pole participate
  if( MYROW( gridPole_ ) == 0 ){

    Real timeTotalSelInvSta, timeTotalSelInvEnd;
    PMatrix<Complex>&          PMloc     = *PMComplexMat_;
    DistSparseMatrix<Complex>& AMat      = shiftComplexMat_;
    DistSparseMatrix<Complex>& AinvMat   = shiftInvComplexMat_;
    // Copy the pattern
    CopyPattern( PatternMat_, AMat );
    blas::Copy( 2*AMat.nnzLocal, AnzvalLocal, 1,
        reinterpret_cast<double*>(AMat.nzvalLocal.Data()), 1 );

    switch (solver) {
      case 0:
        {
          SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;

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

          GetTime( timeTotalSelInvSta );

          luMat.LUstructToPMatrix( PMloc );
        }
        break;
#ifdef WITH_SYMPACK
      case 1:
        {
          symPACK::symPACKMatrix<Complex>& symPACKMat = *symPACKComplexMat_ ;
          symPACK::DistSparseMatrix<Complex> ltAMat;
          if( verbosity >= 2 ){
            statusOFS << "Before ToLowerTriangular." << std::endl;
          }
          Convert(AMat,ltAMat);
          ltAMat.ToLowerTriangular();
          ltAMat.GetLocalGraph().SetSorted(false);
          ltAMat.SortGraph();
          if( verbosity >= 2 ){
            statusOFS << "After ToLowerTriangular." << std::endl;
          }


          Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

          GetTime( timeTotalFactorizationSta );

          // Data redistribution
          if( verbosity >= 2 ){
            statusOFS << "Before Distribute." << std::endl;
          }
          symPACKMat.DistributeMatrix(ltAMat);
          if( verbosity >= 2 ){
            statusOFS << "After Distribute." << std::endl;
          }

          // Numerical factorization
          if( verbosity >= 2 ){
            statusOFS << "Before NumericalFactorize." << std::endl;
          }
          symPACKMat.Factorize();
          // Numerical factorization
          if( verbosity >= 2 ){
            statusOFS << "After NumericalFactorize." << std::endl;
          }

          GetTime( timeTotalFactorizationEnd );

          if( verbosity >= 1 ){
            statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl;
          }

          //make a symmetric matrix out of that....
          GetTime( timeTotalSelInvSta );
          symPACKMatrixToPMatrix( symPACKMat, PMloc );
        }
        break;
#endif
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }

    PMloc.PreSelInv();

    PMloc.SelInv();

    GetTime( timeTotalSelInvEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total selected inversion is " <<
        timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
    }


    Real timePostProcessingSta, timePostProcessingEnd;

    GetTime( timePostProcessingSta );

    //TODO convert to symmAinvMat too
    PMloc.PMatrixToDistSparseMatrix( PatternMat_, AinvMat );

    GetTime( timePostProcessingEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for postprocessing is " <<
        timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
    }

    // Return the data to AinvnzvalLocal
    blas::Copy( 2*AMat.nnzLocal, reinterpret_cast<double*>(AinvMat.nzvalLocal.Data()), 1,
        AinvnzvalLocal, 1 );
  }


  return ;
}         // -----  end of method PPEXSIData::SelInvComplexSymmetricMatrix  -----

void
PPEXSIData::SelInvComplexUnsymmetricMatrix(
    Int               solver,
    double*           AnzvalLocal,
    Int               verbosity,
    double*           AinvnzvalLocal )
{
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadComplexMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  if( isComplexUnsymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeComplexUnsymmetricMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // Only the processor group corresponding to the first pole participate
  if( MYROW( gridPole_ ) == 0 ){

    DistSparseMatrix<Complex>& AMat      = shiftComplexMat_;
    DistSparseMatrix<Complex>& AinvMat   = shiftInvComplexMat_;
    PMatrixUnsym<Complex>&     PMloc     = *PMComplexUnsymMat_;

    // Copy the pattern
    CopyPattern( PatternMat_, AMat );

    blas::Copy( 2*AMat.nnzLocal, AnzvalLocal, 1,
        reinterpret_cast<double*>(AMat.nzvalLocal.Data()), 1 );

    Real timeTotalSelInvSta, timeTotalSelInvEnd;
    switch (solver) {
      case 0:
        {
          SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
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

          GetTime( timeTotalSelInvSta );

          luMat.LUstructToPMatrix( PMloc );
        }
        break;
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }

    PMloc.PreSelInv();

    PMloc.SelInv();

    GetTime( timeTotalSelInvEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for total selected inversion is " <<
        timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
    }


    Real timePostProcessingSta, timePostProcessingEnd;

    GetTime( timePostProcessingSta );

    PMloc.PMatrixToDistSparseMatrix( PatternMat_, AinvMat );

    GetTime( timePostProcessingEnd );

    if( verbosity >= 1 ){
      statusOFS << "Time for postprocessing is " <<
        timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
    }

    // Return the data to AinvnzvalLocal
    blas::Copy( 2*AMat.nnzLocal, reinterpret_cast<double*>(AinvMat.nzvalLocal.Data()), 1,
        AinvnzvalLocal, 1 );
  }


  return ;
}         // -----  end of method PPEXSIData::SelInvComplexUnsymmetricMatrix  -----



void PPEXSIData::CalculateNegativeInertiaReal(
    const std::vector<Real>&       shiftVec,
    std::vector<Real>&             inertiaVec,
    Int                            solver,
    Int                            verbosity ){

  // *********************************************************************
  // Initialize
  // *********************************************************************
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  if( isRealSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeRealSymmetricMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }


  // Rename for convenience.
  // The symbolic information should have already been there.
  DistSparseMatrix<Real>&  HMat     = HRealMat_;
  DistSparseMatrix<Real>&  SMat     = SRealMat_;
  DistSparseMatrix<Real>& AMat      = shiftRealMat_;  // A = H - \lambda  S
  SuperLUMatrix<Real>&    luMat     = *luRealMat_;
  PMatrix<Real>&          PMloc     = *PMRealMat_;

  CopyPattern( PatternMat_, AMat );

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


      Real timeInertiaSta, timeInertiaEnd;
      // *********************************************************************
      // Factorization
      // *********************************************************************
      switch(solver){
        case 0:
          {
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

            GetTime( timeInertiaSta );
            luMat.LUstructToPMatrix( PMloc );
          }
          break;
#ifdef WITH_SYMPACK
        case 1:
          {
            Real timeTotalFactorizationSta, timeTotalFactorizationEnd;
            symPACK::symPACKMatrix<Real>& symPACKMat = *symPACKRealMat_ ;





            symPACK::DistSparseMatrix<Real> ltAMat;
            if( verbosity >= 2 ){
              statusOFS << "Before ToLowerTriangular." << std::endl;
            }
            Convert(AMat,ltAMat);
            ltAMat.ToLowerTriangular();
            if( verbosity >= 2 ){
              statusOFS << "After ToLowerTriangular." << std::endl;
            }

            GetTime( timeTotalFactorizationSta );
            // Data redistribution
            if( verbosity >= 2 ){
              statusOFS << "Before Distribute." << std::endl;
            }
            symPACKMat.DistributeMatrix(ltAMat);
            if( verbosity >= 2 ){
              statusOFS << "After Distribute." << std::endl;
            }

            // Numerical factorization
            if( verbosity >= 2 ){
              statusOFS << "Before NumericalFactorize." << std::endl;
            }
            symPACKMat.Factorize();
            // Numerical factorization
            if( verbosity >= 2 ){
              statusOFS << "After NumericalFactorize." << std::endl;
            }

            GetTime( timeTotalFactorizationEnd );

            if( verbosity >= 1 ){
              statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl;
            }

            GetTime( timeInertiaSta );
            symPACKMatrixToPMatrix( symPACKMat, PMloc );
          }
          break;
#endif
        default:
          ErrorHandling("Unsupported solver.");
          break;
      }

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

  //statusOFS<<"inertiaVecLocal "<<inertiaVecLocal<<std::endl;

  return ;
}         // -----  end of method PPEXSIData::CalculateNegativeInertiaReal -----

void PPEXSIData::CalculateNegativeInertiaComplex(
    const std::vector<Real>&       shiftVec,
    std::vector<Real>&             inertiaVec,
    Int                            solver,
    Int                            verbosity ){

  // *********************************************************************
  // Initialize
  // *********************************************************************
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadComplexMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // Use ComplexSymmetric factorization
  // FIXME This could change when interfacing with LDL^* factorization.
  if( isComplexSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeRealSymmetricMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }


  // Rename for convenience.
  // The symbolic information should have already been there.
  DistSparseMatrix<Complex>&  HMat     = HComplexMat_;
  DistSparseMatrix<Complex>&  SMat     = SComplexMat_;
  DistSparseMatrix<Complex>&  AMat     = shiftComplexMat_;  // A = H - \lambda  S
  SuperLUMatrix<Complex>&    luMat     = *luComplexMat_;
  PMatrix<Complex>&          PMloc     = *PMComplexMat_;

  CopyPattern( PatternMat_, AMat );

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

      Real timeInertiaSta, timeInertiaEnd;
      // *********************************************************************
      // Factorization
      // *********************************************************************
      switch(solver){
        case 0:
          {
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

            GetTime( timeInertiaSta );
            luMat.LUstructToPMatrix( PMloc );
          }
          break;
#ifdef WITH_SYMPACK
        case 1:
          {
            symPACK::symPACKMatrix<Complex>& symPACKMat = *symPACKComplexMat_ ;
            symPACK::DistSparseMatrix<Complex> ltAMat;
            if( verbosity >= 2 ){
              statusOFS << "Before ToLowerTriangular." << std::endl;
            }
            Convert(AMat,ltAMat);
            ltAMat.ToLowerTriangular();
            ltAMat.GetLocalGraph().SetSorted(false);
            ltAMat.SortGraph();
            if( verbosity >= 2 ){
              statusOFS << "After ToLowerTriangular." << std::endl;
            }


            Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

            GetTime( timeTotalFactorizationSta );
            // Data redistribution
            if( verbosity >= 2 ){
              statusOFS << "Before Distribute." << std::endl;
            }
            symPACKMat.DistributeMatrix(ltAMat);
            if( verbosity >= 2 ){
              statusOFS << "After Distribute." << std::endl;
            }

            // Numerical factorization
            if( verbosity >= 2 ){
              statusOFS << "Before NumericalFactorize." << std::endl;
            }
            symPACKMat.Factorize();
            // Numerical factorization
            if( verbosity >= 2 ){
              statusOFS << "After NumericalFactorize." << std::endl;
            }

            GetTime( timeTotalFactorizationEnd );

            if( verbosity >= 1 ){
              statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl;
            }

            GetTime( timeInertiaSta );
            symPACKMatrixToPMatrix( symPACKMat, PMloc );
          }
          break;
#endif
        default:
          ErrorHandling("Unsupported solver.");
          break;
      }


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

  return ;
}         // -----  end of method PPEXSIData::CalculateNegativeInertiaComplex -----

// Main subroutine for the electronic structure calculation
void PPEXSIData::CalculateFermiOperatorReal(
    Int   numPole,
    Real  temperature,
    Real  gap,
    Real  deltaE,
    Real  mu,
    Real  numElectronExact,
    Real  numElectronTolerance,
    Int   solver,
    Int   verbosity,
    Real& numElectron,
    Real& numElectronDrvMu ){

  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  if( isComplexSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeComplexSymmetricMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // *********************************************************************
  // Check the input parameters
  // *********************************************************************
  if( numPole % 2 != 0 ){
    statusOFS << " Warning, numPole is " << numPole <<" not a even number. "<< std::endl << std::flush;
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
  CopyPattern( PatternMat_, AMat );

  CopyPattern( PatternMat_, rhoMat );
  CopyPattern( PatternMat_, rhoDrvMuMat );
  if( isFreeEnergyDensityMatrix )
    CopyPattern( PatternMat_, hmzMat );
  if( isEnergyDensityMatrix )
    CopyPattern( PatternMat_, frcMat );
  if( isDerivativeTMatrix )
    CopyPattern( PatternMat_, rhoDrvTMat );

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
        statusOFS    << "zweightRho       = " << zweightRho_[l] << std::endl;
        statusOFS    << "zweightRhoDrvMu  = " << zweightRhoDrvMu_[l] << std::endl;
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





        Real timeTotalSelInvSta, timeTotalSelInvEnd;
        // *********************************************************************
        // Factorization
        // *********************************************************************
        switch(solver){
          case 0:
            {
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
              GetTime( timeTotalSelInvSta );

              luMat.LUstructToPMatrix( PMloc );
            }
            break;
#ifdef WITH_SYMPACK
          case 1:
            {
              symPACK::symPACKMatrix<Complex>& symPACKMat = *symPACKComplexMat_ ;
              symPACK::DistSparseMatrix<Complex> ltAMat;
              if( verbosity >= 2 ){
                statusOFS << "Before ToLowerTriangular." << std::endl;
              }
              Convert(AMat,ltAMat);
              ltAMat.ToLowerTriangular();
              ltAMat.GetLocalGraph().SetSorted(false);
              ltAMat.SortGraph();
              if( verbosity >= 2 ){
                statusOFS << "After ToLowerTriangular." << std::endl;
              }

              Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

              GetTime( timeTotalFactorizationSta );
              // Data redistribution
              if( verbosity >= 2 ){
                statusOFS << "Before Distribute." << std::endl;
              }
              symPACKMat.DistributeMatrix(ltAMat);
              if( verbosity >= 2 ){
                statusOFS << "After Distribute." << std::endl;
              }

              // Numerical factorization
              if( verbosity >= 2 ){
                statusOFS << "Before NumericalFactorize." << std::endl;
              }
              symPACKMat.Factorize();
              // Numerical factorization
              if( verbosity >= 2 ){
                statusOFS << "After NumericalFactorize." << std::endl;
              }

              GetTime( timeTotalFactorizationEnd );

              if( verbosity >= 1 ){
                statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl;
              }

              GetTime( timeTotalSelInvSta );
              symPACKMatrixToPMatrix( symPACKMat, PMloc );
            }
            break;
#endif
          default:
            ErrorHandling("Unsupported solver.");
            break;
        }



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

        //TODO convert to symmAinvMat too
        PMloc.PMatrixToDistSparseMatrix( PatternMat_, AinvMat );

        if( verbosity >= 2 ){
          statusOFS << "rhoMat.nnzLocal = " << rhoMat.nnzLocal << std::endl;
          statusOFS << "AinvMat.nnzLocal = " << AinvMat.nnzLocal << std::endl;
        }


        // Update the density matrix. The following lines are equivalent to
        //
        //                for( Int i = 0; i < rhoMat.nnzLocal; i++ ){
        //                    rhoMat.nzvalLocal(i) +=
        //                        zweightRho_[l].real() * AinvMat.nzvalLocal(i).imag() +
        //                        zweightRho_[l].imag() * AinvMat.nzvalLocal(i).real();
        //                }
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


  // Compute the energy, and free energy
  {
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
    if( isEnergyDensityMatrix )
    {
      Real local = 0.0;
      if( SMat.size != 0 ){
        local = blas::Dot( SMat.nnzLocal,
            SMat.nzvalLocal.Data(),
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
    if( isFreeEnergyDensityMatrix )
    {
      Real local = 0.0;
      if( SMat.size != 0 ){
        local = blas::Dot( SMat.nnzLocal,
            SMat.nzvalLocal.Data(),
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
      totalFreeEnergy_ += mu * numElectron;
    }
  }
  return ;
}    // -----  end of method PPEXSIData::CalculateFermiOperatorReal  -----



// Main subroutine for the electronic structure calculation with Hermitian Hamiltonian
// and overlap matrices
void PPEXSIData::CalculateFermiOperatorComplex(
    Int   numPole,
    Real  temperature,
    Real  gap,
    Real  deltaE,
    Real  mu,
    Real  numElectronExact,
    Real  numElectronTolerance,
    Int   solver,
    Int   verbosity,
    Real& numElectron,
    Real& numElectronDrvMu,
    Int   method,
    Int   nPoints,
    Real  spin) {

  // These are needed in the point-pole-pexsi parallelization
  Int npPerPoint  = gridPole_->mpisize / nPoints;
  if(npPerPoint < 1) {
    statusOFS << " Error, Point parallelization error, npPerPoint < 1 " << std::endl;
    statusOFS << " mpisize should be at least equal to " << nPoints << std::endl;
  }
  Int myPoint     = gridPole_->mpirank / npPerPoint;
  Int myPointRank = gridPole_->mpirank % npPerPoint;
  Int myRowPoint  = myPointRank / gridPole_->numProcCol;
  MPI_Comm pointColComm, pointRowComm;
  MPI_Comm_split( gridPole_->colComm, myPoint, myPointRank, &pointColComm);
  MPI_Comm_split( gridPole_->colComm, myPointRank, myPoint, &pointRowComm);

  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadComplexMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  if( isComplexUnsymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeComplexUnsymmetricMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // *********************************************************************
  // Check the input parameters
  // *********************************************************************
  if( numPole % 2 != 0 ){
    //ErrorHandling( "Must be even number of poles!" );
    statusOFS << " Warning, numPole is " << numPole <<" not a even number. "<< std::endl << std::flush;
  }

  // *********************************************************************
  // Initialize
  // *********************************************************************
  // Rename for convenience
  DistSparseMatrix<Complex>&  HMat        = HComplexMat_;
  DistSparseMatrix<Complex>&  SMat        = SComplexMat_;

  DistSparseMatrix<Complex>& rhoMat       = rhoComplexMat_;
  DistSparseMatrix<Complex>& rhoDrvMuMat  = rhoDrvMuComplexMat_;
  DistSparseMatrix<Complex>& rhoDrvTMat   = rhoDrvTComplexMat_;
  DistSparseMatrix<Complex>& hmzMat       = freeEnergyDensityComplexMat_;
  DistSparseMatrix<Complex>& frcMat       = energyDensityComplexMat_;

  DistSparseMatrix<Complex>& AMat      = shiftComplexMat_;
  DistSparseMatrix<Complex>& AinvMat   = shiftInvComplexMat_;

  // The symbolic information should have already been there.
  SuperLUMatrix<Complex>& luMat        = *luComplexMat_;
  PMatrixUnsym<Complex>&  PMloc        = *PMComplexUnsymMat_;

  //
  bool isFreeEnergyDensityMatrix = true;
  bool isEnergyDensityMatrix     = true;
  bool isDerivativeTMatrix       = false;

  if( isDerivativeTMatrix == true )
    ErrorHandling("isDerivativeTMatrix has to be false.");

  Real numSpin = spin;

  // Copy the pattern
  CopyPattern( PatternMat_, AMat );
  CopyPattern( PatternMat_, rhoMat );
  fflush(stdout);
  if( isFreeEnergyDensityMatrix )
    CopyPattern( PatternMat_, hmzMat );
  if( isEnergyDensityMatrix )
    CopyPattern( PatternMat_, frcMat );

  // Reinitialize the variables
  SetValue( rhoMat.nzvalLocal, Z_ZERO );
  if( isFreeEnergyDensityMatrix )
    SetValue( hmzMat.nzvalLocal, Z_ZERO );
  if( isEnergyDensityMatrix )
    SetValue( frcMat.nzvalLocal, Z_ZERO );

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
  //
  // This part should be kept the same as CalculateFermiOperatorReal
  //  std::vector<Int>  poleIdx(numPole);
  {
    // Setup a grid from (mu - deltaE, mu + deltaE), and measure
    // the error bound in the L^infty sense on this grid.
    //
    // fdGrid:      Exact Fermi-Dirac function evaluated on xGrid
    // fdPoleGrid:  Fermi-Dirac function using pole expansion
    // evaluated on the grid.
    //    Int numX = 10000;
    //    std::vector<Real>    xGrid( numX );
    //    std::vector<Real>    fdGrid( numX );
    //
    //    Real x0 = mu - deltaE;
    //    Real x1 = mu + deltaE;
    //    Real h  = (x1 - x0) / (numX - 1);
    //    Real ez;
    //    for( Int i = 0; i < numX; i++ ){
    //      xGrid[i]  = x0 + i * h;
    //      if( xGrid[i] - mu >= 0 ){
    //        ez = std::exp(- (xGrid[i] - mu) / temperature );
    //        fdGrid[i] = 2.0 * ez / (1.0 + ez);
    //      }
    //      else{
    //        ez = std::exp((xGrid[i] - mu) / temperature );
    //        fdGrid[i] = 2.0 / (1.0 + ez);
    //      }
    //    }


    numPoleInput = numPole;
    Real tol;
    Int  numPoleSignificant;

    Int poleIter = 0;
    /// test the getPole function. actually, NumPole should be get from here.
    {
      zweightFDM_.resize(numPole);
      zweightEDM_.resize(numPole);

      // first get the real numPole
      poleClass pole;
      bool result;
      // Method 1: Contour integral
      // Method 2: Moussa's pole expansion
      // Method 3: AAA pole expansion
      if(method == 1 || method == 2 || method == 3)
        result = pole.getPole(method, numPole,deltaE/temperature, zshift_, zweightRho_, zweightFDM_, zweightEDM_);
      else{
        statusOFS << " Error, method must be chosen from [ 1, 2, 3] " << std::endl;
        statusOFS << " Method 1: Contour integral" << std::endl;
        statusOFS << " Method 2: Moussa's pole expansion" << std::endl;
        statusOFS << " Method 3: AAA pole expansion" << std::endl;
        statusOFS << " Error, Current method is: " << method << std::endl;
        ErrorHandling("getPole error.");
      }
      if(method == 2)
        isEDMCorrection_ = 1;
      else
        isEDMCorrection_ = 0;

      for(int i = 0; i < zshift_.size(); i++)
      {
        zshift_[i] *= temperature;
        zweightRho_[i] *= temperature;
        zshift_[i] += mu;
        if( isEDMCorrection_ == 0){
          zweightFDM_[i] *= temperature * temperature;
          zweightEDM_[i] *= temperature * temperature;
        }
      }

      if(verbosity >= 2){
        statusOFS << " numPole is "<< numPole << " deltaE  "
          << deltaE/ temperature << std::endl;
        statusOFS << " temperature " << temperature << std::endl;
        statusOFS << " zshift "<< zshift_ << " zweight "
          << zweightRho_ << std::endl;
      }
      if(!result)
      {
        statusOFS << " Eroror getPole wrong method: " << method
          <<" " <<numPole << " " << deltaE/temperature << std::endl;
        statusOFS << " numPoles should be chosen from [ 10, 15, 20, 25, 30, 35, 40] " << std::endl;
        statusOFS << " Please also check whether the temperature is set to be too high or too low " << std::endl;
        ErrorHandling("getPole error.");
      }
      if(numPole != zshift_.size()){
        statusOFS << "  Error: numPole should be changed to "
          << zshift_.size()   << std::endl;
        ErrorHandling("getPole error.");
      }
    }
  }


  // Initialize the number of electrons
  numElectron  = 0.0;
  if( verbosity >= 2 ){
    statusOFS << "zshift" << std::endl << zshift_ << std::endl;
    statusOFS << "zweightRho" << std::endl << zweightRho_ << std::endl;
  }

  // for each pole, perform LDLT factoriation and selected inversion
  Real timePoleSta, timePoleEnd;

  Int numPoleComputed = 0;
  for(Int lidx = 0; lidx < numPole; lidx++){
    if( myRowPoint % numPole == lidx % (gridPole_->numProcRow /nPoints)){

      //Int l = poleIdx[lidx];
      Int l = lidx;

      GetTime( timePoleSta );

      if( verbosity >= 1 ){
        statusOFS << "Pole " << lidx << " processing..." << std::endl;
      }
      if( verbosity >= 2 ){
        statusOFS << "zshift           = " << zshift_[l] << std::endl;
        statusOFS    << "zweightRho       = " << zweightRho_[l] << std::endl;
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

        Real timeTotalSelInvSta, timeTotalSelInvEnd;
        switch (solver) {
          case 0:
            {
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
              GetTime( timeTotalSelInvSta );

              luMat.LUstructToPMatrix( PMloc );
            }
            break;
#ifdef WITH_SYMPACK
          case 1:
            {
              symPACK::symPACKMatrix<Complex>& symPACKMat = *symPACKComplexMat_ ;
              symPACK::DistSparseMatrix<Complex> ltAMat;
              if( verbosity >= 2 ){
                statusOFS << "Before ToLowerTriangular." << std::endl;
              }
              Convert(AMat,ltAMat);
              ltAMat.ToLowerTriangular();
              ltAMat.GetLocalGraph().SetSorted(false);
              ltAMat.SortGraph();
              if( verbosity >= 2 ){
                statusOFS << "After ToLowerTriangular." << std::endl;
              }


              Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

              GetTime( timeTotalFactorizationSta );

              // Data redistribution
              if( verbosity >= 2 ){
                statusOFS << "Before Distribute." << std::endl;
              }
              symPACKMat.DistributeMatrix(ltAMat);
              if( verbosity >= 2 ){
                statusOFS << "After Distribute." << std::endl;
              }

              // Numerical factorization
              if( verbosity >= 2 ){
                statusOFS << "Before NumericalFactorize." << std::endl;
              }
              symPACKMat.Factorize();
              // Numerical factorization
              if( verbosity >= 2 ){
                statusOFS << "After NumericalFactorize." << std::endl;
              }

              GetTime( timeTotalFactorizationEnd );

              if( verbosity >= 1 ){
                statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl;
              }

              //make a symmetric matrix out of that....
              GetTime( timeTotalSelInvSta );
              symPACKMatrixToPMatrix( symPACKMat, PMloc );
            }
            break;
#endif
          default:
            ErrorHandling("Unsupported solver.");
            break;
        }


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

        //TODO convert to symmAinvMat too
        PMloc.PMatrixToDistSparseMatrix( PatternMat_, AinvMat );

        if( verbosity >= 2 ){
          statusOFS << "rhoMat.nnzLocal = " << rhoMat.nnzLocal << std::endl;
          statusOFS << "AinvMat.nnzLocal = " << AinvMat.nnzLocal << std::endl;
        }

        // AinvMat now stores the selected elements of G_l^{T}; now get G_l and store it in transMat
        // using CSCToCSR
        DistSparseMatrix<Complex>& regMat   = AinvMat;
        DistSparseMatrix<Complex>  transMat;

        CSCToCSR( regMat, transMat );
        Complex* ptrReg   = regMat.nzvalLocal.Data();
        Complex* ptrTrans = transMat.nzvalLocal.Data();

        // fac = 1/(2i)
        Complex  fac1 = Complex( 0.0, -0.5 );
        Complex  fac2 = Complex( 0.0, 0.5 );

        // Update the density matrix <==  1/(2i)*w_l * G_l
        blas::Axpy( rhoMat.nnzLocal, fac1*numSpin*zweightRho_[l], transMat.nzvalLocal.Data(), 1,
            rhoMat.nzvalLocal.Data(), 1 );

        // Update the density matrix <==  -1/(2i) * conj(w_l) * G_l^{T}
        Complex* ptrRho = rhoMat.nzvalLocal.Data();
        for( Int g = 0; g < rhoMat.nnzLocal; g++ ){
          ptrRho[g] += fac2 * numSpin * std::conj(zweightRho_[l]) * std::conj(ptrReg[g]);
        }

        // Free energy density matrix
        if( isFreeEnergyDensityMatrix ){
          if( isEDMCorrection_ == 0){
            blas::Axpy( hmzMat.nnzLocal, fac1*numSpin*zweightFDM_[l],
              transMat.nzvalLocal.Data(), 1,
              hmzMat.nzvalLocal.Data(), 1 );

            Complex* ptrFDM = hmzMat.nzvalLocal.Data();
            for( Int g = 0; g < hmzMat.nnzLocal; g++ ){
              ptrFDM[g] += fac2 * numSpin * std::conj(zweightFDM_[l]) * std::conj(ptrReg[g]);
            }

          }
        }

        // Energy density matrix
        if( isEnergyDensityMatrix ){
          if( isEDMCorrection_ == 0){
            blas::Axpy( frcMat.nnzLocal, fac1*zweightEDM_[l]*numSpin,
                transMat.nzvalLocal.Data(), 1,
                frcMat.nzvalLocal.Data(), 1 );

            Complex* ptrEDM = frcMat.nzvalLocal.Data();
            for( Int g = 0; g < frcMat.nnzLocal; g++ ){
              ptrEDM[g] += fac2 * numSpin * std::conj(zweightEDM_[l]) * std::conj(ptrReg[g]);
            }
          }
          else{
            blas::Axpy( frcMat.nnzLocal, fac1*zweightRho_[l]*zshift_[l]*numSpin,
                transMat.nzvalLocal.Data(), 1,
                frcMat.nzvalLocal.Data(), 1 );

            Complex* ptrEDM = frcMat.nzvalLocal.Data();
            for( Int g = 0; g < frcMat.nnzLocal; g++ ){
              ptrEDM[g] += fac2 * numSpin * std::conj(zweightRho_[l]) *std::conj(zshift_[l]) * std::conj(ptrReg[g]);
            }

          }
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
    CpxNumVec nzvalRhoMatLocal = rhoMat.nzvalLocal;
    SetValue( rhoMat.nzvalLocal, Z_ZERO );

    mpi::Allreduce( nzvalRhoMatLocal.Data(), rhoMat.nzvalLocal.Data(),
        //rhoMat.nnzLocal, MPI_SUM, gridPole_->colComm );
      rhoMat.nnzLocal, MPI_SUM, pointColComm);

  }

  // Reduce the free energy density matrix across the processor rows in gridPole_
  if( isFreeEnergyDensityMatrix ){
    CpxNumVec nzvalHmzMatLocal = hmzMat.nzvalLocal;
    SetValue( hmzMat.nzvalLocal, Z_ZERO );

    mpi::Allreduce( nzvalHmzMatLocal.Data(), hmzMat.nzvalLocal.Data(),
      //hmzMat.nnzLocal, MPI_SUM, gridPole_->colComm );
      hmzMat.nnzLocal, MPI_SUM, pointColComm );

  }

  // Reduce the energy density matrix across the processor rows in gridPole_
  if( isEnergyDensityMatrix ){

    CpxNumVec nzvalFrcMatLocal = frcMat.nzvalLocal;
    SetValue( frcMat.nzvalLocal, Z_ZERO );

    mpi::Allreduce( nzvalFrcMatLocal.Data(), frcMat.nzvalLocal.Data(),
        //frcMat.nnzLocal, MPI_SUM, gridPole_->colComm );
      frcMat.nnzLocal, MPI_SUM, pointColComm);

    if( isEDMCorrection_ == 0){
      // for method 1 and 3, add mu * DM to EDM
      blas::Axpy( rhoMat.nnzLocal, mu, rhoMat.nzvalLocal.Data(), 1,
        frcMat.nzvalLocal.Data(), 1 );
    }

  }


  // Compute the number of electrons
  // The number of electrons is computed by Tr[DM*S]
  {
    Complex numElecLocal = Z_ZERO;

    if( SMat.size != 0 ){
      // S is not an identity matrix
      numElecLocal = blas::Dotc( SMat.nnzLocal, SMat.nzvalLocal.Data(),
          1, rhoMat.nzvalLocal.Data(), 1 );
    }
    else{
      // S is an identity matrix
      CpxNumVec& nzval = rhoMat.nzvalLocal;
      for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
        numElecLocal += nzval(diagIdxLocal_[i]);
      }
    } // if ( SMat.size != 0 )
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << std::endl << "numElecLocal = " << numElecLocal << std::endl;
#endif

    Complex numElectronCpx;
    mpi::Allreduce( &numElecLocal, &numElectronCpx, 1, MPI_SUM, gridPole_->rowComm );
    //mpi::Allreduce( &numElecLocal, &numElectronCpx, 1, MPI_SUM, rhoMat.comm );
    //mpi::Allreduce( &numElecLocal, &numElectron, 1, MPI_SUM, gridPole_->rowComm);

    numElectron = numElectronCpx.real();
  }


  // Compute the energy, and free energy
  // NOTE: the free energy is inaccessible if isEDMCorrection_ == 1
  {
    // Energy computed from Tr[H*DM]
    {
      Real local = 0.0;
      local = (blas::Dotc( HComplexMat_.nnzLocal,
            HComplexMat_.nzvalLocal.Data(),
            1, rhoComplexMat_.nzvalLocal.Data(), 1 )).real();
      mpi::Allreduce( &local, &totalEnergyH_, 1, MPI_SUM,
          gridPole_->rowComm );
    }

    if( isEDMCorrection_ == 0)
    {
      // Energy computed from Tr[S*EDM]
      if( isEnergyDensityMatrix )
      {
        Real local = 0.0;
        if( SMat.size != 0 ){
          local = (blas::Dotc( SMat.nnzLocal,
                SMat.nzvalLocal.Data(),
                1, energyDensityComplexMat_.nzvalLocal.Data(), 1 )).real();
        }
        else{
          CpxNumVec& nzval = energyDensityComplexMat_.nzvalLocal;
          for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
            local += nzval(diagIdxLocal_[i]).real();
          }
        }

        mpi::Allreduce( &local, &totalEnergyS_, 1, MPI_SUM,
            gridPole_->rowComm );
      }

      // Free energy
      if( isFreeEnergyDensityMatrix )
      {
        Real local = 0.0;
        if( SMat.size != 0 ){
          local = (blas::Dotc( SMat.nnzLocal,
                SMat.nzvalLocal.Data(),
                1, freeEnergyDensityComplexMat_.nzvalLocal.Data(), 1 )).real();
          }
          else{
            CpxNumVec& nzval = freeEnergyDensityComplexMat_.nzvalLocal;
            for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
              local += nzval(diagIdxLocal_[i]).real();
            }
          }

          mpi::Allreduce( &local, &totalFreeEnergy_, 1, MPI_SUM,
              gridPole_->rowComm );

          // Correction
          totalFreeEnergy_ += mu * numElectron;
      }
    }
  }
  MPI_Comm_free(&pointColComm);
  MPI_Comm_free(&pointRowComm);

  return ;
}    // -----  end of method PPEXSIData::CalculateFermiOperatorComplex   -----



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
    Int        solver,
    Int        symmetricStorage,
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
  switch (solver){
    case 0:
      {
        //Handle SuperLU ordering options
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
        switch (ordering){
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



  if( matrixType != 0 ){
    ErrorHandling("Unsupported matrixType. The variable has to be 0.");
  }

  // Perform symbolic factorization first if required
  if( isSymbolicFactorize == true ){
    if( verbosity >= 1 ){
      statusOFS << "Perform symbolic factorization now." << std::endl;
    }
    if( matrixType == 0 ){
      SymbolicFactorizeRealSymmetricMatrix(
          solver,
          symmetricStorage,
          colPerm,
          numProcSymbFact,
          verbosity );

      SymbolicFactorizeComplexSymmetricMatrix(
          solver,
          symmetricStorage,
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
          ErrorHandling( msg.str().c_str() );
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
              solver,
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
            solver,
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



  return ;
}         // -----  end of method PPEXSIData::DFTDriver  -----


void
PPEXSIData::DFTDriver2 (
    Real       numElectronExact,
    Real       temperature,
    Real       gap,
    Real       deltaE,
    Int        numPole,
    Real       muInertiaTolerance,
    Real       numElectronPEXSITolerance,
    Int        matrixType,
    Int        isSymbolicFactorize,
    Int        solver,
    Int        symmetricStorage,
    Int        ordering,
    Int        numProcSymbFact,
    Int        verbosity,
    Real&      muPEXSI,
    Real&      numElectronPEXSI,
    Real&      muMinInertia,
    Real&      muMaxInertia,
    Int&       numTotalInertiaIter,
    Int        method,
    Int        nPoints,
    Real       spin)
{
  Real timeSta, timeEnd;
  Real timeInertiaSta, timeInertiaEnd;
  Real timePEXSISta, timePEXSIEnd;
  Real timeTotalSta, timeTotalEnd;
  Real timeInertia = 0.0;
  Real timePEXSI   = 0.0;


  // Initial setup
  Real muMin = muMinInertia;
  Real muMax = muMaxInertia;
  //muPEXSI = mu0;

  // parameters.
  // Real K2au = 3.166815E-6;
  // Real T = temperature;
  // Real Beta = 1.0 / temperature;
  Real sigma = 3* temperature; // CHECK CHECK
  Real numElecTol = 1.0e-5;

  // add the points parallelization.
  if(gridPole_->mpisize % nPoints){
    ErrorHandling("nPoints can not divided by MPI used in PEXSI.");
  }
  if(gridPole_->mpisize / nPoints > numPole * gridPole_-> numProcRow * gridPole_-> numProcCol ){
    ErrorHandling("use more points in the PEXSI parallelization.");
  }

  // These are needed in the point-pole-pexsi parallelization
  Int npPerPoint  = gridPole_->mpisize / nPoints;
  Int myPoint     = gridPole_->mpirank / npPerPoint;
  Int myPointRank = gridPole_->mpirank % npPerPoint;

  if( verbosity >= 1 ){
    statusOFS << gridPole_->numProcRow << " numProcCol: " <<  gridPole_->numProcCol << " numPole " << numPole<< std::endl;
    statusOFS << gridPole_->mpisize << " "  << gridPole_->mpirank << " myPoint "  <<  myPointRank<< " nPoints " << nPoints << std::endl;
  }


  // not used in current code.
  Real muLower;
  Real muUpper;


  GetTime( timeTotalSta );

  std::string colPerm;
  switch (solver){
    case 0:
      {
        //Handle SuperLU ordering options
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
        switch (ordering){
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


  if( matrixType != 0 ){
    ErrorHandling("Unsupported matrixType. The variable has to be 0.");
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
          solver,
          symmetricStorage,
          colPerm,
          numProcSymbFact,
          verbosity );

      SymbolicFactorizeComplexSymmetricMatrix(
          solver,
          symmetricStorage,
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


  numTotalInertiaIter = 0;
  if( verbosity >= 1 ){
    statusOFS << std::endl
      << "input (muMin, muMax)   = " << "(" << muMin<< ", " << muMax
      << ")" << std::endl;
  }


  if(muMax - muMin > muInertiaTolerance){

    // do the Inertia Counting
    // Inertia counting loop
    // Hard coded: max inertia counting 10 steps
    const Int maxTotalInertiaIter = 10;

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

      std::vector<Real> NeLower(numShift);
      std::vector<Real> NeUpper(numShift);
      Real hsShift = ( muMax - muMin ) / (numShift - 1);

      Int matrix_size = HRealMat_.size;  // CHECK CHECK
      Int numSpin = spin;
      for(Int l = 0; l < numShift; l ++)
      {
        NeLower[l] = 0.0;
        inertiaVec[l] = 0.0;
        NeUpper[l] = matrix_size * numSpin;
        shiftVec[l] = muMin + l * hsShift;
      }

      if( matrixType == 0 ){
        CalculateNegativeInertiaReal(
            shiftVec,
            inertiaVec,
            solver,
            verbosity );
      }

      // FIXME In the future numSpin should be introduced.
      for( Int l = 0; l < numShift; l++ ){
        inertiaVec[l] *= numSpin;
      }

      Int Idx = (Int) std::ceil( sigma / hsShift ) ;

      for(Int l = Idx; l < numShift; l++)
      {
        NeLower[l]     = 0.5 * ( inertiaVec[l-Idx] + inertiaVec[l] );
        NeUpper[l-Idx] = 0.5 * ( inertiaVec[l-Idx] + inertiaVec[l] );
      }

      Int idxMin = 0;
      Int idxMax = numShift-1;
      for(Int l = 1; l < numShift-1; l++)
      {
        if( ( NeUpper[l] < numElectronExact ) && ( NeUpper[l+1] >= numElectronExact ) )
          idxMin = l;
        if( ( NeLower[l] > numElectronExact ) && ( NeLower[l-1] <= numElectronExact ) )
          idxMax = l;
      }

      if( verbosity >= 1 ){

        statusOFS << std::endl;


        for(Int l = 0; l < numShift; l++)
        {
          statusOFS << std::setiosflags(std::ios::left)
            << "Shift " << std::setw(4) << l
            << " = "
            << std::setiosflags(std::ios::scientific)
            << std::setiosflags(std::ios::showpos)
            << std::setw(LENGTH_DBL_DATA) << shiftVec[l]
            << std::resetiosflags(std::ios::scientific)
            << std::resetiosflags(std::ios::showpos)
            << std::setw(LENGTH_VAR_NAME) << "Inertia    = "
            << std::setw(LENGTH_VAR_DATA) << inertiaVec[l]
            << std::setw(LENGTH_VAR_NAME) << "NeLower  = "
            << std::setw(LENGTH_VAR_DATA) << NeLower[l]
            << std::setw(LENGTH_VAR_NAME) << "NeUpper  = "
            << std::setw(LENGTH_VAR_DATA) << NeUpper[l]
            << std::endl;
        }

        statusOFS << std::endl;
        statusOFS << "idxMin = " << idxMin << ", inertiaVec = " << inertiaVec[idxMin] << std::endl;
        statusOFS << "idxMax = " << idxMax << ", inertiaVec = " << inertiaVec[idxMax] << std::endl;
      }

      muMin = shiftVec[idxMin];
      muMax = shiftVec[idxMax];

      if ( ( ( idxMin == 0 ) && (idxMax == numShift-1 ) ) ||
          ( NeLower[idxMin] == 0 && NeUpper[idxMax] == matrix_size) ||
          muMax - muMin < muInertiaTolerance ) {

        muMinInertia = shiftVec[idxMin];
        muMaxInertia = shiftVec[idxMax];

        if( verbosity >= 1 ){
          statusOFS << std::endl << std::endl
            << "Inertia counting has converged." << std::endl
            << "(muMin, muMax) = ( " << shiftVec[idxMin] << " , " << shiftVec[idxMax] << " ) " << std::endl
            << "(Ne(muMin), Ne(muMax)) = ( " << NeLower[idxMin] << " , " << NeUpper[idxMax]
            << " ) " << std::endl
            << "NeExact = " << numElectronExact << std::endl
            << "Leave the inertia counting "<< std::endl << std::endl;
        }

        break;
      }

      numTotalInertiaIter++;

      if( verbosity >= 1 ){
        statusOFS << std::endl
          << "Inertia Counting " << std::endl
          << "(muMin, muMax)   = " << "(" << muMin<< ", " << muMax
          << ")" << std::endl
          << "numShift           = " << numShift << std::endl;
      }

    } // while (inertiaCount)
  } // do the Inertia Counting

  // PEXSI phase.
  // No option for falling back to inertia counting
  //
  {

    Int numShift = nPoints;//std::max( gridPole_->numProcRow/ numPole, 1);
    if( verbosity >= 1 ) {
      statusOFS << "  >> total Number of Points : " << numShift << std::endl;
      statusOFS << "  >> Pole Expansion Method  : " << method   << std::endl;
      statusOFS << "  >> number of Poles used   : " << numPole  << std::endl;
      statusOFS << std::endl;
    }

    std::vector<Real>  shiftVec( numShift );
    Real hsShift = ( muMax - muMin ) / ( numShift + 1);

    double * NeVec      = new double[numShift];
    double * NeVec_temp = new double[numShift];
    for(Int l = 0; l < numShift; l++)
    {
      shiftVec[l] = muMin + (l+1) * hsShift;
      NeVec[l] = 0.0;
      NeVec_temp[l] = 0.0;
    }

    bool  isConverged = false;

    // call the PEXSI and get the (H - z*S) inverse to get the NeVec[l]

    // communicator that are needed in the point parallelization
    MPI_Comm pointColComm, pointRowComm;
    MPI_Comm_split( gridPole_->colComm, myPoint, myPointRank, &pointColComm);
    MPI_Comm_split( gridPole_->colComm, myPointRank, myPoint, &pointRowComm);

    GetTime( timePEXSISta );

    for(int l = 0; l < numShift ; l++)
    {
      muPEXSI = shiftVec[l];
      Int currentPoint = l % nPoints;

      // point parallelization.
      if(myPoint == currentPoint){

        CalculateFermiOperatorReal3(
            numPole,
            temperature,
            gap,
            deltaE,
            numElectronExact,
            numElectronPEXSITolerance,  // Not used
            solver,
            verbosity,
            muPEXSI,
            numElectronPEXSI,
            method,
            nPoints,
            spin);

        // exact numElectron for the correspoding muPEXSI
        NeVec_temp[l] = numElectronPEXSI;

      }
    }

    // correct the coressponding EDM
    /*
       CalculateEDMCorrectionReal(
       numPole,
       solver,
       verbosity,
       nPoints,
       spin);
    */

    MPI_Allreduce(NeVec_temp, NeVec, numShift, MPI_DOUBLE, MPI_SUM, pointRowComm);

    GetTime( timePEXSIEnd );
    timePEXSI = timePEXSIEnd - timePEXSISta;

    if( verbosity >= 1 ) {
      for(int i = 0; i < numShift; i++){
        statusOFS << " MuPEXSI from Point :" << i << " is  " << shiftVec[i]<< " numElectron " << NeVec[i] << std::endl;
        if( (i < numShift-1 )&& (NeVec[i] > NeVec[i+1]) ) {
          statusOFS << " number of poles is not enough for the calculation, try increase the number of poles. " << std::endl;
        }
      }
    }

    InterpolateDMReal(
        numElectronExact,
        numElectronPEXSI,
        numElectronPEXSITolerance,
        nPoints,
        NeVec,
        muMinInertia,
        muMaxInertia,
        muPEXSI,
        method,
        verbosity);

    MPI_Comm_free(&pointColComm);
    MPI_Comm_free(&pointRowComm);

    delete []NeVec;
    delete []NeVec_temp;
  }

  GetTime( timeTotalEnd );

  if( verbosity == 1 ){
    statusOFS << std::endl
      << "(Point 0) Total number of inertia counts       = " << numTotalInertiaIter << std::endl
      << "(Point 0) Total time for inertia count step    = " << timeInertia << " [s]" << std::endl
      << "(Point 0) Total time for PEXSI step            = " << timePEXSI   << " [s]" << std::endl
      << "(Point 0) Total time for the DFT driver        = " << timeTotalEnd - timeTotalSta   << " [s]" << std::endl
      << std::endl;
  }

  if( 1 ){
    statusOFS << "Final result " << std::endl;
    Print( statusOFS, "mu                          = ", muPEXSI );
    Print( statusOFS, "muMin                       = ", muMinInertia );
    Print( statusOFS, "muMax                       = ", muMaxInertia );
    Print( statusOFS, "Computed number of electron = ", numElectronPEXSI );
    Print( statusOFS, "Exact number of electron    = ", numElectronExact );
    Print( statusOFS, "Total energy (H*DM)         = ", totalEnergyH_ );
    Print( statusOFS, "Total energy (Tr[EDM])      = ", totalEnergyS_ );
    Print( statusOFS, "Total free energy           = ", totalFreeEnergy_ );
    statusOFS << std::endl << std::endl;
  }

  return ;
}         // -----  end of method PPEXSIData::DFTDriver2  -----


void PPEXSIData::CalculateFermiOperatorReal3(
    Int   numPole,
    Real  temperature,
    Real  gap,
    Real  deltaE,
    Real  numElectronExact,
    Real  numElectronTolerance,
    Int   solver,
    Int   verbosity,
    Real& mu,
    Real& numElectron,
    Int   method,
    Int   nPoints,
    Real  spin) {

  // These are needed in the point-pole-pexsi parallelization
  Int npPerPoint  = gridPole_->mpisize / nPoints;
  Int myPoint     = gridPole_->mpirank / npPerPoint;
  Int myPointRank = gridPole_->mpirank % npPerPoint;
  Int myRowPoint  = myPointRank / gridPole_->numProcCol;

  MPI_Comm pointColComm, pointRowComm;
  MPI_Comm_split( gridPole_->colComm, myPoint, myPointRank, &pointColComm);
  MPI_Comm_split( gridPole_->colComm, myPointRank, myPoint, &pointRowComm);

  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  if( isComplexSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeComplexSymmetricMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // *********************************************************************
  // Check the input parameters
  // *********************************************************************
  if( numPole % 2 != 0 ){
    //ErrorHandling( "Must be even number of poles!" );
    statusOFS << " Warning, numPole is " << numPole <<" not a even number. "<< std::endl << std::flush;
  }

  // *********************************************************************
  // Initialize
  // *********************************************************************
  // Rename for convenience
  DistSparseMatrix<Real>&  HMat        = HRealMat_;
  DistSparseMatrix<Real>&  SMat        = SRealMat_;

  DistSparseMatrix<Real>& rhoMat       = rhoRealMat_;
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

  // only calculate FDM and EDM within the EnergyDensity when method == 3
  Real numSpin = spin;

  // Copy the pattern
  CopyPattern( PatternMat_, AMat );
  CopyPattern( PatternMat_, rhoMat );

  if( isEnergyDensityMatrix )
    CopyPattern( PatternMat_, frcMat );

  if( isFreeEnergyDensityMatrix )
    CopyPattern( PatternMat_, hmzMat );

  // Reinitialize the variables
  SetValue( rhoMat.nzvalLocal, 0.0 );

  if( isEnergyDensityMatrix )
    SetValue( frcMat.nzvalLocal, 0.0 );

  if( isFreeEnergyDensityMatrix )
    SetValue( hmzMat.nzvalLocal, 0.0 );

  // FIXME

  // Directly read off the poles from a table.

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
        fdGrid[i] = numSpin * ez / (1.0 + ez);
      }
      else{
        ez = std::exp((xGrid[i] - mu) / temperature );
        fdGrid[i] = numSpin / (1.0 + ez);
      }
    }

    numPoleInput = numPole;
    Real tol;
    Int  numPoleSignificant;

    Int poleIter = 0;
    /// test the getPole function. actually, NumPole should be get from here.
    {
      // first get the real numPole
      poleClass pole;
      bool result;
      zweightFDM_.resize(numPole);
      zweightEDM_.resize(numPole);

      if( method == 1 || method == 2 || method == 3)
        result = pole.getPole(method, numPole,deltaE/temperature, zshift_, zweightRho_, zweightFDM_, zweightEDM_);
      else{
        statusOFS << " Error, method must be chosen from [ 1, 2, 3] " << std::endl;
        statusOFS << " Error, Current method is: " << method << std::endl;
        ErrorHandling("getPole error.");
      }
      if(method == 2)
        isEDMCorrection_ = 1;
      else
        isEDMCorrection_ = 0;

      for(int i = 0; i < zshift_.size(); i++)
      {
        zshift_[i] *= temperature;
        zweightRho_[i] *= temperature;
        zshift_[i] += mu;
        if( isEDMCorrection_ == 0){
          zweightFDM_[i] *= temperature * temperature;
          zweightEDM_[i] *= temperature * temperature;
        }
      }

      if(verbosity >= 2){
        statusOFS << " numPole is "<< numPole << " deltaE  "
          << deltaE/ temperature << std::endl;
        statusOFS << " temperature " << temperature << std::endl;
        statusOFS << " zshift "<< zshift_ << " zweight "
          << zweightRho_ << std::endl;
      }
      if(!result)
      {
        statusOFS << " Eroror getPole wrong method: " << method
          <<" " <<numPole << " " << deltaE/temperature << std::endl;
        statusOFS << " numPoles should be chosen from [ 10, 15, 20, 25, 30, 35, 40] " << std::endl;
        statusOFS << " Please check the temperature setting, if it is too high or too low, there is also potential problem " << std::endl;
        ErrorHandling("getPole error.");
      }
      if(numPole != zshift_.size()){
        statusOFS << "  Error: numPole should be changed to "
          << zshift_.size()   << std::endl;
        ErrorHandling("getPole error.");
      }
    }

  }

  // Initialize the number of electrons
  numElectron  = 0.0;
  if( verbosity >= 2 ){
    statusOFS << "zshift" << std::endl << zshift_ << std::endl;
    statusOFS << "zweightRho" << std::endl << zweightRho_ << std::endl;
  }

  // for each pole, perform LDLT factoriation and selected inversion
  Real timePoleSta, timePoleEnd;

  Int numPoleComputed = 0;
  for(Int lidx = 0; lidx < numPole; lidx++){
    // add another line here to do the
    if( myRowPoint % numPole == lidx % (gridPole_->numProcRow /nPoints)){
      //Int l = poleIdx[lidx];
      Int l = lidx;

      GetTime( timePoleSta );

      if( verbosity >= 1 ){
        statusOFS << "Pole " << lidx << " processing..." << std::endl;
      }
      if( verbosity >= 1 ){
        statusOFS << "zshift           = " << zshift_[l] << std::endl;
        statusOFS    << "zweightRho       = " << zweightRho_[l] << std::endl;
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
        Real timeTotalSelInvSta, timeTotalSelInvEnd;
        switch (solver) {
          case 0:
            {
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
              GetTime( timeTotalSelInvSta );

              luMat.LUstructToPMatrix( PMloc );
            }
            break;
#ifdef WITH_SYMPACK
          case 1:
            {
              symPACK::symPACKMatrix<Complex>& symPACKMat = *symPACKComplexMat_ ;
              symPACK::DistSparseMatrix<Complex> ltAMat;
              if( verbosity >= 2 ){
                statusOFS << "Before ToLowerTriangular." << std::endl;
              }
              Convert(AMat,ltAMat);
              ltAMat.ToLowerTriangular();
              if( verbosity >= 2 ){
                statusOFS << "After ToLowerTriangular." << std::endl;
              }


              Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

              GetTime( timeTotalFactorizationSta );

              // Data redistribution
              if( verbosity >= 2 ){
                statusOFS << "Before Distribute." << std::endl;
              }
              symPACKMat.DistributeMatrix(ltAMat);
              if( verbosity >= 2 ){
                statusOFS << "After Distribute." << std::endl;
              }

              // Numerical factorization
              if( verbosity >= 2 ){
                statusOFS << "Before NumericalFactorize." << std::endl;
              }
              symPACKMat.Factorize();
              // Numerical factorization
              if( verbosity >= 2 ){
                statusOFS << "After NumericalFactorize." << std::endl;
              }

              GetTime( timeTotalFactorizationEnd );

              if( verbosity >= 1 ){
                statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl;
              }

              //make a symmetric matrix out of that....
              GetTime( timeTotalSelInvSta );
              symPACKMatrixToPMatrix( symPACKMat, PMloc );
            }
            break;
#endif
          default:
            ErrorHandling("Unsupported solver.");
            break;
        }

        PMloc.PreSelInv();

        // Main subroutine for selected inversion
        PMloc.SelInv();

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

        //TODO convert to symmAinvMat too
        PMloc.PMatrixToDistSparseMatrix( PatternMat_, AinvMat );

        if( verbosity >= 2 ){
          statusOFS << "rhoMat.nnzLocal = " << rhoMat.nnzLocal << std::endl;
          statusOFS << "AinvMat.nnzLocal = " << AinvMat.nnzLocal << std::endl;
        }


        // Update the density matrix. The following lines are equivalent to
        //
        //                for( Int i = 0; i < rhoMat.nnzLocal; i++ ){
        //                    rhoMat.nzvalLocal(i) +=
        //                        zweightRho_[l].real() * AinvMat.nzvalLocal(i).imag() +
        //                        zweightRho_[l].imag() * AinvMat.nzvalLocal(i).real();
        //                }
        //
        // But done more cache-efficiently with blas.
        Real* AinvMatRealPtr = (Real*)AinvMat.nzvalLocal.Data();
        Real* AinvMatImagPtr = AinvMatRealPtr + 1;
        blas::Axpy( rhoMat.nnzLocal, numSpin*zweightRho_[l].real(), AinvMatImagPtr, 2,
            rhoMat.nzvalLocal.Data(), 1 );
        blas::Axpy( rhoMat.nnzLocal, numSpin*zweightRho_[l].imag(), AinvMatRealPtr, 2,
            rhoMat.nzvalLocal.Data(), 1 );

        // Free energy density matrix
        if( isFreeEnergyDensityMatrix ){
          if( isEDMCorrection_ == 0){
            blas::Axpy( hmzMat.nnzLocal, numSpin*zweightFDM_[l].real(), AinvMatImagPtr, 2,
                hmzMat.nzvalLocal.Data(), 1 );
            blas::Axpy( hmzMat.nnzLocal, numSpin*zweightFDM_[l].imag(), AinvMatRealPtr, 2,
                hmzMat.nzvalLocal.Data(), 1 );
          }
        }


        // Energy density matrix
        if( isEnergyDensityMatrix ){
          if( isEDMCorrection_ == 0){
            blas::Axpy( frcMat.nnzLocal, numSpin*zweightEDM_[l].real(), AinvMatImagPtr, 2,
                frcMat.nzvalLocal.Data(), 1 );
            blas::Axpy( frcMat.nnzLocal, numSpin*zweightEDM_[l].imag(), AinvMatRealPtr, 2,
                frcMat.nzvalLocal.Data(), 1 );
          }
          else{
            Real wzReal = numSpin*(zweightRho_[l].real() * zshift_[l].real()
                - zweightRho_[l].imag() * zshift_[l].imag());
            Real wzImag = numSpin*(zweightRho_[l].real() * zshift_[l].imag()
                + zweightRho_[l].imag() * zshift_[l].real() );
            blas::Axpy( frcMat.nnzLocal, wzReal, AinvMatImagPtr, 2,
                frcMat.nzvalLocal.Data(), 1 );
            blas::Axpy( frcMat.nnzLocal, wzImag, AinvMatRealPtr, 2,
                frcMat.nzvalLocal.Data(), 1 );
          }
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
        rhoMat.nnzLocal, MPI_SUM, pointColComm);

  }

  // Reduce the free energy density matrix across the processor rows in gridPole_
  if( isFreeEnergyDensityMatrix ){
    DblNumVec nzvalHmzMatLocal = hmzMat.nzvalLocal;
    SetValue( hmzMat.nzvalLocal, 0.0 );

    mpi::Allreduce( nzvalHmzMatLocal.Data(), hmzMat.nzvalLocal.Data(),
        hmzMat.nnzLocal, MPI_SUM, pointColComm);

  }


  if( isEnergyDensityMatrix ){
    DblNumVec nzvalFrcMatLocal = frcMat.nzvalLocal;
    SetValue( frcMat.nzvalLocal, 0.0 );

    mpi::Allreduce( nzvalFrcMatLocal.Data(), frcMat.nzvalLocal.Data(),
         frcMat.nnzLocal, MPI_SUM, pointColComm);

    if( isEDMCorrection_ == 0){
      // FIXME
      // EDM correction: add mu * DM for method == 3/1
      blas::Axpy( rhoMat.nnzLocal, mu, rhoMat.nzvalLocal.Data(), 1,
        frcMat.nzvalLocal.Data(), 1 );
    }
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

#if ( _DEBUGlevel_ >= 0 )
    statusOFS << std::endl << "numElecLocal = " << numElecLocal << std::endl;
#endif
    mpi::Allreduce( &numElecLocal, &numElectron, 1, MPI_SUM, gridPole_->rowComm);
  }

  if( isEDMCorrection_ == 0)
  {
    // Energy computed from Tr[S*EDM]
    if( isEnergyDensityMatrix )
    {
      Real local = 0.0;
      if( SMat.size != 0 ){
        local = blas::Dot( SMat.nnzLocal,
            SMat.nzvalLocal.Data(),
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
      //std::cout << " Energy Density Tr[S*FDM] " << totalEnergyS_ << std::endl;
    }


    // Free energy
    if( isFreeEnergyDensityMatrix )
    {
      Real local = 0.0;
      if( SMat.size != 0 ){
        local = blas::Dot( SMat.nnzLocal,
            SMat.nzvalLocal.Data(),
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
      totalFreeEnergy_ += mu * numElectronExact;

      //std::cout << " Free  Energy S Tr[S*FDM] " << totalFreeEnergy_ << std::endl;
    }
  }

  MPI_Comm_free(&pointColComm);
  MPI_Comm_free(&pointRowComm);

  return ;
}    // -----  end of method PPEXSIData::CalculateFermiOperatorReal3  -----

void PPEXSIData::CalculateEDMCorrectionReal(
    Int   numPole,
    Int   solver,
    Int   verbosity,
    Int   nPoints,
    Real  spin) {

  // add the points parallelization.
  /*
     Int npPerPoint  = gridPole_->mpisize / nPoints;
     Int myPoint     = gridPole_->mpirank / npPerPoint;
     Int myPointRank = gridPole_->mpirank % npPerPoint;
     Int myRowPoint  = myPointRank / gridPole_->numProcCol;
     MPI_Comm pointColComm, pointRowComm;
     MPI_Comm_split( gridPole_->colComm, myPoint, myPointRank, &pointColComm);
     MPI_Comm_split( gridPole_->colComm, myPointRank, myPoint, &pointRowComm);
   */

  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadRealMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  if( isComplexSymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeComplexSymmetricMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // *********************************************************************
  // Initialize
  // *********************************************************************
  // Rename for convenience
  DistSparseMatrix<Real>&  SMat        = SRealMat_;
  DistSparseMatrix<Real>& frcMat       = energyDensityRealMat_;
  DistSparseMatrix<Complex>& AMat      = shiftComplexMat_;
  DistSparseMatrix<Complex>& AinvMat   = shiftInvComplexMat_;

  // The symbolic information should have already been there.
  SuperLUMatrix<Complex>& luMat        = *luComplexMat_;
  PMatrix<Complex>&       PMloc        = *PMComplexMat_;

  bool isFreeEnergyDensityMatrix = true;
  bool isEnergyDensityMatrix     = true;
  bool isDerivativeTMatrix       = false;
  Real numSpin = spin;

  // SMat.size == 0, do not need to invert the S matrix.
  // Just accumulate the correction term and return
  if( SMat.size == 0 ) {
    if( isEnergyDensityMatrix ){
      for(Int l = 0; l < numPole; l++){
        for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
          frcMat.nzvalLocal( diagIdxLocal_[i] ) += numSpin*zweightRho_[l].imag();
        }
      }
    }

    Real local = 0.0;

    DblNumVec& nzval = energyDensityRealMat_.nzvalLocal;
    for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
      local += nzval(diagIdxLocal_[i]);
    }
    mpi::Allreduce( &local, &totalEnergyS_, 1, MPI_SUM,
        gridPole_->rowComm );
    return;
  }

  // If S is not Identity matrix
  // Copy the pattern
  CopyPattern( PatternMat_, AMat );

  Int numPoleInput = numPole;
  Real tol;
  Int  numPoleSignificant;
  Int poleIter = 0;

  // for each pole, perform LDLT factoriation and selected inversion
  Real timePoleSta, timePoleEnd;
  if( SMat.size != 0 ){
    // S is not an identity matrix
    for( Int i = 0; i < SMat.nnzLocal; i++ ){
      //AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift_[l] * SMat.nzvalLocal(i);
      AMat.nzvalLocal(i) =  SMat.nzvalLocal(i);
    }

    Real timeTotalSelInvSta, timeTotalSelInvEnd;


    switch (solver) {
      case 0:
        {


          if( verbosity >= 2 ){
            statusOFS << "Sinv Before DistSparseMatrixToSuperMatrixNRloc." << std::endl;
          }
          luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
          if( verbosity >= 2 ){
            statusOFS << "Sinv After DistSparseMatrixToSuperMatrixNRloc." << std::endl;
          }

          Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

          GetTime( timeTotalFactorizationSta );

          // Data redistribution
          if( verbosity >= 2 ){
            statusOFS << "Sinv Before Distribute." << std::endl;
          }
          luMat.Distribute();
          if( verbosity >= 2 ){
            statusOFS << "Sinv After Distribute." << std::endl;
          }

          // Numerical factorization
          if( verbosity >= 2 ){
            statusOFS << "Sinv Before NumericalFactorize." << std::endl;
          }
          luMat.NumericalFactorize();
          if( verbosity >= 2 ){
            statusOFS << "Sinv After NumericalFactorize." << std::endl;
          }
          luMat.DestroyAOnly();

          GetTime( timeTotalFactorizationEnd );

          if( verbosity >= 1 ){
            statusOFS << "Sinv Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl;
          }

          // *********************************************************************
          // Selected inversion
          // *********************************************************************
          GetTime( timeTotalSelInvSta );

          luMat.LUstructToPMatrix( PMloc );


        }
        break;
#ifdef WITH_SYMPACK
      case 1:
        {
          symPACK::symPACKMatrix<Complex>& symPACKMat = *symPACKComplexMat_ ;
          symPACK::DistSparseMatrix<Complex> ltAMat;
          if( verbosity >= 2 ){
            statusOFS << "Before ToLowerTriangular." << std::endl;
          }
          Convert(AMat,ltAMat);
          ltAMat.ToLowerTriangular();
          ltAMat.GetLocalGraph().SetSorted(false);
          ltAMat.SortGraph();
          if( verbosity >= 2 ){
            statusOFS << "After ToLowerTriangular." << std::endl;
          }


          Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

          GetTime( timeTotalFactorizationSta );

          // Data redistribution
          if( verbosity >= 2 ){
            statusOFS << "Before Distribute." << std::endl;
          }
          symPACKMat.DistributeMatrix(ltAMat);
          if( verbosity >= 2 ){
            statusOFS << "After Distribute." << std::endl;
          }

          // Numerical factorization
          if( verbosity >= 2 ){
            statusOFS << "Before NumericalFactorize." << std::endl;
          }
          symPACKMat.Factorize();
          // Numerical factorization
          if( verbosity >= 2 ){
            statusOFS << "After NumericalFactorize." << std::endl;
          }

          GetTime( timeTotalFactorizationEnd );

          if( verbosity >= 1 ){
            statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl;
          }

          //make a symmetric matrix out of that....
          GetTime( timeTotalSelInvSta );
          symPACKMatrixToPMatrix( symPACKMat, PMloc );
        }
        break;
#endif
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }

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
      statusOFS << "Sinv Time for total selected inversion is " <<
        timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
    }

    // *********************************************************************
    // Postprocessing
    // *********************************************************************

    //TODO convert to symmAinvMat too
    PMloc.PMatrixToDistSparseMatrix( PatternMat_, AinvMat );

    if( verbosity >= 1 ){
      statusOFS << "frcMat.nnzLocal = " << frcMat.nnzLocal << std::endl;
      statusOFS << "AinvMat.nnzLocal = " << AinvMat.nnzLocal << std::endl;
    }

  }

  // Accumulate the pre-factor
  // Correction term: \Omega_l * S^{-1}
  Real ImgW = 0.0;
  for(Int l = 0; l < numPole; l++){
    ImgW += numSpin *  zweightRho_[l].imag();
  }

  // Energy density matrix
  if( isEnergyDensityMatrix ){
    for( Int i = 0; i < SMat.nnzLocal; i++ ){
      frcMat.nzvalLocal(i) +=   ImgW * AinvMat.nzvalLocal(i).real();
    }
  }

  // Reduce the energy density matrix across the processor rows in gridPole_
  {
    Real local = 0.0;
    local = blas::Dot( SMat.nnzLocal,
        SMat.nzvalLocal.Data(),
        1, energyDensityRealMat_.nzvalLocal.Data(), 1 );
    mpi::Allreduce( &local, &totalEnergyS_, 1, MPI_SUM,
        gridPole_->rowComm );
  }

  return ;

}  // ---  PPEXSIData::CalculateEDMCorrectionReal  -----

void PPEXSIData::CalculateEDMCorrectionComplex(
    Int   numPole,
    Int   solver,
    Int   verbosity,
    Int   nPoints,
    Real  spin) {

  if(solver==1){
    std::ostringstream msg;
    msg  << std::endl
      << "Solver not supported in this routine." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }


  // add the points parallelization.
  /*
     Int npPerPoint  = gridPole_->mpisize / nPoints;
     Int myPoint     = gridPole_->mpirank / npPerPoint;
     Int myPointRank = gridPole_->mpirank % npPerPoint;
     Int myRowPoint  = myPointRank / gridPole_->numProcCol;

     MPI_Comm pointColComm, pointRowComm;
     MPI_Comm_split( gridPole_->colComm, myPoint, myPointRank, &pointColComm);
     MPI_Comm_split( gridPole_->colComm, myPointRank, myPoint, &pointRowComm);
   */
  if( isMatrixLoaded_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been loaded." << std::endl
      << "Call LoadComplexMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  if( isComplexUnsymmetricSymbolicFactorized_ == false ){
    std::ostringstream msg;
    msg  << std::endl
      << "Matrix has not been factorized symbolically." << std::endl
      << "Call SymbolicFactorizeComplexUnsymmetricMatrix first." << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  // *********************************************************************
  // Initialize
  // *********************************************************************
  // Rename for convenience
  DistSparseMatrix<Complex>&  SMat        = SComplexMat_;
  DistSparseMatrix<Complex>& frcMat       = energyDensityComplexMat_;
  DistSparseMatrix<Complex>& AMat      = shiftComplexMat_;
  DistSparseMatrix<Complex>& AinvMat   = shiftInvComplexMat_;
  // The symbolic information should have already been there.
  SuperLUMatrix<Complex>& luMat        = *luComplexMat_;
  PMatrixUnsym<Complex>&  PMloc        = *PMComplexUnsymMat_;

  // calculate EDM
  bool isFreeEnergyDensityMatrix = true;
  bool isEnergyDensityMatrix     = true;
  bool isDerivativeTMatrix       = false;
  Real numSpin = spin;

  // SMat.size == 0, do not need to invert the S matrix.
  // Just accumulate the correction term and return
  if( SMat.size == 0 ) {
    if( isEnergyDensityMatrix ) {
      for(Int l = 0; l < numPole; l++){
        for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
          frcMat.nzvalLocal( diagIdxLocal_[i] ) += numSpin*zweightRho_[l].imag();
        }
      }
    }
    Real local = 0.0;
    CpxNumVec& nzval = energyDensityComplexMat_.nzvalLocal;
    for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
      local += nzval(diagIdxLocal_[i]).real();
    }
    mpi::Allreduce( &local, &totalEnergyS_, 1, MPI_SUM,
        gridPole_->rowComm );
    statusOFS << " recalculate the totalEnergyS_ " << totalEnergyS_ << std::endl;
    return;
  }

  // Smat.size != 0, need to invert the S matrix
  // Copy the pattern
  CopyPattern( PatternMat_, AMat );

  // for each pole, perform LDLT factoriation and selected inversion
  Real timePoleSta, timePoleEnd;
  if( SMat.size != 0 ){
    // S is not an identity matrix
    for( Int i = 0; i < SMat.nnzLocal; i++ ){
      //AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift_[l] * SMat.nzvalLocal(i);
      AMat.nzvalLocal(i) =  SMat.nzvalLocal(i);
    }

    Real timeTotalSelInvSta, timeTotalSelInvEnd;
    switch (solver) {
      case 0:
        {
          if( verbosity >= 2 ){
            statusOFS << "Sinv Before DistSparseMatrixToSuperMatrixNRloc." << std::endl;
          }
          luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt_ );
          if( verbosity >= 2 ){
            statusOFS << "Sinv After DistSparseMatrixToSuperMatrixNRloc." << std::endl;
          }

          Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

          GetTime( timeTotalFactorizationSta );

          // Data redistribution
          if( verbosity >= 2 ){
            statusOFS << "Sinv Before Distribute." << std::endl;
          }
          luMat.Distribute();
          if( verbosity >= 2 ){
            statusOFS << "Sinv After Distribute." << std::endl;
          }

          // Numerical factorization
          if( verbosity >= 2 ){
            statusOFS << "Sinv Before NumericalFactorize." << std::endl;
          }
          luMat.NumericalFactorize();
          if( verbosity >= 2 ){
            statusOFS << "Sinv After NumericalFactorize." << std::endl;
          }
          luMat.DestroyAOnly();

          GetTime( timeTotalFactorizationEnd );

          if( verbosity >= 1 ){
            statusOFS << "Sinv Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl;
          }

          // *********************************************************************
          // Selected inversion
          // *********************************************************************
          GetTime( timeTotalSelInvSta );

          luMat.LUstructToPMatrix( PMloc );
        }
        break;
      default:
        ErrorHandling("Unsupported solver.");
        break;
    }



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
      statusOFS << "Sinv Time for total selected inversion is " <<
        timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
    }

    // *********************************************************************
    // Postprocessing
    // *********************************************************************

    //TODO convert to symmAinvMat too
    PMloc.PMatrixToDistSparseMatrix( PatternMat_, AinvMat );

    if( verbosity >= 1 ){
      statusOFS << "frcMat.nnzLocal = " << frcMat.nnzLocal << std::endl;
      statusOFS << "AinvMat.nnzLocal = " << AinvMat.nnzLocal << std::endl;
    }

  }

  DistSparseMatrix<Complex>& regMat       = AinvMat;
  DistSparseMatrix<Complex>  transMat;

  // After selected inversion,  AinvMat stores the selected elements of
  // S^{-T}, so transMat stores S^{-1}.
  CSCToCSR( regMat, transMat );

  // Accumulate the pre-factor
  // Correction term: \Omega_l * S^{-1}
  Real ImgW = 0.0;
  for(Int l = 0; l < numPole; l++){
    ImgW += numSpin *  zweightRho_[l].imag();
  }

  // Energy density matrix
  if( isEnergyDensityMatrix ){
    for( Int i = 0; i < SMat.nnzLocal; i++ ){
      frcMat.nzvalLocal(i) +=   ImgW * transMat.nzvalLocal(i);
    }
  }


  {
    Real local = 0.0;
    local = (blas::Dotc( SMat.nnzLocal,
          SMat.nzvalLocal.Data(),
          1, energyDensityComplexMat_.nzvalLocal.Data(), 1 )).real();
    mpi::Allreduce( &local, &totalEnergyS_, 1, MPI_SUM,
        gridPole_->rowComm );
  }
  statusOFS << " recalculate the totalEnergyS_ " << totalEnergyS_ << std::endl;

  return ;

}  // ---  PPEXSIData::CalculateEDMCorrectionComplex  ---

void PPEXSIData::InterpolateDMReal(
    Real                numElectronExact,
    Real                numElectronPEXSI,
    Real                numElectronPEXSITolerance,
    Int                 nPoints,
    Real              * NeVec,
    Real              & muMin,
    Real              & muMax,
    Real              & muPEXSI,
    Int                 method,
    Int                 verbosity){

  DistSparseMatrix<Real>&  SMat        = SRealMat_;

  // get the MPI communicators.
  Int npPerPoint  = gridPole_->mpisize / nPoints;
  if(npPerPoint < 1) {
    statusOFS << " Error, Point parallelization error, npPerPoint < 1 " << std::endl;
  }
  Int myPoint     = gridPole_->mpirank / npPerPoint;
  Int myPointRank = gridPole_->mpirank % npPerPoint;
  Int myRowPoint  = myPointRank / gridPole_->numProcCol;
  MPI_Comm pointColComm, pointRowComm;
  MPI_Comm_split( gridPole_->colComm, myPoint, myPointRank, &pointColComm);
  MPI_Comm_split( gridPole_->colComm, myPointRank, myPoint, &pointRowComm);

  Int l, idxMin, idxMax;
  idxMin = -1;
  idxMax = nPoints;
  Int numShift = nPoints;

  std::vector<Real>  shiftVec( numShift );
  Real hsShift = ( muMax - muMin ) / ( numShift + 1);
  for(Int l = 0; l < numShift; l++)
    shiftVec[l] = muMin + (l+1) * hsShift;

  for(l = 0; l < numShift; l++)
    if( NeVec[l] < numElectronExact - numElectronPEXSITolerance)
    {
      muMin = shiftVec[l];
      idxMin = l;
    }

  for(l = numShift-1; l >= 0; l--)
    if( NeVec[l] > numElectronExact + numElectronPEXSITolerance)
    {
      muMax = shiftVec[l];
      idxMax = l;
    }
  Real mu;
  if(numShift == 1) {
    mu = shiftVec[0];

    // Scale the density matrix. numElectron = numElectronExact in the linear sense
    Real fac = numElectronExact / numElectronPEXSI;
    blas::Scal( rhoRealMat_.nnzLocal, fac, rhoRealMat_.nzvalLocal.Data(), 1 );
    blas::Scal( energyDensityRealMat_.nnzLocal, fac,
        energyDensityRealMat_.nzvalLocal.Data(), 1 );
    numElectronPEXSI = numElectronExact;

    // calculate the total Energy  - check check
    {
      Real local = 0.0;
      local = blas::Dot( HRealMat_.nnzLocal,
          HRealMat_.nzvalLocal.Data(),
          1, rhoRealMat_.nzvalLocal.Data(), 1 );

      mpi::Allreduce( &local, &totalEnergyH_, 1, MPI_SUM,
          gridPole_->rowComm);
    }

    // Energy computed from Tr[S*EDM] check check
    // if( isEnergyDensityMatrix )
    {
      Real local = 0.0;

      if( SMat.size != 0 ){
        local = blas::Dot( SMat.nnzLocal,
            SMat.nzvalLocal.Data(),
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
  }
  else
  {
    if(idxMin == -1) {
      idxMin = 0;
      if(idxMax <= idxMin)
        idxMax = 1;
    }
    if(idxMax == numShift) {
      idxMax = numShift -1;
      if(idxMin >= idxMax)
        idxMin = numShift - 2;
    }

    // check if converged.
    bool isLinearInterplate = true;
    for(Int mu_idx = idxMin; mu_idx <= idxMax; mu_idx++)
    {
      if (std::abs( NeVec[mu_idx] - numElectronExact ) < numElectronPEXSITolerance){
        mu = shiftVec[mu_idx];
        numElectronPEXSI = NeVec[mu_idx];
        isLinearInterplate = false;
        MPI_Bcast(rhoRealMat_.nzvalLocal.Data(), rhoRealMat_.nnzLocal, MPI_DOUBLE, mu_idx, pointRowComm);
        MPI_Bcast(energyDensityRealMat_.nzvalLocal.Data(),
            energyDensityRealMat_.nnzLocal, MPI_DOUBLE, mu_idx, pointRowComm);

        if( verbosity >= 1 ) {
          statusOFS << "PEXSI Converged " <<std::endl;
          Print( statusOFS, "converged idx               = ", mu_idx);
          Print( statusOFS, "converged NeVec[mu_idx]     = ", NeVec[mu_idx]);
          Print( statusOFS, "numElectronTolerance        = ", numElectronPEXSITolerance);
          Print( statusOFS, "Final mu                    = ", mu);
        }
        break;
      }
    }

    if(isLinearInterplate) {
      // Linear interpolation. when not converged.
      numElectronPEXSI = numElectronExact;
      mu = shiftVec[idxMin] + (numElectronExact - NeVec[idxMin])/
        ( NeVec[idxMax] -NeVec[idxMin] ) * (shiftVec[idxMax] - shiftVec[idxMin]);

      if( verbosity >= 1 ) {
        statusOFS << " idxMin = " << idxMin << " idxMax = " << idxMax <<std::endl;
        Print( statusOFS, "idxMin                      = ", idxMin );
        Print( statusOFS, "NeVec[min]                  = ", NeVec[idxMin] );
        Print( statusOFS, "idxMax                      = ", idxMax );
        Print( statusOFS, "NeVec[max]                  = ", NeVec[idxMax] );
        Print( statusOFS, "Final mu                    = ", mu );
      }

      // linear combine density matrix.
      Real facMin = (NeVec[idxMax]-numElectronExact) / (NeVec[idxMax] - NeVec[idxMin]);
      Real facMax = (numElectronExact-NeVec[idxMin]) / (NeVec[idxMax] - NeVec[idxMin]);

      DistSparseMatrix<Real> rhoMatMin ;
      CopyPattern( rhoRealMat_, rhoMatMin );
      rhoMatMin.nzvalLocal = rhoRealMat_.nzvalLocal;
      MPI_Bcast(rhoMatMin.nzvalLocal.Data(), rhoMatMin.nnzLocal, MPI_DOUBLE, idxMin, pointRowComm);

      DistSparseMatrix<Real> rhoMatMax ;
      CopyPattern( rhoRealMat_, rhoMatMax );
      rhoMatMax.nzvalLocal = rhoRealMat_.nzvalLocal;
      MPI_Bcast(rhoMatMax.nzvalLocal.Data(), rhoMatMax.nnzLocal, MPI_DOUBLE, idxMax, pointRowComm);

      SetValue( rhoRealMat_.nzvalLocal, 0.0 );
      blas::Axpy( rhoRealMat_.nnzLocal, facMin, rhoMatMin.nzvalLocal.Data(), 1,
          rhoRealMat_.nzvalLocal.Data(), 1 );

      blas::Axpy( rhoRealMat_.nnzLocal, facMax, rhoMatMax.nzvalLocal.Data(), 1,
          rhoRealMat_.nzvalLocal.Data(), 1 );
      CopyPattern( energyDensityRealMat_, rhoMatMin );
      rhoMatMin.nzvalLocal = energyDensityRealMat_.nzvalLocal;
      MPI_Bcast(rhoMatMin.nzvalLocal.Data(), rhoMatMin.nnzLocal, MPI_DOUBLE, idxMin, pointRowComm);

      CopyPattern( energyDensityRealMat_, rhoMatMax );
      rhoMatMax.nzvalLocal = energyDensityRealMat_.nzvalLocal;
      MPI_Bcast(rhoMatMax.nzvalLocal.Data(), rhoMatMax.nnzLocal, MPI_DOUBLE, idxMax, pointRowComm);

      SetValue( energyDensityRealMat_.nzvalLocal, 0.0 );
      blas::Axpy( energyDensityRealMat_.nnzLocal, facMin, rhoMatMin.nzvalLocal.Data(), 1,
          energyDensityRealMat_.nzvalLocal.Data(), 1 );

      blas::Axpy( energyDensityRealMat_.nnzLocal, facMax, rhoMatMax.nzvalLocal.Data(), 1,
          energyDensityRealMat_.nzvalLocal.Data(), 1 );
    }
    // calculate the totalEnergyH_
    {
      Real local = 0.0;
      local = blas::Dot( HRealMat_.nnzLocal,
          HRealMat_.nzvalLocal.Data(),
          1, rhoRealMat_.nzvalLocal.Data(), 1 );
      mpi::Allreduce( &local, &totalEnergyH_, 1, MPI_SUM,
          gridPole_->rowComm );
    }

    //if( isEnergyDensityMatrix )
    if( isEDMCorrection_ )
    {

      Real local = 0.0;
      if( SMat.size != 0 ){
        local = blas::Dot( SMat.nnzLocal,
            SMat.nzvalLocal.Data(),
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

  }

  // FIXME A placeholder for the FDM -- check check
  // if( isFreeEnergyDensityMatrix )
  if( isEDMCorrection_ )
  {
    DistSparseMatrix<Real>& hmzMat  = freeEnergyDensityRealMat_;
    totalFreeEnergy_                = totalEnergyH_;
    blas::Copy( hmzMat.nnzLocal, rhoRealMat_.nzvalLocal.Data(), 1,
        hmzMat.nzvalLocal.Data(), 1 );
  }

  muPEXSI = mu;

  for(l = 0; l < numShift; l++)
    if( NeVec[l] < numElectronExact - numElectronPEXSITolerance)
      muMin = shiftVec[l];

  for(l = numShift-1; l >= 0; l--)
    if( NeVec[l] > numElectronExact + numElectronPEXSITolerance)
      muMax = shiftVec[l];

  if ( muMax < muMin )
    ErrorHandling( " muMax < muMin , Error occured in the PEXSI ");

  MPI_Comm_free(&pointColComm);
  MPI_Comm_free(&pointRowComm);

} // --- subroutine  interpolateDMReal  ---

void PPEXSIData::InterpolateDMComplex(
    Real                numElectronExact,
    Real                numElectronPEXSI,
    Real                numElectronPEXSITolerance,
    Int                 nPoints,
    Real              * NeVec,
    Real              & muMin,
    Real              & muMax,
    Real              & muPEXSI,
    Int                 method,
    Int                 verbosity){

  DistSparseMatrix<Complex>&  SMat        = SComplexMat_;

  // get the MPI communicators.
  Int npPerPoint  = gridPole_->mpisize / nPoints;
  if(npPerPoint < 1) {
    statusOFS << " Error, Point parallelization error, npPerPoint < 1 " << std::endl;
  }
  Int myPoint     = gridPole_->mpirank / npPerPoint;
  Int myPointRank = gridPole_->mpirank % npPerPoint;
  Int myRowPoint  = myPointRank / gridPole_->numProcCol;
  MPI_Comm pointColComm, pointRowComm;
  MPI_Comm_split( gridPole_->colComm, myPoint, myPointRank, &pointColComm);
  MPI_Comm_split( gridPole_->colComm, myPointRank, myPoint, &pointRowComm);

  Int l, idxMin, idxMax;
  idxMin = -1;
  idxMax = nPoints;
  Int numShift = nPoints;

  std::vector<Real>  shiftVec( numShift );
  Real hsShift = ( muMax - muMin ) / ( numShift + 1);
  for(Int l = 0; l < numShift; l++)
    shiftVec[l] = muMin + (l+1) * hsShift;

  for(l = 0; l < numShift; l++)
    if( NeVec[l] < numElectronExact - numElectronPEXSITolerance)
    {
      muMin = shiftVec[l];
      idxMin = l;
    }

  for(l = numShift-1; l >= 0; l--)
    if( NeVec[l] > numElectronExact + numElectronPEXSITolerance)
    {
      muMax = shiftVec[l];
      idxMax = l;
    }
  Real mu;
  if(numShift == 1) {
    mu = shiftVec[0];

    // Scale the density matrix. numElectron = numElectronExact in the linear sense
    Real fac = numElectronExact / numElectronPEXSI;
    blas::Scal( rhoComplexMat_.nnzLocal, fac, rhoComplexMat_.nzvalLocal.Data(), 1 );
    blas::Scal( energyDensityComplexMat_.nnzLocal, fac,
        energyDensityComplexMat_.nzvalLocal.Data(), 1 );
    numElectronPEXSI = numElectronExact;

    // calculate the total Energy  - check check
    {
      Real local = 0.0;
      local = (blas::Dotc( HComplexMat_.nnzLocal,
            HComplexMat_.nzvalLocal.Data(),
            1, rhoComplexMat_.nzvalLocal.Data(), 1 )).real();

      mpi::Allreduce( &local, &totalEnergyH_, 1, MPI_SUM,
          gridPole_->rowComm);
    }

    // Energy computed from Tr[S*EDM] check check
    // if( isEnergyDensityMatrix )
    {
      Real local = 0.0;

      if( SMat.size != 0 ){
        local = (blas::Dotc( SMat.nnzLocal,
              SMat.nzvalLocal.Data(),
              1, energyDensityComplexMat_.nzvalLocal.Data(), 1 )).real();
      }
      else{
        CpxNumVec& nzval = energyDensityComplexMat_.nzvalLocal;
        for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
          local += nzval(diagIdxLocal_[i]).real();
        }
      }

      mpi::Allreduce( &local, &totalEnergyS_, 1, MPI_SUM,
          gridPole_->rowComm );
    }
  }
  else
  {
    if(idxMin == -1) {
      idxMin = 0;
      if(idxMax <= idxMin)
        idxMax = 1;
    }
    if(idxMax == numShift) {
      idxMax = numShift -1;
      if(idxMin >= idxMax)
        idxMin = numShift - 2;
    }

    // check if converged.
    bool isLinearInterplate = true;
    for(Int mu_idx = idxMin; mu_idx <= idxMax; mu_idx++)
    {
      if (std::abs( NeVec[mu_idx] - numElectronExact ) < numElectronPEXSITolerance){
        mu = shiftVec[mu_idx];
        numElectronPEXSI = NeVec[mu_idx];
        isLinearInterplate = false;
        MPI_Bcast(rhoComplexMat_.nzvalLocal.Data(), rhoComplexMat_.nnzLocal, MPI_DOUBLE_COMPLEX, mu_idx, pointRowComm);
        MPI_Bcast(energyDensityComplexMat_.nzvalLocal.Data(),
            energyDensityComplexMat_.nnzLocal, MPI_DOUBLE_COMPLEX, mu_idx, pointRowComm);

        if( verbosity >= 1 ) {
          statusOFS << "PEXSI Converged " <<std::endl;
          Print( statusOFS, "converged idx               = ", mu_idx);
          Print( statusOFS, "converged NeVec[mu_idx]     = ", NeVec[mu_idx]);
          Print( statusOFS, "numElectronTolerance        = ", numElectronPEXSITolerance);
          Print( statusOFS, "Final mu                    = ", mu);
        }
        break;
      }
    }

    if(isLinearInterplate) {
      // Linear interpolation. when not converged.
      numElectronPEXSI = numElectronExact;
      mu = shiftVec[idxMin] + (numElectronExact - NeVec[idxMin])/
        ( NeVec[idxMax] -NeVec[idxMin] ) * (shiftVec[idxMax] - shiftVec[idxMin]);

      if( verbosity >= 1 ) {
        statusOFS << " idxMin = " << idxMin << " idxMax = " << idxMax <<std::endl;
        Print( statusOFS, "idxMin                      = ", idxMin );
        Print( statusOFS, "NeVec[min]                  = ", NeVec[idxMin] );
        Print( statusOFS, "idxMax                      = ", idxMax );
        Print( statusOFS, "NeVec[max]                  = ", NeVec[idxMax] );
        Print( statusOFS, "Final mu                    = ", mu );
      }

      // linear combine density matrix.
      Real facMin = (NeVec[idxMax]-numElectronExact) / (NeVec[idxMax] - NeVec[idxMin]);
      Real facMax = (numElectronExact-NeVec[idxMin]) / (NeVec[idxMax] - NeVec[idxMin]);

      DistSparseMatrix<Complex> rhoMatMin ;
      CopyPattern( rhoComplexMat_, rhoMatMin );
      rhoMatMin.nzvalLocal = rhoComplexMat_.nzvalLocal;
      MPI_Bcast(rhoMatMin.nzvalLocal.Data(), rhoMatMin.nnzLocal, MPI_DOUBLE, idxMin, pointRowComm);

      DistSparseMatrix<Complex> rhoMatMax ;
      CopyPattern( rhoComplexMat_, rhoMatMax );
      rhoMatMax.nzvalLocal = rhoComplexMat_.nzvalLocal;
      MPI_Bcast(rhoMatMax.nzvalLocal.Data(), rhoMatMax.nnzLocal, MPI_DOUBLE, idxMax, pointRowComm);

      SetValue( rhoComplexMat_.nzvalLocal, Z_ZERO );
      blas::Axpy( rhoComplexMat_.nnzLocal, facMin, rhoMatMin.nzvalLocal.Data(), 1,
          rhoComplexMat_.nzvalLocal.Data(), 1 );

      blas::Axpy( rhoComplexMat_.nnzLocal, facMax, rhoMatMax.nzvalLocal.Data(), 1,
          rhoComplexMat_.nzvalLocal.Data(), 1 );
      CopyPattern( energyDensityComplexMat_, rhoMatMin );
      rhoMatMin.nzvalLocal = energyDensityComplexMat_.nzvalLocal;
      MPI_Bcast(rhoMatMin.nzvalLocal.Data(), rhoMatMin.nnzLocal, MPI_DOUBLE, idxMin, pointRowComm);

      CopyPattern( energyDensityComplexMat_, rhoMatMax );
      rhoMatMax.nzvalLocal = energyDensityComplexMat_.nzvalLocal;
      MPI_Bcast(rhoMatMax.nzvalLocal.Data(), rhoMatMax.nnzLocal, MPI_DOUBLE, idxMax, pointRowComm);

      SetValue( energyDensityComplexMat_.nzvalLocal, Z_ZERO );
      blas::Axpy( energyDensityComplexMat_.nnzLocal, facMin, rhoMatMin.nzvalLocal.Data(), 1,
          energyDensityComplexMat_.nzvalLocal.Data(), 1 );

      blas::Axpy( energyDensityComplexMat_.nnzLocal, facMax, rhoMatMax.nzvalLocal.Data(), 1,
          energyDensityComplexMat_.nzvalLocal.Data(), 1 );
    }
    // calculate the totalEnergyH_
    {
      Real local = 0.0;
      local = (blas::Dotc( HComplexMat_.nnzLocal,
            HComplexMat_.nzvalLocal.Data(),
            1, rhoComplexMat_.nzvalLocal.Data(), 1 )).real();

      mpi::Allreduce( &local, &totalEnergyH_, 1, MPI_SUM,
          gridPole_->rowComm );
    }

    //if( isEnergyDensityMatrix )
    if( isEDMCorrection_ )
    {

      Real local = 0.0;
      if( SMat.size != 0 ){
        local = (blas::Dotc( SMat.nnzLocal,
              SMat.nzvalLocal.Data(),
              1, energyDensityComplexMat_.nzvalLocal.Data(), 1 )).real();
      }
      else{
        CpxNumVec& nzval = energyDensityComplexMat_.nzvalLocal;
        for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
          local += nzval(diagIdxLocal_[i]).real();
        }
      }
      mpi::Allreduce( &local, &totalEnergyS_, 1, MPI_SUM,
          gridPole_->rowComm );
    }

  }

  // FIXME A placeholder for the FDM -- check check
  // if( isFreeEnergyDensityMatrix )
  if( isEDMCorrection_ )
  {
    DistSparseMatrix<Complex>& hmzMat  = freeEnergyDensityComplexMat_;
    totalFreeEnergy_                   = totalEnergyH_;
    blas::Copy( hmzMat.nnzLocal, rhoComplexMat_.nzvalLocal.Data(), 1,
        hmzMat.nzvalLocal.Data(), 1 );
  }

  muPEXSI = mu;

  for(l = 0; l < numShift; l++)
    if( NeVec[l] < numElectronExact - numElectronPEXSITolerance)
      muMin = shiftVec[l];

  for(l = numShift-1; l >= 0; l--)
    if( NeVec[l] > numElectronExact + numElectronPEXSITolerance)
      muMax = shiftVec[l];

  if ( muMax < muMin )
    ErrorHandling( " muMax < muMin , Error occured in the PEXSI ");

  MPI_Comm_free(&pointColComm);
  MPI_Comm_free(&pointRowComm);

}  // ---   interpolateDMComplex  ---


} //  namespace PEXSI
