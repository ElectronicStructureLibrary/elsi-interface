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
/// @file ppexsi_real_old.cpp
/// @brief Real arithmetic routines related to SuperLU_DIST.
/// 
/// This file contains the real arithmetic inertia counting routine and
/// the real arithmietic selected inversion interface.
///
/// This file is listed separately from other PEXSI routines because the
/// inertia count uses real arithmetic and may not be compatible with
/// the complex superlu_dist.
///
/// @date 2013-07-18 
/// @date 2014-01-26 Add selected inversion interface.

#include "ppexsi.hpp"
// NOTE: IMPORTANT: Since some macros in pselinv.hpp and superlu_ddefs.h
// share the same name (such as MYROW, MYCOL), superlu_ddefs.h MUST be
// included AFTER ppexsi.hpp
#include "pexsi/superlu_ddefs.h"
#include "pexsi/Cnames.h"

extern "C"{
void
pdsymbfact(superlu_options_t *options, SuperMatrix *A, 
    ScalePermstruct_t *ScalePermstruct, gridinfo_t *grid,
    LUstruct_t *LUstruct, SuperLUStat_t *stat, int *numProcSymbFact,
    int *info, double *totalMemory, double *maxMemory );
}


namespace PEXSI{

void PPEXSIData::CalculateNegativeInertiaReal(
    const std::vector<Real>&       shiftVec, 
    std::vector<Real>&             inertiaVec,
    const DistSparseMatrix<Real>&  HMat,
    const DistSparseMatrix<Real>&  SMat,
    std::string                    ColPerm,
    Int                            numProcSymbFact
    ){

  // *********************************************************************
  // Initialize
  // *********************************************************************

  // SuperLU Data structure

  SuperMatrix         A;
  superlu_options_t   options;
  ScalePermstruct_t   ScalePermstruct;          
  gridinfo_t*         grid;
  LUstruct_t          LUstruct;
  SOLVEstruct_t       SOLVEstruct;
  SuperLUStat_t       stat;
  Int                 info;

  // Initialize the SuperLU grid
  // Note: cannot use gridSuperLU_ because of the pImpl treatment.
  // Here it assumes that gridSuperLU_ and gridSelInv_ follow the same
  // distribution.

  grid = new gridinfo_t;
  superlu_gridinit( gridSelInv_->comm, gridSelInv_->numProcRow, 
      gridSelInv_->numProcCol, grid );

  // Get the diagonal indices for H and save it n diagIdxLocal_
  {
    Int numColLocal      = HMat.colptrLocal.m() - 1;
    Int numColLocalFirst = HMat.size / gridSelInv_->mpisize;
    Int firstCol         = gridSelInv_->mpirank * numColLocalFirst;

    diagIdxLocal_.clear();

    for( Int j = 0; j < numColLocal; j++ ){
      Int jcol = firstCol + j + 1;
      for( Int i = HMat.colptrLocal(j)-1; 
          i < HMat.colptrLocal(j+1)-1; i++ ){
        Int irow = HMat.rowindLocal(i);
        if( irow == jcol ){
          diagIdxLocal_.push_back( i );
        }
      }
    } // for (j)
  }


  // Set SuperLU options
  set_default_options_dist(&options);

  // The default value of ColPerm uses the default value from SuperLUOptions
  options.Fact              = DOFACT;
  options.RowPerm           = NOROWPERM; // IMPORTANT for symmetric matrices
  options.IterRefine        = NOREFINE;
  options.ParSymbFact       = NO;
  options.Equil             = NO; 
  options.ReplaceTinyPivot  = NO;
  // For output information such as # of nonzeros in L and U
  // and the memory cost, set PrintStat = YES
  options.PrintStat         = NO;
  options.SolveInitialized  = NO;

  if ( ColPerm == "NATURAL" ){
    options.ColPerm = NATURAL;
  } 
  else if( ColPerm == "MMD_AT_PLUS_A" ){
    options.ColPerm = MMD_AT_PLUS_A;
  }
  else if( ColPerm == "METIS_AT_PLUS_A" ){
    options.ColPerm = METIS_AT_PLUS_A;
  }
  else if( ColPerm == "PARMETIS" ){
    options.ColPerm           = PARMETIS;
    options.ParSymbFact       = YES;
  }


  // *********************************************************************
  // Symbolic factorization.  
  // Each numPoleGroup perform independently
  // *********************************************************************
  Real timeSta, timeEnd;
  GetTime( timeSta );
  // Generate the data pattern for the SuperMatrix A
  {
    Int numRowLocal = HMat.colptrLocal.m() - 1;
    Int numRowLocalFirst = HMat.size / gridSelInv_->mpisize;
    Int firstRow = gridSelInv_->mpirank * numRowLocalFirst;

    int_t *colindLocal, *rowptrLocal;
    double *nzvalLocal;
    rowptrLocal = (int_t*)intMalloc_dist(numRowLocal+1);
    colindLocal = (int_t*)intMalloc_dist(HMat.nnzLocal); 
    nzvalLocal  = (double*)doubleMalloc_dist(HMat.nnzLocal);

    std::copy( HMat.colptrLocal.Data(), HMat.colptrLocal.Data() + HMat.colptrLocal.m(),
        rowptrLocal );
    std::copy( HMat.rowindLocal.Data(), HMat.rowindLocal.Data() + HMat.rowindLocal.m(),
        colindLocal );

    // The value is not used here
    std::copy( HMat.nzvalLocal.Data(), HMat.nzvalLocal.Data() + HMat.nzvalLocal.m(),
        nzvalLocal );

    // Important to adjust from FORTRAN convention (1 based) to C convention (0 based) indices
    for(Int i = 0; i < HMat.rowindLocal.m(); i++){
      colindLocal[i]--;
    }

    for(Int i = 0; i < HMat.colptrLocal.m(); i++){
      rowptrLocal[i]--;
    }

    // Construct the distributed matrix according to the SuperLU_DIST format
    dCreate_CompRowLoc_Matrix_dist(&A, HMat.size, HMat.size, HMat.nnzLocal, 
        numRowLocal, firstRow,
        nzvalLocal, colindLocal, rowptrLocal,
        SLU_NR_loc, SLU_D, SLU_GE);
  }
  GetTime( timeEnd );
#if ( _DEBUGlevel_ >= 1 )
  statusOFS << "SuperMatrix is constructed." << std::endl;
  statusOFS << "Time for SuperMatrix construction is " <<
    timeEnd - timeSta << " [s]" << std::endl << std::endl;
#endif
  GetTime( timeSta );
  // Symbolic factorization
  {
    ScalePermstructInit(A.nrow, A.ncol, &ScalePermstruct);
    LUstructInit(A.nrow, A.ncol, &LUstruct);

    PStatInit(&stat);
#if ( _DEBUGlevel_ >= 1 )
    statusOFS << "Before symbfact subroutine." << std::endl;
#endif

    double totalMemory, maxMemory;

    pdsymbfact(&options, &A, &ScalePermstruct, grid, 
        &LUstruct, &stat, &numProcSymbFact, &info,
        &totalMemory, &maxMemory);
    PStatFree(&stat);

#if ( _DEBUGlevel_ >= 0 )
    statusOFS << "Memory cost of symbolic factorization (MB): " << std::endl;
    statusOFS << "Total: " << totalMemory << ", Average: " << 
      totalMemory / ( grid->nprow * grid->npcol )
      << ", Max: " << maxMemory << std::endl << std::endl;
#endif


  }
#if ( _DEBUGlevel_ >= 1 )
  statusOFS << "Symbolic factorization is finished." << std::endl;
#endif

  GetTime( timeEnd );
#if ( _DEBUGlevel_ >= 0 )
  statusOFS << "Time for symbolic factorization is " <<
    timeEnd - timeSta << " [s]" << std::endl << std::endl;
#endif

  // NOTE: A does not need to be destroyed.  The values of A can
  // be reused in later factorization.

  Real timeShiftSta, timeShiftEnd;

  Int numShift = shiftVec.size();
  std::vector<Real>  inertiaVecLocal(numShift);
  inertiaVec.resize(numShift);
  for(Int l = 0; l < numShift; l++){
    inertiaVecLocal[l] = 0.0;
    inertiaVec[l]      = 0.0;
  }

  for(Int l = 0; l < numShift; l++){
    if( gridPole_->mpirank / gridPole_->numProcCol == 
        l % gridPole_->numProcRow ){
      //  MYROW( gridPole_ ) == PROW( l, gridPole_ )

      GetTime( timeShiftSta );

#if ( _DEBUGlevel_ >= 0 )
      statusOFS << "Shift " << l << " = " << shiftVec[l] 
        << " processing..." << std::endl;
#endif

      // *********************************************************************
      // Data conversion
      // *********************************************************************

      {
        NRformat_loc *Astore = (NRformat_loc *) A.Store;
        Astore = (NRformat_loc *) A.Store;
        double* AnzvalLocal  = (double*)(Astore->nzval);
        if( SMat.size != 0 ){
          // S is not an identity matrix
          for( Int i = 0; i < HMat.nnzLocal; i++ ){
            AnzvalLocal[i] = HMat.nzvalLocal[i] - shiftVec[l] * SMat.nzvalLocal[i];
          }
        }
        else{
          // S is an identity matrix
          for( Int i = 0; i < HMat.nnzLocal; i++ ){
            AnzvalLocal[i] = HMat.nzvalLocal[i];
          }

          for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
            AnzvalLocal[ diagIdxLocal_[i] ] -= shiftVec[l];
          }
        } // if (SMat.size != 0 )
      }



      // *********************************************************************
      // Factorization
      // *********************************************************************

      Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

      GetTime( timeTotalFactorizationSta );

      // Data redistribution
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "Before Distribute." << std::endl;
#endif
      {
        Int* perm_c = ScalePermstruct.perm_c;
        NRformat_loc *Astore = (NRformat_loc *) A.Store;
        Int* colind   = Astore->colind;
        Int  nnzLocal = Astore->nnz_loc;

        // Important: recompute the permuted colind
        std::copy( HMat.rowindLocal.Data(), HMat.rowindLocal.Data() +
            HMat.rowindLocal.m(), colind );

        // Important to adjust from FORTRAN convention (1 based) to C convention (0 based) indices
        for(Int i = 0; i < HMat.rowindLocal.m(); i++){
          colind[i]--;
        }

        // Apply column permutation to the original distributed A
        for(Int j = 0; j < nnzLocal; j++)
          colind[j] = perm_c[colind[j]];
        // Distribute Pc*Pr*diag(R)*A*diag(C)*Pc' into L and U storage.  
        // NOTE: the row permutation Pc*Pr is applied internally in the
        // distribution routine. 
        pddistribute(SamePattern_SameRowPerm, A.nrow, 
            &A, &ScalePermstruct, NULL, &LUstruct, grid);
      }
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "After Distribute." << std::endl;
#endif

      // Numerical factorization
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "Before NumericalFactorize." << std::endl;
#endif
      {
        // Estimate the 1-norm
        char norm[1]; *norm = '1';
        double anorm = pdlangs( norm, &A, grid );

        PStatInit(&stat);
        pdgstrf(&options, A.nrow, A.ncol, 
            anorm, &LUstruct, grid, &stat, &info); 
        PStatFree(&stat);
        if( info ){
          std::ostringstream msg;
          msg << "Numerical factorization error, info =  " << info << std::endl;
          ErrorHandling( msg.str().c_str() );
        }
      }
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "After NumericalFactorize." << std::endl;
#endif

      GetTime( timeTotalFactorizationEnd );

#if ( _DEBUGlevel_ >= 0 )
      statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 
#endif

      // *********************************************************************
      // Compute inertia
      // *********************************************************************
      Real timeInertiaSta, timeInertiaEnd;
      GetTime( timeInertiaSta );

      // Compute the negative inertia of the matrix, and save it to
      // inertiaVecLocal[l]
      Real inertiaLocal = 0.0;
      {
        Int  n = A.ncol;
        Int  numSuper = LUstruct.Glu_persist -> supno[ n-1 ] + 1;
        Int *xsup = LUstruct.Glu_persist -> xsup;
        LocalLU_t* Llu   = LUstruct.Llu;
        for( Int jb = 0; jb < CEILING( numSuper, gridSelInv_->numProcCol ); jb++ ){
          Int bnum = GBj( jb, gridSelInv_ );
          if( bnum >= numSuper ) continue;

          Int cnt = 0;                                // Count for the index in LUstruct
          Int cntval = 0;                             // Count for the nonzero values
          const Int* index = Llu->Lrowind_bc_ptr[jb];
          // We only need the information from the diagonal block
          if( index ){ 
            Int  numBlock = index[cnt++];
            Int  lda      = index[cnt++];
            Int  blockIdx = index[cnt++];

            if( blockIdx == bnum ){
              // Diagonal block occurs on this processor.
              // By the convention in SuperLU, the diagonal block occurs
              // at the first block in L.
              Int  numRow   = index[cnt++];
              Int  numCol   = xsup[bnum+1] - xsup[bnum];
              // Check numRow == numCol
              if( numRow != numCol ){
                std::ostringstream msg;
                msg << "This is not a diagonal block." << std::endl 
                  << "blockIdx  = " << blockIdx << std::endl
                  << "numRow    = " << numRow 
                  << ", numCol = " << numCol << std::endl;
                ErrorHandling( msg.str().c_str() );
              }
              NumMat<Real>  nzval( numRow, numCol );

              lapack::Lacpy( 'A', numRow, numCol, 
                  (Real*)(Llu->Lnzval_bc_ptr[jb]), lda, 
                  nzval.Data(), numRow );

              for( Int i = 0; i < numRow; i++ ){
                if( nzval(i, i) < 0 )
                  inertiaLocal++;
              }
            }

          }  // if(index)
        } // for(jb)
      } // Compute the inertia

      mpi::Allreduce( &inertiaLocal, &inertiaVecLocal[l], 1, 
          MPI_SUM, gridPole_->rowComm );
    } // if I am in charge of this shift

  } // for(l)

  // Collect all the negative inertia together
  mpi::Allreduce( &inertiaVecLocal[0], &inertiaVec[0], numShift, 
      MPI_SUM, gridPole_->colComm );

  // Free the data
  Destroy_LU(A.ncol, grid, &LUstruct);
  LUstructFree(&LUstruct); 
  ScalePermstructFree(&ScalePermstruct);
  Destroy_CompRowLoc_Matrix_dist(&A);

  superlu_gridexit(grid);
  delete grid;


  return ;
} 		// -----  end of method PPEXSIData::CalculateNegativeInertiaReal ----- 

//extern "C"
//void PSelInvRealSymmetricInterface ( 
//    int           nrows,                        
//    int           nnz,                          
//    int           nnzLocal,                     
//    int           numColLocal,                  
//    int*          colptrLocal,                  
//    int*          rowindLocal,                  
//    double*       AnzvalLocal,                  
//    int           ordering,                
//    int           npSymbFact,                   
//    MPI_Comm	    comm,                         
//    int           nprow,
//    int           npcol,
//    double*       AinvnzvalLocal,
//    int*          info
//    )
//{
//	Int mpirank, mpisize;
//	MPI_Comm_rank( comm, &mpirank );
//	MPI_Comm_size( comm, &mpisize );
//	Real timeSta, timeEnd;
//	
//
//	// log files
//	std::stringstream  ss;
//	ss << "logPEXSI" << mpirank;
//	// append to previous log files
//	statusOFS.open( ss.str().c_str(), std::ios_base::app );
//		
//
//	try{
//		if( mpisize != nprow * npcol ){
ErrorHandling( " mpisize != nprow * npcol ." );
//		}
//
//    SuperMatrix         A;
//    superlu_options_t   options;
//    ScalePermstruct_t   ScalePermstruct;          
//    gridinfo_t*         grid;
//    LUstruct_t          LUstruct;
//    SOLVEstruct_t       SOLVEstruct;
//    SuperLUStat_t       stat;
//    Int                 info;
//
//    // Initialize the SuperLU grid
//    // Note: cannot use gridSuperLU_ because of the pImpl treatment.
//    // Here it assumes that gridSuperLU_ and gridSelInv_ follow the same
//    // distribution.
//
//    grid = new gridinfo_t;
//    superlu_gridinit( comm, nprow, npcol, grid );
//
//    // Set SuperLU options
//    set_default_options_dist(&options);
//
//    // The default value of ColPerm uses the default value from SuperLUOptions
//    options.Fact              = DOFACT;
//    options.RowPerm           = NOROWPERM; // IMPORTANT for symmetric matrices
//    options.IterRefine        = NOREFINE;
//    options.ParSymbFact       = NO;
//    options.Equil             = NO; 
//    options.ReplaceTinyPivot  = NO;
//    // For output information such as # of nonzeros in L and U
//    // and the memory cost, set PrintStat = YES
//    options.PrintStat         = NO;
//    options.SolveInitialized  = NO;
//
//		switch (ordering){
//			case 0:
//        options.ColPerm           = PARMETIS;
//        options.ParSymbFact       = YES;
//				break;
//			case 1:
//        options.ColPerm           = METIS_AT_PLUS_A;
//				break;
//			case 2:
//        options.ColPerm           = MMD_AT_PLUS_A;
//				break;
//			default:
ErrorHandling("Unsupported ordering strategy.");
//		}
//
//
//		// *********************************************************************
//		// Symbolic factorization 
//		// *********************************************************************
//
//    // Generate the data pattern for the SuperMatrix A
//    {
//      Int numRowLocal = numColLocal;
//      Int numRowLocalFirst = nrows / mpisize;
//      Int firstRow = mpirank * numRowLocalFirst;
//
//      int_t *colindLocal, *rowptrLocal;
//      double *nzvalLocal;
//      rowptrLocal = (int_t*)intMalloc_dist(numRowLocal+1);
//      colindLocal = (int_t*)intMalloc_dist(nnzLocal); 
//      nzvalLocal  = (double*)doubleMalloc_dist(nnzLocal);
//
//      std::copy( colptrLocal, colptrLocal + numColLocal + 1,
//          rowptrLocal );
//      std::copy( rowindLocal, rowindLocal + nnzLocal,
//          colindLocal );
//
//      // Important to adjust from FORTRAN convention (1 based) to C
//      // convention (0 based) indices
//      for(Int i = 0; i < nnzLocal; i++){
//        colindLocal[i]--;
//      }
//
//      for(Int i = 0; i < numColLocal + 1; i++){
//        rowptrLocal[i]--;
//      }
//
//      // Construct the distributed matrix according to the SuperLU_DIST format
//      dCreate_CompRowLoc_Matrix_dist(&A, nrows, nrows, nnzLocal, 
//          numRowLocal, firstRow,
//          nzvalLocal, colindLocal, rowptrLocal,
//          SLU_NR_loc, SLU_D, SLU_GE);
//    }
//
//    GetTime( timeSta );
//    // Symbolic factorization
//    {
//      ScalePermstructInit(A.nrow, A.ncol, &ScalePermstruct);
//      LUstructInit(A.nrow, A.ncol, &LUstruct);
//
//      PStatInit(&stat);
//#if ( _DEBUGlevel_ >= 1 )
//      statusOFS << "Before symbfact subroutine." << std::endl;
//#endif
//
//      double totalMemory, maxMemory;
//
//      pdsymbfact(&options, &A, &ScalePermstruct, grid, 
//          &LUstruct, &stat, &numProcSymbFact, &info,
//          &totalMemory, &maxMemory);
//      PStatFree(&stat);
//
//#if ( _DEBUGlevel_ >= 0 )
//      statusOFS << "Memory cost of symbolic factorization (MB): " << std::endl;
//      statusOFS << "Total: " << totalMemory << ", Average: " << 
//        totalMemory / ( grid->nprow * grid->npcol )
//        << ", Max: " << maxMemory << std::endl << std::endl;
//#endif
//    }
//#if ( _DEBUGlevel_ >= 1 )
//    statusOFS << "Symbolic factorization is finished." << std::endl;
//#endif
//
//    GetTime( timeEnd );
//#if ( _DEBUGlevel_ >= 0 )
//    statusOFS << "Time for symbolic factorization is " <<
//      timeEnd - timeSta << " [s]" << std::endl << std::endl;
//#endif
//		
//    // *********************************************************************
//		// Numerical factorization
//		// *********************************************************************
//    GetTime( timeSta );
//
//    // Numerical factorization
//#if ( _DEBUGlevel_ >= 1 )
//    statusOFS << "Before NumericalFactorize." << std::endl;
//#endif
//    {
//      // Estimate the 1-norm
//      char norm[1]; *norm = '1';
//      double anorm = pdlangs( norm, &A, grid );
//
//      PStatInit(&stat);
//      pdgstrf(&options, A.nrow, A.ncol, 
//          anorm, &LUstruct, grid, &stat, &info); 
//      PStatFree(&stat);
//      if( info ){
//        std::ostringstream msg;
//        msg << "Numerical factorization error, info =  " << info << std::endl;
ErrorHandling( msg.str().c_str() );
//      }
//    }
//#if ( _DEBUGlevel_ >= 1 )
//    statusOFS << "After NumericalFactorize." << std::endl;
//#endif
//
//    GetTime( timeEnd );
//
//#if ( _DEBUGlevel_ >= 0 )
//    statusOFS << "Time for total factorization is " 
//      << timeEnd - timeSta << " [s]" << std::endl; 
//#endif
//
//
//		// *********************************************************************
//		// Selected inversion
//		// *********************************************************************
//
//#if ( _DEBUGlevel_ >= 1 )
//    statusOFS << "After numerical factorization." << std::endl;
//#endif
//
//		GridType g1( comm, nprow, npcol );
//		SuperNodeType super;
//
//		luMat.SymbolicToSuperNode( super );
//		PMatrix PMloc( &g1, &super, &luOpt );
//		luMat.LUstructToPMatrix( PMloc );
//
//
//    // P2p communication version
//		PMloc.ConstructCommunicationPattern();
//
//    // Collective communication version
////		PMloc.ConstructCommunicationPattern_Collectives();
//
//		PMloc.PreSelInv();
//
//		// Main subroutine for selected inversion
//		// Use the broadcast pipelined version of SelInv.
//    
//    // P2p communication version
//    PMloc.SelInv();
//
//    // Collective communication version
////		PMloc.SelInv_Collectives();
//
//#if ( _DEBUGlevel_ >= 1 )
//    statusOFS << "After selected inversion." << std::endl;
//#endif
//
//		DistSparseMatrix<Complex>  AinvMat;       // A^{-1} in DistSparseMatrix format
//
//		PMloc.PMatrixToDistSparseMatrix2( AMat, AinvMat );
//#if ( _DEBUGlevel_ >= 1 )
//    statusOFS << "After conversion to DistSparseMatrix." << std::endl;
//#endif
//
//		// Convert the internal variables to output parameters
//
//		blas::Copy( nnzLocal*2, reinterpret_cast<double*>(AinvMat.nzvalLocal.Data()), 
//				1, AinvnzvalLocal, 1 );
//
//#if ( _DEBUGlevel_ >= 1 )
//    statusOFS << "After copying the matrix of Ainv." << std::endl;
//#endif
//
//
//
//		*info = 0;
//	}
//	catch( std::exception& e ) {
//		statusOFS << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
//			<< std::endl << e.what() << std::endl;
////		std::cerr  << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
////			<< std::endl << e.what() << std::endl;
//		*info = 1;
//	}
//	
//	// Synchronize the info among all processors. 
//	// If any processor gets error message, info = 1
//	Int infoAll = 0;
//	mpi::Allreduce( info, &infoAll, 1, MPI_MAX, comm  );
//	*info = infoAll;
//
//	statusOFS.close();
//	return;
//}  // -----  end of function PSelInvRealSymmetricInterface ----- 


} //  namespace PEXSI

