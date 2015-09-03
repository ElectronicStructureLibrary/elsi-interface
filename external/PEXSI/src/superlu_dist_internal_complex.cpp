/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

Author: Lin Lin and Mathias Jacquelin

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
/// @file superlu_dist_internal_complex.cpp
/// @brief Implementation of internal structures for interfacing with SuperLU_Dist (version 3.0 and later) for complex arithmetic
/// @date 2014-03-17
#include "pexsi/SuperLUGrid.hpp"
#include "pexsi/superlu_dist_internal.hpp"

#include "pexsi/pselinv.hpp"

#include <superlu_zdefs.h>

#include <Cnames.h>
extern "C"{ void
pzsymbfact(superlu_options_t *options, SuperMatrix *A, 
    ScalePermstruct_t *ScalePermstruct, gridinfo_t *grid,
    LUstruct_t *LUstruct, SuperLUStat_t *stat, int *numProcSymbFact,
    int *info, double *totalMemory, double *maxMemory );
}




// SuperLUGrid class
namespace PEXSI{


  class ComplexGridInfo{
    friend class ComplexGridData;
    friend class ComplexSuperLUData_internal;
    protected:
    gridinfo_t          grid;
  };

  ComplexGridData::ComplexGridData(){
    info_ = new ComplexGridInfo;
  }

  ComplexGridData::~ComplexGridData(){
    delete info_;
  }



  ComplexGridData::ComplexGridData(const ComplexGridData & g)
  {
    this->info_ = new ComplexGridInfo;
    this->info_->grid = g.info_->grid;
  }

  ComplexGridData & ComplexGridData::operator = (const ComplexGridData & g)
  {
    //if this is the same object, skip the thing
    if(&g != this){
      delete info_;
      this->info_ = new ComplexGridInfo;
      this->info_->grid = g.info_->grid;
    }

    return *this;
  }

  void ComplexGridData::GridInit( MPI_Comm comm, Int nprow, Int npcol ){
    superlu_gridinit(comm, nprow, npcol, &info_->grid);
  }

  void ComplexGridData::GridExit(  ){
    superlu_gridexit(&info_->grid);
  }

}



// SuperLUData class
namespace PEXSI{

  class ComplexSuperLUData_internal{
    friend class ComplexSuperLUData;
    protected:
    /// @brief SuperLU matrix. 
    SuperMatrix         A;                        

    /// @brief SuperLU options. 
    ///
    /// Note
    /// ----
    ///
    /// It is important to have 
    ///
    /// options.RowPerm           = NOROWPERM;
    /// 
    /// to make sure that symmetric permutation is used.
    ///
    superlu_options_t   options;                  

    /// @brief Saves the permutation vectors.  Only perm_c (permutation of
    /// column as well as rows due to the symmetric permutation) will be used.
    ScalePermstruct_t   ScalePermstruct;          

    /// @brief SuperLU grid structure.
    gridinfo_t*         grid;

    /// @brief Saves the supernodal partition as well as the numerical
    /// values and structures of the L and U structure.
    LUstruct_t          LUstruct;

    /// @brief Used for solve for multivectors.
    SOLVEstruct_t       SOLVEstruct;

    /// @brief SuperLU statistics
    SuperLUStat_t       stat;

    /// @brief Number of processors used for parallel symbolic
    /// factorization and PARMETIS/PT-SCOTCH
    Int                 numProcSymbFact;

    /// @brief SuperLU information
    Int                 info;

    // The following are for consistency checks
    bool                isSuperMatrixAllocated;
    bool                isSuperMatrixFactorized;
    bool                isScalePermstructAllocated;
    bool                isLUstructAllocated;

    Int maxDomains;

    ComplexSuperLUData_internal(const SuperLUGrid<Complex>& g, const SuperLUOptions& opt);
    ~ComplexSuperLUData_internal();
    ComplexSuperLUData_internal(const ComplexSuperLUData_internal& g);
    ComplexSuperLUData_internal & operator = (const ComplexSuperLUData_internal& g);

    void DestroyAOnly();

  };

  ComplexSuperLUData_internal::ComplexSuperLUData_internal(const SuperLUGrid<Complex>& g, const SuperLUOptions& opt){

    isSuperMatrixAllocated     = false;
    isScalePermstructAllocated = false;
    isLUstructAllocated        = false;
    numProcSymbFact            = opt.numProcSymbFact;

    // Options
    set_default_options_dist(&options);

    // The default value of ColPerm uses the default value from SuperLUOptions
    options.Fact              = DOFACT;
    if(opt.RowPerm == "LargeDiag"){
      options.RowPerm         = LargeDiag;
    }
    else{
      options.RowPerm         = NOROWPERM;
    }

    options.IterRefine        = NOREFINE;
    options.ParSymbFact       = NO;
    options.Equil             = NO; 
    options.ReplaceTinyPivot  = YES;
    // For output information such as # of nonzeros in L and U
    // and the memory cost, set PrintStat = YES
    options.PrintStat         = NO;
    options.SolveInitialized  = NO;
    // Necessary to invoke static scheduling of SuperLU
    options.lookahead_etree   = YES;
    options.SymPattern        = YES;

    if ( opt.ColPerm == "NATURAL" ){
      options.ColPerm = NATURAL;
    } 
    else if( opt.ColPerm == "MMD_AT_PLUS_A" ){
      options.ColPerm = MMD_AT_PLUS_A;
    }
    else if( opt.ColPerm == "METIS_AT_PLUS_A" ){
      options.ColPerm = METIS_AT_PLUS_A;
    }
    else if( opt.ColPerm == "PARMETIS" ){
      options.ColPerm           = PARMETIS;
      options.ParSymbFact       = YES;
    }
    else{
      std::ostringstream msg;
      msg << opt.ColPerm << " is not a supported ColPerm type. Try (case sensitive) " << std::endl
        << "NATURAL | MMD_AT_PLUS_A | METIS_AT_PLUS_A | PARMETIS" << std::endl;
#ifdef USE_ABORT
      abort();
#endif
      throw std::runtime_error( msg.str().c_str() );
    }

    // Setup grids
    grid = &(g.ptrData->info_->grid);
  }

  ComplexSuperLUData_internal::~ComplexSuperLUData_internal(){
    if( isLUstructAllocated ){
      Destroy_LU(A.ncol, grid, &LUstruct);
      LUstructFree(&LUstruct); 
    }
    if( isScalePermstructAllocated ){
      ScalePermstructFree(&ScalePermstruct);
    }
    if( options.SolveInitialized ){
      // TODO real arithmetic
      zSolveFinalize(&options, &SOLVEstruct);
    }
    if( isSuperMatrixAllocated ){
      DestroyAOnly();
    }
  }


  ComplexSuperLUData_internal::ComplexSuperLUData_internal(const ComplexSuperLUData_internal& g){
    memcpy(this,&g,sizeof(ComplexSuperLUData_internal));
  }

  ComplexSuperLUData_internal & ComplexSuperLUData_internal::operator = (const ComplexSuperLUData_internal& g){
    if(this!=&g){
      memcpy(this,&g,sizeof(ComplexSuperLUData_internal));
    }

    return *this;
  }




  void ComplexSuperLUData_internal::DestroyAOnly	(  )
  {
#ifndef _RELEASE_
    PushCallStack("ComplexSuperLUData_internal::DestroyAOnly");
#endif
    if( isSuperMatrixAllocated == false ){
#ifdef USE_ABORT
      abort();
#endif
      throw std::logic_error( "SuperMatrix has not been allocated." );
    }
    switch ( A.Stype ){
      case SLU_NC:
        Destroy_CompCol_Matrix_dist(&A);
        break;
      case SLU_NR_loc:
        Destroy_CompRowLoc_Matrix_dist(&A);
        break;
      default:
        std::ostringstream msg;
        msg << "Type " << SLU_NR_loc << " is to be destroyed" << std::endl
          << "This is an unsupported SuperMatrix format to be destroyed." << std::endl;
#ifdef USE_ABORT
        abort();
#endif
        throw std::runtime_error( msg.str().c_str() );
    }
    isSuperMatrixAllocated = false;
#ifndef _RELEASE_
    PopCallStack();
#endif

    return ;
  } 		// -----  end of method ComplexSuperLUData_internal::DestroyAOnly  ----- 






  ComplexSuperLUData::ComplexSuperLUData( const SuperLUGrid<Complex>& g, const SuperLUOptions& opt ){
#ifndef _RELEASE_
    PushCallStack("ComplexSuperLUData::ComplexSuperLUData");
#endif

    ptrData = new ComplexSuperLUData_internal(g,opt);
    if( ptrData == NULL ){
#ifdef USE_ABORT
      abort();
#endif
      throw std::runtime_error( "SuperLUMatrix cannot be allocated." );
    }

#ifndef _RELEASE_
    PopCallStack();
#endif
  }


  ComplexSuperLUData::~ComplexSuperLUData(){
#ifndef _RELEASE_
    PushCallStack("ComplexSuperLUData::~ComplexSuperLUData");
#endif
    delete ptrData;

#ifndef _RELEASE_
    PopCallStack();
#endif
  }



  ComplexSuperLUData::ComplexSuperLUData(const ComplexSuperLUData & g)
  {
#ifndef _RELEASE_
    PushCallStack("ComplexSuperLUData::ComplexSuperLUData");
#endif

    if( g.ptrData == NULL ){
      throw std::runtime_error( "Copied SuperLUMatrix is not allocated." );
    }

    ptrData = new ComplexSuperLUData_internal(*g.ptrData);

    if( ptrData == NULL ){
#ifdef USE_ABORT
      abort();
#endif
      throw std::runtime_error( "SuperLUMatrix cannot be allocated." );
    }

#ifndef _RELEASE_
    PopCallStack();
#endif
  }

  ComplexSuperLUData & ComplexSuperLUData::operator = (const ComplexSuperLUData & g)
  {
    //if this is the same object, skip the thing
    if(&g == this){
      return *this;
    }

    if( g.ptrData == NULL ){
      throw std::runtime_error( "Copied SuperLUMatrix is not allocated." );
    }

    delete ptrData;
    ptrData = new ComplexSuperLUData_internal(*g.ptrData);
    if( ptrData == NULL ){
#ifdef USE_ABORT
      abort();
#endif
      throw std::runtime_error( "SuperLUMatrix cannot be allocated." );
    }

    return *this;
  }



  Int ComplexSuperLUData::m (  ) const	
  {
    return ptrData->A.nrow;
  } 		// -----  end of method ComplexSuperLUData::m  ----- 


  Int ComplexSuperLUData::n (  ) const	
  {
    return ptrData->A.ncol;
  } 		// -----  end of method ComplexSuperLUData::n  ----- 



  void ComplexSuperLUData::DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Complex>& sparseA, const SuperLUOptions & options ){
#ifndef _RELEASE_
    PushCallStack( "ComplexSuperLUData::DistSparseMatrixToSuperMatrixNRloc" );
#endif
    if( ptrData->isSuperMatrixAllocated == true ){
#ifdef USE_ABORT
      abort();
#endif
      throw std::logic_error( "SuperMatrix is already allocated." );
    }
    gridinfo_t* grid = ptrData->grid;

    Int mpirank = grid->iam;
    Int mpisize = grid->nprow * grid->npcol;



    int_t *colindLocal, *rowptrLocal;
    doublecomplex *nzvalLocal;

    Int numRowLocalFirst = sparseA.size / mpisize;
    Int firstRow = mpirank * numRowLocalFirst;
    Int numRowLocal = -1;
    Int nnzLocal = -1;

    numRowLocal = sparseA.colptrLocal.m() - 1;
    nnzLocal = sparseA.nnzLocal;

    colindLocal = (int_t*)intMalloc_dist(sparseA.nnzLocal); 
    nzvalLocal  = (doublecomplex*)doublecomplexMalloc_dist(sparseA.nnzLocal);
    rowptrLocal = (int_t*)intMalloc_dist(numRowLocal+1);

    std::copy( sparseA.colptrLocal.Data(), sparseA.colptrLocal.Data() + sparseA.colptrLocal.m(),
        rowptrLocal );
    std::copy( sparseA.rowindLocal.Data(), sparseA.rowindLocal.Data() + sparseA.rowindLocal.m(),
        colindLocal );
    std::copy( sparseA.nzvalLocal.Data(), sparseA.nzvalLocal.Data() + sparseA.nzvalLocal.m(),
        (Complex*)nzvalLocal );


    // Important to adjust from FORTRAN convention (1 based) to C convention (0 based) indices
    for(Int i = 0; i < nnzLocal; i++){
      colindLocal[i]--;
    }

    for(Int i = 0; i < numRowLocal+1; i++){
      rowptrLocal[i]--;
    }

    // Construct the distributed matrix according to the SuperLU_DIST format
    zCreate_CompRowLoc_Matrix_dist(&ptrData->A, sparseA.size, sparseA.size, nnzLocal, 
        numRowLocal, firstRow,
        nzvalLocal, colindLocal, rowptrLocal,
        SLU_NR_loc, SLU_Z, SLU_GE);

    ptrData->isSuperMatrixAllocated = true;

#ifndef _RELEASE_
    PopCallStack();
#endif

    return;

  } 		// -----  end of method ComplexSuperLUData::DistSparseMatrixToSuperMatrixNRloc ----- 


  void
    ComplexSuperLUData::DestroyAOnly	(  )
    {
#ifndef _RELEASE_
      PushCallStack("ComplexSuperLUData::DestroyAOnly");
#endif
      ptrData->DestroyAOnly();
#ifndef _RELEASE_
      PopCallStack();
#endif
      return ;
    } 		// -----  end of method ComplexSuperLUData::DestroyAOnly  ----- 

  void
    ComplexSuperLUData::SymbolicFactorize	(  )
    {
#ifndef _RELEASE_
      PushCallStack("ComplexSuperLUData::SymbolicFactorize");
#endif
      if( ptrData->isScalePermstructAllocated ){
#ifdef USE_ABORT
        abort();
#endif
        throw std::logic_error( "ScalePermstruct is already allocated." );
      }
      if( ptrData->isLUstructAllocated){
#ifdef USE_ABORT
        abort();
#endif
        throw std::logic_error( "LUstruct is already allocated." );
      }
      //      if(ptrData->options.RowPerm != NOROWPERM ){
      //#ifdef USE_ABORT
      //        abort();
      //#endif
      //        throw std::logic_error( "For PEXSI there must be no row permutation." );
      //      }

      SuperMatrix&  A = ptrData->A;

      ScalePermstructInit(A.nrow, A.ncol, &ptrData->ScalePermstruct);
      LUstructInit(A.nrow, A.ncol, &ptrData->LUstruct);

      PStatInit(&ptrData->stat);
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "Before symbfact subroutine." << std::endl;
#endif

      double totalMemory = 0.0, maxMemory = 0.0;

      pzsymbfact(&ptrData->options, &A, &ptrData->ScalePermstruct, ptrData->grid, 
          &ptrData->LUstruct, &ptrData->stat, &ptrData->numProcSymbFact, &ptrData->info,
          &totalMemory, &maxMemory);
      PStatFree(&ptrData->stat);


      //assert(ptrData->ScalePermstruct.DiagScale == BOTH );


#if ( _DEBUGlevel_ >= 0 )
      statusOFS << "Memory cost of symbolic factorization (MB): " << std::endl;
      statusOFS << "Total: " << totalMemory << ", Average: " << 
        totalMemory / ( ptrData->grid->nprow * ptrData->grid->npcol )
        << ", Max: " << maxMemory << std::endl << std::endl;
#endif



      ptrData->isScalePermstructAllocated = true;
      ptrData->isLUstructAllocated        = true; 

#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method ComplexSuperLUData::SymbolicFactorize  ----- 


  void
    ComplexSuperLUData::Distribute	(  )
    {
#ifndef _RELEASE_
      PushCallStack("ComplexSuperLUData::Distribute");
#endif
      if( ptrData->isScalePermstructAllocated == false ){
#ifdef USE_ABORT
        abort();
#endif
        throw std::logic_error( "ScalePermstruct has not been allocated by SymbolicFactorize." );
      }	
      if( ptrData->isLUstructAllocated == false ){
#ifdef USE_ABORT
        abort();
#endif
        throw std::logic_error( "LUstruct has not been allocated by SymbolicFactorize." );
      }	
      if( ptrData->isSuperMatrixAllocated == false ){
#ifdef USE_ABORT
        abort();
#endif
        throw std::logic_error( "SuperMatrix has not been allocated." );
      }	

      Int* perm_c = ptrData->ScalePermstruct.perm_c;
      NRformat_loc *Astore = (NRformat_loc *) ptrData->A.Store;
      Int* colind = Astore->colind;
      Int  nnzLocal = Astore->nnz_loc;
      // Apply column permutation to the original distributed A
      for(Int j = 0; j < nnzLocal; j++)
        colind[j] = perm_c[colind[j]];
      // Distribute Pc*Pr*diag(R)*A*diag(C)*Pc' into L and U storage.  
      // NOTE: the row permutation Pc*Pr is applied internally in the
      // distribution routine. 
      float dist_mem_use = pzdistribute(SamePattern_SameRowPerm, ptrData->A.nrow, 
          &ptrData->A, &ptrData->ScalePermstruct, NULL, &ptrData->LUstruct, ptrData->grid);

#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method ComplexSuperLUData::Distribute  ----- 


  void
    ComplexSuperLUData::NumericalFactorize	(  )
    {
#ifndef _RELEASE_
      PushCallStack("ComplexSuperLUData::NumericalFactorize");
#endif
      if( !ptrData->isLUstructAllocated ){
#ifdef USE_ABORT
        abort();
#endif
        throw std::logic_error( "LUstruct has not been allocated." );
      }
      // Estimate the 1-norm
      char norm[1]; *norm = '1';
      double anorm = pzlangs( norm, &ptrData->A, ptrData->grid );

      PStatInit(&ptrData->stat);
      pzgstrf(&ptrData->options, ptrData->A.nrow, ptrData->A.ncol, 
          anorm, &ptrData->LUstruct, ptrData->grid, &ptrData->stat, &ptrData->info); 

#ifdef _PRINT_STATS_
      statusOFS<<"******************SUPERLU STATISTICS****************"<<std::endl;
      statusOFS<<"Number of tiny pivots: "<<ptrData->stat.TinyPivots<<std::endl;
      statusOFS<<"Number of look aheads: "<<ptrData->stat.num_look_aheads<<std::endl;
      statusOFS<<"Number of FLOPS during factorization: "<<ptrData->stat.ops[FACT]<<std::endl;
      statusOFS<<"****************************************************"<<std::endl;
#endif




      PStatFree(&ptrData->stat);
      if( ptrData->info != 0 ){
        std::ostringstream msg;
        msg << "Numerical factorization error, info =  " << ptrData->info << std::endl;
#ifdef USE_ABORT
        printf("P%d %s\n",ptrData->grid->iam,msg.str().c_str());
        abort();
#endif
        throw std::runtime_error( msg.str().c_str() );
      }

      //assert(ptrData->ScalePermstruct.DiagScale == BOTH );

      // Prepare for Solve.
      ptrData->options.Fact = FACTORED;

#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method ComplexSuperLUData::NumericalFactorize  ----- 


  void
    ComplexSuperLUData::ConvertNRlocToNC	( ComplexSuperLUData * aptrData )
    {
#ifndef _RELEASE_
      PushCallStack("ComplexSuperLUData::ConvertNRlocToNC");
#endif
      if( !ptrData->isSuperMatrixAllocated ){
#ifdef USE_ABORT
        abort();
#endif
        throw std::runtime_error( "The local SuperMatrix has not been allocated." );
      }
      if( aptrData->ptrData->isSuperMatrixAllocated ){
#ifdef USE_ABORT
        abort();
#endif
        throw std::runtime_error( "The global SuperMatrix has been allocated." );
      }
      // TODO make sure the two grids are the same

      // TODO real arithmetic
      const Int NEED_VALUE = 1;
      pzCompRow_loc_to_CompCol_global(NEED_VALUE, &ptrData->A, ptrData->grid, 
          &aptrData->ptrData->A);

      ptrData->isSuperMatrixAllocated = true;

#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method ComplexSuperLUData::ConvertNRlocToNC  ----- 

  void
    ComplexSuperLUData::MultiplyGlobalMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& bGlobal )
    {
#ifndef _RELEASE_
      PushCallStack("ComplexSuperLUData::MultiplyGlobalMultiVector");
#endif
      char trans[1]; *trans = 'N';
      Int m = xGlobal.m();
      Int nrhs = xGlobal.n();
      bGlobal.Resize( m, nrhs );

      if( ptrData->A.Stype != SLU_NC ){
        std::ostringstream msg;
        msg << "MultiplyGlobalMultiVector requires SLU_NC matrix with type " << SLU_NC << std::endl
          << "The matrix is of type " << ptrData->A.Stype << std::endl
          << "Consider using ConvertNRlocToNC subroutine" << std::endl;
#ifdef USE_ABORT
        abort();
#endif
        throw std::runtime_error( msg.str().c_str() );
      }	
      zFillRHS_dist(trans, nrhs, (doublecomplex*)xGlobal.Data(), m, 
          &ptrData->A, (doublecomplex*) bGlobal.Data(), m);
#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method ComplexSuperLUData::MultiplyGlobalMultiVector  ----- 


  void
    ComplexSuperLUData::DistributeGlobalMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal )
    {
#ifndef _RELEASE_
      PushCallStack("ComplexSuperLUData::DistributeGlobalMultiVector");
#endif
      SuperMatrix& A = ptrData->A;
      if( ptrData->A.Stype != SLU_NR_loc ){
        std::ostringstream msg;
        msg << "DistributeGlobalMultiVector requires SLU_NR_loc matrix with type " << SLU_NR_loc << std::endl
          << "The matrix is of type " << ptrData->A.Stype << std::endl;
#ifdef USE_ABORT
        abort();
#endif
        throw std::runtime_error( msg.str().c_str() );
      }	

      NRformat_loc *Astore = (NRformat_loc *) ptrData->A.Store;

      Int numRowLocal = Astore->m_loc;
      Int firstRow    = Astore->fst_row;
      Int nrhs = xGlobal.n();

      xLocal.Resize( numRowLocal, nrhs );
      SetValue( xLocal, ZERO<Complex>() );
      for( Int j = 0; j < nrhs; j++ ){
        std::copy( xGlobal.VecData(j)+firstRow, xGlobal.VecData(j)+firstRow+numRowLocal,
            xLocal.VecData(j) );
      }

#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method ComplexSuperLUData::DistributeGlobalMultiVector  ----- 


  void ComplexSuperLUData::GatherDistributedMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal )
  {
#ifndef _RELEASE_
    PushCallStack("ComplexSuperLUData::GatherDistributedMultiVector");
#endif
    SuperMatrix& A = ptrData->A;
    if( ptrData->A.Stype != SLU_NR_loc ){
      std::ostringstream msg;
      msg << "GatherDistributedMultiVector requires SLU_NR_loc matrix with type " << SLU_NR_loc << std::endl
        << "The matrix is of type " << ptrData->A.Stype << std::endl;
#ifdef USE_ABORT
      abort();
#endif
      throw std::runtime_error( msg.str().c_str() );
    }	

    NRformat_loc *Astore = (NRformat_loc *) ptrData->A.Store;

    Int numRowLocal = Astore->m_loc;
    Int firstRow    = Astore->fst_row;
    Int nrhs = xGlobal.n();

    Int maxRows = 0;
    Int localRows = xLocal.m();

    MPI_Allreduce(&localRows,&maxRows,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);



    NumMat<Complex> tmpLocal(xGlobal.m(),nrhs);
    SetValue( tmpLocal, ZERO<Complex>() );
    for( Int j = 0; j < nrhs; j++ ){
      std::copy( xLocal.VecData(j), xLocal.VecData(j)+numRowLocal, tmpLocal.VecData(j)+firstRow );
    }

    MPI_Allreduce(tmpLocal.Data(),xGlobal.Data(),2*tmpLocal.m()*tmpLocal.n(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

#ifndef _RELEASE_
    PopCallStack();
#endif

    return ;
  } 		// -----  end of method ComplexSuperLUData::GatherDistributedMultiVector  ----- 


  void
    ComplexSuperLUData::SolveDistMultiVector	( NumMat<Complex>& bLocal, DblNumVec& berr )
    {
#ifndef _RELEASE_
      PushCallStack("ComplexSuperLUData::SolveDistMultiVector");
#endif
      Int nrhs = bLocal.n();
      NRformat_loc *Astore = (NRformat_loc *) ptrData->A.Store;
      Int numRowLocal = Astore->m_loc;
      Int firstRow    = Astore->fst_row;

      berr.Resize( nrhs );

      // TODO Complex arithmetic

      PStatInit(&ptrData->stat);
      pzgssvx(&ptrData->options, &ptrData->A, &ptrData->ScalePermstruct, 
          (doublecomplex*)bLocal.Data(), numRowLocal, nrhs, ptrData->grid,
          &ptrData->LUstruct, &ptrData->SOLVEstruct, berr.Data(), 
          &ptrData->stat, &ptrData->info);
      PStatFree(&ptrData->stat); 

      if ( ptrData->options.SolveInitialized ) {
        zSolveFinalize(&ptrData->options, &ptrData->SOLVEstruct);
      }

      if( ptrData->info ){
        std::ostringstream msg;
        msg << "Numerical solve error, info =  " << ptrData->info << std::endl;
#ifdef USE_ABORT
        abort();
#endif
        throw std::runtime_error( msg.str().c_str() );
      }




#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method ComplexSuperLUData::SolveDistMultiVector  ----- 


  void
    ComplexSuperLUData::CheckErrorDistMultiVector	( NumMat<Complex>& xLocal, NumMat<Complex>& xTrueLocal )
    {
#ifndef _RELEASE_
      PushCallStack("ComplexSuperLUData::CheckErrorDistMultiVector");
#endif
      Int nrhs = xLocal.n();
      NRformat_loc *Astore = (NRformat_loc *) ptrData->A.Store;
      Int numRowLocal = Astore->m_loc;

      pzinf_norm_error(ptrData->grid->iam, numRowLocal, nrhs, 
          (doublecomplex*)xLocal.Data(), numRowLocal, 
          (doublecomplex*)xTrueLocal.Data(), numRowLocal, ptrData->grid);

#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method ComplexSuperLUData::CheckErrorDistMultiVector  ----- 


  void
    ComplexSuperLUData::LUstructToPMatrix	( PMatrix<Complex>& PMloc )
    {
#ifndef _RELEASE_
      PushCallStack("ComplexSuperLUData::LUstructToPMatrix");
#endif
      const LocalLU_t* Llu   = ptrData->LUstruct.Llu;
      const GridType* grid   = PMloc.Grid();
      const SuperNodeType* super = PMloc.SuperNode();
      Int numSuper = PMloc.NumSuper();

      PMloc.ColBlockIdx().clear();
      PMloc.RowBlockIdx().clear();
      PMloc.ColBlockIdx().resize( PMloc.NumLocalBlockCol() );
      PMloc.RowBlockIdx().resize( PMloc.NumLocalBlockRow() );










      // L part   
#ifndef _RELEASE_
      PushCallStack("L part");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "LUstructToPMatrix::L part" << std::endl;
#endif


      for( Int jb = 0; jb < PMloc.NumLocalBlockCol(); jb++ ){
        Int bnum = GBj( jb, grid );
        if( bnum >= numSuper ) continue;

        Int cnt = 0;                                // Count for the index in LUstruct
        Int cntval = 0;                             // Count for the nonzero values
        Int cntidx = 0;                             // Count for the nonzero block indexes
        const Int* index = Llu->Lrowind_bc_ptr[jb];
        if( index ){
          // Not an empty column, start a new column then.
          std::vector<LBlock<Complex> >& Lcol = PMloc.L(jb);
          Lcol.resize( index[cnt++] );


          Int lda = index[cnt++];

          Int savCnt = cnt;
          IntNumVec LIndices(Lcol.size());
          for( Int iblk = 0; iblk < Lcol.size(); iblk++ ){
            Int blockIdx    = index[cnt++];
            Int numRow      = index[cnt++];
            cnt += numRow;
            LIndices[iblk] = blockIdx;
          }
          //      statusOFS<<"Unsorted blockidx for L are: "<<LIndices<<std::endl;
          //sort the array
          IntNumVec sortedIndices(Lcol.size());
          for(Int i = 0; i<sortedIndices.m(); ++i){ sortedIndices[i] = i;}
          ULComparator cmp2(LIndices);
          //sort the row indices (so as to speedup the index lookup
          std::sort(sortedIndices.Data(),sortedIndices.Data()+sortedIndices.m(),cmp2);

          //      statusOFS<<"Sorted indices for L are: "<<sortedIndices<<std::endl;
          //get the inverse perm
          IntNumVec tmp = sortedIndices;
          for(Int i = 0; i<sortedIndices.m(); ++i){ tmp[sortedIndices[i]] = i;}
          sortedIndices = tmp;

          cnt = savCnt;
          for( Int iblk = 0; iblk < Lcol.size(); iblk++ ){
            Int blockIdx    = index[cnt++];
            //determine where to put it
            LBlock<Complex> & LB     = Lcol[sortedIndices[iblk]];
            LB.blockIdx    = blockIdx;

            PMloc.ColBlockIdx(jb).push_back(LB.blockIdx);
            Int LBi = LB.blockIdx / grid->numProcRow; 
            PMloc.RowBlockIdx( LBi ).push_back( bnum );


            LB.numRow      = index[cnt++];
            LB.numCol      = super->superPtr[bnum+1] - super->superPtr[bnum];
            LB.rows        = IntNumVec( LB.numRow, true, const_cast<Int*>(&index[cnt]) );

            IntNumVec rowsPerm(LB.numRow);
            IntNumVec rowsSorted(LB.numRow);
            for(Int i = 0; i<rowsPerm.m(); ++i){ rowsPerm[i] = i;}
            ULComparator cmp(LB.rows);
            //sort the row indices (so as to speedup the index lookup
            std::sort(rowsPerm.Data(),rowsPerm.Data()+rowsPerm.m(),cmp);

            for(Int i = 0; i<LB.rows.m(); ++i){ 
              rowsSorted[i] = LB.rows[rowsPerm[i]];
            }

            LB.rows = rowsSorted;

            LB.nzval.Resize( LB.numRow, LB.numCol );   
            SetValue( LB.nzval, ZERO<Complex>() ); 
            cnt += LB.numRow;

            //sort the nzval
            for(Int j = 0; j<LB.numCol; ++j){
              for(Int i = 0; i<LB.numRow; ++i){
                LB.nzval(i,j) = ((Complex*)(Llu->Lnzval_bc_ptr[jb]+cntval))[rowsPerm[i]+j*lda];
              }
            }

            cntval += LB.numRow;


#if ( _DEBUGlevel_ >= 1 )
            statusOFS 
              << "L part: bnum = " << bnum << ", numBlock = " << Lcol.size()
              << ", blockIdx = " << LB.blockIdx
              << ", numRow = " << LB.numRow 
              << ", numCol = " << LB.numCol << std::endl;
#endif 

          } // for(iblk)


#if ( _DEBUGlevel_ >= 2 )
          statusOFS<<"Real Sorted blockidx for L are: ";
          for( Int iblk = 0; iblk < Lcol.size(); iblk++ ){
            statusOFS<<Lcol[iblk].blockIdx<<" ";
          }
          statusOFS<<std::endl;
#endif


        }  // if(index)
      } // for(jb)
#ifndef _RELEASE_
      PopCallStack();
#endif






      // U part
#ifndef _RELEASE_
      PushCallStack("U part");
#endif
#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "LUstructToPMatrix::U part" << std::endl;
#endif
      for( Int ib = 0; ib < PMloc.NumLocalBlockRow(); ib++ ){
        Int bnum = GBi( ib, grid );
        if( bnum >= numSuper ) continue;

        Int cnt = 0;                                // Count for the index in LUstruct
        Int cntval = 0;                             // Count for the nonzero values
        Int cntidx = 0;                             // Count for the nonzero block indexes
        const Int*    index = Llu->Ufstnz_br_ptr[ib]; 
        const Complex* pval  = reinterpret_cast<const Complex*>(Llu->Unzval_br_ptr[ib]);
        if( index ){ 
          // Not an empty row
          // Compute the number of nonzero columns 

          std::vector<UBlock<Complex> >& Urow = PMloc.U(ib);
          Urow.resize( index[cnt++] );
          cnt = BR_HEADER;


          std::vector<Int> cols;                    //Save the nonzero columns in the current block
          for(Int jblk = 0; jblk < Urow.size(); jblk++ ){
            cols.clear();

            Int blockIdx = index[cnt];
            Int LBj = blockIdx / grid->numProcCol; 
            UBlock<Complex> & UB = Urow[jblk];
            UB.blockIdx = blockIdx;


            PMloc.RowBlockIdx(ib).push_back(UB.blockIdx);
            PMloc.ColBlockIdx( LBj ).push_back( bnum );



            UB.numRow = super->superPtr[bnum+1] - super->superPtr[bnum];
            cnt += UB_DESCRIPTOR;
            //            for( Int j = FirstBlockCol( UB.blockIdx, super ); 
            //                j < FirstBlockCol( UB.blockIdx+1, super ); j++ ){
            //              Int firstRow = index[cnt++];
            //              if( firstRow != FirstBlockCol( bnum+1, super ) )
            //                cols.push_back(j);
            //            }
            //            // Rewind the index
            //            cnt -= super->superPtr[UB.blockIdx+1] - super->superPtr[UB.blockIdx];

            int pos = 0;
            for( Int j = FirstBlockCol( UB.blockIdx, super ); 
                j < FirstBlockCol( UB.blockIdx+1, super ); j++ ){
              Int firstRow = index[cnt++];
              if( firstRow != FirstBlockCol( bnum+1, super ) )
                pos++;
            }
            // Rewind the index
            cnt -= super->superPtr[UB.blockIdx+1] - super->superPtr[UB.blockIdx];
            UB.cols.Resize(pos);
            UB.numCol = pos;


            pos = 0;
            for( Int j = FirstBlockCol( UB.blockIdx, super ); 
                j < FirstBlockCol( UB.blockIdx+1, super ); j++ ){
              Int firstRow = index[cnt++];
              if( firstRow != FirstBlockCol( bnum+1, super ) )
                UB.cols[pos++] = j;
            }
            // Rewind the index
            cnt -= super->superPtr[UB.blockIdx+1] - super->superPtr[UB.blockIdx];


            //            UB.numCol = cols.size();
            //            UB.cols.Resize(cols.size());
            //            std::copy(&cols[0],&cols[0]+cols.size(),&UB.cols[0]);
            //            //UB.cols   = IntNumVec( cols.size(), true, &cols[0] );
            UB.nzval.Resize( UB.numRow, UB.numCol );
            SetValue( UB.nzval, ZERO<Complex>() );

            Int cntcol = 0;
            for( Int j = 0; 
                j < super->superPtr[UB.blockIdx+1] - super->superPtr[UB.blockIdx]; j++ ){
              Int firstRow = index[cnt++];
              if( firstRow != FirstBlockCol( bnum+1, super ) ){
                Int tnrow = FirstBlockCol( bnum+1, super ) - firstRow;
                lapack::Lacpy( 'A', tnrow, 1, &pval[cntval], tnrow,
                    &UB.nzval(firstRow - FirstBlockCol(bnum, super), cntcol),
                    UB.numRow );
                cntcol ++;
                cntval += tnrow;
              }
            } // for( j )

#if ( _DEBUGlevel_ >= 1 )
            statusOFS 
              << "U part: bnum = " << bnum << ", numBlock = " << Urow.size()
              << ", blockIdx = " << UB.blockIdx
              << ", numRow = " << UB.numRow 
              << ", numCol = " << UB.numCol << std::endl;
#endif 

          } // for (jblk)



        } // if( index )

      } // for(ib)
#ifndef _RELEASE_
      PopCallStack();
#endif




      for( Int ib = 0; ib < PMloc.NumLocalBlockRow(); ib++ ){
        std::vector<Int> & rowBlockIdx = PMloc.RowBlockIdx(ib);
        std::sort(rowBlockIdx.begin(),rowBlockIdx.end());
      }

      for( Int jb = 0; jb < PMloc.NumLocalBlockCol(); jb++ ){
        std::vector<Int> & colBlockIdx = PMloc.ColBlockIdx(jb);
        std::sort(colBlockIdx.begin(),colBlockIdx.end());
      }





#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method ComplexSuperLUData::LUstructToPMatrix  ----- 



  void
    ComplexSuperLUData::SymbolicToSuperNode	( SuperNodeType& super )
    {
#ifndef _RELEASE_
      PushCallStack("ComplexSuperLUData::SymbolicToSuperNode");
#endif
      Int n = ptrData->A.ncol;
      // Permutation vector
      Int *perm_c = ptrData->ScalePermstruct.perm_c;
      super.perm.Resize( n );
      std::copy( perm_c, perm_c + n, super.perm.Data() );

      // Construct the inverse permutation vector
      super.permInv.Resize( n );
      for( Int i = 0; i < n; i++ ){
        super.permInv(i) = i;
      }

      std::sort( super.permInv.Data(), super.permInv.Data() + n,
          IndexComp<IntNumVec&>(super.perm) );

      // Row Permutation vector
      Int *perm_r = ptrData->ScalePermstruct.perm_r;
      super.perm_r.Resize( n );
      std::copy( perm_r, perm_r + n, super.perm_r.Data() );

      // Construct the inverse row permutation vector
      super.permInv_r.Resize( n );
      for( Int i = 0; i < n; i++ ){
        super.permInv_r(i) = i;
      }

      std::sort( super.permInv_r.Data(), super.permInv_r.Data() + n,
          IndexComp<IntNumVec&>(super.perm_r) );



      // Supernodal information

      Int *xsup    = ptrData -> LUstruct.Glu_persist -> xsup;
      Int *superno = ptrData -> LUstruct.Glu_persist -> supno;
      Int *etree = ptrData -> LUstruct.etree;

      Int numSuper = superno[ n-1 ] + 1;
      super.superPtr.Resize( numSuper + 1 );
      std::copy( xsup, xsup + numSuper + 1, super.superPtr.Data() );
      super.superIdx.Resize( n );
      std::copy( superno, superno + n, super.superIdx.Data() );
      super.etree.Resize( n );
      std::copy( etree, etree + n, super.etree.Data() );

#ifndef _RELEASE_
      PopCallStack();
#endif

      return ;
    } 		// -----  end of method ComplexSuperLUData::SymbolicToSuperNode  ----- 

}
