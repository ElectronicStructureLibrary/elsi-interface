/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

   Authors: Mathias Jacquelin and Lin Lin
 
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
/// @file SuperLUMatrix_impl.hpp
/// @brief Implementation of the wrapper class for SuperLU internal data structures.
/// @date 2014-03-17
#ifndef _PEXSI_SUPERLUMATRIX_IMPL_HPP_
#define _PEXSI_SUPERLUMATRIX_IMPL_HPP_

namespace PEXSI{

inline SuperLUMatrix<Real>::SuperLUMatrix	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::SuperLUMatrix");
#endif
  ptrData = NULL;
#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix<Real>::SuperLUMatrix  ----- 


inline SuperLUMatrix<Real>::SuperLUMatrix	( const SuperLUGrid<Real>& g, const SuperLUOptions& opt )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::SuperLUMatrix");
#endif
    ptrData = new RealSuperLUData(g,opt);
#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix<Real>::SuperLUMatrix  ----- 

inline SuperLUMatrix<Real>::~SuperLUMatrix	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::~SuperLUMatrix");
#endif
  if( ptrData != NULL )
    delete ptrData;
#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix<Real>::~SuperLUMatrix  ----- 




inline SuperLUMatrix<Real>::SuperLUMatrix(const SuperLUMatrix<Real> & g){
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::SuperLUMatrix");
#endif
    if(g.ptrData==NULL){
      ptrData=NULL;
    }
    else{
      ptrData = new RealSuperLUData(*g.ptrData);
    }
#ifndef _RELEASE_
	PopCallStack();
#endif
}

inline SuperLUMatrix<Real> & SuperLUMatrix<Real>::operator = (const SuperLUMatrix<Real> & g){
    
  if(this!=&g){
    if(ptrData!=NULL){
      delete ptrData;
    }

    if(g.ptrData==NULL){
      ptrData=NULL;
    }
    else{
      ptrData = new RealSuperLUData(*g.ptrData);
    }
  }
    return *this;
}






inline void 
SuperLUMatrix<Real>::Setup ( const SuperLUGrid<Real>& g, const SuperLUOptions& opt )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::Setup");
#endif
//  if( ptrData == NULL ){
    ptrData = new RealSuperLUData(g,opt);
//  }
//  else{
//    #ifdef USE_ABORT
//abort();
//#endif
//throw std::logic_error("SuperLUMatrix has been set up before.");
//  }
#ifndef _RELEASE_
	PopCallStack();
#endif
  return;
} 		// -----  end of method SuperLUMatrix<Real>::Setup  ----- 



inline Int SuperLUMatrix<Real>::m (  ) const	
{
	return ptrData->m();
} 		// -----  end of method SuperLUMatrix<Real>::m  ----- 



inline Int SuperLUMatrix<Real>::n (  ) const	
{
	return ptrData->n();
} 		// -----  end of method SuperLUMatrix<Real>::n  ----- 

inline void
SuperLUMatrix<Real>::DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Real>& sparseA , const SuperLUOptions& opt )
{
#ifndef _RELEASE_
	PushCallStack( "SuperLUMatrix<Real>::DistSparseMatrixToSuperMatrixNRloc" );
#endif
  ptrData->DistSparseMatrixToSuperMatrixNRloc(sparseA, opt );
#ifndef _RELEASE_
	PopCallStack();
#endif
	return;
} 		// -----  end of method SuperLUMatrix<Real>::DistSparseMatrixToSuperMatrixNRloc ----- 


inline void
SuperLUMatrix<Real>::DestroyAOnly	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::DestroyAOnly");
#endif
  ptrData->DestroyAOnly();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::DestroyAOnly  ----- 

inline void
SuperLUMatrix<Real>::SymbolicFactorize	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::SymbolicFactorize");
#endif
  ptrData->SymbolicFactorize();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::SymbolicFactorize  ----- 


inline void
SuperLUMatrix<Real>::Distribute	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::Distribute");
#endif
  ptrData->Distribute();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::Distribute  ----- 


inline void
SuperLUMatrix<Real>::NumericalFactorize	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::NumericalFactorize");
#endif
  ptrData->NumericalFactorize();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::NumericalFactorize  ----- 


inline void
SuperLUMatrix<Real>::ConvertNRlocToNC	( SuperLUMatrix& AGlobal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::ConvertNRlocToNC");
#endif
  ptrData->ConvertNRlocToNC(AGlobal.ptrData);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::ConvertNRlocToNC  ----- 

inline void
SuperLUMatrix<Real>::MultiplyGlobalMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& bGlobal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::MultiplyGlobalMultiVector");
#endif
  ptrData->MultiplyGlobalMultiVector(xGlobal, bGlobal);
#ifndef _RELEASE_
	PopCallStack();
#endif
 
	return ;
} 		// -----  end of method SuperLUMatrix<Real>::MultiplyGlobalMultiVector  ----- 


inline void
SuperLUMatrix<Real>::DistributeGlobalMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& xLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::DistributeGlobalMultiVector");
#endif
  ptrData->DistributeGlobalMultiVector(xGlobal, xLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::DistributeGlobalMultiVector  ----- 


inline void SuperLUMatrix<Real>::GatherDistributedMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& xLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::GatherDistributedMultiVector");
#endif
  ptrData->GatherDistributedMultiVector(xGlobal, xLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::GatherDistributedMultiVector  ----- 


inline void
SuperLUMatrix<Real>::SolveDistMultiVector	( NumMat<Real>& bLocal, DblNumVec& berr )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::SolveDistMultiVector");
#endif
  ptrData->SolveDistMultiVector(bLocal, berr );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::SolveDistMultiVector  ----- 


inline void
SuperLUMatrix<Real>::CheckErrorDistMultiVector	( NumMat<Real>& xLocal, NumMat<Real>& xTrueLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::CheckErrorDistMultiVector");
#endif
  ptrData->CheckErrorDistMultiVector(xLocal, xTrueLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::CheckErrorDistMultiVector  ----- 


inline void
SuperLUMatrix<Real>::LUstructToPMatrix	( PMatrix<Real>& PMloc )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::LUstructToPMatrix");
#endif
  ptrData->LUstructToPMatrix(PMloc);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::LUstructToPMatrix  ----- 



inline void
SuperLUMatrix<Real>::SymbolicToSuperNode	( SuperNodeType& super )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::SymbolicToSuperNode");
#endif
  ptrData->SymbolicToSuperNode(super);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::SymbolicToSuperNode  ----- 

}

namespace PEXSI{

inline SuperLUMatrix<Complex>::SuperLUMatrix	( )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::SuperLUMatrix");
#endif
  ptrData = NULL;
#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix<Complex>::SuperLUMatrix  ----- 


inline SuperLUMatrix<Complex>::SuperLUMatrix	( const SuperLUGrid<Complex>& g, const SuperLUOptions& opt )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::SuperLUMatrix");
#endif
    ptrData = new ComplexSuperLUData(g,opt);
#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix<Complex>::SuperLUMatrix  ----- 

inline SuperLUMatrix<Complex>::~SuperLUMatrix	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::~SuperLUMatrix");
#endif
  if( ptrData != NULL ){
    delete ptrData;
  }
#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix<Complex>::~SuperLUMatrix  ----- 


inline SuperLUMatrix<Complex>::SuperLUMatrix(const SuperLUMatrix<Complex> & g){
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::SuperLUMatrix");
#endif
    if(g.ptrData==NULL){
      ptrData=NULL;
    }
    else{
      ptrData = new ComplexSuperLUData(*g.ptrData);
    }
#ifndef _RELEASE_
	PopCallStack();
#endif
}

inline SuperLUMatrix<Complex> & SuperLUMatrix<Complex>::operator = (const SuperLUMatrix<Complex> & g){
    
  if(this!=&g){
    if(ptrData!=NULL){
      delete ptrData;
    }

    if(g.ptrData==NULL){
      ptrData=NULL;
    }
    else{
      ptrData = new ComplexSuperLUData(*g.ptrData);
    }
  }
    return *this;
}





inline void
SuperLUMatrix<Complex>::Setup ( const SuperLUGrid<Complex>& g, const SuperLUOptions& opt )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::Setup");
#endif
//  if( ptrData == NULL ){
    ptrData = new ComplexSuperLUData(g,opt);
//  }
//  else{
//    #ifdef USE_ABORT
//abort();
//#endif
//throw std::logic_error("SuperLUMatrix has been set up before.");
//  }
#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix<Complex>::Setup  ----- 


inline Int SuperLUMatrix<Complex>::m (  ) const	
{
return ptrData->m();
}		// -----  end of method SuperLUMatrix<Complex>::m  ----- 

inline Int SuperLUMatrix<Complex>::n (  ) const	
{
	return ptrData->n();
} 		// -----  end of method SuperLUMatrix<Complex>::n  ----- 

inline void
SuperLUMatrix<Complex>::DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Complex>& sparseA , const SuperLUOptions& opt)
{
#ifndef _RELEASE_
	PushCallStack( "SuperLUMatrix<Complex>::DistSparseMatrixToSuperMatrixNRloc" );
#endif
  ptrData->DistSparseMatrixToSuperMatrixNRloc(sparseA, opt);
#ifndef _RELEASE_
	PopCallStack();
#endif
	return;
} 		// -----  end of method SuperLUMatrix<Complex>::DistSparseMatrixToSuperMatrixNRloc ----- 

inline void
SuperLUMatrix<Complex>::DestroyAOnly	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::DestroyAOnly");
#endif
  ptrData->DestroyAOnly();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::DestroyAOnly  ----- 

inline void
SuperLUMatrix<Complex>::SymbolicFactorize	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::SymbolicFactorize");
#endif
  ptrData->SymbolicFactorize();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::SymbolicFactorize  ----- 

inline void
SuperLUMatrix<Complex>::Distribute	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::Distribute");
#endif
  ptrData->Distribute();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::Distribute  ----- 

inline void
SuperLUMatrix<Complex>::NumericalFactorize	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::NumericalFactorize");
#endif
  ptrData->NumericalFactorize();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::NumericalFactorize  ----- 

inline void
SuperLUMatrix<Complex>::ConvertNRlocToNC	( SuperLUMatrix& AGlobal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::ConvertNRlocToNC");
#endif
  ptrData->ConvertNRlocToNC(AGlobal.ptrData);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::ConvertNRlocToNC  ----- 

inline void
SuperLUMatrix<Complex>::MultiplyGlobalMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& bGlobal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::MultiplyGlobalMultiVector");
#endif
  ptrData->MultiplyGlobalMultiVector(xGlobal, bGlobal);
#ifndef _RELEASE_
	PopCallStack();
#endif
 
	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::MultiplyGlobalMultiVector  ----- 

inline void
SuperLUMatrix<Complex>::DistributeGlobalMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::DistributeGlobalMultiVector");
#endif
  ptrData->DistributeGlobalMultiVector(xGlobal, xLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::DistributeGlobalMultiVector  ----- 

inline void SuperLUMatrix<Complex>::GatherDistributedMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::GatherDistributedMultiVector");
#endif
  ptrData->GatherDistributedMultiVector(xGlobal, xLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::GatherDistributedMultiVector  ----- 

inline void
SuperLUMatrix<Complex>::SolveDistMultiVector	( NumMat<Complex>& bLocal, DblNumVec& berr )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::SolveDistMultiVector");
#endif
  ptrData->SolveDistMultiVector(bLocal, berr );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::SolveDistMultiVector  ----- 

inline void
SuperLUMatrix<Complex>::CheckErrorDistMultiVector	( NumMat<Complex>& xLocal, NumMat<Complex>& xTrueLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::CheckErrorDistMultiVector");
#endif
  ptrData->CheckErrorDistMultiVector(xLocal, xTrueLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::CheckErrorDistMultiVector  ----- 

inline void
SuperLUMatrix<Complex>::LUstructToPMatrix	( PMatrix<Complex>& PMloc )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::LUstructToPMatrix");
#endif
  ptrData->LUstructToPMatrix(PMloc);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::LUstructToPMatrix  ----- 

inline void
SuperLUMatrix<Complex>::SymbolicToSuperNode	( SuperNodeType& super )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::SymbolicToSuperNode");
#endif
  ptrData->SymbolicToSuperNode(super);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::SymbolicToSuperNode  ----- 

}

#endif //_PEXSI_SUPERLUMATRIX_IMPL_HPP_
