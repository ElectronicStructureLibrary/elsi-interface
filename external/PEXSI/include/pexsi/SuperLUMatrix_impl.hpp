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
  ptrData = NULL;
} 		// -----  end of method SuperLUMatrix<Real>::SuperLUMatrix  -----


inline SuperLUMatrix<Real>::SuperLUMatrix	( const SuperLUGrid<Real>& g, const SuperLUOptions& opt )
{
  ptrData = new RealSuperLUData(g,opt);
} 		// -----  end of method SuperLUMatrix<Real>::SuperLUMatrix  -----

inline SuperLUMatrix<Real>::~SuperLUMatrix	(  )
{
  if( ptrData != NULL )
    delete ptrData;
} 		// -----  end of method SuperLUMatrix<Real>::~SuperLUMatrix  -----




inline SuperLUMatrix<Real>::SuperLUMatrix(const SuperLUMatrix<Real> & g){
  if(g.ptrData==NULL){
    ptrData=NULL;
  }
  else{
    ptrData = new RealSuperLUData(*g.ptrData);
  }
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
  //  if( ptrData == NULL ){
  ptrData = new RealSuperLUData(g,opt);
  //  }
  //  else{
  //ErrorHandling("SuperLUMatrix has been set up before.");
  //  }
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
  ptrData->DistSparseMatrixToSuperMatrixNRloc(sparseA, opt );
  return;
} 		// -----  end of method SuperLUMatrix<Real>::DistSparseMatrixToSuperMatrixNRloc -----


inline void
SuperLUMatrix<Real>::DestroyAOnly	(  )
{
  ptrData->DestroyAOnly();

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::DestroyAOnly  -----

inline void
SuperLUMatrix<Real>::SymbolicFactorize	(  )
{
  ptrData->SymbolicFactorize();

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::SymbolicFactorize  -----


inline void
SuperLUMatrix<Real>::Distribute	(  )
{
  ptrData->Distribute();

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::Distribute  -----


inline void
SuperLUMatrix<Real>::NumericalFactorize	(  )
{
  ptrData->NumericalFactorize();

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::NumericalFactorize  -----


inline void
SuperLUMatrix<Real>::ConvertNRlocToNC	( SuperLUMatrix& AGlobal )
{
  ptrData->ConvertNRlocToNC(AGlobal.ptrData);

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::ConvertNRlocToNC  -----

inline void
SuperLUMatrix<Real>::MultiplyGlobalMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& bGlobal )
{
  ptrData->MultiplyGlobalMultiVector(xGlobal, bGlobal);

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::MultiplyGlobalMultiVector  -----


inline void
SuperLUMatrix<Real>::DistributeGlobalMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& xLocal )
{
  ptrData->DistributeGlobalMultiVector(xGlobal, xLocal );

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::DistributeGlobalMultiVector  -----


inline void SuperLUMatrix<Real>::GatherDistributedMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& xLocal )
{
  ptrData->GatherDistributedMultiVector(xGlobal, xLocal );

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::GatherDistributedMultiVector  -----


inline void
SuperLUMatrix<Real>::SolveDistMultiVector	( NumMat<Real>& bLocal, DblNumVec& berr )
{
  ptrData->SolveDistMultiVector(bLocal, berr );

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::SolveDistMultiVector  -----


inline void
SuperLUMatrix<Real>::CheckErrorDistMultiVector	( NumMat<Real>& xLocal, NumMat<Real>& xTrueLocal )
{
  ptrData->CheckErrorDistMultiVector(xLocal, xTrueLocal );

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::CheckErrorDistMultiVector  -----


inline void
SuperLUMatrix<Real>::LUstructToPMatrix	( PMatrix<Real>& PMloc )
{
  ptrData->LUstructToPMatrix(PMloc);

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::LUstructToPMatrix  -----



inline void
SuperLUMatrix<Real>::SymbolicToSuperNode	( SuperNodeType& super )
{
  ptrData->SymbolicToSuperNode(super);

  return ;
} 		// -----  end of method SuperLUMatrix<Real>::SymbolicToSuperNode  -----

}

namespace PEXSI{

inline SuperLUMatrix<Complex>::SuperLUMatrix	( )
{
  ptrData = NULL;
} 		// -----  end of method SuperLUMatrix<Complex>::SuperLUMatrix  -----


inline SuperLUMatrix<Complex>::SuperLUMatrix	( const SuperLUGrid<Complex>& g, const SuperLUOptions& opt )
{
  ptrData = new ComplexSuperLUData(g,opt);
} 		// -----  end of method SuperLUMatrix<Complex>::SuperLUMatrix  -----

inline SuperLUMatrix<Complex>::~SuperLUMatrix	(  )
{
  if( ptrData != NULL ){
    delete ptrData;
  }
} 		// -----  end of method SuperLUMatrix<Complex>::~SuperLUMatrix  -----


inline SuperLUMatrix<Complex>::SuperLUMatrix(const SuperLUMatrix<Complex> & g){
  if(g.ptrData==NULL){
    ptrData=NULL;
  }
  else{
    ptrData = new ComplexSuperLUData(*g.ptrData);
  }
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
  //  if( ptrData == NULL ){
  ptrData = new ComplexSuperLUData(g,opt);
  //  }
  //  else{
  //ErrorHandling("SuperLUMatrix has been set up before.");
  //  }
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
  ptrData->DistSparseMatrixToSuperMatrixNRloc(sparseA, opt);
  return;
} 		// -----  end of method SuperLUMatrix<Complex>::DistSparseMatrixToSuperMatrixNRloc -----

inline void
SuperLUMatrix<Complex>::DestroyAOnly	(  )
{
  ptrData->DestroyAOnly();

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::DestroyAOnly  -----

inline void
SuperLUMatrix<Complex>::SymbolicFactorize	(  )
{
  ptrData->SymbolicFactorize();

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::SymbolicFactorize  -----

inline void
SuperLUMatrix<Complex>::Distribute	(  )
{
  ptrData->Distribute();

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::Distribute  -----

inline void
SuperLUMatrix<Complex>::NumericalFactorize	(  )
{
  ptrData->NumericalFactorize();

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::NumericalFactorize  -----

inline void
SuperLUMatrix<Complex>::ConvertNRlocToNC	( SuperLUMatrix& AGlobal )
{
  ptrData->ConvertNRlocToNC(AGlobal.ptrData);

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::ConvertNRlocToNC  -----

inline void
SuperLUMatrix<Complex>::MultiplyGlobalMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& bGlobal )
{
  ptrData->MultiplyGlobalMultiVector(xGlobal, bGlobal);

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::MultiplyGlobalMultiVector  -----

inline void
SuperLUMatrix<Complex>::DistributeGlobalMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal )
{
  ptrData->DistributeGlobalMultiVector(xGlobal, xLocal );

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::DistributeGlobalMultiVector  -----

inline void SuperLUMatrix<Complex>::GatherDistributedMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal )
{
  ptrData->GatherDistributedMultiVector(xGlobal, xLocal );

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::GatherDistributedMultiVector  -----

inline void
SuperLUMatrix<Complex>::SolveDistMultiVector	( NumMat<Complex>& bLocal, DblNumVec& berr )
{
  ptrData->SolveDistMultiVector(bLocal, berr );

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::SolveDistMultiVector  -----

inline void
SuperLUMatrix<Complex>::CheckErrorDistMultiVector	( NumMat<Complex>& xLocal, NumMat<Complex>& xTrueLocal )
{
  ptrData->CheckErrorDistMultiVector(xLocal, xTrueLocal );

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::CheckErrorDistMultiVector  -----

inline void
SuperLUMatrix<Complex>::LUstructToPMatrix	( PMatrix<Complex>& PMloc )
{
  ptrData->LUstructToPMatrix(PMloc);

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::LUstructToPMatrix  -----

inline void
SuperLUMatrix<Complex>::SymbolicToSuperNode	( SuperNodeType& super )
{
  ptrData->SymbolicToSuperNode(super);

  return ;
} 		// -----  end of method SuperLUMatrix<Complex>::SymbolicToSuperNode  -----

}

#endif //_PEXSI_SUPERLUMATRIX_IMPL_HPP_
