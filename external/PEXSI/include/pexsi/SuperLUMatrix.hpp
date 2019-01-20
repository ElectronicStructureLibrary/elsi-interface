/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.

Author: Mathias Jacquelin and Lin Lin

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
/// @file SuperLUMatrix.hpp
/// @brief Wrapper class for SuperLU internal data structures.
/// @date 2012-11-12
#ifndef _PEXSI_SUPERLUMATRIX_HPP_
#define _PEXSI_SUPERLUMATRIX_HPP_

// Interface with PSelInv
#include "pexsi/pselinv.hpp"

// Interface with sparse matrix (CSC format)
#include "pexsi/sparse_matrix.hpp"

// Interface with LAPACK
#include "pexsi/lapack.hpp"


#include "pexsi/superlu_dist_internal.hpp"
#include "pexsi/SuperLUGrid.hpp"
//#include "pexsi/SuperLUData.hpp"

namespace PEXSI{


struct SuperNodeType;
template<typename T>  class PMatrix;



/// @class SuperLUMatrix
/// @brief An thin interface to keep the main code insulated from
/// the source code of SuperLU.
///
template<typename T> class SuperLUMatrix{
private:
  /// @struct SuperLUData
  ///
  /// @brief  Data of a matrix in the SuperLU format.
  ///
  /// SuperLUData is only used to define ptrData which is a private member
  /// SuperLUMatrix, and is only defined in superlu_dist_interf.cpp.
  ///
  /// @brief Actual pointer to save the data related to SuperLU.
  void *  ptrData;


public:

  SuperLUMatrix( const SuperLUGrid<T>& g, const SuperLUOptions& opt = SuperLUOptions() ){ }
  ~SuperLUMatrix(){ }
  SuperLUMatrix(const SuperLUMatrix & g){};
  SuperLUMatrix & operator = (const SuperLUMatrix & g){};

  Int m() const {return -1;}

  Int n() const { return -1;}

  /// @brief Convert a distributed sparse matrix in compressed sparse
  /// column format into the SuperLU compressed row format.
  ///
  /// The output is saved in the current %SuperLUMatrix.
  ///
  /// @note
  /// Although LU factorization is used, the routine
  /// assumes that the matrix is strictly symmetric, and therefore the
  /// compressed sparse row (CSR) format, used by SuperLU_DIST, gives
  /// exactly the same matrix as formed by the compresed sparse column
  /// format (CSC).
  ///
  ///
  /// @param[in] sparseA Matrix A saved in distributed compressed
  /// sparse column format.
  ///
  /// @todo
  /// Better way to incorporate both real and complex arithmetic.
  void DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<T>& sparseA , const SuperLUOptions & options){ }

  /// @brief Releases the data in A but keeps other data,
  /// such as LUstruct.
  ///
  /// This allows one to perform factorization of
  /// matrices of the same pattern, such as the option
  ///
  /// fact = SamePattern_SameRowPerm
  ///
  /// in SuperLU_DIST.
  void DestroyAOnly(){ }

  /// @brief Factorizes the superlu matrix symbolically.
  ///
  /// Symbolic factorization contains three steps.
  ///
  /// - Permute the matrix to reduce fill-in.
  /// - Symbolic factorize the matrix.
  /// - Distribute the matrix into 2D block cyclic format.
  ///
  /// This routine is controlled via SuperLUOptions.
  void SymbolicFactorize(){ }

  /// @brief Distribute redistrbutes the SuperMatrix in parallel so that
  /// it is ready for the numerical factorization.
  ///
  /// @todo
  /// Better way to incorporate both real and complex arithmetic.
  void Distribute(){ }

  /// @brief Performs LU factorization numerically.
  ///
  /// The %SuperLUMatrix should have been permuted and distributed.
  ///
  /// @todo
  /// Better way to incorporate both real and complex arithmetic.
  void NumericalFactorize(){ }

  /// @brief ConvertNRlocToNC converts a distributed compressed sparse
  /// row matrix to a global compressed sparse column matrix.
  ///
  /// @param[out] AGlobal
  void ConvertNRlocToNC( SuperLUMatrix<T>& AGlobal ){ }

  /// @brief MultiplyGlobalMultiVector computes b = A * x.
  ///
  /// Both xGlobal and bGlobal are saved globally.
  ///
  /// @param[in] xGlobal
  /// @param[out] bGlobal
  void MultiplyGlobalMultiVector( NumMat<T>& xGlobal, NumMat<T>& bGlobal ){ }

  /// @brief DistributeGlobalMultiVector distributes a global multivector into a
  /// local multivector according to the compressed row format of A.
  ///
  /// @param[in] xGlobal
  /// @param[out] xLocal
  void DistributeGlobalMultiVector( NumMat<T>& xGlobal, NumMat<T>& xLocal ){ }

  void GatherDistributedMultiVector	( NumMat<T>& xGlobal, NumMat<T>& xLocal ){ }

  /// @brief Solve A x = b with b overwritten by x for
  /// distributed multivector.
  ///
  /// @param[in,out] bLocal Right hand side savd in the distributed format.
  /// @param[out] berr The componentwise relative backward error of each solution vector.
  void SolveDistMultiVector( NumMat<T>& bLocal, DblNumVec& berr ){ }


  /// @brief Prints out the error by direct
  /// comparison with the true solution in distributed format.
  ///
  /// @param[in] xLocal The computed solution.
  /// @param[in] xTrueLocal The true solution.
  void CheckErrorDistMultiVector( NumMat<T>& xLocal, NumMat<T>& xTrueLocal ){
  }
  /// @brief LUstructToPMatrix converts the data in LUstruct to PMatrix.
  ///
  /// @paramp[out] PMloc
  void LUstructToPMatrix( PMatrix<T>& PMloc ){ }

  /// @brief SymbolicToSuperNode converts the symbolic information to
  /// SuperNode structure in SelInv.
  ///
  /// @param[out] super
  void SymbolicToSuperNode( SuperNodeType& super ){ }
};

template< > class SuperLUMatrix<Real>{
private:
  /// @struct SuperLUData
  ///
  /// @brief  Data of a matrix in the SuperLU format.
  ///
  /// SuperLUData is only used to define ptrData which is a private member
  /// SuperLUMatrix, and is only defined in superlu_dist_interf.cpp.
  ///
  /// @brief Actual pointer to save the data related to SuperLU.
  RealSuperLUData*  ptrData;

public:

  SuperLUMatrix();
  SuperLUMatrix( const SuperLUGrid<Real>& g, const SuperLUOptions& opt = SuperLUOptions() );
  ~SuperLUMatrix();
  SuperLUMatrix(const SuperLUMatrix & g);
  SuperLUMatrix & operator = (const SuperLUMatrix & g);

  void Setup( const SuperLUGrid<Real>& g, const SuperLUOptions& opt = SuperLUOptions() );

  Int m() const;
  Int n() const;

  /// @brief Convert a distributed sparse matrix in compressed sparse
  /// column format into the SuperLU compressed row format.
  ///
  /// The output is saved in the current %SuperLUMatrix.
  ///
  /// @note
  /// Although LU factorization is used, the routine
  /// assumes that the matrix is strictly symmetric, and therefore the
  /// compressed sparse row (CSR) format, used by SuperLU_DIST, gives
  /// exactly the same matrix as formed by the compresed sparse column
  /// format (CSC).
  ///
  ///
  /// @param[in] sparseA Matrix A saved in distributed compressed
  /// sparse column format.
  ///
  /// @todo
  /// Better way to incorporate both real and complex arithmetic.
  void DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Real>& sparseA , const SuperLUOptions & opt);

  /// @brief Releases the data in A but keeps other data,
  /// such as LUstruct.
  ///
  /// This allows one to perform factorization of
  /// matrices of the same pattern, such as the option
  ///
  /// fact = SamePattern_SameRowPerm
  ///
  /// in SuperLU_DIST.
  void DestroyAOnly();

  /// @brief Factorizes the superlu matrix symbolically.
  ///
  /// Symbolic factorization contains three steps.
  ///
  /// - Permute the matrix to reduce fill-in.
  /// - Symbolic factorize the matrix.
  /// - Distribute the matrix into 2D block cyclic format.
  ///
  /// This routine is controlled via SuperLUOptions.
  void SymbolicFactorize();

  /// @brief Distribute redistrbutes the SuperMatrix in parallel so that
  /// it is ready for the numerical factorization.
  ///
  /// @todo
  /// Better way to incorporate both real and complex arithmetic.
  void Distribute();

  /// @brief Performs LU factorization numerically.
  ///
  /// The %SuperLUMatrix should have been permuted and distributed.
  ///
  /// @todo
  /// Better way to incorporate both real and complex arithmetic.
  void NumericalFactorize();

  /// @brief ConvertNRlocToNC converts a distributed compressed sparse
  /// row matrix to a global compressed sparse column matrix.
  ///
  /// @param[out] AGlobal
  void ConvertNRlocToNC( SuperLUMatrix<Real>& AGlobal );

  /// @brief MultiplyGlobalMultiVector computes b = A * x.
  ///
  /// Both xGlobal and bGlobal are saved globally.
  ///
  /// @param[in] xGlobal
  /// @param[out] bGlobal
  void MultiplyGlobalMultiVector( NumMat<Real>& xGlobal, NumMat<Real>& bGlobal );

  /// @brief DistributeGlobalMultiVector distributes a global multivector into a
  /// local multivector according to the compressed row format of A.
  ///
  /// @param[in] xGlobal
  /// @param[out] xLocal
  void DistributeGlobalMultiVector( NumMat<Real>& xGlobal, NumMat<Real>& xLocal );

  void GatherDistributedMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& xLocal );

  /// @brief Solve A x = b with b overwritten by x for
  /// distributed multivector.
  ///
  /// @param[in,out] bLocal Right hand side savd in the distributed format.
  /// @param[out] berr The componentwise relative backward error of each solution vector.
  void SolveDistMultiVector( NumMat<Real>& bLocal, DblNumVec& berr );

  /// @brief Prints out the error by direct
  /// comparison with the true solution in distributed format.
  ///
  /// @param[in] xLocal The computed solution.
  /// @param[in] xTrueLocal The true solution.
  void CheckErrorDistMultiVector( NumMat<Real>& xLocal, NumMat<Real>& xTrueLocal );

  /// @brief LUstructToPMatrix converts the data in LUstruct to PMatrix.
  ///
  /// @paramp[out] PMloc
  void LUstructToPMatrix( PMatrix<Real>& PMloc );

  /// @brief SymbolicToSuperNode converts the symbolic information to
  /// SuperNode structure in SelInv.
  ///
  /// @param[out] super
  void SymbolicToSuperNode( SuperNodeType& super );
};

template< > class SuperLUMatrix<Complex>{
private:
  /// @struct SuperLUData
  ///
  /// @brief  Data of a matrix in the SuperLU format.
  ///
  /// SuperLUData is only used to define ptrData which is a private member
  /// SuperLUMatrix, and is only defined in superlu_dist_interf.cpp.
  ///
  /// @brief Actual pointer to save the data related to SuperLU.
  ComplexSuperLUData*  ptrData;


public:

  SuperLUMatrix();
  SuperLUMatrix( const SuperLUGrid<Complex>& g, const SuperLUOptions& opt = SuperLUOptions() );
  ~SuperLUMatrix();
  SuperLUMatrix(const SuperLUMatrix & g);
  SuperLUMatrix & operator = (const SuperLUMatrix & g);

  void Setup( const SuperLUGrid<Complex>& g, const SuperLUOptions& opt = SuperLUOptions() );

  Int m() const;
  Int n() const;

  /// @brief Convert a distributed sparse matrix in compressed sparse
  /// column format into the SuperLU compressed row format.
  ///
  /// The output is saved in the current %SuperLUMatrix.
  ///
  /// @note
  /// Although LU factorization is used, the routine
  /// assumes that the matrix is strictly symmetric, and therefore the
  /// compressed sparse row (CSR) format, used by SuperLU_DIST, gives
  /// exactly the same matrix as formed by the compresed sparse column
  /// format (CSC).
  ///
  ///
  /// @param[in] sparseA Matrix A saved in distributed compressed
  /// sparse column format.
  ///
  /// @todo
  /// Better way to incorporate both real and complex arithmetic.
  void DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Complex>& sparseA , const SuperLUOptions & options);

  /// @brief Releases the data in A but keeps other data,
  /// such as LUstruct.
  ///
  /// This allows one to perform factorization of
  /// matrices of the same pattern, such as the option
  ///
  /// fact = SamePattern_SameRowPerm
  ///
  /// in SuperLU_DIST.
  void DestroyAOnly();

  /// @brief Factorizes the superlu matrix symbolically.
  ///
  /// Symbolic factorization contains three steps.
  ///
  /// - Permute the matrix to reduce fill-in.
  /// - Symbolic factorize the matrix.
  /// - Distribute the matrix into 2D block cyclic format.
  ///
  /// This routine is controlled via SuperLUOptions.
  void SymbolicFactorize();

  /// @brief Distribute redistrbutes the SuperMatrix in parallel so that
  /// it is ready for the numerical factorization.
  ///
  /// @todo
  /// Better way to incorporate both real and complex arithmetic.
  void Distribute();

  /// @brief Performs LU factorization numerically.
  ///
  /// The %SuperLUMatrix should have been permuted and distributed.
  ///
  /// @todo
  /// Better way to incorporate both real and complex arithmetic.
  void NumericalFactorize();

  /// @brief ConvertNRlocToNC converts a distributed compressed sparse
  /// row matrix to a global compressed sparse column matrix.
  ///
  /// @param[out] AGlobal
  void ConvertNRlocToNC( SuperLUMatrix<Complex>& AGlobal );

  /// @brief MultiplyGlobalMultiVector computes b = A * x.
  ///
  /// Both xGlobal and bGlobal are saved globally.
  ///
  /// @param[in] xGlobal
  /// @param[out] bGlobal
  void MultiplyGlobalMultiVector( NumMat<Complex>& xGlobal, NumMat<Complex>& bGlobal );

  /// @brief DistributeGlobalMultiVector distributes a global multivector into a
  /// local multivector according to the compressed row format of A.
  ///
  /// @param[in] xGlobal
  /// @param[out] xLocal
  void DistributeGlobalMultiVector( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal );

  void GatherDistributedMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal );

  /// @brief Solve A x = b with b overwritten by x for
  /// distributed multivector.
  ///
  /// @param[in,out] bLocal Right hand side savd in the distributed format.
  /// @param[out] berr The componentwise relative backward error of each solution vector.
  void SolveDistMultiVector( NumMat<Complex>& bLocal, DblNumVec& berr );


  /// @brief Prints out the error by direct
  /// comparison with the true solution in distributed format.
  ///
  /// @param[in] xLocal The computed solution.
  /// @param[in] xTrueLocal The true solution.
  void CheckErrorDistMultiVector( NumMat<Complex>& xLocal, NumMat<Complex>& xTrueLocal );

  /// @brief LUstructToPMatrix converts the data in LUstruct to PMatrix.
  ///
  /// @paramp[out] PMloc
  void LUstructToPMatrix( PMatrix<Complex>& PMloc );

  /// @brief SymbolicToSuperNode converts the symbolic information to
  /// SuperNode structure in SelInv.
  ///
  /// @param[out] super
  void SymbolicToSuperNode( SuperNodeType& super );
};

} // namespace PEXSI

#include "pexsi/SuperLUMatrix_impl.hpp"

#endif // _PEXSI_SUPERLUMATRIX_HPP_

