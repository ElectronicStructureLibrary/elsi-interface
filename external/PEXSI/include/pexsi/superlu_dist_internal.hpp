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
/// @file superlu_dist_internal.hpp
/// @brief Internal structures for interfacing with SuperLU_Dist (version 3.0 and later)
/// @date 2014-03-17

#ifndef _PEXSI_SUPERLU_INTERNAL_HPP_
#define _PEXSI_SUPERLU_INTERNAL_HPP_

#include "pexsi/environment.hpp"
#include <string>
#include "pexsi/sparse_matrix.hpp"
#include "pexsi/NumMat.hpp"
#include "pexsi/NumVec.hpp"
#include "pexsi/SuperLUGrid.hpp"

namespace PEXSI{
/// @struct SuperLUOptions
/// @brief A thin interface for passing parameters to set the SuperLU
/// options.
///
struct SuperLUOptions{
  /// @brief Number of processors for parallel symbolic factorization.
  ///
  /// numProcSymbFact should be a power of 2, and is only useful when
  ///
  /// ColPerm = "PARMETIS"
  Int              numProcSymbFact;

  /// @brief Option of matrixi column permutation strategy.
  ///
  /// The following options of column permutation strategy are available (case
  /// sensitive):
  ///
  /// - "MMD_AT_PLUS_A": Multiple minimum degree ordering. This is the
  /// default option.
  /// - "METIS_AT_PLUS_A": Sequential ordering using METIS. This
  /// requires the usage of METIS package.
  /// - "PARMETIS": Parallel ordering. This requires the usage of
  /// ParMETIS/PT-SCOTCH package.
  /// - "NATURAL": No ordering. This lead to SIGNIFICANTLY higher
  /// computational and storage costs.
  ///
  std::string      ColPerm;

  /// @brief Option of matrix row permutation strategy.
  ///
  /// The following options of row permutation strategy are available (case
  /// sensitive):
  ///
  /// - "LargeDiag": Make the diagonal large relative to off-diagonal elements (MC64).
  /// - "NOROWPERM": No row permutation. This might lead to numerically unstable
  /// factorization, and selected inversion.
  ///
  std::string      RowPerm;

  /// @brief Option whether to equilibrate the system.
  ///
  /// The following options of equilibration strategy are available (case
  /// sensitive):
  ///
  /// - "YES": Scale A's rows and columns to have unit norm.
  /// - "NO": Don't equilibrate the system.
  ///
  std::string      Equil;

  /// @brief Option to specify if matrix is symmetric or not.
  Int              Symmetric;

  /// @brief Option to specify whether selected elements should
  /// be computed in the pattern of the transposed matrix or not.
  /// Note that this has an impact only for unsymmetric matrices.
  Int              Transpose;


  // Member functions to setup the default value
  SuperLUOptions(): numProcSymbFact(0), ColPerm("MMD_AT_PLUS_A"), Symmetric(1), RowPerm("NOROWPERM"), Transpose(0), Equil("NO") {}
};


struct SuperNodeType;
template<typename T> class PMatrix;

class ComplexSuperLUData_internal;
class RealSuperLUData_internal;

class RealSuperLUData{
protected:
  RealSuperLUData_internal * ptrData;
public:
  RealSuperLUData( const SuperLUGrid<Real>& g, const SuperLUOptions& opt );
  ~RealSuperLUData();
  RealSuperLUData(const RealSuperLUData & g);
  RealSuperLUData & operator = (const RealSuperLUData & g);

  Int m() const;
  Int n() const;
  void DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Real>& sparseA , const SuperLUOptions& opt);
  void DestroyAOnly();
  void SymbolicFactorize();
  void Distribute();
  void NumericalFactorize();
  void ConvertNRlocToNC	( RealSuperLUData * aptrData );
  void MultiplyGlobalMultiVector( NumMat<Real>& xGlobal, NumMat<Real>& bGlobal );
  void DistributeGlobalMultiVector( NumMat<Real>& xGlobal, NumMat<Real>& xLocal );
  void GatherDistributedMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& xLocal );
  void SolveDistMultiVector( NumMat<Real>& bLocal, DblNumVec& berr );
  void CheckErrorDistMultiVector( NumMat<Real>& xLocal, NumMat<Real>& xTrueLocal );
  void LUstructToPMatrix( PMatrix<Real>& PMloc );
  void SymbolicToSuperNode( SuperNodeType& super );
};


class ComplexSuperLUData{
protected:
  ComplexSuperLUData_internal * ptrData;
public:
  ComplexSuperLUData( const SuperLUGrid<Complex>& g, const SuperLUOptions& opt );
  ~ComplexSuperLUData();
  ComplexSuperLUData(const ComplexSuperLUData & g);
  ComplexSuperLUData & operator = (const ComplexSuperLUData & g);

  Int m() const;
  Int n() const;
  void DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Complex>& sparseA , const SuperLUOptions& opt);
  void DestroyAOnly();
  void SymbolicFactorize();
  void Distribute();
  void NumericalFactorize();
  void ConvertNRlocToNC	( ComplexSuperLUData * aptrData );
  void MultiplyGlobalMultiVector( NumMat<Complex>& xGlobal, NumMat<Complex>& bGlobal );
  void DistributeGlobalMultiVector( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal );
  void GatherDistributedMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal );
  void SolveDistMultiVector( NumMat<Complex>& bLocal, DblNumVec& berr );
  void CheckErrorDistMultiVector( NumMat<Complex>& xLocal, NumMat<Complex>& xTrueLocal );
  void LUstructToPMatrix( PMatrix<Complex>& PMloc );
  void SymbolicToSuperNode( SuperNodeType& super );
};


}



#endif //_PEXSI_SUPERLU_INTERNAL_HPP_
