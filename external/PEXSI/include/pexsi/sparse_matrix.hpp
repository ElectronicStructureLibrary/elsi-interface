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
/// @file sparse_matrix.hpp
/// @brief Sparse matrix and Distributed sparse matrix in compressed
/// column format.
/// @date 2012-11-10
#ifndef _PEXSI_SPARSE_MATRIX_HPP_
#define _PEXSI_SPARSE_MATRIX_HPP_

#include "pexsi/environment.hpp"
#include "pexsi/NumVec.hpp"

namespace  PEXSI{

/// @struct SparseMatrix
///
/// @brief SparseMatrix describes a sequential sparse matrix saved in
/// compressed sparse column format.
///
/// Note
/// ----
///
/// Since in PEXSI and PPEXSI only symmetric matrix is considered, the
/// compressed sparse row format will also be represented by the
/// compressed sparse column format.
template <class F> struct SparseMatrix{
  Int          size;                            // Matrix dimension
  Int          nnz;                             // Number of nonzeros
  IntNumVec    colptr;                          // Column index pointer
  IntNumVec    rowind;                          // Starting row index pointer
  NumVec<F>    nzval;                           // Nonzero values for the sparse matrix
};

// Commonly used
typedef SparseMatrix<Real>       DblSparseMatrix;
typedef SparseMatrix<Complex>    CpxSparseMatrix;

/// @struct DistSparseMatrix
///
/// @brief DistSparseMatrix describes a Sparse matrix in the compressed
/// sparse column format (CSC) and distributed with column major partition.
///
/// Note
/// ----
///
/// Since in PEXSI and PPEXSI only symmetric matrix is considered, the
/// compressed sparse row format will also be represented by the
/// compressed sparse column format.
///
/// @todo Add the parameter of numColLocal
template <class F> struct DistSparseMatrix{
  /// @brief Matrix dimension.
  Int          size;

  /// @brief Local number of local nonzeros elements on this processor.
  Int          nnzLocal;

  /// @brief Total number of nonzeros elements.
  ///
  /// FIXME: The datatype should be changed to LongInt in the future.
  Int          nnz;

  /// @brief Dimension numColLocal + 1, storing the pointers to the
  /// nonzero row indices and nonzero values in rowptrLocal and
  /// nzvalLocal, respectively.  numColLocal is the number
  /// of local columns saved on this processor. The indices are 1-based
  /// (FORTRAN-convention), i.e.  colptrLocal[0] = 1.
  IntNumVec    colptrLocal;

  /// @brief Dimension nnzLocal, storing the nonzero row indices.
  /// The indices are 1-based (FORTRAN-convention), i.e. the first row
  /// index is 1.
  IntNumVec    rowindLocal;

  /// @brief Dimension nnzLocal, storing the nonzero values.
  NumVec<F>    nzvalLocal;

  /// @brief MPI communicator
  MPI_Comm     comm;

  /// @brief Compute the total number of nonzeros through
  /// MPI_Allreduce
  LongInt      Nnz();

  /// @brief Locally sorts the row indices within every column
  void         SortIndices();

  /// @brief Clear all memory (The matrix becomes a 0-by-0 matrix)
  void Clear();

};

// Commonly used
typedef DistSparseMatrix<Real>       DblDistSparseMatrix;
typedef DistSparseMatrix<Complex>    CpxDistSparseMatrix;



} // namespace PEXSI

#include "pexsi/sparse_matrix_impl.hpp"

#endif // _PEXSI_SPARSE_MATRIX_HPP_
