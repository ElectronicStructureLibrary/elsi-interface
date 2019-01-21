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
/// @file sparse_matrix_impl.hpp
/// @brief Implementation of sparse matrices.
/// @date 2012-11-28
#ifndef _PEXSI_SPARSE_MATRIX_IMPL_HPP_
#define _PEXSI_SPARSE_MATRIX_IMPL_HPP_

#include <iostream>
#include "pexsi/environment.hpp"

using std::ifstream;
using std::ofstream;

#include "pexsi/utility.hpp"


namespace  PEXSI{

struct IndexComparator {
  Int * lookup;
  IndexComparator(Int * v):lookup(v){}
  bool operator() (int i,int j) { return (lookup[i]<lookup[j]);}
};



template <typename F> inline void DistSparseMatrix<F>::SortIndices()
{


  for(Int col = 0; col<colptrLocal.m()-1;++col){
    Int colbeg = colptrLocal[col]-1;
    Int colend = colptrLocal[col+1]-1;

    Int numRow = colend - colbeg;

    std::vector<Int> rowsPerm(numRow);
    std::vector<Int> rowsSorted(numRow);
    std::vector<F> valsSorted(numRow);
    Int * rows = &rowindLocal[colbeg];
    F * vals = &nzvalLocal[colbeg];

    for(Int i = 0; i<rowsPerm.size(); ++i){ rowsPerm[i] = i;}
    IndexComparator cmp(rows);
    //sort the row indices
    std::sort(rowsPerm.begin(),rowsPerm.end(),cmp);

    for(Int i = 0; i<numRow; ++i){ rowsSorted[i] = rows[rowsPerm[i]]; }
    std::copy(rowsSorted.begin(),rowsSorted.end(),rows);

    for(Int i = 0; i<numRow; ++i){ valsSorted[i] = vals[rowsPerm[i]]; }
    std::copy(valsSorted.begin(),valsSorted.end(),vals);
  }




}

template <class F> inline LongInt DistSparseMatrix<F>::Nnz ( )
{
  LongInt nnzLocalLong = nnzLocal;
  LongInt nnz;

  MPI_Allreduce( &nnzLocalLong, &nnz, 1, MPI_LONG_LONG, MPI_SUM,
      comm );


  return nnz;
} 		// -----  end of method DistSparseMatrix<F>::Nnz  -----

template <typename F> inline void DistSparseMatrix<F>::Clear()
{


  size = 0;
  nnzLocal = 0;
  nnz = 0;
  colptrLocal.Clear();
  rowindLocal.Clear();
  nzvalLocal.Clear();







}





} // namespace PEXSI

#endif // _PEXSI_SPARSE_MATRIX_IMPL_HPP_
