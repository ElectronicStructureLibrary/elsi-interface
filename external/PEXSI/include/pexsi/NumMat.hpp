/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.

Authors: Lexing Ying, Mathias Jacquelin and Lin Lin

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
/// @file NumMat.hpp
/// @brief Numerical matrix.
/// @date 2010-09-27
#ifndef _PEXSI_NUMMAT_HPP_
#define _PEXSI_NUMMAT_HPP_

#include "pexsi/environment.hpp"


namespace  PEXSI{
/// @class NumMat
///
/// @brief Numerical matrix.
///
/// NumMat is a portable encapsulation of a pointer to represent a 2D
/// matrix, which can either own (owndata == true) or view (owndata ==
/// false) a piece of data.
template <class F>
  class NumMat
  {
  public:
    /// @brief Whether it owns the data.
    bool owndata_;

    /// @brief The size of the first dimension.
    Int m_;

    /// @brief The size of second dimension.
    Int n_;

    /// @brief The pointer for the actual data.
    F* data_;

    /// @brief Actual size of the allocated buffer. Similar to the capacity of a std::vector
    Int bufsize_;

    /// @brief Helper function allocating the memory pointed by the data_ attribute
    void allocate(F* data=NULL);

    /// @brief Helper function freeing memory pointed by the data_ attribute
    void deallocate();

  public:
    NumMat(Int m=0, Int n=0);
    NumMat(Int m, Int n, bool owndata, F* data);
    NumMat(const NumMat& C);
    ~NumMat();
    NumMat& Copy(const NumMat& C);
    NumMat& operator=(const NumMat& C);
    void Resize(Int m, Int n);
    void Clear();
    const F& operator()(Int i, Int j) const;
    F& operator()(Int i, Int j);
    F* Data() const { return data_; }
    F* VecData(Int j)  const;
    Int m() const { return m_; }
    Int n() const { return n_; }
    Int Size() const {return m_*n_;}
    Int ByteSize() const { return m_*n_*sizeof(F);}
    Int AllocatedSize() const {return bufsize_;}
  };

// Commonly used
typedef NumMat<bool>     BolNumMat;
typedef NumMat<Int>      IntNumMat;
typedef NumMat<Real>     DblNumMat;
typedef NumMat<Complex>  CpxNumMat;

// *********************************************************************
// Utility functions
// *********************************************************************
/// @brief SetValue sets a numerical matrix to a constant val.
template <class F> inline void SetValue(NumMat<F>& M, F val);

/// @brief Energy computes the L2 norm of a matrix (treated as a vector).
template <class F> inline Real Energy(const NumMat<F>& M);

} // namespace PEXSI



#include "pexsi/NumMat_impl.hpp"

#endif // _PEXSI_NUMMAT_HPP_
