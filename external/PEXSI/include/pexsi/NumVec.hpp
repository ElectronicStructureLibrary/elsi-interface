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
/// @file NumVec.hpp
/// @brief  Numerical vector.
/// @date 2010-09-27
#ifndef _PEXSI_NUMVEC_HPP_
#define _PEXSI_NUMVEC_HPP_

#include "pexsi/environment.hpp"

namespace  PEXSI{
/// @class NumVec
///
/// @brief Numerical vector.
///
/// NumVec is a portable encapsulation of a pointer to represent a 1D
/// vector. The main difference between NumVec<F> and std::vector<F> is
/// that NumVec<F> allows the vector to not owning the data, by
/// specifying (owndata_ == false).
template <class F> class NumVec
{
protected:

  /// @brief Helper function allocating the memory pointed by the data_ attribute
  inline void allocate(F* data=NULL);

  /// @brief Helper function freeing memory pointed by the data_ attribute
  inline void deallocate();

  /// @brief The size of the vector.
  Int  m_;
  ///
  /// @brief Whether it owns the data.
  bool owndata_;

  /// @brief The pointer for the actual data.
  F* data_;

  /// @brief The actual storage space allocated
  Int bufsize_;
public:
  NumVec();
  NumVec(Int m);
  NumVec(Int m, bool owndata, F* data);
  NumVec(const NumVec& C);
  ~NumVec();

  NumVec& operator=(const NumVec& C);

  void Resize ( Int m );
  void Clear();

  const F& operator()(Int i) const;
  F& operator()(Int i);
  const F& operator[](Int i) const;
  F& operator[](Int i);

  bool IsOwnData() const { return owndata_; }
  F*   Data() const { return data_; }
  Int  m () const { return m_; }
  Int ByteSize() const { return m_*sizeof(F);}
};

// Commonly used
typedef NumVec<bool>       BolNumVec;
typedef NumVec<Int>        IntNumVec;
typedef NumVec<Real>       DblNumVec;
typedef NumVec<Complex>    CpxNumVec;


// *********************************************************************
// Utility functions
// *********************************************************************
/// @brief SetValue sets a numerical vector to a constant val.
template <class F> inline void SetValue( NumVec<F>& vec, F val );

/// @brief Energy computes the L2 norm of a vector.
template <class F> inline Real Energy( const NumVec<F>& vec );


} // namespace PEXSI

#include "pexsi/NumVec_impl.hpp"

#endif // _PEXSI_NUMVEC_HPP_
