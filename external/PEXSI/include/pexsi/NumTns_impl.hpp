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
/// @file NumTns_impl.hpp
/// @brief Implementation of numerical tensor.
/// @date 2010-09-27
#ifndef _PEXSI_NUMTNS_IMPL_HPP_
#define _PEXSI_NUMTNS_IMPL_HPP_


namespace  PEXSI{

template <class F> NumTns<F>::NumTns(Int m, Int n, Int p): m_(m), n_(n), p_(p), owndata_(true) {
  if(m_>0 && n_>0 && p_>0) { data_ = new F[m_*n_*p_]; if( data_ == NULL ) {
    ErrorHandling("Cannot allocate memory.");}
  } else data_=NULL;
}

template <class F> NumTns<F>::NumTns(Int m, Int n, Int p, bool owndata, F* data): m_(m), n_(n), p_(p), owndata_(owndata) {
  if(owndata_) {
    if(m_>0 && n_>0 && p_>0) { data_ = new F[m_*n_*p_]; if( data_ == NULL ) {
      ErrorHandling("Cannot allocate memory.");}
    } else data_=NULL;
    if(m_>0 && n_>0 && p_>0) { for(Int i=0; i<m_*n_*p_; i++) data_[i] = data[i]; }
  } else {
    data_ = data;
  }
}

template <class F> NumTns<F>::NumTns(const NumTns<F>& C): m_(C.m_), n_(C.n_), p_(C.p_), owndata_(C.owndata_) {
  if(owndata_) {
    if(m_>0 && n_>0 && p_>0) { data_ = new F[m_*n_*p_];
      if( data_ == NULL ) {
        ErrorHandling("Cannot allocate memory.");}
    } else data_=NULL;

    if(m_>0 && n_>0 && p_>0) { for(Int i=0; i<m_*n_*p_; i++) data_[i] = C.data_[i]; }
  } else {
    data_ = C.data_;
  }
}

template <class F> NumTns<F>::~NumTns() {
  if(owndata_) {
    if(m_>0 && n_>0 && p_>0) { delete[] data_; data_ = NULL; }
  }
}

template <class F> NumTns<F>& NumTns<F>::operator=(const NumTns<F>& C) {
  if(owndata_) {
    if(m_>0 && n_>0 && p_>0) { delete[] data_; data_ = NULL; }
  }
  m_ = C.m_; n_=C.n_; p_=C.p_; owndata_=C.owndata_;
  if(owndata_) {
    if(m_>0 && n_>0 && p_>0) { data_ = new F[m_*n_*p_]; if( data_ == NULL ) {
      ErrorHandling("Cannot allocate memory.");} } else data_=NULL;
    if(m_>0 && n_>0 && p_>0) { for(Int i=0; i<m_*n_*p_; i++) data_[i] = C.data_[i]; }
  } else {
    data_ = C.data_;
  }
  return *this;
}

template <class F> void NumTns<F>::Resize(Int m, Int n, Int p)  {
  if( owndata_ == false ){
    ErrorHandling("Tensor being resized must own data.");
  }
  if(m_!=m || n_!=n || p_!=p) {
    if(m_>0 && n_>0 && p_>0) { delete[] data_; data_ = NULL; }
    m_ = m; n_ = n; p_=p;
    if(m_>0 && n_>0 && p_>0) { data_ = new F[m_*n_*p_]; if( data_ == NULL ) {
      ErrorHandling("Cannot allocate memory.");}
    } else data_=NULL;
  }
}

template <class F> const F& NumTns<F>::operator()(Int i, Int j, Int k) const  {
  if( i < 0 || i >= m_ ||
      j < 0 || j >= n_ ||
      k < 0 || k >= p_ ) {
    ErrorHandling( "Index is out of bound." );
  }
  return data_[i+j*m_+k*m_*n_];
}

template <class F> F& NumTns<F>::operator()(Int i, Int j, Int k)  {
  if( i < 0 || i >= m_ ||
      j < 0 || j >= n_ ||
      k < 0 || k >= p_ ) {
    ErrorHandling( "Index is out of bound." );
  }
  return data_[i+j*m_+k*m_*n_];
}


template <class F> F* NumTns<F>::MatData (Int j) const {
  if( j < 0 || j >= p_ ) {
    ErrorHandling( "Index is out of bound." );
  }
  return &(data_[j*m_*n_]);
};

template <class F> F* NumTns<F>::VecData (Int j, Int k) const {
  if( j < 0 || j >= n_ ||
      k < 0 || k >= p_ ) {
    ErrorHandling( "Index is out of bound." );
  }

  return &(data_[k*m_*n_+j*m_]);
};

template <class F> inline void SetValue(NumTns<F>& T, F val)
{
  F *ptr = T.data_;
  for(Int i=0; i < T.m() * T.n() * T.p(); i++) *(ptr++) = val;

  return;
}



template <class F> inline Real Energy(const NumTns<F>& T)
{
  Real sum = 0;

  F *ptr = T.Data();
  for(Int i=0; i < T.m() * T.n() * T.p(); i++)
    sum += abs(ptr[i]) * abs(ptr[i]);

  return sum;
}

} // namespace PEXSI

#endif // _PEXSI_NUMTNS_IMPL_HPP_
