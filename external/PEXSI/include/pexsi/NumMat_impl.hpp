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
/// @file NumMat_impl.hpp
/// @brief Implementation of numerical matrix.
/// @date 2010-09-27
#ifndef _PEXSI_NUMMAT_IMPL_HPP_
#define _PEXSI_NUMMAT_IMPL_HPP_


namespace  PEXSI{

template <class F> inline void NumMat<F>::allocate(F* data) {
  if(owndata_) {
    if(m_>0 && n_>0) { data_ = new F[m_*n_]; if( data_ == NULL ) {
      ErrorHandling("Cannot allocate memory.");}
    } else data_=NULL;
    if(data!=NULL){std::copy(data,data+m_*n_,data_);}
  } else {
    data_ = data;
  }
  bufsize_=m_*n_;
}
template <class F> inline void NumMat<F>::deallocate(){
  if(owndata_) {
    if(bufsize_>0) { delete[] data_; data_ = NULL; bufsize_ = 0; m_=0; n_=0; }
  }
}

template <class F> NumMat<F>::NumMat(Int m, Int n): m_(m), n_(n), owndata_(true) {
  this->allocate();
}

template <class F> NumMat<F>::NumMat(Int m, Int n, bool owndata, F* data): m_(m), n_(n), owndata_(owndata) {
  this->allocate(data);
}

template <class F> NumMat<F>::NumMat(const NumMat& C): m_(C.m_), n_(C.n_), owndata_(C.owndata_) {
  this->allocate(C.data_);
}

template <class F> NumMat<F>::~NumMat() {
  this->deallocate();
}

template <class F> NumMat<F>& NumMat<F>::Copy(const NumMat<F>& C) {
  this->deallocate();
  m_ = C.m_; n_=C.n_; owndata_=C.owndata_;
  this->allocate(C.data_);
  return *this;
}

template <class F> NumMat<F>& NumMat<F>::operator=(const NumMat<F>& C) {
  this->deallocate();
  m_ = C.m_; n_=C.n_; owndata_=C.owndata_;
  this->allocate(C.data_);
  return *this;
}

template <class F> void NumMat<F>::Resize(Int m, Int n)  {
  if( owndata_ == false ){
    ErrorHandling("Matrix being resized must own data.");
  }

  if(m*n > bufsize_) {
    this->deallocate();
    m_ = m; n_ = n;
    this->allocate();
  }
  else{
    m_ = m; n_ = n;
  }
}

template <class F> void NumMat<F>::Clear()  {
  if( owndata_ == false ){
    ErrorHandling("Matrix being cleared must own data.");
  }

  this->deallocate();
  m_ = 0; n_ = 0;
  bufsize_=0;
}





template <class F> const F& NumMat<F>::operator()(Int i, Int j) const  {
  if( i < 0 || i >= m_ ||
      j < 0 || j >= n_ ) {
    ErrorHandling( "Index is out of bound." );
  }
  return data_[i+j*m_];
}

template <class F> F& NumMat<F>::operator()(Int i, Int j)  {
  if( i < 0 || i >= m_ ||
      j < 0 || j >= n_ ) {
    ErrorHandling( "Index is out of bound." );
  }
  return data_[i+j*m_];
}

template <class F> F* NumMat<F>::VecData(Int j)  const
{
  if( j < 0 || j >= n_ ) {
    ErrorHandling( "Index is out of bound." );
  }
  return &(data_[j*m_]);
}


template <class F> inline void SetValue(NumMat<F>& M, F val)
{
  std::fill(M.Data(),M.Data()+M.m()*M.n(),val);
}

template <class F> inline Real Energy(const NumMat<F>& M)
{
  Real sum = 0;
  F *ptr = M.Data();
  for (Int i=0; i < M.m()*M.n(); i++)
    sum += abs(ptr[i]) * abs(ptr[i]);
  return sum;
}


template <class F> inline void
Transpose ( const NumMat<F>& A, NumMat<F>& B )
{
  NumMat<F> tmp;
  F* Adata = A.Data();
  if(&A==&B){
    //inplace transpose
    tmp = A;
    Adata = tmp.Data();
  }

  if( A.m() != B.n() || A.n() != B.m() ){
    B.Resize( A.n(), A.m() );
  }

  F* Bdata = B.Data();
  Int m = A.m(), n = A.n();

//#pragma omp taskloop collapse(2)
  for( Int i = 0; i < m; i++ ){
    for( Int j = 0; j < n; j++ ){
      //      Bdata[ j + n*i ] = F(std::conj(Adata[ i + j*m ]));
      Bdata[ j + n*i ] = Adata[ i + j*m ];
    }
  }

  return ;
}		// -----  end of function Transpose  -----

template <class F> inline void
ConjTranspose ( const NumMat<F>& A, NumMat<F>& B )
{
  NumMat<F> tmp;
  F* Adata = A.Data();
  if(&A==&B){
    //inplace transpose
    tmp = A;
    Adata = tmp.Data();
  }

  if( A.m() != B.n() || A.n() != B.m() ){
    B.Resize( A.n(), A.m() );
  }

  F* Bdata = B.Data();
  Int m = A.m(), n = A.n();

  for( Int i = 0; i < m; i++ ){
    for( Int j = 0; j < n; j++ ){
      //      Bdata[ j + n*i ] = F(std::conj(Adata[ i + j*m ]));
      Bdata[ j + n*i ] = std::conj(Adata[ i + j*m ]);
    }
  }

  return ;
}		// -----  end of function ConjTranspose  -----


template <class F> inline void
Symmetrize( NumMat<F>& A )
{
  if( A.m() != A.n() ){
    ErrorHandling( "The matrix to be symmetrized should be a square matrix." );
  }

  NumMat<F> B;
  Transpose( A, B );

  F* Adata = A.Data();
  F* Bdata = B.Data();

  F  half = (F) 0.5;

  for( Int i = 0; i < A.m() * A.n(); i++ ){
    *Adata = half * (*Adata + *Bdata);
    Adata++; Bdata++;
  }


  return ;
}		// -----  end of function Symmetrize -----


} // namespace PEXSI

#endif // _PEXSI_NUMMAT_IMPL_HPP_
