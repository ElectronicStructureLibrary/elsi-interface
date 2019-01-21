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
/// @file NumVec_impl.hpp
/// @brief Implementation of Numerical Vector.
/// @date 2010-09-27
#ifndef _PEXSI_NUMVEC_IMPL_HPP_
#define _PEXSI_NUMVEC_IMPL_HPP_


namespace  PEXSI{

// Templated form of numerical vectors
//
// The main advantage of this portable NumVec structure is that it can
// either own (owndata == true) or view (owndata == false) a piece of
// data.

template <class F> NumVec<F>::NumVec	() : m_(0), owndata_(true), data_(NULL), bufsize_(0)
{
  m_ = 0;
} 		// -----  end of method NumVec<F>::NumVec  -----


template <class F> NumVec<F>::NumVec	( Int m ) : m_(m), owndata_(true), data_(NULL), bufsize_(0)
{

  this->allocate();

} 		// -----  end of method NumVec<F>::NumVec  -----

template <class F> NumVec<F>::NumVec	( Int m, bool owndata, F* data ) : m_(m), owndata_(owndata), data_(NULL), bufsize_(0)
{

  this->allocate(data);

} 		// -----  end of method NumVec<F>::NumVec  -----

template <class F> NumVec<F>::NumVec	( const NumVec<F>& C ) : m_(C.m_), owndata_(C.owndata_), data_(NULL), bufsize_(0)
{
  this->allocate(C.data_);
} 		// -----  end of method NumVec<F>::NumVec  -----


template < class F > NumVec<F>::~NumVec	(  )
{
  this->deallocate();

} 		// -----  end of method NumVec<F>::~NumVec  -----


template < class F > inline NumVec<F>& NumVec<F>::operator =	( const NumVec& C  )
{
  this->deallocate();
  m_ = C.m_;
  owndata_ = C.owndata_;
  this->allocate(C.data_);

  return *this;
} 		// -----  end of method NumVec<F>::operator=  -----


template < class F > inline void NumVec<F>::Resize	( const Int m )
{
  if( owndata_ == false ){
    ErrorHandling("Vector being resized must own data.");
  }
  if(m > bufsize_) {
    this->deallocate();
    m_ = m;
    this->allocate();
  }
  else{
    m_ = m;
  }
  return ;
} 		// -----  end of method NumVec<F>::Resize  -----

template < class F > inline void NumVec<F>::Clear	( )
{
  if( owndata_ == false ){
    data_ = NULL;
  }
  else{
    this->deallocate();
  }
  return ;
} 		// -----  end of method NumVec<F>::Clear  -----



template <class F> inline F& NumVec<F>::operator()	( Int i )
{
  if( i < 0 || i >= m_ ){
    ErrorHandling( "Index is out of bound." );
  }
  return data_[i];

} 		// -----  end of method NumVec<F>::operator()  -----


template <class F>
inline const F&
NumVec<F>::operator()	( Int i ) const
{
  if( i < 0 || i >= m_ ){
    ErrorHandling( "Index is out of bound." );
  }
  return data_[i];

} 		// -----  end of method NumVec<F>::operator()  -----


template <class F> inline F& NumVec<F>::operator[]	( Int i )
{
  if( i < 0 || i >= m_ ){
    ErrorHandling( "Index is out of bound." );
  }
  return data_[i];
} 		// -----  end of method NumVec<F>::operator[]  -----


template <class F> inline const F& NumVec<F>::operator[]	( Int i ) const
{
  if( i < 0 || i >= m_ ){
    ErrorHandling( "Index is out of bound." );
  }
  return data_[i];

} 		// -----  end of method NumVec<F>::operator[]  -----

template <class F> inline void NumVec<F>::allocate(F* data) {
  if(owndata_) {
    if(m_>0) {
      data_ = new F[m_];
      if( data_ == NULL ) {
        ErrorHandling("Cannot allocate memory.");
      }
    }
    else{
      data_=NULL;
    }
    if(data!=NULL){
      std::copy(data,data+m_,data_);
    }
  }
  else {
    data_ = data;
  }
  bufsize_ = m_;
} 		// -----  end of method NumVec<F>::allocate  -----

template <class F> inline void NumVec<F>::deallocate() {
  if(owndata_) {
    if(bufsize_>0) { delete[] data_; data_ = NULL; bufsize_=0; m_=0; }
  }
} 		// -----  end of method NumVec<F>::deallocate  -----




template <class F> inline void SetValue( NumVec<F>& vec, F val )
{
  std::fill(vec.Data(),vec.Data()+vec.m(),val);
}

template <class F> inline Real Energy( const NumVec<F>& vec )
{
  Real sum = 0;
  for(Int i=0; i<vec.m(); i++){
    sum += std::abs(vec(i)*vec(i));
  }
  return sum;
}


} // namespace PEXSI

#endif // _PEXSI_NUMVEC_IMPL_HPP_
