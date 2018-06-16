/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

Author: Mathias Jacquelin

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
/// @file flops.hpp
/// @brief Flop count utility subroutines.
/// @date 2016-12-09
#ifndef _PEXSI_FLOPS_HPP_ 
#define _PEXSI_FLOPS_HPP_

#include  <stdlib.h>
#include "pexsi/environment.hpp"
#include "pexsi/mpi_interf.hpp"

namespace flops{
//base flops (thx to M. Faverge && LAWN41)  
#define FMULS_AXPY(n_) (n_)
#define FADDS_AXPY(n_) (n_)

#define FMULS_GEMM(m_, n_, k_) ((m_) * (n_) * (k_))
#define FADDS_GEMM(m_, n_, k_) ((m_) * (n_) * (k_))

#define FMULS_TRMM_2(m_, n_) (0.5 * (n_) * (m_) * ((m_)+1))
#define FADDS_TRMM_2(m_, n_) (0.5 * (n_) * (m_) * ((m_)-1))
#define FMULS_TRMM(side_, m_, n_) ( ( (side_) == 'L' ) ? FMULS_TRMM_2((m_), (n_)) : FMULS_TRMM_2((n_), (m_)) )
#define FADDS_TRMM(side_, m_, n_) ( ( (side_) == 'L' ) ? FADDS_TRMM_2((m_), (n_)) : FADDS_TRMM_2((n_), (m_)) )
#define FMULS_TRSM FMULS_TRMM
#define FADDS_TRSM FADDS_TRMM

#define FMULS_GETRI(n_) ( (n_) * ((5. / 6.) + (n_) * ((2. / 3.) * (n_) + 0.5)) )
#define FADDS_GETRI(n_) ( (n_) * ((5. / 6.) + (n_) * ((2. / 3.) * (n_) - 1.5)) )


  template<typename T> bool is_complex(){ return false; }
  template<>  bool is_complex<std::complex<float> >(){ return true;  }
  template<>  bool is_complex<std::complex<double> >(){ return true;  }

  template<typename T>
  double Gemm(int m_, int n_, int k_){
    if (is_complex<T>()){
      return (6. * FMULS_GEMM((double)(m_), (double)(n_), (double)(k_)) + 2.0 * FADDS_GEMM((double)(m_), (double)(n_), (double)(k_)) );
    }
    else{
      return (     FMULS_GEMM((double)(m_), (double)(n_), (double)(k_)) +       FADDS_GEMM((double)(m_), (double)(n_), (double)(k_)) );
    }
  }


  template<typename T>
  double Getri(int n_){
    if (is_complex<T>()){
      return (6. * FMULS_GETRI((double)(n_)) + 2.0 * FADDS_GETRI((double)(n_)) );
    }
    else{
      return (     FMULS_GETRI((double)(n_)) +       FADDS_GETRI((double)(n_)) );
    }
  }

  template<typename T>
  double Axpy(int n_){
    if (is_complex<T>()){
      return (6. * FMULS_AXPY((double)(n_)) + 2.0 * FADDS_AXPY((double)(n_)) );
    }
    else{
      return (     FMULS_AXPY((double)(n_)) +       FADDS_AXPY((double)(n_)) );
    }
  }

  template<typename T>
  double Trsm(char side_, int m_, int n_){
    if (is_complex<T>()){
      return (6. * FMULS_TRSM(side_, (double)(m_), (double)(n_)) + 2.0 * FADDS_TRSM(side_, (double)(m_), (double)(n_)) );
    }
    else{
      return (     FMULS_TRSM(side_, (double)(m_), (double)(n_)) +       FADDS_TRSM(side_, (double)(m_), (double)(n_)) );
    }
  }

} // namespace flops


#endif //_PEXSI_FLOPS_HPP_
