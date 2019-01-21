/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.

Author: Lin Lin and Mathias Jacquelin

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
/// @file utility.hpp
/// @brief Various utility subroutines.
/// @date 2012-09-27
#ifndef _PEXSI_UTILITY_HPP_
#define _PEXSI_UTILITY_HPP_

#include  <stdlib.h>
#include "pexsi/environment.hpp"
#include "pexsi/tinyvec.hpp"
#include "pexsi/NumVec.hpp"
#include "pexsi/NumMat.hpp"
#include "pexsi/NumTns.hpp"
#include "pexsi/sparse_matrix.hpp"
#include "pexsi/mpi_interf.hpp"

#ifdef WITH_SYMPACK
#include <sympack.hpp>
#endif

namespace PEXSI{

// *********************************************************************
// Define constants
// *********************************************************************

// Write format control parameters
const int LENGTH_VAR_NAME = 8;
const int LENGTH_DBL_DATA = 16;
const int LENGTH_INT_DATA = 5;
const int LENGTH_VAR_UNIT = 6;
const int LENGTH_DBL_PREC = 8;
const int LENGTH_VAR_DATA = 16;


// standard case for most serialization/deserialization process.
const std::vector<Int> NO_MASK(1);

// *********************************************************************
// Formatted output stream
// *********************************************************************

// Bool is NOT defined due to ambiguity with Int

inline Int PrintBlock(std::ostream &os, const std::string name){

  os << std::endl<< "*********************************************************************" << std::endl;
  os << name << std::endl;
  os << "*********************************************************************" << std::endl << std::endl;
  return 0;
}

// String
inline Int Print(std::ostream &os, const std::string name) {
  os << std::setiosflags(std::ios::left) << name << std::endl;
  return 0;
};

inline Int Print(std::ostream &os, const char* name) {
  os << std::setiosflags(std::ios::left) << std::string(name) << std::endl;
  return 0;
};

inline Int Print(std::ostream &os, const std::string name, std::string val) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << name
    << std::setw(LENGTH_VAR_DATA) << val
    << std::endl;
  return 0;
};

inline Int Print(std::ostream &os, const std::string name, const char* val) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << name
    << std::setw(LENGTH_VAR_DATA) << std::string(val)
    << std::endl;
  return 0;
};


// Real

// one real number

inline Int Print(std::ostream &os, const std::string name, Real val) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << name
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_DBL_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::endl;
  return 0;
};

inline Int Print(std::ostream &os, const char* name, Real val) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name)
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_DBL_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::endl;
  return 0;
};


inline Int Print(std::ostream &os, const std::string name, Real val, const std::string unit) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << name
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_DBL_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_UNIT) << unit
    << std::endl;
  return 0;
};

inline Int Print(std::ostream &os, const char *name, Real val, const char *unit) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name)
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_DBL_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_UNIT) << std::string(unit)
    << std::endl;
  return 0;
};

// Two real numbers
inline Int Print(std::ostream &os, const std::string name1, Real val1, const std::string unit1,
    const std::string name2, Real val2, const std::string unit2) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << name1
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val1
    << std::setw(LENGTH_VAR_UNIT) << unit1
    << std::setw(LENGTH_VAR_NAME) << name2
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_UNIT) << unit2
    << std::endl;
  return 0;
};

inline Int Print(std::ostream &os, const char *name1, Real val1, const char *unit1,
    char *name2, Real val2, char *unit2) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name1)
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val1
    << std::setw(LENGTH_VAR_UNIT) << std::string(unit1)
    << std::setw(LENGTH_VAR_NAME) << std::string(name2)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_UNIT) << std::string(unit2)
    << std::endl;
  return 0;
};

// Int and Real
inline Int Print(std::ostream &os, const std::string name1, Int val1, const std::string unit1,
    const std::string name2, Real val2, const std::string unit2) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << name1
    << std::setw(LENGTH_INT_DATA) << val1
    << std::setw(LENGTH_VAR_UNIT) << unit1
    << std::setw(LENGTH_VAR_NAME) << name2
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_UNIT) << unit2
    << std::endl;
  return 0;
};

inline Int Print(std::ostream &os, const char *name1, Int val1, const char *unit1,
    char *name2, Real val2, char *unit2) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name1)
    << std::setw(LENGTH_INT_DATA) << val1
    << std::setw(LENGTH_VAR_UNIT) << std::string(unit1)
    << std::setw(LENGTH_VAR_NAME) << std::string(name2)
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_UNIT) << std::string(unit2)
    << std::endl;
  return 0;
};

inline Int Print(std::ostream &os,
    const char *name1, Int val1,
    const char *name2, Real val2,
    char *name3, Real val3) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name1)
    << std::setw(LENGTH_INT_DATA) << val1
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_NAME) << std::string(name2)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
    << std::setw(LENGTH_VAR_NAME) << std::string(name3)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val3
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::endl;
  return 0;
};


inline Int Print(std::ostream &os,
    const char *name1, Int val1,
    const char *name2, Real val2,
    const char *name3, Real val3,
    const char *name4, Real val4 ) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name1)
    << std::setw(LENGTH_INT_DATA) << val1
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_NAME) << std::string(name2)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
    << std::setw(LENGTH_VAR_NAME) << std::string(name3)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val3
    << std::setw(LENGTH_VAR_NAME) << std::string(name4)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val4
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::endl;
  return 0;
};

// Int

// one integer number
inline Int Print(std::ostream &os, std::string name, Int val) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << name
    << std::setw(LENGTH_VAR_DATA) << val
    << std::endl;
  return 0;
};


inline Int Print(std::ostream &os, const char *name, Int val) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name)
    << std::setw(LENGTH_VAR_DATA) << val
    << std::endl;
  return 0;
};


inline Int Print(std::ostream &os, const std::string name, Int val, const std::string unit) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << name
    << std::setw(LENGTH_VAR_DATA) << val
    << std::setw(LENGTH_VAR_UNIT) << unit
    << std::endl;
  return 0;
};


inline Int Print(std::ostream &os, const char* name, Int val, const std::string unit) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name)
    << std::setw(LENGTH_VAR_DATA) << val
    << std::setw(LENGTH_VAR_UNIT) << unit
    << std::endl;
  return 0;
};



// two integer numbers
inline Int Print(std::ostream &os, const std::string name1, Int val1, const std::string unit1,
    const std::string name2, Int val2, const std::string unit2) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << name1
    << std::setw(LENGTH_VAR_DATA) << val1
    << std::setw(LENGTH_VAR_UNIT) << unit1
    << std::setw(LENGTH_VAR_NAME) << name2
    << std::setw(LENGTH_VAR_DATA) << val2
    << std::setw(LENGTH_VAR_UNIT) << unit2
    << std::endl;
  return 0;
};

inline Int Print(std::ostream &os, const char *name1, Int val1, const char *unit1,
    char *name2, Int val2, char *unit2) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name1)
    << std::setw(LENGTH_VAR_DATA) << val1
    << std::setw(LENGTH_VAR_UNIT) << std::string(unit1)
    << std::setw(LENGTH_VAR_NAME) << std::string(name2)
    << std::setw(LENGTH_VAR_DATA) << val2
    << std::setw(LENGTH_VAR_UNIT) << std::string(unit2)
    << std::endl;
  return 0;
};

// Bool

// one boolean number
inline Int Print(std::ostream &os, const std::string name, bool val) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << name;
  if( val == true )
    os << std::setw(LENGTH_VAR_NAME) << "true" << std::endl;
  else
    os << std::setw(LENGTH_VAR_NAME) << "false" << std::endl;
  return 0;
};


inline Int Print(std::ostream &os, const char* name, bool val) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name);
  if( val == true )
    os << std::setw(LENGTH_VAR_NAME) << "true" << std::endl;
  else
    os << std::setw(LENGTH_VAR_NAME) << "false" << std::endl;
  return 0;
};


// Index 3 and Point 3
inline Int Print(std::ostream &os,
    const char *name1, Index3 val ) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name1)
    << std::setw(LENGTH_VAR_DATA) << val[0]
    << std::setw(LENGTH_VAR_DATA) << val[1]
    << std::setw(LENGTH_VAR_DATA) << val[2]
    << std::endl;
  return 0;
};

inline Int Print(std::ostream &os,
    const char *name1, Point3 val ) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name1)
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC) << val[0]
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC) << val[1]
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC) << val[2]
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::endl;
  return 0;
};

inline Int Print(std::ostream &os,
    const char *name1, Int val1,
    const char *name2, Point3 val ) {
  os << std::setiosflags(std::ios::left)
    << std::setw(LENGTH_VAR_NAME) << std::string(name1)
    << std::setw(LENGTH_INT_DATA) << val1
    << std::setw(LENGTH_VAR_NAME) << std::string(name2)
    << std::setiosflags(std::ios::scientific)
    << std::setiosflags(std::ios::showpos)
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC) <<val[0]
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC) <<val[1]
    << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC) <<val[2]
    << std::resetiosflags(std::ios::scientific)
    << std::resetiosflags(std::ios::showpos)
    << std::endl;
  return 0;
};

// *********************************************************************
// Overload << and >> operators for basic data types
// *********************************************************************

// std::vector
template <class F> inline std::ostream& operator<<( std::ostream& os, const std::vector<F>& vec)
{
  os<<vec.size()<<std::endl;
  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  for(Int i=0; i<vec.size(); i++)
    os<<" "<<vec[i];
  os<<std::endl;
  return os;
}

// NumVec
template <class F> inline std::ostream& operator<<( std::ostream& os, const NumVec<F>& vec)
{
  os<<vec.m()<<std::endl;
  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  for(Int i=0; i<vec.m(); i++)
    os<<" "<<vec(i);
  os<<std::endl;
  return os;
}

//template <class F> inline std::istream& operator>>( std::istream& is, NumVec<F>& vec)
//{
//	Int m;  is>>m;  vec.resize(m);
//	for(Int i=0; i<vec.m(); i++)
//		is >> vec(i);
//	return is;
//}

// NumMat
template <class F> inline std::ostream& operator<<( std::ostream& os, const NumMat<F>& mat)
{
  os<<mat.m()<<" "<<mat.n()<<std::endl;
  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  os << std::setprecision(16);
  for(Int i=0; i<mat.m(); i++) {
    for(Int j=0; j<mat.n(); j++)
      os<<" "<<mat(i,j);
    os<<std::endl;
  }
  return os;
}

// NumTns
template <class F> inline std::ostream& operator<<( std::ostream& os, const NumTns<F>& tns)
{
  os<<tns.m()<<" "<<tns.n()<<" "<<tns.p()<<std::endl;
  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  for(Int i=0; i<tns.m(); i++) {
    for(Int j=0; j<tns.n(); j++) {
      for(Int k=0; k<tns.p(); k++) {
        os<<" "<<tns(i,j,k);
      }
      os<<std::endl;
    }
    os<<std::endl;
  }
  return os;
}

// *********************************************************************
// serialize/deserialize for basic types
// More specific serialize/deserialize will be defined in individual
// class files
// *********************************************************************




//template
template <typename T>
inline Int serialize(const T& val, std::ostream& os, const std::vector<Int>& mask)
{
  os.write((char*)&val, sizeof(T));
  return 0;
}

template <typename T>
inline Int deserialize(T& val, std::istream& is, const std::vector<Int>& mask)
{
  is.read((char*)&val, sizeof(T));
  return 0;
}

//bool
inline Int serialize(const bool& val, std::ostream& os, const std::vector<Int>& mask)
{
  os.write((char*)&val, sizeof(bool));
  return 0;
}

inline Int deserialize(bool& val, std::istream& is, const std::vector<Int>& mask)
{
  is.read((char*)&val, sizeof(bool));
  return 0;
}

//char
inline Int serialize(const char& val, std::ostream& os, const std::vector<Int>& mask)
{
  os.write((char*)&val, sizeof(char));
  return 0;
}

inline Int deserialize(char& val, std::istream& is, const std::vector<Int>& mask)
{
  is.read((char*)&val, sizeof(char));
  return 0;
}

inline Int combine(char& val, char& ext)
{
  ErrorHandling( "Combine operation not implemented." );
  return 0;
}

//-------------------
//Int
inline Int serialize(const Int& val, std::ostream& os, const std::vector<Int>& mask)
{
  os.write((char*)&val, sizeof(Int));
  return 0;
}

inline Int deserialize(Int& val, std::istream& is, const std::vector<Int>& mask)
{
  is.read((char*)&val, sizeof(Int));
  return 0;
}

inline Int combine(Int& val, Int& ext)
{
  val += ext;
  return 0;
}

//-------------------
//LongInt
inline Int serialize(const LongInt& val, std::ostream& os, const std::vector<Int>& mask)
{
  os.write((char*)&val, sizeof(LongInt));
  return 0;
}

inline Int deserialize(LongInt& val, std::istream& is, const std::vector<Int>& mask)
{
  is.read((char*)&val, sizeof(LongInt));
  return 0;
}

inline Int combine(LongInt& val, LongInt& ext)
{
  val += ext;
  return 0;
}


//-------------------
//Real
inline Int serialize(const Real& val, std::ostream& os, const std::vector<Int>& mask)
{
  os.write((char*)&val, sizeof(Real));
  return 0;
}

inline Int deserialize(Real& val, std::istream& is, const std::vector<Int>& mask)
{
  is.read((char*)&val, sizeof(Real));
  return 0;
}

inline Int combine(Real& val, Real& ext)
{
  val += ext;
  return 0;
}

//-------------------
//Complex
inline Int serialize(const Complex& val, std::ostream& os, const std::vector<Int>& mask)
{
  os.write((char*)&val, sizeof(Complex));
  return 0;
}

inline Int deserialize(Complex& val, std::istream& is, const std::vector<Int>& mask)
{
  is.read((char*)&val, sizeof(Complex));
  return 0;
}

inline Int combine(Complex& val, Complex& ext)
{
  val += ext;
  return 0;
}

//-------------------
//Index2  TODO
//inline Int serialize(const Index2& val, std::ostream& os, const std::vector<Int>& mask)
//{
//  os.write((char*)&(val[0]), 2*sizeof(Int));
//  return 0;
//}
//
//inline Int deserialize(Index2& val, std::istream& is, const std::vector<Int>& mask)
//{
//  is.read((char*)&(val[0]), 2*sizeof(Int));
//  return 0;
//}
//
//inline Int combine(Index2& val, Index2& ext)
//{
//	ErrorHandling( "Combine operation not implemented." );
//  return 0;
//}

//-------------------
//Point2  TODO
//inline Int serialize(const Point2& val, std::ostream& os, const std::vector<Int>& mask)
//{
//  os.write((char*)&(val[0]), 2*sizeof(Real));
//  return 0;
//}
//
//inline Int deserialize(Point2& val, std::istream& is, const std::vector<Int>& mask)
//{
//  is.read((char*)&(val[0]), 2*sizeof(Real));
//  return 0;
//}
//
//inline Int combine(Point2& val, Point2& ext)
//{
//	ErrorHandling( "Combine operation not implemented." );
//  return 0;
//}

//-------------------
//Index3
inline Int serialize(const Index3& val, std::ostream& os, const std::vector<Int>& mask)
{
  os.write((char*)&(val[0]), 3*sizeof(Int));
  return 0;
}

inline Int deserialize(Index3& val, std::istream& is, const std::vector<Int>& mask)
{
  is.read((char*)&(val[0]), 3*sizeof(Int));
  return 0;
}

inline Int combine(Index3& val, Index3& ext)
{
  ErrorHandling( "Combine operation not implemented." );
  return 0;
}

//-------------------
//Point3
inline Int serialize(const Point3& val, std::ostream& os, const std::vector<Int>& mask)
{
  os.write((char*)&(val[0]), 3*sizeof(Real));
  return 0;
}

inline Int deserialize(Point3& val, std::istream& is, const std::vector<Int>& mask)
{
  is.read((char*)&(val[0]), 3*sizeof(Real));
  return 0;
}

inline Int combine(Point3& val, Point3& ext)
{
  ErrorHandling( "Combine operation not implemented." );
  return 0;
}

//-------------------
//std::vector
template<class T>
Int serialize(const std::vector<T>& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int sz = val.size();
  os.write((char*)&sz, sizeof(Int));
  for(Int k=0; k<sz; k++)
    serialize(val[k], os, mask);
  return 0;
}

template<class T>
Int deserialize(std::vector<T>& val, std::istream& is, const std::vector<Int>& mask)
{
  Int sz;
  is.read((char*)&sz, sizeof(Int));
  val.resize(sz);
  for(Int k=0; k<sz; k++)
    deserialize(val[k], is, mask);
  return 0;
}

template<class T>
Int combine(std::vector<T>& val, std::vector<T>& ext)
{
  ErrorHandling( "Combine operation not implemented." );
  return 0;
}

//-------------------
//std::set
template<class T>
Int serialize(const std::set<T>& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int sz = val.size();
  os.write((char*)&sz, sizeof(Int));
  for(typename std::set<T>::const_iterator mi=val.begin(); mi!=val.end(); mi++)
    serialize((*mi), os, mask);
  return 0;
}

template<class T>
Int deserialize(std::set<T>& val, std::istream& is, const std::vector<Int>& mask)
{
  val.clear();
  Int sz;
  is.read((char*)&sz, sizeof(Int));
  for(Int k=0; k<sz; k++) {
    T t; deserialize(t, is, mask);
    val.insert(t);
  }
  return 0;
}

template<class T>
Int combine(std::set<T>& val, std::set<T>& ext)
{
  ErrorHandling( "Combine operation not implemented." );
  return 0;
}

//-------------------
//std::map
template<class T, class S>
Int serialize(const std::map<T,S>& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int sz = val.size();
  os.write((char*)&sz, sizeof(Int));
  for(typename std::map<T,S>::const_iterator mi=val.begin(); mi!=val.end(); mi++) {
    serialize((*mi).first, os, mask);
    serialize((*mi).second, os, mask);
  }
  return 0;
}

template<class T, class S>
Int deserialize(std::map<T,S>& val, std::istream& is, const std::vector<Int>& mask)
{
  val.clear();
  Int sz;
  is.read((char*)&sz, sizeof(Int));
  for(Int k=0; k<sz; k++) {
    T t;	deserialize(t, is, mask);
    S s;	deserialize(s, is, mask);
    val[t] = s;
  }
  return 0;
}

template<class T, class S>
Int combine(std::map<T,S>& val, std::map<T,S>& ext)
{
  ErrorHandling( "Combine operation not implemented." );
  return 0;
}

//-------------------
//std::pair
template<class T, class S>
Int serialize(const std::pair<T,S>& val, std::ostream& os, const std::vector<Int>& mask)
{
  serialize(val.first, os, mask);
  serialize(val.second, os, mask);
  return 0;
}

template<class T, class S>
Int deserialize(std::pair<T,S>& val, std::istream& is, const std::vector<Int>& mask)
{
  deserialize(val.first, is, mask);
  deserialize(val.second, is, mask);
  return 0;
}

template<class T, class S>
Int combine(std::pair<T,S>& val, std::pair<T,S>& ext)
{
  ErrorHandling( "Combine operation not implemented." );
  return 0;
}

/*
//-------------------
//BolNumVec
inline Int serialize(const BolNumVec& val, std::ostream& os, const std::vector<Int>& mask)
{
Int m = val.m();
os.write((char*)&m, sizeof(Int));
os.write((char*)(val.Data()), m*sizeof(bool));
return 0;
}

inline Int deserialize(BolNumVec& val, std::istream& is, const std::vector<Int>& mask)
{
Int m;
is.read((char*)&m, sizeof(Int));
val.Resize(m);
is.read((char*)(val.Data()), m*sizeof(bool));
return 0;
}

//-------------------
//BolNumMat
inline Int serialize(const BolNumMat& val, std::ostream& os, const std::vector<Int>& mask)
{
Int m = val.m();
Int n = val.n();
os.write((char*)&m, sizeof(Int));
os.write((char*)&n, sizeof(Int));
os.write((char*)(val.Data()), m*n*sizeof(bool));
return 0;
}

inline Int deserialize(BolNumMat& val, std::istream& is, const std::vector<Int>& mask)
{
Int m;
Int n;
is.read((char*)&m, sizeof(Int));
is.read((char*)&n, sizeof(Int));
val.Resize(m,n);
is.read((char*)(val.Data()), m*n*sizeof(bool));
return 0;
}

//-------------------
//BolNumTns
inline Int serialize(const BolNumTns& val, std::ostream& os, const std::vector<Int>& mask)
{
Int m = val.m();  Int n = val.n();  Int p = val.p();
os.write((char*)&m, sizeof(Int));
os.write((char*)&n, sizeof(Int));
os.write((char*)&p, sizeof(Int));
os.write((char*)(val.Data()), m*n*p*sizeof(bool));
return 0;
}

inline Int deserialize(BolNumTns& val, std::istream& is, const std::vector<Int>& mask)
{
Int m,n,p;
is.read((char*)&m, sizeof(Int));
is.read((char*)&n, sizeof(Int));
is.read((char*)&p, sizeof(Int));
val.Resize(m,n,p);
is.read((char*)(val.Data()), m*n*p*sizeof(bool));
return 0;
}
 */

//-------------------
//IntNumVec
inline Int serialize(const IntNumVec& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)(val.Data()), m*sizeof(Int));
  return 0;
}

inline Int deserialize(IntNumVec& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  is.read((char*)&m, sizeof(Int));
  val.Resize(m);
  is.read((char*)(val.Data()), m*sizeof(Int));
  return 0;
}

inline Int combine(IntNumVec& val, IntNumVec& ext)
{
  //val.resize(ext.m());
  assert(val.m()==ext.m());
  for(Int i=0; i<val.m(); i++)    val(i) += ext(i);
  return 0;
}


//-------------------
//IntNumMat
inline Int serialize(const IntNumMat& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  Int n = val.n();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  os.write((char*)(val.Data()), m*n*sizeof(Int));
  return 0;
}

inline Int deserialize(IntNumMat& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  Int n;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  val.Resize(m,n);
  is.read((char*)(val.Data()), m*n*sizeof(Int));
  return 0;
}

inline Int combine(IntNumMat& val, IntNumMat& ext)
{
  //val.resize(ext.m(),ext.n());
  assert(val.m()==ext.m() && val.n()==ext.n());
  for(Int i=0; i<val.m(); i++)
    for(Int j=0; j<val.n(); j++)
      val(i,j) += ext(i,j);
  return 0;
}

//-------------------
//IntNumTns
inline Int serialize(const IntNumTns& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();  Int n = val.n();  Int p = val.p();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  os.write((char*)&p, sizeof(Int));
  os.write((char*)(val.Data()), m*n*p*sizeof(Int));
  return 0;
}

inline Int deserialize(IntNumTns& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m,n,p;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  is.read((char*)&p, sizeof(Int));
  val.Resize(m,n,p);
  is.read((char*)(val.Data()), m*n*p*sizeof(Int));
  return 0;
}

inline Int combine(IntNumTns& val, IntNumTns& ext)
{
  //val.resize(ext.m(),ext.n(),ext.p());
  assert(val.m()==ext.m() && val.n()==ext.n() && val.p()==ext.p());
  for(Int i=0; i<val.m(); i++)
    for(Int j=0; j<val.n(); j++)
      for(Int k=0; k<val.p(); k++)
        val(i,j,k) += ext(i,j,k);
  return 0;
}

//-------------------
//DblNumVec
inline Int serialize(const DblNumVec& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)(val.Data()), m*sizeof(Real));
  return 0;
}

inline Int deserialize(DblNumVec& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  is.read((char*)&m, sizeof(Int));
  val.Resize(m);
  is.read((char*)(val.Data()), m*sizeof(Real));
  return 0;
}

inline Int combine(DblNumVec& val, DblNumVec& ext)
{
  //val.resize(ext.m());
  assert(val.m()==ext.m());
  for(Int i=0; i<val.m(); i++)    val(i) += ext(i);
  return 0;
}

//-------------------
//DblNumMat
inline Int serialize(const DblNumMat& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  Int n = val.n();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  os.write((char*)(val.Data()), m*n*sizeof(Real));
  return 0;
}

inline Int deserialize(DblNumMat& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  Int n;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  val.Resize(m,n);
  is.read((char*)(val.Data()), m*n*sizeof(Real));
  return 0;
}

inline Int combine(DblNumMat& val, DblNumMat& ext)
{
  //val.resize(ext.m(),ext.n());
  assert(val.m()==ext.m() && val.n()==ext.n());
  for(Int i=0; i<val.m(); i++)
    for(Int j=0; j<val.n(); j++)
      val(i,j) += ext(i,j);
  return 0;
}


//-------------------
//DblNumTns
inline Int serialize(const DblNumTns& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();  Int n = val.n();  Int p = val.p();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  os.write((char*)&p, sizeof(Int));
  os.write((char*)(val.Data()), m*n*p*sizeof(Real));
  return 0;
}

inline Int deserialize(DblNumTns& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m,n,p;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  is.read((char*)&p, sizeof(Int));
  val.Resize(m,n,p);
  is.read((char*)(val.Data()), m*n*p*sizeof(Real));
  return 0;
}

inline Int combine(DblNumTns& val, DblNumTns& ext)
{
  //val.resize(ext.m(),ext.n(),ext.p());
  assert(val.m()==ext.m() && val.n()==ext.n() && val.p()==ext.p());
  for(Int i=0; i<val.m(); i++)
    for(Int j=0; j<val.n(); j++)
      for(Int k=0; k<val.p(); k++)
        val(i,j,k) += ext(i,j,k);
  return 0;
}

//-------------------
//CpxNumVec
inline Int serialize(const CpxNumVec& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)(val.Data()), m*sizeof(Complex));
  return 0;
}

inline Int deserialize(CpxNumVec& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  is.read((char*)&m, sizeof(Int));
  val.Resize(m);
  is.read((char*)(val.Data()), m*sizeof(Complex));
  return 0;
}

inline Int combine(CpxNumVec& val, CpxNumVec& ext)
{
  //val.resize(ext.m());
  assert(val.m()==ext.m());
  for(Int i=0; i<val.m(); i++)    val(i) += ext(i);
  return 0;
}

//-------------------
//CpxNumMat
inline Int serialize(const CpxNumMat& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  Int n = val.n();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  os.write((char*)(val.Data()), m*n*sizeof(Complex));
  return 0;
}

inline Int deserialize(CpxNumMat& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  Int n;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  val.Resize(m,n);
  is.read((char*)(val.Data()), m*n*sizeof(Complex));
  return 0;
}

inline Int combine(CpxNumMat& val, CpxNumMat& ext)
{
  //val.resize(ext.m(),ext.n());
  assert(val.m()==ext.m() && val.n()==ext.n());
  for(Int i=0; i<val.m(); i++)
    for(Int j=0; j<val.n(); j++)
      val(i,j) += ext(i,j);
  return 0;
}

//-------------------
//CpxNumTns
inline Int serialize(const CpxNumTns& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();  Int n = val.n();  Int p = val.p();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  os.write((char*)&p, sizeof(Int));
  os.write((char*)(val.Data()), m*n*p*sizeof(Complex));
  return 0;
}

inline Int deserialize(CpxNumTns& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m,n,p;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  is.read((char*)&p, sizeof(Int));
  val.Resize(m,n,p);
  is.read((char*)(val.Data()), m*n*p*sizeof(Complex));
  return 0;
}

inline Int combine(CpxNumTns& val, CpxNumTns& ext)
{
  //val.resize(ext.m(),ext.n(),ext.p());
  assert(val.m()==ext.m() && val.n()==ext.n() && val.p()==ext.p());
  for(Int i=0; i<val.m(); i++)
    for(Int j=0; j<val.n(); j++)
      for(Int k=0; k<val.p(); k++)
        val(i,j,k) += ext(i,j,k);
  return 0;
}

//-------------------
//NumVec
template<class T>
Int inline serialize(const NumVec<T>& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  os.write((char*)&m, sizeof(Int));
  for(Int i=0; i<m; i++)
    serialize(val(i), os, mask);
  return 0;
}

template<class T>
Int inline deserialize(NumVec<T>& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  is.read((char*)&m, sizeof(Int));
  val.Resize(m);
  for(Int i=0; i<m; i++)
    deserialize(val(i), is, mask);
  return 0;
}

template<class T>
Int inline combine(NumVec<T>& val, NumVec<T>& ext)
{
  ErrorHandling( "Combine operation not implemented." );
  return 0;
}

//-------------------
//NumMat
template<class T>
Int inline serialize(const NumMat<T>& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  Int n = val.n();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  for(Int j=0; j<n; j++)
    for(Int i=0; i<m; i++)
      serialize(val(i,j), os, mask);
  return 0;
}
template<class T>
Int inline deserialize(NumMat<T>& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  Int n;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  val.Resize(m,n);
  for(Int j=0; j<n; j++)
    for(Int i=0; i<m; i++)
      deserialize(val(i,j), is, mask);
  return 0;
}

template<class T>
Int inline combine(NumMat<T>& val, NumMat<T>& ext)
{
  ErrorHandling( "Combine operation not implemented." );
  return 0;
}


//-------------------
//NumTns
template<class T>
Int inline serialize(const NumTns<T>& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  Int n = val.n();
  Int p = val.p();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  os.write((char*)&p, sizeof(Int));
  for(Int k=0; k<p; k++)
    for(Int j=0; j<n; j++)
      for(Int i=0; i<m; i++)
        serialize(val(i,j,k), os, mask);
  return 0;
}

template<class T>
Int inline deserialize(NumTns<T>& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  Int n;
  Int p;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  is.read((char*)&p, sizeof(Int));
  val.Resize(m,n,p);
  for(Int k=0; k<p; k++)
    for(Int j=0; j<n; j++)
      for(Int i=0; i<m; i++)
        deserialize(val(i,j,k), is, mask);
  return 0;
}

template<class T>
Int inline combine(NumTns<T>& val, NumTns<T>& ext)
{
  ErrorHandling( "Combine operation not implemented." );
  return 0;
}

//-------------------
//DistSparseMatrix
template<class T>
Int inline serialize(const DistSparseMatrix<T>& val, std::ostream& os, const std::vector<Int>& mask)
{
  serialize( val.size,        os, mask );
  serialize( val.nnz,         os, mask );
  serialize( val.nnzLocal,    os, mask );
  serialize( val.colptrLocal, os, mask );
  serialize( val.rowindLocal, os, mask );
  serialize( val.nzvalLocal,  os, mask );
  // No need to serialize the communicator
  return 0;
}

template<class T>
Int inline deserialize(DistSparseMatrix<T>& val, std::istream& is, const std::vector<Int>& mask)
{
  deserialize( val.size,        is, mask );
  deserialize( val.nnz,         is, mask );
  deserialize( val.nnzLocal,    is, mask );
  deserialize( val.colptrLocal, is, mask );
  deserialize( val.rowindLocal, is, mask );
  deserialize( val.nzvalLocal,  is, mask );
  // No need to deserialize the communicator
  return 0;
}


#ifdef WITH_SYMPACK
Int inline serialize(const symPACK::DistSparseMatrixGraph & val, std::ostream& os, const std::vector<Int>& mask)
{
  serialize( val.size,        os, mask );
  serialize( val.nnz,         os, mask );
  serialize( val.IsExpanded(), os, mask );
  serialize( val.mpirank,     os, mask );
  serialize( val.mpisize,     os, mask );
  serialize( val.baseval,     os, mask );
  serialize( val.keepDiag,    os, mask );
  serialize( val.sorted,      os, mask );
  serialize( val.vertexDist,  os, mask );
  serialize( val.colptr,      os, mask );
  serialize( val.rowind,      os, mask );
  // No need to serialize the communicator
  return 0;
}

Int inline deserialize(symPACK::DistSparseMatrixGraph& val, std::istream& is, const std::vector<Int>& mask)
{
  deserialize( val.size,        is, mask );
  deserialize( val.nnz,         is, mask );
  deserialize( val.expanded, is, mask );
  deserialize( val.mpirank,     is, mask );
  deserialize( val.mpisize,     is, mask );
  deserialize( val.baseval,     is, mask );
  deserialize( val.keepDiag,    is, mask );
  deserialize( val.sorted,      is, mask );
  deserialize( val.vertexDist,  is, mask );
  deserialize( val.colptr,      is, mask );
  deserialize( val.rowind,      is, mask );
  // No need to deserialize the communicator
  return 0;
}

template<class T>
Int inline serialize(const symPACK::DistSparseMatrix<T>& val, std::ostream& os, const std::vector<Int>& mask)
{
  serialize( val.size,            os, mask );
  serialize( val.nnz,             os, mask );
  serialize( val.GetLocalGraph(), os, mask );
  serialize( val.nzvalLocal,      os, mask );
  // No need to serialize the communicator
  return 0;
}

template<class T>
Int inline deserialize(symPACK::DistSparseMatrix<T>& val, std::istream& is, const std::vector<Int>& mask)
{
  deserialize( val.size,            is, mask );
  deserialize( val.nnz,             is, mask );
  deserialize( val.GetLocalGraph(), is, mask );
  deserialize( val.nzvalLocal,      is, mask );
  // No need to deserialize the communicator
  return 0;
}

template<class T>
void inline Convert(const DistSparseMatrix<T>& A, symPACK::DistSparseMatrix<T> & B)
{
    const MPI_Comm & comm = A.comm;
    int mpisize,mpirank;
    MPI_Comm_size(comm,&mpisize);
    MPI_Comm_rank(comm,&mpirank);

    B.comm = comm;
    B.size = A.size;
    B.nnz = A.nnz;

    //Initialize the graph
    symPACK::DistSparseMatrixGraph & BGraph = B.GetLocalGraph();
    //Int colPerProc = A.size / mpisize;
    //BGraph.vertexDist.resize(mpisize+1,colPerProc);
    //BGraph.vertexDist[0] = 1;
    //std::partial_sum(BGraph.vertexDist.begin(),BGraph.vertexDist.end(),BGraph.vertexDist.begin());
    //BGraph.vertexDist.back() = A.size+1;
    BGraph.size = B.size;
    BGraph.nnz = B.nnz;
    BGraph.SetComm(B.comm);
    BGraph.SetExpanded(true);
    BGraph.baseval = 1;
    BGraph.keepDiag = 1;
    BGraph.sorted = 1;
    BGraph.colptr.resize(BGraph.LocalVertexCount()+1);
    assert(A.colptrLocal.m()== BGraph.colptr.size());
    std::copy(A.colptrLocal.Data(),A.colptrLocal.Data()+BGraph.LocalVertexCount()+1,BGraph.colptr.data());

    //Copy nzval
    BGraph.rowind.resize(A.nnzLocal);
    std::copy(A.rowindLocal.Data(),A.rowindLocal.Data()+A.nnzLocal,BGraph.rowind.data());

    B.nzvalLocal.resize(A.nzvalLocal.m());
    std::copy(A.nzvalLocal.Data(),A.nzvalLocal.Data()+A.nzvalLocal.m(),B.nzvalLocal.data());
}


template<class T>
void inline Convert(symPACK::DistSparseMatrix<T>& B, DistSparseMatrix<T> & A)
{
        MPI_Comm & comm = B.comm;
        int mpisize,mpirank;
        MPI_Comm_size(comm,&mpisize);
        MPI_Comm_rank(comm,&mpirank);
        bool wasExpanded = B.GetLocalGraph().IsExpanded();
        if(!wasExpanded && B.nnz>0){
          B.ExpandSymmetric();
          B.SortGraph();
          B.GetLocalGraph().SetBaseval(1);
        }

        //TODO do we have to redistribute when vertexDist is not balanced ?
        Int colPerProc = B.size / mpisize;
        if(mpirank==mpisize-1){ colPerProc = B.size - mpirank*colPerProc;}
        assert(B.GetLocalGraph().LocalVertexCount() == colPerProc);

        A.comm = B.comm;
        A.size = B.size;
        A.nnz = B.nnz;
        A.colptrLocal.Resize(B.GetLocalGraph().colptr.size());
        std::copy(B.GetLocalGraph().colptr.begin(), B.GetLocalGraph().colptr.end(), A.colptrLocal.Data());
        A.rowindLocal.Resize(B.GetLocalGraph().rowind.size());
        std::copy(B.GetLocalGraph().rowind.begin(), B.GetLocalGraph().rowind.end(), A.rowindLocal.Data());
        A.nnzLocal = B.GetLocalGraph().rowind.size();

        A.nzvalLocal.Resize(B.nzvalLocal.size());
        std::copy(B.nzvalLocal.begin(),B.nzvalLocal.end(),A.nzvalLocal.Data());

        if(!wasExpanded && B.nnz>0){
          B.ToLowerTriangular();
        }
}




#endif



template<class T>
Int inline combine(DistSparseMatrix<T>& val, DistSparseMatrix<T>& ext)
{
  ErrorHandling( "Combine operation not implemented." );
  return 0;
}


// *********************************************************************
// Parallel IO functions
// *********************************************************************

Int SeparateRead(std::string name, std::istringstream& is);

Int SeparateWrite(std::string name, std::ostringstream& os);

Int SeparateWriteAscii(std::string name, std::ostringstream& os);

Int SharedRead(std::string name, std::istringstream& is);

Int SharedWrite(std::string name, std::ostringstream& os);

// *********************************************************************
// Ej
// *********************************************************************
inline void IdentityCol( Int col, NumVec<Real>& vec )
{
  for(Int i=0; i<std::min(col,vec.m()); i++)
    vec(i) = 0.0;

  if(col<vec.m())
    vec(col) = 1.0;


  for(Int i=col+1; i<vec.m(); i++)
    vec(i) = 0.0;
}

inline void IdentityCol( Int col, NumVec<Complex>& vec )
{
  for(Int i=0; i<std::min(col,vec.m()); i++)
    vec(i) = Complex(0.0,0.0);

  if(col<vec.m())
    vec(col) = Complex(1.0,0.0);


  for(Int i=col+1; i<vec.m(); i++)
    vec(i) = Complex(0.0,0.0);
}



inline void IdentityCol( IntNumVec & cols, NumMat<Real>& mat )
{
  for(Int j=0;j<cols.m();j++){
    Int col= cols(j);
    for(Int i=0; i<std::min(col,mat.m()); i++)
      mat(i,j) = 0.0;

    if(col<mat.m())
      mat(col,j) = 1.0;


    for(Int i=col+1; i<mat.m(); i++)
      mat(i,j) = 0.0;
  }
}


inline void IdentityCol( IntNumVec & cols, NumMat<Complex>& mat )
{
  for(Int j=0;j<cols.m();j++){
    Int col= cols(j);
    for(Int i=0; i<std::min(col,mat.m()); i++)
      mat(i,j) = Complex(0.0,0.0);

    if(col<mat.m())
      mat(col,j) = Complex(1.0,0.0);


    for(Int i=col+1; i<mat.m(); i++)
      mat(i,j) = Complex(0.0,0.0);
  }
}




// *********************************************************************
// Random numbers
// *********************************************************************
inline void SetRandomSeed(long int seed){
  srand48(seed);
}

inline Real UniformRandom(){
  return (Real)drand48();
}

inline void UniformRandom( NumVec<Real>& vec )
{
  for(Int i=0; i<vec.m(); i++)
    vec(i) = UniformRandom();
}

inline void UniformRandom( NumVec<Complex>& vec )
{
  for(Int i=0; i<vec.m(); i++)
    vec(i) = Complex(UniformRandom(), UniformRandom());
}

inline void UniformRandom( NumMat<Real>& M )
{
  Real *ptr = M.Data();
  for(Int i=0; i < M.m() * M.n(); i++)
    *(ptr++) = UniformRandom();
}

inline void UniformRandom( NumMat<Complex>& M )
{
  Complex *ptr = M.Data();
  for(Int i=0; i < M.m() * M.n(); i++)
    *(ptr++) = Complex(UniformRandom(), UniformRandom());
}


inline void UniformRandom( NumTns<Real>& T )
{
  Real *ptr = T.Data();
  for(Int i=0; i < T.m() * T.n() * T.p(); i++)
    *(ptr++) = UniformRandom();
}

inline void UniformRandom( NumTns<Complex>& T )
{
  Complex *ptr = T.Data();
  for(Int i=0; i < T.m() * T.n() * T.p(); i++)
    *(ptr++) = Complex(UniformRandom(), UniformRandom());
}

// *********************************************************************
// Timing
// *********************************************************************
inline void GetTime(Real&  t){
  t = MPI_Wtime();
}

// *********************************************************************
// Comparator
// *********************************************************************

// Real
inline bool PairLtComparator( const std::pair<Real, Int>& l,
    const std::pair<Real, Int>& r ){
  return l.first < r.first;
}

inline bool PairGtComparator( const std::pair<Real, Int>& l,
    const std::pair<Real, Int>& r ){
  return l.first > r.first;
}

// For sorting with indices
// Example usage:
//   std::sort(val.begin(), val.end(), IndexComp<std::vector<int>&>(indices));
template<class T>
struct IndexComp {
private:
  const T indices_;
public:
  IndexComp (const T indices) : indices_(indices) {}
  bool operator()(const size_t a, const size_t b) const
  { return indices_[a] < indices_[b]; }
};




// *********************************************************************
// Sparse Matrix
// *********************************************************************

template <typename T> void CSCToCSR(DistSparseMatrix<T>& sparseA, DistSparseMatrix<T> & sparseB );

// Real format
void ReadSparseMatrix ( const char* filename, SparseMatrix<Real>& spmat );

void ReadDistSparseMatrix( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );

void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );

void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );

void ReadDistSparseMatrixFormatted( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );

// Complex format
void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm );

void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm );

void ReadDistSparseMatrixFormatted( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm );

template<typename T>
std::string ToMatlabScalar( T val);

template<typename T>
std::string ToMatlabScalar( std::complex<T> val);

template<typename T>
void WriteDistSparseMatrixMatlab(const char * filename, DistSparseMatrix<T> & pspmat, MPI_Comm comm);

template <class F1, class F2>
void
CopyPattern	( const SparseMatrix<F1>& A, SparseMatrix<F2>& B )
{
  B.size        = A.size;
  B.nnz         = A.nnz;
  B.colptr      = A.colptr;
  B.rowind      = A.rowind;
  B.nzval.Resize( A.nnz );
  return ;
}		// -----  end of template function CopyPattern  -----


// TODO Real format
template<typename T>
void
GetDiagonal ( const DistSparseMatrix<T>& A,
    NumVec<T>& diag );


// Functions for DistSparseMatrix

template <class F1, class F2>
void
CopyPattern	( const DistSparseMatrix<F1>& A, DistSparseMatrix<F2>& B )
{
  B.size        = A.size;
  B.nnz         = A.nnz;
  B.nnzLocal    = A.nnzLocal;
  B.colptrLocal = A.colptrLocal;
  B.rowindLocal = A.rowindLocal;
  B.nzvalLocal.Resize( A.nnzLocal );
  B.comm        = A.comm;
  return ;
}		// -----  end of template function CopyPattern  -----

#ifdef WITH_SYMPACK
template <class F1, class F2>
void
CopyPattern	( const symPACK::DistSparseMatrix<F1>& A, symPACK::DistSparseMatrix<F2>& B )
{
  B.size        = A.size;
  B.nnz         = A.nnz;
  const symPACK::DistSparseMatrixGraph & Agraph = A.GetLocalGraph();
  symPACK::DistSparseMatrixGraph & Bgraph = B.GetLocalGraph();
  //copy the graph
  Bgraph = Agraph;
  B.nzvalLocal.resize( Bgraph.LocalEdgeCount() );
  B.comm        = A.comm;
  return ;
}		// -----  end of template function CopyPattern  -----
#endif


/// Perform
/// \f[
///   C = \beta C + \alpha A B
/// \f]

/// @brief Multiply a DistSparseMatrix with a vector that is
/// distributed across all processors participating the operation.
///
/// **Note** It is not assumed that the matrix is symmetric.
///
/// @param[in]  AMat  Distributed sparse matrix in CSC format.
/// @param[in]  B  Dense vector with the same values across all
///                procesors. The size of the vector is the same as
///                the size of the matrix A.
/// @param[in,out] C Dense vector with the same values across all
///                  procesors. The size of the vector is the same as
///                  the size of the matrix A and B.
template <class F>
void
DistSparseMatMultGlobalVec(
    const F                    alpha,
    const DistSparseMatrix<F>& AMat,
    const NumVec<F>&           B,
    const F                    beta,
    NumVec<F>&                 C )
{
  if( B.m() != AMat.size ){
    std::ostringstream msg;
    msg << std::endl
      << "The size of the matrix A is " << AMat.size << std::endl
      << "The size of the vector B is " << B.m() << std::endl
      << "and they must agree with each other.";
    ErrorHandling( msg.str().c_str() );
  }
  if( C.m() != AMat.size ){
    std::ostringstream msg;
    msg << std::endl
      << "The size of the matrix A is " << AMat.size << std::endl
      << "The size of the vector C is " << C.m() << std::endl
      << "and they must agree with each other.";
    ErrorHandling( msg.str().c_str() );
  }

  MPI_Comm  comm = AMat.comm;

  Int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );

  // Initiate a new vector for saving DeltaC = alpha * A * B
  // This is somewhat inefficient implementation but this should not
  // matter at this moment, since the global distribution of a dense
  // vector is not important anyway.
  NumVec<F>  DeltaClocal( AMat.size );
  NumVec<F>  DeltaC( AMat.size );
  SetValue( DeltaClocal, 0.0 );
  SetValue( DeltaC, 0.0 );

  Int numColLocal      = AMat.colptrLocal.m() - 1;
  Int numColLocalFirst = AMat.size / mpisize;
  Int firstCol         = mpirank * numColLocalFirst;

  // NOTE: DistSparseMatrix use 1-based indicies, while access to
  // vectors uses 0-based indicies.
  for( Int j = 0; j < numColLocal; j++ ){
    Int jcol = firstCol + j;
    for( Int i = AMat.colptrLocal(j)-1;
        i < AMat.colptrLocal(j+1)-1; i++ ){
      Int irow = AMat.rowindLocal(i) - 1;
      DeltaClocal( irow ) += AMat.nzvalLocal(i) * B(jcol);
    }
  } // for (j)

  mpi::Allreduce( DeltaClocal.Data(), DeltaC.Data(),
      AMat.size, MPI_SUM, comm );

  for( Int i = 0; i < AMat.size; i++ ){
    C(i) = beta * C(i) + alpha * DeltaC(i);
  }

  return ;
}		// -----  end of template function DistSparseMatMultGlobalVec  -----

// *********************************************************************
// Other numerical routines
// *********************************************************************

// Interpolation
/// @brief Linear interpolates from (x,y) to (xx,yy)
///
/// Note:
///
/// x and xx must be sorted in ascending order.
///
/// if xx[i] < x[0],     yy[i] = y[0]
///    xx[i] > x[end-1], yy[i] = y[end-1]
void
LinearInterpolation (
    const std::vector<Real>& x,
    const std::vector<Real>& y,
    const std::vector<Real>& xx,
    std::vector<Real>& yy );

// Root finding

/// @brief Root finding for monotonic non-decreasing functions.
///
/// Find the solution to \f$y(x) = val\f$ through a set of tabulated
/// values \f$\{(x_i,y_i)\}\f$. Both \f$\{x_i\}\f$ and \f$\{y_i\}\f$
/// follow non-decreasing order. The solution is obtained via linear
/// interpolation.
///
/// This is an internal routine used for searching the chemical potential.
///
/// @param[in] x Dimension: N.
/// @param[in] y Dimension: N.
/// @param[in] val Real scalar.
/// @return Root of \f$y(x) = val\f$.
inline Real
MonotoneRootFinding (
    const std::vector<Real>&  x,
    const std::vector<Real>&  y,
    Real val )
{
  Int numX = x.size();

  if( val <= y[0] || val >= y[numX-1] ){
    std::ostringstream msg;
    msg
      << "The root finding procedure cannot find the solution for y(x)="
      << val << std::endl << "here [min(y) max(y)] = [" << y[0] << ", "
      << y[numX-1] << "]" << std::endl;
    ErrorHandling( msg.str().c_str() );
  }

  std::vector<Real>::const_iterator vi = std::lower_bound(
      y.begin(), y.end(), val );

  Int idx = vi - y.begin();

  Real root;
  if( y[idx] == y[idx-1] ){
    root = x[idx-1];
  }
  else{
    // Linear interpolation
    root = x[idx-1] + ( x[idx] - x[idx-1] ) / ( y[idx] - y[idx-1] )
      * ( val - y[idx-1] );
  }


  return root;
}		// -----  end of function MonotoneRootFinding  -----




} // namespace PEXSI

#include "pexsi/utility_impl.hpp"

#endif // _PEXSI_UTILITY_HPP_
