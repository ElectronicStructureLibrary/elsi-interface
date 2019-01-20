/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.

Authors: Mathias Jacquelin and Lin Lin

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
/// @file SuperLUGrid.hpp
/// @brief SuperLU processor grid.
/// @date 2014-03-17
#ifndef _PEXSI_SUPERLUGRID_HPP_
#define _PEXSI_SUPERLUGRID_HPP_

#include "pexsi/environment.hpp"

//#include "pexsi/superlu_RealGridData.hpp"
//#include "pexsi/superlu_ComplexGridData.hpp"


namespace PEXSI{





class RealGridInfo;
class ComplexGridInfo;

class ComplexGridData{
public:
  ComplexGridInfo * info_;
public:

  ComplexGridData();
  ~ComplexGridData();
  ComplexGridData(const ComplexGridData & g);
  ComplexGridData & operator = (const ComplexGridData & g);

  void GridInit( MPI_Comm comm, Int nprow, Int npcol );
  void GridExit(  );
};


class RealGridData{
public:
  RealGridInfo * info_;
public:

  RealGridData();
  ~RealGridData();
  RealGridData(const RealGridData & g);
  RealGridData & operator = (const RealGridData & g);

  void GridInit( MPI_Comm comm, Int nprow, Int npcol );
  void GridExit(  );
};




/// @class SuperLUGrid
/// @brief A thin interface for the gridinfo_t structure in SuperLU.
template< typename T > class SuperLUGrid{
  /// @brief SuperLUMatrix can have access to the grid information.
public:
  void *        ptrData;
public:
  SuperLUGrid( MPI_Comm comm, int nprow, int npcol ){};
  ~SuperLUGrid(){};

  SuperLUGrid(const SuperLUGrid & g){};
  SuperLUGrid & operator = (const SuperLUGrid & g){};
};


template< > class SuperLUGrid<Real>{
  /// @brief SuperLUMatrix can have access to the grid information.
public:
  RealGridData*        ptrData;
public:
  SuperLUGrid( MPI_Comm comm, int nprow, int npcol );
  ~SuperLUGrid();

  SuperLUGrid(const SuperLUGrid & g);
  SuperLUGrid & operator = (const SuperLUGrid & g);
};


template< > class SuperLUGrid<Complex>{
  /// @brief SuperLUMatrix can have access to the grid information.
public:
  ComplexGridData*        ptrData;
public:
  SuperLUGrid( MPI_Comm comm, int nprow, int npcol );
  ~SuperLUGrid();


  SuperLUGrid(const SuperLUGrid & g);
  SuperLUGrid & operator = (const SuperLUGrid & g);
};

}


#include "pexsi/SuperLUGrid_impl.hpp"

#endif //_PEXSI_SUPERLUGRID_HPP_
