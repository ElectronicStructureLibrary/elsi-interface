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
/// @file SuperLUGrid_impl.hpp
/// @brief Implementation of SuperLU processor grid.
/// @date 2014-03-17
#ifndef _PEXSI_SUPERLUGRID_IMPL_HPP_
#define _PEXSI_SUPERLUGRID_IMPL_HPP_

namespace PEXSI{






inline SuperLUGrid<Real>::SuperLUGrid	( MPI_Comm comm, Int nprow, Int npcol )
{
  ptrData = new RealGridData;
  if( ptrData == NULL ){
    ErrorHandling( "SuperLUGrid cannot be allocated." );
  }
  ptrData->GridInit(comm, nprow, npcol);


  return ;
} 		// -----  end of method SuperLUGrid::SuperLUGrid  -----

inline SuperLUGrid<Real>::~SuperLUGrid	(  )
{
  // NOTE (07/21/2013): Since superlu_gridinit gets a copy of the
  // communicator, it is legal to call superlu_gridexit even if
  // grid->comm is a copy of MPI_COMM_WORLD.

  ptrData->GridExit();

  delete ptrData;

  return ;
} 		// -----  end of method SuperLUGrid::~SuperLUGrid  -----






inline SuperLUGrid<Real>::SuperLUGrid(const SuperLUGrid<Real> & g)
{

  if(g.ptrData == NULL){
    ErrorHandling( "Original SuperLUGrid not allocated." );
  }

  ptrData = new RealGridData(*g.ptrData);
  if( ptrData == NULL ){
    ErrorHandling( "SuperLUGrid cannot be allocated." );
  }


  return ;
} 		// -----  end of method SuperLUGrid::SuperLUGrid  -----

inline SuperLUGrid<Real> & SuperLUGrid<Real>::operator = (const SuperLUGrid<Real> & g){

  if(&g!=this){
    if(this->ptrData!=NULL){
      ptrData->GridExit();
      delete ptrData;
    }

    if(g.ptrData == NULL){
      ErrorHandling( "Original SuperLUGrid not allocated." );
    }
    else{
      ptrData = new RealGridData(*g.ptrData);
      if( ptrData == NULL ){
        ErrorHandling( "SuperLUGrid cannot be allocated." );
      }
    }
  }

  return *this;
}











inline SuperLUGrid<Complex>::SuperLUGrid	( MPI_Comm comm, Int nprow, Int npcol )
{
  ptrData = new ComplexGridData;
  if( ptrData == NULL ){
    ErrorHandling( "SuperLUGrid cannot be allocated." );
  }
  ptrData->GridInit(comm, nprow, npcol);


  return ;
} 		// -----  end of method SuperLUGrid::SuperLUGrid  -----


inline SuperLUGrid<Complex>::~SuperLUGrid	(  )
{
  // NOTE (07/21/2013): Since superlu_gridinit gets a copy of the
  // communicator, it is legal to call superlu_gridexit even if
  // grid->comm is a copy of MPI_COMM_WORLD.

  ptrData->GridExit();

  delete ptrData;

  return ;
} 		// -----  end of method SuperLUGrid::~SuperLUGrid  -----


inline SuperLUGrid<Complex>::SuperLUGrid(const SuperLUGrid<Complex> & g)
{

  if(g.ptrData == NULL){
    ErrorHandling( "Original SuperLUGrid not allocated." );
  }

  ptrData = new ComplexGridData(*g.ptrData);
  if( ptrData == NULL ){
    ErrorHandling( "SuperLUGrid cannot be allocated." );
  }


  return ;
} 		// -----  end of method SuperLUGrid::SuperLUGrid  -----

inline SuperLUGrid<Complex> & SuperLUGrid<Complex>::operator = (const SuperLUGrid<Complex> & g){

  if(&g!=this){
    if(this->ptrData!=NULL){
      ptrData->GridExit();
      delete ptrData;
    }

    if(g.ptrData == NULL){
      ErrorHandling( "Original SuperLUGrid not allocated." );
    }
    else{
      ptrData = new ComplexGridData(*g.ptrData);
      if( ptrData == NULL ){
        ErrorHandling( "SuperLUGrid cannot be allocated." );
      }
    }
  }

  return *this;
}











}


#endif //_PEXSI_SUPERLUGRID_IMPL_HPP_
