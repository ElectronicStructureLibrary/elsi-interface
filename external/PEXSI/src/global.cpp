/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.

Author: Lin Lin

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
#include "pexsi/environment.hpp"
#include <deque>
#ifdef _COREDUMPER_
#include <google/coredumper.h>
//#include <cerrno>
//#include <string.h>
//#include <stdio.h>
#endif

namespace PEXSI{

// *********************************************************************
// IO
// *********************************************************************
std::ofstream  statusOFS;

#ifdef GEMM_PROFILE
std::ofstream  statOFS;
std::deque<int > gemm_stat;
#endif

#if defined(COMM_PROFILE) || defined(COMM_PROFILE_BCAST)
std::ofstream  commOFS;
std::deque<int > comm_stat;
#endif


// *********************************************************************
// Error handling
// *********************************************************************
void ErrorHandling( const char * msg ){
  std::cout<< msg <<std::endl;
  abort();
#ifdef _COREDUMPER_
  int mpirank, mpisize;
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );
  char filename[100];
  sprintf(filename, "core_%d_%d", mpirank, mpisize);

  if( WriteCoreDump(filename) ==0 ) {
    statusOFS << "success: WriteCoreDump to " << filename << std::endl;
  } else {
    statusOFS << "failed:  WriteCoreDump to " << filename << std::endl;
  }
#endif // #ifdef _COREDUMPER_
  throw std::runtime_error( msg );
}

} // namespace PEXSI
