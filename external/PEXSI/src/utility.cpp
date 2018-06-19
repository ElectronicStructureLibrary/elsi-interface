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
/// @file utility.cpp
/// @brief Implementation of the non-templated utility subroutines.
/// @date 2012-09-20
#include "pexsi/utility.hpp"

using namespace std;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;

namespace PEXSI{

// *********************************************************************
// IO functions
// TODO Move this to utility.hpp and make them inline functions
// *********************************************************************
//---------------------------------------------------------
Int SeparateRead(std::string name, std::istringstream& is)
{
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);  //cerr<<filename<<endl;
  ifstream fin(filename);
  if( !fin.good() ){
    ErrorHandling( "File cannot be opened!" );
  }

  is.str( std::string(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>()) );
  fin.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}

//---------------------------------------------------------
Int SeparateWrite(std::string name, std::ostringstream& os)
{
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);
  ofstream fout(filename);
  if( !fout.good() ){
    ErrorHandling( "File cannot be opened!" );
  }
  fout<<os.str();
  fout.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}

//---------------------------------------------------------
Int SharedRead(std::string name, std::istringstream& is)
{
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  vector<char> tmpstr;
  if(mpirank==0) {
    ifstream fin(name.c_str());
    if( !fin.good() ){
      ErrorHandling( "File cannot be opened!" );
    }
    //std::string str(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>());
    //tmpstr.insert(tmpstr.end(), str.begin(), str.end());
    tmpstr.insert(tmpstr.end(), std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>());
    fin.close();
    int size = tmpstr.size();	//cerr<<size<<endl;
    MPI_Bcast((void*)&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast((void*)&(tmpstr[0]), size, MPI_BYTE, 0, MPI_COMM_WORLD);
  } else {
    int size;
    MPI_Bcast((void*)&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    tmpstr.resize(size);
    MPI_Bcast((void*)&(tmpstr[0]), size, MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  is.str( std::string(tmpstr.begin(), tmpstr.end()) );
  //
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}

//---------------------------------------------------------
Int SharedWrite(std::string name, std::ostringstream& os)
{
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  if(mpirank==0) {
    ofstream fout(name.c_str());
    if( !fout.good() ){
      ErrorHandling( "File cannot be opened!" );
    }
    // Corrected by Patrick Seewald 5/30/2015
    fout<<os.str();
    fout.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}


//---------------------------------------------------------
Int SeparateWriteAscii(std::string name, std::ostringstream& os)
{
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);
  ofstream fout(filename, ios::trunc);
  if( !fout.good() ){
    ErrorHandling( "File cannot be opened!" );
  }
  fout<<os.str();
  fout.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}


// *********************************************************************
// Other numerical routines
// *********************************************************************

void
LinearInterpolation (
    const std::vector<Real>& x,
    const std::vector<Real>& y,
    const std::vector<Real>& xx,
    std::vector<Real>& yy )
{
  Int numX  = x.size();
  Int numXX = xx.size();

  for( Int i = 1; i < numX; i++ ){
    if( x[i] <= x[i-1] ){
      //statusOFS<<"x["<<i<<"] = "<<x[i]<<std::endl;
      //statusOFS<<"x["<<i-1<<"] = "<<x[i-1]<<std::endl;
      //statusOFS<<"x"<<x<<std::endl;
      ErrorHandling("x must be sorted strictly ascendingly.");
    }
  }

  for( Int i = 1; i < numXX; i++){
    if( xx[i] < xx[i-1] ){
      //statusOFS<<"xx"<<x<<std::endl;
      ErrorHandling("xx must be sorted ascendingly.");
    }
  }


  yy.resize( numXX );
  std::vector<Real>::const_iterator vi;
  Int ix;
  for( Int i = 0; i < numXX; i++ ){
    if( xx[i] <= x[0] ){
      yy[i] = y[0];
    }
    else if( xx[i] >= x[numX-1] ){
      yy[i] = y[numX-1];
    }
    else{
      vi = std::lower_bound( x.begin(), x.end(), xx[i] );
      ix = vi - x.begin();

      yy[i] = y[ix-1] + (y[ix] - y[ix-1]) / (x[ix] - x[ix-1])
        * (xx[i] - x[ix-1]);
    }
  } // for (i)


  return ;
}		// -----  end of function LinearInterpolation  -----



}  // namespace PEXSI
