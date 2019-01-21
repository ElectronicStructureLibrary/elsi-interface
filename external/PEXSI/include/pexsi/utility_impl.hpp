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
/// @file utility_impl.hpp
/// @brief Implementation of utility subroutines.
/// @date 2012-09-27
#ifndef _PEXSI_UTILITY_IMPL_HPP_
#define _PEXSI_UTILITY_IMPL_HPP_

#include <numeric>

//using namespace std;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;


namespace PEXSI{

// *********************************************************************
// Sparse Matrix
// *********************************************************************

//---------------------------------------------------------
inline void ReadSparseMatrix ( const char* filename, SparseMatrix<Real>& spmat )
{

  // FIXME
  // Binary format
  if( 1 ){
    std::istringstream iss;
    SharedRead( filename, iss );
    deserialize( spmat.size, iss, NO_MASK );
    deserialize( spmat.nnz,  iss, NO_MASK );
    deserialize( spmat.colptr, iss, NO_MASK );
    deserialize( spmat.rowind, iss, NO_MASK );
    deserialize( spmat.nzval, iss, NO_MASK );
  }

  // Ascii format
  if( 0 ) {
    ifstream fin(filename);
    fin >> spmat.size >> spmat.nnz;

    spmat.colptr.Resize( spmat.size+1 );
    spmat.rowind.Resize( spmat.nnz );
    spmat.nzval.Resize ( spmat.nnz );

    for( Int i = 0; i < spmat.size + 1; i++ ){
      fin >> spmat.colptr(i);
    }

    for( Int i = 0; i < spmat.nnz; i++ ){
      fin >> spmat.rowind(i);
    }

    for( Int i = 0; i < spmat.nnz; i++ ){
      fin >> spmat.nzval(i);
    }

    fin.close();
  }

  return ;
}		// -----  end of function ReadSparseMatrix  -----

//---------------------------------------------------------
inline void ReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  int mpirank;  MPI_Comm_rank(comm, &mpirank);
  int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  std::ifstream fin;

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  // Read basic information
  if( mpirank == 0 ){
    fin.open(filename);
    if( !fin.good() ){
      ErrorHandling( "File cannot be opened!" );
    }
    fin.read((char*)&pspmat.size, sizeof(Int));
    fin.read((char*)&pspmat.nnz,  sizeof(Int));
  }

  // FIXME Maybe need LongInt format to read the number of nonzeros later

  MPI_Bcast(&pspmat.size, 1, MPI_INT, 0, comm);
  MPI_Bcast(&pspmat.nnz,  1, MPI_INT, 0, comm);

  // Read colptr
  IntNumVec  colptr(pspmat.size+1);
  if( mpirank == 0 ){
    Int tmp;
    fin.read((char*)&tmp, sizeof(Int));

    if( tmp != pspmat.size+1 ){
      ErrorHandling( "colptr is not of the right size." );
    }

    fin.read((char*)colptr.Data(), sizeof(Int)*tmp);

  }

  MPI_Bcast(colptr.Data(), pspmat.size+1, MPI_INT, 0, comm);

  // Compute the number of columns on each processor
  IntNumVec numColLocalVec(mpisize);
  Int numColLocal, numColFirst;
  numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry
  numColLocal = numColLocalVec[mpirank];

  pspmat.colptrLocal.Resize( numColLocal + 1 );
  for( Int i = 0; i < numColLocal + 1; i++ ){
    pspmat.colptrLocal[i] = colptr[mpirank * numColFirst+i] - colptr[mpirank * numColFirst] + 1;
  }

  // Calculate nnz_loc on each processor
  pspmat.nnzLocal = pspmat.colptrLocal[numColLocal] - pspmat.colptrLocal[0];

  pspmat.rowindLocal.Resize( pspmat.nnzLocal );
  pspmat.nzvalLocal.Resize ( pspmat.nnzLocal );

  // Read and distribute the row indices
  if( mpirank == 0 ){
    Int tmp;
    fin.read((char*)&tmp, sizeof(Int));

    if( tmp != pspmat.nnz ){
      std::ostringstream msg;
      msg
        << "The number of nonzeros in row indices do not match." << std::endl
        << "nnz = " << pspmat.nnz << std::endl
        << "size of row indices = " << tmp << std::endl;
      ErrorHandling( msg.str().c_str() );
    }
    IntNumVec buf;
    Int numRead;
    for( Int ip = 0; ip < mpisize; ip++ ){
      numRead = colptr[ip*numColFirst + numColLocalVec[ip]] -
        colptr[ip*numColFirst];
      buf.Resize(numRead);
      fin.read( (char*)buf.Data(), numRead*sizeof(Int) );

      if( ip > 0 ){
        MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
        MPI_Send(buf.Data(), numRead, MPI_INT, ip, 1, comm);
      }
      else{
        pspmat.rowindLocal = buf;
      }
    }
  }
  else{
    Int numRead;
    MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
    if( numRead != pspmat.nnzLocal ){
      std::ostringstream msg;
      msg << "The number of columns in row indices do not match." << std::endl
        << "numRead  = " << numRead << std::endl
        << "nnzLocal = " << pspmat.nnzLocal << std::endl;
      ErrorHandling( msg.str().c_str() );
    }

    pspmat.rowindLocal.Resize( numRead );
    MPI_Recv( pspmat.rowindLocal.Data(), numRead, MPI_INT, 0, 1, comm, &mpistat );
  }

  // Read and distribute the nonzero values
  if( mpirank == 0 ){
    Int tmp;
    fin.read((char*)&tmp, sizeof(Int));

    if( tmp != pspmat.nnz ){
      std::ostringstream msg;
      msg
        << "The number of nonzeros in values do not match." << std::endl
        << "nnz = " << pspmat.nnz << std::endl
        << "size of values = " << tmp << std::endl;
      ErrorHandling( msg.str().c_str() );
    }
    NumVec<Real> buf;
    Int numRead;
    for( Int ip = 0; ip < mpisize; ip++ ){
      numRead = colptr[ip*numColFirst + numColLocalVec[ip]] -
        colptr[ip*numColFirst];
      buf.Resize(numRead);
      fin.read( (char*)buf.Data(), numRead*sizeof(Real) );

      if( ip > 0 ){
        MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
        MPI_Send(buf.Data(), numRead, MPI_DOUBLE, ip, 1, comm);
      }
      else{
        pspmat.nzvalLocal = buf;
      }
    }
  }
  else{
    Int numRead;
    MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
    if( numRead != pspmat.nnzLocal ){
      std::ostringstream msg;
      msg << "The number of columns in values do not match." << std::endl
        << "numRead  = " << numRead << std::endl
        << "nnzLocal = " << pspmat.nnzLocal << std::endl;
      ErrorHandling( msg.str().c_str() );
    }

    pspmat.nzvalLocal.Resize( numRead );
    MPI_Recv( pspmat.nzvalLocal.Data(), numRead, MPI_DOUBLE, 0, 1, comm, &mpistat );
  }

  // Close the file
  if( mpirank == 0 ){
    fin.close();
  }



  MPI_Barrier( comm );


  return ;
}		// -----  end of function ReadDistSparseMatrix  -----

inline void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  int mpirank;  MPI_Comm_rank(comm, &mpirank);
  int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  Int err = 0;



  int filemode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN;

  MPI_File fout;
  MPI_Status status;



  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fout);

  if (err != MPI_SUCCESS) {
    ErrorHandling( "File cannot be opened!" );
  }

  // FIXME Note that nnz uses the Int data type for consistency of writing / reading
  // Write header
  if( mpirank == 0 ){
    err = MPI_File_write_at(fout, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
    err = MPI_File_write_at(fout, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
  }


  // Compute the number of columns on each processor
  Int numColLocal = pspmat.colptrLocal.m()-1;
  Int numColFirst = pspmat.size / mpisize;
  IntNumVec  colptrChunk(numColLocal+1);

  Int prev_nz = 0;
  MPI_Exscan(&pspmat.nnzLocal, &prev_nz, 1, MPI_INT, MPI_SUM, comm);

  for( Int i = 0; i < numColLocal + 1; i++ ){
    colptrChunk[i] = pspmat.colptrLocal[i] + prev_nz;
  }


  MPI_Datatype memtype, filetype;
  MPI_Aint disps[6];
  int blklens[6];
  MPI_Datatype types[6] = {MPI_INT,MPI_INT, MPI_INT,MPI_INT, MPI_INT,MPI_DOUBLE};

  /* set block lengths (same for both types) */
  blklens[0] = (mpirank==0)?1:0;
  blklens[1] = numColLocal+1;
  blklens[2] = (mpirank==0)?1:0;
  blklens[3] = pspmat.nnzLocal;
  blklens[4] = (mpirank==0)?1:0;
  blklens[5] = pspmat.nnzLocal;




  //Calculate offsets
  MPI_Offset myColPtrOffset, myRowIdxOffset, myNzValOffset;
  myColPtrOffset = 3*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);
  myRowIdxOffset = 3*sizeof(int) + (pspmat.size +1  +  prev_nz)*sizeof(Int);
  myNzValOffset = 4*sizeof(int) + (pspmat.size +1 +  pspmat.nnz)*sizeof(Int)+ prev_nz*sizeof(Real);
  disps[0] = 2*sizeof(int);
  disps[1] = myColPtrOffset;
  disps[2] = myRowIdxOffset;
  disps[3] = sizeof(int)+myRowIdxOffset;
  disps[4] = myNzValOffset;
  disps[5] = sizeof(int)+myNzValOffset;



#if ( _DEBUGlevel_ >= 1 )
  char msg[200];
  char * tmp = msg;
  tmp += sprintf(tmp,"P%d ",mpirank);
  for(int i = 0; i<6; ++i){
    if(i==5)
      tmp += sprintf(tmp, "%d [%d - %d] | ",i,disps[i],disps[i]+blklens[i]*sizeof(double));
    else
      tmp += sprintf(tmp, "%d [%d - %d] | ",i,disps[i],disps[i]+blklens[i]*sizeof(int));
  }
  tmp += sprintf(tmp,"\n");
  printf("%s",msg);
#endif




  MPI_Type_create_struct(6, blklens, disps, types, &filetype);
  MPI_Type_commit(&filetype);

  /* create memory type */
  Int np1 = pspmat.size+1;
  MPI_Address( (void *)&np1,  &disps[0]);
  MPI_Address(colptrChunk.Data(), &disps[1]);
  MPI_Address( (void *)&pspmat.nnz,  &disps[2]);
  MPI_Address((void *)pspmat.rowindLocal.Data(),  &disps[3]);
  MPI_Address( (void *)&pspmat.nnz,  &disps[4]);
  MPI_Address((void *)pspmat.nzvalLocal.Data(),   &disps[5]);

  MPI_Type_create_struct(6, blklens, disps, types, &memtype);
  MPI_Type_commit(&memtype);



  /* set file view */
  err = MPI_File_set_view(fout, 0, MPI_BYTE, filetype, "native",MPI_INFO_NULL);

  /* everyone writes their own row offsets, columns, and
   * data with one big noncontiguous write (in memory and
   * file)
   */
  err = MPI_File_write_all(fout, MPI_BOTTOM, 1, memtype, &status);

  MPI_Type_free(&filetype);
  MPI_Type_free(&memtype);





  MPI_Barrier( comm );

  MPI_File_close(&fout);

  return ;
}		// -----  end of function ParaWriteDistSparseMatrix  -----

inline void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  int mpirank;  MPI_Comm_rank(comm, &mpirank);
  int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  MPI_Datatype type;
  int lens[3];
  MPI_Aint disps[3];
  MPI_Datatype types[3];
  Int err = 0;

  int filemode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;

  MPI_File fin;
  MPI_Status status;

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fin);

  if (err != MPI_SUCCESS) {
    ErrorHandling( "File cannot be opened!" );
  }

  // FIXME Note that nnz uses the Int data type for consistency of writing / reading
  // Read header
  if( mpirank == 0 ){
    err = MPI_File_read_at(fin, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
    err = MPI_File_read_at(fin, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
  }


  /* define a struct that describes all our data */
  lens[0] = 1;
  lens[1] = 1;
  MPI_Address(&pspmat.size, &disps[0]);
  MPI_Address(&pspmat.nnz, &disps[1]);
  types[0] = MPI_INT;
  types[1] = MPI_INT;
  MPI_Type_struct(2, lens, disps, types, &type);
  MPI_Type_commit(&type);


  /* broadcast the header data to everyone */
  MPI_Bcast(MPI_BOTTOM, 1, type, 0, comm);

  MPI_Type_free(&type);

  // Compute the number of columns on each processor
  IntNumVec numColLocalVec(mpisize);
  Int numColLocal, numColFirst;
  numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry
  numColLocal = numColLocalVec[mpirank];
  pspmat.colptrLocal.Resize( numColLocal + 1 );



  MPI_Offset myColPtrOffset = (2 + ((mpirank==0)?0:1) )*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);

  Int np1 = 0;
  lens[0] = (mpirank==0)?1:0;
  lens[1] = numColLocal + 1;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(pspmat.colptrLocal.Data(), &disps[1]);

  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myColPtrOffset, MPI_BOTTOM, 1, type, &status);

  if (err != MPI_SUCCESS) {
    ErrorHandling( "error reading colptr" );
  }
  MPI_Type_free(&type);

  // Calculate nnz_loc on each processor
  pspmat.nnzLocal = pspmat.colptrLocal[numColLocal] - pspmat.colptrLocal[0];


  pspmat.rowindLocal.Resize( pspmat.nnzLocal );
  pspmat.nzvalLocal.Resize ( pspmat.nnzLocal );

  //read rowIdx
  MPI_Offset myRowIdxOffset = (3 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + (pspmat.colptrLocal[0]-1))*sizeof(int);

  lens[0] = (mpirank==0)?1:0;
  lens[1] = pspmat.nnzLocal;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(pspmat.rowindLocal.Data(), &disps[1]);

  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myRowIdxOffset, MPI_BOTTOM, 1, type,&status);

  if (err != MPI_SUCCESS) {
    ErrorHandling( "error reading rowind" );
  }
  MPI_Type_free(&type);


  //read nzval
//  MPI_Offset myNzValOffset = (4 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(Int) + (pspmat.colptrLocal[0]-1)*sizeof(double);
//
//  lens[0] = (mpirank==0)?1:0;
//  lens[1] = pspmat.nnzLocal;
//
//  MPI_Address(&np1, &disps[0]);
//  MPI_Address(pspmat.nzvalLocal.Data(), &disps[1]);
//
//  types[0] = MPI_INT;
//  types[1] = MPI_DOUBLE;
//
//  MPI_Type_create_struct(2, lens, disps, types, &type);
//  MPI_Type_commit(&type);
//
//  err = MPI_File_read_at_all(fin, myNzValOffset, MPI_BOTTOM, 1, type,&status);
//  MPI_Type_free(&type);

  if( mpirank == 0 ){
    MPI_Offset myNzValOffset = (4 )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(int) + (pspmat.colptrLocal[0]-1)*sizeof(double);
    err = MPI_File_read_at(fin, myNzValOffset,(char*)&np1, 1, MPI_INT, &status);
  }

  MPI_Offset myNzValOffset = (5 )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(int) + (pspmat.colptrLocal[0]-1)*sizeof(double);
  err = MPI_File_read_at_all(fin, myNzValOffset, pspmat.nzvalLocal.Data(), pspmat.nnzLocal, MPI_DOUBLE,&status);



  if (err != MPI_SUCCESS) {
    ErrorHandling( "error reading nzval" );
  }



  //convert to local references
  for( Int i = 1; i < numColLocal + 1; i++ ){
    pspmat.colptrLocal[i] = pspmat.colptrLocal[i] -  pspmat.colptrLocal[0] + 1;
  }
  pspmat.colptrLocal[0]=1;

  MPI_Barrier( comm );

  MPI_File_close(&fin);

  return ;
}		// -----  end of function ParaReadDistSparseMatrix  -----

inline void ReadDistSparseMatrixFormatted ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  int mpirank;  MPI_Comm_rank(comm, &mpirank);
  int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  std::ifstream fin;

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  // FIXME Maybe need LongInt format to read the number of nonzeros later

  // Read basic information
  if( mpirank == 0 ){
    fin.open(filename);
    if( !fin.good() ){
      ErrorHandling( "File cannot be opened!" );
    }
    Int dummy;

    fin >> pspmat.size >> dummy ;
    fin >> pspmat.nnz >> dummy;
    // FIXME this is temporary and only applies to 4*4 matrix.
    //	  fin	>> dummy;
  }

  MPI_Bcast(&pspmat.size, 1, MPI_INT, 0, comm);
  MPI_Bcast(&pspmat.nnz,  1, MPI_INT, 0, comm);

  // Read colptr

  IntNumVec  colptr(pspmat.size+1);
  if( mpirank == 0 ){
    Int* ptr = colptr.Data();
    for( Int i = 0; i < pspmat.size+1; i++ )
      fin >> *(ptr++);
  }

  MPI_Bcast(colptr.Data(), pspmat.size+1, MPI_INT, 0, comm);

  // Compute the number of columns on each processor
  IntNumVec numColLocalVec(mpisize);
  Int numColLocal, numColFirst;
  numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry
  numColLocal = numColLocalVec[mpirank];

  pspmat.colptrLocal.Resize( numColLocal + 1 );
  for( Int i = 0; i < numColLocal + 1; i++ ){
    pspmat.colptrLocal[i] = colptr[mpirank * numColFirst+i] - colptr[mpirank * numColFirst] + 1;
  }

  // Calculate nnz_loc on each processor
  pspmat.nnzLocal = pspmat.colptrLocal[numColLocal] - pspmat.colptrLocal[0];

  pspmat.rowindLocal.Resize( pspmat.nnzLocal );
  pspmat.nzvalLocal.Resize ( pspmat.nnzLocal );

  // Read and distribute the row indices
  if( mpirank == 0 ){
    Int tmp;
    IntNumVec buf;
    Int numRead;
    for( Int ip = 0; ip < mpisize; ip++ ){
      numRead = colptr[ip*numColFirst + numColLocalVec[ip]] -
        colptr[ip*numColFirst];
      buf.Resize(numRead);
      Int *ptr = buf.Data();
      for( Int i = 0; i < numRead; i++ ){
        fin >> *(ptr++);
      }
      if( ip > 0 ){
        MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
        MPI_Send(buf.Data(), numRead, MPI_INT, ip, 1, comm);
      }
      else{
        pspmat.rowindLocal = buf;
      }
    }
  }
  else{
    Int numRead;
    MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
    if( numRead != pspmat.nnzLocal ){
      std::ostringstream msg;
      msg << "The number of columns in row indices do not match." << std::endl
        << "numRead  = " << numRead << std::endl
        << "nnzLocal = " << pspmat.nnzLocal << std::endl;
      ErrorHandling( msg.str().c_str() );
    }

    pspmat.rowindLocal.Resize( numRead );
    MPI_Recv( pspmat.rowindLocal.Data(), numRead, MPI_INT, 0, 1, comm, &mpistat );
  }

  //	std::cout << "Proc " << mpirank << " outputs rowindLocal.size() = "
  //		<< pspmat.rowindLocal.m() << endl;


  // Read and distribute the nonzero values
  if( mpirank == 0 ){
    Int tmp;
    NumVec<Real> buf;
    Int numRead;
    for( Int ip = 0; ip < mpisize; ip++ ){
      numRead = colptr[ip*numColFirst + numColLocalVec[ip]] -
        colptr[ip*numColFirst];
      buf.Resize(numRead);
      Real *ptr = buf.Data();
      for( Int i = 0; i < numRead; i++ ){
        fin >> *(ptr++);
      }
      if( ip > 0 ){
        MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
        MPI_Send(buf.Data(), numRead, MPI_DOUBLE, ip, 1, comm);
      }
      else{
        pspmat.nzvalLocal = buf;
      }
    }
  }
  else{
    Int numRead;
    MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
    if( numRead != pspmat.nnzLocal ){
      std::ostringstream msg;
      msg << "The number of columns in values do not match." << std::endl
        << "numRead  = " << numRead << std::endl
        << "nnzLocal = " << pspmat.nnzLocal << std::endl;
      ErrorHandling( msg.str().c_str() );
    }

    pspmat.nzvalLocal.Resize( numRead );
    MPI_Recv( pspmat.nzvalLocal.Data(), numRead, MPI_DOUBLE, 0, 1, comm, &mpistat );
  }

  // Close the file
  if( mpirank == 0 ){
    fin.close();
  }



  MPI_Barrier( comm );


  return ;
}		// -----  end of function ReadDistSparseMatrixFormatted  -----


inline void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm )
{
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  int mpirank;  MPI_Comm_rank(comm, &mpirank);
  int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  MPI_Datatype type;
  int lens[3];
  MPI_Aint disps[3];
  MPI_Datatype types[3];
  Int err = 0;

  int filemode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;

  MPI_File fin;
  MPI_Status status;

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fin);

  if (err != MPI_SUCCESS) {
    ErrorHandling( "File cannot be opened!" );
  }

  // FIXME Note that nnz uses the Int data type for consistency of writing / reading
  // Read header
  if( mpirank == 0 ){
    err = MPI_File_read_at(fin, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
    err = MPI_File_read_at(fin, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
  }


  /* define a struct that describes all our data */
  lens[0] = 1;
  lens[1] = 1;
  MPI_Address(&pspmat.size, &disps[0]);
  MPI_Address(&pspmat.nnz, &disps[1]);
  types[0] = MPI_INT;
  types[1] = MPI_INT;
  MPI_Type_struct(2, lens, disps, types, &type);
  MPI_Type_commit(&type);


  /* broadcast the header data to everyone */
  MPI_Bcast(MPI_BOTTOM, 1, type, 0, comm);

  MPI_Type_free(&type);

  // Compute the number of columns on each processor
  IntNumVec numColLocalVec(mpisize);
  Int numColLocal, numColFirst;
  numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry
  numColLocal = numColLocalVec[mpirank];
  pspmat.colptrLocal.Resize( numColLocal + 1 );



  MPI_Offset myColPtrOffset = (2 + ((mpirank==0)?0:1) )*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);

  Int np1 = 0;
  lens[0] = (mpirank==0)?1:0;
  lens[1] = numColLocal + 1;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(pspmat.colptrLocal.Data(), &disps[1]);

  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myColPtrOffset, MPI_BOTTOM, 1, type, &status);

  if (err != MPI_SUCCESS) {
    ErrorHandling( "error reading colptr" );
  }
  MPI_Type_free(&type);

  // Calculate nnz_loc on each processor
  pspmat.nnzLocal = pspmat.colptrLocal[numColLocal] - pspmat.colptrLocal[0];


  pspmat.rowindLocal.Resize( pspmat.nnzLocal );
  pspmat.nzvalLocal.Resize ( pspmat.nnzLocal );

  //read rowIdx
  MPI_Offset myRowIdxOffset = (3 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + pspmat.colptrLocal[0]-1)*sizeof(Int);

  lens[0] = (mpirank==0)?1:0;
  lens[1] = pspmat.nnzLocal;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(pspmat.rowindLocal.Data(), &disps[1]);

//  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myRowIdxOffset, MPI_BOTTOM, 1, type,&status);

  if (err != MPI_SUCCESS) {
    ErrorHandling( "error reading rowind" );
  }
  MPI_Type_free(&type);


  //read nzval
//  MPI_Offset myNzValOffset = (4 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(int) + (pspmat.colptrLocal[0]-1)*sizeof(Complex);

//  lens[0] = ((mpirank==0)?1:0)*sizeof(Int);
//  lens[1] = pspmat.nnzLocal*sizeof(Complex);
//
//  MPI_Address(&np1, &disps[0]);
//  MPI_Address(pspmat.nzvalLocal.Data(), &disps[1]);
//
//  types[0] = MPI_INT;
//  types[1] = MPI_DOUBLE_COMPLEX;
//  MPI_Type_create_struct(2, lens, disps, types, &type);
//  MPI_Type_commit(&type);
//
//  err = MPI_File_read_at_all(fin, myNzValOffset, MPI_BOTTOM, 1, type,&status);
//
//  MPI_Type_free(&type);

  if( mpirank == 0 ){
    MPI_Offset myNzValOffset = (4 )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(int) + (pspmat.colptrLocal[0]-1)*sizeof(Complex);
    err = MPI_File_read_at(fin, myNzValOffset,(char*)&np1, 1, MPI_INT, &status);
  }

  MPI_Offset myNzValOffset = (5 )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(int) + (pspmat.colptrLocal[0]-1)*sizeof(Complex);
  err = MPI_File_read_at_all(fin, myNzValOffset, pspmat.nzvalLocal.Data(), pspmat.nnzLocal, MPI_DOUBLE_COMPLEX,&status);



  if (err != MPI_SUCCESS) {
    ErrorHandling( "error reading nzval" );
  }



  //convert to local references
  for( Int i = 1; i < numColLocal + 1; i++ ){
    pspmat.colptrLocal[i] = pspmat.colptrLocal[i] -  pspmat.colptrLocal[0] + 1;
  }
  pspmat.colptrLocal[0]=1;

  MPI_Barrier( comm );

  MPI_File_close(&fin);

  return ;
}		// -----  end of function ParaReadDistSparseMatrix  -----

inline void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm )
{
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  int mpirank;  MPI_Comm_rank(comm, &mpirank);
  int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  Int err = 0;



  int filemode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN;

  MPI_File fout;
  MPI_Status status;



  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fout);

  if (err != MPI_SUCCESS) {
    ErrorHandling( "File cannot be opened!" );
  }

  // FIXME Note that nnz uses the Int data type for consistency of writing / reading
  // Write header
  if( mpirank == 0 ){
    err = MPI_File_write_at(fout, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
    err = MPI_File_write_at(fout, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
  }


  // Compute the number of columns on each processor
  Int numColLocal = pspmat.colptrLocal.m()-1;
  Int numColFirst = pspmat.size / mpisize;
  IntNumVec  colptrChunk(numColLocal+1);

  Int prev_nz = 0;
  MPI_Exscan(&pspmat.nnzLocal, &prev_nz, 1, MPI_INT, MPI_SUM, comm);

  for( Int i = 0; i < numColLocal + 1; i++ ){
    colptrChunk[i] = pspmat.colptrLocal[i] + prev_nz;
  }


  MPI_Datatype memtype, filetype;
  MPI_Aint disps[6];
  int blklens[6];
  MPI_Datatype types[6] = {MPI_INT,MPI_INT, MPI_INT,MPI_INT, MPI_INT,MPI_DOUBLE_COMPLEX};

  /* set block lengths (same for both types) */
  blklens[0] = (mpirank==0)?1:0;
  blklens[1] = numColLocal+1;
  blklens[2] = (mpirank==0)?1:0;
  blklens[3] = pspmat.nnzLocal;
  blklens[4] = (mpirank==0)?1:0;
  blklens[5] = pspmat.nnzLocal;




  //Calculate offsets
  MPI_Offset myColPtrOffset, myRowIdxOffset, myNzValOffset;
  myColPtrOffset = 3*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);
  myRowIdxOffset = 3*sizeof(int) + (pspmat.size +1  +  prev_nz)*sizeof(Int);
  myNzValOffset = 4*sizeof(int) + (pspmat.size +1 +  pspmat.nnz)*sizeof(Int)+ prev_nz*sizeof(Complex);
  disps[0] = 2*sizeof(int);
  disps[1] = myColPtrOffset;
  disps[2] = myRowIdxOffset;
  disps[3] = sizeof(int)+myRowIdxOffset;
  disps[4] = myNzValOffset;
  disps[5] = sizeof(int)+myNzValOffset;



#if ( _DEBUGlevel_ >= 1 )
  char msg[200];
  char * tmp = msg;
  tmp += sprintf(tmp,"P%d ",mpirank);
  for(int i = 0; i<6; ++i){
    if(i==5)
      tmp += sprintf(tmp, "%d [%d - %d] | ",i,disps[i],disps[i]+blklens[i]*sizeof(Complex));
    else
      tmp += sprintf(tmp, "%d [%d - %d] | ",i,disps[i],disps[i]+blklens[i]*sizeof(int));
  }
  tmp += sprintf(tmp,"\n");
  printf("%s",msg);
#endif




  MPI_Type_create_struct(6, blklens, disps, types, &filetype);
  MPI_Type_commit(&filetype);

  /* create memory type */
  Int np1 = pspmat.size+1;
  MPI_Address( (void *)&np1,  &disps[0]);
  MPI_Address(colptrChunk.Data(), &disps[1]);
  MPI_Address( (void *)&pspmat.nnz,  &disps[2]);
  MPI_Address((void *)pspmat.rowindLocal.Data(),  &disps[3]);
  MPI_Address( (void *)&pspmat.nnz,  &disps[4]);
  MPI_Address((void *)pspmat.nzvalLocal.Data(),   &disps[5]);

  MPI_Type_create_struct(6, blklens, disps, types, &memtype);
  MPI_Type_commit(&memtype);



  /* set file view */
  err = MPI_File_set_view(fout, 0, MPI_BYTE, filetype, "native",MPI_INFO_NULL);

  /* everyone writes their own row offsets, columns, and
   * data with one big noncontiguous write (in memory and
   * file)
   */
  err = MPI_File_write_all(fout, MPI_BOTTOM, 1, memtype, &status);

  MPI_Type_free(&filetype);
  MPI_Type_free(&memtype);





  MPI_Barrier( comm );

  MPI_File_close(&fout);

  return ;
}		// -----  end of function ParaWriteDistSparseMatrix  -----



inline void ReadDistSparseMatrixFormatted ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm )
{
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  int mpirank;  MPI_Comm_rank(comm, &mpirank);
  int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  std::ifstream fin;

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  // FIXME Maybe need LongInt format to read the number of nonzeros later

  // Read basic information
  if( mpirank == 0 ){
    fin.open(filename);
    if( !fin.good() ){
      ErrorHandling( "File cannot be opened!" );
    }
    Int dummy;

    fin >> pspmat.size >> dummy;
    fin >> pspmat.nnz >> dummy;
    // FIXME this is temporary and only applies to 4*4 matrix.
    //	  fin	>> dummy;
  }

  MPI_Bcast(&pspmat.size, 1, MPI_INT, 0, comm);
  MPI_Bcast(&pspmat.nnz,  1, MPI_INT, 0, comm);

  // Read colptr

  IntNumVec  colptr(pspmat.size+1);
  if( mpirank == 0 ){
    Int* ptr = colptr.Data();
    for( Int i = 0; i < pspmat.size+1; i++ )
      fin >> *(ptr++);
  }

  MPI_Bcast(colptr.Data(), pspmat.size+1, MPI_INT, 0, comm);

  // Compute the number of columns on each processor
  IntNumVec numColLocalVec(mpisize);
  Int numColLocal, numColFirst;
  numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry
  numColLocal = numColLocalVec[mpirank];

  pspmat.colptrLocal.Resize( numColLocal + 1 );
  for( Int i = 0; i < numColLocal + 1; i++ ){
    pspmat.colptrLocal[i] = colptr[mpirank * numColFirst+i] - colptr[mpirank * numColFirst] + 1;
  }

  // Calculate nnz_loc on each processor
  pspmat.nnzLocal = pspmat.colptrLocal[numColLocal] - pspmat.colptrLocal[0];

  pspmat.rowindLocal.Resize( pspmat.nnzLocal );
  pspmat.nzvalLocal.Resize ( pspmat.nnzLocal );

  // Read and distribute the row indices
  if( mpirank == 0 ){
    Int tmp;
    IntNumVec buf;
    Int numRead;
    for( Int ip = 0; ip < mpisize; ip++ ){
      numRead = colptr[ip*numColFirst + numColLocalVec[ip]] -
        colptr[ip*numColFirst];


      buf.Resize(numRead);
      Int *ptr = buf.Data();
      for( Int i = 0; i < numRead; i++ ){
        fin >> *(ptr++);
      }
      if( ip > 0 ){
        MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
        MPI_Send(buf.Data(), numRead, MPI_INT, ip, 1, comm);
      }
      else{
        pspmat.rowindLocal = buf;
      }
    }
  }
  else{
    Int numRead;
    MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
    if( numRead != pspmat.nnzLocal ){
      std::ostringstream msg;
      msg << "The number of columns in row indices do not match." << std::endl
        << "numRead  = " << numRead << std::endl
        << "nnzLocal = " << pspmat.nnzLocal << std::endl;
      ErrorHandling( msg.str().c_str() );
    }

    pspmat.rowindLocal.Resize( numRead );
    MPI_Recv( pspmat.rowindLocal.Data(), numRead, MPI_INT, 0, 1, comm, &mpistat );
  }

  //	std::cout << "Proc " << mpirank << " outputs rowindLocal.size() = "
  //		<< pspmat.rowindLocal.m() << endl;


  // Read and distribute the nonzero values
  if( mpirank == 0 ){
    Int tmp;
    NumVec<Real> buf;
    Int numRead;
    for( Int ip = 0; ip < mpisize; ip++ ){
      numRead = 2* (colptr[ip*numColFirst + numColLocalVec[ip]] -
          colptr[ip*numColFirst]);
      Real *ptr;
      if( ip > 0 ){
        buf.Resize(numRead);
        ptr = buf.Data();
      }
      else{
        ptr = reinterpret_cast<Real*>(pspmat.nzvalLocal.Data());
      }

      //read data in either buf or pspmat.nzvalLocal
      for( Int i = 0; i < numRead; i++ ){
        fin >> *(ptr++);
      }

      if( ip > 0 ){
        MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
        MPI_Send(buf.Data(), numRead, MPI_DOUBLE, ip, 1, comm);
      }
      //			else{
      //        std::copy( buf.Data(), buf.Data() + numRead,  pspmat.nzvalLocal.Data()  );
      //        memcpy( (void*)( pspmat.nzvalLocal.Data() ),
      //            (void*)buf.Data(), sizeof(Real)*numRead );
      //			}
    }
  }
  else{
    Int numRead;
    MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
    if( numRead != 2*pspmat.nnzLocal ){
      std::ostringstream msg;
      msg << "The number of columns in values do not match." << std::endl
        << "numRead    = " << numRead << std::endl
        << "2*nnzLocal = " << 2*pspmat.nnzLocal << std::endl;
      ErrorHandling( msg.str().c_str() );
    }

    pspmat.nzvalLocal.Resize( numRead/2 );
    MPI_Recv( reinterpret_cast<Real*>(pspmat.nzvalLocal.Data()),
        numRead, MPI_DOUBLE, 0, 1, comm, &mpistat );
  }

  // Close the file
  if( mpirank == 0 ){
    fin.close();
  }



  MPI_Barrier( comm );


  return ;
}		// -----  end of function ReadDistSparseMatrixFormatted  -----


template<typename T> inline void GetDiagonal ( const DistSparseMatrix<T>& A, NumVec<T>& diag )
{
  int mpirank, mpisize;
  MPI_Comm_rank( A.comm, &mpirank );
  MPI_Comm_size( A.comm, &mpisize );

  NumVec<T>	 diagLocal( A.size );
  SetValue( diagLocal, ZERO<T>() );
  diag.Resize( A.size );
  SetValue( diag, ZERO<T>() );

  Int numColFirst = A.size / mpisize;
  Int firstCol    = mpirank * numColFirst;
  Int numColLocal = A.colptrLocal.m() - 1;

#if ( _DEBUGlevel_ >= 1 )
  statusOFS << "numColFirst = " << numColFirst << std::endl;
  statusOFS << "A.nzvalLocal.size = " << A.nzvalLocal.m() << std::endl;
  statusOFS << "A.nnzLocal = " << A.nnzLocal << std::endl;
#endif

  // Note that the indices in DistSparseMatrix follows the FORTRAN convention
  for( Int j = 0; j < numColLocal; j++ ){
    Int jcol = j + firstCol + 1;
    Int numRow = A.colptrLocal(j+1) - A.colptrLocal(j);
    const Int* rowPtr = &A.rowindLocal( A.colptrLocal(j) - 1 );
    // NOTE: The rows in DistSparseMatrix are not necessarily ordered.
    // So lower_bound cannot be used here for fast searching. find has to be used.
    const Int* ptr = std::find( rowPtr, rowPtr + numRow, jcol );
    if( ptr == rowPtr + numRow ){
      std::ostringstream msg;
      msg << "Serious problem. Did not find the row corresponding to the column." << std::endl
        << "This happens when j = " << j << ", jcol = " << jcol << ", and the row indices are " << std::endl
        << IntNumVec( numRow, false, const_cast<Int*>(rowPtr) ) << std::endl;
      ErrorHandling( msg.str().c_str() );
    }
    Int diagIdx = ptr - A.rowindLocal.Data();
    diagLocal( jcol - 1 ) = A.nzvalLocal( diagIdx );
  }

  mpi::Allreduce( &diagLocal[0], &diag[0], A.size, MPI_SUM, A.comm );


  return ;
}		// -----  end of function GetDiagonal  -----




template <typename T> void CSCToCSR(DistSparseMatrix<T>& sparseA, DistSparseMatrix<T> & sparseB ){


  int mpirank;
  MPI_Comm_rank(sparseA.comm,&mpirank);

  int mpisize;
  MPI_Comm_size(sparseA.comm,&mpisize);

  Int numRowLocalFirst = sparseA.size / mpisize;
  Int firstRow = mpirank * numRowLocalFirst;
  Int numRowLocal = -1;
  Int nnzLocal = -1;

  sparseB.size = sparseA.size;
  sparseB.nnz = sparseA.nnz;
  sparseB.comm = sparseA.comm;


  LongInt nnz = 0;
  IntNumVec rowindGlobal;
  IntNumVec colptrGlobal;
  //TIMER_START(ToGlobalStructure);
  {
    colptrGlobal.Resize(sparseA.size+1);

    /* Allgatherv for row indices. */
    IntNumVec prevnz(mpisize);
    IntNumVec rcounts(mpisize);
    MPI_Allgather((void*)&sparseA.nnzLocal, 1, MPI_INT, rcounts.Data(), 1, MPI_INT, sparseA.comm);

    prevnz[0] = 0;
    for (Int i = 0; i < mpisize-1; ++i) { prevnz[i+1] = prevnz[i] + rcounts[i]; }

    nnz = 0;
    for (Int i = 0; i < mpisize; ++i) { nnz += rcounts[i]; }
    rowindGlobal.Resize(nnz);

    MPI_Allgatherv(sparseA.rowindLocal.Data(), sparseA.nnzLocal, MPI_INT, rowindGlobal.Data(),rcounts.Data(), prevnz.Data(), MPI_INT, sparseA.comm);

    /* Allgatherv for colptr */
    // Compute the number of columns on each processor
    Int numColFirst = std::max(1,sparseA.size / mpisize);
    SetValue( rcounts, numColFirst );
    rcounts[mpisize-1] = sparseA.size - numColFirst * (mpisize-1);  // Modify the last entry

    IntNumVec rdispls(mpisize);
    rdispls[0] = 0;
    for (Int i = 0; i < mpisize-1; ++i) { rdispls[i+1] = rdispls[i] + rcounts[i]; }

    MPI_Allgatherv(sparseA.colptrLocal.Data(), sparseA.colptrLocal.m()-1, MPI_INT, colptrGlobal.Data(),
        rcounts.Data(), rdispls.Data(), MPI_INT, sparseA.comm);

    /* Recompute column pointers. */
    for (Int p = 1; p < mpisize; p++) {
      Int idx = rdispls[p];
      for (Int j = 0; j < rcounts[p]; ++j) colptrGlobal[idx++] += prevnz[p];
    }
    colptrGlobal(sparseA.size)= nnz+1;
  }
  //TIMER_STOP(ToGlobalStructure);


  IntNumVec rowptrGlobal;
  IntNumVec colindGlobal;
  //Compute Global CSR structure and Local CSR structure
  {
    Int i, j;

    colindGlobal.Resize(nnz);
    rowptrGlobal.Resize(sparseA.size+1);
    IntNumVec row_nnz(sparseA.size);
    SetValue(row_nnz,I_ZERO);

    for (i = 0; i < nnz; i++) {
      Int k = rowindGlobal[i];
      row_nnz[k-1]++;
    }

    rowptrGlobal[0]=1;
    for (i = 1; i <= sparseA.size; i++)
    {
      rowptrGlobal[i] = rowptrGlobal[i-1] + row_nnz[i-1];
      row_nnz[i-1] = 0;
    }

    /*
     *  convert from CSC to CSR.
     *  use row_nnz to keep track of the number of non-zeros
     *  added to each row.
     */
    for (j = 0; j < colptrGlobal.m()-1; j++)
    {
      Int k;
      Int nnz_col;    /* # of non-zeros in column j of A */

      nnz_col = colptrGlobal[j+1] - colptrGlobal[j];

      for (k = colptrGlobal[j]; k < colptrGlobal[j+1]; k++)
      {
        Int i = rowindGlobal[ k - 1 ];               /* row index */
        Int h = rowptrGlobal[i-1] + row_nnz[i-1];    /* non-zero position */

        /* add the non-zero A(i,j) to B */
        colindGlobal[ h - 1 ] = j+1;
        row_nnz[i-1]++;
      }
    }



  }

  //This assertion is only ok for structurally symmetric matrices
  //for(int i=0;i<rowptrGlobal.m();++i){ assert(rowptrGlobal[i]==colptrGlobal[i]); }

  //for(int i=1;i<rowptrGlobal.m();++i){
  //  Int colbeg = rowptrGlobal[i-1];
  //  Int colend = rowptrGlobal[i]-1;
  //  for(Int j=colbeg;j<=colend;++j){
  //    Int row = colindGlobal[j-1];

  //    Int * ptr = std::find(&rowindGlobal[colbeg-1],&rowindGlobal[colend-1]+1,row);
  //  }
  //}

  //compute the local CSR structure
  std::vector<Int > numRowVec(mpisize,numRowLocalFirst);
  numRowVec[mpisize-1] = sparseA.size - (mpisize-1)*numRowLocalFirst;

  numRowLocal = numRowVec[mpirank];
  Int myFirstCol = mpirank*numRowLocalFirst+1;
  Int myLastCol = myFirstCol + numRowLocal -1;

  sparseB.colptrLocal.Resize(numRowLocal+1);
  Int * rowptrLocal = &sparseB.colptrLocal[0];

  //copy local chunk of rowptr
  std::copy(&rowptrGlobal[myFirstCol-1],
      &rowptrGlobal[myLastCol]+1,
      rowptrLocal);

  nnzLocal = rowptrLocal[numRowLocal] - rowptrLocal[0];

  sparseB.nnzLocal = nnzLocal;

  sparseB.rowindLocal.Resize(nnzLocal);
  Int * colindLocal = &sparseB.rowindLocal[0];

  sparseB.nzvalLocal.Resize(nnzLocal);
  T * nzvalLocal  = &sparseB.nzvalLocal[0];

  //copy local chunk of colind
  std::copy(&colindGlobal[rowptrLocal[0]-1],
      &colindGlobal[rowptrLocal[0]-1]+nnzLocal,
      colindLocal);

  //convert rowptr to local references
  for( Int i = 1; i < numRowLocal + 1; i++ ){
    rowptrLocal[i] = rowptrLocal[i] -  rowptrLocal[0] + 1;
  }
  rowptrLocal[0]=1;

  //for(int i=0;i<numRowLocal;++i){
  //  Int col = i+myFirstCol;
  //  Int colbeg = rowptrLocal[i];
  //  Int colend = rowptrLocal[i+1]-1;

  //  for(Int j=colbeg;j<=colend;++j){
  //    Int row = colindLocal[j-1];

  //    const Int * ptr = std::find(&sparseA.rowindLocal[colbeg-1],&sparseA.rowindLocal[colend-1]+1,row);
  //  }

  //}


  {
    //Redistribute the data
    //local pointer to head to rows in local CSR structure
    std::vector<Int> * row_nnz = new std::vector<Int>(numRowLocal,0);

    for(int p =0; p<mpisize;p++){
      std::vector< char > send_buffer;
      std::vector< char > recv_buffer;
      std::vector<int> displs(mpisize+1,0);
      std::vector<int> bufsizes(mpisize,0);

      //parse the Global CSC structure to count
      Int firstCol = p*numRowLocalFirst+1;
      Int lastCol = firstCol + numRowVec[p]-1;
      for(Int col = 1; col<colptrGlobal.m();++col){
        if(col>= firstCol && col<= lastCol){
          Int colbeg = colptrGlobal[col-1];
          Int colend = colptrGlobal[col]-1;
          for(Int i=colbeg;i<=colend;++i){
            Int row = rowindGlobal[i-1];
            //determine where it should be packed ?
            Int p_dest = std::min(mpisize-1,(row-1)/(numRowLocalFirst));
            bufsizes[p_dest]+=sizeof(T);
          }
        }
        else if(col>lastCol){
          break;
        }
      }

      //compute displacement (cum sum)
      std::partial_sum(bufsizes.begin(),bufsizes.end(),displs.begin()+1);

      if(mpirank==p){
        std::vector<int> bufpos(mpisize,0);

        //resize buffers
        Int send_buffer_size = displs[mpisize];
        send_buffer.resize(displs[mpisize]);

        //fill the buffers by parsing local CSC structure
        for(Int j = 1; j<sparseA.colptrLocal.m();++j){
          Int gCol = mpirank*numRowLocalFirst + j;
          Int colbeg = sparseA.colptrLocal[j-1];
          Int colend = sparseA.colptrLocal[j]-1;
          for(Int i=colbeg;i<=colend;++i){
            Int row = sparseA.rowindLocal[i-1];
            //determine where it should be packed ?
            Int p_dest = std::min(mpisize-1,(row-1)/(numRowLocalFirst));
            Int p_displs = displs[p_dest];
            Int p_bufpos = bufpos[p_dest];
            *((T *)&send_buffer.at( displs[p_dest] + bufpos[p_dest] )) = sparseA.nzvalLocal[i-1];
            //advance position
            bufpos[p_dest]+= sizeof(T);
          }
        }
      }



      Int recv_buffer_size = bufsizes[mpirank];
      recv_buffer.resize(bufsizes[mpirank]);

      //scatterv
      MPI_Scatterv(&send_buffer[0], &bufsizes[0], &displs[0], MPI_BYTE, &recv_buffer[0], bufsizes[mpirank], MPI_BYTE, p, sparseA.comm);

      //process data
      Int recv_pos = 0;
      //parse the Global CSC structure to count
      for(Int col = 1; col<colptrGlobal.m();++col){
        if(col>= firstCol && col<= lastCol){
          Int colbeg = colptrGlobal[col-1];
          Int colend = colptrGlobal[col]-1;
          for(Int i=colbeg;i<=colend;++i){
            Int row = rowindGlobal[i-1];
            Int p_dest = std::min(mpisize-1,(row-1)/(numRowLocalFirst));
            if(p_dest==mpirank){
              //compute local CSR coordinates
              Int local_row = row - mpirank*numRowLocalFirst;
              Int h = rowptrLocal[local_row-1] + row_nnz->at(local_row-1);    /* non-zero position */
              *((T*)&nzvalLocal[h-1]) = *((T*)&recv_buffer.at(recv_pos));
              //advance position in CSR buffer
              row_nnz->at(local_row-1)++;
              //advance position in CSC recv_buffer
              recv_pos+=sizeof(T);
            }
          }
        }
        else if(col>lastCol){
          break;
        }
      }
    }
    delete row_nnz;
  }
}



  template<typename T>
  inline std::string ToMatlabScalar( std::complex<T> val){
    std::stringstream s;
    s.precision(std::numeric_limits< std::complex<T> >::max_digits10);
    s.precision(15);
    s<<"complex("<<std::scientific<<std::real(val)<<","<<std::imag(val)<<")";
    return s.str();
  }

  template<typename T>
  inline std::string ToMatlabScalar( T val){
    std::stringstream s;
    s.precision(std::numeric_limits< T >::max_digits10);
    s<<std::scientific<<val;
    return s.str();
  }



template<typename T>
void WriteDistSparseMatrixMatlab(const char * filename, DistSparseMatrix<T> & pspmat, MPI_Comm comm){

  MPI_Barrier( comm );
  int mpirank=0;  MPI_Comm_rank(comm, &mpirank);
  int mpisize=0;  MPI_Comm_size(comm, &mpisize);
  std::string fname (filename);
  std::stringstream sstm;
  sstm<<fname<<"_"<<mpirank<<".m";
  std::ofstream ofile(sstm.str().c_str());

    if( !ofile.good() ){
      ErrorHandling( "File cannot be opened!" );
    }

    Int baseval = 1;

    Int numColFirst = pspmat.size / mpisize;
    Int firstCol    = mpirank * numColFirst;
    Int numColLocal = pspmat.colptrLocal.m() - baseval;

    Int firstLocCol = firstCol;
    Int LocalVertexCount = numColLocal;

    //I
    ofile<<"AinvExp = sparse([";
    for(Int locCol = 0 ; locCol< LocalVertexCount; locCol++){
      Int col = locCol + firstLocCol;
      Int colbeg = pspmat.colptrLocal[locCol]-baseval; //now 0 based
      Int colend = pspmat.colptrLocal[locCol+1]-baseval; // now 0 based

      for(Int pos = colbeg; pos<colend; pos++){
        Int row = pspmat.rowindLocal[pos];
        ofile<<row<<" ";
      }
    }
    ofile<<"],";//<<std::endl;

    //J
    ofile<<"[";
    for(Int locCol = 0 ; locCol< LocalVertexCount; locCol++){
      Int col = locCol + firstLocCol;
      Int colbeg = pspmat.colptrLocal[locCol]-baseval; //now 0 based
      Int colend = pspmat.colptrLocal[locCol+1]-baseval; // now 0 based

      for(Int pos = colbeg; pos<colend; pos++){
        ofile<<col+1<<" ";
      }
    }
    ofile<<"],";//<<std::endl;

    //V
    ofile<<"[";
    ofile.precision(std::numeric_limits< T >::max_digits10);
    for(Int locCol = 0 ; locCol< LocalVertexCount; locCol++){
      Int col = locCol + firstLocCol;
      Int colbeg = pspmat.colptrLocal[locCol]-baseval; //now 0 based
      Int colend = pspmat.colptrLocal[locCol+1]-baseval; // now 0 based

      for(Int pos = colbeg; pos<colend; pos++){
        T val = pspmat.nzvalLocal[pos];
        ofile<<std::scientific<<ToMatlabScalar(val)<<" ";
      }
    }
    ofile<<"]);"<<std::endl;

}

}
#endif //_PEXSI_UTILITY_IMPL_HPP_
