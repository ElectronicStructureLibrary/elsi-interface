/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

Authors: Lin Lin and Mathias Jacquelin

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
/// @file run_pselinv.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @date 2013-04-15
#include  "ppexsi.hpp"

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654


#define INFO(I) info[(I)-1]
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
#define KEEP(I) keep[(I)-1] /* macro s.t. indices match documentation */
#define CNTL(I) cntl[(I)-1] /* macro s.t. indices match documentation */


#include "pexsi/timer.h"

#define _MYCOMPLEX_

#ifdef _MYCOMPLEX_
#define MYSCALAR Complex
#include "zmumps_c.h"
#define MUMPS(a) zmumps_c(a)
#define MUMPS_STRUC_C ZMUMPS_STRUC_C
#else
#define MYSCALAR Real
#include "dmumps_c.h"
#define MUMPS(a) dmumps_c(a)
#define MUMPS_STRUC_C DMUMPS_STRUC_C
#endif


using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout << "Usage" << std::endl << "run_mumps -T [isText] -Inv [doInv] -Fact [0(LU)|1(LLT)|2(LDLT)] -H <Hfile> [-S <Sfile>] -Ord [AMD|...] -rshift [real shift] -ishift [imaginary shift] -Real [0|1]" << std::endl;
}

int main(int argc, char **argv) 
{
  MUMPS_STRUC_C id;

  if( argc < 3 ) {
    Usage();
    return 0;
  }

#if defined(PROFILE) || defined(PMPI)
  TAU_PROFILE_INIT(argc, argv);
#endif

  MPI_Init( &argc, &argv );
  int mpirank, mpisize;
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );


  try{
    MPI_Comm world_comm;

    // *********************************************************************
    // Input parameter
    // *********************************************************************
    std::map<std::string,std::string> options;

    OptionsCreate(argc, argv, options);

    // Default processor number
    Int nprow = 1;
    Int npcol = mpisize;

    Int doFact = 2;
    if( options.find("-Fact") != options.end() ){
      doFact = atoi(options["-Fact"].c_str());
    }

    Int isReal = 1;
    if( options.find("-Real") != options.end() ){ 
      isReal = atoi(options["-Real"].c_str());
    }

    Int doInv = 1;
    if( options.find("-Inv") != options.end() ){
      doInv = atoi(options["-Inv"].c_str());
    }

    Int icntl7 = 7;
    if( options.find("-Ord") != options.end() ){
      if(options["-Ord"] == "AMD"){
        icntl7=0;
      }
      else if(options["-Ord"] == "AMF"){
        icntl7=2;
      }
      else if(options["-Ord"] == "SCOTCH"){
        icntl7=3;
      }
    }

    nprow = mpisize;
    npcol=1;

    //    if( options.find("-r") != options.end() ){
    //      if( options.find("-c") != options.end() ){
    //        nprow= atoi(options["-r"].c_str());
    //        npcol= atoi(options["-c"].c_str());
    //        if(nprow*npcol > mpisize){
    //          ErrorHandling("The number of used processors cannot be higher than the total number of available processors." );
    //        } 
    //      }
    //      else{
    //        ErrorHandling( "When using -r option, -c also needs to be provided." );
    //      }
    //    }
    //    else if( options.find("-c") != options.end() ){
    //      if( options.find("-r") != options.end() ){
    //        nprow= atoi(options["-r"].c_str());
    //        npcol= atoi(options["-c"].c_str());
    //        if(nprow*npcol > mpisize){
    //          ErrorHandling("The number of used processors cannot be higher than the total number of available processors." );
    //        } 
    //      }
    //      else{
    //        ErrorHandling( "When using -c option, -r also needs to be provided." );
    //      }
    //    }

    //Create a communicator with npcol*nprow processors
    MPI_Comm_split(MPI_COMM_WORLD, mpirank<nprow*npcol, mpirank, &world_comm);

    if (mpirank<nprow*npcol){

      MPI_Comm_rank(world_comm, &mpirank );
      MPI_Comm_size(world_comm, &mpisize );
#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
      TAU_PROFILE_SET_CONTEXT(world_comm);
#endif

      stringstream  ss;
      ss << "logTest" << mpirank;
      statusOFS.open( ss.str().c_str() );

      //if( mpisize != nprow * npcol || nprow != npcol ){
      //  ErrorHandling( "nprow == npcol is assumed in this test routine." );
      //}

      if( mpirank == 0 )
        cout << "nprow = " << nprow << ", npcol = " << npcol << endl;

      std::string Hfile, Sfile;
      int isCSC = true;
      if( options.find("-T") != options.end() ){ 
        isCSC= ! atoi(options["-T"].c_str());
      }


      /* load A in coordinates version */
      /* a[k] is as row irn[k] and col jcn[k] */

      IntNumVec irn, jcn;
      //  int * rhs_irn, * rhs_jcn;
      //  struct sparse_matrix_t* A;
      //
      //  int nz_rhs;
      //
      //  double * rhs;
      //  int * irhs_sparse;
      //  int * irhs_ptr;


      int icntl14 = 20;
      if(options.find("-R") != options.end() ){
        icntl14 = atoi(options["-R"].c_str());
      }


      int icntl27 = 8;
      if(options.find("-B") != options.end() ){
        icntl27 = atoi(options["-B"].c_str());
      }

      int diagOnly = 0;
      if(options.find("-Pattern") != options.end() ){
        if(options["-Pattern"]=="DIAG"){
          diagOnly = 1;
        }
      }


      if( options.find("-H") != options.end() ){ 
        Hfile = options["-H"];
      }
      else{
        ErrorHandling("Hfile must be provided.");
      }

      if( options.find("-S") != options.end() ){ 
        Sfile = options["-S"];
      }
      else{
        statusOFS << "-S option is not given. " 
          << "Treat the overlap matrix as an identity matrix." 
          << std::endl << std::endl;
      }

      Real rshift = 0.0, ishift = 0.0;
      if( options.find("-rshift") != options.end() ){ 
        rshift = atof(options["-rshift"].c_str());
      }
      if( options.find("-ishift") != options.end() ){ 
        ishift = atof(options["-ishift"].c_str());
      }


      // *********************************************************************
      // Read input matrix
      // *********************************************************************

      // Setup grid.
      SuperLUGrid<MYSCALAR> g( world_comm, nprow, npcol );
      //      SuperLUGrid<Complex> g1( world_comm, nprow, npcol );

      int      m, n;
      DistSparseMatrix<MYSCALAR>  AMat;
      Real timeSta, timeEnd;


#ifdef _MYCOMPLEX_
      if(!isReal){
        DistSparseMatrix<Complex> HMat;
        DistSparseMatrix<Complex> SMat;
        GetTime( timeSta );
        if(isCSC){
          ParaReadDistSparseMatrix( Hfile.c_str(), HMat, world_comm ); 
        }
        else{
          ReadDistSparseMatrixFormatted( Hfile.c_str(), HMat, world_comm ); 
          ParaWriteDistSparseMatrix( "H.csc", HMat, world_comm ); 
        }

        if( Sfile.empty() ){
          // Set the size to be zero.  This will tell PPEXSI.Solve to treat
          // the overlap matrix as an identity matrix implicitly.
          SMat.size = 0;  
        }
        else{
          if(isCSC){
            ParaReadDistSparseMatrix( Sfile.c_str(), SMat, world_comm ); 
          }
          else{
            ReadDistSparseMatrixFormatted( Sfile.c_str(), SMat, world_comm ); 
            ParaWriteDistSparseMatrix( "S.csc", SMat, world_comm ); 
          }
        }

        GetTime( timeEnd );
        LongInt nnzH = HMat.Nnz();
        if( mpirank == 0 ){
          cout << "Time for reading H and S is " << timeEnd - timeSta << endl;
          cout << "H.size = " << HMat.size << endl;
          cout << "H.nnz  = " << nnzH  << endl;
        }

        // Get the diagonal indices for H and save it n diagIdxLocal_

        std::vector<Int>  diagIdxLocal;
        { 
          Int numColLocal      = HMat.colptrLocal.m() - 1;
          Int numColLocalFirst = HMat.size / mpisize;
          Int firstCol         = mpirank * numColLocalFirst;

          diagIdxLocal.clear();

          for( Int j = 0; j < numColLocal; j++ ){
            Int jcol = firstCol + j + 1;
            for( Int i = HMat.colptrLocal(j)-1; 
                i < HMat.colptrLocal(j+1)-1; i++ ){
              Int irow = HMat.rowindLocal(i);
              if( irow == jcol ){
                diagIdxLocal.push_back( i );
              }
            }
          } // for (j)
        }


        GetTime( timeSta );

        AMat.size          = HMat.size;
        AMat.nnz           = HMat.nnz;
        AMat.nnzLocal      = HMat.nnzLocal;
        AMat.colptrLocal   = HMat.colptrLocal;
        AMat.rowindLocal   = HMat.rowindLocal;
        AMat.nzvalLocal.Resize( HMat.nnzLocal );
        AMat.comm = world_comm;

        MYSCALAR *ptr0 = AMat.nzvalLocal.Data();
        Complex *ptr1 = HMat.nzvalLocal.Data();
        Complex *ptr2 = SMat.nzvalLocal.Data();

#ifdef _MYCOMPLEX_
        Complex zshift = Complex(rshift, ishift);
#else
        Real zshift = Real(rshift);
#endif

        if( SMat.size != 0 ){
          // S is not an identity matrix
          for( Int i = 0; i < HMat.nnzLocal; i++ ){
            AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift * SMat.nzvalLocal(i);
          }
        }
        else{
          // S is an identity matrix
          for( Int i = 0; i < HMat.nnzLocal; i++ ){
            AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
          }

          for( Int i = 0; i < diagIdxLocal.size(); i++ ){
            AMat.nzvalLocal( diagIdxLocal[i] ) -= zshift;
          }
        } // if (SMat.size != 0 )
      }
      else
#endif
      {
        DistSparseMatrix<Real> HMat;
        DistSparseMatrix<Real> SMat;
        Real timeSta, timeEnd;
        GetTime( timeSta );
        if(isCSC){
          ParaReadDistSparseMatrix( Hfile.c_str(), HMat, world_comm ); 
        }
        else{
          ReadDistSparseMatrixFormatted( Hfile.c_str(), HMat, world_comm ); 
          ParaWriteDistSparseMatrix( "H.csc", HMat, world_comm ); 
        }

        if( Sfile.empty() ){
          // Set the size to be zero.  This will tell PPEXSI.Solve to treat
          // the overlap matrix as an identity matrix implicitly.
          SMat.size = 0;  
        }
        else{
          if(isCSC){
            ParaReadDistSparseMatrix( Sfile.c_str(), SMat, world_comm ); 
          }
          else{
            ReadDistSparseMatrixFormatted( Sfile.c_str(), SMat, world_comm ); 
            ParaWriteDistSparseMatrix( "S.csc", SMat, world_comm ); 
          }
        }

        GetTime( timeEnd );
        LongInt nnzH = HMat.Nnz();
        if( mpirank == 0 ){
          cout << "Time for reading H and S is " << timeEnd - timeSta << endl;
          cout << "H.size = " << HMat.size << endl;
          cout << "H.nnz  = " << nnzH  << endl;
        }

        // Get the diagonal indices for H and save it n diagIdxLocal_

        std::vector<Int>  diagIdxLocal;
        { 
          Int numColLocal      = HMat.colptrLocal.m() - 1;
          Int numColLocalFirst = HMat.size / mpisize;
          Int firstCol         = mpirank * numColLocalFirst;

          diagIdxLocal.clear();

          for( Int j = 0; j < numColLocal; j++ ){
            Int jcol = firstCol + j + 1;
            for( Int i = HMat.colptrLocal(j)-1; 
                i < HMat.colptrLocal(j+1)-1; i++ ){
              Int irow = HMat.rowindLocal(i);
              if( irow == jcol ){
                diagIdxLocal.push_back( i );
              }
            }
          } // for (j)
        }


        GetTime( timeSta );

        AMat.size          = HMat.size;
        AMat.nnz           = HMat.nnz;
        AMat.nnzLocal      = HMat.nnzLocal;
        AMat.colptrLocal   = HMat.colptrLocal;
        AMat.rowindLocal   = HMat.rowindLocal;
        AMat.nzvalLocal.Resize( HMat.nnzLocal );
        AMat.comm = world_comm;

        MYSCALAR *ptr0 = AMat.nzvalLocal.Data();
        Real *ptr1 = HMat.nzvalLocal.Data();
        Real *ptr2 = SMat.nzvalLocal.Data();

#ifdef _MYCOMPLEX_
        Complex zshift = Complex(rshift, ishift);
#else
        Real zshift = Real(rshift);
#endif

        if( SMat.size != 0 ){
          // S is not an identity matrix
          for( Int i = 0; i < HMat.nnzLocal; i++ ){
            AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift * SMat.nzvalLocal(i);
          }
        }
        else{
          // S is an identity matrix
          for( Int i = 0; i < HMat.nnzLocal; i++ ){
            AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
          }

          for( Int i = 0; i < diagIdxLocal.size(); i++ ){
            AMat.nzvalLocal( diagIdxLocal[i] ) -= zshift;
          }
        } // if (SMat.size != 0 )
      }












      /*
         DistSparseMatrix<Real> HMat;
         DistSparseMatrix<Real> SMat;
         GetTime( timeSta );
         if(isCSC)
         ParaReadDistSparseMatrix( Hfile.c_str(), HMat, world_comm ); 
         else{
         ReadDistSparseMatrixFormatted( Hfile.c_str(), HMat, world_comm ); 
         ParaWriteDistSparseMatrix( "H.csc", HMat, world_comm ); 
         }

         if( Sfile.empty() ){
      // Set the size to be zero.  This will tell PPEXSI.Solve to treat
      // the overlap matrix as an identity matrix implicitly.
      SMat.size = 0;  
      }
      else{
      if(isCSC)
      ParaReadDistSparseMatrix( Sfile.c_str(), SMat, world_comm ); 
      else{
      ReadDistSparseMatrixFormatted( Sfile.c_str(), SMat, world_comm ); 
      ParaWriteDistSparseMatrix( "S.csc", SMat, world_comm ); 
      }
      }

      GetTime( timeEnd );
      LongInt nnzH = HMat.Nnz();
      if( mpirank == 0 ){
      cout << "Time for reading H and S is " << timeEnd - timeSta << endl;
      cout << "H.size = " << HMat.size << endl;
      cout << "H.nnz  = " << nnzH  << endl;
      }

      // Get the diagonal indices for H and save it n diagIdxLocal_

      std::vector<Int>  diagIdxLocal;
      { 
      Int numColLocal      = HMat.colptrLocal.m() - 1;
      Int numColLocalFirst = HMat.size / mpisize;
      Int firstCol         = mpirank * numColLocalFirst;

      diagIdxLocal.clear();

      for( Int j = 0; j < numColLocal; j++ ){
      Int jcol = firstCol + j + 1;
      for( Int i = HMat.colptrLocal(j)-1; 
      i < HMat.colptrLocal(j+1)-1; i++ ){
      Int irow = HMat.rowindLocal(i);
      if( irow == jcol ){
      diagIdxLocal.push_back( i );
      }
      }
      } // for (j)
      }


      GetTime( timeSta );

      AMat.size          = HMat.size;
      AMat.nnz           = HMat.nnz;
      AMat.nnzLocal      = HMat.nnzLocal;
      AMat.colptrLocal   = HMat.colptrLocal;
      AMat.rowindLocal   = HMat.rowindLocal;
      AMat.nzvalLocal.Resize( HMat.nnzLocal );
      AMat.comm = world_comm;

      MYSCALAR *ptr0 = AMat.nzvalLocal.Data();
      Real *ptr1 = HMat.nzvalLocal.Data();
      Real *ptr2 = SMat.nzvalLocal.Data();

#ifdef _MYCOMPLEX_
      Complex zshift = Complex(rshift, ishift);
#else
      Real zshift = Real(rshift);
#endif

      if( SMat.size != 0 ){
        // S is not an identity matrix
        for( Int i = 0; i < HMat.nnzLocal; i++ ){
          AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift * SMat.nzvalLocal(i);
        }
      }
      else{
        // S is an identity matrix
        for( Int i = 0; i < HMat.nnzLocal; i++ ){
          AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
        }

        for( Int i = 0; i < diagIdxLocal.size(); i++ ){
          AMat.nzvalLocal( diagIdxLocal[i] ) -= zshift;
        }
      } // if (SMat.size != 0 )
      */




        LongInt nnzA = AMat.Nnz();
      if( mpirank == 0 ){
        cout << "nonzero in A (DistSparseMatrix format) = " << nnzA << endl;
      }


      GetTime( timeEnd );

      // Compute the number of columns on each processor
      Int numColFirst;
      numColFirst = AMat.size / mpisize;


      if( mpirank == 0 )
        cout << "Time for constructing the matrix A is " << timeEnd - timeSta << endl;



      //statusOFS<<"AMat local: "<<AMat.nzvalLocal<<endl;

#if ( _DEBUGlevel_ >= 2 )
      statusOFS<<"AMat local: "<<endl;
      for(Int j =0; j<AMat.colptrLocal.m()-1; ++j){
        Int col = j +1 + (mpirank)*numColFirst; 
        for(Int i = AMat.colptrLocal[j]; i<AMat.colptrLocal[j+1]; ++i){
          Int row = AMat.rowindLocal[i-1];
          MYSCALAR val = AMat.nzvalLocal[i-1];
          //          statusOFS<<"("<<row<<","<<col<<")="<<val<<" ";
          statusOFS<<"("<<row<<","<<col<<") ";
        } 
      }
      statusOFS<<endl;
#endif


      IntNumVec irhs_ptr, irhs_sparse;

#ifdef _MYCOMPLEX_
      std::vector<mumps_double_complex> rhs_sparse;
#else
      NumVec<MYSCALAR> rhs_sparse;
#endif

      //MUMPS conversion
      irn.Resize(AMat.nnzLocal);
      jcn.Resize(AMat.nnzLocal);
      //extract the coordinates
      Int index = 0;
      for(Int j =0; j<AMat.colptrLocal.m()-1; ++j){
        Int col = j +1 + (mpirank)*numColFirst; 
        for(Int i = AMat.colptrLocal[j]; i<AMat.colptrLocal[j+1]; ++i){
          Int row = AMat.rowindLocal[i-1];
          irn[index] = row;
          jcn[index] = col;
          index++;
        } 
      }


      //Get the structure of the matrix to compute the sparse rhs on the host
      //      if(!diagOnly){
      IntNumVec numColLocalVec(mpisize);
      IntNumVec displs(mpisize);
      IntNumVec displsColptr(mpisize);

      // Compute the number of columns on each processor
      numColLocalVec.Resize(mpisize);
      Int numColLocal;
      numColFirst = AMat.size / mpisize;
      SetValue( numColLocalVec, numColFirst );
      numColLocalVec[mpisize-1] = AMat.size - numColFirst * (mpisize-1);  // Modify the last entry	
      numColLocal = numColLocalVec[mpirank];

      displsColptr[0]=0;
      for(int p = 1;p<mpisize;++p){
        displsColptr[p]= displsColptr[p-1]+ numColLocalVec[p-1];
      }

      for(int p = 0;p<mpisize;++p){
        numColLocalVec[p]*=sizeof(Int);
        displsColptr[p]*=sizeof(Int);
      }

      //algather nzvalLocal
      IntNumVec nnzLocalVec(mpisize);
      MPI_Gather(&AMat.nnzLocal,sizeof(Int),MPI_BYTE,&nnzLocalVec[0],sizeof(Int),MPI_BYTE,0,world_comm);
      //compute displacements 
      displs[0]=0;
      for(int p = 1;p<mpisize;++p){
        displs[p]= displs[p-1]+ nnzLocalVec[p-1];
      }

      for(int p = 0;p<mpisize;++p){
        nnzLocalVec[p]*=sizeof(Int);
        displs[p]*=sizeof(Int);
      }



      if(mpirank==0){
        irhs_ptr.Resize(AMat.size+1);
        irhs_sparse.Resize(nnzA);
#ifdef _MYCOMPLEX_
        rhs_sparse.resize(nnzA);
        for(int i = 0;i<rhs_sparse.size();++i){rhs_sparse[i].r=0.0;rhs_sparse[i].i=0.0;}
#else
        rhs_sparse.Resize(nnzA);
        SetValue(rhs_sparse,0.0);
#endif
      }



      //do the same for rowind
#if ( _DEBUGlevel_ >= 2 )
      statusOFS<<"RowindLocal : "<<AMat.rowindLocal<<endl;
#endif
      if(mpirank==0){
        MPI_Gatherv(&AMat.rowindLocal[0],AMat.nnzLocal*sizeof(Int),MPI_BYTE,&irhs_sparse[0],&nnzLocalVec[0],&displs[0],MPI_BYTE,0,world_comm);
      }
      else{
        MPI_Gatherv(&AMat.rowindLocal[0],AMat.nnzLocal*sizeof(Int),MPI_BYTE,NULL,NULL,NULL,MPI_BYTE,0,world_comm);
      }

#if ( _DEBUGlevel_ >= 2 )
      if(mpirank==0){
        statusOFS<<"irhs_sparse : "<<irhs_sparse<<endl;
      }
#endif

      for(int p = 0;p<mpisize;++p){
        nnzLocalVec[p]/=sizeof(Int);
        displs[p]/=sizeof(Int);
      }





#if ( _DEBUGlevel_ >= 2 )
      statusOFS<<"Colptr Local: "<<AMat.colptrLocal<<endl;
#endif

      //gatherv the colptrs
      if(mpirank==0){
        MPI_Gatherv(&AMat.colptrLocal[0],numColLocalVec[mpirank],MPI_BYTE,&irhs_ptr[0],&numColLocalVec[0],&displsColptr[0],MPI_BYTE,0,world_comm);
      }
      else{
        MPI_Gatherv(&AMat.colptrLocal[0],numColLocalVec[mpirank],MPI_BYTE,NULL,NULL,NULL,MPI_BYTE,0,world_comm);
      }

      if(mpirank==0){
        Int index = 0;

        for(int p=0;p<mpisize;++p){
          numColLocalVec[p] /= sizeof(Int); 
          displsColptr[p] /= sizeof(Int); 
        }

#if ( _DEBUGlevel_ >= 2 )
        statusOFS<<"NumColLocalVec: "<<numColLocalVec<<endl;
        statusOFS<<"nnzLocalVec: "<<nnzLocalVec<<endl;
        statusOFS<<"irhs_ptr before: "<<irhs_ptr<<endl;
#endif
        for(int p=1;p<mpisize;++p){
          Int * colptr = &irhs_ptr[0] + displsColptr[p];
          Int offset = displs[p]+1; 
          if(numColLocalVec[p]>0){
            for(int j = numColLocalVec[p]-1; j>=0;--j){
              colptr[j] += offset - colptr[0]; 
            }
          }
        }
        irhs_ptr[AMat.size]=nnzA+1;
#if ( _DEBUGlevel_ >= 2 )
        statusOFS<<"irhs_ptr after: "<<irhs_ptr<<endl;
#endif
      }


      if(doFact!=0){
        for(Int j =0; j<irhs_ptr.m()-1; ++j){
          Int col = j+1;
          for(Int i = irhs_ptr[j]; i<irhs_ptr[j+1];++i){
            Int row = irhs_sparse[i-1];
            if(row<col){cerr<<"error at "<<j<<" "<<i-1<<" "<<row<<" vs "<<col<<endl;}
            assert(row>=col);
          }
        }
      }

      //      }
      if(diagOnly){
        rhs_sparse.clear();
        rhs_sparse.resize(AMat.size);

        irhs_ptr.Resize(AMat.size+1);
        for(Int i = 0; i<=AMat.size;++i){ irhs_ptr[i] = i+1; }

        irhs_sparse.Clear();
        irhs_sparse.Resize(AMat.size);
        for(Int i = 0; i<AMat.size;++i){ irhs_sparse[i] = i+1; }
      }

      MPI_Barrier(world_comm);



      /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
      id.job=JOB_INIT; 
      id.par=1; 
      id.sym = doFact;
      id.comm_fortran=world_comm;
      MUMPS(&id);

      id.n = AMat.size; 
      id.nz = nnzA; 
      id.nz_loc = AMat.nnzLocal;
      id.irn_loc = &irn[0];
      id.jcn_loc = &jcn[0];


#if ( _DEBUGlevel_ >= 2 )
      statusOFS<<"Local AMat "<<endl;
      for(Int idx=0; idx<AMat.nnzLocal; ++idx){
        Int row = id.irn_loc[idx];
        Int col = id.jcn_loc[idx];
        MYSCALAR val = AMat.nzvalLocal[idx];
        //          statusOFS<<"("<<row<<","<<col<<")="<<val<<" ";
        statusOFS<<"("<<row<<","<<col<<") ";
      }
      statusOFS<<endl;
#endif



#ifdef _MYCOMPLEX_
      id.a_loc = new mumps_double_complex[AMat.nnzLocal];
      for(Int i =0; i<AMat.nzvalLocal.m();++i){id.a_loc[i].r=std::real(AMat.nzvalLocal[i]);id.a_loc[i].i=std::imag(AMat.nzvalLocal[i]);}
#else
      id.a_loc = &AMat.nzvalLocal[0];
#endif

      id.ICNTL(18)=3;
      id.ICNTL(7)=icntl7;
      id.ICNTL(14) = icntl14;
      //id.KEEP(219)=0;
      //for mumps 5.0beta magic keep for inversion
      //id.KEEP(243)=1;
      //id.KEEP(497)=1;
      //id.KEEP(495)=1;

      GetTime( timeSta );
      id.job=1;
      /* call the mumps package for symbfact*/
      MUMPS(&id);
      GetTime( timeEnd );
      if (mpirank == 0) {
        cout<<"Symbolic factorization done in: "<<timeEnd-timeSta<<"s"<<std::endl;
      }

      GetTime( timeSta );
      id.job=2;
      MUMPS(&id);
      GetTime( timeEnd );
      if (mpirank == 0) {
        cout<<"Factorization done in: "<<timeEnd-timeSta<<"s"<<std::endl;
      }

      if(doInv){


        //proceed with inversion
        if(mpirank==0){
          if(diagOnly){
            id.nz_rhs = AMat.size;
            id.nrhs =AMat.size; 
            id.rhs_sparse = &rhs_sparse[0];
            id.irhs_sparse = &irhs_sparse[0];
            id.irhs_ptr = &irhs_ptr[0];
          }
          else{
            id.nz_rhs = nnzA;
            id.nrhs =AMat.size; 
            id.rhs_sparse = &rhs_sparse[0];
            id.irhs_sparse = &irhs_sparse[0];
            id.irhs_ptr = &irhs_ptr[0];
          }
        }



        GetTime( timeSta );
        id.ICNTL(11)=1;
        /* ask for inversion*/
        id.ICNTL(30)=1;
        //block size for the solve, ICNTL(27) is the parameter to be tweaked
        id.ICNTL(27) = icntl27; 
        /* Call the MUMPS package. */
        id.job=3;
        MUMPS(&id);
        GetTime( timeEnd );
        if (mpirank == 0) {
          cout<<"Inversion done in: "<<timeEnd-timeSta<<"s"<<std::endl;
        }


#ifdef _MYCOMPLEX_
        delete [] id.a_loc;
#endif


        if(!diagOnly){


          //build a distsparsematrix
          DistSparseMatrix<MYSCALAR> Ainv;
          Ainv.colptrLocal = AMat.colptrLocal;
          Ainv.rowindLocal = AMat.rowindLocal;
          Ainv.size = AMat.size;
          Ainv.nnz = AMat.nnz;
          Ainv.nnzLocal = AMat.nnzLocal;
          Ainv.nzvalLocal.Resize(Ainv.nnzLocal);
          Ainv.comm = AMat.comm;



#if ( _DEBUGlevel_ >= 2 )
          if(mpirank==0){
            statusOFS<<"nnzLocalVec: "<<nnzLocalVec<<endl;
            statusOFS<<"displs: "<<displs<<endl;
            statusOFS<<"Ainv global: "<<rhs_sparse.size()<<endl;
            Int curP = 1;
            int index=0;
            for(Int j =0; j<irhs_ptr.m()-1; ++j){
              Int col = j+1;
              for(Int i = irhs_ptr[j]; i<irhs_ptr[j+1];++i){
                Int row = irhs_sparse[i-1];
#ifdef _MYCOMPLEX_
                MYSCALAR val(rhs_sparse[i-1].r,rhs_sparse[i-1].i);
#else
                MYSCALAR val = rhs_sparse[i-1];
#endif
                statusOFS<<"("<<row<<","<<col<<")="<<val<<" ";
                index++;
                if(curP<mpisize){if(index==displs[curP]){++curP;statusOFS<<endl;}}
              }
            }
            statusOFS<<endl;
          }
#endif

#ifdef _MYCOMPLEX_
          std::vector<mumps_double_complex> recvAinv(Ainv.nnzLocal);
#else
          NumVec<MYSCALAR> & recvAinv = Ainv.nzvalLocal;
#endif

          Int sizeElem = sizeof(recvAinv[0]);
          for(int p = 0;p<mpisize;++p){
            nnzLocalVec[p]*=sizeElem;
            displs[p]*=sizeElem;
          }

          //Now do a scatterv from P0 to distribute the values
          if(mpirank==0){
            MPI_Scatterv(&rhs_sparse[0],&nnzLocalVec[0],&displs[0],MPI_BYTE,&recvAinv[0],Ainv.nnzLocal*sizeElem,MPI_BYTE,0,world_comm);
          }
          else{
            MPI_Scatterv(NULL,NULL,NULL,MPI_BYTE,&recvAinv[0],Ainv.nnzLocal*sizeElem,MPI_BYTE,0,world_comm);
          }


#ifdef _MYCOMPLEX_
          for(Int i =0; i<Ainv.nzvalLocal.m();++i){  Ainv.nzvalLocal[i] = MYSCALAR( recvAinv[i].r, recvAinv[i].i);}
          recvAinv.clear();
#endif

#if ( _DEBUGlevel_ >= 2 )
          statusOFS<<"Ainv local: "<<endl;
          for(Int j =0; j<Ainv.colptrLocal.m()-1; ++j){
            Int col = j +1 + (mpirank)*numColFirst; 
            for(Int i = Ainv.colptrLocal[j]; i<Ainv.colptrLocal[j+1];++i){
              Int row = Ainv.rowindLocal[i-1];
#ifdef _MYCOMPLEX_
              MYSCALAR val(Ainv.nzvalLocal[i-1].r,Ainv.nzvalLocal[i-1].i);
#else
              MYSCALAR val = Ainv.nzvalLocal[i-1];
#endif
              statusOFS<<"("<<row<<","<<col<<")="<<val<<" ";
              //            statusOFS<<"("<<row<<","<<col<<") ";
            }
          }
          statusOFS<<endl;
#endif



#ifdef _MYCOMPLEX_
          if(mpirank==0){
            rhs_sparse.clear();
          }
#endif


          ParaWriteDistSparseMatrix( "Ainv_Mumps.csc", Ainv, world_comm ); 
        }

      }

      id.job=JOB_END; MUMPS(&id); /* Terminate instance */

      statusOFS.close();
    }
  }
  catch( std::exception& e )
  {
    std::cerr << "Processor " << mpirank << " caught exception with message: "
      << e.what() << std::endl;
  }

  MPI_Finalize();

  return 0;
}
