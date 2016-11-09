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

#include "pexsi/timer.h"

#define _MYCOMPLEX_

#ifdef _MYCOMPLEX_
#define MYSCALAR Complex
#else
#define MYSCALAR Real
#endif


using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout << "Usage" << std::endl << "run_pselinv -T [isText] -F [doFacto -E [doTriSolve] -Sinv [doSelInv]]  -H <Hfile> -S [Sfile] -colperm [colperm] -r [nprow] -c [npcol] -npsymbfact [npsymbfact] -P [maxpipelinedepth] -SinvBcast [doSelInvBcast] -SinvPipeline [doSelInvPipeline] -SinvHybrid [doSelInvHybrid] -rshift [real shift] -ishift [imaginary shift] -ToDist [doToDist] -Diag [doDiag]" << std::endl;
}

int main(int argc, char **argv) 
{

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

    if( options.find("-r") != options.end() ){
      if( options.find("-c") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
        if(nprow*npcol > mpisize){
          ErrorHandling("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        ErrorHandling( "When using -r option, -c also needs to be provided." );
      }
    }
    else if( options.find("-c") != options.end() ){
      if( options.find("-r") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
        if(nprow*npcol > mpisize){
          ErrorHandling("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        ErrorHandling( "When using -c option, -r also needs to be provided." );
      }
    }

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

#ifdef COMM_PROFILE
      stringstream  ss3;
      ss3 << "comm_stat" << mpirank;
      commOFS.open( ss3.str().c_str());
#endif

      //if( mpisize != nprow * npcol || nprow != npcol ){
      //  ErrorHandling( "nprow == npcol is assumed in this test routine." );
      //}

      if( mpirank == 0 )
        cout << "nprow = " << nprow << ", npcol = " << npcol << endl;


      GridType * g1Ptr;

      g1Ptr = new GridType( world_comm, nprow, npcol );
      GridType &g1 = *g1Ptr;


      //Try the reduction tree
      set<Int> set_ranks;
      for(int i = 1; i<nprow;++i){
        set_ranks.insert(i);
      }

      vector<Int> tree_ranks;
      tree_ranks.push_back(0);
      tree_ranks.insert(tree_ranks.end(),set_ranks.begin(),set_ranks.end());


      Int msgSize = sizeof(MYSCALAR)*nprow*nprow;
#if 1
      TreeBcast2<MYSCALAR> * redDTree;
      double rseed = 0;
      redDTree = TreeBcast2<MYSCALAR>::Create(g1Ptr->colComm,&tree_ranks[0],tree_ranks.size(),msgSize,rseed);

#ifdef COMM_PROFILE
      redDTree->SetGlobalComm(g1Ptr->comm);
#endif

      NumMat<MYSCALAR> buffer(nprow,nprow);
      Int myRow = MYROW(g1Ptr);

      for(int r=0;r<1;++r){

        //if(mpirank==15){i
        //gdb_lock();
        //}

        //bcastUTree2->AllocRecvBuffer();
        if(myRow==0){
          redDTree->SetLocalBuffer(buffer.Data());
          for(int i =0;i<nprow;++i){buffer(myRow,i)=myRow; buffer(i,i)=1;}
          redDTree->SetDataReady(true);
        }
        redDTree->SetTag( 7 );
        //Initialize the tree
        //redDTree->AllocRecvBuffers();
        //Post All Recv requests;
        //       redDTree->PostRecv();



        bool done = redDTree->Progress();


        //try a blocking wait
        redDTree->Wait();

        //      MPI_Barrier(g1Ptr->colComm);

        //root dumps the results
        //if(myRow == 0){
        statusOFS<<"DEBUG"<<endl;
        MYSCALAR * res = redDTree->GetLocalBuffer();
        for(int i =0;i<nprow;++i){
          for(int j =0;j<nprow;++j){
            statusOFS<<res[j*nprow+i]<<" ";
          }
          statusOFS<<endl;
        }
        statusOFS<<endl;
        redDTree->SetLocalBuffer(buffer.Data());
        statusOFS<<buffer<<endl;
        //}

        redDTree->CleanupBuffers();

      }
#else
      TreeReduce<MYSCALAR> * redDTree;
      redDTree = TreeReduce<MYSCALAR>::Create(g1Ptr->colComm,&tree_ranks[0],tree_ranks.size(),msgSize);
#ifdef COMM_PROFILE
      redDTree->SetGlobalComm(g1Ptr->comm);
#endif

      NumMat<MYSCALAR> buffer(nprow,nprow);
      Int myRow = MYROW(g1Ptr);

      for(int r=0;r<1;++r){

        for(int i =0;i<nprow;++i){buffer(myRow,i)=myRow; buffer(i,i)=1;}

        redDTree->SetLocalBuffer(buffer.Data());
        redDTree->SetDataReady(true);
        redDTree->SetTag( 7 );
        //Initialize the tree
        redDTree->AllocRecvBuffers();
        //Post All Recv requests;
        redDTree->PostFirstRecv();




        bool done = redDTree->Progress();


        //try a blocking wait
        redDTree->Wait();

        //      MPI_Barrier(g1Ptr->colComm);

        //root dumps the results
        if(myRow == 0){
          statusOFS<<"DEBUG"<<endl;
          MYSCALAR * res = redDTree->GetLocalBuffer();
          for(int i =0;i<nprow;++i){
            for(int j =0;j<nprow;++j){
              statusOFS<<res[j*nprow+i]<<" ";
            }
            statusOFS<<endl;
          }
          statusOFS<<endl;

          statusOFS<<buffer<<endl;
        }

        redDTree->CleanupBuffers();

      }
#endif


      delete redDTree;

      delete g1Ptr;
#ifdef COMM_PROFILE
      //std::cout<<"DUMPING COMM STATS "<<comm_stat.size()<<" "<<std::endl;
      commOFS<<HEADER_COMM<<std::endl;
      for(auto it = comm_stat.begin(); it!=comm_stat.end(); it+=4){
        commOFS<<LINE_COMM(it)<<std::endl;
      }
#endif


#ifdef COMM_PROFILE
      commOFS.close();
#endif
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
