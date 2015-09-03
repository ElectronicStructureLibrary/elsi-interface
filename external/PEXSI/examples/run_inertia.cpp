/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

   Authors: Lin Lin

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
/// @file run_inertia.cpp
/// @brief Test routine for computing a series of inertia for a matrix
/// stencil. This can be used for estimating the density of states in a
/// local region, when diagonalization is too expensive but
/// factorization is still possible.
/// @date 2014-01-24
#include  "ppexsi.hpp"
#include "pexsi/timer.h"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout << "Usage" << std::endl << "run_inertia -T [isText]  -H <Hfile> -S [Sfile] -ordering [ordering] -r [nprow] -c [npcol] -npsymbfact [npsymbfact]  -muMin [muMin] -muMax [muMax] -numShift [numShift] -npPerShift [npPerShift]" << std::endl;
}

int main(int argc, char **argv) 
{

  if( argc < 3 ) {
    Usage();
    return 0;
  }


  MPI_Init( &argc, &argv );
  int mpirank, mpisize;
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );


  try{
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
          throw std::runtime_error("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        throw std::runtime_error( "When using -r option, -c also needs to be provided." );
      }
    }
    else if( options.find("-c") != options.end() ){
      if( options.find("-r") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
        if(nprow*npcol > mpisize){
          throw std::runtime_error("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        throw std::runtime_error( "When using -c option, -r also needs to be provided." );
      }
    }

		Real muMin;
		if( options.find("-muMin") != options.end() ){
			muMin = std::atof(options["-muMin"].c_str());
		}
		else{
      throw std::logic_error("muMin must be provided.");
		}

		Real muMax;
		if( options.find("-muMax") != options.end() ){
			muMax = std::atof(options["-muMax"].c_str());
		}
		else{
      throw std::logic_error("muMax must be provided.");
		}

		Int numShift;
    if( options.find("-numShift") != options.end() ){
			numShift = std::atof(options["-numShift"].c_str());
		}
		else{
      throw std::logic_error("numShift must be provided.");
		}

		Int npPerShift;
    if( options.find("-npPerShift") != options.end() ){
			npPerShift = std::atoi(options["-npPerShift"].c_str());
		}
		else{
      throw std::logic_error("npPerShift must be provided.");
		}

    std::string Hfile, Sfile;
    int isCSC = true;
    if( options.find("-T") != options.end() ){ 
      isCSC= ! atoi(options["-T"].c_str());
    }

		if( options.find("-H") != options.end() ){ 
			Hfile = options["-H"];
		}
		else{
      throw std::logic_error("Hfile must be provided.");
		}

		if( options.find("-S") != options.end() ){ 
			Sfile = options["-S"];
		}
		else{
      if( mpirank == 0 ){
        std::cout << "-S option is not given. " 
          << "Treat the overlap matrix as an identity matrix." 
          << std::endl << std::endl;
      }
		}


    Int numProcSymbFact;
    if( options.find("-npsymbfact") != options.end() ){ 
      numProcSymbFact = atoi( options["-npsymbfact"].c_str() );
    }
    else{
      if( mpirank == 0 ){
        std::cout << "-npsymbfact option is not given. " 
          << "Use default value (maximum number of procs)." 
          << std::endl << std::endl;
      }
      numProcSymbFact = 0;
    }

    Int ordering;
    if( options.find("-ordering") != options.end() ){ 
      ordering = atoi( options["-ordering"].c_str() );
    }
    else{
      if( mpirank == 0 ){
        std::cout << "-ordering option is not given. " 
          << "Use default value 0 (PARMETIS)." 
          << std::endl << std::endl;
      }
      ordering = 0;
    }
    

		// *********************************************************************
		// Check the input parameters
		// *********************************************************************
		if( mpisize % npPerShift != 0 ){
			throw std::logic_error( "mpisize cannot be divided evenly by npPerShift." );
		}

		if( npPerShift != nprow * npcol  ){
			throw std::runtime_error( "npPerShift should be equal to nprow * npcol." );
		}

    // *********************************************************************
    // Read input matrix
    // *********************************************************************

    Int isProcRead;
    MPI_Comm  readComm;
    Real timeSta, timeEnd;


    /* Split the processors to read matrix */
    if( mpirank < npPerShift )
      isProcRead = 1;
    else
      isProcRead = 0;

    MPI_Comm_split( MPI_COMM_WORLD, isProcRead, mpirank, &readComm );
   
    DistSparseMatrix<Real> HMat;
    DistSparseMatrix<Real> SMat;

    if( isProcRead == 1 ){
      GetTime( timeSta );
      if(isCSC)
        ParaReadDistSparseMatrix( Hfile.c_str(), HMat, readComm ); 
      else{
        ReadDistSparseMatrixFormatted( Hfile.c_str(), HMat, readComm ); 
        ParaWriteDistSparseMatrix( "H.csc", HMat, readComm ); 
      }

      if( Sfile.empty() ){
        // Set the size to be zero.  This will tell PPEXSI.Solve to treat
        // the overlap matrix as an identity matrix implicitly.
        SMat.size = 0;  
      }
      else{
        if(isCSC)
          ParaReadDistSparseMatrix( Sfile.c_str(), SMat, readComm ); 
        else{
          ReadDistSparseMatrixFormatted( Sfile.c_str(), SMat, readComm ); 
          ParaWriteDistSparseMatrix( "S.csc", SMat, readComm ); 
        }
      }

      GetTime( timeEnd );
      LongInt nnzH = HMat.Nnz();
      if( mpirank == 0 ){
        cout << "Time for reading H and S is " << timeEnd - timeSta << endl;
        cout << "H.size = " << HMat.size << endl;
        cout << "H.nnz  = " << nnzH  << endl;
      }
    }

    GetTime( timeSta );

    DblNumVec  shiftList( numShift );
    IntNumVec  inertiaListInt( numShift );

    Int info;

    PPEXSIRawInertiaCountInterface(
        HMat.size,
        HMat.nnz,
        HMat.nnzLocal,
        HMat.colptrLocal.m() - 1,
        HMat.colptrLocal.Data(),
        HMat.rowindLocal.Data(),
        HMat.nzvalLocal.Data(),
        Sfile.empty(),
        SMat.nzvalLocal.Data(),
        muMin,
        muMax,
        numShift,
        ordering,
        npPerShift,
        numProcSymbFact,
        MPI_COMM_WORLD,
        shiftList.Data(),
        inertiaListInt.Data(),
        &info);

    if( info != 0 ){
      if( mpirank == 0 ){
        cout << "Inertia count routine gives info = " << info << ". Exit now." << endl;
      }
      MPI_Finalize();
      return info;
    }

    if( mpirank == 0 ){
      cout << "Time for computing the inertia is " << timeEnd  - timeSta << endl;

      for( Int i = 0; i < numShift; i++ )
        printf( "Shift = %25.15f  inertia = %15d\n", 
            shiftList[i], inertiaListInt[i] / 2 );
    }

  }
  catch( std::exception& e )
  {
    std::cerr << "Processor " << mpirank << " caught exception with message: "
      << e.what() << std::endl;
#ifndef _RELEASE_
    DumpCallStack();
#endif
  }

  MPI_Finalize();

  return 0;
}
