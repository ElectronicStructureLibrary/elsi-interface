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
/// @file run_ksdft.cpp
/// @brief Similar to driver_ksdft_2 but test the CPP interface
/// directly.
///
/// This can be used to test the correctness of new CPP implementation
/// without going through the interfaces.
/// @date 2015-12-07. Initial revision. Derived from run_ppexsi.cpp
#include "ppexsi.hpp"
#include "c_pexsi_interface.h"

using namespace PEXSI;
using namespace std;


void Usage(){
  std::cout 
		<< "run_ksdft" << std::endl;
}

int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

  // PEXSI parameters
  int           isSIdentity;                  


  double        numElectronExact;

  double        muPEXSI;
  double        numElectronPEXSI;
  double        muMinInertia;
  double        muMaxInertia;
  int           numTotalInertiaIter;
  int           numTotalPEXSIIter;
  int           numColLocal;

  char*         Hfile;
  char*         Sfile;
  int           isFormatted;

  PPEXSIPlan    plan;
  PPEXSIOptions options;

  int           nprow, npcol;
  MPI_Comm      readComm;
  int           isProcRead;
  int           outputFileIndex;


//	if( argc < 25 || argc%2 == 0 ) {
//		if( mpirank == 0 ) Usage();
//		MPI_Finalize();
//		return 0;
//	}
			
	
	try{
		// *********************************************************************
		// Input parameter
    // 
    // FIXME: Currently hard coded
		// *********************************************************************

    numElectronExact    = 12.0;
    nprow               = 1;
    npcol               = 1;
    Hfile               = "lap2dr.matrix";
    Sfile               = "";
    isFormatted         = 1;
    isSIdentity         = 1;

    /* Split the processors to read matrix */
    if( mpirank < nprow * npcol )
      isProcRead = 1;
    else
      isProcRead = 0;

    MPI_Comm_split( MPI_COMM_WORLD, isProcRead, mpirank, &readComm );
		
    DistSparseMatrix<Real> HMat;
    DistSparseMatrix<Real> SMat;

    if( isProcRead == 1 ){
      printf("Proc %5d is reading file...\n", mpirank );
      /* Read the matrix head for allocating memory */
      if( isFormatted == 1 ){
        ReadDistSparseMatrixFormatted( Hfile, 
            HMat, readComm ); 
      }
      else{
        ReadDistSparseMatrix( Hfile, HMat, readComm ); 
      }

      numColLocal = HMat.colptrLocal.m() - 1;

      if( mpirank == 0 ){
        printf("On processor 0...\n");
        printf("nrows       = %d\n", HMat.size );
        printf("nnz         = %d\n", HMat.nnz );
        printf("nnzLocal    = %d\n", HMat.nnzLocal );
        printf("numColLocal = %d\n", numColLocal );
      }

      if( isSIdentity == 0 ){
        if( isFormatted == 1 ){
          ReadDistSparseMatrixFormatted( Sfile, SMat, readComm ); 
        }
        else{
          ReadDistSparseMatrix( Sfile, SMat, readComm ); 
        }
      }

      if( mpirank == 0 ){ 
        printf("Finish reading the matrix.\n");
      }
    } // Read the matrix


		// *********************************************************************
		// Check the input parameters
		// *********************************************************************

		// Initialize

    /* Set the outputFileIndex to be the pole index */
    /* The first processor for each pole outputs information */

    if( mpirank % (nprow * npcol) == 0 ){
      outputFileIndex = mpirank / (nprow * npcol);
    }
    else{
      outputFileIndex = -1;
    }

    Int npPerPole = nprow * npcol;
		PPEXSIData pexsi( MPI_COMM_WORLD, nprow, npcol, outputFileIndex );


    /* Step 1. Initialize PEXSI */

    PPEXSISetDefaultOptions( &options );
    options.muMin0 = 0.0;
    options.muMax0 = 0.5;
    options.mu0    = 0.0;
    options.npSymbFact = 1;
    options.ordering = 0;
    options.isInertiaCount = 1;
    options.maxPEXSIIter   = 1;
    options.verbosity = 1;
    options.deltaE   = 20.0;
    options.numPole  = 40;
    options.temperature  = 0.0019; // 300K
    //  options.muInertiaTolerance = 0.0019;
    options.muPEXSISafeGuard  = 0.2; 
    options.numElectronPEXSITolerance = 0.001;
    options.isSymbolicFactorize = 1;


    pexsi.LoadRealSymmetricMatrix(
          HMat.size,                        
          HMat.nnz,                          
          HMat.nnzLocal,                     
          numColLocal,                  
          HMat.colptrLocal.Data(),
          HMat.rowindLocal.Data(),
          HMat.nzvalLocal.Data(),
          isSIdentity,
          SMat.nzvalLocal.Data(),
          options.verbosity );

    // PEXSI Solve
    pexsi.DFTDriver2(
        numElectronExact,
        options.temperature,
        options.gap,
        options.deltaE,
        options.numPole,
        options.isInertiaCount,
        options.muMin0,
        options.muMax0,
        options.mu0,
        options.muInertiaTolerance,
        options.muInertiaExpansion,
        options.numElectronPEXSITolerance,
        options.matrixType,
        options.isSymbolicFactorize,
        options.ordering,
        options.npSymbFact,
        options.verbosity,
        muPEXSI,
        numElectronPEXSI,
        muMinInertia,
        muMaxInertia,
        numTotalInertiaIter );

    // Retrieve data
    if( isProcRead == 1 ){

      if( mpirank == 0 ){
        printf("Output from the main program\n");
        printf("Total energy (H*DM)         = %15.5f\n", pexsi.TotalEnergyH());
        printf("Total energy (S*EDM)        = %15.5f\n", pexsi.TotalEnergyS());
        printf("Total free energy           = %15.5f\n", pexsi.TotalFreeEnergy());
      }
    }

    // No need for clean up
	}
	catch( std::exception& e )
	{
		statusOFS << std::endl << " ERROR!!! Proc " << mpirank << " caught exception with message: "
			<< e.what() << std::endl;
		statusOFS.close();
		statusOFS << std::endl << " ERROR!!! Proc " << mpirank << " caught exception with message: "
			<< e.what() << std::endl;
#ifndef _RELEASE_
		DumpCallStack();
#endif
	}
	
	MPI_Finalize();

	return 0;
}
