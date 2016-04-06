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
/**
 * @file driver_ksdft_2.c
 * @brief Example for using the new driver interface for performing KSDFT
 * calculations.
 *
 * This file is eventually going to be merged with the driver_ksdft.c
 *
 * @date 2015-11-25  Test for DFTDriver2 with updating strategy of pole
 * expansion
 *
 */
#include  <stdio.h>
#include  <stdlib.h>
#include  "c_pexsi_interface.h"

int main(int argc, char **argv) 
{
  int mpirank, mpisize;
  int           nrows;
  int           nnz;                          
  int           nnzLocal;                     
  int           numColLocal;                  
  int*          colptrLocal;                  
  int*          rowindLocal;                  
  double*       HnzvalLocal;                  
  int           isSIdentity;                  
  double*       SnzvalLocal;                  

  double*       DMnzvalLocal;
  double*       EDMnzvalLocal;
  double*       FDMnzvalLocal;
  double*       localDOSnzvalLocal;  

  double        Energy;
  double        eta;

  double        numElectronExact;

  double        muPEXSI;
  double        numElectronPEXSI;
  double        muMinInertia;
  double        muMaxInertia;
  int           numTotalInertiaIter;
  int           numTotalPEXSIIter;

  double        totalEnergyH;
  double        totalEnergyS;
  double        totalFreeEnergy;
  
  char*         Hfile;
  char*         Sfile;
  int           isFormatted;


  PPEXSIPlan    plan;
  PPEXSIOptions options;

  int           i, j;
  int           nprow, npcol;
  MPI_Comm      readComm;
  int           isProcRead;
  int           info;
  int           outputFileIndex;



  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

  /* Below is the data used for the toy matrix */

#if 1
  numElectronExact    = 12.0;
  nprow               = 1;
  npcol               = 1;
  Hfile               = "lap2dr.matrix";
  Sfile               = "";
  isFormatted         = 1;
  isSIdentity         = 1;
#else
#if 1
  numElectronExact    = 7000.0;
  nprow               = 8;
  npcol               = 8;
  Hfile               = "/project/projectdirs/m1027/PEXSI/LU_C_BN_C_1by1/H_LU.csc";
  Sfile               = "";
  isFormatted         = 0;

  isSIdentity         = 1;
  Energy              = 1.0;
  eta                 = 0.001;
#else
  numElectronExact    = 70000.0;
  nprow               = 8;
  npcol               = 8;
  Hfile               = "/project/projectdirs/m1027/PEXSI/DG_Phosphorene_14000/H.csc";
  Sfile               = "";
  isFormatted         = 0;

  isSIdentity         = 1;
  Energy              = 1.0;
  eta                 = 0.001;
#endif
#endif

  /* Split the processors to read matrix */
  if( mpirank < nprow * npcol )
    isProcRead = 1;
  else
    isProcRead = 0;

  MPI_Comm_split( MPI_COMM_WORLD, isProcRead, mpirank, &readComm );

  if( isProcRead == 1 ){
    printf("Proc %5d is reading file...\n", mpirank );
    /* Read the matrix head for allocating memory */
    if( isFormatted == 1 ){
      ReadDistSparseMatrixFormattedHeadInterface(
          Hfile,
          &nrows,
          &nnz,
          &nnzLocal,
          &numColLocal,
          readComm );
    }
    else{
      ReadDistSparseMatrixHeadInterface(
          Hfile,
          &nrows,
          &nnz,
          &nnzLocal,
          &numColLocal,
          readComm );
    }
    
    if( mpirank == 0 ){
      printf("On processor 0...\n");
      printf("nrows       = %d\n", nrows );
      printf("nnz         = %d\n", nnz );
      printf("nnzLocal    = %d\n", nnzLocal );
      printf("numColLocal = %d\n", numColLocal );
    }


    /* Allocate memory visible to processors in the group of readComm */
    colptrLocal             = (int*) malloc( sizeof(int) * (numColLocal+1) );
    rowindLocal             = (int*) malloc( sizeof(int) * nnzLocal );
    HnzvalLocal             = (double*) malloc( sizeof(double) * nnzLocal );
    SnzvalLocal             = (double*) malloc( sizeof(double) * nnzLocal );
    DMnzvalLocal            = (double*) malloc( sizeof(double) * nnzLocal );
    EDMnzvalLocal           = (double*) malloc( sizeof(double) * nnzLocal );
    FDMnzvalLocal           = (double*) malloc( sizeof(double) * nnzLocal );
    localDOSnzvalLocal      = (double*) malloc( sizeof(double) * nnzLocal );

    /* Actually read the matrix */
    if( isFormatted == 1 ){
      ReadDistSparseMatrixFormattedInterface(
          Hfile,
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          HnzvalLocal,
          readComm );
    }
    else{
      ParaReadDistSparseMatrixInterface(
          Hfile,
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          HnzvalLocal,
          readComm );
    }

    if( isSIdentity == 0 ){
      if( isFormatted == 1 ){
        ReadDistSparseMatrixFormattedInterface(
            Hfile,
            nrows,
            nnz,
            nnzLocal,
            numColLocal,
            colptrLocal,
            rowindLocal,
            SnzvalLocal,
            readComm );
      }
      else{
        ParaReadDistSparseMatrixInterface(
            Hfile,
            nrows,
            nnz,
            nnzLocal,
            numColLocal,
            colptrLocal,
            rowindLocal,
            SnzvalLocal,
            readComm );
      }
    }
    
    if( mpirank == 0 ){ 
      printf("Finish reading the matrix.\n");
    }
  }

  /* Call the PEXSI interface */

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
  options.muInertiaTolerance = 0.0120;
  options.muPEXSISafeGuard  = 0.2; 
  options.numElectronPEXSITolerance = 0.001;
  options.isSymbolicFactorize = 1;

  /* Set the outputFileIndex to be the pole index */
  /* The first processor for each pole outputs information */

  if( mpirank % (nprow * npcol) == 0 ){
    outputFileIndex = mpirank / (nprow * npcol);
  }
  else{
    outputFileIndex = -1;
  }

  plan = PPEXSIPlanInitialize( 
      MPI_COMM_WORLD, 
      nprow,
      npcol,
      outputFileIndex, 
      &info );

  PPEXSILoadRealSymmetricHSMatrix( 
      plan, 
      options,
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      isSIdentity,
      SnzvalLocal,
      &info );

  /* Step 2. PEXSI Solve */

  PPEXSIDFTDriver2(
      plan,
      options,
      numElectronExact,
      &muPEXSI,                   
      &numElectronPEXSI,         
      &muMinInertia,              
      &muMaxInertia,             
      &numTotalInertiaIter,   
      &info );

  if( info != 0 ){
    if( mpirank == 0 ){
      printf("PEXSI solve routine gives info = %d. Exit now.\n", info );
    }
    MPI_Finalize();
    return info;
  }


  if( isProcRead == 1 ){
    PPEXSIRetrieveRealSymmetricDFTMatrix(
        plan,
        DMnzvalLocal,
        EDMnzvalLocal,
        FDMnzvalLocal,
        &totalEnergyH,
        &totalEnergyS,
        &totalFreeEnergy,
        &info );

    if( mpirank == 0 ){
      printf("Output from the main program\n");
      printf("Total energy (H*DM)         = %15.5f\n", totalEnergyH);
      printf("Total energy (S*EDM)        = %15.5f\n", totalEnergyS);
      printf("Total free energy           = %15.5f\n", totalFreeEnergy);
    }
  }

  /* Step 3. Solve the problem once again without symbolic factorization */
  {
    if( mpirank == 0 ){
      printf("To test the correctness of the program, solve the problem \n");
      printf("again without symbolic factorization or inertia counting.\n");
    }

    PPEXSILoadRealSymmetricHSMatrix( 
        plan, 
        options,
        nrows,
        nnz,
        nnzLocal,
        numColLocal,
        colptrLocal,
        rowindLocal,
        HnzvalLocal,
        isSIdentity,
        SnzvalLocal,
        &info );

    // No symbolic factorization
    options.muMin0 = muMinInertia;
    options.muMax0 = muMaxInertia;
    options.isInertiaCount = 0;
    options.isSymbolicFactorize = 0;
    // Reuse previous mu to start
    options.mu0 = muPEXSI;

    PPEXSIDFTDriver2(
        plan,
        options,
        numElectronExact,
        &muPEXSI,                   
        &numElectronPEXSI,         
        &muMinInertia,              
        &muMaxInertia,             
        &numTotalInertiaIter,   
        &info );


    if( info != 0 ){
      if( mpirank == 0 ){
        printf("PEXSI solve routine gives info = %d. Exit now.\n", info );
      }
      MPI_Finalize();
      return info;
    }

    if( isProcRead == 1 ){
      PPEXSIRetrieveRealSymmetricDFTMatrix(
          plan,
          DMnzvalLocal,
          EDMnzvalLocal,
          FDMnzvalLocal,
          &totalEnergyH,
          &totalEnergyS,
          &totalFreeEnergy,
          &info );

      if( mpirank == 0 ){
        printf("Output from the main program\n");
        printf("Total energy (H*DM)         = %15.5f\n", totalEnergyH);
        printf("Total energy (S*EDM)        = %15.5f\n", totalEnergyS);
        printf("Total free energy           = %15.5f\n", totalFreeEnergy);
      }
    }
  }

  /* Step 4. Clean up */

  PPEXSIPlanFinalize( 
      plan,
      &info );
  
  if( mpirank == 0 ){ 
    printf("\nAll calculation is finished. Exit the program.\n");
  }


  if( isProcRead == 1 ){
    free( colptrLocal );
    free( rowindLocal );
    free( HnzvalLocal );
    free( SnzvalLocal );
    free( DMnzvalLocal );
    free( EDMnzvalLocal );
    free( FDMnzvalLocal );
    free( localDOSnzvalLocal );
  }

  
  MPI_Comm_free( &readComm );
  MPI_Finalize();

  return 0;
}
