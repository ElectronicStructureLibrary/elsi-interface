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
 * @file driver_fermi_complex.c
 * @brief Example for using the driver interface for performing inertia
 * counting and evaluating the Fermi operator with complex Hermitian
 * matrices.
 *
 * Note: This routine is different from driver_ksdft.c in the sense that
 * it cannot perform searching for mu. However, the routines contained
 * in this file should provide sufficient information for implementing
 * such procedure, such as to be developed in the ELSI interface.
 *
 * @date 2016-09-10  Original version.
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
  double*       RnzvalLocal;                  
  double*       InzvalLocal;                  
  double*       HnzvalLocal;                  
  int           isSIdentity;                  
  double*       SnzvalLocal;                  

  double*       DMnzvalLocal;
  double*       EDMnzvalLocal;
  double*       FDMnzvalLocal;

  int           numShift;
  double*       shiftVec;
  double*       inertiaVec;

  double        Energy;
  double        eta;

  double        numElectronExact;

  double        muPEXSI;
  double        numElectron;
  double        numElectronDrvMu;
  double        muMinInertia;
  double        muMaxInertia;
  int           numTotalInertiaIter;
  int           numTotalPEXSIIter;

  double        totalEnergyH;
  double        totalEnergyS;
  double        totalFreeEnergy;
  
  char*         HfileR;
  char*         HfileI;
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
  HfileR              = "lap2dc_real.matrix";
  HfileI              = "lap2dc_imag.matrix";
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
          HfileR,
          &nrows,
          &nnz,
          &nnzLocal,
          &numColLocal,
          readComm );
    }
//    else{
//      ReadDistSparseMatrixHeadInterface(
//          HfileRj,
//          &nrows,
//          &nnz,
//          &nnzLocal,
//          &numColLocal,
//          readComm );
//    }
    
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
    RnzvalLocal             = (double*) malloc( sizeof(double) * nnzLocal );
    InzvalLocal             = (double*) malloc( sizeof(double) * nnzLocal );
    HnzvalLocal             = (double*) malloc( sizeof(double) * nnzLocal * 2 );
    SnzvalLocal             = (double*) malloc( sizeof(double) * nnzLocal * 2 );
    DMnzvalLocal            = (double*) malloc( sizeof(double) * nnzLocal * 2 );
    EDMnzvalLocal           = (double*) malloc( sizeof(double) * nnzLocal * 2 );
    FDMnzvalLocal           = (double*) malloc( sizeof(double) * nnzLocal * 2 );

    numShift = options.numPole;
    shiftVec                = (double*) malloc( sizeof(double) * numShift );
    inertiaVec              = (double*) malloc( sizeof(double) * numShift );

    /* Actually read the matrix */
    if( isFormatted == 1 ){
      ReadDistSparseMatrixFormattedInterface(
          HfileR,
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          RnzvalLocal,
          readComm );
      
      ReadDistSparseMatrixFormattedInterface(
          HfileI,
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          InzvalLocal,
          readComm );

      for( i = 0; i < nnzLocal; i++ ){
        HnzvalLocal[2*i]   = RnzvalLocal[i];
        HnzvalLocal[2*i+1] = InzvalLocal[i];
      }

    }
//    else{
//      ParaReadDistSparseMatrixInterface(
//          Hfile,
//          nrows,
//          nnz,
//          nnzLocal,
//          numColLocal,
//          colptrLocal,
//          rowindLocal,
//          HnzvalLocal,
//          readComm );
//    }

//    if( isSIdentity == 0 ){
//      if( isFormatted == 1 ){
//        ReadDistSparseMatrixFormattedInterface(
//            Hfile,
//            nrows,
//            nnz,
//            nnzLocal,
//            numColLocal,
//            colptrLocal,
//            rowindLocal,
//            SnzvalLocal,
//            readComm );
//      }
//      else{
//        ParaReadDistSparseMatrixInterface(
//            Hfile,
//            nrows,
//            nnz,
//            nnzLocal,
//            numColLocal,
//            colptrLocal,
//            rowindLocal,
//            SnzvalLocal,
//            readComm );
//      }
//    }
    
    if( mpirank == 0 ){ 
      printf("Finish reading the matrix.\n");
    }
  }

  /* Initialize PEXSI */

  PPEXSISetDefaultOptions( &options );
  options.muMin0 = 0.0;
  options.muMax0 = 0.5;
  options.mu0    = 0.270;
  options.npSymbFact = 1;
  options.ordering = 0;
  options.isInertiaCount = 1;
  options.maxPEXSIIter   = 1;
  options.verbosity = 1;
  options.deltaE   = 20.0;
  options.numPole  = 40;
  options.temperature  = 0.0019; // 300K
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

  PPEXSILoadComplexHSMatrix( 
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

  /* Symbolic factorization */

  // No permutation
  options.rowOrdering = 0;
  PPEXSISymbolicFactorizeComplexSymmetricMatrix( 
      plan,
      options,
      &info );

  // Can change row ordering optionally now.
  // options.rowOrdering = 1;
  PPEXSISymbolicFactorizeComplexUnsymmetricMatrix( 
      plan,
      options,
      HnzvalLocal,
      &info );

  /* Inertia counting */
  numShift = options.numPole;
  for( i = 0; i < numShift; i++ ){
    shiftVec[i] = options.muMin0 + 
      ( options.muMax0 - options.muMin0 ) / numShift * (double)i;
  }


  PPEXSIInertiaCountComplexMatrix( 
      plan,
      options,
      numShift,
      shiftVec,
      inertiaVec,
      &info );

  if( mpirank == 0 ){
    for( i = 0; i < numShift; i++ )
      printf( "Shift = %25.15f  inertia = %25.1f\n", 
          shiftVec[i], inertiaVec[i] );
  }

  /* Evaluate the Fermi operator */
  PPEXSICalculateFermiOperatorComplex(
      plan,
      options,
      options.mu0,
      numElectronExact,
      &numElectron,
      &numElectronDrvMu,
      &info );

  if( mpirank == 0 ){
    printf("numElectron       = %25.15f\n", numElectron);
    printf("numElectronDrvMu  = %25.15f\n", numElectronDrvMu);
  }

  if( info != 0 ){
    if( mpirank == 0 ){
      printf("PEXSI solve routine gives info = %d. Exit now.\n", info );
    }
    MPI_Finalize();
    return info;
  }

  /* Retrieve matrix and energy */

  if( isProcRead == 1 ){
    PPEXSIRetrieveComplexDFTMatrix(
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

  /* Clean up */

  PPEXSIPlanFinalize( 
      plan,
      &info );
  
  if( mpirank == 0 ){ 
    printf("\nAll calculation is finished. Exit the program.\n");
  }


  if( isProcRead == 1 ){
    free( colptrLocal );
    free( rowindLocal );
    free( RnzvalLocal );
    free( InzvalLocal );
    free( HnzvalLocal );
    free( SnzvalLocal );
    free( DMnzvalLocal );
    free( EDMnzvalLocal );
    free( FDMnzvalLocal );
    free( shiftVec );
    free( inertiaVec );
  }

  
  MPI_Comm_free( &readComm );
  MPI_Finalize();

  return 0;
}
