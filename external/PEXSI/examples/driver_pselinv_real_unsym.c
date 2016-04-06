/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

   Authors: Mathias Jacquelin

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
 * @file driver_pselinv_real_unsym.c
 * @brief Example for using the driver interface for parallel selected
 * inversion of a real unsymmetric matrix.
 *
 *
 * @date 2015-01-13 Original version.
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
  double*       AnzvalLocal;
  double*       AinvnzvalLocal;
  int           nprow, npcol;
  int           info;
  char*         Rfile;   

  int           i, j, irow, jcol;
  int           numColLocalFirst, firstCol;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );


  /* Below is the data used for the toy g20 matrix */

  nprow               = 1;
  npcol               = mpisize;
  Rfile               = "big.unsym.matrix";


  /* Read the matrix */
  ReadDistSparseMatrixFormattedHeadInterface(
      Rfile,
      &nrows,
      &nnz,
      &nnzLocal,
      &numColLocal,
      MPI_COMM_WORLD );
      
  if( mpirank == 0 ){
    printf("On processor 0...\n");
    printf("nrows       = %d\n", nrows );
    printf("nnz         = %d\n", nnz );
    printf("nnzLocal    = %d\n", nnzLocal );
    printf("numColLocal = %d\n", numColLocal );
  }

  /* Allocate memory */
  colptrLocal = (int*)malloc( (numColLocal+1) * sizeof(int) );
  rowindLocal = (int*)malloc( nnzLocal * sizeof(int) );
  AnzvalLocal = (double*)malloc( nnzLocal * sizeof(double) );
  AinvnzvalLocal = (double*)malloc( nnzLocal * sizeof(double) );

  /* Read the matrix */
  ReadDistSparseMatrixFormattedInterface(
      Rfile,
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      AnzvalLocal,
      MPI_COMM_WORLD );


  /* Initialize PEXSI */

  PPEXSIOptions  options;
  PPEXSISetDefaultOptions( &options );
  options.npSymbFact = 1;
  options.ordering = 0;
  options.rowOrdering = 1;
  options.symmetric = 0;
  options.transpose = 0;
  options.verbosity = 1;

  PPEXSIPlan   plan;

  plan = PPEXSIPlanInitialize( 
      MPI_COMM_WORLD, 
      nprow,
      npcol,
      mpirank, 
      &info );

  PPEXSILoadRealUnsymmetricHSMatrix( 
      plan, 
      options,
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      AnzvalLocal,
      1,     // S is identity
      NULL,  // S is identity
      &info );

  PPEXSISymbolicFactorizeRealUnsymmetricMatrix( 
      plan,
      options,
      AnzvalLocal,
      &info );

  PPEXSISelInvRealUnsymmetricMatrix (
      plan,
      options,
      AnzvalLocal,
      AinvnzvalLocal,
      &info );



  if( info != 0 ){
    if( mpirank == 0 ){
      printf("PSelInv routine gives info = %d. Exit now.\n", info );
      printf("The error message is in logPEXSI* files.\n" );
    }
    MPI_Finalize();
    return info;
  }


  /* The first processor output the diagonal elements in natural order
   */
  if( mpirank == 0 ){
    numColLocalFirst = nrows / mpisize;
    firstCol         = mpirank * numColLocalFirst;
    for( j = 0; j < numColLocal; j++ ){
      jcol = firstCol + j + 1;
      for( i = colptrLocal[j]-1; 
           i < colptrLocal[j+1]-1; i++ ){
        irow = rowindLocal[i];
        if( irow == jcol ){
          printf("Ainv[%5d,%5d] = %15.10e\n", 
              irow, irow,
              AinvnzvalLocal[i]);
        }
      }
    } // for (j)
  }

  /* Clean up */

  PPEXSIPlanFinalize( 
      plan,
      &info );
  
  if( mpirank == 0 ){ 
    printf("\nAll calculation is finished. Exit the program.\n");
  }



  /* Deallocate memory */
  free( colptrLocal );
  free( rowindLocal );
  free( AnzvalLocal );
  free( AinvnzvalLocal );

  
  MPI_Finalize();

  return 0;
}
