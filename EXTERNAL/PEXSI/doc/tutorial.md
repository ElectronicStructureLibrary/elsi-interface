Tutorial              {#pageTutorial}
========

- @subpage pagePEXSIPlan
- @subpage pagePselinvRealSymmetric
- @subpage pagePselinvComplexSymmetric
- @subpage pagePselinvRealUnsymmetric
- @subpage pagePselinvComplexUnsymmetric
- @subpage pageDFT1
- @subpage pageDFT2

<!-- ************************************************************ -->
@page pagePEXSIPlan Using plans and generating log files
\tableofcontents

%PEXSI is written in C++, and the subroutines cannot directly interface
with other programming languages such as C or FORTRAN.  To solve
this problem, the %PEXSI internal data structure is handled using a
datatype @ref PPEXSIPlan.  The idea and the usage of PPEXSIPlan is similar to
`fftw_plan` in the
[FFTW](http://www.fftw.org/~fftw/fftw3_doc/Using-Plans.html#Using-Plans)
package.

In %PEXSI, a matrix is generally referred to as a "pole". The 
factorization and selected inversion procedure for a pole is computed
in parallel using `numProcRow * numProcCol` processors.
 
When only selected inversion (PSelInv) is used, it is recommended to
set the mpisize of the communicator `comm` to be just `numProcRow * numProcCol`.
 
When %PEXSI is used to evaluate a large number of inverse matrices
such as in the electronic structure calculation, mpisize should be 
`numPole*numProcRow*numProcCol`, where `numPole` inverse matrices
can be processed in parallel.

The output information is controlled by the `outputFileIndex` variable.
For instance, if this index is 1, then the corresponding processor will
output to the file `logPEXSI1`.  If outputFileIndex is negative, then
this processor does NOT output logPEXSI files.

**Note** 

- Each processor must output to a **different** file if outputFileIndex
  is non-negative.  
- When many processors are used, it is **not recommended** for all
  processors to output the log files. This is because the IO takes time
  and can be the bottleneck on many architecture. A good practice is to
  let the master processor output information (generating `logPEXSI0`) or 
  to let the master processor of each pole to output the information.


@code{.c}
#include  "c_pexsi_interface.h"
...
{
  PPEXSIPlan   plan;

  plan = PPEXSIPlanInitialize( 
      comm, 
      numProcRow,
      numProcCol,
      outputFileIndex, 
      &info );

  /* ... Computation using plan ... */

  PPEXSIPlanFinalize(
      plan,
      &info );
} 
@endcode




<!-- ************************************************************ -->
@page pagePselinvRealSymmetric Parallel selected inversion for a real symmetric matrix
\tableofcontents

The parallel selected inversion routine for a real symmetric matrix can
be used as follows. This assumes that the size of `MPI_COMM_WORLD` is
`nprow * npcol`.

@code{.c}
#include  "c_pexsi_interface.h"
...
{
  /* Setup the input matrix in distributed compressed sparse column (CSC) format */ 
  ...;
  /* Initialize PEXSI. 
   * PPEXSIPlan is a handle communicating with the C++ internal data structure */
  PPEXSIPlan   plan;
  
  plan = PPEXSIPlanInitialize( 
      MPI_COMM_WORLD, 
      nprow,
      npcol,
      mpirank, 
      &info );

  /* Tuning parameters of PEXSI. The default options is reasonable to
   * start, and the parameters in options can be changed.  */
  PPEXSIOptions  options;
  PPEXSISetDefaultOptions( &options );


  /* Load the matrix into the internal data structure */
  PPEXSILoadRealSymmetricHSMatrix( 
      plan, 
      options,
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      AnzvalLocal,
      1,     // S is an identity matrix here
      NULL,  // S is an identity matrix here
      &info );

  /* Factorize the matrix symbolically */
  PPEXSISymbolicFactorizeRealSymmetricMatrix( 
      plan,
      options,
      &info );

  /* Main routine for computing selected elements and save into AinvnzvalLocal */
  PPEXSISelInvRealSymmetricMatrix (
      plan,
      options,
      AnzvalLocal,
      AinvnzvalLocal,
      &info );

  ...;
  /* Post processing AinvnzvalLocal */
  ...; 

  PPEXSIPlanFinalize(
      plan,
      &info );
} 
@endcode

This routine computes the selected elements of the matrix 
\f$A^{-1}=(H - z S)^{-1}\f$ in parallel.  The input matrix \f$H\f$
follows the @ref secDistCSC, defined through the variables `colptrLocal`,
`rowindLocal`, `HnzvalLocal`.  The input matrix \f$S\f$ can be omitted if it
is an identity matrix and by setting `isSIdentity=1`. If \f$S\f$ is not
an identity matrix, the nonzero sparsity pattern is assumed to be the
same as the nonzero sparsity pattern of \f$H\f$.  Both `HnzvalLocal` and
`SnzvalLocal` are double precision arrays.  

An example is given in driver_pselinv_real.c, which evaluates the
selected elements of the inverse of the matrix saved in
`examples/lap2dr.matrix`.  See also @ref PPEXSISelInvRealSymmetricMatrix
for detailed information of its usage.




<!-- ************************************************************ -->
@page pagePselinvComplexSymmetric Parallel selected inversion for a complex symmetric matrix
\tableofcontents

The parallel selected inversion routine for a complex symmetric matrix
is very similar to the real symmetric case. An example is given in
driver_pselinv_complex.c. See also @ref PPEXSISelInvComplexSymmetricMatrix
for detailed information of its usage.



<!-- ************************************************************ -->
@page pagePselinvRealUnsymmetric Parallel selected inversion for a real unsymmetric matrix
\tableofcontents

The parallel selected inversion routine for a real unsymmetric matrix can
be used as follows. This assumes that the size of `MPI_COMM_WORLD` is
`nprow * npcol`.

@code{.c}
#include  "c_pexsi_interface.h"
...
{
  /* Setup the input matrix in distributed compressed sparse column (CSC) format */ 
  ...;
  /* Initialize PEXSI. 
   * PPEXSIPlan is a handle communicating with the C++ internal data structure */
  PPEXSIPlan   plan;
  
  plan = PPEXSIPlanInitialize( 
      MPI_COMM_WORLD, 
      nprow,
      npcol,
      mpirank, 
      &info );

  /* Tuning parameters of PEXSI. The default options is reasonable to
   * start, and the parameters in options can be changed.  */
  PPEXSIOptions  options;
  PPEXSISetDefaultOptions( &options );
  

  /* Load the matrix into the internal data structure */
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
      1,     // S is an identity matrix here
      NULL,  // S is an identity matrix here
      &info );

  /* Factorize the matrix symbolically */
  PPEXSISymbolicFactorizeRealUnsymmetricMatrix( 
      plan,
      options,
      &info );

  /* Main routine for computing selected elements and save into AinvnzvalLocal */
  PPEXSISelInvRealUnsymmetricMatrix (
      plan,
      options,
      AnzvalLocal,
      AinvnzvalLocal,
      &info );

  ...;
  /* Post processing AinvnzvalLocal */
  ...; 

  PPEXSIPlanFinalize(
      plan,
      &info );
} 
@endcode

This routine computes the selected elements of the matrix 
\f$A^{-1}=(H - z S)^{-1}\f$ in parallel.  The input matrix \f$H\f$
follows the @ref secDistCSC, defined through the variables `colptrLocal`,
`rowindLocal`, `HnzvalLocal`.  The input matrix \f$S\f$ can be omitted if it
is an identity matrix and by setting `isSIdentity=1`. If \f$S\f$ is not
an identity matrix, the nonzero sparsity pattern is assumed to be the
same as the nonzero sparsity pattern of \f$H\f$.  Both `HnzvalLocal` and
`SnzvalLocal` are double precision arrays.  

An example is given in driver_pselinv_real_unsym.c, which evaluates the
selected elements of the inverse of the matrix saved in
`examples/big.unsym.matrix`.  See also @ref PPEXSISelInvRealUnsymmetricMatrix
for detailed information of its usage.



<!-- ************************************************************ -->
@page pagePselinvComplexUnsymmetric Parallel selected inversion for a complex unsymmetric matrix
\tableofcontents

The parallel selected inversion routine for a complex unsymmetric matrix
is very similar to the real unsymmetric case. An example is given in
driver_pselinv_complex_unsym.c. See also @ref PPEXSISelInvComplexUnsymmetricMatrix
for detailed information of its usage.







<!-- ************************************************************ -->
@page pageDFT1 Solving Kohn-Sham density functional theory: I

The simplest way to use %PEXSI to solve Kohn-Sham density functional
theory is to use the @ref PPEXSIDFTDriver routine. This routine uses
built-in heuristics to obtain values of some parameters in %PEXSI and
provides a relatively small set of adjustable parameters for users to
tune.  This routine estimates the chemical potential self-consistently
using a combined approach of inertia counting procedure and Newton's
iteration through %PEXSI. Some heuristic approach is also implemented in
this routine for dynamic adjustment of the chemical potential and some
stopping criterion.

An example routine is given in driver_ksdft.c, which solves a fake DFT
problem by taking a Hamiltonian matrix from `examples/lap2dr.matrix`.

Here is the structure of the code using the simple driver routine.

@code{.c}
#include  "c_pexsi_interface.h"
...
{
  /* Setup the input matrix in distributed compressed sparse column (CSC) format */ 
  ...;
  /* Initialize PEXSI. 
   * PPEXSIPlan is a handle communicating with the C++ internal data structure */
  PPEXSIPlan   plan;

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

  /* Tuning parameters of PEXSI. See PPEXSIOption for explanation of the
   * parameters */
  PPEXSIOptions  options;
  PPEXSISetDefaultOptions( &options );

  options.numPole  = 60;
  options.temperature  = 0.019; // 3000K
  options.muPEXSISafeGuard  = 0.2; 
  options.numElectronPEXSITolerance = 0.001;

  /* Load the matrix into the internal data structure */
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

  /* Call the simple DFT driver using PEXSI */
  PPEXSIDFTDriver(
      plan,
      options,
      numElectronExact,
      &muPEXSI,                   
      &numElectronPEXSI,         
      &muMinInertia,              
      &muMaxInertia,             
      &numTotalInertiaIter,   
      &numTotalPEXSIIter,   
      &info );

  /* Retrieve the density matrix and other quantities from the plan */

  PPEXSIRetrieveRealSymmetricDFTMatrix(
      plan,
      DMnzvalLocal,
      EDMnzvalLocal,
      FDMnzvalLocal,
      &totalEnergyH,
      &totalEnergyS,
      &totalFreeEnergy,
      &info );

  /* Clean up */
  PPEXSIPlanFinalize(
      plan,
      &info );
} 
@endcode


<!-- ************************************************************ -->
@page pageDFT2 Solving Kohn-Sham density functional theory: II
\tableofcontents

In a DFT calculation, the information of the symbolic factorization can
be reused for different \f$(H,S)\f$ matrix pencil if the sparsity pattern does
not change.  An example routine is given in driver_ksdft.c, which solves
a fake DFT problem by taking a Hamiltonian matrix from
`examples/lap2dr.matrix`.

Here is the structure of the code using the simple driver routine.

@code{.c}
#include  "c_pexsi_interface.h"
...
{
  /* Perform DFT calculation as in the previous note */

  /* Update and obtain another set of H and S */

  /* Solve the problem once again without symbolic factorization */
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

  // No need to perform symbolic factorization 
  options.isSymbolicFactorize = 0;
  // Given a good guess of the chemical potential, no need to perform 
  // inertia counting.
  options.isInertiaCount = 0;
  // Optional update mu0, muMin0, muMax0 in PPEXSIOptions

  PPEXSIDFTDriver(
      plan,
      options,
      numElectronExact,
      &muPEXSI,                   
      &numElectronPEXSI,         
      &muMinInertia,              
      &muMaxInertia,             
      &numTotalInertiaIter,   
      &numTotalPEXSIIter,   
      &info );

  /* Postprocessing */
  
} 
@endcode

@note The built-in heuristics in @ref PPEXSIDFTDriver may not be
optimal. It handles only one \f$(H,S)\f$ pair at a time, and does
not accept multiple matrix pairs \f$\{(H_l,S_l)\}\f$ as in the case of
spin-orbit polarized calculations.  For expert users and developers, it
should be relatively easy to dig into the driver routine, and only use
@ref PEXSI::PPEXSIData::SymbolicFactorizeRealSymmetricMatrix "SymbolicFactorizeRealSymmetricMatrix" 
(for symbolic factorization), 
@ref PEXSI::PPEXSIData::CalculateNegativeInertiaReal "CalculateNegativeInertiaReal" 
(for inertia counting), and
@ref PEXSI::PPEXSIData::CalculateFermiOperatorReal "CalculateFermiOperatorReal" 
(for one-shot %PEXSI calculation) to improve heuristics and extend the
functionalities.
