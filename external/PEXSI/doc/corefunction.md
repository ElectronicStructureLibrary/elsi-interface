Core Functionality            {#pageCoreFunction}
==================

- @subpage pageBasicCore
- @subpage pageDataType
- @subpage pagePole
- @subpage pageFactor
- @subpage pageSelInv
- @subpage pageCCPP
- @subpage pageFORTRAN


<!-- ************************************************************ -->
@page pageBasicCore Basic
\tableofcontents

You should be able to skip most of this section if you only intend to
use the driver routines, described in the @ref pageTutorial section and
in @ref c_pexsi_interface.h.  

In such case for C/C++ programmers, include the interface file:

@code{.c}
#include  "c_pexsi_interface.h"
@endcode

For FORTRAN programmers, there is no interface routines such as
`f_pexsi_interface.F90` yet.  However, the FORTRAN routines can directly
be used.  See @ref pageFORTRAN for more information.

The remaining section is mainly for C++ developers to have more detailed control
of the %PEXSI package.  

For C++ and usage beyond the driver routines, include the following file
@code{.cpp}
#include  "ppexsi.hpp"
@endcode

For developers, 

- For VI/VIM users, PEXSI seems to be best visualized with the following
  options concerning indentation.

~~~~~~~~~~
set tabstop=2       
set shiftwidth=2
set expandtab 
~~~~~~~~~~
 
<!--
- For EMACS users or users of other text editors, **to be added**.
-->


<!-- ************************************************************ -->
@page pageDataType Data type 
\tableofcontents


Basic data type    {#secBasicDataType}
===============

To enhance potential transplantability of the code, some basic data
types are constants are defined in @ref environment.hpp.

The basic data types `int` and `double` are redefined as `Int` and
`Real`, in order to improve compatibility for different architecture
especially on 64-bit machines (**not implemented yet**). 
The 64-bit long integer `int64_t` is also redefined as `LongInt`.

The complex arithmetic is treated using the standard C++ `<complex>`
library.  The complex data type is `std::complex<double>`, and is
redefined as `Complex` in the implementation.  


@code{.cpp}
typedef    int                   Int;
typedef    int64_t               LongInt;
typedef    double                Real;
typedef    std::complex<double>  Complex; 
@endcode

NumVec, NumMat, NumTns             {#secNumStructure}
======================

The design of %PEXSI tries to eliminate as much as possible the direct
usage of pointers. This helps reducing memory leak.  Commonly used
pointers are wrapped into different classes.

The most commonly used are @ref PEXSI::NumVec "NumVec", @ref
PEXSI::NumMat "NumMat", and @ref PEXSI::NumTns "NumTns", which are
wrappers for 1D array (vector), 2D array (matrix) and 3D array (tensor),
respectively.  The arrays are always saved contiguously in memory as a
1D array. Column-major ordering is assumed for arrays of all dimensions,
which makes the arrays directly compatible with BLAS/LAPACK libraries.


These wrapper classes can both actually own an array (by specifying
`owndata_=true`), and just view an array (by specifying `owndata_=false`).
Elements of arrays can be accessed directly as in FORTRAN convention,
such as `A(i)` (NumVec), `A(i,j)` (NumMat), and `A(i,j,k)` (NumTns).  

The underlying pointer can be accessed using the member function `Data()`.

**Commonly used wrapper classes**

`NumVec`:

@code{.cpp}
typedef NumVec<bool>       BolNumVec;
typedef NumVec<Int>        IntNumVec;
typedef NumVec<Real>       DblNumVec;
typedef NumVec<Complex>    CpxNumVec;
@endcode

`NumMat`:

@code{.cpp}
typedef NumMat<bool>       BolNumMat;
typedef NumMat<Int>        IntNumMat;
typedef NumMat<Real>       DblNumMat;
typedef NumMat<Complex>    CpxNumMat;
@endcode

`NumTns`:

@code{.cpp}
typedef NumTns<bool>       BolNumTns;
typedef NumTns<Int>        IntNumTns;
typedef NumTns<Real>       DblNumTns;
typedef NumTns<Complex>    CpxNumTns;
@endcode


Distributed compressed sparse column (CSC) format    {#secDistCSC}
=================================================

We use the Compressed Sparse Column (CSC) format, a.k.a. the Compressed
Column Storage (CCS) format for storing a sparse matrix.  See

http://netlib.org/linalg/html_templates/node92.html

for the explanation of the format.

We adopt the following convention for distributed CSC format for saving a
sparse matrix on distributed memory parallel machines.  We assume the
number of processor is \f$P\f$, the number of rows and columns of the
matrix is \f$N\f$.  The class for distributed memory CSC format matrix
is @ref PEXSI::DistSparseMatrix "DistSparseMatrix".

- `DistSparseMatrix` uses FORTRAN convention (1-based) indices for
  `colptrLocal` and `rowindLocal`, i.e. the first row and the first column indices
  are 1 instead of 0. 
- `mpirank` follows the standard C convention, i.e. the `mpirank` for
  the first processor is 0.
- Each processor holds \f$\lfloor N/P \rfloor\f$ consequentive columns,
  with the exception that the last processor (`mpirank == P-1`) holds a
  all the remaining \f$N - (P-1) \lfloor N/P \rfloor\f$ columns. 
  The first column holds by the i-th processor is \f$ i \lfloor N/P \rfloor\f$.
  The number of columns on each local processor is usually denoted by
  `numColLocal`. 
- `colptrLocal`, which is an integer array of type `IntNumVec` of
  dimension `numColLocal + 1`, stores the pointers to the nonzero row
  indices and nonzero values in `rowptrLocal` and `nzvalLocal`,
  respectively.  
- `rowindLocal`, which is an integer array of type `IntNumVec` of
  dimension `nnzLocal`, stores the nonzero row indices in each column.
- `nzvalLocal`, which is an array of flexible type (usually `Real` or
  `Complex`) `NumVec` of dimension `nnzLocal`, stores the nonzero values
  in each column.


<!-- ************************************************************ -->
@page pagePole Pole expansion   
\tableofcontents

The pole expansion is used to expand Fermi-Dirac functions and other
derived quantities using a number of Green's functions (poles).

> @ref PEXSI::GetPoleDensity "GetPoleDensity" 

Pole expansion for the Fermi-Dirac operator.
This is the most commonly used subroutine for the pole expansion,
and can be used to compute the shifts and weights for calculating
the density matrix, the total energy, and the Hellman-Feynman
force.  This routine obtains the expansion

\f[
  f_{\beta} (z) = \frac{2}{1+e^{\beta z}} \approx 
  \mathrm{Im} \sum_{l=1}^{P} \frac{\omega^{\rho}_l}{z-z_l}
\f]


> @ref PEXSI::GetPoleDensityDrvMu "GetPoleDensityDrvMu" 

Pole expansion for the derivative of the Fermi-Dirac
operator with respect to the chemical potential mu.
This routine can be used to evaluate the derivative of the number
of electrons with respect to the chemical potential for the
Newton step for updating the chemical potential.

Note that \f$f_{\beta}\f$ does not explicitly contain \f$\mu\f$,
so this routine actually computes the expansion

\f[
   -\frac{\partial f_{\beta}}{\partial z} (z) =
   2\beta \frac{e^{\beta z}}{(1+e^{\beta z})^2} 
   \approx \mathrm{Im} \sum_{l=1}^{P}
   \frac{\omega^{\mu}_l}{z-z_l}
\f]


> @ref PEXSI::GetPoleDensityDrvT "GetPoleDensityDrvT" 

Pole expansion for the derivative of the Fermi-Dirac
operator with respect to the temperature T \f$(1/\beta)\f$.

This routine can be used to extrapolate the number of electrons
from a finite temperature calculation to a zero temperature
calculation, using the derivative information.  However, this
functionality is not used anymore in the current version of
%PEXSI.
                                                                
                                                                
\f[
   \frac{\partial f_{\beta}}{\partial (1/\beta)} (z) =
   2 \beta^2 z \frac{e^{\beta z}}{(1+e^{\beta z})^2} 
   \approx \mathrm{Im} \sum_{l=1}^{P}
   \frac{\omega^{T}_l}{z-z_l}
\f]



> @ref PEXSI::GetPoleHelmholtz "GetPoleHelmholtz" 

Pole expansion for the Helmholtz free energy function.

This routine can be used to compute the (Helmholtz) free energy
when finite temperature effect exists. This is especially
important for metallic system and other small gapped systems. 
This routine expands the free energy function

\f[
   f^{\mathcal{F}}_{\beta}(z) = -\frac{2}{\beta} \log 
   (1 + e^{-\beta z}) \approx \mathrm{Im} \sum_{l=1}^{P}
   \frac{\omega^{\mathcal{F}}_l}{z-z_l}
\f]

> @ref PEXSI::GetPoleForce "GetPoleForce" 

This routine can be used to compute the Pulay contribution of the
atomic force in electronic structure calculations.  This term is
especially important when basis set is not complete and changes
with atomic positions.
This routine expands the function used in the energy density matrix.  

\f[
   f^{E}_{\beta}(z) = (z+\mu) f_{\beta}(z) 
   \approx \mathrm{Im} \sum_{l=1}^{P}
   \frac{\omega^{E}_l}{z-z_l}
\f]

Note that when \f$z=H-\mu I\f$, \f$f^{E}_{\beta}(H-\mu I) = H
f_{\beta}(H-\mu I)\f$, and therefore the energy density matrix can be
directly used to compute the band energy without using eigenvalues.


<!-- ************************************************************ -->
@page pageFactor Factorization
\tableofcontents



Procedure for factorization       {#secProcedureFactor}
===========================

Before the selected inversion step, the matrix saved in 
[DistSparseMatrix](@ref PEXSI::DistSparseMatrix) format must first be
factorized.  In principle this can be done with any \f$LDL^T\f$
factorization or \f$LU\f$ factorization routines.  In the current
version of %PEXSI, [SuperLU_DIST
v3.3](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) is used for the
\f$LU\f$ factorization.  

@note
To avoid conflict with other routines in %PEXSI, the SuperLU_DIST
routines are encapsulated in superlu_dist_interf.cpp. Access to
SuperLU_DIST routines are made through a wrapper class
[SuperLUMatrix](@ref PEXSI::SuperLUMatrix).

The basic steps for factorization include:

- Convert a `DistSparseMatrix` into the native format (`SuperMatrix` in
  SuperLU_DIST) of the factorization routine.
  
- Symbolic factorization.

- Numerical factorization.

Related structures and subroutines
----------------------------------

> @ref PEXSI::SuperLUGrid "SuperLUGrid"

A thin interface for the mpi grid strucutre in SuperLU.

> @ref PEXSI::SuperLUOptions "SuperLUOptions"

A thin interface for passing parameters to set the SuperLU
options.

> @ref PEXSI::SuperLUMatrix::DistSparseMatrixToSuperMatrixNRloc "SuperLUMatrix::DistSparseMatrixToSuperMatrixNRloc"

Convert a distributed sparse matrix in compressed sparse
column format into the SuperLU compressed row format.  The output is
saved in the current %SuperLUMatrix.

@note
Although LU factorization is used, the routine
assumes that the matrix is strictly symmetric, and therefore the
compressed sparse row (CSR) format, used by SuperLU_DIST, gives
exactly the same matrix as formed by the compresed sparse column
format (CSC).

> @ref PEXSI::SuperLUMatrix::SymbolicFactorize "SuperLUMatrix::SymbolicFactorize"

This routine factorizes the superlu matrix symbolically.  Symbolic
factorization contains three steps.

- Permute the matrix to reduce fill-in.
- Symbolic factorize the matrix.
- Distribute the matrix into 2D block cyclic format.

This routine is controlled via 
[SuperLUOptions](@ref PEXSI::SuperLUOptions). In particular, the permutation strategy is
controlled by 
[SuperLUOptions::ColPerm](@ref PEXSI::SuperLUOptions::ColPerm).
 
> @ref PEXSI::SuperLUMatrix::NumericalFactorize "SuperLUMatrix::NumericalFactorize"

Performs LU factorization numerically. 


Example
-------


@code{.cpp}
#include "ppexsi.hpp"
{
  ...;
  // Construct AMat
  DistSparseMatrix<Complex>  AMat;
  ...;

  // Setup SuperLU
  SuperLUGrid<Complex> g( comm, nprow, npcol );
  SuperLUOptions luOpt;
  luOpt.ColPerm = "MMD_AT_PLUS_A";
  SuperLUMatrix<Complex> luMat( g );

  // Matrix conversion
  luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt );

  // Symbolic factorization
  luMat.SymbolicFactorize();

  // Numerical factorization
  luMat.NumericalFactorize();

  ...;
}
@endcode

Reuse symbolic factorization      {#secSymbolicReuse}
============================

In SuperLU_DIST, the same symbolic factorization can be reused for
factorizing different matrices.  To reuse the symbolic factorization,
one should follow the steps below.

(After symbolic factorization)
- Destroy the `SuperMatrix`.
- Convert another `DistSparseMatrix` into the native format (`SuperMatrix` in
  SuperLU_DIST) of the factorization routine.
- Redistribute the matrix into 2D block cyclic format.
- Numerical factorization.

Related structures and subroutines
----------------------------------

> @ref PEXSI::SuperLUMatrix::DestroyAOnly "SuperLUMatrix::DestroyAOnly"

Releases the data in A but keeps other data, such as LUstruct. 

This allows one to perform factorization of
matrices of the same pattern, such as the option

`fact = SamePattern_SameRowPerm`

in SuperLU_DIST.

> @ref PEXSI::SuperLUMatrix::Distribute "SuperLUMatrix::Distribute"

Distribute redistrbutes the SuperMatrix in parallel so that it is ready
for the numerical factorization.

Example
-------

@code{.cpp}
#include "ppexsi.hpp"
{
  ...;
  // Construct AMat
  DistSparseMatrix<Complex>  AMat;
  ...;

  // Setup SuperLU
  SuperLUGrid<Complex> g( comm, nprow, npcol );
  SuperLUOptions luOpt;
  luOpt.ColPerm = "MMD_AT_PLUS_A";
  SuperLUMatrix<Complex> luMat( g );

  // Matrix conversion
  luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt );

  // Symbolic factorization
  luMat.SymbolicFactorize();

  // Destroy the SuperMatrix saved in luMat.
  luMat.DestroyAOnly();


  // Construct another matrix BMat with the same sparsity pattern as A.
  DistSparseMatrix<Complex>  BMat; 
  ...;
  // Matrix conversion
  luMat.DistSparseMatrixToSuperMatrixNRloc( BMat, luOpt );
  // Redistribute into 2D block cyclic format.
  luMat.Distribute();


  // Numerical factorization
  luMat.NumericalFactorize();

  ...;
}
@endcode

<!--
Triangular solve and accuracy check    {#secTriangularSolve}
===================================

The triangular solve routines provided by SuperLU_DIST can be used to
check the accuracy of the factorization as well as the selected
inversion.

(After numericl factorization)
- Construct the distributed right hand sides.
- Solve \f$Ax=b\f$. Multiple right hand sides can be solved
  simultaneously.
-->

<!--
Related structures and subroutines
---------------------------------->

<!--
> @ref PEXSI::SuperLUMatrix::SolveDistMultiVector "SuperLUMatrix::SolveDistMultiVector"

Solve A x = b with b overwritten by x for distributed multivector.

> @ref PEXSI::SuperLUMatrix::CheckErrorDistMultiVector "SuperLUMatrix::CheckErrorDistMultiVector"

Print out the error by direct comparison with the true solution in
distributed format.

Example
------->

<!--
The following example performs factorization, solves for a series of
right hand sides and compare the accuracy.

~~~~~~~~~~{.cpp}
#include "ppexsi.hpp"
{
  ...;
  // Construct AMat
  DistSparseMatrix<Complex>  AMat;
  ...;

  // Setup SuperLU
  SuperLUGrid<Complex> g( comm, nprow, npcol );
  SuperLUOptions luOpt;
  luOpt.ColPerm = "MMD_AT_PLUS_A";
  SuperLUMatrix<Complex> luMat( g );

  // Matrix conversion
  luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt );

  // Symbolic factorization
  luMat.SymbolicFactorize();

  // Numerical factorization
  luMat.NumericalFactorize();
  
  // Construct a global matrix (for error checking)
  SuperLUMatrix<Complex> A1( g ), GA( g );
  A1.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt );
  A1.ConvertNRlocToNC( GA );
  
  // Construct the distributed right hand sides and the exact solution.
  CpxNumMat xTrueGlobal(n, nrhs), bGlobal(n, nrhs);
  CpxNumMat xTrueLocal, bLocal;
  UniformRandom( xTrueGlobal );
  GA.MultiplyGlobalMultiVector( xTrueGlobal, bGlobal );
  A1.DistributeGlobalMultiVector( xTrueGlobal, xTrueLocal );
  A1.DistributeGlobalMultiVector( bGlobal,     bLocal );


  // Solve and check the error.
  luMat.SolveDistMultiVector( bLocal, berr );
  luMat.CheckErrorDistMultiVector( bLocal, xTrueLocal );

  ...;
}
~~~~~~~~~~
-->


<!-- ************************************************************ -->
@page pageSelInv Selected Inversion
\tableofcontents



Procedure for Selected Inversion     {#secProcedureSelInv}
================================


After factorizing a [SuperLUMatrix](@ref PEXSI::SuperLUMatrix) luMat (See the [factorization](@ref secProcedureFactor) page for
information on how to perform factorization), the parallel selected inversion can be computed.



@note
To provide a layer of abstraction from the matrix format used during the factorization, the [PMatrix](@ref PEXSI::PMatrix) class is used during the selected inversion.

@note
All major operations of [PMatrix](@ref PEXSI::PMatrix), including the selected inversion, are defined directly as member functions of [PMatrix](@ref PEXSI::PMatrix).

The basic steps for selected inversion are:
- Conversion from [SuperLUMatrix](@ref PEXSI::SuperLUMatrix) to [PMatrix](@ref PEXSI::PMatrix).
- Preparation of communicators and preprocessing.
- Parallel selected inversion.
- Conversion from [PMatrix](@ref PEXSI::PMatrix) back to [DistSparseMatrix](@ref PEXSI::DistSparseMatrix) format.







Related structures and subroutines
----------------------------------

> @ref PEXSI::GridType "GridType"

A thin interface for the mpi grid strucutre in PSelInv.
GridType should be consistent with the grid used by SuperLU.

@note
It is the user's responsibility to enforce the coherence between [SuperLUGrid](@ref PEXSI::SuperLUGrid) and GridType.

> @ref PEXSI::SuperNodeType "SuperNodeType"

A data structure containing the supernodal partitioning of the matrix.

@note
It is the user's responsibility to initialize this data structure after [SuperLUMatrix::SymbolicFactorize](@ref PEXSI::SuperLUMatrix::SymbolicFactorize) has been called.
This is done using the [SuperLUMatrix::SymbolicToSuperNode](@ref PEXSI::SuperLUMatrix::SymbolicToSuperNode] utility routine.

> @ref PEXSI::PMatrix "PMatrix"

PMatrix contains the main data structure and
computational routines for parallel selected inversion.  


> @ref PEXSI::SuperLUMatrix::LUstructToPMatrix "SuperLUMatrix::LUstructToPMatrix"; 

Converts a compressed row format [SuperLUMatrix](@ref PEXSI::SuperLUMatrix) into a PMatrix object, using the compressed column format used by PSelInv.

@note
Although LU factorization is used, the routine
assumes that the matrix is strictly symmetric, and therefore the
compressed sparse row (CSR) format, used by SuperLU_DIST, gives
exactly the same matrix as formed by the compresed sparse column
format (CSC).

> @ref PEXSI::PMatrix::Create "PMatrix::Create"

This static factory routine instantiates the correct PMatrix object type depending on matrix structure.
The matrix structure is specified by the [SuperLUOptions::symmetric](@ref PEXSI::SuperLUOptions::symmetric) attribute of the [SuperLUOptions](@ref PEXSI::SuperLUOptions) data structure.
 

> @ref PEXSI::PMatrix::ConstructCommunicationPattern "PMatrix::ConstructCommunicationPattern" 

This routine creates the MPI_Communicators and communication pattern used later by both PreSelInv and SelInv routines.
The supernodal elimination tree is exploited to add an additional level of parallelism between supernodes.
[ConstructCommunicationPattern_P2p](@ref PEXSI::PMatrix::ConstructCommunicationPattern_P2p) is called by default.

> @ref PEXSI::PMatrix::PreSelInv "PMatrix::PreSelInv" 

  PreSelInv prepares the structure in L and U so that
  SelInv only involves matrix-matrix multiplication.
  
  @note
  PreSelInv assumes that PEXSI::PMatrix::ConstructCommunicationPattern has been executed.

> @ref PEXSI::PMatrix::SelInv "PMatrix::SelInv" 

   SelInv preforms the actual parallel selected inversion.
  [SelInv_P2p](@ref PEXSI::PMatrix::SelInv_P2p) is called by default.

  @note
  SelInv assumes that [PreSelInv](@ref PEXSI::PMatrix::PreSelInv) has been executed.

> @ref PEXSI::PMatrix::PMatrixToDistSparseMatrix "PMatrix::PMatrixToDistSparseMatrix" 

  Converts the PMatrix back to the original [DistSparseMatrix](@ref PEXSI::DistSparseMatrix) format.

Example
-------

@code{.cpp}
#include "ppexsi.hpp"
{
  ...;
  // Construct AMat
  DistSparseMatrix<Complex>  AMat;
  ...;

  /****** NUMERICAL FACTORIZATION ******/
  // Setup SuperLU
  SuperLUGrid<Complex> g( comm, nprow, npcol );
  SuperLUOptions luOpt;
  luOpt.ColPerm = "MMD_AT_PLUS_A";
  luOpt.symmetric = 0;

  SuperLUMatrix<Complex> luMat( g );

  // Matrix conversion
  luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt );

  // Symbolic factorization
  luMat.SymbolicFactorize();


  // Numerical factorization
  luMat.NumericalFactorize();

  /****** SELECTED INVERSION ******/
  GridType gPM( comm, nprow, npcol );

  SuperNodeType super;
  luMat.SymbolicToSuperNode( super );

  PMatrix<Complex> * PMloc = PMatrix<Complex>::Create(&gPM, &super, &luOpt);

  // Conversion to PMatrix
  luMat.LUstructToPMatrix( *PMloc );
   
  //Create the communication pattern
  PMloc->ConstructCommunicationPattern();

  //Prepare for parallel selected inversion 
  PMloc->PreSelInv();

  //Perform the parallel selected inversion
  PMloc->SelInv();

  //Get the result back in DistSparseMatrix format
  DistSparseMatrix<Scalar> Ainv;
  PMloc->PMatrixToDistSparseMatrix( Ainv );

  ...;

  delete PMloc;
}
@endcode






<!-- ************************************************************ -->
@page pageCCPP C/C++ interface
\tableofcontents

The main interface routines are given in @ref c_pexsi_interface.h.  The
routines are callable from C/C++.  

@note C++ users also have the option of directly using the subroutines
provided in @ref ppexsi.cpp.  The usage can be obtained from 
@ref interface.cpp.

<!-- ************************************************************ -->
@page pageFORTRAN FORTRAN interface
\tableofcontents

The FORTRAN interface is based on the ISO_C_BINDING feature, which is
available for FORTRAN 2003 or later.  The usage of FORTRAN interface is
very similar to the C interface as given in the @ref pageTutorial
section. 

In FORTRAN, the PPEXSIPlan data type is `c_intptr_t` (or equivalently
`INTEGER*8`). The naming of the subroutines is  similar to the C
interface as in @ref c_pexsi_interface.h.  All FORTRAN interface
routines are in @ref f_interface.f90.  
For instance, the subroutine @ref PPEXSIPlanInitialize (C/C++) corresponds to
the subroutine @ref f_ppexsi_plan_initialize (FORTRAN).

Example: Parallel selected inversion for a real symmetric matrix

@code{.f90}
integer(c_intptr_t)    :: plan
type(f_ppexsi_options) :: options

! Initialize PEXSI. 
! PPEXSIPlan is a handle communicating with the C++ internal data structure 

! Set the outputFileIndex to be the pole index.
! The first processor for each pole outputs information

if( mod( mpirank, nprow * npcol ) .eq. 0 ) then
  outputFileIndex = mpirank / (nprow * npcol);
else
  outputFileIndex = -1;
endif

plan = f_ppexsi_plan_initialize(&
  MPI_COMM_WORLD,&
  nprow,&
  npcol,&
  outputFileIndex,&
  info )

! Tuning parameters of PEXSI. The default options is reasonable to
! start, and the parameters in options can be changed. 
call f_ppexsi_set_default_options(&
  options )

! Load the matrix into the internal data structure 
call f_ppexsi_load_real_symmetric_hs_matrix(&
      plan,&       
      options,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      colptrLocal,& 
      rowindLocal,&
      HnzvalLocal,&
      1,&
      SnzvalLocal,&
      info ) 

! Factorize the matrix symbolically
call f_ppexsi_symbolic_factorize_real_symmetric_matrix(&
  plan,&
  options,&
  info)

! Main routine for computing selected elements and save into AinvnzvalLocal
call f_ppexsi_selinv_real_symmetric_matrix(& plan,&
  options,&
  AnzvalLocal,&
  AinvnzvalLocal,&
  info)

! Post processing step...

! Release the data saved in the plan
call f_ppexsi_plan_finalize( plan, info )


@endcode


The examples of the FORTRAN interface can be found under `fortran/`
directory, including @ref f_driver_pselinv_real.f90, 
@ref f_driver_pselinv_complex.f90, 
@ref f_driver_ksdft.f90.
