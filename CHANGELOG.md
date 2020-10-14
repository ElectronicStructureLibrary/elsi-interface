# ELSI changelog

## Not released

### PEXSI
* Updated redistributed (PT-)SCOTCH source code to version 6.1.0.

### NTPoly
* Updated redistributed NTPoly source code to version 2.5.1.

### Known issues
* The ELPA code can not be compiled with the NAG Fortran compiler, due to the
  use of GNU extensions in ELPA.
* Depending on the choice of k-points, the complex PEXSI solver may randomly
  fail at the inertia counting stage.

## v2.6.2 (July 2020)

### ELPA
* Fixed a performance regression of the ELPA2 generic kernel.

## v2.6.1 (June 2020)

### PEXSI
* Removed an improper abort from the error handling code of PEXSI.

## v2.6.0 (June 2020)

### ELSI interface
* C compiler and MPI-3 have become mandatory to build ELSI.
* Added an option to choose which sparsity pattern to use when converting input
  dense matrices to the sparse format used by the solver.

### ELPA
* Updated redistributed ELPA source code to version 2020.05.001, which supports
  single-precision calculations, autotuning of runtime parameters, and (NVIDIA)
  GPU acceleration.

### PEXSI
* AAA method has become the default pole expansion method in PEXSI.
* Increased default number of poles from 20 to 30.
* Improved accuracy of pole expansion based on minimax rational approximation.
* Updated redistributed (PT-)SCOTCH source code to version 6.0.9.

### NTPoly
* Updated redistributed NTPoly source code to version 2.5.0.

### SLEPc-SIPs
* Interface compatible with PETSc 3.13 and SLEPc 3.13.

## v2.5.0 (February 2020)

### ELSI interface
* Added utility subroutines to retrieve the internally computed eigenvalues,
  eigenvectors, and occupation numbers when using density matrix solver
  interfaces with an eigensolver.
* Fixed Marzari-Vanderbilt broadening.

### Solvers
* Added support for the Bethe-Salpeter eigensolvers in the BSEPACK library.

### ELPA
* Interface for externally linked ELPA compatible with ELPA 2019.11.

### NTPoly
* Updated redistributed NTPoly source code to version 2.4.0.

### PEXSI
* Updated redistributed SuperLU\_DIST source code to version 6.2.0.
* Added support for computing the electronic entropy via the free energy density
  matrix.

### BSEPACK
* Redistributed source code of BSEPACK 0.1.
* Added parallel BSE eigensolvers PDBSEIG and PZBSEIG.

## v2.4.1 (November 2019)

### ELSI interface
* Fixed energy-weighted density matrix computation when using ELPA and Fermi
  broadening.

### NTPoly
* Updated redistributed NTPoly source code to version 2.3.2.

## v2.4.0 (November 2019)

### ELSI interface
* Fixed density matrix computation when using ELPA and Methfessel-Paxton
  broadening.

### Solvers
* Added support for the tridiagonalization and pentadiagonalization eigensolvers
  implemented in the EigenExa library.
* Added support for the one-stage and two-stage tridiagonalization eigensolvers
  implemented in the MAGMA library.

### EigenExa
* Interface compatible with EigeExa 2.4.
* Added tridiagonalization eigensolver eigen\_s and pentadiagonalization
  eigensolver eigen\_sx.

### SLEPc-SIPs
* Interface compatible with PETSc 3.12 and SLEPc 3.12.

### MAGMA
* Interface compatible with MAGMA 2.5.
* Added one-stage and two-stage eigensolvers.

## v2.3.1 (July 2019)

### SLEPc-SIPs
* Fixed memory leaks in the redistributed SIPs code.
* Interface compatible with PETSc 3.11 and SLEPc 3.11.

## v2.3.0 (June 2019)

### ELSI interface
* Added density matrix extrapolation subroutines for sparse matrices.
* Extended the test suite to increase code coverage.

### ELPA
* Interface for externally linked ELPA compatible with ELPA 2019.05.

### NTPoly
* Updated redistributed NTPoly source code to version 2.3.1.
* Fixed complex matrix conversion from BLACS\_DENSE and GENERIC\_COO to NTPoly.

### PEXSI
* Fixed complex energy-weighted density matrix.
* Added support for linking ELSI against an externally compiled PEXSI.

## v2.2.1 (March 2019)

### ELSI interface
* Fixed LAPACK eigensolver for ill-conditioned overlap matrices.

### PEXSI
* Updated redistributed SuperLU\_DIST source code to version 6.1.1.

## v2.2.0 (February 2019)

### ELSI interface
* Added utility subroutines for geometry optimization and molecular dynamics
  calculations, including reinitialization of ELSI between geometry steps,
  density matrix extrapolation, and Gram-Schmidt orthogonalization of
  eigenvectors.
* Extended the test suite to increase code coverage.

### Matrix formats
* Added arbitrarily distributed coordinate (GENERIC\_COO) format.

### ELPA
* Interface for externally linked ELPA compatible with ELPA 2018.11.
* Fixed single-precision calculations with externally linked ELPA.
* Fixed internal ELPA two-stage real solver with AVX512 kernel.

### NTPoly
* Updated redistributed NTPoly source code to version 2.2.

### OMM
* Fixed libOMM Cholesky flavor with random initial guess.

### PEXSI
* Updated redistributed PEXSI source code to version 1.2.0.
* Updated redistributed SuperLU\_DIST source code to version 6.1.0.

## v2.1.0 (October 2018)

### ELSI interface
* Adopted literature definition of the electronic entropy.
* Added subroutines to query the version number and date stamp of ELSI.

### Solvers
* Added support for the density matrix purification methods implemented in the
  NTPoly library. Implementation of the same methods with dense linear algebra
  has been removed.

### ELPA
* For externally linked ELPA, added options to perform single-precision
  calculations and to automatically tune the internal runtime parameters of the
  solver.
* Interface for externally linked ELPA compatible with ELPA 2018.05.

### NTPoly
* Redistributed source code of NTPoly 2.0.
* Added canonical purification, trace correcting purification, 4th order trace
  resetting purification, and generalized hole-particle canonical purification
  methods.

### PEXSI
* Updated redistributed PEXSI source code to version 1.0.3, which returns the
  complex density matrix and energy-weighted density matrix instead of their
  transpose.

### SLEPc-SIPs
* Updated interface to support PETSc 3.9 and SLEPc 3.9.

## v2.0.2 (June 2018)

### PEXSI
* Updated redistributed PEXSI source code to version 1.0.1, which fixes the
  complex Fermi operator expansion routine.
* Downgraded redistributed (PT-)SCOTCH source code to version 6.0.0, as newer
  versions seem to be incompatible with PEXSI.

## v2.0.1 (June 2018)

### ELSI interface
* Switched to the [semantic versioning scheme](http://semver.org).
* Fixed building ELSI as a shared library with tests enabled.
* Improved stability when calling PBLAS routines pdtran and pztranc.

## v2.0.0 (May 2018)

### ELSI interface
* CMake build system has replaced the Makefile-based system.
* Added support for building ELSI as a shared library.
* Added support for spin channels and k-points.
* Added support for energy-weighted density matrix.
* Added support for electronic entropy calculations.
* Added support for complex sparse matrix formats for eigensolver and density
  matrix solver interfaces.
* Removed optional variables from mutator subroutines.
* Added matrix I/O subroutines using the MPI I/O standard.
* Removed TOMATO dependency for the test suite.
* Added a unified JSON output framework via the FortJSON library.

### Solvers
* Added support for the SLEPc-SIPs solver (PETSc 3.8 and SLEPc 3.8 required).
* Implemented density matrix purification with dense linear algebra operations.

### Matrix formats
* Added 1D block-cyclic compressed sparse column (SIESTA\_CSC) format.

### ELPA
* Updated redistributed ELPA source code to version 2016.11.001.
* Added AVX512 kernel.
* Made the two-stage solver default for all matrix sizes.
* Updated the interface for externally linked ELPA to the AEO version (ELPA
  release 2017.05 or later). GPU acceleration and GPU kernels may be enabled
  through the ELSI interface for externally linked ELPA.

### PEXSI
* Updated redistributed PEXSI source code to version 1.0.0.
* Reduced the default number of poles to 20 without sacrificing accuracy.
* Switched to the PT-SCOTCH library as the default sparse matrix reordering
  software.
* Redistributed SuperLU\_DIST 5.3.0 and (PT-)SCOTCH 6.0.5a libraries. Users may
  still provide their own SuperLU\_DIST library linked against any compatible
  sparse matrix reordering library.
* Removed ParMETIS as a mandatory external dependency for PEXSI.

## v1.0.0 (May 2017)

### Solvers
* ELPA (version 2016.11.001.pre)
* libOMM
* PEXSI (version 0.10.2)

### Matrix formats
* 2D block-cyclic dense (BLACS format)
* 1D block compressed sparse column (PEXSI format)
