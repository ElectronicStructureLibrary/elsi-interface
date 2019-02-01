# ELSI changelog

## Not released

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
* Fixed single precision calculations with an externally linked ELPA.
* Fixed internal ELPA 2-stage real solver with AVX512 kernel.

### NTPoly
* Updated redistributed NTPoly source code to version 2.1.

### OMM
* Fixed libOMM Cholesky flavor with random initial guess.

### PEXSI
* Updated redistributed PEXSI source code to version 1.2.0.
* Updated redistributed SuperLU\_DIST source code to version 6.1.0.

## ELSI v2.1.0 (October 2018)

### ELSI interface
* Adopted literature definition of the electronic entropy.
* Added subroutines to query the version number and date stamp of ELSI.

### Solvers
* Added support for the density matrix purification methods implemented in the
  NTPoly library. Implementation of the same methods with dense linear algebra
  has been removed.

### ELPA
* For externally linked ELPA, added options to perform single precision
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
  complex-valued density matrix and energy density matrix instead of their
  transpose.

### SLEPc-SIPs
* Updated interface to support PETSc 3.9.4 and SLEPc 3.9.2.

## ELSI v2.0.2 (June 2018)

### PEXSI
* Updated redistributed PEXSI source code to version 1.0.1, which fixes the
  complex-valued Fermi operator expansion routine.
* Downgraded redistributed (PT-)SCOTCH source code to version 6.0.0, as newer
  versions seem to be incompatible with PEXSI.

## ELSI v2.0.1 (June 2018)

### ELSI interface
* Switched to the [semantic versioning scheme](http://semver.org).
* Fixed building ELSI as a shared library with tests enabled.
* Improved stability when calling PBLAS routines pdtran and pztranc.

### Known issues
* Depending on the choice of k-points, the complex PEXSI solver may fail at the
  inertia counting stage.

## ELSI v2.0.0 (May 2018)

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
* Added support for the SLEPc-SIPs solver (PETSc 3.8.4 and SLEPc 3.8.3
  required).
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

## ELSI v1.0.0 (May 2017)

### Solvers
* ELPA (version 2016.11.001.pre)
* libOMM
* PEXSI (version 0.10.2)

### Matrix formats
* 2D block-cyclic dense (BLACS format)
* 1D block compressed sparse column (PEXSI format)
