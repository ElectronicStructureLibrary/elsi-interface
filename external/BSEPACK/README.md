BSEPACK
=======

Contributers
------------
### Meiyue Shao and Chao Yang

Scalable Solvers Group
Computational Research Division
Lawrence Berkeley National Laboratory

Last update: version 0.1, August 2016


1. Introduction
---------------
BSEPACK is a parallel ScaLAPACK-style software library for computing all
eigenpairs of definite Bethe--Salpeter Hamiltonian matrices of the form

    H = [        A,        B
          -conj(B), -conj(A) ],

where

    [       A,       B
      conj(B), conj(A) ]

is Hermitian and positive definite.

The source code can be downloaded from the software homepage
<https://sites.google.com/a/lbl.gov/bsepack/>


2. How to Install
-----------------
The library is written in fixed format Fortran 90.  In addition to a Fortran
90/95 compiler, the following libraries are required.

  1. MPI, e.g., OpenMPI or MPICH.
  2. An optimized BLAS library, e.g., ATLAS or OpenBLAS,
     see <http://www.netlib.org/blas/> for a reference implementation.
  3. LAPACK, see <http://www.netlib.org/lapack/>
  4. ScaLAPACK (including BLACS and PBLAS),
     see <http://www.netlib.org/scalapack/>

Follow the instruction below to build the library:

  1. Download and unpack the BSEPACK archive.
  2. Modify the file `make.inc` to setup the compilers, linkers, and
     libraries.  Templates of this file are provided in the directory
     `MAKE_INC/`.
  3. Type `make all` to build the library and test programs.
  4. Run the test programs to check whether the installation is successful.

More detailed instructions regarding installation and testing can be found in
*BSEPACK User's Guide*.


3. Functionalities
------------------
The BSEPACK library provides two main functionalities:

  1. Diagonalizing the definite Bethe--Salpeter Hamiltonian matrix.
  2. Compute or estimate the absorption spectrum.

Solvers based on Tamm--Dancoff approximation (TDA) for both cases are also
provided.  Current release supports both real and complex data types, but only
for double precision.  Support to single precision real and complex data types
are *not* planned in future releases unless strong requests are received.

Examples for computing the eigenvalues and absorption spectrum are provided in
the directory `EXAMPLES/`.  The outputs are valid Matlab/Octave scripts that
can be used directly in Matlab/Octave for postprocessing.  To use your own
matrices, replace the input files (`EXAMPLES/input_real.txt` and/or
`EXAMPLES/input_complex.txt`) by yours.  Detailed usage is described in
*BSEPACK User's Guide*.


4. Comments/Questions/Bug Reports
---------------------------------
Please send your request to <myshao@lbl.gov>.


5. Selected References
----------------------
The dense eigensolver in BSEPACK is implemented based on [1].
The estimation of absorption spectrum uses algorithms in [2,3].

  [1] M. Shao, F. H. da Jornada, C. Yang, J. Deslippe, and S. G. Louie.
      Structure preserving parallel algorithms for solving the
      Bethe--Salpeter eigenvalue problem.
      Linear Algebra Appl., 488:148--167, 2016.
      DOI: 10.1016/j.laa.2015.09.036

  [2] J. Brabec, L. Lin, M. Shao, N. Govind, Y. Saad, C. Yang, and E. G. Ng.
      Efficient algorithms for estimating the absorption spectrum within
      linear response TDDFT.
      J. Chem. Theory Comput., 11(11):5197--5208, 2015.
      DOI: 10.1021/acs.jctc.5b00887

  [3] M. Shao, F. H. da Jornada, L. Lin, C. Yang, J. Deslippe, and S. G. Louie.
      A structure preserving Lanczos algorithm for computing the optical
      absorption spectrum.
      Avaliable as arXiv:1611.02348, 2016.
