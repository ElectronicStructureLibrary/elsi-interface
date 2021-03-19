# ELSI - ELectronic Structure Infrastructure (v2.8.0)

## About

ELSI is a unified software interface designed for electronic structure codes to
connect with various high-performance eigensolvers and density matrix solvers.
For more information, visit the [ELSI interchange](https://elsi-interchange.org)
website.

## Installation

The standard installation of ELSI requires:

* CMake (3.0 or newer)
* Fortran compiler (Fortran 2003)
* C compiler (C99)
* C++ compiler (C++11, optional)
* MPI (MPI-3)
* BLAS, LAPACK, ScaLAPACK
* CUDA (optional)

Installation with recent versions of Cray, GNU, IBM, Intel, and NVIDIA (formerly
PGI) compilers has been tested. For a complete description of the installation
process, please refer to [`./INSTALL.md`](./INSTALL.md).

## More

A User's Guide is available at [`./doc/elsi_manual.pdf`](./doc/elsi_manual.pdf).
For comments, feedback, and suggestions, please
[contact the ELSI team](mailto:elsi-team@duke.edu).

Copyright (c) 2015-2021, the ELSI team. All rights reserved.
