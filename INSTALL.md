# How to install ELSI

## Prerequisites

The installation of ELSI makes use of the CMake software. Minimum requirements:

* CMake (3.0.2 or newer)
* Fortran compiler (Fortran 2003 compliant)
* C compiler (C99 compliant)
* MPI (MPI-3 recommended)
* BLAS, LAPACK, ScaLAPACK (with PBLAS and BLACS)

Enabling the PEXSI solver (highly recommended) requires:

* C++ compiler (C++11 compliant)

Enabling the EigenExa solver requires:

* EigenExa (2.3 or newer)

Enabling the SLEPc-SIPs solver requires:

* SLEPc (3.9 or newer)

Enabling the MAGMA solver requires:

* MAGMA (2.5 or newer)

## Quick Start

We recommend preparing and editing configuration settings in a toolchain file
that can be read by CMake. Edit one of the templates provided in the
[`./toolchains`](./toolchains) directory. Then follow the steps below:

    $ mkdir build
    $ cd build
    $ cmake -DCMAKE_TOOLCHAIN_FILE=YOUR_TOOLCHAIN_FILE ..
      ...
      ...
      -- Generating done
      -- Build files have been written to: /current/dir
    $ make [-j np]
    $ [make install]

`YOUR_TOOLCHAIN_FILE` should be the user's toolchain file. Commands in square
brackets are optional.

## Configuration

### 1) Compilers

CMake automatically detects and sets default compilers. The choices made by
CMake often work, but not necessarily guarantee the optimal performance. In some
cases, the compilers selected by CMake may not be the ones desired by the user.
Therefore, it is mandatory that the user explicitly sets the identification of
compilers in the following keywords:

* `CMAKE_Fortran_COMPILER`
* `CMAKE_C_COMPILER`
* `CMAKE_CXX_COMPILER` (only needed if building ELSI with PEXSI support)

It is highly recommended to specify the compiler flags, in particular the
optimization flags:

* `CMAKE_Fortran_FLAGS`
* `CMAKE_C_FLAGS`
* `CMAKE_CXX_FLAGS`

### 2) Solvers

The ELPA, libOMM, NTPoly, BSEPACK, and PEXSI solver libraries, as well as the
SuperLU\_DIST and PT-SCOTCH libraries (both required by PEXSI) are redistributed
through this ELSI package. Experienced users are encouraged to link the ELSI
interface against externally installed, better optimized solver libraries.
Relevant options are:

* `USE_EXTERNAL_ELPA`
* `USE_EXTERNAL_OMM`
* `USE_EXTERNAL_PEXSI`
* `USE_EXTERNAL_NTPOLY`
* `USE_EXTERNAL_BSEPACK`

All external libraries and include paths should be set via:

* `INC_PATHS`
* `LIB_PATHS`
* `LIBS`

Each of the above is a list separated by ` ` (space) or `;` (semicolon).
`INC_PATHS` and `LIB_PATHS` should be absolute paths. `LIBS` accepts three
formats:

* `-lelpa;-lpexsi;-lblas` (link line style)
* `elpa;pexsi;blas` (name of library)
* `libelpa.a;libpexsi.a;libblas.so` (full name of library)

Please note that in the current version of ELSI, the PEXSI and BSEPACK solvers
are not enabled by default. They may be switched on by `ENABLE_PEXSI` and
`ENABLE_BSEPACK`, respectively. In addition, `ENABLE_SIPS`, `ELABLE_EIGENEXA`,
and `ENABLE_MAGMA` may be used to enable support for the SLEPc, EigenExa, and
MAGMA solvers, respectively. These libraries are not redistributed with ELSI,
thus must be installed separately by the user.

### 3) Tests

Building ELSI test programs may be enabled by `ENABLE_TESTS`. Then, the
compilation may be verified by

    $ make test

or

    $ ctest

The tests may fail to start if launching MPI jobs is prohibited on the user's
working platform.

## Beyond this guide

A complete description of the ELSI build system is available in
[`./doc/elsi_manual.pdf`](./doc/elsi_manual.pdf).

## Troubleshooting

For comments, feedback, and suggestions, please
[contact the ELSI team](mailto:elsi-team@duke.edu).

Copyright (c) 2015-2020, the ELSI team. All rights reserved.
