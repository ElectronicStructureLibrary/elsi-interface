# How to install ELSI

## Prerequisites

The installation of ELSI makes use of the CMake software. Minimum requirements:

* CMake (3.0.2 or newer)
* Fortran compiler (Fortran 2003 compliant)
* C compiler (C99 compliant)
* MPI (MPI-3 recommended)
* BLAS, LAPACK, BLACS, ScaLAPACK

In addition, building the PEXSI solver (highly recommended) requires:

* C++ compiler (C++11 compliant)

By default, solver libraries and their dependencies redistributed within ELSI
will be built. Optionally, they may be substituted by user's optimized versions:

* ELPA (2018.05 or newer)
* libOMM
* PEXSI (1.2 or newer)
* NTPoly (2.2 or newer)

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

The ELPA, libOMM, NTPoly, and PEXSI solver libraries, as well as the
SuperLU\_DIST and PT-SCOTCH libraries are redistributed through this ELSI
package. Experienced users are encouraged to link the ELSI interface against
external, better optimized solver libraries. Relevant options are:

* `USE_EXTERNAL_ELPA`
* `USE_EXTERNAL_OMM`
* `USE_EXTERNAL_PEXSI`
* `USE_EXTERNAL_NTPOLY`

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

Please note that in the current version of ELSI, the PEXSI solver is not enabled
unless `ENABLE_PEXSI` is switched on.

### 3) Shared library

Building ELSI as a shared library may be enabled by `BUILD_SHARED_LIBS`.

### 4) Tests

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

Copyright (c) 2015-2019, the ELSI team. All rights reserved.
