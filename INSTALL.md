# How to install ELSI

## Prerequisites

The installation of ELSI makes use of the CMake software. Minimum requirements:

* CMake (3.0.2 or newer)
* Fortran compiler (Fortran 2003)
* C compiler (C99)
* MPI
* BLAS, LAPACK, BLACS, ScaLAPACK

In addition, building the PEXSI solver (highly recommended) requires:

* C++ compiler (C++11)

By default, solver libraries and their dependencies redistributed within ELSI
will be built. Optionally, they may be substituted by user's optimized versions:

* ELPA (2018.05.001 or newer)
* libOMM
* MatrixSwitch
* NTPoly (2.2 or newer)
* PT-SCOTCH (6.0.0)
* SuperLU\_DIST (5.1.3 or newer)

## CMake basics

The typical workflow of using CMake to build ELSI looks like:

    $ ls
      CMakeLists.txt  external/  src/  test/  ...
    $ mkdir build
    $ cd build
    $ cmake [options] ..
      ...
      ...
      -- Generating done
      -- Build files have been written to: /current/dir
    $ make [-j np]
    $ make install

Whenever CMake is invoked, one of the command line arguments must point to the
path where the top level CMakeLists.txt file exists, hence the `..` in the above
example.

An option may be defined by adding

    -DKeyword=Value

to the command line when invoking CMake. If `Keyword` is of type boolean, its
`Value` can be `ON` or `OFF`. If `Keyword` is a list of libraries or paths, the
items should be separated with `;` (semicolon) or ` ` (space). For example,

    -DCMAKE_INSTALL_PREFIX=/path/to/install/elsi
    -DCMAKE_C_COMPILER=gcc
    -DENABLE_TESTS=OFF
    -DENABLE_PEXSI=ON
    -DINC_PATHS="/path/to/include;/another/path/to/include"
    -DLIBS="library1 library2 library3"

## Configuration

### 1) Compilers

CMake automatically detects and sets default compilers. The choices made by
CMake often work, but not necessarily guarantee the optimal performance. In some
cases, the compilers selected by CMake may not be the ones desired by the user.
Therefore, it is mandatory that the user sets the identification of compilers
explicitly:

    -DCMAKE_Fortran_COMPILER=[YOUR_MPI_FORTRAN_COMPILER]
    -DCMAKE_C_COMPILER=[YOUR_MPI_C_COMPILER]
    -DCMAKE_CXX_COMPILER=[YOUR_MPI_C++_COMPILER]

The C++ compiler is not needed if building ELSI without PEXSI.

It is highly recommended to specify the compiler flags, in particular the
optimization flags:

    -DCMAKE_Fortran_FLAGS=[YOUR_FORTRAN_COMPILE_FLAGS]
    -DCMAKE_C_FLAGS=[YOUR_C_COMPILE_FLAGS]
    -DCMAKE_CXX_FLAGS=[YOUR_C++_COMPILE_FLAGS]

### 2) Solvers

The ELPA, libOMM, NTPoly, and PEXSI solver libraries, as well as the
SuperLU\_DIST and PT-SCOTCH libraries are redistributed through this ELSI
package. Experienced users are encouraged to link the ELSI interface against
external, better optimized solver libraries. Relevant options are:

    -DUSE_EXTERNAL_ELPA=ON
    -DUSE_EXTERNAL_OMM=ON
    -DUSE_EXTERNAL_SUPERLU=ON
    -DUSE_EXTERNAL_NTPOLY=ON

External libraries and their paths should be set via the following keywords:

    -DLIB_PATHS=[ALL_DIRECTORIES_CONTAINING_YOUR_EXTERNAL_LIBRARIES]
    -DINC_PATHS=[ALL_INCLUDE_DIRECTORIES_OF_YOUR_EXTERNAL_LIBRARIES]
    -DLIBS=[NAMES_OF_ALL_YOUR_EXTERNAL_LIBRARIES]

Please note that in the current version of ELSI, the PEXSI solver is not enabled
by default. It may be switched on by specifying:

    -DENABLE_PEXSI=ON

At present, an external version of PEXSI is not officially supported.

### 3) Build targets

By default, a static library will be created as the target of the ELSI
installation. Building ELSI as a shared library may be enabled by:

    -DBUILD_SHARED_LIBS=ON

Building ELSI test programs may be enabled by:

    -DENABLE_TESTS=ON

In either case, linear algebra libraries BLAS, LAPACK, BLACS, and ScaLAPACK
should be present in LIB\_PATHS and LIBS.

If the test programs are turned on, the installation may be verified by

    $ ...
    $ make [-j np]
    $ make test

or

    $ ...
    $ make [-j np]
    $ ctest

Note that the tests may not run if launching MPI jobs is prohibited on the
user's working platform.

### 4) Toolchain files

It is sometimes convenient to edit the settings in a toolchain file that can be
read in by CMake:

    -DCMAKE_TOOLCHAIN_FILE=[YOUR_TOOLCHAIN_FILE]

Example toolchains are provided in the [`./toolchains`](./toolchains) directory,
which the user may use as templates to create new ones.

## Beyond this guide

A complete description of the ELSI build system is available in
[`./doc/elsi_manual.pdf`](./doc/elsi_manual.pdf).

## Troubleshooting

For comments, feedback, and suggestions, please
[contact the ELSI team](mailto:elsi-team@duke.edu).

Copyright (c) 2015-2019, the ELSI team. All rights reserved.