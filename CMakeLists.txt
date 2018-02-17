### TODO ###
# GitLab CI
# pkg-config
# PtScotch support

#  <<<  CMake build overview  >>>
#
#    $ ls
#      CMakeLists.txt external/  src/  test/ ...
#    $ mkdir build
#    $ cd build
#    $ cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/elsi [options ...]
#      ...
#      -- Generating done
#      -- Build files have been written to: /current/dir
#    $ make [-j np]
#    $ make install
#
#  <<<  Required dependencies  >>>
#
#    - CMake 3.0 or newer
#    - Fortran compiler (Fortran 2003)
#    - C compiler (C99)
#
#  <<<  Optional dependencies  >>>
#
#    - ELPA 2017.05.003 or newer      "Needed when internal version is not wanted"
#    - libOMM and MatrixSwitch        "Needed when internal version is not wanted"
#    - C++ compiler (C++11)           "Needed when PEXSI is enabled"
#    - SuperLU_DIST 5.1.2 or newer    "Needed when PEXSI is enabled and internal version is not wanted"
#    - PtScotch 6.0.0 or newer        "Needed when PEXSI is enabled"
#    - BLAS/LAPACK(/BLACS)/ScaLAPACK  "Needed when ELSI tests are enabled"
#
#  <<<  Compilers and flags  >>>
#
#    - CMAKE_Fortran_COMPILER  "Fortran compiler"
#    - CMAKE_Fortran_FLAGS     "Fortran flags"
#    - CMAKE_C_COMPILER        "C compiler"
#    - CMAKE_C_FLAGS           "C flags"
#    - CMAKE_CXX_COMPILER      "C++ compiler"
#    - CMAKE_CXX_FLAGS         "C++ flags"
#
#    * -std=c99, -std=c++11 ect are automatically set
#
#  <<<  Comfiguration options  >>>
#
#    - CMAKE_TOOLCHAIN_FILE  "Use pre-defined toolchain file"
#    - BUILD_SHARED_LIBRARY  "Build ELSI as a shared library (todo)
#    - ENABLE_TESTS          "Build ELSI test programs"
#    - ENABLE_PEXSI          "Enable PEXSI (requires PtScotch library)"
#    - ENABLE_SIPS           "Enable SIPs (todo)"
#    - ENABLE_ELPA_BGQ       "When ELPA library not provided, build internal ELPA BGQ kernel"
#    - ENABLE_ELPA_AVX       "When ELPA library not provided, build internal ELPA AVX kernel"
#    - ENABLE_ELPA_AVX2      "When ELPA library not provided, build internal ELPA AVX2 kernel"
#    - ENABLE_ELPA_AVX512    "When ELPA library not provided, build internal ELPA AVX512 kernel"
#    - ELPA_LIB              "External ELPA libraries as a ";"-separated list"
#    - ELPA_INC              "External ELPA include directories as a ";"-separated list"
#    - OMM_LIB               "External libOMM and MatrixSwitch libraries as a ";"-separated list"
#    - OMM_INC               "External libOMM and MatrixSwitch include directories as a ";"-separated list"
#    - SUPERLU_LIB           "External SUPERLU_DIST libraries as a ";"-separated list"
#    - SUPERLU_INC           "External SUPERLU_DIST include directories as a ";"-separated list"
#    - PTSCOTCH_LIB          "External PtScotch libraries as a ";"-separated list"
#    - PTSCOTCH_INC          "External PtScotch include directories as a ";"-separated list"
#    - MATH_LIB              "External BLAS, LAPACK, (BLACS,) ScaLAPACK libraries as a ";"-separated list"
#    - MATH_INC              "External BLAS, LAPACK, (BLACS,) ScaLAPACK include directories as a ";"-separated list"
#
#    * If ENABLE_TESTS=ON, MATH_LIB must be valid
#    * If ENABLE_PEXSI=ON, PTSCOTCH_LIB and PTSCOTCH_INC must be valid
#    * If and only if both "xxx_LIB" and "xxx_INC" are provided, internal version will not be built
#
#  <<<  Installation  >>>
#
#    - CMAKE_INSTALL_PREFIX      "Location to install ELSI"
#    - CMAKE_INSTALL_BINDIR      "Location within CMAKE_INSTALL_PREFIX to install executables"
#    - CMAKE_INSTALL_LIBDIR      "Location within CMAKE_INSTALL_PREFIX to install libraries"
#    - CMAKE_INSTALL_INCLUDEDIR  "Location within CMAKE_INSTALL_PREFIX to install modules and headers"
#
#    * Defaults are: usr/local, CMAKE_INSTALL_PREFIX/bin, CMAKE_INSTALL_PREFIX/lib, CMAKE_INSTALL_PREFIX/include

CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)

### Project ###
PROJECT(elsi
  VERSION 2018.02
  LANGUAGES Fortran C)
SET(elsi_URL "http://elsi-interchange.org")
SET(elsi_EMAIL "elsi-team@duke.edu")
SET(elsi_LICENSE "BSD 3")
SET(elsi_DESCRIPTION "Electronic Structure Infrastructure")

### Installation paths ###
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include")
SET(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
SET(ARCHIVE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
IF(NOT DEFINED CMAKE_INSTALL_INCLUDEDIR)
  SET(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Location to install module and header files")
ENDIF()
IF(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  SET(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Location to install library files")
ENDIF()
IF(NOT DEFINED CMAKE_INSTALL_BINDIR)
  SET(CMAKE_INSTALL_BINDIR "bin" CACHE PATH "Location to install binary files")
ENDIF()
INCLUDE(color_message)
MESSAGE(STATUS "ELSI will be installed to: ${CMAKE_INSTALL_PREFIX}")

### Options ###
OPTION(ENABLE_TESTS "Build ELSI Fortran test programs" OFF)
OPTION(ENABLE_C_TESTS "Build ELSI C test programs" OFF)
OPTION(ENABLE_PEXSI "Enable PEXSI" OFF)
OPTION(ENABLE_SIPS "Enable SIPs" OFF)
OPTION(ENABLE_ELPA_BGQ "Use ELPA BlueGene Q kernel" OFF)
OPTION(ENABLE_ELPA_AVX "Use ELPA AVX kernel" OFF)
OPTION(ENABLE_ELPA_AVX2 "Use ELPA AVX2 kernel" OFF)
OPTION(ENABLE_ELPA_AVX512 "Use ELPA AVX512 kernel" OFF)

IF(ENABLE_PEXSI)
  ENABLE_LANGUAGE(CXX)
ENDIF()

### MPI ###
FIND_PACKAGE(MPI REQUIRED)
IF(MPI_CXX_FOUND)
  MESSAGE(STATUS "Found ${GREEN}MPI CXX compiler${COLORRESET}: ${MPI_CXX_COMPILER}")
  INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
  SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
  SET(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS})
ELSEIF(ENABLE_PEXSI)
  MESSAGE(FATAL_ERROR "${MAGENTA}MPI CXX compiler not found${COLORRESET}")
ENDIF()
IF(MPI_C_FOUND)
  MESSAGE(STATUS "Found ${GREEN}MPI C compiler${COLORRESET}: ${MPI_C_COMPILER}")
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  SET(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS})
  SET(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS})
ELSEIF(ENABLE_C)
  MESSAGE(FATAL_ERROR "${MAGENTA}MPI C compiler not found${COLORRESET}")
ENDIF()
IF(MPI_Fortran_FOUND)
  MESSAGE(STATUS "Found ${GREEN}MPI Fortran compiler${COLORRESET}: ${MPI_Fortran_COMPILER}")
  INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
  SET(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
  SET(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS})
ELSE()
  MESSAGE(FATAL_ERROR "${MAGENTA}MPI Fortran compiler not found${COLORRESET}")
ENDIF()

### Compilations ###
INCLUDE(elpa)
IF(NOT ELPA_FOUND)
  ADD_SUBDIRECTORY(external/ELPA)
  SET(ELPA_LIB "${CMAKE_BINARY_DIR}/lib/libelpa.a")
ENDIF()

INCLUDE(libomm)
IF(NOT OMM_FOUND)
  ADD_SUBDIRECTORY(external/libOMM)
  SET(OMM_LIB "${CMAKE_BINARY_DIR}/lib/libOMM.a;${CMAKE_BINARY_DIR}/lib/libMatrixSwitch.a")
ENDIF()

IF(ENABLE_PEXSI)
  MESSAGE(STATUS "Enabling PEXSI")

  INCLUDE(superlu)
  IF(NOT SUPERLU_FOUND)
    ADD_SUBDIRECTORY(external/SuperLU_DIST)
    SET(SUPERLU_LIB "${CMAKE_BINARY_DIR}/lib/libsuperlu_dist.a")
  ENDIF()

  ADD_SUBDIRECTORY(external/PEXSI)
  SET(PEXSI_LIB "${CMAKE_BINARY_DIR}/lib/libpexsi.a")
ENDIF()

IF(ENABLE_SIPS)
  INCLUDE(sips)
  MESSAGE(STATUS "Enabling SIPs")
  ADD_SUBDIRECTORY(external/SIPs)
  SET(SIPS_LIB "${CMAKE_BINARY_DIR}/lib/libsips.a")
ENDIF()

ADD_SUBDIRECTORY(src)

IF(ENABLE_TESTS)
  MESSAGE(STATUS "Enabling Fortran test programs")
  IF(ENABLE_C_TESTS)
    MESSAGE(STATUS "Enabling C test programs")
  ENDIF()
  ENABLE_TESTING()
  INCLUDE(mathlib)
  INCLUDE_DIRECTORIES(src)
  ADD_SUBDIRECTORY(test)
ENDIF()
