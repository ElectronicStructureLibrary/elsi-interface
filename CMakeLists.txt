cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

### Project ###
project(elsi
  VERSION 2018.02
  LANGUAGES Fortran C)
set(elsi_URL "http://elsi-interchange.org")
set(elsi_EMAIL "elsi-team@duke.edu")
set(elsi_LICENSE "BSD-3")
set(elsi_DESCRIPTION "Electronic Structure Infrastructure")

### Installation paths ###
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(ARCHIVE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
if(NOT DEFINED CMAKE_INSTALL_INCLUDEDIR)
  set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "Location to install module and header files")
endif()
if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "Location to install library files")
endif()
if(NOT DEFINED CMAKE_INSTALL_BINDIR)
  set(CMAKE_INSTALL_BINDIR "bin" CACHE PATH "Location to install binary files")
endif()

### Options ###
option(ENABLE_TESTS "Build ELSI tests" OFF)
option(ENABLE_PEXSI "Enable PEXSI" OFF)
option(ENABLE_SIPS "Enable SIPs" OFF)
option(ENABLE_ELPA_BGQ "Use ELPA BlueGene Q kernel" OFF)
option(ENABLE_ELPA_AVX "Use ELPA AVX kernel" OFF)
option(ENABLE_ELPA_AVX2 "Use ELPA AVX2 kernel" OFF)
option(ENABLE_ELPA_AVX512 "Use ELPA AVX512 kernel" OFF)

if(ENABLE_PEXSI)
  enable_language(CXX)
endif()

### MPI ###
find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
  message(STATUS "Found MPI CXX compiler: ${MPI_CXX_COMPILER}")
  include_directories(${MPI_CXX_INCLUDE_PATH})
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
#  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS})
#  set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES})
elseif(ENABLE_PEXSI)
  message(FATAL_ERROR "MPI CXX compiler not found")
endif()
if(MPI_C_FOUND)
  message(STATUS "Found MPI C compiler: ${MPI_C_COMPILER}")
  include_directories(${MPI_C_INCLUDE_PATH})
  set(CMAKE_CXX_FLAGS ${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS})
#  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS})
#  set(LIBS ${LIBS} ${MPI_C_LIBRARIES})
elseif(ENABLE_C)
  message(FATAL_ERROR "MPI C compiler not found")
endif()
if(MPI_Fortran_FOUND)
  message(STATUS "Found MPI Fortran compiler: ${MPI_Fortran_COMPILER}")
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
#  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS})
#  set(LIBS ${LIBS} ${MPI_Fortran_LIBRARIES})
else()
  message(FATAL_ERROR "MPI Fortran compiler not found")
endif()

### Compilations ###
include(elpa)
if(NOT EXTERNAL_ELPA)
  add_subdirectory(external/ELPA)
endif()

include(libomm)
if(NOT EXTERNAL_OMM)
  add_subdirectory(external/libOMM)
endif()

if(ENABLE_PEXSI)
  include(pexsi)

  if(NOT EXTERNAL_SUPERLU)
    add_subdirectory(external/PtScotch)
    add_subdirectory(external/SuperLU_DIST)
  endif()

  message(STATUS "Enabling PEXSI")
  add_subdirectory(external/PEXSI)
endif()

if(ENABLE_SIPS)
  include(sips)
  message(STATUS "Enabling SIPs")
  add_subdirectory(external/SIPs)
endif()

add_subdirectory(src)

if(ENABLE_TESTS)
  add_subdirectory(test)
endif()
