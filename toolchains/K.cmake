### K computer ###

SET(CMAKE_Fortran_COMPILER "mpifrtpx" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpifccpx" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpiFCCpx" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-Kfast,-Kparallel,openmp -SCALAPACK -SSL2" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-Kfast,-Kparallel,openmp -SCALAPACK -SSL2" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-Kfast,-Kparallel,openmp -SCALAPACK -SSL2" CACHE STRING "C++ flags")

# PEXSI cannot be compiled due to the lack of C++11 support on K
SET(ENABLE_PEXSI OFF CACHE BOOL "Disable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")

SET(CMAKE_Fortran_MODDIR_FLAG "-M")
