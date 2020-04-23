### PGI ###

SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-O2" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-O1 -c99" CACHE STRING "C flags") # SCOTCH segfaults with "O2"
SET(CMAKE_CXX_FLAGS "-O2 --c++11" CACHE STRING "C++ flags")

SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable Fortran tests")
SET(ENABLE_C_TESTS ON CACHE BOOL "Enable C tests")

SET(LIB_PATHS "/home/wy29/opt/pgi/linux86-64-llvm/19.10/lib/scalapack/scalapack-2.0.2/openmpi-3.1.3/lib /home/wy29/opt/pgi/linux86-64-llvm/19.10/lib" CACHE STRING "External library paths")
SET(LIBS "scalapack lapack blas" CACHE STRING "External libraries")

SET(MPIEXEC_1P "mpirun --mca io romio314 -n 1" CACHE STRING "Command to run serial tests with 1 MPI task")
SET(MPIEXEC_NP "mpirun --mca io romio314 -n 4" CACHE STRING "Command to run parallel tests with multiple MPI tasks")
