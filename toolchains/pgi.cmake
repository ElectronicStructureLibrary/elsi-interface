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

SET(LIB_PATHS "/home/yy244/compilers/mpich-3.2.1/install_pgi19.4/scalapack_libs/scalapack-2.1.0 /home/yy244/compilers/mpich-3.2.1/install_pgi19.4/scalapack_libs/lapack-3.10.0" CACHE STRING "External library paths")
SET(LIBS "scalapack lapack refblas" CACHE STRING "External libraries")

SET(MPIEXEC_1P "mpirun -n 1" CACHE STRING "Command to run serial tests with 1 MPI task")
SET(MPIEXEC_NP "mpirun -n 4" CACHE STRING "Command to run parallel tests with multiple MPI tasks")
