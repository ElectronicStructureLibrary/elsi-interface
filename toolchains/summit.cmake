### Summit ###

SET(CMAKE_Fortran_COMPILER "mpixlf" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "mpixlc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "mpixlC" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-Ofast -qarch=pwr9 -qstrict" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-Ofast -qarch=pwr9 -qstrict -qlanglvl=stdc99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-Ofast -qarch=pwr9 -qstrict -qlanglvl=extended0x" CACHE STRING "C++ flags")
SET(CMAKE_CUDA_FLAGS "-O3 -arch=sm_70" CACHE STRING "CUDA flags")

SET(USE_GPU_CUDA ON CACHE BOOL "Use CUDA-based GPU acceleration in ELPA")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")
SET(ADD_UNDERSCORE OFF CACHE BOOL "Do not suffix C functions with an underscore")

SET(LIB_PATHS "$ENV{OLCF_ESSL_ROOT}/lib64;$ENV{OLCF_CUDA_ROOT}/lib64" CACHE STRING "External library paths")
SET(LIBS "essl;scalapack;lapack;cublas;cudart" CACHE STRING "External libraries")
