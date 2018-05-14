### Theta ###

SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "cc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "CC" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-fast -no-ipo" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-fast -no-ipo -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-fast -no-ipo -std=c++11" CACHE STRING "C++ flags")

SET(ELPA2_KERNEL "AVX512" CACHE STRING "Use ELPA AVX512 kernel")
SET(ENABLE_TESTS ON CACHE BOOL "Enable tests")
SET(ENABLE_PEXSI OFF CACHE BOOL "Enable PEXSI")

# Currently external SuperLU_DIST and PT-SCOTCH libraries are needed on Theta
SET(LIB_PATHS "/home/wzvyu/theta/opt/SuperLU_DIST_5.3.0/lib /home/wzvyu/theta/opt/scotch_6.0.5a/lib" CACHE STRING "External library paths")
SET(INC_PATHS "/home/wzvyu/theta/opt/SuperLU_DIST_5.3.0/SRC /home/wzvyu/theta/opt/scotch_6.0.5a/include" CACHE STRING "External library include paths")
SET(LIBS "superlu_dist ptscotchparmetis ptscotch ptscotcherr scotchmetis scotch scotcherr" CACHE STRING "External libraries")
