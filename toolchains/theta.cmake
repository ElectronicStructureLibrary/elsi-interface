### Theta ###

# test_fortran_11 and test_fortran_13 tests are failing.  Both are real libOMM tests and run to
# completion, but they don't print out "Passed."

SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "MPI Fortran compiler")
SET(CMAKE_C_COMPILER "cc" CACHE STRING "MPI C compiler")
SET(CMAKE_CXX_COMPILER "CC" CACHE STRING "MPI C++ compiler")

SET(CMAKE_Fortran_FLAGS "-fast -no-ipo -g" CACHE STRING "Fortran flags")
SET(CMAKE_C_FLAGS "-fast -no-ipo -g -std=c99" CACHE STRING "C flags")
SET(CMAKE_CXX_FLAGS "-fast -no-ipo -g -std=c++11" CACHE STRING "C++ flags")

SET(ELPA2_KERNEL "AVX512" CACHE STRING "Use ELPA AVX512 kernel")
SET(ENABLE_PEXSI ON CACHE BOOL "Enable PEXSI")
SET(ADD_UNDERSCORE ON CACHE BOOL "Suffix C functions with an underscore")

# If you use testing, you will receive a warning that "Linear algebra libraries not provided".
# This warning message may be ignored: the linear algebra libraries are included in the compiler wrapper.
# You will need to run the test calculations using an interactive job, as the login node does not support
# MPI calculations.
SET(ENABLE_TESTS ON CACHE BOOL "Build ELSI Fortran test programs")
SET(ENABLE_C_TESTS ON CACHE BOOL "Build ELSI C test programs")
SET(MPIEXEC "aprun" CACHE STRING "Name of MPI executable for running tests")
SET(ALWAYS_USE_MPIEXEC ON CACHE BOOL "Always run tests with the MPI executable, even for serial calculations")

SET(PTSCOTCH_DIR "/home/huhn/opt/scotch/6.0.5a" CACHE PATH "PtScotch directory")
SET(PTSCOTCH_LIB "${PTSCOTCH_DIR}/lib/libptscotchparmetis.a;${PTSCOTCH_DIR}/lib/libptscotch.a;${PTSCOTCH_DIR}/lib/libptscotcherr.a;${PTSCOTCH_DIR}/lib/libscotchmetis.a;${PTSCOTCH_DIR}/lib/libscotch.a;${PTSCOTCH_DIR}/lib/libscotcherr.a" CACHE STRING "PtScotch libraries")
SET(PTSCOTCH_INC "${PTSCOTCH_DIR}/include" CACHE STRING "PtScotch include directory")
