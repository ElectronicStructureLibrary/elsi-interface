### Source file ###
LIST(APPEND ntpoly_src
  src/ChebyshevSolversModule.f90
  src/DataTypesModule.f90
  src/DensityMatrixSolversModule.f90
  src/DMatrixModule.f90
  src/EigenBoundsModule.f90
  src/ErrorModule.f90
  src/ExponentialSolversModule.f90
  src/GemmTasksModule.f90
  src/GeometryOptimizationModule.f90
  src/HermiteSolversModule.f90
  src/Interface.f90
  src/InverseSolversModule.f90
  src/LinearSolversModule.f90
  src/LoadBalancerModule.f90
  src/LoggingModule.f90
  src/MatrixMapsModule.f90
  src/MatrixMarketModule.f90
  src/MatrixMemoryPoolModule.f90
  src/NTMPIModule.f90
  src/PermutationModule.f90
  src/PMatrixMemoryPoolModule.f90
  src/PolynomialSolversModule.f90
  src/ProcessGridModule.f90
  src/PSMatrixAlgebraModule.f90
  src/PSMatrixModule.f90
  src/RootSolversModule.f90
  src/SignSolversModule.f90
  src/SMatrixAlgebraModule.f90
  src/SMatrixModule.f90
  src/SolverParametersModule.f90
  src/SquareRootSolversModule.f90
  src/SVectorModule.f90
  src/TimerModule.f90
  src/TrigonometrySolversModule.f90
  src/TripletListModule.f90
  src/TripletModule.f90)

IF(USE_MPI_IALLGATHER)
  LIST(APPEND ntpoly_src src/MatrixReduceModule.f90)
ELSE()
  LIST(APPEND ntpoly_src src/MatrixReduceModuleNoIallgather.f90)
ENDIF()

ADD_LIBRARY(NTPoly ${ntpoly_src})

TARGET_LINK_LIBRARIES(NTPoly PRIVATE ${LIBS})

TARGET_INCLUDE_DIRECTORIES(NTPoly PRIVATE
  ${INC_PATHS}
  ${CMAKE_Fortran_MODULE_DIRECTORY})

INSTALL(TARGETS NTPoly
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
