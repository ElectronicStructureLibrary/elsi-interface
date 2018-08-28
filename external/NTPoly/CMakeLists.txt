### Source file ###
LIST(APPEND ntpoly_src
  src/ChebyshevSolversModule.F90
  src/DataTypesModule.F90
  src/DensityMatrixSolversModule.F90
  src/DMatrixModule.F90
  src/EigenBoundsModule.F90
  src/ErrorModule.F90
  src/ExponentialSolversModule.F90
  src/GemmTasksModule.F90
  src/GeometryOptimizationModule.F90
  src/HermiteSolversModule.F90
  src/Interface.F90
  src/InverseSolversModule.F90
  src/LinearSolversModule.F90
  src/LoadBalancerModule.F90
  src/LoggingModule.F90
  src/MatrixMarketModule.F90
  src/MatrixMemoryPoolModule.F90
  src/MatrixReduceModule.F90
  src/MinimizerSolversModule.F90
  src/PermutationModule.F90
  src/PMatrixMemoryPoolModule.F90
  src/PolynomialSolversModule.F90
  src/ProcessGridModule.F90
  src/PSMatrixAlgebraModule.F90
  src/PSMatrixModule.F90
  src/RootSolversModule.F90
  src/SignSolversModule.F90
  src/SMatrixAlgebraModule.F90
  src/SMatrixModule.F90
  src/SolverParametersModule.F90
  src/SquareRootSolversModule.F90
  src/SVectorModule.F90
  src/TimerModule.F90
  src/TrigonometrySolversModule.F90
  src/TripletListModule.F90
  src/TripletModule.F90)

ADD_LIBRARY(NTPoly ${ntpoly_src})

TARGET_INCLUDE_DIRECTORIES(NTPoly PRIVATE ${INC_PATHS})

INSTALL(TARGETS NTPoly
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
