### Source file ###
LIST(APPEND ntpoly_src
  src/ChebyshevSolversModule.f90
  src/DataTypesModule.f90
  src/DenseMatrixModule.f90
  src/DensityMatrixSolversModule.f90
  src/DistributedMatrixMemoryPoolModule.f90
  src/DistributedSparseMatrixModule.f90
  src/DistributedSparseMatrixAlgebraModule.f90
  src/EigenBoundsModule.f90
  src/ErrorModule.f90
  src/ExponentialSolversModule.f90
  src/FixedSolversModule.f90
  src/GemmTasksModule.f90
  src/GeometryOptimizationModule.f90
  src/HermiteSolversModule.f90
  src/Interface.f90
  src/InverseSolversModule.f90
  src/IterativeSolversModule.f90
  src/LinearSolversModule.f90
  src/LoadBalancerModule.f90
  src/LoggingModule.f90
  src/MatrixGatherModule.f90
  src/MatrixMarketModule.f90
  src/MatrixMemoryPoolModule.f90
  src/MinimizerSolversModule.f90
  src/ParameterConverterModule.f90
  src/PermutationModule.f90
  src/PolynomialSolversModule.f90
  src/ProcessGridModule.f90
  src/RootSolversModule.f90
  src/SignSolversModule.f90
  src/SparseMatrixAlgebraModule.f90
  src/SparseMatrixIteratorModule.f90
  src/SparseMatrixModule.f90
  src/SparseVectorModule.f90
  src/SquareRootSolversModule.f90
  src/TimerModule.f90
  src/TrigonometrySolversModule.f90
  src/TripletListModule.f90
  src/TripletModule.f90)

ADD_LIBRARY(NTPoly ${ntpoly_src})

TARGET_INCLUDE_DIRECTORIES(NTPoly PRIVATE ${INC_PATHS})

INSTALL(TARGETS NTPoly
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
