### SSEIG source files ###
LIST(APPEND sseig_src
  src/SSEIG/dlastd.f
  src/SSEIG/dssmm.f
  src/SSEIG/dssmv.f
  src/SSEIG/dssr2.f
  src/SSEIG/dssr2k.f
  src/SSEIG/dsstd2.f
  src/SSEIG/dsstrd.f
  src/SSEIG/dutil.f
  src/SSEIG/pdlastd.f
  src/SSEIG/pdssmv.f
  src/SSEIG/pdssr2k.f
  src/SSEIG/pdsstd2.f
  src/SSEIG/pdsstrd.f
  src/SSEIG/pdutil.f
  src/SSEIG/pslastd.f
  src/SSEIG/psssmv.f
  src/SSEIG/psssr2k.f
  src/SSEIG/pssstd2.f
  src/SSEIG/pssstrd.f
  src/SSEIG/psutil.f
  src/SSEIG/slastd.f
  src/SSEIG/sssmm.f
  src/SSEIG/sssmv.f
  src/SSEIG/sssr2.f
  src/SSEIG/sssr2k.f
  src/SSEIG/ssstd2.f
  src/SSEIG/ssstrd.f
  src/SSEIG/sutil.f)

### BSEPACK source files ###
LIST(APPEND bsepack_src
  src/BSEPACK/approx_delta.f
  src/BSEPACK/dbsabsp_aux.f
  src/BSEPACK/pdbsabsp.f
  src/BSEPACK/pdbseig.f
  src/BSEPACK/pdbslanczos.f
  src/BSEPACK/pdbssolver1.f
  src/BSEPACK/pdbssolver1_svd.f
  src/BSEPACK/pddot_all.f
  src/BSEPACK/pdembed1.f
  src/BSEPACK/pdembed2.f
  src/BSEPACK/pdlanczos.f
  src/BSEPACK/pdsorteig.f
  src/BSEPACK/pzbsabsp.f
  src/BSEPACK/pzbsaux1.f
  src/BSEPACK/pzbseig.f
  src/BSEPACK/pzbslanczos.f
  src/BSEPACK/pzbssolver1.f
  src/BSEPACK/pzdotc_all.f
  src/BSEPACK/pzembed1.f
  src/BSEPACK/pzembed2.f
  src/BSEPACK/pzlaconj.f
  src/BSEPACK/pzlanczos.f)

ADD_LIBRARY(sseig ${sseig_src})
ADD_LIBRARY(bsepack ${bsepack_src})

IF(ELSI_USE_FIND_PACKAGE_MPI)
  find_package(MPI)
  TARGET_LINK_LIBRARIES(sseig PUBLIC MPI::MPI_C MPI::MPI_Fortran)
  TARGET_LINK_LIBRARIES(bsepack PUBLIC MPI::MPI_C MPI::MPI_Fortran)
ENDIF()

ADD_DEPENDENCIES(bsepack sseig)

TARGET_LINK_LIBRARIES(bsepack PRIVATE sseig ${LIBS})
TARGET_LINK_LIBRARIES(sseig PRIVATE ${LIBS})

IF(CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
  TARGET_COMPILE_OPTIONS(sseig PRIVATE "-qfixed")
  TARGET_COMPILE_OPTIONS(bsepack PRIVATE "-qfixed")
ENDIF()

INSTALL(TARGETS sseig
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
INSTALL(TARGETS bsepack
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
