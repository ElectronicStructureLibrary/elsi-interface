set(UPCXX_PREFIX $ENV{UPCXX_DIR})
if(UPCXX_PREFIX)
  message("UPCXX has been provided, not compiling local version.")

  set(UPCXX_INCLUDE_DIR ${UPCXX_PREFIX}/include ${GASNET_INCLUDE_DIR} ${GASNET_CONDUIT_INCLUDE_DIR})
  set(UPCXX_LIBRARY_PATH ${UPCXX_PREFIX}/lib)

  add_library(libupcxx STATIC IMPORTED)

  set_property(TARGET libupcxx PROPERTY IMPORTED_LOCATION ${UPCXX_PREFIX}/lib/libupcxx.a)

  set(UPCXX_LIBRARIES ${GASNET_LIBRARIES})
  set(UPCXX_LIBRARY_PATH ${UPCXX_LIBRARY_PATH} ${GASNET_LIBRARY_PATH})
  set(UPCXX_DEFINES "-DGASNET_ALLOW_OPTIMIZED_DEBUG")

  add_library(upcxx STATIC IMPORTED)
else()
  set(UPCXX_NAME upcxx)
  set(UPCXX_REPO https://bitbucket.org/upcxx/upcxx.git)

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(UPCXX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
    ExternalProject_Add(${UPCXX_NAME}
      DEPENDS ${GASNET_NAME}
      GIT_REPOSITORY ${UPCXX_REPO}
      UPDATE_COMMAND ""
      INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/upcxx_install
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --with-gasnet=${GASNET_CONDUIT} CC=${MPI_C_COMPILER} CXX=${MPI_CXX_COMPILER} CFLAGS=${CMAKE_C_FLAGS} CXXFLAGS=${UPCXX_CXX_FLAGS}
      )
  else()
    ExternalProject_Add(${UPCXX_NAME}
      DEPENDS ${GASNET_NAME}
      GIT_REPOSITORY ${UPCXX_REPO}
      UPDATE_COMMAND ""
      INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/upcxx_install
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --with-gasnet=${GASNET_CONDUIT} CC=${MPI_C_COMPILER} CXX=${MPI_CXX_COMPILER} CFLAGS=${CMAKE_C_FLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}
      )
  endif()

  ExternalProject_Add_Step(${UPCXX_NAME} bootstrap
    DEPENDEES patch update patch download
    DEPENDERS configure
    COMMAND libtoolize COMMAND autoreconf -fi COMMAND <SOURCE_DIR>/Bootstrap.sh
    WORKING_DIRECTORY <SOURCE_DIR>
    COMMENT "Bootstraping the source directory"
    )

  set(UPCXX_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/upcxx_install/include)
  set(UPCXX_LIBRARY_PATH ${CMAKE_CURRENT_BINARY_DIR}/external/upcxx_install/lib)
  set(UPCXX_LIBRARIES ${GASNET_LIBRARIES})
  set(UPCXX_INCLUDE_DIR ${UPCXX_INCLUDE_DIR} ${GASNET_INCLUDE_DIR} ${GASNET_CONDUIT_INCLUDE_DIR})
  set(UPCXX_LIBRARY_PATH ${UPCXX_LIBRARY_PATH} ${GASNET_LIBRARY_PATH})
  set(UPCXX_DEFINES)

  ExternalProject_Get_Property(${UPCXX_NAME} install_dir)
  add_library(libupcxx STATIC IMPORTED)
  set_property(TARGET libupcxx PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libupcxx.a)
  add_dependencies(libupcxx ${UPCXX_NAME} ${GASNET_NAME})
endif()
