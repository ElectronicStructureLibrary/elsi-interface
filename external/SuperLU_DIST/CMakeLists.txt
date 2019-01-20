### Source files ###
LIST(APPEND superlu_src
  src/colamd.c
  src/comm.c
  src/dbinary_io.c
  src/dcomplex_dist.c
  src/ddistribute.c
  src/dgsequ_dist.c
  src/dlangs_dist.c
  src/dlaqgs_dist.c
  src/dldperm_dist.c
  src/dmach_dist.c
  src/dmemory_dist.c
  src/dmyblas2_dist.c
  src/dreadMM.c
  src/dreadhb.c
  src/dreadrb.c
  src/dreadtriple.c
  src/dreadtriple_noheader.c
  src/dsp_blas2_dist.c
  src/dsp_blas3_dist.c
  src/dutil_dist.c
  src/etree.c 
  src/get_perm_c.c
  src/get_perm_c_parmetis.c
  src/mc64ad_dist.c
  src/memory.c
  src/mmd.c
  src/pdGetDiagU.c
  src/pddistribute.c
  src/pdgsequ.c
  src/pdgsmv.c
  src/pdgsmv_AXglobal.c
  src/pdgsrfs.c
  src/pdgsrfs_ABXglobal.c
  src/pdgssvx.c
  src/pdgssvx_ABglobal.c
  src/pdgstrf.c
  src/pdgstrf2.c
  src/pdgstrs.c
  src/pdgstrs1.c
  src/pdgstrs_Bglobal.c
  src/pdgstrs_lsum.c
  src/pdlangs.c
  src/pdlaqgs.c
  src/pdsymbfact_distdata.c
  src/pdutil.c
  src/psymbfact.c
  src/psymbfact_util.c
  src/pxerr_dist.c
  src/pzGetDiagU.c
  src/pzdistribute.c
  src/pzgsequ.c
  src/pzgsmv.c
  src/pzgsmv_AXglobal.c
  src/pzgsrfs.c
  src/pzgsrfs_ABXglobal.c
  src/pzgssvx.c
  src/pzgssvx_ABglobal.c
  src/pzgstrf.c
  src/pzgstrf2.c
  src/pzgstrs.c
  src/pzgstrs1.c
  src/pzgstrs_lsum.c
  src/pzgstrs_Bglobal.c
  src/pzlangs.c
  src/pzlaqgs.c
  src/pzsymbfact_distdata.c
  src/pzutil.c
  src/smach_dist.c
  src/sp_colorder.c
  src/sp_ienv.c
  src/static_schedule.c
  src/superlu_dist_version.c
  src/superlu_grid.c
  src/superlu_timer.c
  src/symbfact.c
  src/TreeInterface.cpp
  src/util.c
  src/xerr_dist.c
  src/zbinary_io.c
  src/zdistribute.c
  src/zgsequ_dist.c
  src/zlangs_dist.c
  src/zlaqgs_dist.c
  src/zldperm_dist.c
  src/zmemory_dist.c
  src/zmyblas2_dist.c
  src/zreadMM.c
  src/zreadhb.c
  src/zreadrb.c
  src/zreadtriple.c
  src/zreadtriple_noheader.c
  src/zsp_blas2_dist.c
  src/zsp_blas3_dist.c
  src/zutil_dist.c)

ADD_LIBRARY(superlu_dist ${superlu_src})

TARGET_COMPILE_DEFINITIONS(superlu_dist PRIVATE NDEBUG)
TARGET_COMPILE_DEFINITIONS(superlu_dist PRIVATE USE_VENDOR_BLAS)
TARGET_COMPILE_DEFINITIONS(superlu_dist PRIVATE PRNTlevel=0)
TARGET_COMPILE_DEFINITIONS(superlu_dist PRIVATE DEBUGlevel=0)
IF(ADD_UNDERSCORE)
  TARGET_COMPILE_DEFINITIONS(superlu_dist PRIVATE Add_)
ENDIF()

IF(NOT CMAKE_VERSION VERSION_LESS 3.8.2)
  TARGET_COMPILE_FEATURES(superlu_dist PRIVATE c_std_99)
ENDIF()

ADD_DEPENDENCIES(superlu_dist
  ptscotchparmetis
  ptscotch
  ptscotcherr
  scotchmetis
  scotch
  scotcherr)

TARGET_LINK_LIBRARIES(superlu_dist PRIVATE
  ptscotchparmetis
  ptscotch
  ptscotcherr
  scotchmetis
  scotch
  scotcherr
  ${LIBS})

TARGET_INCLUDE_DIRECTORIES(superlu_dist PRIVATE
  src
  ${INC_PATHS}
  ${PROJECT_SOURCE_DIR}/external/SCOTCH/src/libscotch
  ${PROJECT_SOURCE_DIR}/external/SCOTCH/src/libscotchmetis
  ${PROJECT_BINARY_DIR}/generated)

INSTALL(TARGETS superlu_dist
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
