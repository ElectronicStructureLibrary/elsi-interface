INCLUDE_DIRECTORIES(src)

### Header files ###
LIST(APPEND superlu_h
  src/Cnames.h
  src/colamd.h
  src/cublas_utils.h
  src/dcomplex.h
  src/machines.h
  src/psymbfact.h
  src/superlu_ddefs.h
  src/superlu_defs.h
  src/superlu_dist_config.h
  src/superlu_enum_consts.h
  src/superlu_zdefs.h
  src/supermatrix.h
  src/util_dist.h
)

### Source files ###
LIST(APPEND superlu_src
  src/colamd.c
  src/comm.c
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
  src/pdgssvx.c
  src/pdgstrf.c
  src/pdgstrf2.c
  src/pdgstrs.c
  src/pdgstrs1.c
  src/pdgstrs_Bglobal.c
  src/pdgstrs_lsum.c
  src/pdgsrfs_ABXglobal.c
  src/pdgssvx_ABglobal.c
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
  src/pzgssvx.c
  src/pzgssvx_ABglobal.c
  src/pzgstrf.c
  src/pzgstrf2.c
  src/pzgstrs.c
  src/pzgstrs1.c
  src/pzgstrs_lsum.c
  src/pzgstrs_Bglobal.c
  src/pzgsrfs_ABXglobal.c
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
  src/util.c
  src/xerr_dist.c
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
  src/zsp_blas2_dist.c
  src/zsp_blas3_dist.c
  src/zutil_dist.c
)

ADD_LIBRARY(superlu_dist ${superlu_src} ${superlu_h})

TARGET_COMPILE_DEFINITIONS(superlu_dist PRIVATE NDEBUG)
TARGET_COMPILE_DEFINITIONS(superlu_dist PRIVATE USE_VENDOR_BLAS)
TARGET_COMPILE_DEFINITIONS(superlu_dist PRIVATE Add_)
TARGET_COMPILE_DEFINITIONS(superlu_dist PRIVATE PRNTlevel=0)
TARGET_COMPILE_DEFINITIONS(superlu_dist PRIVATE DEBUGlevel=0)
TARGET_COMPILE_FEATURES(superlu_dist PUBLIC c_std_99)

INSTALL(TARGETS superlu_dist
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
)
