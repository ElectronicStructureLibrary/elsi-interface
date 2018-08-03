MACRO(cpp_scotch arg1)
  TARGET_COMPILE_DEFINITIONS(${arg1} PRIVATE SCOTCH_VERSION=6)
  TARGET_COMPILE_DEFINITIONS(${arg1} PRIVATE SCOTCH_RELEASE=0)
  TARGET_COMPILE_DEFINITIONS(${arg1} PRIVATE SCOTCH_PATCHLEVEL=5)
  TARGET_COMPILE_DEFINITIONS(${arg1} PRIVATE COMMON_RANDOM_FIXED_SEED)
  TARGET_COMPILE_DEFINITIONS(${arg1} PRIVATE SCOTCH_RENAME)
  TARGET_COMPILE_DEFINITIONS(${arg1} PRIVATE IDXSIZE32)
  TARGET_COMPILE_DEFINITIONS(${arg1} PRIVATE INTSIZE32)
ENDMACRO()

MACRO(inc_scotch arg1)
  TARGET_INCLUDE_DIRECTORIES(${arg1} PRIVATE
    src/libscotch
    src/libscotchmetis
    ${PROJECT_BINARY_DIR}/generated
    ${INC_PATHS})
ENDMACRO()

SET(SCOTCH_SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/SCOTCH/src")

### Generate headers ###
LIST(APPEND dsize_src src/libscotch/dummysizes.c)
LIST(APPEND pdsize_src src/libscotch/ptdummysizes.c)

ADD_EXECUTABLE(dummysizes ${dsize_src})
ADD_EXECUTABLE(ptdummysizes ${pdsize_src})

cpp_scotch(dummysizes)
cpp_scotch(ptdummysizes)
TARGET_COMPILE_DEFINITIONS(ptdummysizes PRIVATE SCOTCH_PTSCOTCH)

inc_scotch(dummysizes)
inc_scotch(ptdummysizes)

IF(NOT SCOTCH_LAST_RESORT STREQUAL "NOT_USED")
  SEPARATE_ARGUMENTS(SCOTCH_LAST_RESORT UNIX_COMMAND "${SCOTCH_LAST_RESORT}")

  ADD_CUSTOM_COMMAND(OUTPUT ${PROJECT_BINARY_DIR}/generated/scotch.h ${PROJECT_BINARY_DIR}/generated/scotchf.h
    COMMAND ${SCOTCH_LAST_RESORT} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dummysizes ${SCOTCH_SOURCE_DIR}/libscotch/library.h ${PROJECT_BINARY_DIR}/generated/scotch.h
    COMMAND ${SCOTCH_LAST_RESORT} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dummysizes ${SCOTCH_SOURCE_DIR}/libscotch/library_f.h ${PROJECT_BINARY_DIR}/generated/scotchf.h
    DEPENDS dummysizes
    VERBATIM)

  ADD_CUSTOM_COMMAND(OUTPUT ${PROJECT_BINARY_DIR}/generated/ptscotch.h ${PROJECT_BINARY_DIR}/generated/ptscotchf.h
    COMMAND ${SCOTCH_LAST_RESORT} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ptdummysizes ${SCOTCH_SOURCE_DIR}/libscotch/library_pt.h ${PROJECT_BINARY_DIR}/generated/ptscotch.h
    COMMAND ${SCOTCH_LAST_RESORT} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ptdummysizes ${SCOTCH_SOURCE_DIR}/libscotch/library_pt_f.h ${PROJECT_BINARY_DIR}/generated/ptscotchf.h
    DEPENDS ptdummysizes
    VERBATIM)
ELSE()
  ADD_CUSTOM_COMMAND(OUTPUT ${PROJECT_BINARY_DIR}/generated/scotch.h ${PROJECT_BINARY_DIR}/generated/scotchf.h
    COMMAND dummysizes ${SCOTCH_SOURCE_DIR}/libscotch/library.h ${PROJECT_BINARY_DIR}/generated/scotch.h
    COMMAND dummysizes ${SCOTCH_SOURCE_DIR}/libscotch/library_f.h ${PROJECT_BINARY_DIR}/generated/scotchf.h
    DEPENDS dummysizes
    VERBATIM)

  ADD_CUSTOM_COMMAND(OUTPUT ${PROJECT_BINARY_DIR}/generated/ptscotch.h ${PROJECT_BINARY_DIR}/generated/ptscotchf.h
    COMMAND ptdummysizes ${SCOTCH_SOURCE_DIR}/libscotch/library_pt.h ${PROJECT_BINARY_DIR}/generated/ptscotch.h
    COMMAND ptdummysizes ${SCOTCH_SOURCE_DIR}/libscotch/library_pt_f.h ${PROJECT_BINARY_DIR}/generated/ptscotchf.h
    DEPENDS ptdummysizes
    VERBATIM)
ENDIF()

ADD_CUSTOM_COMMAND(OUTPUT ${PROJECT_BINARY_DIR}/generated/parser_ll.c ${PROJECT_BINARY_DIR}/generated/parser_yy.c
  COMMAND bison -pscotchyy -y -b y -d -v -o${PROJECT_BINARY_DIR}/generated/tmp2.c ${SCOTCH_SOURCE_DIR}/libscotch/parser_yy.y
  COMMAND mv ${PROJECT_BINARY_DIR}/generated/tmp2.c ${PROJECT_BINARY_DIR}/generated/parser_yy.c
  COMMAND mv ${PROJECT_BINARY_DIR}/generated/tmp2.h ${PROJECT_BINARY_DIR}/generated/parser_ly.h
  COMMAND flex -Pscotchyy -o${PROJECT_BINARY_DIR}/generated/tmp1.c ${SCOTCH_SOURCE_DIR}/libscotch/parser_ll.l
  COMMAND mv ${PROJECT_BINARY_DIR}/generated/tmp1.c ${PROJECT_BINARY_DIR}/generated/parser_ll.c
  VERBATIM)

### Source files ###
LIST(APPEND scotch_src
  ${PROJECT_BINARY_DIR}/generated/scotch.h
  ${PROJECT_BINARY_DIR}/generated/scotchf.h
  ${PROJECT_BINARY_DIR}/generated/parser_ll.c
  ${PROJECT_BINARY_DIR}/generated/parser_yy.c
  src/libscotch/arch.c
  src/libscotch/arch_build.c
  src/libscotch/arch_cmplt.c
  src/libscotch/arch_cmpltw.c
  src/libscotch/arch_deco.c
  src/libscotch/arch_dist.c
  src/libscotch/arch_hcub.c
  src/libscotch/arch_mesh.c
  src/libscotch/arch_tleaf.c
  src/libscotch/arch_torus.c
  src/libscotch/arch_vcmplt.c
  src/libscotch/arch_vhcub.c
  src/libscotch/bgraph.c
  src/libscotch/bgraph_bipart_bd.c
  src/libscotch/bgraph_bipart_df.c
  src/libscotch/bgraph_bipart_ex.c
  src/libscotch/bgraph_bipart_fm.c
  src/libscotch/bgraph_bipart_gg.c
  src/libscotch/bgraph_bipart_gp.c
  src/libscotch/bgraph_bipart_ml.c
  src/libscotch/bgraph_bipart_st.c
  src/libscotch/bgraph_bipart_zr.c
  src/libscotch/bgraph_check.c
  src/libscotch/bgraph_store.c
  src/libscotch/common.c
  src/libscotch/common_file.c
  src/libscotch/common_file_compress.c
  src/libscotch/common_file_uncompress.c
  src/libscotch/common_integer.c
  src/libscotch/common_memory.c
  src/libscotch/common_string.c
  src/libscotch/common_stub.c
  src/libscotch/common_thread.c
  src/libscotch/fibo.c
  src/libscotch/gain.c
  src/libscotch/geom.c
  src/libscotch/graph.c
  src/libscotch/graph_base.c
  src/libscotch/graph_band.c
  src/libscotch/graph_check.c
  src/libscotch/graph_coarsen.c
  src/libscotch/graph_induce.c
  src/libscotch/graph_io.c
  src/libscotch/graph_io_chac.c
  src/libscotch/graph_io_habo.c
  src/libscotch/graph_io_mmkt.c
  src/libscotch/graph_io_scot.c
  src/libscotch/graph_list.c
  src/libscotch/graph_match.c
  src/libscotch/hall_order_hd.c
  src/libscotch/hall_order_hf.c
  src/libscotch/hall_order_hx.c
  src/libscotch/hgraph.c
  src/libscotch/hgraph_check.c
  src/libscotch/hgraph_induce.c
  src/libscotch/hgraph_order_bl.c
  src/libscotch/hgraph_order_cp.c
  src/libscotch/hgraph_order_gp.c
  src/libscotch/hgraph_order_hd.c
  src/libscotch/hgraph_order_hf.c
  src/libscotch/hgraph_order_hx.c
  src/libscotch/hgraph_order_kp.c
  src/libscotch/hgraph_order_nd.c
  src/libscotch/hgraph_order_si.c
  src/libscotch/hgraph_order_st.c
  src/libscotch/hmesh.c
  src/libscotch/hmesh_check.c
  src/libscotch/hmesh_hgraph.c
  src/libscotch/hmesh_induce.c
  src/libscotch/hmesh_mesh.c
  src/libscotch/hmesh_order_bl.c
  src/libscotch/hmesh_order_cp.c
  src/libscotch/hmesh_order_gr.c
  src/libscotch/hmesh_order_gp.c
  src/libscotch/hmesh_order_hd.c
  src/libscotch/hmesh_order_hf.c
  src/libscotch/hmesh_order_hx.c
  src/libscotch/hmesh_order_nd.c
  src/libscotch/hmesh_order_si.c
  src/libscotch/hmesh_order_st.c
  src/libscotch/kgraph.c
  src/libscotch/kgraph_band.c
  src/libscotch/kgraph_check.c
  src/libscotch/kgraph_map_bd.c
  src/libscotch/kgraph_map_cp.c
  src/libscotch/kgraph_map_df.c
  src/libscotch/kgraph_map_ex.c
  src/libscotch/kgraph_map_fm.c
  src/libscotch/kgraph_map_ml.c
  src/libscotch/kgraph_map_rb.c
  src/libscotch/kgraph_map_rb_map.c
  src/libscotch/kgraph_map_rb_part.c
  src/libscotch/kgraph_map_st.c
  src/libscotch/kgraph_store.c
  src/libscotch/library_arch.c
  src/libscotch/library_arch_f.c
  src/libscotch/library_arch_build.c
  src/libscotch/library_arch_build_f.c
  src/libscotch/library_common_f.c
  src/libscotch/library_geom.c
  src/libscotch/library_geom_f.c
  src/libscotch/library_graph.c
  src/libscotch/library_graph_f.c
  src/libscotch/library_graph_base.c
  src/libscotch/library_graph_base_f.c
  src/libscotch/library_graph_check.c
  src/libscotch/library_graph_check_f.c
  src/libscotch/library_graph_coarsen.c
  src/libscotch/library_graph_coarsen_f.c
  src/libscotch/library_graph_color.c
  src/libscotch/library_graph_color_f.c
  src/libscotch/library_graph_io_chac.c
  src/libscotch/library_graph_io_chac_f.c
  src/libscotch/library_graph_io_habo.c
  src/libscotch/library_graph_io_habo_f.c
  src/libscotch/library_graph_io_mmkt.c
  src/libscotch/library_graph_io_mmkt_f.c
  src/libscotch/library_graph_io_scot.c
  src/libscotch/library_graph_io_scot_f.c
  src/libscotch/library_graph_map.c
  src/libscotch/library_graph_map_f.c
  src/libscotch/library_graph_map_io.c
  src/libscotch/library_graph_map_io_f.c
  src/libscotch/library_graph_map_view.c
  src/libscotch/library_graph_map_view_f.c
  src/libscotch/library_graph_order.c
  src/libscotch/library_graph_order_f.c
  src/libscotch/library_graph_part_ovl.c
  src/libscotch/library_graph_part_ovl_f.c
  src/libscotch/library_mapping.c
  src/libscotch/library_mesh.c
  src/libscotch/library_mesh_f.c
  src/libscotch/library_mesh_graph.c
  src/libscotch/library_mesh_graph_f.c
  src/libscotch/library_mesh_io_habo.c
  src/libscotch/library_mesh_io_habo_f.c
  src/libscotch/library_mesh_io_scot.c
  src/libscotch/library_mesh_io_scot_f.c
  src/libscotch/library_mesh_order.c
  src/libscotch/library_mesh_order_f.c
  src/libscotch/library_order.c
  src/libscotch/library_parser.c
  src/libscotch/library_parser_f.c
  src/libscotch/library_random.c
  src/libscotch/library_random_f.c
  src/libscotch/library_strat.c
  src/libscotch/library_version.c
  src/libscotch/library_version_f.c
  src/libscotch/mapping.c
  src/libscotch/mapping_io.c
  src/libscotch/mesh.c
  src/libscotch/mesh_check.c
  src/libscotch/mesh_coarsen.c
  src/libscotch/mesh_graph.c
  src/libscotch/mesh_induce_sepa.c
  src/libscotch/mesh_io.c
  src/libscotch/mesh_io_habo.c
  src/libscotch/mesh_io_scot.c
  src/libscotch/order.c
  src/libscotch/order_check.c
  src/libscotch/order_io.c
  src/libscotch/parser.c
  src/libscotch/vgraph.c
  src/libscotch/vgraph_check.c
  src/libscotch/vgraph_separate_bd.c
  src/libscotch/vgraph_separate_df.c
  src/libscotch/vgraph_separate_es.c
  src/libscotch/vgraph_separate_fm.c
  src/libscotch/vgraph_separate_gg.c
  src/libscotch/vgraph_separate_gp.c
  src/libscotch/vgraph_separate_ml.c
  src/libscotch/vgraph_separate_st.c
  src/libscotch/vgraph_separate_th.c
  src/libscotch/vgraph_separate_vw.c
  src/libscotch/vgraph_separate_zr.c
  src/libscotch/vgraph_store.c
  src/libscotch/vmesh.c
  src/libscotch/vmesh_check.c
  src/libscotch/vmesh_separate_fm.c
  src/libscotch/vmesh_separate_gg.c
  src/libscotch/vmesh_separate_gr.c
  src/libscotch/vmesh_separate_ml.c
  src/libscotch/vmesh_separate_zr.c
  src/libscotch/vmesh_separate_st.c
  src/libscotch/vmesh_store.c
  src/libscotch/wgraph.c
  src/libscotch/wgraph_check.c
  src/libscotch/wgraph_part_fm.c
  src/libscotch/wgraph_part_gg.c
  src/libscotch/wgraph_part_gp.c
  src/libscotch/wgraph_part_ml.c
  src/libscotch/wgraph_part_rb.c
  src/libscotch/wgraph_part_st.c
  src/libscotch/wgraph_part_zr.c
  src/libscotch/wgraph_store.c)

LIST(APPEND ptscotch_src
  ${PROJECT_BINARY_DIR}/generated/ptscotch.h
  ${PROJECT_BINARY_DIR}/generated/ptscotchf.h
  src/libscotch/bdgraph.c
  src/libscotch/bdgraph_bipart_bd.c
  src/libscotch/bdgraph_bipart_df.c
  src/libscotch/bdgraph_bipart_ex.c
  src/libscotch/bdgraph_bipart_ml.c
  src/libscotch/bdgraph_bipart_sq.c
  src/libscotch/bdgraph_bipart_st.c
  src/libscotch/bdgraph_bipart_zr.c
  src/libscotch/bdgraph_check.c
  src/libscotch/bdgraph_gather_all.c
  src/libscotch/bdgraph_store.c
  src/libscotch/comm.c
  src/libscotch/dgraph.c
  src/libscotch/dgraph_allreduce.c
  src/libscotch/dgraph_band.c
  src/libscotch/dgraph_build.c
  src/libscotch/dgraph_build_grid3d.c
  src/libscotch/dgraph_build_hcub.c
  src/libscotch/dgraph_check.c
  src/libscotch/dgraph_coarsen.c
  src/libscotch/dgraph_fold.c
  src/libscotch/dgraph_fold_comm.c
  src/libscotch/dgraph_fold_dup.c
  src/libscotch/dgraph_gather.c
  src/libscotch/dgraph_gather_all.c
  src/libscotch/dgraph_ghst.c
  src/libscotch/dgraph_halo.c
  src/libscotch/dgraph_induce.c
  src/libscotch/dgraph_io_load.c
  src/libscotch/dgraph_io_save.c
  src/libscotch/dgraph_match.c
  src/libscotch/dgraph_match_sync_coll.c
  src/libscotch/dgraph_match_sync_ptop.c
  src/libscotch/dgraph_match_check.c
  src/libscotch/dgraph_redist.c
  src/libscotch/dgraph_scatter.c
  src/libscotch/dgraph_view.c
  src/libscotch/dmapping.c
  src/libscotch/dmapping_io.c
  src/libscotch/dorder.c
  src/libscotch/dorder_gather.c
  src/libscotch/dorder_io.c
  src/libscotch/dorder_io_block.c
  src/libscotch/dorder_io_tree.c
  src/libscotch/dorder_perm.c
  src/libscotch/dorder_tree_dist.c
  src/libscotch/hdgraph.c
  src/libscotch/hdgraph_check.c
  src/libscotch/hdgraph_fold.c
  src/libscotch/hdgraph_gather.c
  src/libscotch/hdgraph_induce.c
  src/libscotch/hdgraph_order_nd.c
  src/libscotch/hdgraph_order_si.c
  src/libscotch/hdgraph_order_sq.c
  src/libscotch/hdgraph_order_st.c
  src/libscotch/kdgraph.c
  src/libscotch/kdgraph_gather.c
  src/libscotch/kdgraph_map_rb.c
  src/libscotch/kdgraph_map_rb_map.c
  src/libscotch/kdgraph_map_rb_part.c
  src/libscotch/kdgraph_map_st.c
  src/libscotch/library_dgraph.c
  src/libscotch/library_dgraph_f.c
  src/libscotch/library_dgraph_band.c
  src/libscotch/library_dgraph_band_f.c
  src/libscotch/library_dgraph_build.c
  src/libscotch/library_dgraph_build_f.c
  src/libscotch/library_dgraph_build_grid3d.c
  src/libscotch/library_dgraph_build_grid3d_f.c
  src/libscotch/library_dgraph_check.c
  src/libscotch/library_dgraph_check_f.c
  src/libscotch/library_dgraph_coarsen.c
  src/libscotch/library_dgraph_coarsen_f.c
  src/libscotch/library_dgraph_gather.c
  src/libscotch/library_dgraph_gather_f.c
  src/libscotch/library_dgraph_grow.c
  src/libscotch/library_dgraph_halo.c
  src/libscotch/library_dgraph_halo_f.c
  src/libscotch/library_dgraph_induce.c
  src/libscotch/library_dgraph_induce_f.c
  src/libscotch/library_dgraph_io_load.c
  src/libscotch/library_dgraph_io_load_f.c
  src/libscotch/library_dgraph_io_save.c
  src/libscotch/library_dgraph_io_save_f.c
  src/libscotch/library_dgraph_map.c
  src/libscotch/library_dgraph_map_f.c
  src/libscotch/library_dgraph_map_view.c
  src/libscotch/library_dgraph_map_view_f.c
  src/libscotch/library_dgraph_order.c
  src/libscotch/library_dgraph_order_f.c
  src/libscotch/library_dgraph_order_gather.c
  src/libscotch/library_dgraph_order_gather_f.c
  src/libscotch/library_dgraph_order_io.c
  src/libscotch/library_dgraph_order_io_f.c
  src/libscotch/library_dgraph_order_io_block.c
  src/libscotch/library_dgraph_order_io_block_f.c
  src/libscotch/library_dgraph_order_perm.c
  src/libscotch/library_dgraph_order_perm_f.c
  src/libscotch/library_dgraph_order_tree_dist.c
  src/libscotch/library_dgraph_order_tree_dist_f.c
  src/libscotch/library_dgraph_redist.c
  src/libscotch/library_dgraph_redist_f.c
  src/libscotch/library_dgraph_scatter.c
  src/libscotch/library_dgraph_scatter_f.c
  src/libscotch/library_dgraph_stat.c
  src/libscotch/library_dgraph_stat_f.c
  src/libscotch/library_dmapping.c
  src/libscotch/library_dorder.c
  src/libscotch/vdgraph.c
  src/libscotch/vdgraph_check.c
  src/libscotch/vdgraph_gather_all.c
  src/libscotch/vdgraph_separate_bd.c
  src/libscotch/vdgraph_separate_df.c
  src/libscotch/vdgraph_separate_ml.c
  src/libscotch/vdgraph_separate_sq.c
  src/libscotch/vdgraph_separate_st.c
  src/libscotch/vdgraph_separate_zr.c
  src/libscotch/vdgraph_store.c)

SET_SOURCE_FILES_PROPERTIES(${PROJECT_BINARY_DIR}/generated/scotch.h
  PROPERTIES GENERATED TRUE)

SET_SOURCE_FILES_PROPERTIES(${PROJECT_BINARY_DIR}/generated/scotchf.h
  PROPERTIES GENERATED TRUE)

SET_SOURCE_FILES_PROPERTIES(${PROJECT_BINARY_DIR}/generated/ptscotch.h
  PROPERTIES GENERATED TRUE)

SET_SOURCE_FILES_PROPERTIES(${PROJECT_BINARY_DIR}/generated/ptscotchf.h
  PROPERTIES GENERATED TRUE)

SET_SOURCE_FILES_PROPERTIES(${PROJECT_BINARY_DIR}/generated/parser_ll.c
  PROPERTIES GENERATED TRUE)

SET_SOURCE_FILES_PROPERTIES(${PROJECT_BINARY_DIR}/generated/parser_yy.c
  PROPERTIES GENERATED TRUE)

SET_SOURCE_FILES_PROPERTIES(src/libscotch/common.c
  PROPERTIES COMPILE_DEFINITIONS SCOTCH_COMMON_RENAME)

SET_SOURCE_FILES_PROPERTIES(src/libscotch/common_file.c
  PROPERTIES COMPILE_DEFINITIONS SCOTCH_COMMON_RENAME)

SET_SOURCE_FILES_PROPERTIES(src/libscotch/common_file_compress.c
  PROPERTIES COMPILE_DEFINITIONS SCOTCH_COMMON_RENAME)

SET_SOURCE_FILES_PROPERTIES(src/libscotch/common_file_uncompress.c
  PROPERTIES COMPILE_DEFINITIONS SCOTCH_COMMON_RENAME)

SET_SOURCE_FILES_PROPERTIES(src/libscotch/common_integer.c
  PROPERTIES COMPILE_DEFINITIONS SCOTCH_COMMON_RENAME)

SET_SOURCE_FILES_PROPERTIES(src/libscotch/common_memory.c
  PROPERTIES COMPILE_DEFINITIONS SCOTCH_COMMON_RENAME)

SET_SOURCE_FILES_PROPERTIES(src/libscotch/common_string.c
  PROPERTIES COMPILE_DEFINITIONS SCOTCH_COMMON_RENAME)

SET_SOURCE_FILES_PROPERTIES(src/libscotch/common_stub.c
  PROPERTIES COMPILE_DEFINITIONS SCOTCH_COMMON_RENAME)

SET_SOURCE_FILES_PROPERTIES(src/libscotch/common_thread.c
  PROPERTIES COMPILE_DEFINITIONS SCOTCH_COMMON_RENAME)

LIST(APPEND scotcherr_src src/libscotch/library_error.c)

LIST(APPEND ptscotcherr_src src/libscotch/library_error_pt.c)

LIST(APPEND scotchmetis_src
  src/libscotchmetis/metis_graph_order.c
  src/libscotchmetis/metis_graph_order_f.c
  src/libscotchmetis/metis_graph_part.c
  src/libscotchmetis/metis_graph_part_f.c)

LIST(APPEND ptscotchparmetis_src
  src/libscotchmetis/parmetis_dgraph_order.c
  src/libscotchmetis/parmetis_dgraph_order_f.c
  src/libscotchmetis/parmetis_dgraph_part.c
  src/libscotchmetis/parmetis_dgraph_part_f.c)

ADD_LIBRARY(scotch ${scotch_src})
ADD_LIBRARY(scotcherr ${scotcherr_src})
ADD_LIBRARY(scotchmetis ${scotchmetis_src})
ADD_LIBRARY(ptscotch ${ptscotch_src})
ADD_LIBRARY(ptscotcherr ${ptscotcherr_src})
ADD_LIBRARY(ptscotchparmetis ${ptscotchparmetis_src})

ADD_DEPENDENCIES(scotcherr scotch)
ADD_DEPENDENCIES(scotchmetis scotch)
ADD_DEPENDENCIES(ptscotch scotch)
ADD_DEPENDENCIES(ptscotcherr ptscotch)
ADD_DEPENDENCIES(ptscotchparmetis ptscotch)

cpp_scotch(scotch)
cpp_scotch(scotcherr)
cpp_scotch(ptscotch)
cpp_scotch(ptscotcherr)
TARGET_COMPILE_DEFINITIONS(ptscotch PRIVATE SCOTCH_PTSCOTCH)
TARGET_COMPILE_DEFINITIONS(ptscotcherr PRIVATE SCOTCH_PTSCOTCH)
TARGET_COMPILE_DEFINITIONS(scotchmetis PRIVATE SCOTCHLIB=ptscotch)
TARGET_COMPILE_DEFINITIONS(ptscotchparmetis PRIVATE SCOTCH_PTSCOTCH)
TARGET_COMPILE_DEFINITIONS(ptscotchparmetis PRIVATE SCOTCHLIB=ptscotch)

IF(NOT CMAKE_VERSION VERSION_LESS 3.8.2)
  TARGET_COMPILE_FEATURES(scotch PRIVATE c_std_99)
  TARGET_COMPILE_FEATURES(scotcherr PRIVATE c_std_99)
  TARGET_COMPILE_FEATURES(scotchmetis PRIVATE c_std_99)
  TARGET_COMPILE_FEATURES(ptscotch PRIVATE c_std_99)
  TARGET_COMPILE_FEATURES(ptscotcherr PRIVATE c_std_99)
  TARGET_COMPILE_FEATURES(ptscotchparmetis PRIVATE c_std_99)
ENDIF()

inc_scotch(scotch)
inc_scotch(scotcherr)
inc_scotch(scotchmetis)
inc_scotch(ptscotch)
inc_scotch(ptscotcherr)
inc_scotch(ptscotchparmetis)

INSTALL(TARGETS scotch
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})

INSTALL(TARGETS scotcherr
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})

INSTALL(TARGETS scotchmetis
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})

INSTALL(TARGETS ptscotch
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})

INSTALL(TARGETS ptscotcherr
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})

INSTALL(TARGETS ptscotchparmetis
  EXPORT elsiConfig
  DESTINATION ${CMAKE_INSTALL_LIBDIR})
