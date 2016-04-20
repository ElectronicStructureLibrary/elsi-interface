
LIBPEXSI=$(top_builddir)/external/PEXSI/src/libpexsi.a


LDADD = $(LIBPEXSI) $(LIBSPARMETISLIB) $(SUPERLU_LIBS)


thisincludedir = $(top_builddir)/external/PEXSI/external/include

#AM_CFLAGS = $(LIBOMM) -I$(LIBMATSWITCH)/. -I$(thisincludedir)
#AM_CFLAGS   = -I$(thisincludedir)
AM_CFLAGS   =  -I$(thisincludedir) -I$(top_builddir)/external/PEXSI/include -I$(top_builddir)/external/include
AM_CXXFLAGS =  -I$(thisincludedir) -I$(top_builddir)/external/PEXSI/include -I$(top_builddir)/external/include

noinst_PROGRAMS = driver_ksdft driver_ksdft_2 driver_pselinv_real driver_pselinv_complex \
                  driver_pselinv_complex_unsym driver_pselinv_real_unsym \
		  run_ksdft run_inertia run_pselinv run_pselinv_unsym run_test_reduce 

driver_ksdft_SOURCE = driver_ksdft.c

driver_ksdft2_SOURCE = driver_ksdft_2.c

driver_pselinv_real_SOURCE = driver_pselinv_real.c

driver_pselinv_complex_SOURCE = driver_pselinv_complex.c

driver_pselinv_complex_unsym_SOURCE = driver_pselinv_complex_unsym.c

driver_pselinv_real_unsym_SOURCE = driver_pselinv_real_unsym.c

run_ksdf_SOURCE = run_ksdft.cpp

run_inertia_SOURCE = run_inertia.cpp

run_pselinv_SOURCE = run_pselinv.cpp

run_pselinv_unsym_SOURCE = run_pselinv_unsym.cpp

run_test_reduce_SOURCE = run_test_reduce.cpp


