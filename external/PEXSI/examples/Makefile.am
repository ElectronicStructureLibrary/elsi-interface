if !HAS_LIBPARMETIS
LIBPARMETIS=$(top_builddir)/external/lib/libparmetis.a
LIBMETIS=$(top_builddir)/external/lib/libmetis.a
else
LIBPARMETIS=$(LIBSPARMETISLIB)
LIBMETIS=$(LIBSPARMETISLIB)
endif

if !HAS_SUPERLU
LIBSUPERLU=$(top_builddir)/external/lib/libsuperlu.a
else
LIBSUPERLU=$(SUPERLU_LIBS)
endif

#this has to change
LIBPEXSI=$(top_builddir)/external/PEXSI/src/libpexsi.a

LDADD = $(LIBPEXSI) $(LIBPARMETIS) $(LIBMETIS) $(LIBSUPERLU) $(SCALAPACK_FCFLAGS)

thisincludedir = $(top_builddir)/external/PEXSI/external/include

#AM_CFLAGS = $(LIBOMM) -I$(LIBMATSWITCH)/. -I$(thisincludedir)
#AM_CFLAGS   = -I$(thisincludedir)
AM_FCFLAGS =  -I$(thisincludedir) -I$(top_builddir)/external/PEXSI/include -I$(top_builddir)/external/PEXSI/src

AM_CFLAGS   =  -I$(thisincludedir) -I$(top_builddir)/external/PEXSI/include -I$(top_builddir)/external/include -I$(top_builddir)/external/PEXSI/src

AM_CXXFLAGS = -I$(thisincludedir) -I$(top_builddir)/external/PEXSI/include -I$(top_builddir)/external/include -I$(top_builddir)/external/PEXSI/src

noinst_PROGRAMS = driver_ksdft driver_pselinv_real \
                  driver_pselinv_complex \
                  run_pselinv run_ppexsi_old 
# run_ppexsi  
#run_inertia f_ppexsi

driver_ksdft_SOURCES = driver_ksdft.c

driver_pselinv_real_SOURCES = driver_pselinv_real.c

driver_pselinv_complex_SOURCES = driver_pselinv_complex.c

#run_inertia_SOURCES = run_inertia.cpp

run_pselinv_SOURCES = run_pselinv.cpp

#run_ppexsi_SOURCES = run_ppexsi.cpp

run_ppexsi_old_SOURCES = run_ppexsi_old.cpp

#f_ppexsi_SOURCES = f_ppexsi.f90
