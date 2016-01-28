
#SUBDIRS = MatrixSwitch-0.1.2

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
LIBPEXSI=$(top_builddir)/external/PEXSI/src/libpexsi.a

#LIBS  = ${PEXSI_LIB} ${DSUPERLU_LIB} ${PAR_ND_LIB} ${SEQ_ND_LIB} ${SCALAPACK_LIB} ${BLACS_LIB} ${FORTRAN_LIB}


#LIBOMM = $(LIBOMM_FLAGS: -D%=-WF,-D%)


#lib_LIBRARIES = libpexsi.a 
bin_PROGRAMS = f_driver_ksdft f_driver_pselinv_real f_driver_pselinv_complex

LDADD = $(LIBPEXSI)  $(LIBPARMETIS) $(LIBMETIS) $(LIBSUPERLU)  $(SCALAPACK_FCFLAGS)


thisincludedir = $(top_builddir)/external/PEXSI/external/include

#AM_CFLAGS = $(LIBOMM) -I$(LIBMATSWITCH)/. -I$(thisincludedir)
AM_CFLAGS   = -I$(thisincludedir)
AM_FCFLAGS =  -I$(thisincludedir) -I$(top_builddir)/external/PEXSI/include -I$(top_builddir)/external/PEXSI/src -lstdc++
#AM_CXXFLAGS = -lstdc++ -std=c++11 -I$(thisincludedir) -I$(top_builddir)/external/PEXSI/include

#include_HEADERS = f_ppexsi_interface.mod


f_driver_ksdft_SOURCES = f_driver_ksdft.f90 #f_interface.f90

f_driver_pselinv_real_SOURCES = f_driver_pselinv_real.f90 # f_interface.f90

f_driver_pselinv_complex_SOURCES = f_driver_pselinv_complex.f90  #f_interface.f90



