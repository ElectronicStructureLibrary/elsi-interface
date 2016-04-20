
lib_LIBRARIES = libpexsi.a 

nobase_include_HEADERS= f_ppexsi_interface.mod

LDADD = $(LIBSPARMETISLIB) $(SUPERLU_LIBS) $(SCALAPACK_FCFLAGS)


#thisincludedir = $(top_builddir)/external/PEXSI/external/include

#AM_CFLAGS = $(LIBOMM) -I$(LIBMATSWITCH)/. -I$(thisincludedir)
#AM_CFLAGS   = -I$(thisincludedir)
AM_CFLAGS   = -I$(top_builddir)/external/PEXSI/include -I$(top_builddir)/external/include
AM_CXXFLAGS = -I$(top_builddir)/external/PEXSI/include -I$(top_builddir)/external/include

libpexsi_a_SOURCES = pdsymbfact.c pzsymbfact.c get_perm_c_parmetis.c \
                     interface.cpp ppexsi.cpp pole.cpp TreeBcast.cpp \
		     superlu_dist_internal_complex.cpp superlu_dist_internal_real.cpp \
		     mpi_interf.cpp lapack.cpp blas.cpp utility.cpp global.cpp timer.cpp \
                     f_interface.f90 
#interface_old.cpp

