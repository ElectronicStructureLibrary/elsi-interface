


LIBOMMLIB=$(top_builddir)/external/libOMM/src/libOMM.a

if HAVE_OPENMP
LIBELPA=$(top_builddir)/external/ELPA/.libs/libelpa_openmp.a
else
LIBELPA=$(top_builddir)/external/ELPA/.libs/libelpa.a
endif

LIBELSI=$(top_builddir)/ELSI/src/libelsi.a

#this has to change
LIBPEXSI=$(top_builddir)/external/PEXSI/src/libpexsi.a

#LIBOMM = $(LIBOMM_FLAGS: -D%=-WF,-D%)

MATRIXSWITCH = $(top_builddir)/external/libOMM/src/MatrixSwitch-0.1.2/src
LIBMATRIXSWITCH = $(MATRIXSWITCH)/libmatrixswitch.a
ELPAMODULES = $(top_builddir)/external/ELPA/modules

lib_LIBRARIES = libelsi.a 

LDADD = $(LIBOMMLIB) $(LIBMATRIXSWITCH) $(LIBSPARMETISLIB) $(SUPERLU_LIBS) $(HDF5_FLIBS) $(HDF5_LIBS)  $(LIBELSI)


thisincludedir = $(top_builddir)/external/PEXSI/external/include

#AM_CFLAGS = $(LIBOMM) -I$(LIBMATSWITCH)/. -I$(thisincludedir)
#AM_CFLAGS   = -I$(thisincludedir)
#AM_CXXFLAGS = -lstdc++ -std=c++11 -I$(thisincludedir) -I$(top_builddir)/external/PEXSI/include

AM_FCFLAGS =  -I$(thisincludedir) -I$(top_builddir)/external/PEXSI/include -I$(top_builddir)/external/PEXSI/src -lstdc++ $(HDF5_FFLAGS) -I$(MATRIXSWITCH) -I$(ELPAMODULES)

libelsi_a_SOURCES = elsi_dimensions.F90 elsi_timings.f90 elsi_matrix_conversion.f90 \
                    elsi_hdf5_tools.F90 elsi_mpi_tools.f90 elsi_interface.f90

