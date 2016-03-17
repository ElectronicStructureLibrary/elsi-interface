


#this has to change

LIBPEXSI=$(top_builddir)/external/PEXSI/src/libpexsi.a

LIBOMMLIB=$(top_builddir)/external/libOMM/src/libOMM.a

if HAVE_OPENMP
LIBELPA=$(top_builddir)/external/ELPA/.libs/libelpa_openmp.a
else
LIBELPA=$(top_builddir)/external/ELPA/.libs/libelpa.a
endif

LIBELSI=$(top_builddir)/ELSI/src/libelsi.a

ELSIMODULES = $(top_builddir)/ELSI/src


#LIBOMM = $(LIBOMM_FLAGS: -D%=-WF,-D%)

MATRIXSWITCH = $(top_builddir)/external/libOMM/src/MatrixSwitch-0.1.2/src
LIBMATRIXSWITCH = $(MATRIXSWITCH)/libmatrixswitch.a
ELPAMODULES = $(top_builddir)/external/ELPA/modules


LDADD = $(LIBOMMLIB) $(LIBMATRIXSWITCH) $(SCALAPACK_FCFLAGS)

thisincludedir = $(top_builddir)/external/include

AM_FCFLAGS =  -I$(thisincludedir) $(HDF5_FFLAGS) -I$(MATRIXSWITCH)

noinst_PROGRAMS = example1.x example2.x example3.x example4.x example5.x example6.x example7.x

example1_x_SOURCES = example1.F90

example2_x_SOURCES = example2.F90

example3_x_SOURCES = example3.F90 

example4_x_SOURCES = example4.F90

example5_x_SOURCES = example5.F90

example6_x_SOURCES = example6.F90

example7_x_SOURCES = example7.F90 


