


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


LDADD = $(LIBELSI) $(LIBELPA) $(LIBOMMLIB) $(LIBMATRIXSWITCH) $(LIBPEXSI) $(LIBSUPERLU) $(LIBPARMETIS) $(LIBMETIS) $(HDF5_FLIBS) $(HDF5_LIBS)  $(SCALAPACK_FCFLAGS)

thisincludedir = $(top_builddir)/external/include

AM_FCFLAGS =  -I$(thisincludedir) -I$(top_builddir)/external/PEXSI/include $(HDF5_FFLAGS) -I$(MATRIXSWITCH) -I$(ELPAMODULES) -I$(ELSIMODULES) -I$(top_builddir)/external/PEXSI/src

noinst_PROGRAMS = construct_and_write_elpa.x \
                  set_and_convert_pexsi.x \
                  set_and_write_elpa.x \
                  read_and_write_elpa.x \
                  read_and_solve_elpa.x \
                  construct_and_write_omm.x \
                  set_and_write_omm.x \
                  read_and_write_omm.x \
                  read_and_write_pexsi.x \
                  read_and_solve_omm.x \
                  read_and_solve_pexsi.x


construct_and_write_elpa_x_SOURCES = construct_and_write_elpa.f90

set_and_convert_pexsi_x_SOURCES = set_and_convert_pexsi.f90

set_and_write_elpa_x_SOURCES = set_and_write_elpa.f90

read_and_write_elpa_x_SOURCES = read_and_write_elpa.f90

read_and_solve_elpa_x_SOURCES = read_and_solve_elpa.f90

construct_and_write_omm_x_SOURCES = construct_and_write_omm.f90

set_and_write_omm_x_SOURCES = set_and_write_omm.f90

read_and_write_omm_x_SOURCES = read_and_write_omm.f90

read_and_write_pexsi_x_SOURCES = read_and_write_pexsi.f90

read_and_solve_omm_x_SOURCES = read_and_solve_omm.f90

read_and_solve_pexsi_x_SOURCES = read_and_solve_pexsi.f90


