
SUBDIRS = MatrixSwitch-0.1.2

MATSWITCH=$(top_builddir)/external/libOMM/src/MatrixSwitch-0.1.2/src
LIBMATSWITCH=$(MATSWITCH)/libmatrixswitch.a


AM_FCFLAGS = $(LIBOMM) -I$(MATSWITCH)/.


lib_LIBRARIES = libOMM.a

LDADD = $(LIBMATSWITCH)

libOMM_a_SOURCES = omm_params.F90 omm_rand.F90 omm_quartic.F90 omm_ops.F90 omm.F90 omm_wrappers.F90

