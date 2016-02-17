
LIBOMM = $(LIBOMM_FLAGS: -D%=-WF,-D%)

#AM_CPPFLAGS = @LIBOMM_FLAGS@
AM_FCFLAGS = $(LIBOMM)

#$(DEFS:-D%=-WF,-D%) 

#bin_PROGRAMS = toy
#
#toy_SOURCES = toy.F90

#noinst_LIBRARIES = libmatrixswitch.a
lib_LIBRARIES = libmatrixswitch.a

libmatrixswitch_a_SOURCES=MatrixSwitch.F90
