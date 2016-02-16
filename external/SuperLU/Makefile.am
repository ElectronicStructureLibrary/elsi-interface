# This must match definition of AC_CONFIG_MACRO_DIR in configure.ac
#ACLOCAL_AMFLAGS = -I config

PEXSI_DIR = $(top_builddir)/external/PEXSI

EXTDIR = $(top_builddir)/external

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


SUPERLUDIR= $(EXTDIR)/SuperLU/SuperLU
.PHONY: clean SuperLU

all: SuperLU

if !HAS_SUPERLU
SuperLU:
	@echo Start building SuperLU... 
	-cp $(EXTDIR)/SuperLU/make.inc.superlu $(SUPERLUDIR)/make.inc
	$(MAKE) $(AM_MAKEFLAGS) -C $(SUPERLUDIR) MPICC="$(MPICC)" MPIFC="$(MPIF77)" PEXSI_DIR="$(abs_srcdir)/../.." FLIBS="$(LIBS) $(FCLIBS)" AR=$(AR) RANLIB=$(RANLIB) THIS_DIR="$(abs_srcdir)/../.."
	-cp $(SUPERLUDIR)/SRC/*.h $(EXTDIR)/include/.
	@echo SuperLU installed.

else
SuperLU:
endif


#	-cp $(SUPERLUDIR)/lib/libsuperlu.a $(EXTDIR)/lib/.

clean:
	$(MAKE) -C $(SUPERLUDIR) clean
	rm -f $(EXTDIR)/lib/libsuperlu.a

