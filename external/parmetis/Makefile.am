# This must match definition of AC_CONFIG_MACRO_DIR in configure.ac
#ACLOCAL_AMFLAGS = -I config

PEXSI_DIR = $(top_builddir)/external/PEXSI

EXTDIR = $(top_builddir)/external

MPIFC = $(MPIF77)

if !HAS_LIBPARMETIS
LIBPARMETIS=$(top_builddir)/external/lib/libparmetis.a
LIBMETIS=$(top_builddir)/external/lib/libmetis.a
else
LIBPARMETIS=$(LIBSPARMETISLIB)
LIBMETIS=$(LIBSPARMETISLIB)
endif

CPU=$(shell uname -m | sed "s/\\ /_/g")
SYS=$(shell uname -s)

if HAVE_OPENMP
OMP=1
endif

.PHONY: clean parMetis

#install: parMetis
	
all: parMetis

libraries: parMetis

if !HAS_LIBPARMETIS
parMetis:
	mkdir -p $(EXTDIR)/lib
	mkdir -p $(EXTDIR)/include
	@echo Start building parMetis...
	$(MAKE) -C $(EXTDIR)/parmetis/parMetis $(AM_MAKEFLAGS) config prefix=$(EXTDIR) cc=$(MPICC) cxx=$(MPICXX) cputype=$(CPU) systype=$(SYS) openmp=$(OMP)
	$(MAKE)  $(AM_MAKEFLAGS) -C $(EXTDIR)/parmetis/parMetis/build/$(SYS)-$(CPU)
	$(MAKE)  $(AM_MAKEFLAGS) -C $(EXTDIR)/parmetis/parMetis/build/$(SYS)-$(CPU)  install
	@echo parMetis installed.
	@echo Start installing parMetis...
	mkdir -p $(EXTDIR)/include
	mkdir -p $(EXTDIR)/lib
	-cp `find $(EXTDIR)/parmetis/parMetis/include -name parmetis.h` $(EXTDIR)/include/.
	-cp `find $(EXTDIR)/parmetis/parMetis/metis/include -name metis.h` $(EXTDIR)/include/.
	-cp `find $(EXTDIR)/parmetis/parMetis/build -name libparmetis.a` $(EXTDIR)/lib/.
	-cp `find $(EXTDIR)/parmetis/parMetis/build -name libmetis.a` $(EXTDIR)/lib/.
	@echo parMetis installed.
else
parMetis:
endif



clean:
	$(MAKE) -C $(EXTDIR)/parmetis/parMetis $(AM_MAKEFLAGS) clean
	rm -f ../include/* ../lib/*
	rm -f $(EXTDIR)/lib/libmetis.a $(EXTDIR)/lib/libparmetis.a $(EXTDIR)/lib/libsuperlu.a

