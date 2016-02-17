# This must match definition of AC_CONFIG_MACRO_DIR in configure.ac
#ACLOCAL_AMFLAGS = -I config

PEXSI_DIR = $(top_builddir)/external/PEXSI
MPIFC = $(MPIF77)

PARMETISLIB  = -L$(PEXSI_DIR)/external/lib -lparmetis

.PHONY: clean

all: parMetis SuperLU

SuperLU: parMetis
	@echo Start building SuperLU... $(MPIFC) $(PEXSI_DIR)/external/SuperLU
	mkdir -p $(PEXSI_DIR)/external/SuperLU
	cp fetch-superlu.sh $(PEXSI_DIR)/external/SuperLU/fetch.sh
	($(am__cd) $(PEXSI_DIR)/external/SuperLU && ./fetch.sh) || exit 1;
	cp make.inc.superlu $(PEXSI_DIR)/external/SuperLU/make.inc
	$(am__cd) $(PEXSI_DIR)/external/SuperLU && $(MAKE) $(AM_MAKEFLAGS) MPICC="$(MPICC)" MPIFC="$(MPIF77)" PEXSI_DIR="$(abs_srcdir)/.." FLIBS="$(FCLIBS)" AR=$(AR) RANLIB=$(RANLIB)
	-cp $(PEXSI_DIR)/external/SuperLU/SRC/*.h $(PEXSI_DIR)/external/include/.
	-cp $(PEXSI_DIR)/external/SuperLU/lib/libsuperlu.a $(PEXSI_DIR)/external/lib/.
	-cp $(PEXSI_DIR)/external/lib/libsuperlu.a $(PEXSI_DIR)/lib/.
	@echo SuperLU installed.

parMetis:
	@echo Start building parMetis...
	mkdir -p parMetis
	cp fetch-parmetis.sh $(PEXSI_DIR)/external/parMetis/fetch.sh
	($(am__cd)  $(PEXSI_DIR)/external/parMetis && ./fetch.sh ) || exit 1;
	($(am__cd) $(PEXSI_DIR)/external/parMetis && $(MAKE) $(AM_MAKEFLAGS) config prefix=$(PEXSI_DIR)/external) || exit 1;
	($(am__cd) $(PEXSI_DIR)/external/parMetis/build/Linux-x86_64 && make && $(MAKE) $(AM_MAKEFLAGS) install) || exit 1;
	mkdir -p include
	mkdir -p lib
	-cp `find $(PEXSI_DIR)/external/parMetis/include -name parmetis.h` $(PEXSI_DIR)/external/include/.
	-cp `find $(PEXSI_DIR)/external/parMetis/metis/include -name metis.h` $(PEXSI_DIR)/external/include/.
	-cp `find $(PEXSI_DIR)/external/parMetis/build -name libparmetis.a` $(PEXSI_DIR)/external/lib/.
	-cp `find $(PEXSI_DIR)/external/parMetis/build -name libmetis.a` $(PEXSI_DIR)/external/lib/.
	($(am__cd) $(PEXSI_DIR) && mkdir -p lib ) || exit 1;
	-cp $(PEXSI_DIR)/external/lib/libparmetis.a $(PEXSI_DIR)/lib/.
	-cp $(PEXSI_DIR)/external/lib/libmetis.a $(PEXSI_DIR)/lib/.
	@echo parMetis installed.


clean-local:
	rm -rf parMetis SuperLU bin include lib
	rm -f ../lib/libmetis.a ../lib/libparmetis.a ../lib/libsuperlu.a

