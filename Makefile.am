# This must match definition of AC_CONFIG_MACRO_DIR in configure.ac
ACLOCAL_AMFLAGS = -I config

SUBDIRS = external ELSI

.PHONY: external libparmetis libsuperlu pexsi_external external MatrixSwitch libOMM elpa elsi

.NOTPARALLEL: external pexsi_external

elsi-example: libelsi
	@echo Start building ELSI examples
	$(MAKE) -C ELSI/example $(AM_MAKEFLAGS)
	@echo ELSI examples installed.

elsi: libOMM libpexsi elpa
	@echo Start building ELSI
	$(MAKE) -C ELSI/src  $(AM_MAKEFLAGS)
	@echo ELSI installed.

#cd external/ELPA && ${MAKE} $(AM_MAKEFLAGS) # ; ${MAKE} install
elpa:
	@echo Start building ELPA...
	 ${MAKE} $(AM_MAKEFLAGS) -C external/ELPA
	@echo ELPA installed.

libOMM: MatrixSwitch
	$(MAKE) -C external/libOMM/src

MatrixSwitch:
	$(MAKE) -C external/libOMM/src/MatrixSwitch-0.1.2/src

pexsi: external
	$(MAKE) -C external/PEXSI/src 
	cp external/PEXSI/src/f_ppexsi_interface.mod  external/PEXSI/include

external: pexsi_external

pexsi_external: libparmetis libsuperlu

if !HAS_SUPERLU
libsuperlu: 
	$(MAKE) -C external/SuperLU SuperLU
else
libsuperlu:
endif

if !HAS_LIBPARMETIS
libparmetis: 
	$(MAKE) -C external/parmetis parMetis
else
libparmetis:
endif

