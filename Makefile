-include make.sys
export

ELPA_HOME = external/ELPA

OMM_HOME = external/libOMM

#PEXSI_HOME = external/PEXSI

all: elpa libomm interface

elpa:
	@echo Start building ELPA...
	cd external/ELPA && ${MAKE} -f Makefile.elpa; ${MAKE} -f Makefile.elpa install
	@echo ELPA installed.
libomm:
	@echo Start building libOMM...
	cd external/libOMM/src && ${MAKE}
	cd external/libOMM/examples && ${MAKE}
	@echo libOMM installed.

pexsi:
	@echo Start building PEXSI...
#	cd external/PEXSI && ${MAKE}
	@echo PEXSI installed.

interface:
	@echo Start building ELSI-interface...
	cd ELSI && ${MAKE}
	@echo ELSI-interface installed.

install:
	cd ELSI && ${MAKE} install
	@echo ELSI library created.

check:
	cd ELSI && ${MAKE} check

clean:
	cd ELSI && ${MAKE} clean
	cd external/ELPA && ${MAKE} -f Makefile.elpa clean
	cd external/libOMM/src && ${MAKE} clean
	cd external/libOMM/src/MatrixSwitch-0.1.2/src && ${MAKE} clean
	cd external/libOMM/examples && ${MAKE} clean
#	cd external/PEXSI && ${MAKE} cleanall
