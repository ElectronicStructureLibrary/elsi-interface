include make.sys
export

all: elpa libomm pexsi elsi

elpa:
	@echo Start building ELPA...
	cd $(ELPA_DIR) && ${MAKE} -f Makefile.elpa; ${MAKE} -f Makefile.elpa install
	@echo ELPA installed.
libomm:
	@echo Start building libOMM...
	cd $(libOMM_DIR)/src && ${MAKE}
	cd $(libOMM_DIR)/examples && ${MAKE}
	@echo libOMM installed.

pexsi:
	@echo Start building PEXSI...
	cd $(PEXSI_DIR)/src && ${MAKE}
	-cp $(PEXSI_DIR)/src/f_ppexsi_interface.mod $(PEXSI_DIR)/include
	cd $(PEXSI_DIR)/fortran && ${MAKE} 
	@echo PEXSI installed.

elsi:
	@echo Start building ELSI...
	cd ELSI && ${MAKE}
	@echo ELSI installed.

external: pexsi_external

pexsi_external:
	cd $(PEXSI_DIR)/external && ${MAKE}

install:
	cd ELSI && ${MAKE} install
	@echo ELSI library created.

check:
	cd ELSI && ${MAKE} check

clean:
	cd ELSI && ${MAKE} clean
	cd $(ELPA_DIR) && ${MAKE} -f Makefile.elpa clean
	cd $(libOMM_DIR)/src && ${MAKE} clean
	cd $(libOMM_DIR)/src/MatrixSwitch-0.1.2/src && ${MAKE} clean
	cd $(libOMM_DIR)/examples && ${MAKE} clean
	cd $(PEXSI_DIR)/src && ${MAKE} cleanall
	cd $(PEXSI_DIR)/fortran && ${MAKE} cleanall
	rm -f $(PEXSI_DIR)/include/f_ppexsi_interface.mod
	rm -f $(PEXSI_DIR)/src/f_ppexsi_interface.mod

cleanall: clean
	cd $(PEXSI_DIR)/external && ${MAKE} clean
