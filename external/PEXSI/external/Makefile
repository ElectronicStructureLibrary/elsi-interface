CC = $(MPICC)

all: parMetis SuperLU

parMetis:
	@echo Start building parMetis...
	mkdir -p parMetis
	cp fetch-parmetis.sh $(PEXSI_DIR)/external/parMetis/fetch.sh
	cd $(PEXSI_DIR)/external/parMetis && ./fetch.sh
	cd $(PEXSI_DIR)/external/parMetis && make config prefix=$(PEXSI_DIR)/external
	cd $(PEXSI_DIR)/external/parMetis/build/Linux-x86_64 && make && make install
	mkdir -p include
	mkdir -p lib
	-cp `find $(PEXSI_DIR)/external/parMetis/include -name parmetis.h` $(PEXSI_DIR)/external/include/.
	-cp `find $(PEXSI_DIR)/external/parMetis/metis/include -name metis.h` $(PEXSI_DIR)/external/include/.
	-cp `find $(PEXSI_DIR)/external/parMetis/build -name libparmetis.a` $(PEXSI_DIR)/external/lib/.
	-cp `find $(PEXSI_DIR)/external/parMetis/build -name libmetis.a` $(PEXSI_DIR)/external/lib/.
	@echo parMetis installed.

SuperLU:
	@echo Start building SuperLU...
	mkdir -p SuperLU
	cp fetch-superlu.sh $(PEXSI_DIR)/external/SuperLU/fetch.sh
	cd $(PEXSI_DIR)/external/SuperLU && ./fetch.sh
	cp make.inc.superlu SuperLU/make.inc
	cd $(PEXSI_DIR)/external/SuperLU && make
	-cp $(PEXSI_DIR)/external/SuperLU/SRC/*.h $(PEXSI_DIR)/external/include/.
	-cp $(PEXSI_DIR)/external/SuperLU/lib/libsuperlu_dist_3.3.a $(PEXSI_DIR)/external/lib/.
	@echo SuperLU installed.

clean:
	rm -rf parMetis SuperLU bin include lib
