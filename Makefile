###############################
# TOP-LEVEL MAKEFILE FOR ELSI #
###############################

RELEASE_DATE = 2017-05-27

# A "make.sys" file is mandatory
include make.sys

# Pre-compiled external libraries
ifneq ($(strip $(EXTERNAL_ELPA)),yes)
  ALL_OBJ  += elpa
  INST_OBJ += installelpa
  ELPA_OMM  = elpa
  ELSI_ELPA = elsi_elpa.o
else
  ELSI_ELPA = elsi_aeo.o
  ELSI_ELPA_API = elsi_elpa_api.o
endif

ifneq ($(strip $(EXTERNAL_OMM)),yes)
  ALL_OBJ  += omm
  INST_OBJ += installomm
endif

ifeq ($(strip $(DISABLE_CXX)),yes)
  DISABLE_PEXSI = yes
endif

ifneq ($(strip $(DISABLE_PEXSI)),yes)
  ifneq ($(strip $(EXTERNAL_PEXSI)),yes)
    ALL_OBJ  += pexsi
    INST_OBJ += installpexsi
  endif
else
    $(info ====================================)
    $(info = PEXSI disabled by DISABLE_PEXSI. =)
    $(info ====================================)
endif

ifeq ($(strip $(ENABLE_SIPS)),yes)
  ifneq ($(strip $(EXTERNAL_SIPS)),yes)
    ALL_OBJ  += sips
    INST_OBJ += installsips
  endif
endif

ifneq ($(strip $(PTSCOTCH_LIB)),)
  ORDERING_LIB = $(PTSCOTCH_LIB)
else
  ORDERING_LIB = $(PARMETIS_LIB) $(METIS_LIB)
endif

# Default archive tools
ARCHIVE      ?= ar
ARCHIVEFLAGS ?= cr
RANLIB       ?= ranlib

# Default installation directories
THIS_DIR   = $(shell pwd)
PREFIX    ?= $(THIS_DIR)
LIB_DIR    = $(PREFIX)/lib
INC_DIR    = $(PREFIX)/include
BUILD_DIR ?= $(THIS_DIR)/build
ELSI_DIR   = $(THIS_DIR)/src/ELSI
ELSI_LIB   = -L$(BUILD_DIR)/ELSI -lelsi
ELSI_INC   = -I$(BUILD_DIR)/ELSI

# Default external libraries
ELPA_DIR  ?= $(THIS_DIR)/src/ELPA
ELPA_LIB  ?= -L$(BUILD_DIR)/ELPA -lelpa
ELPA_INC  ?= -I$(BUILD_DIR)/ELPA
OMM_DIR   ?= $(THIS_DIR)/src/libOMM
OMM_LIB   ?= -L$(BUILD_DIR)/OMM -lOMM -lMatrixSwitch
OMM_INC   ?= -I$(BUILD_DIR)/OMM
PEXSI_DIR ?= $(THIS_DIR)/src/PEXSI
PEXSI_LIB ?= -L$(BUILD_DIR)/PEXSI -lpexsi
PEXSI_LIB += $(SUPERLU_LIB) $(ORDERING_LIB)
PEXSI_INC ?= -I$(BUILD_DIR)/PEXSI
SIPS_DIR  ?= $(THIS_DIR)/src/SIPs
SIPS_LIB  ?= -L$(BUILD_DIR)/SIPS -lqetsc
SIPS_LIB  += $(SLEPC_LIB) $(PETSC_LIB)
SIPS_INC  += $(SLEPC_INC) $(PETSC_INC)
SIPS_INC  += -I$(BUILD_DIR)/SIPS

# Default compiler settings
FFLAGS_I   ?= $(FFLAGS)
FFLAGS_E   ?= $(FFLAGS)
CFLAGS_E   ?= $(CFLAGS)
FFLAGS_O   ?= $(FFLAGS)
FFLAGS_P   ?= $(FFLAGS)
CFLAGS_P   ?= $(CFLAGS)
CXXFLAGS_P ?= $(CXXFLAGS)
LINKER     ?= $(MPIFC)
FLINKER    ?= $(LINKER)
CLINKER    ?= $(LINKER)

# Default architecture settings
ELPA_OMP     ?= no
ELPA_GPU     ?= no
ELPA2_KERNEL ?= Generic

# MPI for "make check"
MPI_EXEC ?= mpirun
MPI_SIZE ?= 4

export

LIBS = $(ELSI_LIB) $(OMM_LIB) $(ELPA_LIB)
INCS = $(ELSI_INC) $(OMM_INC) $(ELPA_INC)

ifneq ($(strip $(DISABLE_PEXSI)),yes)
  LIBS += $(PEXSI_LIB)
  INCS += $(PEXSI_INC)
else
  STUBS += stub_pexsi.o
endif

ifeq ($(strip $(ENABLE_SIPS)),yes)
  LIBS += $(SIPS_LIB)
  INCS += $(SIPS_INC)
else
  STUBS += stub_sips.o
endif

ifeq ($(strip $(ENABLE_CHESS)),yes)
  LIBS += $(CHESS_LIB)
  INCS += $(CHESS_INC)
else
  STUBS += stub_chess.o
endif

.PHONY: all elpa omm pexsi sips elsi \
        create_dir install installelpa installomm installpexsi installsips \
        clean cleanelsi cleanelpa cleanomm cleanpexsi cleansips \
        test testc check checkc

all: $(ALL_OBJ) elsi

elpa:
	@echo ==========================
	@echo = Start building ELPA... =
	@echo ==========================
	mkdir -p $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/ELPA
	cd $(BUILD_DIR)/ELPA && $(MAKE) -f $(ELPA_DIR)/Makefile.elsi
	@echo ===================
	@echo = ELPA installed. =
	@echo ===================

omm: $(ELPA_OMM)
	@echo ============================
	@echo = Start building libOMM... =
	@echo ============================
	mkdir -p $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/OMM
	cd $(BUILD_DIR)/OMM && $(MAKE) -f $(OMM_DIR)/Makefile.elsi
	@echo =====================
	@echo = libOMM installed. =
	@echo =====================

pexsi:
	@echo ===========================
	@echo = Start building PEXSI... =
	@echo ===========================
	mkdir -p $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/PEXSI
	cd $(BUILD_DIR)/PEXSI && $(MAKE) -f $(PEXSI_DIR)/Makefile.elsi
	@echo ====================
	@echo = PEXSI installed. =
	@echo ====================

sips:
	@echo ==========================
	@echo = Start building SIPs... =
	@echo ==========================
	mkdir -p $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/SIPS
	cd $(BUILD_DIR)/SIPS && $(MAKE) -f $(SIPS_DIR)/Makefile.elsi
	@echo ===================
	@echo = SIPs installed. =
	@echo ===================

elsi: $(ALL_OBJ)
	@echo ==========================
	@echo = Start building ELSI... =
	@echo ==========================
	mkdir -p $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/ELSI
	cd $(BUILD_DIR)/ELSI && $(MAKE) -f $(ELSI_DIR)/Makefile.elsi
	@echo ===============================
	@echo = ELSI compiled successfully. =
	@echo ===============================

installelpa: create_dir
	cp $(BUILD_DIR)/ELPA/*.mod $(INC_DIR)
	cp $(BUILD_DIR)/ELPA/*.a $(LIB_DIR)

installomm: create_dir
	cp $(BUILD_DIR)/OMM/*.mod $(INC_DIR)
	cp $(BUILD_DIR)/OMM/*.a $(LIB_DIR)

installpexsi: create_dir
	cp $(BUILD_DIR)/PEXSI/*.mod $(INC_DIR)
	cp $(BUILD_DIR)/PEXSI/*.a $(LIB_DIR)

installsips: create_dir
	cp $(BUILD_DIR)/SIPS/*.mod $(INC_DIR)
	cp $(BUILD_DIR)/SIPS/*.a $(LIB_DIR)

install: create_dir $(INST_OBJ)
	cp $(BUILD_DIR)/ELSI/*.mod $(INC_DIR)
	cp $(BUILD_DIR)/ELSI/*.h $(INC_DIR)
	cp $(BUILD_DIR)/ELSI/*.a $(LIB_DIR)
	@echo ================================
	@echo = ELSI installed successfully. =
	@echo ================================

create_dir:
	mkdir -p $(INC_DIR)
	mkdir -p $(LIB_DIR)

test: elsi
	cd $(BUILD_DIR)/ELSI && $(MAKE) -f $(ELSI_DIR)/test/Makefile.elsi test
	@echo =========
	@echo = Done. =
	@echo =========

testc: elsi
	cd $(BUILD_DIR)/ELSI && $(MAKE) -f $(ELSI_DIR)/test/Makefile.elsi testc
	@echo =========
	@echo = Done. =
	@echo =========

check: test
	@echo ========================================
	@echo = Running ELSI Fortran test programs.. =
	@echo ========================================
	cd $(BUILD_DIR)/ELSI && $(MAKE) -f $(ELSI_DIR)/test/Makefile.elsi checkf
	@echo =========
	@echo = Done. =
	@echo =========

checkc: testc
	@echo ==================================
	@echo = Running ELSI C test programs.. =
	@echo ==================================
	cd $(BUILD_DIR)/ELSI && $(MAKE) -f $(ELSI_DIR)/test/Makefile.elsi checkc
	@echo =========
	@echo = Done. =
	@echo =========

clean:
	@echo ====================
	@echo = Removing ELSI... =
	@echo ====================
	rm -rf $(INC_DIR)
	rm -rf $(LIB_DIR)
	rm -rf $(BUILD_DIR)
	@echo =========
	@echo = Done. =
	@echo =========

cleanelsi:
	@echo ====================
	@echo = Removing ELSI... =
	@echo ====================
	rm -rf $(BUILD_DIR)/ELSI

cleanelpa:
	@echo ====================
	@echo = Removing ELPA... =
	@echo ====================
	rm -rf $(BUILD_DIR)/ELPA

cleanomm:
	@echo ======================
	@echo = Removing libOMM... =
	@echo ======================
	rm -rf $(BUILD_DIR)/OMM

cleanpexsi:
	@echo =====================
	@echo = Removing PEXSI... =
	@echo =====================
	rm -rf $(BUILD_DIR)/PEXSI

cleansips:
	@echo ====================
	@echo = Removing SIPs... =
	@echo ====================
	rm -rf $(BUILD_DIR)/SIPS
