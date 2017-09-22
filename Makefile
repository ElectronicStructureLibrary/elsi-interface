###############################
# TOP-LEVEL MAKEFILE FOR ELSI #
###############################

# A "make.sys" file is mandatory
include make.sys

# Pre-compiled external libraries
ifneq ($(strip $(EXTERNAL_ELPA)),yes)
  ALL_OBJ += elpa
  CLEAN_OBJ += cleanelpa
  ELPA_OMM = elpa
endif

ifneq ($(strip $(EXTERNAL_OMM)),yes)
  ALL_OBJ += omm
  CLEAN_OBJ += cleanomm
endif

ifeq ($(strip $(DISABLE_CXX)),yes)
  DISABLE_PEXSI = yes
endif

ifneq ($(strip $(EXTERNAL_PEXSI)),yes)
  ifneq ($(strip $(DISABLE_PEXSI)),yes)
    ALL_OBJ += pexsi
    CLEAN_OBJ += cleanpexsi
  else
    $(info ====================================)
    $(info = PEXSI disabled by DISABLE_PEXSI. =)
    $(info ====================================)
    STUBS += stub_pexsi.o
  endif
endif

# Default archive tools
ARCHIVE      ?= ar
ARCHIVEFLAGS ?= cr
RANLIB       ?= ranlib

# Default external libraries
THIS_DIR   = $(shell pwd)
ELSI_DIR   = $(THIS_DIR)/src/ELSI

BIN_DIR    = $(THIS_DIR)/bin
LIB_DIR    = $(THIS_DIR)/lib
INC_DIR    = $(THIS_DIR)/include
BUILD_DIR  = $(THIS_DIR)/build

ELPA_DIR  ?= $(THIS_DIR)/src/ELPA
ELPA_INC  ?= -I$(INC_DIR)
ELPA_LIB  ?= -L$(LIB_DIR) -lelpa

OMM_DIR   ?= $(THIS_DIR)/src/libOMM
OMM_INC   ?= -I$(INC_DIR)
OMM_LIB   ?= -L$(LIB_DIR) -lOMM -lMatrixSwitch

PEXSI_DIR ?= $(THIS_DIR)/src/PEXSI
PEXSI_INC ?= -I$(INC_DIR)
PEXSI_LIB ?= -L$(LIB_DIR) -lpexsi $(SUPERLU_LIB) $(PARMETIS_LIB) $(METIS_LIB)

# Default compiler settings
FFLAGS_I   ?= $(FFLAGS) $(SCALAPACK_INC)
FFLAGS_E   ?= $(FFLAGS) $(SCALAPACK_INC)
CFLAGS_E   ?= $(CFLAGS) $(SCALAPACK_INC)
FFLAGS_O   ?= $(FFLAGS) $(SCALAPACK_INC)
FFLAGS_P   ?= $(FFLAGS) $(SCALAPACK_INC)
CFLAGS_P   ?= $(CFLAGS) $(SCALAPACK_INC)
CXXFLAGS_P ?= $(CXXFLAGS) $(SCALAPACK_INC)
LINKER     ?= $(MPIFC)
LDFLAGS    ?= $(FFLAGS_I)

# Default architecture settings
ELPA_OMP      ?= no
ELPA_GPU      ?= no
ELPA2_KERNEL  ?= Generic

# Name of MPI executable
MPI_EXEC ?= mpirun

# Create C interfaces
C_INTERFACE ?= yes
ifeq ($(strip $(C_INTERFACE)),yes)
  C_BINDING += test_dm_complex_c.x
  C_BINDING += test_dm_real_c.x
  C_BINDING += test_ev_complex_c.x
  C_BINDING += test_ev_real_c.x
  C_BINDING += test_standard_ev_real_c.x
endif

export

LIBS = $(ELPA_LIB) $(OMM_LIB)
INCS = $(ELPA_INC) $(OMM_INC) -I$(INC_DIR)

ifneq ($(strip $(DISABLE_PEXSI)),yes)
  LIBS += $(PEXSI_LIB)
  INCS += $(PEXSI_INC)
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

.PHONY: all elpa omm pexsi elsi install check checkc clean cleanelsi cleanelpa cleanomm cleanpexsi

all: $(ALL_OBJ) elsi

elpa:
	@echo ==========================
	@echo = Start building ELPA... =
	@echo ==========================
	mkdir -p $(INC_DIR)
	mkdir -p $(LIB_DIR)
	cd $(ELPA_DIR) && $(MAKE) -f Makefile.elsi install
	@echo ===================
	@echo = ELPA installed. =
	@echo ===================

omm: $(ELPA_OMM)
	@echo ============================
	@echo = Start building libOMM... =
	@echo ============================
	mkdir -p $(INC_DIR)
	mkdir -p $(LIB_DIR)
	cd $(OMM_DIR) && $(MAKE) -f Makefile.elsi install
	@echo =====================
	@echo = libOMM installed. =
	@echo =====================

pexsi:
	@echo ===========================
	@echo = Start building PEXSI... =
	@echo ===========================
	mkdir -p $(INC_DIR)
	mkdir -p $(LIB_DIR)
	cd $(PEXSI_DIR)/src && $(MAKE) -f Makefile.elsi install
	@echo ====================
	@echo = PEXSI installed. =
	@echo ====================

elsi: $(ALL_OBJ)
	@echo ==========================
	@echo = Start building ELSI... =
	@echo ==========================
	mkdir -p $(INC_DIR)
	mkdir -p $(LIB_DIR)
	mkdir -p $(BIN_DIR)
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && $(MAKE) -f $(ELSI_DIR)/Makefile.elsi
	@echo ===============================
	@echo = ELSI compiled successfully. =
	@echo ===============================

install:
	cd $(BUILD_DIR) && $(MAKE) -f $(ELSI_DIR)/Makefile.elsi install
	@echo ================================
	@echo = ELSI installed successfully. =
	@echo ================================

check:
	@echo ========================================
	@echo = Running ELSI Fortran test programs.. =
	@echo ========================================
	cd $(BUILD_DIR) && $(MAKE) -f $(ELSI_DIR)/Makefile.elsi check
	@echo ========================================
	@echo = ELSI Fortran test programs finished. =
	@echo ========================================

checkc:
	@echo ==================================
	@echo = Running ELSI C test programs.. =
	@echo ==================================
	cd $(BUILD_DIR) && $(MAKE) -f $(ELSI_DIR)/Makefile.elsi checkc
	@echo ==================================
	@echo = ELSI C test programs finished. =
	@echo ==================================

clean: $(CLEAN_OBJ) cleanelsi

cleanelsi:
	@echo ====================
	@echo = Removing ELSI... =
	@echo ====================
	rm -rf $(INC_DIR)
	rm -rf $(LIB_DIR)
	rm -rf $(BIN_DIR)
	rm -rf $(BUILD_DIR)
	@echo =========
	@echo = Done. =
	@echo =========

cleanelpa:
	@echo ====================
	@echo = Removing ELPA... =
	@echo ====================
	rm -f $(ELPA_DIR)/*.o
	rm -f $(ELPA_DIR)/*.mod
	rm -f $(ELPA_DIR)/*.a

cleanomm:
	@echo ======================
	@echo = Removing libOMM... =
	@echo ======================
	rm -f $(OMM_DIR)/*.o
	rm -f $(OMM_DIR)/*.mod
	rm -f $(OMM_DIR)/*.a

cleanpexsi:
	@echo =====================
	@echo = Removing PEXSI... =
	@echo =====================
	rm -f $(PEXSI_DIR)/src/*.o
	rm -f $(PEXSI_DIR)/src/*.d
	rm -f $(PEXSI_DIR)/src/*.d.*
	rm -f $(PEXSI_DIR)/src/*.mod
	rm -f $(PEXSI_DIR)/src/*.a
