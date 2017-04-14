###############################
# TOP-LEVEL MAKEFILE FOR ELSI #
###############################

# A "make.sys" file is mandatory
include make.sys

# Default archive tools
ARCHIVE      ?= ar
ARCHIVEFLAGS ?= cr
RANLIB       ?= ranlib

# Default external libraries
THIS_DIR ?= $(shell pwd)

BIN_DIR   = $(THIS_DIR)/bin
LIB_DIR   = $(THIS_DIR)/lib
INC_DIR   = $(THIS_DIR)/include
BUILD_DIR = $(THIS_DIR)/build

ELSI_DIR  = $(THIS_DIR)/src/ELSI

ELPA_DIR  = $(THIS_DIR)/src/ELPA
ELPA_LIB  = -L$(LIB_DIR) -lelpa \
            -I$(INC_DIR)

OMM_DIR   = $(THIS_DIR)/src/libOMM
OMM_LIB   = -I$(INC_DIR) \
            -L$(LIB_DIR)/lib -lOMM -lMatrixSwitch -lpspblas -ltomato

PEXSI_DIR = $(THIS_DIR)/src/PEXSI
PEXSI_LIB = -I$(INC_DIR)/include \
            -L$(LIB_DIR)/lib -lpexsi \
            $(SUPERLU_LIB) $(PARMETIS_LIB) $(METIS_LIB)

# Compiler settings
FFLAGS_I   ?= $(FFLAGS) $(SCALAPACK_FCFLAGS)
FFLAGS_E   ?= $(FFLAGS) $(SCALAPACK_FCFLAGS)
CFLAGS_E   ?= $(CFLAGS) $(SCALAPACK_FCFLAGS)
FFLAGS_O   ?= $(FFLAGS) $(SCALAPACK_FCFLAGS)
FFLAGS_P   ?= $(FFLAGS) $(SCALAPACK_FCFLAGS)
CFLAGS_P   ?= $(CFLAGS) $(SCALAPACK_FCFLAGS)
CXXFLAGS_P ?= $(CXXFLAGS) $(SCALAPACK_FCFLAGS)
LINKER     ?= $(MPIFC)
LDFLAGS    ?= $(FFLAGS_I)

# Default architecture settings
ELPA2_KERNEL ?= Generic
OPENMP       ?= no

# Name of MPI executable
MPI_EXEC ?= mpirun

# Create C interfaces
C_INTERFACE ?= yes
ifeq ($(strip $(C_INTERFACE)),yes)
  C_BINDING = elsi_c_interface.o
endif

export

# Use stub if PEXSI is disabled in make.sys
ALL_OBJ   = elpa omm
CLEAN_OBJ = cleanelpa cleanomm
LIBS      = $(ELPA_LIB) $(OMM_LIB)

ifeq ($(strip $(DISABLE_CXX)),yes)
  $(info ===================================================)
  $(info = PEXSI disabled by user's choice of DISABLE_CXX. =)
  $(info ===================================================)
  STUBS += stub_pexsi.o
else
  ALL_OBJ += pexsi
  CLEAN_OBJ += cleanpexsi
  LIBS += $(PEXSI_LIB)
endif

ifeq ($(strip $(ENABLE_SIPS)),yes)
  LIBS += $(SIPS_LIB)
else
  STUBS += stub_sips.o
endif

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

omm: elpa
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
	cd $(ELSI_DIR) && $(MAKE) -f Makefile.elsi
	@echo ===============================
	@echo = ELSI compiled successfully. =
	@echo ===============================

install:
	rm -f $(ELSI_DIR)/*log*
	cd $(ELSI_DIR) && $(MAKE) -f Makefile.elsi install
	@echo ======================================
	@echo = ELSI library created successfully. =
	@echo ======================================

check:
	@echo ================================
	@echo = Running ELSI test programs.. =
	@echo ================================
	cd $(ELSI_DIR) && $(MAKE) -f Makefile.elsi check
	@echo ================================
	@echo = ELSI test programs finished. =
	@echo ================================

clean: $(CLEAN_OBJ) cleanelsi

cleanelsi:
	@echo ====================
	@echo = Removing ELSI... =
	@echo ====================
	cd $(ELSI_DIR) && $(MAKE) -f Makefile.elsi clean
	rm -rf $(INC_DIR) $(LIB_DIR) $(BIN_DIR) $(BUILD_DIR)
	@echo =========
	@echo = Done. =
	@echo =========

cleanelpa:
	@echo ====================
	@echo = Removing ELPA... =
	@echo ====================
	cd $(ELPA_DIR) && $(MAKE) -f Makefile.elsi clean

cleanomm:
	@echo ======================
	@echo = Removing libOMM... =
	@echo ======================
	cd $(OMM_DIR) && $(MAKE) -f Makefile.elsi clean

cleanpexsi:
	@echo =====================
	@echo = Removing PEXSI... =
	@echo =====================
	cd $(PEXSI_DIR)/src && $(MAKE) -f Makefile.elsi clean
