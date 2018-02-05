###############################
# TOP-LEVEL MAKEFILE FOR ELSI #
###############################

RELEASE_DATE = 2017-05-27

# A "make.sys" file is mandatory
include make.sys

# Pre-compiled external libraries
ifneq ($(strip $(EXTERNAL_ELPA)),yes)
  ALL_OBJ   += elpa
  CLEAN_OBJ += cleanelpa
  ELPA_OMM   = elpa
  ELSI_ELPA  = elsi_elpa.o
else
  ELSI_ELPA = elsi_aeo.o
endif

ifneq ($(strip $(EXTERNAL_OMM)),yes)
  ALL_OBJ   += omm
  CLEAN_OBJ += cleanomm
endif

ifeq ($(strip $(DISABLE_CXX)),yes)
  DISABLE_PEXSI = yes
endif

ifneq ($(strip $(DISABLE_PEXSI)),yes)
  ifneq ($(strip $(EXTERNAL_PEXSI)),yes)
    ALL_OBJ   += pexsi
    CLEAN_OBJ += cleanpexsi
  endif
else
    $(info ====================================)
    $(info = PEXSI disabled by DISABLE_PEXSI. =)
    $(info ====================================)
endif

ifeq ($(strip $(ENABLE_SIPS)),yes)
  ifneq ($(strip $(EXTERNAL_SIPS)),yes)
    ALL_OBJ   += sips
    CLEAN_OBJ += cleansips
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
ELSI_DIR   = $(THIS_DIR)/src/ELSI
LIB_DIR    = $(THIS_DIR)/lib
INC_DIR    = $(THIS_DIR)/include
BUILD_DIR  = $(THIS_DIR)/build

# Default external libraries
ELPA_DIR  ?= $(THIS_DIR)/src/ELPA
ELPA_LIB  ?= -L$(LIB_DIR) -lelpa
OMM_DIR   ?= $(THIS_DIR)/src/libOMM
OMM_LIB   ?= -L$(LIB_DIR) -lOMM -lMatrixSwitch
PEXSI_DIR ?= $(THIS_DIR)/src/PEXSI
PEXSI_LIB ?= -L$(LIB_DIR) -lpexsi
PEXSI_LIB += $(SUPERLU_LIB) $(ORDERING_LIB)
SIPS_DIR  ?= $(THIS_DIR)/src/SIPs
SIPS_LIB  ?= -L$(LIB_DIR) -lqetsc
SIPS_LIB  += $(SLEPC_LIB) $(PETSC_LIB)
SIPS_INC  += $(SLEPC_INC) $(PETSC_INC)

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

# Create C interfaces
C_INTERFACE ?= no
ifeq ($(strip $(C_INTERFACE)),yes)
  C_BINDING += test_dm_complex_c.x
  C_BINDING += test_dm_real_c.x
  C_BINDING += test_ev_complex_c.x
  C_BINDING += test_ev_real_c.x
  C_BINDING += test_standard_ev_real_c.x
endif

export

LIBS = $(OMM_LIB) $(ELPA_LIB)
INCS = $(OMM_INC) $(ELPA_INC) -I$(INC_DIR)

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

.PHONY: all elpa omm pexsi sips elsi install check checkc clean cleanelsi cleanelpa cleanomm cleanpexsi cleansips

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

sips:
	@echo ==========================
	@echo = Start building SIPs... =
	@echo ==========================
	mkdir -p $(INC_DIR)
	mkdir -p $(LIB_DIR)
	cd $(SIPS_DIR) && $(MAKE) -f Makefile.elsi install
	@echo ===================
	@echo = SIPs installed. =
	@echo ===================

elsi: $(ALL_OBJ)
	@echo ==========================
	@echo = Start building ELSI... =
	@echo ==========================
	mkdir -p $(INC_DIR)
	mkdir -p $(LIB_DIR)
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

cleansips:
	@echo ====================
	@echo = Removing SIPs... =
	@echo ====================
	rm -f $(SIPS_DIR)/*.o
	rm -f $(SIPS_DIR)/*.mod
	rm -f $(SIPS_DIR)/*.a
