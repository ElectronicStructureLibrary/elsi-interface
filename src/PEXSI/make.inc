CC      = $(MPICC)
CXX     = $(MPICXX)
FC      = $(MPIFC)
LOADER  = $(MPICXX)

CP      = cp
RM      = rm
RMFLAGS = -f

PEXSI_INCLUDE = -I$(PEXSI_DIR)/include
INCLUDES      = $(PEXSI_INCLUDE) $(SUPERLU_INC)

CPP_LIB       = $(CXX_LDFLAGS)
PEXSI_LIB     = $(PEXSI_DIR)/src/libpexsi.a

COMPILE_DEF  = -DRELEASE

PAR_ND_LIB = $(PARMETIS_LIB)
SEQ_ND_LIB = $(METIS_LIB)

LIBS = $(PEXSI_LIB) $(SUPERLU_LIB) $(PARMETIS_LIB) $(METIS_LIB) $(SCALAPACK_LIB)

COMPILE_DEF += -DAdd_

CFLAGS    = $(CFLAGS_P) $(INCLUDES)
FFLAGS    = $(FFLAGS_P) $(INCLUDES)
CXXFLAGS  = $(CXXFLAGS_P) $(CPPFLAG) $(INCLUDES)
CCDEFS    = $(COMPILE_DEF)
CPPDEFS   = $(COMPILE_DEF)

# Generate auto-dependencies
%.d: %.c
	@set -e; rm -f $@; \
	$(CC) -M $(CCDEFS) $(CFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CXX) -M $(CPPDEFS) $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$
