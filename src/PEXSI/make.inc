CC      = $(MPICC)
CXX     = $(MPICXX)
FC      = $(MPIFC)

PEXSI_INCLUDE = -I$(PEXSI_DIR)/include
INCLUDES      = $(PEXSI_INCLUDE) $(SUPERLU_INC)

PEXSI_LIB     = $(PEXSI_DIR)/src/libpexsi.a

COMPILE_DEF  = -DRELEASE -DAdd_

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
