CC  = $(MPICC)
CXX = $(MPICXX)
FC  = $(MPIFC)

ifeq ($(strip $(PEXSI_ADD_UNDERSCORE)),no)
  COMPILE_DEF = -DRELEASE
else
  COMPILE_DEF = -DRELEASE -DAdd_
endif

CFLAGS   = $(CFLAGS_P) -I$(PEXSI_DIR)/include $(SUPERLU_INC)
FFLAGS   = $(FFLAGS_P) -I$(PEXSI_DIR)/include $(SUPERLU_INC)
CXXFLAGS = $(CXXFLAGS_P) -I$(PEXSI_DIR)/include $(SUPERLU_INC)
CCDEFS   = $(COMPILE_DEF)
CPPDEFS  = $(COMPILE_DEF)

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
