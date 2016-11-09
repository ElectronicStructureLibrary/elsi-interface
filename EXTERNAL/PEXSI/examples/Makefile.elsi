# NOTE: This Makefile does NOT support auto-dependency for the .h files.
# If the header files are changed, do "make clean" first.

include ../make.inc

SRCS_CPP  = run_pselinv.cpp run_pselinv_unsym.cpp \
 						run_ksdft.cpp run_fermi.cpp run_fermi_complex.cpp
SRCS_C    = driver_pselinv_real.c driver_pselinv_complex.c driver_pselinv_real_unsym.c \
					 	driver_pselinv_complex_unsym.c driver_ksdft.c driver_ksdft_2.c \
						driver_fermi_complex.c

all: driver_pselinv_real driver_pselinv_real_unsym \
	driver_pselinv_complex driver_pselinv_complex_unsym \
	driver_ksdft driver_ksdft_2 driver_fermi_complex


ifeq (${USE_SYMPACK}, 1)
run_pselinv_sympack: run_pselinv_sympack.o ${PEXSI_LIB} ../include/pexsi/sympack_interf_impl.hpp 
	($(LOADER) -o $@_${SUFFIX} run_pselinv_sympack.o  $(LOADOPTS) )

endif

# Below are interfaces open to all users
driver_ksdft_2: driver_ksdft_2.o ${PEXSI_LIB}
	($(LOADER) -o $@_${SUFFIX} driver_ksdft_2.o  $(LOADOPTS) )

driver_ksdft: driver_ksdft.o ${PEXSI_LIB}
	($(LOADER) -o $@_${SUFFIX} driver_ksdft.o  $(LOADOPTS) )

driver_pselinv_real: driver_pselinv_real.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} driver_pselinv_real.o  $(LOADOPTS) )

driver_pselinv_complex: driver_pselinv_complex.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} driver_pselinv_complex.o  $(LOADOPTS) )

driver_pselinv_real_unsym: driver_pselinv_real_unsym.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} driver_pselinv_real_unsym.o  $(LOADOPTS) )

driver_pselinv_complex_unsym: driver_pselinv_complex_unsym.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} driver_pselinv_complex_unsym.o  $(LOADOPTS) )

driver_fermi_complex: driver_fermi_complex.o ${PEXSI_LIB}
	($(LOADER) -o $@_${SUFFIX} driver_fermi_complex.o  $(LOADOPTS) )


# Below are routines mainly for developers
run_fermi: run_fermi.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} run_fermi.o  $(LOADOPTS) )

run_fermi_complex: run_fermi_complex.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} run_fermi_complex.o  $(LOADOPTS) )

run_ksdft: run_ksdft.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} run_ksdft.o  $(LOADOPTS) )

run_inertia: run_inertia.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} run_inertia.o  $(LOADOPTS) )

run_pselinv: run_pselinv.o ${PEXSI_LIB} ../include/pexsi/*.hpp 
	($(LOADER) -o $@_${SUFFIX} run_pselinv.o  $(LOADOPTS) )

run_pselinv_unsym: run_pselinv_unsym.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} run_pselinv_unsym.o  $(LOADOPTS) )


run_test_reduce: run_test_reduce.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} run_test_reduce.o  $(LOADOPTS) )

run_test_reduce.o: run_test_reduce.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

run_inertia.o: run_inertia.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

run_pselinv.o: run_pselinv.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

run_pselinv_unsym.o: run_pselinv_unsym.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

run_ppexsi: run_ppexsi.o ${PEXSI_LIB} 
	($(LOADER) -o $@ run_ppexsi.o  $(LOADOPTS) )

run_ppexsi_old: run_ppexsi_old.o ${PEXSI_LIB} 
	($(LOADER) -o $@ run_ppexsi_old.o  $(LOADOPTS) )

run_ppexsi_old.o: run_ppexsi_old.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

run_superlu.o: run_superlu.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 
run_superlu: run_superlu.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} run_superlu.o  $(LOADOPTS) )

OBJS = ${SRCS_CPP:.cpp=.o} ${SRCS_C:.c=.o} 
DEPS = ${SRCS_CPP:.cpp=.d} ${SRCS_C:.c=.d}
EXES = ${SRCS_CPP:.cpp=} ${SRCS_C:.c=} 

# Compilation replacement rules

%.o: %.c
	${CC} -c ${CFLAGS} ${CCDEFS} $< 
%.o: %.cpp
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 
%.o: %.f
	${FC} -c ${FFLAGS} $<
%.o: %.F
	${FC} -c ${FFLAGS} $<


-include ${DEPS}

${PEXSI_LIB}:
	(cd ${PEXSI_DIR}/src; make all)

cleanall:
	rm -f ${EXES} ${OBJS} ${DEPS} *.d.* *.o *.mod
