# NOTE: This Makefile does NOT support auto-dependency for the .h files.
# If the header files are changed, do "make clean" first.

include ../make.inc

SRCS_F90  = f_driver_pselinv_real.f90 f_driver_pselinv_complex.f90 \
						f_driver_pselinv_real_unsym.f90 f_driver_pselinv_complex_unsym.f90 \
						f_driver_ksdft.f90 f_interface.f90

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

all: f_driver_ksdft f_driver_pselinv_real f_driver_pselinv_complex

f_driver_ksdft: f_interface.o f_driver_ksdft.o ${PEXSI_LIB} 
	($(FC) -o $@_${SUFFIX} f_driver_ksdft.o  $(FLOADOPTS) )

f_driver_ksdft.o: f_driver_ksdft.f90
	  ${FC} -c ${FFLAGS} $<


f_driver_pselinv_real: f_interface.o f_driver_pselinv_real.o ${PEXSI_LIB} 
	($(FC) -o $@_${SUFFIX} f_driver_pselinv_real.o  $(FLOADOPTS) )

f_driver_pselinv_real.o: f_driver_pselinv_real.f90
	  ${FC} -c ${FFLAGS} $<


f_driver_pselinv_complex: f_interface.o f_driver_pselinv_complex.o ${PEXSI_LIB} 
	($(FC) -o $@_${SUFFIX} f_driver_pselinv_complex.o  $(FLOADOPTS) )

f_driver_pselinv_complex.o: f_driver_pselinv_complex.f90
	  ${FC} -c ${FFLAGS} $<

f_driver_pselinv_real_unsym: f_interface.o f_driver_pselinv_real_unsym.o ${PEXSI_LIB} 
	($(FC) -o $@_${SUFFIX} f_driver_pselinv_real_unsym.o  $(FLOADOPTS) )

f_driver_pselinv_real_unsym.o: f_driver_pselinv_real_unsym.f90
	  ${FC} -c ${FFLAGS} $<


f_driver_pselinv_complex_unsym: f_interface.o f_driver_pselinv_complex_unsym.o ${PEXSI_LIB} 
	($(FC) -o $@_${SUFFIX} f_driver_pselinv_complex_unsym.o  $(FLOADOPTS) )

f_driver_pselinv_complex_unsym.o: f_driver_pselinv_complex_unsym.f90
	  ${FC} -c ${FFLAGS} $<

f_interface.o: f_interface.f90
	  ${FC} -c ${FFLAGS} $<

f_ppexsi_old: f_ppexsi_old.o ${PEXSI_LIB} 
	  ($(FC) -o $@_${SUFFIX} f_ppexsi_old.o  $(FLOADOPTS) )

f_ppexsi_old.o: f_ppexsi_old.f90
	  ${FC} -c ${FFLAGS} $<

-include ${DEPS}

${PEXSI_LIB}:
	(cd ${PEXSI_DIR}/src; make all)

cleanall:
	rm -f ${EXES} ${OBJS} ${DEPS} *.d.* *.o *.mod
