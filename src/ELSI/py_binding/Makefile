F2PY     = f2py-2.7
MPIFC    = mpif90-mpich-mp
FFLAGS   = -O3 -ffree-line-length-none 
ELSI_DIR = ../../..

#METIS_LIB    = -L/home/wph6/opt/ParMETIS/4.0.3/build/Linux-x86_64/libmetis -lmetis
#PARMETIS_LIB = -L/home/wph6/opt/ParMETIS/4.0.3/build/Linux-x86_64/libparmetis -lparmetis
#SUPERLU_LIB  = -L/home/wph6/opt/SuperLU_DIST/5.1.3/lib -lsuperlu_dist
#SUPERLU_INC  = -I/home/wph6/opt/SuperLU_DIST/5.1.3/SRC

SCALAPACK_LIB = -L/opt/local/lib -lscalapack -llapack -lblas 

.PHONY: all clean test 

all:  elsi_python.so tomato_python.so elsi_scalapack_python.so

elsi_python.so: elsi_python.f90
	$(F2PY) --f90exec=$(MPIFC) --f90flags="$(FFLAGS)" -m elsi_python -c elsi_python.f90 -I$(ELSI_DIR)/include $(SUPERLU_INC) $(SCALAPACK_INC) -L$(ELSI_DIR)/lib -lelsi -lelpa -lOMM -lMatrixSwitch -lpspblas -ltomato $(SCALAPACK_LIB)

tomato_python.so: tomato_python.f90
	$(F2PY) --f90exec=$(MPIFC) --f90flags="$(FFLAGS)" -m tomato_python -c tomato_python.f90 -I$(ELSI_DIR)/include -L$(ELSI_DIR)/lib -lMatrixSwitch -lpspblas -ltomato $(SCALAPACK_LIB)

elsi_scalapack_python.so: elsi_scalapack_python.f90
	$(F2PY) --f90exec=$(MPIFC) --f90flags="$(FFLAGS)" -m elsi_scalapack_python -c elsi_scalapack_python.f90 -I$(ELSI_DIR)/include -L$(ELSI_DIR)/lib -lMatrixSwitch -lpspblas $(SCALAPACK_LIB)

clean:
	rm -f elsi_python.so tomato_python.so elsi_scalapack_python.so

test:  elsi_python.so
	python -c "import elsi_python"
