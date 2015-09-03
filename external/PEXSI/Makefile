include make.inc

all: lib examples 

lib: pexsi_lib

examples: pexsi_examples

pexsi_lib:
	cd src && ${MAKE}

examples: pexsi_examples
  
pexsi_examples: pexsi_lib
	cd examples && ${MAKE}

clean:
	cd src && ${MAKE} clean

cleanall:
	cd src && ${MAKE} cleanall
	cd examples && ${MAKE} cleanall
