-include make.inc

all: lib 

lib: pexsi_lib

pexsi_lib:
	cd src && ${MAKE}

#examples: pexsi_examples
  
#pexsi_examples: pexsi_lib
#	cd examples && ${MAKE}

clean:
	cd src && ${MAKE} cleanall
	cd examples && ${MAKE} cleanall
