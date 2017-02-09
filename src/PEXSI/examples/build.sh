#! /bin/bash
cd ../src; make cleanall; make -j; cd ../examples; rm $1; rm $1.o; make $1 
