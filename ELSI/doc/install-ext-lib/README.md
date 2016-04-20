
Parmetis, metis and SuperLU installation helper
===============================================
Author: v4m4

The two scripts in this folder:
fetch-parmetis.sh
fetch-superlu.sh

are meant to help to install external libraries 
requiered to compile PEXSI. They should work with 
Bash-shell.

You should modify these three variables to make them work.
EXTDIR=/tmp/install
MPICC=mpicc
MPICXX=mpicxx


I
You have also to define in make.inc.superlu the location 
where BLAS subroutines are installed.


For example:
BLASLIB      = -L${MKLROOT}/lib/intel64 \
               -I${MKLROOT}/include/intel64/lp64 -lmkl_intel_lp64 \
               -lmkl_sequential -lmkl_core -lpthread

