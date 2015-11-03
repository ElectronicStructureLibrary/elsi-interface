#    This file is part of ELPA.
#
#    The ELPA library was originally created by the ELPA consortium, 
#    consisting of the following organizations:
#
#    - Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG), 
#    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
#      Informatik,
#    - Technische Universität München, Lehrstuhl für Informatik mit
#      Schwerpunkt Wissenschaftliches Rechnen , 
#    - Fritz-Haber-Institut, Berlin, Abt. Theorie, 
#    - Max-Plack-Institut für Mathematik in den Naturwissenschaftrn, 
#      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition, 
#      and  
#    - IBM Deutschland GmbH
#
#
#    More information can be found here:
#    http://elpa.rzg.mpg.de/
#
#    ELPA is free software: you can redistribute it and/or modify
#    it under the terms of the version 3 of the license of the 
#    GNU Lesser General Public License as published by the Free 
#    Software Foundation.
#
#    ELPA is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
#
#    ELPA reflects a substantial effort on the part of the original
#    ELPA consortium, and we ask you to respect the spirit of the
#    license that we chose: i.e., please contribute any changes you
#    may have back to the original ELPA library distribution, and keep
#    any derivatives of ELPA under the same license that we chose for
#    the original distribution, the GNU Lesser General Public License.

# ------------------------------------------------------------------------------
# Please set the variables below according to your system!
# ------------------------------------------------------------------------------

# It is strongly advised, not to use this Makefile, but build ELPA with the 
# configure autotools. If you want to do it anyway, here are some hints


# Settings for Intel Fortran (Linux) with OpenMP and GENERIC_KERNELS:
#

MPICC = mpicxx
FLAGS =  -O3 -funsafe-loop-optimizations \
         -funsafe-math-optimizations -ftree-vect-loop-version \
         -ftree-vectorize
# BLAS
BLAS = -L${MKLROOT}/lib/intel64 -lmkl_blas95_lp64 \
       -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
       -lpthread -lm
# LAPACK
LAPACK = -L${MKLROOT}/lib/intel64 -lmkl_lapack95_lp64
# SCALAPACK
SCALAPACK = -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 \
            -lmkl_blacs_intelmpi_lp64
#ELPA
ELPA_DIR = /home/bl126/devel/ELSI/external/ELPA
ELPA = -L$(ELPA_DIR)/lib -lelpa 

#FORTRAN RUNTIME Environment
FLIB = -L/opt/intel/lib/intel64 -lifcore -limf -lsvml

#LIBS
LIBS = ${ELPA} ${SCALAPACK} ${LAPACK} ${BLAS} ${FLIB}


all: test_real_c_version

test_real_c_version:
	$(MPICC) $(FLAGS) -o test_real_c_version test_real_c_version.c ${LIBS}

clean:
	rm -f *.o test_real_c_version
