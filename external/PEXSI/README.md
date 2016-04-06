PEXSI: Pole EXpansion and Selected Inversion 
============================================

For full documentation, including instructions for installation please see

http://www.pexsi.org

Installation       
============



Build ParMETIS
--------------

Download ParMETIS (latest version 4.0.3) from

http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz

Follow the installation step to install ParMETIS.

@attention After untar the ParMETIS package, in Install.txt

    Edit the file metis/include/metis.h and specify the width (32 or
    64 bits) of the elementary data type used in ParMetis (and
    METIS). This is controled by the IDXTYPEWIDTH constant.

    For now, on a 32 bit architecture you can only specify a width
    of 32, whereas for a 64 bit architecture you can specify a width
    of either 32 or 64 bits.

@attention In our experience for most cases, the following setup work
fine.

    #define IDXTYPEWIDTH 32



Build SuperLU_DIST
------------------

**SuperLU_DIST v4.3 starting from PEXSI v0.9.2**


Download SuperLU_DIST (latest version 4.3) from

http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_dist_4.3.tar.gz

Follow the installation step to install SuperLU_DIST.

@attention Our experience shows that on some machines it may be better
to build SuperLU_DIST with -O2 option than the more aggresive
optimization options provided by vendors.

@attention In SuperLU_DIST v4.3, some functions conflict when both real
and complex arithmetic factorization is needed. This can be temporarily
solved by adding  `-Wl,--allow-multiple-definition` in the linking
option.

@attention In SuperLU_DIST v4.3, there could be some excessive outputs.
This can be removed by going to the SRC/ directory of superlu, and
comment out the line starting with `printf(".. dQuery_Space` in
dmemory_dist.c. Do the same thing for the line starting with
`printf(".. zQuery_Space` in zmemory_dist.c.

@attention Please note that the number of processors for symbolic
factorization cannot be too large when PARMETIS is used together with
SuperLU. The exact number of processors for symbolic factorization is
unfortunately a **magic parameter**.  See @ref pageFAQ.

**SuperLU_DIST v3.3 for PEXSI v0.9.0 and before**

Download SuperLU_DIST (latest version 3.3) from

http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_dist_3.3.tar.gz


(Optional) Build PT-Scotch
--------------------------

On some machines, ParMETIS may only allow to use a relatively small
number of processors for the matrix permutation. In such circumstance, a
workaround can be to use PT-Scotch, which can be downloaded from
(latest version 6.0.0)

https://gforge.inria.fr/frs/download.php/31831/scotch_6.0.0.tar.gz

Follow the installation step to install PT-Scotch.

In INSTALL.TXT, pay special attention to the following
sections in order to compile PT-Scotch correctly.

    2.3) Integer size issues
    2.5) Threads issues


PT-Scotch is also METIS-Compatible.  See the following section in
INSTALL.TXT for more information.

    2.9) MeTiS compatibility library

Edit make.inc
-------------

Configuration of PEXSI is controlled by a single `make.inc` file.
Examples of the `make.inc` file are given under the `config/` directory.

Find `make.inc` with the most similar architecture, and copy to the main
PEXSI directory (using Edison for example, the latest Intel computer
at NERSC, a CRAY X30 machine).  `${PEXSI_DIR}` stands for the main
directory of PEXSI.

    cd ${PEXSI_DIR}
    cp config/make.inc.CRAY_XC30.intel make.inc

Edit the variables in make.inc. 
    
    PEXSI_DIR     = Main directory for PEXSI
    DSUPERLU_DIR  = Main directory for SuperLU_DIST
    PARMETIS_DIR  = Main directory for ParMETIS 
    PTSCOTCH_DIR  = Main directory for PT-Scotch

Edit the compiler options, for instance

    CC           = cc
    CXX          = CC
    FC           = ftn
    LOADER       = CC

@note 

- Starting from PEXSI v0.8.0, `-std=c++11` is required in
`CXXFLAGS`. 

- Starting from PEXSI v0.9.2, `-std=c99` is required to be compatible
  with SuperLU_DIST v4.3.

- For **FORTRAN** users, `CPP_LIB=-lstdc++ -lmpi -lmpi_cxx` is often needed.
Check this if there is link error.

- PEXSI can be compiled using `debug` or `release` mode in
by the variable `COMPILE_MODE` in `make.inc`.  This variable mainly controls the
compiling flag `-DRELEASE`.  The `debug` mode introduces tracing of call
stacks at all levels of functions, and may significantly slow down the
code.  For production runs, use `release` mode.

- The `USE_PROFILE` options is for internal test purpose. Usually
set this to 0.

Build the PEXSI library
------------------------

If make.inc is configured correctly,
    
    make 
    make install

Should build the PEXSI library under the `build` directory ready to be
used in an external package.  If
examples are needed (not necessary if you use PEXSI in an external
package), type 

    make all

which will generate C examples in `examples/` directory and FORTRAN examples in
`fortran/` directory, respectively.

For more information on the examples, see @ref pageTutorial.

Tests
-----

After example files are compiled, go to the `examples/` directory, and

    examples$ mpirun -n 1 ./driver_pselinv_complex_(suffix)

should return the diagonal of the matrix
\f[
  (A + i I)^{-1}
\f]
saved on the 0-th processor, where \f$A\f$ is the five-point
discretization of a Laplacian operator on a 2D domain.  The result can
be compared with `examples/driver_pselinv_complex.out` to check the
correctness of the result. For more examples see @ref pageTutorial.

