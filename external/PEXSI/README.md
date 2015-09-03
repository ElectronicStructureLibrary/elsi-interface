PEXSI: Pole EXpansion and Selected Inversion 
============================================

For full documentation, including instructions for installation please see

http://www.pexsi.org

Installation       
============



Build ParMETIS
--------------

Download ParMETIS (latest version 4.0.2) from

http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.2.tar.gz

Follow the installation step to install ParMETIS.

After untar the ParMETIS package, in Install.txt

    Edit the file metis/include/metis.h and specify the width (32 or
    64 bits) of the elementary data type used in ParMetis (and
    METIS). This is controled by the IDXTYPEWIDTH constant.

    For now, on a 32 bit architecture you can only specify a width
    of 32, whereas for a 64 bit architecture you can specify a width
    of either 32 or 64 bits.

In our experience for most cases, the following setup work
fine.

    #define IDXTYPEWIDTH 32


Build SuperLU_DIST
------------------

Download SuperLU_DIST (latest version 3.3) from

http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_dist_3.3.tar.gz

Follow the installation step to install SuperLU_DIST.

Our experience shows that on some machines it may be better
to build SuperLU_DIST with -O2 option than the more aggresive
optimization options provided by vendors.

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
at NERSC).  `${PEXSI_DIR}` stands for the main directory of PEXSI.

    cd ${PEXSI_DIR}
    cp config/make.inc.edison.intel make.inc

Edit the variables in make.inc. 
    
    PEXSI_DIR     = Main directory for PEXSI
    DSUPERLU_DIR  = Main directory for SuperLU_DIST
    METIS_DIR     = Main directory for METIS
    PARMETIS_DIR  = Main directory for ParMETIS 
    PTSCOTCH_DIR  = Main directory for PT-Scotch


PEXSI can be compiled using `debug` or `release` mode in
by the variable `COMPILE_MODE` in `make.inc`.  This variable mainly controls the
compiling flag `-DRELEASE`.  The `debug` mode introduces tracing of call
stacks at all levels of functions, and may significantly slow down the
code.  For production runs, use `release` mode.

The `*.profile` configuration files are for debugging purpose and
can be ignored.

Build the PEXSI library
------------------------

If make.inc is configured correctly,
    
    cd ${PEXSI_DIR}
    cd src
    make

should produce `libpexsi_(suffix).a` under `src/`.

Build examples
--------------

After `libpexsi_(suffix).a` is built, all driver routines are readily to be
compiled.  For example, the selected inversion for a complex matrix has
the test routine

    cd ${PEXSI_DIR}
    cd examples
    make driver_pselinv_complex

should produce `driver_pselinv_complex`, which can be executed with MPI.

For more information on the examples, see @ref pageTutorial.

Tests
-----

After driver_pselinv_complex is compiled, 

    examples$ mpirun -n 1 ./driver_pselinv_complex

should return the diagonal of the matrix
\f[
  (A + i I)^{-1}
\f]
saved on the 0-th processor, where \f$A\f$ is the five-point
discretization of a Laplacian operator on a 2D domain.  The result can
be compared with `examples/driver_pselinv_complex.out` to check the
correctness of the result.
