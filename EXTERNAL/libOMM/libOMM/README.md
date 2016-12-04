libOMM
======

Author
------

Fabiano Corsetti, Imperial College London, UK

*   Email: fabiano.corsetti08 {at} imperial.ac.uk
*   Homepage: <http://www.cmth.ph.ic.ac.uk/people/f.corsetti/>

Description
-----------

libOMM solves the Kohn-Sham equation as a generalized eigenvalue problem for a fixed Hamiltonian. It implements the orbital minimization method (OMM), which works within a density matrix formalism. The basic strategy of the OMM is to find the set of Wannier functions (WFs) describing the occupied subspace by direct unconstrained minimization of an appropriately-constructed functional. The density matrix can then be calculated from the WFs. The solver is usually employed within an outer self-consistency (SCF) cycle. Therefore, the WFs resulting from one SCF iteration can be saved and then re-used as the initial guess for the next iteration.

Installation
------------

1.  Enter the `src` directory.
2.  Copy `make.inc.example` to `make.inc` and modify it to suit your needs. `MSLIBPATH` should point to the MatrixSwitch directory (default in `make.inc.example` is for the version included in the distribution). MatrixSwitch should be compiled with the `-DCONV` flag. Available options for `FPPFLAGS` are:
    * `-DMPI`: enable MPI parallel routines
    * `-DLAP`: enable LAPACK routines (currently necessary for preconditioning/Cholesky factorization)
    * `-DSLAP`: enable ScaLAPACK routines (requires `-DMPI`)
    * `-DNORAND`: fixed seed for the random number generator. Enable for testing purposes.
    * `-DCBIND`: use ISO_C_BINDING for LOGICAL inputs in the wrapper interfaces. Enable for linking to C.
3.  Type `make`.

Testing
-------

The `examples` directory contains a number of small programs that make use of libOMM with MatrixSwitch. These can be useful both for testing the installation and for learning how to use the library. To compile them:

1.   Enter the `examples` directory.
2.   Copy `make.inc.example` to `make.inc` and modify it to suit your needs. Be aware that `make.inc` in the `src` directory will also be used.
3.   Type `make`.

Each example contains a header explaining what the program does and providing sample output to compare against.

Publications
------------

The algorithms and implementation are described in: F. Corsetti, *The orbital minimization method for electronic structure calculations with finite-range atomic basis sets*, Comput. Phys. Commun. **185**, 873 (2014). <http://dx.doi.org/10.1016/j.cpc.2013.12.008>

Documentation
-------------

A complete documentation is maintained at: <http://esl.cecam.org/libOMM>. Also see the examples in the `examples` directory.

License
-------

Copyright &copy; 2014, Fabiano Corsetti

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1.   Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2.   Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
