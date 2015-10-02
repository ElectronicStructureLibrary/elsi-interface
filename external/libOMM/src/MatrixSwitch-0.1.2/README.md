MatrixSwitch
============

Version
-------

0.1.2

Date
----

02.06.2015

Author
------

Fabiano Corsetti, CIC nanoGUNE, Donostia-San Sebastián, Spain

*   Email: f.corsetti {at} nanogune.eu
*   Homepage: <http://www.nanogune.eu/fabiano-corsetti/>

Description
-----------

MatrixSwitch is a module which acts as an intermediary interface layer between high-level routines for physics-related algorithms and low-level routines dealing with matrix storage and manipulation. This allows the high-level routines to be written in a way which is physically transparent, and enables them to switch seamlessly between different software implementations of the matrix operations.

Installation
------------

1.  Enter the `src` directory.
2.  Copy `make.inc.example` to `make.inc` and modify it to suit your needs. Available options for `FPPFLAGS` are:
    * `-DMPI`: enable MPI parallel routines
    * `-DLAP`: enable LAPACK routines
    * `-DSLAP`: enable ScaLAPACK routines (requires MPI)
    * `-DCONV`: enable automatic conversion of scalar types (real/complex) to agree with matrix definitions (real/complex). Note that conversions from complex to real will simply discard the imaginary part.
3.  Type `make`.

Testing
-------

The `examples` directory contains a number of small programs that make use of MatrixSwitch. These can be useful both for testing the installation and for learning how to use the library. To compile them:

1.   Enter the `examples` directory.
2.   Copy `make.inc.example` to `make.inc` and modify it to suit your needs. Be aware that `make.inc` in the `src` directory will also be used.
3.   Type `make`.

Each example contains a header explaining what the program does and providing sample output to compare against.

Usage
-----

`MatrixSwitch` is a module that you can `use` in Fortran routines. Note that both the `.a` and `.mod` files need to be available. An example compilation command for a code using MatrixSwitch is: `gfortran MyCode.f90 /path/to/MatrixSwitch-x.y.z/src/MatrixSwitch.a -I/path/to/MatrixSwitch-x.y.z/src/ -llapack -lblas`

Documentation
-------------

A complete documentation is maintained at: <http://esl.cecam.org/mediawiki/index.php/MatrixSwitch>. Also see the examples in the `examples` directory.

License
-------

Copyright &copy; 2014, Fabiano Corsetti

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1.   Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2.   Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
