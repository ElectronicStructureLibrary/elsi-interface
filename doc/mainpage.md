Welcome to the API documentation of ELSI {#mainpage}
========================================

ELSI provides and enhances scalable, open-source software library solutions for electronic structure calculations in materials science, condensed matter physics, chemistry, and many other fields. ELSI focuses on methods that solve or circumvent eigenproblems in electronic structure theory. For more information, please visit the [ELSI Interchange](https://elsi-interchange.org) website, or contact us via [email](mailto:elsi-team@duke.edu).

The eigensolvers and density matrix solvers supported in ELSI include:

* Shared-memory eigensolvers

    * [LAPACK](https://www.netlib.org/lapack): One-stage and two-stage tridiagonalization-based dense eigensolvers.

    * [MAGMA](https://icl.utk.edu/magma): GPU-accelerated one-stage and two-stage tridiagonalization-based dense eigensolvers.

* Distributed-memory eigensolvers

    * [ELPA](https://elpa.mpcdf.mpg.de): One-stage and two-stage tridiagonalization-based dense eigensolvers.

    * [EigenEXA](https://www.r-ccs.riken.jp/labs/lpnctrt/en/projects/eigenexa): One-stage tridiagonalization-based and pentadiagonalization-based dense eigensolvers.

    * [SLEPc](https://slepc.upv.es): Sparse eigensolver based on parallel spectrum slicing.

* Distributed-memory density matrix solvers

    * [libOMM](https://esl.cecam.org/LibOMM): Orbital minimization method based on dense linear algebra.

    * [PEXSI](https://pexsi.org): Pole expansion and selected inversion based on sparse linear algebra.

    * [NTPoly](https://william-dawson.github.io/NTPoly): Linear scaling density matrix purification based on sparse linear algebra.

The design of ELSI focues on portability from laptop-type computers all the way up to the most efficient massively parallel supercomputers and new architectures. Work is in progress to support additional solver libraries, providing electronic structure code developers and users with a flexible, customizable choice of solution for the central algebraic problems in large-scale electronic structure simulations.

The ELSI software is adopted by electronic structure code projects such as [DFTB+](https://www.dftbplus.org), DGDFT, [FHI-aims](https://aimsclub.fhi-berlin.mpg.de), and [SIESTA](https://departments.icmab.es/leem/siesta).
