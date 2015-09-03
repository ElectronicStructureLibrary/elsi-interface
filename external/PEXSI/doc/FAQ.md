Frequently asked questions {#pageFAQ}
==========================

General questions
-----------------

**Q: How do I know whether %PEXSI works for my application?**

A: %PEXSI may not necessarily be faster than the diagonalization method or
other competitive methods.  The simplest way to see whether %PEXSI
brings acceleration for your applications is to use %PEXSI to compute
the selected elements of the inverse for a typical matrix from your
applications.  See @ref pagePselinvRealSymmetric and @ref driver_pselinv_real.c for how to do this.

**Q: Can I just use the selected inversion routine?**

A: The parallel selected inversion (PSelInv) is a standalone routine.
See @ref pagePselinvRealSymmetric for an example.


**Q: Does %PEXSI accelerate dense matrix computation?**

A: No.  The acceleration is based on the sparsity of the LU factor or
the Cholesky factor.  %PEXSI should not be fast if the matrix is dense
or nearly dense.

**Q: Does %PEXSI work for asymmetric matrices?**

A: Currently %PEXSI only works for (real or complex) symmetric matrices.
We plan to release the version that supports asymmetric matrices in the
near future.  If you have applications in mind and %PEXSI may bring
significant acceleration to your application, please contact us as in
@ref pageTrouble, with some description of your matrix and application.

**Q: How to control the amount of output information and which processor
outputs the log file?**

See @ref pagePEXSIPlan. Also the `verbosity` option in @ref PPEXSIOptions
control the amount of information.

Installation
------------

**Q: The FORTRAN routine cannot compile.**

A: For FORTRAN users, `CPP_LIB=-lstdc++ -lmpi -lmpi_cxx` is often
needed.  Check this first if there is link error.


Performance
-----------

**Q: What if %PEXSI is numerically unstable or inaccurate for my
application?**

A: If you are using the pole expansion, the expansion converges
exponentially with respect to the number of poles.  So first increase
the number of poles until the error saturates.  After this step, error
only comes from the selected inversion phase. The selected inversion is
in principle an exact method, but may possibly suffer from numerical
instability issue due to the round off erros.  This is possible due to
the lack of dynamic pivoting strategies.  If this problem persists,
please contact us as in @ref pageTrouble, with some description of your
matrix and application.

**Q: When using ParMETIS/PT-Scotch, I got segmentaiton fault in the
factorization phase.**

A: We have observed that when ParMETIS/PT-Scotch is used, the number of
processors for symbolic factorization (`npSymbFact`) cannot be larger
than a magic number depending on the matrix and the machine.  This may
to be related to the parallel symbolic factorization routine in
SuperLU\_DIST.  If this problem happens, try to reduce `npSymbFact` to
a small number (such as 4 or 16), or even try to use the sequential
symbolic factorization if feasible.
