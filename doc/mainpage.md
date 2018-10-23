Welcome to the documentation of ELSI Interface           {#mainpage}
==============================================

This is the source code documentation of the ELSI interface software. For more information, please visit the <a href="http://elsi-interchange.org">ELSI Interchange</a> website, or contact us via email: elsi-team@duke.edu

ELSI provides and enhances scalable, open-source software library solutions for electronic structure theory for simulations in materials science, condensed matter physics, chemistry, molecular biochemistry, and many other fields. ELSI focuses particularly on methods that solve or circumvent eigenvalue problems in the self-consistent field cycle of density-functional theory, the key bottleneck “cubic scaling wall” of electronic structure theory. In addition, the same infrastructure will help overcome other challenging eigenvalue problems as well.

We place particular emphasis on efficient support across platforms, from laptop-type computers all the way to the most efficient available massively parallel computers and new architectures (GPU and manycore processors).

The libraries supported and enhanced in ELSI are:

* <a href="http://elpa.mpcdf.mpg.de">ELPA</a> (massively parallel dense eigensolvers)

* <a href="http://william-dawson.github.io/NTPoly">libOMM</a> (orbital minimization method)

* <a href="http://esl.cecam.org/LibOMM">libOMM</a> (density matrix purification)

* <a href="http://pexsi.org">PEXSI</a> (pole expansion and selected inversion)

* <a href="http://keceli.github.io/SLEPc-SIPs">SLEPc-SIPs</a> (sparse eigensolver based on shift-and-invert spectral transformations)

We expect that related strong libraries will be supported by ELSI in the future, providing electronic structure codes and users with a flexible, customizable choice of solutions for the central algebraic problems that are often the bottleneck of our simulations.
